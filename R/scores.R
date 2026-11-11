#' Evaluate the frozen maximum standardized score at candidate target values
#'
#' @param fit A fitted lifeconformal object, or its stored prediction state.
#' @param newdata Baseline data with one row per subject.
#' @param values Numeric target values in their original occupation-time units.
#'   Supply an n-by-d matrix/data frame with named target columns. Columns are
#'   reordered to the training target order. A numeric vector is accepted for
#'   one target, or a named target vector for one subject. A one-target matrix
#'   need not have a column name. Values must be finite; this score evaluator
#'   does not solve the separate joint-support membership problem.
#' @param initial_state Optional common or subject-specific initial-state labels.
#' @returns A numeric vector with one maximum absolute standardized deviation
#'   per baseline row. Subject order is preserved, including zero or one row.
#' @details Score scales include the training-frozen variance floor and relative
#'   target widths. No calibration histories, target refitting or model training
#'   are needed. Candidate values are not interpreted as already standardized.
#' @noRd
evaluate_life_score <- function(fit, newdata, values, initial_state = NULL) {
  moments <- predict_life_moments(fit, newdata, initial_state)
  n <- nrow(moments$mean)
  targets <- colnames(fit$problem$rewards)
  values <- align_life_score_values(values, n, targets)
  score <- numeric(n)
  for (j in seq_along(targets)) {
    score <- pmax(score, abs((values[, j] - moments$mean[, j]) / moments$scale[, j]))
  }
  unname(score)
}

#' Inspect the baseline outcome-model distribution of the frozen score
#'
#' @inheritParams evaluate_life_score
#' @param q Finite, nonnegative, strictly increasing score thresholds.
#' @returns A list with n-by-length(q) probability, the unchanged ordered q,
#'   numerical metadata, and lower/upper numerical probability bounds for
#'   deterministic Markov calculations. Simulation has no such deterministic
#'   bounds, so lower and upper are NULL in that branch. The result describes
#'   a working outcome distribution, not a calibrated or guaranteed coverage.
#' @details Always uses fit$outcome, independently of any separate augmentation
#'   fit. Monotone binary Markov rewards use a single native batch; other Markov
#'   rewards use deterministic joint-grid propagation with the fitted accuracy
#'   and resource controls. There is no Markov simulation fallback.
#'
#' Flexible models generate one joint target bank per profile, retaining all
#'   fitted history features. The same bank supplies frozen baseline moments
#'   and all requested score CDF values. This matches the stored profile-keyed
#'   moment evaluator and avoids duplicate simulations. Batches of at most 64
#'   profiles bound working target-bank storage independently of sample size.
#'   More thresholds do not trigger additional continuation draws.
#' @noRd
predict_life_coverage <- function(fit, newdata, q, initial_state = NULL) {
  if (!is.numeric(q) || is.complex(q) || !is.null(dim(q)) || !length(q) ||
      anyNA(q) || any(!is.finite(q)) || any(q < 0) || any(diff(q) <= 0)) {
    stop("q must be finite, nonnegative and strictly increasing.", call. = FALSE)
  }
  q <- unname(q)
  if (is.matrix(newdata)) newdata <- preprocessing_data_frame(newdata)
  if (!is.data.frame(newdata)) {
    stop("newdata must be a baseline data frame, one row per subject.", call. = FALSE)
  }
  initial <- prediction_initial_states(fit, newdata, initial_state)
  view <- apply_preprocessor(newdata, fit$preprocessor, "outcome")
  n <- nrow(newdata)
  d <- ncol(fit$problem$rewards)
  if (!n) {
    probability <- matrix(0, n, length(q))
    return(list(probability = probability, q = q, lower = probability,
                upper = probability,
                numerical = list(method = "empty_prediction",
                  working_model = fit$outcome$type, profiles = 0L,
                  maximum_probability_bound_width = 0,
                  deterministic_probability_bounds = TRUE,
                  approximation = "not_evaluated", coalesced = FALSE)))
  }
  if (fit$outcome$type != "markov_cox") {
    return(predict_flexible_score_coverage(fit, view$x, initial, q))
  }
  moments <- predict_markov_cox_moments(fit$outcome, view$x,
                                        fit$problem$rewards, initial)
  moments <- exact_deterministic_moments(moments, fit, initial)
  scale <- life_score_scale(moments$variance, fit)
  jumps <- predict_transition_jumps(fit$outcome, view$x)
  accrued <- matrix(0, n, d)
  orientation <- orient_markov_rewards(fit$problem$rewards, jumps$edges)
  if (!is.null(orientation)) {
    transformed <- orient_markov_query(orientation, moments$mean, accrued,
                                        rep(0, n), fit$problem$tau)
    result <- markov_monotone_batch(
      jumps$jump_times, jumps$base_increments, jumps$edges, jumps$relative_risk,
      transformed$rewards, fit$problem$tau, transformed$center, scale, q,
      subject_index = seq_len(n), query_times = rep(0, n), query_state = initial,
      accrued = transformed$accrued, weight = rep(1, n), flip = transformed$flip,
      original_center = transformed$original_center,
      original_accrued = transformed$original_accrued
    )
    result$diagnostics$working_model <- fit$outcome$type
    result$diagnostics$complemented_targets <- names(transformed$flip)[transformed$flip]
    result$diagnostics$maximum_probability_bound_width <- 0
    result$diagnostics$deterministic_probability_bounds <- TRUE
    result$diagnostics$approximation <- "floating_point_only"
    result$diagnostics$coalesced <- FALSE
    return(list(probability = result$contributions, q = q,
                lower = result$contributions, upper = result$contributions,
                numerical = result$diagnostics))
  }
  profiles <- life_exact_profile_groups(cbind(jumps$relative_risk, initial,
                                              moments$mean, scale))
  probability <- lower <- upper <- matrix(0, length(profiles$first), length(q))
  diagnostics <- vector("list", length(profiles$first))
  for (g in seq_along(profiles$first)) {
    i <- profiles$first[g]
    result <- with_life_grid_context(markov_conditional_coverage_grid(
      jumps$jump_times, NULL, fit$problem$rewards, fit$problem$tau,
      query_times = 0, query_state = initial[i], accrued = accrued[i, , drop = FALSE],
      center = moments$mean[i, , drop = FALSE], scale = scale[i, , drop = FALSE], q = q,
      side = "right", grid_size = fit$control$grid_size,
      max_nodes = fit$control$max_nodes,
      tolerance = fit$control$probability_tolerance, strict = TRUE,
      kernel_factors = life_markov_kernel_factors(jumps, i, fit$problem$absorbing)),
      profile_index = i, stage = "working_cdf",
      subject_id = if ("id" %in% names(newdata)) newdata$id[i] else NULL)
    probability[g, ] <- result$probability[1L, ]
    lower[g, ] <- result$lower[1L, ]
    upper[g, ] <- result$upper[1L, ]
    diagnostics[[g]] <- result$diagnostics
  }
  coalesced <- any(vapply(diagnostics, function(x) x$coalesced_state_layers > 0,
                          logical(1)))
  if (length(profiles$first) < n) {
    probability <- probability[profiles$group, , drop = FALSE]
    lower <- lower[profiles$group, , drop = FALSE]
    upper <- upper[profiles$group, , drop = FALSE]
    diagnostics <- diagnostics[profiles$group]
  }
  list(probability = probability, q = q, lower = lower, upper = upper,
       numerical = list(method = "deterministic_joint_grid_envelopes",
         working_model = fit$outcome$type, profiles = n,
         evaluated_profiles = length(profiles$first),
         reused_profiles = n - length(profiles$first),
         maximum_probability_bound_width = max(upper - lower),
         deterministic_probability_bounds = TRUE,
         approximation = if (coalesced) "bounded_joint_grid" else "floating_point_only",
         coalesced = coalesced,
         tolerance_met = TRUE, profile_diagnostics = diagnostics))
}

#' Align target values without silent recycling or positional target guessing
#' @noRd
align_life_score_values <- function(values, n, targets) {
  d <- length(targets)
  if (is.numeric(values) && !is.complex(values) && is.null(dim(values))) {
    if (d == 1L) {
      if (length(values) != n) {
        stop("values must have one value per baseline row for a scalar target.",
             call. = FALSE)
      }
      values <- matrix(values, n, 1L, dimnames = list(NULL, targets))
    } else if (n == 1L && length(values) == d && !is.null(names(values))) {
      values <- matrix(values, 1L, d, dimnames = list(NULL, names(values)))
    } else {
      stop("For multiple targets supply an n-by-d matrix with target column names.",
           call. = FALSE)
    }
  }
  if (is.data.frame(values)) values <- as.matrix(values)
  if (!is.matrix(values) || !is.numeric(values) || is.complex(values) ||
      !identical(dim(values), c(n, d)) || anyNA(values) || any(!is.finite(values))) {
    stop("values must be a finite numeric n-by-d target matrix matching baseline rows.",
         call. = FALSE)
  }
  labels <- colnames(values)
  if (is.null(labels) && d == 1L) {
    colnames(values) <- targets
  } else if (is.null(labels) || anyNA(labels) || anyDuplicated(labels) ||
             !setequal(labels, targets)) {
    stop("values must name every trained target exactly once: ",
         paste(targets, collapse = ", "), ".", call. = FALSE)
  }
  values[, targets, drop = FALSE]
}

#' Apply the frozen variance floor and named relative score widths
#' @noRd
life_score_scale <- function(variance, fit) {
  if (any(variance > 0 & variance < .Machine$double.xmin)) {
    stop("Score variance underflows at the current time unit; rescale time before fitting.",
         call. = FALSE)
  }
  # Apply the floor to standard deviations: its square can underflow even
  # when a deterministic target's requested positive width is representable.
  scale <- pmax(sqrt(variance), fit$control$scale_floor * fit$problem$tau)
  if (!is.null(fit$problem$score_weights)) {
    scale <- sweep(scale, 2L,
                    fit$problem$score_weights[colnames(fit$problem$rewards)], `*`)
  }
  if (any(!is.finite(scale)) || any(scale <= 0)) {
    stop("Score scales overflow or underflow; check the time unit and relative widths.",
         call. = FALSE)
  }
  scale
}

#' Reuse baseline flexible target banks for moments and all score thresholds
#' @noRd
predict_flexible_score_coverage <- function(fit, x, initial, q) {
  adapter <- fit$outcome
  n <- nrow(x)
  d <- ncol(fit$problem$rewards)
  model_d <- ncol(adapter$problem$rewards)
  selected <- match(colnames(fit$problem$rewards), colnames(adapter$problem$rewards))
  if (anyNA(selected)) {
    stop("The fitted outcome model lacks a requested target.", call. = FALSE)
  }
  B <- fit$control$draws
  probability <- matrix(0, n, length(q))
  for (start in seq.int(1L, n, by = 64L)) {
    rows <- seq.int(start, min(n, start + 63L))
    nr <- length(rows)
    draws <- simulate_flexible_targets(adapter, x[rows, , drop = FALSE],
      list(time = rep(0, nr), state = initial[rows], duration = rep(0, nr),
           previous_transitions = rep(0, nr), accrued = matrix(0, nr, model_d)), B = B)
    draws <- draws[, , selected, drop = FALSE]
    center <- variance <- matrix(0, nr, d,
                                  dimnames = list(NULL, colnames(fit$problem$rewards)))
    for (j in seq_len(d)) {
      values <- matrix(draws[, , j], nr, B)
      moments <- flex_empirical_moments(values)
      center[, j] <- moments$mean
      variance[, j] <- moments$variance
    }
    scale <- life_score_scale(variance, fit)
    probability[rows, ] <- flexible_score_cdf(draws, center, scale, q)
  }
  list(probability = probability, q = q, lower = NULL, upper = NULL,
       numerical = list(method = "profile_keyed_continuation_score_cdf",
         working_model = adapter$type, profiles = n, draws = B, seed = adapter$seed,
         deterministic_probability_bounds = FALSE,
         maximum_probability_bound_width = NA_real_,
         approximation = "monte_carlo", batch_size = 64L))
}
