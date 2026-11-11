#' Calibrate a trained prediction rule on independent life histories
#'
#' Estimate the coverage curve for a frozen trained score and choose the first
#' candidate threshold reaching the rank-adjusted coverage cutoff. The default
#' augmented estimator uses partial histories; `method = "ipw"` uses only
#' observed coverage determinations and censoring weights.
#'
#' @param object A fitted object returned by [life_train()].
#' @param data Independent calibration histories in the training formula format.
#' @param alpha Desired joint miscoverage probability, default 0.10.
#' @param alpha_components Optional named vector allocating miscoverage to
#'   every target. Positive entries must sum to at most `alpha`. Each marginal
#'   score is calibrated at its allocated level; their intersection gives the
#'   corresponding joint region by the union bound. NULL uses one joint score.
#' @param method `"aipw"` (default) or `"ipw"` uses the coverage determining
#'   time. `"ipw_full"` weights at absorption or the horizon using the fitted
#'   censoring model. `"complete"` calibrates unweighted scores on histories
#'   observed through absorption or the horizon.
#' @param q Strictly increasing nonnegative finite candidate thresholds. NULL
#'   uses a fixed grid with fine spacing near usual standardized scores and
#'   additional values up to the full-support bound implied by `scale_floor`.
#'   Choose an explicit grid before inspecting calibration outcomes; selecting
#'   it adaptively from those outcomes is not the fixed-grid calibration rule.
#' @param id Identifier column or vector. The stored column rule is reused.
#' @param initial_state Optional calibration initial states; the training rule
#'   is reused when it can be applied to new data.
#' @param check_overlap Logical. By default identifiers shared with training
#'   cause an error. Set FALSE only for independent samples whose identifiers
#'   were numbered separately. Distinct labels alone do not establish independence.
#' @param control Optional numerical overrides of the training controls. Model
#'   fitting and score scale controls remain fixed. For a three-state process
#'   with one absorbing state, two targets ordered as time in the first living
#'   state and total time alive, and Markov-Cox augmentation,
#'   `markov_occupation = "condition_on_alive"` can use scalar occupation
#'   bounds conditional on eligible death times. The default `"general"` uses
#'   the existing joint solver. Unsupported profiles and ambiguous terminal
#'   boundaries retain the general calculation. `calibration_accuracy`
#'   defaults to `"auto"`: Markov augmentation uses `"threshold"` accuracy;
#'   flexible augmentation retains its Monte Carlo evaluation. Explicit
#'   `"probability"` requires every evaluated Markov conditional probability
#'   to meet `probability_tolerance`. `"threshold"` instead
#'   refines only until signed numerical bounds resolve the first crossing
#'   over the complete fixed `q` grid. This mode starts at `grid_size = 16`
#'   unless `grid_size` was explicitly supplied at training or calibration;
#'   a calibration override takes precedence. Direct working-CDF
#'   prediction retains its trained grid and probability-accuracy contract.
#'   Threshold mode supports Markov augmentation and IPW; flexible Monte Carlo
#'   augmentation has no deterministic numerical brackets and is unsupported.
#'   Threshold mode permits up to seven grid doublings per profile, stopping
#'   as soon as the crossing is resolved; `max_nodes` remains enforced in both modes.
#' @param x,response,time,status,start,stop,from,to,event Vector input as in
#'   [life_train()]. The input representation and predictor names must match
#'   the stored schema.
#'
#' @details
#' For `n` calibration subjects, the cutoff is
#' `ceiling((n + 1) * (1 - alpha)) / n`. The first crossing of the raw coverage
#' curve is retained, including when the curve is nonmonotone. Decisions compare
#' raw contribution sums with the integer rank, avoiding a second rounding
#' through ratios. In threshold mode, displayed estimates and numerical bounds
#' divide their raw sums by `n`. A displayed mean can still equal the rounded
#' cutoff at an earlier candidate whose raw sum is below the rank. An unattainable
#' rank or absent crossing returns the full target support.
#'
#' The IPW estimator uses the censoring survival left limit at the observed
#' stopping endpoint. Augmentation uses through-time histories and strictly
#' earlier censoring jumps, with coherent inverse-survival increments. Signed
#' contributions are not clipped. Coverage guarantees require the assumptions
#' and nuisance-model conditions of the method; a fitted working model does
#' not itself give a universal finite-sample coverage guarantee.
#'
#' Coordinate envelopes give valid conservative coverage certificates for
#' multiple targets. They need not give the earliest possible certificate when
#' the joint support has additional constraints. Numerical conditional coverage
#' approximations are recorded separately from these structural certificates.
#'
#' Threshold accuracy compares summed signed bounds against the integer
#' conformal rank. Every earlier candidate must be excluded before a crossing
#' is resolved; no monotonicity is assumed. It can return wide probability
#' brackets away from the selected threshold. Diagnostics record whether the
#' originally requested probability tolerance was met separately from whether
#' the threshold was resolved. An ambiguous threshold after the refinement or
#' node budget is exhausted raises `lifeconformal_threshold_ambiguity`, with
#' the complete curve bounds in the condition. It does not silently return
#' full support. These numerical bounds cover discretization error and exclude
#' floating-point rounding; they are not statistical coverage guarantees.
#'
#' Eligible two-target threshold calculations can start with marginal
#' probability bounds. If these do not resolve the cutoff, selected profiles
#' receive an ordinary joint calculation at the same grid before any grid
#' doubling. Stage diagnostics distinguish marginal, ordinary joint and exact
#' monotone calculations. The representative full curve remains approximate;
#' the signed numerical certificate controls the selected cutoff.
#'
#' @returns A `lifeconformal_calibration` containing the frozen fit, selected
#'   threshold, raw coverage curve, actual cutoff, method, sample size,
#'   full-support flag, and diagnostics. Curve columns `numerical_lower` and
#'   `numerical_upper` bound deterministic numerical error after signed
#'   aggregation; they are not statistical confidence limits. Flexible-model
#'   Monte Carlo calculations do not supply those deterministic bounds.
#'   Use [predict()] for new baseline data.
#' @export
#' @examples
#' histories <- life_example_data(120, seed = 1)
#' fit <- life_train(LifeHist(time, status) ~ age,
#'                   data = histories[histories$id <= 80, ], absorbing = 2,
#'                   targets = list(alive = 0), tau = 5)
#' cal <- life_calibrate(fit, data = histories[histories$id > 80, ])
#' predict(cal, newdata = data.frame(age = c(50, 70)))
life_calibrate <- function(object, data = NULL, alpha = 0.10,
                            method = c("aipw", "ipw", "ipw_full", "complete"), q = NULL,
                            alpha_components = NULL,
                            id = NULL, initial_state = NULL,
                            check_overlap = TRUE, control = list(),
                            x = NULL, response = NULL, time = NULL,
                            status = NULL, start = NULL, stop = NULL,
                            from = NULL, to = NULL, event = NULL) {
  call <- match.call()
  if (!inherits(object, "lifeconformal_fit"))
    base::stop("object must be returned by life_train().", call. = FALSE)
  # Preserve the formerly unambiguous abbreviations for determining-time IPW.
  if (length(method) == 1L && method %in% c("i", "ip")) method <- "ipw"
  method <- match.arg(method)
  if (!is.numeric(alpha) || is.complex(alpha) || !is.null(dim(alpha)) ||
      length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1)
    base::stop("alpha must be strictly between zero and one.", call. = FALSE)
  if (!is.logical(check_overlap) || length(check_overlap) != 1L ||
      is.na(check_overlap))
    base::stop("check_overlap must be TRUE or FALSE.", call. = FALSE)
  controls <- resolve_calibration_controls(object, control)
  if (!is.null(alpha_components)) {
    result <- component_calibration(object, alpha_components, alpha,
      calibration_args = list(data = data, method = method, q = q, id = id,
        initial_state = initial_state, check_overlap = check_overlap,
        control = control, x = x, response = response, time = time,
        status = status, start = start, stop = stop, from = from, to = to,
        event = event))
    result$call <- call
    return(result)
  }
  context <- prepare_life_calibration_context(object, data, id, initial_state,
    check_overlap, x, response, time, status, start, stop, from, to, event)
  calibrate_life_prepared(object, context, alpha, method, q, controls, call)
}

#' Evaluate augmentation at a common endpoint with exact structural overrides
#' @noRd
evaluate_augmented_curve <- function(fit, history, x, censor_x, moments, cdt, q,
                                      controls, queries = NULL, rank = NULL) {
  orientation <- if (fit$augmentation$type == "markov_cox")
    orient_markov_rewards(fit$problem$rewards, fit$augmentation$edges) else NULL
  if (!is.null(orientation) &&
      observed_monotone_histories(history, orientation$rewards)) {
    return(evaluate_monotone_augmented_curve(fit, history, x, censor_x,
                                              moments, cdt, q, orientation, queries))
  }
  if (identical(controls$calibration_accuracy, "threshold")) {
    requested_alive <- identical(controls$markov_occupation, "condition_on_alive")
    if (requested_alive && alive_occupation_supported(fit$problem$rewards, fit$problem$absorbing)) {
      conditioned <- tryCatch(alive_threshold_curve(fit, history, x, censor_x,
        moments, cdt, q, controls, queries, rank), lifeconformal_grid_budget=function(e) NULL,
        lifeconformal_grid_precision=function(e) NULL)
      if (!is.null(conditioned)) return(conditioned)
      controls$markov_occupation <- "general"
    }
    if (is.null(queries)) queries <- prepare_calibration_queries(fit, history, censor_x, cdt)
    result <- evaluate_threshold_augmented_curve(fit, history, x, moments, cdt, q,
      controls, queries, rank)
    if (requested_alive) result$diagnostics$occupation_engine <- "general_fallback"
    return(result)
  }
  n <- nrow(history$subjects)
  result <- matrix(0, n, length(q))
  lower_result <- upper_result <- if (fit$augmentation$type == "markov_cox") result else NULL
  if (is.null(queries)) queries <- prepare_calibration_queries(fit, history, censor_x, cdt)
  prepared <- queries$prepared
  if (fit$augmentation$type != "markov_cox" && is.null(prepared$flexible) &&
      !identical(colnames(fit$problem$rewards),
                 colnames(fit$augmentation$problem$rewards))) {
    prepared$flexible <- prepare_occupation_history(history, fit$augmentation$problem)
  }
  numerical <- c(list(conditional_coverage_evaluated = TRUE,
                     max_probability_bound_width = 0, queried_snapshots = 0L,
                     structural_overrides = 0L), queries$diagnostics)
  if (identical(controls$markov_occupation, "condition_on_alive"))
    numerical$alive_profile_counts <- c(window=0L, conditioned=0L, fallback=0L, terminal_boundary=0L)
  for (i in seq_len(n)) {
    first <- if (i == 1L) 1L else queries$endpoint_index[i - 1L] + 1L
    rows <- seq.int(first, queries$endpoint_index[i])
    times <- queries$time[rows]
    weight <- queries$weight[rows]
    snapshot <- list(time = times, state = queries$snapshot$state[rows],
      accrued = queries$snapshot$accrued[rows, , drop = FALSE])
    h <- length(times)
    if (all(cdt$resolved[i, ]) && times[1L] >= max(cdt$time[i, ])) {
      # Every censoring jump is at or after determination: the completed
      # common-endpoint expression telescopes exactly to known membership.
      value <- as.numeric(cdt$membership[i, ])
      result[i, ] <- value
      if (!is.null(lower_result))
        lower_result[i, ] <- upper_result[i, ] <- value
      numerical$structural_overrides <- numerical$structural_overrides + h * length(q)
      next
    }
    certified <- outer(times, cdt$time[i, ], `>=`) &
      matrix(rep(cdt$resolved[i, ], each = h), h)
    known <- matrix(rep(as.numeric(cdt$membership[i, ]), each = h), h)
    known[!certified] <- NA_real_
    missing <- rowSums(is.na(known))
    need <- which(missing > 0L)
    values <- known
    lower_values <- upper_values <- known
    streamed <- NULL
    if (length(need)) {
      if (!is.null(lower_result) && is.null(prepared$markov_jumps))
        prepared$markov_jumps <- calibration_markov_jumps(fit, x, prepared)
      evaluated <- with_life_grid_context(conditional_life_coverage(fit, x[i, , drop = FALSE],
        snapshot, need, moments$mean[i, ], moments$scale[i, ], q, controls,
        prepared = prepared, subject_index = i,
        known = known[need, , drop = FALSE],
        query_weights = if (!is.null(lower_result)) weight[need] else NULL),
        profile_index = i, stage = "calibration", subject_id = history$subjects$id[i])
      if (!is.null(numerical$alive_profile_counts)) {
        engine <- evaluated$diagnostics$method
        kind <- if (identical(engine, "alive_time_window")) "window" else
          if (identical(engine, "conditioned_on_alive_time")) "conditioned" else "fallback"
        numerical$alive_profile_counts[kind] <- numerical$alive_profile_counts[kind] + 1L
        if (identical(evaluated$diagnostics$alive_fallback, "terminal atom boundary"))
          numerical$alive_profile_counts["terminal_boundary"] <- numerical$alive_profile_counts["terminal_boundary"]+1L
      }
      streamed <- evaluated$reduced
      if (is.null(streamed)) {
        values[need, ] <- evaluated$probability
        if (!is.null(lower_result)) {
          lower_values[need, ] <- evaluated$lower
          upper_values[need, ] <- evaluated$upper
        }
      }
      numerical$max_probability_bound_width <- max(
        numerical$max_probability_bound_width, evaluated$bound_width)
      numerical$queried_snapshots <- numerical$queried_snapshots + length(need)
    }
    numerical$structural_overrides <- numerical$structural_overrides +
      length(known) - sum(missing)
    if (!is.null(streamed)) {
      # The native pass compensates all signed unresolved-query sums and
      # applies structural completion before weighting. Fully certified rows
      # need no distribution calculation. Their membership is constant in
      # time, so sum their weights once with R's extended-precision sum where
      # available. Adding the two contributions only once
      # avoids retaining probability and bracket matrices for every query.
      certain <- which(missing == 0L)
      fixed <- if (length(certain))
        known[certain[1L], ] * sum(weight[certain]) else numeric(length(q))
      value <- streamed$probability + fixed
      low <- streamed$lower + fixed
      high <- streamed$upper + fixed
    } else {
      exact <- !is.na(known)
      values[exact] <- known[exact]
      lower_values[exact] <- upper_values[exact] <- known[exact]
      value <- values[h, ] * weight[h]
      if (!is.null(lower_result)) {
        low <- lower_values[h, ] * weight[h]
        high <- upper_values[h, ] * weight[h]
      }
      if (h > 1L) {
        value <- value + colSums(values[-h, , drop = FALSE] * weight[-h])
        if (!is.null(lower_result)) {
          low <- low + colSums(upper_values[-h, , drop = FALSE] * weight[-h])
          high <- high + colSums(lower_values[-h, , drop = FALSE] * weight[-h])
        }
      }
    }
    origin <- cdt$resolved[i, ] & cdt$time[i, ] == 0
    value[origin] <- as.numeric(cdt$membership[i, origin])
    result[i, ] <- value
    if (!is.null(lower_result)) {
      low[origin] <- high[origin] <- value[origin]
      lower_result[i, ] <- low
      upper_result[i, ] <- high
    }
  }
  if (!is.null(lower_result)) {
    numerical$curve_lower_sum <- colSums(lower_result)
    numerical$curve_upper_sum <- colSums(upper_result)
    numerical$curve_lower <- numerical$curve_lower_sum / n
    numerical$curve_upper <- numerical$curve_upper_sum / n
    numerical$max_curve_bound_width <- max(numerical$curve_upper - numerical$curve_lower)
  }
  numerical$conditional_coverage_evaluated <- numerical$queried_snapshots > 0L
  list(contributions = result, diagnostics = numerical)
}

#' Evaluate all thresholds together for conditional target predictions
#' @noRd
conditional_life_coverage <- function(fit, x, snapshot, need, center, scale, q,
                                       controls, prepared, subject_index, known,
                                       query_weights = NULL, grid_refine = TRUE,
                                       use_marginal = FALSE) {
  adapter <- fit$augmentation
  if (adapter$type == "markov_cox") {
    jumps <- prepared$markov_jumps
    profile <- subject_index
    if (is.null(jumps)) {
      jumps <- predict_transition_jumps(adapter, x)
      profile <- 1L
    }
    k <- length(fit$problem$states)
    h <- length(need)
    orientation <- orient_markov_rewards(fit$problem$rewards, jumps$edges)
    transformed <- if (!is.null(orientation)) orient_markov_query(orientation,
      matrix(rep(center, each = h), h), snapshot$accrued[need, , drop = FALSE],
      snapshot$time[need], fit$problem$tau) else NULL
    coherent <- !is.null(transformed)
    if (coherent && any(orientation$flip)) {
      active <- orientation$rewards[snapshot$state[need], , drop = FALSE] == 1
      coherent <- all(transformed$accrued[active] ==
                        matrix(snapshot$time[need], h, ncol(active))[active]) &&
        all(transformed$original_accrued[, orientation$flip, drop = FALSE][
          active[, orientation$flip, drop = FALSE]] == 0)
    }
    if (coherent) {
      kernels <- lapply(seq_along(jumps$jump_times), function(j) {
        inc <- matrix(0, k, k)
        inc[cbind(jumps$edges$from, jumps$edges$to)] <-
          jumps$base_increments[j, ] * jumps$relative_risk[profile, ]
        make_markov_kernel(inc, absorbing = fit$problem$absorbing)$kernel
      })
      evaluated <- markov_conditional_coverage_monotone(jumps$jump_times, kernels,
        transformed$rewards, fit$problem$tau, snapshot$time[need],
        snapshot$state[need], transformed$accrued, transformed$center,
        matrix(rep(scale, each = h), h), q, side = "right",
        flip = transformed$flip, original_center = transformed$original_center,
        original_accrued = transformed$original_accrued)
      return(list(probability = evaluated$probability, lower = evaluated$probability,
                   upper = evaluated$probability, bound_width = 0,
                   diagnostics = evaluated$diagnostics))
    }
    grid_evaluator <- if (isTRUE(use_marginal))
      markov_conditional_coverage_marginal else markov_conditional_coverage_grid
    alive <- identical(controls$markov_occupation, "condition_on_alive") &&
      alive_occupation_supported(fit$problem$rewards, fit$problem$absorbing)
    args <- list(jump_times=jumps$jump_times, kernels=NULL,
      rewards=fit$problem$rewards, tau=fit$problem$tau, query_times=snapshot$time[need],
      query_state=snapshot$state[need], accrued=snapshot$accrued[need, , drop = FALSE],
      center=matrix(rep(center, each=h), h), scale=matrix(rep(scale, each=h), h), q=q,
      side = "right", grid_size = controls$grid_size, max_nodes = controls$max_nodes,
      tolerance = controls$probability_tolerance, strict = grid_refine,
      max_refinements = if (grid_refine) 4L else 0L,
      # Threshold certification can use a smaller numerical partition while
      # retaining every mass and envelope. The signed bounds, not the point
      # approximation, must resolve the full candidate-grid crossing.
      exact_limit = if (grid_refine) 256L else 32L,
      point_mass_threshold = if (grid_refine) .001 else .1,
      query_known = known, query_weights = query_weights,
      kernel_factors = life_markov_kernel_factors(jumps, profile, fit$problem$absorbing))
    evaluated <- if (alive)
      alive_conditional_coverage(args, grid_evaluator, screen=length(q)>8L) else
      do.call(grid_evaluator,args)
    return(list(probability = evaluated$probability, lower = evaluated$lower,
                 upper = evaluated$upper, reduced = evaluated$reduced,
                 bound_width = evaluated$diagnostics$maximum_probability_error,
                 diagnostics = evaluated$diagnostics))
  }
  if (!is.null(prepared$flexible)) {
    snapshot <- evaluate_occupation_history(prepared$flexible,
      rep(subject_index, length(need)), snapshot$time[need], side = "right")
    need <- seq_along(need)
    prepared <- prepared$flexible
  }
  snaps <- flexible_history_snapshots(prepared, subject_index, snapshot, need)
  draws <- simulate_flexible_targets(adapter,
    x[rep(1L, length(need)), , drop = FALSE], snaps, B = controls$draws)
  selected <- match(colnames(fit$problem$rewards), colnames(adapter$problem$rewards))
  draws <- draws[, , selected, drop = FALSE]
  h <- length(need)
  probability <- flexible_score_cdf(draws, matrix(rep(center, each = h), h),
                                    matrix(rep(scale, each = h), h), q)
  list(probability = probability, bound_width = NA_real_)
}

#' Recover predictable duration and transition counts for flexible continuations
#' @noRd
flexible_history_snapshots <- function(prepared, subject_index, snapshot, need) {
  rows <- prepared$event_rows[[subject_index]]
  event_times <- prepared$intervals$stop[rows]
  previous <- findInterval(snapshot$time[need], event_times)
  entry <- c(0, event_times)[previous + 1L]
  list(time = snapshot$time[need], state = snapshot$state[need],
       accrued = snapshot$accrued[need, , drop = FALSE],
       duration = snapshot$time[need] - entry, previous_transitions = previous)
}
