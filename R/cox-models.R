#' Prepare subject-level censoring risk records
#'
#' @param history A validated `life_history_data` object.
#' @param x Finite numeric subject-by-predictor matrix, with unique column
#'   names. A zero-column matrix specifies an intercept-only model.
#' @param tau Positive finite prediction horizon.
#' @returns A list with `data` (subject index, time, censoring event), aligned
#'   `x`, and counts of subjects excluded because their path was complete at
#'   the origin. A genuine censoring observation at zero is retained.
#' @details Follow-up ends at the earliest of reported observation, absorption,
#'   and the horizon. Only actual censoring strictly before the horizon is a
#'   censoring event. Absorption and horizon restriction introduce no
#'   censoring events. Preparation uses O(n p) storage and O(n) indexing.
#' @noRd
make_censoring_data <- function(history, x, tau) {
  .cox_check_history_design(history, x, tau)
  subjects <- history$subjects
  absorption <- subjects$absorption_time
  absorption[is.na(absorption)] <- Inf
  keep <- absorption > 0
  time <- pmin(subjects$observation_end, absorption, tau)
  event <- subjects$censoring & subjects$observation_end < tau &
    subjects$observation_end <= absorption
  list(data = data.frame(subject_index = subjects$subject_index[keep],
                         time = time[keep], event = as.integer(event[keep])),
       x = x[keep, , drop = FALSE],
       excluded_origin_absorption = sum(!keep), n_subjects = nrow(subjects))
}

#' Fit a frozen Cox censoring adapter
#'
#' @inheritParams make_censoring_data
#' @param options Named list. `ties` is `"breslow"` (default) or `"efron"`;
#'   `singular` is `"drop"` (warn and set aliased coefficients to zero) or
#'   `"error"`; `iter.max` defaults to 50 and `eps` to 1e-9. Censoring also
#'   accepts `survival_floor`, a value in (0, 1]. Its default is the smallest
#'   normalized positive double, a computational guard rather than a
#'   statistical positivity threshold. Public controls may set a larger floor.
#' @returns A serializable `life_censor_cox` adapter with its centered Cox
#'   baseline step law, saved backend fit, controls, and diagnostics.
#' @details Uses `survival::coxph()` when regression coefficients need fitting,
#'   followed by a direct centered baseline
#'   extraction on the observed event times. Predictions use stored numeric state; they do
#'   not reconstruct a model frame or call `survfit()` per time or subject.
#'   Near but unequal times remain distinct (`timefix = FALSE`). Fits with no
#'   censoring events return survival one and do not invent a failure. Fits
#'   reporting nonconvergence or potentially infinite coefficients fail with
#'   a diagnostic instead of silently supplying an unstable working model.
#' @noRd
fit_censor_cox <- function(history, x, tau, options = list()) {
  controls <- .cox_resolve_options(options, censoring = TRUE)
  risk <- make_censoring_data(history, x, tau)
  model <- .cox_fit_numeric(risk$data$time, risk$data$event, risk$x,
                            controls, label = "Censoring Cox model")
  structure(list(type = "cox", backend = "survival",
                 backend_version = as.character(utils::packageVersion("survival")),
                 model = model, tau = tau, controls = controls,
                 follow_up_support = c(0, model$diagnostics$maximum_risk_time),
                 predictor_names = colnames(x), n_predictors = ncol(x),
                 jump_times = model$jump_times,
                 temporal_law = "right_continuous_exponential_cumulative_hazard",
                 diagnostics = c(model$diagnostics,
                   list(excluded_origin_absorption = risk$excluded_origin_absorption,
                        n_subjects = risk$n_subjects,
                        administrative_censoring_events = 0L))),
            class = "life_censor_cox")
}

#' Predict aligned left or right censoring survival
#'
#' @param adapter A fitted `life_censor_cox` adapter.
#' @param x New numeric predictor matrix with the frozen column order.
#' @param time Nonnegative query time(s) no greater than the fitted horizon.
#'   A scalar is recycled; otherwise there is one value per query.
#' @param side `"right"` includes a jump at the query time; `"left"` excludes
#'   it. In particular, the left value at zero is one.
#' @param subject_index Optional predictor-row index for each query, allowing
#'   many time queries to share each once-computed linear predictor. With NULL,
#'   one query per predictor row is expected.
#' @returns A numeric vector in query order. The same survival floor is applied
#'   to the entire step law before either side or any increment is evaluated.
#' @details Linear predictors are computed once per input profile; ordered
#'   baseline lookup uses `findInterval()`. Working storage is O(n + m), with
#'   n profiles and m queries. The log cumulative-hazard calculation avoids
#'   overflowing relative risks for extreme new censoring profiles.
#' @noRd
predict_censor_survival <- function(adapter, x, time,
                                    side = c("right", "left"),
                                    subject_index = NULL) {
  if (!inherits(adapter, "life_censor_cox")) {
    stop("`adapter` must be a life_censor_cox model.", call. = FALSE)
  }
  side <- match.arg(side)
  query <- .cox_censor_queries(adapter, x, time, subject_index)
  lp <- .cox_linear_predictor(adapter$model, x)
  index <- findInterval(query$time, adapter$jump_times)
  if (side == "left" && length(adapter$jump_times)) {
    equal <- index > 0L
    equal[equal] <- adapter$jump_times[index[equal]] == query$time[equal]
    index <- index - equal
  }
  cumulative <- c(0, adapter$model$cumulative_hazard)[index + 1L]
  exp(.cox_log_survival(cumulative, lp[query$subject_index],
                        adapter$controls$survival_floor))
}

#' Predict coherent censoring jump increments
#'
#' @inheritParams predict_censor_survival
#' @returns A data frame in query order with `subject_index`, `time`,
#'   `g_left`, `g_right`, `hazard_increment` and `inverse_increment`.
#'   Queries away from fitted jumps have zero increments. The hazard jump is
#'   `1 - g_right/g_left`, not the Cox cumulative-hazard increment. The inverse
#'   jump is `1/g_right - 1/g_left`. Both derive from the floored survival law.
#' @details This routine evaluates only the supplied profile/time pairs; it
#'   does not form a full subject-by-jump cross product. Repeated subject
#'   indices reuse linear predictors. Increments use `expm1()` to retain
#'   precision for small jumps, including jumps crossing the survival floor.
#' @noRd
predict_censor_jumps <- function(adapter, x, time, subject_index = NULL) {
  if (!inherits(adapter, "life_censor_cox")) {
    stop("`adapter` must be a life_censor_cox model.", call. = FALSE)
  }
  query <- .cox_censor_queries(adapter, x, time, subject_index)
  lp <- .cox_linear_predictor(adapter$model, x)[query$subject_index]
  right <- findInterval(query$time, adapter$jump_times)
  equal <- right > 0L
  if (length(adapter$jump_times)) {
    equal[equal] <- adapter$jump_times[right[equal]] == query$time[equal]
  }
  left <- right - equal
  cumulative <- c(0, adapter$model$cumulative_hazard)
  log_left <- .cox_log_survival(cumulative[left + 1L], lp,
                                adapter$controls$survival_floor)
  log_right <- .cox_log_survival(cumulative[right + 1L], lp,
                                 adapter$controls$survival_floor)
  g_right <- exp(log_right)
  hazard <- -expm1(log_right - log_left)
  data.frame(subject_index = query$subject_index, time = query$time,
             g_left = exp(log_left), g_right = g_right,
             hazard_increment = hazard, inverse_increment = hazard / g_right)
}

#' Prepare the start-stop risk set for one transition
#'
#' @inheritParams make_censoring_data
#' @param transition Two internal state indices, origin followed by destination.
#' @returns A list with counting-process `data` and aligned `x`. Other exits
#'   from the origin are censoring observations for this cause-specific fit.
#' @details Uses all observed risk through actual follow-up, with time since
#'   the common prediction origin. Recurrent visits create separate (start,
#'   stop] intervals. Zero-duration terminal records carry no risk. The horizon
#'   validates the prediction problem but does not truncate outcome fitting.
#' @noRd
make_transition_risk_data <- function(history, transition, x, tau) {
  .cox_check_history_design(history, x, tau)
  k <- length(history$schema$state_labels)
  if (!is.numeric(transition) || length(transition) != 2L ||
      anyNA(transition) || any(!is.finite(transition)) ||
      any(transition != floor(transition) | transition < 1 | transition > k) ||
      transition[1L] == transition[2L]) {
    stop("`transition` must contain two distinct valid state indices.", call. = FALSE)
  }
  intervals <- history$intervals
  rows <- which(intervals$from == transition[1L] & intervals$stop > intervals$start)
  .cox_transition_view(intervals, rows, transition[2L], x, tau)
}

#' Fit transition-specific Cox models on a frozen empirical graph
#'
#' @inheritParams make_censoring_data
#' @inheritParams fit_censor_cox
#' @param edges Optional data frame with internal `from` and `to` indices.
#'   It must contain exactly the observed training pairs, in the desired
#'   frozen order. Extra graph metadata columns are ignored. NULL uses order
#'   of first appearance. Pairs observed only after tau remain in the graph
#'   with zero fitted increments over the prediction horizon, although their
#'   regression coefficients are estimated from their full observed risk sets.
#' @returns A serializable `life_markov_cox` adapter. `jump_times`,
#'   `base_increments` (jump by edge), and `edges` share one edge order.
#'   `models` retain the centered baseline and coefficient state for each edge.
#' @details No empirical absence is treated as a structural impossibility.
#'   Each origin's risk-row index is built once using all actual observed
#'   follow-up. Centered baseline increments are extracted only through tau.
#'   This function supplies intensity
#'   increments, not a multistate `survfit()` transition-probability curve.
#'   The deterministic Markov calculations normalize rows whose sums exceed
#'   one after the profile risks are combined with these increments.
#' @noRd
fit_markov_cox <- function(history, x, tau, options = list(), edges = NULL) {
  .cox_check_history_design(history, x, tau)
  controls <- .cox_resolve_options(options, censoring = FALSE)
  intervals <- history$intervals
  observed <- unique(intervals[intervals$event, c("from", "to"), drop = FALSE])
  rownames(observed) <- NULL
  if (is.null(edges)) edges <- observed
  edges <- .cox_validate_edges(edges, observed)
  edge_names <- paste(edges$from, edges$to, sep = "->")
  valid_rows <- which(intervals$stop > intervals$start)
  origin_rows <- split(valid_rows, intervals$from[valid_rows])
  models <- vector("list", nrow(edges))
  names(models) <- edge_names
  for (e in seq_len(nrow(edges))) {
    rows <- origin_rows[[as.character(edges$from[e])]]
    if (is.null(rows)) rows <- integer()
    risk <- .cox_transition_view(intervals, rows, edges$to[e], x, tau)
    displayed_edge <- paste(history$schema$state_labels[edges$from[e]],
                             history$schema$state_labels[edges$to[e]], sep = "->")
    sparse <- ncol(risk$x) > 0L && sum(risk$data$event) < controls$min_events_per_transition
    design <- if (sparse) risk$x[, FALSE, drop = FALSE] else risk$x
    models[[e]] <- .cox_fit_numeric(risk$data$stop, risk$data$event, design,
                                    controls, start = risk$data$start,
                                    label = paste("Transition", displayed_edge, "Cox model"),
                                    baseline_horizon = tau)
    if (sparse) {
      models[[e]]$coefficients <- stats::setNames(numeric(ncol(x)), colnames(x))
      models[[e]]$center <- models[[e]]$coefficients
      models[[e]]$diagnostics$coefficient_fit <- "sparse_intercept_only"
    }
  }
  jump_times <- sort(unique(unlist(lapply(models, `[[`, "jump_times"),
                                   use.names = FALSE)))
  if (is.null(jump_times)) jump_times <- numeric()
  base_increments <- matrix(0, length(jump_times), nrow(edges),
                             dimnames = list(NULL, edge_names))
  for (e in seq_along(models)) {
    base_increments[match(models[[e]]$jump_times, jump_times), e] <-
      models[[e]]$increments
  }
  center <- if (ncol(x)) colMeans(x) else numeric()
  coefficients <- matrix(0, ncol(x), nrow(edges),
                          dimnames = list(colnames(x), edge_names))
  center_offsets <- numeric(nrow(edges))
  for (e in seq_along(models)) {
    coefficients[, e] <- models[[e]]$coefficients
    center_offsets[e] <- sum((center - models[[e]]$center) * models[[e]]$coefficients)
  }
  structure(list(type = "markov_cox", backend = "survival",
                 backend_version = as.character(utils::packageVersion("survival")),
                 tau = tau, edges = edges, jump_times = jump_times,
                 base_increments = base_increments, models = models,
                 coefficients = coefficients, center = center,
                 center_offsets = center_offsets, controls = controls,
                 predictor_names = colnames(x), n_predictors = ncol(x),
                 absorbing = history$schema$absorbing,
                 state_labels = history$schema$state_labels,
                 temporal_law = "deterministic_jump_markov",
                 diagnostics = list(ties = controls$ties, timefix = FALSE,
                   n_edges = nrow(edges), n_jumps = length(jump_times),
                   zero_event_edges = edge_names[vapply(models,
                     function(m) m$diagnostics$n_events == 0L, logical(1))],
                   zero_prediction_event_edges = edge_names[vapply(models,
                     function(m) m$diagnostics$n_prediction_events == 0L, logical(1))],
                   by_edge = lapply(models, `[[`, "diagnostics"))),
            class = "life_markov_cox")
}

#' Predict the factored Markov intensity jumps for a profile batch
#'
#' @param adapter A `life_markov_cox` adapter.
#' @param x New numeric profile-by-predictor matrix in frozen column order.
#' @param times Optional ordered time window, supplied as two endpoints in
#'   [0, tau]. Only jumps strictly after its first and at or before its second
#'   endpoint are returned. NULL returns the complete fitted grid.
#' @returns A list containing `jump_times`, shared `base_increments`, `edges`,
#'   and profile-by-edge `relative_risk`. No profile kernel lists are built.
#' @details Matrix multiplication computes all edge linear predictors at once
#'   after subtracting a common training center. Edge-specific center offsets
#'   retain the reference covariates used in the baseline fits. Working storage
#'   is O(n p + n E + J E), not O(n J K squared). Unrepresentable positive
#'   relative risks fail explicitly; they are not silently clipped.
#' @noRd
predict_transition_jumps <- function(adapter, x, times = NULL) {
  if (!inherits(adapter, "life_markov_cox")) {
    stop("`adapter` must be a life_markov_cox model.", call. = FALSE)
  }
  .cox_check_design(x, adapter$n_predictors, adapter$predictor_names)
  used <- if (ncol(x)) rowSums(abs(adapter$coefficients)) > 0 else logical()
  lp <- if (any(used)) {
    sweep(x[, used, drop = FALSE], 2L, adapter$center[used], "-") %*%
      adapter$coefficients[used, , drop = FALSE]
  } else matrix(0, nrow(x), nrow(adapter$edges))
  if (ncol(lp)) lp <- sweep(lp, 2L, adapter$center_offsets, "+")
  risks <- exp(lp)
  if (any(!is.finite(risks)) || any(risks <= 0)) {
    stop("New profile relative risks are outside the representable Cox range.",
         call. = FALSE)
  }
  dimnames(risks) <- list(rownames(x), colnames(adapter$base_increments))
  keep <- rep(TRUE, length(adapter$jump_times))
  if (!is.null(times)) {
    if (!is.numeric(times) || length(times) != 2L || any(!is.finite(times)) ||
        times[1L] < 0 || times[2L] > adapter$tau || times[1L] > times[2L]) {
      stop("`times` must be an ordered two-endpoint window within [0, tau].",
           call. = FALSE)
    }
    keep <- adapter$jump_times > times[1L] & adapter$jump_times <= times[2L]
  }
  list(jump_times = adapter$jump_times[keep],
       base_increments = adapter$base_increments[keep, , drop = FALSE],
       edges = adapter$edges, relative_risk = risks)
}

#' Compute fitted Markov-Cox target moments without path simulation
#'
#' @inheritParams predict_transition_jumps
#' @param rewards Numeric state-by-target reward matrix.
#' @param initial_state One internal state index or one per profile.
#' @param accrued Optional common or profile-specific accrued target rewards.
#' @param keep_state_moments Also return moments for every possible initial state.
#' @returns The output of [markov_reward_moments_batch()], including kernel
#'   adjustment diagnostics. This is a baseline moment interface; it does not
#'   implement a conditional target distribution for partial histories.
#' @noRd
predict_markov_cox_moments <- function(adapter, x, rewards, initial_state,
                                       accrued = NULL, keep_state_moments = FALSE) {
  if (!is.matrix(rewards) || nrow(rewards) != length(adapter$state_labels)) {
    stop("`rewards` must have one row per fitted state.", call. = FALSE)
  }
  jumps <- predict_transition_jumps(adapter, x)
  markov_reward_moments_batch(jumps$jump_times, jumps$base_increments,
                               jumps$edges, jumps$relative_risk, rewards,
                               adapter$tau, initial_state, accrued,
                               adapter$absorbing, keep_state_moments)
}

#' Validate numeric Cox designs and the frozen predictor schema
#' @noRd
.cox_check_design <- function(x, p = NULL, predictor_names = NULL,
                               allow_empty_rows = FALSE) {
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) ||
      any(!is.finite(x)) || (!allow_empty_rows && nrow(x) < 1L)) {
    stop("`x` must be a finite numeric subject-by-predictor matrix.", call. = FALSE)
  }
  if (ncol(x) && (is.null(colnames(x)) || anyNA(colnames(x)) ||
      any(!nzchar(colnames(x))) || anyDuplicated(colnames(x)))) {
    stop("`x` must have unique nonempty predictor column names.", call. = FALSE)
  }
  if (!is.null(p) && (ncol(x) != p ||
      (p > 0L && !identical(colnames(x), predictor_names)))) {
    stop("`x` does not match the frozen Cox predictor columns and order.",
         call. = FALSE)
  }
  invisible(x)
}

#' Validate history/design alignment before derived risk construction
#' @noRd
.cox_check_history_design <- function(history, x, tau) {
  if (!inherits(history, "life_history_data")) {
    stop("`history` must be a validated life_history_data object.", call. = FALSE)
  }
  .cox_check_design(x)
  if (nrow(x) != nrow(history$subjects)) {
    stop("`x` must have one row per canonical subject.", call. = FALSE)
  }
  if (!is.numeric(tau) || is.complex(tau) || length(tau) != 1L ||
      !is.finite(tau) || tau <= 0) {
    stop("`tau` must be a finite positive scalar.", call. = FALSE)
  }
  invisible(NULL)
}

#' Resolve explicit Cox controls without passing arbitrary backend arguments
#' @noRd
.cox_resolve_options <- function(options, censoring) {
  defaults <- list(ties = "breslow", singular = "drop", iter.max = 50L, eps = 1e-9)
  if (censoring) {
    defaults$survival_floor <- .Machine$double.xmin
  } else {
    defaults$min_events_per_transition <- 0L
    defaults$separation_action <- "error"
  }
  if (!is.list(options) || (length(options) && (is.null(names(options)) ||
      anyNA(names(options)) || any(!nzchar(names(options))) ||
      anyDuplicated(names(options)) || any(!names(options) %in% names(defaults))))) {
    stop("Cox options must be a uniquely named list of supported controls: ",
         paste(names(defaults), collapse = ", "), ".", call. = FALSE)
  }
  defaults[names(options)] <- options
  if (!censoring) {
    k <- defaults$min_events_per_transition
    if (!is.numeric(k) || is.complex(k) || !is.null(dim(k)) || length(k) != 1L ||
        !is.finite(k) || k < 0 || k != floor(k))
      stop("min_events_per_transition must be a nonnegative whole number.", call. = FALSE)
  }
  valid_choice <- function(value, choices) is.character(value) &&
    length(value) == 1L && !is.na(value) && value %in% choices
  if (!censoring && !valid_choice(defaults$separation_action, c("error", "baseline")))
    stop("Markov Cox separation_action must be 'error' or 'baseline'.", call. = FALSE)
  if (!valid_choice(defaults$ties, c("breslow", "efron")) ||
      !valid_choice(defaults$singular, c("drop", "error"))) {
    stop("Cox `ties` must be breslow/efron and `singular` drop/error.", call. = FALSE)
  }
  if (!is.numeric(defaults$iter.max) || length(defaults$iter.max) != 1L ||
      !is.finite(defaults$iter.max) || defaults$iter.max < 1 ||
      defaults$iter.max != floor(defaults$iter.max) ||
      defaults$iter.max > .Machine$integer.max) {
    stop("Cox `iter.max` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(defaults$eps) || length(defaults$eps) != 1L ||
      !is.finite(defaults$eps) || defaults$eps <= 0 || defaults$eps >= 1) {
    stop("Cox `eps` must be strictly between zero and one.", call. = FALSE)
  }
  if (censoring && (!is.numeric(defaults$survival_floor) ||
      length(defaults$survival_floor) != 1L || !is.finite(defaults$survival_floor) ||
      defaults$survival_floor < .Machine$double.xmin || defaults$survival_floor > 1)) {
    stop("Cox `survival_floor` must be between .Machine$double.xmin and one.",
         call. = FALSE)
  }
  defaults$iter.max <- as.integer(defaults$iter.max)
  defaults
}

#' Fit one numeric Cox risk view and extract a centered step baseline once
#' @noRd
.cox_fit_numeric <- function(stop, event, x, controls, start = NULL, label,
                              baseline_horizon = Inf) {
  p <- ncol(x)
  n <- nrow(x)
  n_events <- sum(event)
  zero <- stats::setNames(numeric(p), colnames(x))
  diagnostics <- list(n_risk_records = n, n_events = n_events,
                       n_prediction_events = sum(event[stop <= baseline_horizon]),
                       maximum_risk_time = if (n) max(stop) else 0,
                       no_events = n_events == 0L, singular_columns = character(),
                       warnings = character(), converged = TRUE,
                       iterations = 0L, ties = controls$ties, timefix = FALSE)
  if (!n_events) {
    return(list(fit = NULL, coefficients = zero,
                 center = if (n && p) colMeans(x) else zero,
                 jump_times = numeric(), increments = numeric(),
                 cumulative_hazard = numeric(), diagnostics = diagnostics))
  }
  # An intercept-only Cox model has no regression coefficients to optimize.
  # A single risk record likewise cannot identify any covariate effect. The
  # existing observed-time sweep gives the complete baseline in either case,
  # without asking the backend to compute irrelevant concordance statistics
  # (which also lose their matrix dimensions for n = 1 in older survival).
  if (!p || n == 1L) {
    if (p) {
      diagnostics$singular_columns <- colnames(x)
      detail <- paste(colnames(x), collapse = ", ")
      if (controls$singular == "error")
        stop(label, " has unidentifiable predictor columns: ", detail, ".",
             call. = FALSE)
      warning(label, " has unidentifiable predictor columns; coefficients set to zero: ",
              detail, ".", call. = FALSE)
    }
    baseline <- .cox_centered_baseline(stop, event, numeric(n), start,
                                       controls$ties, baseline_horizon)
    cumulative <- cumsum(baseline$increments)
    if (any(!is.finite(cumulative)))
      stop(label, " produced an unrepresentable centered cumulative hazard.", call. = FALSE)
    diagnostics$baseline_method <- "observed_time_log_risk_sweep"
    diagnostics$coefficient_fit <- if (p) "single_record_singular" else "intercept_only"
    return(list(fit = NULL, coefficients = zero,
                 center = if (p) stats::setNames(as.numeric(x[1L, ]), colnames(x)) else zero,
                 jump_times = baseline$time, increments = baseline$increments,
                 cumulative_hazard = cumulative, diagnostics = diagnostics))
  }
  risk_frame <- data.frame(.stop = stop, .event = event)
  if (!is.null(start)) risk_frame$.start <- start
  internal_names <- paste0(".x", seq_len(p))
  if (p) risk_frame[internal_names] <- as.data.frame(x)
  response <- if (is.null(start)) "survival::Surv(.stop, .event)" else
    "survival::Surv(.start, .stop, .event)"
  formula <- stats::as.formula(paste(response, "~", if (p)
    paste(internal_names, collapse = " + ") else "1"), env = asNamespace("survival"))
  backend_warnings <- character()
  fit <- withCallingHandlers(
    survival::coxph(formula, data = risk_frame, ties = controls$ties,
      singular.ok = TRUE, model = TRUE, x = TRUE, y = TRUE, nocenter = NULL,
      na.action = stats::na.fail,
      control = survival::coxph.control(iter.max = controls$iter.max,
                                        eps = controls$eps, timefix = FALSE)),
    warning = function(w) {
      backend_warnings <<- c(backend_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (length(backend_warnings) && any(grepl(
      "did not converge|failed to converge|infinite|ran out of iterations",
      backend_warnings, ignore.case = TRUE))) {
    # The counting-process backend reports "beta"; the right-censored
    # backend reports "coefficient". Only the converged-likelihood warning
    # diagnoses separation. Other warnings, including iteration exhaustion,
    # retain the usual error behavior.
    separation <- "^Loglik converged before variable [[:digit:], ]+; (beta|coefficient) may be infinite"
    if (identical(controls$separation_action, "baseline") &&
        all(grepl(separation, backend_warnings))) {
      model <- .cox_fit_numeric(stop, event, x[, FALSE, drop = FALSE], controls,
                                start, label, baseline_horizon)
      model$coefficients <- model$center <- zero
      model$diagnostics$coefficient_fit <- "separation_baseline_only"
      model$diagnostics$warnings <- unique(backend_warnings)
      return(model)
    }
    stop(label, " failed: ", paste(unique(backend_warnings), collapse = "; "),
         call. = FALSE)
  }
  # coxph deliberately suppresses its usual nonconvergence warning when
  # iter.max is one (a backend debugging facility). Reaching the requested
  # cap is therefore checked independently before accepting a fitted adapter.
  if (p && length(fit$iter) && any(fit$iter >= controls$iter.max)) {
    stop(label, " reached the iteration limit; convergence was not established.",
         call. = FALSE)
  }
  coefficients <- stats::coef(fit)
  if (is.null(coefficients)) coefficients <- numeric()
  singular <- is.na(coefficients)
  if (any(singular)) {
    diagnostics$singular_columns <- colnames(x)[singular]
    detail <- paste(diagnostics$singular_columns, collapse = ", ")
    if (controls$singular == "error") {
      stop(label, " has unidentifiable predictor columns: ", detail, ".", call. = FALSE)
    }
    warning(label, " has unidentifiable predictor columns; coefficients set to zero: ",
            detail, ".", call. = FALSE)
    coefficients[singular] <- 0
  }
  if (any(!is.finite(coefficients)) || any(!is.finite(fit$loglik))) {
    stop(label, " produced nonfinite estimates.", call. = FALSE)
  }
  if (length(backend_warnings)) {
    warning(label, ": ", paste(unique(backend_warnings), collapse = "; "), call. = FALSE)
  }
  baseline <- .cox_centered_baseline(stop, event, fit$linear.predictors,
                                      start, controls$ties, baseline_horizon)
  increments <- baseline$increments
  cumulative <- cumsum(increments)
  if (any(!is.finite(cumulative)))
    stop(label, " produced an unrepresentable centered cumulative hazard.", call. = FALSE)
  center <- if (p) fit$means else numeric()
  names(coefficients) <- names(center) <- colnames(x)
  diagnostics$warnings <- unique(backend_warnings)
  diagnostics$iterations <- fit$iter
  diagnostics$baseline_method <- "observed_time_log_risk_sweep"
  list(fit = fit, coefficients = coefficients, center = center,
       jump_times = baseline$time, increments = increments,
       cumulative_hazard = cumulative, diagnostics = diagnostics)
}

#' Extract Cox baseline jumps without altering observed event times
#'
#' @param stop,event,linear_predictors Aligned risk-record endpoints, event
#'   indicators, and centered fitted log relative risks.
#' @param start NULL for right-censored records, or counting-process entry times.
#' @param ties Breslow or Efron tie handling, matching the coefficient fit.
#' @param horizon Extract only event jumps at or before this time; coefficient
#'   fitting and the supplied risk records may extend beyond the horizon.
#' @returns Actual distinct event times and positive centered hazard increments.
#' @details
#' Risk membership is exactly start < t <= stop. An active-risk binary sum tree
#' stores log-sums of positive relative risks. Updating a leaf recomputes its
#' ancestors from their children; it never subtracts two large cumulative sums.
#' This avoids cancellation when a high-risk record leaves, and exponentiation
#' overflow when finite log risks span a wide range. Genuine unequal times are
#' never rounded, jittered, interpolated, or merged.
#'
#' The sweep costs O(n log n + J) after sorting n risk records for J event times,
#' with O(n + J) storage. Right-censored data, or intervals all entered before
#' the first queried event, use a linear log-tail-sum pass after sorting instead
#' of tree updates. Intervals ending between consecutive event times are skipped
#' when they never enter an evaluated risk set. Efron corrections visit tied event records
#' once. No information matrices or subject-by-time risk matrices are computed.
#'
#' This direct extraction also avoids the near-tie sentinel construction in
#' older survival::agsurv() versions, where half an event-time ulp added to the
#' last entry time can round back to that same entry time and corrupt risk sums.
#' @noRd
.cox_centered_baseline <- function(stop, event, linear_predictors, start = NULL,
                                    ties = "breslow", horizon = Inf) {
  event_rows <- which(event == 1L & stop <= horizon)
  if (!length(event_rows)) return(list(time = numeric(), increments = numeric()))
  n <- length(stop)
  if (is.null(linear_predictors)) linear_predictors <- numeric(n)
  if (length(linear_predictors) != n || any(!is.finite(linear_predictors)))
    stop("Cox risk-record linear predictors must be finite.", call. = FALSE)
  time <- sort(unique(stop[event_rows]))
  counts <- tabulate(match(stop[event_rows], time), nbins = length(time))
  exit_order <- order(stop)
  if (is.null(start) || all(start < time[1L])) {
    tail_log_risk <- rep(-Inf, n + 1L)
    for (i in seq.int(n, 1L)) {
      a <- linear_predictors[exit_order[i]]; b <- tail_log_risk[i + 1L]
      tail_log_risk[i] <- if (b == -Inf) a else
        max(a, b) + log1p(exp(-abs(a - b)))
    }
    log_risk <- tail_log_risk[findInterval(time, stop[exit_order], left.open = TRUE) + 1L]
  } else {
    leaves <- 2^ceiling(log2(n))
    tree <- rep(-Inf, 2 * leaves)
    exit_next <- 1L
    entry_order <- order(start)
    entry_next <- 1L
    log_risk <- numeric(length(time))
    for (j in seq_along(time)) {
      u <- time[j]
      while (exit_next <= n && stop[exit_order[exit_next]] < u) {
        node <- leaves + exit_order[exit_next] - 1L
        if (tree[node] != -Inf) {
          tree[node] <- -Inf
          while (node > 1) {
            node <- node %/% 2L
            a <- tree[2 * node]; b <- tree[2 * node + 1L]
            tree[node] <- if (a == -Inf) b else if (b == -Inf) a else
              max(a, b) + log1p(exp(-abs(a - b)))
          }
        }
        exit_next <- exit_next + 1L
      }
      while (entry_next <= n && start[entry_order[entry_next]] < u) {
        row <- entry_order[entry_next]
        if (stop[row] >= u) {
          node <- leaves + row - 1L
          tree[node] <- linear_predictors[row]
          while (node > 1) {
            node <- node %/% 2L
            a <- tree[2 * node]; b <- tree[2 * node + 1L]
            tree[node] <- if (a == -Inf) b else if (b == -Inf) a else
              max(a, b) + log1p(exp(-abs(a - b)))
          }
        }
        entry_next <- entry_next + 1L
      }
      log_risk[j] <- tree[1L]
    }
  }
  if (any(!is.finite(log_risk)))
    stop("A Cox event has no finite positive fitted risk set.", call. = FALSE)
  log_increment <- log(counts) - log_risk
  if (ties == "efron" && any(counts > 1L)) {
    event_rows <- event_rows[order(stop[event_rows])]
    last <- cumsum(counts)
    for (j in which(counts > 1L)) {
      d <- counts[j]
      values <- linear_predictors[event_rows[seq.int(last[j] - d + 1L, last[j])]]
      maximum <- max(values)
      event_log_risk <- maximum + log(sum(exp(values - maximum)))
      difference <- event_log_risk - log_risk[j]
      if (difference > 64 * .Machine$double.eps *
          max(1, abs(event_log_risk), abs(log_risk[j])))
        stop("Tied-event risk exceeds its Cox risk set.", call. = FALSE)
      fraction <- exp(min(0, difference))
      denominator <- (1 - fraction) + (d - seq.int(0L, d - 1L)) / d * fraction
      log_increment[j] <- log(sum(1 / denominator)) - log_risk[j]
    }
  }
  increments <- exp(log_increment)
  if (any(!is.finite(increments)) || any(increments <= 0))
    stop("Centered Cox baseline jumps are outside the representable range.", call. = FALSE)
  list(time = time, increments = increments)
}

#' Evaluate a saved centered Cox linear predictor without formula reconstruction
#' @noRd
.cox_linear_predictor <- function(model, x) {
  used <- model$coefficients != 0
  if (!any(used)) return(numeric(nrow(x)))
  lp <- as.vector(sweep(x[, used, drop = FALSE], 2L, model$center[used], "-") %*%
                    model$coefficients[used])
  if (any(!is.finite(lp))) {
    stop("New profile Cox linear predictors are nonfinite.", call. = FALSE)
  }
  lp
}

#' Evaluate the floored Cox survival law on the log scale
#' @noRd
.cox_log_survival <- function(cumulative, lp, floor) {
  result <- numeric(length(cumulative))
  use <- cumulative > 0
  if (any(use)) {
    lower <- log(floor)
    result[use] <- pmax(lower, -exp(pmin(log(cumulative[use]) + lp[use],
                                         log(-lower))))
  }
  result
}

#' Validate aligned censoring queries and permit explicit repeated profile indices
#' @noRd
.cox_censor_queries <- function(adapter, x, time, subject_index) {
  if (!inherits(adapter, c("life_censor_cox", "life_censor_forest"))) {
    stop("`adapter` must be a supported censoring model.", call. = FALSE)
  }
  .cox_check_design(x, adapter$n_predictors, adapter$predictor_names,
                     allow_empty_rows = TRUE)
  if (is.null(subject_index)) subject_index <- seq_len(nrow(x))
  if (!is.numeric(subject_index) || is.complex(subject_index) ||
      !is.null(dim(subject_index)) ||
      any(!is.finite(subject_index)) || any(subject_index != floor(subject_index) |
        subject_index < 1 | subject_index > nrow(x))) {
    stop("`subject_index` must be a numeric vector indexing rows of `x`.", call. = FALSE)
  }
  if (!is.numeric(time) || is.complex(time) || !is.null(dim(time)) ||
      any(!is.finite(time)) ||
      any(time < 0 | time > adapter$tau) ||
      !(length(time) %in% c(1L, length(subject_index)))) {
    stop("`time` must be a finite numeric scalar or vector with one value per queried subject, within [0, tau].",
         call. = FALSE)
  }
  if (length(time) == 1L) time <- rep(time, length(subject_index))
  list(subject_index = subject_index, time = time)
}

#' Build one cause-specific risk view from preindexed origin intervals
#' @noRd
.cox_transition_view <- function(intervals, rows, destination, x, tau) {
  selected <- intervals[rows, , drop = FALSE]
  event <- selected$event & !is.na(selected$to) & selected$to == destination
  list(data = data.frame(subject_index = selected$subject_index,
                         start = selected$start, stop = selected$stop,
                         event = as.integer(event), original_row = selected$original_row),
       x = x[selected$subject_index, , drop = FALSE])
}

#' Match the specified graph to empirical training transitions
#' @noRd
.cox_validate_edges <- function(edges, observed) {
  if (!is.data.frame(edges) || anyDuplicated(names(edges)) ||
      !all(c("from", "to") %in% names(edges))) {
    stop("`edges` must be a data frame containing unique from/to columns.", call. = FALSE)
  }
  edges <- edges[c("from", "to")]
  if (!all(vapply(edges, is.numeric, logical(1))) || anyNA(edges) ||
      anyDuplicated(edges) || nrow(edges) != nrow(observed) ||
      !setequal(paste(edges$from, edges$to), paste(observed$from, observed$to))) {
    stop("`edges` must contain exactly the distinct observed training transitions.",
         call. = FALSE)
  }
  rownames(edges) <- NULL
  edges
}
