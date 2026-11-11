#' Calibrate named marginal scores under an allocated joint error budget
#'
#' @param object Trained lifeconformal_fit with its full frozen score model.
#' @param alpha_components Numeric error budgets, named for every target.
#'   Names may arrive in any order and are aligned with the stored targets.
#' @param alpha Requested joint error budget. Allocations must be positive and
#'   their sum must not exceed alpha, apart from floating-point summation error.
#' @param calibration_args Named arguments forwarded to life_calibrate(), such
#'   as data, vector inputs, method, q and control. Do not include object,
#'   alpha or alpha_components.
#' @returns A lifeconformal_calibration retaining the original fitted object
#'   and named marginal thresholds. The full_support flag is TRUE only when
#'   every component fell back; component_full_support and
#'   any_component_full_support retain partial fallback information. The curve
#'   has target/q/estimate columns. Retained contributions are a named list of
#'   subject-by-threshold matrices, rather than a sum across target components.
#' @details Each temporary problem selects one target while reusing all fitted
#'   outcome, censoring and augmentation models, original preprocessing,
#'   horizon and scale floor. Thus it calibrates the original marginal
#'   standardized deviations; it does not fit new scores. Flexible adapters
#'   retain their original multivariate history features and model problem.
#'   Their common prediction interface selects requested target outputs only
#'   after constructing full model histories.
#'
#' Parsing, frozen baseline moments, censoring curves and signed endpoint
#'   queries are shared within one call. Each target retains its own score,
#'   determination times, coverage curve and allocated cutoff. No cache is
#'   retained across calls. Flexible AIPW components also share continuation
#'   draws in bounded query chunks; signed sums retain their original order,
#'   and full contribution matrices are allocated only when requested on this
#'   path. Marginal regions are intersected with the original joint
#'   occupation support at prediction. The union bound uses the sum of the
#'   allocated errors under each marginal method's applicable conditions.
#' @noRd
component_calibration <- function(object, alpha_components, alpha,
                                   calibration_args = list()) {
  if (!inherits(object, "lifeconformal_fit")) {
    stop("object must be returned by life_train().", call. = FALSE)
  }
  if (!identical(object$problem$score, "standardized_max")) {
    stop("Component allocation requires the frozen standardized score family.",
         call. = FALSE)
  }
  target_names <- colnames(object$problem$rewards)
  budgets <- validate_component_allocation(alpha_components, alpha, target_names)
  if (!is.list(calibration_args) || (length(calibration_args) &&
      (is.null(names(calibration_args)) || anyNA(names(calibration_args)) ||
       any(!nzchar(names(calibration_args))) ||
       anyDuplicated(names(calibration_args))))) {
    stop("calibration_args must be a list with unique, nonempty names.", call. = FALSE)
  }
  if (any(names(calibration_args) %in% c("object", "alpha", "alpha_components"))) {
    stop("calibration_args must not replace object, alpha or alpha_components.",
         call. = FALSE)
  }
  unknown <- setdiff(names(calibration_args), names(formals(life_calibrate)))
  if (length(unknown)) stop("Unknown calibration argument(s): ",
                            paste(unknown, collapse = ", "), ".", call. = FALSE)
  method <- if (is.null(calibration_args$method)) "aipw" else
    match.arg(calibration_args$method, c("aipw", "ipw", "ipw_full", "complete"))
  controls <- resolve_calibration_controls(object,
    if (is.null(calibration_args$control)) list() else calibration_args$control)
  context_args <- intersect(names(calibration_args),
                            names(formals(prepare_life_calibration_context)))
  context <- do.call(prepare_life_calibration_context,
                     c(list(object = object), calibration_args[context_args]))
  scalar_objects <- lapply(seq_along(target_names), function(j)
    subset_component_problem(object, j))
  grids <- lapply(scalar_objects, resolve_calibration_grid, q = calibration_args$q)
  n <- nrow(context$history$subjects)
  active <- which(ceiling((n + 1) * (1 - budgets)) <= n)
  if (length(active)) {
    moment_object <- if (length(active) == length(target_names)) object else
      subset_component_problem(object, active)
    context <- prepare_life_calibration_predictions(context, object, method,
                                                     moment_object)
  }
  shared <- NULL
  if (method == "aipw" && object$augmentation$type != "markov_cox" &&
      length(active) > 1L && controls$calibration_accuracy != "threshold") {
    precision <- calibration_sum_precision()
    if (!is.null(precision)) shared <- evaluate_flexible_component_curves(
      object, context, scalar_objects, grids, active, controls, precision)
  }
  results <- vector("list", length(target_names))
  names(results) <- target_names
  for (j in seq_along(target_names)) {
    scalar <- scalar_objects[[j]]
    component_context <- subset_life_calibration_context(context, scalar)
    results[[j]] <- calibrate_life_prepared(scalar, component_context,
      unname(budgets[j]), method, grids[[j]], controls, match.call(),
      evaluated = if (is.null(shared)) NULL else shared[[j]])
  }
  first <- results[[1L]]
  extract_number <- function(field) {
    setNames(vapply(results, function(z) as.numeric(z[[field]]), numeric(1)),
              target_names)
  }
  full_support <- setNames(vapply(results, function(z) z$full_support, logical(1)),
                            target_names)
  reasons <- setNames(vapply(results, function(z) z$reason, character(1)), target_names)
  curves <- lapply(seq_along(results), function(j) {
    data.frame(target = rep(target_names[j], nrow(results[[j]]$curve)),
               results[[j]]$curve, row.names = NULL)
  })
  curve <- do.call(rbind, curves)
  rownames(curve) <- NULL
  component_diagnostics <- lapply(results, `[[`, "diagnostics")
  diagnostics <- list(
    overlap_checked = first$diagnostics$overlap_checked,
    overlapping_labels = first$diagnostics$overlapping_labels,
    allocated_alpha = sum(budgets),
    components = component_diagnostics,
    numerical = list(components = lapply(component_diagnostics, `[[`, "numerical")),
    preparation = "one shared history and nuisance preparation; target-specific calibration curves"
  )
  structure(list(fit = object, alpha = alpha, alpha_components = budgets,
                 allocated_alpha = sum(budgets), method = first$method,
                 n = first$n, threshold = extract_number("threshold"),
                 full_support = all(full_support),
                 component_full_support = full_support,
                 any_component_full_support = any(full_support),
                 reason = "component_allocation", component_reason = reasons,
                 cutoff = extract_number("cutoff"), rank = extract_number("rank"),
                 curve = curve, diagnostics = diagnostics, control = first$control,
                 contributions = if (first$control$keep_contributions)
                   lapply(results, `[[`, "contributions") else NULL,
                 call = match.call(), object_version = 1L),
            class = "lifeconformal_calibration")
}

#' Align and validate a complete marginal allocation of the joint error budget
#' @noRd
validate_component_allocation <- function(alpha_components, alpha, target_names) {
  if (!is.numeric(alpha) || is.complex(alpha) || length(alpha) != 1L ||
      is.na(alpha) || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be strictly between zero and one.", call. = FALSE)
  }
  if (!is.numeric(alpha_components) || is.complex(alpha_components) ||
      !is.null(dim(alpha_components)) || !length(alpha_components) ||
      is.null(names(alpha_components)) || anyNA(names(alpha_components)) ||
      any(!nzchar(names(alpha_components))) || anyDuplicated(names(alpha_components))) {
    stop("alpha_components must be a numeric vector with unique, nonempty target names.",
         call. = FALSE)
  }
  if (!setequal(names(alpha_components), target_names)) {
    stop("alpha_components must name every target exactly once: ",
         paste(target_names, collapse = ", "), ".", call. = FALSE)
  }
  if (anyNA(alpha_components) || any(!is.finite(alpha_components)) ||
      any(alpha_components <= 0 | alpha_components >= 1)) {
    stop("Each alpha_components allocation must lie strictly between zero and one.",
         call. = FALSE)
  }
  # Allow only summation-scale rounding, not a user-level overspend tolerance.
  tolerance <- 8 * .Machine$double.eps * alpha
  if (sum(alpha_components) > alpha + tolerance) {
    stop("The sum of alpha_components must not exceed alpha.", call. = FALSE)
  }
  alpha_components[target_names]
}

#' Select a marginal problem while preserving every fitted model and feature
#' @noRd
subset_component_problem <- function(object, target_index) {
  problem <- object$problem
  targets <- compile_targets(problem$targets$sets[target_index],
                              problem$states, problem$tau)
  support <- problem$structural_support
  support$targets <- targets
  support$rewards <- targets$rewards
  for (name in c("reward_min", "reward_max", "lower", "upper"))
    support[[name]] <- support[[name]][, target_index, drop = FALSE]
  support$joint_support <- targets$support
  problem$targets <- targets
  problem$A <- targets$A
  problem$rewards <- targets$rewards
  if (!is.null(problem$score_weights)) {
    problem$score_weights <- problem$score_weights[targets$names]
  }
  problem$structural_support <- support
  object$problem <- problem
  # Internal orchestration uses fitted model helpers. Remove the original
  # joint callable wrappers from this temporary object to prevent accidental
  # evaluation against the wrong target definition. The returned calibration
  # object retains the untouched original fit and its callable functions.
  object$functions <- NULL
  object
}

#' Resolve calibration controls without changing the fitted score or censoring rule
#' @noRd
resolve_calibration_controls <- function(object, control) {
  # Validate the override list before merging, so duplicate or unnamed entries
  # cannot disappear through list assignment.
  validated <- resolve_life_controls(control)
  control <- validated[names(control)]
  if (any(names(control) %in% c("scale_floor", "survival_floor")))
    stop("Calibration cannot change trained scale_floor or survival_floor.", call. = FALSE)
  controls <- object$control
  controls[names(control)] <- control
  # Choose the numerical task before choosing its starting resolution.
  # Existing saved fits retain their explicit mode; only 'auto' dispatches.
  if (identical(controls$calibration_accuracy, "auto"))
    controls$calibration_accuracy <- if (object$augmentation$type == "markov_cox")
      "threshold" else "probability"
  trained_grid <- "grid_size" %in% object$control_overrides ||
    (!"control_overrides" %in% names(object) && !identical(object$control$grid_size, 128L))
  if (identical(controls$calibration_accuracy, "threshold") &&
      !"grid_size" %in% names(control) && !trained_grid) controls$grid_size <- 16L
  resolve_life_controls(controls)
}

#' Parse an independent calibration sample once, without evaluating any model
#' @noRd
prepare_life_calibration_context <- function(object, data = NULL, id = NULL,
    initial_state = NULL, check_overlap = TRUE, x = NULL, response = NULL,
    time = NULL, status = NULL, start = NULL, stop = NULL, from = NULL,
    to = NULL, event = NULL) {
  if (!is.logical(check_overlap) || length(check_overlap) != 1L ||
      is.na(check_overlap))
    base::stop("check_overlap must be TRUE or FALSE.", call. = FALSE)
  formula <- if (!is.null(data)) object$schema$formula else NULL
  response <- resolve_public_response(formula, response, time, status, start,
                                       stop, from, to, event)
  input <- resolve_input(data = data, x = x, id = id, response = response,
                          schema = object$schema)
  history <- normalize_histories(input,
    initial_state = stored_initial_state(object$schema, initial_state),
    status_map = object$schema$status_map,
    absorbing = object$problem$absorbing_labels, states = object$problem$states)
  # Normalization has already validated every path. Only the frozen problem
  # boundary and declared structural restrictions remain to be checked here.
  check_problem_history(history, object$problem$states)
  check_structural_transitions(history, object$problem$structural_support)
  overlap <- intersect(as.character(history$subjects$id),
                         as.character(object$training_ids))
  if (check_overlap && length(overlap))
    base::stop("Calibration identifiers overlap training: ",
               paste(utils::head(overlap), collapse = ", "),
               ". Supply an independent sample or set check_overlap = FALSE only for independently numbered IDs.",
               call. = FALSE)
  list(history = history,
       diagnostics = list(overlap_checked = check_overlap,
                          overlapping_labels = length(overlap), numerical = list()),
       moments = NULL, outcome_view = NULL, censor_view = NULL,
       queries = NULL, censor_curves = NULL)
}

#' Evaluate target-independent calibration work once within the current call
#' @noRd
prepare_life_calibration_predictions <- function(context, object, method,
                                                   moment_object = object) {
  history <- context$history
  context$outcome_view <- apply_preprocessor(history$baseline,
                                               object$preprocessor, "outcome")
  context$moments <- predict_life_moments(moment_object, history$baseline,
    object$problem$states[history$subjects$initial_state],
    outcome_view = context$outcome_view)
  if (method == "complete") return(context)
  context$censor_view <- apply_preprocessor(history$baseline,
                                            object$preprocessor, "censor")
  if (method == "aipw") {
    # Endpoint and inverse-survival increments depend on the observed history,
    # not the target, threshold, or first coverage-determination time.
    absorption <- history$subjects$absorption_time
    absorption[is.na(absorption)] <- Inf
    context$queries <- prepare_calibration_queries(object, history,
      context$censor_view$x,
      list(endpoint = pmin(history$subjects$observation_end, absorption,
                            object$problem$tau)))
  } else if (object$censor$type == "survival_forest") {
    context$censor_curves <- censor_forest_curves(object$censor,
                                                   context$censor_view$x)
  }
  context
}

#' Select target coordinates without mutating shared history or model features
#' @noRd
subset_life_calibration_context <- function(context, object) {
  targets <- colnames(object$problem$rewards)
  if (!is.null(context$moments)) {
    selected <- match(targets, colnames(context$moments$mean))
    if (anyNA(selected)) {
      # Only rank-infeasible components lack moments; their core returns before
      # prediction. They must not inherit an unrelated target's moments.
      context$moments <- NULL
    } else for (name in c("mean", "variance", "second_moment", "scale")) {
      if (!is.null(context$moments[[name]]))
        context$moments[[name]] <- context$moments[[name]][, selected, drop = FALSE]
    }
  }
  if (!is.null(context$queries)) {
    prepared <- context$queries$prepared
    selected <- match(targets, colnames(prepared$rewards))
    context$queries$snapshot$accrued <-
      context$queries$snapshot$accrued[, selected, drop = FALSE]
    context$queries$prepared$prefix <- prepared$prefix[, selected, drop = FALSE]
    context$queries$prepared$rewards <- prepared$rewards[, selected, drop = FALSE]
    if (object$augmentation$type != "markov_cox" &&
        !identical(targets, colnames(object$augmentation$problem$rewards)))
      context$queries$prepared$flexible <- prepared
  }
  context
}

#' Resolve the trained, outcome-independent candidate threshold grid
#' @noRd
resolve_calibration_grid <- function(object, q) {
  if (is.null(q)) {
    weights <- object$problem$score_weights
    maximum <- 1 / (object$control$scale_floor *
                      if (is.null(weights)) 1 else min(weights))
    if (!is.finite(maximum) || maximum > 1e6)
      stop("The full-support score bound is too large; supply an explicit q grid.",
           call. = FALSE)
    q <- unique(c(seq(0, 4, by = 0.05), seq(4.1, 10, by = 0.1),
                  seq(11, max(11, ceiling(maximum)), length.out = 90L)))
  }
  if (!is.numeric(q) || is.complex(q) || !is.null(dim(q)) || !length(q) ||
      anyNA(q) || any(!is.finite(q)) || any(q < 0) || any(diff(q) <= 0))
    stop("q must be finite, nonnegative and strictly increasing.", call. = FALSE)
  q
}

#' Calibrate one frozen score using an immutable, within-call sample context
#' @param evaluated Optional within-call flexible component evaluation with
#'   completed raw sums, means, determination arrays and diagnostics.
#' @noRd
calibrate_life_prepared <- function(object, context, alpha, method, q, controls,
                                     call = NULL, evaluated = NULL) {
  q <- resolve_calibration_grid(object, q)
  history <- context$history
  n <- nrow(history$subjects)
  diagnostics <- context$diagnostics
  contributions <- NULL
  full_history <- !is.na(history$subjects$absorption_time) |
    history$subjects$observation_end >= object$problem$tau
  if (method == "complete" && !all(full_history))
    stop("Complete calibration requires follow-up through absorption or the horizon for every subject.", call. = FALSE)
  rank <- ceiling((as.double(n) + 1) * (1 - alpha))
  if (rank > n) {
    inversion <- invert_coverage_curve(numeric(), numeric(), n, alpha)
    diagnostics$skipped <- "rank_requires_full_support"
    diagnostics$numerical <- list(conditional_coverage_evaluated = FALSE)
  } else {
    if (identical(controls$calibration_accuracy, "threshold") && method == "aipw" &&
        object$augmentation$type != "markov_cox")
      stop("calibration_accuracy = 'threshold' requires Markov augmentation; flexible Monte Carlo continuation does not provide deterministic numerical bounds.",
           call. = FALSE)
    if (is.null(context$moments))
      context <- prepare_life_calibration_predictions(context, object, method)
    moments <- context$moments
    cdt <- if (is.null(evaluated))
      compute_cdt(history, object$problem, moments$mean, moments$scale, q) else
        evaluated$cdt
    if (method == "ipw_full") {
      cdt$time[] <- rep(cdt$endpoint, length(q))
      cdt$resolved[] <- rep(full_history, length(q))
      cdt$membership[!cdt$resolved] <- NA
    }
    if (method == "complete") {
      contributions <- matrix(as.double(cdt$membership), n, length(q))
      diagnostics$numerical <- list(conditional_coverage_evaluated = FALSE)
    } else if (method != "aipw") {
      times <- as.vector(cdt$time)
      indices <- rep(seq_len(n), times = length(q))
      if (is.null(context$censor_curves)) {
        g <- predict_life_censor(object$censor, context$censor_view$x,
                                  times, "left", indices)
      } else {
        curves <- context$censor_curves
        knot <- findInterval(times, curves$jump_times, left.open = TRUE)
        g <- rep(1, length(times))
        use <- knot > 0L
        g[use] <- curves$survival[cbind(indices[use], knot[use])]
      }
      contributions <- matrix(compute_ipw_contributions(
        as.vector(cdt$membership), as.vector(cdt$resolved), g), n, length(q))
      diagnostics$numerical <- list(conditional_coverage_evaluated = FALSE)
    } else {
      if (is.null(evaluated)) evaluated <- evaluate_augmented_curve(object, history,
          context$outcome_view$x, context$censor_view$x, moments, cdt, q, controls,
          queries = context$queries, rank = rank)
      contributions <- evaluated$contributions
      diagnostics$numerical <- evaluated$diagnostics
    }
    estimate_sum <- if (is.null(evaluated$estimate_sum)) colSums(contributions) else
      evaluated$estimate_sum
    # Threshold estimates and bounds display the same raw-sum normalization.
    # A separately rounded mean can otherwise lie one ULP outside its bounds.
    estimate <- if (identical(controls$calibration_accuracy, "threshold"))
      estimate_sum / n else if (is.null(evaluated$estimate)) colMeans(contributions) else
        evaluated$estimate
    diagnostics$numerical$curve_estimate_sum <- estimate_sum
    if (is.null(diagnostics$numerical$curve_lower_sum) &&
        (method != "aipw" || object$augmentation$type == "markov_cox")) {
      diagnostics$numerical$curve_lower_sum <- estimate_sum
      diagnostics$numerical$curve_upper_sum <- estimate_sum
      diagnostics$numerical$curve_lower <- estimate
      diagnostics$numerical$curve_upper <- estimate
      diagnostics$numerical$max_curve_bound_width <- 0
    }
    if (identical(controls$calibration_accuracy, "threshold") &&
        !identical(diagnostics$numerical$accuracy_mode, "threshold")) {
      # IPW and exact monotone Markov evaluations also use the integer-rank
      # comparison, without a ratio comparison at a rounded tie.
      total <- estimate_sum
      certificate <- calibration_threshold_certificate(total, total, rank)
      diagnostics$numerical[names(certificate)] <- certificate
      diagnostics$numerical$accuracy_mode <- "threshold"
      diagnostics$numerical$curve_lower <- total / n
      diagnostics$numerical$curve_upper <- total / n
      diagnostics$numerical$curve_estimate_sum <- total
      diagnostics$numerical$curve_lower_sum <- total
      diagnostics$numerical$curve_upper_sum <- total
      diagnostics$numerical$max_curve_bound_width <- 0
      diagnostics$numerical$bound_scope <-
        "numerical discretization, excluding floating-point rounding error"
    }
    inversion <- invert_coverage_curve(q, estimate, n, alpha,
      estimate_sum = estimate_sum)
    if (!is.null(diagnostics$numerical[["curve_lower"]])) {
      inversion$curve$numerical_lower <- diagnostics$numerical$curve_lower
      inversion$curve$numerical_upper <- diagnostics$numerical$curve_upper
      certificate <- calibration_threshold_certificate(
        diagnostics$numerical$curve_lower_sum,
        diagnostics$numerical$curve_upper_sum, rank)
      point_index <- if (inversion$full_support) Inf else inversion$index
      certificate$threshold_stable <- certificate$threshold_stable &&
        certificate$possible_index == point_index
      diagnostics$numerical[names(certificate)] <- certificate
    }
    diagnostics$determination <- cdt$diagnostics
    diagnostics$resolved_fraction <- colMeans(cdt$resolved)
  }
  if (is.null(inversion$curve$numerical_lower)) {
    exact <- method != "aipw" || object$augmentation$type == "markov_cox"
    bound <- if (exact) inversion$curve$estimate else
      rep(NA_real_, nrow(inversion$curve))
    inversion$curve$numerical_lower <- bound
    inversion$curve$numerical_upper <- bound
    diagnostics$numerical$max_curve_bound_width <- if (exact) 0 else NA_real_
    diagnostics$numerical$threshold_stable <- if (exact) TRUE else NA
  }
  diagnostics$numerical$accuracy_mode <- controls$calibration_accuracy
  if (is.null(diagnostics$numerical$probability_tolerance_met)) {
    width <- diagnostics$numerical$max_probability_bound_width
    diagnostics$numerical$probability_tolerance_met <- if (is.null(width))
      if (method != "aipw" || object$augmentation$type == "markov_cox") TRUE else NA else
        width <= controls$probability_tolerance
  }
  structure(list(fit = object, alpha = alpha, method = method, n = n,
    threshold = inversion$threshold, full_support = inversion$full_support,
    reason = inversion$reason, cutoff = inversion$criterion, rank = inversion$rank,
    curve = inversion$curve, diagnostics = diagnostics, control = controls,
    contributions = if (controls$keep_contributions) contributions else NULL,
    call = call, object_version = 1L), class = "lifeconformal_calibration")
}
