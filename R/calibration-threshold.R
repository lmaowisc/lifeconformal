#' Resolve a first crossing from complete signed numerical sum bounds
#' @noRd
calibration_threshold_certificate <- function(lower_sum, upper_sum, rank) {
  if (!is.numeric(lower_sum) || !is.numeric(upper_sum) ||
      is.complex(lower_sum) || is.complex(upper_sum) ||
      !is.null(dim(lower_sum)) || !is.null(dim(upper_sum)) ||
      !length(lower_sum) || length(lower_sum) != length(upper_sum) ||
      anyNA(lower_sum) || anyNA(upper_sum) || any(!is.finite(lower_sum)) ||
      any(!is.finite(upper_sum)) || any(lower_sum > upper_sum) ||
      !is.numeric(rank) || is.complex(rank) || !is.null(dim(rank)) || length(rank) != 1L ||
      !is.finite(rank) || rank < 1 || rank != floor(rank))
    stop("Threshold certification requires finite ordered sum bounds and a positive integer rank.",
         call. = FALSE)
  first <- function(x) { at <- which(x); if (length(at)) at[1L] else Inf }
  possible <- first(upper_sum >= rank)
  certain <- first(lower_sum >= rank)
  list(possible_index = possible, certain_index = certain,
       threshold_stable = possible == certain,
       full_support = is.infinite(possible))
}

#' Evaluate one immutable recurrent profile for a threshold refinement pass
#' @noRd
evaluate_threshold_profile <- function(fit, profile, x, center, scale, q,
                                        controls, prepared, subject_index,
                                        use_marginal = FALSE) {
  times <- profile$snapshot$time
  h <- length(times)
  membership <- as.numeric(profile$membership)
  if (all(profile$resolved) && times[1L] >= max(profile$determination))
    return(list(probability = membership, lower = membership, upper = membership,
      bound_width = 0, queried_snapshots = 0L, structural_overrides = h * length(q),
      diagnostics = list()))
  certified <- outer(times, profile$determination, `>=`) &
    matrix(rep(profile$resolved, each = h), h)
  known <- matrix(rep(membership, each = h), h)
  known[!certified] <- NA_real_
  need <- which(rowSums(is.na(known)) > 0L)
  certain <- which(rowSums(is.na(known)) == 0L)
  fixed <- if (length(certain)) membership * sum(profile$weight[certain]) else
    numeric(length(q))
  if (length(need)) {
    evaluated <- with_life_grid_context(conditional_life_coverage(fit, x,
      profile$snapshot, need, center, scale, q, controls, prepared,
      subject_index, known[need, , drop = FALSE],
      query_weights = profile$weight[need], grid_refine = FALSE,
      use_marginal = use_marginal),
      subject_index, "calibration", profile$id)
    if (!is.null(evaluated$reduced)) {
      value <- evaluated$reduced$probability + fixed
      lower <- evaluated$reduced$lower + fixed
      upper <- evaluated$reduced$upper + fixed
    } else {
      # A history can become coherent for a monotone continuation after the
      # global batched monotone route was unavailable. Preserve exact overrides.
      p <- evaluated$probability; low <- evaluated$lower; high <- evaluated$upper
      k <- known[need, , drop = FALSE]; exact <- !is.na(k)
      p[exact] <- low[exact] <- high[exact] <- k[exact]
      weight <- profile$weight[need]
      negative <- weight < 0
      lower_matrix <- low; upper_matrix <- high
      lower_matrix[negative, ] <- high[negative, , drop = FALSE]
      upper_matrix[negative, ] <- low[negative, , drop = FALSE]
      value <- colSums(p * weight) + fixed
      lower <- colSums(lower_matrix * weight) + fixed
      upper <- colSums(upper_matrix * weight) + fixed
    }
  } else {
    evaluated <- list(bound_width = 0, diagnostics = list())
    value <- lower <- upper <- fixed
  }
  origin <- profile$resolved & profile$determination == 0
  value[origin] <- lower[origin] <- upper[origin] <- membership[origin]
  list(probability = value, lower = lower, upper = upper,
       bound_width = evaluated$bound_width, queried_snapshots = length(need),
       structural_overrides = sum(certified), diagnostics = evaluated$diagnostics)
}

#' Resolve calibration using signed bounds, refining only consequential profiles
#'
#' The bounds enclose discretization uncertainty; ordinary floating-point
#' arithmetic is not a directed-rounding interval enclosure. The candidate grid,
#' nuisance predictions, query snapshots, CDT certificates and score remain fixed.
#' @noRd
evaluate_threshold_augmented_curve <- function(fit, history, x, moments, cdt, q,
                                                controls, queries, rank,
                                                max_refinements = 7L) {
  n <- nrow(history$subjects)
  profiles <- lapply(seq_len(n), function(i) {
    first <- if (i == 1L) 1L else queries$endpoint_index[i - 1L] + 1L
    rows <- seq.int(first, queries$endpoint_index[i])
    list(snapshot = list(time = queries$time[rows], state = queries$snapshot$state[rows],
      accrued = queries$snapshot$accrued[rows, , drop = FALSE]),
      weight = queries$weight[rows], determination = cdt$time[i, ],
      resolved = cdt$resolved[i, ], membership = cdt$membership[i, ],
      id = history$subjects$id[i])
  })
  prepared <- queries$prepared
  # Risks and jump factors are shared across every refinement, never saved in fit.
  unresolved <- vapply(profiles, function(z)
    !all(z$resolved) || z$snapshot$time[1L] < max(z$determination), logical(1))
  if (any(unresolved))
    prepared$markov_jumps <- calibration_markov_jumps(fit, x, prepared)
  contributions <- lower <- upper <- matrix(0, n, length(q))
  widths <- numeric(n); levels <- evaluations <- queried <- overrides <- integer(n)
  details <- vector("list", n); blocked <- rep(FALSE, n); budget_errors <- list()
  marginal_active <- rep(FALSE, n)
  marginal_evaluations <- ordinary_evaluations <- scalar_attempts <- integer(n)
  marginal_fallback_attempts <- marginal_fallbacks <- integer(n)
  marginal_eligible <- marginal_threshold_eligible(fit, controls)
  evaluate <- function(i, level, use_marginal = FALSE) {
    local_controls <- controls
    local_controls$grid_size <- controls$grid_size * 2^level
    evaluate_threshold_profile(fit, profiles[[i]], x[i, , drop = FALSE],
      moments$mean[i, ], moments$scale[i, ], q, local_controls, prepared, i,
      use_marginal = use_marginal)
  }
  for (i in seq_len(n)) {
    value <- evaluate(i, 0L, marginal_eligible && unresolved[i])
    contributions[i, ] <- value$probability
    lower[i, ] <- value$lower; upper[i, ] <- value$upper
    widths[i] <- value$bound_width; queried[i] <- value$queried_snapshots
    overrides[i] <- value$structural_overrides; details[[i]] <- value$diagnostics
    evaluations[i] <- as.integer(value$queried_snapshots > 0L)
    marginal_active[i] <- isTRUE(value$diagnostics$marginal_applied)
    if (marginal_active[i]) marginal_evaluations[i] <- evaluations[i] else
      ordinary_evaluations[i] <- evaluations[i]
    scalar_attempts[i] <- if (is.null(value$diagnostics$marginal_scalar_attempts)) 0L else
      value$diagnostics$marginal_scalar_attempts
  }
  profile_stages <- function() {
    monotone <- vapply(details, function(z)
      identical(z$method, "deterministic_exact_monotone_clocks"), logical(1))
    ifelse(queried == 0L, "structural", ifelse(marginal_active, "marginal",
      ifelse(monotone, "monotone", "joint")))
  }
  profile_grid_sizes <- function(include_unqueried = FALSE) {
    result <- controls$grid_size * 2^levels
    if (!include_unqueried) result[queried == 0L] <- NA_real_
    if (marginal_eligible)
      result[profile_stages() %in% c("marginal", "monotone")] <- NA_real_
    result
  }
  stage_diagnostics <- function() {
    if (!marginal_eligible) return(list())
    list(marginal_active = marginal_active,
      marginal_evaluations = marginal_evaluations,
      marginal_scalar_attempts = scalar_attempts,
      applied_marginal_scalar_passes = 2L * sum(marginal_evaluations),
      ordinary_evaluations = ordinary_evaluations,
      marginal_fallback_attempts = marginal_fallback_attempts,
      marginal_fallbacks = marginal_fallbacks,
      profile_stages = profile_stages(),
      profile_marginal_grid_sizes = ifelse(marginal_active,
        controls$grid_size, NA_real_))
  }
  rounds <- 0L
  repeat {
    lower_sum <- colSums(lower); upper_sum <- colSums(upper)
    certificate <- calibration_threshold_certificate(lower_sum, upper_sum, rank)
    if (certificate$threshold_stable) break
    focus <- certificate$possible_index
    uncertainty <- upper[, focus] - lower[, focus]
    eligible <- which(uncertainty > 0 & (marginal_active | levels < max_refinements) & !blocked)
    if (!length(eligible)) {
      diagnostics <- c(certificate, list(rank = rank, q = q,
        curve_lower = lower_sum / n, curve_upper = upper_sum / n,
        lower_sum = lower_sum, upper_sum = upper_sum,
        initial_grid_size = controls$grid_size,
        grid_size = profile_grid_sizes(include_unqueried = TRUE),
        max_nodes = controls$max_nodes, max_refinements = max_refinements,
        refinements = levels, profile_diagnostics = details,
        budget_errors = budget_errors,
        requested_probability_tolerance = controls$probability_tolerance,
        max_probability_bound_width = max(widths),
        probability_tolerance_met = max(widths) <= controls$probability_tolerance),
        stage_diagnostics())
      stop(structure(list(message = paste0("The calibration threshold remains ambiguous within ",
        "the numerical discretization bounds. Increase control = list(grid_size = ...) ",
        "or max_nodes in life_calibrate(), or use calibration_accuracy = 'probability' ",
        "to request conditional-probability accuracy instead; that option ",
        "does not require a numerically resolved threshold."), call = NULL,
        diagnostics = diagnostics),
        class = c("lifeconformal_threshold_ambiguity", "error", "condition")))
    }
    # A bounded batch avoids re-summing every profile after each tiny update.
    batch <- eligible[order(uncertainty[eligible], decreasing = TRUE)]
    batch <- utils::head(batch, max(1L, ceiling(length(eligible) / 8)))
    for (i in batch) {
      from_marginal <- marginal_active[i]
      next_level <- if (from_marginal) levels[i] else levels[i] + 1L
      value <- tryCatch(evaluate(i, next_level),
        lifeconformal_grid_budget = function(e) e)
      evaluations[i] <- evaluations[i] + 1L
      ordinary_evaluations[i] <- ordinary_evaluations[i] + 1L
      if (from_marginal) marginal_fallback_attempts[i] <- marginal_fallback_attempts[i] + 1L
      if (inherits(value, "lifeconformal_grid_budget")) {
        blocked[i] <- TRUE
        budget_errors[[as.character(i)]] <- value
        next
      }
      levels[i] <- next_level
      if (from_marginal) marginal_fallbacks[i] <- marginal_fallbacks[i] + 1L
      marginal_active[i] <- FALSE
      contributions[i, ] <- value$probability
      lower[i, ] <- value$lower; upper[i, ] <- value$upper
      widths[i] <- value$bound_width; details[[i]] <- value$diagnostics
    }
    rounds <- rounds + 1L
  }
  point_sum <- colSums(contributions)
  point_crossing <- which(point_sum >= rank)
  point_index <- if (length(point_crossing)) point_crossing[1L] else Inf
  if (point_index != certificate$possible_index)
    stop("Floating-point aggregation disagrees with the resolved threshold bounds.", call. = FALSE)
  numerical <- c(certificate, list(accuracy_mode = "threshold",
    bound_scope = "numerical discretization, excluding floating-point rounding error",
    conditional_coverage_evaluated = any(queried > 0L),
    queried_snapshots = sum(queried), structural_overrides = sum(overrides),
    max_probability_bound_width = max(widths),
    requested_probability_tolerance = controls$probability_tolerance,
    probability_tolerance_met = max(widths) <= controls$probability_tolerance,
    curve_lower = lower_sum / n, curve_upper = upper_sum / n,
    curve_lower_sum = lower_sum, curve_upper_sum = upper_sum,
    curve_estimate_sum = point_sum,
    max_curve_bound_width = max((upper_sum - lower_sum) / n),
    initial_grid_size = controls$grid_size,
    profiles = n, evaluated_profiles = sum(queried > 0L),
    profile_grid_sizes = profile_grid_sizes(),
    profile_refinements = levels, profile_evaluations = sum(evaluations),
    refinement_rounds = rounds, profile_diagnostics = details,
    blocked_profiles = which(blocked), budget_errors = budget_errors,
    contribution_payload_bytes = 24 * n * length(q)), stage_diagnostics(), queries$diagnostics)
  list(contributions = contributions, diagnostics = numerical)
}
