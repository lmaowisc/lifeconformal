#' Extract one profile's factored Markov kernel inputs without jump matrices
#' @noRd
life_markov_kernel_factors <- function(jumps, profile, absorbing) {
  list(base_increments = jumps$base_increments,
       relative_risk = as.double(jumps$relative_risk[profile, ]),
       from = as.double(jumps$edges$from), to = as.double(jumps$edges$to),
       absorbing = as.double(absorbing))
}

#' Reuse immutable augmentation factors across marginal calibration curves
#'
#' The cache belongs to one prepared calibration sample and is never attached
#' to a fitted or returned object. Target subsets share the same augmentation
#' adapter and baseline design matrix, so their jump factors are identical.
#' @noRd
calibration_markov_jumps <- function(fit, x, prepared) {
  cache <- prepared$augmentation_cache
  if (is.environment(cache)) {
    if (!exists("markov_jumps", envir = cache, inherits = FALSE))
      cache$markov_jumps <- predict_transition_jumps(fit$augmentation, x)
    return(cache$markov_jumps)
  }
  predict_transition_jumps(fit$augmentation, x)
}

#' Group exactly equal numeric profiles without rounding floating-point values
#' @noRd
life_exact_profile_groups <- function(values) {
  # Match each coordinate before constructing integer keys. Decimal formatting
  # of the original doubles could incorrectly merge distinct fitted risks.
  key <- rep("", nrow(values))
  for (j in seq_len(ncol(values)))
    key <- paste0(key, ":", match(values[, j], values[, j]))
  group <- match(key, unique(key))
  list(group = group, first = which(!duplicated(group)))
}

#' Attach a compact public profile context to typed numerical failures
#' @noRd
with_life_grid_context <- function(code, profile_index, stage, subject_id = NULL) {
  annotate <- function(error) {
    error$diagnostics$profile_index <- profile_index
    error$diagnostics$stage <- stage
    label <- paste(if (stage == "calibration") "Calibration subject" else
      "Prediction profile", profile_index)
    if (!is.null(subject_id) && length(subject_id) == 1L && !is.na(subject_id)) {
      error$diagnostics$subject_id <- subject_id
      text <- encodeString(as.character(subject_id), quote = "'")
      if (nchar(text) > 84L) text <- paste0(substr(text, 1L, 80L), "...'")
      label <- paste0(label, " (", text, ")")
    }
    error$message <- paste0(label, ": ", conditionMessage(error))
    stop(error)
  }
  tryCatch(force(code), lifeconformal_grid_precision = annotate,
             lifeconformal_grid_budget = annotate)
}

#' Prepare a common-endpoint signed query representation once per calibration
#' @noRd
prepare_calibration_queries <- function(fit, history, censor_x, cdt) {
  n <- nrow(history$subjects)
  jumps <- fit$censor$jump_times
  endpoint <- cdt$endpoint
  count <- findInterval(endpoint, jumps)
  equal <- count > 0L
  equal[equal] <- jumps[count[equal]] == endpoint[equal]
  count <- count - equal
  curves <- if (fit$censor$type == "survival_forest") {
    censor_forest_curves(fit$censor, censor_x)$survival
  } else if (length(jumps)) {
    matrix(predict_life_censor(fit$censor, censor_x, rep(jumps, each = n),
             "right", rep(seq_len(n), length(jumps))), n, length(jumps))
  } else matrix(numeric(), n, 0L)
  # Endpoints use the same frozen step curves as the augmentation increments.
  # Looking up the last strictly earlier knot also avoids a second forest
  # prediction, which would recompute the entire survival matrix.
  g_end <- rep(1, n)
  after_jump <- count > 0L
  g_end[after_jump] <- curves[cbind(which(after_jump), count[after_jump])]
  times <- weights <- vector("list", n)
  for (i in seq_len(n)) {
    j <- seq_len(count[i])
    times[[i]] <- c(jumps[j], endpoint[i])
    if (length(j)) {
      gr <- curves[i, j]
      gl <- c(1, gr[-length(gr)])
      weights[[i]] <- c(-(gl - gr) / gl / gr, 1 / g_end[i])
    } else weights[[i]] <- 1 / g_end[i]
  }
  subject_index <- rep(seq_len(n), count + 1L)
  times <- unlist(times, use.names = FALSE)
  weights <- unlist(weights, use.names = FALSE)
  nonzero <- weights != 0
  subject_index <- subject_index[nonzero]
  times <- times[nonzero]
  weights <- weights[nonzero]
  prepared <- prepare_occupation_history(history, fit$problem)
  prepared$augmentation_cache <- new.env(parent = emptyenv())
  snapshots <- evaluate_occupation_history(prepared, subject_index, times,
                                            side = "right")
  list(subject_index = subject_index, time = times, weight = weights,
       snapshot = snapshots, prepared = prepared,
       endpoint_index = cumsum(tabulate(subject_index, nbins = n)),
       diagnostics = list(n_queries = length(times),
         max_endpoint_weight = max(1 / g_end),
         floored_endpoint_fraction = mean(g_end <= fit$control$survival_floor)))
}

#' Check whether observed histories support exact monotone-clock batching
#' @noRd
observed_monotone_histories <- function(history, rewards) {
  rows <- which(history$intervals$event)
  if (!length(rows)) return(TRUE)
  all(rewards[history$intervals$to[rows], , drop = FALSE] <=
        rewards[history$intervals$from[rows], , drop = FALSE])
}

#' Evaluate the exact monotone-target augmented curve in one compiled pass
#' @noRd
evaluate_monotone_augmented_curve <- function(fit, history, x, censor_x,
                                               moments, cdt, q, orientation = NULL,
                                               queries = NULL) {
  if (is.null(queries)) queries <- prepare_calibration_queries(fit, history, censor_x, cdt)
  snapshot <- queries$snapshot
  if (is.null(orientation))
    orientation <- orient_markov_rewards(fit$problem$rewards, fit$augmentation$edges)
  transformed <- orient_markov_query(orientation, moments$mean, snapshot$accrued,
                                      queries$time, fit$problem$tau)
  active <- transformed$rewards[snapshot$state, , drop = FALSE] == 1
  # This equality follows from the verified nonincreasing observed reward
  # paths: a target still active at t has accrued exactly t since the origin.
  # Canonicalization removes only floating summation roundoff of that identity.
  for (j in seq_len(ncol(snapshot$accrued))) {
    transformed$accrued[active[, j], j] <- queries$time[active[, j]]
    if (any(transformed$flip) && !transformed$flip[j])
      transformed$original_accrued[active[, j], j] <- queries$time[active[, j]]
  }
  if (!any(transformed$flip)) transformed$original_accrued <- transformed$accrued
  last_certificate <- apply(cdt$time, 1L, max)
  last_certificate[rowSums(cdt$resolved) < length(q)] <- Inf
  queried_snapshots <- sum(queries$time < last_certificate[queries$subject_index])
  if (!queried_snapshots) {
    contribution <- matrix(as.numeric(cdt$membership), nrow(cdt$membership), length(q))
    return(list(contributions = contribution,
      diagnostics = c(list(conditional_coverage_evaluated = FALSE,
                            method = "structurally_resolved",
                            queried_snapshots = 0L,
                            max_probability_bound_width = 0), queries$diagnostics)))
  }
  jumps <- calibration_markov_jumps(fit, x, queries$prepared)
  evaluated <- markov_monotone_batch(jumps$jump_times, jumps$base_increments,
    jumps$edges, jumps$relative_risk, transformed$rewards, fit$problem$tau,
    transformed$center, moments$scale, q, queries$subject_index, queries$time,
    snapshot$state, transformed$accrued, queries$weight,
    determination_time = cdt$time, resolved = cdt$resolved,
    membership = cdt$membership, flip = transformed$flip,
    original_center = transformed$original_center,
    original_accrued = transformed$original_accrued)
  origin <- cdt$resolved & cdt$time == 0
  evaluated$contributions[origin] <- as.numeric(cdt$membership[origin])
  evaluated$diagnostics$origin_exact_overrides <- sum(origin)
  evaluated$diagnostics$complemented_targets <- names(transformed$flip)[transformed$flip]
  list(contributions = evaluated$contributions,
       diagnostics = c(list(conditional_coverage_evaluated = TRUE,
                            method = "exact_batched_monotone_rewards",
                            queried_snapshots = queried_snapshots,
                            max_probability_bound_width = 0),
                       queries$diagnostics, evaluated$diagnostics))
}
