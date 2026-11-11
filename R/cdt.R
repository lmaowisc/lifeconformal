#' Prepare structural coverage certification
#'
#' @param history Validated observed histories.
#' @param problem Compiled target and structural-support specification.
#' @param center,scale Subject-by-target matrices defining the frozen maximum
#'   standardized absolute-deviation score. Scales must be strictly positive.
#' @returns A prepared object referencing the histories and fixed score arrays.
#' @details The reachable-state coordinate envelope is conservative for joint
#'   targets. Its certificates are valid stopping certificates, but may be later
#'   than the earliest determination available from exact joint constraints.
#'   No fitted transition probabilities enter the calculation.
#' @noRd
prepare_cdt_sweep <- function(history, problem, center, scale) {
  if (!inherits(history, "life_history_data") || !inherits(problem, "life_problem")) {
    stop("CDT computation requires validated histories and a compiled problem.",
         call. = FALSE)
  }
  if (!identical(as.character(history$schema$state_labels),
                 as.character(problem$states))) {
    stop("CDT history and problem state dictionaries differ.", call. = FALSE)
  }
  n <- nrow(history$subjects)
  d <- ncol(problem$rewards)
  check_array <- function(x, name) {
    if (is.numeric(x) && is.null(dim(x)) && d == 1L && length(x) == n) {
      x <- matrix(x, n, 1L)
    }
    if (!is.matrix(x) || !is.numeric(x) || is.complex(x) ||
        !identical(dim(x), c(n, d)) || any(!is.finite(x))) {
      stop(name, " must be a finite subject-by-target numeric matrix.", call. = FALSE)
    }
    if (!is.null(colnames(x)) &&
        !identical(colnames(x), colnames(problem$rewards))) {
      stop(name, " target columns differ from the compiled target order.", call. = FALSE)
    }
    x
  }
  center <- check_array(center, "center")
  scale <- check_array(scale, "scale")
  if (any(scale <= 0)) stop("CDT scales must be strictly positive.", call. = FALSE)
  structure(list(history = history, problem = problem, center = center,
                 scale = scale, n = n, d = d), class = "life_cdt_prepared")
}

#' Compute an outer occupation envelope and score bounds at one snapshot
#' @noRd
cdt_box_bounds <- function(accrued, state, time, problem, center, scale,
                           endpoint_open = TRUE) {
  remain <- max(0, problem$tau - time)
  minimum <- problem$structural_support$reward_min[state, ]
  maximum <- problem$structural_support$reward_max[state, ]
  reward <- problem$rewards[state, ]
  # Form the complete binary occupation totals before score subtraction.
  # Rounding tau-time first can move an inclusive boundary by one ulp even
  # when the accumulated and remaining occupation cover the whole horizon.
  lower <- markov_protected_total(problem$tau, time, minimum, accrued)
  upper <- markov_protected_total(problem$tau, time, maximum, accrued)
  # A right-continuous finite-transition path spends a positive initial
  # sojourn in its current state. An extreme requiring a different reward is
  # therefore a nonattained infimum/supremum. Do not apply this fact before a
  # transition at the evaluation endpoint: through-time history uses its new
  # state. At the horizon every coordinate is the observed accrued value.
  lower_open <- endpoint_open && remain > 0
  upper_open <- lower_open
  lower_open <- lower_open & reward != minimum
  upper_open <- upper_open & reward != maximum
  lower_distance <- (lower - center) / scale
  upper_distance <- (center - upper) / scale
  low_score <- max(0, lower_distance, upper_distance)
  high_score <- max(abs((lower - center) / scale),
                    abs((upper - center) / scale))
  low_open <- any(lower_open & lower_distance == low_score) ||
    any(upper_open & upper_distance == low_score)
  list(lower = lower, upper = upper, lower_open = lower_open,
       upper_open = upper_open, low_score = low_score,
       high_score = high_score, low_open = low_open)
}

#' Classify all thresholds using one structural snapshot
#'
#' @returns Logical resolved/membership vectors, coordinate envelopes, score
#'   bounds and a conservative-envelope flag. Unresolved membership is NA.
#' @noRd
classify_coverage_box <- function(accrued, state, time, problem, center, scale, q) {
  b <- cdt_box_bounds(accrued, state, time, problem, center, scale)
  inside <- q >= b$high_score
  outside <- q < b$low_score | (q == b$low_score & b$low_open)
  membership <- rep(NA, length(q))
  membership[outside] <- FALSE
  membership[inside] <- TRUE
  c(b, list(resolved = inside | outside, membership = membership, outer = TRUE))
}

#' Count sorted thresholds below, or at, a score boundary
#' @noRd
cdt_count <- function(q, value, inclusive) {
  left <- 0L
  right <- length(q) + 1L
  while (right - left > 1L) {
    middle <- as.integer(floor((left + right) / 2))
    if (q[middle] < value || (inclusive && q[middle] == value)) {
      left <- middle
    } else {
      right <- middle
    }
  }
  left
}

#' Find affine within-sojourn certificate roots
#'
#' @param segment Cached coordinate bounds, current reward and reachable reward
#'   extrema, interval endpoints and target scales for a single sojourn.
#' @param q Thresholds known to be certified by the segment endpoint.
#' @param membership TRUE for inclusion, FALSE for exclusion.
#' @returns Candidate times obtained analytically, one per threshold.
#' @noRd
candidate_determination_times <- function(segment, q, membership) {
  if (!length(q)) return(numeric())
  z <- segment
  elapsed <- if (membership) rep(0, length(q)) else rep(Inf, length(q))
  for (j in seq_along(z$reward)) {
    rise <- z$reward[j] - z$minimum[j]
    fall <- z$maximum[j] - z$reward[j]
    if (membership) {
      required_low <- z$center[j] - q * z$scale[j] - z$lower[j]
      required_high <- z$upper[j] - z$center[j] - q * z$scale[j]
      # Compare stationary endpoints on the score scale. At equality, q
      # times scale can round below the distance used to construct q.
      t_low <- if (rise > 0) required_low / rise else
        ifelse((z$center[j] - z$lower[j]) / z$scale[j] <= q, 0, Inf)
      t_high <- if (fall > 0) required_high / fall else
        ifelse((z$upper[j] - z$center[j]) / z$scale[j] <= q, 0, Inf)
      elapsed <- pmax(elapsed, t_low, t_high)
    } else {
      t_low <- if (rise > 0)
        (z$center[j] + q * z$scale[j] - z$lower[j]) / rise else
          rep(Inf, length(q))
      t_high <- if (fall > 0)
        (z$upper[j] - z$center[j] + q * z$scale[j]) / fall else
          rep(Inf, length(q))
      elapsed <- pmin(elapsed, t_low, t_high)
    }
  }
  z$start + pmax(0, elapsed)
}

#' Advance a computed root by one representable floating-point value
#' @useDynLib lifeconformal, .registration = TRUE
#' @noRd
cdt_nextafter <- function(x) {
  .Call("life_nextafter", as.double(x), PACKAGE = "lifeconformal")
}

#' Recheck analytic roots in vector batches
#' @noRd
cdt_verify_roots <- function(roots, segment, accrued, q, membership, tau,
                              event_at_end) {
  z <- segment
  valid_time <- is.finite(roots) & roots >= z$start & roots <= z$end
  verify_time <- roots
  verify_time[!valid_time] <- z$end
  remain <- pmax(0, tau - verify_time)
  can_open <- remain > 0 & !(event_at_end & verify_time == z$end)
  certified <- rep(membership, length(roots))
  for (j in seq_along(z$reward)) {
    # Complete the original binary totals, rather than advancing already
    # rounded bounds. The latter can erase an inclusive endpoint at a root.
    lower <- if (z$reward[j] == z$minimum[j]) rep(z$lower[j], length(roots)) else
      vapply(verify_time, function(u) markov_sum3(u, z$start, accrued[j]), numeric(1))
    upper <- if (z$reward[j] == z$maximum[j]) rep(z$upper[j], length(roots)) else
      vapply(verify_time, function(u) markov_sum3(tau, u, accrued[j]), numeric(1))
    dl <- (lower - z$center[j]) / z$scale[j]
    du <- (z$center[j] - upper) / z$scale[j]
    if (membership) {
      certified <- certified & pmax(abs(dl), abs(du)) <= q
    } else {
      certified <- certified | dl > q | du > q |
        (dl == q & can_open & z$reward[j] != z$minimum[j]) |
        (du == q & can_open & z$reward[j] != z$maximum[j])
    }
  }
  valid_time & certified
}

#' Sweep sorted score thresholds through observed sojourns once
#'
#' @param prepared An object returned by prepare_cdt_sweep().
#' @param q Strictly increasing nonnegative numeric thresholds; Inf is allowed.
#' @returns A life_cdt object. time is the observed stopping endpoint X_q,
#'   resolved is R_q, and membership is known only when resolved. certificate_time
#'   is NA for unresolved observations. Matrix rows follow subjects and columns
#'   follow q. The result reports conservative coordinate-envelope certification
#'   and the number of floating-point root checks requiring a later certificate.
#' @details Each history is scanned once. On each segment, binary searches find
#'   the newly certified threshold prefix and suffix, and their affine roots are
#'   solved in batches. Every threshold is removed from the unresolved block
#'   once. With a bounded number of floating-point refinements, work is
#'   O(m (d + log L) + L d) per subject, with O(L + d) working
#'   memory, apart from the requested subject-by-threshold output matrices.
#'
#'   A transition at a segment endpoint is processed before applying open-bound
#'   certificates at that time. The region includes its boundary. A root that
#'   cannot be verified in floating-point arithmetic is advanced through at most
#'   eight adjacent representable double values, with each candidate rechecked.
#'   Remaining uncertain roots are bracketed by the certified segment endpoint
#'   and refined by up to 64 bisections, retaining the certified upper endpoint.
#'   All numerical later certificates are recorded; no user-scale epsilon is added.
#' @noRd
sweep_cdt_thresholds <- function(prepared, q) {
  if (!inherits(prepared, "life_cdt_prepared")) {
    stop("prepared must come from prepare_cdt_sweep().", call. = FALSE)
  }
  if (!is.numeric(q) || is.complex(q) || !is.null(dim(q)) || !length(q) ||
      anyNA(q) || any(q < 0) || anyDuplicated(q) || any(diff(q) <= 0)) {
    stop("q must be a strictly increasing nonnegative numeric vector.", call. = FALSE)
  }
  answer <- .Call("life_cdt_sweep", prepared, as.double(q), PACKAGE = "lifeconformal")
  answer$q <- q
  answer$subject_id <- prepared$history$subjects$id
  answer$outer <- TRUE
  answer$certificate <- "reachable_coordinate_envelope"
  class(answer) <- "life_cdt"
  answer
}

#' R implementation retained to verify the compiled shared sweep
#' @noRd
sweep_cdt_thresholds_R <- function(prepared, q) {
  if (!inherits(prepared, "life_cdt_prepared")) {
    stop("prepared must come from prepare_cdt_sweep().", call. = FALSE)
  }
  if (!is.numeric(q) || is.complex(q) || !is.null(dim(q)) || !length(q) ||
      anyNA(q) || any(q < 0) || anyDuplicated(q) || any(diff(q) <= 0)) {
    stop("q must be a strictly increasing nonnegative numeric vector.", call. = FALSE)
  }
  n <- prepared$n
  L <- length(q)
  p <- prepared$problem
  h <- prepared$history
  interval_stop <- h$intervals$stop
  interval_event <- h$intervals$event
  interval_to <- h$intervals$to
  interval_censor <- h$intervals$terminal_censor
  endpoint <- pmin(h$subjects$observation_end, p$tau)
  time <- matrix(endpoint, n, L)
  membership <- matrix(NA, n, L)
  certificate_time <- matrix(NA_real_, n, L)
  later <- matrix(FALSE, n, L)
  segment_count <- 0L
  bound_count <- 0L
  root_count <- 0L
  fallbacks <- 0L
  root_adjustments <- 0L
  root_refinements <- 0L
  max_adjustment <- 0
  for (i in seq_len(n)) {
    mu <- prepared$center[i, ]
    sc <- prepared$scale[i, ]
    current <- h$subjects$initial_state[i]
    initial_reward <- p$rewards[current, ]
    constant_reward <- rep(TRUE, prepared$d)
    accrued <- numeric(prepared$d)
    t <- 0
    lo <- 1L
    hi <- L
    record <- function(indices, when, value, fallback = FALSE) {
      if (!length(indices)) return(invisible(NULL))
      time[i, indices] <<- when
      membership[i, indices] <<- value
      certificate_time[i, indices] <<- when
      later[i, indices] <<- fallback
      invisible(NULL)
    }
    classify_now <- function(b, when) {
      bound_count <<- bound_count + 1L
      last_out <- min(hi, cdt_count(q, b$low_score, b$low_open))
      if (last_out >= lo) {
        record(seq.int(lo, last_out), when, FALSE)
        lo <<- last_out + 1L
      }
      first_in <- max(lo, cdt_count(q, b$high_score, FALSE) + 1L)
      if (first_in <= hi) {
        record(seq.int(first_in, hi), when, TRUE)
        hi <<- first_in - 1L
      }
      invisible(NULL)
    }
    start_bounds <- cdt_box_bounds(accrued, current, t, p, mu, sc)
    classify_now(start_bounds, t)
    if (lo > hi) next
    for (r in h$row_index[[i]]) {
      if (lo > hi || t >= p$tau) break
      end <- min(interval_stop[r], p$tau)
      event_at_end <- interval_event[r] && interval_stop[r] <= p$tau
      duration <- end - t
      reward <- p$rewards[current, ]
      end_accrued <- accrued + duration * reward
      constant_reward <- constant_reward & (duration == 0 | reward == initial_reward)
      end_accrued[constant_reward] <- end * initial_reward[constant_reward]
      end_bounds <- cdt_box_bounds(end_accrued, current, end, p, mu, sc,
                                   endpoint_open = !event_at_end)
      bound_count <- bound_count + 1L
      segment_count <- segment_count + 1L
      segment <- list(start = t, end = end, lower = start_bounds$lower,
                      upper = start_bounds$upper, reward = reward,
                      minimum = p$structural_support$reward_min[current, ],
                      maximum = p$structural_support$reward_max[current, ],
                      center = mu, scale = sc)
      last_out <- min(hi, cdt_count(q, end_bounds$low_score, end_bounds$low_open))
      first_in <- max(lo, cdt_count(q, end_bounds$high_score, FALSE) + 1L)
      finish <- function(indices, value) {
        if (!length(indices)) return(invisible(NULL))
        roots <- candidate_determination_times(segment, q[indices], value)
        root_count <<- root_count + length(indices)
        fallback <- rep(FALSE, length(indices))
        # Verify newly removed thresholds in target-wise vector batches;
        # there is no R-level loop over thresholds or bound recomputation for
        # the unresolved grid at every history row.
        initial_roots <- roots
        ok <- cdt_verify_roots(roots, segment, accrued, q[indices], value,
                                p$tau, event_at_end)
        for (step in seq_len(8L)) {
          need <- which(!ok & is.finite(roots) & roots >= t & roots < end)
          if (!length(need)) break
          roots[need] <- cdt_nextafter(roots[need])
          bound_count <<- bound_count + length(need)
          ok[need] <- cdt_verify_roots(roots[need], segment, accrued,
                                      q[indices[need]], value, p$tau, event_at_end)
        }
        adjacent_adjusted <- ok & roots != initial_roots
        refine <- which(!ok & is.finite(roots) & roots >= t & roots < end)
        if (length(refine)) {
          upper <- rep(end, length(refine))
          bound_count <<- bound_count + length(refine)
          upper_ok <- cdt_verify_roots(upper, segment, accrued,
                                       q[indices[refine]], value, p$tau, event_at_end)
          refine <- refine[upper_ok]
          upper <- upper[upper_ok]
          if (length(refine)) {
            lower <- roots[refine]
            for (iteration in seq_len(64L)) {
              middle <- lower + (upper - lower) / 2
              moving <- middle > lower & middle < upper
              if (!any(moving)) break
              bound_count <<- bound_count + length(refine)
              valid <- cdt_verify_roots(middle, segment, accrued,
                                        q[indices[refine]], value, p$tau, event_at_end)
              upper[moving & valid] <- middle[moving & valid]
              lower[moving & !valid] <- middle[moving & !valid]
            }
            roots[refine] <- upper
            ok[refine] <- TRUE
            root_refinements <<- root_refinements + length(refine)
          }
        }
        fallback <- !ok
        roots[fallback] <- end
        adjusted <- roots != initial_roots
        finite_adjustment <- adjusted & is.finite(initial_roots)
        if (any(finite_adjustment)) {
          max_adjustment <<- max(max_adjustment,
                                  roots[finite_adjustment] - initial_roots[finite_adjustment])
        }
        root_adjustments <<- root_adjustments + sum(adjacent_adjusted)
        bound_count <<- bound_count + length(indices)
        fallbacks <<- fallbacks + sum(fallback)
        record(indices, roots, value, adjusted | fallback)
        invisible(NULL)
      }
      if (last_out >= lo) {
        finish(seq.int(lo, last_out), FALSE)
        lo <- last_out + 1L
      }
      first_in <- max(lo, first_in)
      if (first_in <= hi) {
        finish(seq.int(first_in, hi), TRUE)
        hi <- first_in - 1L
      }
      accrued <- end_accrued
      t <- end
      if (event_at_end) current <- interval_to[r]
      start_bounds <- cdt_box_bounds(accrued, current, t, p, mu, sc)
      if (lo <= hi) classify_now(start_bounds, t)
      if (t >= p$tau || current %in% p$absorbing || interval_censor[r]) break
    }
  }
  resolved <- !is.na(membership)
  structure(list(q = q, time = time, resolved = resolved,
                 membership = membership, certificate_time = certificate_time,
                 numerical_later_certificate = later, endpoint = endpoint,
                 subject_id = h$subjects$id, outer = TRUE,
                 certificate = "reachable_coordinate_envelope",
                 diagnostics = list(segments = segment_count,
                                    bound_evaluations = bound_count,
                                    analytic_roots = root_count,
                                    numerical_later_certificates = sum(later),
                                    adjacent_float_adjustments = root_adjustments,
                                    bracket_refinements = root_refinements,
                                    endpoint_fallbacks = fallbacks,
                                    max_root_delay = max_adjustment,
                                    complete_subjects = sum(rowSums(resolved) == L))),
            class = "life_cdt")
}

#' Find structural coverage certificates for fixed score regions
#' @noRd
compute_cdt <- function(history, problem, center, scale, q) {
  sweep_cdt_thresholds(prepare_cdt_sweep(history, problem, center, scale), q)
}

#' Independent one-threshold reference for coverage certification
#'
#' @details This reference scans each threshold's complete history separately,
#'   evaluates all affine candidate roots, and checks certificates directly.
#'   It is for small-case numerical verification, not production calibration.
#' @noRd
compute_cdt_reference <- function(history, problem, center, scale, q) {
  prepared <- prepare_cdt_sweep(history, problem, center, scale)
  # Validate q once using the production boundary, without a duplicate scan.
  if (!is.numeric(q) || !length(q) || anyNA(q) || any(q < 0) ||
      any(diff(q) <= 0)) stop("Invalid reference thresholds.", call. = FALSE)
  n <- prepared$n
  times <- matrix(pmin(history$subjects$observation_end, problem$tau), n, length(q))
  membership <- matrix(NA, n, length(q))
  for (i in seq_len(n)) for (k in seq_along(q)) {
    a <- numeric(prepared$d)
    state <- history$subjects$initial_state[i]
    initial_reward <- problem$rewards[state, ]
    constant_reward <- rep(TRUE, prepared$d)
    t <- 0
    mu <- prepared$center[i, ]
    sc <- prepared$scale[i, ]
    classify <- function(accrued, s, when, open = TRUE) {
      b <- cdt_box_bounds(accrued, s, when, problem, mu, sc, open)
      if (q[k] >= b$high_score) return(TRUE)
      if (q[k] < b$low_score || (q[k] == b$low_score && b$low_open)) return(FALSE)
      NA
    }
    value <- classify(a, state, t)
    if (!is.na(value)) {
      times[i, k] <- t
      membership[i, k] <- value
      next
    }
    for (r in history$row_index[[i]]) {
      row <- history$intervals[r, ]
      end <- min(row$stop, problem$tau)
      event <- row$event && row$stop <= problem$tau
      reward <- problem$rewards[state, ]
      bounds <- cdt_box_bounds(a, state, t, problem, mu, sc)
      rise <- reward - problem$structural_support$reward_min[state, ]
      fall <- problem$structural_support$reward_max[state, ] - reward
      root_low_in <- ifelse(rise > 0, (mu - q[k] * sc - bounds$lower) / rise,
                            ifelse((mu - bounds$lower) / sc <= q[k], 0, Inf))
      root_up_in <- ifelse(fall > 0, (bounds$upper - mu - q[k] * sc) / fall,
                           ifelse((bounds$upper - mu) / sc <= q[k], 0, Inf))
      root_low_out <- ifelse(rise > 0, (mu + q[k] * sc - bounds$lower) / rise, Inf)
      root_up_out <- ifelse(fall > 0, (bounds$upper - mu + q[k] * sc) / fall, Inf)
      candidates <- sort(unique(c(t + max(0, root_low_in, root_up_in),
                                   t + pmax(0, root_low_out),
                                   t + pmax(0, root_up_out), end)))
      candidates <- candidates[is.finite(candidates) & candidates >= t & candidates <= end]
      for (u in candidates) {
        candidate_accrued <- markov_protected_total(u, t, reward, a)
        value <- classify(candidate_accrued, state, u, !(event && u == end))
        if (!is.na(value)) {
          times[i, k] <- u
          membership[i, k] <- value
          break
        }
      }
      if (!is.na(membership[i, k])) break
      a <- a + (end - t) * reward
      constant_reward <- constant_reward & (end == t | reward == initial_reward)
      a[constant_reward] <- end * initial_reward[constant_reward]
      t <- end
      if (event) state <- row$to
      value <- classify(a, state, t)
      if (!is.na(value)) {
        times[i, k] <- t
        membership[i, k] <- value
        break
      }
      if (t >= problem$tau || state %in% problem$absorbing || row$terminal_censor) break
    }
  }
  list(time = times, resolved = !is.na(membership), membership = membership)
}
