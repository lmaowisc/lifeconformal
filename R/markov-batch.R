#' Compute baseline Markov reward moments for a batch of Cox profiles
#'
#' Internal deterministic moment engine for a common jump grid with
#' transition-specific baseline increments and time-constant relative risks.
#' It applies the same stochastic kernels and stable total-variance recursion
#' as [markov_reward_moments()], without constructing a kernel list or running
#' an R-level time loop separately for every profile.
#'
#' @param jump_times Strictly increasing positive numeric vector, at most tau.
#' @param base_increments Numeric jump-by-edge matrix of nonnegative baseline
#'   intensity increments. A zero-row matrix is allowed for no jumps.
#' @param edges Data frame or matrix with numeric integer `from` and `to`
#'   columns giving distinct directed edges in the order of the increment
#'   and relative-risk columns. State indices start at one; diagonal edges
#'   are not allowed.
#' @param relative_risk Numeric profile-by-edge matrix of strictly positive
#'   finite relative risks. At least one profile is required. Named columns
#'   must agree in order with named columns of `base_increments`.
#' @param rewards Finite numeric state-by-target reward matrix.
#' @param tau Finite positive common prediction horizon.
#' @param initial_state One integer state index, or one per profile.
#' @param accrued Optional common target vector or profile-by-target matrix
#'   of already accrued rewards. It is added after propagation, leaving the
#'   future-reward variances unchanged. NULL means zero.
#' @param absorbing Integer indices of declared absorbing states. Their
#'   outgoing baseline increments must all be zero.
#' @param keep_state_moments Logical; also return profile-by-state-by-target
#'   arrays for the mean and variance at the origin.
#'
#' @returns A list with profile-by-target `mean`, `variance`, `second_moment`,
#'   and `diagnostics`. When requested, `state_mean` and `state_variance`
#'   contain all origin states. State means include the supplied accrued
#'   rewards; the user of that internal option must choose compatible histories.
#' @details
#' This routine computes moments at the prediction origin; it does not supply
#' dynamic history queries or a conditional target distribution. At each jump,
#' outgoing increments are multiplied by the profile's relative risks. Rows
#' with a sum above one are rescaled to sum to one.
#' The diagonal is the remaining probability. Every origin uses the same
#' pre-update destination moments when several origins change simultaneously.
#' Positive transitions that underflow to zero during multiplication or row
#' normalization are rejected explicitly; supplied zero increments and
#' representable subnormal probabilities remain valid.
#'
#' States are shifted in time only when touched by an active edge, or when
#' the output is assembled. Per-state time anchors avoid repeatedly shifting
#' every state, and avoid subtracting nearly equal horizon-sized offsets when
#' the true occupation time is tiny. Mixtures use nonnegative weighted
#' centered squares, rather than subtracting two raw second moments.
#'
#' R-level loops run over active jumps and edges, not profiles. Arithmetic is
#' O(N d sum_j (E_j + K_j)), where E_j is the active edge count and K_j is the
#' number of touched states. Validation additionally reads the supplied J-by-E
#' and N-by-E matrices once. Working storage is O(N K d + N E), excluding
#' inputs and optional output arrays; no N-by-J-by-K kernel array is formed.
#' @noRd
markov_reward_moments_batch <- function(jump_times, base_increments, edges,
                                        relative_risk, rewards, tau,
                                        initial_state, accrued = NULL,
                                        absorbing = integer(0),
                                        keep_state_moments = FALSE) {
  spec <- validate_markov_batch(jump_times, base_increments, edges,
                                 relative_risk, rewards, tau, initial_state,
                                 accrued, absorbing, keep_state_moments)
  n <- nrow(relative_risk)
  k <- nrow(rewards)
  d <- ncol(rewards)
  nd <- n * d
  means <- rep(list(numeric(nd)), k)
  variances <- means
  no_underflow <- rep(FALSE, nd)
  variance_underflow <- rep(list(no_underflow), k)
  anchors <- rep(tau, k)
  # Target blocks contain consecutive profiles; profile weights therefore
  # recycle across target blocks without expanding to an N-by-d matrix.
  expanded_rewards <- lapply(seq_len(k), function(s) rep(rewards[s, ], each = n))
  has_reward <- rowSums(abs(rewards)) > 0
  risks <- lapply(seq_len(ncol(relative_risk)), function(e) relative_risk[, e])
  adjusted_by_profile <- integer(n)
  adjusted_by_state <- integer(k)
  largest_outgoing <- 0
  active_rows <- which(rowSums(base_increments > 0) > 0)
  active_edges <- which(colSums(base_increments > 0) > 0)
  can_vary <- markov_reward_variability(rewards, spec$from[active_edges],
                                        spec$to[active_edges])
  expanded_variability <- lapply(seq_len(k), function(s) rep(can_vary[s, ], each = n))

  for (j in rev(active_rows)) {
    u <- jump_times[j]
    active <- which(base_increments[j, ] > 0)
    origins <- unique(spec$from[active])
    touched <- unique(c(origins, spec$to[active]))
    for (s in touched) {
      if (has_reward[s] && anchors[s] != u) {
        means[[s]] <- means[[s]] + (anchors[s] - u) * expanded_rewards[[s]]
      }
      anchors[s] <- u
    }
    updated_means <- means
    updated_variances <- variances
    updated_underflow <- variance_underflow
    for (s in origins) {
      out <- active[spec$from[active] == s]
      weights <- lapply(out, function(e) {
        value <- risks[[e]] * base_increments[j, e]
        if (any(value == 0))
          stop("Positive Markov intensity increments underflow; check the fitted hazards and prediction profiles.",
               call. = FALSE)
        value
      })
      total <- weights[[1L]]
      if (length(weights) > 1L) {
        for (h in seq.int(2L, length(weights))) total <- total + weights[[h]]
      }
      if (any(!is.finite(total))) {
        stop("Outgoing Markov increments overflow in the profile batch.", call. = FALSE)
      }
      largest_outgoing <- max(largest_outgoing, total)
      adjusted <- total > 1
      if (any(adjusted)) {
        scale <- pmax(1, total)
        weights <- lapply(weights, function(w) {
          value <- w / scale
          if (any(w > 0 & value == 0))
            stop("Positive Markov transition probabilities underflow during row normalization; check the fitted intensity increments.",
                 call. = FALSE)
          value
        })
        adjusted_by_profile <- adjusted_by_profile + adjusted
        adjusted_by_state[s] <- adjusted_by_state[s] + sum(adjusted)
        total <- weights[[1L]]
        if (length(weights) > 1L) {
          for (h in seq.int(2L, length(weights))) total <- total + weights[[h]]
        }
      }
      stay <- pmax(0, 1 - total)
      mass <- stay + total
      if (any(mass != 1)) {
        stay <- stay / mass
        weights <- lapply(weights, function(w) w / mass)
      }
      center <- means[[s]] * stay
      for (h in seq_along(out)) {
        center <- center + means[[spec$to[out[h]]]] * weights[[h]]
      }
      # Scale before squaring: a zero or tiny mixture probability must not
      # turn an otherwise finite weighted variance into 0 * Inf.
      mixed_variance <- variances[[s]] * stay +
        ((means[[s]] - center) * sqrt(stay))^2
      for (h in seq_along(out)) {
        destination <- spec$to[out[h]]
        mixed_variance <- mixed_variance + variances[[destination]] * weights[[h]] +
          ((means[[destination]] - center) * sqrt(weights[[h]]))^2
      }
      if (any(expanded_variability[[s]] & mixed_variance < .Machine$double.xmin)) {
        # A zero variance may be a true point mass or a lost tiny spread.
        # Check actual positive-probability branches before rejecting units.
        lower <- rep(Inf, nd)
        upper <- rep(-Inf, nd)
        positive_variance <- rep(FALSE, nd)
        for (h in seq_len(length(out) + 1L)) {
          destination <- if (h == 1L) s else spec$to[out[h - 1L]]
          weight <- if (h == 1L) stay else weights[[h - 1L]]
          branch_active <- rep(weight > 0, d)
          lower[branch_active] <- pmin(lower[branch_active], means[[destination]][branch_active])
          upper[branch_active] <- pmax(upper[branch_active], means[[destination]][branch_active])
          positive_variance[branch_active] <- positive_variance[branch_active] |
            variances[[destination]][branch_active] > 0 |
            variance_underflow[[destination]][branch_active]
        }
        updated_underflow[[s]] <- mixed_variance < .Machine$double.xmin &
          (upper > lower | positive_variance)
      } else {
        updated_underflow[[s]] <- no_underflow
      }
      updated_means[[s]] <- center
      updated_variances[[s]] <- mixed_variance
    }
    means <- updated_means
    variances <- updated_variances
    variance_underflow <- updated_underflow
  }

  needed <- if (keep_state_moments) seq_len(k) else unique(spec$initial_state)
  for (s in needed) {
    if (has_reward[s] && anchors[s] != 0) {
      means[[s]] <- means[[s]] + anchors[s] * expanded_rewards[[s]]
    }
  }
  out_mean <- matrix(0, n, d, dimnames = spec$output_names)
  out_variance <- out_mean
  for (s in unique(spec$initial_state)) {
    selected <- which(spec$initial_state == s)
    positions <- as.vector(outer(selected, (seq_len(d) - 1L) * n, "+"))
    out_mean[positions] <- means[[s]][positions]
    out_variance[positions] <- variances[[s]][positions]
    check_markov_variance_scale(out_variance[positions],
                                variance_underflow[[s]][positions])
  }
  out_mean <- out_mean + spec$accrued
  second <- out_variance + out_mean^2
  if (any(!is.finite(out_mean)) || any(!is.finite(second)) ||
      any(!is.finite(out_variance))) {
    stop("Batched Markov moments overflow; check reward and time units.", call. = FALSE)
  }
  result <- list(mean = out_mean, variance = out_variance, second_moment = second,
                 diagnostics = list(method = "deterministic_batched_cox_moments",
                                    n_profiles = n, active_jump_count = length(active_rows),
                                    adjusted_by_profile = adjusted_by_profile,
                                    adjusted_by_state = adjusted_by_state,
                                    largest_outgoing_sum = largest_outgoing))
  if (keep_state_moments) {
    state_mean <- array(0, c(n, k, d),
                         dimnames = list(rownames(relative_risk), rownames(rewards),
                                         colnames(rewards)))
    state_variance <- state_mean
    for (s in seq_len(k)) {
      check_markov_variance_scale(variances[[s]], variance_underflow[[s]])
      state_mean[, s, ] <- matrix(means[[s]], n, d) + spec$accrued
      state_variance[, s, ] <- matrix(variances[[s]], n, d)
    }
    if (any(!is.finite(state_mean)) || any(!is.finite(state_variance))) {
      stop("Batched state moments overflow; check reward and time units.", call. = FALSE)
    }
    result$state_mean <- state_mean
    result$state_variance <- state_variance
  }
  result
}

#' Validate the shared-grid batched Markov moment inputs once
#' @noRd
validate_markov_batch <- function(jump_times, base_increments, edges,
                                   relative_risk, rewards, tau, initial_state,
                                   accrued, absorbing, keep_state_moments) {
  real_matrix <- function(x) is.matrix(x) && is.numeric(x) &&
    !is.complex(x) && all(is.finite(x))
  if (!real_matrix(rewards) || any(dim(rewards) < 1L)) {
    stop("`rewards` must be a finite numeric state-by-target matrix.", call. = FALSE)
  }
  k <- nrow(rewards)
  d <- ncol(rewards)
  if (!is.numeric(tau) || is.complex(tau) || !is.null(dim(tau)) ||
      length(tau) != 1L || !is.finite(tau) || tau <= 0) {
    stop("`tau` must be a finite positive scalar.", call. = FALSE)
  }
  if (!is.numeric(jump_times) || is.complex(jump_times) ||
      !is.null(dim(jump_times)) || any(!is.finite(jump_times)) ||
      any(jump_times <= 0 | jump_times > tau) || any(diff(jump_times) <= 0)) {
    stop("`jump_times` must be strictly increasing in (0, tau].", call. = FALSE)
  }
  if (!(is.data.frame(edges) || is.matrix(edges)) ||
      is.null(colnames(edges)) || !all(c("from", "to") %in% colnames(edges)) ||
      anyDuplicated(colnames(edges))) {
    stop("`edges` must have unique numeric from/to columns.", call. = FALSE)
  }
  from <- edges[, "from"]
  to <- edges[, "to"]
  valid_indices <- function(x, size = NULL) {
    is.numeric(x) && !is.complex(x) && is.null(dim(x)) &&
      (is.null(size) || length(x) == size) && all(is.finite(x)) &&
      all(x == floor(x) & x >= 1 & x <= k)
  }
  if (!valid_indices(from) || !valid_indices(to, length(from)) ||
      any(from == to) || anyDuplicated(data.frame(from, to))) {
    stop("`edges` must contain distinct off-diagonal state-index pairs.", call. = FALSE)
  }
  e <- length(from)
  if (!real_matrix(base_increments) ||
      !identical(dim(base_increments), c(length(jump_times), e)) ||
      any(base_increments < 0)) {
    stop("`base_increments` must be a nonnegative jump-by-edge matrix.", call. = FALSE)
  }
  if (!real_matrix(relative_risk) || nrow(relative_risk) < 1L ||
      ncol(relative_risk) != e || any(relative_risk <= 0)) {
    stop("`relative_risk` must be a positive finite profile-by-edge matrix.", call. = FALSE)
  }
  if (!is.null(colnames(base_increments)) && !is.null(colnames(relative_risk)) &&
      !identical(colnames(base_increments), colnames(relative_risk))) {
    stop("Named increment and risk columns must have the same edge order.", call. = FALSE)
  }
  n <- nrow(relative_risk)
  if (!valid_indices(initial_state) || !length(initial_state) ||
      !(length(initial_state) %in% c(1L, n))) {
    stop("`initial_state` must be one valid state index or one per profile.", call. = FALSE)
  }
  if (length(initial_state) == 1L) initial_state <- rep(initial_state, n)
  if (!valid_indices(absorbing) || anyDuplicated(absorbing)) {
    stop("`absorbing` must contain distinct valid state indices.", call. = FALSE)
  }
  if (length(absorbing) && any(base_increments[, from %in% absorbing, drop = FALSE] > 0)) {
    stop("Declared absorbing states have outgoing baseline increments.", call. = FALSE)
  }
  for (labels in dimnames(rewards)) {
    if (!is.null(labels) &&
        (anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels))) {
      stop("Reward dimension names must be distinct and nonempty.", call. = FALSE)
    }
  }
  output_names <- list(rownames(relative_risk), colnames(rewards))
  if (is.null(accrued)) {
    accrued <- matrix(0, n, d)
  } else if (is.numeric(accrued) && !is.complex(accrued) && is.null(dim(accrued)) &&
             length(accrued) == d && all(is.finite(accrued))) {
    if (!is.null(names(accrued)) && !is.null(colnames(rewards)) &&
        !identical(names(accrued), colnames(rewards))) {
      stop("Named accrued values must match the target order.", call. = FALSE)
    }
    accrued <- matrix(rep(accrued, each = n), n, d)
  } else if (!real_matrix(accrued) || !identical(dim(accrued), c(n, d))) {
    stop("`accrued` must be a finite target vector or profile-by-target matrix.", call. = FALSE)
  }
  if ((!is.null(colnames(accrued)) && !is.null(colnames(rewards)) &&
       !identical(colnames(accrued), colnames(rewards))) ||
      (!is.null(rownames(accrued)) && !is.null(rownames(relative_risk)) &&
       !identical(rownames(accrued), rownames(relative_risk)))) {
    stop("Named accrued matrix dimensions must match the profile and target order.", call. = FALSE)
  }
  dimnames(accrued) <- output_names
  if (!is.logical(keep_state_moments) || length(keep_state_moments) != 1L ||
      !is.null(dim(keep_state_moments)) || is.na(keep_state_moments)) {
    stop("`keep_state_moments` must be TRUE or FALSE.", call. = FALSE)
  }
  list(from = as.integer(from), to = as.integer(to),
       initial_state = as.integer(initial_state), accrued = accrued,
       output_names = output_names)
}
