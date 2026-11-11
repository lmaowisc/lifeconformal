#' Detect binary targets that are irreversible along fitted transitions
#'
#' This is a working-model computation check, not a structural assertion about
#' possible patient histories. Cycles are allowed within equal-reward classes.
#' @param rewards Numeric state-by-target binary reward matrix.
#' @param edges A two-column matrix/data frame with one-based from/to indices.
#' @returns TRUE if no supplied edge increases any target reward.
#' @noRd
markov_monotone_eligible <- function(rewards, edges) {
  if (!is.matrix(rewards) || !is.numeric(rewards) || is.complex(rewards) ||
      any(dim(rewards) < 1L) || anyNA(rewards) ||
      any(!rewards %in% c(0, 1))) return(FALSE)
  if ((!is.matrix(edges) && !is.data.frame(edges)) || ncol(edges) != 2L) {
    stop("edges must have two columns containing state indices.", call. = FALSE)
  }
  # Base R may give a zero-row edge matrix logical storage. With no edges,
  # monotonicity is vacuous and there are no state indices to validate.
  if (!nrow(edges)) return(TRUE)
  from <- edges[, 1L]
  to <- edges[, 2L]
  if (!is.numeric(from) || !is.numeric(to) || is.complex(from) || is.complex(to) ||
      anyNA(from) || anyNA(to) ||
      any(from != floor(from) | from < 1 | from > nrow(rewards)) ||
      any(to != floor(to) | to < 1 | to > nrow(rewards))) {
    stop("edges must contain valid one-based state indices.", call. = FALSE)
  }
  !any(rewards[to, , drop = FALSE] > rewards[from, , drop = FALSE])
}

#' Exact conditional coverage for monotone binary Markov rewards
#'
#' A binary reward that never increases along a positive-probability fitted
#' edge is a first-exit clock: while its current reward is one, its final value
#' is accrued reward plus the time until exit (or the remaining horizon).
#' Consequently a backward state recursion suffices, without an occupation
#' grid, atom enumeration, or simulation. A transition that switches off a
#' target is retained only for thresholds covering its exit clock. The terminal
#' value checks the horizon atom for targets still active in each state.
#'
#' Queries share one pass whenever their active-target deadline offsets and
#' scales agree. For ordinary origin-based histories, accrued reward equals the
#' query time for every still-active target, so all queries of a baseline
#' profile share a single pass. Targets already inactive are tested against
#' their known accrued value separately. Arbitrary finite accrued values are
#' also supported by grouping distinct effective deadline specifications.
#'
#' @param jump_times,kernels,rewards,tau See markov_reward_moments(). Rewards
#'   must be binary and nonincreasing along every positive-probability edge.
#' @param query_times Conditioning times in [0, tau].
#' @param query_state One-based current state per query.
#' @param accrued,center,scale Query-by-target matrices; a target vector is
#'   accepted for one query. Scales must be positive.
#' @param q Finite nonnegative strictly increasing score thresholds.
#' @param side One common left/right convention or one per query. Right means
#'   history through a jump and therefore excludes that jump from the future.
#' @param flip Logical target vector identifying complemented rewards.
#' @param original_center,original_accrued Original target-unit query matrices
#'   for complemented targets. Direct original-unit score evaluation preserves
#'   inclusive atom comparisons; complemented calls require active accrued
#'   clocks to equal their conditioning time.
#' @returns A list with identical query-by-threshold `probability`, `lower`,
#'   and `upper` matrices and exact-method diagnostics. Numerical rounding is
#'   the only approximation. No unobserved working transition is used as a
#'   structural coverage-determination certificate.
#' @details
#' With G distinct deadline groups, J jumps, E fitted edges, K states, H queries,
#'   d targets and L thresholds, work is O(G J (K + E) L + H d L), with rolling
#'   O(K L) recursion storage plus the requested H-by-L output. Each edge uses
#'   one score cutoff shared by all thresholds. Groups and event/query indices
#'   are prepared once. Native batching can use this routine as a reference.
#' @noRd
markov_conditional_coverage_monotone <- function(
    jump_times, kernels, rewards, tau, query_times, query_state,
    accrued, center, scale, q, side = "right", flip = NULL,
    original_center = NULL, original_accrued = NULL) {
  spec <- prepare_markov_recursion(jump_times, kernels, rewards, tau,
                                    0, "left", NULL)
  prepared <- prepare_monotone_queries(query_times, query_state, accrued,
                                       center, scale, q, side, rewards, tau)
  h <- length(query_times)
  d <- ncol(rewards)
  k <- nrow(rewards)
  kernels <- spec$kernels
  jump_times <- spec$jump_times
  used <- matrix(FALSE, k, k)
  for (kernel in kernels) used <- used | kernel > 0
  edges <- which(used, arr.ind = TRUE)
  if (!markov_monotone_eligible(rewards, edges)) {
    stop("Exact clock recursion requires binary rewards that never increase along fitted edges.",
         call. = FALSE)
  }
  accrued <- prepared$accrued
  center <- prepared$center
  scale <- prepared$scale
  side <- prepared$side
  current_rewards <- rewards[query_state, , drop = FALSE]
  if (is.null(flip)) flip <- rep(FALSE, d)
  if (!is.logical(flip) || !is.null(dim(flip)) || length(flip) != d || anyNA(flip))
    stop("flip must contain one logical value per target.", call. = FALSE)
  canonical <- any(flip)
  if (canonical) {
    if (is.null(original_center) || is.null(original_accrued))
      stop("Complemented clocks need original_center and original_accrued.", call. = FALSE)
    original <- prepare_monotone_queries(query_times, query_state, original_accrued,
                                         original_center, scale, q, side, rewards, tau)
    original_center <- original$center
    original_accrued <- original$accrued
    active <- current_rewards == 1
    if (any(accrued[active] != matrix(query_times, h, d)[active]))
      stop("Complemented clock queries require active accrued rewards equal to time.",
           call. = FALSE)
    if (any(original_accrued[, flip, drop = FALSE][active[, flip, drop = FALSE]] != 0))
      stop("Active complemented clocks require zero original accrued reward.", call. = FALSE)
  }
  effective_center <- if (canonical) original_center else center
  # An arbitrary observed accrued value is added to each future exit clock
  # before centering. Absorbing it into an effective center reassociates the
  # arithmetic and can change an inclusive atom comparison by one ulp.
  offset_targets <- current_rewards == 1 & accrued != matrix(query_times, h, d)
  offset_time <- offset_accrued <- matrix(0, h, d)
  offset_time[offset_targets] <- matrix(query_times, h, d)[offset_targets]
  offset_accrued[offset_targets] <- accrued[offset_targets]
  fixed_score <- matrix(abs((accrued - center) / scale), h, d)
  if (canonical) {
    fixed_value <- original_accrued
    for (a in seq_len(h)) fixed_value[a, flip] <-
      markov_protected_total(tau, query_times[a], rep(1, sum(flip)),
                             original_accrued[a, flip])
    fixed_score <- abs((fixed_value - original_center) / scale)
  }
  fixed_score[current_rewards == 1] <- 0
  fixed_score <- if (d == 1L) as.vector(fixed_score) else
    do.call(pmax, unname(as.data.frame(fixed_score)))
  probability <- outer(fixed_score, q, "<=") * 1
  active_rows <- which(rowSums(current_rewards) > 0 & fixed_score <= max(q))
  if (!length(active_rows)) {
    return(monotone_coverage_result(probability, h, length(jump_times), 0L,
                                    spec$kernel_row_error))
  }
  # Hexadecimal doubles form exact keys: decimal printing must not collapse
  # nearby but distinct centers or scales into a falsely shared pass.
  parameters <- cbind(effective_center, scale, offset_time, offset_accrued,
                       offset_targets * 1)[active_rows, , drop = FALSE]
  keys <- do.call(paste, c(lapply(seq_len(ncol(parameters)), function(j)
    sprintf("%a", parameters[, j])), sep = ":"))
  groups <- split(active_rows, match(keys, unique(keys)))
  bins <- findInterval(query_times, jump_times)
  at_jump <- match(query_times, jump_times, nomatch = 0L)
  before <- side == "left" & at_jump > 0L
  bins[before] <- bins[before] - 1L
  # Each edge's switched-off targets are immutable across times and groups.
  off_diagonal <- which(used & row(used) != col(used), arr.ind = TRUE)
  drops <- lapply(seq_len(nrow(off_diagonal)), function(e)
    which(rewards[off_diagonal[e, 1L], ] > rewards[off_diagonal[e, 2L], ]))
  for (rows in groups) {
    answer <- monotone_clock_pass(jump_times, kernels, rewards, tau,
                                  effective_center[rows[1L], ],
                                  scale[rows[1L], ], q, bins[rows],
                                  query_state[rows], off_diagonal, drops, flip,
                                  offset_targets[rows[1L], ],
                                  offset_time[rows[1L], ], offset_accrued[rows[1L], ])
    probability[rows, ] <- probability[rows, , drop = FALSE] * answer
  }
  monotone_coverage_result(probability, h, length(jump_times), length(groups),
                            spec$kernel_row_error)
}

#' Validate query arrays without changing their dimensions
#' @noRd
prepare_monotone_queries <- function(query_times, query_state, accrued,
                                      center, scale, q, side, rewards, tau) {
  numeric_vector <- function(x) is.numeric(x) && !is.complex(x) &&
    is.null(dim(x)) && all(is.finite(x))
  h <- length(query_times)
  d <- ncol(rewards)
  if (!numeric_vector(query_times) || !h ||
      any(query_times < 0 | query_times > tau)) {
    stop("query_times must be nonempty and lie in [0, tau].", call. = FALSE)
  }
  if (!numeric_vector(query_state) || length(query_state) != h ||
      any(query_state != floor(query_state) | query_state < 1 |
          query_state > nrow(rewards))) {
    stop("Supply one valid query_state per query.", call. = FALSE)
  }
  if (!numeric_vector(q) || !length(q) || any(q < 0) || any(diff(q) <= 0)) {
    stop("q must be finite, nonnegative and strictly increasing.", call. = FALSE)
  }
  if (length(side) == 1L) side <- rep(side, h)
  if (!is.character(side) || length(side) != h || anyNA(side) ||
      any(!side %in% c("left", "right"))) {
    stop("side must contain one left/right value per query or one common value.",
         call. = FALSE)
  }
  as_query_matrix <- function(x, name, positive = FALSE) {
    if (h == 1L && numeric_vector(x) && length(x) == d) x <- matrix(x, 1L)
    if (!is.matrix(x) || !is.numeric(x) || is.complex(x) ||
        !identical(dim(x), c(h, d)) || any(!is.finite(x)) ||
        (positive && any(x <= 0))) {
      stop(name, " must be a finite query-by-target matrix",
           if (positive) " with positive entries." else ".", call. = FALSE)
    }
    x
  }
  list(accrued = as_query_matrix(accrued, "accrued"),
       center = as_query_matrix(center, "center"),
       scale = as_query_matrix(scale, "scale", TRUE), side = side)
}

#' Backward substochastic recursion with shared ordered threshold cutoffs
#' @noRd
monotone_clock_pass <- function(jump_times, kernels, rewards, tau, center,
                                 scale, q, query_bins, query_state, edges, drops,
                                 flip = rep(FALSE, ncol(rewards)),
                                 offset_targets = rep(FALSE, ncol(rewards)),
                                 offset_time = numeric(ncol(rewards)),
                                 offset_accrued = numeric(ncol(rewards))) {
  k <- nrow(rewards)
  terminal_value <- ifelse(flip, 0, tau)
  for (j in which(offset_targets)) terminal_value[j] <-
    markov_protected_total(tau, offset_time[j], 1, offset_accrued[j])
  threshold <- abs((terminal_value - center) / scale)
  terminal_score <- vapply(seq_len(k), function(s) {
    active <- rewards[s, ] == 1
    if (any(active)) max(threshold[active]) else 0
  }, numeric(1))
  value <- outer(terminal_score, q, "<=") * 1
  result <- matrix(0, length(query_bins), length(q))
  order <- order(query_bins, decreasing = TRUE)
  cursor <- 1L
  emit <- function(bin) {
    start <- cursor
    while (cursor <= length(order) && query_bins[order[cursor]] == bin) {
      cursor <<- cursor + 1L
    }
    if (cursor > start) {
      idx <- order[seq.int(start, cursor - 1L)]
      result[idx, ] <<- value[query_state[idx], , drop = FALSE]
    }
  }
  emit(length(jump_times))
  needed_jumps <- if (min(query_bins) < length(jump_times))
    seq.int(length(jump_times), min(query_bins) + 1L) else integer()
  for (j in needed_jumps) {
    kernel <- kernels[[j]]
    next_value <- value * diag(kernel)
    for (e in seq_len(nrow(edges))) {
      source <- edges[e, 1L]
      destination <- edges[e, 2L]
      p <- kernel[source, destination]
      if (p == 0) next
      changed <- drops[[e]]
      if (!length(changed)) {
        next_value[source, ] <- next_value[source, ] + p * value[destination, ]
      } else {
        original_value <- ifelse(flip[changed], tau - jump_times[j], jump_times[j])
        shifted <- offset_targets[changed]
        for (a in which(shifted)) original_value[a] <-
          markov_protected_total(jump_times[j], offset_time[changed[a]], 1,
                                 offset_accrued[changed[a]])
        edge_score <- max(abs((original_value - center[changed]) / scale[changed]))
        cutoff <- findInterval(edge_score, q, left.open = TRUE)
        if (cutoff < length(q)) {
          keep <- seq.int(cutoff + 1L, length(q))
          next_value[source, keep] <- next_value[source, keep] +
            p * value[destination, keep]
        }
      }
    }
    value <- next_value
    emit(j - 1L)
  }
  result
}

#' Assemble exact-clock results with explicit numerical metadata
#' @noRd
monotone_coverage_result <- function(probability, h, j, groups, row_error) {
  if (any(!is.finite(probability)) || any(probability < -1e-12 | probability > 1 + 1e-12)) {
    stop("Exact Markov clock recursion produced invalid probabilities.", call. = FALSE)
  }
  probability[] <- pmin(1, pmax(0, probability))
  list(probability = probability, lower = probability, upper = probability,
       diagnostics = list(method = "deterministic_exact_monotone_clocks",
                          query_count = h, jump_count = j,
                          deadline_groups = groups, maximum_probability_error = 0,
                          tolerance_met = TRUE, kernel_row_error = row_error,
                          approximation = "floating_point_only"))
}
