#' Stream exact monotone Markov coverage into signed subject contributions
#'
#' Native backward recursion for binary rewards that never increase along
#' fitted edges. Each baseline profile shares one state-by-threshold recursion
#' across its through-time history queries. No target grid, path simulation,
#' or query-by-threshold matrix is constructed.
#'
#' @param jump_times,base_increments,edges,relative_risk,rewards,tau Shared
#'   Cox/Markov representation as in markov_reward_moments_batch().
#' @param center,scale Finite subject-by-target matrices; scale is positive.
#' @param q Increasing finite nonnegative score thresholds.
#' @param subject_index,query_times,query_state Subject index, time, and current
#'   state for each query. All queries condition through the event time.
#' @param accrued Query-by-target matrix. A currently active binary reward
#'   must have accrued exactly the query time. An incompatible observed history
#'   requires the more general conditional-distribution method.
#' @param weight Finite signed query weights. Common-endpoint augmentation uses
#'   the endpoint inverse survival and negative earlier inverse-survival jumps.
#' @param determination_time,resolved,membership Optional subject-by-threshold
#'   structural certificates. At or after a resolved certificate time, its known
#'   membership overrides the working-model probability exactly.
#' @param flip Optional logical target vector indicating complemented rewards.
#' @param original_center,original_accrued Original target-unit matrices when
#'   any reward is complemented. Scores are evaluated directly in these units
#'   to preserve inclusive atom comparisons under floating-point subtraction.
#' @returns Subject-by-threshold contributions and diagnostics. Working storage
#'   is O(K L + H + E), excluding input arrays and the n-by-L returned matrix.
#'   The result is exact for the fitted jump kernels apart from floating-point
#'   arithmetic. Signed accumulation uses compensated summation. Zero-weight
#'   queries do not extend the backward recursion. Model row rescaling is
#'   counted in diagnostics. Repeated query accumulation and structural
#'   certificates use contiguous per-subject storage before the final result
#'   row is copied into R's column-major matrix.
#'   Positive transition probabilities that underflow during multiplication
#'   or row normalization are rejected instead of treated as structural zeros.
#' @noRd
markov_monotone_batch <- function(jump_times, base_increments, edges,
                                  relative_risk, rewards, tau, center, scale, q,
                                  subject_index, query_times, query_state,
                                  accrued, weight, determination_time = NULL,
                                  resolved = NULL, membership = NULL,
                                  flip = NULL, original_center = NULL,
                                  original_accrued = NULL) {
  n <- nrow(relative_risk)
  if (is.null(n) || n < 1L) stop("relative_risk must have subject rows.", call. = FALSE)
  validate_markov_batch(jump_times, base_increments, edges, relative_risk,
                        rewards, tau, rep(1L, n), NULL, integer(0), FALSE)
  edges <- as.matrix(edges[, c("from", "to"), drop = FALSE])
  if (!markov_monotone_eligible(rewards, edges))
    stop("Exact batched clocks require nonincreasing binary rewards.", call. = FALSE)
  d <- ncol(rewards)
  h <- length(query_times)
  matrix_ok <- function(x, nr, nc) is.matrix(x) && is.numeric(x) &&
    !is.complex(x) && identical(dim(x), c(nr, nc)) && all(is.finite(x))
  if (!matrix_ok(center, n, d) || !matrix_ok(scale, n, d) || any(scale <= 0))
    stop("center and scale must be finite subject-by-target matrices with positive scales.",
         call. = FALSE)
  vector_ok <- function(x, size) is.numeric(x) && !is.complex(x) &&
    is.null(dim(x)) && length(x) == size && all(is.finite(x))
  if (!vector_ok(q, length(q)) || !length(q) || any(q < 0) || any(diff(q) <= 0))
    stop("q must be finite, nonnegative and strictly increasing.", call. = FALSE)
  if (!vector_ok(subject_index, h) || any(subject_index != floor(subject_index) |
      subject_index < 1 | subject_index > n) || !vector_ok(query_times, h) ||
      any(query_times < 0 | query_times > tau) || !vector_ok(query_state, h) ||
      any(query_state != floor(query_state) | query_state < 1 | query_state > nrow(rewards)) ||
      !vector_ok(weight, h) || !matrix_ok(accrued, h, d))
    stop("Invalid subject indices, times, states, accrued rewards or query weights.", call. = FALSE)
  active <- rewards[query_state, , drop = FALSE] == 1
  expected <- matrix(query_times, h, d)
  if (any(accrued[active] != expected[active])) {
    stop(structure(list(message = "Active monotone rewards must equal query time; use the general conditional method for this history.",
                        call = NULL),
                   class = c("lifeconformal_monotone_history", "error", "condition")))
  }
  if (is.null(flip)) flip <- rep(FALSE, d)
  if (!is.logical(flip) || !is.null(dim(flip)) || length(flip) != d || anyNA(flip))
    stop("flip must be one nonmissing logical value per target.", call. = FALSE)
  if (any(flip) && (is.null(original_center) || is.null(original_accrued)))
    stop("Complemented clocks require original_center and original_accrued.", call. = FALSE)
  if (is.null(original_center)) original_center <- center
  if (is.null(original_accrued)) original_accrued <- accrued
  if (!matrix_ok(original_center, n, d) || !matrix_ok(original_accrued, h, d))
    stop("Original score matrices must match their oriented dimensions.", call. = FALSE)
  if (any(flip) && any(original_accrued[, flip, drop = FALSE][active[, flip, drop = FALSE]] != 0))
    stop("An active complemented clock must have zero original accrued reward.", call. = FALSE)
  if (any(original_center[, !flip, drop = FALSE] != center[, !flip, drop = FALSE]) ||
      any(original_accrued[, !flip, drop = FALSE] != accrued[, !flip, drop = FALSE]) ||
      (any(flip) &&
       (any(center[, flip, drop = FALSE] != tau - original_center[, flip, drop = FALSE]) ||
        any(accrued[, flip, drop = FALSE] !=
            matrix(query_times, h, sum(flip)) - original_accrued[, flip, drop = FALSE]))))
    stop("Original and oriented target arrays are inconsistent.", call. = FALSE)
  l <- length(q)
  supplied <- !vapply(list(determination_time, resolved, membership), is.null, logical(1))
  if (any(supplied) && !all(supplied))
    stop("Supply determination_time, resolved and membership together.", call. = FALSE)
  if (!any(supplied)) {
    determination_time <- matrix(0, n, l)
    resolved <- membership <- matrix(FALSE, n, l)
  } else {
    if (!is.matrix(determination_time) || !is.numeric(determination_time) ||
        is.complex(determination_time) || !identical(dim(determination_time), c(n, l)) ||
        !is.matrix(resolved) || !is.logical(resolved) || anyNA(resolved) ||
        !identical(dim(resolved), c(n, l)) || !is.matrix(membership) ||
        !is.logical(membership) || !identical(dim(membership), c(n, l)) ||
        anyNA(membership[resolved]) || any(!is.finite(determination_time[resolved])) ||
        any(determination_time[resolved] < 0 | determination_time[resolved] > tau))
      stop("Invalid subject-by-threshold structural certificates.", call. = FALSE)
  }
  as_double_matrix <- function(x) { storage.mode(x) <- "double"; x }
  storage.mode(edges) <- "integer"
  spec <- list(jump_times = as.double(jump_times),
               base_increments = as_double_matrix(base_increments), edges = edges,
               relative_risk = as_double_matrix(relative_risk),
               rewards = as_double_matrix(rewards), tau = as.double(tau),
               center = as_double_matrix(center), scale = as_double_matrix(scale), q = as.double(q),
               subject_index = as.integer(subject_index), query_times = as.double(query_times),
               query_state = as.integer(query_state), accrued = as_double_matrix(accrued),
               weight = as.double(weight), determination_time = as_double_matrix(determination_time),
               resolved = resolved, membership = membership, flip = flip,
               original_center = as_double_matrix(original_center),
               original_accrued = as_double_matrix(original_accrued))
  out <- .Call("life_markov_monotone_batch", spec, PACKAGE = "lifeconformal")
  list(contributions = out$contributions,
       diagnostics = list(method = "deterministic_native_monotone_batch",
                          subjects = n, queries = h, jumps = length(jump_times),
                          thresholds = l, adjusted_rows = out$adjusted_rows,
                          maximum_probability_error = 0, tolerance_met = TRUE))
}
