#' Orient monotone binary rewards as nonincreasing clocks
#'
#' A binary target that only increases along fitted edges is complemented.
#' Its occupation reward becomes `1 - reward`, its total becomes `tau - Y`,
#' and its score center becomes `tau - center`. Absolute standardized scores
#' are algebraically unchanged. A numerical evaluator preserving inclusive
#' atom comparisons should evaluate reflected scores in the original target
#' units; reassociating floating-point subtraction can change a boundary by one
#' unit in the last place. Targets that both increase and decrease require the general
#' deterministic distribution engine.
#'
#' @param rewards Numeric state-by-target reward matrix.
#' @param edges Data frame or matrix with one-based `from` and `to` columns.
#' @returns A `life_markov_orientation` containing oriented `rewards` and a
#'   logical target vector `flip`, or NULL when the rewards are not binary or
#'   a target changes in both directions. Constant targets retain their original
#'   orientation. This checks the fitted graph, not structural feasibility.
#' @noRd
orient_markov_rewards <- function(rewards, edges) {
  if (!is.matrix(rewards) || !is.numeric(rewards) || is.complex(rewards) ||
      any(dim(rewards) < 1L) || anyNA(rewards) ||
      any(!rewards %in% c(0, 1))) return(NULL)
  if ((!is.matrix(edges) && !is.data.frame(edges)) || is.null(colnames(edges)) ||
      !all(c("from", "to") %in% colnames(edges)))
    stop("edges must have named from and to columns.", call. = FALSE)
  from <- edges[, "from"]
  to <- edges[, "to"]
  valid <- function(x) is.numeric(x) && !is.complex(x) && is.null(dim(x)) &&
    all(is.finite(x)) && all(x == floor(x) & x >= 1 & x <= nrow(rewards))
  if (!valid(from) || !valid(to) || length(from) != length(to))
    stop("edges must contain valid one-based state indices.", call. = FALSE)
  delta <- rewards[to, , drop = FALSE] - rewards[from, , drop = FALSE]
  increases <- colSums(delta > 0) > 0
  decreases <- colSums(delta < 0) > 0
  if (any(increases & decreases)) return(NULL)
  flip <- increases
  names(flip) <- colnames(rewards)
  oriented <- rewards
  if (any(flip)) oriented[, flip] <- 1 - oriented[, flip, drop = FALSE]
  structure(list(rewards = oriented, flip = flip),
            class = "life_markov_orientation")
}

#' Apply a binary reward orientation to score and history arrays
#'
#' @param orientation Object returned by orient_markov_rewards().
#' @param center Finite baseline-profile-by-target score centers.
#' @param accrued Optional finite query-by-target accumulated rewards.
#' @param time Conditioning time per accrued row. A common scalar is accepted.
#' @param tau Positive common horizon.
#' @returns Oriented rewards, centers, accrued rewards and the flip vector,
#'   together with `original_center` and `original_accrued` for canonical score
#'   evaluation at floating-point atom boundaries.
#'   Baseline and query matrices may have different row counts. Scale values
#'   and score thresholds are unchanged and are not copied into this result.
#' @details
#' An increasing target's accrued complement is `time - accrued`, because the
#'   target clock began at zero. This operation does not certify that an observed
#'   path follows the working graph. Verify the observed reward path separately
#'   before using the exact native clock eligibility condition.
#' @noRd
orient_markov_query <- function(orientation, center, accrued = NULL,
                                 time = NULL, tau) {
  if (!is.list(orientation) || !is.matrix(orientation$rewards) ||
      !is.numeric(orientation$rewards) || is.complex(orientation$rewards) ||
      any(dim(orientation$rewards) < 1L) || anyNA(orientation$rewards) ||
      any(!orientation$rewards %in% c(0, 1)) ||
      !is.logical(orientation$flip) || anyNA(orientation$flip) ||
      length(orientation$flip) != ncol(orientation$rewards))
    stop("Supply a valid binary reward orientation.", call. = FALSE)
  d <- ncol(orientation$rewards)
  matrix_ok <- function(x) is.matrix(x) && is.numeric(x) && !is.complex(x) &&
    ncol(x) == d && all(is.finite(x))
  if (!matrix_ok(center) || !is.numeric(tau) || is.complex(tau) ||
      length(tau) != 1L || !is.finite(tau) || tau <= 0)
    stop("Supply finite score centers and a positive finite horizon.", call. = FALSE)
  original_center <- center
  original_accrued <- accrued
  flip <- orientation$flip
  if (any(flip)) center[, flip] <- tau - center[, flip, drop = FALSE]
  if (!is.null(accrued)) {
    if (!matrix_ok(accrued) || !is.numeric(time) || is.complex(time) ||
        !is.null(dim(time)) || !(length(time) %in% c(1L, nrow(accrued))) ||
        any(!is.finite(time)) || any(time < 0 | time > tau))
      stop("Supply finite accrued rewards and matching conditioning times in [0,tau].",
           call. = FALSE)
    if (length(time) == 1L) time <- rep(time, nrow(accrued))
    if (any(flip)) {
      accrued[, flip] <- matrix(time, nrow(accrued), sum(flip)) -
        accrued[, flip, drop = FALSE]
    }
  } else if (!is.null(time)) {
    stop("time is used only when accrued rewards are supplied.", call. = FALSE)
  }
  if (any(!is.finite(center)) || (!is.null(accrued) && any(!is.finite(accrued))))
    stop("Oriented reward values overflow; check time units.", call. = FALSE)
  list(rewards = orientation$rewards, flip = flip, center = center,
       accrued = accrued, original_center = original_center,
       original_accrued = original_accrued)
}
