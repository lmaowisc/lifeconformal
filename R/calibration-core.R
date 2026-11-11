#' Compute stopped inverse-probability calibration contributions
#'
#' This reference kernel implements `R * I / G(X_q-)` for one candidate
#' threshold. It does not evaluate or require conditional coverage models.
#'
#' @param observed_membership Logical or zero/one vector giving observed region
#'   membership. Missing values are allowed only where `resolved` is false.
#' @param resolved Logical vector indicating that coverage was determined within
#'   the observed history under the declared observation convention.
#' @param g_left Numeric vector of inclusive observation probabilities
#'   `G(X_q-)`, in `(0, 1]`, in the same subject order.
#' @returns A numeric contribution vector, in input order. Unresolved entries
#'   are zero, without multiplying missing membership values by zero.
#' @details The caller supplies valid coverage-determining times and their
#'   observed-determination indicators. This kernel cannot infer these from
#'   probabilities. Its time and working
#'   memory costs are linear in the number of subjects.
#' @noRd
compute_ipw_contributions <- function(observed_membership, resolved, g_left) {
  n <- .calibration_check_endpoints(observed_membership, resolved, g_left)
  contribution <- numeric(n)
  contribution[resolved] <- observed_membership[resolved] / g_left[resolved]
  names(contribution) <- names(observed_membership)
  contribution
}

#' Compute the stopped step-censoring augmented contribution
#'
#' Uses conditional coverage given the history through the evaluation time.
#'
#' @inheritParams compute_ipw_contributions
#' @param stop_time Nonnegative finite vector of stopped observation times
#'   `X_q`, in the same subject order.
#' @param m_stop Numeric vector of conditional coverage probabilities at the
#'   stopped time. Missing values are permitted for resolved subjects, whose
#'   exact membership is used instead.
#' @param jumps Either `NULL` for no censoring jumps or a data frame with columns
#'   `subject_index`, `time`, `g_left`, `g_right`, and `m_plus`. Subject indices
#'   refer to the endpoint-vector order. The survival columns give the two
#'   limits of the same fitted survival function; `m_plus` is conditional
#'   coverage given history through that time. Include all censoring jumps
#'   strictly before each subject's stopped time. Rows at or after the stopped
#'   time are ignored, so in particular a jump at the horizon is excluded.
#' @returns A numeric vector, in subject order, of
#'   `[R I + (1-R) m_stop] / G(X_q-)` minus the sum of
#'   `m_plus * [1/G(u) - 1/G(u-)]` over `u < X_q`. Values are not clipped.
#' @details Survival values must come from one coherent fitted function,
#'   including any positivity floor applied upstream. This kernel verifies
#'   each jump's direction and range; construction of a globally coherent
#'   survival function is the censoring adapter's responsibility. The caller
#'   applies exact membership overrides after coverage determination. The
#'   arithmetic uses one pass through endpoint vectors and jump rows; working
#'   memory is linear in their combined length. It does not fit models or
#'   perform conditional-distribution prediction.
#' @noRd
compute_aipw_contributions <- function(observed_membership, resolved,
                                       stop_time, g_left, m_stop,
                                       jumps = NULL) {
  n <- .calibration_check_endpoints(observed_membership, resolved, g_left)
  .calibration_check_numeric(stop_time, "stop_time", n)
  if (any(stop_time < 0)) {
    stop("`stop_time` must be nonnegative.", call. = FALSE)
  }
  .calibration_check_probability(m_stop, "m_stop", n, allow_na = TRUE)
  if (anyNA(m_stop[!resolved])) {
    stop("`m_stop` is required for every unresolved subject.", call. = FALSE)
  }

  numerator <- as.numeric(observed_membership)
  numerator[!resolved] <- m_stop[!resolved]
  contribution <- numerator / g_left

  if (!is.null(jumps)) {
    required <- c("subject_index", "time", "g_left", "g_right", "m_plus")
    if (!is.data.frame(jumps) || anyDuplicated(names(jumps)) ||
        !all(required %in% names(jumps))) {
      stop("`jumps` must be a data frame with unique columns including ",
           "subject_index, time, g_left, g_right, and m_plus.", call. = FALSE)
    }
    nj <- nrow(jumps)
    .calibration_check_numeric(jumps$subject_index, "jumps$subject_index", nj)
    index <- jumps$subject_index
    if (any(index != floor(index) | index < 1 | index > n)) {
      stop("`jumps$subject_index` must index the endpoint vectors.",
           call. = FALSE)
    }
    .calibration_check_numeric(jumps$time, "jumps$time", nj)
    if (any(jumps$time < 0)) {
      stop("`jumps$time` must be nonnegative.", call. = FALSE)
    }

    # Endpoint ties are excluded before requiring the nuisance values. No
    # conditional coverage prediction is needed for these unused rows.
    keep <- jumps$time < stop_time[index]
    if (any(keep)) {
      used <- jumps[keep, required, drop = FALSE]
      index <- used$subject_index
      if (anyDuplicated(used[c("subject_index", "time")])) {
        stop("`jumps` contains duplicate subject/time rows before stopping.",
             call. = FALSE)
      }
      .calibration_check_survival(used$g_left, "jumps$g_left", nrow(used))
      .calibration_check_survival(used$g_right, "jumps$g_right", nrow(used))
      .calibration_check_probability(used$m_plus, "jumps$m_plus", nrow(used))
      if (any(used$g_right > used$g_left)) {
        stop("Censoring survival must not increase at a jump.", call. = FALSE)
      }

      # This form avoids cancellation of two nearly equal inverse values and
      # avoids forming g_left * g_right, which may underflow unnecessarily.
      inverse_increment <- (used$g_left - used$g_right) /
        used$g_left / used$g_right
      weighted_increment <- used$m_plus * inverse_increment
      correction <- rowsum(weighted_increment, index, reorder = FALSE)
      subject_index <- as.integer(rownames(correction))
      contribution[subject_index] <- contribution[subject_index] -
        correction[, 1L]
    }
  }

  if (any(!is.finite(contribution))) {
    stop("Augmented contributions overflowed; check censoring positivity ",
         "and jump coherence.", call. = FALSE)
  }
  names(contribution) <- names(observed_membership)
  contribution
}

#' Invert an ordered empirical coverage curve
#'
#' Selects the first threshold at which estimated coverage reaches the cutoff.
#' The raw curve may be nonmonotone for either calibration method.
#'
#' @param q Strictly increasing finite numeric vector of candidate thresholds.
#' @param estimate Finite numeric vector of raw empirical coverage estimates,
#'   one per threshold. Signed augmented values are retained unchanged.
#' @param n Positive integer number of independent calibration subjects.
#' @param alpha Desired miscoverage probability, strictly between zero and one.
#' @param estimate_sum Optional unscaled contribution sums, one per candidate.
#'   When supplied, compare these directly with the integer rank; the displayed
#'   mean estimates remain unchanged. NULL retains the mean-based internal API.
#' @returns A list with selected `threshold`, one-based `index`, `full_support`
#'   flag, `reason`, conformal `rank`, `criterion`, `n`, `alpha`, and the raw
#'   `curve` data frame. Full support uses threshold `Inf` and index `NA`.
#' @details Set `k = ceiling((n + 1) * (1 - alpha))`. If `k = n + 1`, return
#'   full support. Otherwise select the first candidate with estimate at least
#'   `k/n`; if none exists, return full support. An empty curve is allowed only
#'   when the rank already requires full support. The scan uses linear time;
#'   the raw input curve is preserved without isotonic fitting or bisection.
#' @noRd
invert_coverage_curve <- function(q, estimate, n, alpha, estimate_sum = NULL) {
  .calibration_check_numeric(n, "n", 1L)
  if (n < 1 || n != floor(n) || n > .Machine$integer.max) {
    stop("`n` must be a positive integer calibration sample size.",
         call. = FALSE)
  }
  .calibration_check_numeric(alpha, "alpha", 1L)
  if (alpha <= 0 || alpha >= 1) {
    stop("`alpha` must lie strictly between zero and one.", call. = FALSE)
  }
  .calibration_check_numeric(q, "q")
  .calibration_check_numeric(estimate, "estimate", length(q))
  if (!is.null(estimate_sum))
    .calibration_check_numeric(estimate_sum, "estimate_sum", length(q))
  if (length(q) > 1L && any(diff(q) <= 0)) {
    stop("`q` must be strictly increasing.", call. = FALSE)
  }

  rank <- ceiling((as.double(n) + 1) * (1 - alpha))
  criterion <- rank / n
  full_support <- rank > n
  reason <- "rank_exceeds_sample"
  index <- NA_integer_
  if (!full_support) {
    if (!length(q)) {
      stop("A nonempty threshold grid is required when the rank is at most n.",
           call. = FALSE)
    }
    crossing <- if (is.null(estimate_sum)) which(estimate >= criterion) else
      which(estimate_sum >= rank)
    full_support <- !length(crossing)
    reason <- if (full_support) "no_crossing" else "crossing"
    if (!full_support) {
      index <- crossing[1L]
    }
  }

  list(
    threshold = if (full_support) Inf else q[index],
    index = index,
    full_support = full_support,
    reason = reason,
    rank = rank,
    criterion = criterion,
    n = as.integer(n),
    alpha = alpha,
    curve = data.frame(q = q, estimate = estimate)
  )
}

#' Validate aligned contribution endpoint inputs
#' @noRd
.calibration_check_endpoints <- function(observed_membership, resolved,
                                         g_left) {
  if (!is.logical(resolved) || !is.null(dim(resolved)) ||
      !length(resolved) || anyNA(resolved)) {
    stop("`resolved` must be a nonempty logical vector without missing values.",
         call. = FALSE)
  }
  n <- length(resolved)
  if (!(is.logical(observed_membership) || is.numeric(observed_membership)) ||
      is.complex(observed_membership) ||
      !is.null(dim(observed_membership)) || length(observed_membership) != n) {
    stop("`observed_membership` must be a logical or zero/one vector aligned ",
         "with `resolved`.", call. = FALSE)
  }
  known <- !is.na(observed_membership)
  if (any(!observed_membership[known] %in% c(0, 1)) ||
      anyNA(observed_membership[resolved])) {
    stop("Membership must be zero or one, and observed for resolved subjects.",
         call. = FALSE)
  }
  .calibration_check_survival(g_left, "g_left", n)
  n
}

#' Validate a numeric vector without recycling
#' @noRd
.calibration_check_numeric <- function(x, label, n = NULL, allow_na = FALSE) {
  valid <- is.numeric(x) && !is.complex(x) && is.null(dim(x)) &&
    (is.null(n) || length(x) == n)
  if (!valid || any(!is.finite(x) & !(allow_na & is.na(x)))) {
    stop("`", label, "` must be a finite numeric vector",
         if (is.null(n)) "." else paste0(" of length ", n, "."),
         call. = FALSE)
  }
  invisible(NULL)
}

#' Validate probability values for the calibration kernel
#' @noRd
.calibration_check_probability <- function(x, label, n, allow_na = FALSE) {
  .calibration_check_numeric(x, label, n, allow_na)
  if (any(x < 0 | x > 1, na.rm = allow_na)) {
    stop("`", label, "` must lie in [0, 1].", call. = FALSE)
  }
  invisible(NULL)
}

#' Validate positive survival and representable inverse weights
#' @noRd
.calibration_check_survival <- function(x, label, n) {
  .calibration_check_numeric(x, label, n)
  if (any(x <= 0 | x > 1) || any(!is.finite(1 / x))) {
    stop("`", label, "` must lie in (0, 1] with finite inverse weights.",
         call. = FALSE)
  }
  invisible(NULL)
}
