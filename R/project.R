#' Project a joint prediction region onto linear summaries
#'
#' Obtain prediction intervals for sums, differences, or other linear
#' combinations of the calibrated target vector. The calculation retains the
#' shared occupation-time constraints, including overlapping targets.
#'
#' @param object A `lifeconformal_prediction` returned by [predict()].
#' @param weights A finite named numeric vector for one linear combination, or
#'   a numeric matrix whose rows define combinations and whose named columns
#'   identify targets. Omitted targets receive coefficient zero. Matrix row
#'   names name the derived summaries; unnamed rows use `summary1`, `summary2`,
#'   and so on. A vector produces one summary named `summary1`. Positive,
#'   negative, and zero coefficients are permitted. Unknown or duplicated target
#'   names are rejected; named columns are reordered to the trained target order.
#'
#' @returns A data frame with columns `id`, `summary`, `lower`, `upper`, and
#'   `status`, in prediction-row then summary order. Attributes include `weights`
#'   (the full aligned coefficient matrix), `alpha`, `tau`,
#'   `inherits_joint_coverage = TRUE`, `conservative_outer_support`, and numerical
#'   diagnostics. An empty region has status `"empty_region"`, missing endpoints,
#'   and emits a warning of class `lifeconformal_empty_projection`.
#'
#' @details
#' For occupation fractions `v`, this function projects the stored closed
#' region defined by `v >= 0`, `sum(v) = 1`, and
#' `lower / tau <= A %*% v <= upper / tau`. Only states structurally reachable
#' from the subject's initial state are retained. `A` maps state occupations to
#' the original targets. Minimizing and maximizing each requested linear
#' combination over this common region respects its joint constraints.
#'
#' If the original target vector belongs to the joint prediction region, all
#' requested linear summaries belong to these intervals simultaneously. No new
#' calibration or multiple-testing adjustment is needed. This implication
#' inherits the original region's coverage assumptions and nominal `alpha`;
#' projection does not strengthen its statistical guarantee.
#'
#' Structural reachability supplies a conservative closed support envelope.
#' Consequently the result projects that implemented envelope, not necessarily
#' the exact set of attainable paths. Open path endpoints are enclosed.
#'
#' Nontrivial intersections use the suggested package `lpSolve`; install it to
#' use those projections. Full-support and single-point regions are evaluated
#' analytically. Occupations are scaled to fractions before optimization. Dual
#' bounds are corrected for residual constraint violations and widened for
#' floating-point arithmetic, including signed coefficient sums and time-unit
#' conversions. Ill-conditioned combinations may require conservative widening.
#' Solver failure is reported instead of replacing
#' the joint projection by a sum of marginal intervals.
#'
#' @seealso [predict.lifeconformal_fit()], [life_calibrate()]
#' @export
#' @examples
#' if (requireNamespace("lpSolve", quietly = TRUE)) {
#'   dat <- life_example_data(80, seed = 11)
#'   fit <- life_train(LifeHist(time, status) ~ 1, dat[1:40, ],
#'     absorbing = 2, targets = list(alive = 0, dead = 2), tau = 5)
#'   cal <- life_calibrate(fit, dat[41:80, ], method = "ipw")
#'   region <- predict(cal, data.frame(id = "new"))
#'   project_life(region, c(alive = 1, dead = 1))
#'   contrasts <- rbind(alive_minus_dead = c(alive = 1, dead = -1),
#'                      total = c(alive = 1, dead = 1))
#'   project_life(region, contrasts)
#' }
project_life <- function(object, weights) {
  specification <- prepare_life_projection(object, weights)
  s <- specification
  n <- s$n
  J <- nrow(s$weights)
  answer <- data.frame(id = rep(s$id, each = J),
                       summary = rep(rownames(s$weights), times = n),
                       lower = rep(NA_real_, n * J), upper = rep(NA_real_, n * J),
                       status = rep("ok", n * J), stringsAsFactors = FALSE)
  cache <- new.env(parent = emptyenv())
  templates <- vector("list", length(s$reachable))
  calls <- 0L
  unique_regions <- 0L
  max_residual <- 0
  max_padding <- 0
  max_scaling_padding <- 0
  for (i in seq_len(n)) {
    initial <- s$initial[i]
    if (is.null(templates[[initial]])) {
      A <- s$A[, s$reachable[[initial]], drop = FALSE]
      # States with the same target incidence yield the same projected
      # occupation variable. Aggregate them before constructing any LP.
      A <- A[, !duplicated(t(A)), drop = FALSE]
      coefficients <- projection_coefficients(s$weights, A)
      templates[[initial]] <- c(list(A = A), coefficients, list(
                                    minimum = apply(A, 1L, min),
                                    maximum = apply(A, 1L, max)))
    }
    template <- templates[[initial]]
    lower <- projection_round_scale(s$lower[i, ], s$tau, "lower", divide = TRUE)
    upper <- projection_round_scale(s$upper[i, ], s$tau, "upper", divide = TRUE)
    key <- paste(c(initial, sprintf("%a", s$lower[i, ]),
                   sprintf("%a", s$upper[i, ])), collapse = "|")
    if (exists(key, envir = cache, inherits = FALSE)) {
      projected <- get(key, envir = cache, inherits = FALSE)
    } else {
      projected <- if (any(s$lower[i, ] > s$upper[i, ]))
        empty_life_projection(J) else
        project_occupation_polytope(template, lower, upper)
      assign(key, projected, envir = cache)
      calls <- calls + projected$solver_calls
      unique_regions <- unique_regions + 1L
      max_residual <- max(max_residual, projected$dual_residual)
      max_padding <- max(max_padding, projected$roundoff_padding)
    }
    rows <- (i - 1L) * J + seq_len(J)
    if (projected$empty) answer$status[rows] <- "empty_region" else {
      answer$lower[rows] <- projection_round_scale(projected$lower, s$tau, "lower")
      answer$upper[rows] <- projection_round_scale(projected$upper, s$tau, "upper")
      max_scaling_padding <- max(max_scaling_padding,
        abs(answer$lower[rows] - s$tau * projected$lower),
        abs(answer$upper[rows] - s$tau * projected$upper))
    }
  }
  if (any(!is.finite(answer$lower[answer$status == "ok"])) ||
      any(!is.finite(answer$upper[answer$status == "ok"]))) {
    stop("The projected endpoints overflow; rescale the supplied weights.", call. = FALSE)
  }
  attributes(answer)[c("weights", "alpha", "tau", "inherits_joint_coverage",
                       "conservative_outer_support", "diagnostics")] <-
    list(s$weights, object$alpha, s$tau, TRUE, isTRUE(s$outer),
         list(solver_calls = calls, unique_regions = unique_regions,
              max_dual_residual = max_residual,
              max_roundoff_padding = s$tau * max_padding + max_scaling_padding,
              support = "closed_reachable_state_simplex_intersect_box"))
  if (any(answer$status == "empty_region")) {
    ids <- unique(as.character(answer$id[answer$status == "empty_region"]))
    warning(structure(list(message = paste0("The joint prediction region is empty for ",
      length(ids), " subject label(s): ", paste(utils::head(ids, 5L), collapse = ", "),
      ". Derived intervals have missing endpoints and status 'empty_region'."),
      call = NULL), class = c("lifeconformal_empty_projection", "warning", "condition")))
  }
  answer
}

#' Validate and align a joint prediction region and projection weights
#' @noRd
prepare_life_projection <- function(object, weights) {
  if (!inherits(object, "lifeconformal_prediction") ||
      !is.matrix(object$center) || !is.numeric(object$center)) {
    stop("object must be a lifeconformal prediction region.", call. = FALSE)
  }
  target <- colnames(object$center)
  n <- nrow(object$center)
  d <- ncol(object$center)
  if (!d || is.null(target) || anyNA(target) || any(!nzchar(target)) || anyDuplicated(target))
    stop("The prediction object has an invalid target dictionary.", call. = FALSE)
  check_names <- function(x) !is.null(x) && !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x)
  if (is.numeric(weights) && !is.complex(weights) && is.null(dim(weights))) {
    if (!length(weights) || !check_names(names(weights)))
      stop("Vector weights require unique nonempty target names.", call. = FALSE)
    weights <- matrix(weights, 1L, dimnames = list("summary1", names(weights)))
  }
  if (!is.matrix(weights) || !is.numeric(weights) || is.complex(weights) ||
      !nrow(weights) || !ncol(weights) || any(!is.finite(weights)) ||
      !check_names(colnames(weights)))
    stop("weights must be a finite numeric vector or matrix with named target columns.", call. = FALSE)
  if (any(!colnames(weights) %in% target))
    stop("Unknown target in weights: ", paste(setdiff(colnames(weights), target), collapse = ", "),
         ".", call. = FALSE)
  if (is.null(rownames(weights))) rownames(weights) <- paste0("summary", seq_len(nrow(weights)))
  if (!check_names(rownames(weights)))
    stop("Summary row names must be unique and nonempty.", call. = FALSE)
  aligned <- matrix(0, nrow(weights), d, dimnames = list(rownames(weights), target))
  aligned[, match(colnames(weights), target)] <- weights
  support <- object$region$support
  A <- support$joint_support$A
  tau <- object$tau
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0 ||
      !is.matrix(A) || !is.numeric(A) || any(!is.finite(A)) ||
      any(!A %in% c(0, 1)) || nrow(A) != d ||
      !ncol(A) || !check_names(rownames(A)) || !setequal(rownames(A), target))
    stop("The prediction object lacks a valid joint occupation support.", call. = FALSE)
  A <- A[match(target, rownames(A)), , drop = FALSE]
  if (!identical(colnames(A), as.character(support$states)))
    stop("The joint support state columns differ from its state dictionary.", call. = FALSE)
  if (!is.numeric(object$alpha) || length(object$alpha) != 1L ||
      !is.finite(object$alpha) || object$alpha <= 0 || object$alpha >= 1)
    stop("The prediction object has an invalid nominal miscoverage level.", call. = FALSE)
  initial <- object$region$initial_state
  reachable <- support$reachable
  if (!is.numeric(initial) || length(initial) != n || anyNA(initial) ||
      any(initial != floor(initial)) || any(initial < 1 | initial > ncol(A)) ||
      !is.list(reachable) || length(reachable) != ncol(A))
    stop("The prediction object has invalid initial-state reachability.", call. = FALSE)
  for (k in unique(initial)) {
    r <- reachable[[k]]
    if (!is.numeric(r) || !length(r) || anyNA(r) || any(r != floor(r)) ||
        any(r < 1 | r > ncol(A)) || anyDuplicated(r))
      stop("The prediction object has invalid reachable-state indices.", call. = FALSE)
  }
  b <- object$bounds
  if (!is.data.frame(b) || !all(c("id", "target", "lower", "upper") %in% names(b)) ||
      nrow(b) != n * d || !is.numeric(b$lower) || !is.numeric(b$upper) ||
      any(!is.finite(b$lower)) || any(!is.finite(b$upper)))
    stop("The prediction object has invalid coordinate bounds.", call. = FALSE)
  lower <- upper <- matrix(0, n, d, dimnames = list(NULL, target))
  id <- b$id[if (n) seq.int(1L, n * d, by = d) else integer()]
  for (i in seq_len(n)) {
    rows <- (i - 1L) * d + seq_len(d)
    if (anyNA(b$target[rows]) || anyDuplicated(b$target[rows]) ||
        !setequal(b$target[rows], target) ||
        !identical(as.character(b$id[rows]), rep(as.character(id[i]), d)))
      stop("Prediction bounds must retain one complete target block per subject.", call. = FALSE)
    order <- match(target, b$target[rows])
    lower[i, ] <- b$lower[rows[order]]
    upper[i, ] <- b$upper[rows[order]]
  }
  list(n = n, weights = aligned, A = A, tau = tau, lower = lower, upper = upper,
       id = id, reachable = reachable, initial = initial, outer = support$outer)
}

#' Project a closed occupation polytope without discarding shared constraints
#' @noRd
project_occupation_polytope <- function(template, lower, upper) {
  A <- template$A
  coefficients <- template$coefficients
  J <- nrow(coefficients)
  empty_result <- function(calls = 0L) empty_life_projection(J, calls)
  if (any(lower > upper) || any(lower > template$maximum) ||
      any(upper < template$minimum)) return(empty_result())
  minimum <- apply(template$coefficient_lower, 1L, min)
  maximum <- apply(template$coefficient_upper, 1L, max)
  if (all(lower <= template$minimum) && all(upper >= template$maximum)) {
    return(list(lower = minimum, upper = maximum, empty = FALSE,
                 solver_calls = 0L, dual_residual = 0,
                 roundoff_padding = template$coefficient_padding))
  }
  if (!requireNamespace("lpSolve", quietly = TRUE))
    stop("This joint-region projection requires package 'lpSolve'; install it with install.packages('lpSolve').",
         call. = FALSE)
  K <- ncol(A)
  # All constraints use <=. Two rows express the simplex equality; this
  # makes dual multipliers nonnegative and permits residual correction.
  M <- rbind(A, -A, rep(1, K), rep(-1, K))
  rhs <- c(upper, -lower, 1, -1)
  primal <- lpSolve::lp("min", rep(0, K), M, rep("<=", nrow(M)), rhs)
  if (primal$status == 2L) return(empty_result(1L))
  if (primal$status != 0L)
    stop("The projection feasibility solver failed (status ", primal$status, ").", call. = FALSE)
  residual <- max(0, M %*% primal$solution - rhs, -primal$solution)
  if (!is.finite(residual) || residual > 1e-7)
    stop("The projection solver returned an insufficiently feasible region point.", call. = FALSE)
  result_lower <- result_upper <- numeric(J)
  calls <- 1L
  max_residual <- 0
  max_padding <- template$coefficient_padding
  for (j in seq_len(J)) {
    c <- coefficients[j, ]
    # If the rounded coefficients coincide, their outward envelopes still
    # contain every exact signed sum. Feasibility plus those extrema suffices.
    if (all(c == c[1L])) {
      result_lower[j] <- minimum[j]
      result_upper[j] <- maximum[j]
      next
    }
    lo <- projection_dual_bound(-template$coefficient_lower[j, ], M, rhs)
    hi <- projection_dual_bound(template$coefficient_upper[j, ], M, rhs)
    calls <- calls + 2L
    result_lower[j] <- max(minimum[j], -lo$value)
    result_upper[j] <- min(maximum[j], hi$value)
    max_residual <- max(max_residual, lo$residual, hi$residual)
    max_padding <- max(max_padding, lo$padding, hi$padding)
  }
  if (any(result_lower > result_upper)) return(empty_result(calls))
  list(lower = result_lower, upper = result_upper, empty = FALSE,
       solver_calls = calls, dual_residual = max_residual, roundoff_padding = max_padding)
}

#' Obtain a conservative upper objective bound from LP dual multipliers
#' @noRd
projection_dual_bound <- function(objective, M, rhs) {
  dual <- lpSolve::lp("min", rhs, t(M), rep(">=", ncol(M)), objective)
  if (dual$status != 0L)
    stop("The projection bound solver failed (status ", dual$status, ").", call. = FALSE)
  lambda <- pmax(0, dual$solution)
  product <- as.vector(crossprod(M, lambda))
  residual <- max(0, objective - product)
  value <- sum(rhs * lambda) + residual
  # For feasible simplex v, objective'v <= rhs'lambda +
  # max(0, objective - M'lambda). The correction is valid even if the solver's
  # dual coefficients have a small feasibility residual. Add ordinary dot-product
  # roundoff allowances before returning the upper endpoint.
  eps <- .Machine$double.eps
  operations <- nrow(M) + 8
  gamma <- operations * eps / (1 - operations * eps)
  padding <- gamma * (1 + sum(abs(rhs * lambda)) + max(abs(objective)) +
    max(colSums(abs(M) * lambda)))
  if (!is.finite(value) || !is.finite(padding))
    stop("Numerical overflow in the projection dual bound; rescale the weights.", call. = FALSE)
  list(value = value + padding, residual = residual, padding = padding)
}

#' Return an explicit empty projection without a reversed numeric interval
#' @noRd
empty_life_projection <- function(n, calls = 0L) {
  list(lower = rep(NA_real_, n), upper = rep(NA_real_, n), empty = TRUE,
       solver_calls = calls, dual_residual = 0, roundoff_padding = 0)
}

#' Enclose one rounded time-unit conversion in the required outward direction
#' @noRd
projection_round_scale <- function(x, tau, side, divide = FALSE) {
  result <- if (divide) x / tau else x * tau
  if (any(!is.finite(result)))
    stop("The projected endpoints overflow; rescale the supplied weights.", call. = FALSE)
  exact <- x == 0 | tau == 1 | (divide & x == tau)
  if (side == "lower") result[!exact] <- -cdt_nextafter(-result[!exact])
  else result[!exact] <- cdt_nextafter(result[!exact])
  result
}

#' Preserve signed occupation coefficients with exact floating expansions
#' @noRd
projection_coefficients <- function(weights, A) {
  value <- lower <- upper <- matrix(0, nrow(weights), ncol(A))
  for (k in seq_len(ncol(A))) {
    selected <- which(A[, k] != 0)
    for (j in seq_len(nrow(weights))) {
      coefficient <- projection_sum_interval(weights[j, selected])
      value[j, k] <- coefficient[1L]
      lower[j, k] <- coefficient[2L]
      upper[j, k] <- coefficient[3L]
    }
  }
  list(coefficients = value, coefficient_lower = lower, coefficient_upper = upper,
       coefficient_padding = max(value - lower, upper - value))
}

#' Add finite doubles with an exact residual, without extended precision
#' @noRd
projection_two_sum <- function(a, b) {
  value <- a + b
  if (!is.finite(value))
    stop("The linear combination overflows; rescale the supplied weights.", call. = FALSE)
  b_virtual <- value - a
  error <- (a - (value - b_virtual)) + (b - b_virtual)
  c(value, error)
}

#' Enclose the mathematical sum of the supplied floating-point coefficients
#' @noRd
projection_sum_interval <- function(x) {
  # Nonoverlapping expansion terms retain cancellation residuals exactly.
  # FastTwoSum is error-free for |a| >= |b| with IEEE gradual underflow;
  # overflow is rejected. No platform-specific long-double sum is assumed.
  partials <- numeric(length(x))
  used <- 0L
  for (a in x) {
    next_used <- 0L
    for (index in seq_len(used)) {
      b <- partials[index]
      if (abs(a) < abs(b)) {
        temporary <- a
        a <- b
        b <- temporary
      }
      value <- a + b
      if (!is.finite(value))
        stop("The linear combination overflows; rescale the supplied weights.", call. = FALSE)
      residual <- b - (value - a)
      if (residual != 0) {
        next_used <- next_used + 1L
        partials[next_used] <- residual
      }
      a <- value
    }
    used <- next_used + 1L
    partials[used] <- a
  }
  value <- lower <- upper <- 0
  for (term in partials[seq_len(used)]) {
    lo <- projection_two_sum(lower, term)
    hi <- projection_two_sum(upper, term)
    lower <- if (lo[2L] < 0) -cdt_nextafter(-lo[1L]) else lo[1L]
    upper <- if (hi[2L] > 0) cdt_nextafter(hi[1L]) else hi[1L]
    value <- value + term
  }
  c(value, lower, upper)
}
