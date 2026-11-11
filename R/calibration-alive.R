#' Eligibility for conditioning a nested occupation time on time alive
#' @noRd
alive_occupation_supported <- function(rewards, absorbing) {
  if (!identical(dim(rewards), c(3L,2L)) || length(absorbing)!=1L ||
      any(!rewards %in% c(0,1))) return(FALSE)
  key <- rewards[,1L] + 2L*rewards[,2L]
  identical(sort(unname(as.integer(key))), c(0L,2L,3L)) && key[absorbing]==0
}

#' Bound a two-clock event by the scalar home clock under a terminal-time gate
#'
#' The first reward identifies one living state. The second covers both living
#' states until death. A live history has accrued its calendar time in the second
#' clock. The accepted death dates and survival through the horizon therefore
#' depend only on the second score coordinate. A backward recursion conditions
#' the fitted stochastic kernels on those terminal dates. The existing scalar
#' enclosure for the first reward then gives a joint enclosure on multiplication
#' by the conditioning probability. On a large candidate block, a cheaper
#' screen intersects marginal Frechet bounds with the maximum probability of
#' a feasible death outside the required alive interval.
#' @noRd
alive_conditional_coverage <- function(args, baseline, screen) {
  original <- args
  fallback <- function(reason) {
    out <- do.call(baseline, original)
    out$diagnostics$alive_fallback <- reason
    out
  }
  h <- length(args$query_times)
  q <- args$q
  times <- args$jump_times
  factor <- args$kernel_factors
  if (h > 1000000/length(q) || !alive_occupation_supported(args$rewards,factor$absorbing))
    return(fallback("unsupported rewards or size"))
  # The first occurrence of a state in a training sample need not follow its
  # clinical label. Re-index the two transient states and death jointly.
  order_states <- match(c(3,2,0), args$rewards[,1L]+2*args$rewards[,2L])
  args$rewards <- args$rewards[order_states,,drop=FALSE]
  args$query_state <- match(args$query_state,order_states)
  factor$from <- as.double(match(factor$from,order_states))
  factor$to <- as.double(match(factor$to,order_states))
  factor$absorbing <- 3
  args$kernel_factors <- factor
  live <- args$query_state != 3L
  if (any(args$accrued[live,2] != args$query_times[live]) ||
      any(args$center[,2] != args$center[1,2]) || any(args$scale[,2] != args$scale[1,2]))
    return(fallback("unsupported profile"))
  active <- rowSums(factor$base_increments[,factor$to == 3,drop=FALSE]) > 0
  dates <- c(times, args$tau)
  terminal <- c(active, TRUE)
  distance <- abs((dates - args$center[1,2]) / args$scale[1,2])
  eps <- 128 * .Machine$double.eps * (length(times) + 2)
  # The general joint recursion accumulates alive time in binary64. Within a
  # jump-count-scaled neighborhood of an atom, use that original tie convention.
  proximity <- eps * (1 + (args$tau + abs(args$center[1,2])) / args$scale[1,2] + max(q))
  if (!all(is.finite(c(distance, proximity))) ||
      any(vapply(q,function(z) any(abs(distance[terminal] - z) <= proximity), logical(1))))
    return(fallback("terminal atom boundary"))
  native <- list(jump_times = as.double(times), query_times = as.double(args$query_times),
    query_state = as.double(args$query_state), tau = as.double(args$tau),
    kernel_factors = factor)
  scalar <- args
  scalar$query_known <- scalar$query_weights <- NULL
  scalar$strict <- FALSE
  scalar$max_refinements <- 0L
  scalar$exact_limit <- 32L
  scalar$point_mass_threshold <- .1
  coordinate <- function(j) {
    z <- scalar
    for (name in c("rewards", "accrued", "center", "scale"))
      z[[name]] <- args[[name]][, j, drop = FALSE]
    z
  }
  home <- coordinate(1L)
  value <- tryCatch({
    if (screen) {
      first <- do.call(markov_conditional_coverage_grid, home)
      second <- do.call(markov_conditional_coverage_grid, coordinate(2L))
      spec <- c(native, args[c("center","scale","accrued","q")])
      spec$q <- as.double(q)
      missed <- .Call("life_markov_alive_conditional", spec, FALSE, PACKAGE="lifeconformal")
      if (is.null(missed)) return(fallback("terminal mass error"))
      low <- pmax(first$lower + second$lower - 1, first$lower - missed)
      high <- pmin(first$upper, second$upper)
      point <- (low + high)/2
    } else {
      low <- point <- high <- matrix(0, h, length(q))
      allowed <- outer(distance, q, "<=")
      # Acceptance masks are nested. Equal numbers of allowed active atoms
      # therefore identify identical masks on every positive death transition.
      for (columns in split(seq_along(q), colSums(allowed[terminal,,drop=FALSE]))) {
        gated <- .Call("life_markov_alive_conditional",
          c(native, list(allow=allowed[,columns[1L]])), TRUE, PACKAGE="lifeconformal")
        success <- gated$h
        if (any(!is.finite(success)) || any(success < -eps | success > 1 + eps))
          return(fallback("terminal probability error"))
        if (any(!live)) {
          score <- apply(abs((args$accrued[!live,,drop=FALSE] -
            args$center[!live,,drop=FALSE]) / args$scale[!live,,drop=FALSE]), 1L, max)
          point[!live,columns] <- low[!live,columns] <- high[!live,columns] <-
            outer(score, q[columns], "<=")
        }
        if (!any(success[live] > 0)) next
        home$kernels <- gated$kernels
        home$kernel_factors <- NULL
        home$q <- q[columns]
        dist <- do.call(markov_conditional_coverage_grid, home)
        low[live,columns] <- dist$lower[live,,drop=FALSE] * success[live]
        high[live,columns] <- dist$upper[live,,drop=FALSE] * success[live]
        point[live,columns] <- dist$probability[live,,drop=FALSE] * success[live]
      }
    }
    if (any(!is.finite(c(low,point,high))) || any(low > high + eps) ||
        any(low < -1-eps | high < -eps | high > 1+eps)) return(fallback("unordered enclosure"))
    low <- pmax(low - eps, 0)
    high <- pmin(high + eps, 1)
    if (any(low > high)) return(fallback("roundoff envelope"))
    point <- pmax(low,pmin(high,point))
    known <- args$query_known
    fixed <- !is.na(known)
    point[fixed] <- low[fixed] <- high[fixed] <- known[fixed]
    width <- max(high-low)
    if (isTRUE(args$strict) && width > args$tolerance) return(fallback("probability tolerance"))
    diagnostics <- list(method=if (screen) "alive_time_window" else "conditioned_on_alive_time",
      maximum_probability_error=width, query_count=h, jump_count=length(times))
    if (is.null(args$query_weights))
      return(list(probability=point,lower=low,upper=high,reduced=NULL,diagnostics=diagnostics))
    reduction <- .Call("life_markov_marginal_reduce", list(probability=point,lower=low,
      upper=high,weights=as.double(args$query_weights),
      order=order(args$query_times,decreasing=TRUE,method="radix")), PACKAGE="lifeconformal")
    list(probability=NULL,lower=NULL,upper=NULL,reduced=reduction,diagnostics=diagnostics)
  }, lifeconformal_grid_budget=function(e) NULL, lifeconformal_grid_precision=function(e) NULL)
  if (is.null(value)) fallback("scalar numerical budget") else value
}

#' Certify a first crossing using the fitted terminal-time law
#'
#' Every finite candidate is initially enclosed. Only candidates whose signed
#' sums can still change the first crossing receive a finer conditional pass.
#' A failed certificate returns control to the established general solver.
#' @noRd
alive_threshold_curve <- function(fit, history, x, censor_x, moments, cdt, q,
                                  controls, queries, rank) {
  n <- nrow(history$subjects)
  original_size <- controls$grid_size
  requested_tolerance <- controls$probability_tolerance
  controls$calibration_accuracy <- "probability"
  controls$probability_tolerance <- 1
  evaluate <- function(indices, size) {
    part <- cdt
    part$q <- q[indices]
    for (key in c("time","resolved","membership"))
      part[[key]] <- cdt[[key]][, indices, drop = FALSE]
    controls$grid_size <- size
    evaluate_augmented_curve(fit, history, x, censor_x, moments, part,
      q[indices], controls, queries=queries, rank=rank)
  }
  if (is.null(queries)) queries <- prepare_calibration_queries(fit, history, censor_x, cdt)
  queries$prepared$markov_jumps <- calibration_markov_jumps(fit,x,queries$prepared)
  initial <- max(128L, original_size)
  out <- evaluate(seq_along(q), initial)
  d <- out$diagnostics
  low <- d$curve_lower_sum
  high <- d$curve_upper_sum
  stages <- list()
  counts <- d$alive_profile_counts
  for (size in unique(pmax(initial, c(1024L,4096L,8192L,16384L)))) {
    # Resolve the first possible crossing, with up to seven nearby candidates
    # sharing the same terminal gates and subject-level occupation distributions.
    repeat {
      cert <- calibration_threshold_certificate(low,high,rank)
      if (cert$threshold_stable) break
      unknown <- which(seq_along(q) < cert$certain_index & low < rank & high >= rank)
      idx <- utils::head(unknown, 8L)
      z <- evaluate(idx, size)
      d2 <- z$diagnostics
      lo <- pmax(low[idx],d2$curve_lower_sum)
      hi <- pmin(high[idx],d2$curve_upper_sum)
      if (any(lo > hi)) return(NULL)
      low[idx] <- lo
      high[idx] <- hi
      out$contributions[,idx] <- z$contributions
      counts <- counts + d2$alive_profile_counts
      stages[[length(stages)+1L]] <- list(grid_size=size,candidates=idx)
      if (high[idx[1L]] >= rank && low[idx[1L]] < rank) break
    }
    if (cert$threshold_stable) break
  }
  point <- colSums(out$contributions)
  # Intersected passes may have different representatives. Widen to include
  # the reported contribution sum, then certify its first crossing again.
  low <- pmin(low,point)
  high <- pmax(high,point)
  cert <- calibration_threshold_certificate(low,high,rank)
  crossing <- which(point >= rank)
  at <- if(length(crossing)) crossing[1L] else Inf
  # Points and brackets travel together. A disagreement or unresolved numerical
  # tie delegates to the original public algorithm with the same fixed inputs.
  if (!cert$threshold_stable || at != cert$possible_index) return(NULL)
  # At floating-point proximity to the actual rank defer to the original
  # arithmetic even when the conditioned numerical brackets appear decisive.
  rank_guard <- 128 * .Machine$double.eps *
    (1 + rank + max(colSums(abs(out$contributions))))
  if (any(rank - high[seq_along(q) < at] <= rank_guard) ||
      (is.finite(at) && low[at] - rank <= rank_guard)) return(NULL)
  d[names(cert)] <- cert
  d$curve_lower_sum <- low
  d$curve_upper_sum <- high
  d$curve_lower <- low/n
  d$curve_upper <- high/n
  d$max_curve_bound_width <- max((high-low)/n)
  d$accuracy_mode <- "threshold"
  d$initial_grid_size <- initial
  d$requested_probability_tolerance <- requested_tolerance
  d$bound_scope <- "numerical discretization with terminal-conditioning roundoff padding"
  d$occupation_engine <- "condition_on_alive"
  d$alive_profile_counts <- counts
  d$alive_refinements <- stages
  out$diagnostics <- d
  out
}
