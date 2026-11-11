#' Evaluate a scalar marginal using a coherent first-exit clock when eligible
#'
#' Only the initial two-target marginal threshold path calls this dispatcher.
#' Binary rewards must never increase on a fitted edge, all queries must use
#' through-time conditioning, active accrued rewards must equal query times
#' exactly, and the profile must share one center and scale. Unsupported inputs
#' retain the ordinary grid and its validation. The fitted kernel constructor
#' is shared with the native general grid; no R kernel matrix per jump is made.
#'
#' The query output, each rolling state buffer and kernel ledger each have a
#' one-million-cell cap. The exact-clock result has equal numerical brackets,
#' subject to the same floating-point scope as the existing monotone engine.
#' A clock has no occupation grid: its scalar grid diagnostics are NA, while
#' the enclosing marginal profile retains its requested partition context.
#' @inheritParams markov_conditional_coverage_grid
#' @param ... Internal scalar-grid controls, checked before dispatch.
#' @returns The scalar probability-grid contract, with method-specific diagnostics.
#' @noRd
markov_marginal_scalar_coverage <- function(jump_times,kernels,rewards,tau,query_times,query_state,
    accrued,center,scale,q,...) {
  original <- markov_conditional_coverage_grid
  max_cells <- max_kernel_cells <- 1000000L
  numeric_vector <- function(x) is.numeric(x) && !is.complex(x) &&
    is.null(dim(x)) && all(is.finite(x))
  extra <- list(...)
  args <- c(list(jump_times=jump_times,kernels=kernels,rewards=rewards,tau=tau,
    query_times=query_times,query_state=query_state,accrued=accrued,
    center=center,scale=scale,q=q),extra)
  h <- length(query_times); l <- length(q); jnum <- length(jump_times)
  f <- extra$kernel_factors
  numeric_matrix <- function(x,nr,nc) is.matrix(x) && is.numeric(x) &&
    !is.complex(x) && identical(dim(x),c(as.integer(nr),as.integer(nc))) && all(is.finite(x))
  allowed <- c('side','grid_size','max_nodes','exact_limit','point_mass_threshold',
    'tolerance','max_refinements','strict','query_transform','reducer',
    'query_known','query_weights','kernel_factors')
  supported <- is.matrix(rewards) && is.numeric(rewards) && !is.complex(rewards) &&
    ncol(rewards)==1L && nrow(rewards)>0L && all(rewards %in% c(0,1)) &&
    (is.null(rownames(rewards)) || (!anyNA(rownames(rewards)) &&
      all(nzchar(rownames(rewards))) && !anyDuplicated(rownames(rewards)))) &&
    h>0L && l>0L && length(names(extra))==length(extra) &&
    !anyDuplicated(names(extra)) && all(names(extra) %in% allowed) &&
    identical(extra$strict,FALSE) && identical(as.double(extra$max_refinements),0) &&
    identical(as.double(extra$grid_size),16) && identical(as.double(extra$exact_limit),32) &&
    identical(as.double(extra$point_mass_threshold),.1) && identical(extra$side,'right') &&
    is.null(extra$query_transform) && is.null(extra$reducer) &&
    is.null(extra$query_known) && is.null(extra$query_weights) &&
    numeric_vector(tau) && length(tau)==1L && tau>0 &&
    numeric_vector(query_times) && all(query_times>=0 & query_times<=tau) &&
    numeric_vector(jump_times) && all(jump_times>0 & jump_times<=tau) && all(diff(jump_times)>0) &&
    numeric_vector(q) && all(q>=0) && all(diff(q)>0) &&
    numeric_vector(extra$tolerance) && length(extra$tolerance)==1L &&
    extra$tolerance>=0 && extra$tolerance<=1 &&
    numeric_vector(extra$max_nodes) && length(extra$max_nodes)==1L &&
    extra$max_nodes>=nrow(rewards) && extra$max_nodes==floor(extra$max_nodes) &&
    extra$max_nodes<=.Machine$integer.max &&
    numeric_matrix(accrued,h,1L) && numeric_matrix(center,h,1L) && numeric_matrix(scale,h,1L) &&
    all(scale>0) && all(center==center[1L]) && all(scale==scale[1L]) &&
    numeric_vector(query_state) && length(query_state)==h &&
    all(query_state>=1 & query_state<=nrow(rewards) & query_state==floor(query_state)) &&
    all(accrued[rewards[query_state,1L]==1,1L]==query_times[rewards[query_state,1L]==1]) &&
    is.null(kernels) && is.list(f) &&
    all(c('base_increments','relative_risk','from','to','absorbing') %in% names(f))
  if(!supported) return(do.call(original,args))
  # An overflowing support endpoint can still leave a particular reachable
  # target finite. Let the ordinary grid decide that case rather than treating
  # Inf as a clock cutoff or rejecting every profile with that support bound.
  endpoints <- c(center[1L]/scale[1L],(tau-center[1L])/scale[1L])
  inactive <- rewards[query_state,1L]==0
  if(any(!is.finite(endpoints)) ||
     any(!is.finite((accrued[inactive,1L]-center[1L])/scale[1L])))
    return(do.call(original,args))
  k <- nrow(rewards); e <- length(f$relative_risk)
  factors_valid <- e<=k*(k-1) && length(f$absorbing)<=k &&
    all(vapply(f[c('relative_risk','from','to','absorbing')],numeric_vector,logical(1))) &&
    length(f$from)==e && length(f$to)==e && numeric_matrix(f$base_increments,jnum,e) &&
    all(f$relative_risk>0) && all(f$base_increments>=0) &&
    all(f$from>=1 & f$from<=k & f$from==floor(f$from)) &&
    all(f$to>=1 & f$to<=k & f$to==floor(f$to)) && all(f$from!=f$to) &&
    all(f$absorbing>=1 & f$absorbing<=k & f$absorbing==floor(f$absorbing)) &&
    !anyDuplicated(f$absorbing) && !anyDuplicated(paste(f$from,f$to)) &&
    all(rewards[f$to,1L]<=rewards[f$from,1L])
  # The ordinary validator retains responsibility for unsupported/malformed
  # factors. All memory guards precede the new native allocation.
  if(!factors_valid || h>floor(max_cells/l) || k>floor(max_cells/l) ||
     k>floor(max_kernel_cells/k) || (jnum>0 && k*k>floor(max_kernel_cells/jnum)))
    return(do.call(original,args))
  f <- lapply(f,function(x){if(is.numeric(x)) storage.mode(x)<-'double';x})
  value <- .Call("life_markov_monotone_queries",list(jump_times=as.double(jump_times),
    kernel_factors=f,rewards=as.double(rewards),tau=as.double(tau),
    query_times=as.double(query_times),query_state=as.double(query_state),
    accrued=as.double(accrued),center=as.double(center[1L]),scale=as.double(scale[1L]),
    q=as.double(q),max_cells=as.double(max_cells),max_kernel_cells=as.double(max_kernel_cells),
    inspect_kernels=FALSE), PACKAGE="lifeconformal")
  p <- value$probability
  diagnostics <- list(method='deterministic_exact_monotone_clocks',engine='native_scalar_monotone_queries',
    maximum_probability_error=0,tolerance_met=TRUE,query_count=h,jump_count=jnum,
    peak_nodes=0,traversed_nodes=0,grid_size=NA_integer_,exact_limit=NA_integer_,
    point_mass_threshold=NA_real_,kernel_row_error=value$kernel_row_error,
    adjusted_kernel_rows=value$adjusted_rows,clock_state_updates=value$state_updates,
    monotone_query_applied=TRUE)
  list(probability=p,lower=p,upper=p,reduced=NULL,diagnostics=diagnostics)
}
