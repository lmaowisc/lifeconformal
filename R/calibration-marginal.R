#' Identify threshold profiles eligible for two scalar marginal passes
#' @noRd
marginal_threshold_eligible <- function(fit, controls) {
  identical(fit$augmentation$type, "markov_cox") &&
    !is.complex(fit$problem$rewards) && ncol(fit$problem$rewards) == 2L &&
    all(fit$problem$rewards %in% c(0, 1)) &&
    length(controls$grid_size) == 1L && controls$grid_size == 16
}

#' Bound an initial two-target threshold query using scalar marginal laws
#'
#' Marginal bounds give the Frechet bounds L = max(0, L1 + L2 - 1) and
#' U = min(U1, U2), with a midpoint representative clipped into [L, U].
#' Joint known memberships override these values before signed reduction.
#' The threshold coordinator separately selects same-grid joint fallback;
#' this helper does not change the candidate grid or refine either marginal.
#'
#' Temporary query matrices are limited to one million cells. Unsupported
#' inputs retain the original grid validator. Signed aggregation uses the
#' existing native query order and binary64 products and corrections.
#' @noRd
markov_conditional_coverage_marginal <- function(jump_times, kernels, rewards, tau,
                                                query_times, query_state, accrued,
                                                center, scale, q, ...,
                                                max_cells = 1000000L) {
  stopifnot(is.numeric(max_cells), length(max_cells) == 1L,
    is.finite(max_cells), max_cells >= 1, max_cells <= 1000000,
    max_cells == floor(max_cells))
  original_grid <- markov_conditional_coverage_grid
  extra <- list(...)
  args <- c(list(jump_times = jump_times, kernels = kernels, rewards = rewards,
    tau = tau, query_times = query_times, query_state = query_state,
    accrued = accrued, center = center, scale = scale, q = q), extra)
  supported <- is.matrix(rewards) && is.numeric(rewards) &&
    !is.complex(rewards) && ncol(rewards) == 2L &&
    all(rewards %in% c(0, 1)) && length(query_times) > 0L &&
    identical(extra$strict, FALSE) && identical(as.double(extra$max_refinements), 0) &&
    length(extra$grid_size) == 1L && extra$grid_size == 16 &&
    identical(as.double(extra$exact_limit), 32) &&
    identical(as.double(extra$point_mass_threshold), .1) &&
    identical(extra$side, "right") && is.null(extra$query_transform) &&
    is.null(extra$reducer)
  h <- length(query_times)
  length_q <- length(q)
  query_shape <- all(vapply(list(accrued, center, scale), function(x)
    is.matrix(x) && identical(dim(x), c(h, 2L)), logical(1)))
  known <- args$query_known
  weights <- args$query_weights
  known_valid <- is.null(known) || (is.matrix(known) && is.numeric(known) &&
    !is.complex(known) && identical(dim(known), c(h, length_q)) &&
    all(is.na(known) | known %in% c(0, 1)))
  weights_valid <- is.null(weights) || (is.numeric(weights) && !is.complex(weights) &&
    is.null(dim(weights)) && length(weights) == h && all(is.finite(weights)))
  # Valid one-query vectors and malformed full-query controls return to the
  # original validator before component slicing can reinterpret their shape.
  if (!supported || !query_shape || !known_valid || !weights_valid ||
      length_q == 0L || h > floor(max_cells / length_q))
    return(do.call(original_grid, args))
  marginal_args <- args
  marginal_args$query_known <- NULL
  marginal_args$query_weights <- NULL
  low <- matrix(0, h, length_q)
  high <- matrix(1, h, length_q)
  details <- vector("list", 2L)
  for (j in 1:2) {
    local <- marginal_args
    for (field in c("rewards", "accrued", "center", "scale"))
      local[[field]] <- marginal_args[[field]][, j, drop = FALSE]
    value <- tryCatch(do.call(markov_marginal_scalar_coverage, local),
      lifeconformal_grid_budget = function(e) e)
    if (inherits(value, "lifeconformal_grid_budget")) {
      ordinary <- do.call(original_grid, args)
      ordinary$diagnostics$marginal_scalar_attempts <- j
      ordinary$diagnostics$marginal_budget_fallback <- list(component = j, error = value)
      return(ordinary)
    }
    low <- low + value$lower
    high <- pmin(high, value$upper)
    details[[j]] <- value$diagnostics
  }
  low <- pmax(low - 1, 0)
  # Preserve the explicit widening for a subtraction-rounding reversal.
  low <- pmin(low, high)
  point <- pmax(low, pmin(high, low + (high - low) / 2))
  if (!is.null(known)) {
    certified <- !is.na(known)
    point[certified] <- low[certified] <- high[certified] <- known[certified]
  }
  width <- max(high - low)
  tolerance <- if ("tolerance" %in% names(extra)) extra$tolerance else
    eval(formals(original_grid)$tolerance, envir = environment(original_grid))
  diagnostics <- list(method = "deterministic_marginal_frechet_grid",
    maximum_probability_error = width, tolerance_met = width <= tolerance,
    query_count = h, jump_count = length(jump_times),
    peak_nodes = max(vapply(details, `[[`, numeric(1), "peak_nodes")),
    traversed_nodes = sum(vapply(details, `[[`, numeric(1), "traversed_nodes")),
    grid_size = extra$grid_size, point_mass_threshold = extra$point_mass_threshold,
    exact_limit = extra$exact_limit, marginal_diagnostics = details,
    marginal_applied = TRUE, marginal_scalar_attempts = 2L)
  if (is.null(weights))
    return(list(probability = point, lower = low, upper = high,
      reduced = NULL, diagnostics = diagnostics))
  # Descending-time, stable-input-index order and native binary64 reduction
  # match the existing streamed grid; R colSums can change a signed crossing.
  order <- order(query_times, decreasing = TRUE, method = "radix")
  reduced <- .Call("life_markov_marginal_reduce", list(probability = point,
    lower = low, upper = high, weights = as.double(weights),
    order = as.integer(order)), PACKAGE = "lifeconformal")
  list(probability = NULL, lower = NULL, upper = NULL,
    reduced = reduced, diagnostics = diagnostics)
}
