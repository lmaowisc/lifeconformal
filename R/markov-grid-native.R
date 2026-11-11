#' Native sparse Markov grid pass with the reference numerical guarantees
#'
#' This is the compiled implementation of one `markov_grid_pass()` refinement
#' pass. It retains the R implementation as an independent reference and uses
#' the same state-relative cells, protected atoms, moment/envelope bounds,
#' query-side convention, transformation hook and streaming reducer. Numerical
#' tolerance selection and refinement remain in the calling coordinator.
#'
#' Cell membership uses compact direct lookup or an open-addressed numeric-key
#' hash, with first-occurrence order and per-call workspace reuse. Exact point
#' nodes of mass at least `point_mass_threshold` retain separate coordinate groups; this is
#' re-evaluated at each layer and does not permanently protect descendants.
#' Weighted moments use
#' two passes, avoiding cancellation of raw second moments. Inclusive CDFs use
#' threshold buckets rather than sorting every node score; all thresholds
#' reuse each node's standardized coordinates and moments. No probability is
#' intentionally discarded, and no simulated Markov paths are used.
#'
#' @inheritParams markov_grid_pass
#' @param query_known Optional query-by-threshold matrix, with NA for unresolved
#'   membership and zero/one for valid structural certificates. Known entries
#'   skip native CDF and moment-bound work, then receive exact probabilities.
#' @param query_weights Optional finite weight per query. Returns streamed
#'   signed sums in `reduced`, including lower/upper bounds with the correct
#'   weight signs. Cannot be combined with an R reducer.
#' @param kernel_factors Optional factored Cox jump specification used instead
#'   of `kernels`: shared `base_increments`, one profile's `relative_risk`,
#'   `from`, `to`, and `absorbing` state indices. The native engine constructs
#'   and validates the same normalized kernels without an R matrix per jump.
#' @returns The same list contract as markov_grid_pass, with a native-engine
#'   diagnostic. R callback reducers initialize afresh on every call.
#' @noRd
markov_grid_pass_native <- function(jump_times, kernels, rewards, tau,
                                     query_times, query_state, accrued, center,
                                     scale, q, side, step, origin, max_nodes,
                                     exact_limit, query_transform = NULL,
                                     reducer = NULL, query_known = NULL,
                                     query_weights = NULL, kernel_factors = NULL,
                                     point_mass_threshold = 0.001) {
  numeric_matrix <- function(x) {
    if (!is.matrix(x) || !is.numeric(x) || is.complex(x)) {
      stop("Native grid matrices must be numeric matrices.", call. = FALSE)
    }
    storage.mode(x) <- "double"
    x
  }
  numeric_vector <- function(x, name) {
    if (!is.numeric(x) || is.complex(x) || !is.null(dim(x))) {
      stop(name, " must be a numeric vector.", call. = FALSE)
    }
    as.double(x)
  }
  if (is.null(kernel_factors)) {
    if (!is.list(kernels)) stop("kernels must be a list.", call. = FALSE)
  } else {
    if (!is.null(kernels))
      stop("Supply kernels or kernel_factors, not both.", call. = FALSE)
    required <- c("base_increments", "relative_risk", "from", "to", "absorbing")
    if (!is.list(kernel_factors) || !all(required %in% names(kernel_factors)))
      stop("kernel_factors must supply base_increments, relative_risk, from, to and absorbing.",
           call. = FALSE)
    kernel_factors$base_increments <- numeric_matrix(kernel_factors$base_increments)
    for (name in required[-1L])
      kernel_factors[[name]] <- numeric_vector(kernel_factors[[name]], name)
  }
  if (!is.character(side)) stop("side must contain left/right values.", call. = FALSE)
  if (length(side) == 1L) side <- rep(side, length(query_times))
  if (anyNA(side) || any(!side %in% c("right", "left"))) {
    stop("side must contain left/right values.", call. = FALSE)
  }
  if (!is.null(query_weights) && !is.null(reducer)) {
    stop("Use query_weights or reducer, not both.", call. = FALSE)
  }
  if (!is.null(query_known)) query_known <- numeric_matrix(query_known)
  tryCatch(
    .Call("life_markov_grid_native", list(
      jump_times = numeric_vector(jump_times, "jump_times"),
      kernels = if (is.null(kernel_factors)) lapply(kernels, numeric_matrix) else NULL,
      kernel_factors = kernel_factors, rewards = numeric_matrix(rewards),
      tau = numeric_vector(tau, "tau"), query_times = numeric_vector(query_times, "query_times"),
      query_state = numeric_vector(query_state, "query_state"), accrued = numeric_matrix(accrued),
      center = numeric_matrix(center), scale = numeric_matrix(scale),
      q = numeric_vector(q, "q"), side = as.integer(side == "left"),
      step = numeric_vector(step, "step"), origin = numeric_vector(origin, "origin"),
      max_nodes = numeric_vector(max_nodes, "max_nodes"),
      exact_limit = numeric_vector(exact_limit, "exact_limit"),
      point_mass_threshold = numeric_vector(point_mass_threshold, "point_mass_threshold"),
      query_transform = query_transform,
      reducer = reducer, query_known = query_known,
      query_weights = if (is.null(query_weights)) NULL else
        numeric_vector(query_weights, "query_weights")), PACKAGE = "lifeconformal"),
    error = function(e) {
      if (startsWith(conditionMessage(e), "Native Markov grid node budget exceeded")) {
        stop_markov_grid_budget(max_nodes)
      }
      stop(e)
    })
}
