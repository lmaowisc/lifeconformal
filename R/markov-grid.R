#' Evaluate deterministic conditional Markov coverage on a sparse reward grid
#'
#' All history queries for one baseline profile share a backward distribution
#' pass. Between fitted jumps, reward nodes translate exactly; at jumps their
#' joint masses mix with the supplied transition kernel. Only when a mixture
#' exceeds `exact_limit` are nearby nodes coalesced on a joint reward grid.
#' Every coalesced node retains lower and upper coordinate envelopes containing
#' all of its underlying rewards. These envelopes give lower and upper bounds
#' on the working-model coverage probability, including at distribution atoms.
#'
#' @param jump_times,kernels,rewards,tau See [markov_reward_moments()].
#' @param query_times Numeric conditioning times in `[0, tau]`.
#' @param query_state One-based current state for each query.
#' @param accrued,center,scale Finite query-by-target matrices. `scale` is
#'   strictly positive. A target vector is accepted when there is one query.
#' @param q Finite nonnegative increasing score thresholds.
#' @param side A common `"right"`/`"left"` side, or one per query. Through-time
#'   histories use `"right"`, excluding the kernel at the conditioning time.
#' @param grid_size Number of cells per nonconstant target coordinate, either
#'   one common value or one per target. Refinement doubles these counts.
#' @param max_nodes Maximum nodes across all states at a completed time layer.
#'   A larger intermediate mixture is also checked against this resource bound.
#' @param exact_limit Maximum mixture size retained without grid coalescing.
#' @param point_mass_threshold Minimum mass for retaining an additional exact
#'   point separately at each coalescence, from 0.001 through 2. Values above
#'   one disable only this optional retention; intrinsic no-change atoms are
#'   always preserved. This changes the partition, with all probability and
#'   support envelopes retained, rather than discarding small probabilities.
#' @param tolerance Maximum permitted coverage-probability bracket width.
#' @param max_refinements Number of additional passes at doubled resolution;
#'   the default permits four further doublings while retaining `max_nodes`.
#'   At most 20 are allowed; explicit node budgets remain active at each pass.
#' @param strict If TRUE, failure to meet `tolerance` raises a typed numerical
#'   precision error. FALSE returns the bracket and failure flag explicitly.
#' @param query_transform Optional internal function `(index, value)` returning
#'   updated probability/lower/upper vectors for one query. Calibration uses
#'   this hook for valid structural membership overrides before assessing
#'   numerical precision. Modeled extreme probabilities are not overrides.
#' @param reducer Optional internal list of functions `initialize()`,
#'   `add(state, index, probability, lower, upper)`, and `finalize(state)`.
#'   It aggregates queries without allocating query-by-threshold result matrices.
#'   A fresh state is initialized on every refinement pass; callbacks must not
#'   modify external accumulation state shared across passes.
#' @param query_known Optional query-by-threshold matrix of NA for unresolved
#'   membership or zero/one for valid structural certificates. The native pass
#'   skips distribution-query work for known entries.
#' @param query_weights Optional signed weight per query for native streaming.
#'   Cannot be combined with `reducer`. `reduced` contains signed probability,
#'   lower-bound and upper-bound sums; negative weights reverse bracket sides.
#'
#' @returns A list with query-by-threshold `probability`, `lower`, `upper`, and
#'   `diagnostics`. `probability` is the CDF of weighted representative nodes,
#'   projected into the probability brackets. It is approximate when nodes
#'   have merged. Exact within-node component variances tighten envelope bounds
#'   using one-sided variance bounds. The bracket width bounds the absolute
#'   error against the fitted kernel distribution, apart from floating
#'   arithmetic. This is a numerical bound, not a statistical confidence bound.
#'   With a reducer, the three full matrices are NULL and `reduced` contains
#'   its finalized output; the same error diagnostics are still returned.
#' @details
#' No nonzero probability is intentionally discarded, no target independence
#' is assumed, and no simulation fallback is used. The no-future-transition
#' atom in each reward-equivalent class of states is preserved separately.
#' Transitions between states with identical target rewards do not blur this
#' atom. Other atoms remain exact until coalescing; their probability remains
#' in the envelopes thereafter. A model
#' probability of zero or one is never a structural coverage certificate.
#'
#' Grids use state-relative coordinates `future_reward + time * reward(state)`.
#' A no-transition sojourn therefore leaves grid coordinates unchanged, avoiding
#' repeated interpolation of diagonal stay mass at every fitted jump.
#' Grid resolution is measured in total reward units over the full horizon.
#' Recursive coalescing can widen envelopes beyond one cell; the reported
#' probability brackets account for that accumulation. Refinement is controlled
#' by these brackets rather than just agreement of adjacent approximations.
#' Time is proportional to propagated nodes and queried nodes, with sorting
#' for coalescing and weighted CDFs; it is not constant in subjects or queries.
#' Rolling state distributions avoid storing every time layer.
#' @noRd
markov_conditional_coverage_grid <- function(
    jump_times, kernels, rewards, tau, query_times, query_state,
    accrued, center, scale, q, side = "right", grid_size = 64L,
    max_nodes = 100000L, exact_limit = 256L, tolerance = 0.01,
    max_refinements = 4L, strict = TRUE, query_transform = NULL, reducer = NULL,
    query_known = NULL, query_weights = NULL, kernel_factors = NULL,
    point_mass_threshold = 0.001) {
  if (is.null(kernel_factors)) {
    spec <- prepare_markov_recursion(jump_times, kernels, rewards, tau,
                                      0, "left", NULL)
  } else {
    # Validate common state/reward/horizon fields without constructing J R
    # matrices. Native validation also checks every factored edge and jump.
    prepare_markov_recursion(numeric(), list(), rewards, tau, 0, "left", NULL)
    if (!is.null(kernels))
      stop("Supply kernels or kernel_factors, not both.", call. = FALSE)
    spec <- list(jump_times = jump_times, kernels = NULL, kernel_row_error = 0)
  }
  h <- length(query_times)
  d <- ncol(rewards)
  k <- nrow(rewards)
  numeric_vector <- function(x) is.numeric(x) && !is.complex(x) &&
    is.null(dim(x)) && all(is.finite(x))
  integer_scalar <- function(x, lower) numeric_vector(x) && length(x) == 1L &&
    x >= lower && x == floor(x) && x <= .Machine$integer.max
  if (!numeric_vector(query_times) || !h ||
      any(query_times < 0 | query_times > tau)) {
    stop("`query_times` must be nonempty and lie in [0, tau].", call. = FALSE)
  }
  if (!numeric_vector(query_state) || length(query_state) != h ||
      any(query_state != floor(query_state) | query_state < 1 | query_state > k)) {
    stop("Supply one valid `query_state` per query.", call. = FALSE)
  }
  if (length(side) == 1L) side <- rep(side, h)
  if (!is.character(side) || length(side) != h || anyNA(side) ||
      any(!side %in% c("left", "right"))) {
    stop("`side` must be 'left' or 'right', one common value or one per query.",
         call. = FALSE)
  }
  query_matrix <- function(x, name, positive = FALSE) {
    if (h == 1L && numeric_vector(x) && length(x) == d) x <- matrix(x, 1L)
    if (!is.matrix(x) || !is.numeric(x) || is.complex(x) ||
        !identical(dim(x), c(h, d)) || any(!is.finite(x)) ||
        (positive && any(x <= 0))) {
      stop("`", name, "` must be a finite query-by-target matrix",
           if (positive) " with positive entries." else ".", call. = FALSE)
    }
    x
  }
  accrued <- query_matrix(accrued, "accrued")
  center <- query_matrix(center, "center")
  scale <- query_matrix(scale, "scale", TRUE)
  if (!numeric_vector(q) || !length(q) || any(q < 0) || any(diff(q) <= 0)) {
    stop("`q` must be finite, nonnegative and strictly increasing.", call. = FALSE)
  }
  if (length(grid_size) == 1L) grid_size <- rep(grid_size, d)
  if (!numeric_vector(grid_size) || length(grid_size) != d ||
      any(grid_size < 1 | grid_size != floor(grid_size) |
          grid_size > .Machine$integer.max)) {
    stop("`grid_size` must contain positive integer cell counts.", call. = FALSE)
  }
  if (!integer_scalar(max_nodes, 1))
    stop("`max_nodes` must be a positive integer node budget.", call. = FALSE)
  if (max_nodes < k) stop_markov_grid_budget(max_nodes, minimum_nodes = k)
  if (!integer_scalar(exact_limit, 1))
    stop("`exact_limit` must be a positive integer.", call. = FALSE)
  if (!numeric_vector(point_mass_threshold) || length(point_mass_threshold) != 1L ||
      point_mass_threshold < .001 || point_mass_threshold > 2)
    stop("`point_mass_threshold` must be a finite number from 0.001 through 2.",
         call. = FALSE)
  if (!integer_scalar(max_refinements, 0) || max_refinements > 20)
    stop("`max_refinements` must be an integer from zero to 20.", call. = FALSE)
  if (!numeric_vector(tolerance) || length(tolerance) != 1L ||
      tolerance < 0 || tolerance > 1 || !is.logical(strict) ||
      length(strict) != 1L || is.na(strict)) {
    stop("`tolerance` must be in [0, 1] and `strict` must be TRUE or FALSE.",
         call. = FALSE)
  }
  if (!is.null(query_transform) && !is.function(query_transform)) {
    stop("`query_transform` must be NULL or a function.", call. = FALSE)
  }
  if (!is.null(reducer) && (!is.list(reducer) ||
      !all(c("initialize", "add", "finalize") %in% names(reducer)) ||
      !all(vapply(reducer[c("initialize", "add", "finalize")],
                   is.function, logical(1))))) {
    stop("`reducer` must supply initialize, add, and finalize functions.", call. = FALSE)
  }
  range <- apply(rewards, 2L, function(x) diff(range(x))) * tau
  origin <- apply(rewards, 2L, function(x) min(0, x)) * tau
  if (any(!is.finite(range)) || any(!is.finite(origin))) {
    stop("Reward grid ranges overflow; check reward and time units.", call. = FALSE)
  }
  passes <- vector("list", max_refinements + 1L)
  for (refinement in seq.int(0L, max_refinements)) {
    current_grid <- grid_size * 2^refinement
    if (any(!is.finite(current_grid)) || any(current_grid > .Machine$integer.max)) {
      stop("Refined grid counts exceed the supported integer range.", call. = FALSE)
    }
    step <- ifelse(range == 0, Inf, range / current_grid)
    result <- markov_grid_pass_native(spec$jump_times, spec$kernels, rewards, tau,
                               query_times, query_state, accrued, center,
                               scale, q, side, step, origin, max_nodes,
                               exact_limit, query_transform, reducer,
                               query_known, query_weights, kernel_factors,
                               point_mass_threshold)
    width <- result$diagnostics$maximum_probability_error
    passes[[refinement + 1L]] <- list(grid_size = current_grid,
                                     maximum_width = width,
                                     peak_nodes = result$diagnostics$peak_nodes)
    if (width <= tolerance) break
  }
  result$diagnostics$grid_size <- current_grid
  result$diagnostics$initial_grid_size <- grid_size
  result$diagnostics$max_nodes <- max_nodes
  result$diagnostics$exact_limit <- exact_limit
  result$diagnostics$maximum_probability_error <- width
  result$diagnostics$tolerance <- tolerance
  result$diagnostics$tolerance_met <- width <= tolerance
  result$diagnostics$refinements <- refinement
  result$diagnostics$passes <- passes[seq_len(refinement + 1L)]
  result$diagnostics$kernel_row_error <- if (is.null(kernel_factors))
    spec$kernel_row_error else result$diagnostics$native_kernel_row_error
  if (strict && width > tolerance) {
    message <- paste0("Deterministic Markov coverage has probability bracket width ",
                      signif(width, 4), ", exceeding tolerance ", tolerance,
                      ". Increase `grid_size` via `control = list(grid_size = ...)` ",
                      "in life_train() for working-CDF prediction, or life_calibrate() ",
                      "for calibration. The requested probability_tolerance remains ",
                      "enforced. Numerical diagnostics are attached to this error.")
    stop(structure(list(message = message, call = NULL,
                        diagnostics = result$diagnostics),
                   class = c("lifeconformal_grid_precision", "error", "condition")))
  }
  result
}

#' Run one deterministic sparse-grid pass for one Markov profile
#' @noRd
markov_grid_pass <- function(jump_times, kernels, rewards, tau, query_times,
                              query_state, accrued, center, scale, q, side,
                              step, origin, max_nodes, exact_limit,
                              query_transform = NULL, reducer = NULL,
                              point_mass_threshold = 0.001) {
  k <- nrow(rewards)
  d <- ncol(rewards)
  h <- length(query_times)
  probability <- lower <- upper <- if (is.null(reducer)) matrix(0, h, length(q)) else NULL
  reduction_state <- if (!is.null(reducer)) reducer$initialize() else NULL
  maximum_width <- 0
  record_query <- function(a, value) {
    if (!is.null(query_transform)) {
      value <- query_transform(a, value)
      valid <- is.list(value) && all(c("probability", "lower", "upper") %in% names(value))
      if (valid) valid <- all(vapply(value[c("probability", "lower", "upper")],
        function(x) is.numeric(x) && !is.complex(x) && is.null(dim(x)) &&
          length(x) == length(q) && all(is.finite(x)), logical(1)))
      if (!valid || any(value$lower < -1e-12) || any(value$upper > 1 + 1e-12) ||
          any(value$probability < value$lower - 1e-12) ||
          any(value$probability > value$upper + 1e-12)) {
        stop("`query_transform` returned invalid probability brackets.", call. = FALSE)
      }
    }
    maximum_width <<- max(maximum_width, value$upper - value$lower)
    if (is.null(reducer)) {
      probability[a, ] <<- value$probability
      lower[a, ] <<- value$lower
      upper[a, ] <<- value$upper
    } else {
      reduction_state <<- reducer$add(reduction_state, a, value$probability,
                                      value$lower, value$upper)
    }
    invisible(NULL)
  }
  distributions <- lapply(seq_len(k), function(s) {
    value <- matrix(0, 1L, d)
    list(value = value, lower = value, upper = value, variance = value,
         variance_reliable = matrix(TRUE, 1L, d),
         probability = 1, protected = TRUE, constant = matrix(TRUE, 1L, d))
  })
  times <- sort(unique(c(jump_times[jump_times >= min(query_times)], query_times)),
                decreasing = TRUE)
  jump_index <- match(times, jump_times)
  query_index <- split(seq_len(h), match(query_times, times))
  anchor <- tau
  peak_nodes <- k
  coalesced_layers <- 0L
  maximum_mass_error <- 0
  traversed_nodes <- 0
  query_value <- function(a) {
    state <- query_state[a]
    nodes <- distributions[[state]]
    shift <- (anchor - query_times[a]) * rewards[state, ]
    canonical <- (tau - query_times[a]) * rewards[state, ]
    canonical <- matrix(canonical, nrow(nodes$value), d, byrow = TRUE)
    for (field in c("value", "lower", "upper")) {
      nodes[[field]] <- sweep(nodes[[field]], 2L, shift, "+")
      nodes[[field]][nodes$constant] <- canonical[nodes$constant]
    }
    total <- markov_protected_total(tau, query_times[a], rewards[state, ], accrued[a, ])
    markov_grid_cdf(nodes, accrued[a, ], center[a, ], scale[a, ], q,
                    protected_target = total)
  }
  for (i in seq_along(times)) {
    time <- times[i]
    duration <- anchor - time
    j <- jump_index[i]
    # Query-only times translate a temporary view, never the stored law.
    # This keeps ordinary atoms independent of other requested histories.
    if (!is.na(j) && duration > 0) {
      for (s in seq_len(k)) {
        shift <- duration * rewards[s, ]
        canonical <- matrix((tau - time) * rewards[s, ],
          nrow(distributions[[s]]$value), d, byrow = TRUE)
        constant <- distributions[[s]]$constant
        for (field in c("value", "lower", "upper")) {
          distributions[[s]][[field]] <-
            sweep(distributions[[s]][[field]], 2L, shift, "+")
          distributions[[s]][[field]][constant] <- canonical[constant]
        }
        distributions[[s]]$variance[constant] <- 0
        distributions[[s]]$variance_reliable[constant] <- TRUE
      }
    }
    if (!is.na(j)) anchor <- time
    queries <- query_index[[as.character(i)]]
    right <- queries[side[queries] == "right"]
    if (length(right)) {
      for (a in right) {
        value <- query_value(a)
        record_query(a, value)
      }
    }
    if (!is.na(j)) {
      kernel <- kernels[[j]]
      next_layer <- distributions
      for (s in seq_len(k)) {
        destinations <- which(kernel[s, ] > 0)
        if (length(destinations) == 1L && destinations == s) next
        pieces <- lapply(destinations, function(destination) {
          z <- distributions[[destination]]
          p <- z$probability * kernel[s, destination]
          if (any(p == 0 & z$probability > 0)) {
            stop("A positive Markov grid mass underflowed to zero.", call. = FALSE)
          }
          z$probability <- p
          # A transition exactly at tau contributes no further occupation.
          z$protected <- z$protected &
            (time == tau || all(rewards[destination, ] == rewards[s, ]))
          z$constant <- sweep(z$constant, 2L,
            time == tau | rewards[destination, ] == rewards[s, ], "&")
          z
        })
        count <- sum(vapply(pieces, function(z) length(z$probability), integer(1)))
        if (count > max_nodes) stop_markov_grid_budget(max_nodes)
        mixed <- list(value = do.call(rbind, lapply(pieces, `[[`, "value")),
                      lower = do.call(rbind, lapply(pieces, `[[`, "lower")),
                      upper = do.call(rbind, lapply(pieces, `[[`, "upper")),
                      variance = do.call(rbind, lapply(pieces, `[[`, "variance")),
                      variance_reliable = do.call(rbind, lapply(pieces, `[[`, "variance_reliable")),
                      probability = unlist(lapply(pieces, `[[`, "probability"),
                                           use.names = FALSE),
                      protected = unlist(lapply(pieces, `[[`, "protected"),
                                         use.names = FALSE),
                      constant = do.call(rbind, lapply(pieces, `[[`, "constant")))
        traversed_nodes <- traversed_nodes + count
        if (count > exact_limit) {
          mixed <- coalesce_markov_grid(mixed, step, origin - time * rewards[s, ],
                                       point_mass_threshold)
          coalesced_layers <- coalesced_layers + 1L
        }
        next_layer[[s]] <- mixed
      }
      distributions <- next_layer
      counts <- vapply(distributions, function(z) length(z$probability), integer(1))
      peak_nodes <- max(peak_nodes, sum(counts))
      if (sum(counts) > max_nodes) stop_markov_grid_budget(max_nodes)
      error <- max(vapply(distributions, function(z)
        abs(sum(z$probability) - 1), numeric(1)))
      maximum_mass_error <- max(maximum_mass_error, error)
      if (!is.finite(error) || error > 1e-10) {
        stop("Markov grid probabilities lost unit mass.", call. = FALSE)
      }
    }
    left <- queries[side[queries] == "left"]
    if (length(left)) {
      for (a in left) {
        value <- query_value(a)
        record_query(a, value)
      }
    }
  }
  list(probability = probability, lower = lower, upper = upper,
       reduced = if (!is.null(reducer)) reducer$finalize(reduction_state) else NULL,
       diagnostics = list(method = "deterministic_joint_grid_envelopes",
                          peak_nodes = peak_nodes,
                          coalesced_state_layers = coalesced_layers,
                          traversed_nodes = traversed_nodes,
                          mass_error = maximum_mass_error,
                          maximum_probability_error = maximum_width,
                          query_count = h, jump_count = length(jump_times)))
}

#' Coalesce nearby joint reward nodes while retaining their exact envelopes
#' @noRd
coalesce_markov_grid <- function(nodes, step, origin, point_mass_threshold = 0.001) {
  n <- length(nodes$probability)
  d <- ncol(nodes$value)
  if (is.null(nodes$constant)) nodes$constant <- matrix(nodes$protected, n, d)
  # Mixed-radix double keys avoid formatting every node into a string. Use
  # them only within the exactly represented integer range; higher-dimensional
  # sparse grids fall back to delimited coordinate keys.
  cells <- floor(sweep(sweep(nodes$value, 2L, origin, "-"), 2L, step, "/"))
  cells[, is.infinite(step)] <- 0
  minima <- apply(cells, 2L, min)
  spans <- 2 * (apply(cells, 2L, max) - minima + 1)
  numeric_keys <- all(is.finite(spans)) && prod(spans) < 2^52
  if (numeric_keys) {
    strides <- c(1, cumprod(spans))[seq_len(d)]
    key <- as.vector((2 * sweep(cells, 2L, minima, "-") +
      nodes$constant * 1) %*% strides)
  } else {
    key <- do.call(paste, c(as.data.frame(cells),
      as.data.frame(nodes$constant * 1), sep = ":"))
  }
  protected <- which(nodes$protected)
  # All protected components for this origin have exactly the same reward:
  # the process remained within its reward-equivalent state class. Merge them
  # with one another, but never with a neighboring grid cell.
  key[protected] <- if (numeric_keys) max(key) + 1 else "atom"
  # A sizeable exact point node must not be smeared into a neighboring cell.
  # This numerical partition rule does not discard mass or alter envelopes.
  # Retest its current mass at every layer; unlike intrinsic no-change atoms,
  # descendants do not inherit permanent protection. A unit-mass layer has at
  # most 1 / point_mass_threshold such point nodes; the explicit node budget still
  # applies to the complete distribution and intermediate mixtures.
  point <- which(!nodes$protected & nodes$probability >= point_mass_threshold &
                   rowSums(nodes$lower != nodes$upper) == 0L)
  if (length(point)) {
    coordinates <- nodes$value[point, , drop = FALSE]
    coordinates[coordinates == 0] <- 0 # Match native equality for signed zero.
    point_key <- do.call(paste, c(lapply(as.data.frame(coordinates),
      function(x) sprintf("%a", x)), sep = ":"))
    proof_key <- do.call(paste, c(as.data.frame(nodes$constant[point, , drop = FALSE] * 1), sep = ""))
    point_key <- paste(point_key, proof_key, sep = "/")
    if (numeric_keys) key[point] <- max(key) + match(point_key, unique(point_key))
    else key[point] <- paste0("point:", point_key)
  }
  group <- match(key, unique(key))
  size <- max(group)
  if (size == n) return(nodes)
  mass <- as.vector(rowsum(nodes$probability, group, reorder = FALSE))
  value <- rowsum(nodes$value * nodes$probability, group, reorder = FALSE) / mass
  lower <- upper <- matrix(0, size, d)
  for (j in seq_len(d)) {
    ord <- order(group, nodes$lower[, j])
    keep <- !duplicated(group[ord])
    lower[group[ord[keep]], j] <- nodes$lower[ord[keep], j]
    ord <- order(group, -nodes$upper[, j])
    keep <- !duplicated(group[ord])
    upper[group[ord[keep]], j] <- nodes$upper[ord[keep], j]
  }
  # Roundoff in a weighted mean must not put its representative outside its
  # envelope. This adjustment does not narrow the mathematical envelope.
  value <- pmax(lower, pmin(upper, value))
  delta <- nodes$value - value[group, , drop = FALSE]
  squared <- delta^2
  summand <- nodes$variance + squared
  weighted <- summand * nodes$probability
  # A zero obtained by squaring a nonzero tiny displacement is not a known
  # point mass. Propagate an explicit reliability flag, including subnormal
  # moments whose relative rounding error may invalidate a Cantelli bound.
  unreliable <- !nodes$variance_reliable |
    (delta != 0 & squared < .Machine$double.xmin) |
    (summand > 0 & weighted < .Machine$double.xmin)
  variance <- rowsum(weighted, group, reorder = FALSE) / mass
  reliable <- rowsum(unreliable * 1, group, reorder = FALSE) == 0 &
    (variance == 0 | variance >= .Machine$double.xmin)
  is_protected <- rep(FALSE, size)
  is_protected[group[protected]] <- TRUE
  constant <- rowsum((!nodes$constant) * 1, group, reorder = FALSE) == 0
  variance[constant] <- 0
  reliable[constant] <- TRUE
  list(value = value, lower = lower, upper = upper, variance = variance,
       variance_reliable = reliable,
       probability = mass, protected = is_protected, constant = constant)
}

#' Query all thresholds from representative and envelope score distributions
#' @noRd
markov_grid_cdf <- function(nodes, accrued, center, scale, q, protected_target = NULL) {
  standardized <- function(x) {
    total <- sweep(x, 2L, accrued, "+")
    if (!is.null(protected_target)) {
      constant <- nodes$constant
      if (is.null(constant)) constant <- matrix(nodes$protected, nrow(x), ncol(x))
      canonical <- matrix(protected_target, nrow(x), ncol(x), byrow = TRUE)
      total[constant] <- canonical[constant]
    }
    sweep(sweep(total, 2L, center, "-"), 2L, scale, "/")
  }
  value <- standardized(nodes$value)
  low <- standardized(nodes$lower)
  high <- standardized(nodes$upper)
  if (any(!is.finite(value)) || any(!is.finite(low)) || any(!is.finite(high))) {
    stop("Standardized reward scores overflow; check centers and scales.", call. = FALSE)
  }
  row_max <- function(x) if (ncol(x) == 1L) as.vector(x) else
    do.call(pmax, unname(as.data.frame(x)))
  score <- row_max(abs(value))
  minimum <- row_max(pmax(low, -high, 0))
  maximum <- row_max(pmax(abs(low), abs(high)))
  probability <- weighted_markov_cdf(score, nodes$probability, q)
  lower <- weighted_markov_cdf(maximum, nodes$probability, q)
  upper <- weighted_markov_cdf(minimum, nodes$probability, q)
  # Each node's representative and component variances are mixture moments.
  # Cantelli bounds applied to
  # each face, together with a union bound, tighten broad support envelopes.
  # Only ambiguous q/node pairs are evaluated; exact atom queries take the
  # weighted-CDF route above without these additional operations.
  ambiguous_q <- which(upper - lower > 1e-14)
  if (length(ambiguous_q)) {
    if (any(!is.finite(nodes$variance)) || any(nodes$variance < 0)) {
      stop("Within-node reward variances overflow; check reward units.", call. = FALSE)
    }
    standard_deviation <- sweep(sqrt(nodes$variance), 2L, scale, "/")
    unreliable_scale <- nodes$variance > 0 &
      standard_deviation < .Machine$double.xmin
    standard_deviation[!nodes$variance_reliable | unreliable_scale] <- Inf
    for (a in ambiguous_q) {
      idx <- which(minimum <= q[a] & maximum > q[a])
      means <- value[idx, , drop = FALSE]
      deviations <- standard_deviation[idx, , drop = FALSE]
      distances_lower <- means + q[a]
      distances_upper <- q[a] - means
      # At distance zero, a positive-variance component gives no useful bound;
      # zero variance is a known point and observes the inclusive boundary.
      failure <- matrix(0, nrow(means), ncol(means))
      inside_upper <- rep(1, nrow(means))
      for (distances in list(distances_lower, distances_upper)) {
        positive <- distances >= 0
        ratio <- matrix(0, nrow(means), ncol(means))
        nonzero <- deviations > 0
        ratio[nonzero] <- 1 / (1 +
          (distances[nonzero] / deviations[nonzero])^2)
        failure <- failure + ifelse(positive, ratio, 1)
        face_upper <- ifelse(positive, 1, ratio)
        if (ncol(face_upper) == 1L) inside_upper <- pmin(inside_upper, face_upper[, 1L]) else
          inside_upper <- pmin(inside_upper, do.call(pmin, unname(as.data.frame(face_upper))))
      }
      inside_lower <- pmax(0, 1 - rowSums(failure))
      mass <- nodes$probability[idx]
      lower[a] <- lower[a] + sum(mass * inside_lower)
      upper[a] <- upper[a] - sum(mass * (1 - inside_upper))
    }
  }
  # The representative CDF need not lie between the moment-based brackets:
  # a point at its mean can overstate mass at an atom. Project it into the
  # certified interval, preserving its explicit numerical error bound.
  probability <- pmax(lower, pmin(upper, probability))
  list(probability = probability, lower = lower, upper = upper)
}

#' Evaluate a weighted inclusive empirical CDF without a node-by-q matrix
#' @noRd
weighted_markov_cdf <- function(value, probability, q) {
  ord <- order(value)
  sums <- c(0, cumsum(probability[ord]))
  pmin(1, pmax(0, sums[findInterval(q, value[ord]) + 1L]))
}

#' Signal an explicit deterministic grid resource-budget failure
#' @noRd
stop_markov_grid_budget <- function(max_nodes, minimum_nodes = NULL) {
  requirement <- if (is.null(minimum_nodes)) "" else
    paste0(" At least ", minimum_nodes, " nodes are needed for the fitted states.")
  stop(structure(list(message = paste0("Deterministic Markov grid exceeds ",
                                       "`max_nodes = ", max_nodes,
                                       "`.", requirement,
                                       " Increase `max_nodes` via `control = list(max_nodes = ...)` ",
                                       "in life_train() for working-CDF prediction, or ",
                                       "life_calibrate() for calibration. A larger budget ",
                                       "uses more memory; the requested probability_tolerance ",
                                       "remains enforced."), call = NULL,
                        diagnostics = list(max_nodes = max_nodes,
                                           minimum_nodes = minimum_nodes)),
                 class = c("lifeconformal_grid_budget", "error", "condition")))
}
