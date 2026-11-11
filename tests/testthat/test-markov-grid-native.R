native_grid_inputs <- function(jumps = 5L, dimensions = 2L) {
  times <- seq_len(jumps) / (jumps + 1)
  kernel <- matrix(c(.65, .25, .10, .20, .65, .15, 0, 0, 1), 3, byrow = TRUE)
  rewards <- rbind(c(1, 1), c(0, 1), c(0, 0))
  if (dimensions == 1L) rewards <- rewards[, 1, drop = FALSE]
  query_times <- c(0, times[1], times[1], .71, 1)
  h <- length(query_times)
  list(jump_times = times, kernels = rep(list(kernel), jumps), rewards = rewards,
       tau = 1, query_times = query_times, query_state = c(1, 1, 2, 2, 3),
       accrued = matrix(rep(query_times / 3, dimensions), h, dimensions),
       center = matrix(rep(c(.4, .7)[seq_len(dimensions)], each = h), h, dimensions),
       scale = matrix(.3, h, dimensions), q = seq(0, 4, .0625),
       side = c("right", "left", "right", "right", "left"),
       step = rep(.15, dimensions), origin = rep(0, dimensions),
       max_nodes = 100000L, exact_limit = 2L)
}

test_that("native recurrent passes preserve the R envelopes and moments", {
  for (dimension in 1:2) {
    args <- native_grid_inputs(dimensions = dimension)
    native <- do.call(markov_grid_pass_native, args)
    reference <- do.call(markov_grid_pass, args)
    for (name in c("probability", "lower", "upper")) {
      expect_equal(native[[name]], reference[[name]], tolerance = 1e-12)
    }
    expect_equal(native$diagnostics$peak_nodes, reference$diagnostics$peak_nodes)
    expect_equal(native$diagnostics$coalesced_state_layers,
                 reference$diagnostics$coalesced_state_layers)
    expect_equal(native$diagnostics$traversed_nodes, reference$diagnostics$traversed_nodes)
    expect_equal(native$diagnostics$maximum_probability_error,
                 reference$diagnostics$maximum_probability_error, tolerance = 1e-12)
    expect_identical(native$diagnostics$engine, "native_cpp")
  }
})

test_that("native probability brackets cover independently enumerated paths", {
  set.seed(4107)
  for (scenario in seq_len(6)) {
    args <- native_grid_inputs(jumps = 4L)
    args$kernels <- lapply(args$jump_times, function(t) {
      kernel <- matrix(stats::runif(9, .05, 1), 3)
      kernel <- kernel / rowSums(kernel)
      kernel[3, ] <- c(0, 0, 1)
      kernel
    })
    args$step <- rep(stats::runif(1, .15, .35), 2)
    native <- do.call(markov_grid_pass_native, args)
    for (i in seq_along(args$query_times)) {
      atoms <- propagate_markov_rewards(args$jump_times, args$kernels, args$rewards,
        args$tau, start = args$query_times[i], side = args$side[i],
        accrued = args$accrued[i, ], state = args$query_state[i])$distributions[[1L]]
      standardized <- abs(sweep(sweep(atoms$values, 2, args$center[i, ], "-"),
                                 2, args$scale[i, ], "/"))
      score <- apply(standardized, 1, max)
      exact <- vapply(args$q, function(q) sum(atoms$probability[score <= q]), numeric(1))
      expect_true(all(native$lower[i, ] <= exact + 1e-12))
      expect_true(all(native$upper[i, ] >= exact - 1e-12))
    }
    expect_true(all(native$probability >= native$lower - 1e-12))
    expect_true(all(native$probability <= native$upper + 1e-12))
  }
})

test_that("native exact atoms retain boundaries and reward-equivalent protection", {
  args <- native_grid_inputs(jumps = 2L)
  args$rewards <- matrix(c(1, 1, 0), 3, 1)
  args$accrued <- matrix(0, 5, 1)
  args$center <- matrix(1, 5, 1)
  args$scale <- matrix(1, 5, 1)
  args$step <- 1
  args$origin <- 0
  args$q <- c(0, 1)
  got <- do.call(markov_grid_pass_native, args)
  # Within-class changes preserve the exact alive-through-tau atom.
  expected <- .65 * .9 + .25 * .85
  expect_equal(got$lower[1, 1], expected)
  expect_equal(got$upper[1, 1], expected)
  args$exact_limit <- 100000L
  args$q <- seq(0, 2, 1/6)
  native <- do.call(markov_grid_pass_native, args)
  reference <- do.call(markov_grid_pass, args)
  expect_equal(native$probability, reference$probability, tolerance = 1e-13)
  expect_identical(native$lower, native$upper)
})

test_that("native grouping supports constant and high-dimensional coordinates", {
  args <- native_grid_inputs(jumps = 3L)
  args$rewards <- cbind(args$rewards, constant = 1, -args$rewards)
  args$accrued <- cbind(args$accrued, args$query_times, -args$accrued)
  args$center <- cbind(args$center, 1, -args$center)
  args$scale <- matrix(.3, 5, 5)
  # Product of occupied coordinate spans exceeds the safe mixed-radix range.
  args$step <- c(1e-7, 1e-7, Inf, 1e-7, 1e-7)
  args$origin <- c(0, 0, 0, -1, -1)
  native <- do.call(markov_grid_pass_native, args)
  reference <- do.call(markov_grid_pass, args)
  expect_equal(native$probability, reference$probability, tolerance = 1e-12)
  expect_equal(native$lower, reference$lower, tolerance = 1e-12)
  expect_equal(native$upper, reference$upper, tolerance = 1e-12)
  expect_equal(native$diagnostics$peak_nodes, reference$diagnostics$peak_nodes)
})

test_that("native empty jump grids preserve one-row and one-threshold shapes", {
  args <- native_grid_inputs(dimensions = 1L)
  args$jump_times <- numeric()
  args$kernels <- list()
  args$query_times <- .5
  args$query_state <- 1L
  args$accrued <- matrix(.5, 1, 1)
  args$center <- matrix(1, 1, 1)
  args$scale <- matrix(1, 1, 1)
  args$q <- 0
  args$side <- "right"
  native <- do.call(markov_grid_pass_native, args)
  expect_identical(dim(native$probability), c(1L, 1L))
  expect_equal(native$probability[1, 1], 1)
  expect_equal(native$diagnostics$maximum_probability_error, 0)
})

test_that("native streaming and structural overrides match full results", {
  args <- native_grid_inputs()
  weights <- c(1, -.3, -.2, -.1, .8)
  calls <- 0L
  reducer <- list(
    initialize = function() { calls <<- calls + 1L; matrix(0, 3, length(args$q)) },
    add = function(state, index, probability, lower, upper) {
      w <- weights[index]
      state[1, ] <- state[1, ] + w * probability
      state[2, ] <- state[2, ] + w * if (w >= 0) lower else upper
      state[3, ] <- state[3, ] + w * if (w >= 0) upper else lower
      state
    }, finalize = identity)
  override <- function(index, value) {
    if (index == 1L) value$probability[] <- value$lower[] <- value$upper[] <- 1
    value
  }
  full <- do.call(markov_grid_pass_native, c(args, list(query_transform = override)))
  streamed <- do.call(markov_grid_pass_native, c(args,
    list(query_transform = override, reducer = reducer)))
  reference <- do.call(markov_grid_pass, c(args,
    list(query_transform = override, reducer = reducer)))
  expect_equal(streamed$reduced, reference$reduced, tolerance = 1e-12)
  expect_equal(streamed$reduced[1, ], as.vector(weights %*% full$probability), tolerance = 1e-12)
  expect_null(streamed$probability)
  expect_null(streamed$lower)
  expect_null(streamed$upper)
  expect_equal(calls, 2L)
  expect_equal(streamed$diagnostics$maximum_probability_error,
               full$diagnostics$maximum_probability_error)
  expect_true(all(full$probability[1, ] == 1))
  expect_error(do.call(markov_grid_pass_native, c(args,
    list(query_transform = function(index, value) list(probability = 2)))), "invalid probability")
  expect_error(do.call(markov_grid_pass_native, c(args,
    list(query_transform = function(index, value) stop("intentional callback failure")))),
    "intentional callback failure")
  # Callback errors must release native state and leave the next call usable.
  expect_equal(do.call(markov_grid_pass_native, args)$probability,
               do.call(markov_grid_pass, args)$probability, tolerance = 1e-12)
})

test_that("native grid budgets and malformed numerical inputs fail safely", {
  args <- native_grid_inputs()
  invalid <- args; invalid$max_nodes <- 3L
  expect_error(do.call(markov_grid_pass_native, invalid), class = "lifeconformal_grid_budget")
  invalid <- args; invalid$query_state[1] <- 0
  expect_error(do.call(markov_grid_pass_native, invalid), "Invalid native")
  invalid <- args; invalid$scale[1, 1] <- 0
  expect_error(do.call(markov_grid_pass_native, invalid), "positive")
  invalid <- args; invalid$q <- c(1, 0)
  expect_error(do.call(markov_grid_pass_native, invalid), "increasing")
  invalid <- args; invalid$step[1] <- 0
  expect_error(do.call(markov_grid_pass_native, invalid), "cell widths")
  invalid <- args; invalid$kernels[[1]][1, 1] <- -1
  expect_error(do.call(markov_grid_pass_native, invalid), "nonnegative")
  invalid <- args; invalid$query_times <- as.complex(invalid$query_times)
  expect_error(do.call(markov_grid_pass_native, invalid), "numeric vector")
  invalid <- args; invalid$accrued <- matrix(0, 1, 2)
  expect_error(do.call(markov_grid_pass_native, invalid), "dimensions")
})

test_that("structured certificates and signed streaming match callback results", {
  args <- native_grid_inputs()
  known <- matrix(NA_real_, length(args$query_times), length(args$q))
  known[1, ] <- 1
  known[2, 1:15] <- 0
  known[3:5, 50:65] <- 1
  weights <- c(2, -.7, -.2, -.1, .3)
  override <- function(index, value) {
    fixed <- !is.na(known[index, ])
    for (name in c("probability", "lower", "upper")) value[[name]][fixed] <- known[index, fixed]
    value
  }
  callback <- do.call(markov_grid_pass_native, c(args, list(query_transform = override)))
  direct <- do.call(markov_grid_pass_native, c(args, list(query_known = known)))
  for (name in c("probability", "lower", "upper")) {
    expect_equal(direct[[name]], callback[[name]], tolerance = 1e-12)
  }
  weighted <- do.call(markov_grid_pass_native,
                     c(args, list(query_known = known, query_weights = weights)))
  expect_equal(weighted$reduced$probability, as.vector(weights %*% direct$probability),
               tolerance = 1e-12)
  positive <- weights >= 0
  lower <- colSums(direct$lower[positive, , drop = FALSE] * weights[positive]) +
    colSums(direct$upper[!positive, , drop = FALSE] * weights[!positive])
  upper <- colSums(direct$upper[positive, , drop = FALSE] * weights[positive]) +
    colSums(direct$lower[!positive, , drop = FALSE] * weights[!positive])
  expect_equal(weighted$reduced$lower, lower, tolerance = 1e-12)
  expect_equal(weighted$reduced$upper, upper, tolerance = 1e-12)
  expect_null(weighted$probability)
  expect_equal(weighted$diagnostics$maximum_probability_error,
               direct$diagnostics$maximum_probability_error, tolerance = 1e-12)
  known[1, 1] <- 2
  expect_error(do.call(markov_grid_pass_native, c(args, list(query_known = known))),
               "NA, zero, or one")
  expect_error(do.call(markov_grid_pass_native, c(args, list(query_weights = NA_real_))),
               "finite")
  expect_error(do.call(markov_grid_pass_native, c(args,
    list(query_known = matrix(1, nrow(known), ncol(known)),
         query_weights = rep(1e308, nrow(known))))), "aggregation overflowed")
})

test_that("R conditions and callback allocations unwind native state safely", {
  args <- native_grid_inputs(jumps = 2L)
  callback <- function(index, value) {
    gc()
    if (index == 2L) {
      stop(structure(list(message = "native callback condition", call = NULL),
                     class = c("life_grid_callback_test", "error", "condition")))
    }
    value
  }
  expect_error(do.call(markov_grid_pass_native, c(args, list(query_transform = callback))),
               class = "life_grid_callback_test")
  expect_equal(do.call(markov_grid_pass_native, args)$probability,
               do.call(markov_grid_pass, args)$probability, tolerance = 1e-12)
})
