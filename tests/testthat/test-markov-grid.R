grid_test_spec <- function() {
  list(times = c(1, 2, 3),
       kernels = list(matrix(c(.6, .4, 0, .3, .4, .3, 0, 0, 1), 3, byrow = TRUE),
                      matrix(c(.5, .3, .2, .2, .5, .3, 0, 0, 1), 3, byrow = TRUE),
                      matrix(c(.8, .1, .1, .2, .5, .3, 0, 0, 1), 3, byrow = TRUE)),
       rewards = rbind(c(1, 1), c(0, 1), c(0, 0)), tau = 4)
}

coalescing_grid_test_spec <- function() {
  z <- grid_test_spec()
  # Keep genuinely distinct low-mass paths below the exact-point retention
  # threshold, so coarse passes must merge uncertain target values.
  rare <- 1e-4
  kernel <- matrix(c(1 - 2 * rare, rare, rare,
                     rare, 1 - 2 * rare, rare,
                     0, 0, 1), 3, byrow = TRUE)
  z$kernels <- rep(list(kernel), length(z$times))
  z
}

test_that("exact grid mode matches atom enumeration with query sides", {
  z <- grid_test_spec()
  times <- c(0, 1, 1, 1.5, 2, 3, 4)
  sides <- c("right", "left", "right", "right", "left", "right", "left")
  state <- c(1, 2, 1, 2, 1, 3, 2)
  accrued <- cbind(times / 3, times / 2)
  center <- matrix(c(2, 3), length(times), 2, byrow = TRUE)
  scale <- matrix(c(1, .7), length(times), 2, byrow = TRUE)
  q <- seq(0, 5, .125)
  got <- markov_conditional_coverage_grid(z$times, z$kernels, z$rewards,
    z$tau, times, state, accrued, center, scale, q, side = sides,
    exact_limit = 10000, tolerance = 0)
  reference <- t(vapply(seq_along(times), function(i) {
    atoms <- propagate_markov_rewards(z$times, z$kernels, z$rewards,
      z$tau, start = times[i], side = sides[i], accrued = accrued[i, ],
      state = state[i])$distributions[[1]]
    score <- apply(abs(sweep(sweep(atoms$values, 2, center[i, ], "-"),
                              2, scale[i, ], "/")), 1, max)
    vapply(q, function(x) sum(atoms$probability[score <= x]), numeric(1))
  }, numeric(length(q))))
  expect_equal(got$probability, reference, tolerance = 1e-13)
  expect_equal(got$lower, reference, tolerance = 1e-13)
  expect_equal(got$upper, reference, tolerance = 1e-13)
  expect_true(got$diagnostics$tolerance_met)
  expect_equal(got$diagnostics$maximum_probability_error, 0)
})

test_that("coalesced joint envelopes contain the exact CDF", {
  z <- coalescing_grid_test_spec()
  q <- sort(unique(c(seq(0, 5, .025), sqrt(2))))
  exact <- markov_conditional_coverage_grid(z$times, z$kernels, z$rewards,
    z$tau, 0, 1, c(0, 0), c(1.8, 2.7), c(1, .8), q,
    exact_limit = 10000, tolerance = 0)
  coarse <- markov_conditional_coverage_grid(z$times, z$kernels, z$rewards,
    z$tau, 0, 1, c(0, 0), c(1.8, 2.7), c(1, .8), q,
    grid_size = 2, exact_limit = 2, tolerance = 1, max_refinements = 0)
  expect_true(coarse$diagnostics$coalesced_state_layers > 0)
  expect_gt(coarse$diagnostics$maximum_probability_error, 0)
  expect_true(all(coarse$lower <= exact$probability + 1e-13))
  expect_true(all(coarse$upper >= exact$probability - 1e-13))
  expect_true(all(coarse$lower <= coarse$probability + 1e-13))
  expect_true(all(coarse$upper >= coarse$probability - 1e-13))
  expect_true(all(diff(coarse$probability[1, ]) >= 0))
  expect_lt(coarse$diagnostics$mass_error, 1e-13)
  fine <- markov_conditional_coverage_grid(z$times, z$kernels, z$rewards,
    z$tau, 0, 1, c(0, 0), c(1.8, 2.7), c(1, .8), q,
    grid_size = 64, exact_limit = 2, tolerance = 0, max_refinements = 0)
  expect_equal(fine$probability, exact$probability, tolerance = 1e-13)
})

test_that("deterministic and no-transition atoms remain exact", {
  rewards <- matrix(c(1, 0), 2, 1)
  kernel <- matrix(c(.6, .4, 0, 1), 2, byrow = TRUE)
  got <- markov_conditional_coverage_grid(c(1, 2), list(kernel, kernel),
    rewards, 3, c(0, 0, 2), c(1, 2, 1), matrix(0, 3, 1),
    matrix(3, 3, 1), matrix(1, 3, 1), c(0, 1, 2, 3),
    grid_size = 1, exact_limit = 1, tolerance = 1, max_refinements = 0)
  expect_equal(got$lower[1, 1], .36)
  expect_equal(got$upper[1, 1], .36)
  expect_equal(got$probability[2, ], c(0, 0, 0, 1))
  empty <- markov_conditional_coverage_grid(numeric(), list(), rewards, 3,
    1, 1, 1, 3, 1, c(0, 1), tolerance = 0)
  expect_equal(as.vector(empty$probability), c(1, 1))
  expect_equal(empty$lower, empty$upper)
})

test_that("grid refinement and failures are explicit", {
  z <- coalescing_grid_test_spec()
  args <- list(z$times, z$kernels, z$rewards, z$tau,
               0, 1, c(0, 0), c(1.8, 2.7), c(1, .8), seq(0, 5, .1))
  expect_error(do.call(markov_conditional_coverage_grid,
    c(args, list(grid_size = 1, exact_limit = 1, tolerance = 0,
                 max_refinements = 0))), class = "lifeconformal_grid_precision")
  got <- do.call(markov_conditional_coverage_grid,
    c(args, list(grid_size = 1, exact_limit = 1, tolerance = 0,
                 max_refinements = 0, strict = FALSE)))
  expect_false(got$diagnostics$tolerance_met)
  expect_gt(got$diagnostics$maximum_probability_error, 0)
  refined <- do.call(markov_conditional_coverage_grid,
    c(args, list(grid_size = 1, exact_limit = 1, tolerance = 0,
                 max_refinements = 8)))
  expect_true(refined$diagnostics$tolerance_met)
  expect_gt(refined$diagnostics$refinements, 0)
  exact <- do.call(markov_conditional_coverage_grid,
    c(args, list(exact_limit = 10000, tolerance = 0)))
  expect_equal(refined$probability, exact$probability, tolerance = 1e-13)
  expect_true(all(got$lower <= exact$probability + 1e-13))
  expect_true(all(got$upper >= exact$probability - 1e-13))
  expect_error(do.call(markov_conditional_coverage_grid,
    c(args, list(max_nodes = 3))), class = "lifeconformal_grid_budget")
})

test_that("one profile's batched query pass matches separate calls", {
  z <- grid_test_spec()
  times <- c(3, 0, 1.2, 1.2)
  states <- c(1, 2, 3, 1)
  accrued <- cbind(times / 2, times / 2)
  centers <- matrix(c(2, 3), 4, 2, byrow = TRUE)
  scales <- matrix(1, 4, 2)
  q <- seq(0, 4, .25)
  many <- markov_conditional_coverage_grid(z$times, z$kernels, z$rewards,
    z$tau, times, states, accrued, centers, scales, q,
    grid_size = 2, exact_limit = 2, tolerance = 1)
  # Different query times can create different floating additions but not a
  # different exact distribution; use exact mode to compare those recursions.
  exact <- markov_conditional_coverage_grid(z$times, z$kernels, z$rewards,
    z$tau, times, states, accrued, centers, scales, q,
    exact_limit = 10000, tolerance = 0)
  for (i in seq_along(times)) {
    one <- markov_conditional_coverage_grid(z$times, z$kernels, z$rewards,
      z$tau, times[i], states[i], accrued[i, ], centers[i, ], scales[i, ], q,
      exact_limit = 10000, tolerance = 0)
    expect_equal(exact$probability[i, ], as.vector(one$probability),
                 tolerance = 1e-13)
    expect_true(all(many$lower[i, ] <= one$probability + 1e-13))
    expect_true(all(many$upper[i, ] >= one$probability - 1e-13))
  }
})

test_that("constant and negative rewards and joint dependence are retained", {
  rewards <- rbind(c(-1, 2), c(1, 2))
  kernel <- matrix(c(.5, .5, 0, 1), 2, byrow = TRUE)
  got <- markov_conditional_coverage_grid(1, list(kernel), rewards, 2,
    0, 1, c(0, 0), c(-1, 4), c(1, 1), c(0, 1, 2),
    exact_limit = 1, grid_size = 2, tolerance = 0)
  expect_equal(as.vector(got$probability), c(0, 1, 1))
  z <- grid_test_spec()
  clone <- unserialize(serialize(z, NULL))
  set.seed(512)
  before <- .Random.seed
  do.call(markov_conditional_coverage_grid,
    list(clone$times, clone$kernels, clone$rewards, clone$tau,
         0, 1, c(0, 0), c(2, 3), c(1, 1), c(0, 1, 2), tolerance = 0))
  expect_identical(.Random.seed, before)
})

test_that("grid inputs and inclusive weighted CDF ties are checked", {
  expect_equal(weighted_markov_cdf(c(2, 1, 1), c(.2, .3, .5),
                                  c(0, 1, 1.5, 2)), c(0, .8, .8, 1))
  z <- grid_test_spec()
  args <- list(z$times, z$kernels, z$rewards, z$tau,
               0, 1, c(0, 0), c(2, 3), c(1, 1), c(0, 1))
  invalid <- args; invalid[[9]] <- c(0, 1)
  expect_error(do.call(markov_conditional_coverage_grid, invalid), "scale")
  invalid <- args; invalid[[5]] <- 5
  expect_error(do.call(markov_conditional_coverage_grid, invalid), "query_times")
  invalid <- args; invalid[[6]] <- 0
  expect_error(do.call(markov_conditional_coverage_grid, invalid), "query_state")
  invalid <- args; invalid[[10]] <- c(1, 0)
  expect_error(do.call(markov_conditional_coverage_grid, invalid), "q")
})

test_that("reward-equivalent switching preserves boundary atom probabilities", {
  rewards <- matrix(c(1, 1, 0), 3, 1)
  kernel <- matrix(c(.7, .2, .1, .3, .6, .1, 0, 0, 1), 3, byrow = TRUE)
  times <- seq_len(8) / 8
  got <- markov_conditional_coverage_grid(times, rep(list(kernel), 8),
    rewards, 1.5, 0, 1, 0, 1.5, 1, c(0, .1, 1),
    exact_limit = 1, grid_size = 2, tolerance = 1, max_refinements = 0)
  expect_equal(got$lower[1, 1], .9^8, tolerance = 1e-13)
  expect_equal(got$upper[1, 1], .9^8, tolerance = 1e-13)
})

test_that("moment and support probability brackets hold across varied kernels", {
  # Fixed deterministic probabilities avoid coupling this reference check to
  # RNG behavior while testing heterogeneous simultaneous transition kernels.
  for (b in seq_len(8)) {
    kernels <- lapply(seq_len(4), function(j) {
      p <- matrix(abs(sin(seq_len(9) + b * j)), 3, 3)
      p / rowSums(p)
    })
    rewards <- rbind(c(-1, 1), c(.25, 1), c(1, 0))
    args <- list(c(.3, .7, 1.1, 1.6), kernels, rewards, 2,
                 c(0, .5), c(1, 2), matrix(c(0, .1, 0, .3), 2),
                 matrix(c(.3, .4, 1, 1.2), 2), matrix(c(.7, .5, 1, 1), 2),
                 seq(0, 4, .025))
    reference <- do.call(markov_conditional_coverage_grid,
      c(args, list(exact_limit = 10000, tolerance = 0)))
    coarse <- do.call(markov_conditional_coverage_grid,
      c(args, list(grid_size = 3, exact_limit = 1, tolerance = 1,
                   max_refinements = 0)))
    expect_true(all(coarse$lower <= reference$probability + 1e-12))
    expect_true(all(coarse$upper >= reference$probability - 1e-12))
    expect_true(all(coarse$probability >= coarse$lower - 1e-12))
    expect_true(all(coarse$probability <= coarse$upper + 1e-12))
    expect_true(all(apply(coarse$probability, 1, diff) >= -1e-12))
  }
})

test_that("streamed reducers equal full matrices and restart on refinement", {
  z <- coalescing_grid_test_spec()
  q <- seq(0, 5, .1)
  args <- list(z$times, z$kernels, z$rewards, z$tau,
               c(0, 1), c(1, 2), matrix(0, 2, 2),
               matrix(c(1.8, 2.7), 2, 2, byrow = TRUE),
               matrix(c(1, .8), 2, 2, byrow = TRUE), q,
               grid_size = 1, exact_limit = 1, tolerance = 0, max_refinements = 8)
  weights <- c(1, -.5)
  init_count <- 0L
  reducer <- list(
    initialize = function() { init_count <<- init_count + 1L; numeric(length(q)) },
    add = function(state, index, probability, lower, upper)
      state + weights[index] * probability,
    finalize = identity)
  full <- do.call(markov_conditional_coverage_grid, args)
  streamed <- do.call(markov_conditional_coverage_grid, c(args, list(reducer = reducer)))
  expect_equal(streamed$reduced, as.vector(weights %*% full$probability))
  expect_null(streamed$probability)
  expect_null(streamed$lower)
  expect_null(streamed$upper)
  expect_equal(init_count, streamed$diagnostics$refinements + 1L)
  expect_gt(init_count, 1L)
  exact_args <- args
  exact_args$exact_limit <- 10000
  exact <- do.call(markov_conditional_coverage_grid, exact_args)
  expect_equal(full$probability, exact$probability, tolerance = 1e-13)
  expect_true(all(full$lower <= exact$probability + 1e-13))
  expect_true(all(full$upper >= exact$probability - 1e-13))
})

test_that("structural override hook precedes precision assessment", {
  z <- coalescing_grid_test_spec()
  override <- function(index, value) {
    value$probability[] <- value$lower[] <- value$upper[] <- 1
    value
  }
  args <- list(z$times, z$kernels, z$rewards, z$tau, 0, 1, c(0, 0),
    c(1.8, 2.7), c(1, .8), seq(0, 5, .1), grid_size = 1, exact_limit = 1,
    tolerance = 0, max_refinements = 0)
  unresolved <- do.call(markov_conditional_coverage_grid,
    c(args, list(strict = FALSE)))
  expect_gt(unresolved$diagnostics$maximum_probability_error, 0)
  expect_false(unresolved$diagnostics$tolerance_met)
  got <- do.call(markov_conditional_coverage_grid,
    c(args, list(query_transform = override)))
  expect_true(got$diagnostics$tolerance_met)
  expect_true(all(got$probability == 1))
  expect_error(markov_conditional_coverage_grid(z$times, z$kernels, z$rewards,
    z$tau, 0, 1, c(0, 0), c(1.8, 2.7), c(1, .8), seq(0, 5, .1),
    query_transform = function(index, value) list(probability = 2)),
    "invalid probability")
})
