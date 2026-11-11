monotone_atom_reference <- function(times, kernels, rewards, tau, query_times,
                                    states, accrued, center, scale, q, side) {
  t(vapply(seq_along(query_times), function(i) {
    atoms <- propagate_markov_rewards(times, kernels, rewards, tau,
      start = query_times[i], side = side[i], accrued = accrued[i, ],
      state = states[i])$distributions[[1L]]
    standardized <- abs(sweep(sweep(atoms$values, 2L, center[i, ], "-"),
                               2L, scale[i, ], "/"))
    score <- if (ncol(rewards) == 1L) as.vector(standardized) else
      apply(standardized, 1L, max)
    vapply(q, function(x) sum(atoms$probability[score <= x]), numeric(1))
  }, numeric(length(q))))
}

test_that("exact clocks retain event atoms and both conditioning sides", {
  kernels <- list(matrix(c(.6, .4, 0, 1), 2, byrow = TRUE),
                  matrix(c(.5, .5, 0, 1), 2, byrow = TRUE))
  rewards <- matrix(c(1, 0), 2, 1)
  got <- markov_conditional_coverage_monotone(c(1, 2), kernels, rewards, 3,
    c(0, 1, 1, 2, 3), c(1, 1, 1, 2, 1), matrix(c(0, 1, 1, 1, 3), 5, 1),
    matrix(2, 5, 1), matrix(1, 5, 1), c(0, .5, 1),
    side = c("right", "left", "right", "right", "left"))
  expect_equal(got$probability, rbind(c(.3, .3, 1), c(.3, .3, 1),
                                     c(.5, .5, 1), c(0, 0, 1), c(0, 0, 1)))
  expect_identical(got$probability, got$lower)
  expect_identical(got$probability, got$upper)
  expect_equal(got$diagnostics$deadline_groups, 1L)
  expect_equal(got$diagnostics$maximum_probability_error, 0)
})

test_that("joint progressive targets agree with exact path enumeration", {
  times <- c(1, 2, 3, 4)
  kernels <- list(matrix(c(.5, .3, .2, 0, .6, .4, 0, 0, 1), 3, byrow = TRUE),
                  matrix(c(.8, .1, .1, 0, .7, .3, 0, 0, 1), 3, byrow = TRUE),
                  matrix(c(.7, .2, .1, 0, .8, .2, 0, 0, 1), 3, byrow = TRUE),
                  matrix(c(.6, .2, .2, 0, .5, .5, 0, 0, 1), 3, byrow = TRUE))
  rewards <- rbind(c(1, 1), c(0, 1), c(0, 0))
  qt <- c(0, 1, 1, 1.5, 2, 2, 3.5, 4, 5)
  states <- c(1, 1, 2, 2, 1, 2, 3, 1, 2)
  side <- c("right", "left", "right", "right", "left", "right", "right", "left", "right")
  accrued <- cbind(c(0, 1, 1, 1, 2, 1, 1, 4, 2),
                    c(0, 1, 1, 1.5, 2, 2, 2, 4, 5))
  center <- matrix(c(2, 3), length(qt), 2, byrow = TRUE)
  scale <- matrix(c(1, 2), length(qt), 2, byrow = TRUE)
  q <- seq(0, 5, .125)
  got <- markov_conditional_coverage_monotone(times, kernels, rewards, 5,
    qt, states, accrued, center, scale, q, side)
  expected <- monotone_atom_reference(times, kernels, rewards, 5, qt, states,
                                      accrued, center, scale, q, side)
  expect_equal(got$probability, expected, tolerance = 1e-13)
  expect_equal(got$diagnostics$deadline_groups, 1L)
  for (i in seq_along(qt)) {
    one <- markov_conditional_coverage_monotone(times, kernels, rewards, 5,
      qt[i], states[i], accrued[i, ], center[i, ], scale[i, ], q, side[i])
    expect_equal(as.vector(one$probability), as.vector(got$probability[i, ]),
                 tolerance = 1e-13)
  }
})

test_that("reward-equivalent recurrent classes retain exact clock reduction", {
  kernels <- list(matrix(c(.5, .3, .2, .2, .5, .3, 0, 0, 1), 3, byrow = TRUE),
                  matrix(c(.6, .2, .2, .1, .7, .2, 0, 0, 1), 3, byrow = TRUE))
  rewards <- matrix(c(1, 1, 0), 3, 1)
  expect_true(markov_monotone_eligible(rewards, rbind(c(1, 2), c(2, 1), c(2, 3))))
  got <- markov_conditional_coverage_monotone(c(1, 2), kernels, rewards, 3,
    0, 1, 0, 3, 1, c(0, 1, 2, 3))
  # Staying anywhere in the alive class through both kernels retains tau atom.
  expect_equal(got$probability[1, 1], .5 * .8 + .3 * .8)
  reference <- monotone_atom_reference(c(1, 2), kernels, rewards, 3,
    0, 1, matrix(0, 1, 1), matrix(3, 1, 1), matrix(1, 1, 1), c(0, 1, 2, 3), "right")
  expect_equal(got$probability, reference)
})

test_that("general query accrual and scale differences get separate exact groups", {
  kernels <- list(matrix(c(.5, .5, 0, 1), 2, byrow = TRUE))
  rewards <- matrix(c(1, 0), 2, 1)
  qt <- c(0, .5, .5, 1, 1.5)
  states <- c(1, 1, 1, 1, 2)
  accrued <- matrix(c(0, .25, 0, 1, .25), 5, 1)
  center <- matrix(c(1, 1, 2, 1.5, 1), 5, 1)
  scale <- matrix(c(1, 1, 2, 1, 1), 5, 1)
  side <- c("right", "right", "right", "left", "right")
  q <- seq(0, 3, .125)
  got <- markov_conditional_coverage_monotone(1, kernels, rewards, 2,
    qt, states, accrued, center, scale, q, side)
  reference <- monotone_atom_reference(1, kernels, rewards, 2, qt, states,
                                      accrued, center, scale, q, side)
  expect_equal(got$probability, reference, tolerance = 1e-13)
  expect_equal(got$diagnostics$deadline_groups, 4L)
  perm <- c(4, 1, 5, 2, 3)
  reordered <- markov_conditional_coverage_monotone(1, kernels, rewards, 2,
    qt[perm], states[perm], accrued[perm, , drop = FALSE], center[perm, , drop = FALSE],
    scale[perm, , drop = FALSE], q, side[perm])
  expect_equal(reordered$probability, got$probability[perm, , drop = FALSE])
})

test_that("no jumps, final jumps and fully inactive targets are exact", {
  rewards <- rbind(c(1, 1), c(0, 1), c(0, 0))
  none <- markov_conditional_coverage_monotone(numeric(), list(), rewards, 3,
    c(0, 1, 2), 1:3, rbind(c(0, 0), c(.5, 1), c(.5, 1.5)),
    matrix(0, 3, 2), matrix(1, 3, 2), c(0, 1, 1.5, 3))
  expect_equal(none$probability, rbind(c(0, 0, 0, 1), c(0, 0, 0, 1), c(0, 0, 1, 1)))
  kernel <- matrix(c(.5, .2, .3, 0, .5, .5, 0, 0, 1), 3, byrow = TRUE)
  final <- markov_conditional_coverage_monotone(3, list(kernel), rewards, 3,
    c(0, 3, 3), c(1, 1, 1), rbind(c(0, 0), c(3, 3), c(3, 3)),
    matrix(3, 3, 2), matrix(1, 3, 2), c(0, 1), c("right", "left", "right"))
  expect_equal(final$probability, matrix(1, 3, 2))
  inactive <- markov_conditional_coverage_monotone(3, list(kernel), rewards, 3,
    2, 3, c(1, 2), c(1, 2), c(1, 1), 0)
  expect_identical(dim(inactive$probability), c(1L, 1L))
  expect_equal(inactive$probability[1, 1], 1)
  expect_equal(inactive$diagnostics$deadline_groups, 0L)
})

test_that("random small monotone systems match independent exact distributions", {
  set.seed(507)
  for (rep in seq_len(12)) {
    times <- seq_len(4)
    kernels <- lapply(times, function(i) {
      a <- stats::runif(3)
      a <- a / sum(a)
      b <- stats::runif(2)
      unname(rbind(a, c(0, b / sum(b)), c(0, 0, 1)))
    })
    rewards <- rbind(c(1, 1), c(0, 1), c(0, 0))
    qt <- c(0, 1, 1.5, 2, 3, 4)
    states <- sample(1:3, length(qt), replace = TRUE)
    accrued <- matrix(stats::runif(length(qt) * 2, 0, .75), length(qt), 2)
    center <- matrix(stats::runif(length(qt) * 2, 1, 4), length(qt), 2)
    scale <- matrix(stats::runif(length(qt) * 2, .5, 1.5), length(qt), 2)
    side <- sample(c("left", "right"), length(qt), replace = TRUE)
    q <- seq(0, 5, length.out = 17)
    got <- markov_conditional_coverage_monotone(times, kernels, rewards, 5,
      qt, states, accrued, center, scale, q, side)
    reference <- monotone_atom_reference(times, kernels, rewards, 5, qt, states,
                                        accrued, center, scale, q, side)
    expect_equal(got$probability, reference, tolerance = 1e-12)
  }
})

test_that("nonmonotone targets and invalid query contracts fail explicitly", {
  recurrent <- matrix(c(.5, .5, .5, .5), 2)
  rewards <- matrix(c(1, 0), 2, 1)
  expect_false(markov_monotone_eligible(rewards, rbind(c(1, 2), c(2, 1))))
  expect_false(markov_monotone_eligible(matrix(c(.5, 0), 2), matrix(c(1, 2), 1)))
  expect_error(markov_conditional_coverage_monotone(1, list(recurrent), rewards, 2,
    0, 1, 0, 1, 1, 0), "never increase")
  args <- list(numeric(), list(), rewards, 2, 0, 1, 0, 1, 1, 0)
  invalid <- args; invalid[[9]] <- 0
  expect_error(do.call(markov_conditional_coverage_monotone, invalid), "scale")
  invalid <- args; invalid[[5]] <- 3
  expect_error(do.call(markov_conditional_coverage_monotone, invalid), "query_times")
  invalid <- args; invalid[[6]] <- 3
  expect_error(do.call(markov_conditional_coverage_monotone, invalid), "query_state")
  invalid <- args; invalid[[10]] <- c(1, 0)
  expect_error(do.call(markov_conditional_coverage_monotone, invalid), "q must")
  expect_error(do.call(markov_conditional_coverage_monotone, c(args, list(side = "before"))),
               "side must")
})
