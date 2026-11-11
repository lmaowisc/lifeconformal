test_that("protected totals preserve inclusive half-ULP atoms", {
  tau <- 1 + .Machine$double.eps
  time <- .Machine$double.eps / 2
  tiny <- .Machine$double.xmin * .Machine$double.eps
  expect_identical(markov_protected_total(tau, time, c(1, 0), c(time, tiny)),
                   c(tau, tiny))
  expect_identical(markov_protected_total(tau, time, 1, tiny), tau)
  # The exact halfway result rounds to the even lower endpoint.
  expect_identical(markov_protected_total(tau, time, 1, 0), 1)
  rewards <- matrix(c(1, 0), 2)
  atoms <- propagate_markov_rewards(rewards = rewards, tau = tau,
                                    start = time, accrued = time)
  expect_identical(unname(atoms$distributions[[1]]$values[1, 1]), tau)
  expect_identical(unname(atoms$distributions[[2]]$values[1, 1]), time)
  moments <- markov_reward_moments(rewards = rewards, tau = tau,
                                   start = time, accrued = time)
  expect_identical(unname(moments$mean[, 1]), c(tau, time))
  expect_identical(as.vector(moments$variance), c(0, 0))
  terminal <- list(matrix(c(0, 1, 1, 0), 2))
  atoms <- propagate_markov_rewards(tau, terminal, rewards, tau,
                                    start = time, accrued = time)
  moments <- markov_reward_moments(tau, terminal, rewards, tau,
                                   start = time, accrued = time)
  expect_identical(unname(atoms$distributions[[1]]$values[1, 1]), tau)
  expect_identical(unname(moments$mean[1, 1]), tau)
  expect_identical(as.vector(moments$variance), c(0, 0))
})

test_that("no-change component provenance survives reward-equivalent jumps", {
  times <- c(.1, .3, .6, 1.2, 2.1, 2.7, 3.2)
  kernel <- matrix(c(.4, .6, .7, .3), 2, byrow = TRUE)
  kernels <- rep(list(kernel), length(times))
  for (unit in c(1, 1e-155)) {
    # The first component is known constant even when the second changes.
    rewards <- cbind(c(1, 1), c(1, 0))
    constant <- markov_reward_moments(times * unit, kernels,
      rewards[, 1, drop = FALSE], 3.3 * unit)
    expect_identical(as.vector(constant$mean), rep(3.3 * unit, 2))
    expect_identical(as.vector(constant$variance), c(0, 0))
  }
  rewards <- cbind(c(1, 1), c(1, 0))
  moments <- markov_reward_moments(times, kernels, rewards, 3.3)
  expect_identical(unname(moments$mean[, 1]), c(3.3, 3.3))
  expect_identical(unname(moments$variance[, 1]), c(0, 0))
  expect_true(all(moments$variance[, 2] > 0))
  atoms <- propagate_markov_rewards(times, kernels, rewards, 3.3)
  for (distribution in atoms$distributions) {
    expect_true(all(distribution$values[, 1] == 3.3))
    expect_named(distribution, c("values", "probability"))
    expect_equal(sum(distribution$probability), 1)
  }
  # A last instantaneous departure earns no occupation in its destination.
  last <- propagate_markov_rewards(1, list(kernel), matrix(c(1, 0), 2), 1)
  expect_equal(nrow(last$distributions[[1]]$values), 1)
  expect_equal(last$distributions[[1]]$values[1, 1], 1)
  expect_equal(last$distributions[[1]]$probability, 1)
})

test_that("exact clocks preserve offset and complemented deterministic totals", {
  tau <- 1 + .Machine$double.eps
  time <- .Machine$double.eps / 2
  tiny <- .Machine$double.xmin * .Machine$double.eps
  direct <- markov_conditional_coverage_monotone(numeric(), list(), matrix(1),
    tau, time, 1, tiny, tau, 1, 0)
  expect_identical(as.vector(direct$probability), 1)
  # A known exit uses the same accurate duration-plus-accrued clock.
  exit <- markov_conditional_coverage_monotone(tau,
    list(matrix(c(0, 0, 1, 1), 2)), matrix(c(1, 0)),
    2, time, 1, tiny, tau, 1, 0)
  expect_identical(as.vector(exit$probability), 1)
  complement <- markov_conditional_coverage_monotone(numeric(), list(),
    matrix(c(1, 0)), tau, time, 2, 0, 0, 1, 0,
    flip = TRUE, original_center = tau, original_accrued = time)
  expect_identical(as.vector(complement$probability), 1)
  batch <- markov_monotone_batch(numeric(), matrix(numeric(), 0, 0),
    data.frame(from = integer(), to = integer()), matrix(numeric(), 1, 0),
    matrix(c(1, 0)), tau, matrix(0), matrix(1), 0,
    1, time, 2, matrix(0), 1, flip = TRUE,
    original_center = matrix(tau), original_accrued = matrix(time))
  expect_identical(as.vector(batch$contributions), 1)
})
