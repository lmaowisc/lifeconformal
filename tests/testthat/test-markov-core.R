test_that("intensity jumps give the specified adjusted stochastic kernel", {
  increments <- matrix(c(0, .8, .7, .2, 0, .1, 0, 0, 0), 3, byrow = TRUE)
  result <- make_markov_kernel(increments, absorbing = 3L)
  expect_equal(result$kernel[1, ], c(0, 8 / 15, 7 / 15))
  expect_equal(result$kernel[2, ], c(.2, .7, .1))
  expect_equal(result$kernel[3, ], c(0, 0, 1))
  expect_equal(rowSums(result$kernel), rep(1, 3))
  expect_identical(result$adjusted_rows, 1L)
  expect_equal(result$original_row_sums, c(1.5, .3, 0))
  expect_equal(result$row_scale, c(2 / 3, 1, 1))
  expect_error(make_markov_kernel(increments, absorbing = 1L), "absorbing")
  allowed <- increments > 0
  allowed[1, 2] <- FALSE
  expect_error(make_markov_kernel(increments, allowed = allowed), "prohibited")
  expect_error(make_markov_kernel(diag(2)), "zero diagonal")
  expect_error(make_markov_kernel(matrix(-1, 2, 2)), "nonnegative")
})

test_that("a single jump has the exact two-point occupation distribution", {
  p <- matrix(c(.7, .3, 0, 1), 2, byrow = TRUE)
  rewards <- rbind(home = c(home = 1, dead = 0),
                   dead = c(home = 0, dead = 1))
  moments <- markov_reward_moments(2, list(p), rewards, tau = 5)
  expect_equal(unname(moments$mean), rbind(c(4.1, .9), c(0, 5)))
  expect_equal(unname(moments$variance), rbind(c(1.89, 1.89), c(0, 0)))
  expect_equal(unname(moments$second_moment), rbind(c(18.7, 2.7), c(0, 25)))
  atoms <- propagate_markov_rewards(2, list(p), rewards, tau = 5)
  expect_named(atoms$distributions, c("home", "dead"))
  home <- atoms$distributions$home
  expect_equal(unname(home$values), rbind(c(5, 0), c(2, 3)))
  expect_equal(home$probability, c(.7, .3))
  expect_equal(rowSums(home$values), c(5, 5))
  expect_equal(atoms$distributions$dead$values,
               matrix(c(0, 5), 1, dimnames = list(NULL, c("home", "dead"))))
  expect_equal(atoms$distributions$dead$probability, 1)
})

test_that("query sides, accrued rewards and final intervals are respected", {
  p <- matrix(c(.7, .3, 0, 1), 2, byrow = TRUE)
  rewards <- matrix(c(1, 0), 2, dimnames = list(NULL, "alive"))
  left <- markov_reward_moments(2, list(p), rewards, 5,
                                start = 2, side = "left", accrued = 2)
  right <- markov_reward_moments(2, list(p), rewards, 5,
                                 start = 2, side = "right", accrued = 2)
  expect_equal(unname(left$mean[1, 1]), 4.1)
  expect_equal(unname(left$variance[1, 1]), 1.89)
  expect_equal(unname(right$mean[1, 1]), 5)
  expect_equal(unname(right$variance[1, 1]), 0)
  expect_equal(right$diagnostics$jump_count, 0)
  atoms <- propagate_markov_rewards(2, list(p), rewards, 5,
                                    start = 2, side = "left", accrued = 2,
                                    state = 1)
  expect_equal(as.vector(atoms$distributions[[1]]$values), c(5, 2))
  expect_equal(atoms$distributions[[1]]$probability, c(.7, .3))
  expect_equal(atoms$state, 1L)
  at_horizon <- markov_reward_moments(5, list(p), rewards, 5,
                                      start = 5, side = "left", accrued = 3)
  expect_equal(as.vector(at_horizon$mean), c(3, 3))
  expect_equal(as.vector(at_horizon$variance), c(0, 0))
  baseline <- markov_reward_moments(5, list(p), rewards, 5)
  expect_equal(as.vector(baseline$mean), c(5, 0))
  expect_equal(as.vector(baseline$variance), c(0, 0))
})

test_that("no jumps and one-state models remain deterministic", {
  rewards <- matrix(c(1, 2), 1, dimnames = list("only", c("a", "b")))
  moments <- markov_reward_moments(rewards = rewards, tau = 7,
                                   start = 2, accrued = c(3, 4))
  expect_equal(moments$mean, matrix(c(8, 14), 1, dimnames = dimnames(rewards)))
  expect_equal(moments$variance, matrix(0, 1, 2, dimnames = dimnames(rewards)))
  atoms <- propagate_markov_rewards(rewards = rewards, tau = 7,
                                    start = 2, accrued = c(3, 4))
  expect_equal(atoms$distributions[[1]]$values,
               matrix(c(8, 14), 1, dimnames = list(NULL, c("a", "b"))))
  expect_equal(atoms$distributions[[1]]$probability, 1)
})

test_that("an interior query ignores past jumps and shortens the next sojourn", {
  first <- matrix(c(0, 1, 0, 1), 2, byrow = TRUE)
  second <- matrix(c(.5, .5, 0, 1), 2, byrow = TRUE)
  rewards <- matrix(c(1, 0), 2)
  moments <- markov_reward_moments(c(1, 3), list(first, second), rewards, 5,
                                   start = 2, accrued = .75)
  expect_equal(moments$mean[1, 1], 2.75)
  expect_equal(moments$variance[1, 1], 1)
  atoms <- propagate_markov_rewards(c(1, 3), list(first, second), rewards, 5,
                                    start = 2, accrued = .75, state = 1)
  expect_equal(as.vector(atoms$distributions[[1]]$values), c(3.75, 1.75))
  expect_equal(atoms$distributions[[1]]$probability, c(.5, .5))
})

test_that("backward moments and atoms agree with independent path enumeration", {
  kernels <- list(matrix(c(.6, .3, .1, .2, .6, .2, 0, 0, 1), 3, byrow = TRUE),
                  matrix(c(.5, .2, .3, .1, .6, .3, 0, 0, 1), 3, byrow = TRUE))
  rewards <- rbind(c(1, 1), c(0, 1), c(0, 0))
  times <- c(1, 3)
  tau <- 5
  moments <- markov_reward_moments(times, kernels, rewards, tau)
  propagated <- propagate_markov_rewards(times, kernels, rewards, tau)
  # Enumerate the two post-jump states, accumulating rewards forward.
  paths <- expand.grid(after_first = 1:3, after_second = 1:3)
  for (initial in 1:3) {
    probability <- numeric(nrow(paths))
    values <- matrix(0, nrow(paths), 2)
    for (i in seq_len(nrow(paths))) {
      a <- paths$after_first[i]
      b <- paths$after_second[i]
      probability[i] <- kernels[[1]][initial, a] * kernels[[2]][a, b]
      values[i, ] <- rewards[initial, ] + 2 * rewards[a, ] + 2 * rewards[b, ]
    }
    reference_mean <- colSums(values * probability)
    reference_var <- colSums(sweep(values, 2, reference_mean)^2 * probability)
    expect_equal(as.vector(moments$mean[initial, ]), reference_mean)
    expect_equal(as.vector(moments$variance[initial, ]), reference_var)
    distribution <- propagated$distributions[[initial]]
    expect_equal(colSums(distribution$values * distribution$probability),
                 reference_mean)
    expect_equal(colSums(sweep(distribution$values, 2, reference_mean)^2 *
                          distribution$probability), reference_var)
    for (threshold in c(0, 1, 2, 3, 4, 5)) {
      reference_cdf <- sum(probability[values[, 1] <= threshold & values[, 2] <= 4])
      actual_cdf <- sum(distribution$probability[
        distribution$values[, 1] <= threshold & distribution$values[, 2] <= 4])
      expect_equal(actual_cdf, reference_cdf)
    }
  }
})

test_that("moment variances are stable under a large accrued offset", {
  p <- matrix(c(.7, .3, 0, 1), 2, byrow = TRUE)
  result <- markov_reward_moments(2, list(p), matrix(c(1, 0), 2), 5,
                                  accrued = 1e12)
  expect_equal(result$variance[1, 1], 1.89)
  expect_equal(result$variance[2, 1], 0)
})

test_that("atom merging is exact and the resource limit is enforced", {
  p <- matrix(c(.5, .5, .5, .5), 2)
  same <- propagate_markov_rewards(c(1, 2, 3), rep(list(p), 3),
                                   matrix(c(1, 1), 2), tau = 4, max_atoms = 2)
  expect_equal(same$diagnostics$peak_atoms, 2)
  expect_equal(nrow(same$distributions[[1]]$values), 1)
  expect_equal(same$distributions[[1]]$values[1, 1], 4)
  expect_equal(same$distributions[[1]]$probability, 1)
  near <- propagate_markov_rewards(1, list(p), matrix(c(1, 1 + 1e-8), 2), 2)
  expect_equal(nrow(near$distributions[[1]]$values), 2)
  expect_error(propagate_markov_rewards(1, list(p), matrix(c(1, 0), 2), 2,
                                        max_atoms = 3),
               class = "lifeconformal_atom_budget")
  expect_error(propagate_markov_rewards(rewards = matrix(c(1, 0), 2), tau = 2,
                                        max_atoms = 1),
               class = "lifeconformal_atom_budget")
})

test_that("invalid grids, kernels and query arguments fail explicitly", {
  p <- diag(2)
  rewards <- matrix(c(1, 0), 2)
  expect_error(markov_reward_moments(c(2, 1), list(p, p), rewards, 5),
               "strictly increasing")
  expect_error(markov_reward_moments(0, list(p), rewards, 5), "positive")
  expect_error(markov_reward_moments(6, list(p), rewards, 5), "tau")
  expect_error(markov_reward_moments(1, list(p, p), rewards, 5), "one matrix")
  expect_error(markov_reward_moments(1, list(p / 2), rewards, 5), "sum to one")
  expect_error(markov_reward_moments(1, list(diag(3)), rewards, 5), "state-by-state")
  expect_error(markov_reward_moments(1, list(p), rewards, 5, start = -1), "start")
  expect_error(markov_reward_moments(1, list(p), rewards, 5, accrued = c(1, 2)),
               "one value per target")
  expect_error(propagate_markov_rewards(1, list(p), rewards, 5, state = 3),
               "state indices")
  expect_error(propagate_markov_rewards(1, list(p), rewards, 5, max_atoms = 1.5),
               "positive integer")
  expect_error(markov_reward_moments(matrix(c(1, 3, 2, 4), 2),
                                     rep(list(p), 4), rewards, 5),
               "strictly increasing")
  expect_error(markov_reward_moments(1, list(p), rewards, 5, accrued = matrix(0)),
               "one value per target")
  expect_error(propagate_markov_rewards(1, list(p), rewards, 5, state = matrix(1)),
               "state indices")
  expect_error(propagate_markov_rewards(1, list(p), rewards, 5, max_atoms = matrix(10)),
               "positive integer")
  named_rewards <- rewards
  rownames(named_rewards) <- c("home", "dead")
  dimnames(p) <- list(c("dead", "home"), c("dead", "home"))
  expect_error(markov_reward_moments(1, list(p), named_rewards, 5), "state order")
})

test_that("deterministic propagation leaves the random stream unchanged", {
  set.seed(109)
  before <- .Random.seed
  p <- matrix(c(.7, .3, 0, 1), 2, byrow = TRUE)
  rewards <- matrix(c(1, 0), 2)
  markov_reward_moments(2, list(p), rewards, 5)
  result <- propagate_markov_rewards(2, list(p), rewards, 5)
  expect_identical(.Random.seed, before)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(result, path)
  expect_identical(readRDS(path), result)
})
