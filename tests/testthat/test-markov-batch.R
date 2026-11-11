batch_reference <- function(times, base, edges, risk, rewards, tau, initial,
                             accrued = NULL) {
  n <- nrow(risk)
  k <- nrow(rewards)
  d <- ncol(rewards)
  if (length(initial) == 1L) initial <- rep(initial, n)
  if (is.null(accrued)) accrued <- matrix(0, n, d)
  if (is.null(dim(accrued))) accrued <- matrix(rep(accrued, each = n), n, d)
  mean <- variance <- matrix(0, n, d)
  state_mean <- state_variance <- array(0, c(n, k, d))
  for (i in seq_len(n)) {
    kernels <- lapply(seq_along(times), function(j) {
      increments <- matrix(0, k, k)
      increments[cbind(edges[, "from"], edges[, "to"])] <- base[j, ] * risk[i, ]
      make_markov_kernel(increments)$kernel
    })
    reference <- markov_reward_moments(times, kernels, rewards, tau,
                                       accrued = accrued[i, ])
    mean[i, ] <- reference$mean[initial[i], ]
    variance[i, ] <- reference$variance[initial[i], ]
    state_mean[i, , ] <- reference$mean
    state_variance[i, , ] <- reference$variance
  }
  list(mean = mean, variance = variance,
       state_mean = state_mean, state_variance = state_variance)
}

test_that("batch moments match separate reference profiles with simultaneous exits", {
  edges <- cbind(from = c(1, 1, 2, 2), to = c(2, 3, 1, 3))
  times <- c(1, 2, 3, 6, 7)
  base <- matrix(c(.2, .6, .1, 0, .2, .5, .5, .1, 0, .1,
                   .2, .9, .4, 0, .3, .1, .5, .2, 0, .1), 5, 4)
  risk <- matrix(c(.4, 1, 1.8, .8, 1.1, .9, 3, .5, 1.1, 1.3, .8, 1.7), 3, 4)
  rewards <- cbind(home = c(1, 0, 0), alive = c(1, 1, 0))
  rownames(rewards) <- c("home", "hospital", "death")
  rownames(risk) <- c("A", "B", "C")
  colnames(base) <- colnames(risk) <- c("12", "13", "21", "23")
  accrued <- cbind(home = c(1, 2, 0), alive = c(3, 5, 0))
  initial <- c(1, 2, 3)
  actual <- markov_reward_moments_batch(times, base, edges, risk, rewards, 8,
                                        initial, accrued, absorbing = 3,
                                        keep_state_moments = TRUE)
  expected <- batch_reference(times, base, edges, risk, rewards, 8, initial, accrued)
  expect_equal(unname(actual$mean), expected$mean, tolerance = 1e-12)
  expect_equal(unname(actual$variance), expected$variance, tolerance = 1e-12)
  expect_equal(unname(actual$state_mean), expected$state_mean, tolerance = 1e-12)
  expect_equal(unname(actual$state_variance), expected$state_variance, tolerance = 1e-12)
  expect_equal(actual$second_moment, actual$variance + actual$mean^2)
  expect_equal(dim(actual$mean), c(3L, 2L))
  expect_equal(dimnames(actual$mean), list(rownames(risk), colnames(rewards)))
  expect_true(sum(actual$diagnostics$adjusted_by_profile) > 0)
  expect_equal(sum(actual$diagnostics$adjusted_by_profile),
                sum(actual$diagnostics$adjusted_by_state))
})

test_that("sparse jump traversal and repeated initial states preserve exact moments", {
  edges <- data.frame(from = c(1, 1, 2, 2), to = c(2, 3, 1, 3))
  times <- seq_len(12)
  base <- matrix(0, 12, 4)
  base[cbind(c(1, 4, 6, 10), c(1, 3, 4, 2))] <- c(.2, .5, .8, .4)
  risk <- matrix(seq(.4, 1.9, length.out = 20), 5, 4)
  rewards <- cbind(c(1, 0, 0), c(1, 1, 0), c(-1, 2, 4))
  initial <- rep(1, 5)
  actual <- markov_reward_moments_batch(times, base, edges, risk, rewards, 13,
                                        initial, accrued = c(1e6, 2e6, 1e5))
  expected <- batch_reference(times, base, edges, risk, rewards, 13, initial,
                               c(1e6, 2e6, 1e5))
  expect_equal(unname(actual$mean), expected$mean, tolerance = 1e-12)
  expect_equal(unname(actual$variance), expected$variance, tolerance = 1e-12)
  expect_equal(actual$diagnostics$active_jump_count, 4L)
})

test_that("state updates are simultaneous rather than in-place", {
  edges <- cbind(from = c(1, 2), to = c(2, 1))
  rewards <- matrix(c(1, 0), 2, 1)
  actual <- markov_reward_moments_batch(1, matrix(c(1, 1), 1), edges,
                                        matrix(1, 2, 2), rewards, 3, c(1, 2),
                                        keep_state_moments = TRUE)
  expect_equal(drop(actual$mean), c(1, 2))
  expect_equal(drop(actual$variance), c(0, 0))
  expect_equal(drop(actual$state_mean[1, , ]), c(1, 2))
})

test_that("lazy shifts retain tiny occupation times at a huge horizon", {
  edges <- cbind(from = 1, to = 2)
  actual <- markov_reward_moments_batch(1e-6, matrix(1, 1), edges,
                                        matrix(1, 1), matrix(c(1, 0), 2),
                                        1e12, 1, absorbing = 2)
  expect_equal(drop(actual$mean), 1e-6, tolerance = 1e-15)
  expect_equal(drop(actual$variance), 0)
  enormous <- markov_reward_moments_batch(1, matrix(1, 1), edges,
                                          matrix(1, 1), matrix(c(1, 0), 2),
                                          1e200, 1, absorbing = 2)
  expect_equal(drop(enormous$mean), 1)
  expect_equal(drop(enormous$variance), 0)
})

test_that("no jumps, no edges, single dimensions and final-interval rewards work", {
  edges <- matrix(integer(), 0, 2, dimnames = list(NULL, c("from", "to")))
  rewards <- matrix(2, 1, 1, dimnames = list("absorbed", "target"))
  actual <- markov_reward_moments_batch(numeric(), matrix(numeric(), 0, 0),
                                        edges, matrix(numeric(), 1, 0), rewards,
                                        4, 1, accrued = 3, absorbing = 1,
                                        keep_state_moments = TRUE)
  expect_equal(unname(actual$mean), matrix(11, 1, 1))
  expect_equal(unname(actual$variance), matrix(0, 1, 1))
  expect_equal(dim(actual$state_mean), c(1L, 1L, 1L))
  expect_equal(drop(actual$second_moment), c(target = 121))
  expect_equal(actual$diagnostics$active_jump_count, 0L)
  zero_rows <- markov_reward_moments_batch(c(1, 2), matrix(numeric(), 2, 0),
                                           edges, matrix(numeric(), 2, 0), rewards,
                                           4, c(1, 1), keep_state_moments = TRUE)
  expect_equal(unname(zero_rows$mean), matrix(8, 2, 1))
})

test_that("a jump at tau cannot change a time-integrated reward", {
  edges <- cbind(from = 1, to = 2)
  actual <- markov_reward_moments_batch(4, matrix(1, 1), edges,
                                        matrix(2, 1), matrix(c(1, 0), 2), 4, 1)
  expect_equal(drop(actual$mean), 4)
  expect_equal(drop(actual$variance), 0)
})

test_that("batch numerical and structural input errors are explicit", {
  edges <- cbind(from = c(1, 2), to = c(2, 3))
  base <- matrix(c(.2, .3), 1)
  risk <- matrix(1, 2, 2)
  rewards <- matrix(c(1, 1, 0), 3)
  run <- function(...) {
    arguments <- list(jump_times = 1, base_increments = base, edges = edges,
                       relative_risk = risk, rewards = rewards, tau = 3,
                       initial_state = c(1, 1))
    modifications <- list(...)
    arguments[names(modifications)] <- modifications
    do.call(markov_reward_moments_batch, arguments)
  }
  expect_error(run(relative_risk = matrix(0, 2, 2)), "positive finite")
  expect_error(run(relative_risk = matrix(1e308, 2, 2),
                   base_increments = matrix(1e308, 1, 2)), "overflow")
  expect_error(run(edges = cbind(from = c(1, 1), to = c(2, 2))), "distinct")
  expect_error(run(edges = cbind(from = c(1, 2), to = c(1, 3))), "off-diagonal")
  expect_error(run(base_increments = -base), "nonnegative")
  expect_error(run(jump_times = matrix(1, 1)), "strictly increasing")
  expect_error(run(initial_state = c(1, 1, 1)), "one per profile")
  expect_error(run(absorbing = 2), "absorbing states")
  expect_error(run(accrued = matrix(1, 1, 2)), "profile-by-target")
  expect_error(run(keep_state_moments = NA), "TRUE or FALSE")
  colnames(base) <- c("12", "23")
  colnames(risk) <- c("23", "12")
  expect_error(run(base_increments = base, relative_risk = risk), "same edge order")
})

test_that("batch moments do not modify the random stream", {
  set.seed(135)
  before <- .Random.seed
  markov_reward_moments_batch(1, matrix(.2, 1), cbind(from = 1, to = 2),
                               matrix(1, 3, 1), matrix(c(1, 0), 2), 3, 1)
  expect_identical(.Random.seed, before)
})
