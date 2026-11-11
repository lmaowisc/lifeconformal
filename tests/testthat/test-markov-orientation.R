test_that("binary reward directions are detected componentwise", {
  rewards <- cbind(alive = c(1, 1, 0, 0), death_a = c(0, 0, 1, 0),
                   death_b = c(0, 0, 0, 1), always = rep(1, 4), never = rep(0, 4))
  edges <- data.frame(from = c(1, 1, 2, 2), to = c(2, 3, 3, 4))
  got <- orient_markov_rewards(rewards, edges)
  expect_equal(got$flip, c(alive = FALSE, death_a = TRUE, death_b = TRUE,
                          always = FALSE, never = FALSE))
  expect_true(markov_monotone_eligible(got$rewards, as.matrix(edges)))
  expect_equal(got$rewards[, "death_a"], 1 - rewards[, "death_a"])
  expect_equal(got$rewards[, "alive"], rewards[, "alive"])
  expect_false(markov_monotone_eligible(rewards, as.matrix(edges)))
  cycle <- rbind(edges, data.frame(from = 3, to = 1))
  expect_null(orient_markov_rewards(rewards, cycle))
  invalid <- rewards; invalid[1, 1] <- -.5
  expect_null(orient_markov_rewards(invalid, edges))
  empty <- edges[FALSE, ]
  expect_false(any(orient_markov_rewards(rewards, empty)$flip))
})

test_that("mixed alive and absorbing occupations preserve exact joint coverage", {
  rewards <- cbind(alive = c(1, 1, 0, 0), death_a = c(0, 0, 1, 0),
                   death_b = c(0, 0, 0, 1), always = rep(1, 4))
  edges <- data.frame(from = c(1, 1, 2, 2), to = c(2, 3, 3, 4))
  jump_times <- c(.5, 1.5, 2.5)
  kernels <- lapply(seq_along(jump_times), function(j) {
    increments <- matrix(0, 4, 4)
    increments[as.matrix(edges)] <- c(.2, .1, .3, .15) * j / 2
    make_markov_kernel(increments)$kernel
  })
  time <- c(0, 1, 2, 3)
  state <- c(1, 2, 3, 4)
  accrued <- rbind(c(0, 0, 0, 0), c(1, 0, 0, 1),
                    c(1.5, .5, 0, 2), c(2.5, 0, .5, 3))
  center <- matrix(c(2.5, 1, .5, 4), 4, 4, byrow = TRUE)
  scale <- matrix(c(.8, .6, .7, 1), 4, 4, byrow = TRUE)
  q <- seq(0, 5, .125)
  orientation <- orient_markov_rewards(rewards, edges)
  transformed <- orient_markov_query(orientation, center, accrued, time, tau = 4)
  clock <- markov_conditional_coverage_monotone(jump_times, kernels,
    transformed$rewards, 4, time, state, transformed$accrued,
    transformed$center, scale, q, flip = transformed$flip,
    original_center = transformed$original_center,
    original_accrued = transformed$original_accrued)
  for (i in seq_along(time)) {
    atoms <- propagate_markov_rewards(jump_times, kernels, rewards, 4,
      start = time[i], state = state[i], accrued = accrued[i, ])$distributions[[1]]
    scores <- apply(abs(sweep(sweep(atoms$values, 2, center[i, ], "-"),
                              2, scale[i, ], "/")), 1, max)
    truth <- vapply(q, function(x) sum(atoms$probability[scores <= x]), numeric(1))
    expect_equal(clock$probability[i, ], truth, tolerance = 1e-13)
    # The same transformation preserves every joint atom's score, including
    # zero/terminal atoms and overlapping/complementary target components.
    reflected <- atoms$values
    reflected[, orientation$flip] <- 4 - reflected[, orientation$flip]
    score2 <- apply(abs(sweep(sweep(reflected, 2, transformed$center[i, ], "-"),
                              2, scale[i, ], "/")), 1, max)
    expect_equal(scores, score2, tolerance = 1e-13)
  }
})

test_that("canonical complemented scores preserve inclusive floating atom ties", {
  rewards <- matrix(c(0, 1), 2, 1)
  edges <- data.frame(from = 1, to = 2)
  orientation <- orient_markov_rewards(rewards, edges)
  original_center <- matrix(.1, 1, 1)
  original_accrued <- matrix(0, 1, 1)
  transformed <- orient_markov_query(orientation, original_center,
    original_accrued, 0, tau = 1)
  kernel <- matrix(c(.5, .5, 0, 1), 2, byrow = TRUE)
  q <- c(.6, .7, .8)
  clock <- markov_conditional_coverage_monotone(.2, list(kernel), transformed$rewards,
    1, 0, 1, transformed$accrued, transformed$center, matrix(1, 1, 1), q,
    flip = transformed$flip, original_center = original_center,
    original_accrued = original_accrued)
  native <- markov_monotone_batch(.2, matrix(.5, 1, 1), edges, matrix(1, 1, 1),
    transformed$rewards, 1, transformed$center, matrix(1, 1, 1), q,
    1, 0, 1, transformed$accrued, 1, flip = transformed$flip,
    original_center = original_center, original_accrued = original_accrued)
  atoms <- propagate_markov_rewards(.2, list(kernel), rewards, 1,
    state = 1)$distributions[[1]]
  score <- abs(atoms$values[, 1] - .1)
  truth <- vapply(q, function(x) sum(atoms$probability[score <= x]), numeric(1))
  expect_equal(truth, c(.5, .5, 1))
  expect_equal(as.vector(clock$probability), truth)
  expect_equal(as.vector(native$contributions), truth)
  # A post-entry query must retain its original accrued value; calculating
  # time - (time - original_accrued) need not reproduce the same double.
  time <- .3
  original_accrued <- matrix(time - .2, 1, 1)
  transformed <- orient_markov_query(orientation, original_center,
    original_accrued, time, tau = 1)
  # Exact binary64-input oracle: float(Fraction(1) - Fraction(.3) +
  # Fraction(.3-.2)) is .8 (0x1.999999999999ap-1). Rounding 1-time
  # first gives the previous double and moves this inclusive boundary.
  truth_target <- .8
  q <- sort(unique(c(.6, .7, abs(truth_target - .1), .8)))
  native <- markov_monotone_batch(.2, matrix(.5, 1, 1), edges, matrix(1, 1, 1),
    transformed$rewards, 1, transformed$center, matrix(1, 1, 1), q,
    1, time, 2, transformed$accrued, 1, flip = transformed$flip,
    original_center = original_center, original_accrued = original_accrued)
  expect_equal(as.vector(native$contributions), as.numeric(abs(truth_target - .1) <= q))
})

test_that("orientation keeps active observed clocks coherent", {
  rewards <- cbind(alive = c(1, 1, 0), dead = c(0, 0, 1))
  edges <- data.frame(from = c(1, 1, 2), to = c(2, 3, 3))
  orientation <- orient_markov_rewards(rewards, edges)
  times <- c(0, 1, 2.5, 4)
  states <- c(1, 2, 3, 3)
  accrued <- cbind(alive = c(0, 1, 2.5, 2.5), dead = c(0, 0, 0, 1.5))
  got <- orient_markov_query(orientation, matrix(c(2, 3), 1),
                              accrued, times, tau = 5)
  active <- got$rewards[states, , drop = FALSE] == 1
  expect_equal(got$accrued[active], matrix(times, 4, 2)[active])
  expect_equal(got$center, matrix(c(2, 2), 1))
  expect_equal(got$accrued[, 1], got$accrued[, 2])
  expect_equal(nrow(got$center), 1L)
  expect_equal(nrow(got$accrued), 4L)
  baseline <- orient_markov_query(orientation, matrix(c(2, 3), 1), tau = 5)
  expect_null(baseline$accrued)
})

test_that("orientation validation handles empty and singleton matrices", {
  rewards <- matrix(c(0, 1), 2, 1, dimnames = list(c("a", "b"), "entered"))
  edges <- data.frame(from = 1, to = 2)
  orientation <- orient_markov_rewards(rewards, edges)
  one <- orient_markov_query(orientation, matrix(1, 1, 1), matrix(.25, 1, 1), .5, 2)
  expect_equal(one$accrued, matrix(.25, 1, 1))
  expect_equal(one$rewards, matrix(c(1, 0), 2, 1, dimnames = dimnames(rewards)))
  empty <- orient_markov_query(orientation, matrix(numeric(), 0, 1),
    matrix(numeric(), 0, 1), numeric(), 2)
  expect_equal(dim(empty$accrued), c(0L, 1L))
  expect_error(orient_markov_rewards(rewards, data.frame(from = 0, to = 2)), "indices")
  expect_error(orient_markov_query(orientation, matrix(1, 1), tau = 0), "horizon")
  expect_error(orient_markov_query(orientation, matrix(1, 1), time = 1, tau = 2), "only")
})
