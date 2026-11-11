grid_query_subset <- function(args, i) {
  for (name in c("query_times", "query_state", "side")) args[[name]] <- args[[name]][i]
  for (name in c("accrued", "center", "scale")) args[[name]] <- args[[name]][i, , drop = FALSE]
  args
}

test_that("protected occupation totals retain cancellation and half-even residuals", {
  eps <- .Machine$double.eps
  tiny <- .Machine$double.xmin * eps
  expect_identical(markov_protected_total(1 + eps, eps / 2, 1, eps / 2), 1 + eps)
  expect_identical(markov_protected_total(1 + eps, eps / 2, 1, 0), 1)
  expect_identical(markov_protected_total(1 + eps, eps / 2, 1, eps / 4), 1 + eps)
  expect_identical(markov_protected_total(2 * tiny, tiny, 1, tiny), 2 * tiny)
  expect_identical(markov_protected_total(.Machine$double.xmax,
    .Machine$double.xmax, 1, .Machine$double.xmax), .Machine$double.xmax)
  expect_identical(markov_protected_total(1 + eps, eps / 2, c(0, 1),
    c(0, eps / 2)), c(0, 1 + eps))
  expect_identical(markov_protected_total(1.1, .2, .5, .1), .1 + (1.1 - .2) * .5)

  times <- c(0, eps / 2, 1 + eps)
  args <- list(jump_times = numeric(), kernels = list(),
    rewards = matrix(c(1, 0), 2L, 1L), tau = 1 + eps,
    query_times = times, query_state = rep(1L, 3L), accrued = matrix(times),
    center = matrix(1 + eps, 3L, 1L), scale = matrix(1, 3L, 1L), q = 0,
    side = rep("right", 3L), step = .1, origin = 0,
    max_nodes = 1000L, exact_limit = 256L)
  for (engine in list(markov_grid_pass, markov_grid_pass_native)) {
    result <- do.call(engine, args)
    expect_identical(result$probability, matrix(1, 3L, 1L))
    expect_identical(result$lower, result$probability)
    expect_identical(result$upper, result$probability)
    at_horizon <- args
    at_horizon$jump_times <- args$tau
    at_horizon$kernels <- list(matrix(c(0, 1, 1, 0), 2L, 2L))
    terminal <- do.call(engine, at_horizon)
    expect_identical(terminal$probability, result$probability)
    expect_identical(terminal$lower, result$lower)
    expect_identical(terminal$upper, result$upper)
  }
})

test_that("query cuts preserve intrinsic no-change atoms at non-dyadic times", {
  for (tau in c(1.1, 3.3)) {
    times <- sort(unique(c(seq(0, tau, by = .1), tau)))
    n <- length(times)
    args <- list(jump_times = numeric(), kernels = list(),
      rewards = matrix(c(1, 0), 2L, 1L), tau = tau,
      query_times = times, query_state = rep(1L, n), accrued = matrix(times),
      center = matrix(tau, n, 1L), scale = matrix(1, n, 1L), q = 0,
      side = rep("right", n), step = .1, origin = 0,
      max_nodes = 1000L, exact_limit = 256L)
    for (engine in list(markov_grid_pass, markov_grid_pass_native)) {
      full <- do.call(engine, args)
      expect_identical(full$probability, matrix(1, n, 1L))
      expect_identical(full$lower, full$probability)
      expect_identical(full$upper, full$probability)
      permutation <- rev(seq_len(n))
      reordered <- do.call(engine, grid_query_subset(args, permutation))
      expect_identical(reordered$probability, full$probability[permutation, , drop = FALSE])
    }
  }
})

test_that("ordinary recurrent atoms are independent of other query times", {
  tau <- 1.1
  jumps <- c(.1, .2, .4, .7)
  swap <- matrix(c(0, 1, 1, 0), 2L, 2L)
  times <- sort(unique(c(seq(0, tau, by = .03), jumps, tau)))
  # This observed path alternates states at every fitted jump. Its elapsed
  # reward is computed chronologically, independently of either grid engine.
  past <- function(t, side) {
    included <- jumps[jumps < t | (jumps == t & side == "right")]
    boundaries <- c(0, included, t)
    state <- if (length(included) %% 2L) 2L else 1L
    reward <- sum(diff(boundaries) * rep(c(1, 0), length.out = length(boundaries) - 1L))
    c(state = state, reward = reward)
  }
  queries <- rep(times, each = 2L)
  sides <- rep(c("right", "left"), length(times))
  history <- vapply(seq_along(queries), function(i) past(queries[i], sides[i]), numeric(2))
  n <- length(queries)
  args <- list(jump_times = jumps, kernels = rep(list(swap), length(jumps)),
    rewards = matrix(c(1, 0), 2L, 1L), tau = tau,
    query_times = queries, query_state = history["state", ],
    accrued = matrix(history["reward", ]), center = matrix(.7, n, 1L),
    scale = matrix(1, n, 1L), q = c(0, .Machine$double.eps, .1, 2),
    side = sides, step = .1, origin = 0, max_nodes = 1000L, exact_limit = 256L)
  for (engine in list(markov_grid_pass, markov_grid_pass_native)) {
    full <- do.call(engine, args)
    singles <- lapply(seq_len(n), function(i) do.call(engine, grid_query_subset(args, i)))
    for (field in c("probability", "lower", "upper")) {
      expect_identical(full[[field]], do.call(rbind, lapply(singles, `[[`, field)))
      permutation <- rev(seq_len(n))
      reordered <- do.call(engine, grid_query_subset(args, permutation))
      expect_identical(reordered[[field]], full[[field]][permutation, , drop = FALSE])
    }
  }
})

test_that("reward-equivalent fitted jumps preserve the canonical intrinsic atom", {
  tau <- 3.3
  jumps <- seq(.1, 3.2, by = .1)
  kernel <- rbind(c(.4, .6, 0), c(.7, .3, 0), c(0, 0, 1))
  times <- sort(unique(c(seq(0, tau, by = .07), jumps, tau)))
  n <- length(times)
  args <- list(jump_times = jumps, kernels = rep(list(kernel), length(jumps)),
    rewards = matrix(c(1, 1, 0), 3L, 1L), tau = tau,
    query_times = times, query_state = rep(1L, n), accrued = matrix(times),
    center = matrix(tau, n, 1L), scale = matrix(1, n, 1L), q = 0,
    side = rep("right", n), step = .1, origin = 0,
    max_nodes = 1000L, exact_limit = 2L)
  for (engine in list(markov_grid_pass, markov_grid_pass_native)) {
    full <- do.call(engine, args)
    expect_equal(full$probability, matrix(1, n, 1L), tolerance = 1e-14)
    expect_equal(full$lower, full$probability, tolerance = 1e-14)
    expect_equal(full$upper, full$probability, tolerance = 1e-14)
  }
})

test_that("joint nodes preserve constant components through mixing and coalescence", {
  eps <- .Machine$double.eps
  tau <- 1 + eps
  time <- eps / 2
  rewards <- rbind(c(1, 1), c(1, 0), c(0, 0))
  for (hospital_probability in c(1, .5)) {
    kernel <- rbind(c(0, hospital_probability, 1 - hospital_probability),
      c(0, 1, 0), c(0, 0, 1))
    args <- list(jump_times = .5, kernels = list(kernel), rewards = rewards,
      tau = tau, query_times = time, query_state = 1,
      accrued = matrix(c(time, time), 1L, 2L),
      center = matrix(c(tau, .5), 1L, 2L), scale = matrix(1, 1L, 2L),
      q = 0, side = "right", step = c(10, 10), origin = c(0, 0),
      max_nodes = 1000L, exact_limit = 256L)
    for (engine in list(markov_grid_pass, markov_grid_pass_native)) {
      exact <- do.call(engine, args)
      expect_identical(exact$probability[1L, 1L], hospital_probability)
      expect_identical(exact$lower, exact$upper)
      # The two paths have different alive-time proofs. A large common grid
      # cell must not erase the canonical endpoint of the hospital path.
      coarse_args <- args
      coarse_args$exact_limit <- 1L
      coarse_args$point_mass_threshold <- 1
      coarse <- do.call(engine, coarse_args)
      expect_lte(coarse$lower[1L, 1L], hospital_probability)
      expect_gte(coarse$upper[1L, 1L], hospital_probability)
    }
  }
})
