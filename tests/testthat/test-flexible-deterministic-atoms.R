flexible_atom_adapter <- function(tau) {
  structure(list(type = "transition_neural", seed = 47L, n_predictors = 0L,
    x_names = NULL, edges = data.frame(from = 1L, to = 2L),
    fits = list(list(constant = 0)),
    options = list(step_size = tau, max_events = 5L),
    problem = list(tau = tau, absorbing = 2L,
      rewards = matrix(c(0, 1, 1, 1, 0, 0), 2L, 3L,
        dimnames = list(NULL, c("terminal", "total", "zero"))))),
    class = "life_flexible_outcome")
}

test_that("deterministic flexible totals retain a half-ulp occupation atom", {
  tau <- 1 + .Machine$double.eps
  time <- .Machine$double.eps / 2
  expect_false(identical(time + (tau - time), tau))
  adapter <- flexible_atom_adapter(tau)
  snapshots <- list(time = c(time, time), state = c(2L, 1L),
    duration = c(time, time), previous_transitions = c(0L, 0L),
    accrued = rbind(c(time, time, 0), c(0, time, 0)))
  set.seed(271)
  rng <- .Random.seed
  stored <- serialize(adapter, NULL)
  draws <- simulate_flexible_targets(adapter, matrix(0, 2L, 0L), snapshots, B = 7L)
  # The terminal target varies by state: its known-terminal shortcut must be
  # accurate independently of the all-state constant-reward reset below it.
  expect_identical(draws[1L, , "terminal"], rep(tau, 7L))
  expect_identical(draws[2L, , "terminal"], rep(0, 7L))
  expect_identical(as.numeric(draws[, , "total"]), rep(tau, 14L))
  expect_identical(as.numeric(draws[, , "zero"]), rep(0, 14L))
  expect_identical(.Random.seed, rng)
  expect_identical(serialize(adapter, NULL), stored)
  expect_identical(draws, simulate_flexible_targets(adapter,
    matrix(0, 2L, 0L), snapshots, B = 7L))
  expect_identical(.Random.seed, rng)
})

test_that("terminal flexible totals leave an absent caller RNG absent", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  if (had_seed) rm(".Random.seed", envir = .GlobalEnv)
  tau <- 1 + .Machine$double.eps
  time <- .Machine$double.eps / 2
  adapter <- flexible_atom_adapter(tau)
  draws <- simulate_flexible_targets(adapter, matrix(0, 1L, 0L),
    list(time = time, state = 2L, duration = time, previous_transitions = 0L,
      accrued = matrix(c(time, time, 0), 1L, 3L)), B = 3L)
  expect_identical(draws[1L, , "terminal"], rep(tau, 3L))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})
