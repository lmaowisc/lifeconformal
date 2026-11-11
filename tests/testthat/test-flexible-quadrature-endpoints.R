flexible_endpoint_fixture <- function(start, stop, event = TRUE) {
  history <- normalize_histories(resolve_input(id = c(1L, 1L),
    response = LifeHist(c(start, stop), c(1L, if (event) 2L else 0L))),
    absorbing = 2L)
  list(history = history, x = matrix(0, 1L, 0L),
    problem = compile_problem(history, list(alive = c(0L, 1L)), tau = 2))
}

test_that("rounded terminal cell counts retain all exposure and transitions", {
  start <- .4 + .001 * sin(35)
  stop <- start + .6
  expect_identical(ceiling((stop - start) / .1), 7)
  expect_identical(start + 6 * .1, stop)
  for (event in c(TRUE, FALSE)) {
    z <- flexible_endpoint_fixture(start, stop, event)
    d <- make_flexible_transition_data(z$history, z$x, z$problem, .1,
      max_rows = 10L)
    row <- which(d$from == match(1L, z$problem$states))
    expect_length(row, 6L)
    boundaries <- c(start + (0:5) * .1, stop)
    expect_identical(d$exposure[row], diff(boundaries))
    expect_equal(sum(d$exposure[row]), stop - start, tolerance = 1e-15)
    expect_equal(sum(d$exposure), stop, tolerance = 1e-15)
    expect_true(all(d$exposure > 0))
    expect_identical(d$features[row, 1L], start + (0:5) * .1 + .05)
    expect_identical(sum(d$to > 0L), if (event) 2L else 1L)
    expect_identical(tail(d$to, 1L), if (event) 3L else 0L)
    expect_error(make_flexible_transition_data(z$history, z$x, z$problem,
      .1, max_rows = 9L), "budget")
  }
})

test_that("positive one-ulp final cells are not removed by a tolerance", {
  start <- .4 + .001 * sin(35)
  boundary <- start + 6 * .1
  stop <- boundary + 2^-53
  expect_gt(stop, boundary)
  z <- flexible_endpoint_fixture(start, stop)
  d <- make_flexible_transition_data(z$history, z$x, z$problem, .1)
  row <- which(d$from == match(1L, z$problem$states))
  expect_length(row, 7L)
  expect_identical(tail(d$exposure, 1L), 2^-53)
  expect_identical(tail(d$to, 1L), 3L)
  expect_true(all(d$exposure > 0))
  expect_equal(sum(d$exposure[row]), stop - start, tolerance = 1e-15)
})

test_that("unresolvable interior cells still report a small step size", {
  # Direct private fixture avoids allocating the immense pre-entry sojourn.
  z <- flexible_endpoint_fixture(.4, 1)
  z$history$intervals <- z$history$intervals[2L, , drop = FALSE]
  z$history$intervals$start <- 2^52
  z$history$intervals$stop <- 2^52 + 1
  z$history$row_index <- list(1L)
  expect_error(make_flexible_transition_data(z$history, z$x, z$problem, .1),
    "step_size is too small for the numerical resolution")
})

test_that("ordinary shifted decimal intervals train a public transition forest", {
  skip_if_not_installed("rpart")
  training <- do.call(rbind, lapply(seq_len(48L), function(i) {
    a <- .4 + (i %% 5) / 7 + .001 * sin(i)
    b <- a + .6 + (i %% 7) / 9
    if (i %% 3 == 0) {
      data.frame(id = i, time = a / 2, status = 0L, x = sin(i * .7))
    } else {
      data.frame(id = i, time = c(a, b),
        status = c(1L, if (i %% 3 == 1) 0L else 2L), x = sin(i * .7))
    }
  }))
  set.seed(17)
  rng <- .Random.seed
  fit <- life_train(LifeHist(time, status) ~ x, data = training, absorbing = 2,
    targets = list(event_free = 0, alive = c(0, 1)), tau = 3,
    allowed_transitions = data.frame(from = c(0, 1), to = c(1, 2)),
    augmentation_model = "transition_forest",
    augmentation_options = list(num_trees = 3L, step_size = .1),
    control = list(draws = 8L), seed = 12)
  expect_s3_class(fit, "lifeconformal_fit")
  expect_identical(.Random.seed, rng)
  expect_identical(fit$augmentation$options$step_size, .1)
})
