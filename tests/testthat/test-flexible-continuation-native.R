forest_block_fixture <- function() {
  compiled <- list(roots = c(1L, 4L, 7L, 10L),
    feature = c(1L, 0L, 0L, 2L, 0L, 0L, 1L, 0L, 0L, 0L),
    left = c(2L, 0L, 0L, 5L, 0L, 0L, 8L, 0L, 0L, 0L),
    right = c(3L, 0L, 0L, 6L, 0L, 0L, 9L, 0L, 0L, 0L),
    direction = c(-1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L),
    cut = c(0, 0, 0, .5, 0, 0, 0, 0, 0, 0),
    value = c(0, .1, .3, 0, .7, .9, 0, .2, .05, .6), n_features = 6L,
    tree_start = c(1L, 4L, 7L, 10L), tree_end = c(3L, 6L, 9L, 10L))
  compiled$lookup <- compile_forest_lookup(compiled)
  compiled$additive <- compile_forest_additive(compiled)
  structure(list(type = "transition_forest", fits = rep(list(list(constant = NULL,
    compiled = compiled)), 4L), edges = data.frame(from = c(1L, 1L, 2L, 2L), to = c(2L, 3L, 1L, 3L)),
    problem = list(tau = 4, rewards = matrix(c(1, 0, 0, 1, 1, 1), 3L, 2L,
      dimnames = list(NULL, c("before", "total"))), absorbing = 3L),
    options = list(step_size = .1, max_events = 100L), x_names = "x", n_predictors = 1L,
    seed = 37L), class = "life_flexible_outcome")
}

forest_block_reference <- function(adapter, x, snapshots, B) {
  testthat::with_mocked_bindings(simulate_flexible_targets(adapter, x, snapshots, B = B),
    flex_forest_block_precision = function(...) NULL, .package = "lifeconformal")
}

test_that("native forest blocks preserve seeded banks, targets and input state", {
  adapter <- forest_block_fixture()
  precision <- flex_forest_block_precision()
  skip_if(is.null(precision), "Runtime accumulator format uses the R fallback")
  expect_true(.Call(life_forest_block_supported, adapter, 16L * 31L, 31L, precision))
  x <- matrix(rep(c(-1L, 0L, 1L), 11L), 33L, dimnames = list(NULL, "x"))
  time <- rep(c(0, .213, .714), 11L)
  snapshots <- list(time = time, state = rep(c(1L, 2L, 1L), 11L),
    duration = rep(c(0, .11, .19), 11L), previous_transitions = rep(c(0L, 1L, 2L), 11L),
    accrued = cbind(time * .5, time))
  set.seed(179); rng <- .Random.seed; stored <- serialize(adapter, NULL)
  stored_v2 <- serialize(adapter, NULL, version = 2)
  expected <- forest_block_reference(adapter, x, snapshots, 31L)
  actual <- simulate_flexible_targets(adapter, x, snapshots, B = 31L)
  expect_identical(actual, expected)
  expect_identical(.Random.seed, rng)
  expect_identical(serialize(adapter, NULL), stored)
  expect_identical(serialize(adapter, NULL, version = 2), stored_v2)
  expect_identical(actual, simulate_flexible_targets(unserialize(stored), x, snapshots, B = 31L))
  reversed <- lapply(snapshots, function(z) if (is.matrix(z)) z[33:1, , drop = FALSE] else z[33:1])
  expect_identical(as.numeric(actual[33:1, , , drop = FALSE]), as.numeric(
    simulate_flexible_targets(adapter, x[33:1, , drop = FALSE], reversed, B = 31L)))
  expect_identical(as.numeric(actual[, 1:7, , drop = FALSE]), as.numeric(
    simulate_flexible_targets(adapter, x, snapshots, B = 7L)))
  expect_identical(as.numeric(actual[, , "total"]), rep(adapter$problem$tau, 33L * 31L))
  one <- lapply(snapshots, function(z) if (is.matrix(z)) z[1L, , drop = FALSE] else z[1L])
  expect_identical(simulate_flexible_targets(adapter, x[1L, , drop = FALSE], one, B = 1L),
    forest_block_reference(adapter, x[1L, , drop = FALSE], one, 1L))
})

forest_component_fixture <- function() {
  records <- function(prefix, n) do.call(rbind, lapply(seq_len(n), function(i) {
    first <- .4 + (i %% 5L) / 7 + .001 * sin(i)
    last <- first + .6 + (i %% 7L) / 9
    if (i %% 3L == 0L) data.frame(id = paste0(prefix, i), time = first / 2, status = 0L, x = sin(i * .7))
    else data.frame(id = paste0(prefix, i), time = c(first, last),
      status = c(1L, if (i %% 3L == 1L) 0L else 2L), x = sin(i * .7))
  }))
  fit <- life_train(LifeHist(time, status) ~ x, records("train", 24L), absorbing = 2,
    targets = list(before = 0, alive = c(0, 1), total = c(0, 1, 2)), tau = 3,
    allowed_transitions = data.frame(from = c(0, 1), to = c(1, 2)),
    augmentation_model = "transition_forest",
    augmentation_options = list(num_trees = 3L, step_size = .13),
    control = list(draws = 16L, scale_floor = .05), seed = 12L)
  list(fit = fit, data = records("cal", 12L))
}

test_that("bounded component sharing uses full augmentation histories in native blocks", {
  skip_if_not_installed("rpart")
  skip_if(is.null(flex_forest_block_precision()) || is.null(calibration_sum_precision()),
    "Runtime accumulator format uses the R fallback")
  fixture <- forest_component_fixture()
  before <- serialize(fixture$fit, NULL, version = 2)
  set.seed(817); rng <- .Random.seed
  original_state <- flex_forest_block_state; original_simulate <- simulate_flexible_targets
  original_curves <- evaluate_flexible_component_curves
  calls <- 0L; targets <- integer(); queries <- integer()
  run <- function() life_calibrate(fixture$fit, fixture$data, alpha = .5,
    alpha_components = c(before = .2, alive = .2, total = .001), q = c(0, .5, 1, 100),
    control = list(keep_contributions = TRUE))
  candidate <- with_mocked_bindings(run(),
    flex_forest_block_state = function(...) { calls <<- calls + 1L; original_state(...) },
    simulate_flexible_targets = function(adapter, x, snapshots, ...) {
      targets <<- c(targets, ncol(adapter$problem$rewards)); queries <<- c(queries, nrow(x))
      original_simulate(adapter, x, snapshots, ...)
    }, evaluate_flexible_component_curves = function(...)
      original_curves(..., query_batch = 50L, bank_budget = 1024), .package = "lifeconformal")
  reference <- with_mocked_bindings(run(),
    flex_forest_block_precision = function(...) NULL,
    evaluate_flexible_component_curves = function(...)
      original_curves(..., query_batch = 50L, bank_budget = 1024), .package = "lifeconformal")
  candidate$call <- reference$call <- NULL
  expect_identical(candidate, reference)
  expect_gt(calls, 0L)
  expect_true(all(targets == 3L))
  expect_true(all(queries <= 2L))
  expect_identical(serialize(fixture$fit, NULL, version = 2), before)
  expect_identical(.Random.seed, rng)
  expect_identical(fixture$fit$outcome$type, "markov_cox")
  expect_identical(fixture$fit$augmentation$type, "transition_forest")
})

test_that("forest dispatch preserves resource and arithmetic fallbacks", {
  adapter <- forest_block_fixture(); precision <- flex_forest_block_precision()
  skip_if(is.null(precision), "Runtime accumulator format uses the R fallback")
  expect_false(.Call(life_forest_block_supported, adapter, .Machine$integer.max, 1L, precision))
  malformed <- adapter; malformed$fits[[1L]]$compiled$roots <- 999999L
  # An ineligible workspace must fall back before inspecting unused metadata.
  expect_false(.Call(life_forest_block_supported, malformed, .Machine$integer.max, 1L, precision))
  old <- adapter; old$fits[[1L]]$compiled$additive <- NULL
  expect_false(.Call(life_forest_block_supported, old, 7L, 7L, precision))
  x <- matrix(0, 1L, 1L, dimnames = list(NULL, "x"))
  snapshot <- list(time = .2, state = 1L, duration = .2, previous_transitions = 0L,
    accrued = matrix(c(.1, .2), 1L, 2L))
  expect_identical(simulate_flexible_targets(old, x, snapshot, B = 7L),
    forest_block_reference(old, x, snapshot, 7L))
  expect_null(flex_forest_block_precision(c(53L, 64L, 1024L, 10000L),
    machine = list(double.digits = 53L, double.max.exp = 1024L,
      longdouble.digits = 64L, longdouble.max.exp = 16384L), wide = TRUE))
})

test_that("native forest errors preserve models and caller RNG", {
  adapter <- forest_block_fixture(); adapter$options$max_events <- 0L
  for (e in seq_along(adapter$fits)) adapter$fits[[e]]$compiled$additive$value <-
    adapter$fits[[e]]$compiled$additive$value * 1e8
  x <- matrix(0, 1L, 1L, dimnames = list(NULL, "x"))
  snapshot <- list(time = 0, state = 1L, duration = 0, previous_transitions = 0L,
    accrued = matrix(0, 1L, 2L))
  set.seed(171); rng <- .Random.seed; stored <- serialize(adapter, NULL)
  expect_error(simulate_flexible_targets(adapter, x, snapshot, B = 7L), "max_events")
  expect_error(forest_block_reference(adapter, x, snapshot, 7L), "max_events")
  expect_identical(.Random.seed, rng)
  expect_identical(serialize(adapter, NULL), stored)
})

test_that("native blocks finish terminal arrivals before another bank", {
  adapter <- forest_block_fixture(); adapter$fits <- list(list(constant = 1))
  adapter$edges <- data.frame(from = 1L, to = 3L)
  adapter$options$step_size <- 1; adapter$problem$tau <- 100
  precision <- flex_forest_block_precision()
  skip_if(is.null(precision), "Runtime accumulator format uses the R fallback")
  path <- flex_forest_block_state(31, 1L, 0, 31, 32, 0,
    matrix(c(31, 31), 1L, 2L), 0L, 31L)
  bank <- matrix(0, 1L, 64L); bank[1L, 63L] <- exp(-.5); bank[1L, 64L] <- .5
  controls <- list(draws = 1L, precision = precision, steps = 1L, max_iter = 200L, trace = 1L)
  result <- .Call(life_forest_block_advance, adapter, matrix(0L, 1L), 1L, path, bank, controls)
  expect_identical(result$path$iteration, 32L)
  expect_identical(result$path$state, 3L)
  expect_identical(result$path$tt, 31.5)
  # A subsequent block sees terminal completion without advancing iteration.
  controls$steps <- 32L
  complete <- .Call(life_forest_block_advance, adapter, matrix(0L, 1L), 1L,
    result$path, bank, controls)
  expect_identical(complete$path$iteration, 32L)
  expect_identical(complete$path$tt, 100)
  expect_length(complete$trace, 0L)
  expect_identical(as.numeric(complete$path$accrued), c(31.5, 100))
  bad <- path; bad$count <- .Machine$integer.max
  expect_error(.Call(life_forest_block_advance, adapter, matrix(0L, 1L), 1L, bad,
    bank, modifyList(controls, list(steps = 1L))), "path state")
})
