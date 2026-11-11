component_sharing_data <- function(prefix, n) {
  do.call(rbind, lapply(seq_len(n), function(i) {
    first <- .4 + (i %% 5L) / 7 + .001 * sin(i)
    last <- first + .6 + (i %% 7L) / 9
    if (i %% 3L == 0L) {
      times <- first / 2; status <- 0L
    } else {
      times <- c(first, last)
      status <- c(1L, if (i %% 3L == 1L) 0L else 2L)
    }
    data.frame(id = paste0(prefix, i), time = times, status = status, x = sin(i * .7))
  }))
}

component_sharing_fixture <- function(augmentation = "transition_forest") {
  if (identical(augmentation, "transition_forest")) skip_if_not_installed("rpart")
  options <- switch(augmentation,
    transition_forest = list(num_trees = 3L, step_size = .13),
    transition_neural = list(hidden = 2L, decay = .1, maxit = 1000L, step_size = .13),
    list())
  fit <- life_train(LifeHist(time, status) ~ x,
    data = component_sharing_data("train", 24L), absorbing = 2,
    targets = list(event_free = 0, alive = c(0, 1), total = c(0, 1, 2)), tau = 3,
    allowed_transitions = data.frame(from = c(0, 1), to = c(1, 2)),
    score_weights = c(event_free = .7, alive = 1.3, total = 3),
    augmentation_model = if (identical(augmentation, "markov_cox")) NULL else augmentation,
    augmentation_options = options, control = list(draws = 16L, scale_floor = .05), seed = 12L)
  list(fit = fit, data = component_sharing_data("cal", 12L))
}

component_sharing_without_call <- function(value) {
  value$call <- NULL
  value
}

component_sharing_compare <- function(fixture, keep = TRUE, q = NULL,
    budgets = c(event_free = .15, alive = .15, total = .15), method = "aipw") {
  before <- serialize(fixture$fit, NULL, version = 2)
  withr::local_seed(817)
  rng <- .Random.seed
  run <- function() life_calibrate(fixture$fit, data = fixture$data,
    alpha = .5, alpha_components = budgets, method = method, q = q,
    control = list(keep_contributions = keep))
  candidate <- run()
  expect_identical(.Random.seed, rng)
  reference <- with_mocked_bindings(run(),
    calibration_sum_precision = function(...) NULL, .package = "lifeconformal")
  expect_identical(component_sharing_without_call(candidate), component_sharing_without_call(reference))
  expect_identical(serialize(fixture$fit, NULL, version = 2), before)
  expect_identical(.Random.seed, rng)
  candidate
}

test_that("shared forest and neural components exactly preserve public results", {
  skip_if(is.null(calibration_sum_precision()), "R and native accumulation formats differ")
  for (augmentation in c("transition_forest", "transition_neural")) {
    fixture <- component_sharing_fixture(augmentation)
    # Replaying a version-2 saved fit must not introduce persistent sample state.
    fixture$fit <- unserialize(serialize(fixture$fit, NULL, version = 2))
    original <- evaluate_flexible_component_curves
    calls <- 0L
    results <- with_mocked_bindings({
      full <- component_sharing_compare(fixture, keep = TRUE)
      partial <- component_sharing_compare(fixture, keep = FALSE,
        budgets = c(event_free = .18, alive = .001, total = .12))
      partial_kept <- component_sharing_compare(fixture, keep = TRUE,
        budgets = c(event_free = .18, alive = .001, total = .12))
      list(full = full, partial = partial, partial_kept = partial_kept)
    }, evaluate_flexible_component_curves = function(...) {
      calls <<- calls + 1L
      original(...)
    }, .package = "lifeconformal")
    expect_identical(calls, 3L)
    full <- results$full; partial <- results$partial
    expect_false(identical(full$curve$q[full$curve$target == "event_free"],
      full$curve$q[full$curve$target == "alive"]))
    expect_lt(sum(full$curve$target == "total"), sum(full$curve$target == "event_free"))
    expect_true(all(full$contributions$total == 1))
    for (target in names(full$contributions))
      expect_identical(full$diagnostics$components[[target]]$numerical$curve_estimate_sum,
        colSums(full$contributions[[target]]))
    expect_null(partial$contributions)
    expect_null(results$partial_kept$contributions$alive)
    expect_identical(results$partial_kept$curve, partial$curve)
    expect_identical(partial$component_reason[["alive"]], "rank_exceeds_sample")
    expect_identical(unique(partial$curve$target), c("event_free", "total"))
  }
})

test_that("bounded flexible banks preserve results and release state on failures", {
  skip_if(is.null(calibration_sum_precision()), "R and native accumulation formats differ")
  fixture <- component_sharing_fixture()
  native_count <- function() .Call("life_calibration_sum_info", PACKAGE = "lifeconformal")[[7L]]
  baseline <- native_count()
  evaluator <- evaluate_flexible_component_curves
  simulator <- simulate_flexible_targets
  sizes <- integer(); inside_shared <- FALSE
  bounded <- with_mocked_bindings(component_sharing_compare(fixture,
    q = c(0, .5, 1, 100), keep = FALSE),
    evaluate_flexible_component_curves = function(...) {
      inside_shared <<- TRUE
      on.exit(inside_shared <<- FALSE)
      # With 16 draws, 3 targets and 4 thresholds, 1024 bytes caps each bank
      # at two queries, independently of the total number of subject queries.
      evaluator(..., query_batch = 50L, bank_budget = 1024)
    }, simulate_flexible_targets = function(adapter, x, snapshots, ...) {
      if (inside_shared) sizes <<- c(sizes, nrow(x))
      simulator(adapter, x, snapshots, ...)
    }, .package = "lifeconformal")
  expect_gt(length(sizes), 0L)
  expect_true(all(sizes <= 2L))
  expect_true(any(sizes == 2L))
  expect_null(bounded$contributions)
  expect_identical(native_count(), baseline)
  before <- serialize(fixture$fit, NULL, version = 2)
  withr::local_seed(901)
  rng <- .Random.seed
  expect_error(with_mocked_bindings(
    life_calibrate(fixture$fit, fixture$data, alpha = .5,
      alpha_components = c(event_free = .15, alive = .15, total = .15), q = c(0, .5, 1)),
    simulate_flexible_targets = function(...) stop("injected continuation failure"),
    .package = "lifeconformal"), "injected continuation failure")
  expect_identical(native_count(), baseline)
  expect_identical(serialize(fixture$fit, NULL, version = 2), before)
  expect_identical(.Random.seed, rng)
  # A fresh successful call after the error cannot inherit partial sums.
  component_sharing_compare(fixture, q = c(0, .5, 1, 100), keep = TRUE)
  expect_identical(native_count(), baseline)
})

test_that("fully known component curves need no flexible continuation", {
  skip_if(is.null(calibration_sum_precision()), "R and native accumulation formats differ")
  fixture <- component_sharing_fixture()
  value <- with_mocked_bindings(component_sharing_compare(fixture, q = 10000),
    simulate_flexible_targets = function(...) stop("unexpected continuation"),
    .package = "lifeconformal")
  expect_true(all(value$curve$estimate == 1))
  for (diagnostics in value$diagnostics$components) {
    expect_identical(diagnostics$numerical$queried_snapshots, 0L)
    expect_false(diagnostics$numerical$conditional_coverage_evaluated)
  }
})

test_that("single active, IPW and Markov component paths bypass sharing", {
  fixture <- component_sharing_fixture()
  markov <- component_sharing_fixture("markov_cox")
  values <- with_mocked_bindings({
    one <- component_sharing_compare(fixture, q = c(0, 1, 100),
      budgets = c(event_free = .2, alive = .001, total = .001))
    ipw <- component_sharing_compare(fixture, q = c(0, 1, 100), method = "ipw")
    deterministic <- component_sharing_compare(markov, q = c(0, 1, 100))
    list(one = one, ipw = ipw, deterministic = deterministic)
  }, evaluate_flexible_component_curves = function(...) stop("unexpected shared path"),
    .package = "lifeconformal")
  expect_identical(unique(values$one$curve$target), "event_free")
  expect_true(all(values$one$component_reason[c("alive", "total")] == "rank_exceeds_sample"))
  expect_identical(values$ipw$method, "ipw")
  expect_identical(values$deterministic$fit$augmentation$type, "markov_cox")
})
