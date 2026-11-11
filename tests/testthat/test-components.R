component_complete_data <- function(prefix, n, calibration = FALSE) {
  i <- seq_len(n)
  first <- if (calibration) 0.5 + i %% 5 else 1 + i %% 4
  second <- first + 1 + (if (calibration) i %% 3 else (i * 3) %% 5)
  data.frame(id = rep(paste0(prefix, i), each = 2L),
             time = as.vector(rbind(first, second)),
             status = rep(c(1, 2), n))
}

component_training_fit <- function() {
  life_train(LifeHist(time, status) ~ 1,
              data = component_complete_data("train", 24),
              absorbing = 2, targets = list(event_free = 0, alive = c(0, 1)),
              tau = 12)
}

test_that("component allocation calibrates the original marginal scores without refitting", {
  fit <- component_training_fit()
  calibration <- component_complete_data("cal", 12, calibration = TRUE)
  truth <- matrix(calibration$time, nrow = 12, byrow = TRUE,
                   dimnames = list(NULL, c("event_free", "alive")))
  baseline <- data.frame(row.names = seq_len(12))
  moments <- predict(fit, baseline)
  scores <- abs(truth - moments$mean) / moments$scale
  q <- sort(unique(c(0, as.vector(scores), 100)))
  before <- serialize(fit, NULL)
  local_mocked_bindings(
    fit_life_censor = function(...) stop("Unexpected censoring refit"),
    fit_life_outcome = function(...) stop("Unexpected outcome refit"),
    fit_preprocessor = function(...) stop("Unexpected preprocessing refit"),
    .package = "lifeconformal"
  )
  allocations <- c(alive = 0.15, event_free = 0.20)
  for (method in c("ipw", "aipw")) {
    cal <- component_calibration(fit, allocations, alpha = 0.4,
      calibration_args = list(data = calibration, method = method, q = q,
                               control = list(keep_contributions = TRUE)))
    expect_identical(cal$fit, fit)
    expect_identical(cal$alpha_components, c(event_free = 0.20, alive = 0.15))
    expect_equal(cal$allocated_alpha, 0.35)
    expect_equal(cal$alpha, 0.4)
    expect_identical(names(cal$threshold), c("event_free", "alive"))
    for (j in seq_len(2)) {
      rank <- ceiling(13 * (1 - cal$alpha_components[j]))
      expect_equal(unname(cal$threshold[j]), unname(sort(scores[, j])[rank]))
      expect_equal(unname(cal$rank[j]), unname(rank))
      expect_equal(unname(cal$contributions[[j]]), unname(1 * outer(scores[, j], q, "<=")))
      curve <- cal$curve[cal$curve$target == colnames(scores)[j], ]
      expect_equal(curve$estimate, colMeans(1 * outer(scores[, j], q, "<=")))
    }
    expect_false(cal$full_support)
    expect_false(cal$any_component_full_support)
    expect_identical(cal$component_full_support, c(event_free = FALSE, alive = FALSE))
    prediction <- predict(cal, data.frame(id = "new"))
    expected_lower <- pmax(0, moments$mean[1, ] - cal$threshold * moments$scale[1, ])
    expected_upper <- pmin(12, moments$mean[1, ] + cal$threshold * moments$scale[1, ])
    expect_equal(prediction$bounds$lower, unname(expected_lower))
    expect_equal(prediction$bounds$upper, unname(expected_upper))
    expect_identical(prediction$region$support$joint_support$A, fit$problem$A)
    expect_identical(serialize(fit, NULL), before)
  }
})

test_that("partial and complete marginal fallbacks retain distinct meanings", {
  fit <- component_training_fit()
  calibration <- component_complete_data("cal", 12, calibration = TRUE)
  partial <- component_calibration(fit, c(event_free = 0.25, alive = 0.01), 0.3,
    list(data = calibration, method = "ipw", q = c(0, 1, 2, 100)))
  expect_identical(partial$component_full_support, c(event_free = FALSE, alive = TRUE))
  expect_false(partial$full_support)
  expect_true(partial$any_component_full_support)
  expect_true(is.finite(partial$threshold["event_free"]))
  expect_equal(unname(partial$threshold["alive"]), Inf)
  pred <- predict(partial, data.frame(id = "new"))
  expect_equal(pred$bounds$lower[pred$bounds$target == "alive"], 0)
  expect_equal(pred$bounds$upper[pred$bounds$target == "alive"], 12)
  expect_identical(partial$component_reason["alive"], c(alive = "rank_exceeds_sample"))
  expect_identical(unique(partial$curve$target), "event_free")

  full <- component_calibration(fit, c(event_free = 0.01, alive = 0.01), 0.1,
    list(data = calibration, method = "aipw"))
  expect_true(full$full_support)
  expect_true(all(full$component_full_support))
  expect_equal(nrow(full$curve), 0)
  expect_identical(names(full$curve),
                   c("target", "q", "estimate", "numerical_lower", "numerical_upper"))
  full_pred <- predict(full, data.frame(id = "new"))
  expect_equal(full_pred$bounds$lower, c(0, 0))
  expect_equal(full_pred$bounds$upper, c(12, 12))
})

test_that("allocations must cover every target and obey the joint budget", {
  targets <- c("event_free", "alive")
  check <- function(x, alpha = 0.3) validate_component_allocation(x, alpha, targets)
  expect_identical(check(c(alive = 0.1, event_free = 0.2)),
                   c(event_free = 0.2, alive = 0.1))
  expect_error(check(c(0.1, 0.2)), "target names")
  expect_error(check(c(event_free = 0.1)), "every target")
  expect_error(check(c(event_free = 0.1, unknown = 0.1)), "every target")
  expect_error(check(c(alive = 0.1, alive = 0.1)), "unique")
  expect_error(check(c(event_free = 0, alive = 0.1)), "strictly between")
  expect_error(check(c(event_free = -0.1, alive = 0.1)), "strictly between")
  expect_error(check(c(event_free = NA_real_, alive = 0.1)), "strictly between")
  expect_error(check(c(event_free = Inf, alive = 0.1)), "strictly between")
  expect_error(check(c(event_free = 0.2, alive = 0.2)), "sum.*not exceed")
  expect_error(check(c(event_free = 0.2, alive = 0.1 + 1e-10)), "sum.*not exceed")
  expect_error(check(c(event_free = 0.1, alive = 0.1), alpha = 0), "alpha")
  fit <- component_training_fit()
  expect_error(component_calibration(fit, c(event_free = 0.1, alive = 0.1), 0.3,
                                      list(alpha_components = NULL)), "must not replace")
  expect_error(component_calibration(fit, c(event_free = 0.1, alive = 0.1), 0.3,
                                      list("unnamed")), "unique, nonempty names")
})

test_that("temporary component problems preserve fitted model history features", {
  train <- data.frame(id = paste0("train", 1:12), time = 3, status = 0)
  calibration <- data.frame(id = paste0("cal", 1:12), time = 2, status = 0)
  fit <- life_train(LifeHist(time, status) ~ 1, data = train,
                    absorbing = 1, targets = list(alive = 0, dead = 1), tau = 10,
                    outcome_model = "transition_neural", seed = 7,
                    control = list(draws = 16L))
  scalar <- subset_component_problem(fit, 2L)
  expect_identical(scalar$outcome, fit$outcome)
  expect_identical(scalar$augmentation, fit$augmentation)
  expect_identical(colnames(scalar$problem$rewards), "dead")
  expect_identical(colnames(scalar$outcome$problem$rewards), c("alive", "dead"))
  expect_identical(scalar$control, fit$control)
  expect_identical(scalar$training_ids, fit$training_ids)
  baseline <- data.frame(row.names = 1:2)
  expect_equal(predict_life_moments(scalar, baseline)$mean,
               predict(fit, baseline)$mean[, "dead", drop = FALSE])
  cal <- component_calibration(fit, c(alive = 0.15, dead = 0.15), 0.3,
    list(data = calibration, method = "aipw", q = c(0, 100),
           control = list(keep_contributions = TRUE)))
  expect_equal(unname(cal$threshold), c(0, 0))
  expect_equal(cal$curve$estimate, rep(1, 4))
  expect_equal(cal$diagnostics$components$alive$numerical$queried_snapshots, 12)
  expect_equal(cal$diagnostics$components$dead$numerical$queried_snapshots, 12)
  expect_identical(cal$fit$outcome$problem$rewards, fit$problem$rewards)
})

test_that("public component allocation matches the internal fixed-score orchestration", {
  fit <- component_training_fit()
  data <- component_complete_data("cal", 12, calibration = TRUE)
  allocation <- c(alive = 0.15, event_free = 0.20)
  public <- life_calibrate(fit, data = data, alpha = 0.4,
                           alpha_components = allocation, method = "ipw",
                           q = c(0, 0.5, 1, 2, 4))
  internal <- component_calibration(fit, allocation, 0.4,
    list(data = data, method = "ipw", q = c(0, 0.5, 1, 2, 4)))
  expect_identical(public$threshold, internal$threshold)
  expect_identical(public$curve, internal$curve)
  expect_identical(public$component_full_support, internal$component_full_support)
  newdata <- data.frame(id = "new")
  expect_identical(predict(public, newdata)$bounds, predict(internal, newdata)$bounds)
})

component_censored_data <- function(prefix, n) {
  records <- lapply(seq_len(n), function(i) {
    first <- 0.5 + (i %% 7L) / 4 + 0.01 * sin(i * .31)
    second <- first + 1 + (i %% 5L) / 3 + 0.01 * cos(i * .73)
    pattern <- as.character(i %% 4L)
    times <- switch(pattern, "0" = first / 2, "1" = c(first, second),
                    "2" = second, "3" = c(first, second))
    status <- switch(pattern, "0" = 0L, "1" = c(1L, 0L),
                     "2" = 2L, "3" = c(1L, 2L))
    data.frame(id = paste0(prefix, i), time = times, status = status,
               x = sin(i * .91))
  })
  do.call(rbind, records)
}

test_that("shared censored preparation equals separate frozen marginal calibrations", {
  fit <- life_train(LifeHist(time, status) ~ x,
    data = component_censored_data("train", 48L), absorbing = 2,
    targets = list(event_free = 0, alive = c(0, 1)), tau = 6,
    score_weights = c(alive = 1, event_free = .5))
  data <- component_censored_data("cal", 24L)
  budgets <- c(event_free = .12, alive = .15)
  before <- serialize(fit, NULL)
  q <- c(0, .5, 1, 2, 5, 100)
  for (method in c("ipw", "aipw")) {
    shared <- life_calibrate(fit, data, alpha = .3,
      alpha_components = budgets, method = method, q = q,
      control = list(keep_contributions = TRUE))
    for (j in seq_along(budgets)) {
      separate <- life_calibrate(subset_component_problem(fit, j), data,
        alpha = unname(budgets[j]), method = method, q = q,
        control = list(keep_contributions = TRUE))
      expect_identical(unname(shared$threshold[j]), separate$threshold)
      expect_identical(shared$contributions[[j]], separate$contributions)
      curve <- shared$curve[shared$curve$target == names(budgets)[j], -1L]
      rownames(curve) <- NULL
      expect_identical(curve, separate$curve)
      expect_identical(shared$diagnostics$components[[j]], separate$diagnostics)
    }
    # Reusing an earlier calibration context across calls would miss this
    # change in observed follow-up and produce the wrong contribution matrix.
    shorter <- data
    shorter$time <- shorter$time * .7
    changed <- life_calibrate(fit, shorter, alpha = .3,
      alpha_components = budgets, method = method, q = q,
      control = list(keep_contributions = TRUE))
    expect_false(identical(shared$contributions, changed$contributions))
  }
  expect_identical(serialize(fit, NULL), before)
})

test_that("component parsing, baseline predictions and endpoint queries occur once", {
  fit <- component_training_fit()
  data <- component_censored_data("cal", 24L)
  normalize_original <- normalize_histories
  moments_original <- predict_life_moments
  queries_original <- prepare_calibration_queries
  counts <- c(parse = 0L, moments = 0L, queries = 0L)
  local_mocked_bindings(
    normalize_histories = function(...) {
      counts["parse"] <<- counts["parse"] + 1L
      normalize_original(...)
    },
    predict_life_moments = function(...) {
      counts["moments"] <<- counts["moments"] + 1L
      moments_original(...)
    },
    prepare_calibration_queries = function(...) {
      counts["queries"] <<- counts["queries"] + 1L
      queries_original(...)
    }, .package = "lifeconformal")
  for (method in c("ipw", "aipw")) {
    counts[] <- 0L
    life_calibrate(fit, data, alpha = .3,
      alpha_components = c(event_free = .15, alive = .15),
      method = method, q = c(0, 1, 4, 100))
    expect_identical(counts, c(parse = 1L, moments = 1L,
                               queries = as.integer(method == "aipw")))
  }
})

test_that("rank-infeasible components skip all nuisance prediction", {
  fit <- component_training_fit()
  data <- component_censored_data("cal", 4L)
  local_mocked_bindings(
    predict_life_moments = function(...) stop("Unneeded baseline prediction"),
    prepare_calibration_queries = function(...) stop("Unneeded endpoint plan"),
    predict_life_censor = function(...) stop("Unneeded censoring prediction"),
    .package = "lifeconformal")
  cal <- life_calibrate(fit, data, alpha = .1,
    alpha_components = c(event_free = .01, alive = .01))
  expect_true(cal$full_support)
  expect_equal(nrow(cal$curve), 0)
  expect_error(life_calibrate(fit, data, alpha = .4,
    alpha_components = c(event_free = .2, alive = .2), q = c(1, 0)),
    "strictly increasing")
})

test_that("forest censoring curves are shared across all component scores", {
  skip_if_not_installed("ranger")
  fit <- life_train(LifeHist(time, status) ~ x,
    data = component_censored_data("train", 48L), absorbing = 2,
    targets = list(event_free = 0, alive = c(0, 1)), tau = 6,
    censor_model = "survival_forest",
    censor_options = list(num_trees = 5L, min_node_size = 3L), seed = 8L)
  data <- component_censored_data("cal", 24L)
  original <- censor_forest_curves
  calls <- 0L
  local_mocked_bindings(censor_forest_curves = function(...) {
    calls <<- calls + 1L
    original(...)
  }, .package = "lifeconformal")
  for (method in c("ipw", "aipw")) {
    calls <- 0L
    shared <- life_calibrate(fit, data, alpha = .3,
      alpha_components = c(event_free = .15, alive = .15),
      method = method, q = c(0, .5, 1, 2, 5, 100),
      control = list(keep_contributions = TRUE))
    expect_identical(calls, 1L)
    for (j in 1:2) {
      separate <- life_calibrate(subset_component_problem(fit, j), data,
        alpha = .15, method = method, q = c(0, .5, 1, 2, 5, 100),
        control = list(keep_contributions = TRUE))
      expect_identical(shared$contributions[[j]], separate$contributions)
      expect_identical(unname(shared$threshold[j]), separate$threshold)
    }
  }
})
