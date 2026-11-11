flex_options_workflow_data <- function(prefix = "train", n = 18L) {
  i <- seq_len(n)
  data.frame(id = paste0(prefix, i), time = .5 + ((i * 7) %% 13) / 5,
    status = ifelse(i %% 3 == 0, 0L, 2L), z = ((i * 5) %% 11) / 10,
    w = ((i * 3) %% 7) / 7)
}

flex_options_workflow_fit <- function(data = flex_options_workflow_data(),
    outcome_model = "transition_neural", outcome_options = NULL,
    censor_model = "survival_forest", censor_options = list(num_trees = 3L),
    augmentation_model = NULL, augmentation_options = list(), seed = 0) {
  if (is.null(outcome_options)) outcome_options <- if (outcome_model == "transition_forest")
    list(num_trees = 3L, step_size = 1, minsplit = 4L, minbucket = 2L) else
    list(hidden = 1L, decay = .2, maxit = 1000L, step_size = 1)
  life_train(LifeHist(time, status) ~ z, data = data, absorbing = 2,
    targets = list(alive = 0), tau = 2.8, censor_formula = ~ w,
    censor_model = censor_model, censor_options = censor_options,
    outcome_model = outcome_model, outcome_options = outcome_options,
    augmentation_model = augmentation_model, augmentation_options = augmentation_options,
    control = list(draws = 16L), seed = seed)
}

test_that("flexible option errors identify unsupported and malformed settings", {
  for (model in c("transition_forest", "transition_neural")) {
    expect_error(flex_outcome_options(list(typo = 1), list(tau = 2), model), "Unknown.*typo")
    expect_error(flex_outcome_options(list(step_size = 1, step_size = 2),
      list(tau = 2), model), "unique, nonempty names")
    expect_error(flex_outcome_options(list(1), list(tau = 2), model), "unique, nonempty names")
    for (value in list(NULL, NA_real_, NaN, Inf, 1 + 1i, "1", c(1, 2), matrix(1))) {
      expect_error(flex_outcome_options(list(step_size = value), list(tau = 2), model),
        "step_size must be a finite positive number")
    }
  }
  for (name in c("hidden", "maxit", "decay", "reltol", "max_training_rows", "max_events")) {
    expect_error(flex_outcome_options(stats::setNames(list(NULL), name),
      list(tau = 2), "transition_neural"), name)
  }
  for (name in c("num_trees", "max_depth", "minsplit", "minbucket", "cp")) {
    expect_error(flex_outcome_options(stats::setNames(list(NULL), name),
      list(tau = 2), "transition_forest"), name)
  }
  expect_null(flex_outcome_options(list(mtry = NULL), list(tau = 2),
    "transition_forest")$mtry)
  expect_error(flex_outcome_options(list(num_trees = 1.5), list(tau = 2),
    "transition_forest"), "num_trees must be an integer")
  expect_error(flex_outcome_options(list(max_depth = 31), list(tau = 2),
    "transition_forest"), "max_depth cannot exceed 30")
})

test_that("censoring forest options are validated even without censoring events", {
  d <- flex_options_workflow_data()
  d$status[] <- 2L
  expect_error(flex_options_workflow_fit(data = d, censor_options = list(mtry = 2L)),
    "mtry exceeds the predictor count")
  for (value in list(NULL, NA_real_, 0, 1 + 1i, matrix(2))) {
    expect_error(flex_options_workflow_fit(data = d,
      censor_options = list(num_trees = value)), "num_trees must be a finite positive number")
  }
  expect_error(flex_options_workflow_fit(data = d,
    censor_options = list(mtry = 1L, mtry = 1L)), "unique, nonempty names")
  expect_error(flex_options_workflow_fit(data = d,
    censor_options = list(seed = 1L)), "Unknown.*seed")
  fitted <- flex_options_workflow_fit(data = d)
  expect_identical(fitted$censor$options$mtry, 1)
  expect_false(fitted$censor$diagnostics$backend_fitted)
  expect_equal(fitted$functions$censor_survival(data.frame(w = c(.1, .2)), c(1, 2)), c(1, 1))
})

test_that("flexible models keep predictor roles and summaries understandable", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("rpart")
  for (outcome in c("transition_forest", "transition_neural")) {
    augmentation <- if (outcome == "transition_forest") "transition_neural" else "transition_forest"
    augmentation_options <- if (augmentation == "transition_neural")
      list(hidden = 1L, decay = .3, maxit = 1000L, step_size = .8) else
      list(num_trees = 2L, step_size = .8, minsplit = 4L, minbucket = 2L)
    fit <- flex_options_workflow_fit(outcome_model = outcome,
      augmentation_model = augmentation, augmentation_options = augmentation_options)
    newdata <- data.frame(z = c(.1, .8))
    expect_true(all(is.finite(predict(fit, newdata)$mean)))
    coverage <- predict(fit, newdata, type = "coverage", q = c(0, 1, 2))
    expect_identical(coverage$numerical$approximation, "monte_carlo")
    expect_identical(coverage$numerical$deterministic_probability_bounds, FALSE)
    expect_identical(coverage$numerical$maximum_probability_bound_width, NA_real_)
    expect_true(all(is.finite(fit$functions$censor_survival(data.frame(w = c(.2, .6)), 1))))
    expect_error(predict(fit, data.frame(w = 1)), "Missing baseline column.*z")
    expect_error(predict(fit, data.frame(z = NA_real_)), "finite, nonmissing")
    expect_error(predict(fit, data.frame(z = 1 + 1i)), "must be numeric")
    expect_error(fit$functions$censor_survival(newdata, 1), "Missing baseline column.*w")
    details <- summary(fit)
    expect_identical(details$model_options$outcome, fit$outcome$options)
    expect_identical(details$model_options$censoring$mtry, 1)
    expect_identical(details$model_options$augmentation, fit$augmentation$options)
    expect_identical(details$model_diagnostics$outcome$resolved_options, fit$outcome$options)
    expect_identical(details$model_diagnostics$augmentation$seed, fit$augmentation$seed)
    expect_identical(details$model_diagnostics$censoring$backend_fitted, TRUE)
    caldata <- flex_options_workflow_data("cal", 12L)
    caldata$time <- caldata$time * .93
    for (method in c("ipw", "aipw")) {
      cal <- life_calibrate(fit, data = caldata, method = method,
        alpha = .25, q = c(0, 1, 2, 4), control = list(draws = 20L))
      expect_identical(cal$fit$control$draws, 16L)
      expect_true(all(is.finite(predict(cal, newdata)$bounds$lower)))
    }
  }
})

test_that("fixed public seeds separate censoring and outcome random streams", {
  skip_if_not_installed("ranger")
  d <- flex_options_workflow_data()
  set.seed(81); before <- .Random.seed
  cox <- flex_options_workflow_fit(censor_model = "cox", censor_options = list())
  forest <- flex_options_workflow_fit()
  expect_identical(cox$outcome$seed, forest$outcome$seed)
  expect_identical(predict(cox, d["z"])$mean, predict(forest, d["z"])$mean)
  expect_identical(.Random.seed, before)
  extra <- flex_options_workflow_fit(augmentation_model = "transition_neural",
    augmentation_options = list(hidden = 1L, decay = .3, maxit = 1000L, step_size = 1))
  expect_identical(extra$outcome$seed, forest$outcome$seed)
  expect_identical(predict(extra, d["z"])$mean, predict(forest, d["z"])$mean)
  expect_identical(.Random.seed, before)
  expect_error(flex_options_workflow_fit(outcome_options = list(seed = 25L)),
    "Unknown.*seed")
  for (seed in list(NA_real_, NaN, Inf, 1 + 1i, "1", c(1, 2), matrix(1))) {
    expect_error(flex_with_seed(seed, 1), "seed must be an integer")
  }
})

test_that("an unseeded deterministic fit leaves the caller stream unchanged", {
  d <- flex_options_workflow_data()
  set.seed(70); before <- .Random.seed
  fit <- life_train(LifeHist(time, status) ~ 1, data = d, absorbing = 2,
    targets = list(alive = 0), tau = 2.8, seed = NULL)
  expect_identical(.Random.seed, before)
  expect_null(fit$fitting_seeds)
})

test_that("neural convergence errors report the original state labels", {
  expect_error(flex_options_workflow_fit(censor_model = "cox", censor_options = list(),
    outcome_options = list(hidden = 1L, maxit = 1L, step_size = 1)),
    "edge 0 -> 2; inspect model/data scaling")
})
