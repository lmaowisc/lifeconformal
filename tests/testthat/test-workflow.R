workflow_survival_data <- function(prefix = "train", censored = FALSE) {
  if (!censored) {
    return(data.frame(id = paste0(prefix, seq_len(12)),
                      time = rep(c(1, 2, 4, 7, 9, 12), 2), status = 1))
  }
  i <- seq_len(48)
  data.frame(id = paste0(prefix, i), time = ((i * 5) %% 17 + 1) / 2,
             status = as.integer(!i %% 5 %in% c(0, 1)),
             age = 40 + (i * 7 %% 31), z = (i * 11 %% 19) / 10,
             treatment = factor(rep(c("A", "B"), 24)))
}

workflow_recurrent_data <- function(prefix = "train") {
  records <- vector("list", 12L)
  for (i in seq_len(12L)) {
    initial <- (i - 1L) %% 3L
    if (initial == 0L) {
      time <- c(2, 4, 7)
      status <- c(1, 7, 2)
    } else if (initial == 1L) {
      time <- c(1, 3, 6)
      status <- c(7, 1, 2)
    } else {
      time <- 11
      status <- 0
    }
    records[[i]] <- data.frame(id = paste0(prefix, i), time = time,
                               status = status, baseline_state = initial)
  }
  do.call(rbind, records)
}

test_that("complete-data survival recovers the exact split-conformal rank for both methods", {
  training <- workflow_survival_data()
  calibration <- data.frame(id = paste0("cal", seq_len(12)),
                            time = c(1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15),
                            status = 1)
  fit <- life_train(LifeHist(time, status) ~ 1, data = training,
                    absorbing = 1, targets = list(alive = 0), tau = 10)
  newdata <- data.frame(row.names = seq_len(12))
  moments <- predict(fit, newdata)
  expected_mean <- mean(pmin(training$time, 10))
  expected_variance <- mean((pmin(training$time, 10) - expected_mean)^2)
  expect_equal(unname(moments$mean[, 1]), rep(expected_mean, 12))
  expect_equal(unname(moments$variance[, 1]), rep(expected_variance, 12))
  expect_identical(fit$functions$censor_survival(newdata, seq(0, 10, length.out = 12)),
                   rep(1, 12))
  scores <- abs(pmin(calibration$time, 10) - moments$mean[, 1]) / moments$scale[, 1]
  q <- sort(unique(c(0, scores, 2)))
  membership <- 1 * outer(scores, q, "<=")
  cutoff <- ceiling((nrow(calibration) + 1) * 0.75)
  expected_threshold <- unname(sort(scores)[cutoff])
  for (method in c("ipw", "aipw")) {
    cal <- life_calibrate(fit, data = calibration, alpha = 0.25,
                          method = method, q = q,
                          control = list(keep_contributions = TRUE))
    expect_s3_class(cal, "lifeconformal_calibration")
    expect_identical(cal$method, method)
    expect_equal(unname(cal$contributions), unname(membership))
    expect_equal(cal$curve$estimate, colMeans(membership))
    expect_equal(cal$threshold, expected_threshold)
    expect_equal(cal$cutoff, cutoff / nrow(calibration))
    expect_false(cal$full_support)
    region <- predict(cal, data.frame(id = c("p1", "p2")))
    expect_s3_class(region, "lifeconformal_prediction")
    expect_identical(region$bounds$id, c("p1", "p2"))
    expect_identical(region$bounds$target, c("alive", "alive"))
    expect_equal(region$bounds$lower,
                 rep(max(0, expected_mean - expected_threshold * sqrt(expected_variance)), 2))
    expect_equal(region$bounds$upper,
                 rep(min(10, expected_mean + expected_threshold * sqrt(expected_variance)), 2))
    if (method == "aipw") expect_equal(cal$diagnostics$numerical$queried_snapshots, 0)
  }
})

test_that("formula and vector workflows agree on covariate fits and calibration", {
  train <- workflow_survival_data(censored = TRUE)
  caldata <- workflow_survival_data("cal", censored = TRUE)[seq_len(20), ]
  formula_fit <- life_train(LifeHist(time, status) ~ age + treatment, data = train,
                            absorbing = 1, targets = list(alive = 0), tau = 8)
  vector_fit <- life_train(x = train[c("age", "treatment")], id = train$id,
                           time = train$time, status = train$status,
                           absorbing = 1, targets = list(alive = 0), tau = 8)
  newdata <- train[c(8, 2, 19), c("age", "treatment")]
  expect_equal(predict(formula_fit, newdata)$mean, predict(vector_fit, newdata)$mean)
  expect_equal(predict(formula_fit, newdata)$variance, predict(vector_fit, newdata)$variance)
  expect_equal(formula_fit$functions$censor_survival(newdata, c(1, 3, 7)),
               vector_fit$functions$censor_survival(newdata, c(1, 3, 7)))
  q <- c(0, 0.4, 0.8, 1.4, 4)
  a <- life_calibrate(formula_fit, data = caldata, method = "ipw", alpha = 0.2, q = q,
                       control = list(keep_contributions = TRUE))
  b <- life_calibrate(vector_fit, x = caldata[c("treatment", "age")], id = caldata$id,
                       response = LifeHist(caldata$time, caldata$status),
                       method = "ipw", alpha = 0.2, q = q,
                       control = list(keep_contributions = TRUE))
  expect_equal(a$contributions, b$contributions)
  expect_equal(a$curve, b$curve)
  expect_equal(a$threshold, b$threshold)
  expect_equal(predict(a, newdata)$bounds, predict(b, newdata)$bounds)
})

test_that("event-status and explicit transitions produce equivalent analyses", {
  train <- workflow_survival_data()
  caldata <- workflow_survival_data("cal")
  a <- life_train(LifeHist(time, status) ~ 1, data = train,
                  absorbing = 1, targets = list(alive = 0), tau = 10)
  b <- life_train(id = train$id, start = rep(0, nrow(train)), stop = train$time,
                  from = rep(0, nrow(train)), to = rep(1, nrow(train)),
                  event = rep(1, nrow(train)),
                  absorbing = 1, targets = list(alive = 0), tau = 10)
  newdata <- data.frame(id = c("p1", "p2"))
  expect_equal(predict(a, newdata)$mean,
               predict(b, newdata, initial_state = 0)$mean)
  expect_error(predict(b, newdata), "Supply initial_state")
  ca <- life_calibrate(a, data = caldata, alpha = 0.25, method = "ipw", q = c(0, 1, 2))
  cb <- life_calibrate(b, id = caldata$id, start = rep(0, nrow(caldata)),
                       stop = caldata$time, from = rep(0, nrow(caldata)),
                       to = rep(1, nrow(caldata)), event = rep(1, nrow(caldata)),
                       alpha = 0.25, method = "ipw", q = c(0, 1, 2))
  expect_equal(ca$curve, cb$curve)
  expect_equal(predict(ca, newdata)$bounds,
               predict(cb, newdata, initial_state = 0)$bounds)
})

test_that("separate censoring predictors remain frozen and are required only where used", {
  train <- workflow_survival_data(censored = TRUE)
  caldata <- workflow_survival_data("cal", censored = TRUE)[seq_len(20), ]
  fit <- life_train(LifeHist(time, status) ~ age + treatment, data = train,
                    censor_formula = ~ z, absorbing = 1,
                    targets = list(alive = 0), tau = 8)
  expect_identical(fit$censor$predictor_names, "z")
  expect_identical(fit$outcome$predictor_names, c("age", "treatmentB"))
  expect_setequal(fit$schema$predictor_names, c("age", "treatment", "z"))
  expect_silent(predict(fit, train[seq_len(2), c("age", "treatment")]))
  expect_silent(fit$functions$censor_survival(train[seq_len(2), "z", drop = FALSE],
                                             c(2, 5)))
  expected <- fit$functions$censor_survival(train[seq_len(2), "z", drop = FALSE],
                                            c(2, 5))
  expect_equal(fit$functions$censor_survival(train[seq_len(2), ], c(2, 5)), expected)
  expect_error(life_calibrate(fit, data = caldata[setdiff(names(caldata), "z")],
                               method = "ipw", alpha = 0.2), "z")
  before <- serialize(fit, NULL)
  cal <- life_calibrate(fit, data = caldata, method = "aipw", alpha = 0.2,
                        q = c(0, 0.4, 0.8, 1.4, 4))
  expect_s3_class(cal, "lifeconformal_calibration")
  expect_true(all(is.finite(cal$curve$estimate)))
  expect_identical(serialize(fit, NULL), before)
  expect_equal(cal$fit$censor$model$coefficients, fit$censor$model$coefficients)
})

test_that("recurrent returns and subject-specific initial states survive the complete workflow", {
  train <- workflow_recurrent_data()
  caldata <- workflow_recurrent_data("cal")
  map <- c("1" = 1, "2" = 2, "7" = 0)
  fit <- life_train(LifeHist(time, status) ~ 1, data = train,
                    initial_state = ~ baseline_state, status_map = map,
                    absorbing = 2,
                    targets = list(home = 0, alive = c(0, 1), dead = 2), tau = 8)
  expect_identical(fit$schema$initial_state_rule$type, "column")
  expect_equal(fit$diagnostics$n_subjects, 12)
  expect_equal(fit$censor$diagnostics$excluded_origin_absorption, 4)
  labels <- fit$problem$states
  pairs <- paste(labels[fit$problem$working_graph$edges$from],
                  labels[fit$problem$working_graph$edges$to])
  expect_true("1 0" %in% pairs)
  newdata <- data.frame(id = c("home", "hospital", "dead"), baseline_state = c(0, 1, 2))
  moments <- predict(fit, newdata)
  expect_equal(unname(moments$mean[3, ]), c(0, 0, 8))
  expect_equal(unname(moments$variance[3, ]), c(0, 0, 0))
  expect_equal(unname(moments$scale[3, ]), rep(0.08, 3))
  expect_error(predict(fit, newdata["id"]), "baseline_state")
  expect_error(predict(fit, newdata, initial_state = c(0, 1, 7)), "Unknown initial state")
  for (method in c("aipw", "ipw")) {
    cal <- life_calibrate(fit, data = caldata, alpha = 0.25, method = method,
                          q = c(0, 0.5, 1, 2, 4))
    expect_equal(cal$n, 12)
    expect_equal(cal$diagnostics$resolved_fraction, rep(1, 5))
    pred <- predict(cal, newdata)
    expect_equal(pred$bounds$lower[pred$bounds$id == "dead"], c(0, 0, 8))
    expect_equal(pred$bounds$upper[pred$bounds$id == "dead"], c(0, 0, 8))
    expect_identical(pred$region$support$joint_support$A, fit$problem$A)
  }
})

test_that("vector initial-state rules require a new specification for new subjects", {
  train <- workflow_recurrent_data()
  caldata <- workflow_recurrent_data("cal")
  fit <- life_train(id = train$id, time = train$time, status = train$status,
                    initial_state = train$baseline_state,
                    status_map = c("1" = 1, "2" = 2, "7" = 0),
                    absorbing = 2, targets = list(alive = c(0, 1)), tau = 8)
  expect_error(predict(fit, data.frame(id = "new")), "Supply initial_state")
  expect_error(life_calibrate(fit, id = caldata$id, time = caldata$time,
                               status = caldata$status, alpha = 0.25), "Supply initial_state")
  cal <- life_calibrate(fit, id = caldata$id, time = caldata$time,
                        status = caldata$status, initial_state = caldata$baseline_state,
                        alpha = 0.25, method = "ipw", q = c(0, 1, 2, 4))
  expect_s3_class(predict(cal, data.frame(id = "new"), initial_state = 2),
                   "lifeconformal_prediction")
})

test_that("no-event outcome fits retain a positive scale and conservative support", {
  train <- data.frame(id = paste0("train", 1:12), time = rep(12, 12), status = 0)
  caldata <- transform(train, id = paste0("cal", 1:12))
  fit <- life_train(LifeHist(time, status) ~ 1, data = train,
                    absorbing = 1, targets = list(alive = 0), tau = 10)
  expect_equal(nrow(fit$problem$working_graph$edges), 0)
  expect_true(fit$problem$structural_support$adjacency["0", "1"])
  moments <- predict(fit, data.frame(id = "new"))
  expect_equal(unname(moments$mean), matrix(10, 1, 1))
  expect_equal(unname(moments$variance), matrix(0, 1, 1))
  expect_equal(unname(moments$scale), matrix(0.1, 1, 1))
  for (method in c("ipw", "aipw")) {
    cal <- life_calibrate(fit, data = caldata, alpha = 0.25, method = method,
                          q = c(0, 1, 100))
    expect_equal(cal$threshold, 0)
    expect_equal(cal$curve$estimate, rep(1, 3))
    expect_equal(predict(cal, data.frame(id = "new"))$bounds$lower, 10)
    expect_equal(predict(cal, data.frame(id = "new"))$bounds$upper, 10)
  }
})

test_that("rank fallback returns full support and obvious training reuse is rejected", {
  train <- workflow_survival_data()
  fit <- life_train(LifeHist(time, status) ~ 1, data = train,
                    absorbing = 1, targets = list(alive = 0), tau = 10)
  expect_error(life_calibrate(fit, data = train), "overlap training")
  independent <- workflow_survival_data("cal")[1:3, ]
  for (method in c("ipw", "aipw")) {
    cal <- life_calibrate(fit, data = independent, method = method, alpha = 0.1)
    expect_true(cal$full_support)
    expect_identical(cal$reason, "rank_exceeds_sample")
    expect_identical(cal$diagnostics$skipped, "rank_requires_full_support")
    expect_equal(cal$threshold, Inf)
    bounds <- predict(cal, data.frame(id = "new"))$bounds
    expect_equal(bounds$lower, 0)
    expect_equal(bounds$upper, 10)
    expect_equal(nrow(cal$curve), 0)
  }
  bypass <- life_calibrate(fit, data = train, alpha = 0.25, method = "ipw",
                           check_overlap = FALSE, q = c(0, 2))
  expect_false(bypass$diagnostics$overlap_checked)
  expect_equal(bypass$diagnostics$overlapping_labels, 12)
})

test_that("working-model certainty never substitutes for observed structural determination", {
  train <- data.frame(id = paste0("train", 1:12), time = 3, status = 0)
  caldata <- data.frame(id = paste0("cal", 1:12), time = 2, status = 0)
  fit <- life_train(LifeHist(time, status) ~ 1, data = train,
                    absorbing = 1, targets = list(alive = 0), tau = 10)
  # The empirical outcome model has no transitions, but a future event remains
  # structurally possible. At q=0, no censored subject can certify alive=10.
  ipw <- life_calibrate(fit, data = caldata, method = "ipw", alpha = 0.25,
                        q = c(0, 100), control = list(keep_contributions = TRUE))
  aipw <- life_calibrate(fit, data = caldata, method = "aipw", alpha = 0.25,
                         q = c(0, 100), control = list(keep_contributions = TRUE))
  expect_equal(ipw$diagnostics$resolved_fraction, c(0, 1))
  expect_equal(aipw$diagnostics$resolved_fraction, c(0, 1))
  expect_equal(ipw$curve$estimate, c(0, 1))
  expect_equal(aipw$curve$estimate, c(1, 1))
  expect_equal(ipw$threshold, 100)
  expect_equal(aipw$threshold, 0)
  expect_false(ipw$diagnostics$numerical$conditional_coverage_evaluated)
  expect_equal(aipw$diagnostics$numerical$queried_snapshots, 12)
  expect_equal(aipw$diagnostics$numerical$max_probability_bound_width, 0)
  empty <- predict(aipw, data.frame(id = character()))
  expect_equal(nrow(empty$bounds), 0)
  expect_identical(dim(empty$center), c(0L, 1L))
})

test_that("IPW does not invoke conditional coverage or augmentation engines", {
  train <- workflow_survival_data(censored = TRUE)
  caldata <- workflow_survival_data("cal", censored = TRUE)[seq_len(20), ]
  fit <- life_train(LifeHist(time, status) ~ age, data = train,
                    absorbing = 1, targets = list(alive = 0), tau = 8)
  local_mocked_bindings(
    evaluate_augmented_curve = function(...) stop("Unexpected augmentation evaluation"),
    conditional_life_coverage = function(...) stop("Unexpected conditional coverage"),
    markov_conditional_coverage_grid = function(...) stop("Unexpected reward distribution"),
    .package = "lifeconformal"
  )
  cal <- life_calibrate(fit, data = caldata, alpha = 0.2, method = "ipw",
                        q = c(0, 0.5, 1, 2, 4))
  expect_false(cal$diagnostics$numerical$conditional_coverage_evaluated)
  expect_s3_class(cal, "lifeconformal_calibration")
})

test_that("saved fitted functions and calibrated predictions reproduce their results", {
  train <- workflow_survival_data(censored = TRUE)
  caldata <- workflow_survival_data("cal", censored = TRUE)[seq_len(20), ]
  fit <- life_train(LifeHist(time, status) ~ age + treatment, data = train,
                    absorbing = 1, targets = list(alive = 0), tau = 8)
  cal <- life_calibrate(fit, data = caldata, alpha = 0.2, method = "ipw",
                        q = c(0, 0.5, 1, 2, 4))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(list(fit = fit, cal = cal), path)
  restored <- readRDS(path)
  newdata <- train[c(2, 7, 8), c("id", "age", "treatment")]
  expect_identical(restored$fit$functions$moments(newdata), fit$functions$moments(newdata))
  expect_identical(restored$fit$functions$censor_survival(newdata, c(0, 3, 8), "left"),
                   fit$functions$censor_survival(newdata, c(0, 3, 8), "left"))
  expect_identical(predict(restored$cal, newdata), predict(cal, newdata))
  expect_equal(fit$functions$moments(newdata), predict(fit, newdata))
  reordered <- predict(restored$cal, newdata)$bounds[c(3, 1, 2), ]
  rownames(reordered) <- NULL
  expect_identical(predict(restored$cal, newdata[c(3, 1, 2), ])$bounds, reordered)
})

test_that("public prediction and calibration reject incompatible changes", {
  train <- workflow_survival_data(censored = TRUE)
  caldata <- workflow_survival_data("cal", censored = TRUE)[seq_len(20), ]
  fit <- life_train(LifeHist(time, status) ~ treatment, data = train,
                    absorbing = 1, targets = list(alive = 0), tau = 8)
  expect_error(predict(fit, data.frame(treatment = "unseen")), "level")
  expect_error(predict(fit, data.frame(treatment = "A"), unused = TRUE), "Unused")
  expect_error(life_calibrate(fit, data = caldata, control = list(scale_floor = 0.1)),
               "cannot change trained")
  expect_error(life_calibrate(fit, data = caldata, control = list(survival_floor = 0.1)),
               "cannot change trained")
  expect_error(life_calibrate(fit, data = caldata, q = c(1, 0)), "strictly increasing")
  expect_error(life_calibrate(fit, data = caldata, q = c(0, 0)), "strictly increasing")
  expect_error(life_calibrate(fit, data = caldata, alpha = 0), "strictly between")
  bad <- caldata
  bad$status[1] <- 9
  expect_error(life_calibrate(fit, data = bad), "fixed dictionary")
  empty <- data.frame(treatment = factor(character(), levels = levels(train$treatment)))
  expect_identical(dim(predict(fit, empty)$mean), c(0L, 1L))
})
