test_that("complete and endpoint-weighted calibration agree with direct sums", {
  train <- data.frame(id = 1:12, time = c(.1, .2, .3, .4, .5, .6, rep(1, 6)),
                      status = c(0, 1, 0, 1, 0, 1, rep(0, 6)))
  fit <- life_train(LifeHist(time, status) ~ 1, train, absorbing = 1,
                    tau = 1, targets = list(alive = 0))
  full <- data.frame(id = 20:25, time = c(.15, .35, .7, 1, 1, 1),
                     status = c(1, 1, 1, 0, 0, 0))
  censored <- full
  censored$time[6] <- .25
  q <- c(0, .5, 1, 2, 5, 100)
  score <- predict(fit, full["id"], type = "score",
                   values = matrix(full$time, ncol = 1, dimnames = list(NULL, "alive")))
  truth <- 1 * outer(as.double(score), q, "<=")
  cal <- function(data, method) life_calibrate(fit, data, q = q, alpha = .3,
    method = method, control = list(keep_contributions = TRUE))
  expect_equal(unname(cal(full, "complete")$contributions), unname(truth))
  expect_error(cal(censored, "complete"), "through absorption or the horizon")
  expect_identical(cal(censored, "ip")$method, "ipw")
  # The three censoring jumps have 12, 10, and 8 training patients at risk.
  weight <- exp(as.vector(outer(full$time, c(.1, .3, .5), ">") %*% (1 / c(12, 10, 8))))
  weight[6] <- 0
  expect_equal(unname(cal(censored, "ipw_full")$contributions),
               unname(truth * weight), tolerance = 1e-12)
  for (method in c("complete", "ipw_full", "ipw", "aipw")) {
    data <- if (method == "complete") full else censored
    curve <- life_calibrate(fit, data, q = q, alpha = .3, method = method,
      control = list(calibration_accuracy = "probability", keep_contributions = TRUE))
    expect_identical(curve$curve$estimate, colMeans(curve$contributions))
    expect_identical(curve$curve$numerical_lower, curve$curve$estimate)
    expect_identical(curve$curve$numerical_upper, curve$curve$estimate)
    expect_identical(curve$diagnostics$numerical$curve_lower_sum, colSums(curve$contributions))
  }
})

test_that("the prespecified sparse-transition rule fits a baseline-only transition", {
  i <- seq_len(60)
  d <- data.frame(id = i, z = sin(i), time = (i %% 19 + 1) / 20,
                  status = as.integer(i %% 3L == 0L))
  train <- function(formula, cutoff) life_train(formula, d, absorbing = 1,
    tau = 1, targets = list(alive = 0),
    outcome_options = list(min_events_per_transition = cutoff))
  sparse <- train(LifeHist(time, status) ~ z, 21L)
  unadjusted <- train(LifeHist(time, status) ~ 1, 0L)
  x <- data.frame(z = c(-1, 1))
  score <- function(fit) as.double(predict(fit, x, type = "score",
    values = matrix(.7, 2L, 1L, dimnames = list(NULL, "alive"))))
  expect_equal(score(sparse), score(unadjusted), tolerance = 1e-12)
  expect_error(train(LifeHist(time, status) ~ z, -1), "nonnegative")
})

test_that("only a separated Cox transition uses the optional baseline refit", {
  i <- seq_len(100)
  entry <- data.frame(id = i, time = .1 + ((7*i) %% 23)/100, status = 1L, z = i %% 2L)
  exit_time <- .6 + ((11*i) %% 31)/100
  exits <- function(death) data.frame(id = i, time = ifelse(death, exit_time, 1),
                                      status = ifelse(death, 2L, 0L), z = i %% 2L)
  separated <- rbind(entry, exits(i %% 2L == 1L))
  ordinary <- rbind(entry, exits(i %% 3L == 1L))
  fit <- function(data, options = list(), formula = LifeHist(time, status) ~ z)
    life_train(formula, data[order(data$id, data$time), ], absorbing = 2, tau = 1,
               targets = list(alive = c(0, 1)), censor_formula = ~1,
               outcome_options = options)
  expect_error(fit(separated), "Transition 1->2.*may be infinite")
  policy <- list(min_events_per_transition = 12L, separation_action = "baseline")
  adjusted <- fit(separated, policy)
  unadjusted <- fit(separated, formula = LifeHist(time, status) ~ 1)
  before <- fit(ordinary)
  after <- fit(ordinary, policy)
  expect_identical(adjusted$outcome$coefficients[, 1], before$outcome$coefficients[, 1])
  expect_equal(unname(adjusted$outcome$coefficients[, 2]), 0)
  edge <- adjusted$outcome$models[[2]]
  expect_identical(edge$increments, unadjusted$outcome$models[[2]]$increments)
  expect_identical(edge$diagnostics$coefficient_fit, "separation_baseline_only")
  expect_match(edge$diagnostics$warnings, "may be infinite")
  expect_identical(before$outcome$base_increments, after$outcome$base_increments)
  expect_identical(before$outcome$coefficients, after$outcome$coefficients)
  expect_error(fit(ordinary, c(policy, list(iter.max = 1L))), "iteration limit|did not converge")
  expect_error(fit(ordinary, list(separation_action = "ignore")), "separation_action")
  expect_error(.cox_resolve_options(policy, censoring = TRUE), "supported controls")
  censor_separation <- data.frame(id = i, time = ifelse(i %% 2L, exit_time, 1),
                                  status = 0L, z = i %% 2L)
  expect_error(life_train(LifeHist(time, status) ~ z, censor_separation,
    absorbing = 1, tau = 1, targets = list(alive = 0), censor_formula = ~z,
    outcome_options = policy), "Censoring Cox model failed.*may be infinite")
})
