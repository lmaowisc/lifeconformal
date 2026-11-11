test_that("malformed controls fail before entering numerical models", {
  expect_error(resolve_life_controls(setNames(list(1), NA_character_)),
               "unique, nonempty names")
  expect_error(resolve_life_controls(list(grid_size = 1.5)), "integer")
  expect_error(resolve_life_controls(list(probability_tolerance = 0)), "in")
  expect_error(resolve_life_controls(list(unknown = 1)), "Unknown control")
  for (key in c("scale_floor", "survival_floor", "probability_tolerance",
                "grid_size", "max_nodes", "draws")) {
    for (value in list(1i, matrix(.5), numeric(), NA_real_))
      expect_error(resolve_life_controls(setNames(list(value), key)), key)
  }
  expect_error(resolve_life_controls(list(keep_contributions = matrix(TRUE))),
    "keep_contributions")
})

test_that("explicit censor floors agree with stored public diagnostics", {
  data <- life_example_data(60, seed = 371)
  fit <- life_train(LifeHist(time, status) ~ 1, data = data,
                    absorbing = 2, targets = list(alive = 0), tau = 5,
                    censor_options = list(survival_floor = 0.95))
  expect_equal(fit$control$survival_floor, 0.95)
  expect_equal(fit$censor$controls$survival_floor, 0.95)
  g <- fit$functions$censor_survival(data.frame(id = "new"), 5)
  expect_gte(g, 0.95)
  calibration <- life_example_data(30, seed = 372)
  calibration$id <- paste0("cal", calibration$id)
  result <- life_calibrate(fit, calibration, q = c(0, 1, 2, 4))
  expect_equal(result$control$survival_floor, 0.95)
  expect_true(all(is.finite(result$curve$estimate)))
})
