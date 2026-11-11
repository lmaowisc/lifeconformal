output_method_fixture <- function() {
  data <- life_example_data(72, seed = 411)
  fit <- life_train(LifeHist(time, status) ~ age, data = data[data$id <= 48, ],
    absorbing = 2, targets = list(alive = 0, dead = 2), tau = 5)
  calibration <- data[data$id > 48, ]
  joint <- life_calibrate(fit, data = calibration, method = "ipw", alpha = .25,
                         q = c(0, .5, 1, 2, 100))
  partial <- life_calibrate(fit, data = calibration, method = "ipw", alpha = .25,
    alpha_components = c(alive = .2, dead = .001), q = c(0, .5, 1, 2, 100))
  full <- life_calibrate(fit, data = calibration, method = "ipw", alpha = .001)
  list(fit = fit, joint = joint, partial = partial, full = full)
}

output_pdf <- function(code) {
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file, width = 8, height = 6)
  on.exit({ grDevices::dev.off(); unlink(file) }, add = TRUE)
  force(code)
}

test_that("summaries retain effective model settings and component meaning", {
  z <- output_method_fixture()
  fit <- summary(z$fit)
  expect_named(fit$model_options, c("outcome", "censoring", "augmentation"))
  expect_named(fit$model_diagnostics, c("outcome", "censoring", "augmentation"))
  expect_identical(fit$model_options$outcome, z$fit$outcome$controls)
  expect_identical(fit$model_options$censoring, z$fit$censor$controls)
  expect_identical(fit$model_diagnostics$augmentation, z$fit$augmentation$diagnostics)
  partial <- summary(z$partial)
  expect_identical(partial$alpha_components, z$partial$alpha_components)
  expect_identical(partial$allocated_alpha, sum(z$partial$alpha_components))
  expect_identical(partial$component_full_support, c(alive = FALSE, dead = TRUE))
  expect_identical(partial$component_reason, z$partial$component_reason)
  expect_true(partial$any_component_full_support)
  expect_identical(partial$rank, z$partial$rank)
  expect_identical(summary(z$full)$reason, "rank_exceeds_sample")
  prediction <- predict(z$partial, data.frame(id = c("A", "B"), age = c(45, 65)))
  shown <- summary(prediction)
  expect_identical(shown$threshold, z$partial$threshold)
  expect_identical(shown$alpha_components, z$partial$alpha_components)
  expect_identical(shown$component_full_support, z$partial$component_full_support)
  expect_identical(shown$coordinate_bounds, "conservative_floating_score_envelope")
})

test_that("print and data-frame outputs distinguish joint regions and fallbacks", {
  z <- output_method_fixture()
  expect_output(result <- print(z$fit), "Life-history training fit")
  expect_identical(result, z$fit)
  expect_output(print(z$partial), "Allocated total miscoverage")
  expect_output(print(z$partial), "rank_exceeds_sample")
  expect_output(print(z$full), "requested rank exceeds")
  prediction <- predict(z$partial, data.frame(id = factor(c("A", "B")), age = c(45, 65)))
  expect_output(result <- print(prediction), "Full support for: dead")
  expect_output(print(prediction), "joint region")
  expect_identical(result, prediction)
  table <- as.data.frame(prediction, row.names = paste0("r", 1:4))
  expect_identical(names(table), c("id", "target", "lower", "upper"))
  expect_identical(rownames(table), paste0("r", 1:4))
  expect_identical(table$id, rep(factor(c("A", "B")), each = 2L))
  expect_identical(table$lower, prediction$bounds$lower)
  full <- predict(z$full, data.frame(age = 50))
  expect_output(print(full), "Full target support returned")
  empty <- predict(z$joint, data.frame(age = numeric()))
  expect_equal(nrow(as.data.frame(empty)), 0L)
  expect_identical(summary(empty)$n, 0L)
  expect_error(plot(empty), "no prediction rows")
})

test_that("calibration plots accept ordinary label, type and style overrides", {
  z <- output_method_fixture()
  output_pdf({
  for (object in list(z$joint, z$partial, z$full)) {
    expect_no_warning(result <- plot(object, xlab = "Score cutoff", ylab = "Coverage",
      main = "Custom title", type = "p", col = "navy", pch = 16, lwd = 2))
    expect_identical(result, object)
  }
  expect_no_warning(plot(z$partial, col = c("navy", "orange"), lty = c(1, 3),
                         lwd = c(1, 2), pch = c(16, 17), type = "b", legend = FALSE))
  expect_error(plot(z$joint, legend = NA), "legend must")
  expect_error(plot(z$joint, "unnamed"), "unique, nonempty names")
  expect_error(do.call(plot, c(list(x = z$joint), list(xlab = "a", xlab = "b"))),
               "unique, nonempty names")
  # The actual cutoff must remain visible even outside the observed curve range.
  scalar <- z$joint
  scalar$curve$estimate <- rep(-.5, nrow(scalar$curve))
  expect_no_warning(plot(scalar))
  limits <- graphics::par("usr")[3:4]
  expect_lte(limits[1], -.5)
  expect_gte(limits[2], scalar$cutoff)
  })
})

test_that("prediction plot styling reaches intervals and suppressed axes stay suppressed", {
  z <- output_method_fixture()
  prediction <- predict(z$joint, data.frame(id = c("A", "B"), age = c(45, 65)))
  output_pdf({
  expect_no_warning(result <- plot(prediction, xlab = "Years", ylab = "Patient",
    main = "Custom intervals", type = "p", pch = 16, col = "navy", lwd = 3))
  expect_identical(result, prediction)
  segments <- NULL
  axes <- 0L
  local_mocked_bindings(
    segments = function(...) { segments <<- list(...) },
    axis = function(...) { axes <<- axes + 1L }, .package = "graphics")
  expect_no_warning(plot(prediction, col = "orange", lty = 3, lwd = 4, axes = FALSE))
  expect_identical(segments$col, "orange")
  expect_identical(segments$lty, 3)
  expect_identical(segments$lwd, 4)
  expect_identical(axes, 0L)
  })
})

test_that("singleton and zero-threshold predictions remain plottable", {
  z <- output_method_fixture()
  z$joint$threshold <- 0
  zero <- predict(z$joint, data.frame(age = 50))
  dead <- predict(z$joint, data.frame(age = 50), initial_state = 2)
  output_pdf({
  expect_no_warning(plot(zero, xlab = "Years"))
  points <- NULL
  local_mocked_bindings(points = function(...) { points <<- list(...) },
                         .package = "graphics")
  expect_no_warning(plot(dead))
  expect_identical(points[[1]], c(0, 5))
  expect_identical(points$pch, 16)
  })
})

test_that("default horizontal labels fit margins and explicit label size is honored", {
  z <- output_method_fixture()
  prediction <- predict(z$joint, data.frame(id = "A long subject identifier", age = 50))
  output_pdf({
    label_size <- NULL
    original_margins <- graphics::par("mar")
    original_label_size <- graphics::par("cex.axis")
    local_mocked_bindings(axis = function(side, ...) {
      if (side == 2) label_size <<- list(...)$cex.axis
    }, .package = "graphics")
    plot(prediction)
    expect_lt(label_size, graphics::par("cex.axis"))
    plot(prediction, cex.axis = .7)
    expect_identical(label_size, .7)
    expect_identical(graphics::par("mar"), original_margins)
    expect_identical(graphics::par("cex.axis"), original_label_size)
  })
})
