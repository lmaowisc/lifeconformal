rounding_region <- function(center, scale, threshold, tau = 10,
                             reward_min = 0, reward_max = 1) {
  if (is.null(dim(center))) center <- matrix(center, ncol = 1L)
  if (is.null(dim(scale))) scale <- matrix(scale, nrow(center), ncol(center))
  targets <- paste0("target", seq_len(ncol(center)))
  colnames(center) <- colnames(scale) <- targets
  support <- list(reward_min = matrix(rep(reward_min, length.out = ncol(center)), 1L),
                  reward_max = matrix(rep(reward_max, length.out = ncol(center)), 1L))
  object <- structure(list(fit = list(problem = list(tau = tau, structural_support = support)),
    threshold = threshold, alpha = .1, method = "ipw", full_support = all(is.infinite(threshold)),
    diagnostics = list(numerical = list())), class = "lifeconformal_calibration")
  local_mocked_bindings(predict_life_moments = function(...) {
    list(mean = center, scale = scale, initial_state = rep(1L, nrow(center)))
  })
  predict(object, data.frame(id = seq_len(nrow(center))))
}

rounding_neighbors <- function(x, steps = 8L) {
  values <- x
  up <- down <- x
  for (i in seq_len(steps)) {
    up <- cdt_nextafter(up)
    down <- -cdt_nextafter(-down)
    values <- c(values, up, down)
  }
  unique(values[is.finite(values)])
}

test_that("prediction envelopes include rounded score equality beyond both naive endpoints", {
  # Hexadecimal literals preserve the intended binary64 fixture on platforms
  # with different precision in decimal-to-double conversion.
  center <- 0x1.e3b614aap-1
  scale <- 0x1.9ef4feab8f5c3p-1
  q <- 0x1.1bf20eec8p+1
  naive_upper <- center + q * scale
  value <- cdt_nextafter(naive_upper)
  expect_gt(value, naive_upper)
  expect_identical(abs((value - center) / scale), q)
  region <- rounding_region(center, scale, q)
  expect_gte(region$bounds$upper, value)
  expect_identical(region$region$coordinate_bounds,
                   "conservative_floating_score_envelope")

  center <- center + 2
  naive_lower <- center - q * scale
  value <- -cdt_nextafter(-naive_lower)
  expect_lt(value, naive_lower)
  expect_gt(value, 0)
  expect_identical(abs((value - center) / scale), q)
  region <- rounding_region(center, scale, q)
  expect_lte(region$bounds$lower, value)
})

test_that("zero cutoffs include division-underflow score atoms", {
  # A nonzero distance divided by a very large frozen relative scale can
  # round to zero. The algebraic width q*scale is zero and would omit it.
  value <- 1e-20
  expect_identical(value / 1e307, 0)
  region <- rounding_region(0, 1e307, 0)
  expect_gte(region$bounds$upper, value)
  expect_identical(region$bounds$lower, 0)
  negative_side <- rounding_region(1e-20, 1e307, 0)
  expect_identical(abs((0 - 1e-20) / 1e307), 0)
  expect_identical(negative_side$bounds$lower, 0)
})

test_that("outward bounds contain admitted neighboring values across numeric scales", {
  settings <- expand.grid(center = c(0, 1e-200, .5, 2, 1e200),
    scale = c(.Machine$double.xmin, 1e-200, .8, 1e200, 1e307),
    q = c(0, .Machine$double.xmin, .125, 2.218324532965198, 1e200))
  missed <- 0L
  admitted <- 0L
  for (i in seq_len(nrow(settings))) {
    z <- settings[i, ]
    candidates <- unique(c(rounding_neighbors(z$center - z$q * z$scale),
                           rounding_neighbors(z$center + z$q * z$scale),
                           z$center, 0, 1e-320, 1e-20, 1e200, 1e308))
    candidates <- candidates[candidates >= 0 & candidates <= 1e308]
    score <- abs((candidates - z$center) / z$scale)
    covered <- candidates[score <= z$q]
    region <- rounding_region(z$center, z$scale, z$q, tau = 1e308)
    missed <- missed + sum(covered < region$bounds$lower | covered > region$bounds$upper)
    admitted <- admitted + length(covered)
  }
  expect_gt(admitted, 500)
  expect_identical(missed, 0L)
})

test_that("full support, deterministic coordinates and empty rows remain exact", {
  region <- rounding_region(matrix(c(2, 3, 4, 5), 2L),
                             matrix(c(.5, 1, 2, 3), 2L), c(Inf, 0))
  expect_identical(region$bounds$lower[region$bounds$target == "target1"], c(0, 0))
  expect_identical(region$bounds$upper[region$bounds$target == "target1"], c(10, 10))
  expect_true(all(is.finite(region$bounds$lower)))
  expect_true(all(is.finite(region$bounds$upper)))
  deterministic <- rounding_region(matrix(c(0, 10), 1L), matrix(c(1, 1), 1L),
                                    c(0, 0), reward_min = c(0, 1), reward_max = c(0, 1))
  expect_identical(deterministic$bounds$lower, c(0, 10))
  expect_identical(deterministic$bounds$upper, c(0, 10))
  empty <- rounding_region(matrix(numeric(), 0L, 2L), matrix(numeric(), 0L, 2L), c(0, Inf))
  expect_equal(nrow(empty$bounds), 0L)
  expect_identical(dim(empty$center), c(0L, 2L))
})


test_that("public calibration and prediction retain a computed-score boundary", {
  train <- data.frame(id = paste0("train", 1:8),
    time = c(.1, .2, .4, .8, 1.2, 1.8, 2.5, 4), status = 2)
  fit <- life_train(LifeHist(time, status) ~ 1, data = train,
    absorbing = 2, targets = list(alive = 0), tau = 5)
  patient <- data.frame(id = "new")
  moments <- predict(fit, patient)
  # The fitted moments can vary by an ULP across supported math libraries.
  # Select a genuine omitted neighbor for those actual moments, retaining
  # the strict boundary predicate rather than weakening it to a tolerance.
  candidates <- rounding_neighbors(0x1.59a6046f29689p+1, steps = 128L)
  scores <- abs((candidates - moments$mean[1]) / moments$scale[1])
  witnesses <- candidates[candidates > moments$mean[1] + scores * moments$scale[1]]
  expect_gt(length(witnesses), 0L)
  value <- witnesses[1L]
  q <- predict(fit, patient, type = "score", values = value)
  # Establish that ordinary mean + q*scale would omit this admitted target.
  expect_gt(value, moments$mean[1] + q * moments$scale[1])
  calibration <- data.frame(id = paste0("cal", 1:12), time = value, status = 2)
  for (method in c("ipw", "aipw")) {
    cal <- life_calibrate(fit, data = calibration, alpha = .25, method = method,
                          q = c(0, q, q + 1))
    expect_identical(cal$threshold, q)
    region <- predict(cal, patient)
    expect_lte(region$bounds$lower, value)
    expect_gte(region$bounds$upper, value)
    expect_identical(region$threshold, cal$threshold)
    expect_identical(region$center, moments$mean)
    expect_identical(region$scale, moments$scale)
  }
})
