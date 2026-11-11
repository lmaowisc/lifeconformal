calibration_sum_test_call <- function(operation, ...) {
  .Call(paste0("life_calibration_sum_", operation), ..., PACKAGE = "lifeconformal")
}

with_calibration_sum_test <- function(ncols, code,
    precision = calibration_sum_precision()) {
  skip_if(is.null(precision), "R and native accumulation formats differ")
  state <- calibration_sum_test_call("create", as.integer(ncols), precision)
  on.exit(calibration_sum_test_call("close", state), add = TRUE)
  code(state)
}

test_that("equal mantissas with a wider R exponent range retain the legacy path", {
  # The constructor selects double when its mantissa count matches. A wider
  # exponent range in R's long double therefore still requires fallback.
  expect_null(calibration_sum_precision(
    native = c(53L, 53L, 1024L, 16384L, 8L, 16L, 0L),
    machine = list(longdouble.digits = 53L, longdouble.max.exp = 16384L),
    wide = TRUE))
})

test_that("chunked calibration sums retain R's product and sum rounding", {
  set.seed(618204)
  fixtures <- list(
    cancellation = list(values = cbind(rep(0, 7), rep(1, 7), seq(0, 1, length.out = 7)),
      weights = c(2^53, 1, -2^53, 2^53, -2^53, -2^-50, 2^-51)),
    product = list(values = matrix(c(.5 + 2^-31, 1), 2), weights = c(2 - 2^-29, -1)),
    mean = list(values = matrix(c(1, 1, .Machine$double.eps), 3), weights = rep(1, 3)),
    overflow = list(values = matrix(1, 4, 1), weights = rep(.Machine$double.xmax, 4)),
    random = list(values = matrix(runif(155), 31, 5),
      weights = runif(31, -1, 1) * 2^sample(-100:100, 31, replace = TRUE)))
  for (fixture in fixtures) {
    x <- fixture$values; w <- fixture$weights
    before <- serialize(fixture, NULL)
    expected_sum <- unname(colSums(x * w))
    expected_mean <- unname(colMeans(x * w))
    with_calibration_sum_test(ncol(x), function(state) {
      calibration_sum_test_call("add", state, x, w)
      expect_identical(calibration_sum_test_call("value", state, 1), expected_sum)
      expect_identical(calibration_sum_test_call("value", state, as.double(nrow(x))), expected_mean)
      for (size in c(1L, 3L)) {
        calibration_sum_test_call("reset", state)
        for (begin in seq.int(1L, nrow(x), by = size)) {
          rows <- seq.int(begin, min(nrow(x), begin + size - 1L))
          calibration_sum_test_call("add", state, x[rows, , drop = FALSE], w[rows])
        }
        expect_identical(calibration_sum_test_call("value", state, 1), expected_sum)
        expect_identical(calibration_sum_test_call("value", state, as.double(nrow(x))), expected_mean)
      }
      retained <- calibration_sum_test_call("value", state, 1)
      calibration_sum_test_call("reset", state)
      expect_identical(calibration_sum_test_call("value", state, 1), rep(0, ncol(x)))
      expect_identical(retained, expected_sum)
    })
    expect_identical(serialize(fixture, NULL), before)
  }
  expect_identical(colSums(fixtures$product$values * fixtures$product$weights), 0)
  if (isTRUE(capabilities("long.double")) &&
      .Machine$longdouble.digits > .Machine$double.digits) {
    # Dividing an already rounded sum loses a bit required by colMeans here.
    expect_false(identical(colMeans(fixtures$mean$values),
      colSums(fixtures$mean$values) / 3))
  }
})

test_that("subject endpoint addition precedes global extended accumulation", {
  with_calibration_sum_test(1L, function(query) {
    calibration_sum_test_call("add", query, matrix(1, 1, 1), -1)
    calibration_sum_test_call("add", query, matrix(1, 1, 1), -2^-54)
    # Keeping this negative-query sum extended through endpoint addition would
    # incorrectly return -2^-54 on an extended-precision platform.
    endpoint <- 1
    contribution <- endpoint + calibration_sum_test_call("value", query, 1)
    expect_identical(contribution, 0)
    expected <- endpoint + colSums(matrix(1, 2, 1) * c(-1, -2^-54))
    expect_identical(contribution, expected)
    subjects <- cbind(c(contribution, 1, 1 + .Machine$double.eps), c(.25, .5, .75))
    with_calibration_sum_test(ncol(subjects), function(global) {
      for (i in seq_len(nrow(subjects)))
        calibration_sum_test_call("add", global, subjects[i, , drop = FALSE], 1)
      expect_identical(calibration_sum_test_call("value", global, 1), unname(colSums(subjects)))
      expect_identical(calibration_sum_test_call("value", global, 3), unname(colMeans(subjects)))
    })
  })
})

test_that("double accumulation fallback follows sequential binary64 arithmetic", {
  x <- cbind(rep(1, 7), c(0, .5, 1, .75, .25, 1, 0))
  weights <- c(2^53, 1, -2^53, 4, -4, 2^-53, -2^-54)
  products <- x * weights
  expected <- vapply(seq_len(ncol(x)), function(j)
    Reduce(`+`, products[, j], init = 0), numeric(1))
  with_calibration_sum_test(ncol(x), function(state) {
    for (i in seq_len(nrow(x)))
      calibration_sum_test_call("add", state, x[i, , drop = FALSE], weights[i])
    expect_identical(calibration_sum_test_call("value", state, 1), expected)
    expect_identical(calibration_sum_test_call("value", state, 7), expected / 7)
  }, precision = as.integer(.Machine$double.digits))
})

test_that("empty calibration accumulators preserve column-reduction conventions", {
  with_calibration_sum_test(2L, function(state) {
    x <- matrix(numeric(), 0, 2)
    calibration_sum_test_call("add", state, x, numeric())
    expect_identical(calibration_sum_test_call("value", state, 1), colSums(x))
    expect_identical(calibration_sum_test_call("value", state, 0), colMeans(x))
  })
  with_calibration_sum_test(0L, function(state) {
    calibration_sum_test_call("add", state, matrix(numeric(), 3, 0), rep(1, 3))
    expect_identical(calibration_sum_test_call("value", state, 1), numeric())
    expect_identical(calibration_sum_test_call("value", state, 3), numeric())
  })
})

test_that("native calibration sums reject malformed inputs without partial updates", {
  precision <- calibration_sum_precision()
  skip_if(is.null(precision), "R and native accumulation formats differ")
  for (ncols in list(2, -1L, NA_integer_, c(1L, 2L), matrix(2L)))
    expect_error(calibration_sum_test_call("create", ncols, precision), "ncols")
  unsupported <- max(calibration_sum_test_call("info")[1:2]) + 1L
  for (digits in list(53, NA_integer_, unsupported, c(53L, 53L)))
    expect_error(calibration_sum_test_call("create", 2L, digits), "precision")
  with_calibration_sum_test(2L, function(state) {
    calibration_sum_test_call("add", state, matrix(c(1, 2, 3, 4), 2), c(1, -1))
    before <- calibration_sum_test_call("value", state, 1)
    for (values in list(matrix(1L, 2, 2), rep(1, 4), matrix(1, 2, 3),
      matrix(c(1, 2, 3, NA_real_), 2), matrix(c(1, 2, 3, Inf), 2))) {
      expect_error(calibration_sum_test_call("add", state, values, c(1, -1)), "values")
      expect_identical(calibration_sum_test_call("value", state, 1), before)
    }
    for (weights in list(c(1L, 2L), 1, matrix(c(1, 2)), c(1, NA_real_), c(1, Inf))) {
      expect_error(calibration_sum_test_call("add", state, matrix(1, 2, 2), weights), "row_weights")
      expect_identical(calibration_sum_test_call("value", state, 1), before)
    }
    expect_error(calibration_sum_test_call("add", state,
      matrix(c(0, 0, .Machine$double.xmax, .Machine$double.xmax), 2), c(2, 2)), "products")
    expect_identical(calibration_sum_test_call("value", state, 1), before)
    for (divisor in list(1L, c(1, 2), NA_real_, Inf, -1, 0))
      expect_error(calibration_sum_test_call("value", state, divisor), "divisor")
    foreign <- getNativeSymbolInfo("life_calibration_sum_create", PACKAGE = "lifeconformal")$address
    stale <- unserialize(serialize(state, NULL))
    for (pointer in list(NULL, 1, methods::new("externalptr"), foreign, stale)) {
      expect_error(calibration_sum_test_call("add", pointer, matrix(0, 1, 2), 1), "live")
      expect_error(calibration_sum_test_call("value", pointer, 1), "live")
      expect_error(calibration_sum_test_call("reset", pointer), "live")
    }
  })
})

test_that("calibration sum cleanup releases states on error and garbage collection", {
  precision <- calibration_sum_precision()
  skip_if(is.null(precision), "R and native accumulation formats differ")
  gc()
  baseline <- calibration_sum_test_call("info")[[7L]]
  expect_error(with_calibration_sum_test(2L, function(state) {
    calibration_sum_test_call("add", state, matrix(c(1, 2, 3, NA_real_), 2), c(1, 1))
  }), "values")
  expect_identical(calibration_sum_test_call("info")[[7L]], baseline)
  state <- calibration_sum_test_call("create", 2L, precision)
  calibration_sum_test_call("close", state)
  expect_error(calibration_sum_test_call("value", state, 1), "live")
  expect_error(calibration_sum_test_call("close", state), "live")
  local({
    states <- lapply(seq_len(3L), function(i) calibration_sum_test_call("create", i, precision))
    expect_length(states, 3L)
    expect_identical(calibration_sum_test_call("info")[[7L]], baseline + 3L)
  })
  gc(); gc()
  expect_identical(calibration_sum_test_call("info")[[7L]], baseline)
})
