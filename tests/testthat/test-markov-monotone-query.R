monotone_query_fixture <- function() {
  list(jump_times = c(.25, .5), kernels = NULL,
    rewards = matrix(c(1, 0), 2, 1), tau = 1,
    query_times = c(0, .25, .5, 1), query_state = rep(1, 4),
    accrued = matrix(c(0, .25, .5, 1), 4, 1),
    center = matrix(0, 4, 1), scale = matrix(1, 4, 1),
    q = c(0, .25 - 2^-55, .25, .25 + 2^-54, .5 - 2^-54, .5, .5 + 2^-53, 1),
    side = "right", grid_size = 16L, max_nodes = 100000L, exact_limit = 32L,
    point_mass_threshold = .1, tolerance = .01, max_refinements = 0L, strict = FALSE,
    kernel_factors = list(base_increments = matrix(c(.25, .5), 2, 1),
      relative_risk = 1, from = 1, to = 2, absorbing = 2))
}

monotone_query_spec <- function(a) {
  list(jump_times = as.double(a$jump_times), rewards = as.double(a$rewards),
    tau = as.double(a$tau), kernel_factors = lapply(a$kernel_factors, function(x) {
      if (is.numeric(x)) storage.mode(x) <- "double"
      x
    }), query_times = as.double(a$query_times), query_state = as.double(a$query_state),
    accrued = as.double(a$accrued), center = as.double(a$center[1]),
    scale = as.double(a$scale[1]), q = as.double(a$q), max_cells = 1e6,
    max_kernel_cells = 1e6, inspect_kernels = FALSE)
}

monotone_query_native <- function(x) {
  .Call("life_markov_monotone_queries", x, PACKAGE = "lifeconformal")
}

test_that("installed scalar clock preserves exact exit atoms and through-time queries", {
  a <- monotone_query_fixture()
  expected <- rbind(c(0, 0, .25, .25, .25, .625, .625, 1),
    c(0, 0, 0, 0, 0, .5, .5, 1), c(rep(0, 7), 1), c(rep(0, 7), 1))
  raw <- monotone_query_native(monotone_query_spec(a))
  value <- do.call(markov_marginal_scalar_coverage, a)
  expect_identical(raw$probability, expected)
  expect_identical(value$probability, expected)
  expect_identical(value$lower, value$upper)
  expect_identical(value$lower, expected)
  expect_true(value$diagnostics$monotone_query_applied)
  expect_true(is.na(value$diagnostics$grid_size))
  expect_identical(value$diagnostics$peak_nodes, 0)
  for (order in list(c(4, 1, 3, 2), c(1, 1, 4, 2, 3))) {
    b <- a
    for (f in c("query_times", "query_state")) b[[f]] <- a[[f]][order]
    for (f in c("accrued", "center", "scale")) b[[f]] <- a[[f]][order, , drop = FALSE]
    expect_identical(do.call(markov_marginal_scalar_coverage, b)$probability,
      expected[order, , drop = FALSE])
  }
  for (i in seq_along(a$query_times)) {
    b <- a
    for (f in c("query_times", "query_state")) b[[f]] <- a[[f]][i]
    for (f in c("accrued", "center", "scale")) b[[f]] <- a[[f]][i, , drop = FALSE]
    expect_identical(do.call(markov_marginal_scalar_coverage, b)$probability,
      expected[i, , drop = FALSE])
  }
  a$jump_times <- .3
  a$tau <- 1 + .Machine$double.eps
  a$query_times <- c(.Machine$double.eps / 2, .3, a$tau)
  a$query_state <- rep(1, 3)
  a$accrued <- matrix(a$query_times, 3, 1)
  a$center <- matrix(.3, 3, 1)
  a$scale <- matrix(1, 3, 1)
  a$q <- c(0, 1)
  a$kernel_factors$base_increments <- matrix(1, 1, 1)
  expect_identical(do.call(markov_marginal_scalar_coverage, a)$probability,
    rbind(c(1, 1), c(0, 1), c(0, 1)))
  a$jump_times <- a$tau
  a$center[] <- a$tau
  expect_identical(do.call(markov_marginal_scalar_coverage, a)$probability, matrix(1, 3, 2))
})

test_that("scalar clocks share generic-grid kernel construction under row rescaling", {
  a <- monotone_query_fixture()
  a$rewards <- matrix(c(1, 1, 0), 3, 1)
  a$kernel_factors <- list(base_increments = matrix(c(.2, .3, .7, .1, .4, .5, .2, .1), 2, 4),
    relative_risk = c(2.2, 1.3, 5.1, 2.9), from = c(1, 1, 2, 2),
    to = c(2, 3, 1, 3), absorbing = 3)
  for (multiplier in c(0, .1, 1, 10)) {
    b <- a
    b$kernel_factors$base_increments <- a$kernel_factors$base_increments * multiplier
    f <- b$kernel_factors
    kernels <- lapply(seq_along(b$jump_times), function(j) {
      inc <- matrix(0, 3, 3)
      inc[cbind(f$from, f$to)] <- f$base_increments[j, ] * f$relative_risk
      make_markov_kernel(inc, 3)$kernel
    })
    kernels <- prepare_markov_recursion(b$jump_times, kernels, b$rewards, b$tau,
      0, "left", NULL)$kernels
    spec <- monotone_query_spec(b)
    spec$inspect_kernels <- TRUE
    native <- monotone_query_native(spec)
    expect_identical(native$kernel_storage, as.double(unlist(kernels)))
    clock_args <- b[c("jump_times", "rewards", "tau", "query_times", "query_state",
      "accrued", "center", "scale", "q", "side")]
    clock_args$kernels <- kernels
    oracle <- do.call(markov_conditional_coverage_monotone, clock_args)
    expect_equal(native$probability, oracle$probability, tolerance = 3e-15)
    b$exact_limit <- 100000L
    exact <- do.call(markov_conditional_coverage_grid, b)
    expect_equal(native$probability, exact$probability, tolerance = 3e-15)
  }
})

test_that("unsupported scalar clock inputs preserve the original validator", {
  base <- monotone_query_fixture()
  for (kind in c("incoherent", "left", "center", "scale", "increasing", "strict",
                 "grid", "duplicate", "unnamed", "known", "weights", "missing_tolerance", "state_labels")) {
    a <- base
    if (kind == "incoherent") a$accrued[2, 1] <- a$query_times[2] + .Machine$double.eps
    if (kind == "left") a$side <- "left"
    if (kind == "center") a$center[2, 1] <- .1
    if (kind == "scale") a$scale[2, 1] <- 2
    if (kind == "increasing") a$rewards[, 1] <- c(0, 1)
    if (kind == "strict") a$strict <- TRUE
    if (kind == "grid") a$grid_size <- 32L
    if (kind == "duplicate") a <- c(a, list(tolerance = .02))
    if (kind == "unnamed") a <- c(a, list(.02))
    if (kind == "known") a$query_known <- matrix(NA_real_, 4, length(a$q))
    if (kind == "weights") a$query_weights <- c(-.5, 0, -.5, 2)
    if (kind == "missing_tolerance") a$tolerance <- NULL
    if (kind == "state_labels") rownames(a$rewards) <- c("duplicate", "duplicate")
    value <- tryCatch(do.call(markov_marginal_scalar_coverage, a), error = identity)
    expected <- tryCatch(do.call(markov_conditional_coverage_grid, a), error = identity)
    if (inherits(expected, "error")) {
      expect_s3_class(value, "error")
      expect_identical(conditionMessage(value), conditionMessage(expected))
    } else expect_identical(value, expected)
  }
})

test_that("native scalar clock caps precede allocation and errors unwind safely", {
  a <- monotone_query_fixture()
  original <- monotone_query_spec(a)
  expected <- monotone_query_native(original)
  for (kind in c("cells", "kernels", "empty_query", "empty_q", "state", "active_clock",
                 "risk_count", "from_count", "to_count", "absorbing_count", "underflow")) {
    x <- original
    if (kind == "cells") x$max_cells <- 1
    if (kind == "kernels") x$max_kernel_cells <- 1
    if (kind == "empty_query") x$query_times <- numeric()
    if (kind == "empty_q") x$q <- numeric()
    if (kind == "state") x$query_state[1] <- NA_real_
    if (kind == "active_clock") x$accrued[1] <- .1
    if (kind %in% c("risk_count", "from_count", "to_count", "absorbing_count")) {
      x$jump_times <- numeric()
      x$kernel_factors$base_increments <- matrix(numeric(), 0, 1)
      field <- switch(kind, risk_count = "relative_risk", from_count = "from",
        to_count = "to", absorbing_count = "absorbing")
      x$kernel_factors[[field]] <- rep(1, 3)
    }
    if (kind == "underflow") {
      x$kernel_factors$base_increments[1, 1] <- .Machine$double.xmin * .Machine$double.eps
      x$kernel_factors$relative_risk <- .1
    }
    expect_error(monotone_query_native(x))
    expect_identical(monotone_query_native(original), expected)
  }
  x <- original
  x$query_times <- x$query_state <- x$accrued <- rep(0, 1000)
  x$query_state[] <- 1
  x$q <- seq(0, 2, length.out = 1000)
  expect_identical(dim(monotone_query_native(x)$probability), c(1000L, 1000L))
  x$query_times <- c(x$query_times, 0)
  x$query_state <- c(x$query_state, 1)
  x$accrued <- c(x$accrued, 0)
  expect_error(monotone_query_native(x), "cell cap")
})

test_that("unrepresentable clock scores retain the ordinary grid's outcome", {
  base <- monotone_query_fixture()
  base$jump_times <- numeric()
  base$kernel_factors$base_increments <- matrix(numeric(), 0, 1)
  base$query_times <- 0
  base$query_state <- 1
  base$accrued <- base$center <- matrix(0, 1, 1)
  base$scale <- matrix(1, 1, 1)
  base$q <- c(0, 1, 2, 3)
  valid <- monotone_query_native(monotone_query_spec(base))
  for (kind in c("active_endpoint", "inactive_accrued", "inactive_scale")) {
    a <- base
    if (kind == "active_endpoint") {
      a$tau <- .Machine$double.xmax
      a$center[] <- -.Machine$double.xmax
      a$scale[] <- .Machine$double.xmax
    } else {
      a$query_state <- 2
      if (kind == "inactive_accrued") {
        a$accrued[] <- -.Machine$double.xmax
        a$center[] <- a$scale[] <- .Machine$double.xmax
      } else {
        a$accrued[] <- .25
        a$scale[] <- .Machine$double.xmin * .Machine$double.eps
      }
    }
    got <- tryCatch(do.call(markov_marginal_scalar_coverage, a), error = identity)
    expected <- tryCatch(do.call(markov_conditional_coverage_grid, a), error = identity)
    expect_s3_class(got, "error")
    expect_identical(conditionMessage(got), conditionMessage(expected))
    expect_error(monotone_query_native(monotone_query_spec(a)), "scores overflow")
    expect_identical(monotone_query_native(monotone_query_spec(base)), valid)
  }
  # A global support bound can overflow while the actual deterministic target
  # has finite score zero. Fallback must return that ordinary result, not error.
  base$center[] <- base$tau
  base$scale[] <- .Machine$double.xmin * .Machine$double.eps
  got <- do.call(markov_marginal_scalar_coverage, base)
  expect_identical(got, do.call(markov_conditional_coverage_grid, base))
  expect_identical(got$probability, matrix(1, 1, 4))
  expect_false(isTRUE(got$diagnostics$monotone_query_applied))
})

test_that("plain marginal uses only the eligible clock and retains joint signed reduction", {
  a <- monotone_query_fixture()
  a$rewards <- cbind(c(1, 0, 0), c(1, 1, 0))
  a$accrued <- cbind(c(0, .1, .25, .5), a$accrued)
  a$center <- cbind(rep(.3, 4), a$center)
  a$scale <- cbind(rep(1, 4), a$scale)
  a$kernel_factors <- list(base_increments = matrix(c(.2, .3, .1, .2, .4, .1, .2, .1), 2, 4),
    relative_risk = rep(1, 4), from = c(1, 1, 2, 2), to = c(2, 3, 1, 3), absorbing = 3)
  value <- do.call(markov_conditional_coverage_marginal, a)
  details <- value$diagnostics$marginal_diagnostics
  expect_false(isTRUE(details[[1]]$monotone_query_applied))
  expect_true(details[[2]]$monotone_query_applied)
  expect_identical(details[[2]]$grid_size, NA_integer_)
  expect_identical(value$diagnostics$grid_size, 16L)
  a$query_known <- matrix(NA_real_, 4, length(a$q))
  a$query_known[2, ] <- rep(c(0, 1), length.out = length(a$q))
  value <- do.call(markov_conditional_coverage_marginal, a)
  expect_identical(value$lower[2, ], a$query_known[2, ])
  expect_identical(value$upper[2, ], a$query_known[2, ])
  weights <- c(-2^53, 1, 2^53, -1)
  got <- do.call(markov_conditional_coverage_marginal, c(a, list(query_weights = weights)))
  oracle <- markov_grid_pass_native(numeric(), list(), matrix(0, 1, 2), 1,
    a$query_times, rep(1, 4), matrix(0, 4, 2), matrix(0, 4, 2), matrix(1, 4, 2),
    a$q, "right", c(1, 1), c(0, 0), 100L, 32L, query_weights = weights,
    query_transform = function(index, ignored)
      lapply(value[c("probability", "lower", "upper")], function(z) unname(z[index, ])))
  expect_identical(got$reduced, oracle$reduced)
})
