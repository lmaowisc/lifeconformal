marginal_threshold_fixture <- function() {
  list(jump_times = numeric(), kernels = list(),
    rewards = matrix(c(1, 0, 0, 1), 2, 2), tau = 1,
    query_times = c(0, .5), query_state = c(1L, 1L),
    accrued = matrix(c(0, .5, 0, 0), 2, 2), center = matrix(0, 2, 2),
    scale = matrix(1, 2, 2), q = c(0, .5, 1), side = "right",
    grid_size = 16L, max_nodes = 1000L, exact_limit = 32L,
    point_mass_threshold = .1, tolerance = .01, strict = FALSE,
    max_refinements = 0L, query_weights = c(-.5, 1.5))
}

test_that("marginal bounds preserve native signed arithmetic and inclusive atoms", {
  a <- marginal_threshold_fixture()
  value <- do.call(markov_conditional_coverage_marginal, a)
  joint <- do.call(markov_conditional_coverage_grid, a)
  expect_identical(value$reduced, joint$reduced)
  expect_true(value$diagnostics$marginal_applied)
  expect_identical(value$diagnostics$marginal_scalar_attempts, 2L)
  a$tau <- 1 + .Machine$double.eps
  a$query_times <- rep(.Machine$double.eps / 2, 2)
  a$accrued[, 1] <- a$query_times
  a$center[, 1] <- a$tau
  a$q <- c(0, .1, 1)
  expect_identical(do.call(markov_conditional_coverage_marginal, a)$reduced,
    do.call(markov_conditional_coverage_grid, a)$reduced)
  a$tolerance <- NULL
  expect_length(do.call(markov_conditional_coverage_marginal, a)$diagnostics$tolerance_met, 1)
})

test_that("marginal guards retain ordinary validation and capped streaming", {
  a <- marginal_threshold_fixture()
  original <- do.call(markov_conditional_coverage_grid, a)
  expect_identical(do.call(markov_conditional_coverage_marginal,
    c(a, list(max_cells = 1L))), original)
  for (name in c("strict", "max_refinements", "side", "exact_limit",
                 "point_mass_threshold", "grid_size")) {
    changed <- a
    changed[[name]] <- switch(name, strict = TRUE, max_refinements = 1L,
      side = "left", exact_limit = 256L, point_mass_threshold = .001, grid_size = 32L)
    expect_identical(do.call(markov_conditional_coverage_marginal, changed),
      do.call(markov_conditional_coverage_grid, changed))
  }
  a$query_known <- matrix(.5, 2, 3)
  ordinary <- tryCatch(do.call(markov_conditional_coverage_grid, a), error = identity)
  marginal <- tryCatch(do.call(markov_conditional_coverage_marginal, a), error = identity)
  expect_s3_class(marginal, "error")
  expect_identical(conditionMessage(marginal), conditionMessage(ordinary))
  a <- marginal_threshold_fixture()
  a$rewards <- a$rewards + 0i
  expect_error(do.call(markov_conditional_coverage_marginal, a))
})

test_that("joint known memberships precede signed marginal reduction", {
  a <- marginal_threshold_fixture()
  a$query_known <- matrix(c(0, 1, NA, NA, 1, 0), 2, 3)
  weighted <- do.call(markov_conditional_coverage_marginal, a)
  weights <- a$query_weights
  a$query_weights <- NULL
  matrix_value <- do.call(markov_conditional_coverage_marginal, a)
  exact <- !is.na(a$query_known)
  expect_identical(matrix_value$probability[exact], a$query_known[exact])
  expect_identical(matrix_value$lower[exact], a$query_known[exact])
  expect_identical(matrix_value$upper[exact], a$query_known[exact])
  oracle <- markov_grid_pass_native(numeric(), list(), matrix(0, 1, 2), 1,
    a$query_times, c(1, 1), matrix(0, 2, 2), matrix(0, 2, 2), matrix(1, 2, 2),
    a$q, "right", c(1, 1), c(0, 0), 100L, 32L, query_weights = weights,
    query_transform = function(index, ignored)
      lapply(matrix_value[c("probability", "lower", "upper")], function(z) unname(z[index, ])))
  expect_identical(weighted$reduced, oracle$reduced)
})

test_that("native marginal reduction retains the default-floor rank boundary", {
  gr <- c(1 - 2^-16, 2^-20)
  gl <- c(1, gr[1])
  weights <- c(-(gl - gr) / gl / gr, 1 / gr[2])
  probability <- matrix(1, 3, 1)
  result <- .Call("life_markov_marginal_reduce", list(probability = probability,
    lower = probability, upper = probability, weights = weights, order = c(3L, 2L, 1L)),
    PACKAGE = "lifeconformal")
  oracle <- markov_grid_pass_native(numeric(), list(), matrix(0, 1, 1), 1,
    c(0, .25, .5), c(1, 1, 1), matrix(0, 3, 1), matrix(0, 3, 1), matrix(1, 3, 1),
    0, "right", 1, 0, 100L, 32L, query_weights = weights)
  expect_identical(result, oracle$reduced)
  expect_lt(result$probability, 1)
  expect_identical(as.numeric(colSums(probability * weights)), 1)
  bad <- list(probability = probability, lower = probability, upper = probability,
    weights = weights, order = c(1L, 1L, 3L))
  expect_error(.Call("life_markov_marginal_reduce", bad, PACKAGE = "lifeconformal"), "permutation")
  big <- matrix(numeric(), 0L, 1000001L)
  bad <- list(probability = big, lower = big, upper = big, weights = numeric(), order = integer())
  expect_error(.Call("life_markov_marginal_reduce", bad, PACKAGE = "lifeconformal"), "budget")
  expect_identical(.Call("life_markov_marginal_reduce", list(probability = probability,
    lower = probability, upper = probability, weights = weights, order = c(3L, 2L, 1L)),
    PACKAGE = "lifeconformal"), result)
})

test_that("scalar budget fallback retains attempted work and other errors propagate", {
  a <- marginal_threshold_fixture()
  original <- markov_conditional_coverage_grid
  expected <- do.call(original, a)
  local_mocked_bindings(markov_conditional_coverage_grid = function(...) {
    args <- list(...)
    if (ncol(args$rewards) == 1L) stop_markov_grid_budget(args$max_nodes)
    do.call(original, args)
  }, .package = "lifeconformal")
  value <- do.call(markov_conditional_coverage_marginal, a)
  expect_identical(value$reduced, expected$reduced)
  expect_identical(value$diagnostics$marginal_scalar_attempts, 1L)
  expect_s3_class(value$diagnostics$marginal_budget_fallback$error, "lifeconformal_grid_budget")
  local_mocked_bindings(markov_conditional_coverage_grid = function(...)
    stop("positive mass underflow fixture"), .package = "lifeconformal")
  expect_error(do.call(markov_conditional_coverage_marginal, a), "positive mass underflow fixture")
})

run_marginal_stage_fixture <- function(values, max_refinements = 0L) {
  seen <- character()
  local_mocked_bindings(calibration_markov_jumps = function(...) list(),
    evaluate_threshold_profile = function(fit, profile, x, center, scale, q,
      controls, prepared, subject_index, use_marginal = FALSE) {
      stage <- if (use_marginal) "marginal" else paste0("joint", controls$grid_size)
      seen <<- c(seen, stage)
      value <- values[[stage]]
      if (inherits(value, "condition")) stop(value)
      list(probability = value$lower + (value$upper - value$lower) / 2,
        lower = value$lower, upper = value$upper,
        bound_width = max(value$upper - value$lower), queried_snapshots = 1L,
        structural_overrides = 0L, diagnostics = if (!is.null(value$diagnostics))
          value$diagnostics else list(marginal_applied = use_marginal,
            marginal_scalar_attempts = if (use_marginal) 2L else 0L))
    }, .package = "lifeconformal")
  fit <- list(augmentation = list(type = "markov_cox"),
    problem = list(rewards = rbind(c(0, 1), c(1, 1))))
  queries <- list(time = .5, snapshot = list(state = 1L, accrued = matrix(0, 1, 2)),
    weight = 1, endpoint_index = 1L, prepared = list(), diagnostics = list())
  cdt <- list(time = matrix(Inf, 1, 3), resolved = matrix(FALSE, 1, 3),
    membership = matrix(NA_real_, 1, 3))
  value <- tryCatch(evaluate_threshold_augmented_curve(fit,
    list(subjects = data.frame(id = "stage")), matrix(0, 1, 0),
    list(mean = matrix(0, 1, 2), scale = matrix(1, 1, 2)), cdt, 0:2,
    list(grid_size = 16L, max_nodes = 1000L, probability_tolerance = .01),
    queries, rank = 1, max_refinements = max_refinements), error = identity)
  list(value = value, stages = seen)
}

test_that("same-grid fallback precedes doublings and reports its own stage", {
  marginal <- list(lower = c(0, 0, 1), upper = c(1, 1, 1))
  exact <- list(lower = c(0, 1, 1), upper = c(0, 1, 1))
  x <- run_marginal_stage_fixture(list(marginal = marginal, joint16 = exact))
  expect_identical(x$stages, c("marginal", "joint16"))
  d <- x$value$diagnostics
  expect_identical(d$profile_refinements, 0L)
  expect_identical(d$marginal_fallback_attempts, 1L)
  expect_identical(d$marginal_fallbacks, 1L)
  expect_identical(d$ordinary_evaluations, 1L)
  expect_identical(d$profile_stages, "joint")
  expect_equal(d$profile_grid_sizes, 16)
  expect_equal(d$possible_index, 2)
  x <- run_marginal_stage_fixture(list(marginal = exact))
  expect_identical(x$stages, "marginal")
  expect_identical(x$value$diagnostics$profile_stages, "marginal")
  expect_true(is.na(x$value$diagnostics$profile_grid_sizes))
  expect_equal(x$value$diagnostics$profile_marginal_grid_sizes, 16)
})

test_that("ambiguity retains marginal evidence after a failed fallback", {
  marginal <- list(lower = c(0, 0, 1), upper = c(1, 1, 1))
  budget <- structure(list(message = "node budget fixture", call = NULL, diagnostics = list()),
    class = c("lifeconformal_grid_budget", "error", "condition"))
  x <- run_marginal_stage_fixture(list(marginal = marginal, joint16 = budget))
  expect_s3_class(x$value, "lifeconformal_threshold_ambiguity")
  d <- x$value$diagnostics
  expect_true(d$marginal_active)
  expect_true(is.na(d$grid_size))
  expect_identical(d$profile_stages, "marginal")
  expect_identical(d$marginal_fallback_attempts, 1L)
  expect_identical(d$marginal_fallbacks, 0L)
  expect_length(d$budget_errors, 1)
  expect_identical(d$marginal_scalar_attempts, 2L)
})

test_that("nonnested passes recompute every candidate and exact monotone stages stay distinct", {
  values <- list(marginal = list(lower = c(0, 0, 1), upper = c(1, 1, 1)),
    joint16 = list(lower = c(.1, .4, 1), upper = c(.9, 1, 1)),
    joint32 = list(lower = c(0, .5, 1), upper = c(1, .9, 1)))
  x <- run_marginal_stage_fixture(values, max_refinements = 1L)
  expect_identical(x$stages, c("marginal", "joint16", "joint32"))
  expect_s3_class(x$value, "lifeconformal_threshold_ambiguity")
  expect_equal(x$value$diagnostics$possible_index, 1)
  expect_equal(x$value$diagnostics$certain_index, 3)
  exact <- list(lower = c(0, 1, 1), upper = c(0, 1, 1),
    diagnostics = list(method = "deterministic_exact_monotone_clocks"))
  x <- run_marginal_stage_fixture(list(marginal = exact))
  expect_identical(x$value$diagnostics$profile_stages, "monotone")
  expect_true(is.na(x$value$diagnostics$profile_grid_sizes))
  expect_false(x$value$diagnostics$marginal_active)
  expect_identical(x$value$diagnostics$ordinary_evaluations, 1L)
})
