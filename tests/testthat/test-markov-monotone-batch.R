native_monotone_fixture <- function() {
  list(jump_times = c(.5, 1, 2, 3),
       base_increments = rbind(c(.15, .05, .2), c(.2, .1, .15),
                               c(.1, .05, .1), c(.15, .1, .1)),
       edges = data.frame(from = c(1, 1, 2), to = c(2, 3, 3)),
       relative_risk = rbind(c(1, 1, 1), c(.5, 2, 1.5), c(3, 4, 2)),
       rewards = rbind(c(1, 1), c(0, 1), c(0, 0)), tau = 4,
       center = rbind(c(2, 3), c(1, 2.5), c(2.5, 3)),
       scale = rbind(c(.7, 1), c(1, .8), c(.8, .5)),
       q = seq(0, 5, .125))
}

native_monotone_reference <- function(z, subject_index, query_times,
                                        query_state, accrued, weight) {
  out <- matrix(0, nrow(z$relative_risk), length(z$q))
  for (i in unique(subject_index)) {
    indices <- which(subject_index == i)
    kernels <- lapply(seq_along(z$jump_times), function(j) {
      increments <- matrix(0, nrow(z$rewards), nrow(z$rewards))
      increments[as.matrix(z$edges)] <- z$base_increments[j, ] * z$relative_risk[i, ]
      make_markov_kernel(increments)$kernel
    })
    fit <- markov_conditional_coverage_monotone(z$jump_times, kernels, z$rewards,
      z$tau, query_times[indices], query_state[indices], accrued[indices, , drop = FALSE],
      z$center[rep(i, length(indices)), , drop = FALSE],
      z$scale[rep(i, length(indices)), , drop = FALSE], z$q)
    out[i, ] <- colSums(fit$probability * weight[indices])
  }
  out
}

test_that("native batched monotone contribution equals independent query reference", {
  z <- native_monotone_fixture()
  query <- list(subject_index = c(2, 1, 3, 1, 2, 3, 1),
                query_times = c(.5, 0, 1, 2, 3, 4, 1),
                query_state = c(1, 1, 2, 2, 3, 3, 1),
                accrued = rbind(c(.5, .5), c(0, 0), c(.5, 1), c(1, 2),
                                c(.5, 2), c(1, 3), c(1, 1)),
                weight = c(-.2, 1.5, -.3, -.4, 1.2, 1.3, -.1))
  got <- do.call(markov_monotone_batch, c(z, query))
  reference <- do.call(native_monotone_reference, c(list(z), query))
  expect_equal(got$contributions, reference, tolerance = 2e-13)
  expect_equal(got$diagnostics$maximum_probability_error, 0)
  # Isolating each query verifies through-time boundaries, not just a weighted
  # total in which errors could cancel.
  for (a in seq_along(query$weight)) {
    one <- lapply(query, function(x) if (is.matrix(x)) x[a, , drop = FALSE] else x[a])
    one$weight <- 1
    got <- do.call(markov_monotone_batch, c(z, one))
    ref <- do.call(native_monotone_reference, c(list(z), one))
    expect_equal(got$contributions, ref, tolerance = 2e-13)
  }
})

test_that("native certificates override exactly at determining time", {
  z <- native_monotone_fixture()
  query <- list(subject_index = c(1, 1), query_times = c(.5, 1),
                query_state = c(1, 1), accrued = rbind(c(.5, .5), c(1, 1)),
                weight = c(-.25, 1.25))
  raw_first <- do.call(markov_monotone_batch, c(z, list(subject_index = 1,
    query_times = .5, query_state = 1, accrued = matrix(c(.5, .5), 1), weight = -.25)))
  time <- matrix(1, 3, length(z$q))
  resolved <- matrix(FALSE, 3, length(z$q)); resolved[1, ] <- TRUE
  membership <- matrix(NA, 3, length(z$q)); membership[1, ] <- z$q >= 2
  got <- do.call(markov_monotone_batch, c(z, query,
    list(determination_time = time, resolved = resolved, membership = membership)))
  expect_equal(got$contributions[1, ], raw_first$contributions[1, ] + 1.25 * (z$q >= 2))
})

test_that("native monotone clocks include reward-equivalent recurrent transitions", {
  z <- native_monotone_fixture()
  z$rewards <- matrix(c(1, 1, 0), 3)
  z$edges <- rbind(z$edges, data.frame(from = 2, to = 1))
  z$base_increments <- cbind(z$base_increments, .1)
  z$relative_risk <- cbind(z$relative_risk, c(1, 2, 3))
  z$center <- matrix(c(2, 3, 2.5), 3)
  z$scale <- matrix(1, 3, 1)
  query <- list(subject_index = c(1, 2, 3), query_times = c(0, 1, 2),
                query_state = c(1, 2, 1), accrued = matrix(c(0, 1, 2), 3),
                weight = c(1, -.5, 2))
  expect_equal(do.call(markov_monotone_batch, c(z, query))$contributions,
               do.call(native_monotone_reference, c(list(z), query)), tolerance = 2e-13)
})

test_that("native row adjustment and no-jump or empty-query cases are valid", {
  z <- native_monotone_fixture()
  z$base_increments <- z$base_increments * 10
  query <- list(subject_index = 3, query_times = 0, query_state = 1,
                accrued = matrix(c(0, 0), 1), weight = 1)
  got <- do.call(markov_monotone_batch, c(z, query))
  expect_gt(got$diagnostics$adjusted_rows, 0)
  expect_equal(got$contributions, do.call(native_monotone_reference, c(list(z), query)),
               tolerance = 2e-13)
  z$jump_times <- numeric(); z$base_increments <- matrix(numeric(), 0, 3)
  expect_equal(do.call(markov_monotone_batch, c(z, query))$contributions,
               do.call(native_monotone_reference, c(list(z), query)))
  query <- list(subject_index = numeric(), query_times = numeric(),
                query_state = numeric(), accrued = matrix(numeric(), 0, 2), weight = numeric())
  expect_equal(do.call(markov_monotone_batch, c(z, query))$contributions,
               matrix(0, 3, length(z$q)))
})

test_that("native entry point and wrapper reject malformed inputs safely", {
  expect_error(.Call("life_markov_monotone_batch", list(), PACKAGE = "lifeconformal"),
               "named list")
  z <- native_monotone_fixture()
  query <- list(subject_index = 1, query_times = 1, query_state = 1,
                accrued = matrix(c(.5, 1), 1), weight = 1)
  expect_error(do.call(markov_monotone_batch, c(z, query)),
               class = "lifeconformal_monotone_history")
  query$accrued[] <- 1
  query$subject_index <- 4
  expect_error(do.call(markov_monotone_batch, c(z, query)), "Invalid subject")
  query$subject_index <- 1
  z$rewards[2, 1] <- .5
  expect_error(do.call(markov_monotone_batch, c(z, query)), "binary")
})

test_that("compensated signed accumulation retains small residuals", {
  z <- native_monotone_fixture()
  z$jump_times <- numeric(); z$base_increments <- matrix(numeric(), 0, 3)
  z$q <- c(10, 20)
  got <- do.call(markov_monotone_batch, c(z, list(subject_index = c(1, 1, 1, 2),
    query_times = c(3, 2, 1, 0), query_state = c(1, 1, 1, 1),
    accrued = cbind(c(3, 2, 1, 0), c(3, 2, 1, 0)),
    weight = c(1e16, 1, -1e16, 0))))
  expect_equal(got$contributions[1, ], c(1, 1))
  expect_equal(got$contributions[2:3, ], matrix(0, 2, 2))
})
