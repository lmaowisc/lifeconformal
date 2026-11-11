test_that("IPW uses observed determination and never multiplies missing values", {
  value <- compute_ipw_contributions(
    c(a = 1, b = NA, c = 0, d = 1),
    c(TRUE, FALSE, TRUE, TRUE),
    c(0.5, 0.2, 0.4, 1)
  )
  expect_identical(value, c(a = 2, b = 0, c = 0, d = 1))
  expect_identical(
    compute_ipw_contributions(c(NA, NA), c(FALSE, FALSE), c(1, 0.8)),
    c(0, 0)
  )
  expect_identical(
    names(formals(compute_ipw_contributions)),
    c("observed_membership", "resolved", "g_left")
  )
})

test_that("IPW input validation prevents recycling and unobserved membership", {
  expect_error(compute_ipw_contributions(c(1, 0), c(TRUE, TRUE), 1), "length 2")
  expect_error(compute_ipw_contributions(NA, TRUE, 1), "observed for resolved")
  expect_error(compute_ipw_contributions(0.5, TRUE, 1), "zero or one")
  expect_error(compute_ipw_contributions(1, 1, 1), "logical vector")
  expect_error(compute_ipw_contributions(1, NA, 1), "logical vector")
  expect_error(compute_ipw_contributions(1, TRUE, 0), "finite inverse")
  expect_error(compute_ipw_contributions(1, TRUE, 1.1), "finite inverse")
  expect_error(compute_ipw_contributions(1, TRUE, 1e-320), "finite inverse")
  expect_error(compute_ipw_contributions(1, TRUE, NA_real_), "finite numeric")
  expect_error(compute_ipw_contributions(1, TRUE, 1 + 0i), "finite numeric")
  expect_error(compute_ipw_contributions(matrix(1), TRUE, 1), "aligned")
  expect_error(compute_ipw_contributions(numeric(), logical(), numeric()),
               "nonempty")
})

test_that("step AIPW agrees with a hand calculation and retains signed values", {
  jumps <- data.frame(
    subject_index = c(1, 1, 2, 3),
    time = c(1, 2, 1, 1),
    g_left = c(1, 0.8, 1, 1),
    g_right = c(0.8, 0.5, 0.8, 0.6),
    m_plus = c(0.3, 0.7, 0.4, 0.9)
  )
  result <- compute_aipw_contributions(
    c(1, NA, 0), c(TRUE, FALSE, TRUE), c(3, 2, 4),
    c(0.5, 0.8, 0.6), c(NA, 0.25, NA), jumps
  )
  expect_equal(result, c(1.4, 0.2125, -0.6), tolerance = 1e-14)
  expect_gt(result[1], 1)
  expect_lt(result[3], 0)
  expect_equal(
    compute_aipw_contributions(
      c(1, NA, 0), c(TRUE, FALSE, TRUE), c(3, 2, 4),
      c(0.5, 0.8, 0.6), c(NA, 0.25, NA), jumps[c(4, 2, 3, 1), ]
    ), result
  )
})

test_that("endpoint jumps are excluded before nuisance evaluation", {
  jumps <- data.frame(
    subject_index = c(1, 1, 1), time = c(1, 2, 3),
    g_left = c(1, NA, NA), g_right = c(0.8, NA, NA),
    m_plus = c(0.4, NA, NA)
  )
  expect_equal(
    compute_aipw_contributions(NA_real_, FALSE, 2, 0.8, 0.25, jumps),
    0.2125
  )
  # Censoring at the prediction origin contains no strictly earlier jump.
  expect_equal(
    compute_aipw_contributions(
      NA_real_, FALSE, 0, 1, 0.25,
      transform(jumps[1, ], time = 0, g_left = NA, g_right = NA, m_plus = NA)
    ),
    0.25
  )
  # A jump at zero does enter the sum for a later stopped time.
  expect_equal(
    compute_aipw_contributions(
      NA_real_, FALSE, 1, 0.8, 0.25, transform(jumps[1, ], time = 0)
    ),
    0.2125
  )
})

test_that("full-support AIPW telescopes under coherent survival", {
  jumps <- data.frame(
    subject_index = c(1, 1, 2), time = c(1, 2, 0),
    g_left = c(1, 0.8, 1), g_right = c(0.8, 0.5, 0.7),
    m_plus = c(1, 1, 1)
  )
  expect_equal(
    compute_aipw_contributions(
      c(1, NA), c(TRUE, FALSE), c(3, 1), c(0.5, 0.7), c(NA, 1), jumps
    ), c(1, 1), tolerance = 1e-14
  )
})

test_that("stopped and common-endpoint contributions agree after determination", {
  stopped <- compute_aipw_contributions(
    1, TRUE, 2, 0.8, NA_real_,
    data.frame(subject_index = 1, time = 1, g_left = 1,
               g_right = 0.8, m_plus = 0.4)
  )
  extended <- compute_aipw_contributions(
    1, TRUE, 4, 0.4, NA_real_,
    data.frame(subject_index = c(1, 1, 1), time = c(1, 2, 3),
               g_left = c(1, 0.8, 0.6), g_right = c(0.8, 0.6, 0.4),
               m_plus = c(0.4, 1, 1))
  )
  expect_equal(stopped, extended, tolerance = 1e-14)
})

test_that("complete data give observed membership for both kernels", {
  membership <- c(1, 0, 1, 1, 0)
  resolved <- rep(TRUE, length(membership))
  g <- rep(1, length(membership))
  expect_identical(compute_ipw_contributions(membership, resolved, g), membership)
  expect_identical(
    compute_aipw_contributions(membership, resolved, 1:5, g,
                              rep(NA_real_, length(membership))),
    membership
  )
})

test_that("exact finite-support expectation checks isolate nuisance conditions", {
  # Independent censoring at time 1; binary outcome revealed at time 2.
  # Enumerate outcome/censoring combinations instead of Monte Carlo sampling.
  p <- 0.35
  censor_probability <- 0.4
  probability <- c(p * censor_probability, (1 - p) * censor_probability,
                   p * (1 - censor_probability),
                   (1 - p) * (1 - censor_probability))
  mean_contribution <- function(working_censor_probability, working_outcome) {
    survival <- 1 - working_censor_probability
    membership <- c(NA, NA, 1, 0)
    resolved <- c(FALSE, FALSE, TRUE, TRUE)
    g <- c(1, 1, survival, survival)
    jumps <- data.frame(
      subject_index = c(3, 4), time = c(1, 1), g_left = c(1, 1),
      g_right = c(survival, survival),
      m_plus = rep(working_outcome, 2)
    )
    c(
      ipw = sum(probability * compute_ipw_contributions(membership, resolved, g)),
      aipw = sum(probability * compute_aipw_contributions(
        membership, resolved, c(1, 1, 2, 2), g,
        c(working_outcome, working_outcome, NA, NA), jumps
      ))
    )
  }
  expect_equal(unname(mean_contribution(censor_probability, 0.8)), c(p, p))
  wrong_censor <- mean_contribution(0.2, p)
  expect_equal(unname(wrong_censor["aipw"]), p)
  expect_equal(unname(wrong_censor["ipw"]), p * 0.6 / 0.8)
  expect_equal(unname(mean_contribution(0.2, 0.8)["aipw"]), 0.4625)
})

test_that("AIPW validates jump shape, bounds and conditional probabilities", {
  j <- data.frame(subject_index = 1, time = 1, g_left = 1,
                  g_right = 0.8, m_plus = 0.4)
  eval_jumps <- function(jumps = j) {
    compute_aipw_contributions(NA_real_, FALSE, 2, 0.8, 0.3, jumps)
  }
  expect_error(eval_jumps(j[-1]), "data frame")
  expect_error(eval_jumps(transform(j, subject_index = 2)), "index the endpoint")
  expect_error(eval_jumps(transform(j, subject_index = 1.5)), "index the endpoint")
  expect_error(eval_jumps(transform(j, time = -1)), "nonnegative")
  expect_error(eval_jumps(transform(j, g_left = 0.7)), "must not increase")
  expect_error(eval_jumps(transform(j, g_right = 0)), "finite inverse")
  expect_error(eval_jumps(transform(j, m_plus = 1.1)), "in \\[0, 1\\]")
  expect_error(eval_jumps(transform(j, m_plus = NA_real_)), "finite numeric")
  expect_error(eval_jumps(rbind(j, j)), "duplicate subject/time")
  expect_error(
    compute_aipw_contributions(NA_real_, FALSE, 2, 0.8, NA_real_, j),
    "required for every unresolved"
  )
  expect_error(compute_aipw_contributions(1, TRUE, -1, 1, NA_real_), "nonnegative")
})

test_that("inversion takes the first crossing of the raw nonmonotone curve", {
  estimate <- c(0.4, 0.9, 0.7, 1.1)
  result <- invert_coverage_curve(1:4, estimate, n = 9, alpha = 0.2)
  expect_identical(result$threshold, 2L)
  expect_identical(result$index, 2L)
  expect_false(result$full_support)
  expect_identical(result$reason, "crossing")
  expect_equal(result$rank, 8)
  expect_equal(result$criterion, 8 / 9)
  expect_identical(result$curve$estimate, estimate)
})

test_that("inversion handles equality, noncrossing and the boundary rank", {
  result <- invert_coverage_curve(0:2, c(0.2, 8 / 9, 0.1), 9, 0.2)
  expect_identical(result$threshold, 1L)
  result <- invert_coverage_curve(0:2, c(0.1, 0.5, 0.7), 9, 0.2)
  expect_identical(result$threshold, Inf)
  expect_true(result$full_support)
  expect_identical(result$reason, "no_crossing")
  expect_identical(result$index, NA_integer_)
  result <- invert_coverage_curve(0:2, c(0.2, 1, 2), 9, 0.05)
  expect_true(result$full_support)
  expect_equal(result$rank, 10)
  expect_identical(result$reason, "rank_exceeds_sample")
  expect_true(invert_coverage_curve(numeric(), numeric(), 9, 0.05)$full_support)
})

test_that("complete-data inversion recovers conformal order statistics and ties", {
  scores <- 1:99
  result <- invert_coverage_curve(scores, scores / 99, 99, 0.1)
  expect_equal(result$threshold, 90)
  scores <- c(0, 0, 1, 2)
  q <- unique(scores)
  estimate <- vapply(q, function(threshold) mean(scores <= threshold), 0)
  result <- invert_coverage_curve(q, estimate, length(scores), 0.4)
  expect_equal(result$rank, 3)
  expect_equal(result$threshold, 1)
})

test_that("inversion rejects malformed grids, sample sizes and levels", {
  expect_error(invert_coverage_curve(c(1, 1), c(0, 1), 9, 0.2), "strictly increasing")
  expect_error(invert_coverage_curve(c(2, 1), c(0, 1), 9, 0.2), "strictly increasing")
  expect_error(invert_coverage_curve(c(1, Inf), c(0, 1), 9, 0.2), "finite numeric")
  expect_error(invert_coverage_curve(1:2, NA_real_, 9, 0.2), "length 2")
  expect_error(invert_coverage_curve(1:2, c(0, NaN), 9, 0.2), "finite numeric")
  expect_error(invert_coverage_curve(numeric(), numeric(), 9, 0.2), "nonempty")
  expect_error(invert_coverage_curve(1, 1, 0, 0.2), "positive integer")
  expect_error(invert_coverage_curve(1, 1, 2.5, 0.2), "positive integer")
  expect_error(invert_coverage_curve(1, 1, 9, 0), "strictly between")
  expect_error(invert_coverage_curve(1, 1, 9, 1), "strictly between")
})
