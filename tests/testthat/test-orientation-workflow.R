orientation_workflow_data <- function(prefix, n = 48L) {
  i <- seq_len(n)
  data.frame(id = paste0(prefix, i), time = ((i * 7L) %% 19L + 1L) / 3,
               status = ifelse(i %% 4L == 0L, 0L, 2L))
}

test_that("increasing and mixed targets use exact public working coverage", {
  training <- orientation_workflow_data("train")
  fit <- life_train(LifeHist(time, status) ~ 1, data = training,
    absorbing = 2, targets = list(alive = 0, dead = 2), tau = 5)
  baseline <- data.frame(id = "new")
  moments <- predict(fit, baseline)
  expect_equal(sum(moments$mean), 5, tolerance = 1e-13)
  x <- apply_preprocessor(baseline, fit$preprocessor, "outcome")$x
  jumps <- predict_transition_jumps(fit$outcome, x)
  kernels <- lapply(seq_along(jumps$jump_times), function(j) {
    increment <- matrix(0, 2, 2)
    increment[cbind(jumps$edges$from, jumps$edges$to)] <-
      jumps$base_increments[j, ] * jumps$relative_risk[1, ]
    make_markov_kernel(increment, absorbing = fit$problem$absorbing)$kernel
  })
  # In this two-state model an exit at u gives exactly (u, tau-u).
  # Calculate its probabilities independently from the product survival.
  # Repeatedly adding many sojourn lengths is not a canonical atom-boundary
  # reference: floating reassociation can move an atom by one ulp.
  hazards <- vapply(kernels, function(kernel) kernel[1, 2], numeric(1))
  survival <- cumprod(c(1, 1 - hazards))
  atoms <- list(values = cbind(alive = c(jumps$jump_times, 5),
                               dead = 5 - c(jumps$jump_times, 5)),
                probability = c(survival[seq_along(hazards)] * hazards,
                                tail(survival, 1)))
  scores <- predict(fit, data.frame(row.names = seq_len(nrow(atoms$values))),
                     type = "score", values = atoms$values)
  q <- sort(unique(c(0, scores, 3)))
  truth <- vapply(q, function(x) sum(atoms$probability[scores <= x]), numeric(1))
  result <- predict(fit, baseline, type = "coverage", q = q)
  expect_equal(as.vector(result$probability), truth, tolerance = 1e-13)
  expect_identical(result$numerical$method, "deterministic_native_monotone_batch")
  expect_identical(result$numerical$complemented_targets, "dead")
  expect_equal(result$lower, result$upper)
  # An initially absorbing subject already accumulates the increasing target.
  absorbed <- predict(fit, baseline, initial_state = 2, type = "coverage", q = c(0, 1))
  expect_equal(as.vector(absorbed$probability), c(1, 1))
  constant <- life_train(LifeHist(time, status) ~ 1, data = training,
    absorbing = 2, targets = list(total = c(0, 2)), tau = 5)
  expect_identical(as.vector(predict(constant, baseline)$mean), 5)
  expect_equal(as.vector(predict(constant, baseline, type = "coverage", q = 0)$probability), 1)
})

test_that("public increasing-target AIPW matches direct original-unit conditional atoms", {
  training <- orientation_workflow_data("train")
  calibration <- orientation_workflow_data("cal", 36L)
  calibration$time <- calibration$time * .8
  fit <- life_train(LifeHist(time, status) ~ 1, data = training,
    absorbing = 2, targets = list(dead = 2), tau = 5)
  q <- c(0, .5, 1, 1.5, 2, 3)
  result <- life_calibrate(fit, data = calibration, alpha = .2, q = q,
                            control = list(keep_contributions = TRUE))
  expect_identical(result$diagnostics$numerical$method, "exact_batched_monotone_rewards")
  expect_identical(result$diagnostics$numerical$complemented_targets, "dead")
  # Direct original-unit atom integration checks the public dispatch and its
  # signed common-endpoint weights without using a reflected clock evaluator.
  input <- resolve_input(data = calibration, schema = fit$schema)
  history <- normalize_histories(input, absorbing = fit$problem$absorbing_labels,
    states = fit$problem$states)
  view <- apply_preprocessor(history$baseline, fit$preprocessor, "outcome")
  censor_view <- apply_preprocessor(history$baseline, fit$preprocessor, "censor")
  moments <- predict_life_moments(fit, history$baseline, initial_state = 0)
  cdt <- compute_cdt(history, fit$problem, moments$mean, moments$scale, q)
  queries <- prepare_calibration_queries(fit, history, censor_view$x, cdt)
  jumps <- predict_transition_jumps(fit$augmentation, view$x[1, , drop = FALSE])
  kernels <- lapply(seq_along(jumps$jump_times), function(j) {
    increment <- matrix(0, 2, 2)
    increment[cbind(jumps$edges$from, jumps$edges$to)] <- jumps$base_increments[j, ]
    make_markov_kernel(increment, absorbing = fit$problem$absorbing)$kernel
  })
  direct <- matrix(0, nrow(history$subjects), length(q))
  for (a in seq_along(queries$time)) {
    i <- queries$subject_index[a]
    atoms <- propagate_markov_rewards(jumps$jump_times, kernels, fit$problem$rewards,
      5, start = queries$time[a], state = queries$snapshot$state[a],
      accrued = as.vector(queries$snapshot$accrued[a, ]))$distributions[[1]]
    scores <- abs((atoms$values[, 1] - moments$mean[i, 1]) / moments$scale[i, 1])
    probability <- vapply(q, function(x) sum(atoms$probability[scores <= x]), numeric(1))
    known <- cdt$resolved[i, ] & queries$time[a] >= cdt$time[i, ]
    probability[known] <- as.numeric(cdt$membership[i, known])
    direct[i, ] <- direct[i, ] + queries$weight[a] * probability
  }
  origin <- cdt$resolved & cdt$time == 0
  direct[origin] <- as.numeric(cdt$membership[origin])
  expect_equal(result$contributions, direct, tolerance = 2e-12)
  expect_equal(result$curve$estimate, colMeans(direct), tolerance = 2e-12)
  expect_equal(result$curve$numerical_lower, result$curve$numerical_upper)
  regions <- predict(result, data.frame(id = c("a", "b")))
  expect_identical(regions$bounds$target, c("dead", "dead"))
  expect_true(all(regions$bounds$lower >= 0 & regions$bounds$upper <= 5))
})
