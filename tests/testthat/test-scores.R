score_survival_fixture <- function(targets = list(alive = 0), score_weights = NULL,
                                    outcome_model = "markov_cox", outcome_options = list(),
                                    draws = 32L) {
  data <- data.frame(id = paste0("t", seq_len(12)),
                     time = rep(c(1, 2, 4, 7, 9, 12), 2), status = 1)
  life_train(LifeHist(time, status) ~ 1, data = data, absorbing = 1,
              targets = targets, tau = 10, score_weights = score_weights,
              outcome_model = outcome_model, outcome_options = outcome_options,
              control = list(draws = draws), seed = 23)
}

test_that("scores use original target units and frozen weighted scale parameters", {
  fit <- score_survival_fixture(list(first = 0, second = 0),
                                 c(second = 1, first = 0.5))
  data <- data.frame(id = c("a", "b"))
  moments <- predict_life_moments(fit, data)
  expect_identical(fit$problem$score_weights, c(first = 0.5, second = 1))
  expect_equal(moments$scale[, "first"], 0.5 * moments$scale[, "second"])
  values <- data.frame(second = c(4, 8), first = c(1, 7))
  expected <- pmax(abs((values$first - moments$mean[, "first"]) /
                         moments$scale[, "first"]),
                    abs((values$second - moments$mean[, "second"]) /
                         moments$scale[, "second"]))
  expect_equal(evaluate_life_score(fit, data, values), unname(expected))
  expect_identical(evaluate_life_score(fit, data, values),
                   evaluate_life_score(fit, data, as.matrix(values[c("first", "second")])))
  expect_equal(evaluate_life_score(fit, data[1, , drop = FALSE],
                                    c(second = 4, first = 1)), unname(expected[1]))
  expect_identical(evaluate_life_score(fit, data, moments$mean), c(0, 0))
  component <- subset_component_problem(fit, 1L)
  expect_identical(component$problem$score_weights, c(first = 0.5))
  expect_equal(predict_life_moments(component, data)$scale,
               moments$scale[, "first", drop = FALSE])
})

test_that("scalar, empty and deterministic scores preserve their dimensions", {
  fit <- score_survival_fixture()
  data <- data.frame(id = c("a", "b"))
  expect_equal(evaluate_life_score(fit, data, c(2, 4)),
               evaluate_life_score(fit, data, matrix(c(2, 4), 2L)))
  expect_length(evaluate_life_score(fit, data[1, , drop = FALSE], 2), 1)
  expect_identical(evaluate_life_score(fit, data.frame(id = character()), numeric()), numeric())
  # An initially absorbed path has target zero; its nonzero scale is the frozen
  # floor, not a newly estimated sample standard deviation of candidate values.
  expect_equal(evaluate_life_score(fit, data, c(0, 0.1), initial_state = 1), c(0, 1))
  two <- score_survival_fixture(list(alive = 0, dead = 1))
  empty <- matrix(numeric(), 0L, 2L, dimnames = list(NULL, c("dead", "alive")))
  expect_identical(evaluate_life_score(two, data.frame(id = character()), empty), numeric())
})

test_that("score inputs cannot recycle rows or guess unnamed multivariate targets", {
  fit <- score_survival_fixture(list(alive = 0, dead = 1))
  data <- data.frame(id = c("a", "b"))
  expect_error(evaluate_life_score(fit, data, matrix(1, 2, 2)), "name every trained target")
  expect_error(evaluate_life_score(fit, data, c(alive = 1, dead = 2)), "n-by-d matrix")
  expect_error(evaluate_life_score(fit, data, data.frame(alive = 1, dead = 2)), "baseline rows")
  expect_error(evaluate_life_score(fit, data, data.frame(alive = c(1, 2), wrong = c(2, 3))),
               "name every trained target")
  duplicate <- matrix(1, 2, 2, dimnames = list(NULL, c("alive", "alive")))
  expect_error(evaluate_life_score(fit, data, duplicate), "name every trained target")
  bad <- data.frame(alive = c(1, NA_real_), dead = c(2, 3))
  expect_error(evaluate_life_score(fit, data, bad), "finite numeric")
  bad$alive <- c(1, Inf)
  expect_error(evaluate_life_score(fit, data, bad), "finite numeric")
  bad$alive <- factor(c("one", "two"))
  expect_error(evaluate_life_score(fit, data, bad), "finite numeric")
  scalar <- score_survival_fixture()
  expect_error(evaluate_life_score(scalar, data, 1), "one value per baseline row")
})

test_that("baseline monotone working coverage matches empirical survival atoms at equality", {
  fit <- score_survival_fixture()
  truth <- rep(c(1, 2, 4, 7, 9, 10), 2)
  baseline <- data.frame(id = paste0("p", seq_along(truth)))
  scores <- evaluate_life_score(fit, baseline, truth)
  q <- sort(unique(c(0, scores, 100)))
  expected <- vapply(q, function(x) mean(scores <= x), numeric(1))
  newdata <- data.frame(id = c("new1", "new2"))
  result <- predict_life_coverage(fit, newdata, q)
  expect_identical(result$q, q)
  expect_identical(dim(result$probability), c(2L, length(q)))
  expect_equal(result$probability, rbind(expected, expected), ignore_attr = TRUE,
               tolerance = 1e-13)
  expect_identical(result$lower, result$probability)
  expect_identical(result$upper, result$probability)
  expect_identical(result$numerical$method, "deterministic_native_monotone_batch")
  expect_identical(result$numerical$working_model, "markov_cox")
  expect_true(all(diff(result$probability[1, ]) >= 0))
  expect_equal(result$probability[, length(q)], c(1, 1))
  absorbed <- predict_life_coverage(fit, newdata, c(0, 1), initial_state = 1)
  expect_identical(absorbed$probability, matrix(1, 2, 2))
  # Inspection is the score/outcome working CDF, even with a separate or invalid
  # augmentation object: augmentation does not define these probabilities.
  fit$augmentation <- list(type = "not_an_outcome_model")
  local_mocked_bindings(conditional_life_coverage = function(...) stop("Augmentation used"),
                         .package = "lifeconformal")
  expect_identical(predict_life_coverage(fit, newdata, q)$probability, result$probability)
})

test_that("general Markov working coverage matches exact recurrent joint atoms", {
  i <- seq_len(6)
  records <- lapply(i, function(j) {
    data.frame(id = paste0("t", j), time = if (j <= 3) c(1, 2, 4) else c(1, 3, 5),
                 status = c(1, 7, 2))
  })
  fit <- life_train(LifeHist(time, status) ~ 1, data = do.call(rbind, records),
                    absorbing = 2, status_map = c("1" = 1, "7" = 0, "2" = 2),
                    targets = list(home = 0, alive = c(0, 1)), tau = 5,
                    score_weights = c(home = 0.5, alive = 1))
  baseline <- data.frame(id = "new")
  x <- apply_preprocessor(baseline, fit$preprocessor, "outcome")$x
  jumps <- predict_transition_jumps(fit$outcome, x)
  expect_false(markov_monotone_eligible(fit$problem$rewards, jumps$edges))
  kernels <- lapply(seq_along(jumps$jump_times), function(j) {
    increment <- matrix(0, 3, 3)
    increment[cbind(jumps$edges$from, jumps$edges$to)] <- jumps$base_increments[j, ]
    make_markov_kernel(increment, absorbing = fit$problem$absorbing)$kernel
  })
  atoms <- propagate_markov_rewards(jumps$jump_times, kernels, fit$problem$rewards,
                                     5, state = match("0", as.character(fit$problem$states)))$distributions[[1L]]
  scores <- evaluate_life_score(fit, data.frame(row.names = seq_len(nrow(atoms$values))),
                                 atoms$values)
  q <- sort(unique(c(0, scores, 100)))
  expected <- vapply(q, function(x) sum(atoms$probability[scores <= x]), numeric(1))
  result <- predict_life_coverage(fit, baseline, q)
  expect_equal(as.vector(result$probability), expected, tolerance = 1e-13)
  expect_equal(result$lower, result$upper, tolerance = 1e-13)
  expect_identical(result$numerical$method, "deterministic_joint_grid_envelopes")
  expect_equal(result$numerical$maximum_probability_bound_width, 0)
})

test_that("coverage validation and empty queries retain the public output contract", {
  fit <- score_survival_fixture()
  data <- data.frame(id = "new")
  for (q in list(numeric(), c(1, 0), c(0, 0), -1, NA_real_, Inf, "1", matrix(1))) {
    expect_error(predict_life_coverage(fit, data, q), "q must be finite")
  }
  one <- predict_life_coverage(fit, data, 0)
  expect_identical(dim(one$probability), c(1L, 1L))
  empty <- predict_life_coverage(fit, data.frame(id = character()), c(0, 1))
  expect_identical(dim(empty$probability), c(0L, 2L))
  expect_identical(empty$q, c(0, 1))
  expect_identical(empty$numerical$method, "empty_prediction")
})

test_that("empty working graphs are monotone regardless of empty matrix storage type", {
  rewards <- matrix(c(1, 0), 2, 1)
  logical_edges <- matrix(logical(), 0L, 2L)
  expect_true(markov_monotone_eligible(rewards, logical_edges))
  expect_true(markov_monotone_eligible(rewards,
    data.frame(from = integer(), to = integer())))
  expect_error(markov_monotone_eligible(rewards, matrix(logical(), 0L, 3L)), "two columns")
  expect_error(markov_monotone_eligible(rewards, matrix(TRUE, 1L, 2L)), "state indices")
  train <- data.frame(id = paste0("t", seq_len(12)), time = 12, status = 0)
  fit <- life_train(LifeHist(time, status) ~ 1, data = train,
                    absorbing = 1, targets = list(alive = 0), tau = 10)
  result <- predict_life_coverage(fit, data.frame(id = "new"), c(0, 0.5, 1))
  expect_identical(result$probability, matrix(1, 1L, 3L))
  expect_identical(result$lower, result$upper)
})

test_that("public score and working-coverage interfaces agree with stored callable functions", {
  fit <- score_survival_fixture(list(first = 0, second = 0),
                                 c(first = 0.5, second = 1))
  data <- data.frame(id = c("a", "b"))
  values <- data.frame(second = c(5, 8), first = c(2, 7))
  q <- c(0, 0.5, 1, 2, 100)
  score <- evaluate_life_score(fit, data, values)
  coverage <- predict_life_coverage(fit, data, q)
  expect_identical(predict(fit, data, type = "score", values = values), score)
  expect_identical(fit$functions$score(data, values), score)
  expect_identical(predict(fit, data, type = "coverage", q = q), coverage)
  expect_identical(fit$functions$coverage(data, q), coverage)
  expect_identical(predict(fit, data), predict(fit, data, type = "moments"))
  expect_error(predict(fit, data, type = "score"), "finite numeric")
  expect_error(predict(fit, data, type = "coverage"), "q must be finite")
  expect_error(predict(fit, data, type = "moments", q = q), "q is used only")
  expect_error(predict(fit, data, type = "coverage", values = values, q = q), "values is used only")
  restored <- unserialize(serialize(fit, NULL))
  expect_identical(restored$functions$score(data, values), score)
  expect_identical(restored$functions$coverage(data, q), coverage)
})

test_that("forest coverage reuses one joint bank and preserves seed and row-order reproducibility", {
  skip_if_not_installed("rpart")
  fit <- score_survival_fixture(list(alive = 0, dead = 1), c(alive = 0.5, dead = 1),
    outcome_model = "transition_forest",
    outcome_options = list(num_trees = 3L, step_size = 2, minsplit = 4L, minbucket = 2L))
  baseline <- data.frame(id = c("a", "b", "c"))
  initial <- c(0, 1, 0)
  q <- c(0, 0.25, 0.5, 1, 2, 100)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(981)
  before <- .Random.seed
  original <- simulate_flexible_targets
  calls <- 0L
  local_mocked_bindings(simulate_flexible_targets = function(...) {
    calls <<- calls + 1L
    original(...)
  }, .package = "lifeconformal")
  result <- predict_life_coverage(fit, baseline, q, initial_state = initial)
  expect_equal(calls, 1)
  expect_identical(.Random.seed, before)
  expect_null(result$lower)
  expect_null(result$upper)
  expect_false(result$numerical$deterministic_probability_bounds)
  expect_equal(result$numerical$draws, 32)
  expect_identical(result$numerical$working_model, "transition_forest")
  expect_identical(result$probability[1, ], result$probability[3, ])
  expect_equal(result$probability[2, ], rep(1, length(q)))
  reordered <- predict_life_coverage(fit, baseline[c(3, 1, 2), , drop = FALSE], q,
                                     initial_state = initial[c(3, 1, 2)])
  expect_identical(reordered$probability, result$probability[c(3, 1, 2), , drop = FALSE])
  moments <- predict_life_moments(fit, baseline, initial)
  x <- apply_preprocessor(baseline, fit$preprocessor, "outcome")$x
  state <- match(as.character(initial), as.character(fit$problem$states))
  draws <- original(fit$outcome, x, list(time = rep(0, 3), state = state,
    duration = rep(0, 3), previous_transitions = rep(0, 3), accrued = matrix(0, 3, 2)), B = 32)
  oracle <- flexible_score_cdf(draws, moments$mean, moments$scale, q)
  expect_identical(result$probability, oracle)
  restored <- unserialize(serialize(fit, NULL))
  expect_identical(predict_life_coverage(restored, baseline, q, initial), result)
})
