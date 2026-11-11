flexible_fixture <- function(n = 16L) {
  id <- rep(seq_len(n), each = 3L)
  first <- 1 + (seq_len(n) %% 4L) / 5
  second <- first + 1.2 + (seq_len(n) %% 3L) / 5
  last <- 4 + (seq_len(n) %% 2L) / 3
  stop <- as.vector(rbind(first, second, last))
  start <- as.vector(rbind(0, first, second))
  from <- rep(c(0, 1, 0), n)
  event <- rep(1L, 3L * n)
  to <- rep(c(1, 0, 2), n)
  censor <- seq_len(n) %% 3L == 0L
  event[3L * which(censor)] <- 0L; to[3L * which(censor)] <- NA
  x <- matrix(seq_len(n) / n, n, 1L, dimnames = list(NULL, "z"))
  raw <- resolve_input(id = id, x = x[id, , drop = FALSE],
    response = LifeHist(start = start, stop = stop, from = from, to = to, event = event))
  history <- normalize_histories(raw, absorbing = 2)
  problem <- compile_problem(history, list(home = 0, alive = c(0, 1)), tau = 5)
  list(history = history, x = x, problem = problem)
}

test_that("flexible quadrature respects recurrent exposure, events and past history", {
  z <- flexible_fixture(4)
  d <- make_flexible_transition_data(z$history, z$x, z$problem, step_size = 0.4)
  expect_equal(sum(d$exposure), sum(z$history$intervals$stop - z$history$intervals$start))
  expect_equal(sum(d$to > 0), sum(z$history$intervals$event))
  expect_true(all(d$exposure > 0 & d$exposure <= 0.4 + 1e-14))
  # Third sojourn is a recurrent visit; it carries two prior transitions and
  # all occupation accumulated before returning to the initial state.
  recurrent <- which(d$subject_index == 1L & d$from == 1L & d$features[, 4L] == 2)
  expect_true(length(recurrent) > 0)
  expect_true(all(d$features[recurrent, 5L] > 1))
  expect_error(make_flexible_transition_data(z$history, z$x, z$problem,
                                            step_size = .001, max_rows = 5), "budget")
})

test_that("neural likelihood gradient agrees with numerical differentiation", {
  design <- cbind(1, c(-1, 0, 1), c(.5, -.2, .8))
  par <- seq(-.3, .3, length.out = 9)
  f <- function(p) flex_neural_objective(p, design, c(.2, 1, .5), c(1, 0, 1), 2, .03)
  analytic <- flex_neural_objective(par, design, c(.2, 1, .5), c(1, 0, 1), 2, .03, TRUE)
  numeric <- vapply(seq_along(par), function(j) {
    p <- par; p[j] <- p[j] + 1e-6
    m <- par; m[j] <- m[j] - 1e-6
    (f(p) - f(m)) / 2e-6
  }, numeric(1))
  expect_equal(analytic, numeric, tolerance = 1e-7)
  expect_true(all(is.finite(flex_softplus(c(-1000, 0, 1000)))))
  # Subdividing identical predictor/exposure pieces must not change the loss
  # or its penalty strength merely because the numerical row count increased.
  divided <- rep(seq_len(3L), each = 2L)
  expect_equal(flex_neural_objective(par, design, c(.2, 1, .5), c(1, 0, 1),
                                    2, .03, normalizer = 3),
    flex_neural_objective(par, design[divided, ], c(.1, .1, .5, .5, .25, .25),
                          c(0, 1, 0, 0, 0, 1), 2, .03, normalizer = 3))
})

test_that("neural models fit recurrent intensities and survive serialization", {
  z <- flexible_fixture()
  set.seed(123); old <- .Random.seed
  fit <- fit_transition_neural(z$history, z$x, z$problem,
    options = list(hidden = 2L, maxit = 1000L, step_size = 1), seed = 47)
  expect_identical(.Random.seed, old)
  rates <- predict_transition_rates(fit, z$x, rep(1, nrow(z$x)),
    rep(1L, nrow(z$x)), rep(1, nrow(z$x)), rep(0L, nrow(z$x)),
    matrix(1, nrow(z$x), 2L))
  expect_true(all(is.finite(rates) & rates >= 0))
  expect_true(all(rates[, fit$edges$from != 1L, drop = FALSE] == 0))
  restored <- unserialize(serialize(fit, NULL))
  expect_equal(predict_transition_rates(restored, z$x, rep(1, nrow(z$x)),
    rep(1L, nrow(z$x)), rep(1, nrow(z$x)), rep(0L, nrow(z$x)),
    matrix(1, nrow(z$x), 2L)), rates)
  expect_error(fit_transition_neural(z$history, z$x, z$problem,
                                     list(unknown = 1)), "Unknown")
})

test_that("survival-tree ensembles bootstrap subjects and preserve the graph", {
  skip_if_not_installed("rpart")
  z <- flexible_fixture()
  set.seed(123); old <- .Random.seed
  fit <- fit_transition_forest(z$history, z$x, z$problem,
    options = list(num_trees = 4, step_size = 1, minsplit = 4, minbucket = 2), seed = 23)
  expect_identical(.Random.seed, old)
  expect_identical(fit$diagnostics$resampling, "whole subjects")
  rates <- predict_transition_rates(fit, z$x, rep(1, nrow(z$x)),
    rep(1L, nrow(z$x)), rep(1, nrow(z$x)), rep(0L, nrow(z$x)),
    matrix(1, nrow(z$x), 2L))
  expect_true(all(is.finite(rates) & rates >= 0))
  expect_true(all(rates[, fit$edges$from != 1L, drop = FALSE] == 0))
  expect_equal(predict_transition_rates(unserialize(serialize(fit, NULL)),
    z$x, rep(1, nrow(z$x)), rep(1L, nrow(z$x)), rep(1, nrow(z$x)),
    rep(0L, nrow(z$x)), matrix(1, nrow(z$x), 2L)), rates)
})

constant_flexible_adapter <- function(rate = .5) {
  structure(list(type = "transition_neural", fits = list(list(constant = rate)),
    edges = data.frame(from = 1L, to = 2L),
    problem = list(tau = 3, rewards = matrix(c(1, 0, 1, 1), 2, 2,
      dimnames = list(NULL, c("before", "total"))), absorbing = 2L),
    options = list(step_size = .3, max_events = 5L), x_names = "z",
    n_predictors = 1L, seed = 45L), class = "life_flexible_outcome")
}

test_that("flexible simulation matches exponential rewards and known continuation", {
  fit <- constant_flexible_adapter()
  x <- matrix(0, 1, 1, dimnames = list(NULL, "z"))
  snap <- list(time = 0, state = 1L, duration = 0, previous_transitions = 0L,
               accrued = matrix(0, 1, 2))
  set.seed(82); old <- .Random.seed
  draws <- simulate_flexible_targets(fit, x, snap, B = 10000)
  expect_identical(.Random.seed, old)
  expect_equal(mean(draws[1, , 1]), (1 - exp(-.5 * 3)) / .5, tolerance = .025)
  expect_equal(unname(draws[1, , 2]), rep(3, 10000), tolerance = 1e-12)
  snap <- list(time = 1, state = 2L, duration = .2, previous_transitions = 1L,
               accrued = matrix(c(.8, 1), 1, 2))
  known <- simulate_flexible_targets(fit, x, snap, B = 3)
  expect_equal(unname(known[1, , 1]), rep(.8, 3))
  expect_equal(unname(known[1, , 2]), rep(3, 3))
})

test_that("continuation streams are stable across ordering, batching, and draw counts", {
  fit <- constant_flexible_adapter()
  fit$options$step_size <- .05 # Exercise a second fixed-size stream block.
  x <- matrix(c(0, 1), 2, 1, dimnames = list(NULL, "z"))
  snap <- list(time = c(0, .4), state = c(1L, 1L), duration = c(0, .4),
               previous_transitions = c(0L, 0L), accrued = rbind(c(0, 0), c(.4, .4)))
  a <- simulate_flexible_targets(fit, x, snap, B = 20)
  sr <- lapply(snap, function(v) if (is.matrix(v)) v[2:1, , drop = FALSE] else v[2:1])
  b <- simulate_flexible_targets(fit, x[2:1, , drop = FALSE], sr, B = 20)
  expect_equal(unname(a[1, , ]), unname(b[2, , ]))
  one <- lapply(snap, function(v) if (is.matrix(v)) v[1, , drop = FALSE] else v[1])
  c <- simulate_flexible_targets(fit, x[1, , drop = FALSE], one, B = 30)
  expect_equal(unname(a[1, , ]), unname(c[1, 1:20, ]))
  probs <- flexible_score_cdf(a, matrix(1, 2, 2), matrix(1, 2, 2), c(0, 1, 2, Inf))
  expect_true(all(probs[, -1] >= probs[, -4]))
  expect_equal(probs[, 4], c(1, 1))
})

test_that("flexible draw budget raises errors rather than truncating paths", {
  fit <- constant_flexible_adapter(100)
  fit$problem$absorbing <- integer(0)
  fit$edges <- data.frame(from = c(1L, 2L), to = c(2L, 1L))
  fit$fits <- list(list(constant = 100), list(constant = 100))
  fit$options$max_events <- 2L
  expect_error(simulate_flexible_targets(fit,
    matrix(0, 1, 1, dimnames = list(NULL, "z")),
    list(time = 0, state = 1L, duration = 0, previous_transitions = 0L,
         accrued = matrix(0, 1, 2)), B = 10), "max_events")
})

test_that("forest censoring matches backend survival with coherent sides and floor", {
  skip_if_not_installed("ranger")
  z <- flexible_fixture(30)
  fit <- fit_censor_forest(z$history, z$x, tau = 5,
    options = list(num_trees = 12L, min_node_size = 3L, survival_floor = .9), seed = 91)
  set.seed(61); before_prediction <- .Random.seed
  curves <- censor_forest_curves(fit, z$x[1:2, , drop = FALSE])
  expect_true(identical(.Random.seed, before_prediction))
  expect_true(all(curves$survival >= .9 & curves$survival <= 1))
  t <- curves$jump_times[length(curves$jump_times)]
  left <- predict_censor_forest(fit, z$x[1:2, , drop = FALSE], t, "left")
  right <- predict_censor_forest(fit, z$x[1:2, , drop = FALSE], t, "right")
  expect_true(all(left >= right))
  expect_equal(right, curves$survival[, ncol(curves$survival)])
  expect_equal(predict_censor_forest(unserialize(serialize(fit, NULL)),
    z$x[1:2, , drop = FALSE], t, "right"), right)
  expect_equal(predict_censor_forest(fit, z$x[1:2, , drop = FALSE], 0), c(1, 1))
  expect_error(predict_censor_forest(fit, z$x, 6), "within")
  jumps <- predict_censor_forest_jumps(fit, z$x[1:2, , drop = FALSE], c(t, 0, t),
                                      subject_index = c(2L, 1L, 1L))
  expect_equal(jumps$g_right, c(right[2L], 1, right[1L]))
  expect_equal(jumps$hazard_increment, 1 - jumps$g_right / jumps$g_left)
  expect_equal(jumps$inverse_increment, 1 / jumps$g_right - 1 / jumps$g_left)
})

test_that("no censoring events gives unit survival without fitting a forest", {
  z <- flexible_fixture(2)
  fit <- fit_censor_forest(z$history, z$x, tau = 5, seed = 2)
  expect_null(fit$model)
  expect_equal(predict_censor_forest(fit, z$x, c(2, 5)), c(1, 1))
})

test_that("censoring forest distinguishes origin censoring from initial absorption", {
  skip_if_not_installed("ranger")
  h <- normalize_histories(resolve_input(id = 1:4,
    response = LifeHist(c(0, 2, 3, 20), c(0, 0, 2, 0))),
    initial_state = c(0, 0, 0, 2), absorbing = 2)
  x <- matrix(0, 4L, 0L)
  fit <- fit_censor_forest(h, x, 5, list(num_trees = 3, min_node_size = 2), seed = 17)
  expect_equal(fit$diagnostics$n_censor_events, 2)
  expect_equal(fit$diagnostics$excluded_origin_absorption, 1)
  expect_equal(predict_censor_forest(fit, x, 0, "left"), rep(1, 4))
  expect_true(all(predict_censor_forest(fit, x, 0, "right") <= 1))
  p <- compile_problem(h, list(alive = 0), 5)
  risk <- make_flexible_transition_data(h, x, p, 1)
  expect_false(any(risk$subject_index == 4L))
})

test_that("explicit simulation seeds freeze the RNG algorithm and restore it", {
  old_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
  RNGkind("Mersenne-Twister")
  a <- flex_with_seed(21, stats::runif(5))
  RNGkind("L'Ecuyer-CMRG")
  set.seed(31); old <- .Random.seed
  b <- flex_with_seed(21, stats::runif(5))
  expect_identical(a, b)
  expect_identical(RNGkind()[1L], "L'Ecuyer-CMRG")
  expect_identical(.Random.seed, old)
})
