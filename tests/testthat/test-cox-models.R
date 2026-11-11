cox_test_history <- function(id, time, status, initial_state = NULL,
                              absorbing = 2, status_map = NULL) {
  normalize_histories(resolve_input(id = id, response = LifeHist(time, status)),
                        initial_state = initial_state, absorbing = absorbing,
                        status_map = status_map)
}

cox_test_fixture <- function() {
  i <- seq_len(96)
  time <- 1 + (i * 17L) %% 29L / 3
  status <- ifelse(i %% 3L == 0L, 2, 0)
  x <- cbind(z = sin(i * 1.3), w = cos(i * 2.1))
  list(history = cox_test_history(i, time, status), x = x, tau = 8)
}

test_that("censoring risk truncates at absorption and horizon without fake failures", {
  history <- cox_test_history(c(1, 2, 3, 4, 5, 6), c(0, 3, 7, 8, 10, 12),
                              c(0, 2, 0, 0, 0, 0),
                              initial_state = c(0, 0, 0, 0, 0, 2))
  x <- matrix(numeric(), 6, 0)
  risk <- make_censoring_data(history, x, 8)
  expect_equal(risk$data$subject_index, 1:5)
  expect_equal(risk$data$time, c(0, 3, 7, 8, 8))
  expect_equal(risk$data$event, c(1, 0, 1, 0, 0))
  expect_equal(risk$excluded_origin_absorption, 1)
  expect_equal(dim(risk$x), c(5, 0))
  adapter <- fit_censor_cox(history, x, 8)
  expect_equal(adapter$jump_times, c(0, 7))
  expect_equal(adapter$model$increments, c(1 / 5, 1 / 3))
  expect_equal(predict_censor_survival(adapter, x, 0, "left"), rep(1, 6))
  expect_equal(predict_censor_survival(adapter, x, 0, "right"), rep(exp(-1 / 5), 6))
  expect_equal(predict_censor_jumps(adapter, x, 8)$hazard_increment, rep(0, 6))
})

test_that("no-event and no-risk fits retain deterministic survival one", {
  h <- cox_test_history(1:3, c(3, 5, 7), c(2, 0, 0))
  x <- cbind(z = 1:3)
  f <- fit_censor_cox(h, x, 5)
  expect_true(f$diagnostics$no_events)
  expect_null(f$model$fit)
  expect_length(f$jump_times, 0)
  expect_equal(predict_censor_survival(f, x, c(0, 2, 5)), rep(1, 3))
  expect_equal(predict_censor_jumps(f, x, c(0, 2, 5))$inverse_increment, rep(0, 3))
  absorbed <- cox_test_history(1:2, c(4, 9), c(0, 0), initial_state = 2)
  g <- fit_censor_cox(absorbed, matrix(numeric(), 2, 0), 5)
  expect_equal(g$diagnostics$n_risk_records, 0)
  expect_equal(g$diagnostics$excluded_origin_absorption, 2)
  expect_equal(predict_censor_survival(g, matrix(numeric(), 2, 0), 5), c(1, 1))
})

test_that("near but unequal observed jump times are not merged by the backend", {
  times <- c(1, 1 + 1e-12, 2, 3)
  h <- cox_test_history(1:4, times, rep(0, 4), absorbing = integer(0))
  f <- fit_censor_cox(h, matrix(numeric(), 4, 0), 4)
  expect_identical(f$jump_times, times)
  expect_equal(f$model$increments, 1 / (4:1))
  expect_equal(predict_censor_survival(f, matrix(numeric(), 1, 0),
                                       1 + 0.5e-12), exp(-1 / 4))
})

test_that("Cox censoring predictions agree with independent survival curves", {
  fixture <- cox_test_fixture()
  risk <- make_censoring_data(fixture$history, fixture$x, fixture$tau)
  d <- data.frame(time = risk$data$time, event = risk$data$event,
                  z = risk$x[, "z"], w = risk$x[, "w"])
  for (ties in c("breslow", "efron")) {
    f <- fit_censor_cox(fixture$history, fixture$x, fixture$tau, list(ties = ties))
    reference <- survival::coxph(survival::Surv(time, event) ~ z + w,
                                  data = d, ties = ties, nocenter = NULL,
                                  control = survival::coxph.control(timefix = FALSE))
    new <- data.frame(z = c(-1, 0.2, 0.8), w = c(0.5, -0.2, 1))
    curves <- survival::survfit(reference, newdata = new, se.fit = FALSE)
    time <- c(0, 2.5, 7.75)
    index <- findInterval(time, curves$time)
    expected <- vapply(seq_along(time), function(i) {
      if (!index[i]) 1 else curves$surv[index[i], i]
    }, numeric(1))
    expect_equal(predict_censor_survival(f, as.matrix(new), time), expected,
                  tolerance = 1e-12)
    expect_equal(unname(f$model$coefficients), unname(stats::coef(reference)),
                  tolerance = 1e-12)
    expect_identical(f$controls$ties, ties)
  }
})

test_that("censoring floor preserves left/right and discrete jump coherence", {
  h <- cox_test_history(1:5, 1:5, rep(0, 5), absorbing = integer(0))
  x <- matrix(numeric(), 5, 0)
  f <- fit_censor_cox(h, x, 6, list(survival_floor = 0.5))
  jumps <- predict_censor_jumps(f, x, 1:5)
  expect_equal(jumps$g_left, c(1, head(jumps$g_right, -1)))
  expect_equal(jumps$g_right, pmax(0.5, exp(-cumsum(1 / (5:1)))))
  expect_equal(jumps$hazard_increment, 1 - jumps$g_right / jumps$g_left)
  expect_equal(jumps$inverse_increment, 1 / jumps$g_right - 1 / jumps$g_left)
  expect_equal(sum(jumps$inverse_increment), 1 / tail(jumps$g_right, 1) - 1)
  expect_equal(jumps$hazard_increment[1], -expm1(-1 / 5))
  expect_false(isTRUE(all.equal(jumps$hazard_increment[1], 1 / 5)))
  expect_equal(predict_censor_survival(f, x, 1:5, "left"), jumps$g_left)
  expect_equal(predict_censor_survival(f, x, 1:5), jumps$g_right)
  repeated <- predict_censor_jumps(f, x[1, , drop = FALSE], 1:5, rep(1L, 5))
  expect_equal(repeated$g_right, jumps$g_right)
  expect_equal(predict_censor_jumps(f, x, rep(3.5, 5))$inverse_increment, rep(0, 5))
  all_one <- fit_censor_cox(h, x, 6, list(survival_floor = 1))
  expect_equal(predict_censor_jumps(all_one, x, 1:5)$g_right, rep(1, 5))
  # A frozen model with an extreme finite linear predictor still has coherent
  # finite weights at the computational floor; exp(lp) is never formed.
  extreme <- f
  extreme$n_predictors <- 1L
  extreme$predictor_names <- "z"
  extreme$model$coefficients <- c(z = 1000)
  extreme$model$center <- c(z = 0)
  extreme$controls$survival_floor <- .Machine$double.xmin
  tail <- predict_censor_jumps(extreme, cbind(z = 100), 1)
  expect_true(is.finite(tail$inverse_increment))
  expect_equal(tail$hazard_increment, 1)
  expect_equal(tail$g_right, .Machine$double.xmin, tolerance = 1e-14)
})

test_that("cause-specific risk views retain recurrent entries and calendar time", {
  h <- cox_test_history(c(1, 1, 1, 1, 2, 2, 3),
                         c(2, 4, 6, 9, 3, 8, 10), c(1, 3, 1, 2, 1, 0, 0),
                         status_map = c("1" = 1, "3" = 0, "2" = 2))
  x <- cbind(z = c(7, 8, 9))
  origin <- match(0, h$schema$state_labels)
  dest <- match(1, h$schema$state_labels)
  risk <- make_transition_risk_data(h, c(origin, dest), x, 8)
  expect_equal(risk$data$start, c(0, 4, 0, 0))
  expect_equal(risk$data$stop, c(2, 6, 3, 10))
  expect_equal(risk$data$event, c(1, 1, 1, 0))
  expect_equal(risk$x[, 1], c(7, 7, 8, 9))
  f <- fit_markov_cox(h, matrix(numeric(), 3, 0), 8)
  edge <- which(f$edges$from == origin & f$edges$to == dest)
  expect_equal(f$models[[edge]]$jump_times, c(2, 3, 6))
  expect_equal(f$models[[edge]]$increments, c(1 / 3, 1 / 2, 1 / 2))
  expect_length(f$diagnostics$zero_event_edges, 0)
  expect_equal(f$diagnostics$zero_prediction_event_edges,
                 paste(match(1, h$schema$state_labels), match(2, h$schema$state_labels),
                       sep = "->"))
})

test_that("transition Cox fits agree with separate start-stop survival models", {
  fixture <- cox_test_fixture()
  f <- fit_markov_cox(fixture$history, fixture$x, fixture$tau)
  expect_equal(nrow(f$edges), 1)
  risk <- make_transition_risk_data(fixture$history,
                                    unlist(f$edges[1, ]), fixture$x, fixture$tau)
  d <- cbind(risk$data, as.data.frame(risk$x))
  reference <- survival::coxph(survival::Surv(start, stop, event) ~ z + w,
                                data = d, ties = "breslow", nocenter = NULL,
                                control = survival::coxph.control(timefix = FALSE))
  b <- survival::basehaz(reference, centered = TRUE)
  positive <- c(TRUE, diff(b$hazard) > 0) & b$hazard > 0 & b$time <= fixture$tau
  expect_equal(f$jump_times, b$time[positive])
  expect_equal(as.vector(f$base_increments), diff(c(0, b$hazard[positive])),
                 tolerance = 1e-12)
  new <- fixture$x[c(3, 5, 7), , drop = FALSE]
  predicted <- predict_transition_jumps(f, new)
  expected <- exp(sweep(new, 2, reference$means, "-") %*% stats::coef(reference))
  expect_equal(unname(predicted$relative_risk), unname(expected), tolerance = 1e-12)
  expect_equal(predict_transition_jumps(f, new, c(2, 5))$jump_times,
                 f$jump_times[f$jump_times > 2 & f$jump_times <= 5])
})

test_that("outcome coefficients use all actual follow-up and prediction jumps stop at tau", {
  fixture <- cox_test_fixture()
  f <- fit_markov_cox(fixture$history, fixture$x, fixture$tau)
  risk <- make_transition_risk_data(fixture$history,
                                    unlist(f$edges[1, ]), fixture$x, fixture$tau)
  expect_gt(max(risk$data$stop), fixture$tau)
  expect_gt(sum(risk$data$event), sum(risk$data$event[risk$data$stop <= fixture$tau]))
  d <- cbind(risk$data, as.data.frame(risk$x))
  full <- survival::coxph(survival::Surv(start, stop, event) ~ z + w, data = d,
                          ties = "breslow", nocenter = NULL,
                          control = survival::coxph.control(timefix = FALSE))
  truncated <- d
  truncated$event <- truncated$event * (truncated$stop <= fixture$tau)
  truncated$stop <- pmin(truncated$stop, fixture$tau)
  truncated <- survival::coxph(survival::Surv(start, stop, event) ~ z + w,
                                data = truncated, ties = "breslow", nocenter = NULL,
                                control = survival::coxph.control(timefix = FALSE))
  expect_equal(unname(f$models[[1]]$coefficients), unname(stats::coef(full)), tolerance = 1e-12)
  expect_gt(max(abs(stats::coef(full) - stats::coef(truncated))), 1e-5)
  expect_true(all(f$jump_times <= fixture$tau))
  expect_equal(f$models[[1]]$diagnostics$n_events, sum(d$event))
  censor <- make_censoring_data(fixture$history, fixture$x, fixture$tau)
  expect_lte(max(censor$data$time), fixture$tau)
  expect_false(any(censor$data$time == fixture$tau & censor$data$event == 1))
})

cox_exact_tie_fixture <- function() {
  records <- lapply(seq_len(48L), function(i) {
    first <- 0.5 + (i %% 7L) / 4
    second <- first + 1 + (i %% 5L) / 3
    pattern <- as.character(i %% 4L)
    times <- switch(pattern, "0" = first / 2, "1" = c(first, second),
                      "2" = second, "3" = c(first, second))
    status <- switch(pattern, "0" = 0L, "1" = c(1L, 0L),
                       "2" = 2L, "3" = c(1L, 2L))
    data.frame(id = paste0("train", i), time = times, status = status, x = sin(i * .91))
  })
  do.call(rbind, records)
}

test_that("the original unjittered near-tie fixture has correct Cox risk denominators", {
  data <- cox_exact_tie_fixture()
  raw <- resolve_input(LifeHist(time, status) ~ x, data = data)
  history <- normalize_histories(raw, absorbing = 2)
  x <- matrix(history$baseline$x, ncol = 1, dimnames = list(NULL, "x"))
  for (ties in c("breslow", "efron")) {
    expect_no_warning(f <- fit_markov_cox(history, x, 6, list(ties = ties)))
    for (edge in seq_len(nrow(f$edges))) {
      model <- f$models[[edge]]
      risk <- make_transition_risk_data(history, unlist(f$edges[edge, ]), x, 6)
      lp <- model$fit$linear.predictors
      times <- sort(unique(risk$data$stop[risk$data$event == 1]))
      reference <- vapply(times, function(t) {
        active <- risk$data$start < t & risk$data$stop >= t
        deaths <- risk$data$event == 1 & risk$data$stop == t
        total <- sum(exp(lp[active])); d <- sum(deaths)
        if (ties == "breslow") d / total else
          sum(1 / (total - (seq_len(d) - 1) / d * sum(exp(lp[deaths]))))
      }, numeric(1))
      detail <- survival::coxph.detail(model$fit)
      expect_identical(model$jump_times, times)
      expect_equal(model$increments, reference, tolerance = 1e-13)
      expect_equal(model$increments, detail$hazard, tolerance = 1e-13)
      expect_true(all(model$increments > 0))
    }
  }
})

test_that("baseline sweep preserves strict entry sides, ties, and extreme log risks", {
  start <- c(0, 0, 0, 1, 1.5, 2)
  stop <- c(2, 2, 3, 3, 3, 4)
  event <- c(1, 1, 1, 1, 0, 1)
  lp <- log(c(2, 3, 4, 5, 6, 7))
  for (ties in c("breslow", "efron")) {
    got <- .cox_centered_baseline(stop, event, lp, start, ties)
    reference <- vapply(got$time, function(t) {
      active <- start < t & t <= stop
      deaths <- event == 1 & stop == t
      d <- sum(deaths); r <- sum(exp(lp[active]))
      if (ties == "breslow") d / r else
        sum(1 / (r - (seq_len(d) - 1) / d * sum(exp(lp[deaths]))))
    }, numeric(1))
    expect_equal(got$increments, reference, tolerance = 1e-14)
  }
  extreme <- .cox_centered_baseline(c(1, 3, 3, 4), c(0, 1, 1, 1),
    c(700, -700, -710, -705), c(0, 0, 2, 2))
  expect_equal(log(extreme$increments), c(log(2) + 700 - log(1 + exp(-10) + exp(-5)), 705),
                 tolerance = 1e-13)
  origin <- .cox_centered_baseline(c(0, 2), c(1, 1), c(700, -700))
  expect_equal(log(origin$increments), c(-700, 700), tolerance = 1e-13)
  boundary <- .cox_centered_baseline(c(2, 3), c(1, 1), c(0, 0), c(0, 2))
  expect_equal(boundary$increments, c(1, 1))
})

test_that("transition diagnostics display the original state labels", {
  data <- data.frame(id = 1:8, time = 1:8, status = rep(2, 8))
  h <- normalize_histories(resolve_input(LifeHist(time, status) ~ 1, data = data),
                            initial_state = "well", absorbing = "dead",
                            status_map = c("2" = "dead"))
  expect_warning(fit_markov_cox(h, cbind(constant = rep(1, 8)), 10),
                   "Transition well->dead.*unidentifiable")
})

test_that("empirical edges remain frozen including after-horizon transitions", {
  h <- cox_test_history(c(1, 1, 2, 3), c(12, 15, 3, 8), c(1, 0, 2, 0))
  x <- matrix(numeric(), 3, 0)
  f <- fit_markov_cox(h, x, 10)
  expect_equal(nrow(f$edges), 2)
  expect_equal(unname(f$base_increments[, 1]), 0)
  reordered <- f$edges[2:1, ]
  g <- fit_markov_cox(h, x, 10, edges = reordered)
  rownames(reordered) <- NULL
  expect_equal(g$edges, reordered)
  expect_equal(unname(g$base_increments), unname(f$base_increments[, 2:1, drop = FALSE]))
  expect_error(fit_markov_cox(h, x, 10, edges = f$edges[1, , drop = FALSE]),
                 "exactly the distinct observed")
  no_edges <- cox_test_history(1:2, c(4, 8), c(0, 0), absorbing = 2)
  empty <- fit_markov_cox(no_edges, matrix(numeric(), 2, 0), 10)
  expect_equal(dim(empty$base_increments), c(0, 0))
  expect_equal(dim(predict_transition_jumps(empty, matrix(numeric(), 3, 0))$relative_risk),
                 c(3, 0))
})

test_that("fitted batch moments match independent per-profile Markov recursion", {
  fixture <- cox_test_fixture()
  f <- fit_markov_cox(fixture$history, fixture$x, fixture$tau)
  rewards <- matrix(c(1, 0), 2, 1, dimnames = list(NULL, "alive"))
  x <- fixture$x[c(2, 10, 50), , drop = FALSE]
  result <- predict_markov_cox_moments(f, x, rewards, 1)
  jumps <- predict_transition_jumps(f, x)
  for (i in seq_len(nrow(x))) {
    kernels <- lapply(seq_along(jumps$jump_times), function(j) {
      increment <- matrix(0, 2, 2)
      increment[cbind(jumps$edges$from, jumps$edges$to)] <-
        jumps$base_increments[j, ] * jumps$relative_risk[i, ]
      make_markov_kernel(increment, absorbing = f$absorbing)$kernel
    })
    reference <- markov_reward_moments(jumps$jump_times, kernels, rewards, f$tau)
    expect_equal(unname(result$mean[i, ]), unname(reference$mean[1, ]), tolerance = 1e-12)
    expect_equal(unname(result$variance[i, ]), unname(reference$variance[1, ]),
                   tolerance = 1e-12)
  }
})

test_that("centering survives large predictor origins and saved-object reloads", {
  fixture <- cox_test_fixture()
  shifted <- fixture$x
  shifted[, 1] <- shifted[, 1] + 1e6
  f <- fit_censor_cox(fixture$history, shifted, fixture$tau)
  unshifted <- fit_censor_cox(fixture$history, fixture$x, fixture$tau)
  expect_equal(predict_censor_survival(f, shifted, 4),
                 predict_censor_survival(unshifted, fixture$x, 4), tolerance = 1e-9)
  g <- fit_markov_cox(fixture$history, shifted, fixture$tau)
  filename <- tempfile(fileext = ".rds")
  saveRDS(list(censor = f, outcome = g, x = shifted), filename)
  reread <- readRDS(filename)
  expect_equal(predict_censor_survival(reread$censor, reread$x, 4),
                 predict_censor_survival(f, shifted, 4))
  expect_equal(predict_transition_jumps(reread$outcome, reread$x),
                 predict_transition_jumps(g, shifted))
  expect_identical(environment(f$model$fit$terms), asNamespace("survival"))
  # A genuinely fresh process verifies that extraction and prediction do not
  # depend on local risk-frame or formula objects left in the fitting call.
  functions_file <- tempfile(fileext = ".R")
  script <- tempfile(fileext = ".R")
  output <- tempfile(fileext = ".rds")
  dump(c("predict_censor_survival", "predict_transition_jumps",
          ".cox_censor_queries", ".cox_check_design", ".cox_linear_predictor",
          ".cox_log_survival"), functions_file,
        envir = environment(predict_censor_survival))
  writeLines(c("args <- commandArgs(TRUE)", "source(args[1])",
                "a <- readRDS(args[2])",
                "result <- list(survival = predict_censor_survival(a$censor, a$x, 4),",
                " jumps = predict_transition_jumps(a$outcome, a$x),",
                " baseline = survival::basehaz(a$censor$model$fit, centered=TRUE))",
                "saveRDS(result, args[3])"), script)
  status <- system2(file.path(R.home("bin"), "Rscript"),
                      c("--vanilla", shQuote(script), shQuote(functions_file),
                        shQuote(filename), shQuote(output)), stdout = TRUE, stderr = TRUE)
  expect_null(attr(status, "status"), info = paste(status, collapse = "\n"))
  fresh <- readRDS(output)
  expect_equal(fresh$survival, predict_censor_survival(f, shifted, 4))
  expect_equal(fresh$jumps, predict_transition_jumps(g, shifted))
  expect_equal(fresh$baseline, survival::basehaz(f$model$fit, centered = TRUE))
})

test_that("singular, unstable and incompatible Cox inputs fail or diagnose explicitly", {
  fixture <- cox_test_fixture()
  collinear <- cbind(fixture$x, duplicate = fixture$x[, 1])
  expect_warning(f <- fit_censor_cox(fixture$history, collinear, fixture$tau),
                   "unidentifiable.*duplicate")
  expect_equal(f$diagnostics$singular_columns, "duplicate")
  expect_equal(unname(f$model$coefficients["duplicate"]), 0)
  constant <- cbind(z = rep(9, nrow(fixture$x)))
  expect_warning(constant_fit <- fit_censor_cox(fixture$history, constant, fixture$tau),
                   "unidentifiable.*z")
  intercept <- fit_censor_cox(fixture$history,
                              matrix(numeric(), nrow(fixture$x), 0), fixture$tau)
  expect_equal(predict_censor_survival(constant_fit, constant, 5),
                 predict_censor_survival(intercept,
                                         matrix(numeric(), nrow(fixture$x), 0), 5))
  expect_error(fit_censor_cox(fixture$history, collinear, fixture$tau,
                              list(singular = "error")), "unidentifiable")
  expect_error(fit_censor_cox(fixture$history, fixture$x, fixture$tau,
                              list(iter.max = 1L)), "iteration limit|did not converge")
  separated <- cox_test_history(1:10, 1:10, rep(0, 10), absorbing = integer(0))
  expect_error(fit_censor_cox(separated, cbind(z = 1:10), 11),
                 "failed.*(infinite|did not converge)")
  expect_error(fit_censor_cox(fixture$history, fixture$x, fixture$tau,
                              list(typo = 1)), "supported controls")
  expect_error(predict_censor_survival(f, collinear[, 3:1], 2), "frozen")
  expect_error(predict_censor_survival(f, collinear, 9), "within")
  expect_error(predict_censor_jumps(f, collinear, 2, subject_index = 0), "index")
  expect_error(make_censoring_data(fixture$history, fixture$x[-1, ], 5), "one row")
})
