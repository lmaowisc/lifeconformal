backend_workflow_data <- function(prefix, n, seed) {
  flex_with_seed(seed, {
    z <- stats::runif(n, -.8, .8)
    t1 <- stats::runif(n, .4, 1.2)
    t2 <- t1 + stats::runif(n, .5, 1.1)
    t3 <- t2 + stats::runif(n, 1, 2.5)
    records <- lapply(seq_len(n), function(i) {
      pattern <- i %% 4L
      time <- switch(as.character(pattern),
        "0" = t1[i] / 2, "1" = c(t1[i], (t1[i] + t2[i]) / 2),
        "2" = c(t1[i], t2[i], t3[i]), "3" = c(t1[i], t2[i], t3[i]))
      status <- switch(as.character(pattern),
        "0" = 0, "1" = c(1, 0), "2" = c(1, 7, 0), "3" = c(1, 7, 2))
      data.frame(id = paste0(prefix, i), time = time, status = status, z = z[i])
    })
    do.call(rbind, records)
  })
}

backend_workflow_options <- function(outcome) {
  switch(outcome, markov_cox = list(),
    transition_forest = list(num_trees = 3L, step_size = .75,
                             minsplit = 5L, minbucket = 2L, max_depth = 4L),
    transition_neural = list(hidden = 2L, decay = .1, step_size = .75,
                             maxit = 1000L))
}

backend_workflow_fit <- function(train, censor, outcome,
                                  targets = list(alive = c(0, 1)), ...) {
  life_train(LifeHist(time, status) ~ z, data = train, tau = 4.5,
    absorbing = 2, status_map = c("1" = 1, "7" = 0, "2" = 2),
    targets = targets, censor_model = censor, outcome_model = outcome,
    censor_options = if (censor == "cox") list() else
      list(num_trees = 5L, min_node_size = 3L),
    outcome_options = backend_workflow_options(outcome),
    control = list(draws = 32L, keep_contributions = TRUE), seed = 902L, ...)
}

backend_forest_stopped_reference <- function(fit, data, q) {
  history <- normalize_histories(resolve_input(data = data, schema = fit$schema),
    initial_state = 0, status_map = fit$schema$status_map, absorbing = 2,
    states = fit$problem$states)
  x <- apply_preprocessor(history$baseline, fit$preprocessor, "outcome")$x
  censor_x <- apply_preprocessor(history$baseline, fit$preprocessor, "censor")$x
  moments <- predict(fit, history$baseline)
  cdt <- compute_cdt(history, fit$problem, moments$mean, moments$scale, q)
  prepared <- prepare_occupation_history(history, fit$problem)
  raw <- flex_with_seed(89L, stats::predict(fit$censor$model,
    data = as.data.frame(censor_x), num.threads = 1L, seed = 89L))
  # This lookup uses ranger's unfiltered grid directly, independently of the
  # adapter's reduced censor-event grid and the calibration query table.
  g <- function(i, time, side) {
    index <- findInterval(time, raw$unique.death.times, left.open = side == "left")
    out <- rep(1, length(time)); use <- index > 0L
    out[use] <- pmax(raw$survival[i, index[use]], fit$control$survival_floor)
    out
  }
  n <- nrow(history$subjects)
  ipw <- aipw <- matrix(0, n, length(q))
  for (i in seq_len(n)) {
    jumps <- predict_transition_jumps(fit$outcome, x[i, , drop = FALSE])
    k <- length(fit$problem$states)
    kernels <- lapply(seq_along(jumps$jump_times), function(j) {
      increment <- matrix(0, k, k)
      increment[cbind(jumps$edges$from, jumps$edges$to)] <-
        jumps$base_increments[j, ] * jumps$relative_risk[1L, ]
      make_markov_kernel(increment, absorbing = fit$problem$absorbing)$kernel
    })
    times <- sort(unique(c(cdt$time[i, ],
      fit$censor$jump_times[fit$censor$jump_times < max(cdt$time[i, ])])))
    snap <- evaluate_occupation_history(prepared, rep(i, length(times)), times, "right")
    h <- length(times)
    probability <- markov_conditional_coverage_monotone(jumps$jump_times, kernels,
      fit$problem$rewards, fit$problem$tau, times, snap$state, snap$accrued,
      matrix(rep(moments$mean[i, ], each = h), h),
      matrix(rep(moments$scale[i, ], each = h), h), q)$probability
    certified <- outer(times, cdt$time[i, ], ">=") &
      matrix(rep(cdt$resolved[i, ], each = h), h)
    exact <- matrix(rep(as.numeric(cdt$membership[i, ]), each = h), h)
    probability[certified] <- exact[certified]
    for (j in seq_along(q)) {
      stop_time <- cdt$time[i, j]
      endpoint <- match(stop_time, times)
      used <- which(times %in% fit$censor$jump_times & times < stop_time)
      increment <- data.frame(subject_index = rep(1L, length(used)), time = times[used],
        g_left = g(i, times[used], "left"), g_right = g(i, times[used], "right"),
        m_plus = probability[used, j])
      gl <- g(i, stop_time, "left")
      ipw[i, j] <- compute_ipw_contributions(cdt$membership[i, j], cdt$resolved[i, j], gl)
      aipw[i, j] <- compute_aipw_contributions(cdt$membership[i, j],
        cdt$resolved[i, j], stop_time, gl, probability[endpoint, j], increment)
    }
  }
  list(ipw = ipw, aipw = aipw)
}

test_that("all censoring and outcome backend pairs complete both public workflows", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("rpart")
  train <- backend_workflow_data("train", 32L, 1357L)
  calibration <- backend_workflow_data("cal", 8L, 2468L)
  patients <- data.frame(id = c("new_a", "new_b"), z = c(-.2, .3))
  set.seed(83); before_rng <- .Random.seed
  for (censor in c("cox", "survival_forest")) {
    for (outcome in c("markov_cox", "transition_forest", "transition_neural")) {
      label <- paste(censor, outcome)
      fit <- backend_workflow_fit(train, censor, outcome)
      expect_identical(.Random.seed, before_rng, info = label)
      expect_equal(nrow(fit$problem$working_graph$edges), 3L, info = label)
      moments <- predict(fit, patients)
      expect_true(all(is.finite(moments$mean) & moments$mean >= 0 &
                        moments$mean <= 4.5), info = label)
      expect_true(all(moments$variance >= 0 & moments$scale > 0), info = label)
      frozen <- serialize(fit, NULL)
      for (method in c("ipw", "aipw")) {
        cal <- life_calibrate(fit, data = calibration, method = method,
                              alpha = .3, q = c(0, .75, 1.5, 4))
        expect_identical(.Random.seed, before_rng, info = paste(label, method))
        expect_true(identical(serialize(fit, NULL), frozen), info = paste(label, method))
        expect_true(all(is.finite(cal$curve$estimate)), info = paste(label, method))
        expect_equal(dim(cal$contributions), c(8L, 4L), info = paste(label, method))
        expect_equal(cal$curve$estimate, colMeans(cal$contributions),
                     info = paste(label, method))
        region <- predict(cal, patients)
        expect_s3_class(region, "lifeconformal_prediction")
        expect_true(all(region$bounds$lower >= 0 & region$bounds$upper <= 4.5 &
                          region$bounds$lower <= region$bounds$upper),
                    info = paste(label, method))
        expect_identical(.Random.seed, before_rng, info = paste(label, method))
        if (method == "ipw") {
          expect_false(cal$diagnostics$numerical$conditional_coverage_evaluated)
        } else {
          # These independent calibration histories include censoring before
          # outcome determination, so augmentation must actually be exercised.
          expect_true(cal$diagnostics$numerical$queried_snapshots > 0,
                      info = paste(label, method))
        }
      }
    }
  }
})

test_that("flexible component allocation preserves the full history feature model", {
  skip_if_not_installed("rpart")
  train <- backend_workflow_data("train", 32L, 1357L)
  calibration <- backend_workflow_data("cal", 12L, 2468L)
  targets <- list(home = 0, alive = c(0, 1))
  patients <- data.frame(id = c("new_a", "new_b"), z = c(-.2, .3))
  for (outcome in c("transition_forest", "transition_neural")) {
    fit <- backend_workflow_fit(train, "cox", outcome, targets)
    before <- predict(fit, patients)
    frozen <- serialize(fit, NULL)
    joint <- life_calibrate(fit, data = calibration, alpha = .4,
      alpha_components = c(home = .2, alive = .2), method = "aipw",
      q = c(0, .75, 1.5, 4), control = list(draws = 48L))
    expect_true(identical(serialize(fit, NULL), frozen))
    expect_identical(predict(fit, patients), before)
    expect_equal(joint$alpha_components, c(home = .2, alive = .2))
    expect_true(all(is.finite(joint$threshold) | is.infinite(joint$threshold)))
    expect_equal(predict(joint, patients)$center, before$mean)
    expect_equal(predict(joint, patients)$scale, before$scale)
    bounds <- predict(joint, patients)$bounds
    expect_equal(bounds$target, rep(c("home", "alive"), 2L))
    expect_true(all(bounds$lower >= 0 & bounds$upper <= 4.5 & bounds$lower <= bounds$upper))
  }
})

test_that("forest weighting agrees with direct backend survival and stopped augmentation", {
  skip_if_not_installed("ranger")
  train <- backend_workflow_data("train", 24L, 1357L)
  calibration <- backend_workflow_data("cal", 4L, 2468L)
  fit <- backend_workflow_fit(train, "survival_forest", "markov_cox")
  q <- c(0, .75, 1.5, 4)
  reference <- backend_forest_stopped_reference(fit, calibration, q)
  expect_equal(fit$censor$jump_times,
               sort(unique(train$time[train$status == 0 & train$time < 4.5])))
  for (method in c("ipw", "aipw")) {
    cal <- life_calibrate(fit, data = calibration, method = method,
                          alpha = .4, q = q)
    expect_equal(unname(cal$contributions), unname(reference[[method]]), tolerance = 1e-10)
  }
})

test_that("a separate flexible augmentation backend leaves the trained score unchanged", {
  skip_if_not_installed("rpart")
  train <- backend_workflow_data("train", 32L, 1357L)
  calibration <- backend_workflow_data("cal", 8L, 2468L)
  fit <- backend_workflow_fit(train, "cox", "transition_forest",
    targets = list(home = 0, alive = c(0, 1)),
    augmentation_model = "transition_neural",
    augmentation_options = backend_workflow_options("transition_neural"))
  expect_false(fit$augmentation_reused)
  expect_identical(fit$outcome$type, "transition_forest")
  expect_identical(fit$augmentation$type, "transition_neural")
  patients <- data.frame(z = c(-.2, .3))
  before <- predict(fit, patients)
  frozen <- serialize(fit, NULL)
  a <- life_calibrate(fit, data = calibration, alpha = .3, q = c(0, .75, 1.5, 4),
                      control = list(draws = 48L))
  expect_true(identical(serialize(fit, NULL), frozen))
  expect_equal(predict(a, patients)$center, before$mean)
  expect_equal(predict(a, patients)$scale, before$scale)
  expect_gt(a$diagnostics$numerical$queried_snapshots, 0)
  # Repeating augmentation with its fixed profile streams reproduces the curve.
  b <- life_calibrate(fit, data = calibration, alpha = .3, q = c(0, .75, 1.5, 4),
                      control = list(draws = 48L))
  expect_identical(a$curve, b$curve)
})

test_that("full flexible fits and calibrations reproduce after serialization", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("rpart")
  train <- backend_workflow_data("train", 32L, 1357L)
  calibration <- backend_workflow_data("cal", 8L, 2468L)
  fit <- backend_workflow_fit(train, "survival_forest", "transition_forest",
                              targets = list(home = 0, alive = c(0, 1)))
  cal <- life_calibrate(fit, data = calibration, alpha = .3, q = c(0, .75, 1.5, 4))
  patients <- data.frame(id = c("a", "b"), z = c(-.2, .3))
  expected <- predict(cal, patients)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(cal, path)
  restored <- readRDS(path)
  old_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
  RNGkind("L'Ecuyer-CMRG"); set.seed(12)
  before <- .Random.seed
  repeated_fit <- backend_workflow_fit(train, "survival_forest", "transition_forest",
                                        targets = list(home = 0, alive = c(0, 1)))
  expect_identical(.Random.seed, before)
  expect_identical(predict(repeated_fit, patients), predict(fit, patients))
  expect_identical(predict(restored, patients), expected)
  expect_identical(.Random.seed, before)
  expect_identical(restored$fit$functions$moments(patients), predict(fit, patients))
})

test_that("an installed saved flexible workflow runs in a fresh R process", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("rpart")
  package_path <- getNamespaceInfo(asNamespace("lifeconformal"), "path")
  skip_if_not(file.exists(file.path(package_path, "Meta", "package.rds")),
              "Fresh-process library loading is exercised on the installed package.")
  train <- backend_workflow_data("train", 24L, 1357L)
  calibration <- backend_workflow_data("cal", 8L, 2468L)
  fit <- backend_workflow_fit(train, "survival_forest", "transition_neural")
  cal <- life_calibrate(fit, data = calibration, alpha = .3, q = c(0, .75, 1.5, 4))
  patients <- data.frame(id = c("a", "b"), z = c(-.2, .3))
  bundle <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  log <- tempfile(fileext = ".txt")
  on.exit(unlink(c(bundle, script, log)), add = TRUE)
  saveRDS(list(calibration = cal, patients = patients,
               expected = predict(cal, patients), moments = predict(fit, patients)), bundle)
  writeLines(c(
    sprintf("library(lifeconformal, lib.loc = %s)", deparse(dirname(package_path))),
    sprintf("saved <- readRDS(%s)", deparse(bundle)),
    "RNGkind(\"L'Ecuyer-CMRG\"); set.seed(107)",
    "before <- .Random.seed",
    "stopifnot(identical(predict(saved$calibration, saved$patients), saved$expected))",
    "stopifnot(identical(saved$calibration$fit$functions$moments(saved$patients), saved$moments))",
    "stopifnot(identical(.Random.seed, before))"), script)
  status <- system2(file.path(R.home("bin"), "Rscript"),
                     c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  expect_equal(status, 0L, info = paste(readLines(log, warn = FALSE), collapse = "\n"))
})
