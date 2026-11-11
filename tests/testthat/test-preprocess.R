test_that("formula and explicit matrix routes give aligned subject-level designs", {
  baseline <- data.frame(age = c(51, 66, 43, 74),
                         treatment = factor(c("B", "A", "B", "A"), levels = c("A", "B")))
  formula_fit <- fit_preprocessor(baseline, response_not_evaluated() ~ age + treatment)
  vector_fit <- fit_preprocessor(baseline)
  expect_identical(formula_fit$outcome$x, vector_fit$outcome$x)
  expect_identical(formula_fit$outcome, formula_fit$censor)
  expect_identical(formula_fit$preprocessor$outcome, formula_fit$preprocessor$censor)
  expect_named(formula_fit$outcome$baseline, c("age", "treatment"))
  expect_equal(unname(formula_fit$outcome$x[, "age"]), baseline$age)
  expect_identical(apply_preprocessor(baseline[2:1], formula_fit$preprocessor)$x,
                   formula_fit$outcome$x)
  explicit <- matrix(seq_len(12), nrow = 4)
  matrix_fit <- fit_preprocessor(explicit)
  expect_identical(colnames(matrix_fit$outcome$x), c("x1", "x2", "x3"))
  expect_equal(unname(matrix_fit$outcome$x[, 2]), as.numeric(explicit[, 2]))
})

test_that("censoring uses its own frozen predictor selection", {
  baseline <- data.frame(age = 40:49, treatment = rep(c("A", "B"), 5),
                         censor_score = c(3:11, 20))
  fit <- fit_preprocessor(baseline, ~ age + treatment, ~ scale(censor_score))
  expect_named(fit$outcome$baseline, c("age", "treatment"))
  expect_named(fit$censor$baseline, "censor_score")
  expect_named(fit$baseline, c("age", "treatment", "censor_score"))
  new <- data.frame(censor_score = 20)
  expect_equal(as.vector(apply_preprocessor(new, fit$preprocessor, "censor")$x),
               (20 - mean(baseline$censor_score)) / stats::sd(baseline$censor_score))
  expect_error(apply_preprocessor(new, fit$preprocessor), "Missing baseline column")
  expect_error(fit_preprocessor(baseline, ~age, censor_score ~age), "one-sided")
})

test_that("scale and nested transformations use training constants", {
  baseline <- data.frame(x = c(1, 2, 4, 7, 12))
  new <- data.frame(x = c(10, 100))
  fit <- fit_preprocessor(baseline, ~ scale(x) + I(scale(log(x))^2))
  expected <- cbind((new$x - mean(baseline$x)) / stats::sd(baseline$x),
                    ((log(new$x) - mean(log(baseline$x))) / stats::sd(log(baseline$x)))^2)
  actual <- apply_preprocessor(new, fit$preprocessor)$x
  expect_equal(as.vector(actual), as.vector(expected), tolerance = 1e-13)
  expect_equal(as.vector(apply_preprocessor(new[1, , drop = FALSE], fit$preprocessor)$x),
               as.vector(actual[1, , drop = FALSE]))
  expect_equal(as.vector(apply_preprocessor(new[2, , drop = FALSE], fit$preprocessor)$x),
               as.vector(actual[2, , drop = FALSE]))
  rms <- fit_preprocessor(baseline, ~ scale(x, center = FALSE))
  expected_rms <- new$x / sqrt(sum(baseline$x^2) / (nrow(baseline) - 1))
  expect_equal(as.vector(apply_preprocessor(new, rms$preprocessor)$x), expected_rms)
  fixed <- fit_preprocessor(baseline, ~ scale(x, scale = 5, center = 2))
  expect_equal(as.vector(apply_preprocessor(new, fixed$preprocessor)$x), (new$x - 2) / 5)
  unscaled <- fit_preprocessor(baseline, ~ scale(x, scale = FALSE, center = FALSE))
  expect_equal(as.vector(apply_preprocessor(new, unscaled$preprocessor)$x), new$x)
})

test_that("polynomial and spline coefficients and knots are frozen", {
  baseline <- data.frame(x = seq(0, 12, length.out = 25))
  new <- data.frame(x = c(0.25, 7.3, 11.5))
  fit <- fit_preprocessor(baseline, ~ poly(x, 3) + splines::ns(x, df = 4) +
                           splines::bs(x, df = 5))
  expected <- cbind(stats::predict(stats::poly(baseline$x, 3), new$x),
                    stats::predict(splines::ns(baseline$x, df = 4), new$x),
                    stats::predict(splines::bs(baseline$x, df = 5), new$x))
  actual <- apply_preprocessor(new, fit$preprocessor)$x
  expect_equal(as.vector(actual), as.vector(expected), tolerance = 1e-12)
  expect_equal(as.vector(apply_preprocessor(new[1, , drop = FALSE], fit$preprocessor)$x),
               as.vector(actual[1, , drop = FALSE]), tolerance = 1e-12)
  nested <- fit_preprocessor(baseline, ~ I(poly(x, 2)^2) + scale(splines::ns(x, 3)))
  train_spline <- splines::ns(baseline$x, 3)
  expected_nested <- cbind(stats::predict(stats::poly(baseline$x, 2), new$x)^2,
                           scale(stats::predict(train_spline, new$x),
                                 center = colMeans(train_spline),
                                 scale = apply(train_spline, 2, stats::sd)))
  expect_equal(as.vector(apply_preprocessor(new, nested$preprocessor)$x),
               as.vector(expected_nested), tolerance = 1e-12)
  simple <- fit_preprocessor(baseline, ~ poly(x, 2, simple = TRUE))
  expect_equal(as.vector(apply_preprocessor(new, simple$preprocessor)$x),
               as.vector(stats::predict(stats::poly(baseline$x, 2), new$x)),
               tolerance = 1e-12)
  raw <- fit_preprocessor(baseline, ~ poly(x, 2, raw = TRUE))
  expect_equal(as.vector(apply_preprocessor(new, raw$preprocessor)$x),
               as.vector(cbind(new$x, new$x^2)))
  reordered <- fit_preprocessor(baseline, ~splines::ns(df = 4, x = x) +
                                 scale(center = FALSE, x = x))
  expected_reordered <- cbind(stats::predict(splines::ns(baseline$x, 4), new$x),
                              new$x / sqrt(sum(baseline$x^2) / (nrow(baseline) - 1)))
  expect_equal(as.vector(apply_preprocessor(new, reordered$preprocessor)$x),
               as.vector(expected_reordered), tolerance = 1e-12)
  explicit_null <- fit_preprocessor(baseline, ~poly(x, 2, coefs = NULL))
  expect_equal(as.vector(apply_preprocessor(new, explicit_null$preprocessor)$x),
               as.vector(stats::predict(stats::poly(baseline$x, 2), new$x)))
})

test_that("factors, interactions and contrasts retain training definitions", {
  previous <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(previous), add = TRUE)
  baseline <- data.frame(x = 1:6, group = rep(c("red", "blue", "green"), 2),
                         grade = ordered(rep(c("low", "high"), 3), levels = c("low", "high")),
                         flag = c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE))
  fit <- fit_preprocessor(baseline, ~x * group + grade + flag)
  expect_true(is.factor(fit$outcome$baseline$group))
  expect_true(is.ordered(fit$outcome$baseline$grade))
  options(contrasts = c("contr.treatment", "contr.helmert"))
  new <- baseline[c(2, 3), ]
  new$group <- factor(new$group, levels = rev(levels(fit$outcome$baseline$group)))
  new$grade <- as.character(new$grade)
  expect_equal(as.vector(apply_preprocessor(new, fit$preprocessor)$x),
               as.vector(fit$outcome$x[c(2, 3), , drop = FALSE]))
  custom <- matrix(c(2, -1, -1, 0, 1, -1), 3, 2,
                    dimnames = list(NULL, c("first", "second")))
  specified <- fit_preprocessor(baseline, ~ group,
                                model_options = list(contrasts = list(group = custom)))
  expect_identical(specified$preprocessor$outcome$contrasts$group,
                   `rownames<-`(custom, levels(specified$outcome$baseline$group)))
  expect_equal(as.vector(apply_preprocessor(baseline, specified$preprocessor)$x),
               as.vector(specified$outcome$x))
  full <- fit_preprocessor(baseline, ~ 0 + group)
  expect_equal(ncol(full$outcome$x), 3L)
  expect_true(all(rowSums(full$outcome$x) == 1))
})

test_that("intercept-only and empty prediction designs retain dimensions", {
  empty <- data.frame(row.names = paste0("s", seq_len(4)))
  fit <- fit_preprocessor(empty, ~1)
  expect_identical(dim(fit$outcome$x), c(4L, 0L))
  expect_identical(dim(apply_preprocessor(empty[1, , drop = FALSE], fit$preprocessor)$x),
                   c(1L, 0L))
  expect_identical(dim(apply_preprocessor(empty[FALSE, , drop = FALSE], fit$preprocessor)$x),
                   c(0L, 0L))
  expect_identical(dim(fit_preprocessor(matrix(numeric(), 4L, 0L))$outcome$x), c(4L, 0L))
  baseline <- data.frame(x = 1:6, g = factor(rep(c("A", "B"), 3)))
  fit <- fit_preprocessor(baseline, ~poly(x, 2) + scale(x) + g)
  zero <- apply_preprocessor(baseline[FALSE, , drop = FALSE], fit$preprocessor)$x
  expect_identical(dim(zero), c(0L, ncol(fit$outcome$x)))
  expect_identical(colnames(zero), colnames(fit$outcome$x))
  expect_identical(attr(zero, "assign"), attr(fit$outcome$x, "assign"))
  expect_identical(attr(zero, "contrasts"), attr(fit$outcome$x, "contrasts"))
  one <- fit_preprocessor(data.frame(x = 6), ~x)
  expect_identical(dim(one$outcome$x), c(1L, 1L))
  expect_error(fit_preprocessor(baseline[FALSE, ]), "training subject")
})

test_that("missing values, unseen levels and incompatible raw types fail", {
  baseline <- data.frame(x = 1:4, g = factor(c("A", "B", "A", "B")))
  fit <- fit_preprocessor(baseline)
  expect_error(apply_preprocessor(data.frame(x = 1, g = "C"), fit$preprocessor),
               "unseen level.*C")
  expect_error(apply_preprocessor(data.frame(x = NA_real_, g = "A"), fit$preprocessor),
               "finite, nonmissing")
  expect_error(apply_preprocessor(data.frame(x = Inf, g = "A"), fit$preprocessor),
               "finite, nonmissing")
  expect_error(apply_preprocessor(data.frame(x = "1", g = "A"), fit$preprocessor),
               "must be numeric")
  expect_error(apply_preprocessor(data.frame(x = 1, g = 1), fit$preprocessor),
               "character or factor")
  expect_error(fit_preprocessor(data.frame(x = c(1, NA))), "finite, nonmissing")
  expect_error(fit_preprocessor(data.frame(x = c(1, 0)), ~log(x)), "finite")
  expect_error(fit_preprocessor(baseline, ~I(1)), "one row per subject")
  expect_error(fit_preprocessor(baseline, ~I(c(x, x))), "one row per subject")
  expect_error(fit_preprocessor(data.frame(x = c(1, 1)), ~scale(x)), "missing|finite")
  expect_error(fit_preprocessor(data.frame(g = c("A", "A"))), "two declared levels")
  expect_error(fit_preprocessor(data.frame(x = Sys.Date())), "must be numeric")
  expect_error(fit_preprocessor(data.frame(x = 1:4), ~x,
                                model_options = list(unknown = 1)), "only contrasts")
  expect_error(fit_preprocessor(baseline, model_options = list(contrasts = list(x = "contr.sum"))),
               "noncategorical")
})

test_that("unsupported sample-dependent or external transformations are explicit", {
  baseline <- data.frame(x = 1:5, y = 6:10)
  k <- 2
  custom <- function(x) x - mean(x)
  expect_error(fit_preprocessor(baseline, ~I(x - mean(x))), "Unsupported formula")
  expect_error(fit_preprocessor(baseline, ~custom(x)), "Unsupported formula")
  expect_error(fit_preprocessor(baseline, ~I(x/k)), "constants must be literals")
  expect_error(fit_preprocessor(baseline, ~stats::rnorm(x)), "Unsupported formula")
  expect_error(fit_preprocessor(baseline, ~factor(x)), "Unsupported formula")
  expect_error(fit_preprocessor(baseline, ~poly(x, y, degree = 2)), "literal constants")
  expect_error(fit_preprocessor(baseline, ~scale(x, center = y)), "literal constants")
  expect_error(fit_preprocessor(baseline, ~.), "explicitly")
  log <- function(x) rep(-999, length(x))
  fit <- fit_preprocessor(baseline, ~log(x))
  expect_equal(as.vector(fit$outcome$x), base::log(baseline$x))
  expect_identical(environment(fit$preprocessor$outcome$terms), baseenv())
})

test_that("saved preprocessing is independent of training scope and options", {
  fit <- local({
    accidental_training_state <- new.env(parent = emptyenv())
    accidental_training_state$large_data <- seq_len(100000)
    baseline <- data.frame(x = seq(1, 10, length.out = 30),
                           group = factor(rep(c("A", "B", "C"), 10)))
    fit_preprocessor(baseline, ~poly(x, 2) + scale(x) + splines::ns(x, 3) + group)
  })
  new <- data.frame(x = c(2.5, 8), group = c("B", "C"))
  before <- apply_preprocessor(new, fit$preprocessor)
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file), add = TRUE)
  saveRDS(fit$preprocessor, file)
  after <- readRDS(file)
  expect_identical(apply_preprocessor(new, after), before)
  expect_lt(as.numeric(object.size(after)), 100000)
  expect_identical(environment(after$outcome$terms), baseenv())
  expect_null(after$training_data)
})
