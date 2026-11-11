test_that("LifeHist preserves either response representation", {
  a <- LifeHist(c(20, 35, 80), c(1, 2, 0))
  expect_s3_class(a, "life_response")
  expect_identical(a$format, "event")
  expect_identical(a$data$time, c(20, 35, 80))
  b <- LifeHist(start = c(0, 20, 35), stop = c(20, 35, 80),
                from = c(0, 1, 0), to = c(1, 0, NA), event = c(1, 1, 0))
  expect_identical(b$format, "transition")
  expect_identical(b$data$from, c(0, 1, 0))
  expect_error(LifeHist(), "either")
  expect_error(LifeHist(1, 0, start = 0), "either")
  expect_error(LifeHist(time = 1), "status")
  expect_error(LifeHist(1:2, 0), "same positive length")
  expect_error(LifeHist(matrix(1:2), c(0, 0)), "vectors")
})

test_that("formula and vector extraction preserve identical selected baselines", {
  d <- data.frame(id = c("a", "a", "b"), time = c(2, 4, 7),
                  status = c(1, 0, 0), age = c(60, 60, 70),
                  group = factor(c("B", "B", "A")), unused = c(NA, 1, 2))
  a <- resolve_input(LifeHist(time, status) ~ age + group, data = d)
  b <- resolve_input(response = LifeHist(d$time, d$status), id = d$id,
                      x = d[c("age", "group")])
  expect_identical(a$response, b$response)
  expect_identical(a$baseline, b$baseline)
  expect_identical(a$id, b$id)
  expect_identical(a$schema$predictor_names, c("age", "group"))
  expect_s3_class(a$baseline$group, "factor")
  expect_error(resolve_input(LifeHist(time, status) ~ ., data = d), "explicitly")
  expect_error(resolve_input(LifeHist(time, status) ~ age, data = d,
                              x = d["age"]), "combine")
  expect_error(resolve_input(LifeHist(time, status) ~ other, data = d), "other")
  expect_error(resolve_input(response = LifeHist(d$time, d$status), id = "a"),
               "one value per record")
})

test_that("intercept-only and matrix input have explicit baseline dimensions", {
  d <- data.frame(id = 1:2, time = 1:2, status = 0)
  a <- resolve_input(LifeHist(time, status) ~ 1, data = d)
  expect_identical(dim(a$baseline), c(2L, 0L))
  b <- resolve_input(response = LifeHist(d$time, d$status), id = d$id,
                      x = matrix(c(30, 40), ncol = 1))
  expect_identical(names(b$baseline), "x1")
  expect_error(resolve_input(response = LifeHist(d$time, d$status), id = d$id,
                              x = matrix(1, ncol = 1)), "one row per record")
})

test_that("required values fail rather than being omitted", {
  response <- LifeHist(c(1, 2), c(0, 0))
  expect_error(resolve_input(response = response, id = c(1, NA)), "nonmissing")
  expect_error(resolve_input(response = response, id = c("a", "")), "nonempty")
  expect_error(resolve_input(response = response, id = 1:2,
                              x = data.frame(age = c(30, NA))), "age")
  expect_error(resolve_input(response = response, id = 1:2,
                              x = data.frame(age = c(30, Inf))), "age")
})

test_that("a stored formula mapping is reusable on another sample", {
  d <- data.frame(id = 1:2, time = 1:2, status = 0, age = c(40, 50))
  a <- resolve_input(LifeHist(time, status) ~ age, data = d)
  d$id <- 3:4
  b <- resolve_input(data = d, schema = a$schema)
  expect_identical(b$id, 3:4)
  expect_identical(b$baseline, d["age"])
  vector_input <- resolve_input(response = LifeHist(d$time, d$status), id = d$id,
                                 x = d["age"], schema = a$schema)
  expect_identical(vector_input$baseline, b$baseline)
  expect_error(resolve_input(response = LifeHist(d$time, d$status), id = d$id,
                              x = data.frame(other = c(40, 50)), schema = a$schema),
               "stored predictor schema")
})
