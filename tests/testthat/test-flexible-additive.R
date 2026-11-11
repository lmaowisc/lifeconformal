additive_forest_fixture <- function() {
  m <- list(roots=c(1L,4L,7L,10L),
    feature=c(1L,0L,0L,2L,0L,0L,1L,0L,0L,0L),
    left=c(2L,0L,0L,5L,0L,0L,8L,0L,0L,0L),
    right=c(3L,0L,0L,6L,0L,0L,9L,0L,0L,0L),
    direction=c(-1L,0L,0L,1L,0L,0L,1L,0L,0L,0L),
    cut=c(0,0,0,1,0,0,0,0,0,0),
    value=c(0,.1,.3,0,.7,.9,0,1e-10,1e-11,.6), n_features=3L,
    tree_start=c(1L,4L,7L,10L),tree_end=c(3L,6L,9L,10L))
  m$lookup <- compile_forest_lookup(m)
  m$additive <- compile_forest_additive(m)
  m
}

test_that("additive tables preserve raw cuts and the finite ensemble", {
  m <- additive_forest_fixture()
  z <- sort(unique(c(-1e100, 0, 1, 1e100,
    -cdt_nextafter(-c(0,1)), cdt_nextafter(c(0,1)))))
  x <- as.matrix(expand.grid(z,z,z)); storage.mode(x) <- "double"
  expected <- (ifelse(x[,1] < 0, .1, .3) +
    ifelse(x[,2] < 1, .9, .7) + ifelse(x[,1] < 0, 1e-11, 1e-10) + .6) / 4
  rates <- function(model, features = x, rows = seq_len(nrow(features)))
    .Call(life_forest_rates, model, features, as.integer(rows))
  expect_length(m$additive$feature1, 3L)
  expect_equal(rates(m), expected, tolerance = 1e-14)
  rows <- c(20:1,1L,1L,5L)
  expect_equal(rates(m,rows=rows), expected[rows], tolerance = 1e-14)
  expect_identical(rates(m,rows=integer()), numeric())
  expect_identical(rates(m,rows=1L), rates(m)[1L])
  expect_identical(rates(unserialize(serialize(m,NULL))),rates(m))
  expect_identical(rates(m,x[rev(seq_len(nrow(x))),,drop=FALSE]),rev(rates(m)))
  for (scale in c(1e-200,1e200)) {
    scaled <- m; scaled$value <- m$value*scale
    scaled$additive <- compile_forest_additive(scaled)
    reference <- scaled; reference$lookup <- reference$additive <- NULL
    actual <- rates(scaled); direct <- rates(reference)
    expect_true(all(abs(actual-direct) <=
      8*length(m$roots)*.Machine$double.eps*pmax(abs(direct),.Machine$double.xmin)))
  }
})

test_that("optional additive table budgets fall back safely", {
  m <- additive_forest_fixture()
  expect_null(compile_forest_additive(m,max_bytes=1))
  expect_null(compile_forest_additive(m,max_cells=1L))
  expect_lte(as.numeric(object.size(m$additive)),16*1024^2)
  absent <- m; absent$lookup <- NULL
  expect_null(compile_forest_additive(absent))
  incomplete <- m; incomplete$lookup$complete <- FALSE
  expect_null(compile_forest_additive(incomplete))
  tree <- list(constant=2)
  no_tables <- compile_flexible_forest(list(tree),c("f1","f2"),table_budget=0)
  expect_null(no_tables$lookup)
  expect_null(no_tables$additive)
  expect_identical(.Call(life_forest_rates,no_tables,matrix(1:8,4,2),1:4),rep(2,4))
})

test_that("additive tables reject unsafe dimensions and values", {
  m <- additive_forest_fixture(); x <- matrix(c(-1,0,2,1,2,3),2,3,byrow=TRUE)
  rate <- function(model) .Call(life_forest_rates,model,x,1:2)
  bad <- m; bad$additive <- 1L
  expect_error(rate(bad),"Invalid compiled")
  bad <- m; bad$additive$size1 <- integer()
  expect_error(rate(bad),"dimensions")
  bad <- m; bad$additive$trees <- NA_integer_
  expect_error(rate(bad),"dimensions")
  for (field in c("feature1","feature2","size1","size2","offset")) {
    bad <- m; bad$additive[[field]][1] <- .Machine$integer.max
    expect_error(rate(bad),"table")
    bad$additive[[field]][1] <- NA_integer_
    expect_error(rate(bad),"table")
  }
  for (value in c(Inf,NaN,-1)) {
    bad <- m; bad$additive$value[] <- value
    expect_error(rate(bad),"rate")
  }
})
