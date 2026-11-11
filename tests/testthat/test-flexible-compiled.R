test_that("compiled numeric trees preserve backend split directions and equalities", {
  skip_if_not_installed("rpart")
  x <- seq_len(16)
  make_tree <- function(event) {
    frame <- data.frame(f1 = x, exposure = 1, event = event)
    model <- rpart::rpart(cbind(exposure, event) ~ f1, data = frame,
      method = "poisson", control = rpart::rpart.control(minsplit = 2,
        minbucket = 1, cp = 0, maxdepth = 1, xval = 0,
        maxsurrogate = 0, usesurrogate = 0))
    list(constant = NULL, model = model, columns = "f1")
  }
  trees <- list(make_tree(rep(c(1, 9), each = 8)),
    make_tree(rep(c(9, 1), each = 8)), make_tree(rep(3, 16)),
    list(constant = 0, model = NULL, columns = "f1"))
  directions <- vapply(trees[1:2], function(x) x$model$splits[1, "ncat"], numeric(1))
  expect_setequal(directions, c(-1, 1))
  cuts <- vapply(trees[1:2], function(x) x$model$splits[1, "index"], numeric(1))
  values <- unique(c(-1e100, x, -cdt_nextafter(-cuts), cuts, cdt_nextafter(cuts), 1e100))
  features <- matrix(values, ncol = 1, dimnames = list(NULL, "f1"))
  reference <- numeric(length(values))
  for (tree in trees) reference <- reference + if (is.null(tree$constant))
    as.numeric(predict(tree$model, as.data.frame(features), type = "vector")) else tree$constant
  reference <- reference / length(trees)
  compiled <- compile_flexible_forest(trees, "f1")
  expect_identical(.Call(life_forest_rates, compiled, features, seq_along(values)), reference)
  rows <- c(length(values), 2L, 2L, 1L)
  expect_identical(.Call(life_forest_rates, compiled, features, rows), reference[rows])
  expect_identical(.Call(life_forest_rates, compiled, features, integer()), numeric())
  expect_identical(.Call(life_forest_rates, unserialize(serialize(compiled, NULL)),
                        features, seq_along(values)), reference)
  cache <- new.env(parent = emptyenv())
  cache$nodes <- matrix(0L, length(values), length(trees))
  cached <- function(f, rows = seq_along(values), path = seq_along(values))
    .Call(life_forest_rates_cached, compiled, f, rows, cache, path)
  expect_identical(cached(features), reference)
  expect_true(all(cache$nodes > 0L))
  expect_identical(cached(features), reference)
  # Reuse the same draw slots with values on both sides of every cut, including
  # exact equality and adjacent representable values. Invalid leaves must miss.
  reverse <- rev(seq_along(values))
  expect_identical(cached(features[reverse, , drop = FALSE]), reference[reverse])
  expect_identical(cached(features, rows), reference[rows])
  old_leaf <- cache$nodes[1L, 1L]
  cache$nodes[1L, 1L] <- -1L
  expect_error(cached(features), "Invalid leaf")
  cache$nodes[1L, 1L] <- old_leaf
  bad_box <- compiled; bad_box$box_start[old_leaf] <- -1L
  expect_error(.Call(life_forest_rates_cached, bad_box, features,
    seq_along(values), cache, seq_along(values)), "leaf bounds")
  bad_box <- compiled; bad_box$box_count[old_leaf] <- -1L
  expect_error(.Call(life_forest_rates_cached, bad_box, features,
    seq_along(values), cache, seq_along(values)), "leaf bounds")
  bad_cache <- new.env(parent = emptyenv())
  bad_cache$nodes <- matrix(0L, length(values), length(trees) + 1L)
  expect_error(.Call(life_forest_rates_cached, compiled, features,
    seq_along(values), bad_cache, seq_along(values)), "cache dimensions")
})

test_that("compiled forests reject malformed traversal state safely", {
  leaf <- compile_flexible_forest(list(list(constant = 2)), c("f1", "f2"))
  features <- matrix(1:8, 4, 2)
  expect_identical(.Call(life_forest_rates, leaf, features, 1:4), rep(2, 4))
  expect_error(.Call(life_forest_rates, leaf, features, NA_integer_), "row index")
  expect_error(.Call(life_forest_rates, leaf, features, 0L), "row index")
  expect_error(.Call(life_forest_rates, leaf, features, 5L), "row index")
  bad <- leaf; bad$roots <- 2L
  expect_error(.Call(life_forest_rates, bad, features, 1L), "child index")
  bad <- leaf; bad$value <- Inf
  expect_error(.Call(life_forest_rates, bad, features, 1L), "terminal.*rate")
  bad <- leaf; bad$feature <- 1L; bad$direction <- -1L
  bad$cut <- 0; bad$left <- bad$right <- 1L
  expect_error(.Call(life_forest_rates, bad, features, 1L), "cycle")
  expect_error(.Call(life_forest_rates, leaf, features[, 1, drop = FALSE], 1L),
               "features.*match")
})

test_that("neural optimizer callbacks share exact calculations without changing fits", {
  design <- cbind(1, seq(-1, 1, length.out = 8), rep(c(-.5, .5), 4))
  exposure <- seq(.1, .8, length.out = 8)
  count <- c(0, 1, 0, 0, 1, 0, 0, 1)
  par <- seq(-.4, .4, length.out = 9)
  evaluator <- flex_neural_evaluator(design, exposure, count, 2L, .03, 8)
  objective <- function(p, gradient = FALSE)
    flex_neural_objective(p, design, exposure, count, 2L, .03, gradient, normalizer = 8)
  expect_identical(evaluator$fn(par), objective(par))
  expect_identical(evaluator$gr(par), objective(par, TRUE))
  expect_identical(evaluator$count(), 1L)
  changed <- par; changed[1L] <- changed[1L] + .01
  expect_identical(evaluator$gr(changed), objective(changed, TRUE))
  expect_identical(evaluator$fn(changed), objective(changed))
  expect_identical(evaluator$count(), 2L)
  expect_identical(evaluator$fn(par), objective(par))
  expect_identical(evaluator$count(), 3L)
  reference <- stats::optim(par, objective, gr = function(p) objective(p, TRUE),
    method = "BFGS", control = list(maxit = 1000, reltol = 1e-8))
  cached <- stats::optim(par, evaluator$fn, gr = evaluator$gr,
    method = "BFGS", control = list(maxit = 1000, reltol = 1e-8))
  expect_identical(cached, reference)
  expect_lt(evaluator$count() - 3L, sum(cached$counts))
})

compiled_forest_fixture <- function() {
  n <- 24L
  first <- .4 + seq_len(n) / 100
  second <- first + .3 + seq_len(n) %% 3 / 10
  last <- 2 + seq_len(n) %% 4 / 10
  data <- data.frame(id = rep(seq_len(n), each = 3L),
    start = c(rbind(0, first, second)), stop = c(rbind(first, second, last)),
    from = rep(c(0, 1, 0), n), to = rep(c(1, 0, 2), n), event = 1,
    z = rep(seq_len(n) / n, each = 3L), w = rep(rep(c(-1, 1), n / 2), each = 3L))
  life_train(LifeHist(start = start, stop = stop, from = from, to = to,
    event = event) ~ z + w, data, absorbing = 2,
    targets = list(home = 0, alive = c(0, 1)), tau = 3,
    outcome_model = "transition_forest", outcome_options = list(
      num_trees = 8, step_size = .25, minsplit = 4, minbucket = 2, max_depth = 4),
    control = list(draws = 37), seed = 103)
}

test_that("fitted compiled forests match saved-backend rate and continuation predictions", {
  skip_if_not_installed("rpart")
  fit <- compiled_forest_fixture()
  adapter <- fit$outcome
  legacy <- adapter
  for (e in seq_along(legacy$fits)) legacy$fits[[e]]$compiled <- NULL
  n <- 101L
  x <- cbind(z = seq(-2, 3, length.out = n), w = rep(c(-1, 1), length.out = n))
  time <- seq(0, 3, length.out = n)
  state <- rep(1:3, length.out = n)
  duration <- time / 2
  previous <- rep(0:4, length.out = n)
  accrued <- cbind(time / 3, time)
  reference <- predict_transition_rates(legacy, x, time, state, duration, previous, accrued)
  # Additive split tables preserve the ensemble while regrouping nonnegative
  # floating-point sums; allow only roundoff against the backend tree order.
  expect_equal(predict_transition_rates(adapter, x, time, state, duration, previous, accrued),
                   reference, tolerance = 1e-13)
  restored <- unserialize(serialize(adapter, NULL))
  expect_identical(predict_transition_rates(restored, x, time, state, duration, previous, accrued),
    predict_transition_rates(adapter, x, time, state, duration, previous, accrued))
  take <- c(1L, 7L, 18L, 65L)
  snapshots <- list(time = time[take], state = state[take], duration = duration[take],
    previous_transitions = previous[take], accrued = accrued[take, , drop = FALSE])
  set.seed(619); before <- .Random.seed
  a <- simulate_flexible_targets(adapter, x[take, , drop = FALSE], snapshots, B = 37)
  b <- simulate_flexible_targets(legacy, x[take, , drop = FALSE], snapshots, B = 37)
  expect_equal(a, b, tolerance = 1e-12)
  expect_identical(.Random.seed, before)
  expect_identical(simulate_flexible_targets(restored, x[take, , drop = FALSE], snapshots, B = 37), a)
  # A saved version without native traversal arrays retains its backend route.
  old_fit <- fit; old_fit$outcome <- old_fit$augmentation <- legacy
  profiles <- data.frame(z = c(.1, .5, .9), w = c(-1, 1, -1))
  expect_equal(predict(fit, profiles, initial_state = 0)$mean,
                   predict(old_fit, profiles, initial_state = 0)$mean, tolerance = 1e-12)
  expect_equal(predict(fit, profiles, initial_state = 0)$variance,
                   predict(old_fit, profiles, initial_state = 0)$variance, tolerance = 1e-12)
})


test_that("native neural standardization preserves exact sweep arithmetic", {
  f <- matrix(seq(-1e6, 1e6, length.out = 70), 10, 7)
  center <- seq(-.3, .8, length.out = 7)
  scale <- c(1, .001, 3, 1e-5, 1e5, 7, .9)
  rows <- c(10L, 2L, 2L, 1L, 7L)
  expected <- cbind(1, sweep(sweep(f[rows, , drop = FALSE], 2L, center), 2L, scale, "/"))
  expect_identical(.Call(life_flexible_design, f, rows, center, scale), expected)
  f <- matrix(1:70, 10, 7)
  expected <- cbind(1, sweep(sweep(f[rows, , drop = FALSE], 2L, center), 2L, scale, "/"))
  expect_identical(.Call(life_flexible_design, f, rows, center, scale), expected)
  expect_identical(dim(.Call(life_flexible_design, f, integer(), center, scale)), c(0L, 8L))
  expect_error(.Call(life_flexible_design, f, NA_integer_, center, scale), "row index")
  expect_error(.Call(life_flexible_design, f, rows, center[-1L], scale), "dimensions")
  scale[1L] <- 0
  expect_error(.Call(life_flexible_design, f, rows, center, scale), "standardization")
})


test_that("simulation private evaluation retains generated-feature overflow guards", {
  skip_if_not_installed("rpart")
  fit <- compiled_forest_fixture()$outcome
  expect_error(predict_transition_rates(fit, matrix(c(.1, 1), 1,
    dimnames = list(NULL, c("z", "w"))), 1, 1L, .5, 0L,
    matrix(c(Inf, 1), 1), .validated = TRUE), "history features")
})


test_that("large integer budgets disable optional leaf caches before allocating", {
  skip_if_not_installed("rpart")
  adapter <- compiled_forest_fixture()$outcome
  for (e in seq_along(adapter$fits)) if (!is.null(adapter$fits[[e]]$compiled))
    adapter$fits[[e]]$compiled$roots <- rep(adapter$fits[[e]]$compiled$roots[1L], 4000L)
  simulator <- simulate_flexible_targets
  context <- new.env(parent = environment(simulator))
  # Stop at the existing uniform-bank allocation. This exercises the budget
  # decision without allocating a huge bank or fitting 4,000 separate trees.
  context$matrix <- function(data, nrow, ncol, ...) {
    context$dimensions <- c(nrow, ncol)
    stop("uniform bank allocation boundary", call. = FALSE)
  }
  environment(simulator) <- context
  x <- matrix(c(.1, 1), 1, dimnames = list(NULL, c("z", "w")))
  snapshots <- list(time = 0, state = 1L, duration = 0,
    previous_transitions = 0L, accrued = matrix(0, 1, 2))
  expect_error(simulator(adapter, x, snapshots, B = 1000000L),
    "uniform bank allocation boundary")
  expect_identical(context$dimensions, c(1000000L, 64L))
})

test_that("native leaf caches require a local binding and evaluate it safely", {
  compiled <- compile_flexible_forest(list(list(constant=2)),c("f1","f2"))
  x <- matrix(1:8,4,2)
  rate <- function(cache) .Call(life_forest_rates_cached,compiled,x,1:4,cache,1:4)
  local <- new.env(parent=emptyenv()); local$nodes <- matrix(0L,4,1)
  expect_identical(rate(local),rep(2,4))
  expect_identical(local$nodes,matrix(1L,4,1))
  expect_identical(rate(local),rep(2,4))
  missing <- new.env(parent=emptyenv())
  expect_error(rate(missing),"missing local nodes")
  inherited <- new.env(parent=local)
  expect_error(rate(inherited),"missing local nodes")
  delayed <- new.env(parent=emptyenv()); counter <- new.env(parent=emptyenv()); counter$n <- 0L
  delayedAssign("nodes",{counter$n <- counter$n+1L; matrix(0L,4,1)},
    eval.env=environment(),assign.env=delayed)
  expect_identical(rate(delayed),rep(2,4))
  expect_identical(rate(delayed),rep(2,4))
  expect_identical(counter$n,1L)
  active <- new.env(parent=emptyenv())
  makeActiveBinding("nodes",function() matrix(0L,4,1),active)
  expect_identical(rate(active),rep(2,4))
  invalid <- new.env(parent=emptyenv()); invalid$nodes <- NULL
  expect_error(rate(invalid),"cache dimensions")
})
