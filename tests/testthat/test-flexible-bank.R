forest_bank_fixture <- function() {
  structure(list(type = "transition_forest", n_predictors = 1L, x_names = "x", seed = 37L,
    fits = list(list(constant = .25), list(constant = .25)),
    edges = data.frame(from = c(1L, 2L), to = c(2L, 1L)),
    problem = list(tau = 8, rewards = matrix(c(1, 0, 1, 1), 2L, 2L,
      dimnames = list(NULL, c("before", "total"))), absorbing = integer()),
    options = list(step_size = .125, max_events = 100L)), class = "life_flexible_outcome")
}

forest_bank_case <- function(B = 3L, block = 1L, adapter = forest_bank_fixture(), audit = 1L) {
  precision <- flex_forest_block_precision()
  skip_if(is.null(precision), "Runtime accumulator format uses the R fallback")
  group <- rep(1:2, each = B)
  token <- .Call(life_forest_bank_create, adapter, group, as.integer(B), precision,
    as.integer(block), 64 * 1024^2, audit)
  finished <- FALSE
  on.exit(if (!finished) .Call(life_forest_bank_release, token), add = TRUE)
  expected <- matrix(0, B * 2L, 64L)
  for (ii in 1:2) {
    seed <- flex_profile_seed(100L + ii, block)
    old <- flex_with_seed(seed, {
      values <- matrix(stats::runif(B * 64L), B, byrow = TRUE)
      list(values = values, seed = .Random.seed, kind = RNGkind())
    })
    actual <- flex_with_seed(seed, {
      .Call(life_forest_bank_fill, token, ii)
      list(seed = .Random.seed, kind = RNGkind())
    })
    expect_identical(actual$seed, old$seed)
    expect_identical(actual$kind, old$kind)
    expected[(ii - 1L) * B + seq_len(B), ] <- old$values
  }
  .Call(life_forest_bank_seal, token)
  finished <- TRUE
  list(adapter = adapter, group = group, B = B, token = token, bank = expected,
    x = matrix(c(0L, 1L), 2L, 1L),
    path = flex_forest_block_state(rep(.25, B * 2L), rep(1:2, each = B),
      rep(.125, B * 2L), rep(1, B * 2L), rep(.375, B * 2L), rep(0, B * 2L),
      cbind(rep(.125, B * 2L), rep(.25, B * 2L)), integer(B * 2L), (block - 1L) * 32L),
    controls = list(draws = as.integer(B), precision = precision, steps = 32L,
      max_iter = 1000L, trace = as.integer(B <= 3L)))
}

test_that("owned banks preserve full draws, post-fill RNG and native traces", {
  withr::local_seed(849)
  rng <- .Random.seed; kind <- RNGkind()
  for (B in c(3L, 1000L)) {
    z <- forest_bank_case(B)
    on.exit(.Call(life_forest_bank_release, z$token), add = TRUE)
    info <- .Call(life_forest_bank_inspect, z$token)
    expect_equal(info[["completed"]], 2)
    expect_lte(info[["reserved_payload"]], 64 * 1024^2)
    expect_gte(info[["reserved_payload"]], 2 * length(z$bank) * 8)
    copy <- .Call(life_forest_bank_copy, z$token)
    expect_true(identical(copy, z$bank, num.eq = FALSE))
    expect_error(.Call(life_forest_bank_copy, z$token), "already consumed")
    copy[1L] <- NA_real_
    stored <- serialize(list(z$adapter, z$x, z$group, z$path), NULL, version = 2)
    expected <- .Call(life_forest_block_advance, z$adapter, z$x, z$group,
      z$path, z$bank, z$controls)
    gc()
    actual <- .Call(life_forest_bank_advance, z$adapter, z$x, z$group,
      z$path, z$token, z$controls)
    expect_true(identical(actual, expected, num.eq = FALSE))
    expect_identical(serialize(list(z$adapter, z$x, z$group, z$path), NULL, version = 2), stored)
    expect_identical(.Call(life_forest_bank_live), 0L)
    expect_error(.Call(life_forest_bank_advance, z$adapter, z$x, z$group,
      z$path, z$token, z$controls), "Released")
    expect_identical(.Random.seed, rng)
    expect_identical(RNGkind(), kind)
  }
})

test_that("bank identity and path failures consume their owner", {
  for (which in c("group", "block", "draws", "model", "path", "dimensions")) {
    z <- forest_bank_case(1L)
    on.exit(.Call(life_forest_bank_release, z$token), add = TRUE)
    if (which == "group") z$group <- c(1L, 3L)
    if (which == "block") z$path$iteration <- 32L
    if (which == "draws") z$controls$draws <- 2L
    if (which == "model") z$adapter <- unserialize(serialize(z$adapter, NULL))
    if (which == "path") z$path$count[1L] <- 1L
    if (which == "dimensions") z$x <- matrix(0, 2L, 2L)
    expect_error(.Call(life_forest_bank_advance, z$adapter, z$x, z$group,
      z$path, z$token, z$controls), "bank|continuation")
    expect_identical(.Call(life_forest_bank_live), 0L)
    expect_error(.Call(life_forest_bank_seal, z$token), "Released")
  }
  z <- forest_bank_case(1L)
  restored <- unserialize(serialize(z$token, NULL, version = 2))
  expect_error(.Call(life_forest_bank_advance, z$adapter, z$x, z$group,
    z$path, restored, z$controls), "Released")
  .Call(life_forest_bank_release, z$token)
  expect_error(.Call(life_forest_bank_advance, z$adapter, z$x, z$group,
    z$path, z$bank, z$controls), "Foreign")
})

test_that("incomplete banks, helper errors and unreachable handles release safely", {
  adapter <- forest_bank_fixture(); precision <- flex_forest_block_precision()
  skip_if(is.null(precision), "Runtime accumulator format uses the R fallback")
  token <- .Call(life_forest_bank_create, adapter, 1L, 1L, precision, 1L, 64 * 1024^2, 0L)
  expect_error(.Call(life_forest_bank_seal, token), "not complete")
  expect_error(.Call(life_forest_bank_copy, token), "unavailable")
  alias <- token
  .Call(life_forest_bank_release, token)
  .Call(life_forest_bank_release, alias)
  expect_error(.Call(life_forest_bank_fill, alias, 1L), "Released")
  withr::local_seed(739); rng <- .Random.seed
  expect_error(flex_forest_bank_generate(adapter, c(1L, 2L), 1L, precision,
    c(101L, NA_integer_), 1L), "seed")
  expect_identical(.Random.seed, rng)
  expect_identical(.Call(life_forest_bank_live), 0L)
  invisible(local({
    abandoned <- .Call(life_forest_bank_create, adapter, 1L, 1L, precision, 1L, 64 * 1024^2, 0L)
    expect_identical(.Call(life_forest_bank_live), 1L)
    NULL
  }))
  gc()
  expect_identical(.Call(life_forest_bank_live), 0L)
})

test_that("generation restores an absent seed and the caller's RNG kind", {
  adapter <- forest_bank_fixture(); precision <- flex_forest_block_precision()
  skip_if(is.null(precision), "Runtime accumulator format uses the R fallback")
  old_kind <- RNGkind(); had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv) else
      if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  rm(".Random.seed", envir = .GlobalEnv)
  kind <- RNGkind()
  token <- flex_forest_bank_generate(adapter, c(1L, 2L), 1L, precision, c(101L, 102L), 1L)
  .Call(life_forest_bank_release, token)
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
  expect_identical(RNGkind(), kind)
})

test_that("resource fallback preserves complete matrix-route targets and seeds", {
  adapter <- forest_bank_fixture(); precision <- flex_forest_block_precision()
  skip_if(is.null(precision), "Runtime accumulator format uses the R fallback")
  expect_false(flex_forest_bank_supported(adapter, rep(1L, 1001L), 1001L, precision))
  expect_false(flex_forest_bank_supported(adapter, 1:17, 1L, precision))
  expect_false(flex_forest_bank_supported(adapter, 1L, 1L, precision, budget = 0))
  expect_error(.Call(life_forest_bank_create, adapter, 1L, 1L, precision, 1L, 0, 0L), "unsupported")
  x <- matrix(rep(0:1, length.out = 33L), 33L, 1L, dimnames = list(NULL, "x"))
  snapshot <- list(time = rep(.25, 33L), state = rep(1:2, length.out = 33L),
    duration = rep(.125, 33L), previous_transitions = integer(33L),
    accrued = cbind(rep(.125, 33L), rep(.25, 33L)))
  withr::local_seed(739); rng <- .Random.seed; stored <- serialize(adapter, NULL, version = 2)
  generate <- flex_forest_bank_generate; calls <- 0L
  actual <- with_mocked_bindings(simulate_flexible_targets(adapter, x, snapshot, B = 7L),
    flex_forest_bank_generate = function(...) { calls <<- calls + 1L; generate(...) },
    .package = "lifeconformal")
  expected <- with_mocked_bindings(simulate_flexible_targets(adapter, x, snapshot, B = 7L),
    flex_forest_bank_supported = function(...) FALSE, .package = "lifeconformal")
  expect_gt(calls, 0L)
  expect_true(identical(actual, expected, num.eq = FALSE))
  expect_identical(serialize(adapter, NULL, version = 2), stored)
  expect_identical(.Random.seed, rng)
  expect_identical(.Call(life_forest_bank_live), 0L)
  # An explicit larger draw request uses the original route, never a smaller B.
  adapter$problem$tau <- .5
  one <- lapply(snapshot, function(z) if (is.matrix(z)) z[1L, , drop = FALSE] else z[1L])
  oversized <- with_mocked_bindings(simulate_flexible_targets(adapter, x[1L, , drop = FALSE], one, B = 1001L),
    flex_forest_bank_generate = function(...) stop("unexpected generated bank"), .package = "lifeconformal")
  expect_identical(dim(oversized), c(1L, 1001L, 2L))
  expect_identical(.Random.seed, rng)
})
