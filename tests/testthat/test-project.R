projection_fixture <- function(lower, upper, targets = list(home = 0, hospital = 1, alive = c(0, 1)),
                                states = c(0, 1, 2), initial = 1L, tau = 10,
                                allowed = NULL) {
  d <- length(targets)
  lower <- matrix(lower, ncol = d, byrow = TRUE)
  upper <- matrix(upper, ncol = d, byrow = TRUE)
  n <- nrow(lower)
  support <- build_structural_support(states, absorbing = 2,
    allowed_transitions = allowed, targets = targets, tau = tau)
  center <- matrix(0, n, d, dimnames = list(NULL, names(targets)))
  bounds <- data.frame(id = rep(seq_len(n), each = d),
                       target = rep(names(targets), n),
                       lower = as.vector(t(lower)), upper = as.vector(t(upper)))
  structure(list(bounds = bounds, center = center, tau = tau, alpha = .1,
                 region = list(support = support, initial_state = rep(initial, length.out = n))),
            class = "lifeconformal_prediction")
}

test_that("sums use shared occupation constraints instead of marginal interval sums", {
  skip_if_not_installed("lpSolve")
  p <- projection_fixture(c(2, 3, 7), c(8, 9, 9))
  z <- project_life(p, c(home = 1, hospital = 1))
  expect_equal(z$lower, 7, tolerance = 1e-12)
  expect_equal(z$upper, 9, tolerance = 1e-12)
  expect_identical(z$status, "ok")
  expect_equal(attr(z, "alpha"), .1)
  expect_true(attr(z, "inherits_joint_coverage"))
  expect_true(attr(z, "conservative_outer_support"))
  expect_identical(colnames(attr(z, "weights")), c("home", "hospital", "alive"))
  expect_equal(as.numeric(attr(z, "weights")), c(1, 1, 0))
})

test_that("signed summaries preserve nesting and simultaneous interpretations", {
  skip_if_not_installed("lpSolve")
  p <- projection_fixture(c(2, 3, 7), c(8, 9, 9))
  w <- rbind(home_minus_hospital = c(home = 1, hospital = -1, alive = 0),
              alive_minus_home = c(home = -1, hospital = 0, alive = 1),
              zero = c(home = 0, hospital = 0, alive = 0))
  z <- project_life(p, w)
  expect_identical(z$summary, rownames(w))
  expect_equal(z$lower, c(-5, 3, 0), tolerance = 1e-12)
  expect_equal(z$upper, c(3, 7, 0), tolerance = 1e-12)
  expect_equal(attr(z, "diagnostics")$solver_calls, 5)
  w_reorder <- w[, c("alive", "hospital", "home"), drop = FALSE]
  expect_identical(project_life(p, w_reorder), z)
  nested <- projection_fixture(c(2, 7), c(8, 9), targets = list(home = 0, alive = c(0, 1)))
  contrast <- project_life(nested, c(alive = 1, home = -1))
  expect_equal(contrast$lower, 0, tolerance = 1e-12)
  expect_equal(contrast$upper, 7, tolerance = 1e-12)
})

test_that("full support and singleton states use analytic projections", {
  p <- projection_fixture(c(0, 0, 0), c(10, 10, 10))
  w <- rbind(alive = c(home = 1, hospital = 1), difference = c(home = -1, hospital = 1))
  z <- project_life(p, w)
  expect_equal(z$lower, c(0, -10))
  expect_equal(z$upper, c(10, 10))
  expect_equal(attr(z, "diagnostics")$solver_calls, 0L)
  absorbed <- projection_fixture(c(0, 10), c(0, 10),
                                  targets = list(alive = c(0, 1), dead = 2), initial = 3L)
  a <- project_life(absorbed, c(alive = 1, dead = -2))
  expect_equal(a$lower, -20)
  expect_equal(a$upper, -20)
  expect_equal(attr(a, "diagnostics")$solver_calls, 0L)
  restricted <- projection_fixture(c(0, 0, 0), c(10, 10, 10),
                                    allowed = data.frame(from = 0, to = 2))
  h <- project_life(restricted, c(hospital = 1))
  expect_equal(h$lower, 0)
  expect_equal(h$upper, 0)
})

test_that("empty joint intersections are explicit and never reversed intervals", {
  skip_if_not_installed("lpSolve")
  p <- projection_fixture(c(7, 7, 0), c(8, 8, 10))
  expect_warning(z <- project_life(p, c(home = 1, hospital = 1)),
                   class = "lifeconformal_empty_projection")
  expect_identical(z$status, "empty_region")
  expect_true(is.na(z$lower) && is.na(z$upper))
  p$bounds$lower[1] <- 9
  expect_warning(z <- project_life(p, c(home = 1)), class = "lifeconformal_empty_projection")
  expect_equal(attr(z, "diagnostics")$solver_calls, 0L)
})

test_that("projection shapes survive one target, one row and zero rows", {
  skip_if_not_installed("lpSolve")
  p <- projection_fixture(2, 6, targets = list(alive = c(0, 1)))
  z <- project_life(p, c(alive = 2))
  expect_equal(z$lower, 4, tolerance = 1e-12)
  expect_equal(z$upper, 12, tolerance = 1e-12)
  p$center <- p$center[FALSE, , drop = FALSE]
  p$bounds <- p$bounds[FALSE, ]
  p$region$initial_state <- integer()
  empty <- project_life(p, matrix(c(1, -1), 2, dimnames = list(c("pos", "neg"), "alive")))
  expect_s3_class(empty, "data.frame")
  expect_equal(nrow(empty), 0)
  expect_equal(dim(attr(empty, "weights")), c(2L, 1L))
  expect_equal(attr(empty, "alpha"), .1)
})

test_that("identical regions reuse LP solutions within a projection call", {
  skip_if_not_installed("lpSolve")
  p <- projection_fixture(rep(c(2, 3, 7), 20), rep(c(8, 9, 9), 20))
  z <- project_life(p, c(home = 1, hospital = 1))
  expect_equal(nrow(z), 20L)
  expect_equal(attr(z, "diagnostics")$unique_regions, 1L)
  expect_equal(attr(z, "diagnostics")$solver_calls, 3L)
  expect_equal(z$lower, rep(7, 20), tolerance = 1e-12)
  expect_identical(unserialize(serialize(z, NULL)), z)
})

test_that("weights and stored region alignment fail with useful errors", {
  p <- projection_fixture(c(0, 0, 0), c(10, 10, 10))
  expect_error(project_life(p, c(1, 2)), "target names")
  expect_error(project_life(p, c(home = 1, home = 2)), "unique")
  expect_error(project_life(p, c(unknown = 1)), "Unknown target")
  expect_error(project_life(p, c(home = NA_real_)), "finite numeric")
  expect_error(project_life(p, matrix(c(1, 2), 1)), "named target")
  expect_error(project_life(p, matrix(c(1, 2), 2,
                 dimnames = list(c("same", "same"), "home"))), "Summary row names")
  p$bounds$target[2] <- "home"
  expect_error(project_life(p, c(home = 1)), "complete target block")
})

test_that("the public toy-data workflow yields total-time projections", {
  skip_if_not_installed("lpSolve")
  dat <- life_example_data(80, seed = 11)
  fit <- life_train(LifeHist(time, status) ~ 1, dat[1:40, ],
                     absorbing = 2, targets = list(alive = 0, dead = 2), tau = 5)
  cal <- life_calibrate(fit, dat[41:80, ], method = "ipw")
  region <- predict(cal, data.frame(id = "new"))
  projected <- project_life(region, c(alive = 1, dead = 1))
  expect_equal(projected$lower, 5, tolerance = 1e-12)
  expect_equal(projected$upper, 5, tolerance = 1e-12)
  expect_identical(projected$id, "new")
})

test_that("LP projections match independent two-dimensional vertex enumeration", {
  skip_if_not_installed("lpSolve")
  set.seed(543)
  for (iteration in seq_len(20)) {
    point <- stats::runif(3)
    point <- point / sum(point)
    gamma <- c(point[1:2], sum(point[1:2]))
    lower <- pmax(0, gamma - stats::runif(3, .02, .3))
    upper <- pmin(1, gamma + stats::runif(3, .02, .3))
    p <- projection_fixture(10 * lower, 10 * upper)
    weights <- stats::setNames(sample(-3:3, 3, replace = TRUE), c("home", "hospital", "alive"))
    z <- project_life(p, weights)
    A2 <- rbind(c(1, 0), c(0, 1), c(1, 1))
    constraints <- rbind(A2, -A2, c(-1, 0), c(0, -1), c(1, 1))
    rhs <- c(upper, -lower, 0, 0, 1)
    vertices <- matrix(numeric(), 0, 2)
    for (pair in utils::combn(seq_len(nrow(constraints)), 2, simplify = FALSE)) {
      B <- constraints[pair, , drop = FALSE]
      if (abs(det(B)) < 1e-12) next
      v <- as.vector(solve(B, rhs[pair]))
      if (all(constraints %*% v <= rhs + 1e-12)) vertices <- rbind(vertices, v)
    }
    objective <- as.vector(10 * vertices %*% c(weights[1] + weights[3], weights[2] + weights[3]))
    expect_equal(z$lower, min(objective), tolerance = 1e-10)
    expect_equal(z$upper, max(objective), tolerance = 1e-10)
    expect_true(z$lower <= min(objective) + 1e-12)
    expect_true(z$upper >= max(objective) - 1e-12)
  }
})

test_that("conservative reachability projections do not claim attainable path support", {
  skip_if_not_installed("lpSolve")
  support <- build_structural_support(c(0, 1, 2), absorbing = c(1, 2),
    allowed_transitions = data.frame(from = c(0, 0), to = c(1, 2)),
    targets = list(first = 1, second = 2), tau = 10)
  p <- structure(list(bounds = data.frame(id = 1, target = c("first", "second"),
                                           lower = c(3, 3), upper = c(8, 8)),
                      center = matrix(0, 1, 2, dimnames = list(NULL, c("first", "second"))),
                      tau = 10, alpha = .1,
                      region = list(support = support, initial_state = 1L)),
                 class = "lifeconformal_prediction")
  z <- project_life(p, c(first = 1, second = 1))
  expect_equal(z$lower, 6, tolerance = 1e-12)
  expect_equal(z$upper, 10, tolerance = 1e-12)
  expect_true(attr(z, "conservative_outer_support"))
})

test_that("signed coefficient cancellation cannot erase a reachable summary", {
  targets <- list(a = c(0, 1, 2), b = 0, c = c(0, 1, 2))
  p <- projection_fixture(c(1, 0, 1), c(1, 1, 1), targets = targets, tau = 1)
  # Mathematically the identical a/c targets cancel, leaving exactly b.
  # A direct BLAS product can round both state coefficients to zero.
  z <- project_life(p, c(a = 1e16, b = 1, c = -1e16))
  expect_equal(z$lower, 0, tolerance = 0)
  expect_equal(z$upper, 1, tolerance = 0)
  expect_equal(attr(z, "diagnostics")$solver_calls, 0L)
  if (requireNamespace("lpSolve", quietly = TRUE)) {
    restricted <- projection_fixture(c(1, .25, 1), c(1, .75, 1),
                                      targets = targets, tau = 1)
    interval <- project_life(restricted, c(a = 1e16, b = 1, c = -1e16))
    expect_lte(interval$lower, .25)
    expect_gte(interval$upper, .75)
    expect_equal(interval$lower, .25, tolerance = 1e-12)
    expect_equal(interval$upper, .75, tolerance = 1e-12)
  }
})

test_that("coefficient envelopes retain exact residuals and unrepresentable sums", {
  # These references use exact binary powers and algebra, independent of R's
  # platform-dependent long-double implementation of sum().
  for (x in list(c(2^54, 1, -2^54), c(-2^54, 2^54, 1),
                c(1, 2^54, -2^54), c(2^54, -2^54, 1)))
    expect_equal(projection_sum_interval(x), c(1, 1, 1), tolerance = 0)
  positive <- projection_sum_interval(c(1, 2^-53))
  expect_equal(positive, c(1, 1, cdt_nextafter(1)), tolerance = 0)
  negative <- projection_sum_interval(c(-1, -2^-53))
  expect_equal(negative, c(-1, -cdt_nextafter(1), -1), tolerance = 0)
  subnormal <- .Machine$double.xmin * .Machine$double.eps
  expect_equal(projection_sum_interval(c(.Machine$double.xmin,
    -.Machine$double.xmin, subnormal)), rep(subnormal, 3L), tolerance = 0)
})

test_that("analytic and LP extrema preserve outward coefficient enclosures", {
  p <- projection_fixture(c(0, 0), c(1, 1),
    targets = list(a = 0, b = 0), tau = 1)
  z <- project_life(p, c(a = 1, b = 2^-53))
  expect_equal(z$lower, 0, tolerance = 0)
  expect_equal(z$upper, cdt_nextafter(1), tolerance = 0)
  if (requireNamespace("lpSolve", quietly = TRUE)) {
    singleton <- projection_fixture(c(1, 1), c(1, 1),
      targets = list(a = 0, b = 0), tau = 1)
    bounded <- project_life(singleton, c(a = 1, b = 2^-53))
    expect_lte(bounded$lower, 1)
    expect_gte(bounded$upper, cdt_nextafter(1))
  }
})

test_that("time-unit conversions widen outward and retain visibly empty boxes", {
  skip_if_not_installed("lpSolve")
  p <- projection_fixture(.2, .2, targets = list(alive = 0), tau = .3)
  z <- project_life(p, c(alive = 1))
  expect_lte(z$lower, .2)
  expect_gte(z$upper, .2)
  expect_equal(z$lower, .2, tolerance = 1e-12)
  expect_equal(z$upper, .2, tolerance = 1e-12)
  p$bounds$lower <- cdt_nextafter(.2)
  expect_warning(empty <- project_life(p, c(alive = 1)),
                   class = "lifeconformal_empty_projection")
  expect_identical(empty$status, "empty_region")
})
