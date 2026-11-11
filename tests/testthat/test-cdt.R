cdt_fixture <- function(id = 1, time = 8, status = 0, absorbing = 2,
                        initial_state = 0, targets = list(home = 0), tau = 10,
                        states = c(0, 1, 2), status_map = NULL,
                        allowed_transitions = NULL) {
  h <- normalize_histories(resolve_input(id = id, response = LifeHist(time, status)),
                            initial_state = initial_state, absorbing = absorbing,
                            states = states, status_map = status_map)
  p <- compile_problem(h, targets, tau, allowed_transitions = allowed_transitions)
  list(history = h, problem = p)
}

test_that("affine roots certify inside and outside during a sojourn", {
  f <- cdt_fixture()
  outside <- compute_cdt(f$history, f$problem, 2, 1, c(0, 1, 4, 7, 8, Inf))
  expect_equal(as.numeric(outside$time), c(2, 3, 6, 8, 0, 0))
  expect_identical(as.logical(outside$membership), c(FALSE, FALSE, FALSE, NA, TRUE, TRUE))
  inside <- compute_cdt(f$history, f$problem, 8, 1, c(0, 1, 2, 5, 8))
  expect_equal(as.numeric(inside$time), c(8, 8, 6, 3, 0))
  expect_identical(as.logical(inside$membership), c(FALSE, NA, TRUE, TRUE, TRUE))
  expect_equal(outside$diagnostics$numerical_later_certificates, 0)
  expect_true(outside$outer)
})

test_that("open extrema and observed transition endpoints use different histories", {
  # Staying home through t=3 excludes the exact home-time boundary 3 because
  # a future exit occurs strictly after 3. Observing death at 3 instead makes
  # the complete target equal 3, so the closed region includes that value.
  censored <- cdt_fixture(time = 3)
  dead <- cdt_fixture(time = 3, status = 2)
  a <- compute_cdt(censored$history, censored$problem, 2, 1, 1)
  b <- compute_cdt(dead$history, dead$problem, 2, 1, 1)
  expect_identical(a$membership[1, 1], FALSE)
  expect_identical(b$membership[1, 1], TRUE)
  expect_equal(a$time, matrix(3))
  expect_equal(b$time, matrix(3))
  # At the horizon the target is known; there is no open extreme left.
  horizon <- cdt_fixture(time = 10)
  c <- compute_cdt(horizon$history, horizon$problem, 9, 1, 1)
  expect_identical(c$membership[1, 1], TRUE)
  expect_equal(c$time[1, 1], 8)
  d <- compute_cdt(horizon$history, horizon$problem, 10, 1, 0)
  expect_identical(d$membership[1, 1], TRUE)
  expect_equal(d$time[1, 1], 10)
})

test_that("multiple targets retain conservative joint support certification", {
  f <- cdt_fixture(targets = list(home = 0, hospital = 1))
  a <- compute_cdt(f$history, f$problem, matrix(c(8, 8), 1), matrix(1, 1, 2), c(0, 1, 2, 8))
  # Both times cannot be 8 over a horizon of 10, but coordinate bounds alone
  # need not certify that at zero. This honest delay must be labelled outer.
  expect_true(a$outer)
  expect_gt(a$time[1, 1], 0)
  expect_identical(a$membership[1, 4], TRUE)
  expect_identical(a$membership[1, 1:3], c(FALSE, FALSE, FALSE))
})

test_that("initial and later absorption finish reward targets without future observation", {
  initial <- cdt_fixture(time = 1, initial_state = 2,
                         targets = list(dead = 2, alive = c(0, 1)))
  a <- compute_cdt(initial$history, initial$problem, matrix(c(10, 0), 1), matrix(1, 1, 2), c(0, 1))
  expect_equal(a$time, matrix(0, 1, 2))
  expect_true(all(a$membership))
  later <- cdt_fixture(time = 4, status = 2, targets = list(dead = 2))
  b <- compute_cdt(later$history, later$problem, 6, 1, c(0, 1))
  expect_true(all(b$resolved))
  expect_equal(b$time, matrix(4, 1, 2))
  expect_true(all(b$membership))
})

test_that("zero follow-up and same-time transition/censoring are respected", {
  f <- cdt_fixture(time = 0)
  a <- compute_cdt(f$history, f$problem, 5, 1, c(0, 5))
  expect_equal(a$time, matrix(0, 1, 2))
  expect_identical(a$membership[1, ], c(NA, TRUE))
  tied <- cdt_fixture(id = c(1, 1), time = c(3, 3), status = c(1, 0))
  b <- compute_cdt(tied$history, tied$problem, 2, 1, 1)
  expect_false(b$resolved[1, 1])
  expect_equal(b$time[1, 1], 3)
})

test_that("shared sweeps agree with a separate threshold-by-threshold reference", {
  f <- cdt_fixture(id = c(1, 2, 1, 3, 2, 1, 3),
                   time = c(2, 3, 5, 4, 9, 8, 8),
                   status = c(1, 1, 2, 1, 0, 0, 3),
                   status_map = c("1" = 1, "2" = 0, "3" = 2),
                   targets = list(home = 0, alive = c(0, 1)))
  mu <- matrix(c(2, 6, 5, 7, 3, 9), 3, 2, byrow = TRUE)
  sc <- matrix(c(1, 2, 2, 1, 1, 1), 3, 2, byrow = TRUE)
  q <- seq(0, 10, by = .125)
  a <- compute_cdt(f$history, f$problem, mu, sc, q)
  b <- compute_cdt_reference(f$history, f$problem, mu, sc, q)
  expect_equal(a$time, b$time)
  expect_identical(a$resolved, b$resolved)
  expect_identical(a$membership, b$membership)
  expect_lte(a$diagnostics$segments, nrow(f$history$intervals))
  expect_lte(a$diagnostics$analytic_roots, length(a$time))
  expect_identical(unserialize(serialize(a, NULL)), a)
})

test_that("future path changes cannot alter already observed certificates", {
  short <- cdt_fixture(time = 5)
  long <- cdt_fixture(id = c(1, 1), time = c(6, 9), status = c(1, 0))
  a <- compute_cdt(short$history, short$problem, 2, 1, c(0, 1, 2))
  b <- compute_cdt(long$history, long$problem, 2, 1, c(0, 1, 2))
  expect_equal(a$time, b$time)
  expect_identical(a$membership, b$membership)
})

test_that("CDT validates score alignment and sorted thresholds", {
  f <- cdt_fixture()
  expect_error(compute_cdt(f$history, f$problem, 2, 0, 1), "strictly positive")
  expect_error(compute_cdt(f$history, f$problem, 2, 1, c(2, 1)), "strictly increasing")
  expect_error(compute_cdt(f$history, f$problem, 2, 1, c(Inf, Inf)), "strictly increasing")
  expect_error(compute_cdt(f$history, f$problem, 2, 1, -1), "nonnegative")
  expect_error(compute_cdt(f$history, f$problem, matrix(2, dimnames = list(NULL, "wrong")), 1, 1),
               "target columns")
  expect_error(compute_cdt(f$history, f$problem, c(2, 3), 1, 1), "subject-by-target")
})

test_that("representable-root correction is deterministic and explicitly recorded", {
  next_one <- cdt_nextafter(1)
  expect_gt(next_one, 1)
  expect_equal(next_one - 1, .Machine$double.eps)
  expect_gt(cdt_nextafter(0), 0)
  expect_identical(cdt_nextafter(Inf), Inf)
  f <- cdt_fixture(time = 8.123)
  q <- seq(0, 4, length.out = 501)
  set.seed(17)
  before <- .Random.seed
  a <- compute_cdt(f$history, f$problem, 2.12345, 1.23456, q)
  expect_identical(.Random.seed, before)
  expect_gt(a$diagnostics$numerical_later_certificates, 0)
  expect_lte(a$diagnostics$max_root_delay, 1e-13)
  expect_equal(a$diagnostics$endpoint_fallbacks, 0)
  expect_equal(sum(a$numerical_later_certificate),
                a$diagnostics$numerical_later_certificates)
  # Outside certificates equal the analytic crossing up to the documented
  # representable arithmetic adjustment, not a later observation-time grid.
  outside <- which(a$membership[1, ] %in% FALSE)
  expect_equal(a$time[1, outside], 2.12345 + q[outside] * 1.23456,
                tolerance = 1e-13)
})

test_that("certificates agree with complete trajectory membership in generated cases", {
  set.seed(8104)
  for (run in seq_len(30)) {
    tim <- sort(sample(1:9, 3))
    f <- cdt_fixture(id = rep(1, 4), time = c(tim, 10), status = c(1, 2, 1, 0),
                     status_map = c("1" = 1, "2" = 0, "3" = 2),
                     targets = list(home = 0, alive = c(0, 1)))
    mu <- matrix(sample(0:10, 2), 1)
    sc <- matrix(sample(1:3, 2, replace = TRUE), 1)
    q <- seq(0, 10, .25)
    a <- compute_cdt(f$history, f$problem, mu, sc, q)
    b <- compute_cdt_reference(f$history, f$problem, mu, sc, q)
    truth <- c(tim[1] + tim[3] - tim[2], 10)
    expected <- max(abs((truth - mu[1, ]) / sc[1, ])) <= q
    expect_true(all(a$resolved))
    expect_identical(a$membership[1, ], expected)
    expect_equal(a$time, b$time, tolerance = 1e-12)
    expect_identical(a$membership, b$membership)
  }
})

test_that("compiled and R shared sweeps agree for distinct decimal score profiles", {
  set.seed(912)
  n <- 12L
  tim <- as.vector(vapply(seq_len(n), function(i) c(1, 3, 5, 8) + runif(1, 0, .5),
                           numeric(4)))
  f <- cdt_fixture(id = rep(seq_len(n), each = 4), time = tim,
                     status = rep(c(1, 2, 1, 0), n),
                     status_map = c("1" = 1, "2" = 0, "3" = 2),
                     targets = list(home = 0, alive = c(0, 1)))
  mu <- cbind(runif(n, 2, 8), runif(n, 4, 10))
  sc <- cbind(runif(n, 1, 3), runif(n, 1, 3))
  q <- seq(0, 6, length.out = 101)
  prepared <- prepare_cdt_sweep(f$history, f$problem, mu, sc)
  native <- sweep_cdt_thresholds(prepared, q)
  reference <- sweep_cdt_thresholds_R(prepared, q)
  expect_identical(native$time, reference$time)
  expect_identical(native$membership, reference$membership)
  expect_identical(native$resolved, reference$resolved)
  expect_identical(native$numerical_later_certificate,
                    reference$numerical_later_certificate)
  expect_identical(native$certificate_time, reference$certificate_time)
  expect_lte(native$diagnostics$segments, length(tim))
  expect_equal(native$diagnostics$endpoint_fallbacks, 0)
})

test_that("several absorbing states and structural restrictions retain labels", {
  f <- cdt_fixture(id = c("a", "b"), time = c(3, 7), status = c(2, 3),
                     absorbing = c(2, 3), states = c(0, 1, 2, 3),
                     targets = list(terminal = c(2, 3)))
  a <- compute_cdt(f$history, f$problem, matrix(c(7, 3), 2), matrix(1, 2), c(0, 1))
  expect_true(all(a$membership))
  expect_equal(a$time[, 1], c(3, 7))
  expect_identical(a$subject_id, c("a", "b"))
  constrained <- cdt_fixture(time = 8, absorbing = 2,
                              allowed_transitions = data.frame(from = 0, to = 2),
                              targets = list(hospital = 1))
  b <- compute_cdt(constrained$history, constrained$problem, 0, 1, 0)
  expect_equal(b$time[1, 1], 0)
  expect_true(b$membership[1, 1])
})

test_that("inclusive stationary score boundaries preserve stopping times on prefixes", {
  # The observed home-only path has feasible home times [t,1]. Its
  # complementary hospital target has feasible times [0,1-t]. At this
  # q, reconstructing the fixed endpoint via q*scale loses one ulp.
  center <- 0x1.e5580b1c51ac3p-1
  scale <- 0x1.49acae903a6cdp-3
  q <- (1 - center) / scale
  expect_lt(q * scale, 1 - center)
  ends <- c(.8, .925, .95, 1)
  grid <- c(-cdt_nextafter(-q), q, cdt_nextafter(q))
  for (home in c(TRUE, FALSE)) {
    f <- cdt_fixture(id = seq_along(ends), time = ends, status = rep(0, length(ends)),
      targets = if (home) list(home = 0) else list(hospital = 1), tau = 1)
    mu <- if (home) center else 1 - center
    inside <- function(t) {
      limits <- if (home) c(t, 1) else c(0, 1 - t)
      max(abs((limits - mu) / scale)) <= q
    }
    # Independent bisection of the actual standardized score on this path.
    lo <- 0
    hi <- 1
    repeat {
      mid <- lo + (hi - lo) / 2
      if (mid <= lo || mid >= hi) break
      if (inside(mid)) hi <- mid else lo <- mid
    }
    expect_false(inside(lo))
    expect_true(inside(hi))
    expect_lt(hi, ends[2])
    prepared <- prepare_cdt_sweep(f$history, f$problem, rep(mu, 4), rep(scale, 4))
    native <- sweep_cdt_thresholds(prepared, grid)
    r <- sweep_cdt_thresholds_R(prepared, grid)
    ref <- compute_cdt_reference(f$history, f$problem, rep(mu, 4), rep(scale, 4), grid)
    expect_identical(native$time, r$time)
    expect_identical(native$membership, r$membership)
    expect_identical(native$membership[, 2], c(NA, TRUE, TRUE, TRUE))
    expect_identical(ref$membership[, 2], native$membership[, 2])
    expect_equal(native$time[, 2], pmin(ends, hi), tolerance = 1e-15)
    expect_equal(ref$time[, 2], pmin(ends, hi), tolerance = 1e-15)
    expect_identical(native$time[2:4, 2], rep(native$time[2, 2], 3))
    expect_false(any(native$membership[, 1] %in% TRUE))
    expect_true(all(native$membership[2:4, 3]))
    expect_equal(native$diagnostics$endpoint_fallbacks, 0)
  }
})
