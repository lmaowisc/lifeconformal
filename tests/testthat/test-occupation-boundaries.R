occupation_boundary_fixture <- function(kind = c("absorbing_alive", "recurrent", "terminal_death")) {
  kind <- match.arg(kind)
  tau <- 1 + .Machine$double.eps
  time <- .Machine$double.eps / 2
  if (kind == "absorbing_alive") {
    history <- normalize_histories(resolve_input(id = 1L,
      response = LifeHist(time, 1L)), initial_state = 0L,
      absorbing = c(1L, 2L), states = 0:2)
    allowed <- data.frame(from = c(0L, 0L), to = c(1L, 2L))
  } else {
    times <- if (kind == "recurrent") c(time, tau, tau) else c(time, tau)
    status <- if (kind == "recurrent") c(1L, 3L, 0L) else c(1L, 2L)
    history <- normalize_histories(resolve_input(id = rep(1L, length(times)),
      response = LifeHist(times, status)), initial_state = 0L,
      absorbing = 2L, states = 0:2, status_map = c("1" = 1L, "2" = 2L, "3" = 0L))
    allowed <- NULL
  }
  problem <- compile_problem(history, list(alive = c(0L, 1L)), tau,
    allowed_transitions = allowed)
  list(history = history, problem = problem, tau = tau, time = time)
}

occupation_boundary_sweeps <- function(fixture, center, q) {
  prepared <- prepare_cdt_sweep(fixture$history, fixture$problem,
    matrix(center), matrix(1))
  list(native = sweep_cdt_thresholds(prepared, q),
    R = sweep_cdt_thresholds_R(prepared, q),
    separate_reference = compute_cdt_reference(fixture$history, fixture$problem,
      matrix(center), matrix(1), q))
}

test_that("occupation completion preserves complete binary reward totals", {
  tau <- 1 + .Machine$double.eps
  time <- .Machine$double.eps / 2
  reward <- matrix(rep(c(1L, 0L), 3L), 3L, 2L, byrow = TRUE)
  accrued <- cbind(c(time, 1, 0), c(0, 0, 0))
  value <- occupation_complete(rep(tau, 3L), c(time, 1, 0), reward, accrued)
  expect_identical(dim(value), c(3L, 2L))
  expect_identical(as.numeric(value), c(rep(tau, 3L), rep(0, 3L)))
  # Integer inputs use the same public internal wrapper and retain dimensions.
  integer_value <- occupation_complete(c(2L, 3L), c(1L, 2L),
    matrix(1L, 2L, 1L), matrix(c(1L, 2L), 2L, 1L))
  expect_identical(as.numeric(integer_value), c(2, 3))
  expect_identical(dim(integer_value), c(2L, 1L))
  expect_error(occupation_complete(tau, c(time, time), matrix(1, 2L), matrix(time, 2L)))
  expect_error(occupation_complete(c(tau, tau), c(time, time), matrix(1, 2L), matrix(time, 1L, 2L)))
})

test_that("native occupation completion enforces aligned double-array inputs", {
  native <- function(end, start, reward, accrued)
    .Call("life_occupation_complete", end, start, reward, accrued,
      PACKAGE = "lifeconformal")
  valid <- list(end = 2, start = 1, reward = matrix(1), accrued = matrix(1))
  for (field in names(valid)) {
    bad <- valid
    if (is.matrix(bad[[field]])) storage.mode(bad[[field]]) <- "integer" else
      bad[[field]] <- as.integer(bad[[field]])
    expect_error(do.call(native, bad), "double arrays")
  }
  expect_error(native(2, 1, matrix(c(1, 0), 2L, 1L),
    matrix(c(1, 0), 1L, 2L)), "matching dimensions")
  expect_error(native(c(2, 3), 1, matrix(1), matrix(1)), "one entry per matrix row")
  expect_error(native(2, matrix(1), matrix(1), matrix(1)), "one entry per matrix row")
  expect_error(native(2, 1, 1, matrix(1)), "matching dimensions")
  for (times in list(c(Inf, 1), c(2, Inf), c(NA_real_, 1), c(2, NA_real_),
                     c(2, -1), c(1, 2))) {
    expect_error(native(times[1L], times[2L], matrix(1), matrix(1)),
      "finite with 0 <= start <= end")
  }
  for (reward in c(NA_real_, Inf, -.5, .5, 2))
    expect_error(native(2, 1, matrix(reward), matrix(1)), "rewards must be binary")
  for (accrued in c(NA_real_, Inf, -.5))
    expect_error(native(2, 1, matrix(1), matrix(accrued)), "accrued values finite and nonnegative")

  reward <- matrix(c(1, 0, 0, 1), 2L,
    dimnames = list(c("first", "second"), c("one", "two")))
  accrued <- matrix(c(1, 0, 0, 1), 2L,
    dimnames = list(c("A", "B"), c("alive", "hospital")))
  arguments <- list(end = c(2, 3), start = c(1, 1), reward = reward, accrued = accrued)
  before <- serialize(arguments, NULL)
  value <- do.call(native, arguments)
  expect_identical(dimnames(value), dimnames(accrued))
  expect_identical(serialize(arguments, NULL), before)
  empty_rows <- matrix(numeric(), 0L, 2L, dimnames = list(NULL, c("alive", "home")))
  expect_identical(native(numeric(), numeric(), empty_rows, empty_rows), empty_rows)
  empty_targets <- matrix(numeric(), 2L, 0L, dimnames = list(c("A", "B"), NULL))
  expect_identical(native(c(1, 2), c(0, 0), empty_targets, empty_targets), empty_targets)
})

test_that("zero duration preserves accrued occupation despite a changed reward", {
  tau <- 1 + .Machine$double.eps
  accrued <- rbind(c(tau, 0), c(0, tau), c(0, 0))
  reward <- rbind(c(0L, 1L), c(1L, 0L), c(1L, 0L))
  value <- occupation_complete(c(tau, tau, 0), c(tau, tau, 0), reward, accrued)
  expect_identical(as.numeric(value), as.numeric(accrued))
  expect_identical(dim(value), dim(accrued))
  f <- occupation_boundary_fixture("terminal_death")
  prepared <- prepare_occupation_history(f$history, f$problem)
  expect_identical(as.numeric(prepared$prefix[, 1L]), c(f$time, f$tau))
  snapshot <- evaluate_occupation_history(prepared, c(1L, 1L),
    c(f$tau, f$tau), side = c("left", "right"))
  expect_identical(as.numeric(snapshot$accrued[, 1L]), c(f$tau, f$tau))
  expect_identical(snapshot$state, c(2L, 3L))
  for (value in occupation_boundary_sweeps(f, f$tau, 0)) {
    expect_identical(as.logical(value$membership), TRUE)
    expect_identical(as.numeric(value$time), f$tau)
    expect_identical(as.logical(value$resolved), TRUE)
  }
  # A zero-duration censoring row may use the newly entered reward class.
  # It must not invalidate the all-alive prefix established through the jump.
  history <- normalize_histories(resolve_input(id = rep(1L, 3L),
    response = LifeHist(c(f$time, f$tau, f$tau), c(1L, 3L, 0L))),
    initial_state = 0L, absorbing = 2L, states = 0:3)
  problem <- compile_problem(history, list(alive = c(0L, 1L)), f$tau)
  tied <- prepare_occupation_history(history, problem)
  expect_identical(as.numeric(tied$prefix[, 1L]), c(f$time, f$tau, f$tau))
  expect_identical(as.numeric(tied$rewards[tied$intervals$from, 1L]), c(1, 1, 0))
  expect_identical(as.numeric(evaluate_occupation_history(tied, 1L, f$tau)$accrued), f$tau)
})

test_that("reward-equivalent transitions preserve exact all-alive history prefixes", {
  f <- occupation_boundary_fixture("recurrent")
  prepared <- prepare_occupation_history(f$history, f$problem)
  expect_identical(as.numeric(prepared$prefix[, 1L]), c(f$time, f$tau, f$tau))
  snapshot <- evaluate_occupation_history(prepared, rep(1L, 4L),
    c(f$time, 1, f$tau, f$tau), side = c("right", "right", "left", "right"))
  expect_identical(as.numeric(snapshot$accrued[, 1L]), c(f$time, 1, f$tau, f$tau))
  for (value in occupation_boundary_sweeps(f, f$tau, 0)) {
    expect_identical(as.logical(value$membership), TRUE)
    expect_identical(as.numeric(value$time), f$tau)
    expect_identical(as.logical(value$resolved), TRUE)
  }
  q <- c(1, f$tau, 1 + 2 * .Machine$double.eps)
  for (value in occupation_boundary_sweeps(f, 0, q)) {
    expect_identical(as.logical(value$membership), c(FALSE, TRUE, TRUE))
    expect_identical(as.numeric(value$time), c(1, 0, 0))
  }
})

test_that("absorbing alive paths retain inclusive atoms and adjacent exclusions", {
  f <- occupation_boundary_fixture("absorbing_alive")
  prepared <- prepare_occupation_history(f$history, f$problem)
  snapshot <- evaluate_occupation_history(prepared, c(1L, 1L), c(f$time, f$tau))
  expect_identical(as.numeric(snapshot$accrued[, 1L]), c(f$time, f$tau))
  box <- classify_coverage_box(f$time, 2L, f$time, f$problem,
    center = f$tau, scale = 1, q = 0)
  expect_identical(as.numeric(box$lower), f$tau)
  expect_identical(as.numeric(box$upper), f$tau)
  expect_identical(box$membership, TRUE)
  for (value in occupation_boundary_sweeps(f, f$tau, 0)) {
    expect_identical(as.logical(value$membership), TRUE)
    expect_identical(as.numeric(value$time), f$time)
    if (!is.null(value$certificate_time))
      expect_identical(as.numeric(value$certificate_time), f$time)
  }
  q <- c(1, f$tau, 1 + 2 * .Machine$double.eps)
  for (value in occupation_boundary_sweeps(f, 0, q)) {
    expect_identical(as.logical(value$membership), c(FALSE, TRUE, TRUE))
    expect_identical(as.numeric(value$time), c(f$time, 0, 0))
    expect_identical(as.logical(value$resolved), rep(TRUE, 3L))
  }
})

test_that("within-sojourn root verification recomputes exact occupation at the horizon", {
  f <- occupation_boundary_fixture("recurrent")
  initial <- cdt_box_bounds(f$time, 2L, f$time, f$problem,
    center = f$tau, scale = 1)
  segment <- list(start = f$time, end = f$tau, lower = initial$lower,
    upper = initial$upper, reward = 1, minimum = 0, maximum = 1,
    center = f$tau, scale = 1)
  roots <- c(1, f$tau)
  for (event_at_end in c(FALSE, TRUE)) {
    expect_identical(cdt_verify_roots(roots, segment, accrued = f$time,
      q = c(0, 0), membership = TRUE, tau = f$tau,
      event_at_end = event_at_end), c(FALSE, TRUE))
  }
  # The preceding representable time is still a feasible interval, and the
  # exact horizon is an included atom. A broad numeric tolerance would hide it.
  before <- classify_coverage_box(1, 2L, 1, f$problem, f$tau, 1, 0)
  at <- classify_coverage_box(f$tau, 2L, f$tau, f$problem, f$tau, 1, 0)
  expect_false(before$resolved)
  expect_identical(before$membership, NA)
  expect_true(at$resolved)
  expect_identical(at$membership, TRUE)
})
