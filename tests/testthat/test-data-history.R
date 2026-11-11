make_event_input <- function(id, time, status, x = NULL) {
  resolve_input(id = id, response = LifeHist(time, status), x = x)
}

make_transition_input <- function(id, start, stop, from, to, event, x = NULL) {
  resolve_input(id = id,
                response = LifeHist(start = start, stop = stop, from = from,
                                    to = to, event = event), x = x)
}

test_that("recurrent event and explicit inputs describe the same canonical paths", {
  ids <- c("a", "a", "a", "b", "c", "c")
  x <- data.frame(age = c(60, 60, 60, 55, 70, 70))
  e <- make_event_input(ids, c(20, 35, 80, 50, 15, 70), c(1, 2, 0, 0, 1, 3), x)
  t <- make_transition_input(ids, c(0, 20, 35, 0, 0, 15),
                             c(20, 35, 80, 50, 15, 70), c(0, 1, 0, 0, 0, 1),
                             c(1, 0, NA, NA, 1, 2), c(1, 1, 0, 0, 1, 1), x)
  a <- normalize_histories(e, status_map = c("1" = 1, "2" = 0, "3" = 2),
                            absorbing = 2)
  b <- normalize_histories(t, absorbing = 2)
  expect_identical(a$intervals, b$intervals)
  expect_identical(a$subjects, b$subjects)
  expect_identical(a$baseline, b$baseline)
  expect_identical(a$row_index, list(1:3, 4L, 5:6))
  expect_identical(a$subjects$observation_end, c(80, 50, 70))
  expect_identical(a$subjects$censoring, c(TRUE, TRUE, FALSE))
  expect_identical(a$diagnostics$n_transitions, 4L)
  expect_identical(a$diagnostics$rows_omitted, 0L)
})

test_that("subject indexing preserves interleaved records and typed baselines", {
  e <- make_event_input(c("b", "a", "b", "a"), c(1, 2, 3, 4), c(1, 1, 0, 0),
                        data.frame(group = factor(c("B", "A", "B", "A"))))
  a <- normalize_histories(e, absorbing = integer(0))
  expect_identical(a$subjects$id, c("b", "a"))
  expect_identical(a$intervals$original_row, 1:4)
  expect_identical(a$row_index, list(c(1L, 3L), c(2L, 4L)))
  expect_identical(a$intervals$start, c(0, 0, 1, 2))
  expect_identical(as.character(a$baseline$group), c("B", "A"))
  expect_s3_class(a$baseline$group, "factor")
})

test_that("state and event coding remain separate with several absorbing states", {
  e <- make_event_input(c(1, 1, 1, 2), c(2, 3, 9, 5), c(1, 2, 3, 4))
  a <- normalize_histories(e, initial_state = "home",
                            status_map = c("1" = "hospital", "2" = "home",
                                           "3" = "death", "4" = "discharge"),
                            absorbing = c("death", "discharge"))
  labels <- a$schema$state_labels
  expect_identical(labels[a$subjects$absorption_state], c("death", "discharge"))
  expect_identical(labels[a$intervals$from], c("home", "hospital", "home", "home"))
  expect_identical(a$subjects$absorption_time, c(9, 5))
  expect_error(normalize_histories(e, status_map = c("0" = 0, "1" = 1), absorbing = 3),
               "nonzero")
  expect_error(normalize_histories(e, status_map = c("1" = 1), absorbing = 3),
               "Unknown event code.*subject '1'.*row 2")
})

test_that("initial states can be inferred, common, or repeated by subject", {
  d <- data.frame(id = c(1, 1, 2), time = c(2, 4, 3), status = c(2, 0, 0),
                  initial = c(1, 1, 0), age = c(40, 40, 50))
  raw <- resolve_input(LifeHist(time, status) ~ age, data = d)
  a <- normalize_histories(raw, initial_state = ~ initial, absorbing = integer(0))
  expect_equal(a$schema$state_labels[a$subjects$initial_state], c(1, 0))
  expect_identical(a$schema$initial_state_rule$type, "column")
  expect_error(normalize_histories(raw, initial_state = c(1, 0, 0), absorbing = 2),
               "Initial state changes.*row 2")
  t <- make_transition_input(1, 0, 4, 2, NA_real_, 0)
  b <- normalize_histories(t, absorbing = integer(0))
  expect_equal(b$schema$state_labels[b$subjects$initial_state], 2)
  expect_error(normalize_histories(t, initial_state = 0, absorbing = integer(0)),
               "Origin state")
})

test_that("the only admissible tied records are transition then terminal censoring", {
  e <- make_event_input(c(1, 1), c(5, 5), c(1, 0))
  a <- normalize_histories(e, absorbing = integer(0))
  t <- make_transition_input(c(1, 1), c(0, 5), c(5, 5), c(0, 1), c(1, NA), c(1, 0))
  b <- normalize_histories(t, absorbing = integer(0))
  expect_identical(a$intervals, b$intervals)
  expect_identical(a$subjects$censoring, TRUE)
  expect_error(normalize_histories(make_event_input(c(1, 1, 1), c(5, 5, 6), c(1, 2, 0)),
                                   absorbing = integer(0)), "positive sojourn")
  expect_error(normalize_histories(make_event_input(c(1, 1), c(5, 5), c(0, 1)),
                                   absorbing = integer(0)), "after censoring")
})

test_that("origin censoring and origin absorption retain distinct meanings", {
  a <- normalize_histories(make_event_input(1, 0, 0), absorbing = integer(0))
  b <- normalize_histories(make_transition_input(1, 0, 0, 0, NA_real_, 0),
                            absorbing = integer(0))
  expect_identical(a$intervals, b$intervals)
  expect_identical(a$subjects$censoring, TRUE)
  ordinary_na <- normalize_histories(make_transition_input(1, 0, 0, 0, NA, 0),
                                     absorbing = integer(0))
  expect_identical(ordinary_na$intervals, b$intervals)
  expect_identical(ordinary_na$schema$state_labels, 0)
  expect_identical(a$diagnostics$n_origin_censored, 1L)
  expect_error(normalize_histories(make_event_input(c(1, 1), c(0, 2), c(1, 0)),
                                   absorbing = integer(0)), "after zero")
  c <- normalize_histories(make_event_input(1, 100, 0), initial_state = 2,
                            absorbing = 2)
  expect_identical(c$subjects$observation_end, 100)
  expect_identical(c$subjects$path_completion, 0)
  expect_identical(c$subjects$censoring, FALSE)
  expect_identical(c$intervals$terminal_censor, FALSE)
})

test_that("invalid time order and explicit gaps fail with identifiable records", {
  expect_error(normalize_histories(make_event_input(c("a", "a"), c(5, 4), c(1, 0)),
                                   absorbing = 2), "chronological.*subject 'a'.*row 2")
  expect_error(normalize_histories(make_transition_input(c(1, 1), c(0, 6), c(5, 9),
                                                         c(0, 1), c(1, NA), c(1, 0)),
                                   absorbing = 2), "contiguous.*row 2")
  expect_error(normalize_histories(make_transition_input(1, 1, 5, 0, NA_real_, 0),
                                   absorbing = 2), "start at zero")
  expect_error(normalize_histories(make_event_input(1, Inf, 0), absorbing = 2),
               "finite, nonnegative")
  expect_error(normalize_histories(make_event_input(1, -1, 0), absorbing = 2),
               "finite, nonnegative")
})

test_that("terminal and absorption rules are explicit", {
  expect_error(normalize_histories(make_event_input(1, 2, 1), absorbing = 2),
               "final censoring")
  expect_error(normalize_histories(make_event_input(c(1, 1), c(2, 3), c(2, 0)),
                                   absorbing = 2), "after censoring or absorbing")
  expect_error(normalize_histories(make_event_input(1, 2, 1), initial_state = 2,
                                   absorbing = 2), "leaves a declared absorbing")
  expect_error(normalize_histories(make_event_input(c(1, 1), c(2, 3), c(1, 0)),
                                   initial_state = 1, absorbing = integer(0)),
               "different state")
  expect_error(normalize_histories(make_event_input(1, 2, 0)), "Declare absorbing")
})

test_that("a fixed state dictionary is enforced without learning from new data", {
  a <- normalize_histories(make_event_input(1, 3, 2), absorbing = 2,
                            states = c(0, 1, 2, 3))
  expect_equal(a$schema$state_labels, c(0, 1, 2, 3))
  expect_error(normalize_histories(make_event_input(1, 3, 4), absorbing = 2,
                                   states = c(0, 1, 2)), "fixed dictionary: 4")
  expect_error(normalize_histories(make_event_input(1, 3, 1 + 1e-15),
                                   initial_state = 1, absorbing = 1 + 1e-15),
               "indistinguishable text")
})

test_that("factor mappings preserve their event-code names", {
  map <- factor(c("hospital", "death"))
  names(map) <- c("7", "9")
  a <- normalize_histories(make_event_input(c(1, 1), c(2, 5), c(7, 9)),
                            initial_state = "home", status_map = map,
                            absorbing = "death")
  expect_identical(a$schema$state_labels[a$intervals$to], c("hospital", "death"))
})

test_that("baseline changes fail while factor types survive normalization", {
  raw <- make_event_input(c("a", "a"), c(2, 3), c(1, 0),
                          data.frame(age = c(30, 31)))
  expect_error(normalize_histories(raw, absorbing = 2), "age.*subject 'a'.*row 2")
})

test_that("canonical histories survive serialization", {
  a <- normalize_histories(make_event_input(c(1, 1, 2), c(2, 5, 9), c(1, 0, 2)),
                            absorbing = 2)
  expect_identical(unserialize(serialize(a, NULL)), a)
})
