problem_event_history <- function(id, time, status, initial_state = 0,
                                   absorbing = 2, states = c(0, 1, 2),
                                   status_map = c("1" = 1, "2" = 0, "3" = 2)) {
  normalize_histories(
    resolve_input(id = id, response = LifeHist(time, status)),
    initial_state = initial_state, absorbing = absorbing,
    states = states, status_map = status_map
  )
}

test_that("targets use actual labels and preserve joint occupation structure", {
  targets <- compile_targets(
    list(home = c(0, 0), alive = c(7, 0), total = c(7, 0, 3)),
    states = c(7, 0, 3), tau = 10
  )
  expected <- rbind(home = c(0, 1, 0), alive = c(1, 1, 0), total = c(1, 1, 1))
  colnames(expected) <- c("7", "0", "3")
  expect_identical(targets$A, expected)
  expect_identical(targets$rewards, t(expected))
  expect_identical(targets$sets$home, 0)
  expect_equal(targets$support$lower, c(home = 0, alive = 0, total = 10))
  expect_equal(targets$support$upper, c(home = 10, alive = 10, total = 10))
  # A single common state-duration vector retains home <= alive and total=10.
  v <- c(2, 3, 5)
  z <- drop(targets$support$A %*% v)
  expect_equal(z, c(home = 3, alive = 5, total = 10))
  expect_identical(targets$support$state_time_total, 10)

  single <- compile_targets(list(stay = factor("home")), "home", 4)
  expect_identical(dim(single$A), c(1L, 1L))
  expect_equal(single$rewards[1, 1], 1)
  expect_equal(single$support$lower, c(stay = 4))
  expect_equal(single$support$upper, c(stay = 4))
})

test_that("invalid target definitions and horizons cannot silently change a problem", {
  expect_error(compile_targets(list(0), c(0, 1), 5), "unique, nonempty names")
  expect_error(compile_targets(list(a = 0, a = 1), c(0, 1), 5), "unique")
  expect_error(compile_targets(setNames(list(0), " "), c(0, 1), 5), "nonempty names")
  expect_error(compile_targets(list(home = numeric()), c(0, 1), 5), "nonempty")
  expect_error(compile_targets(list(home = 0), c(1, 2), 5), "unknown state label.*0")
  expect_error(compile_targets(list(home = 1), c(7, 0, 3), 5), "unknown state label.*1")
  expect_error(compile_targets(list(home = NA_real_), c(0, 1), 5), "missing")
  expect_error(compile_targets(list(home = TRUE), c(0, 1), 5), "state labels")
  expect_error(compile_targets(list(home = 0), c(0, 0), 5), "unique labels")
  for (tau in list(0, -1, Inf, NA_real_, numeric(), c(1, 2), "5")) {
    expect_error(compile_targets(list(home = 0), c(0, 1), tau), "tau")
  }
})

test_that("the working graph counts recurrent events and observed origin risk once", {
  h <- problem_event_history(
    c("a", "b", "a", "b", "a", "c"), c(2, 1, 5, 4, 9, 7),
    c(1, 1, 2, 3, 0, 0), initial_state = c(0, 0, 0, 0, 0, 2),
    states = c(2, 0, 1, 4)
  )
  graph <- build_working_graph(h)
  expect_identical(graph$edges$from, c(2L, 3L, 3L))
  expect_identical(graph$edges$to, c(3L, 2L, 1L))
  expect_identical(graph$edges$n_events, c(2L, 1L, 1L))
  expect_equal(graph$origin_risk_time, c("2" = 0, "0" = 7, "1" = 6, "4" = 0))
  expect_equal(graph$edges$origin_risk_time, c(7, 6, 6))
  expect_identical(graph$edges$from_label, c(0, 1, 1))
  expect_identical(graph$edges$to_label, c(1, 0, 2))
  expect_identical(graph$absorbing, 1L)
  expect_false(graph$adjacency["4", "0"])
  expect_equal(sum(graph$adjacency), 3)
  expect_identical(h$subjects$observation_end, c(9, 4, 7))
  expect_error(build_working_graph(h, c(0, 1, 2, 4)), "dictionaries differ")
})

test_that("no-event working graphs stay empty without declaring absorption", {
  h <- problem_event_history(c("a", "b"), c(0, 5), c(0, 0), absorbing = integer(0))
  graph <- build_working_graph(h)
  expect_equal(nrow(graph$edges), 0)
  expect_identical(graph$edges$from, integer())
  expect_equal(sum(graph$adjacency), 0)
  expect_equal(graph$origin_risk_time, c("0" = 5, "1" = 0, "2" = 0))
  expect_length(graph$absorbing, 0)
  support <- build_structural_support(c(0, 1, 2), integer(0),
                                      targets = list(home = 0), tau = 10)
  expect_true(all(support$reachability))
})

test_that("default support allows unobserved returns and respects multiple absorbing states", {
  support <- build_structural_support(
    c(0, 1, 2, 3), absorbing = c(2, 3),
    targets = list(home = 0, alive = c(0, 1), dead = c(2, 3)), tau = 10
  )
  expect_false(support$specified)
  expect_false(any(diag(support$adjacency)))
  expect_true(support$adjacency["1", "0"])
  expect_false(any(support$adjacency[c("2", "3"), ]))
  expect_true(all(support$reachability[c("0", "1"), ]))
  expect_equal(support$reachable[["2"]], 3L)
  expect_equal(support$reachable[["3"]], 4L)
  expect_equal(unname(support$reward_min["1", ]), c(0, 0, 0))
  expect_equal(unname(support$reward_max["1", ]), c(1, 1, 1))
  expect_equal(unname(support$lower["2", ]), c(0, 0, 10))
  expect_identical(support$lower["2", ], support$upper["2", ])
  expect_true(support$outer)
  expect_identical(support$bound_attainment, "not_assumed")
})

test_that("declared reachability handles cycles, isolated states and empty graphs", {
  edges <- data.frame(from = c(0, 1, 1), to = c(1, 0, 2))
  support <- build_structural_support(
    c(2, 0, 1, 4), 2, edges, list(home = 0, alive = c(0, 1)), 10
  )
  expect_true(support$specified)
  expect_true(all(support$reachability["0", c("2", "0", "1")]))
  expect_false(support$reachability["0", "4"])
  expect_equal(support$reachable[["4"]], 4L)
  # A structurally isolated state need not be declared biologically absorbing.
  expect_identical(support$absorbing, 1L)
  expect_equal(unname(support$lower["4", ]), c(0, 0))
  expect_equal(unname(support$upper["4", ]), c(0, 0))

  empty <- build_structural_support(
    c("a", "b"), character(), data.frame(from = character(), to = character()),
    list(a = "a", b = "b"), 3
  )
  expect_equal(unname(empty$reachability), diag(TRUE, 2))
  expect_equal(unname(empty$lower), diag(3, 2))
  expect_identical(empty$lower, empty$upper)

  matrix_graph <- support$adjacency[c("4", "1", "2", "0"),
                                     c("1", "0", "4", "2")]
  from_matrix <- build_structural_support(
    c(2, 0, 1, 4), 2, matrix_graph, support$targets, 10
  )
  expect_identical(from_matrix$adjacency, support$adjacency)
  expect_identical(from_matrix$reachability, support$reachability)
  from_edges <- build_structural_support(c(2, 0, 1, 4), 2,
                                         unname(as.matrix(edges)), support$targets, 10)
  expect_identical(from_edges$adjacency, support$adjacency)
})

test_that("structural graph validation prevents label and absorption mistakes", {
  make_support <- function(graph, absorbing = 2) {
    build_structural_support(c(0, 1, 2), absorbing, graph, list(home = 0), 10)
  }
  expect_error(make_support(data.frame(from = 0, to = 7)), "unknown state")
  expect_error(make_support(data.frame(from = 0, to = 0)), "self-transitions")
  expect_error(make_support(data.frame(from = 2, to = 1)), "absorbing")
  expect_error(make_support(data.frame(from = c(0, 0), to = c(1, 1))), "duplicate")
  expect_error(make_support(data.frame(a = 0, b = 1)), "columns from and to")
  expect_error(make_support(matrix(FALSE, 2, 2)), "K by K")
  expect_error(make_support(diag(TRUE, 3)), "self-transitions")
  expect_error(make_support(NULL, c(2, 2)), "duplicate")
  expect_error(make_support(NULL, 7), "absent from states")
  expect_error(make_support(NULL, NULL), "Declare absorbing")
  named <- matrix(FALSE, 3, 3, dimnames = list(c("0", "1", "1"), c("0", "1", "2")))
  expect_error(make_support(named), "row and column names")
  compiled <- compile_targets(list(home = 0), c(0, 1, 2), 10)
  expect_error(build_structural_support(c(0, 1, 2), 2, targets = compiled, tau = 9),
               "different state dictionary or horizon")
})

test_that("a frozen problem separates fitted edges from structural knowledge", {
  train <- problem_event_history(c("train", "train"), c(2, 5), c(1, 0))
  problem <- compile_problem(train, list(home = 0, alive = c(0, 1)), 10)
  expect_identical(problem$states, c(0, 1, 2))
  expect_identical(problem$state_index, c("0" = 1L, "1" = 2L, "2" = 3L))
  expect_identical(problem$A, problem$targets$A)
  expect_identical(problem$rewards, problem$targets$rewards)
  expect_false(problem$working_graph$adjacency["0", "2"])
  expect_true(problem$structural_support$adjacency["0", "2"])
  cal <- problem_event_history("cal", 3, 3)
  expect_s3_class(prepare_occupation_history(cal, problem), "life_occupation_history")
  expect_identical(unserialize(serialize(problem, NULL)), problem)

  restricted <- compile_problem(train, list(home = 0), 10,
                                allowed_transitions = data.frame(from = 0, to = 1))
  expect_error(prepare_occupation_history(cal, restricted),
               "structural support.*subject 'cal'.*row 1")
  expect_error(compile_problem(train, list(home = 0), 10,
                               allowed_transitions = data.frame(from = 0, to = 2)),
               "structural support.*subject 'train'.*row 1")
  expect_error(compile_problem(train, list(home = 0), 10, score = NA_character_),
               "score-family")
})

test_that("accrued occupations match a recurrent path and absorbing continuation", {
  h <- problem_event_history(
    c("a", "b", "a", "c", "a"), c(2, 6, 5, 12, 8), c(1, 0, 2, 0, 3),
    initial_state = c(0, 0, 0, 2, 0)
  )
  problem <- compile_problem(h, list(home = 0, hospital = 1, alive = c(0, 1),
                                    dead = 2, total = c(0, 1, 2)), 10)
  prepared <- prepare_occupation_history(h, problem)
  times <- c(0, 2, 4, 5, 7, 8, 10, 20)
  snapshots <- evaluate_occupation_history(prepared, rep(1L, length(times)), times)
  expected <- cbind(home = c(0, 2, 2, 2, 4, 5, 5, 5),
                     hospital = c(0, 0, 2, 3, 3, 3, 3, 3),
                     alive = c(0, 2, 4, 5, 7, 8, 8, 8),
                     dead = c(0, 0, 0, 0, 0, 0, 2, 2),
                     total = c(0, 2, 4, 5, 7, 8, 10, 10))
  expect_identical(snapshots$accrued, expected)
  expect_identical(problem$states[snapshots$state], c(0, 1, 1, 0, 0, 2, 2, 2))
  expect_identical(snapshots$requested_time, times)
  expect_identical(snapshots$time, pmin(times, 10))
  expect_identical(prepared$intervals$original_row, seq_len(5))
  initial_absorbed <- evaluate_occupation_history(prepared, rep(3L, 3), c(0, 7, 20))
  expect_equal(initial_absorbed$accrued[, "dead"], c(0, 7, 10))
  expect_equal(initial_absorbed$accrued[, "alive"], c(0, 0, 0))
  expect_true(all(initial_absorbed$absorbed))
  expect_equal(problem$working_graph$origin_risk_time["2"], c("2" = 0))
})

test_that("transition-side states differ while accrued occupation is continuous", {
  h <- problem_event_history(rep("a", 3), c(2, 5, 8), c(1, 2, 3))
  problem <- compile_problem(h, list(home = 0, alive = c(0, 1), dead = 2), 10)
  prepared <- prepare_occupation_history(h, problem)
  right <- evaluate_occupation_history(prepared, rep(1L, 4), c(0, 2, 5, 8), "right")
  left <- evaluate_occupation_history(prepared, rep(1L, 4), c(0, 2, 5, 8), "left")
  expect_identical(right$accrued, left$accrued)
  expect_identical(problem$states[right$state], c(0, 1, 0, 2))
  expect_identical(problem$states[left$state], c(0, 0, 1, 0))
  expect_false(left$absorbed[4])
  expect_true(right$absorbed[4])
  mixed <- evaluate_occupation_history(prepared, rep(1L, 4), c(5, 2, 5, 2),
                                       c("left", "right", "right", "left"))
  expect_identical(problem$states[mixed$state], c(1, 1, 0, 0))
  expect_identical(mixed$accrued, right$accrued[c(3, 2, 3, 2), , drop = FALSE])
})

test_that("censoring boundaries never lead to unobserved extrapolation", {
  h <- problem_event_history(c("a", "a", "b"), c(2, 2, 0), c(1, 0, 0))
  problem <- compile_problem(h, list(home = 0, hospital = 1), 5)
  prepared <- prepare_occupation_history(h, problem)
  right <- evaluate_occupation_history(prepared, c(1, 2), c(2, 0))
  left <- evaluate_occupation_history(prepared, c(1, 2), c(2, 0), "left")
  expect_identical(problem$states[right$state], c(1, 0))
  expect_identical(problem$states[left$state], c(0, 0))
  expect_identical(right$accrued, left$accrued)
  expect_equal(unname(right$accrued), rbind(c(2, 0), c(0, 0)))
  expect_error(evaluate_occupation_history(prepared, 1L, 2.1),
               "exceeds observed censoring.*subject 'a'")
  expect_error(evaluate_occupation_history(prepared, 2L, 0.1),
               "exceeds observed censoring.*subject 'b'")
  expect_error(evaluate_occupation_history(prepared, 1L, 100), "censoring")

  later <- problem_event_history("later", 20, 0)
  short <- compile_problem(later, list(home = 0), 5)
  at_horizon <- evaluate_occupation_history(prepare_occupation_history(later, short),
                                            1L, 100)
  expect_equal(unname(at_horizon$accrued[1, 1]), 5)
  expect_equal(later$subjects$observation_end, 20)
})

test_that("horizon clipping keeps a jump at tau but excludes its later reward", {
  h <- problem_event_history(rep("a", 3), c(2, 5, 8), c(1, 2, 3))
  problem <- compile_problem(h, list(home = 0, hospital = 1, dead = 2), 5)
  prepared <- prepare_occupation_history(h, problem)
  right <- evaluate_occupation_history(prepared, 1L, 12)
  left <- evaluate_occupation_history(prepared, 1L, 12, "left")
  expect_equal(unname(right$accrued), matrix(c(2, 3, 0), 1))
  expect_identical(right$accrued, left$accrued)
  expect_equal(problem$states[right$state], 0)
  expect_equal(problem$states[left$state], 1)
  expect_true(right$at_horizon)
  expect_equal(problem$working_graph$edges$n_events, c(1L, 1L, 1L))
})

test_that("snapshots cannot reveal events occurring after their requested side", {
  h <- problem_event_history(c("a", "b", "b", "a", "b"),
                             c(2, 2, 6, 7, 9), c(1, 1, 2, 3, 0))
  problem <- compile_problem(h, list(home = 0, hospital = 1, dead = 2), 10)
  prepared <- prepare_occupation_history(h, problem)
  for (side in c("right", "left")) {
    a <- evaluate_occupation_history(prepared, rep(1L, 4), c(0, 1, 2, 4), side)
    b <- evaluate_occupation_history(prepared, rep(2L, 4), c(0, 1, 2, 4), side)
    expect_identical(a$accrued, b$accrued)
    expect_identical(a$state, b$state)
    expect_identical(a$absorbed, b$absorbed)
  }
  reloaded <- unserialize(serialize(prepared, NULL))
  expect_identical(evaluate_occupation_history(reloaded, c(2, 1), c(4, 1)),
                   evaluate_occupation_history(prepared, c(2, 1), c(4, 1)))
})

test_that("structural envelopes contain an exact complete-path target at earlier times", {
  h <- problem_event_history(rep("a", 3), c(2, 5, 8), c(1, 2, 3))
  p <- compile_problem(h, list(home = 0, hospital = 1, alive = c(0, 1), dead = 2), 10)
  prepared <- prepare_occupation_history(h, p)
  q <- evaluate_occupation_history(prepared, rep(1L, 8), c(0, 1, 2, 4, 5, 7, 8, 10))
  lower <- q$accrued + (10 - q$time) * p$structural_support$reward_min[q$state, ]
  upper <- q$accrued + (10 - q$time) * p$structural_support$reward_max[q$state, ]
  truth <- q$accrued[rep(8L, 8L), ]
  expect_true(all(lower <= truth))
  expect_true(all(upper >= truth))
  expect_identical(lower[7, ], truth[7, ])
  expect_identical(upper[7, ], truth[7, ])
})

test_that("occupation query validation preserves shapes and never recycles identifiers", {
  h <- problem_event_history("a", 5, 0)
  p <- compile_problem(h, list(home = 0), 5)
  prepared <- prepare_occupation_history(h, p)
  one <- evaluate_occupation_history(prepared, 1L, 2)
  expect_identical(dim(one$accrued), c(1L, 1L))
  empty <- evaluate_occupation_history(prepared, integer(), numeric())
  expect_identical(dim(empty$accrued), c(0L, 1L))
  expect_length(empty$state, 0)
  expect_error(evaluate_occupation_history(prepared, 1L, c(1, 2)), "one valid")
  expect_error(evaluate_occupation_history(prepared, 0L, 1), "valid internal")
  expect_error(evaluate_occupation_history(prepared, 1.5, 1), "valid internal")
  expect_error(evaluate_occupation_history(prepared, 1e100, 1), "valid internal")
  expect_error(evaluate_occupation_history(prepared, 1L, Inf), "finite")
  expect_error(evaluate_occupation_history(prepared, 1L, -1), "nonnegative")
  expect_error(evaluate_occupation_history(prepared, 1L, 1, "before"), "side")
  expect_error(evaluate_occupation_history(prepared, c(1, 1), c(1, 2),
                                          c("right", "left", "right")), "side")
})

test_that("a single initially absorbing state has deterministic full-horizon occupation", {
  h <- problem_event_history("a", 20, 0, initial_state = "dead", absorbing = "dead",
                             states = "dead", status_map = NULL)
  p <- compile_problem(h, list(dead = "dead"), 10)
  expect_equal(p$working_graph$origin_risk_time, c(dead = 0))
  expect_equal(nrow(p$working_graph$edges), 0)
  expect_identical(dim(p$structural_support$reachability), c(1L, 1L))
  expect_true(p$structural_support$reachability[1, 1])
  expect_identical(p$structural_support$lower, p$structural_support$upper)
  q <- evaluate_occupation_history(prepare_occupation_history(h, p),
                                    rep(1L, 3), c(0, 2, 10))
  expect_equal(unname(q$accrued), matrix(c(0, 2, 10), 3L))
  expect_identical(q$state, rep(1L, 3))
})
