#' Compile occupation targets using actual state labels
#'
#' @param targets Nonempty named list of nonempty sets of state labels.
#' @param states Ordered dictionary of actual state labels, not integer indices.
#' @param tau Finite, strictly positive common prediction horizon.
#' @returns A life_targets object with target-by-state incidence A, its
#'   state-by-target transpose rewards, and a conservative joint support
#'   representation. In that representation z = A %*% v, v >= 0 and
#'   sum(v) = tau. It retains overlap and nesting of the targets.
#' @details A state may appear in several targets. Repeated labels within one
#'   target are removed because a target is a set. All numerical arrays retain
#'   two dimensions when there is only one state or target.
#' @noRd
compile_targets <- function(targets, states, tau) {
  states <- check_problem_states(states)
  check_problem_tau(tau)
  if (!is.list(targets) || is.data.frame(targets) || !length(targets) ||
      is.null(names(targets)) || anyNA(names(targets)) ||
      any(!nzchar(trimws(names(targets)))) || anyDuplicated(names(targets))) {
    stop("targets must be a nonempty list with unique, nonempty names.",
         call. = FALSE)
  }
  keys <- as.character(states)
  indices <- vector("list", length(targets))
  names(indices) <- names(targets)
  A <- matrix(0, length(targets), length(states),
              dimnames = list(names(targets), keys))
  for (j in seq_along(targets)) {
    label <- paste0("Target '", names(targets)[j], "'")
    values <- state_labels(targets[[j]], label)
    if (!length(values)) stop(label, " must be nonempty.", call. = FALSE)
    index <- match(as.character(values), keys)
    if (anyNA(index)) {
      stop(label, " has unknown state label(s): ",
           paste(unique(as.character(values[is.na(index)])), collapse = ", "),
           ".", call. = FALSE)
    }
    indices[[j]] <- unique(index)
    targets[[j]] <- states[indices[[j]]]
    A[j, indices[[j]]] <- 1
  }
  rewards <- t(A)
  lower <- tau * apply(A, 1L, min)
  upper <- tau * apply(A, 1L, max)
  support <- list(type = "state_occupation_simplex_image", A = A,
                  state_time_total = tau, lower = lower, upper = upper,
                  outer = TRUE, closed = TRUE)
  structure(list(states = states, state_keys = keys, names = names(targets),
                 sets = targets, indices = indices, n_targets = length(targets),
                 A = A, rewards = rewards, tau = tau, support = support),
            class = "life_targets")
}

#' Summarize the training-observed working transition graph
#'
#' @param history Validated canonical life_history_data.
#' @param states The unchanged canonical state dictionary.
#' @returns A life_working_graph containing edges in first-observed order,
#'   event counts and total origin risk duration. Edges use internal indices;
#'   the dictionary and label columns retain the corresponding actual labels.
#' @details Risk durations use actual observed follow-up, before any
#'   horizon-specific fitting view. An initially absorbed subject contributes
#'   zero risk even when its reported follow-up endpoint is later than zero.
#'   Aggregation is linear in the number of records, apart from constructing
#'   the small state adjacency matrix. Missing edges are working-model zeros,
#'   never declarations of absorption or structural impossibility.
#' @noRd
build_working_graph <- function(history, states = history$schema$state_labels) {
  check_problem_history(history, states)
  d <- history$intervals
  keys <- as.character(states)
  K <- length(states)
  duration <- d$stop - d$start
  duration[d$from %in% history$schema$absorbing] <- 0
  origin_risk_time <- setNames(numeric(K), keys)
  totals <- rowsum(duration, d$from, reorder = FALSE)
  origin_risk_time[as.integer(rownames(totals))] <- totals[, 1L]

  event_rows <- which(d$event)
  # Double keys avoid integer overflow without materializing a K-by-K count
  # array. The state dictionary already supplies consecutive integer indices.
  key <- (as.double(d$from[event_rows]) - 1) * K + d$to[event_rows]
  first <- !duplicated(key)
  edge_rows <- event_rows[first]
  edges <- data.frame(edge_index = seq_along(edge_rows),
                      from = d$from[edge_rows], to = d$to[edge_rows],
                      n_events = tabulate(match(key, key[first]),
                                          nbins = length(edge_rows)))
  edges$origin_risk_time <- unname(origin_risk_time[edges$from])
  edges$from_label <- states[edges$from]
  edges$to_label <- states[edges$to]
  adjacency <- matrix(FALSE, K, K, dimnames = list(keys, keys))
  if (nrow(edges)) adjacency[cbind(edges$from, edges$to)] <- TRUE
  structure(list(states = states, edges = edges, adjacency = adjacency,
                 origin_risk_time = origin_risk_time,
                 absorbing = history$schema$absorbing,
                 scope = "all_observed_training_transitions"),
            class = "life_working_graph")
}

#' Compile structural reachability and conservative occupation support
#'
#' @param states Ordered dictionary of actual state labels.
#' @param absorbing Explicit absorbing-state labels; integer(0) means none.
#' @param allowed_transitions NULL for every off-diagonal transition except
#'   exits from absorbing states, or a from/to data frame or two-column matrix
#'   of actual labels. A logical adjacency matrix is also accepted, in state
#'   order when unnamed, or reordered using its row and column labels.
#' @param targets Named target sets or a compatible life_targets object.
#' @param tau Common prediction horizon.
#' @returns A life_structural_support with structural adjacency, reflexive
#'   reachability, reachable-state lists and state-by-target reward_min and
#'   reward_max. The same limits multiplied by tau give lower and upper bounds
#'   at the origin. The joint future support is the image of nonnegative
#'   remaining-state durations summing to remaining time, restricted to
#'   reachable states.
#' @details Bounds are closed conservative envelopes, not assertions that an
#'   endpoint is attained by a path with strictly positive sojourns. For a
#'   snapshot (state k, time t, accrued a), coordinate bounds are
#'   a + (tau-t)*reward_min[k,] and a + (tau-t)*reward_max[k,]. Identical reward
#'   bounds are exact for that component. Empirical working edges are never
#'   consulted. Reachability uses one adjacency traversal per state.
#' @noRd
build_structural_support <- function(states, absorbing, allowed_transitions = NULL,
                                     targets, tau) {
  states <- check_problem_states(states)
  targets <- resolve_problem_targets(targets, states, tau)
  if (missing(absorbing) || is.null(absorbing)) {
    stop("Declare absorbing states explicitly; use integer(0) for none.",
         call. = FALSE)
  }
  absorbing <- state_labels(absorbing, "absorbing", allow_empty = TRUE)
  if (anyDuplicated(as.character(absorbing))) {
    stop("absorbing contains duplicate state labels.", call. = FALSE)
  }
  keys <- as.character(states)
  absorbing_index <- match(as.character(absorbing), keys)
  if (anyNA(absorbing_index)) {
    stop("absorbing contains labels absent from states.", call. = FALSE)
  }
  K <- length(states)
  specified <- !is.null(allowed_transitions)
  if (!specified) {
    adjacency <- matrix(TRUE, K, K, dimnames = list(keys, keys))
    diag(adjacency) <- FALSE
    adjacency[absorbing_index, ] <- FALSE
  } else {
    adjacency <- compile_structural_adjacency(allowed_transitions, states)
    if (any(adjacency[absorbing_index, , drop = FALSE])) {
      stop("allowed_transitions permits an exit from a declared absorbing state.",
           call. = FALSE)
    }
  }
  adjacent <- lapply(seq_len(K), function(k) unname(which(adjacency[k, ])))
  reachability <- matrix(FALSE, K, K, dimnames = list(keys, keys))
  if (!specified) {
    reachability <- adjacency
    diag(reachability) <- TRUE
  } else {
    for (k in seq_len(K)) {
      visited <- logical(K)
      visited[k] <- TRUE
      queue <- integer(K)
      queue[1L] <- k
      head <- tail <- 1L
      while (head <= tail) {
        next_states <- adjacent[[queue[head]]]
        next_states <- next_states[!visited[next_states]]
        if (length(next_states)) {
          visited[next_states] <- TRUE
          queue[seq.int(tail + 1L, tail + length(next_states))] <- next_states
          tail <- tail + length(next_states)
        }
        head <- head + 1L
      }
      reachability[k, ] <- visited
    }
  }
  reachable <- lapply(seq_len(K), function(k) unname(which(reachability[k, ])))
  names(reachable) <- keys
  reward_min <- reward_max <- matrix(0, K, targets$n_targets,
                                     dimnames = list(keys, targets$names))
  for (k in seq_len(K)) {
    values <- targets$rewards[reachable[[k]], , drop = FALSE]
    reward_min[k, ] <- apply(values, 2L, min)
    reward_max[k, ] <- apply(values, 2L, max)
  }
  edge_indices <- which(adjacency, arr.ind = TRUE)
  edges <- data.frame(from = edge_indices[, 1L], to = edge_indices[, 2L],
                      row.names = NULL)
  structure(list(states = states, absorbing = absorbing_index,
                 specified = specified, adjacency = adjacency,
                 edges = edges, adjacent = adjacent,
                 reachable = reachable, reachability = reachability,
                 targets = targets, rewards = targets$rewards,
                 reward_min = reward_min, reward_max = reward_max,
                 lower = tau * reward_min, upper = tau * reward_max,
                 tau = tau, outer = TRUE, closed = TRUE,
                 bound_attainment = "not_assumed",
                 joint_support = targets$support),
            class = "life_structural_support")
}

#' Freeze the state, target and graph definition of a training problem
#'
#' @param history Validated canonical training histories.
#' @param targets Named occupation target sets, using actual state labels.
#' @param tau Common prediction horizon.
#' @param allowed_transitions Optional declared structural graph.
#' @param score Fixed score family identifier retained for later construction.
#' @returns A life_problem whose working graph and structural support have
#'   separate fields. Numerical state indices always refer to states in their
#'   original, frozen order. Training records contradicting structural support
#'   are rejected with original subject and row identifiers.
#' @noRd
compile_problem <- function(history, targets, tau, allowed_transitions = NULL,
                             score = "standardized_max") {
  check_problem_history(history, history$schema$state_labels)
  if (!is.character(score) || length(score) != 1L || is.na(score) ||
      !nzchar(score)) {
    stop("score must be one nonempty score-family identifier.", call. = FALSE)
  }
  states <- history$schema$state_labels
  compiled <- compile_targets(targets, states, tau)
  support <- build_structural_support(
    states, states[history$schema$absorbing], allowed_transitions, compiled, tau
  )
  check_structural_transitions(history, support)
  structure(list(states = states,
                 state_index = setNames(seq_along(states), as.character(states)),
                 initial_state_rule = history$schema$initial_state_rule,
                 absorbing = history$schema$absorbing,
                 absorbing_labels = states[history$schema$absorbing],
                 targets = compiled, A = compiled$A, rewards = compiled$rewards,
                 tau = tau, score = score,
                 working_graph = build_working_graph(history, states),
                 structural_support = support),
            class = "life_problem")
}

#' Prepare reusable occupation prefixes for observed-history queries
#'
#' @param history Validated canonical histories using the frozen dictionary.
#' @param problem Compiled life_problem.
#' @returns A life_occupation_history containing one target prefix per record
#'   and event-row references per subject. It stores no copied path prefixes.
#' @details Preparation costs O(m*d) time and storage for m records and d
#'   targets. The interval table remains in supplied order, including when
#'   subjects are interleaved. Its rewards cover actual observations; evaluation
#'   is clipped to the frozen horizon. Structural contradictions are checked
#'   once here, rather than inside every requested query.
#' @noRd
prepare_occupation_history <- function(history, problem) {
  if (!inherits(problem, "life_problem")) {
    stop("problem must be returned by compile_problem().", call. = FALSE)
  }
  check_problem_history(history, problem$states)
  if (!setequal(history$schema$absorbing, problem$absorbing)) {
    stop("The history and problem absorbing-state declarations differ.",
         call. = FALSE)
  }
  check_structural_transitions(history, problem$structural_support)
  d <- history$intervals
  rewards <- problem$rewards
  contribution <- (d$stop - d$start) * rewards[d$from, , drop = FALSE]
  prefix <- matrix(0, nrow(d), ncol(rewards),
                    dimnames = list(NULL, colnames(rewards)))
  event_rows <- vector("list", nrow(history$subjects))
  for (i in seq_along(history$row_index)) {
    rows <- history$row_index[[i]]
    for (j in seq_len(ncol(rewards))) {
      prefix[rows, j] <- cumsum(contribution[rows, j])
      # Until a positive-duration reward change, occupation is exactly the
      # elapsed time times its initial binary reward. Duration subtraction
      # followed by cumulative addition must not erase this invariant.
      initial_reward <- rewards[history$subjects$initial_state[i], j]
      changed <- d$stop[rows] > d$start[rows] &
        rewards[d$from[rows], j] != initial_reward
      constant <- cumsum(changed) == 0L
      prefix[rows[constant], j] <- d$stop[rows[constant]] * initial_reward
    }
    event_rows[[i]] <- rows[d$event[rows]]
  }
  structure(list(subjects = history$subjects, intervals = d,
                 row_index = history$row_index, event_rows = event_rows,
                 prefix = prefix, rewards = rewards, tau = problem$tau,
                 absorbing = problem$absorbing, states = problem$states),
            class = "life_occupation_history")
}

#' Evaluate accrued occupations and current states at aligned history queries
#'
#' @param prepared Result of prepare_occupation_history().
#' @param subject_index Internal subject indices, one per query.
#' @param time Finite nonnegative requested times, one per query. Values beyond
#'   tau are evaluated at tau. Scalar recycling is not performed.
#' @param side Through-time state (right) or state just before a transition
#'   (left). A scalar side applies to every query; an aligned side vector is
#'   also accepted. Accrued occupations are continuous at transition times.
#' @returns A life_history_snapshot list with aligned subject_index, id,
#'   requested_time, clipped time, side, state, accrued matrix, absorbed and
#'   at_horizon. An empty query returns correctly shaped empty fields.
#' @details Queries at censoring use the last observed state, including an
#'   event followed by censoring at that time. Queries beyond censoring are
#'   rejected rather than extrapolated. Known absorbing trajectories continue
#'   through tau, even when absorption occurs at zero. Event-prefix lookup
#'   costs O(log(m_i) + d) per arbitrary query, with subjects grouped once.
#'   Later observed events cannot affect an earlier snapshot.
#' @noRd
evaluate_occupation_history <- function(prepared, subject_index, time,
                                        side = c("right", "left")) {
  if (!inherits(prepared, "life_occupation_history")) {
    stop("prepared must be returned by prepare_occupation_history().", call. = FALSE)
  }
  n <- length(time)
  if (!is.numeric(subject_index) || !is.null(dim(subject_index)) ||
      is.complex(subject_index) || length(subject_index) != n ||
      anyNA(subject_index) || any(!is.finite(subject_index)) ||
      any(subject_index != floor(subject_index)) ||
      any(subject_index < 1 | subject_index > nrow(prepared$subjects))) {
    stop("subject_index must contain one valid internal subject index per query.",
         call. = FALSE)
  }
  if (!is.null(dim(time))) stop("time must be a vector.", call. = FALSE)
  check_history_times(time, "time")
  if (missing(side)) side <- "right"
  if (!is.character(side) || anyNA(side) ||
      !(length(side) == 1L || length(side) == n) ||
      any(!side %in% c("right", "left"))) {
    stop("side must be 'right' or 'left', common or one value per query.",
         call. = FALSE)
  }
  if (length(side) == 1L) side <- rep(side, n)
  subject_index <- as.integer(subject_index)
  subjects <- prepared$subjects
  requested_time <- time
  time <- pmin(time, prepared$tau)
  unavailable <- subjects$censoring[subject_index] &
    time > subjects$observation_end[subject_index]
  if (any(unavailable)) {
    q <- which(unavailable)[1L]
    stop("History query exceeds observed censoring for subject '",
         as.character(subjects$id[subject_index[q]]), "' (query ", q, ").",
         call. = FALSE)
  }
  state <- subjects$initial_state[subject_index]
  accrued <- matrix(0, n, ncol(prepared$rewards),
                     dimnames = list(NULL, colnames(prepared$rewards)))
  groups <- split(seq_len(n), subject_index)
  d <- prepared$intervals
  for (rows in groups) {
    i <- subject_index[rows[1L]]
    events <- prepared$event_rows[[i]]
    jumps <- d$stop[events]
    index <- findInterval(time[rows], jumps)
    had_event <- index > 0L
    query_rows <- rows[had_event]
    event <- events[index[had_event]]
    previous_time <- numeric(length(rows))
    if (any(had_event)) {
      accrued[query_rows, ] <- prepared$prefix[event, , drop = FALSE]
      previous_time[had_event] <- d$stop[event]
      state[query_rows] <- d$to[event]
    }
    accrued[rows, ] <- occupation_complete(time[rows], previous_time,
      prepared$rewards[state[rows], , drop = FALSE],
      accrued[rows, , drop = FALSE])
    # Occupation is continuous. Only the conditioning state changes sides;
    # using the same accrual avoids subtraction at an exact transition tie.
    before_jump <- had_event & side[rows] == "left"
    before_jump[had_event] <- before_jump[had_event] &
      time[rows[had_event]] == jumps[index[had_event]]
    if (any(before_jump)) {
      state[rows[before_jump]] <- d$from[events[index[before_jump]]]
    }
  }
  structure(list(subject_index = subject_index, id = subjects$id[subject_index],
                 requested_time = requested_time, time = time, side = side,
                 state = state, accrued = accrued,
                 absorbed = state %in% prepared$absorbing,
                 at_horizon = time == prepared$tau),
            class = "life_history_snapshot")
}

#' Validate the state dictionary without changing its order or labels
#' @noRd
check_problem_states <- function(states) {
  states <- state_labels(states, "states")
  if (!length(states) || anyDuplicated(as.character(states))) {
    stop("states must be a nonempty dictionary of unique labels.", call. = FALSE)
  }
  states
}

#' Validate a common occupation horizon
#' @noRd
check_problem_tau <- function(tau) {
  if (!is.numeric(tau) || is.complex(tau) || !is.null(dim(tau)) ||
      length(tau) != 1L || is.na(tau) || !is.finite(tau) || tau <= 0) {
    stop("tau must be one finite, strictly positive numeric horizon.", call. = FALSE)
  }
  invisible(NULL)
}

#' Check the immutable history dictionary at an internal boundary
#' @noRd
check_problem_history <- function(history, states) {
  if (!inherits(history, "life_history_data")) {
    stop("history must be a validated life_history_data object.", call. = FALSE)
  }
  states <- check_problem_states(states)
  if (!identical(history$schema$state_keys, as.character(states))) {
    stop("The history and problem state dictionaries differ.", call. = FALSE)
  }
  invisible(NULL)
}

#' Reuse compiled targets only for the identical dictionary and horizon
#' @noRd
resolve_problem_targets <- function(targets, states, tau) {
  check_problem_tau(tau)
  if (!inherits(targets, "life_targets")) {
    return(compile_targets(targets, states, tau))
  }
  if (!identical(targets$state_keys, as.character(states)) ||
      targets$tau != tau) {
    stop("Compiled targets use a different state dictionary or horizon.", call. = FALSE)
  }
  targets
}

#' Map a declared structural graph from actual labels to adjacency indices
#' @noRd
compile_structural_adjacency <- function(allowed_transitions, states) {
  K <- length(states)
  keys <- as.character(states)
  graph <- allowed_transitions
  if (is.matrix(graph) && is.logical(graph)) {
    if (!identical(dim(graph), c(K, K)) || anyNA(graph)) {
      stop("A logical allowed_transitions matrix must be K by K without missing values.",
           call. = FALSE)
    }
    labels <- dimnames(graph)
    if (!is.null(labels)) {
      if (is.null(labels[[1L]]) || is.null(labels[[2L]]) ||
          anyDuplicated(labels[[1L]]) || anyDuplicated(labels[[2L]]) ||
          !setequal(labels[[1L]], keys) || !setequal(labels[[2L]], keys)) {
        stop("Adjacency row and column names must match the state dictionary.",
             call. = FALSE)
      }
      graph <- graph[keys, keys, drop = FALSE]
    }
    dimnames(graph) <- list(keys, keys)
  } else {
    if (!(is.data.frame(graph) || is.matrix(graph)) || ncol(graph) != 2L) {
      stop("allowed_transitions must be a from/to table or a logical adjacency matrix.",
           call. = FALSE)
    }
    if (is.data.frame(graph) || !is.null(colnames(graph))) {
      if (!setequal(colnames(graph), c("from", "to"))) {
        stop("An allowed_transitions edge table must have columns from and to.",
             call. = FALSE)
      }
      graph <- graph[, c("from", "to"), drop = FALSE]
    }
    from <- state_labels(graph[, 1L], "allowed_transitions from", allow_empty = TRUE)
    to <- state_labels(graph[, 2L], "allowed_transitions to", allow_empty = TRUE)
    from <- match(as.character(from), keys)
    to <- match(as.character(to), keys)
    if (anyNA(from) || anyNA(to)) {
      stop("allowed_transitions contains an unknown state label.", call. = FALSE)
    }
    if (anyDuplicated((as.double(from) - 1) * K + to)) {
      stop("allowed_transitions contains duplicate edges.", call. = FALSE)
    }
    graph <- matrix(FALSE, K, K, dimnames = list(keys, keys))
    if (length(from)) graph[cbind(from, to)] <- TRUE
  }
  if (any(diag(graph))) {
    stop("allowed_transitions must omit self-transitions.", call. = FALSE)
  }
  graph
}

#' Reject observed transitions contradicted by declared structural information
#' @noRd
check_structural_transitions <- function(history, support) {
  d <- history$intervals
  rows <- which(d$event)
  allowed <- support$adjacency[cbind(d$from[rows], d$to[rows])]
  if (any(!allowed)) {
    j <- rows[which(!allowed)[1L]]
    history_error("Observed transition contradicts declared structural support",
                  history$subjects$id[d$subject_index[j]], d$original_row[j])
  }
  invisible(NULL)
}
