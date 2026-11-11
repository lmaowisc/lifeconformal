#' Normalize observed histories to indexed intervals
#'
#' @param input A life_input object returned by resolve_input().
#' @param initial_state Common state, one value per record, one-sided column
#'   formula, or NULL to use state zero for event data and infer explicit origins.
#' @param status_map Named vector mapping nonzero event codes to state labels.
#' @param absorbing Declared absorbing-state labels; an empty vector means none.
#' @param states Optional fixed state dictionary used for later samples.
#' @returns A life_history_data object with subjects, intervals, baseline,
#'   row_index, schema, diagnostics, and a NULL x pending preprocessing.
#' @details A sole censoring observation at time zero is allowed. Genuine
#'   transitions are strictly later than zero. A positive-time transition may be
#'   followed immediately by terminal censoring. A subject initially absorbed
#'   may have one follow-up record; its path is complete at zero, so the record
#'   does not create censoring risk. Actual observation endpoints are retained.
#' @noRd
normalize_histories <- function(input, initial_state = NULL, status_map = NULL,
                                absorbing, states = NULL) {
  if (!inherits(input, "life_input")) {
    stop("input must be returned by resolve_input().", call. = FALSE)
  }
  if (missing(absorbing) || is.null(absorbing)) {
    stop("Declare absorbing states explicitly; use integer(0) for none.",
         call. = FALSE)
  }
  state_spec <- resolve_states(input, initial_state, status_map, absorbing, states)
  history <- if (input$response$format == "event") {
    events_to_intervals(input, state_spec)
  } else {
    transitions_to_intervals(input, state_spec)
  }
  validate_histories(history)
}

#' Validate and normalize state labels
#' @noRd
state_labels <- function(x, name, allow_missing = FALSE, allow_empty = FALSE,
                          id = NULL) {
  if (!length(x) && allow_empty) return(character())
  if (allow_missing && is.logical(x) && all(is.na(x))) {
    return(rep(NA_character_, length(x)))
  }
  if (!is.atomic(x) || !is.null(dim(x)) || is.complex(x) ||
      !(is.numeric(x) || is.character(x) || is.factor(x))) {
    stop(name, " must contain numeric, character or factor state labels.",
         call. = FALSE)
  }
  invalid <- if (allow_missing) rep(FALSE, length(x)) else is.na(x)
  if (is.numeric(x)) invalid <- invalid | (!is.finite(x) & !is.na(x))
  if (any(invalid)) {
    message <- paste0(name, " contains missing or nonfinite labels")
    row <- which(invalid)[1L]
    if (length(id) == length(x)) history_error(message, id[row], row)
    stop(message, ".", call. = FALSE)
  }
  if (is.factor(x)) {
    old_names <- names(x)
    x <- as.character(x)
    names(x) <- old_names
  }
  if (any(!nzchar(as.character(x)[!is.na(x)]))) {
    row <- which(!is.na(x) & !nzchar(as.character(x)))[1L]
    message <- paste0(name, " contains an empty label")
    if (length(id) == length(x)) history_error(message, id[row], row)
    stop(message, ".", call. = FALSE)
  }
  if (is.numeric(x) && anyDuplicated(as.character(unique(x[!is.na(x)])))) {
    stop(name, " has distinct numeric codes with indistinguishable text labels; use character codes.",
         call. = FALSE)
  }
  x
}

#' Resolve labels and initial states without imposing a state ordering
#' @noRd
resolve_states <- function(raw, initial_state = NULL, status_map = NULL,
                           absorbing, states = NULL) {
  d <- raw$response$data
  n <- nrow(d)
  ids <- unique(raw$id)
  subject_index <- match(raw$id, ids)
  first <- match(seq_along(ids), subject_index)
  row_index <- unname(split(seq_len(n), factor(subject_index,
                                              levels = seq_along(ids))))
  from <- if (raw$response$format == "transition")
    state_labels(d$from, "from", id = raw$id) else NULL
  initial_rule <- list(type = "common", value = 0)
  if (inherits(initial_state, "formula")) {
    if (length(initial_state) != 2L || !is.symbol(initial_state[[2L]])) {
      stop("initial_state must be a one-sided column formula with data.",
           call. = FALSE)
    }
    column <- as.character(initial_state[[2L]])
    if (is.null(raw$source_data)) {
      stop("The initial_state column rule requires data containing '", column,
           "'; for vector input, supply initial_state repeated once per record.",
           call. = FALSE)
    }
    if (!column %in% names(raw$source_data)) {
      stop("Initial-state column not found: ", column, ".", call. = FALSE)
    }
    initial_state <- raw$source_data[[column]]
    initial_rule <- list(type = "column", column = column)
  } else if (is.null(initial_state)) {
    if (raw$response$format == "transition") {
      initial_state <- from[first]
      initial_state <- initial_state[subject_index]
      initial_rule <- list(type = "first_from")
    } else {
      initial_state <- 0
    }
  } else if (length(initial_state) == 1L) {
    initial_rule <- list(type = "common", value = initial_state)
  } else {
    initial_rule <- list(type = "vector")
  }
  initial_state <- state_labels(initial_state, "initial_state", id = raw$id)
  if (length(initial_state) == 1L) initial_state <- rep(initial_state, n)
  if (length(initial_state) != n) {
    stop("initial_state must be common or have one value per record.", call. = FALSE)
  }
  if (any(initial_state != initial_state[first][subject_index])) {
    j <- which(initial_state != initial_state[first][subject_index])[1L]
    history_error("Initial state changes within subject", raw$id[j], j)
  }
  absorbing <- state_labels(absorbing, "absorbing", allow_empty = TRUE)
  if (anyDuplicated(as.character(absorbing))) {
    stop("absorbing contains duplicate state labels.", call. = FALSE)
  }
  if (raw$response$format == "event") {
    status <- state_labels(d$status, "status", id = raw$id)
    event <- as.character(status) != "0"
    destination <- status
    destination[!event] <- NA
    if (!is.null(status_map)) {
      if (!length(status_map) || is.null(names(status_map)) ||
          anyNA(names(status_map)) || any(!nzchar(names(status_map))) ||
          anyDuplicated(names(status_map)) || "0" %in% names(status_map)) {
        stop("status_map needs unique, nonempty, nonzero event-code names.",
             call. = FALSE)
      }
      status_map <- state_labels(status_map, "status_map")
      key <- match(as.character(status), names(status_map))
      if (any(event & is.na(key))) {
        j <- which(event & is.na(key))[1L]
        history_error(paste0("Unknown event code '", status[j], "'"), raw$id[j], j)
      }
      destination <- unname(status_map[key])
      destination[!event] <- NA
    }
    from <- NULL
  } else {
    if (!is.null(status_map)) {
      stop("status_map applies only to event-status input.", call. = FALSE)
    }
    if (!(is.numeric(d$event) || is.logical(d$event))) {
      stop("event must contain only zero and one (or FALSE and TRUE).", call. = FALSE)
    }
    invalid_event <- is.na(d$event) | !d$event %in% c(0, 1)
    if (any(invalid_event)) {
      row <- which(invalid_event)[1L]
      history_error("event must contain only zero and one (or FALSE and TRUE)",
                    raw$id[row], row)
    }
    event <- as.logical(d$event)
    destination <- state_labels(d$to, "to", allow_missing = TRUE, id = raw$id)
    invalid_to <- (event & is.na(destination)) | (!event & !is.na(destination))
    if (any(invalid_to)) {
      row <- which(invalid_to)[1L]
      history_error("to must be present for transitions and missing for terminal censoring",
                    raw$id[row], row)
    }
  }
  sources <- list(initial_state[first], destination[event], from,
                  unname(status_map), absorbing)
  sources <- sources[lengths(sources) > 0L]
  # Character encodings are also used for status-map matching. Fail rather than
  # silently merging two distinct numeric labels with the same printed value.
  numeric_sources <- Filter(is.numeric, c(sources, list(states)))
  numeric_labels <- unique(unlist(numeric_sources, use.names = FALSE))
  numeric_labels <- numeric_labels[!is.na(numeric_labels)]
  if (anyDuplicated(as.character(numeric_labels))) {
    stop("Distinct numeric states have indistinguishable text labels; use distinct character labels.",
         call. = FALSE)
  }
  found <- do.call(c, sources)
  found <- unique(found[!is.na(found)])
  if (is.null(states)) {
    states <- found
  } else {
    states <- state_labels(states, "states")
    if (!length(states) || anyDuplicated(as.character(states))) {
      stop("states must be a nonempty dictionary of unique labels.", call. = FALSE)
    }
    absent <- setdiff(as.character(found), as.character(states))
    if (length(absent)) {
      message <- paste0("State label(s) absent from the fixed dictionary: ",
                         paste(absent, collapse = ", "))
      fields <- list(initial_state = initial_state, from = from)
      fields[[if (raw$response$format == "event") "status/status_map" else "to"]] <- destination
      for (name in names(fields)) {
        invalid <- as.character(fields[[name]]) %in% absent
        if (any(invalid)) {
          row <- which(invalid)[1L]
          history_error(paste0(message, "; check ", name), raw$id[row], row)
        }
      }
      stop(message, ".", call. = FALSE)
    }
  }
  keys <- as.character(states)
  list(labels = states, keys = keys, event = event,
       initial = match(as.character(initial_state[first]), keys),
       from = match(as.character(from), keys),
       to = match(as.character(destination), keys),
       absorbing = match(as.character(absorbing), keys), status_map = status_map,
       initial_rule = initial_rule, ids = ids, subject_index = subject_index,
       first = first, row_index = row_index)
}

#' Convert destination events to origin-state intervals
#' @noRd
events_to_intervals <- function(raw, state_spec) {
  d <- raw$response$data
  check_history_times(d$time, "time", raw$id)
  n <- nrow(d)
  start <- numeric(n)
  from <- integer(n)
  for (i in seq_along(state_spec$row_index)) {
    rows <- state_spec$row_index[[i]]
    prior_rows <- rows[-length(rows)]
    start[rows] <- c(0, d$time[prior_rows])
    from[rows] <- c(state_spec$initial[i], state_spec$to[prior_rows])
  }
  intervals <- data.frame(subject_index = state_spec$subject_index,
                          start = start, stop = d$time, from = from,
                          to = state_spec$to, event = state_spec$event,
                          terminal_censor = !state_spec$event,
                          original_row = seq_len(n))
  assemble_history(raw, state_spec, intervals)
}

#' Normalize already explicit origin-state intervals
#' @noRd
transitions_to_intervals <- function(raw, state_spec) {
  d <- raw$response$data
  check_history_times(d$start, "start", raw$id)
  check_history_times(d$stop, "stop", raw$id)
  intervals <- data.frame(subject_index = state_spec$subject_index,
                          start = d$start, stop = d$stop, from = state_spec$from,
                          to = state_spec$to, event = state_spec$event,
                          terminal_censor = !state_spec$event,
                          original_row = seq_len(nrow(d)))
  assemble_history(raw, state_spec, intervals)
}

#' Assemble subject and interval tables while preserving original row order
#' @noRd
assemble_history <- function(raw, state_spec, intervals) {
  index <- state_spec$subject_index
  first <- state_spec$first
  last <- vapply(state_spec$row_index, function(z) z[length(z)], integer(1))
  initially_absorbed <- state_spec$initial %in% state_spec$absorbing
  absorbed <- initially_absorbed |
    (intervals$event[last] & intervals$to[last] %in% state_spec$absorbing)
  absorption_time <- rep(NA_real_, length(first))
  absorption_time[absorbed] <- intervals$stop[last[absorbed]]
  absorption_time[initially_absorbed] <- 0
  absorption_state <- rep(NA_integer_, length(first))
  absorption_state[absorbed] <- intervals$to[last[absorbed]]
  absorption_state[initially_absorbed] <- state_spec$initial[initially_absorbed]
  intervals$terminal_censor[initially_absorbed[index]] <- FALSE
  subjects <- data.frame(subject_index = seq_along(first), id = state_spec$ids,
                         initial_state = state_spec$initial,
                         observation_end = intervals$stop[last],
                         censoring = !absorbed,
                         absorption_state = absorption_state,
                         absorption_time = absorption_time)
  # Complete paths at origin are distinct from reported observation endpoints.
  subjects$path_completion <- absorption_time
  for (name in names(raw$baseline)) {
    value <- raw$baseline[[name]]
    different <- value != value[first][index]
    if (any(different)) {
      j <- which(different)[1L]
      history_error(paste0("Baseline column '", name, "' changes within subject"),
                    raw$id[j], j)
    }
  }
  baseline <- raw$baseline[first, , drop = FALSE]
  rownames(baseline) <- NULL
  schema <- raw$schema
  schema$state_labels <- state_spec$labels
  schema$state_keys <- state_spec$keys
  schema$absorbing <- state_spec$absorbing
  schema$status_map <- state_spec$status_map
  schema$initial_state_rule <- state_spec$initial_rule
  schema$tie_convention <- "transition_before_terminal_censoring"
  schema$origin_convention <- "no_transition_at_zero; sole_censoring_at_zero_allowed"
  structure(list(subjects = subjects, intervals = intervals, baseline = baseline,
                 x = NULL, row_index = state_spec$row_index, schema = schema,
                 diagnostics = NULL), class = "life_history_data")
}

#' Check canonical paths once at the data boundary
#' @noRd
validate_histories <- function(history, problem = NULL) {
  if (!inherits(history, "life_history_data")) {
    stop("Expected a life_history_data object.", call. = FALSE)
  }
  d <- history$intervals
  subjects <- history$subjects
  absorbing <- history$schema$absorbing
  for (i in seq_len(nrow(subjects))) {
    rows <- history$row_index[[i]]
    current <- subjects$initial_state[i]
    previous_time <- 0
    ended <- FALSE
    initial_absorbed <- current %in% absorbing
    for (j in rows) {
      fail <- function(message) history_error(message, subjects$id[i], d$original_row[j])
      if (ended) fail(paste0("Record appears after censoring or absorbing entry; ",
                            "that entry must be the subject's last record"))
      if (d$start[j] != previous_time) {
        fail(paste0("Intervals must start at zero and remain contiguous in supplied order; ",
          "start is ", format(d$start[j], digits = 17), ", expected ",
          format(previous_time, digits = 17)))
      }
      if (d$stop[j] < d$start[j]) fail(paste0("Times are not chronological; ",
        if (history$schema$input_format == "event") "time" else "stop",
        " is ", format(d$stop[j], digits = 17), ", before ",
        format(d$start[j], digits = 17)))
      if (is.na(d$from[j]) || d$from[j] != current) fail(paste0(
        "Origin state breaks path continuity; from must equal the current state '",
        history$schema$state_labels[current], "'"))
      if (d$event[j]) {
        if (current %in% absorbing) fail("Transition leaves a declared absorbing state")
        if (d$stop[j] <= d$start[j]) {
          fail("Genuine transitions must occur after zero and have positive sojourn duration")
        }
        if (d$to[j] == current) fail("A transition must enter a different state")
        current <- d$to[j]
        ended <- current %in% absorbing
      } else {
        zero_at_origin <- length(rows) == 1L && d$stop[j] == 0
        terminal_tie <- j != rows[1L] && d$stop[j] == previous_time
        if (d$stop[j] == d$start[j] && !zero_at_origin && !terminal_tie) {
          fail("Zero-duration intervals require terminal censoring after a transition")
        }
        ended <- TRUE
      }
      previous_time <- d$stop[j]
    }
    if (initial_absorbed && length(rows) != 1L) {
      history_error("An initially absorbed subject may have only one follow-up record",
                    subjects$id[i], rows[1L])
    }
    if (!ended) {
      history_error("A nonabsorbing subject needs a final censoring record (status = 0 or event = 0)",
                    subjects$id[i], rows[length(rows)])
    }
  }
  if (!is.null(problem) && !identical(history$schema$state_keys,
                                     as.character(problem$states))) {
    stop("The history and problem state dictionaries differ.", call. = FALSE)
  }
  history$diagnostics <- list(n_subjects = nrow(subjects), n_records = nrow(d),
                              n_transitions = sum(d$event),
                              n_censored = sum(subjects$censoring),
                              n_absorbed = sum(!subjects$censoring),
                              n_origin_censored = sum(subjects$censoring &
                                                        subjects$observation_end == 0),
                              n_initially_absorbed = sum(subjects$absorption_time == 0,
                                                        na.rm = TRUE),
                              rows_omitted = 0L)
  history
}

#' Check finite numeric elapsed times
#' @noRd
check_history_times <- function(x, name, id = NULL) {
  message <- paste0(name, " must contain finite, nonnegative numeric elapsed times")
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)))
    stop(message, ".", call. = FALSE)
  invalid <- is.na(x) | !is.finite(x) | x < 0
  if (any(invalid)) {
    row <- which(invalid)[1L]
    if (length(id) == length(x)) history_error(message, id[row], row)
    stop(message, ".", call. = FALSE)
  }
  invisible(NULL)
}

#' Report a path error using the original subject and row
#' @noRd
history_error <- function(message, id, row) {
  stop(message, " (subject '", as.character(id), "', row ", row, ").", call. = FALSE)
}
