#' Validate shared numerical and fitting controls
#' @noRd
resolve_life_controls <- function(control = list()) {
  defaults <- list(scale_floor = 0.01, survival_floor = 1e-8,
                   grid_size = 128L, max_nodes = 100000L,
                   probability_tolerance = 0.01, calibration_accuracy = "auto",
                   draws = 1000L, keep_contributions = FALSE,
                   markov_occupation = "general")
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || anyNA(names(control)) || any(!nzchar(names(control))) ||
       anyDuplicated(names(control))))) {
    stop("control must be a list with unique, nonempty names.", call. = FALSE)
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) stop("Unknown control(s): ",
                            paste(unknown, collapse = ", "), ".", call. = FALSE)
  defaults[names(control)] <- control
  if (!is.character(defaults$markov_occupation) || length(defaults$markov_occupation) != 1L ||
      !is.null(dim(defaults$markov_occupation)) || is.na(defaults$markov_occupation) ||
      !defaults$markov_occupation %in% c("general", "condition_on_alive"))
    stop("markov_occupation must be 'general' or 'condition_on_alive'.", call. = FALSE)
  if (!is.character(defaults$calibration_accuracy) ||
      length(defaults$calibration_accuracy) != 1L ||
      !is.null(dim(defaults$calibration_accuracy)) ||
      is.na(defaults$calibration_accuracy) ||
      !defaults$calibration_accuracy %in% c("auto", "probability", "threshold"))
    stop("calibration_accuracy must be 'auto', 'probability', or 'threshold'.", call. = FALSE)
  defaults$calibration_accuracy <- unname(defaults$calibration_accuracy)
  for (key in c("scale_floor", "survival_floor", "probability_tolerance")) {
    z <- defaults[[key]]
    if (!is.numeric(z) || is.complex(z) || !is.null(dim(z)) ||
        length(z) != 1L || !is.finite(z) ||
        z <= 0 || z > 1) stop(key, " must be in (0, 1].", call. = FALSE)
  }
  for (key in c("grid_size", "max_nodes", "draws")) {
    z <- defaults[[key]]
    if (!is.numeric(z) || is.complex(z) || !is.null(dim(z)) ||
        length(z) != 1L || !is.finite(z) ||
        z < 2 || z != floor(z) || z > .Machine$integer.max) {
      stop(key, " must be an integer of at least two.", call. = FALSE)
    }
    defaults[[key]] <- as.integer(z)
  }
  if (!is.logical(defaults$keep_contributions) ||
      !is.null(dim(defaults$keep_contributions)) ||
      length(defaults$keep_contributions) != 1L ||
      is.na(defaults$keep_contributions)) {
    stop("keep_contributions must be TRUE or FALSE.", call. = FALSE)
  }
  defaults
}

#' Evaluate randomized fitting without changing the caller's explicit seed
#' @noRd
with_life_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) {
    stop("seed must be NULL or a nonnegative integer.", call. = FALSE)
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  })
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(as.integer(seed))
  force(code)
}

#' Construct a common response for the public vector interface
#' @noRd
resolve_public_response <- function(formula, response, time, status, start,
                                    stop, from, to, event) {
  supplied <- any(!vapply(list(time, status, start, stop, from, to, event),
                          is.null, logical(1)))
  if (!is.null(formula)) {
    if (!is.null(response) || supplied)
      stop("Do not combine a formula with response vectors.", call. = FALSE)
    return(NULL)
  }
  if (!is.null(response)) {
    if (supplied) stop("Use response or individual response vectors, not both.",
                       call. = FALSE)
    return(response)
  }
  LifeHist(time = time, status = status, start = start, stop = stop,
           from = from, to = to, event = event)
}

#' Resolve a previously stored initial-state rule
#' @noRd
stored_initial_state <- function(schema, supplied = NULL) {
  if (!is.null(supplied)) return(supplied)
  rule <- schema$initial_state_rule
  if (is.null(rule)) rule <- schema$initial_rule
  if (rule$type == "common") return(rule$value)
  if (rule$type == "column")
    return(stats::as.formula(call("~", as.name(rule$column)), env = baseenv()))
  if (rule$type == "first_from") return(NULL)
  stop("Supply initial_state: training used a subject-specific state vector.",
       call. = FALSE)
}
