#' Describe observed life histories
#'
#' Construct the response used in a life-history formula. Supply either event
#' times and destination-status codes, or explicit start--stop transitions.
#' Subject identifiers and initial/absorbing states are supplied separately to
#' the analysis. This constructor preserves records; it does not sort them.
#'
#' @param time Finite, nonnegative numeric times since the prediction origin,
#'   one per event or final
#'   observation record. All times use the same unit.
#' @param status Destination-event codes. Zero denotes the end of observation;
#'   every other code denotes entry into a state. A separate status map is needed
#'   when a nonzero code represents a return to state zero.
#' @param start,stop Numeric interval endpoints for explicit transitions.
#' @param from,to Origin and destination state labels. For a terminal censoring
#'   interval, `to` must be missing.
#' @param event Logical or numeric zero/one transition indicator. Zero denotes
#'   the terminal observation interval.
#'
#' @details
#' Genuine transitions must occur after time zero. A subject censored at the
#' origin can have a single record at zero. A transition followed by censoring
#' at the same positive time is represented in that order. Other simultaneous
#' transitions are ambiguous and are rejected when histories are validated.
#'
#' `LifeHist()` checks the representation and vector lengths. Subject-level
#' chronology, state continuity, and baseline covariates are validated by the
#' analysis data parser. Numeric, character, or factor state codes are accepted;
#' factor levels are treated as labels, not integer codes.
#' In an analysis formula, response variables must be columns of `data`.
#' Supply response vectors directly to the analysis for vectors stored elsewhere.
#' Invalid records are rejected before model fitting. Errors identify the
#' offending field and, when available, the original subject and row.
#'
#' @returns A `life_response` list containing `format` and `data`. The data frame
#'   contains the supplied response vectors with their original row order.
#' @export
#' @examples
#' LifeHist(time = c(20, 35, 80), status = c(1, 2, 0))
#' LifeHist(start = c(0, 20, 35), stop = c(20, 35, 80),
#'          from = c(0, 1, 0), to = c(1, 0, NA), event = c(1, 1, 0))
LifeHist <- function(time = NULL, status = NULL, start = NULL, stop = NULL,
                     from = NULL, to = NULL, event = NULL) {
  event_input <- !is.null(time) || !is.null(status)
  transition_input <- any(!vapply(list(start, stop, from, to, event),
                                  is.null, logical(1)))
  if (event_input == transition_input) {
    base::stop("Supply either time/status or start/stop/from/to/event, not both.",
         call. = FALSE)
  }
  values <- if (event_input) list(time = time, status = status) else
    list(start = start, stop = stop, from = from, to = to, event = event)
  absent <- names(values)[vapply(values, is.null, logical(1))]
  if (length(absent)) {
    base::stop("Missing response field(s): ", paste(absent, collapse = ", "), ".",
         call. = FALSE)
  }
  vector_ok <- vapply(values, function(z) is.atomic(z) && is.null(dim(z)),
                      logical(1))
  if (any(!vector_ok)) {
    base::stop("Response fields must be vectors: ",
         paste(names(values)[!vector_ok], collapse = ", "), ".", call. = FALSE)
  }
  lengths <- vapply(values, length, integer(1))
  if (!lengths[1L] || any(lengths != lengths[1L])) {
    base::stop("Response vectors must have the same positive length; recycling is not allowed. Received ",
         paste(paste0(names(lengths), " = ", lengths), collapse = ", "), ".",
         call. = FALSE)
  }
  structure(list(format = if (event_input) "event" else "transition",
                 data = as.data.frame(values, stringsAsFactors = FALSE)),
            class = "life_response")
}

#' Check a response container before interpreting its rows
#' @noRd
validate_life_response <- function(response) {
  invalid <- function(detail) stop("Invalid LifeHist() response: ", detail,
    " Reconstruct it with LifeHist(time, status) or ",
    "LifeHist(start = start, stop = stop, from = from, to = to, event = event).",
    call. = FALSE)
  if (!is.list(response) || !is.character(response$format) ||
      length(response$format) != 1L || is.na(response$format) ||
      !response$format %in% c("event", "transition"))
    invalid("format must be 'event' or 'transition'.")
  d <- response$data
  if (!is.data.frame(d) || !nrow(d)) invalid("data must be a nonempty data frame.")
  if (anyNA(names(d)) || anyDuplicated(names(d)) || any(!nzchar(names(d))))
    invalid("response columns must have unique, nonempty names.")
  fields <- if (response$format == "event") c("time", "status") else
    c("start", "stop", "from", "to", "event")
  absent <- setdiff(fields, names(d))
  if (length(absent)) invalid(paste0("Missing response field(s): ",
                                    paste(absent, collapse = ", "), "."))
  valid <- vapply(fields, function(name) {
    z <- d[[name]]
    is.atomic(z) && is.null(dim(z)) && length(z) == nrow(d)
  }, logical(1))
  if (any(!valid)) invalid(paste0("Response fields must be vectors with one value per row: ",
                                  paste(fields[!valid], collapse = ", "), "."))
  invisible(NULL)
}

#' Extract response and baseline data without fitting transformations
#'
#' @param formula A two-sided life-history formula, or NULL for vector input.
#' @param data A data frame for formula input.
#' @param x Repeated baseline data for vector input, as a data frame or matrix.
#' @param id Identifier column name for formula input, or an identifier vector.
#' @param response A life_response object for vector input.
#' @param schema An existing schema supplying the formula mapping for new data.
#' @param extra_predictors Additional raw baseline column names required by a
#'   separate censoring formula. These are stored with the predictor schema.
#' @returns A life_input object. Baseline variables remain in their original
#'   types; training-only preprocessing constructs the eventual design matrix.
#' @noRd
resolve_input <- function(formula = NULL, data = NULL, x = NULL, id = NULL,
                          response = NULL, schema = NULL,
                          extra_predictors = character()) {
  if (!is.null(schema)) {
    if (!is.null(formula)) {
      stop("Use the stored formula mapping when a schema is supplied.", call. = FALSE)
    }
    if (is.null(response) && !is.null(schema$formula)) formula <- schema$formula
    if (is.null(id) && !is.null(schema$id_column)) id <- schema$id_column
    extra_predictors <- schema$extra_predictors
  }
  if (!is.null(formula)) {
    if (!inherits(formula, "formula") || length(formula) != 3L) {
      stop("formula must be a two-sided formula with a LifeHist response.",
           call. = FALSE)
    }
    if (!is.data.frame(data) || !nrow(data)) {
      stop("Formula input requires a nonempty data frame.", call. = FALSE)
    }
    if (is.null(names(data)) || anyNA(names(data)) || anyDuplicated(names(data)) ||
        any(!nzchar(names(data)))) {
      stop("Data columns must have unique, nonempty names.", call. = FALSE)
    }
    # Normalize the container before column selection. A data-frame subclass
    # may give one-index `[` row semantics or retain unrequested columns.
    # Column classes and row order remain available to the shared validation.
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    if (!is.null(x) || !is.null(response)) {
      stop("Do not combine formula input with x or a response object.", call. = FALSE)
    }
    response_columns <- all.vars(formula[[2L]])
    missing_response <- setdiff(response_columns, names(data))
    if (length(missing_response)) {
      stop("Missing response column(s): ", paste(missing_response, collapse = ", "),
           ". Use data columns in the formula, or supply response vectors directly.",
           call. = FALSE)
    }
    response <- tryCatch(
      eval(formula[[2L]], envir = data, enclos = environment(formula)),
      error = function(e) stop("Cannot construct the formula response: ",
        conditionMessage(e), " Use a LifeHist() response with the required data columns.",
        call. = FALSE))
    predictors <- unique(c(all.vars(formula[[3L]]), extra_predictors))
    if ("." %in% predictors) {
      stop("Select baseline covariates explicitly; '.' is not supported.", call. = FALSE)
    }
    absent <- setdiff(predictors, names(data))
    if (length(absent)) {
      stop("Missing baseline column(s): ", paste(absent, collapse = ", "), ".",
           call. = FALSE)
    }
    x <- data[predictors]
    if (is.null(id)) id <- "id"
    id_column <- if (is.character(id) && length(id) == 1L) id else NULL
    if (!is.null(id_column)) {
      if (!id_column %in% names(data)) {
        stop("Identifier column not found: ", id_column, ".", call. = FALSE)
      }
      id <- data[[id_column]]
    }
  } else {
    predictors <- NULL
    id_column <- NULL
    if (!is.null(data)) {
      stop("Supply a formula with data, or use response, id and x vectors.",
           call. = FALSE)
    }
  }
  if (!inherits(response, "life_response")) {
    stop("The response must be constructed by LifeHist().", call. = FALSE)
  }
  validate_life_response(response)
  n <- nrow(response$data)
  if (!is.null(data) && nrow(data) != n) {
    stop("The response must retain every row of data.", call. = FALSE)
  }
  if (length(id) != n || !is.atomic(id) || !is.null(dim(id)) ||
      is.complex(id) || !(is.numeric(id) || is.character(id) || is.factor(id))) {
    stop("id must be a numeric, character or factor vector with one value per record.",
         call. = FALSE)
  }
  invalid_id <- is.na(id) | !nzchar(as.character(id))
  if (is.numeric(id)) invalid_id <- invalid_id | !is.finite(id)
  if (any(invalid_id)) {
    stop("Subject identifiers in id must be nonmissing, finite and nonempty",
         " (row ", which(invalid_id)[1L], ").", call. = FALSE)
  }
  if (is.null(x)) x <- data.frame(row.names = seq_len(n))
  if (!is.data.frame(x) && !is.matrix(x)) {
    stop("x must be a data frame or matrix with one row per record.", call. = FALSE)
  }
  if (nrow(x) != n) {
    stop("Baseline data must have one row per record; recycling is not allowed.",
         call. = FALSE)
  }
  if (is.matrix(x) && is.null(colnames(x))) {
    colnames(x) <- paste0("x", seq_len(ncol(x)))
  }
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (anyNA(names(x)) || anyDuplicated(names(x)) || any(!nzchar(names(x)))) {
    stop("Baseline columns must have unique, nonempty names.", call. = FALSE)
  }
  if (!is.null(schema)) {
    if (!setequal(names(x), schema$predictor_names)) {
      stop("Baseline column names must match the stored predictor schema.", call. = FALSE)
    }
    x <- x[schema$predictor_names]
  }
  for (name in names(x)) {
    value <- x[[name]]
    if (!is.atomic(value) || !is.null(dim(value)) || is.complex(value)) {
      stop("Baseline column '", name,
           "' must be a vector without missing or nonfinite values.", call. = FALSE)
    }
    supported <- is.factor(value) || is.character(value) ||
      ((is.numeric(value) || is.logical(value)) && !is.object(value))
    if (!supported) stop("Baseline column '", name,
      "' must be numeric, logical, character, or a factor; convert other classes explicitly.",
      call. = FALSE)
    invalid <- is.na(value)
    if (is.numeric(value)) invalid <- invalid | !is.finite(value)
    if (any(invalid)) {
      row <- which(invalid)[1L]
      history_error(paste0("Baseline column '", name,
        "' must be a vector without missing or nonfinite values"), id[row], row)
    }
  }
  if (!is.null(schema) && !identical(response$format, schema$input_format)) {
    stop("The response format differs from the stored input schema.", call. = FALSE)
  }
  structure(list(response = response, id = id, baseline = x, source_data = data,
                 schema = list(input_format = response$format, formula = formula,
                               predictor_names = names(x), id_column = id_column,
                               extra_predictors = extra_predictors)),
            class = "life_input")
}
