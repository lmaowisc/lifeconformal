#' Fit subject-level predictor mappings on training data
#'
#' This boundary receives one verified baseline row per subject. Formula
#' environments are discarded: predictors must be columns and constants must
#' be literals. Supported functions are deterministic, pointwise base-R
#' transformations and `poly()`, `splines::ns()`, `splines::bs()`, and `scale()`.
#' Polynomial coefficients, spline knots and scaling constants are fitted here,
#' including when nested in another supported expression, and never recomputed
#' on calibration or prediction rows. Other sample-dependent transformations
#' must be prepared externally using training-fitted constants.
#'
#' @param baseline Subject-level data frame or matrix. An unnamed matrix receives
#'   names x1, x2, and so on. Character predictors are categorical; their training
#'   levels become factors. Other supported types are numeric, logical, and
#'   factors. Required values must be finite and nonmissing.
#' @param formula_spec One- or two-sided formula. Only its right-hand side is
#'   used. NULL includes the explicitly supplied baseline columns, or uses an
#'   intercept-only specification when there are none. Dot expansion is not
#'   supported.
#' @param censor_formula Optional one-sided formula for censoring predictors;
#'   NULL reuses the outcome mapping without refitting it.
#' @param model_options Preprocessing options only. The optional `contrasts`
#'   entry is a named list of numeric contrast matrices or standard stats
#'   contrast names. Matrices are frozen even when the global contrast option
#'   later changes. Backend fitting options belong to the model adapters.
#' @returns A list with `preprocessor` (class life_preprocessor), `baseline`
#'   (typed union of required columns), and `outcome`/`censor` predictor views.
#'   Each view has `baseline` and numeric matrix `x`. The conventional intercept
#'   column is removed from `x` because the survival models fit their own baseline
#'   law; intercept-only inputs therefore produce an n-by-zero matrix. Formula
#'   intercept settings are retained in the terms and factor coding. The frozen
#'   object stores no training rows or caller environment.
#' @noRd
fit_preprocessor <- function(baseline, formula_spec = NULL,
                             censor_formula = NULL, model_options = list()) {
  baseline <- preprocessing_data_frame(baseline)
  if (!nrow(baseline)) {
    stop("Preprocessing requires at least one training subject.", call. = FALSE)
  }
  if (!is.list(model_options) ||
      (length(model_options) && (is.null(names(model_options)) ||
       anyNA(names(model_options)) || anyDuplicated(names(model_options)) ||
       any(!nzchar(names(model_options))))) ||
      any(!names(model_options) %in% "contrasts")) {
    stop("Preprocessing model_options must be a named list containing only contrasts.",
         call. = FALSE)
  }
  contrast_options <- model_options$contrasts
  if (!is.null(contrast_options) &&
      (!is.list(contrast_options) || (length(contrast_options) &&
       is.null(names(contrast_options))) ||
       anyNA(names(contrast_options)) || anyDuplicated(names(contrast_options)) ||
       any(!nzchar(names(contrast_options))))) {
    stop("contrasts must be a named list with unique, nonempty names.", call. = FALSE)
  }
  outcome_formula <- preprocessing_formula(formula_spec, names(baseline))
  if (!is.null(censor_formula) &&
      (!inherits(censor_formula, "formula") || length(censor_formula) != 2L)) {
    stop("censor_formula must be a one-sided formula.", call. = FALSE)
  }
  censor_formula <- if (is.null(censor_formula)) outcome_formula else
    preprocessing_formula(censor_formula, names(baseline))
  required <- unique(c(all.vars(outcome_formula), all.vars(censor_formula)))
  schema <- fit_baseline_schema(baseline, required)
  baseline <- apply_baseline_schema(baseline, schema)
  outcome <- fit_predictor_mapping(baseline, outcome_formula, contrast_options)
  censor <- if (identical(outcome_formula, censor_formula)) outcome else
    fit_predictor_mapping(baseline, censor_formula, contrast_options)
  used_contrasts <- union(names(outcome$mapping$contrasts),
                          names(censor$mapping$contrasts))
  if (length(setdiff(names(contrast_options), used_contrasts))) {
    stop("Contrast options refer to unknown or noncategorical term(s): ",
         paste(setdiff(names(contrast_options), used_contrasts), collapse = ", "),
         ".", call. = FALSE)
  }
  preprocessor <- structure(
    list(version = 1L, required_variables = required, baseline_schema = schema,
         outcome = outcome$mapping, censor = censor$mapping),
    class = "life_preprocessor")
  list(preprocessor = preprocessor, baseline = baseline,
       outcome = outcome$view, censor = censor$view)
}

#' Apply a frozen mapping without learning from the new sample
#' @param baseline One baseline row per new subject, with named columns.
#' @param preprocessor The life_preprocessor returned by fit_preprocessor.
#' @param purpose Select the outcome or censoring mapping. Only the variables
#'   required for that purpose need to be supplied.
#' @returns A life_predictor_view containing typed `baseline` and numeric `x`.
#'   Zero-row and one-row predictions retain the training column dimensions.
#' @noRd
apply_preprocessor <- function(baseline, preprocessor,
                               purpose = c("outcome", "censor")) {
  purpose <- match.arg(purpose)
  if (!inherits(preprocessor, "life_preprocessor") ||
      !identical(preprocessor$version, 1L)) {
    stop("preprocessor must be a supported life_preprocessor object.", call. = FALSE)
  }
  mapping <- preprocessor[[purpose]]
  baseline <- apply_baseline_schema(
    preprocessing_data_frame(baseline),
    preprocessor$baseline_schema[mapping$required_variables])
  if (!nrow(baseline)) {
    x <- matrix(numeric(), 0L, length(mapping$columns),
                dimnames = list(rownames(baseline), mapping$columns))
    attr(x, "assign") <- mapping$assign
    if (length(mapping$contrasts)) attr(x, "contrasts") <- mapping$contrasts
  } else {
    frame <- tryCatch(stats::model.frame(mapping$terms, data = baseline,
                                na.action = stats::na.fail,
                                xlev = mapping$xlevels, drop.unused.levels = FALSE),
      error = function(e) stop("Cannot evaluate baseline formula ",
        paste(deparse(stats::formula(mapping$terms)), collapse = " "), ": ",
        conditionMessage(e), ". Check the transformation inputs and domains.", call. = FALSE))
    if (nrow(frame) != nrow(baseline)) {
      stop("Formula transformations must retain one row per subject; use ~1 for no covariates.",
           call. = FALSE)
    }
    x <- preprocessing_matrix(mapping$terms, frame, mapping$contrasts)
    if (!identical(colnames(x), mapping$columns)) {
      stop("The frozen predictor mapping produced incompatible design columns.",
           call. = FALSE)
    }
  }
  structure(list(baseline = baseline, x = x), class = "life_predictor_view")
}

#' Convert an explicit predictor container without changing its row count
#' @noRd
preprocessing_data_frame <- function(baseline) {
  if (!is.data.frame(baseline) && !is.matrix(baseline)) {
    stop("baseline must be a subject-level data frame or matrix.", call. = FALSE)
  }
  if (is.matrix(baseline) && ncol(baseline) && is.null(colnames(baseline))) {
    colnames(baseline) <- paste0("x", seq_len(ncol(baseline)))
  }
  baseline <- as.data.frame(baseline, stringsAsFactors = FALSE)
  if (anyNA(names(baseline)) || anyDuplicated(names(baseline)) ||
      any(!nzchar(names(baseline)))) {
    stop("Baseline columns must have unique, nonempty names.", call. = FALSE)
  }
  baseline
}

#' Remove the response and all dependencies on the caller's environment
#' @noRd
preprocessing_formula <- function(formula_spec, column_names) {
  if (is.null(formula_spec)) {
    rhs <- if (!length(column_names)) 1 else
      Reduce(function(a, b) call("+", a, b), lapply(column_names, as.name))
  } else {
    if (!inherits(formula_spec, "formula") ||
        !length(formula_spec) %in% c(2L, 3L)) {
      stop("formula_spec must be a one- or two-sided formula.", call. = FALSE)
    }
    rhs <- formula_spec[[length(formula_spec)]]
  }
  if ("." %in% all.vars(rhs)) {
    stop("Select baseline covariates explicitly; '.' is not supported.", call. = FALSE)
  }
  rhs <- validate_predictor_expression(rhs)
  stats::as.formula(call("~", rhs), env = baseenv())
}

#' Validate and canonicalize the supported, serializable formula language
#' @noRd
validate_predictor_expression <- function(expr) {
  if (!is.call(expr)) {
    return(expr)
  }
  head <- expr[[1L]]
  namespace <- NULL
  if (is.call(head) && identical(head[[1L]], quote(`::`)) &&
      length(head) == 3L) {
    namespace <- as.character(head[[2L]])
    fun <- as.character(head[[3L]])
  } else if (is.symbol(head)) {
    fun <- as.character(head)
  } else {
    stop("Formula functions must be explicitly named supported transformations.",
         call. = FALSE)
  }
  base_functions <- c("+", "-", "*", "/", ":", "^", "(", "I", "c", "log",
                      "log1p", "log2", "log10", "exp", "expm1", "sqrt", "abs",
                      "sign", "sin", "cos", "tan", "asin", "acos", "atan",
                      "ceiling", "floor", "trunc", "round", "signif", "pmin",
                      "pmax", "as.numeric", "as.double", "as.integer", "scale",
                      "==", "!=", ">", ">=", "<", "<=", "!", "&", "|", "%in%")
  expected <- if (fun %in% base_functions) "base" else if (fun == "poly") "stats" else
    if (fun %in% c("ns", "bs")) "splines" else NA_character_
  if (is.na(expected) || (!is.null(namespace) && namespace != expected)) {
    stop("Unsupported formula transformation '", fun,
         "'. Use supported pointwise transforms, poly(), splines::ns(), ",
         "splines::bs(), or scale(); fit other transformations on training data externally.",
         call. = FALSE)
  }
  if (expected != "base") expr[[1L]] <- call("::", as.name(expected), as.name(fun))
  if (length(expr) > 1L) {
    for (i in seq.int(2L, length(expr))) {
      expr[i] <- list(validate_predictor_expression(expr[[i]]))
    }
  }
  expr
}

#' Record raw predictor types and training factor levels
#' @noRd
fit_baseline_schema <- function(baseline, required) {
  absent <- setdiff(required, names(baseline))
  if (length(absent)) {
    stop("Missing baseline column(s): ", paste(absent, collapse = ", "),
         ". Formula constants must be literals, not variables from an external environment.",
         call. = FALSE)
  }
  result <- lapply(required, function(name) {
    value <- baseline[[name]]
    if (is.factor(value)) {
      list(type = "factor", levels = levels(value), ordered = is.ordered(value),
           contrasts = attr(value, "contrasts"))
    } else if (is.character(value)) {
      list(type = "factor", levels = sort(unique(value)), ordered = FALSE,
           contrasts = NULL)
    } else if (is.logical(value) && !is.object(value)) {
      list(type = "logical")
    } else if (is.numeric(value) && !is.complex(value) && !is.object(value)) {
      list(type = "numeric")
    } else {
      stop("Baseline column '", name,
           "' must be numeric, logical, character, or a factor; convert other classes explicitly.",
           call. = FALSE)
    }
  })
  names(result) <- required
  result
}

#' Validate and align raw predictors using their training type schema
#' @noRd
apply_baseline_schema <- function(baseline, schema) {
  absent <- setdiff(names(schema), names(baseline))
  if (length(absent)) {
    stop("Missing baseline column(s): ", paste(absent, collapse = ", "), ".",
         call. = FALSE)
  }
  baseline <- baseline[names(schema)]
  for (name in names(schema)) {
    value <- baseline[[name]]
    rule <- schema[[name]]
    if (!is.atomic(value) || !is.null(dim(value))) {
      stop("Baseline column '", name, "' must contain finite, nonmissing values.",
           call. = FALSE)
    }
    invalid <- is.na(value)
    if (is.numeric(value)) invalid <- invalid | !is.finite(value)
    if (any(invalid)) stop("Baseline column '", name,
      "' must contain finite, nonmissing values (baseline row ",
      which(invalid)[1L], ").", call. = FALSE)
    if (rule$type == "factor") {
      if (!is.character(value) && !is.factor(value)) {
        stop("Baseline column '", name, "' must be character or factor.", call. = FALSE)
      }
      value <- as.character(value)
      novel <- unique(value[!value %in% rule$levels])
      if (length(novel)) {
        stop("Baseline column '", name, "' has unseen level(s): ",
             paste(utils::head(novel, 5L), collapse = ", "),
             if (length(novel) > 5L) ", ..." else "",
             " (baseline row ", which(!value %in% rule$levels)[1L],
             "); use levels seen in training: ",
             paste(utils::head(rule$levels, 5L), collapse = ", "),
             if (length(rule$levels) > 5L) ", ..." else "", ".", call. = FALSE)
      }
      value <- factor(value, levels = rule$levels, ordered = rule$ordered)
      if (!is.null(rule$contrasts)) attr(value, "contrasts") <- rule$contrasts
    } else if (rule$type == "logical") {
      if (!is.logical(value) || is.object(value)) {
        stop("Baseline column '", name, "' must be logical.", call. = FALSE)
      }
    } else if (!is.numeric(value) || is.complex(value) || is.object(value)) {
      stop("Baseline column '", name, "' must be numeric.", call. = FALSE)
    }
    baseline[[name]] <- value
  }
  baseline
}

#' Freeze training-dependent transformations recursively in prediction calls
#' @noRd
freeze_predictor_expression <- function(expr, baseline) {
  if (!is.call(expr)) return(expr)
  if (length(expr) > 1L) {
    for (i in seq.int(2L, length(expr))) {
      expr[i] <- list(freeze_predictor_expression(expr[[i]], baseline))
    }
  }
  head <- expr[[1L]]
  fun <- if (is.call(head)) as.character(head[[3L]]) else as.character(head)
  if (fun %in% c("poly", "ns", "bs", "scale")) {
    definition <- switch(fun, poly = stats::poly, ns = splines::ns,
                          bs = splines::bs, scale = base::scale.default)
    # makepredictcall.ns/bs retain the first data argument, so normalize named
    # arguments before freezing even when the user writes ns(df = 3, x = age).
    expr <- match.call(definition, expr)
    # These functions learn only from their first data argument. A parameter
    # depending on other rows could otherwise change the map on a new sample.
    args <- as.list(expr)[-1L]
    argument_names <- names(args)
    if (is.null(argument_names)) argument_names <- rep("", length(args))
    data_argument <- if ("x" %in% argument_names) match("x", argument_names) else 1L
    parameters <- args[-data_argument]
    if (any(vapply(parameters, function(z) length(all.vars(z)) > 0L, logical(1)))) {
      stop("Parameters of ", fun, "() must be literal constants; ",
           "use separate univariate polynomial terms for interactions.", call. = FALSE)
    }
  }
  if (fun %in% c("poly", "ns", "bs")) {
    fitting_call <- expr
    # poly(simple = TRUE) deliberately discards its coefficients. Retain them
    # for makepredictcall while preserving the caller's requested output shape.
    if (fun == "poly") fitting_call$simple <- FALSE
    value <- eval(fitting_call, baseline, baseenv())
    expr <- stats::makepredictcall(value, expr)
  } else if (fun == "scale") {
    value <- eval(expr, baseline, baseenv())
    args <- expr
    center <- attr(value, "scaled:center")
    scale <- attr(value, "scaled:scale")
    expr <- as.call(list(quote(base::scale), args$x,
                        center = if (is.null(center)) FALSE else center,
                        scale = if (is.null(scale)) FALSE else scale))
  }
  expr
}

#' Fit one formula mapping and create its initial predictor view
#' @noRd
fit_predictor_mapping <- function(baseline, formula, contrast_options) {
  required <- all.vars(formula)
  baseline <- baseline[required]
  terms <- stats::terms(formula)
  predvars <- attr(terms, "variables")
  if (length(predvars) > 1L) {
    for (i in seq.int(2L, length(predvars))) {
      predvars[[i]] <- freeze_predictor_expression(predvars[[i]], baseline)
    }
  }
  attr(terms, "predvars") <- predvars
  frame <- tryCatch(stats::model.frame(terms, data = baseline, na.action = stats::na.fail,
                              drop.unused.levels = FALSE),
    error = function(e) stop("Cannot evaluate baseline formula ",
      paste(deparse(stats::formula(terms)), collapse = " "), ": ",
      conditionMessage(e), ". Check the transformation inputs and domains.", call. = FALSE))
  if (nrow(frame) != nrow(baseline)) {
    stop("Formula transformations must retain one row per subject; use ~1 for no covariates.",
         call. = FALSE)
  }
  terms <- attr(frame, "terms")
  environment(terms) <- baseenv()
  categorical <- vapply(frame, function(z) is.factor(z) || is.logical(z), logical(1))
  contrasts <- lapply(names(frame)[categorical], function(name) {
    frozen_contrast_matrix(frame[[name]], contrast_options[[name]], name)
  })
  names(contrasts) <- names(frame)[categorical]
  xlevels <- lapply(frame[vapply(frame, is.factor, logical(1))], levels)
  x <- preprocessing_matrix(terms, frame, contrasts)
  mapping <- list(terms = terms, required_variables = required,
                  xlevels = xlevels, contrasts = contrasts,
                  columns = colnames(x), assign = attr(x, "assign"))
  list(mapping = mapping,
       view = structure(list(baseline = baseline, x = x), class = "life_predictor_view"))
}

#' Resolve factor coding to a matrix independent of search paths and options
#' @noRd
frozen_contrast_matrix <- function(value, option, name) {
  if (is.logical(value)) value <- factor(value, levels = c(FALSE, TRUE))
  levels <- levels(value)
  if (length(levels) < 2L) {
    stop("Categorical predictor '", name,
         "' needs at least two declared levels; remove a constant predictor.",
         call. = FALSE)
  }
  if (is.null(option)) option <- attr(value, "contrasts")
  if (is.null(option)) option <- getOption("contrasts")[[1L + is.ordered(value)]]
  if (is.character(option) && length(option) == 1L &&
      option %in% c("contr.treatment", "contr.sum", "contr.helmert", "contr.poly", "contr.SAS")) {
    option <- getExportedValue("stats", option)(levels, contrasts = TRUE)
  }
  if (!is.matrix(option) || !is.numeric(option) || any(!is.finite(option)) ||
      nrow(option) != length(levels) || ncol(option) < 1L ||
      ncol(option) >= length(levels)) {
    stop("Contrasts for '", name,
         "' must be a finite matrix with one row per level and fewer columns, ",
         "or a standard stats contrast name.", call. = FALSE)
  }
  rownames(option) <- levels
  option
}

#' Construct a finite numeric matrix with an explicitly absent Cox intercept
#' @noRd
preprocessing_matrix <- function(terms, frame, contrasts) {
  x <- stats::model.matrix(terms, frame,
                           contrasts.arg = if (length(contrasts)) contrasts else NULL)
  keep <- colnames(x) != "(Intercept)"
  assign <- attr(x, "assign")[keep]
  x <- x[, keep, drop = FALSE]
  storage.mode(x) <- "double"
  attr(x, "assign") <- assign
  if (length(contrasts)) attr(x, "contrasts") <- contrasts
  if (any(!is.finite(x))) {
    invalid <- which(!is.finite(x), arr.ind = TRUE)[1L, ]
    stop("Formula transformations must produce finite design columns; column '",
         colnames(x)[invalid[2L]], "' is nonfinite at baseline row ", invalid[1L],
         ". Check the transformation inputs and domains.", call. = FALSE)
  }
  if (nrow(x) != nrow(frame) || anyDuplicated(colnames(x))) {
    stop("Formula transformations must produce finite, uniquely named design columns without dropping rows.",
         call. = FALSE)
  }
  x
}
