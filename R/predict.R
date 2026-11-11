#' Resolve baseline initial states for prediction
#' @noRd
prediction_initial_states <- function(fit, newdata, initial_state = NULL) {
  n <- nrow(newdata)
  rule <- fit$schema$initial_state_rule
  if (is.null(initial_state)) {
    if (rule$type == "common") initial_state <- rule$value
    else if (rule$type == "column") {
      if (!rule$column %in% names(newdata))
        stop("New data require initial-state column '", rule$column, "'.",
             call. = FALSE)
      initial_state <- newdata[[rule$column]]
    } else stop("Supply initial_state for these new subjects.", call. = FALSE)
  }
  initial_state <- state_labels(initial_state, "initial_state", allow_empty = n == 0L)
  if (length(initial_state) == 1L) initial_state <- rep(initial_state, n)
  if (length(initial_state) != n)
    stop("initial_state must be common or have one value per new subject.",
         call. = FALSE)
  index <- match(as.character(initial_state), as.character(fit$problem$states))
  if (anyNA(index)) stop("Unknown initial state in new data.", call. = FALSE)
  index
}

#' Predict target moments and scales at baseline
#' @noRd
predict_life_moments <- function(fit, newdata, initial_state = NULL,
                                 outcome_view = NULL) {
  if (is.matrix(newdata)) newdata <- preprocessing_data_frame(newdata)
  if (!is.data.frame(newdata))
    stop("newdata must be a baseline data frame, one row per subject.", call. = FALSE)
  initial <- prediction_initial_states(fit, newdata, initial_state)
  view <- if (is.null(outcome_view))
    apply_preprocessor(newdata, fit$preprocessor, "outcome") else outcome_view
  if (!nrow(newdata)) {
    empty <- matrix(numeric(), 0L, ncol(fit$problem$rewards),
                     dimnames = list(NULL, colnames(fit$problem$rewards)))
    return(list(mean = empty, variance = empty, scale = empty,
                initial_state = initial, diagnostics = list()))
  }
  moments <- if (fit$outcome$type == "markov_cox") {
    predict_markov_cox_moments(fit$outcome, view$x, fit$problem$rewards, initial)
  } else flexible_target_moments(fit$outcome, view$x, initial,
                                  B = fit$control$draws)
  if (fit$outcome$type != "markov_cox") {
    selected <- match(colnames(fit$problem$rewards), colnames(moments$mean))
    for (name in c("mean", "variance", "second_moment"))
      moments[[name]] <- moments[[name]][, selected, drop = FALSE]
  }
  moments <- exact_deterministic_moments(moments, fit, initial)
  moments$scale <- life_score_scale(moments$variance, fit)
  if (any(!is.finite(moments$scale)) || any(moments$scale <= 0))
    stop("Score scales overflow or underflow; check the time unit and relative widths.", call. = FALSE)
  moments$initial_state <- initial
  moments
}

#' Preserve centers that are fixed by structural occupation support
#' @noRd
exact_deterministic_moments <- function(moments, fit, initial) {
  support <- fit$problem$structural_support
  lower <- support$reward_min[initial, , drop = FALSE]
  upper <- support$reward_max[initial, , drop = FALSE]
  known <- lower == upper
  if (any(known)) {
    # Direct multiplication avoids accumulating interval-rounding error into
    # a deterministic center, which could change an inclusive zero-score atom.
    exact <- fit$problem$tau * lower[known]
    moments$mean[known] <- exact
    moments$variance[known] <- 0
    moments$second_moment[known] <- exact^2
  }
  moments
}

#' Evaluate either censoring adapter with an explicit side
#' @noRd
predict_life_censor <- function(adapter, x, time, side, subject_index = NULL) {
  if (adapter$type == "cox")
    return(predict_censor_survival(adapter, x, time, side, subject_index))
  predict_censor_forest(adapter, x, time, side, subject_index)
}

#' Predict baseline target moments or calibrated prediction regions
#'
#' @param object A trained or calibrated lifeconformal object.
#' @param newdata Baseline covariate data frame or matrix with one row per new
#'   subject. Required raw column names and factor levels are those used during
#'   training; named columns may be reordered. An unnamed matrix receives names
#'   `x1`, `x2`, and so on, matching training with an unnamed `x` matrix. Use a
#'   data frame for mixed numeric and categorical predictors.
#' @param initial_state Optional common initial-state label or one per subject.
#'   The stored common or column-based rule is used when omitted.
#' @param type For a trained object, return `"moments"` (default), `"score"`
#'   for supplied target values, or `"coverage"` for the working outcome-model
#'   score distribution. Working probabilities are not calibrated coverage.
#' @param values Target values for `type = "score"`, with one row per subject
#'   and named columns matching the trained targets.
#' @param q Ordered score thresholds for `type = "coverage"`.
#' @param ... Reserved for future extensions. Unused arguments are rejected.
#'
#' @returns For a fitted object, `type = "moments"` gives a list of target
#'   `mean`, `variance`, and `scale` matrices; `"score"` gives a numeric score
#'   per subject; `"coverage"` gives a list with a probability matrix, thresholds,
#'   and numerical diagnostics for the working outcome distribution.
#'   Coverage diagnostics identify `approximation` and
#'   `maximum_probability_bound_width`; exact Markov calculations report zero,
#'   bounded grid calculations report their largest bracket width, and flexible
#'   Monte Carlo calculations report NA for deterministic probability bounds.
#'   For a calibrated object, a `lifeconformal_prediction` containing
#'   a `bounds` table, the joint region and its support constraints. Coordinate
#'   intervals enclose projections; joint coverage concerns the entire target
#'   vector in the retained region.
#' @details Predictions use frozen preprocessing and fitted models. They do not
#'   update model coefficients, factor levels, or calibration thresholds.
#'   Occupation-time bounds are intersected with conservative structural bounds
#'   for the initial state. Joint simplex-image support is retained in the
#'   returned object because overlapping targets need not vary independently.
#'   Returned coordinate bounds are conservative numerical envelopes of the
#'   floating-point score region. Directed outward rounding includes every
#'   representable target whose computed score is at most the cutoff, including
#'   division-underflow cases at a zero cutoff. Structural support is then
#'   applied exactly; deterministic coordinates remain singletons.
#' @export
predict.lifeconformal_fit <- function(object, newdata, initial_state = NULL,
                                       type = c("moments", "score", "coverage"),
                                       values = NULL, q = NULL, ...) {
  check_empty_dots(...)
  type <- match.arg(type)
  if (type != "score" && !is.null(values))
    stop("values is used only with type = 'score'.", call. = FALSE)
  if (type != "coverage" && !is.null(q))
    stop("q is used only with type = 'coverage'.", call. = FALSE)
  switch(type, moments = predict_life_moments(object, newdata, initial_state),
         score = evaluate_life_score(object, newdata, values, initial_state),
         coverage = predict_life_coverage(object, newdata, q, initial_state))
}

#' @rdname predict.lifeconformal_fit
#' @export
predict.lifeconformal_calibration <- function(object, newdata,
                                              initial_state = NULL, ...) {
  check_empty_dots(...)
  if (is.matrix(newdata)) newdata <- preprocessing_data_frame(newdata)
  fit <- object$fit
  moments <- predict_life_moments(fit, newdata, initial_state)
  n <- nrow(moments$mean)
  d <- ncol(moments$mean)
  support <- fit$problem$structural_support
  lower_support <- fit$problem$tau *
    support$reward_min[moments$initial_state, , drop = FALSE]
  upper_support <- fit$problem$tau *
    support$reward_max[moments$initial_state, , drop = FALSE]
  # Invert the computed abs((value - mean) / scale) conservatively. If its
  # rounded value is <= finite q, the unrounded quotient is < nextafter(q, +Inf).
  # Round scale*q_up upward, then allow one further representable value for
  # rounding in the original subtraction. Finally round the addition outward.
  # This also covers q=0 when a nonzero quotient underflows to zero. Inf gives
  # full coordinate support. No tolerance in data units or arbitrary epsilon
  # is used, and the fitted score/threshold remain unchanged.
  q_up <- cdt_nextafter(rep(object$threshold, length.out = d))
  width <- sweep(moments$scale, 2L, q_up, `*`)
  width[] <- cdt_nextafter(cdt_nextafter(width))
  lower <- moments$mean - width
  upper <- moments$mean + width
  lower[] <- -cdt_nextafter(-lower)
  upper[] <- cdt_nextafter(upper)
  lower <- pmax(lower_support, lower)
  upper <- pmin(upper_support, upper)
  # R's pmin/pmax preserve the dimensions of the first matrix argument.
  ids <- if ("id" %in% names(newdata)) newdata$id else seq_len(n)
  bounds <- data.frame(id = rep(ids, each = d),
                       target = rep(colnames(moments$mean), times = n),
                       lower = as.vector(t(lower)), upper = as.vector(t(upper)),
                       stringsAsFactors = FALSE)
  structure(list(bounds = bounds, center = moments$mean, scale = moments$scale,
                 threshold = object$threshold, tau = fit$problem$tau,
                 alpha = object$alpha, method = object$method,
                 full_support = object$full_support,
                 component_full_support = object$component_full_support,
                 alpha_components = object$alpha_components,
                 region = list(type = "standardized_max_intersect_support",
                               coordinate_bounds = "conservative_floating_score_envelope",
                               support = support,
                               initial_state = moments$initial_state),
                 numerical = object$diagnostics$numerical),
            class = "lifeconformal_prediction")
}

#' @noRd
check_empty_dots <- function(...) {
  if (length(list(...))) stop("Unused arguments in ... .", call. = FALSE)
}

#' Print and summarize lifeconformal analyses
#' @param x,object A lifeconformal analysis object.
#' @param ... Reserved for future extensions.
#' @returns Print methods invisibly return their input. Summary methods return
#'   a compact list of model, target, calibration, and diagnostic information.
#' @export
print.lifeconformal_fit <- function(x, ...) {
  cat("Life-history training fit\n")
  cat(length(x$training_ids), "subjects;", length(x$problem$states), "states;",
      ncol(x$problem$rewards), "targets; horizon", x$problem$tau, "\n")
  cat("Outcome:", x$outcome$type, " | Censoring:", x$censor$type, "\n")
  invisible(x)
}

#' @rdname print.lifeconformal_fit
#' @export
summary.lifeconformal_fit <- function(object, ...) {
  adapters <- list(outcome = object$outcome, censoring = object$censor,
                   augmentation = object$augmentation)
  list(n = length(object$training_ids), states = object$problem$states,
       initial_state = object$schema$initial_state_rule,
       absorbing = object$problem$absorbing_labels,
       targets = object$problem$targets$sets, tau = object$problem$tau,
       outcome = object$outcome$type, censoring = object$censor$type,
       augmentation = object$augmentation$type,
       augmentation_reused = object$augmentation_reused,
       model_options = lapply(adapters, function(adapter)
         if (!is.null(adapter$options)) adapter$options else adapter$controls),
       model_diagnostics = lapply(adapters, `[[`, "diagnostics"),
       control = object$control, diagnostics = object$diagnostics)
}

#' @rdname print.lifeconformal_fit
#' @export
print.lifeconformal_calibration <- function(x, ...) {
  cat("Calibrated life-history prediction regions\n")
  cat(x$n, "calibration subjects;", x$method, "; nominal coverage",
      format(1 - x$alpha), "\n")
  if (!is.null(x$alpha_components)) {
    print(data.frame(target = names(x$threshold), alpha = unname(x$alpha_components),
      threshold = unname(x$threshold), cutoff = unname(x$cutoff),
      full_support = unname(x$component_full_support),
      reason = unname(x$component_reason)),
      row.names = FALSE)
    cat("Allocated total miscoverage:", format(sum(x$alpha_components)), "\n")
  } else {
    cat("Threshold:", format(x$threshold),
        if (x$full_support) "(full target support)" else "", "\n")
    cat("Rank-adjusted cutoff:", format(x$cutoff), "\n")
    if (x$full_support) cat("Reason:",
      if (identical(x$reason, "rank_exceeds_sample"))
        "the requested rank exceeds the calibration sample size."
      else "no candidate threshold reached the cutoff.", "\n")
  }
  width <- x$diagnostics$numerical$max_curve_bound_width
  if (!is.null(width) && is.finite(width) && width > 0) {
    cat("Numerical coverage-curve bound width:", format(width), "\n")
    stable <- x$diagnostics$numerical$threshold_stable
    if (!is.null(stable) && !is.na(stable))
      cat("Threshold accuracy:", if (stable) "resolved" else "unresolved", "\n")
    met <- x$diagnostics$numerical$probability_tolerance_met
    if (identical(met, FALSE))
      cat("Conditional probability tolerance: not met; inspect numerical bounds.\n")
  }
  invisible(x)
}

#' @rdname print.lifeconformal_fit
#' @export
summary.lifeconformal_calibration <- function(object, ...) {
  list(n = object$n, alpha = object$alpha, method = object$method,
       alpha_components = object$alpha_components,
       allocated_alpha = object$allocated_alpha,
       threshold = object$threshold, full_support = object$full_support,
       component_full_support = object$component_full_support,
       any_component_full_support = object$any_component_full_support,
       component_reason = object$component_reason,
       rank = object$rank, reason = object$reason,
       cutoff = object$cutoff, curve = object$curve,
       diagnostics = object$diagnostics)
}

#' @rdname print.lifeconformal_fit
#' @export
print.lifeconformal_prediction <- function(x, ...) {
  cat("Life-history prediction bounds; nominal joint coverage",
      format(1 - x$alpha), "\n")
  if (isTRUE(x$full_support)) cat("Full target support returned.\n")
  else if (any(x$component_full_support))
    cat("Full support for:", paste(names(x$component_full_support)[x$component_full_support],
                                    collapse = ", "), "\n")
  if (ncol(x$center) > 1L)
    cat("Coordinate bounds enclose a joint region with occupation-time constraints.\n")
  print(x$bounds, row.names = FALSE)
  invisible(x)
}

#' @rdname print.lifeconformal_fit
#' @export
summary.lifeconformal_prediction <- function(object, ...) {
  list(n = nrow(object$center), targets = colnames(object$center),
       tau = object$tau, alpha = object$alpha, method = object$method,
       threshold = object$threshold, alpha_components = object$alpha_components,
       full_support = object$full_support,
       component_full_support = object$component_full_support,
       coordinate_bounds = object$region$coordinate_bounds, bounds = object$bounds,
       joint_support = object$region$support$joint_support)
}

#' Extract prediction bounds as a data frame
#' @param x A lifeconformal prediction object.
#' @param row.names Optional row names for the returned table.
#' @param optional Logical, passed to [as.data.frame()].
#' @param ... Additional arguments to the data-frame method.
#' @returns A data frame with subject identifiers, target names, lower and
#'   upper bounds. The full prediction object retains the joint constraints.
#' @export
as.data.frame.lifeconformal_prediction <- function(x, row.names = NULL,
                                                    optional = FALSE, ...) {
  as.data.frame(x$bounds, row.names = row.names, optional = optional, ...)
}

#' Plot prediction intervals by subject and target
#' @param x A lifeconformal prediction object.
#' @param ... Named arguments to [graphics::plot.default()], including `xlab`,
#'   `ylab`, `main`, axis limits and graphical parameters. `col`, `lty` and
#'   `lwd` also control interval segments. Intervals are always drawn; `type`
#'   optionally adds points or lines at their lower endpoints.
#' @returns The prediction object, invisibly.
#' @details The default vertical axis labels each subject and target. Set
#'   `axes = FALSE` or `yaxt = "n"` to suppress these labels, or `yaxt = "s"`
#'   for the usual numeric axis. Joint predictions retain their support
#'   constraints; displayed coordinate intervals can include combinations that
#'   are outside that joint support. Singleton and zero-cutoff intervals receive
#'   a point marker. Default horizontal labels are scaled to fit the current
#'   left margin; `cex.axis` overrides this adjustment. Empty predictions have
#'   no rows to plot.
#' @export
plot.lifeconformal_prediction <- function(x, ...) {
  bounds <- x$bounds
  if (!nrow(bounds)) stop("There are no prediction rows to plot.", call. = FALSE)
  y <- rev(seq_len(nrow(bounds)))
  dots <- list(...)
  args <- life_plot_arguments(list(x = bounds$lower, y = y,
    xlim = range(c(bounds$lower, bounds$upper)), ylim = range(c(0.5, y + 0.5)),
    type = "n", xlab = "Occupation time", ylab = "", yaxt = "n",
    col = graphics::par("fg"), lty = 1, lwd = 2), dots)
  do.call(graphics::plot.default, args)
  graphics::segments(bounds$lower, y, bounds$upper, y,
                      col = args$col, lty = args$lty, lwd = args$lwd)
  singleton <- bounds$lower == bounds$upper
  if (!is.null(x$threshold)) singleton <- singleton |
    rep(rep(x$threshold, length.out = ncol(x$center)) == 0, times = nrow(x$center))
  if (any(singleton)) graphics::points(
    (bounds$lower[singleton] / 2 + bounds$upper[singleton] / 2), y[singleton],
    pch = 16, col = rep(args$col, length.out = nrow(bounds))[singleton],
    cex = if (is.null(args$cex)) 1 else rep(args$cex, length.out = nrow(bounds))[singleton])
  if (is.null(dots$yaxt) && !identical(dots$axes, FALSE)) {
    axis_args <- list(side = 2, at = y,
      labels = paste(bounds$id, bounds$target, sep = ": "),
      las = if (is.null(dots$las)) 1 else dots$las)
    keep <- intersect(names(dots), c("cex.axis", "col.axis", "font.axis"))
    axis_args[keep] <- dots[keep]
    if (is.null(dots$cex.axis) && axis_args$las %in% c(1, 2)) {
      size <- graphics::par("cex.axis")
      available <- graphics::par("mai")[2] -
        (graphics::par("mgp")[2] + .2) * graphics::par("csi") * graphics::par("mex")
      widest <- max(graphics::strwidth(axis_args$labels, units = "inches", cex = size,
        font = if (is.null(axis_args$font.axis)) graphics::par("font.axis") else axis_args$font.axis))
      if (available > 0 && widest > available)
        axis_args$cex.axis <- size * available / widest
    }
    do.call(graphics::axis, axis_args)
  }
  invisible(x)
}

#' Plot a calibration coverage curve
#' @param x A lifeconformal calibration object.
#' @param ... Named arguments to [graphics::plot.default()], including labels,
#'   title, axis limits and graphical parameters. For component calibration,
#'   `col`, `lty`, `lwd`, `pch` and `cex` are recycled in target order.
#' @param legend Logical; show the target legend for component calibration.
#' @returns The input object, invisibly.
#' @details The raw curve is shown, including nonmonotone values and values
#'   outside `[0, 1]`. The dashed line is the actual rank-adjusted cutoff.
#'   Default vertical limits include the cutoff. Components without a curve
#'   because their rank exceeds the sample size are labeled as full support.
#'   If no curve was evaluated, the plot explains the full-support fallback.
#' @export
plot.lifeconformal_calibration <- function(x, ..., legend = TRUE) {
  if (!is.logical(legend) || length(legend) != 1L || is.na(legend))
    stop("legend must be TRUE or FALSE.", call. = FALSE)
  dots <- list(...)
  if (!nrow(x$curve)) {
    args <- life_plot_arguments(list(x = c(0, 1), y = c(0, 1), type = "n",
      axes = FALSE, xlab = "", ylab = "", main = "Full target support"), dots)
    do.call(graphics::plot.default, args)
    usr <- graphics::par("usr")
    graphics::text(mean(usr[1:2]), mean(usr[3:4]),
      paste(strwrap("The requested rank exceeds the calibration sample size.", width = 45),
            collapse = "\n"), col = if (is.null(args$col)) graphics::par("fg") else args$col,
      cex = if (is.null(args$cex)) .9 else args$cex)
    return(invisible(x))
  }
  if (!is.null(x$alpha_components)) {
    target <- names(x$alpha_components)
    args <- life_plot_arguments(list(x = range(x$curve$q),
      y = range(c(x$curve$estimate, x$cutoff)), type = "l", xlab = "Threshold",
      ylab = "Estimated coverage", col = seq_along(target), lty = 1, lwd = 1), dots)
    canvas <- args
    canvas$type <- "n"
    do.call(graphics::plot.default, canvas)
    colors <- rep(args$col, length.out = length(target))
    widths <- rep(args$lwd, length.out = length(target))
    types <- rep(args$lty, length.out = length(target))
    has_curve <- logical(length(target))
    for (j in seq_along(target)) {
      rows <- x$curve$target == target[j]
      has_curve[j] <- any(rows)
      if (has_curve[j]) {
        line_args <- list(x = x$curve$q[rows], y = x$curve$estimate[rows],
                          type = args$type, col = colors[j], lty = types[j], lwd = widths[j])
        for (name in intersect(names(args), c("pch", "cex", "bg")))
          line_args[[name]] <- rep(args[[name]], length.out = length(target))[j]
        do.call(graphics::lines, line_args)
      }
      graphics::abline(h = x$cutoff[j], lty = 2, col = colors[j], lwd = widths[j])
    }
    if (legend) {
      labels <- paste0(target, ifelse(x$component_full_support, " (full support)", ""))
      point_type <- args$type %in% c("p", "b", "o")
      pch <- if (point_type) rep(if (is.null(args$pch)) 1 else args$pch,
                                length.out = length(target)) else rep(NA, length(target))
      pch[!has_curve] <- NA
      graphics::legend("bottomright", legend = labels, col = colors,
        lty = ifelse(has_curve & args$type == "p", NA, ifelse(has_curve, types, 2)),
        lwd = widths, pch = pch)
    }
    return(invisible(x))
  }
  args <- life_plot_arguments(list(x = x$curve$q, y = x$curve$estimate,
    ylim = range(c(x$curve$estimate, x$cutoff)), type = "l", xlab = "Threshold",
    ylab = "Estimated coverage", col = graphics::par("fg"), lwd = 1), dots)
  do.call(graphics::plot.default, args)
  graphics::abline(h = x$cutoff, lty = 2, col = args$col, lwd = args$lwd)
  invisible(x)
}

#' Merge graphical defaults without duplicate actual arguments
#' @noRd
life_plot_arguments <- function(defaults, dots) {
  if (length(dots) && (is.null(names(dots)) || anyNA(names(dots)) ||
      any(!nzchar(names(dots))) || anyDuplicated(names(dots))))
    stop("Plot arguments in ... must have unique, nonempty names.", call. = FALSE)
  if (any(names(dots) %in% c("x", "y")))
    stop("Plot coordinates are supplied by the lifeconformal object.", call. = FALSE)
  defaults[names(dots)] <- dots
  defaults
}
