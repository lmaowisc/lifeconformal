#' Evaluate code with an optional local random seed
#' @noRd
flex_with_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  if (length(seed) != 1L || !is.numeric(seed) || !is.null(dim(seed)) || !is.finite(seed) ||
      seed < 1 || seed > .Machine$integer.max || seed != floor(seed))
    stop("seed must be an integer between 1 and .Machine$integer.max.", call. = FALSE)
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(as.integer(seed))
  force(code)
}

#' Validate a flexible-model option list
#' @noRd
flex_options <- function(options, defaults) {
  if (!is.list(options) || (length(options) &&
      (is.null(names(options)) || anyNA(names(options)) || any(!nzchar(names(options))) ||
       anyDuplicated(names(options))))) stop("Model options must be a list with unique, nonempty names.", call. = FALSE)
  unknown <- setdiff(names(options), names(defaults))
  if (length(unknown)) stop("Unknown flexible-model option(s): ",
                            paste(unknown, collapse = ", "), ".", call. = FALSE)
  utils::modifyList(defaults, options, keep.null = TRUE)
}

#' Validate a numeric scalar control
#' @noRd
flex_positive <- function(x, name, integer = FALSE, zero = FALSE) {
  if (length(x) != 1L || !is.numeric(x) || is.complex(x) ||
      !is.null(dim(x)) || !is.finite(x) ||
      if (zero) x < 0 else x <= 0)
    stop(name, " must be a finite ", if (zero) "nonnegative" else "positive",
         " number.", call. = FALSE)
  if (integer && (x != floor(x) || x > .Machine$integer.max))
    stop(name, " must be an integer.", call. = FALSE)
  invisible(x)
}

#' Validate the explicit baseline design for a flexible adapter
#' @noRd
flex_x <- function(x, n = nrow(x), names = NULL) {
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) || nrow(x) != n ||
      any(!is.finite(x))) stop("x must be a finite numeric matrix with one row per subject.", call. = FALSE)
  if (!is.null(names) && !identical(colnames(x), names))
    stop("Baseline design columns differ from training.", call. = FALSE)
  x
}

#' Fit a random survival forest for censoring
#'
#' The response is one subject-level time/event pair; transition records are
#' never treated as independent censoring observations. Every censoring-event
#' knot is retained. Predictions use the backend survival function itself and its
#' left/right values, with a single coherent positive floor.
#' @noRd
fit_censor_forest <- function(history, x, tau, options = list(), seed = NULL) {
  x <- flex_x(x, nrow(history$subjects))
  flex_positive(tau, "tau")
  o <- flex_options(options, list(num_trees = 200L, mtry = NULL,
    min_node_size = 10L, max_depth = 0L, num_threads = 1L, survival_floor = 1e-8))
  for (nm in c("num_trees", "min_node_size", "num_threads"))
    flex_positive(o[[nm]], nm, integer = TRUE)
  flex_positive(o$max_depth, "max_depth", integer = TRUE, zero = TRUE)
  flex_positive(o$survival_floor, "survival_floor")
  if (o$survival_floor > 1) stop("survival_floor cannot exceed one.", call. = FALSE)
  risk <- make_censoring_data(history, x, tau)
  end <- risk$data$time
  event <- risk$data$event
  predictor_count <- max(1L, ncol(x))
  if (is.null(o$mtry)) o$mtry <- max(1L, floor(sqrt(predictor_count)))
  flex_positive(o$mtry, "mtry", integer = TRUE)
  if (o$mtry > predictor_count) stop("mtry exceeds the predictor count.", call. = FALSE)
  stored_seed <- flex_with_seed(seed, sample.int(.Machine$integer.max, 1L))
  model <- NULL
  if (any(event)) {
    if (!requireNamespace("ranger", quietly = TRUE))
      stop("Install the optional 'ranger' package for censor_model = 'survival_forest'.", call. = FALSE)
    xx <- risk$x
    if (!ncol(xx)) xx <- matrix(0, nrow(xx), 1L, dimnames = list(NULL, ".constant"))
    if (is.null(colnames(xx))) colnames(xx) <- paste0("x", seq_len(ncol(xx)))
    model <- flex_with_seed(stored_seed, ranger::ranger(
      x = as.data.frame(xx), y = survival::Surv(end, event),
      num.trees = o$num_trees, mtry = o$mtry, min.node.size = o$min_node_size,
      max.depth = o$max_depth, num.threads = o$num_threads,
      splitrule = "logrank", importance = "none", oob.error = FALSE,
      write.forest = TRUE, seed = stored_seed, verbose = FALSE))
    # The backend sees force(code) from the local RNG wrapper as its call.
    # Its source reference retains a mutable development-source environment;
    # prediction needs only the forest. The public training call is saved by
    # life_train, so discard this misleading and potentially large backend call.
    model$call <- NULL
  }
  structure(list(type = "survival_forest", model = model, tau = tau,
    x_names = colnames(x), predictor_names = colnames(x), n_predictors = ncol(x),
    options = o, controls = o, seed = stored_seed,
    # Materialize sort's lazy wrapper once, so lookups cannot change the saved
    # object's serialized representation by populating an internal vector cache.
    jump_times = if (is.null(model)) numeric() else sort(unique(end[event > 0L])) + 0,
    diagnostics = list(n_subjects = risk$n_subjects, n_risk_subjects = nrow(risk$data),
      excluded_origin_absorption = risk$excluded_origin_absorption,
      n_censor_events = sum(event), backend = "ranger",
      backend_fitted = !is.null(model), resolved_options = o, seed = stored_seed,
      grid = "all censoring-event times; no administrative horizon jump",
      survival_floor = o$survival_floor)),
    class = "life_censor_forest")
}

#' Obtain batched coherent survival curves from a censoring forest
#' @noRd
censor_forest_curves <- function(adapter, x) {
  x <- flex_x(x, names = adapter$x_names)
  if (ncol(x) != adapter$n_predictors) stop("Wrong predictor count.", call. = FALSE)
  if (!nrow(x) || is.null(adapter$model)) return(list(
    jump_times = adapter$jump_times,
    survival = matrix(1, nrow(x), length(adapter$jump_times)), floored = 0L))
  if (!requireNamespace("ranger", quietly = TRUE)) stop("Install 'ranger' to predict this saved forest.", call. = FALSE)
  if (!ncol(x)) x <- matrix(0, nrow(x), 1L, dimnames = list(NULL, ".constant"))
  if (is.null(colnames(x))) colnames(x) <- paste0("x", seq_len(ncol(x)))
  # ranger otherwise samples a prediction seed from R even for deterministic
  # survival curves. Freeze its seed and preserve the caller's RNG state.
  pred <- flex_with_seed(adapter$seed,
    stats::predict(adapter$model, data = as.data.frame(x),
                   num.threads = adapter$options$num_threads, seed = adapter$seed))
  g <- matrix(pred$survival, nrow = nrow(x))
  times <- pred$unique.death.times
  if (ncol(g) != length(times) || any(!is.finite(g)) || any(g < 0 | g > 1) ||
      (ncol(g) > 1L && any(g[, -1L, drop = FALSE] -
                            g[, -ncol(g), drop = FALSE] > 1e-12)))
    stop("The forest returned an invalid survival curve.", call. = FALSE)
  below <- g < adapter$options$survival_floor
  g[below] <- adapter$options$survival_floor
  selected <- match(adapter$jump_times, times)
  if (anyNA(selected)) stop("The forest omitted a requested censoring-event knot.", call. = FALSE)
  omitted <- setdiff(seq_along(times), selected)
  if (length(omitted)) {
    previous <- cbind(1, g)[, omitted, drop = FALSE]
    if (any(abs(previous - g[, omitted, drop = FALSE]) > 1e-12))
      stop("The forest returned a survival jump outside censoring-event support.", call. = FALSE)
  }
  list(jump_times = adapter$jump_times, survival = g[, selected, drop = FALSE],
       floored = sum(below[, selected, drop = FALSE]))
}

#' Evaluate one censor-survival value per subject and query time
#' @noRd
predict_censor_forest <- function(adapter, x, time, side = c("left", "right"),
                                  subject_index = NULL) {
  side <- match.arg(side)
  x <- flex_x(x, names = adapter$x_names)
  query <- .cox_censor_queries(adapter, x, time, subject_index)
  curves <- censor_forest_curves(adapter, x)
  j <- findInterval(query$time, curves$jump_times, left.open = side == "left")
  answer <- rep(1, length(query$time))
  use <- j > 0L
  answer[use] <- curves$survival[cbind(query$subject_index[use], j[use])]
  answer
}

#' Derive both forest censoring increments from one set of fitted curves
#' @noRd
predict_censor_forest_jumps <- function(adapter, x, time, subject_index = NULL) {
  query <- .cox_censor_queries(adapter, x, time, subject_index)
  curves <- censor_forest_curves(adapter, x)
  lookup <- function(left) {
    j <- findInterval(query$time, curves$jump_times, left.open = left)
    out <- rep(1, length(j)); use <- j > 0L
    out[use] <- curves$survival[cbind(query$subject_index[use], j[use])]
    out
  }
  gl <- lookup(TRUE); gr <- lookup(FALSE)
  data.frame(subject_index = query$subject_index, time = query$time,
    g_left = gl, g_right = gr, hazard_increment = (gl - gr) / gl,
    inverse_increment = (gl - gr) / gl / gr)
}

#' Construct predictable features for flexible transition intensities
#' @noRd
flex_features <- function(x, time, duration, previous_transitions, accrued) {
  f <- cbind(x, time, duration, previous_transitions, accrued)
  colnames(f) <- paste0("f", seq_len(ncol(f)))
  f
}

#' Build recurrent-transition likelihood quadrature rows
#'
#' Every positive observed risk interval, including follow-up beyond the target
#' horizon, contributes its exact exposure and transition count. Numerical rate
#' cells are anchored to state entry and have width step_size. The rate in a cell
#' uses the midpoint history under continued residence in the origin state; this
#' deterministic extension uses no future observed event or censoring time. A
#' final partial cell retains its actual exposure and the same fixed midpoint.
#' Thus the exposure/count likelihood is exact for the declared piecewise-
#' constant working intensity. Subject indices support whole-subject resampling.
#' @noRd
make_flexible_transition_data <- function(history, x, problem, step_size,
                                           max_rows = 1000000L) {
  x <- flex_x(x, nrow(history$subjects))
  flex_positive(step_size, "step_size")
  flex_positive(max_rows, "max_rows", integer = TRUE)
  z <- history$intervals
  lengths <- z$stop - z$start
  lengths[z$from %in% problem$absorbing] <- 0
  pieces <- ceiling(lengths / step_size)
  # Division can round an exact cell count upward. Remove only a terminal
  # cell whose anchored left endpoint is already at the observed endpoint;
  # a representable positive final exposure must remain, however small.
  last <- which(pieces > 1L & is.finite(pieces))
  last <- last[z$start[last] + (pieces[last] - 1L) * step_size >= z$stop[last]]
  pieces[last] <- pieces[last] - 1L
  if (sum(pieces) > max_rows) stop("Flexible transition quadrature exceeds max_training_rows; increase step_size or the explicit budget.", call. = FALSE)
  n <- sum(pieces)
  p <- ncol(x) + 3L + ncol(problem$rewards)
  feature <- matrix(0, n, p, dimnames = list(NULL, paste0("f", seq_len(p))))
  exposure <- numeric(n); origin <- destination <- subject <- integer(n)
  offset <- 0L
  for (i in seq_len(nrow(history$subjects))) {
    accrued <- numeric(ncol(problem$rewards)); count <- 0L
    for (j in history$row_index[[i]]) {
      np <- pieces[j]
      if (np > 0L) {
        left <- z$start[j] + (seq_len(np) - 1L) * step_size
        # Share each interior boundary between neighboring cells rather than
        # recomputing it by a differently associated floating-point addition.
        right <- c(left[-1L], z$stop[j])
        if (any(right <= left))
          stop("step_size is too small for the numerical resolution of the observed times.", call. = FALSE)
        mid <- left + step_size / 2
        rows <- offset + seq_len(np)
        aa <- matrix(accrued, np, length(accrued), byrow = TRUE) +
          outer(mid - z$start[j], problem$rewards[z$from[j], ])
        feature[rows, ] <- flex_features(x[rep(i, np), , drop = FALSE], mid,
                                         mid - z$start[j], rep(count, np), aa)
        exposure[rows] <- right - left
        origin[rows] <- z$from[j]; subject[rows] <- i
        if (z$event[j]) destination[rows[length(rows)]] <- z$to[j]
        offset <- offset + np
        accrued <- accrued + lengths[j] * problem$rewards[z$from[j], ]
      }
      if (z$event[j]) count <- count + 1L
    }
  }
  list(features = feature, exposure = exposure, from = origin, to = destination,
       subject_index = subject, n_subjects = nrow(history$subjects),
       step_size = step_size,
       training_end = if (any(lengths > 0)) max(z$stop[lengths > 0]) else 0,
       approximation = "state-entry-anchored intensity cells with midpoint features")
}

#' Resolve common flexible-outcome settings
#' @noRd
flex_outcome_options <- function(options, problem, type) {
  common <- list(step_size = problem$tau / 50, max_training_rows = 1000000L,
                 max_events = 1000L)
  if (type == "transition_forest") common <- c(common,
    list(num_trees = 100L, max_depth = 8L, minsplit = 20L, minbucket = 7L,
         mtry = NULL, cp = 0.001))
  else common <- c(common, list(hidden = 4L, decay = 0.01, maxit = 1000L,
                                reltol = 1e-7))
  o <- flex_options(options, common)
  for (nm in names(o)) {
    if (nm == "mtry" && is.null(o[[nm]])) next
    flex_positive(o[[nm]], nm,
    integer = nm %in% c("max_training_rows", "max_events", "num_trees",
                       "max_depth", "minsplit", "minbucket", "mtry", "hidden", "maxit"),
    zero = nm %in% c("decay", "cp"))
  }
  if (type == "transition_forest" && o$max_depth > 30L) stop("max_depth cannot exceed 30.", call. = FALSE)
  o
}

#' Fit subject-resampled survival trees to recurrent transition intensities
#'
#' Each tree fits a Poisson exposure/count likelihood to one transition's risk
#' intervals. Competing exits contribute their exposure and no event for that
#' transition. Whole subjects are bootstrapped; all their repeated visits retain
#' the same multiplicity. Random feature subsets are selected once per tree.
#' @noRd
fit_transition_forest <- function(history, x, problem, options = list(), seed = NULL) {
  if (!requireNamespace("rpart", quietly = TRUE)) stop("Install the optional 'rpart' package for outcome_model = 'transition_forest'.", call. = FALSE)
  o <- flex_outcome_options(options, problem, "transition_forest")
  data <- make_flexible_transition_data(history, x, problem, o$step_size, o$max_training_rows)
  edges <- problem$working_graph$edges[, c("from", "to"), drop = FALSE]
  if (is.null(o$mtry)) o$mtry <- max(1L, floor(sqrt(ncol(data$features))))
  if (o$mtry > ncol(data$features)) stop("mtry exceeds the feature count.", call. = FALSE)
  stored_seed <- flex_with_seed(seed, sample.int(.Machine$integer.max, 1L))
  # Optional exact split tables share a 64 MiB budget across all transitions.
  # A large graph keeps checked tree traversal after exhausting this budget.
  remaining_table_budget <- 64 * 1024^2
  fits <- flex_with_seed(stored_seed, lapply(seq_len(nrow(edges)), function(e) {
    rows <- which(data$from == edges$from[e])
    count <- as.integer(data$to[rows] == edges$to[e])
    if (!sum(count)) return(list(constant = 0, trees = NULL))
    f <- as.data.frame(data$features[rows, , drop = FALSE])
    f$.exposure <- data$exposure[rows]; f$.event <- count
    groups <- data$subject_index[rows]
    trees <- lapply(seq_len(o$num_trees), function(b) {
      mult <- tabulate(sample.int(data$n_subjects, data$n_subjects, replace = TRUE),
                       nbins = data$n_subjects)
      weights <- mult[groups]
      use <- which(weights > 0L)
      selected <- sample.int(ncol(data$features), o$mtry)
      columns <- colnames(data$features)[selected]
      if (!length(use) || !sum(count[use] * weights[use]))
        return(list(constant = 0, model = NULL, columns = columns))
      form <- stats::reformulate(columns, response = "cbind(.exposure, .event)")
      fit <- rpart::rpart(form, data = f[use, , drop = FALSE], weights = weights[use],
        method = "poisson", control = rpart::rpart.control(cp = o$cp,
          minsplit = o$minsplit, minbucket = o$minbucket, maxdepth = o$max_depth,
          xval = 0L, maxsurrogate = 0L, usesurrogate = 0L), model = FALSE, x = FALSE, y = FALSE)
      # Formula evaluation is unnecessary after fitting. Remove its temporary
      # training environment; prediction uses stored terms and numeric columns.
      environment(fit$terms) <- baseenv()
      fit$call <- NULL
      # New-data prediction does not need the leaf index of every training
      # quadrature row, often the largest component of a small fitted tree.
      fit$where <- NULL
      fit$variable.importance <- NULL
      list(constant = NULL, model = fit, columns = columns)
    })
    compiled <- compile_flexible_forest(trees, colnames(data$features),
      table_budget = min(24 * 1024^2, remaining_table_budget))
    used <- sum(vapply(compiled[c("lookup", "additive")], function(z)
      if (is.null(z)) 0 else as.numeric(utils::object.size(z)), numeric(1)))
    remaining_table_budget <<- max(0, remaining_table_budget - used)
    list(constant = NULL, trees = trees, compiled = compiled)
  }))
  structure(list(type = "transition_forest", fits = fits, edges = edges,
    problem = problem, options = o, x_names = colnames(x), n_predictors = ncol(x),
    seed = stored_seed, diagnostics = list(n_quadrature_rows = nrow(data$features),
      quadrature = data$approximation, resampling = "whole subjects",
      training_end = data$training_end, fitting_follow_up = "full observed X",
      backend = "rpart Poisson exposure/count trees", step_size = o$step_size,
      prediction_engine = "additive split-bin tables with checked traversal fallback",
      resolved_options = o, seed = stored_seed)),
    class = "life_flexible_outcome")
}

#' Store direct numeric tree traversal arrays once after fitting
#'
#' The transition learner uses complete numeric features and disables
#' surrogates. A split's parity determines which child receives values below
#' its cut; equality takes the other child, matching predict.rpart. Competitor
#' splits remain in the backend fit but are not used for prediction.
#' @noRd
compile_flexible_forest <- function(trees, feature_names, table_budget = 24 * 1024^2) {
  counts <- vapply(trees, function(tree)
    if (is.null(tree$constant)) nrow(tree$model$frame) else 1L, integer(1))
  total <- sum(counts)
  if (!length(trees) || !is.finite(total) || total > .Machine$integer.max)
    stop("Transition forest exceeds the supported node count.", call. = FALSE)
  result <- list(roots = integer(length(trees)), feature = integer(total),
    left = integer(total), right = integer(total), direction = integer(total),
    cut = numeric(total), value = numeric(total), n_features = length(feature_names))
  offset <- 0L
  for (b in seq_along(trees)) {
    tree <- trees[[b]]
    rows <- offset + seq_len(counts[b])
    if (!is.null(tree$constant)) {
      result$roots[b] <- rows[1L]
      result$value[rows] <- tree$constant
    } else {
      frame <- tree$model$frame
      node <- as.numeric(row.names(frame))
      root <- match(1, node)
      branch <- which(frame$var != "<leaf>")
      if (is.na(root) || anyDuplicated(node) || any(!is.finite(frame$yval)) ||
          any(frame$yval < 0))
        stop("The fitted transition tree has invalid nodes or rates.", call. = FALSE)
      result$roots[b] <- offset + root
      result$value[rows] <- frame$yval
      if (length(branch)) {
        first <- 1L + c(0L, cumsum((frame$var != "<leaf>") +
          frame$ncompete + frame$nsurrogate))[-(nrow(frame) + 1L)]
        split <- tree$model$splits[first[branch], , drop = FALSE]
        feature <- match(frame$var[branch], feature_names)
        left <- match(2 * node[branch], node)
        right <- match(2 * node[branch] + 1, node)
        if (anyNA(c(feature, left, right)) ||
            any(!split[, "ncat"] %in% c(-1, 1)) ||
            any(!is.finite(split[, "index"])))
          stop("The transition forest requires complete numeric tree splits.", call. = FALSE)
        at <- rows[branch]
        result$feature[at] <- feature
        result$left[at] <- offset + left
        result$right[at] <- offset + right
        result$direction[at] <- as.integer(split[, "ncat"])
        result$cut[at] <- split[, "index"]
      }
    }
    offset <- offset + counts[b]
  }
  # A leaf's half-open numeric box is exactly the intersection of its splits.
  # Sparse bounds permit simulation-local reuse without changing any tree cut.
  result$tree_start <- cumsum(c(1L, counts[-length(counts)]))
  result$tree_end <- cumsum(counts)
  result$box_start <- result$box_count <- integer(total)
  boxes <- vector("list", total)
  visit <- function(node, lower, upper) {
    column <- result$feature[node]
    if (!column) {
      used <- which(is.finite(lower) | is.finite(upper))
      result$box_count[node] <<- length(used)
      boxes[[node]] <<- list(feature = used, lower = lower[used], upper = upper[used])
      return(invisible(NULL))
    }
    below <- if (result$direction[node] == -1L) result$left[node] else result$right[node]
    above <- if (result$direction[node] == -1L) result$right[node] else result$left[node]
    cut <- result$cut[node]
    below_upper <- upper; below_upper[column] <- min(upper[column], cut)
    visit(below, lower, below_upper)
    above_lower <- lower; above_lower[column] <- max(lower[column], cut)
    visit(above, above_lower, upper)
  }
  for (root in result$roots)
    visit(root, rep(-Inf, length(feature_names)), rep(Inf, length(feature_names)))
  result$box_start <- c(0L, cumsum(result$box_count))[-(total + 1L)]
  result$box_feature <- as.integer(unlist(lapply(boxes, `[[`, "feature"), use.names = FALSE))
  result$box_lower <- as.numeric(unlist(lapply(boxes, `[[`, "lower"), use.names = FALSE))
  result$box_upper <- as.numeric(unlist(lapply(boxes, `[[`, "upper"), use.names = FALSE))
  result$lookup <- compile_forest_lookup(result, max_bytes = min(8 * 1024^2, table_budget))
  remaining <- table_budget - if (is.null(result$lookup)) 0 else as.numeric(utils::object.size(result$lookup))
  result$additive <- compile_forest_additive(result, max_bytes = min(16 * 1024^2, remaining))
  result
}

#' Compile bounded exact leaf tables for trees using one or two features
#'
#' Raw split cuts define half-open bins. Tables store original leaf indices,
#' never rounded rates; predictions still sum trees in their original order.
#' Large/high-dimensional trees retain checked traversal. The stored optional
#' lookup has an 8 MiB cap, and each transient tree grid has at most 65,536 cells.
#' @noRd
compile_forest_lookup <- function(compiled, max_bytes = 8 * 1024^2,
                                  max_cells = 65536L) {
  trees <- length(compiled$roots); nodes <- length(compiled$feature)
  # Bound optional preprocessing before constructing lists of split cuts.
  if (!trees || 32 * nodes + 32 * trees > max_bytes) return(NULL)
  parts <- vector("list", trees)
  for (b in seq_len(trees)) {
    at <- seq.int(compiled$tree_start[b], compiled$tree_end[b])
    at <- at[compiled$feature[at] > 0L]
    variables <- sort(unique(compiled$feature[at]))
    if (!length(variables) || length(variables) > 2L) next
    cuts <- lapply(variables, function(v) sort(unique(compiled$cut[at[compiled$feature[at] == v]])))
    size <- lengths(cuts) + 1L
    cells <- prod(size)
    if (!is.finite(cells) || cells > max_cells) next
    parts[[b]] <- list(variables = variables, cuts = cuts, size = size, cells = cells)
  }
  variables <- sort(unique(as.integer(unlist(lapply(parts, `[[`, "variables"), use.names = FALSE))))
  k <- length(variables)
  if (!k) return(NULL)
  global <- lapply(variables, function(v) sort(unique(as.numeric(unlist(lapply(parts, function(part) {
    where <- match(v, part$variables)
    if (is.na(where)) NULL else part$cuts[[where]]
  }), use.names = FALSE)))))
  counts <- lengths(global)
  used_bytes <- 4096 + 8 * sum(counts) + 12 * k + 28 * trees
  if (used_bytes > max_bytes) return(NULL)
  result <- list(variables = variables,
    cut_start = as.integer(c(0, cumsum(counts))[-(k + 1L)]),
    cut_count = as.integer(counts), cuts = as.numeric(unlist(global, use.names = FALSE)),
    feature1 = integer(trees), feature2 = integer(trees),
    offset = rep.int(-1L, trees), size1 = rep.int(1L, trees), size2 = rep.int(1L, trees),
    map_start1 = integer(trees), map_start2 = integer(trees))
  leaves <- maps <- vector("list", trees)
  leaf_offset <- map_offset <- 0L
  for (b in seq_len(trees)) {
    part <- parts[[b]]
    if (is.null(part)) next
    slots <- match(part$variables, variables)
    map_sizes <- counts[slots] + 1L
    required <- 4 * (part$cells + sum(map_sizes))
    if (used_bytes + required > max_bytes) next
    n1 <- part$size[1L]
    n2 <- if (length(part$size) == 2L) part$size[2L] else 1L
    x1 <- rep(c(-Inf, part$cuts[[1L]]), times = n2)
    x2 <- if (n2 > 1L) rep(c(-Inf, part$cuts[[2L]]), each = n1) else NULL
    leaf <- rep.int(compiled$roots[b], part$cells)
    steps <- 0L
    repeat {
      active <- which(compiled$feature[leaf] != 0L)
      if (!length(active)) break
      steps <- steps + 1L
      if (steps > nodes) stop("Invalid cycle while compiling forest lookup.", call. = FALSE)
      at <- leaf[active]
      value <- x1[active]
      if (!is.null(x2)) {
        second <- compiled$feature[at] == part$variables[2L]
        value[second] <- x2[active[second]]
      }
      below <- value < compiled$cut[at]
      left <- ifelse(compiled$direction[at] == -1L, below, !below)
      leaf[active] <- ifelse(left, compiled$left[at], compiled$right[at])
    }
    map1 <- c(0L, findInterval(global[[slots[1L]]], part$cuts[[1L]]))
    map2 <- if (length(slots) == 2L)
      c(0L, findInterval(global[[slots[2L]]], part$cuts[[2L]])) else integer()
    result$feature1[b] <- slots[1L]
    result$feature2[b] <- if (length(slots) == 2L) slots[2L] else 0L
    result$offset[b] <- leaf_offset
    result$size1[b] <- n1; result$size2[b] <- n2
    result$map_start1[b] <- map_offset
    result$map_start2[b] <- map_offset + length(map1)
    leaves[[b]] <- as.integer(leaf); maps[[b]] <- as.integer(c(map1, map2))
    leaf_offset <- leaf_offset + length(leaf)
    map_offset <- map_offset + length(maps[[b]])
    used_bytes <- used_bytes + required
  }
  if (!any(result$offset >= 0L)) return(NULL)
  result$leaf <- as.integer(unlist(leaves, use.names = FALSE))
  result$map <- as.integer(unlist(maps, use.names = FALSE))
  result$complete <- all(result$offset >= 0L | compiled$feature[compiled$roots] == 0L)
  result$storage_bytes <- 0
  result$storage_bytes <- as.numeric(utils::object.size(result))
  if (result$storage_bytes > max_bytes) return(NULL)
  result
}

#' Stable softplus rate link and derivative
#' @noRd
flex_softplus <- function(z) pmax(z, 0) + log1p(exp(-abs(z))) + 1e-10

#' Neural intensity objective with its analytic gradient
#' @noRd
flex_neural_objective <- function(par, design, exposure, count, hidden, decay,
                                  gradient = FALSE, normalizer = sum(exposure)) {
  p <- ncol(design); nw <- p * hidden
  w <- matrix(par[seq_len(nw)], p, hidden)
  v <- par[nw + seq_len(hidden + 1L)]
  h <- tanh(design %*% w)
  hb <- cbind(1, h)
  eta <- drop(hb %*% v)
  rate <- flex_softplus(eta)
  if (!gradient) return(sum(exposure * rate - count * log(rate)) / normalizer +
                        decay * sum(par[-(nw + 1L)]^2) / 2)
  d <- (exposure - count / rate) * stats::plogis(eta) / normalizer
  gv <- drop(crossprod(hb, d))
  gw <- crossprod(design, (1 - h^2) * outer(d, v[-1L]))
  penalty <- decay * par; penalty[nw + 1L] <- 0
  c(gw, gv) + penalty
}

#' Share exact hidden-layer calculations between optimizer callbacks
#'
#' This cache exists only for one transition fit. Line-search values at new
#' parameter vectors are recomputed; objective and gradient requests at the
#' identical vector share their activations and rates. No cache is saved in
#' the fitted adapter, and optimizer settings/evaluation counts are unchanged.
#' @noRd
flex_neural_evaluator <- function(design, exposure, count, hidden, decay,
                                   normalizer) {
  p <- ncol(design); nw <- p * hidden
  last <- NULL; current <- NULL; derivative <- NULL
  core_evaluations <- 0L
  evaluate <- function(par, gradient = FALSE) {
    if (!identical(par, last)) {
      w <- matrix(par[seq_len(nw)], p, hidden)
      v <- par[nw + seq_len(hidden + 1L)]
      h <- tanh(design %*% w)
      hb <- cbind(1, h)
      eta <- drop(hb %*% v)
      rate <- flex_softplus(eta)
      current <<- list(h = h, hb = hb, eta = eta, rate = rate, v = v,
        value = sum(exposure * rate - count * log(rate)) / normalizer +
          decay * sum(par[-(nw + 1L)]^2) / 2)
      # Force a private numeric copy, independent of optimizer argument reuse.
      last <<- par + 0
      derivative <<- NULL
      core_evaluations <<- core_evaluations + 1L
    }
    if (!gradient) return(current$value)
    if (is.null(derivative)) {
      d <- (exposure - count / current$rate) * stats::plogis(current$eta) / normalizer
      gv <- drop(crossprod(current$hb, d))
      gw <- crossprod(design, (1 - current$h^2) * outer(d, current$v[-1L]))
      penalty <- decay * par; penalty[nw + 1L] <- 0
      derivative <<- c(gw, gv) + penalty
    }
    derivative
  }
  list(fn = function(par) evaluate(par), gr = function(par) evaluate(par, TRUE),
       count = function() core_evaluations)
}

#' Fit a neural counting-process intensity model
#'
#' A one-hidden-layer tanh network with a positive softplus rate is fitted to
#' each transition's exposure/count likelihood using an analytic gradient and
#' BFGS. Predictable features and their training standardization are frozen.
#' This is a recurrent intensity likelihood, not a right-censored-row classifier.
#' @noRd
fit_transition_neural <- function(history, x, problem, options = list(), seed = NULL) {
  o <- flex_outcome_options(options, problem, "transition_neural")
  data <- make_flexible_transition_data(history, x, problem, o$step_size, o$max_training_rows)
  edges <- problem$working_graph$edges[, c("from", "to"), drop = FALSE]
  stored_seed <- flex_with_seed(seed, sample.int(.Machine$integer.max, 1L))
  fits <- flex_with_seed(stored_seed, lapply(seq_len(nrow(edges)), function(e) {
    rows <- which(data$from == edges$from[e])
    count <- as.integer(data$to[rows] == edges$to[e])
    if (!sum(count)) return(list(constant = 0))
    f <- data$features[rows, , drop = FALSE]
    weights <- data$exposure[rows] / sum(data$exposure[rows])
    center <- colSums(f * weights)
    scale <- sqrt(colSums(sweep(f, 2L, center)^2 * weights)); scale[scale == 0] <- 1
    design <- cbind(1, sweep(sweep(f, 2L, center), 2L, scale, "/"))
    nw <- ncol(design) * o$hidden
    par <- stats::rnorm(nw + o$hidden + 1L, sd = 0.05)
    rate <- sum(count) / sum(data$exposure[rows])
    par[nw + 1L] <- if (rate > 30) rate else log(expm1(rate))
    evaluator <- flex_neural_evaluator(design, data$exposure[rows], count,
      hidden = o$hidden, decay = o$decay, normalizer = data$n_subjects)
    fit <- stats::optim(par, evaluator$fn, gr = evaluator$gr, method = "BFGS",
      control = list(maxit = o$maxit, reltol = o$reltol))
    if (fit$convergence != 0L || any(!is.finite(fit$par)))
      stop("Neural transition fit did not converge for edge ",
           as.character(problem$states[edges$from[e]]), " -> ",
           as.character(problem$states[edges$to[e]]),
           "; inspect model/data scaling and consider increasing maxit or decay.",
           call. = FALSE)
    list(constant = NULL, w = matrix(fit$par[seq_len(nw)], ncol(design), o$hidden),
      v = fit$par[nw + seq_len(o$hidden + 1L)], center = center, scale = scale,
      convergence = fit$convergence, objective = fit$value, evaluations = fit$counts,
      shared_core_evaluations = evaluator$count())
  }))
  structure(list(type = "transition_neural", fits = fits, edges = edges,
    problem = problem, options = o, x_names = colnames(x), n_predictors = ncol(x),
    seed = stored_seed, diagnostics = list(n_quadrature_rows = nrow(data$features),
      quadrature = data$approximation, backend = "tanh/softplus intensity network; stats::optim BFGS",
      training_end = data$training_end, fitting_follow_up = "full observed X",
      likelihood_normalization = "per training subject", step_size = o$step_size,
      resolved_options = o, seed = stored_seed)),
    class = "life_flexible_outcome")
}

#' Validate a flexible continuation snapshot without predicting rates
#' @noRd
flex_validate_snapshot <- function(adapter, x, time, state, duration,
                                   previous_transitions, accrued,
                                   within_horizon = FALSE) {
  if (!inherits(adapter, "life_flexible_outcome") ||
      !adapter$type %in% c("transition_forest", "transition_neural"))
    stop("A flexible transition model is required; Markov models use deterministic recursions.", call. = FALSE)
  x <- flex_x(x, names = adapter$x_names); n <- nrow(x)
  if (ncol(x) != adapter$n_predictors || any(lengths(list(time, state, duration,
      previous_transitions)) != n) || !is.matrix(accrued) ||
      nrow(accrued) != n || ncol(accrued) != ncol(adapter$problem$rewards) ||
      any(!is.finite(c(time, state, duration, previous_transitions, accrued))) ||
      any(time < 0 | duration < 0 | duration > time) ||
      (within_horizon && any(time > adapter$problem$tau)) ||
      any(state < 1 | state > nrow(adapter$problem$rewards) | state != floor(state)) ||
      any(previous_transitions < 0 | previous_transitions != floor(previous_transitions)))
    stop("Invalid flexible-model history features.", call. = FALSE)
  invisible(n)
}

#' Evaluate flexible transition intensities in a batch of partial histories
#' @noRd
predict_transition_rates <- function(adapter, x, time, state, duration,
                                      previous_transitions, accrued,
                                      .validated = FALSE, .cache = NULL,
                                      .cache_rows = NULL) {
  if (!.validated)
    flex_validate_snapshot(adapter, x, time, state, duration, previous_transitions, accrued)
  n <- nrow(x)
  rates <- matrix(0, n, nrow(adapter$edges))
  if (!n || !ncol(rates)) return(rates)
  f <- flex_features(x, time, duration, previous_transitions, accrued)
  # Public inputs were checked at simulation entry; arithmetic can still
  # overflow later when rewards or times approach floating-point limits.
  if (.validated && any(!is.finite(f)))
    stop("Invalid flexible-model history features.", call. = FALSE)
  for (e in seq_len(nrow(adapter$edges))) {
    rows <- which(state == adapter$edges$from[e] & !(state %in% adapter$problem$absorbing))
    if (!length(rows)) next
    fit <- adapter$fits[[e]]
    if (!is.null(fit$constant)) { rates[rows, e] <- fit$constant; next }
    if (adapter$type == "transition_forest") {
      if (!is.null(fit$compiled)) {
        rates[rows, e] <- if (is.null(.cache[[e]]))
          .Call(life_forest_rates, fit$compiled, f, rows) else
          .Call(life_forest_rates_cached, fit$compiled, f, rows, .cache[[e]], .cache_rows)
      } else {
        # Older serialized adapters retain their supported backend route.
        if (!requireNamespace("rpart", quietly = TRUE))
          stop("Install 'rpart' to predict this saved transition forest.", call. = FALSE)
        ff <- as.data.frame(f[rows, , drop = FALSE])
        total <- numeric(length(rows))
        for (tree in fit$trees) total <- total + if (!is.null(tree$constant)) tree$constant else
          as.numeric(stats::predict(tree$model, newdata = ff, type = "vector"))
        rates[rows, e] <- total / length(fit$trees)
      }
    } else {
      design <- .Call(life_flexible_design, f, rows, fit$center, fit$scale)
      rates[rows, e] <- flex_softplus(drop(cbind(1, tanh(design %*% fit$w)) %*% fit$v))
    }
  }
  if (any(!is.finite(rates)) || any(rates < 0)) stop("Nonfinite flexible transition rate.", call. = FALSE)
  rates
}

#' Assign stable simulation streams to an explicit numerical profile
#' @noRd
flex_profile_seed <- function(seed, profile) {
  bytes <- charToRaw(paste(sprintf("%a", as.numeric(profile)), collapse = ","))
  value <- as.double(seed)
  for (b in as.integer(bytes)) value <- (value * 131 + b) %% 2147483646
  as.integer(value + 1)
}

#' Simulate only future rewards from flexible history-dependent models
#'
#' Intensities are frozen on state-entry-anchored cells of width step_size,
#' evaluated at each cell's midpoint under continued residence. A conditional
#' query inside a cell retains its original boundary and rate, so it conditions
#' the same fitted temporal law rather than starting a new numerical clock.
#' A sampled departure is continuous within the remaining cell; its destination
#' is sampled from competing intensity ratios.
#' The clock, sojourn duration, prior count and rewards update after every event.
#' Streams depend on the numerical profile and snapshot, not row order or chunk
#' size. Uniform banks are reused over score thresholds. The explicit event
#' budget raises an error rather than truncating an explosive path silently.
#' @noRd
simulate_flexible_targets <- function(adapter, x, snapshots, B = 256L,
                                      seed = adapter$seed) {
  owned_bank <- NULL
  on.exit(if (!is.null(owned_bank)) .Call(life_forest_bank_release, owned_bank), add = TRUE)
  x <- flex_x(x, names = adapter$x_names); n <- nrow(x)
  flex_positive(B, "B", integer = TRUE)
  flex_positive(seed, "seed", integer = TRUE)
  required <- c("time", "state", "duration", "previous_transitions", "accrued")
  if (!all(required %in% names(snapshots))) stop("Incomplete flexible continuation snapshot.", call. = FALSE)
  # The same checker enforces the full numerical snapshot contract once.
  flex_validate_snapshot(adapter, x, snapshots$time, snapshots$state,
    snapshots$duration, snapshots$previous_transitions, snapshots$accrued,
    within_horizon = TRUE)
  d <- ncol(adapter$problem$rewards); tau <- adapter$problem$tau
  answer <- array(0, c(n, B, d), dimnames = list(NULL, NULL, colnames(adapter$problem$rewards)))
  if (!n) return(answer)
  step <- adapter$options$step_size
  max_iter <- ceiling(tau / step) + adapter$options$max_events + 2L
  if (!is.finite(max_iter) || max_iter > .Machine$integer.max / (2 * B))
    stop("Continuation uniform-bank budget is too large.", call. = FALSE)
  # Working-graph sinks have deterministic future rewards under this adapter.
  # This shortcut is a model calculation, never a structural CDT certificate.
  model_terminal <- union(adapter$problem$absorbing,
    setdiff(seq_len(nrow(adapter$problem$rewards)), adapter$edges$from))
  known <- which(snapshots$time == tau | snapshots$state %in% model_terminal)
  for (i in known) {
    value <- markov_protected_total(tau, snapshots$time[i],
      adapter$problem$rewards[snapshots$state[i], ], snapshots$accrued[i, ])
    answer[i, , ] <- matrix(value, B, d, byrow = TRUE)
  }
  unknown <- setdiff(seq_len(n), known)
  # Fixed-size independent stream blocks avoid generating the maximum-event
  # budget's random numbers for paths that finish early. A block key includes
  # the numerical profile and block number, never the batch or draw count.
  uniform_steps <- 32L
  cache_trees <- if (adapter$type == "transition_forest")
    vapply(adapter$fits, function(fit) if (!is.null(fit$compiled$box_start) &&
      !isTRUE(fit$compiled$lookup$complete))
      length(fit$compiled$roots) else 0L, integer(1)) else integer(length(adapter$fits))
  # Disable optional caching if even one profile would need over 32 MiB.
  # This changes only lookup reuse, never draws, tree count or the fitted law.
  if (4 * B * sum(cache_trees) > 32 * 1024^2) cache_trees[] <- 0L
  # Uniform banks plus optional private leaf indices share a 64 MiB batch
  # budget. One profile is the irreducible unit; no fitted/user object is mutated.
  bytes_per_profile <- B * (2 * uniform_steps * 8 + sum(cache_trees) * 4)
  batch_size <- max(1L, min(16L, floor(64 * 1024^2 / bytes_per_profile)))
  # Complete additive forests can advance one existing uniform block in native
  # code. Eligibility is decided before advancing any path. Old/cached forests,
  # unsupported arithmetic and oversized workspaces retain this same R loop.
  native_forest <- FALSE
  if (length(unknown) && adapter$type == "transition_forest" && !any(cache_trees > 0L)) {
    native_n <- B * min(batch_size, length(unknown))
    integer_overflow <- is.integer(snapshots$previous_transitions) &&
      any(snapshots$previous_transitions > .Machine$integer.max - adapter$options$max_events)
    if (!integer_overflow && native_n <= .Machine$integer.max) {
      forest_precision <- flex_forest_block_precision()
      native_forest <- !is.null(forest_precision) && .Call(life_forest_block_supported,
        adapter, as.integer(native_n), as.integer(B), forest_precision)
    }
  }
  for (batch_start in seq.int(1L, max(1L, length(unknown)), by = batch_size)) {
    if (!length(unknown)) break
    ids <- unknown[seq.int(batch_start, min(length(unknown), batch_start + batch_size - 1L))]
    group <- rep(ids, each = B)
    generated_bank <- native_forest && flex_forest_bank_supported(adapter, group, B, forest_precision)
    leaf_cache <- lapply(cache_trees, function(trees) {
      if (!trees) return(NULL)
      cache <- new.env(parent = emptyenv())
      cache$nodes <- matrix(0L, length(group), trees)
      cache
    })
    tt <- snapshots$time[group]; state <- snapshots$state[group]
    duration <- snapshots$duration[group]
    entry <- tt - duration
    cell <- floor(duration / step)
    cell_end <- entry + (cell + 1) * step
    at_boundary <- cell_end <= tt
    if (any(at_boundary)) {
      cell[at_boundary] <- cell[at_boundary] + 1
      cell_end[at_boundary] <- entry[at_boundary] + (cell[at_boundary] + 1) * step
    }
    if (any(cell_end <= tt))
      stop("step_size is too small for the numerical resolution of query times.", call. = FALSE)
    previous <- snapshots$previous_transitions[group]
    accrued <- snapshots$accrued[group, , drop = FALSE]
    count <- integer(length(group))
    bank <- if (generated_bank) NULL else matrix(0, length(group), 2L * uniform_steps)
    profile_seeds <- integer(length(ids))
    for (ii in seq_along(ids)) {
      i <- ids[ii]
      profile <- c(x[i, ], snapshots$time[i], snapshots$state[i], snapshots$duration[i],
                   snapshots$previous_transitions[i], snapshots$accrued[i, ])
      profile_seeds[ii] <- flex_profile_seed(seed, profile)
    }
    iteration <- 0L
    repeat {
      absorbed <- state %in% model_terminal & tt < tau
      if (any(absorbed)) {
        aa <- which(absorbed)
        accrued[aa, ] <- accrued[aa, , drop = FALSE] +
          adapter$problem$rewards[state[aa], , drop = FALSE] * (tau - tt[aa])
        tt[aa] <- tau
      }
      active <- which(tt < tau)
      if (!length(active)) break
      iteration <- iteration + 1L
      if (iteration > max_iter) stop("Flexible continuation iteration budget exceeded.", call. = FALSE)
      uniform_column <- (iteration - 1L) %% uniform_steps + 1L
      if (uniform_column == 1L) {
        block <- (iteration - 1L) %/% uniform_steps + 1L
        if (generated_bank) {
          bank <- flex_forest_bank_generate(adapter, group, B, forest_precision, profile_seeds, block)
          owned_bank <- bank
        } else for (ii in seq_along(ids)) {
          rows <- (ii - 1L) * B + seq_len(B)
          bank[rows, ] <- flex_with_seed(flex_profile_seed(profile_seeds[ii], block),
            matrix(stats::runif(B * 2L * uniform_steps), B, byrow = TRUE))
        }
      }
      if (native_forest) {
        if (generated_bank) {
          native_result <- .Call(life_forest_bank_advance,
            adapter, x, as.integer(group),
            flex_forest_block_state(tt, state, entry, cell, cell_end, previous,
              accrued, count, iteration - 1L), bank,
            list(draws = as.integer(B), precision = forest_precision, steps = uniform_steps,
              max_iter = as.integer(max_iter), trace = 0L))
        } else {
          native_result <- .Call(life_forest_block_advance,
            adapter, x, as.integer(group),
            flex_forest_block_state(tt, state, entry, cell, cell_end, previous,
              accrued, count, iteration - 1L), bank,
            list(draws = as.integer(B), precision = forest_precision, steps = uniform_steps,
              max_iter = as.integer(max_iter), trace = 0L))
        }
        owned_bank <- NULL
        if (is.null(native_result))
          stop("Forest continuation eligibility changed after preflight.", call. = FALSE)
        # The native result preserves the validated state shapes. On return,
        # this loop flushes terminal paths before considering another bank.
        list2env(native_result$path, envir = environment())
        next
      }
      h <- pmin(cell_end[active], tau) - tt[active]
      if (any(h <= 0))
        stop("The flexible intensity-cell clock failed to advance.", call. = FALSE)
      reward <- adapter$problem$rewards[state[active], , drop = FALSE]
      midpoint <- entry[active] + (cell[active] + .5) * step
      shift <- midpoint - tt[active]
      # Until its first simulated departure, every draw from a given snapshot
      # has the same state-entry clock and accumulated history. Predict that
      # exact feature row once, then expand its rate; after departure every
      # path is evaluated individually. Random banks and event order are intact.
      representative <- seq_along(active)
      fresh <- which(count[active] == 0L)
      if (length(fresh)) representative[fresh] <-
        fresh[match(group[active[fresh]], group[active[fresh]])]
      selected <- which(representative == seq_along(active))
      rate_rows <- active[selected]
      rates <- predict_transition_rates(adapter, x[group[rate_rows], , drop = FALSE],
        midpoint[selected], state[rate_rows], (cell[rate_rows] + .5) * step,
        previous[rate_rows], accrued[rate_rows, , drop = FALSE] +
          reward[selected, , drop = FALSE] * shift[selected],
        .validated = TRUE, .cache = leaf_cache, .cache_rows = rate_rows)
      if (length(selected) != length(active))
        rates <- rates[match(representative, selected), , drop = FALSE]
      total <- rowSums(rates)
      if (any(!is.finite(total)))
        stop("Total outgoing transition intensity exceeds the numerical range.", call. = FALSE)
      wait <- rep(Inf, length(active)); positive <- total > 0
      wait[positive] <- -log(bank[active[positive], 2L * uniform_column - 1L]) / total[positive]
      event <- wait < h
      elapsed <- pmin(wait, h)
      accrued[active, ] <- accrued[active, , drop = FALSE] + reward * elapsed
      advanced <- tt[active] + elapsed
      if (any(event & advanced <= tt[active]))
        stop("A simulated waiting time is below the numerical resolution of the time scale.", call. = FALSE)
      tt[active] <- advanced
      no_event <- active[!event]
      tt[no_event] <- pmin(cell_end[no_event], tau)
      cell[no_event] <- cell[no_event] + 1
      cell_end[no_event] <- entry[no_event] + (cell[no_event] + 1) * step
      if (any(event)) {
        ae <- active[event]
        count[ae] <- count[ae] + 1L
        if (any(count[ae] > adapter$options$max_events))
          stop("Flexible continuation exceeded max_events; inspect fitted intensities or increase the explicit budget.", call. = FALSE)
        rr <- rates[event, , drop = FALSE]
        threshold <- bank[ae, 2L * uniform_column] * total[event]
        destination <- integer(length(ae))
        cumulative <- numeric(length(ae))
        for (e in seq_len(ncol(rr))) {
          cumulative <- cumulative + rr[, e]
          hit <- destination == 0L & threshold < cumulative
          destination[hit] <- adapter$edges$to[e]
        }
        if (any(destination == 0L)) stop("Failed competing-transition probability lookup.", call. = FALSE)
        state[ae] <- destination; entry[ae] <- tt[ae]; cell[ae] <- 0
        cell_end[ae] <- entry[ae] + step
        previous[ae] <- previous[ae] + 1L
      }
    }
    for (ii in seq_along(ids)) answer[ids[ii], , ] <-
      accrued[(ii - 1L) * B + seq_len(B), , drop = FALSE]
  }
  # A reward constant in every state has an exact deterministic future value;
  # preserve that atom rather than accumulate pathwise floating-point drift.
  constant <- which(vapply(seq_len(d), function(j)
    all(adapter$problem$rewards[, j] == adapter$problem$rewards[1L, j]), logical(1)))
  for (j in constant) for (i in seq_len(n)) answer[i, , j] <-
    markov_protected_total(tau, snapshots$time[i],
      adapter$problem$rewards[1L, j], snapshots$accrued[i, j])
  attr(answer, "diagnostics") <- list(method = "continuous departures on state-entry-anchored intensity cells",
    step_size = step, draws = B, seed = seed, profile_keyed_streams = TRUE)
  answer
}

#' Compute empirical moments without hiding underflow as determinism
#' @noRd
flex_empirical_moments <- function(values) {
  mean <- rowMeans(values)
  variance <- rowMeans((values - mean)^2)
  if (any(variance < .Machine$double.xmin &
          (variance > 0 | rowSums(values != values[, 1L]) > 0))) {
    stop("Score variance underflows at the current time unit; rescale time before fitting.",
         call. = FALSE)
  }
  list(mean = mean, variance = variance)
}

#' Compute frozen baseline moments from flexible continuation banks
#' @noRd
flexible_target_moments <- function(adapter, x, initial_state, B = 256L) {
  flex_x(x, names = adapter$x_names)
  flex_positive(B, "B", integer = TRUE)
  n <- nrow(x); d <- ncol(adapter$problem$rewards)
  if (length(initial_state) == 1L) initial_state <- rep(initial_state, n)
  if (length(initial_state) != n) stop("Supply one initial state per profile.", call. = FALSE)
  mean <- variance <- matrix(0, n, d, dimnames = list(NULL, colnames(adapter$problem$rewards)))
  diagnostics <- NULL
  # Moments need only a bounded group of target banks at once, not the full
  # new-sample by draws by targets array used by conditional score diagnostics.
  for (start in seq.int(1L, max(1L, n), by = 64L)) {
    if (!n) break
    rows <- seq.int(start, min(n, start + 63L)); nr <- length(rows)
    draws <- simulate_flexible_targets(adapter, x[rows, , drop = FALSE],
      list(time = rep(0, nr), state = initial_state[rows], duration = rep(0, nr),
           previous_transitions = rep(0, nr), accrued = matrix(0, nr, d)), B = B)
    for (j in seq_len(d)) {
      z <- matrix(draws[, , j], nr, B)
      moments <- flex_empirical_moments(z)
      mean[rows, j] <- moments$mean
      variance[rows, j] <- moments$variance
    }
    diagnostics <- attr(draws, "diagnostics")
  }
  list(mean = mean, variance = variance, second_moment = variance + mean^2,
       diagnostics = diagnostics)
}

#' Reuse joint continuation draws across standardized maximum-score thresholds
#' @noRd
flexible_score_cdf <- function(draws, center, scale, q) {
  dims <- dim(draws)
  if (length(dims) != 3L || !identical(dim(center), c(dims[1L], dims[3L])) ||
      !identical(dim(scale), dim(center)) || any(!is.finite(c(draws, center, scale))) ||
      any(scale <= 0) || !is.numeric(q) || is.complex(q) || anyNA(q) || any(q < 0) ||
      is.unsorted(q, strictly = TRUE)) stop("Invalid joint draws or score specification.", call. = FALSE)
  n <- dims[1L]; b <- dims[2L]; d <- dims[3L]
  result <- matrix(0, n, length(q))
  for (i in seq_len(n)) {
    scores <- numeric(b)
    for (j in seq_len(d)) scores <- pmax(scores, abs((draws[i, , j] - center[i, j]) / scale[i, j]))
    result[i, ] <- findInterval(q, sort(scores)) / b
  }
  result
}

#' Compile bounded additive tables for the finite tree ensemble
#'
#' Tables use raw numeric split bins. Regrouping nonnegative rates changes
#' arithmetic association only; trees, split cuts, leaf rates and their mean
#' remain unchanged. Tables share a bounded budget across fitted transitions.
#' @noRd
compile_forest_additive <- function(compiled, max_bytes = 16 * 1024^2,
                                    max_cells = 1048576L) {
  lookup <- compiled$lookup
  if (is.null(lookup) || !isTRUE(lookup$complete)) return(NULL)
  trees <- length(compiled$roots)
  first <- lookup$feature1; second <- lookup$feature2
  constants <- compiled$feature[compiled$roots] == 0L
  first[constants] <- second[constants] <- 0L
  keys <- paste(first, second, sep = ':')
  order <- unique(keys)
  groups <- match(keys, order)
  starts <- match(order, keys)
  f1 <- first[starts]; f2 <- second[starts]
  counts <- lookup$cut_count
  size1 <- size2 <- rep.int(1L, length(order))
  size1[f1 > 0L] <- counts[f1[f1 > 0L]] + 1L
  size2[f2 > 0L] <- counts[f2[f2 > 0L]] + 1L
  cells <- as.double(size1) * size2
  if (any(cells > max_cells) || 4096 + 8 * sum(cells) + 24 * length(order) > max_bytes)
    return(NULL)
  offsets <- c(0, cumsum(cells))[-(length(cells) + 1L)]
  table <- numeric(sum(cells))
  for (b in seq_len(trees)) {
    group <- groups[b]; n1 <- size1[group]; n2 <- size2[group]
    where <- offsets[group] + seq_len(cells[group])
    if (constants[b]) {
      table[where] <- table[where] + compiled$value[compiled$roots[b]]
    } else {
      a <- lookup$map[lookup$map_start1[b] + seq_len(n1)]
      d <- if (second[b]) lookup$map[lookup$map_start2[b] + seq_len(n2)] else 0L
      at <- lookup$offset[b] + as.vector(outer(a, lookup$size1[b] * d, '+')) + 1L
      table[where] <- table[where] + compiled$value[lookup$leaf[at]]
    }
  }
  if (any(!is.finite(table))) return(NULL)
  answer <- list(feature1 = as.integer(f1), feature2 = as.integer(f2),
    size1 = as.integer(size1), size2 = as.integer(size2),
    offset = as.integer(offsets), value = table, trees = as.integer(trees))
  if (as.numeric(utils::object.size(answer)) > max_bytes) return(NULL)
  answer
}
