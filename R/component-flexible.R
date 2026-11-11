#' Match native calibration accumulation to the running R arithmetic
#' @param native,machine,wide Native and R floating-point formats. Defaults use
#'   the current runtime; explicit formats support portable contract checks.
#' @returns The matching mantissa precision, or NULL to retain the ordinary
#'   calibration path on a platform with a different native number format.
#' @noRd
calibration_sum_precision <- function(
    native = .Call("life_calibration_sum_info", PACKAGE = "lifeconformal"),
    machine = .Machine, wide = isTRUE(capabilities("long.double"))) {
  digits <- if (wide) machine$longdouble.digits else machine$double.digits
  exponent <- if (wide) machine$longdouble.max.exp else machine$double.max.exp
  if (length(digits) != 1L || length(exponent) != 1L) return(NULL)
  # Check the branch selected by the native constructor. Equal mantissas do
  # not imply equal exponent ranges on every conforming floating-point ABI.
  at <- if (digits == native[1L]) c(1L, 3L) else c(2L, 4L)
  if (!identical(as.integer(c(digits, exponent)), native[at])) return(NULL)
  as.integer(digits)
}

#' Share flexible continuations across marginal calibration curves
#'
#' @param object Full trained fit, retaining the augmentation model's targets.
#' @param context Shared history, nuisance moments and signed calibration queries.
#' @param scalar_objects,grids Temporary marginal fits and their fixed grids.
#' @param active Components whose allocated finite-sample ranks are attainable.
#' @param controls Resolved calibration controls.
#' @param precision Native accumulation precision validated against running R.
#' @param query_batch Maximum queries per continuation bank.
#' @param bank_budget Budget used to reduce batch size for large draw or grid
#'   dimensions. One query's required outputs are the minimum working size.
#' @returns Per-component contributions when requested, means and raw sums,
#'   structural determination arrays, and the ordinary augmentation diagnostics.
#' @details Subjects and negative query rows retain their original order. Each
#'   double weighted product enters the same extended sum used by R colSums.
#'   That sum is rounded to double before the endpoint term is added. Whole
#'   curve sums and means use separate final rounding, preserving integer-rank
#'   inversion. Simulation banks and probability scratch are discarded after
#'   each chunk; full contribution matrices are retained only when requested.
#'   Existing dense structural determination arrays remain component-specific.
#' @noRd
evaluate_flexible_component_curves <- function(object, context, scalar_objects,
    grids, active, controls, precision, query_batch = 16L,
    bank_budget = 32 * 1024^2) {
  states <- list()
  on.exit(for (state in states)
    .Call("life_calibration_sum_close", state, PACKAGE = "lifeconformal"), add = TRUE)
  create <- function(n) {
    state <- .Call("life_calibration_sum_create", as.integer(n), precision,
      PACKAGE = "lifeconformal")
    # A constructor interrupted before this assignment retains its native
    # finalizer; every state received here is released immediately on exit.
    states[[length(states) + 1L]] <<- state
    state
  }
  add <- function(state, value, weight)
    .Call("life_calibration_sum_add", state, value, as.double(weight),
      PACKAGE = "lifeconformal")
  total <- function(state, divisor = 1)
    .Call("life_calibration_sum_value", state, as.double(divisor),
      PACKAGE = "lifeconformal")
  reset <- function(state)
    .Call("life_calibration_sum_reset", state, PACKAGE = "lifeconformal")
  history <- context$history
  n <- nrow(history$subjects); d <- length(scalar_objects)
  target_names <- colnames(object$problem$rewards)
  contributions <- if (controls$keep_contributions) lapply(seq_len(d), function(j)
    if (j %in% active) matrix(0, n, length(grids[[j]])) else NULL) else NULL
  estimates <- sums <- lapply(seq_len(d), function(j) numeric())
  diagnostics <- vector("list", d)
  if (length(active)) {
    columns <- match(target_names, colnames(context$moments$mean))
    cdt <- lapply(seq_len(d), function(j) if (j %in% active)
      compute_cdt(history, scalar_objects[[j]]$problem,
        context$moments$mean[, columns[j], drop = FALSE],
        context$moments$scale[, columns[j], drop = FALSE], grids[[j]]) else NULL)
    queries <- context$queries; prepared <- queries$prepared
    selected <- match(target_names, colnames(object$augmentation$problem$rewards))
    stopifnot(!anyNA(selected))
    if (!identical(target_names, colnames(object$augmentation$problem$rewards)))
      prepared <- prepare_occupation_history(history, object$augmentation$problem)
    curve_sum <- query_sum <- vector("list", d)
    for (j in active) {
      curve_sum[[j]] <- create(length(grids[[j]]))
      query_sum[[j]] <- create(length(grids[[j]]))
      diagnostics[[j]] <- c(list(conditional_coverage_evaluated = FALSE,
        max_probability_bound_width = 0, queried_snapshots = 0L,
        structural_overrides = 0L), queries$diagnostics)
    }
    # Bound bank output and query-by-grid scratch independently of total H.
    # One query's required outputs remain the lower bound on working memory.
    per_query <- 8 * controls$draws * ncol(object$augmentation$problem$rewards) +
      32 * max(lengths(grids[active]))
    batch <- as.integer(max(1, min(query_batch, floor(bank_budget / per_query))))
    for (i in seq_len(n)) {
      first <- if (i == 1L) 1L else queries$endpoint_index[i - 1L] + 1L
      rows <- seq.int(first, queries$endpoint_index[i]); h <- length(rows)
      times <- queries$time[rows]; weight <- queries$weight[rows]
      last <- vapply(active, function(j) if (all(cdt[[j]]$resolved[i, ]))
        max(cdt[[j]]$time[i, ]) else Inf, numeric(1))
      fast <- active[times[1L] >= last]
      working <- setdiff(active, fast)
      endpoints <- vector("list", d)
      for (j in fast)
        diagnostics[[j]]$structural_overrides <- diagnostics[[j]]$structural_overrides +
          h * length(grids[[j]])
      for (j in working) reset(query_sum[[j]])
      if (length(working)) for (begin in seq.int(1L, h, by = batch)) {
        chunk <- seq.int(begin, min(h, begin + batch - 1L))
        needed <- chunk[times[chunk] < max(last)]
        draws <- NULL
        if (length(needed)) {
          snapshot <- evaluate_occupation_history(prepared, rep(i, length(needed)),
            times[needed], side = "right")
          snaps <- flexible_history_snapshots(prepared, i, snapshot, seq_along(needed))
          draws <- simulate_flexible_targets(object$augmentation,
            context$outcome_view$x[rep(i, length(needed)), , drop = FALSE],
            snaps, B = controls$draws)
        }
        for (j in working) {
          m <- length(chunk); L <- length(grids[[j]])
          certified <- outer(times[chunk], cdt[[j]]$time[i, ], `>=`) &
            matrix(rep(cdt[[j]]$resolved[i, ], each = m), m, L)
          values <- matrix(rep(as.numeric(cdt[[j]]$membership[i, ]), each = m), m, L)
          need <- which(rowSums(!certified) > 0L)
          if (length(need)) {
            probabilities <- flexible_score_cdf(
              draws[match(chunk[need], needed), , selected[j], drop = FALSE],
              matrix(context$moments$mean[i, columns[j]], length(need), 1L),
              matrix(context$moments$scale[i, columns[j]], length(need), 1L), grids[[j]])
            known <- certified[need, , drop = FALSE]
            probabilities[known] <- values[need, , drop = FALSE][known]
            values[need, ] <- probabilities
            diagnostics[[j]]$queried_snapshots <- diagnostics[[j]]$queried_snapshots + length(need)
          }
          diagnostics[[j]]$structural_overrides <- diagnostics[[j]]$structural_overrides +
            as.double(sum(certified))
          negative <- chunk < h
          if (any(negative)) add(query_sum[[j]], values[negative, , drop = FALSE], weight[chunk[negative]])
          if (any(!negative)) endpoints[[j]] <- values[!negative, ] * weight[h]
        }
      }
      for (j in active) {
        value <- if (j %in% fast) as.numeric(cdt[[j]]$membership[i, ]) else
          endpoints[[j]] + total(query_sum[[j]])
        origin <- cdt[[j]]$resolved[i, ] & cdt[[j]]$time[i, ] == 0
        value[origin] <- as.numeric(cdt[[j]]$membership[i, origin])
        if (controls$keep_contributions) contributions[[j]][i, ] <- value
        add(curve_sum[[j]], matrix(value, 1L), 1)
      }
    }
    for (j in active) {
      sums[[j]] <- total(curve_sum[[j]])
      estimates[[j]] <- total(curve_sum[[j]], n)
      diagnostics[[j]]$conditional_coverage_evaluated <- diagnostics[[j]]$queried_snapshots > 0L
      if (diagnostics[[j]]$conditional_coverage_evaluated)
        diagnostics[[j]]$max_probability_bound_width <- NA_real_
    }
  }
  lapply(seq_len(d), function(j) if (j %in% active)
    list(contributions = if (controls$keep_contributions) contributions[[j]] else NULL,
      estimate = estimates[[j]], estimate_sum = sums[[j]], cdt = cdt[[j]],
      diagnostics = diagnostics[[j]]) else NULL)
}
