#' Construct a stochastic kernel from transition-intensity jumps
#'
#' Off-diagonal intensity increments
#' give transition probabilities. If their sum in a row exceeds one, divide
#' that row by its sum and record the adjustment. The diagonal is the remaining
#' probability. Structural zeros and absorbing states are checked explicitly.
#' If row normalization would turn a positive transition into zero by
#' underflow, the calculation stops with an explicit numerical diagnostic.
#'
#' @param increments Numeric square matrix of nonnegative off-diagonal
#'   intensity increments, with a zero diagonal.
#' @param absorbing Integer indices of declared absorbing states.
#' @param allowed Optional logical matrix marking allowed off-diagonal
#'   transitions. A positive increment at a prohibited entry is an error.
#' @returns A list with `kernel`, `adjusted_rows`, `original_row_sums`, and
#'   `row_scale` (one unless a row was adjusted).
#' @noRd
make_markov_kernel <- function(increments, absorbing = integer(0),
                               allowed = NULL) {
  if (!is.matrix(increments) || !is.numeric(increments) || is.complex(increments) ||
      nrow(increments) < 1L || nrow(increments) != ncol(increments) ||
      any(!is.finite(increments)) || any(increments < 0) ||
      any(diag(increments) != 0)) {
    stop("`increments` must be a finite, nonnegative square matrix with a zero diagonal.",
         call. = FALSE)
  }
  k <- nrow(increments)
  check_markov_state_indices(absorbing, k, "absorbing", allow_empty = TRUE)
  if (!is.null(allowed)) {
    if (!is.matrix(allowed) || !is.logical(allowed) || anyNA(allowed) ||
        !identical(dim(allowed), dim(increments))) {
      stop("`allowed` must be a logical matrix matching `increments`.",
           call. = FALSE)
    }
    if (any(increments[!allowed] > 0)) {
      stop("Positive intensity increments occur at prohibited transitions.",
           call. = FALSE)
    }
  }
  original_sum <- rowSums(increments)
  if (any(!is.finite(original_sum))) {
    stop("Intensity row sums overflow; rescale the fitted input before use.",
         call. = FALSE)
  }
  if (length(absorbing) && any(original_sum[absorbing] > 0)) {
    stop("Declared absorbing states have outgoing intensity increments.",
         call. = FALSE)
  }
  adjusted <- which(original_sum > 1)
  row_scale <- rep(1, k)
  row_scale[adjusted] <- 1 / original_sum[adjusted]
  kernel <- increments * row_scale
  if (any(increments > 0 & kernel == 0))
    stop("Positive Markov transition probabilities underflow during row normalization; check the fitted intensity increments.",
         call. = FALSE)
  diag(kernel) <- pmax(0, 1 - rowSums(kernel))
  list(kernel = kernel, adjusted_rows = adjusted,
       original_row_sums = original_sum, row_scale = row_scale)
}

#' Compute component moments of future Markov occupation rewards
#'
#' Uses deterministic backward mixing at fitted jump times and reward shifts
#' between jumps. A stable total-variance
#' recursion avoids subtracting two large, almost equal raw moments. Only
#' component variances are computed; joint dependence is retained by
#' [propagate_markov_rewards()] when a distribution is needed.
#'
#' @param jump_times Strictly increasing positive jump times through `tau`.
#' @param kernels List of square stochastic matrices in jump-time order.
#'   Unnamed matrices use the state order in `rewards`. If a kernel is named,
#'   both dimensions must use the same state labels in that order; named kernels
#'   must agree with one another and with named rows of `rewards`.
#' @param rewards Numeric state-by-target matrix of rewards per time unit.
#' @param tau Finite positive prediction horizon.
#' @param start Finite query time in `[0, tau]`.
#' @param side Query history side. `"right"` conditions on the state through
#'   `start`, excluding a kernel at that time. `"left"` includes that kernel.
#' @param accrued Known accumulated target vector, added after propagation.
#'   Defaults to a zero vector of length `ncol(rewards)`.
#' @returns A list with state-by-target `mean`, `variance`, `second_moment`,
#'   and numerical `diagnostics`. States refer to the selected query side.
#' @details Moment arrays use rolling storage of `O(K d)`. This reference
#'   engine validates and scans dense kernels in `O(J K^2)` work. Additional
#'   mixing work uses positive entries and is proportional to the number of
#'   nonzero entries times the number of targets. Kernel validation may allocate
#'   normalized copies, using up to `O(J K^2)` additional storage. Input matrices whose row sums
#'   differ from one by at most `1e-12` are normalized to remove roundoff;
#'   diagnostics retain the maximum such deviation. Other invalid kernels fail.
#' @noRd
markov_reward_moments <- function(jump_times = numeric(), kernels = list(),
                                  rewards, tau, start = 0,
                                  side = c("right", "left"), accrued = NULL) {
  spec <- prepare_markov_recursion(jump_times, kernels, rewards, tau, start,
                                   side, accrued)
  k <- nrow(rewards)
  d <- ncol(rewards)
  n_jump <- length(spec$jump_times)
  last_time <- if (n_jump) spec$jump_times[n_jump] else start
  means <- (tau - last_time) * rewards
  variances <- matrix(0, k, d, dimnames = dimnames(rewards))
  constant <- matrix(TRUE, k, d)

  if (n_jump) {
    for (j in rev(seq_len(n_jump))) {
      mixed_mean <- matrix(0, k, d, dimnames = dimnames(rewards))
      mixed_var <- mixed_mean
      mixed_constant <- matrix(FALSE, k, d)
      kernel <- spec$kernels[[j]]
      for (state in seq_len(k)) {
        destinations <- which(kernel[state, ] > 0)
        weights <- kernel[state, destinations]
        next_mean <- means[destinations, , drop = FALSE]
        known <- if (spec$jump_times[j] == tau) rep(TRUE, d) else
          vapply(seq_len(d), function(target)
            all(constant[destinations, target] &
                  rewards[destinations, target] == rewards[state, target]), logical(1))
        mixed_constant[state, ] <- known
        center <- colSums(next_mean * weights)
        # Equal-reward positive paths have a known deterministic component.
        # Weighted sums of identical doubles must not create fake variance.
        center[known] <- (tau - spec$jump_times[j]) * rewards[state, known]
        centered <- sweep(next_mean, 2L, center, "-")
        mixed_mean[state, ] <- center
        mixed_var[state, ] <- colSums(
          (variances[destinations, , drop = FALSE] + centered^2) * weights
        )
        mixed_var[state, known] <- 0
        if (any(mixed_var[state, ] < .Machine$double.xmin)) {
          lower <- upper <- next_mean[1L, ]
          positive_variance <- variances[destinations[1L], ] > 0
          if (length(destinations) > 1L) for (a in seq.int(2L, length(destinations))) {
            lower <- pmin(lower, next_mean[a, ])
            upper <- pmax(upper, next_mean[a, ])
            positive_variance <- positive_variance | variances[destinations[a], ] > 0
          }
          check_markov_variance_scale(mixed_var[state, ],
                                      (upper > lower | positive_variance) & !known)
        }
      }
      previous_time <- if (j > 1L) spec$jump_times[j - 1L] else start
      means <- mixed_mean + (spec$jump_times[j] - previous_time) * rewards
      variances <- mixed_var
      constant <- mixed_constant
      canonical <- (tau - previous_time) * rewards
      means[constant] <- canonical[constant]
      if (any(!is.finite(means)) || any(!is.finite(variances))) {
        stop("Markov reward moments overflow at a backward step.", call. = FALSE)
      }
    }
  }
  means <- sweep(means, 2L, spec$accrued, "+")
  for (state in seq_len(k)) if (any(constant[state, ])) {
    known <- constant[state, ]
    means[state, known] <- markov_protected_total(tau, start,
      rewards[state, known], spec$accrued[known])
  }
  second <- variances + means^2
  if (any(!is.finite(means)) || any(!is.finite(second))) {
    stop("Markov reward moments overflow; check reward and time units.",
         call. = FALSE)
  }
  list(mean = means, variance = variances, second_moment = second,
       diagnostics = list(method = "deterministic_backward_moments",
                          jump_count = n_jump, start = start, side = spec$side,
                          kernel_row_error = spec$kernel_row_error))
}

#' Reject unrepresentable stochastic variances without rejecting point masses
#' @noRd
check_markov_variance_scale <- function(variance, stochastic) {
  if (any(stochastic & variance < .Machine$double.xmin)) {
    stop("Markov reward variance underflows at the current time unit; rescale time before fitting.",
         call. = FALSE)
  }
  invisible(NULL)
}

#' Detect which target rewards can change along the fitted working graph
#' @noRd
markov_reward_variability <- function(rewards, from, to) {
  varying <- matrix(FALSE, nrow(rewards), ncol(rewards))
  repeat {
    previous <- varying
    for (e in seq_along(from)) {
      varying[from[e], ] <- varying[from[e], ] | varying[to[e], ] |
        rewards[from[e], ] != rewards[to[e], ]
    }
    if (identical(varying, previous)) return(varying)
  }
}

#' Propagate an exact finite Markov reward distribution
#'
#' Backward kernel mixtures and deterministic reward translations preserve
#' the joint target vector,
#' including atoms from no further transition and from absorbing trajectories.
#' Exactly equal floating-point target vectors are merged; nearby vectors are
#' kept distinct. No nonzero probability is pruned, and no simulation is used.
#'
#' @inheritParams markov_reward_moments
#' @param state Optional single integer state index. `NULL` returns a
#'   distribution for every possible query state.
#' @param max_atoms Positive integer budget for the sum of atom counts across
#'   states in each rolling time layer. Exceeding it raises a
#'   `lifeconformal_atom_budget` error rather than changing algorithms silently.
#' @returns A list with `distributions`, `state`, and `diagnostics`. Each
#'   distribution has an atom-by-target matrix `values` and a probability
#'   vector `probability`. Distributions and matrix columns retain input names.
#' @details This routine is an exact small-problem engine. It is not the
#'   scalable accrued-coordinate-grid approximation planned for large problems.
#'   Work is proportional to the atoms traversed across positive kernel edges,
#'   including canonical-key construction; exact atom counts can grow rapidly.
#'   Two rolling layers and per-state lookup tables bound storage in terms of
#'   `max_atoms`. The maximum uncorrected output mass error is returned.
#' @noRd
propagate_markov_rewards <- function(jump_times = numeric(), kernels = list(),
                                     rewards, tau, start = 0,
                                     side = c("right", "left"), accrued = NULL,
                                     state = NULL, max_atoms = 100000L) {
  spec <- prepare_markov_recursion(jump_times, kernels, rewards, tau, start,
                                   side, accrued)
  k <- nrow(rewards)
  d <- ncol(rewards)
  if (!is.null(state)) {
    check_markov_state_indices(state, k, "state", allow_empty = FALSE)
    if (length(state) != 1L) {
      stop("`state` must be a single state index or NULL.", call. = FALSE)
    }
  }
  if (!is.numeric(max_atoms) || is.complex(max_atoms) ||
      !is.null(dim(max_atoms)) || length(max_atoms) != 1L ||
      !is.finite(max_atoms) || max_atoms < 1 ||
      max_atoms != floor(max_atoms) || max_atoms > .Machine$integer.max) {
    stop("`max_atoms` must be a positive integer budget.", call. = FALSE)
  }
  if (max_atoms < k) stop_markov_atom_budget(max_atoms)
  n_jump <- length(spec$jump_times)
  last_time <- if (n_jump) spec$jump_times[n_jump] else start
  distributions <- lapply(seq_len(k), function(s) {
    list(values = matrix((tau - last_time) * rewards[s, ], nrow = 1L,
                         dimnames = list(NULL, colnames(rewards))),
         probability = 1, constant = matrix(TRUE, 1L, d))
  })
  peak_atoms <- k
  traversed_atoms <- 0
  max_mass_error <- 0

  if (n_jump) {
    for (j in rev(seq_len(n_jump))) {
      previous_time <- if (j > 1L) spec$jump_times[j - 1L] else start
      duration <- spec$jump_times[j] - previous_time
      next_layer <- vector("list", k)
      layer_atoms <- 0L
      kernel <- spec$kernels[[j]]
      for (s in seq_len(k)) {
        lookup <- new.env(hash = TRUE, parent = emptyenv())
        capacity <- min(1024L, max_atoms - layer_atoms)
        value_rows <- vector("list", capacity)
        constant_rows <- vector("list", capacity)
        probability <- numeric(capacity)
        count <- 0L
        shift <- duration * rewards[s, ]
        for (destination in which(kernel[s, ] > 0)) {
          source <- distributions[[destination]]
          for (atom in seq_along(source$probability)) {
            value <- source$values[atom, ] + shift
            # A destination entered at tau earns no occupation. This numerical
            # proof is about positive-duration rewards, not a departure flag.
            known <- source$constant[atom, ] &
              (spec$jump_times[j] == tau | rewards[destination, ] == rewards[s, ])
            value[known] <- (tau - previous_time) * rewards[s, known]
            if (any(!is.finite(value))) {
              stop("Markov reward atoms overflow; check reward and time units.",
                   call. = FALSE)
            }
            value[value == 0] <- 0  # Canonicalize signed zero only.
            mass <- kernel[s, destination] * source$probability[atom]
            if (mass == 0) {
              stop("A positive Markov atom probability underflowed to zero.",
                   call. = FALSE)
            }
            # Keep provenance until accrued totals are evaluated. Equal
            # rounded coordinates need not have the same no-change proof.
            key <- paste(paste(as.integer(known), collapse = ""),
                         paste(sprintf("%a", value), collapse = ":"), sep = "/")
            index <- lookup[[key]]
            if (is.null(index)) {
              if (layer_atoms + count >= max_atoms) {
                stop_markov_atom_budget(max_atoms)
              }
              count <- count + 1L
              if (count > capacity) {
                capacity <- min(max_atoms - layer_atoms,
                                max(count, 2 * capacity))
                length(value_rows) <- capacity
                length(constant_rows) <- capacity
                length(probability) <- capacity
              }
              value_rows[[count]] <- value
              constant_rows[[count]] <- known
              probability[count] <- mass
              lookup[[key]] <- count
            } else {
              probability[index] <- probability[index] + mass
            }
            traversed_atoms <- traversed_atoms + 1
          }
        }
        values <- do.call(rbind, value_rows[seq_len(count)])
        dimnames(values) <- list(NULL, colnames(rewards))
        probability <- probability[seq_len(count)]
        mass_error <- abs(sum(probability) - 1)
        max_mass_error <- max(max_mass_error, mass_error)
        if (!is.finite(mass_error) || mass_error > 1e-10) {
          stop("Markov atom probabilities lost unit mass during propagation.",
               call. = FALSE)
        }
        next_layer[[s]] <- list(values = values, probability = probability,
                               constant = do.call(rbind, constant_rows[seq_len(count)]))
        layer_atoms <- layer_atoms + count
      }
      distributions <- next_layer
      peak_atoms <- max(peak_atoms, layer_atoms)
    }
  }
  selected_states <- if (is.null(state)) seq_len(k) else as.integer(state)
  distributions <- distributions[selected_states]
  distributions <- lapply(seq_along(distributions), function(a) {
    distribution <- distributions[[a]]
    distribution$values <- sweep(distribution$values, 2L, spec$accrued, "+")
    canonical <- markov_protected_total(tau, start, rewards[selected_states[a], ],
                                        spec$accrued)
    for (target in seq_len(d)) distribution$values[distribution$constant[, target], target] <-
      canonical[target]
    if (any(!is.finite(distribution$values))) {
      stop("Accrued Markov reward atoms overflow.", call. = FALSE)
    }
    # Restore the public atom representation after the internal provenance
    # split. Only exactly equal final binary64 coordinates are combined.
    keys <- apply(distribution$values, 1L, function(value)
      paste(sprintf("%a", value), collapse = ":"))
    group <- match(keys, unique(keys))
    first <- !duplicated(group)
    mass <- as.vector(rowsum(distribution$probability, group, reorder = FALSE))
    list(values = distribution$values[first, , drop = FALSE], probability = mass)
  })
  names(distributions) <- rownames(rewards)[selected_states]
  list(distributions = distributions, state = selected_states,
       diagnostics = list(method = "deterministic_exact_atoms",
                          jump_count = n_jump, start = start, side = spec$side,
                          peak_atoms = peak_atoms,
                          traversed_atoms = traversed_atoms,
                          max_atoms = max_atoms,
                          mass_error = max_mass_error,
                          kernel_row_error = spec$kernel_row_error))
}

#' Validate and slice a fitted Markov jump representation
#' @inheritParams markov_reward_moments
#' @returns Validated kernels and query-specific numerical settings.
#' @noRd
prepare_markov_recursion <- function(jump_times, kernels, rewards, tau, start,
                                     side, accrued) {
  if (!is.matrix(rewards) || !is.numeric(rewards) || is.complex(rewards) ||
      any(dim(rewards) < 1L) || any(!is.finite(rewards))) {
    stop("`rewards` must be a finite numeric state-by-target matrix.",
         call. = FALSE)
  }
  if (!is.numeric(tau) || is.complex(tau) || !is.null(dim(tau)) || length(tau) != 1L ||
      !is.finite(tau) || tau <= 0 ||
      !is.numeric(start) || is.complex(start) || !is.null(dim(start)) ||
      length(start) != 1L || !is.finite(start) ||
      start < 0 || start > tau) {
    stop("`tau` must be positive and finite, and `start` must lie in [0, tau].",
         call. = FALSE)
  }
  if (!is.numeric(jump_times) || is.complex(jump_times) ||
      !is.null(dim(jump_times)) || any(!is.finite(jump_times)) ||
      any(jump_times <= 0) || any(jump_times > tau) ||
      any(diff(jump_times) <= 0)) {
    stop("`jump_times` must be strictly increasing, positive and no greater than tau.",
         call. = FALSE)
  }
  if (!is.list(kernels) || !is.null(dim(kernels)) || length(kernels) != length(jump_times)) {
    stop("`kernels` must contain one matrix per jump time.", call. = FALSE)
  }
  side <- match.arg(side, c("right", "left"))
  k <- nrow(rewards)
  state_names <- rownames(rewards)
  if (!is.null(state_names) &&
      (anyNA(state_names) || any(!nzchar(state_names)) || anyDuplicated(state_names))) {
    stop("Named reward rows must use distinct, nonempty state labels.", call. = FALSE)
  }
  max_row_error <- 0
  for (j in seq_along(kernels)) {
    kernel <- kernels[[j]]
    if (!is.matrix(kernel) || !is.numeric(kernel) || is.complex(kernel) ||
        !identical(dim(kernel), c(k, k)) || any(!is.finite(kernel)) ||
        any(kernel < 0)) {
      stop("Each Markov kernel must be a finite nonnegative state-by-state matrix.",
           call. = FALSE)
    }
    if (!is.null(dimnames(kernel))) {
      labels <- rownames(kernel)
      if (is.null(labels) || !identical(labels, colnames(kernel)) ||
          anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels)) {
        stop("Named kernels must use identical, distinct state labels on both dimensions.",
             call. = FALSE)
      }
      if (!is.null(state_names) && !identical(labels, state_names)) {
        stop("Kernel state labels must match the reward and other kernel state order.",
             call. = FALSE)
      }
      state_names <- labels
    }
    row_sums <- rowSums(kernel)
    row_error <- max(abs(row_sums - 1))
    if (row_error > 1e-12) {
      stop("Each Markov kernel row must sum to one (tolerance 1e-12).",
           call. = FALSE)
    }
    max_row_error <- max(max_row_error, row_error)
    if (row_error > 0) kernels[[j]] <- kernel / row_sums
  }
  if (is.null(accrued)) accrued <- numeric(ncol(rewards))
  if (!is.numeric(accrued) || is.complex(accrued) ||
      !is.null(dim(accrued)) || length(accrued) != ncol(rewards) ||
      any(!is.finite(accrued))) {
    stop("`accrued` must be a finite numeric vector with one value per target.",
         call. = FALSE)
  }
  keep <- if (side == "right") jump_times > start else jump_times >= start
  list(jump_times = jump_times[keep], kernels = kernels[keep], accrued = accrued,
       side = side, kernel_row_error = max_row_error)
}

#' Validate internal one-based state indices
#' @param x Numeric indices.
#' @param k Number of states.
#' @param name Argument name used in an error message.
#' @param allow_empty Whether a zero-length vector is accepted.
#' @returns Invisibly `NULL` on success.
#' @noRd
check_markov_state_indices <- function(x, k, name, allow_empty) {
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) ||
      (!allow_empty && !length(x)) || any(!is.finite(x)) ||
      any(x != floor(x)) || any(x < 1 | x > k) || anyDuplicated(x)) {
    stop(sprintf("`%s` must contain distinct integer state indices from 1 to %d.",
                 name, k), call. = FALSE)
  }
  invisible(NULL)
}

#' Stop an exact atom calculation at its declared resource budget
#' @param max_atoms Atom budget.
#' @returns Does not return; raises a typed error.
#' @noRd
stop_markov_atom_budget <- function(max_atoms) {
  message <- paste0("Exact Markov propagation exceeds `max_atoms = ",
                    max_atoms, "`. Increase this explicit budget or use a ",
                    "validated deterministic grid engine when available.")
  stop(structure(list(message = message, call = NULL),
                 class = c("lifeconformal_atom_budget", "error", "condition")))
}
