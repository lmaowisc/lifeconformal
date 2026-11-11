#' Match forest continuation sums to the running R arithmetic
#'
#' The native loop uses the same extended accumulator as rowSums, followed by
#' one double conversion. Unsupported compiler/runtime formats retain the R
#' loop; no fitted-model or numerical control is changed.
#' @param native,machine,wide Native and runtime floating-point format details.
#' @returns Matching mantissa precision, or NULL for the ordinary R loop.
#' @noRd
flex_forest_block_precision <- function(
    native = .Call(life_forest_block_info), machine = .Machine,
    wide = isTRUE(capabilities("long.double"))) {
  digits <- if (wide) machine$longdouble.digits else machine$double.digits
  exponent <- if (wide) machine$longdouble.max.exp else machine$double.max.exp
  if (length(digits) != 1L || length(exponent) != 1L) return(NULL)
  at <- if (identical(as.integer(digits), native[1L])) c(1L, 3L) else c(2L, 4L)
  if (!identical(as.integer(c(digits, exponent)), native[at])) return(NULL)
  as.integer(digits)
}

#' Prepare private path state for one continuation uniform block
#'
#' The native call duplicates these vectors and retains no pointer afterward.
#' Original, converted and returned state are included in its workspace bound.
#' Baseline/model arrays remain read-only inputs and are not copied here.
#' @param tt,state,entry,cell,cell_end,previous,accrued,count,iteration Current
#'   private simulation arrays and absolute iteration number.
#' @returns A typed state list for the native block evaluator.
#' @noRd
flex_forest_block_state <- function(tt, state, entry, cell, cell_end, previous,
                                    accrued, count, iteration) {
  storage.mode(accrued) <- "double"
  list(tt = as.double(tt), state = as.integer(state), entry = as.double(entry),
    cell = as.double(cell), cell_end = as.double(cell_end), previous = as.double(previous),
    accrued = accrued, count = as.integer(count), iteration = as.integer(iteration))
}

#' Check generated-bank resources before drawing any uniforms
#'
#' The native bound includes the existing continuation workspace, old/new bank
#' overlap, stored group identities and metadata. Unsupported sizes retain the
#' ordinary matrix route without changing batch size, draws or seeds.
#' @param adapter Frozen flexible outcome adapter.
#' @param group Full B-row groups of numerical profile indices.
#' @param B Number of draws per profile.
#' @param precision Validated native accumulator precision.
#' @param budget Maximum call-local live workspace in bytes, not total RSS.
#' @returns Whether a generated bank fits the supported format and budget.
#' @noRd
flex_forest_bank_supported <- function(adapter, group, B, precision,
                                       budget = 64 * 1024^2) {
  .Call(life_forest_bank_supported, adapter, as.integer(group), as.integer(B),
    precision, as.double(budget), 0L)
}

#' Generate an owned uniform bank in the original profile seed scopes
#'
#' Every profile is filled on refresh, including profiles whose paths finished
#' in an earlier block. The native consumer releases the sealed token on return
#' or error; this helper releases incomplete generation. Tokens never enter a
#' fit, saved object or public result.
#' @param adapter,group,B,precision,budget As in flex_forest_bank_supported.
#' @param profile_seeds Existing numerical-profile seeds, in batch order.
#' @param block Absolute 32-step bank index.
#' @returns A complete owned bank token for one native advancement call.
#' @noRd
flex_forest_bank_generate <- function(adapter, group, B, precision, profile_seeds,
                                      block, budget = 64 * 1024^2) {
  token <- .Call(life_forest_bank_create, adapter, as.integer(group), as.integer(B),
    precision, as.integer(block), as.double(budget), 0L)
  complete <- FALSE
  on.exit(if (!complete) .Call(life_forest_bank_release, token), add = TRUE)
  for (ii in seq_along(profile_seeds))
    flex_with_seed(flex_profile_seed(profile_seeds[ii], block),
      .Call(life_forest_bank_fill, token, as.integer(ii)))
  .Call(life_forest_bank_seal, token)
  complete <- TRUE
  token
}
