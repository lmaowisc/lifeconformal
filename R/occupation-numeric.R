#' Complete binary occupations at aligned history query times
#'
#' @param end,start Finite query and preceding event times, one per matrix row.
#' @param reward,accrued Aligned query-by-target binary rewards and accrued
#'   occupations. The returned matrix preserves accrued's dimensions and names.
#' @details The binary total is rounded after adding the accrued occupation,
#'   avoiding cancellation from a separately rounded elapsed interval.
#' @noRd
occupation_complete <- function(end, start, reward, accrued) {
  storage.mode(reward) <- "double"
  storage.mode(accrued) <- "double"
  .Call("life_occupation_complete", as.double(end), as.double(start), reward,
        accrued, PACKAGE = "lifeconformal")
}
