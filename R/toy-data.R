#' Generate synthetic histories for learning the package workflow
#'
#' Create a small, reproducible event-status data set with one baseline age
#' and one possible transition, from state zero to absorbing state two.
#' These are artificial observations for examples, not clinical data or a
#' simulation study establishing statistical performance.
#'
#' @param n Positive integer number of independent synthetic subjects.
#' @param seed Optional nonnegative integer random seed, default 1. An explicit
#'   seed preserves the caller's random-number state. NULL uses and advances
#'   the caller's random-number stream.
#'
#' @returns A data frame with one row per subject and four columns:
#'   \describe{
#'     \item{id}{Integer subject identifier, from 1 to n.}
#'     \item{time}{Observed follow-up time, in arbitrary time units.}
#'     \item{status}{Two for entry into the absorbing state, zero for censoring
#'       or the administrative end of follow-up.}
#'     \item{age}{Synthetic baseline age, uniformly distributed from 40 to 80.}
#'   }
#' @details
#' All subjects start in state zero. Conditional on age, the transition time
#' and censoring time are independent exponentials with rates
#' `0.14 * exp(0.018 * (age - 60))` and
#' `0.04 * exp(0.008 * (age - 60))`, respectively. Administrative follow-up ends
#' at time eight. For example, use `tau = 5`, `absorbing = 2`, and
#' `targets = list(alive = 0)` to predict time in state zero through time five.
#'
#' Split the generated subjects into separate training and calibration samples.
#' Repeated calls restart the identifier sequence; separately generated samples
#' require identifiers that distinguish their subjects. This helper simulates
#' example data only. Markov-Cox predictions in [life_train()] and
#' [life_calibrate()] use deterministic recursions.
#'
#' @seealso [LifeHist()], [life_train()], [life_calibrate()]
#' @export
#' @examples
#' example_data <- life_example_data(n = 6, seed = 1)
#' example_data
#' identical(example_data, life_example_data(n = 6, seed = 1))
life_example_data <- function(n = 200, seed = 1) {
  if (!is.numeric(n) || is.complex(n) || length(n) != 1L ||
      !is.null(dim(n)) || !is.finite(n) || n < 1 || n != floor(n) ||
      n > .Machine$integer.max) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }
  n <- as.integer(n)
  with_life_seed(seed, {
    age <- stats::runif(n, min = 40, max = 80)
    event_time <- stats::rexp(n, rate = 0.14 * exp(0.018 * (age - 60)))
    censor_time <- stats::rexp(n, rate = 0.04 * exp(0.008 * (age - 60)))
    observed <- pmin(event_time, censor_time, 8)
    event <- event_time <= censor_time & event_time <= 8
    data.frame(id = seq_len(n), time = observed,
               status = ifelse(event, 2L, 0L), age = age)
  })
}
