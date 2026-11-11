#' @details
#' The workflow separates score training and nuisance-model fitting from
#' calibration on an independent sample. The first version uses sample splitting;
#' cross-fitting is deferred. Use [LifeHist()] to describe observed histories,
#' [life_train()] to fit working models, [life_calibrate()] on independent
#' subjects, and [predict()] to obtain new baseline prediction regions.
#' Markov calculations use deterministic reward recursions. Flexible transition
#' models use frozen continuation streams and report their numerical settings.
#' @importFrom stats setNames predict
#'
#' @examples
#' getNamespaceVersion("lifeconformal")
#'
#' @keywords internal
"_PACKAGE"
