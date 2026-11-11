#' Train prediction scores and nuisance models for life histories
#'
#' Fit the working outcome and censoring models on training subjects. The
#' returned object freezes preprocessing, state definitions, target definitions,
#' and model fits for subsequent independent calibration with [life_calibrate()].
#'
#' @param formula A formula such as `LifeHist(time, status) ~ age + treatment`.
#'   Use `~ 1` on the right for an intercept-only model.
#' @param data A data frame of chronologically ordered records within subjects.
#' @param absorbing Explicit absorbing-state labels. Use `integer(0)` for none.
#' @param targets Named list of state sets whose occupation times are predicted,
#'   for example `list(event_free = 0, alive = c(0, 1))`.
#' @param tau Positive common prediction horizon in the units of the time data.
#' @param id Identifier column name for formula input (default `"id"`), or a
#'   vector with one identifier per record for vector input. Identifiers must
#'   identify the same subject across training and calibration samples.
#' @param initial_state Initial state, either common, a vector with one entry
#'   per record, or a one-sided column formula. NULL uses state zero for event
#'   records and the first `from` state for explicit intervals.
#' @param status_map Named vector mapping nonzero status codes to state labels.
#'   Status zero always means censoring. Use a nonzero code to record a return
#'   to state zero, for example `c("1" = 1, "2" = 2, "3" = 0)`.
#' @param censor_model Censoring model: `"cox"` or `"survival_forest"`.
#' @param outcome_model Outcome model: `"markov_cox"`, `"transition_forest"`, or
#'   `"transition_neural"`.
#' @param censor_formula Optional one-sided formula selecting censoring
#'   predictors. NULL reuses the outcome predictors and transformations.
#' @param censor_options,outcome_options Named lists of model fitting options.
#'   See [life_model_options] for supported names, defaults and backend requirements.
#' @param augmentation_model NULL reuses the fitted outcome model. Alternatively
#'   select one of the outcome model names to fit a separate augmentation model
#'   on these training subjects, leaving the score unchanged.
#' @param augmentation_options Options for a separate augmentation fit.
#' @param allowed_transitions Optional structural transition specification.
#'   Empirically observed transitions define the fitted graph; they do not by
#'   themselves rule out other transitions when certifying coverage.
#' @param score_weights Optional named positive relative widths, one per target.
#'   NULL uses one for every target. A smaller value gives that target a
#'   narrower interval at a fixed common threshold. Weights are frozen during
#'   training and are distinct from component-specific error allocations.
#' @param control Named list of numerical controls. `scale_floor` (default 0.01)
#'   floors the standard deviation at `scale_floor * tau` before multiplication
#'   by `score_weights`. `survival_floor` (1e-8) truncates fitted censoring survival coherently.
#'   `grid_size` (128), `max_nodes` (100000), and `probability_tolerance` (0.01)
#'   control deterministic occupation calculations. `calibration_accuracy`
#'   (`"auto"`) selects threshold accuracy for Markov calibration and the
#'   probability/Monte Carlo path for flexible augmentation. Explicit
#'   `"probability"` and `"threshold"` modes are described in [life_calibrate()]. Working-CDF
#'   prediction always enforces `probability_tolerance`. `draws` (1000) controls
#'   flexible-model continuation draws. `keep_contributions` (FALSE) retains
#'   subject-level calibration contributions when requested.
#' @param seed Optional integer seed for randomized models. An explicit seed
#'   preserves the caller's random-number state and initializes separate
#'   censoring, outcome and augmentation fitting streams. Changing one model
#'   role then leaves the other roles' random draws unchanged. NULL uses the
#'   caller's stream.
#' @param x Baseline covariate matrix or data frame, repeated unchanged within
#'   subject, for vector input.
#' @param response A [LifeHist()] response for vector input. Alternatively supply
#'   its individual fields below.
#' @param time,status,start,stop,from,to,event Response vectors; see [LifeHist()].
#'
#' @details
#' Every subject contributes a complete sequence of observed transitions and
#' either an absorbing entry or a final censoring record. All covariates are
#' baseline covariates and must remain unchanged within a subject. No records
#' are silently sorted or omitted. This version does not perform sample
#' splitting or cross-fitting.
#'
#' Markov-Cox target moments and conditional coverage use deterministic reward
#' recursions. The flexible outcome models use history-dependent transition
#' models and continuation draws. The score model determines candidate regions;
#' calibration estimates coverage using separately fitted nuisance functions.
#' Reusing an outcome fit for these two roles saves fitting work.
#'
#' @returns A `lifeconformal_fit` containing `problem`, `schema`, `preprocessor`,
#'   fitted `censor`, `outcome`, and `augmentation` models, numerical `control`,
#'   `functions` for evaluating fitted predictions, and training diagnostics.
#'   Objects can be saved with [saveRDS()] and restored with [readRDS()].
#' @export
#' @examples
#' histories <- life_example_data(80, seed = 1)
#' fit <- life_train(LifeHist(time, status) ~ age, data = histories,
#'                   absorbing = 2, targets = list(alive = 0), tau = 5)
#' fit
#' predict(fit, newdata = data.frame(age = c(50, 70)))$mean
life_train <- function(formula = NULL, data = NULL, absorbing, targets, tau,
                        id = NULL, initial_state = NULL, status_map = NULL,
                        censor_model = c("cox", "survival_forest"),
                        outcome_model = c("markov_cox", "transition_forest",
                                          "transition_neural"),
                        censor_formula = NULL, censor_options = list(),
                        outcome_options = list(), augmentation_model = NULL,
                        augmentation_options = list(),
                        allowed_transitions = NULL, score_weights = NULL, control = list(),
                        seed = NULL, x = NULL, response = NULL,
                        time = NULL, status = NULL, start = NULL, stop = NULL,
                        from = NULL, to = NULL, event = NULL) {
  call <- match.call()
  if (missing(absorbing) || missing(targets) || missing(tau))
    base::stop("Supply absorbing, targets, and tau explicitly.", call. = FALSE)
  censor_model <- match.arg(censor_model)
  outcome_model <- match.arg(outcome_model)
  if (!is.null(augmentation_model)) {
    augmentation_model <- match.arg(augmentation_model,
      c("markov_cox", "transition_forest", "transition_neural"))
  } else if (length(augmentation_options)) {
    base::stop("augmentation_options requires a separate augmentation_model.",
               call. = FALSE)
  }
  control_overrides <- names(control)
  control <- resolve_life_controls(control)
  if (!is.null(censor_formula) &&
      (!inherits(censor_formula, "formula") || length(censor_formula) != 2L))
    base::stop("censor_formula must be a one-sided formula.", call. = FALSE)
  response <- resolve_public_response(formula, response, time, status, start,
                                       stop, from, to, event)
  input <- resolve_input(formula, data, x, id, response,
                         extra_predictors = all.vars(censor_formula))
  history <- normalize_histories(input, initial_state, status_map, absorbing)
  if (!is.null(history$schema$formula))
    environment(history$schema$formula) <- environment(LifeHist)
  problem <- compile_problem(history, targets, tau,
                              allowed_transitions = allowed_transitions)
  target_names <- colnames(problem$rewards)
  if (is.null(score_weights)) score_weights <- stats::setNames(rep(1, length(target_names)), target_names)
  if (!is.numeric(score_weights) || is.complex(score_weights) ||
      is.null(names(score_weights)) || anyNA(names(score_weights)) ||
      anyDuplicated(names(score_weights)) ||
      !setequal(names(score_weights), target_names) ||
      any(!is.finite(score_weights)) || any(score_weights <= 0))
    base::stop("score_weights must give one named positive finite width per target.", call. = FALSE)
  problem$score_weights <- score_weights[target_names]
  prep <- fit_preprocessor(history$baseline, formula, censor_formula)
  if (!is.list(censor_options) || !is.list(outcome_options) ||
      !is.list(augmentation_options))
    base::stop("Model options must be named lists.", call. = FALSE)
  if (is.null(censor_options$survival_floor))
    censor_options$survival_floor <- control$survival_floor
  else {
    control$survival_floor <- censor_options$survival_floor
    control <- resolve_life_controls(control)
  }
  models <- with_life_seed(seed, {
    # Allocate role streams before any backend can consume random numbers.
    # With no explicit public seed, retain the caller-stream convention and
    # avoid touching the RNG for wholly deterministic fitting.
    fitting_seeds <- if (is.null(seed)) NULL else stats::setNames(
      as.list(sample.int(.Machine$integer.max, 3L)),
      c("censor", "outcome", "augmentation"))
    censor <- with_life_seed(fitting_seeds$censor,
      fit_life_censor(history, prep$censor$x, problem,
                     censor_model, censor_options))
    outcome <- with_life_seed(fitting_seeds$outcome,
      fit_life_outcome(history, prep$outcome$x, problem,
                      outcome_model, outcome_options))
    augmentation <- if (is.null(augmentation_model)) outcome else
      with_life_seed(fitting_seeds$augmentation,
        fit_life_outcome(history, prep$outcome$x, problem,
                         augmentation_model, augmentation_options))
    list(censor = censor, outcome = outcome, augmentation = augmentation,
         fitting_seeds = fitting_seeds)
  })
  fit <- structure(c(list(problem = problem, schema = history$schema,
                          preprocessor = prep$preprocessor), models,
                     list(control = control, control_overrides = control_overrides, seed = seed,
                          training_ids = history$subjects$id,
                          diagnostics = history$diagnostics,
                          augmentation_reused = is.null(augmentation_model),
                          call = call, object_version = 1L,
                          package_version = "0.1.0")),
                   class = "lifeconformal_fit")
  fit$functions <- make_fitted_functions(fit)
  fit
}

#' Dispatch model training at the adapter boundary
#' @noRd
fit_life_censor <- function(history, x, problem, type, options) {
  switch(type, cox = fit_censor_cox(history, x, problem$tau, options),
         survival_forest = fit_censor_forest(history, x, problem$tau, options))
}

#' @rdname fit_life_censor
#' @noRd
fit_life_outcome <- function(history, x, problem, type, options) {
  switch(type,
         markov_cox = fit_markov_cox(history, x, problem$tau, options,
                                     edges = problem$working_graph$edges),
         transition_forest = fit_transition_forest(history, x, problem, options),
         transition_neural = fit_transition_neural(history, x, problem, options))
}

#' Capture only the fitted prediction state in callable interfaces
#' @noRd
make_fitted_functions <- function(fit) {
  # Exclude training identifiers, the call and the function list itself from
  # closure environments. Shared R objects are not copied until modified.
  env <- new.env(parent = environment(make_fitted_functions))
  env$s <- fit[c("problem", "schema", "preprocessor", "censor", "outcome",
                 "augmentation", "control", "seed")]
  eval(quote(list(
    moments = function(newdata, initial_state = NULL)
      predict_life_moments(s, newdata, initial_state),
    score = function(newdata, values, initial_state = NULL)
      evaluate_life_score(s, newdata, values, initial_state),
    coverage = function(newdata, q, initial_state = NULL)
      predict_life_coverage(s, newdata, q, initial_state),
    censor_survival = function(newdata, time, side = c("left", "right")) {
      v <- apply_preprocessor(newdata, s$preprocessor, "censor")
      predict_life_censor(s$censor, v$x, time, match.arg(side))
    })), envir = env)
}
