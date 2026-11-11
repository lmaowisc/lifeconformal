#' Working models and fitting options
#'
#' Supported model choices and named fitting options for [life_train()].
#' Model fitting and tuning belong to the training sample. Calibration reuses
#' the fitted models and the frozen score on independent subjects.
#'
#' @section Model roles:
#' `censor_model` selects `"cox"` (default) or `"survival_forest"`.
#' `outcome_model` selects `"markov_cox"` (default), `"transition_forest"`, or
#' `"transition_neural"`. These choices are independent. Supply fitting options
#' in `censor_options` and `outcome_options`. `augmentation_model = NULL`
#' reuses the outcome fit; a separate outcome family uses `augmentation_options`
#' with the same option names as that family. Separate augmentation uses the
#' outcome predictor mapping. `censor_formula` can select different baseline
#' predictors for censoring. All model options require unique, nonempty names;
#' unknown options are rejected. Use the public `seed` argument for seeds.
#'
#' @section Cox fitting:
#' Censoring and Markov-Cox outcome fitting support `ties = "breslow"`
#' (`"efron"` is also available), `singular = "drop"` (warn and set aliased
#' coefficients to zero; `"error"` rejects them), `iter.max = 50L`, and
#' `eps = 1e-9`. Nonconvergence is an error. Distinct observed times are
#' preserved with `timefix = FALSE`. Outcome fits use common-origin start-stop
#' risk intervals through actual observed follow-up, including after `tau`.
#' The prediction horizon restricts the target, not outcome fitting exposure.
#'
#' Censoring fits use one derived record per subject, ending at observation end,
#' absorption or `tau`, whichever is first. Only censoring strictly before `tau`
#' is an event. Both censoring families accept `survival_floor`, inherited from
#' `control$survival_floor` (default `1e-8`) unless explicitly supplied in
#' `censor_options`. The floor is applied to one coherent survival distribution
#' before left/right limits and increments are calculated. Increasing it changes
#' the estimator; it does not establish censoring positivity.
#'
#' For Markov transitions, `min_events_per_transition` (default zero) fits
#' a baseline-only Cox model to each edge below the specified event count.
#' By default, a Cox separation warning is an error. With
#' `separation_action = "baseline"`, a transition is refitted without
#' covariates when its partial likelihood converges before a coefficient
#' and the backend warns that this coefficient may be infinite. The
#' transition's diagnostics record the refit and the original warning.
#' Other failures retain their usual behavior.
#' @section Censoring forests:
#' `"survival_forest"` uses optional package \pkg{ranger}. Its options are
#' `num_trees = 200L`, `mtry = NULL`, `min_node_size = 10L`, `max_depth = 0L`,
#' `num_threads = 1L`, and the censoring survival floor described above.
#' Zero maximum depth leaves depth unrestricted. NULL `mtry` resolves to the
#' rounded-down square root of the design-column count, at least one;
#' intercept-only fits use a constant feature. All censoring-event knots are
#' retained. No-event data return survival one without fitting a forest.
#'
#' @section Flexible outcome fitting:
#' Both flexible families use baseline predictors plus predictable time,
#' duration in state, previous transition count and accrued target rewards.
#' Common options are `step_size = tau / 50`, `max_training_rows = 1000000L`
#' and `max_events = 1000L`. Intensities are constant within cells anchored at
#' state entry, using midpoint features under continued residence in that state.
#' Observed endpoints determine exposure, not the cell midpoint. Conditional
#' continuation inside an existing cell retains that cell's rate until its
#' boundary. Resource limits raise errors instead of silently truncating paths.
#'
#' `"transition_forest"` uses optional package \pkg{rpart} for Poisson
#' exposure/count trees. Additional defaults are `num_trees = 100L`,
#' `max_depth = 8L` (at most 30), `minsplit = 20L`, `minbucket = 7L`,
#' `mtry = NULL`, and `cp = 0.001`. NULL `mtry` resolves to the rounded-down
#' square root of the full feature count, at least one. Subjects are bootstrapped
#' as whole units; feature subsets are sampled once per tree. These intensity
#' forests differ from the right-censored survival forest used for censoring.
#'
#' `"transition_neural"` uses tanh hidden units, a softplus intensity and
#' [stats::optim()] with an analytic gradient. Additional defaults are
#' `hidden = 4L`, `decay = 0.01`, `maxit = 1000L`, and `reltol = 1e-7`.
#' The quadratic penalty excludes the output intercept; loss is normalized per
#' training subject. Nonfinite fits and nonconvergence are errors. No separate
#' neural-network package is needed.
#'
#' @section Numerical accuracy and reproducibility:
#' Markov predictions use deterministic recursions. Flexible outcomes use stored
#' profile-specific continuation streams. `control$draws` defaults to 1000;
#' increasing it examines Monte Carlo precision, whereas reducing `step_size`
#' examines time approximation. Refit and recalibrate when changing the score's
#' numerical specification. A calibration-only `draws` override changes
#' augmentation precision while retaining the trained score's draw budget.
#'
#' The score standard deviation is floored at `control$scale_floor * tau`
#' before multiplication by the target's `score_weights` entry. See
#' [life_calibrate()] for threshold versus probability accuracy and numerical
#' ambiguity. Deterministic probability brackets cover discretization error,
#' excluding floating-point rounding; flexible Monte Carlo calculations do not
#' supply those bounds. Numerical accuracy is separate from statistical coverage.
#'
#' An explicit training seed preserves the caller's random-number kind/state
#' and initializes separate fitting streams for the three model roles.
#' `seed = NULL` uses the caller's stream for randomized fitting. Save fitted
#' objects with [saveRDS()]; required optional backend packages must also be
#' installed in a fresh session. Inspect resolved options with
#' `summary(fit)$model_options` and diagnostics with
#' `summary(fit)$model_diagnostics`. Retain the analysis script and [sessionInfo()].
#'
#' @name life_model_options
#' @seealso [life_train()], [life_calibrate()], [predict.lifeconformal_fit()]
NULL
