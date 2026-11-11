# lifeconformal 0.1.0

* Added `life_train()`, `life_calibrate()`, and S3 prediction, print, summary,
  and plotting methods for independent training and calibration subjects.
* Added formula and vector interfaces for ordered event-status and explicit
  transition histories, with initial-state rules, absorbing states, and
  separate structural and empirical transition graphs.
* Added frozen predictor preprocessing, Cox censoring regression and
  transition-specific Markov-Cox outcome fitting. Optional forest censoring,
  transition forests, and neural transition adapters support flexible models.
* Added IPW and augmented calibration, raw first-crossing threshold inversion,
  component error allocation, and full-support fallback diagnostics.
* Added deterministic Markov reward moments, exact monotone-target coverage,
  bounded occupation-grid calculations, and optimized structural coverage
  certificates. Markov calculations do not use path simulation.
* Added `project_life()` for linear summaries of the stored joint prediction
  region, retaining overlapping-target and reachable-occupation constraints.
* Added reproducible synthetic example data, an executable workflow vignette,
  reference documentation, and regression tests.
* Explicit fitting seeds use independent streams for censoring, outcome and
  augmentation roles. Fitted summaries expose effective model controls and
  diagnostics. Numerical failures identify supported recovery controls.
* Formula responses require data columns, and matrix prediction inputs retain
  the trained predictor naming rule. Standard plot labels, limits and graphical
  parameters can be customized through the S3 methods.
* History queries no longer alter propagated Markov distributions. Accurate
  binary occupation totals preserve deterministic boundary values in grid,
  monotone, moment, flexible-prediction and structural stopping calculations.
  Joint grids retain each component's no-change proof through coalescing.
* Valid flexible-model intervals no longer create an artificial empty terminal
  cell from rounded cell counts. Native forest caches use portable public R
  interfaces compatible with the declared minimum R version.
* AIPW calibration with multiple active flexible components shares continuation
  draws in bounded query batches and avoids full contribution matrices unless
  requested. Ordered reductions preserve curves, raw-sum threshold decisions
  and diagnostics. High-draw flexible calibration remains computationally
  expensive despite the measured improvement.
* Added a package citation, a maintainer bug-report contact, and complete
  installed model-option documentation.
