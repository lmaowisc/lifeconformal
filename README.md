# lifeconformal

Predict occupation times from censored multistate life histories. Train a fixed
score and nuisance models, calibrate on independent subjects, and obtain scalar
intervals or joint prediction regions for new patients.

**Preprint.** Lu Mao (2026). Conformal prediction for censored life history data.
[arXiv:2611.15724](https://arxiv.org/abs/2611.15724).

## Installation

Install the package from GitHub with R 4.1.0 or later. Installation from source
requires a C/C++ toolchain. To build the tutorial, also install Pandoc or use
the copy bundled with RStudio.

```r
install.packages("remotes")
remotes::install_github("lmaowisc/lifeconformal", dependencies = TRUE,
                        build_vignettes = TRUE)
```

## Example

```r
library(lifeconformal)                         # Load the package
histories <- life_example_data(120, seed = 1)  # Simulate 120 patient histories
training <- histories[histories$id <= 80, ]    # Use 80 patients for training
calibration <- histories[histories$id > 80, ]  # Reserve 40 patients for calibration

fit <- life_train(               # Fit outcome and censoring models
  LifeHist(time, status) ~ age,  # Model the life history by baseline age
  data = training,               # Use the training histories
  absorbing = 2,                 # State 2 (death) is absorbing
  targets = list(alive = 0),     # Predict time alive in state 0
  tau = 5)                       # Predict through five time units

cal <- life_calibrate(fit, data = calibration)       # Augmented calibration (default)
predict(cal, newdata = data.frame(age = c(50, 70)))  # 90% prediction intervals at ages 50 and 70

cal_ipw <- life_calibrate(fit, data = calibration,  # Reuse the same fit and patients
                          method = "ipw")           # Inverse-probability weighting alone
```

Each row records a transition or the end of observation. In the example, state
0 denotes the initial state, state 2 is absorbing, and the target is time
spent in state 0 through time 5. Training and calibration use separate patients.

The [package tutorial](vignettes/lifeconformal.Rmd) covers joint prediction for
multiple state occupation times, alternative working models, and calibration
options. After installation with vignettes, open it with
`vignette("lifeconformal")`. Use `citation("lifeconformal")` for the package
citation and [GitHub Issues](https://github.com/lmaowisc/lifeconformal/issues)
to report a problem.
