# lifeconformal

Predict occupation times from censored multistate life histories. Train a fixed
score and nuisance models, calibrate on independent subjects, and obtain scalar
intervals or joint prediction regions for new patients.

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
library(lifeconformal)
histories <- life_example_data(120, seed = 1)
training <- histories[histories$id <= 80, ]
calibration <- histories[histories$id > 80, ]
fit <- life_train(
  LifeHist(time, status) ~ age,
  data = training, absorbing = 2,
  targets = list(alive = 0), tau = 5
)
cal <- life_calibrate(fit, data = calibration)
predict(cal, newdata = data.frame(age = c(50, 70)))

# Reuse the trained score for inverse weighting alone.
cal_ipw <- life_calibrate(fit, data = calibration, method = "ipw")
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
