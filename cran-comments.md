# CRAN submission comments -- pecanr 0.3.0

## Test environments
* Local: macOS [VERSION], R 4.5.3 (aarch64-apple-darwin20)
* devtools::check_win_devel(): Windows Server, R-devel -- [FILL IN RESULTS]
* rhub::rhub_check(): linux (ubuntu-latest) -- [FILL IN RESULTS]
* rhub::rhub_check(): windows (windows-latest) -- [FILL IN RESULTS]
* rhub::rhub_check(): atlas (Fedora Linux) -- [FILL IN RESULTS]

## R CMD CHECK results
0 errors | 0 warnings | 0 notes

## Changes since last submission
* New function `eta2p_omnibus()` computes a single factor-level partial
  eta-squared for a multi-level factor or multi-df interaction, corresponding
  to the omnibus (multi-df) test of that effect rather than to individual
  contrasts.
* `eta2p()` and `batch_eta2p()` now accept factor predictors directly (by
  variable name or coefficient name), reading predictor variances and factor
  random slopes from the model design matrix and model frame; manual recoding
  of factors to numeric is no longer required.
* `eta2p()` gains an optional parametric-bootstrap confidence interval
  (`ci = TRUE`, with `ci_level`, `n_boot`, and `seed` arguments) via
  `lme4::bootMer()`.
* `eta2p()` and `batch_eta2p()` gain a `partial_predictors` argument for
  unique (semipartial) variance under correlated predictors (residualizing the
  focal predictor on the others); defaults to `FALSE`.
* Interaction effect sizes are now computed from the mean-centered constituent
  predictors, making the interaction effect size invariant to the location of
  its constituent predictors. For models with uncentered continuous
  interaction components, interaction effect sizes may differ from those
  produced by earlier versions.
* The default for `verbose` in `eta2p()` is now `FALSE`, matching
  `batch_eta2p()`.
* Fixed a bug in which operative effect sizes (`operative = TRUE`) could
  misclassify factor grouping variables as within/between, because the focal
  effect's coefficient name (e.g. `"Age75"`) was matched directly against data
  columns; coefficient names are now resolved to their underlying variable
  before classification.
* Fixed a bug in which the random slope of an interaction term (e.g. a random
  slope on `x1:x2`) was silently omitted from the error denominator when its
  random-effect name did not match the fixed-effect product column.

## Downstream dependencies
There are no downstream dependencies.
