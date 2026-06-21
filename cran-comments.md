# CRAN submission comments -- pecanr 0.3.0

## Test environments
* Local: macOS Tahoe 26.2, R 4.5.3 (aarch64-apple-darwin20) -- 0 errors, 0 warnings, 0 notes
* devtools::check_win_devel(): Windows Server, R-devel -- 0 errors, 0 warnings, 1 note
* rhub::rhub_check(): linux (ubuntu-latest), R-devel -- 0 errors, 0 warnings, 0 notes
* rhub::rhub_check(): macos-arm64 (macos-latest), R-* -- 0 errors, 0 warnings, 0 notes
* rhub::rhub_check(): windows (windows-latest), R-* -- 0 errors, 0 warnings, 0 notes

## R CMD CHECK results
Local and R-hub (linux, macos-arm64, windows): 0 errors | 0 warnings | 0 notes.
win-builder R-devel: 0 errors | 0 warnings | 1 note (see Notes).

## Notes
* On win-builder R-devel, the CRAN incoming feasibility check flagged:
  "Possibly misspelled words in DESCRIPTION: df, semipartial". These are not
  misspellings: "df" abbreviates degrees of freedom, and "semipartial" is the
  standard statistical term for the variance-explained measure the package
  computes. This note did not appear on the local or R-hub checks.

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
* `eta2p()` and `batch_eta2p()` gain a `partial_predictors` argument
  controlling the numerator for single (non-interaction) predictors. By default
  (`TRUE`) the variance attributed to a predictor is its unique (semipartial)
  variance -- equivalently its total variance times its tolerance -- so it
  reflects only the predictor's unique contribution under correlation. Set
  `partial_predictors = FALSE` for the previous total-variance numerator. For
  centered, orthogonal designs the two are identical, so correlated-predictor
  models are the only ones affected.
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
