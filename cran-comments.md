# CRAN submission comments -- pecanr 0.4.0

This is a minor release that adds three features to `eta2p()` and removes one
function (see Changes below). Existing `eta2p()` estimates are unchanged.

## Test environments
* Local: macOS, R 4.5.3 (aarch64-apple-darwin20) -- 0 errors, 0 warnings, 0 notes
  (an "unable to verify current time" note appears only when the check machine
  is offline and is environmental).
* win-builder (R-devel) and R-hub (linux, macos-arm64, windows): to be run
  before submission.

## R CMD check results
Local: 0 errors | 0 warnings | 0 notes, excluding the environmental
"unable to verify current time" note on an offline machine.

## Notes

* As in prior releases, CRAN's incoming feasibility check may flag "Possibly
  misspelled words in DESCRIPTION: df, semipartial". These are not
  misspellings: "df" abbreviates degrees of freedom, and "semipartial" is the
  standard statistical term for the variance-explained measure the package
  computes.

## Changes since last release

This release adds features and removes one function. Existing `eta2p()`
estimates (general, operative, and interaction effects) are unchanged, verified
against 0.3.1.

* `eta2p()` gains a `slope` argument: supply a hypothetical fixed-effect
  coefficient (e.g. `slope = 0.03`) and the effect size is computed for that
  slope, holding the predictor variance and the error structure fixed from the
  model. Useful for sensitivity and power analyses (e.g. alongside simr).
* The returned `eta2p_lmm` object now also reports the numerator ingredients:
  `coefficient`, `var_x` (predictor variance), `var_x_unique`, and `tolerance`.
* Added concise `print()` and comprehensive `summary()` methods for
  `eta2p_lmm` objects.
* Removed `eta2p_omnibus()` pending a revised numerator definition (the current
  variance-of-the-summed-contribution numerator overstates the effect when a
  covariate is correlated with the factor).

## Downstream dependencies

There are no downstream dependencies.
