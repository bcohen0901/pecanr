# CRAN submission comments -- pecanr 0.3.1

This is a patch release fixing two bugs in the error-variance computation
that silently inflated effect sizes for common model specifications (see
Changes below). We are submitting it promptly after discovering the bugs.

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

* On win-builder R-devel, the 0.3.0 CRAN incoming feasibility check flagged:
  "Possibly misspelled words in DESCRIPTION: df, semipartial". These are not
  misspellings: "df" abbreviates degrees of freedom, and "semipartial" is the
  standard statistical term for the variance-explained measure the package
  computes. The DESCRIPTION text is unchanged in this release, so the same
  note may appear again.

## Changes since last submission

This release contains only bug fixes and their regression tests; there are no
API changes.

* Fixed a bug in which random-effect terms stored under numeric-suffixed
  group names were silently dropped from the error variance. lme4 stores
  multiple random-effect terms for the same grouping factor -- as produced by
  the double-bar syntax `(x || g)` or by explicit `(1 | g) + (0 + x | g)`
  terms -- under group names `"g"`, `"g.1"`, `"g.2"`, etc. The error-variance
  calculators matched only the base name, so components under suffixed names
  (typically random slopes) were omitted from the denominator, inflating
  eta2p. Affected all three design calculators (crossed, nested, mixed) and
  both general and operative effect sizes.
* Fixed a related bug for nested-syntax grouping terms such as
  `(1 | photo:model)`: `as.data.frame(VarCorr())` reports the combined group
  name dot-separated on current lme4 (e.g. `"photo.model"`), so the previous
  colon-based fallback never matched and the component was silently dropped
  in nested and mixed designs.
* Added a regression test suite (`tests/testthat/test-double-bar.R`) that
  verifies the error variance against a direct hand computation from
  `VarCorr()` across crossed, nested, and mixed designs for the affected
  syntaxes, and guards against false group-name matches.

## Downstream dependencies

There are no downstream dependencies.
