# CRAN submission comments -- pecanr 0.2.0

## Test environments
* Local: macOS [VERSION], R 4.5.0 (aarch64-apple-darwin20)
* devtools::check_win_devel(): Windows Server, R-devel -- 0 errors, 0 warnings, 1 note
* rhub::rhub_check(): linux (ubuntu-latest) -- 0 errors, 0 warnings, 0 notes
* rhub::rhub_check(): windows (windows-latest) -- 0 errors, 0 warnings, 0 notes
* rhub::rhub_check(): atlas (Fedora Linux) -- 0 errors, 0 warnings, 0 notes

## R CMD CHECK results
0 errors | 0 warnings | 1 note

## Notes
* "unable to verify current time" -- this is a system-level network issue on the
  test machine and is unrelated to the package. It appears consistently across
  all platforms and cannot be resolved by the package author.

## Changes since last submission
* Added `design = "mixed"` to `eta2p()` and `batch_eta2p()` for models with
  both crossed and nested random effects simultaneously (e.g., photos nested
  within models, crossed with participants). Includes a new internal helper
  `calc_error_mixed()` and `detect_within_between_mixed()`.
* Fixed bug in operative effect size calculation for crossed designs: 
  `detect_within_between()` previously used hardcoded `$subj`/`$item` keys
  which caused intercept variances to be silently omitted from the operative
  denominator. Keys are now indexed by actual variable name.
* Changed behavior for operative effect sizes with 3+ crossed factors: 
  third and higher factors are now correctly gated on within/between status
  rather than always being included. This is a minor breaking change for
  users with 3+ crossed factors using `operative = TRUE`.
* `batch_eta2p()` output columns for within/between status are now named
  `within_<varname>` (e.g., `within_participant`, `within_item`) rather than
  the hardcoded `within_subj`/`within_item`. This is a minor breaking change
  for code that references those columns by name.

## Downstream dependencies
This is a new package. There are no existing downstream dependencies.
