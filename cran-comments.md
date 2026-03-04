# CRAN submission comments -- pecanr 0.1.2

## Test environments
* Local: macOS 26.2, R 4.5.0 (aarch64-apple-darwin20)
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
In response to CRAN reviewer feedback:
* Wrapped package/software names in single quotes in Description
  (e.g. 'lme4')
* Added references in Description field in the required
  authors (year) <doi:...> format
* Replaced \dontrun{} with \donttest{} throughout examples
* Added Correll, Mellinger, and Pedersen (2022) <doi:10.3758/s13428-021-01687-2>
  to references and DESCRIPTION

## Downstream dependencies
This is a new package. There are no existing downstream dependencies.
