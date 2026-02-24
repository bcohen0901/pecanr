# CRAN submission comments — pecan 0.1.0

## Test environments

* Local: macOS 26.2, R 4.5.0 (aarch64-apple-darwin20)
* GitHub Actions: ubuntu-latest (release, devel), windows-latest (release), macOS-latest (release)
* devtools::check_win_devel(): Windows Server, R-devel
* rhub::check_for_cran(): Ubuntu Linux, Fedora Linux

## R CMD CHECK results

0 errors | 0 warnings | 0 notes

## Notes

* "unable to verify current time" — this is a system-level network issue on the
  test machine and is unrelated to the package. It appears consistently across
  all platforms and cannot be resolved by the package author.

## Downstream dependencies

This is a new package. There are no existing downstream dependencies.
