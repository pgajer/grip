## Submission

This is an update to `grip` that unifies topology-first and edge-length-metric
layouts under `grip()` and gives `trace.grip()` the same metric selection. The
short-lived `weighted.grip()` and `trace.weighted.grip()` compatibility entry
points have been removed. Documentation and tests now use the unified API.

## Test environments

* macOS Tahoe 26.6.1 (Apple Silicon), R-devel 4.7.0
  (2026-06-24 r90190), Apple clang, C++17
* GitHub Actions: R-devel, R-release, and R-oldrel-1

## R CMD check results

0 errors | 0 warnings | 1 note

The local `R CMD check --as-cran` completed successfully under R-devel. The
single note is environment-only: the locally installed HTML Tidy is too old
for optional HTML manual validation.

The local test suite reported 2,205 passes and three skips because the
suggested package `DT` was not installed. The GitHub Actions
R-devel/R-release/R-oldrel-1 matrix, where suggested dependencies are
installed, completed successfully for the version 0.1.3 package code.

## Notes

The words flagged by the incoming spell check are intentional: Gajer and
Kobourov are author surnames, KK is the standard abbreviation for
Kamada-Kawai, and multiscale is a technical term.
