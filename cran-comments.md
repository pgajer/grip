## Submission

This is a deliberately breaking 0.2.0 update. It removes 40 long-form
compatibility aliases deprecated in the 0.1.x series and makes 53 low-level
graph-family construction helpers internal, reducing the export count from 194
to 101. The retained graph-family functions continue to return the same
documented graph bundles. A complete replacement map is available in
`help("grip-0.2-migration")`.

CRAN reported no reverse dependencies for `grip` when checked on 2026-08-21.

## Test environments

* macOS Tahoe 26.6.1 (Apple Silicon), R-devel 4.7.0
  (2026-06-24 r90190), Apple clang, C++17
* GitHub Actions: R-devel, R-release, and R-oldrel-1

## R CMD check results

0 errors | 0 warnings | 1 note

The local `R CMD check --as-cran` completed successfully under R-devel on
2026-08-29. The sole note is environment-only: the locally installed HTML Tidy
is too old for optional HTML manual validation.

The package-check test suite reported 2,193 passes, no failures or warnings,
and six skips for optional Shiny tests. A direct development test run reported
2,203 passes with no failures, warnings, or skips.

## Notes

The words flagged by the incoming spell check are intentional: Gajer and
Kobourov are author surnames, KK is the standard abbreviation for
Kamada-Kawai, and multiscale is a technical term.
