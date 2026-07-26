## Resubmission

This is a resubmission of `grip`, an R package implementing the GRIP
algorithm and related weighted and geodesic methods for 2D and 3D graph
layout.

Per CRAN feedback, I replaced the unnecessary `\dontrun{}` wrappers and
rewrote the runnable examples so they execute quickly.

Since the previous submission, version 0.1.1 adds weighted, geodesic-MDS,
geodesic-KK, landmark, MISF, and interactive layout workflows. Compiled
parallel work is capped at two threads, including automatic thread selection.
Non-package utilities and machine-specific paths are not installed.

## Test environments

* macOS Tahoe 26.3.1 (Apple Silicon), R 4.6.1
* Apple clang 17.0.0, C++17, clean user compiler settings via
  `R_MAKEVARS_USER=/dev/null`
* macOS Tahoe 26.3.1 (Apple Silicon), R-devel 4.7.0 (2026-06-24 r90190)

## R CMD check results

0 errors | 0 warnings | 2 notes

The clean submission-style check was run with R-devel on the exact source
archive `grip_0.1.1.tar.gz`. All 2,194 test expectations passed with no skips.
The test suite is self-contained and does not depend on another source
checkout.

The two notes were the expected `New submission` incoming-feasibility note and
a local-tool note that recent HTML Tidy and V8 were unavailable for optional
HTML manual validation.

The package name was also checked against the current CRAN and Bioconductor
package indexes; no case-insensitive `grip` match was found.

## Notes

This is a revised new-package submission responding to prior CRAN feedback.
