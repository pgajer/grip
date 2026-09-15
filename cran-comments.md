## Development validation, 2026-09-14

This is a development validation record for grip 0.2.0.9001, not submission
comments for an upload. CRAN currently publishes grip 0.2.0 (2026-08-31).
The recommended next release is 0.3.0; the development version is unchanged.

The package now includes two additional installed HTML vignettes: a
function guide and synthetic graph-family layout examples. The API catalog
covers all 104 explicit exports and three S3 registrations, including the
one overlapping name. This documentation work changes no solver behavior.
Previously pending changes include stress-minimizing metric MDS, a separate
classical-scaling entry point, coordinate/surface scores, and corrected
shortest-path symmetry. The MDS migration help documents the behavior change.

## Current test environment and results

* macOS Tahoe 26.6.1, Apple Silicon (aarch64-apple-darwin23).
* R-devel 4.7.0, 2026-06-24 r90190: a dated development snapshot.
* Package compilation: Apple clang 21.0.0, C++17, MacOSX26.5 SDK.
* `make check-clean`: fresh source build followed by
  `R_MAKEVARS_USER=/dev/null R CMD check --as-cran` on its tarball.

Final result: 0 errors, 0 warnings, 2 notes.

* Expected development-version note: version contains large components
  (0.2.0.9001). Select a release version before submission.
* Environment note: local HTML Tidy is too old for the optional HTML
  manual validation.

The tarball check recorded 2,308 passing test assertions, no failures or
warnings, and 15 skips: six because Shiny could not load in the test
subprocess and nine because smacof could not load there. All declared
Suggests are installed locally. A separate direct run of the metric-MDS
test file against this installed candidate passed 58 assertions without
failures, warnings, or skips, using smacof 2.1-7. This does not remove the
skips from the CRAN-style result; the optional Shiny tests were not rerun.

All six vignettes built and rebuilt, examples passed, PDF manual generation
passed, and both new vignettes' HTML, Rmd sources, and extracted R code were
verified in the tarball and installed vignette index. The installed new
vignette sources exactly match the maintained files. Documentation and
README regeneration produced no additional changes. Local pkgdown previews
and reference-link checks passed; no website was published.

The published reverse-suggesting package geosmooth 0.1.0 was source-reviewed.
Its weighted-GRIP-plus-edge-KK and classical-MDS-plus-edge-KK integration
helpers passed focused checks using the installed grip candidate. This is
not a full reverse-dependency check.

Final tarball SHA-256:
`d37cc31e4609d2f295c474acac5ab32fa7cfb20361be8de99b10c46c32c6402a`.

## Before submission

Choose the release version, build and validate that exact archive on current
R release and R-devel platforms, including Windows and Linux, and complete
reverse-dependency checks with optional dependencies available. Review the
optional-package test skips. No current-candidate Windows, Linux, Win-builder,
or R-hub result is claimed here. Confirm submission status with the maintainer
before uploading: local records contain no active submission marker, but do
not establish the state of CRAN's private queue. Replace this development
record with release-specific comments. Prior example, thread-limit, portable
filename, deterministic ordering, and native allocation fixes are preserved.
