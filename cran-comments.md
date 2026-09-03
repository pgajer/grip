## Post-publication validation, 2026-09-03

This records maintenance validation after publication of grip 0.2.0, not a
new submission. No upload is planned and the published version tag is unchanged.

Changes remove local DESCRIPTION-file path attributes from benchmark session
metadata, clarify provenance and thread-control documentation, and consistently
exclude optional generated images from source builds. Numerical benchmark
results and the public API are unchanged. NEWS separates these development
changes from the published 0.2.0 release.

## Test environment and results

macOS Tahoe 26.6.1 (Apple Silicon), R-devel 4.7.0
(2026-06-24 r90190), Apple clang, C++17. This is a dated local R-devel snapshot,
not a claim to have tested the latest R-devel revision.

Two full checks of the same fresh tarball, including a repeat with
`R_MAKEVARS_USER=/dev/null`, each reported:

0 errors | 1 warning | 1 note

* Warning: version 0.2.0 is already published on CRAN. Expected for this
  retrospective check; this version must not be uploaded again as an update.
* Note: the local HTML Tidy is too old for HTML manual validation.

Each tarball check reported 2,255 passing test assertions, no failures or
warnings, and six optional Shiny skips. A direct development test run covered
the optional tests with no failures, warnings, or skips. Documentation and
README regeneration produced no unexplained changes.

The corrected optional grip example in the local gflow vignette was executed
with grip 0.2.0, dgraphs 0.2.0, and ivue 0.1.0 and produced finite 250-by-3
coordinates and an ivue visualization.

Tarball SHA-256:
`51eb82c7d5982858c6e39942ae9b81140b9ceb017cd0d8e59e73320c5f314c3f`.

## Before any future submission

Choose an appropriate new version, rebuild and recheck its exact tarball on
current release/development platforms, review downstream compatibility, and
replace this retrospective record with the submission-specific results.
