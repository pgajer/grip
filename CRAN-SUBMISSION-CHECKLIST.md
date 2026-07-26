# CRAN Submission Checklist for `grip`

## Ready now

- [x] Build the source tarball with `R CMD build .`
- [x] Run local `R CMD check --as-cran` on `grip_0.1.1.tar.gz`
- [x] Re-run `R CMD check --as-cran` in a clean local build environment with
  `R_MAKEVARS_USER=/dev/null`
- [x] Confirm vignette build succeeds
- [x] Confirm tests pass locally: 2,194 pass, 0 fail, 0 warn, 0 skips
- [x] Confirm URL checks pass with `urlchecker::url_check()`
- [x] Draft `cran-comments.md`
- [x] Add Linux GitHub Actions check prep
- [x] Add Win-builder submission helper script
- [x] Polish the `Description` field in `DESCRIPTION`
- [x] Add short runnable Rd examples
- [x] Confirm `pkgdown::check_pkgdown()` reports no problems
- [x] Confirm package-native parallel work is capped at two threads
- [x] Remove development-only installed files and machine-specific paths
- [x] Remove cross-repository test dependencies and retain self-contained
  parity coverage
- [x] Confirm the maintainer email is active
- [x] Confirm Peter Gajer is the sole package author
- [x] Confirm the provenance and release suitability of bundled data and assets
- [x] Confirm GPL-3 is the intended package license
- [x] Confirm no current Bioconductor 3.23 package is named `grip`
- [x] Confirm the exact tarball is 4,102,926 bytes with SHA-256
  `f0a76ddb95854a237aab2f9b59e996a9e9145a5adf0093894f552709ef5f3523`

## Recommended before upload

- [x] Submit the exact archive to Win-builder `R-release`, `R-devel`, and
  `R-oldrelease`
- [x] Confirm the Win-builder `R-release` result reports 0 errors, 0 warnings,
  and only the expected new-submission note
- [ ] Confirm the Win-builder `R-devel` and `R-oldrelease` result emails report
  no errors, warnings, or unexpected notes
- [x] Complete the local R-devel check: 0 errors, 0 warnings, and two
  environment/submission-status notes (`New submission`; local HTML
  Tidy/V8 unavailable)
- [x] Commit and push the release changes
- [x] Confirm the GitHub Actions `R-devel`, `R-release`, and `R-oldrel-1`
  matrix and pkgdown workflow pass on the release commit
- [x] Run and confirm R-hub Linux, Windows, and macOS checks
- [x] Verify the exact tarball to upload is `grip_0.1.1.tar.gz`

## Submission steps

- [ ] Submit through the CRAN web form at `https://cran.r-project.org/submit.html`
- [ ] Confirm the submission email
- [ ] Include the text from `cran-comments.md` in the submission comments field
- [ ] Do not resubmit while this version is still pending unless CRAN asks for changes

## Useful commands

- `R CMD build .`
- `env R_MAKEVARS_USER=/dev/null R CMD check --as-cran grip_0.1.1.tar.gz`
- `make check-clean`
- `make winbuilder-release`
- `make winbuilder-devel`
