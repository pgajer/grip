# CRAN Submission Checklist for `grip`

## Ready now

- [x] Build the source tarball with `R CMD build .`
- [x] Run local `R CMD check --as-cran` on `grip_0.1.0.tar.gz`
- [x] Re-run `R CMD check --as-cran` in a clean local build environment with
  `R_MAKEVARS_USER=/dev/null`
- [x] Confirm vignette build succeeds
- [x] Confirm tests pass locally
- [x] Confirm URL checks pass with `urlchecker::url_check()`
- [x] Draft `cran-comments.md`
- [x] Add Linux GitHub Actions check prep
- [x] Add Win-builder submission helper script
- [x] Polish the `Description` field in `DESCRIPTION`
- [x] Add short runnable Rd examples

## Recommended before upload

- [ ] Run Win-builder checks for `R-release` and `R-devel`
- [ ] Trigger the Linux GitHub Actions workflow and confirm it passes
- [ ] Optionally run one extra non-macOS check, such as R-hub
- [ ] Verify the exact tarball to upload is `grip_0.1.0.tar.gz`

## Submission steps

- [ ] Submit through the CRAN web form at `https://cran.r-project.org/submit.html`
- [ ] Confirm the submission email
- [ ] Include the text from `cran-comments.md` in the submission comments field
- [ ] Do not resubmit while this version is still pending unless CRAN asks for changes

## Useful commands

- `R CMD build .`
- `env R_MAKEVARS_USER=/dev/null R CMD check --as-cran grip_0.1.0.tar.gz`
- `make check-clean`
- `make winbuilder-release`
- `make winbuilder-devel`
