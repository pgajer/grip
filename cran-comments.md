## Resubmission

This is a resubmission of `grip`, an R package implementing the GRIP
algorithm for 2D and 3D graph layout.

Per CRAN feedback, I replaced the unnecessary `\dontrun{}` wrappers and
rewrote the runnable examples so they execute quickly.

## Test environments

* macOS Tahoe 26.3.1 (Apple Silicon), R 4.5.2
* macOS Tahoe 26.3.1 (Apple Silicon), R 4.5.2, clean user compiler settings
  via `R_MAKEVARS_USER=/dev/null`

## R CMD check results

0 errors | 0 warnings | 2 notes

The clean submission-style check was run with:

* `env R_MAKEVARS_USER=/dev/null R CMD build .`
* `env R_MAKEVARS_USER=/dev/null R CMD check --as-cran grip_0.1.0.tar.gz`

## Notes

1. `This is a new submission.`

   This is expected because `grip` has not yet been accepted by CRAN; this
   upload is a revised resubmission responding to CRAN feedback.
2. `HTML Tidy` is not recent enough on the local machine, so HTML validation
   was skipped during `R CMD check --as-cran`:

   `Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.`

   This is an environment-specific note rather than a package issue.
