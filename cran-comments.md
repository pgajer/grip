## Submission

This is a new submission of `grip`, an R package implementing the GRIP
algorithm for 2D and 3D graph layout.

## Test environments

* macOS 26.3.1 (Apple Silicon), R 4.5.2
* macOS 26.3.1 (Apple Silicon), R 4.5.2, clean user build flags via
  `R_MAKEVARS_USER=/dev/null`

## R CMD check results

0 errors | 0 warnings | 2 notes

## Notes

1. `This is a new submission.`
2. `HTML Tidy` is not recent enough on the local machine, so HTML validation
   was skipped during `R CMD check --as-cran`:

   `Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.`

   This is an environment-specific note rather than a package issue.

## Additional checks

* `urlchecker::url_check()` reports that all URLs are correct.
* Package tests pass locally.
