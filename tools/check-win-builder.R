#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

target <- if (length(args) >= 1L) args[[1L]] else "release"
manual <- !("--no-manual" %in% args)
email <- Sys.getenv("WIN_BUILDER_EMAIL", unset = "")

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("The 'devtools' package is required to submit to Win-builder.")
}

submit <- switch(
  target,
  release = devtools::check_win_release,
  devel = devtools::check_win_devel,
  oldrelease = devtools::check_win_oldrelease,
  stop("Target must be one of: release, devel, oldrelease")
)

email_label <- if (nzchar(email)) email else "the Maintainer address in DESCRIPTION"

message("Submitting package to Win-builder for target: ", target)
message("Results will be emailed to ", email_label, ".")
message("Pass WIN_BUILDER_EMAIL to override the default recipient.")

submit(
  pkg = ".",
  manual = manual,
  email = if (nzchar(email)) email else NULL,
  quiet = FALSE
)
