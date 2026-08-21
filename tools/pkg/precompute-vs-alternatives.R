#!/usr/bin/env Rscript

## Compatibility entry point for repository workflows. The distributed and
## submission-facing generator lives under inst/scripts so that it is included
## in the package source and installed package.

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (!length(file_arg)) {
  stop("Could not determine the path to this script.")
}

script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
source(file.path(repo_root, "inst", "scripts", "precompute-vs-alternatives.R"))
