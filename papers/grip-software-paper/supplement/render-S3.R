#!/usr/bin/env Rscript
lib <- Sys.getenv("GRIP_RJOURNAL_PACKAGE_LIBRARY")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
file_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)[1]
source_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
paper_dir <- dirname(source_dir)
args <- commandArgs(TRUE)
output <- if (length(args)) args[1] else file.path(paper_dir, "build", "supplement")
dir.create(output, recursive=TRUE, showWarnings=FALSE)
output <- normalizePath(output)
rmarkdown::render(file.path(source_dir, "S3-controlled-examples.Rmd"),
  output_dir=output, intermediates_dir=output, knit_root_dir=paper_dir,
  clean=FALSE, quiet=FALSE, envir=new.env(parent=globalenv()))
