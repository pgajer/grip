#!/usr/bin/env Rscript
# Run after make build. Isolate a complete dependency library: R CMD check's
# restricted-library setup assumes packages under .Library have all transitive
# dependencies there too, which is not true for some split local installations.
args <- commandArgs(TRUE)
if (length(args) != 1L || !file.exists(args[[1L]])) {
  stop("Usage: Rscript tools/pkg/check-full.R path/to/grip_version.tar.gz")
}
lib <- tempfile("grip-full-library-")
dir.create(lib)
packages <- installed.packages()
packages <- packages[!duplicated(packages[, "Package"]), , drop = FALSE]
packages <- packages[is.na(packages[, "Priority"]), , drop = FALSE]
for (i in seq_len(nrow(packages))) {
  from <- file.path(packages[i, "LibPath"], packages[i, "Package"])
  to <- file.path(lib, packages[i, "Package"])
  ok <- if (.Platform$OS.type == "windows") file.copy(from, lib, recursive = TRUE) else file.symlink(from, to)
  if (!ok) stop("Could not stage dependency: ", from)
}
Sys.setenv(R_LIBS = lib, GRIP_FULL_DEPENDENCIES = "true", RGL_USE_NULL = "true",
           R_MAKEVARS_USER = if (.Platform$OS.type == "windows") "NUL" else "/dev/null")
cat("Full-dependency library:", lib, "\n")
print(packages[, c("Package", "Version", "LibPath"), drop = FALSE])
status <- system2(file.path(R.home("bin"), "R"),
                  c("CMD", "check", "--as-cran", shQuote(args[[1L]])))
unlink(lib, recursive = TRUE)
quit(status = status)
