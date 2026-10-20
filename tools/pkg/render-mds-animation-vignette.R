#!/usr/bin/env Rscript
# Run after an optimized native build; never trigger a debug compilation.
options(rgl.useNULL = TRUE)
pkgload::load_all("../ivue", recompile = FALSE, quiet = TRUE)
pkgload::load_all(".", recompile = FALSE, quiet = TRUE)
output <- file.path(getwd(), "output", "vignette-previews", "doc")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
rmarkdown::render("vignettes/metric-mds-animation.Rmd", output_dir = output,
                  envir = new.env(parent = globalenv()), quiet = FALSE)
