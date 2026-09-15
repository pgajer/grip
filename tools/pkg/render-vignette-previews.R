#!/usr/bin/env Rscript
# Standalone installed-vignette previews, with local help links. No site publish.
root <- normalizePath(".")
out <- file.path(root, "output", "vignette-previews")
dir.create(file.path(out, "doc"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "html"), recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(root, quiet = TRUE, export_all = FALSE, helpers = FALSE)
for (input in list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE)) {
  rmarkdown::render(input, output_dir = file.path(out, "doc"),
                    envir = new.env(parent = globalenv()), quiet = TRUE)
}
# Mirror installed help, including method pages, next to the previews.
for (input in list.files("man", pattern = "[.]Rd$", full.names = TRUE)) {
  tools::Rd2HTML(input, out = file.path(out, "html", sub("[.]Rd$", ".html", basename(input))),
                 package = "grip")
}
file.copy(file.path(R.home("doc"), "html", "R.css"), file.path(out, "html", "R.css"), overwrite = TRUE)
message("Rendered package vignette previews in ", file.path(out, "doc"))
