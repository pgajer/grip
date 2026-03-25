#!/usr/bin/env Rscript

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to render README.Rmd.")
}

if (!file.exists("README.Rmd")) {
  stop("README.Rmd was not found in the current working directory.")
}

args <- commandArgs(trailingOnly = TRUE)
render_html <- "--html" %in% args || "--all" %in% args

rmarkdown::render(
  input = "README.Rmd",
  output_format = "github_document",
  output_file = "README.md",
  quiet = TRUE,
  envir = new.env(parent = globalenv())
)

if (isTRUE(render_html)) {
  rmarkdown::render(
    input = "README.Rmd",
    output_format = "html_document",
    output_file = "README.html",
    quiet = TRUE,
    envir = new.env(parent = globalenv())
  )
} else if (file.exists("README.html")) {
  unlink("README.html")
}

message("Rendered README.md from README.Rmd")
