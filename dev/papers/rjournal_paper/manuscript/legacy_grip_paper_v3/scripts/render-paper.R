#!/usr/bin/env Rscript

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to render grip-paper-v3.Rmd.")
}

if (!requireNamespace("rjtools", quietly = TRUE)) {
  stop("Package 'rjtools' is required to render the R Journal paper.")
}

input_file <- "grip-paper-v3.Rmd"
if (!file.exists(input_file)) {
  stop(input_file, " was not found in the current working directory.")
}

args <- commandArgs(trailingOnly = TRUE)
render_pdf <- !("--html-only" %in% args)
render_html <- "--html" %in% args || "--all" %in% args
timestamped <- !("--no-timestamp" %in% args)

timestamp_suffix <- format(Sys.time(), "%Y%m%d_%H%M%S")
base_name <- tools::file_path_sans_ext(input_file)

output_name <- function(ext) {
  if (timestamped) {
    sprintf("%s_%s.%s", base_name, timestamp_suffix, ext)
  } else {
    sprintf("%s.%s", base_name, ext)
  }
}

stable_name <- function(ext) {
  sprintf("%s.%s", base_name, ext)
}

render_one <- function(output_format, ext) {
  stable_out <- stable_name(ext)
  rmarkdown::render(
    input = input_file,
    output_format = output_format,
    quiet = FALSE,
    envir = new.env(parent = globalenv())
  )
  if (!timestamped) {
    return(stable_out)
  }

  stamped_out <- output_name(ext)
  ok <- file.copy(stable_out, stamped_out, overwrite = TRUE)
  if (!isTRUE(ok)) {
    stop("Failed to copy ", stable_out, " to ", stamped_out)
  }
  stamped_out
}

outputs <- character()

if (isTRUE(render_pdf)) {
  outputs <- c(outputs, render_one("rjtools::rjournal_pdf_article", "pdf"))
}

if (isTRUE(render_html)) {
  outputs <- c(outputs, render_one("rjtools::rjournal_article", "html"))
}

message("Rendered: ", paste(outputs, collapse = ", "))
