#!/usr/bin/env Rscript

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to render grip-software-paper.Rmd.")
}

if (!requireNamespace("rjtools", quietly = TRUE)) {
  stop("Package 'rjtools' is required to render the R Journal paper.")
}

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (!length(file_arg)) {
  stop("Could not determine the path to this script.")
}

script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1L])))
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)
input_file <- "grip-software-paper.Rmd"
if (!file.exists(input_file)) {
  stop(input_file, " was not found in the current working directory.")
}

args <- commandArgs(trailingOnly = TRUE)
render_pdf <- !("--html-only" %in% args)
render_html <- "--html" %in% args || "--all" %in% args
timestamped <- !("--no-timestamp" %in% args)

timestamp_suffix <- format(Sys.time(), "%Y%m%d_%H%M%S", tz = "America/New_York")
base_name <- tools::file_path_sans_ext(input_file)
build_dir <- file.path(getwd(), "build")
build_info_tex <- file.path(build_dir, "manuscript_build_info.tex")

escape_tex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x, perl = TRUE)
  x
}

write_build_info <- function() {
  dir.create(build_dir, recursive = TRUE, showWarnings = FALSE)

  git_version <- tryCatch(
    system2("git", c("-C", repo_root, "describe", "--tags", "--always", "--dirty"), stdout = TRUE, stderr = FALSE),
    error = function(e) "unversioned"
  )
  if (!length(git_version)) git_version <- "unversioned"

  git_build_number <- tryCatch(
    system2("git", c("-C", repo_root, "rev-list", "--count", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) "NA"
  )
  if (!length(git_build_number)) git_build_number <- "NA"

  git_commit <- tryCatch(
    system2("git", c("-C", repo_root, "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) "unknown"
  )
  if (!length(git_commit)) git_commit <- "unknown"

  build_datetime <- format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S %Z",
    tz = "America/New_York"
  )

  lines <- c(
    sprintf("\\renewcommand{\\manuscriptversion}{%s}", escape_tex(git_version[[1L]])),
    sprintf("\\renewcommand{\\manuscriptbuildnumber}{%s}", escape_tex(git_build_number[[1L]])),
    sprintf("\\renewcommand{\\manuscriptcommit}{%s}", escape_tex(git_commit[[1L]])),
    sprintf("\\renewcommand{\\manuscriptbuilddatetime}{%s}", escape_tex(build_datetime))
  )
  writeLines(lines, build_info_tex, useBytes = TRUE)
}

write_build_info()

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
