#!/usr/bin/env Rscript

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to render grip-paper.Rmd.")
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

paper_candidates <- list(
  list(
    dir = file.path(
      repo_root,
      "dev",
      "papers",
      "rjournal_paper",
      "manuscript",
      "legacy_grip_paper_v3"
    ),
    input = "grip-paper-v3.Rmd"
  ),
  list(
    dir = file.path(
      repo_root,
      "dev",
      "papers",
      "rjournal_paper",
      "manuscript",
      "legacy_r_journal_drafts"
    ),
    input = "grip-paper.Rmd"
  )
)

selected <- NULL
for (candidate in paper_candidates) {
  if (file.exists(file.path(candidate$dir, candidate$input))) {
    selected <- candidate
    break
  }
}

if (is.null(selected)) {
  stop("Could not locate an R Journal manuscript source under dev/papers/rjournal_paper/manuscript/.")
}

paper_dir <- normalizePath(selected$dir, mustWork = TRUE)
input_file <- selected$input
build_dir <- file.path(paper_dir, "build")
build_info_tex <- file.path(build_dir, "manuscript_build_info.tex")

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(paper_dir)

args <- commandArgs(trailingOnly = TRUE)
render_pdf <- !("--html-only" %in% args)
render_html <- "--html" %in% args || "--all" %in% args
timestamped <- !("--no-timestamp" %in% args)

timestamp_suffix <- format(Sys.time(), "%Y%m%d_%H%M%S")
base_name <- tools::file_path_sans_ext(input_file)

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

  build_datetime <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

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
    return(file.path(paper_dir, stable_out))
  }

  stamped_out <- output_name(ext)
  if (file.exists(stamped_out)) {
    unlink(stamped_out)
  }

  ok <- file.rename(stable_out, stamped_out)
  if (!isTRUE(ok)) {
    ok <- file.copy(stable_out, stamped_out, overwrite = TRUE)
    if (!isTRUE(ok)) {
      stop("Failed to move ", stable_out, " to ", stamped_out)
    }
    unlink(stable_out)
  }

  file.path(paper_dir, stamped_out)
}

outputs <- character()

if (isTRUE(render_pdf)) {
  outputs <- c(outputs, render_one("rjtools::rjournal_pdf_article", "pdf"))
}

if (isTRUE(render_html)) {
  outputs <- c(outputs, render_one("rjtools::rjournal_article", "html"))
}

message("Rendered: ", paste(outputs, collapse = ", "))
