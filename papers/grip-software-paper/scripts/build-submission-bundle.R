#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (!length(file_arg)) {
  stop("Could not determine the path to this script.")
}

script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
paper_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
repo_root <- normalizePath(file.path(paper_dir, "..", ".."), mustWork = TRUE)
description <- read.dcf(file.path(repo_root, "DESCRIPTION"))
package_version <- unname(description[1L, "Version"])

bundle_name <- paste0("gajer-ravel-grip-rjournal-", package_version)
submission_root <- file.path(repo_root, "output", "rjournal_paper", "submission")
bundle_dir <- file.path(submission_root, bundle_name)
zip_path <- file.path(submission_root, paste0(bundle_name, ".zip"))

if (dir.exists(bundle_dir)) {
  unlink(bundle_dir, recursive = TRUE, force = TRUE)
}
if (file.exists(zip_path)) {
  unlink(zip_path, force = TRUE)
}
dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  "grip-software-paper.Rmd",
  "grip-software-paper.R",
  "grip-software-paper.tex",
  "grip-software-paper.bib",
  "grip-software-paper.pdf",
  "grip-software-paper.html",
  "RJournal.sty",
  "RJwrapper.tex",
  "_Rpackages.txt",
  "citation_verification.html"
)

source_paths <- file.path(paper_dir, required_files)
missing <- required_files[!file.exists(source_paths)]
if (length(missing)) {
  stop(
    "Build the stable PDF and HTML before bundling. Missing: ",
    paste(missing, collapse = ", ")
  )
}

copied <- file.copy(source_paths, bundle_dir, overwrite = TRUE)
if (!all(copied)) {
  stop("Failed to copy one or more required manuscript files.")
}

copy_tree <- function(relative_path) {
  source <- file.path(paper_dir, relative_path)
  if (!dir.exists(source)) {
    stop("Missing required directory: ", source)
  }
  target_parent <- file.path(bundle_dir, dirname(relative_path))
  dir.create(target_parent, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(
    source,
    target_parent,
    recursive = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )
  if (!isTRUE(ok)) {
    stop("Failed to copy required directory: ", relative_path)
  }
}

copy_tree("reproducibility")
copy_tree("grip-software-paper_files")

motivation_dir <- file.path(bundle_dir, "motivation-letter")
dir.create(motivation_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.copy(
  file.path(paper_dir, "motivation-letter.md"),
  motivation_dir,
  overwrite = TRUE
)) {
  stop("Failed to copy the motivation letter.")
}

submission_scripts <- file.path(bundle_dir, "scripts")
dir.create(submission_scripts, recursive = TRUE, showWarnings = FALSE)
if (!file.copy(
  file.path(paper_dir, "scripts", "render-paper.R"),
  submission_scripts,
  overwrite = TRUE
)) {
  stop("Failed to copy the manuscript rendering script.")
}

build_info_source <- file.path(paper_dir, "build", "manuscript_build_info.tex")
if (file.exists(build_info_source)) {
  build_target <- file.path(bundle_dir, "build")
  dir.create(build_target, recursive = TRUE, showWarnings = FALSE)
  copied_build_info <- file.copy(
    build_info_source,
    build_target,
    overwrite = TRUE
  )
  if (!isTRUE(copied_build_info)) {
    stop("Failed to copy manuscript build information.")
  }
}

bundle_files <- list.files(
  bundle_dir,
  recursive = TRUE,
  all.files = TRUE,
  no.. = TRUE
)
bundle_files <- sort(c(bundle_files, "SUBMISSION_MANIFEST.txt"))
writeLines(
  c(
    paste0("Bundle: ", bundle_name),
    paste0("Package version: ", package_version),
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", tz = "America/New_York")),
    "",
    "Files:",
    bundle_files
  ),
  file.path(bundle_dir, "SUBMISSION_MANIFEST.txt"),
  useBytes = TRUE
)

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(submission_root)
status <- system2(
  "zip",
  c("-q", "-r", "-X", basename(zip_path), bundle_name)
)
if (!identical(status, 0L) || !file.exists(zip_path)) {
  stop("Failed to create submission ZIP: ", zip_path)
}

message("Submission directory: ", bundle_dir)
message("Submission ZIP: ", zip_path)
