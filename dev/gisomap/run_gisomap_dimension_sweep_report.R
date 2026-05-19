#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
script_dir <- if (!is.na(script_file) && nzchar(script_file)) dirname(normalizePath(script_file)) else getwd()
repo <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
if (!file.exists(file.path(repo, "DESCRIPTION"))) {
  repo <- normalizePath(getwd(), mustWork = TRUE)
}
setwd(repo)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("pkgload is required to run this report script")
}

pkgload::load_all(repo, quiet = TRUE)
source(file.path(repo, "dev", "gisomap", "gisomap_dimension_sweep_helper.R"))

dims <- 2:5
sweeps <- gisomap.dimension.sweep.examples(
  dims = dims,
  initializers = c("metric_mds", "weighted_grip"),
  edge_kk_max_iter = 25L,
  weighted_rounds = 8L,
  seed = 101L,
  verbose = TRUE
)

out_dir <- file.path(repo, "output", "gisomap")
report_file <- file.path(out_dir, "gisomap_dimension_sweep_report_2026-05-19.html")
results_file <- file.path(out_dir, "gisomap_dimension_sweep_results_2026-05-19.csv")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(gisomap.dimension.bind.results(sweeps), results_file, row.names = FALSE)
gisomap.dimension.report(
  sweeps = sweeps,
  file = report_file,
  title = "GISOMAP Dimension Sweep Report",
  subtitle = paste0(
    "MDS -> edge-KK and weighted-GRIP -> edge-KK over dim = ",
    min(dims), ":", max(dims)
  )
)

message("Wrote results: ", results_file)
message("Wrote report: ", report_file)
