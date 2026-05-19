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

fixtures <- gisomap.dimension.known.fixtures(
  gflow_repo = "/Users/pgajer/current_projects/gflow"
)
sweeps <- gisomap.dimension.sweep.fixtures(
  fixtures = fixtures,
  initializers = c("metric_mds", "weighted_grip"),
  edge_kk_max_iter = 80L,
  weighted_rounds = 12L,
  seed = 511L,
  verbose = TRUE
)

out_dir <- file.path(repo, "output", "gisomap")
report_file <- file.path(out_dir, "gisomap_known_dimension_report_2026-05-19.html")
results_file <- file.path(out_dir, "gisomap_known_dimension_results_2026-05-19.csv")
summary_file <- file.path(out_dir, "gisomap_known_dimension_summary_2026-05-19.csv")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
results <- gisomap.dimension.bind.results(sweeps)
summary <- gisomap.dimension.expected_summary(results)
utils::write.csv(results, results_file, row.names = FALSE)
utils::write.csv(summary, summary_file, row.names = FALSE)
gisomap.dimension.report(
  sweeps = sweeps,
  file = report_file,
  title = "GISOMAP Known-Dimension Diagnostic Report",
  subtitle = "Quadform 2D/3D and Euclidean 4D/5D/6D fixtures with expected dimension markers"
)

message("Wrote results: ", results_file)
message("Wrote summary: ", summary_file)
message("Wrote report: ", report_file)
