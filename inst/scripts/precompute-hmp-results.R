#!/usr/bin/env Rscript

## Recompute the HMP/U01 candidate summaries and layouts used by the real-data
## vignette (not the current R Journal paper's uncoarsened UMB-HMP example).
##
## Set GRIP_HMP_RESULTS_OUTPUT to choose the output file. In a grip source
## checkout the default replaces the bundled package artifact. When run from an
## installed package or the journal supplement, the default is
## hmp_vignette_results.rds in the current working directory.

script_path <- if (!is.null(sys.frames()[[1]]$ofile)) {
  normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE)
} else {
  normalizePath(
    "inst/scripts/precompute-hmp-results.R",
    winslash = "/",
    mustWork = FALSE
  )
}

source(file.path(dirname(script_path), "portable-session-info.R"))

source_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  winslash = "/",
  mustWork = FALSE
)
description_path <- file.path(source_root, "DESCRIPTION")
is_source_checkout <- file.exists(description_path) &&
  identical(
    unname(read.dcf(description_path, fields = "Package")[[1L]]),
    "grip"
  )

data_rds <- Sys.getenv("GRIP_HMP_DATA_RDS", unset = "")

if (nzchar(data_rds)) {
  if (!file.exists(data_rds)) {
    stop("GRIP_HMP_DATA_RDS does not exist: ", data_rds)
  }
  suppressPackageStartupMessages(library(grip))
  hmp.u01.gc.coarse <- readRDS(data_rds)
} else if (is_source_checkout) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Package 'pkgload' is required when running from a source checkout.")
  }
  pkgload::load_all(
    source_root,
    quiet = TRUE,
    export_all = FALSE,
    helpers = FALSE
  )
  load(file.path(source_root, "data", "hmp.u01.gc.coarse.rda"))
} else {
  suppressPackageStartupMessages(library(grip))
  data("hmp.u01.gc.coarse", package = "grip", envir = environment())
}

default_output <- if (is_source_checkout) {
  file.path(
    source_root,
    "inst",
    "extdata",
    "hmp_u01_gc_coarse",
    "vignette_results.rds"
  )
} else {
  file.path(getwd(), "hmp_vignette_results.rds")
}
out.path <- Sys.getenv("GRIP_HMP_RESULTS_OUTPUT", unset = default_output)

vertex_data <- hmp.u01.gc.coarse$vertex_data
adj_list <- hmp.u01.gc.coarse$adj_list
weight_list <- hmp.u01.gc.coarse$weight_list
n <- length(adj_list)

preset_cmp <- compare.layouts(
  adj_list = adj_list,
  weight_list = weight_list,
  n = n,
  dim = 3L,
  candidates = c("default", "tree", "torus"),
  clusters = vertex_data$cst,
  seeds = 1L,
  sample.size.stress = 1200L,
  sample.size.nonedge = 3000L,
  return.layouts = TRUE
)

local_search_cmp <- compare.layouts(
  adj_list = adj_list,
  weight_list = weight_list,
  n = n,
  dim = 3L,
  search = list(
    candidate.prefix = "hmp.local",
    placement = "barycenter",
    rounds = 192L,
    final_rounds = c(224L, 288L),
    num_init = 30L,
    num_nbrs = 10L,
    r = 0.10,
    s = 2.0,
    repulsion_factor = c(0.75, 1.25)
  ),
  clusters = vertex_data$cst,
  seeds = 1L,
  sample.size.stress = 1200L,
  sample.size.nonedge = 3000L,
  return.layouts = TRUE
)

top_preset <- as.character(preset_cmp$summary$candidate[[1L]])
top_search_names <- as.character(
  utils::head(local_search_cmp$summary$candidate, 3L)
)

layout_store <- list(
  preset = lapply(
    c("default", "tree", "torus"),
    function(nm) preset_cmp$layouts[[nm]][["1"]]
  ),
  local = lapply(
    top_search_names,
    function(nm) local_search_cmp$layouts[[nm]][["1"]]
  )
)
names(layout_store$preset) <- c("default", "tree", "torus")
names(layout_store$local) <- top_search_names

out <- list(
  graph_label = "HMP/U01 coarsened giant component",
  graph_info = hmp.u01.gc.coarse$graph_info,
  preset_summary = preset_cmp$summary,
  preset_runs = preset_cmp$runs,
  local_search_summary = local_search_cmp$summary,
  local_search_runs = local_search_cmp$runs,
  top_preset = top_preset,
  top_local_candidates = top_search_names,
  layouts = layout_store,
  session_info = portable_session_info()
)

dir.create(dirname(out.path), recursive = TRUE, showWarnings = FALSE)
saveRDS(out, file = out.path, compress = "xz")
message("Saved HMP results to ", normalizePath(out.path, mustWork = TRUE))
