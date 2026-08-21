script_path <- if (!is.null(sys.frames()[[1]]$ofile)) {
  normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE)
} else {
  normalizePath("data-raw/hmp_u01_gc_coarse_vignette_results.R", winslash = "/", mustWork = FALSE)
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)

suppressPackageStartupMessages({
  library(pkgload)
})

pkgload::load_all(repo_root, quiet = TRUE, export_all = FALSE, helpers = FALSE)
load(file.path(repo_root, "data", "hmp.u01.gc.coarse.rda"))

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
top_search_names <- as.character(utils::head(local_search_cmp$summary$candidate, 3L))

layout_store <- list(
  preset = lapply(c("default", "tree", "torus"), function(nm) preset_cmp$layouts[[nm]][["1"]]),
  local = lapply(top_search_names, function(nm) local_search_cmp$layouts[[nm]][["1"]])
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
  layouts = layout_store
)

saveRDS(
  out,
  file = file.path(repo_root, "inst", "extdata", "hmp_u01_gc_coarse", "vignette_results.rds"),
  compress = "xz"
)
