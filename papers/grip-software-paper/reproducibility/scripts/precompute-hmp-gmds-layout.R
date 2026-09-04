#!/usr/bin/env Rscript

## Precompute the weighted metric-MDS -> edge-KK layout used for the UMB-HMP
## real-data example in the R Journal paper. The all-pairs distance matrix is
## needed only during metric MDS and is deliberately omitted from the result.

if (!requireNamespace("grip", quietly = TRUE)) {
  stop("Install grip >= 0.2.0 before running this script.")
}
if (utils::packageVersion("grip") < "0.2.0") {
  stop("grip >= 0.2.0 is required.")
}

file.arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(file.arg)) stop("Run this file with Rscript.")
script.dir <- dirname(normalizePath(sub("^--file=", "", file.arg[[1L]])))
reproducibility.dir <- normalizePath(file.path(script.dir, ".."), mustWork = TRUE)

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) {
  args[[1L]]
} else {
  file.path(reproducibility.dir, "precomputed", "hmp_gmds", "hmp_gmds_layout.rds")
}

data("hmp.gc", package = "grip", envir = environment())
n <- length(hmp.gc$adj_list)

message("Computing weighted metric-MDS coordinates for ", n, " vertices...")
metric.time <- system.time({
  metric.fit <- grip::metric.mds(
    adj_list = hmp.gc$adj_list,
    weight_list = hmp.gc$weight_list,
    n = n,
    dim = 2L,
    diagnostics = FALSE
  )
})
metric.coords <- metric.fit$coords
rm(metric.fit)
invisible(gc())

message("Refining the metric-MDS coordinates with edge-KK...")
edge.prepared <- grip::prepare.edge.kk(
  adj_list = hmp.gc$adj_list,
  weight_list = hmp.gc$weight_list,
  n = n
)
metric.diagnostics <- grip::score.gmds(
  metric.coords,
  prepared = edge.prepared,
  scale_mode = "profiled"
)
edge.time <- system.time({
  edge.fit <- grip::edge.kk(
    coords = metric.coords,
    prepared = edge.prepared,
    dim = 2L,
    scale_mode = "profiled",
    max_iter = 50L,
    return_trace = FALSE,
    diagnostics = TRUE,
    seed = 1L,
    engine = "cpp"
  )
})

result <- list(
  n = n,
  m = hmp.gc$graph_info$edge_count,
  metric_mds = metric.coords,
  metric_mds_edge_kk = edge.fit$coords,
  metric_mds_diagnostics = metric.diagnostics,
  edge_kk_diagnostics = edge.fit$diagnostics,
  settings = list(
    graph_metric = "weighted shortest-path distance",
    dimension = 2L,
    initializer = "metric MDS",
    refinement = "edge-KK",
    edge_kk_max_iter = 50L,
    edge_kk_scale_mode = "profiled",
    seed = 1L
  ),
  elapsed_seconds = c(
    metric_mds = unname(metric.time[["elapsed"]]),
    edge_kk = unname(edge.time[["elapsed"]])
  ),
  package_version = as.character(utils::packageVersion("grip")),
  r_version = R.version.string
)

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
saveRDS(result, output, compress = "xz")
message("Saved ", normalizePath(output, winslash = "/", mustWork = TRUE))
print(result$elapsed_seconds)
print(result$metric_mds_diagnostics)
print(result$edge_kk_diagnostics)
