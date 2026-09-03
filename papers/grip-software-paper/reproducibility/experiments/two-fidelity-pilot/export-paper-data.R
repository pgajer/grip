#!/usr/bin/env Rscript
# Compact publication inputs; run after check-results.R and plot-results.R.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
input <- normalizePath(args[1], mustWork = TRUE)
output <- args[2]
lib <- Sys.getenv("GRIP_RJOURNAL_PACKAGE_LIBRARY")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
validation <- read.csv(file.path(input, "reference-validation.csv"))
fits <- lapply(1:5, function(r) readRDS(file.path(input,
  sprintf("fit-%02d-k%d.rds", r, validation$best_k[validation$replicate == r]))))
loss <- vapply(fits, function(f) f$scores$xg_error[1], numeric(1))
representative <- order(loss, seq_along(loss))[3]
f <- fits[[representative]]
nearest <- function(target, candidates = seq_len(nrow(f$coords))) {
  candidates[which.min(rowSums(sweep(f$coords[candidates, 1:2], 2, target, "-")^2))]
}
ends <- sort(c(nearest(c(-.8, 0)), nearest(c(.8, 0))))
p <- f$prepared
pid <- which(p$pair_matrix[, 1] == ends[1] & p$pair_matrix[, 2] == ends[2])
idx <- seq.int(p$flat_pair_edge_offsets[pid] + 1L, p$flat_pair_edge_offsets[pid + 1L])
route <- cbind(p$flat_edge_u[idx] + 1L, p$flat_edge_v[idx] + 1L)
graph_cases <- readRDS(file.path(input, sprintf("graphs-%02d.rds", f$replicate)))
sources <- readRDS(file.path(input, sprintf("cloud-%02d.rds", f$replicate)))$sources
control_ends <- c(nearest(c(-.8, 0), sources), nearest(c(.8, 0)))
controls <- lapply(c(3L, f$k), function(k) {
  g <- graph_cases$graphs[[as.character(k)]]
  ig <- igraph::graph_from_edgelist(g$edges, directed = FALSE)
  path <- as.integer(igraph::shortest_paths(ig, from = control_ends[1],
    to = control_ends[2], weights = g$weights)$vpath[[1]])
  list(k = k, edges = g$edges, weights = g$weights, bridges = g$bridges,
       mst_edges = g$mst_edges, route = path,
       distance = g$distances[control_ends[1], control_ends[2]])
})
paper <- list(schema_version = 1L, representative = representative,
  coords = f$coords, candidates = f$candidates, edges = f$graph$edges,
  weights = f$graph$weights, ends = ends, route = route, controls = controls,
  control_ends = control_ends, k = f$k,
  scores = do.call(rbind, lapply(fits, `[[`, "scores")),
  audit = do.call(rbind, lapply(fits, function(f) f$audit$scores)),
  curves = read.csv(file.path(input, "calibration-curves.csv")),
  validation = validation, status = read.csv(file.path(input, "fit-status.csv")),
  timings = read.csv(file.path(input, "layout-timings.csv")),
  independent_validation = read.csv(file.path(input, "independent-validation.csv")),
  provenance = list(source = "five-cloud pilot; see adjacent reproduction scripts",
                    selected_by = "median within-cloud minimum X-to-G loss",
                    sample_size = 1000L, embedding_dimension = 3L,
                    reference_pairs = 119744L, layout_pairs = 499500L,
                    export_session = sessionInfo()))
stopifnot(nrow(paper$scores) == 15L, nrow(paper$curves) == 390L,
          max(abs(paper$coords[, 3] - .8 * (paper$coords[, 1]^2 - paper$coords[, 2]^2))) < 1e-14,
          all(vapply(paper$candidates, ncol, integer(1)) == 3L),
          isTRUE(all.equal(paper$scores, read.csv(file.path(input, "layout-scores.csv")),
                           tolerance = 1e-12, check.attributes = FALSE)))
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
saveRDS(paper, output, compress = "xz")
stopifnot(identical(readRDS(output), paper))
cat(output, "\n", file.info(output)$size, "bytes\n")
