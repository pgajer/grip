#!/usr/bin/env Rscript
# Check the teaching script without complicating its four-step interface.
file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
here <- dirname(normalizePath(file))
code <- parse(file.path(here, "single-saddle-ivue.R"))

# Override only top-level settings in memory; never rewrite the script.
run.example <- function(settings, sample.only = FALSE) {
  env <- new.env(parent = globalenv())
  widgets <- list()
  for (expr in code) {
    assigned <- if (is.call(expr) && identical(expr[[1]], as.name("<-"))) {
      as.character(expr[[2]])
    } else ""
    if (sample.only && identical(assigned, "colors")) break
    if (assigned %in% names(settings)) expr[[3]] <- settings[[assigned]]
    value <- eval(expr, env)
    if (inherits(value, "htmlwidget")) widgets[[length(widgets) + 1L]] <- value
  }
  list(objects = env, widgets = widgets)
}

# The same seed reproduces the actual published sample, not just its equation.
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
cloud <- run.example(list(seed = 2211005), sample.only = TRUE)$objects$X
paper <- readRDS(file.path(here, "../../precomputed/two-fidelity-saddle.rds"))
stopifnot(paper$representative == 5, isTRUE(all.equal(cloud, paper$coords,
                                                  tolerance = 1e-14)))
cat("Published cloud-5 sampling reproduced.\n")

# Evaluate every plotting call without opening a browser or native window.
options(rgl.useNULL = TRUE)
result <- run.example(list(n = 80, max.iter = 10))
e <- result$objects
stopifnot(identical(names(e$graphs), as.character(3:20)),
          identical(names(e$mds), names(e$graphs)),
          identical(names(e$mds.edge.kk), names(e$graphs)),
          length(result$widgets) == 7L,
          max(abs(e$X[, 3] - e$C * (e$X[, 1]^2 - e$X[, 2]^2))) == 0)
for (w in result$widgets) {
  stopifnot(nrow(attr(w, "ivue")$X) == e$n)
}
for (key in names(e$graphs)) {
  g <- e$graphs[[key]]
  ig <- igraph::graph_from_edgelist(g$edge_matrix, directed = FALSE)
  input.lengths <- sqrt(rowSums((e$X[g$edge_matrix[, 1], ] -
                                  e$X[g$edge_matrix[, 2], ])^2))
  stopifnot(igraph::vcount(ig) == e$n, igraph::is_connected(ig),
            g$n_mst_edges_added == g$n_components_before - 1L,
            max(abs(input.lengths - g$edge_weight)) < 1e-12)
  for (z in list(e$mds[[key]], e$mds.edge.kk[[key]])) {
    stopifnot(identical(dim(z), c(80L, 3L)), all(is.finite(z)))
  }
}
cat("All 18 graphs, MDS/edge-KK fits, and seven ivue calls passed on n=80.\n")

# Verify that removing full-path preparation preserves the fitting problem.
for (key in c("3", "10", "20")) {
  g <- e$graphs[[key]]
  p <- grip::prepare.geodesic.kk(g$edge_matrix, n = e$n,
                               edge_weights = g$edge_weight, tie_mode = "single")
  full.mds <- grip::classical.mds(prepared = p, dim = 3, diagnostics = FALSE)$coords
  stopifnot(max(abs(as.vector(dist(full.mds)) -
                    as.vector(dist(e$mds[[key]])))) < 1e-8)
  full.kk <- grip::edge.kk(
    coords = e$mds[[key]], prepared = p, dim = 3, max_iter = e$max.iter,
    stiffness_method = "density", density_mix_schedule = c(0, .25, .5, .75, 1),
    scale_mode = "profiled", edge_length_epsilon = 0,
    diagnostics = FALSE, return_trace = FALSE, seed = e$seed
  )$coords
  stopifnot(max(abs(full.kk - e$mds.edge.kk[[key]])) < 1e-8)
}
cat("Lightweight fits agree with full-path preparation at k=3,10,20.\n")
