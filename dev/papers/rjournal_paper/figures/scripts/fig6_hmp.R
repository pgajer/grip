#!/usr/bin/env Rscript
## Figure 6: HMP microbial network (1,828 vertices) colored by CST.
## Three panels: default, tree preset, torus preset.
##
## NOTE: This script uses the bundled precomputed vignette results if
## available. If not, it runs the full comparison (takes several minutes).

library(grip)

figure_output_path <- function(file.name) {
  args <- commandArgs(trailingOnly = FALSE)
  file.arg <- grep("^--file=", args, value = TRUE)
  script.dir <- if (length(file.arg)) {
    dirname(normalizePath(sub("^--file=", "", file.arg[1L])))
  } else {
    getwd()
  }
  repo.root <- normalizePath(file.path(script.dir, "../../../../.."))
  out.dir <- file.path(repo.root, "output", "rjournal_paper", "figures")
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
  file.path(out.dir, file.name)
}

data(hmp.u01.gc.coarse)

## --- build edge matrix from adjacency list ---
edge.matrix.from.adj <- function(adj.list) {
  edges <- list(); idx <- 0L
  for (u in seq_along(adj.list)) {
    nbrs <- adj.list[[u]][adj.list[[u]] > u]
    for (v in nbrs) { idx <- idx + 1L; edges[[idx]] <- c(u, v) }
  }
  do.call(rbind, edges)
}
hmp.edges <- edge.matrix.from.adj(hmp.u01.gc.coarse$adj_list)

## --- CST coloring ---
hmp.cst <- hmp.u01.gc.coarse$vertex_data$cst
cst.levels <- sort(unique(hmp.cst))
cst.cols <- setNames(
  grDevices::hcl.colors(length(cst.levels), "Dark 3"),
  cst.levels
)
hmp.cst.col <- cst.cols[hmp.cst]

## --- try precomputed results first ---
hmp.results <- tryCatch({
  path <- system.file("extdata", "hmp_u01_gc_coarse",
                       "vignette_results.rds", package = "grip")
  if (nzchar(path)) readRDS(path) else NULL
}, error = function(e) NULL)

if (!is.null(hmp.results)) {
  layouts <- list(
    default = hmp.results$layouts$preset$default,
    tree    = hmp.results$layouts$preset$tree,
    torus   = hmp.results$layouts$preset$torus
  )
} else {
  cat("Precomputed results not found; computing layouts (this may take minutes)...\n")
  hmp.n <- length(hmp.u01.gc.coarse$adj_list)

  layouts <- list(
    default = grip.layout(
      adj_list = hmp.u01.gc.coarse$adj_list,
      weight_list = hmp.u01.gc.coarse$weight_list,
      n = hmp.n, dim = 3, seed = 1
    ),
    tree = grip.layout(
      adj_list = hmp.u01.gc.coarse$adj_list,
      weight_list = hmp.u01.gc.coarse$weight_list,
      n = hmp.n, dim = 3, preset = "tree", seed = 1
    ),
    torus = grip.layout(
      adj_list = hmp.u01.gc.coarse$adj_list,
      weight_list = hmp.u01.gc.coarse$weight_list,
      n = hmp.n, dim = 3, preset = "torus", seed = 1
    )
  )
}

## --- figure ---
out.file <- figure_output_path("fig6_hmp.png")
png(out.file, width = 2800, height = 900, res = 200)
op <- par(mfrow = c(1, 3), mar = c(1, 1, 3, 1))
grip.plot(layouts$default, hmp.edges,
          projection = "ortho", vertex.col = hmp.cst.col, main = "default")
grip.plot(layouts$tree, hmp.edges,
          projection = "ortho", vertex.col = hmp.cst.col, main = "tree preset")
grip.plot(layouts$torus, hmp.edges,
          projection = "ortho", vertex.col = hmp.cst.col, main = "torus preset")
par(op)
dev.off()
cat("Wrote ", out.file, "\n", sep = "")
