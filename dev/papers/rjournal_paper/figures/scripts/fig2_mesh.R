#!/usr/bin/env Rscript
## Figure 2: 8x8 mesh laid out with the mesh preset in 2D.

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

edges <- edges.mesh(8, 8)
coords <- grip.layout(edges, n = 64, dim = 2, preset = "mesh", seed = 1)

out.file <- figure_output_path("fig2_mesh.png")
png(out.file, width = 1000, height = 1000, res = 200)
grip.plot(coords, edges,
          pch = 16, cex = 0.6,
          axes = FALSE, xlab = "", ylab = "", frame.plot = FALSE,
          main = "8x8 mesh (mesh preset)")
dev.off()
cat("Wrote ", out.file, "\n", sep = "")
