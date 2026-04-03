#!/usr/bin/env Rscript
## Figure 3: 8x12 torus in 3D with the torus preset, orthographic projection.

library(grip)

figure_output_path <- function(file.name) {
  args <- commandArgs(trailingOnly = FALSE)
  file.arg <- grep("^--file=", args, value = TRUE)
  script.dir <- if (length(file.arg)) {
    dirname(normalizePath(sub("^--file=", "", file.arg[1L])))
  } else {
    getwd()
  }
  out.dir <- file.path(normalizePath(file.path(script.dir, "..")), "figures")
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
  file.path(out.dir, file.name)
}

edges <- edges.torus(8, 12)
coords <- grip.layout(edges, n = max(edges), dim = 3,
                       preset = "torus", seed = 3)

out.file <- figure_output_path("fig3_torus.png")
png(out.file, width = 1100, height = 1000, res = 200)
grip.plot(coords, edges, projection = "ortho", main = "Torus (8x12)")
dev.off()
cat("Wrote ", out.file, "\n", sep = "")
