#!/usr/bin/env Rscript
## Figure 4: Sierpinski carpet level 3, default vs carpet preset side by side.

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

edges <- edges.sierpinski.carpet(3)
n <- max(edges)

carpet.default <- grip.layout(edges, n = n, dim = 2, seed = 12)
carpet.preset  <- grip.layout(edges, n = n, dim = 2, preset = "carpet", seed = 12)

out.file <- figure_output_path("fig4_carpet.png")
png(out.file, width = 2200, height = 1000, res = 200)
op <- par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
grip.plot(carpet.default, edges,
          main = "Carpet (default)",
          pch = 16, cex = 0.3,
          axes = FALSE, xlab = "", ylab = "", frame.plot = FALSE)
grip.plot(carpet.preset, edges,
          main = "Carpet (carpet preset)",
          pch = 16, cex = 0.3,
          axes = FALSE, xlab = "", ylab = "", frame.plot = FALSE)
par(op)
dev.off()
cat("Wrote ", out.file, "\n", sep = "")
