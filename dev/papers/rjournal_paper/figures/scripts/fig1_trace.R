#!/usr/bin/env Rscript
## Figure 1: Multiscale trace of a 6x6 mesh layout
##
## Shows four stages of the GRIP refinement process from coarsest
## initialization through progressive refinement to the final layout.

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

edges <- edges.mesh(6, 6)
tr <- grip.layout.trace(
  edges, n = 36, dim = 2,
  preset = "mesh", seed = 31,
  trace = "level"
)

frame.idx <- unique(c(
  1L,
  max(2L, floor(length(tr$frames) / 3)),
  max(3L, floor(2 * length(tr$frames) / 3)),
  length(tr$frames)
))

out.file <- figure_output_path("fig1_trace.png")
png(out.file, width = 2400, height = 700, res = 200)
op <- par(mfrow = c(1, 4), mar = c(1, 1, 2.5, 1))
labels <- c("Coarsest start", "Early refinement",
            "Late refinement", "Final layout")
for (i in seq_along(frame.idx)) {
  coords <- tr$frames[[frame.idx[i]]]
  active <- complete.cases(coords[, 1:2])
  xy <- coords[active, 1:2, drop = FALSE]
  active.edges <- edges[active[edges[, 1]] & active[edges[, 2]], , drop = FALSE]

  xlim <- range(xy[, 1]); ylim <- range(xy[, 2])
  xpad <- max(0.08 * diff(xlim), 0.2)
  ypad <- max(0.08 * diff(ylim), 0.2)

  plot(xy[, 1], xy[, 2], type = "n", asp = 1, axes = FALSE,
       xlab = "", ylab = "",
       xlim = xlim + c(-xpad, xpad),
       ylim = ylim + c(-ypad, ypad),
       main = labels[i])
  if (nrow(active.edges) > 0)
    segments(coords[active.edges[, 1], 1], coords[active.edges[, 1], 2],
             coords[active.edges[, 2], 1], coords[active.edges[, 2], 2],
             col = "gray80")
  points(xy[, 1], xy[, 2], pch = 16, cex = 0.6)
}
par(op)
dev.off()
cat("Wrote ", out.file, "\n", sep = "")
