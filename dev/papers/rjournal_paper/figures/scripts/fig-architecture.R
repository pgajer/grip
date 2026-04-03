## ---- Figure: Package architecture diagram (Section 3) -----------------------
##
## Three-layer architecture schematic drawn with base R graphics:
##   Layer 1: R front-end (APIs, presets, result packaging)
##   Layer 2: C++ back-end (MIS filtration, force-directed, Dijkstra/BFS)
##   Layer 3: Workflow layer (scoring, comparison, diagnostics, generators)
##
## This is a new figure -- no vignette equivalent exists.
##
## Output: figures/architecture.pdf

source("00-common.R")

save_pdf("architecture.pdf", width = 7.5, height = 5.5, {

  op <- par(mar = c(0.5, 0.5, 0.5, 0.5), bg = "white")

  plot.new()
  plot.window(xlim = c(0, 10), ylim = c(0, 8))

  ## Colours
  col_r    <- "#4A90D9"
  col_cpp  <- "#D35400"
  col_wf   <- "#27AE60"
  col_bg_r <- "#D6E9F8"
  col_bg_cpp <- "#FDEBD0"
  col_bg_wf  <- "#D5F5E3"
  col_border <- "gray30"

  ## ---- Layer boxes ----

  ## R front-end (top)
  rect(0.5, 5.5, 9.5, 7.5, col = col_bg_r, border = col_border, lwd = 1.5)
  text(5, 7.15, "R Front-end", font = 2, cex = 1.1, col = col_r)
  text(5, 6.65, "grip.layout()   grip.layout.weighted()   grip.layout.trace()", cex = 0.72)
  text(5, 6.25, "Preset resolution   Disconnected-component handling   Result packaging", cex = 0.65, col = "gray40")
  text(5, 5.85, "Edge-list / adjacency-list input   seed reproducibility", cex = 0.65, col = "gray40")

  ## C++ back-end (middle)
  rect(0.5, 2.8, 9.5, 5.1, col = col_bg_cpp, border = col_border, lwd = 1.5)
  text(5, 4.75, "C++ Back-end (Rcpp)", font = 2, cex = 1.1, col = col_cpp)
  text(5, 4.25, "MIS filtration   Coarse-to-fine placement   Force-directed refinement", cex = 0.72)
  text(5, 3.80, "BFS / Dijkstra shortest paths   Global repulsion   LGKK engine", cex = 0.72)
  text(5, 3.35, "Weighted MIS filtration   Weighted neighborhood caches", cex = 0.65, col = "gray40")
  text(5, 2.95, "Adaptive temperature control   Insertion-anchor strategies", cex = 0.65, col = "gray40")

  ## Workflow layer (bottom)
  rect(0.5, 0.3, 9.5, 2.4, col = col_bg_wf, border = col_border, lwd = 1.5)
  text(5, 2.05, "Workflow Layer", font = 2, cex = 1.1, col = col_wf)
  text(5, 1.55, "grip.score.layout()   grip.compare.layouts()   Trace diagnostics", cex = 0.72)
  text(5, 1.10, "Graph family generators (29 families)   Precomputed benchmarks", cex = 0.72)
  text(5, 0.65, "Composite scoring   Parameter grid search   Procrustes stability", cex = 0.65, col = "gray40")

  ## ---- Arrows between layers ----
  arrows(3, 5.5, 3, 5.15, length = 0.12, lwd = 2, col = "gray50")
  arrows(7, 5.5, 7, 5.15, length = 0.12, lwd = 2, col = "gray50")
  text(2.2, 5.32, "Rcpp", cex = 0.65, col = "gray50", font = 3)

  ## Workflow layer connects to both
  arrows(2, 2.8, 2, 2.45, length = 0.1, lwd = 1.5, col = "gray50", code = 3)
  arrows(8, 5.5, 8, 7.5,  length = 0.0, lwd = 0, col = NA)  ## invisible spacer

  ## Workflow <-> R front-end (scores any layout)
  arrows(8.5, 2.4, 8.5, 5.5, length = 0.1, lwd = 1.5, col = "gray50",
         code = 3, lty = 2)
  text(9.3, 3.95, "scores any\nlayout source", cex = 0.55, col = "gray50",
       font = 3)

  par(op)
})

cat("Wrote:", file.path(fig_dir, "architecture.pdf"), "\n")
