## ---- Figure: Comparison with igraph/graphlayouts (Section 8) ----------------
##
## Uses precomputed benchmark results from inst/extdata/vs_alternatives/.
## Shows side-by-side layout panels for one or two benchmark graphs (mesh and
## Sierpinski carpet), plus a quality-metric bar chart.
##
## Adapted from grip-vs-alternatives.Rmd.
##
## Output: figures/comparison-panel.pdf, figures/comparison-metrics.pdf

source("00-common.R")

## Load precomputed benchmark results
res_path <- system.file(
  "extdata", "vs_alternatives", "benchmark_results.rds",
  package = "grip"
)
if (!nzchar(res_path)) {
  stop("Precomputed benchmark results not found. ",
       "Run inst/scripts/precompute-vs-alternatives.R first.")
}
res <- readRDS(res_path)

## ---- Figure A: Visual comparison on 12x12 mesh ----------------------------

m <- res$mesh

save_pdf("comparison-panel-mesh.pdf", width = 12, height = 7.5, {
  par(mfrow = c(2, 3), mar = c(1, 1, 3, 1), bg = "white")

  plot_layout_simple(m$layouts$fr,           m$edges, "FR (igraph)",           cex = 0.5)
  plot_layout_simple(m$layouts$kk,           m$edges, "KK (igraph)",           cex = 0.5)
  plot_layout_simple(m$layouts$drl,          m$edges, "DrL (igraph)",          cex = 0.5)
  plot_layout_simple(m$layouts$stress,       m$edges, "Stress (graphlayouts)", cex = 0.5)
  plot_layout_simple(m$layouts$grip.default, m$edges, "grip default",          cex = 0.5)
  plot_layout_simple(m$layouts$grip.mesh,    m$edges, "grip mesh preset",      cex = 0.5)
})

cat("Wrote:", file.path(fig_dir, "comparison-panel-mesh.pdf"), "\n")

## ---- Figure B: Visual comparison on level-4 Sierpinski carpet ---------------

c4 <- res$carpet

save_pdf("comparison-panel-carpet.pdf", width = 12, height = 7.5, {
  par(mfrow = c(2, 3), mar = c(1, 1, 3, 1), bg = "white")

  plot_layout_simple(c4$layouts$fr,           c4$edges, "FR (igraph)",           cex = 0.3)
  plot_layout_simple(c4$layouts$kk,           c4$edges, "KK (igraph)",           cex = 0.3)
  plot_layout_simple(c4$layouts$drl,          c4$edges, "DrL (igraph)",          cex = 0.3)
  plot_layout_simple(c4$layouts$stress,       c4$edges, "Stress (graphlayouts)", cex = 0.3)
  plot_layout_simple(c4$layouts$grip.default, c4$edges, "grip default",          cex = 0.3)
  plot_layout_simple(c4$layouts$grip.carpet,  c4$edges, "grip carpet preset",    cex = 0.3)
})

cat("Wrote:", file.path(fig_dir, "comparison-panel-carpet.pdf"), "\n")

## ---- Figure C: Quality metric comparison bar chart -------------------------
##
## Compact summary: sampled stress for mesh + carpet across all methods.

save_pdf("comparison-metrics.pdf", width = 9, height = 5, {
  op <- par(mfrow = c(1, 2), mar = c(8, 4, 3, 1), bg = "white")

  ## Mesh stress
  mesh_scores <- res$mesh$scores
  mesh_stress <- mesh_scores$sampled.stress
  names(mesh_stress) <- mesh_scores$method
  barplot(mesh_stress, las = 2, col = "#4A90D9",
          main = "12x12 Mesh: Sampled Stress",
          ylab = "Sampled stress (lower is better)",
          cex.names = 0.7)

  ## Carpet stress
  carpet_scores <- res$carpet$scores
  carpet_stress <- carpet_scores$sampled.stress
  names(carpet_stress) <- carpet_scores$method
  barplot(carpet_stress, las = 2, col = "#D35400",
          main = "Sierpinski Carpet: Sampled Stress",
          ylab = "Sampled stress (lower is better)",
          cex.names = 0.7)

  par(op)
})

cat("Wrote:", file.path(fig_dir, "comparison-metrics.pdf"), "\n")
