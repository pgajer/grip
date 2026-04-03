## ---- Figure: Combinatorial vs Weighted GRIP triptych (Section 7.1) ----------
##
## Three-panel figure: target surface geometry | combinatorial GRIP | weighted GRIP
## on a 5x5 saddle mesh. Adapted from weighted-grip-intro.Rmd.
##
## Output: figures/triptych-saddle.pdf

source("00-common.R")

## Create the weighted saddle mesh (uses fallback if installed grip is older)
surface_mesh <- make_saddle_mesh(5, 5, amplitude = 0.9)

## Combinatorial GRIP (ignores edge weights)
coords_comb <- grip.layout(
  surface_mesh$edges,
  n      = surface_mesh$n,
  dim    = 3,
  preset = "mesh",
  seed   = 1
)

## Weighted GRIP (respects edge-length geometry)
coords_wt <- grip.layout.weighted(
  surface_mesh$edges,
  n            = surface_mesh$n,
  edge_weights = surface_mesh$edge_weights,
  dim          = 3,
  preset       = "mesh",
  seed         = 1
)

## Compute geodesic KK scores for the caption/table (if available)
if (exists("grip.prepare.geodesic.kk", mode = "function")) {
  gkk_prepared <- grip.prepare.geodesic.kk(
    surface_mesh$edges,
    n            = surface_mesh$n,
    edge_weights = surface_mesh$edge_weights
  )

  score_comb <- grip.score.geodesic.kk(coords_comb, prepared = gkk_prepared)
  score_wt   <- grip.score.geodesic.kk(coords_wt,   prepared = gkk_prepared)

  cat("\n--- Geodesic KK scores ---\n")
  cat("Combinatorial GRIP  weighted RMSE:", score_comb$gkk.weighted.rmse, "\n")
  cat("Weighted GRIP       weighted RMSE:", score_wt$gkk.weighted.rmse, "\n")
} else {
  cat("grip.prepare.geodesic.kk not available; skipping GKK scoring.\n")
}

## Draw figure
save_pdf("triptych-saddle.pdf", width = 12.8, height = 4.2, {
  plot_layout_triptych(
    list(
      surface_mesh$coords_surface,
      coords_comb,
      coords_wt
    ),
    edges       = surface_mesh$edges,
    titles      = c("Target geometry", "Combinatorial GRIP", "Weighted GRIP"),
    projection  = "ortho",
    vertex_cols = c("#666666", "black", "#1F3B73")
  )
})

cat("Wrote:", file.path(fig_dir, "triptych-saddle.pdf"), "\n")
