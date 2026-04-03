## ---- Figure: 2D versus 3D weighted layout (Section 7.4) ---------------------
##
## Side-by-side weighted GRIP on a saddle mesh in 2D and 3D, with a score
## table showing the quality difference. Adapted from weighted-grip-intro.Rmd.
##
## Output: figures/2d-vs-3d.pdf

source("00-common.R")

## Build weighted saddle mesh (uses fallback if installed grip is older)
surface_mesh <- make_saddle_mesh(5, 5, amplitude = 0.9)

## Weighted GRIP in 2D
coords_2d <- grip.layout.weighted(
  surface_mesh$edges,
  n            = surface_mesh$n,
  edge_weights = surface_mesh$edge_weights,
  dim          = 2,
  preset       = "mesh",
  seed         = 2
)

## Weighted GRIP in 3D
coords_3d <- grip.layout.weighted(
  surface_mesh$edges,
  n            = surface_mesh$n,
  edge_weights = surface_mesh$edge_weights,
  dim          = 3,
  preset       = "mesh",
  seed         = 2
)

## Score both with geodesic KK (if available)
if (exists("grip.prepare.geodesic.kk", mode = "function")) {
  gkk_prepared <- grip.prepare.geodesic.kk(
    surface_mesh$edges,
    n            = surface_mesh$n,
    edge_weights = surface_mesh$edge_weights
  )

  score_2d <- grip.score.geodesic.kk(coords_2d, prepared = gkk_prepared)
  score_3d <- grip.score.geodesic.kk(coords_3d, prepared = gkk_prepared)

  cat("\n--- Dimension comparison (geodesic KK scores) ---\n")
  cat("2D  weighted RMSE:", score_2d$gkk.weighted.rmse,
      " mean rel error:", score_2d$gkk.mean.rel.path.error, "\n")
  cat("3D  weighted RMSE:", score_3d$gkk.weighted.rmse,
      " mean rel error:", score_3d$gkk.mean.rel.path.error, "\n")
} else {
  cat("grip.prepare.geodesic.kk not available; skipping GKK scoring.\n")
}

## Draw figure
save_pdf("2d-vs-3d.pdf", width = 10.2, height = 4.3, {
  op <- par(mfrow = c(1, 2), mar = c(1.2, 1.2, 3, 1.2), bg = "white")

  grip.plot(
    coords_2d,
    surface_mesh$edges,
    main       = "Weighted GRIP in 2D",
    vertex.col = "black",
    edge.col   = "gray82"
  )

  grip.plot(
    coords_3d,
    surface_mesh$edges,
    projection = "ortho",
    main       = "Weighted GRIP in 3D",
    vertex.col = "#1F3B73",
    edge.col   = "gray82"
  )

  par(op)
})

cat("Wrote:", file.path(fig_dir, "2d-vs-3d.pdf"), "\n")
