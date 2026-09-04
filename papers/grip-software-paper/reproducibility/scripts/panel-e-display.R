# Display-only preparation for Figure 6E and the subsection 3.2 example.
# Scoring must use the unaligned coordinates, not this returned configuration.
panel_e_display <- function(pilot, coords) {
  reference <- readRDS("reproducibility/precomputed/saddle-reference-diagnostics.rds")
  cloud <- reference$clouds[[as.character(pilot$representative)]]
  X <- pilot$coords
  stopifnot(isTRUE(all.equal(X, cloud$coords)), identical(dim(coords), dim(X)))
  centered <- sweep(coords, 2, colMeans(coords), "-")
  target <- sweep(X, 2, colMeans(X), "-")
  fit <- svd(crossprod(centered, target))
  aligned <- sweep((sum(fit$d) / sum(centered^2)) *
    centered %*% fit$u %*% t(fit$v), 2, colMeans(X), "+")
  # The common bounds are those used for all three published panels.
  bounds <- apply(rbind(do.call(rbind, cloud$aligned),
    c(-1, -1, -.8), c(1, 1, .8)), 2, range)
  bounds[1, ] <- bounds[1, ] - .05
  bounds[2, ] <- bounds[2, ] + .05
  corners <- as.matrix(expand.grid(bounds[, 1], bounds[, 2], bounds[, 3]))
  common <- ivue::layer3D.callback(function(ctx) {
    rgl::points3d(corners, col = "white", alpha = 0, size = 1)
  })
  colors <- grDevices::colorRampPalette(
    c("#173D65", "#86AFC4", "#D9B18B", "#8E4921"))(100)[
      pmin(100, pmax(1, 1 + floor((X[, 1] + 1) * 49.5)))]
  list(coords = aligned, triangles = cloud$triangles,
       colors = colors, bounds = common)
}
