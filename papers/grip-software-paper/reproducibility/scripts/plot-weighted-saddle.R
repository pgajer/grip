# Common renderer for the manuscript and the 10x10 / 15x15 comparison figures.
# All coordinates are three-dimensional; only the display is projected to 2D.

weighted_saddle_projection <- function(saddle) {
  target <- saddle$target_coords
  align <- function(z) {
    centered <- sweep(z, 2L, colMeans(z), "-")
    reference <- sweep(target, 2L, colMeans(target), "-")
    fit <- svd(crossprod(centered, reference))
    scale <- sum(fit$d) / sum(centered^2)
    sweep(scale * centered %*% fit$u %*% t(fit$v), 2L, colMeans(target), "+")
  }
  z <- list(
    "Target saddle" = target,
    "Metric MDS" = align(saddle$layouts[["Metric MDS"]]),
    "Metric MDS + edge-KK" = align(saddle$layouts[["Metric MDS + edge-KK"]]),
    "Weighted KK (igraph)" = align(saddle$layouts[["Weighted KK (igraph)"]]),
    "Weighted GRIP" = align(saddle$layouts[["Weighted GRIP"]]),
    "Weighted GRIP + edge-KK" = align(saddle$layouts[["Weighted GRIP + edge-KK"]])
  )
  lapply(z, grip::project.3d, azimuth = 35, elevation = 22)
}

weighted_saddle_limits <- function(saddles) {
  projected <- lapply(saddles, weighted_saddle_projection)
  all <- do.call(rbind, unlist(projected, recursive = FALSE))
  padded <- function(x) range(x) + c(-1, 1) * 0.05 * diff(range(x))
  list(x = padded(all[, 1L]), y = padded(all[, 2L]))
}

plot_weighted_saddle <- function(saddle, limits = NULL) {
  if (is.null(limits)) limits <- weighted_saddle_limits(list(saddle))
  projected <- weighted_saddle_projection(saddle)
  edges <- saddle$edges
  colors <- c("#666666", "#684A88", "#684A88", "#3A6B35", "#1F3B73", "#1F3B73")
  old <- par(mfrow = c(2, 3), mar = c(1, 1, 3, 1), cex = 1.1)
  on.exit(par(old))
  for (i in seq_along(projected)) {
    z <- projected[[i]]
    plot(z[, 1L], z[, 2L], type = "n", asp = 1, axes = FALSE,
         xlab = "", ylab = "", xlim = limits$x, ylim = limits$y,
         main = names(projected)[[i]])
    segments(z[edges[, 1L], 1L], z[edges[, 1L], 2L],
             z[edges[, 2L], 1L], z[edges[, 2L], 2L], col = "gray65", lwd = 1.2)
    points(z[, 1L], z[, 2L], pch = 16, cex = 0.35, col = colors[[i]])
  }
  invisible(limits)
}
