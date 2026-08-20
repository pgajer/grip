## ---- common helpers for all paper figure scripts ----------------------------
##
## Source this file at the top of every figure script:
##   source("00-common.R")
##
## It loads grip, sets the output directory, and defines shared helpers.

library(grip)

## Generated figures live outside dev/, which contains source files only.
script_dir <- dirname(normalizePath(sys.frame(1)$ofile %||% "."))
repo_root <- normalizePath(file.path(script_dir, "../../../../.."))
fig_dir <- file.path(repo_root, "output", "rjournal_paper", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## Convenience: save a PDF figure at R-Journal-friendly dimensions
save_pdf <- function(filename, width, height, expr) {
  path <- file.path(fig_dir, filename)
  pdf(path, width = width, height = height, family = "Helvetica")
  on.exit(dev.off(), add = TRUE)
  force(expr)
  invisible(path)
}

## Triptych helper (adapted from vignettes)
plot_layout_triptych <- function(coords_list,
                                 edges,
                                 titles,
                                 projection = NULL,
                                 vertex_cols = rep("black", length(coords_list)),
                                 edge_col = "gray82") {
  op <- par(
    mfrow = c(1, length(coords_list)),
    mar   = c(1.2, 1.2, 3, 1.2),
    bg    = "white"
  )
  on.exit(par(op), add = TRUE)

  for (i in seq_along(coords_list)) {
    grip.plot(
      coords_list[[i]], edges,
      projection = projection,
      main       = titles[[i]],
      vertex.col = vertex_cols[[i]],
      edge.col   = edge_col
    )
  }
}

## Side-by-side pair helper
plot_layout_pair <- function(coords_left, coords_right, edges, titles,
                             projection = NULL,
                             vertex_cols = c("black", "black"),
                             edge_col = "gray82") {
  op <- par(mfrow = c(1, 2), mar = c(1.2, 1.2, 3, 1.2), bg = "white")
  on.exit(par(op), add = TRUE)

  grip.plot(coords_left,  edges, projection = projection,
            main = titles[[1]], vertex.col = vertex_cols[[1]], edge.col = edge_col)
  grip.plot(coords_right, edges, projection = projection,
            main = titles[[2]], vertex.col = vertex_cols[[2]], edge.col = edge_col)
}

## Trace frame helper (adapted from grip-examples.Rmd)
plot_trace_frame <- function(coords, edges, main = "") {
  active <- stats::complete.cases(coords[, 1:2, drop = FALSE])
  xy <- coords[active, 1:2, drop = FALSE]

  xlim <- range(xy[, 1])
  ylim <- range(xy[, 2])
  xpad <- 0.08 * diff(xlim); if (!is.finite(xpad) || xpad == 0) xpad <- 0.2
  ypad <- 0.08 * diff(ylim); if (!is.finite(ypad) || ypad == 0) ypad <- 0.2

  plot(xy[, 1], xy[, 2], type = "n", asp = 1, axes = FALSE,
       xlab = "", ylab = "",
       xlim = xlim + c(-xpad, xpad),
       ylim = ylim + c(-ypad, ypad),
       main = main)

  active_edges <- edges[active[edges[, 1]] & active[edges[, 2]], , drop = FALSE]
  if (nrow(active_edges) > 0) {
    segments(coords[active_edges[, 1], 1], coords[active_edges[, 1], 2],
             coords[active_edges[, 2], 1], coords[active_edges[, 2], 2],
             col = "gray82")
  }
  points(xy[, 1], xy[, 2], pch = 16, cex = 0.55, col = "black")
}

## Simple layout plotter for comparison panels
plot_layout_simple <- function(coords, edges, title,
                               vertex_col = "#1F3B73",
                               edge_col = "gray80",
                               cex = 0.5) {
  plot(coords[, 1], coords[, 2],
       asp = 1, pch = 16, cex = cex,
       col = vertex_col, xlab = "", ylab = "",
       axes = FALSE, main = title)
  if (!is.null(edges) && nrow(edges) > 0) {
    segments(coords[edges[, 1], 1], coords[edges[, 1], 2],
             coords[edges[, 2], 1], coords[edges[, 2], 2],
             col = edge_col)
  }
}

## Null-coalescing operator if not available
`%||%` <- function(a, b) if (is.null(a)) b else a

## Fallback for mesh.surface.graph if the installed grip version is older
## than the source tree. Builds a mesh, lifts it onto a saddle surface,
## and computes edge weights from Euclidean distances in the embedding.
make_saddle_mesh <- function(h, w, amplitude = 0.9) {

  ## Try the real function first
  if (exists("mesh.surface.graph", mode = "function")) {
    return(mesh.surface.graph(h, w, surface = "saddle", amplitude = amplitude))
  }

  cat("mesh.surface.graph not found; using inline fallback.\n")

  edges <- edges.mesh(h, w)
  n <- h * w

  ## Parametric grid coordinates in [0,1]x[0,1]
  u <- rep(seq(0, 1, length.out = w), each = h)
  v <- rep(seq(0, 1, length.out = h), times = w)

  ## Saddle surface: z = amplitude * (u^2 - v^2), centered
  uc <- u - 0.5
  vc <- v - 0.5
  x <- uc
  y <- vc
  z <- amplitude * (uc^2 - vc^2)
  coords_surface <- cbind(x, y, z)

  ## Edge weights = Euclidean distance in 3D embedding, median-normalized
  d <- sqrt(rowSums((coords_surface[edges[, 1], ] - coords_surface[edges[, 2], ])^2))
  d <- d / median(d)

  list(
    edges          = edges,
    n              = n,
    edge_weights   = d,
    coords_surface = coords_surface
  )
}

cat("Common helpers loaded. fig_dir =", fig_dir, "\n")
