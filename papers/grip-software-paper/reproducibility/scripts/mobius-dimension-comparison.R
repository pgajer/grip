#!/usr/bin/env Rscript
# Modern illustration inspired by Gajer et al. (2004), Section 3.6.
# Not a reproduction of the historical implementation, graphs, or projection.
# Source this file for the manuscript's 300-vertex pair, or run it with an
# output-directory argument to regenerate all six panels of Supplement S2.

mobius_graph <- function(longitudinal, transverse) {
  stopifnot(longitudinal >= 3L, transverse >= 2L)
  id <- function(i, j) (i - 1L) * transverse + j
  edges <- list()
  boundary <- logical()
  for (i in seq_len(longitudinal)) {
    for (j in seq_len(transverse)) {
      next.id <- if (i < longitudinal) id(i + 1L, j) else id(1L, transverse + 1L - j)
      edges[[length(edges) + 1L]] <- c(id(i, j), next.id)
      boundary <- c(boundary, j %in% c(1L, transverse))
      if (j < transverse) {
        edges[[length(edges) + 1L]] <- c(id(i, j), id(i, j + 1L))
        boundary <- c(boundary, FALSE)
      }
    }
  }
  edges <- do.call(rbind, edges)
  storage.mode(edges) <- "integer"
  # Reference used ONLY for a common display orientation, never for layout.
  u <- rep(2 * pi * (seq_len(longitudinal) - 1L) / longitudinal, each = transverse)
  half.width <- (transverse - 1L) * pi / longitudinal
  v <- rep(seq(-half.width, half.width, length.out = transverse), longitudinal)
  reference <- cbind((1 + v * cos(u / 2)) * cos(u),
                     (1 + v * cos(u / 2)) * sin(u), v * sin(u / 2))
  stopifnot(!any(edges[, 1L] == edges[, 2L]),
            nrow(unique(t(apply(edges, 1L, sort)))) == nrow(edges),
            nrow(edges) == longitudinal * (2L * transverse - 1L))
  list(edges = edges, boundary = boundary, reference = reference,
       longitudinal = longitudinal, transverse = transverse)
}

mobius_center_scale <- function(x) {
  x <- scale(x, center = TRUE, scale = FALSE)
  x / sqrt(mean(rowSums(x^2)))
}

mobius_orient <- function(x, reference) {
  x <- mobius_center_scale(x)
  fit <- svd(crossprod(x, mobius_center_scale(reference)))
  x %*% fit$u %*% t(fit$v)
}

mobius_settings <- function() {
  list(rounds = 160L, final_rounds = 256L, num_init = 24L,
       num_nbrs = 24L, seed = 1L)
}

mobius_compute <- function(longitudinal, transverse) {
  if (!requireNamespace("grip", quietly = TRUE) ||
      utils::packageVersion("grip") < "0.2.0") {
    stop("The Mobius comparison requires grip 0.2.0 or later.")
  }
  graph <- mobius_graph(longitudinal, transverse)
  n <- nrow(graph$reference)
  message("Computing n = ", n, ", dimensions 3 and 4; seed = 1")
  layouts <- lapply(3:4, function(d) do.call(grip::weighted.grip.nd,
    c(list(edges = graph$edges, n = n,
           edge_weights = rep(1, nrow(graph$edges)), dim = d), mobius_settings())))
  stopifnot(all(vapply(layouts, function(z) all(is.finite(z)), logical(1))))
  pca <- stats::prcomp(layouts[[2L]], center = TRUE, scale. = FALSE)
  projected <- pca$x[, 1:3, drop = FALSE]
  retained <- sum(pca$sdev[1:3]^2) / sum(pca$sdev^2)
  edge.lengths <- function(z) sqrt(rowSums(
    (z[graph$edges[, 1L], , drop = FALSE] -
       z[graph$edges[, 2L], , drop = FALSE])^2))
  # These ratios precede display normalization. Orthogonal projection cannot
  # increase a length, but high retained variance need not preserve local edges.
  projection.edge.ratio <- edge.lengths(projected) / edge.lengths(layouts[[2L]])
  stopifnot(all(is.finite(projection.edge.ratio)),
            all(projection.edge.ratio >= 0 & projection.edge.ratio <= 1 + 1e-10))
  list(graph = graph, layouts = layouts,
       projected = projected, projection = pca$rotation[, 1:3, drop = FALSE],
       retained.variance = retained, projection.edge.ratio = projection.edge.ratio,
       display = list(mobius_orient(layouts[[1L]], graph$reference),
                      mobius_orient(projected, graph$reference)))
}

mobius_panel <- function(result, method, show_n = TRUE, show_variance = TRUE,
                         plot_limit = 1.5) {
  xy <- grip::project.3d(result$display[[method]], azimuth = 35, elevation = 22)
  graph <- result$graph
  stopifnot(all(is.finite(xy)), max(abs(xy)) < plot_limit)
  plot(NA, xlim = c(-plot_limit, plot_limit),
       ylim = c(-plot_limit, plot_limit), asp = 1,
       axes = FALSE, xlab = "", ylab = "")
  label <- if (method == 1L) "Direct 3D" else "4D projected to 3D"
  if (show_n) label <- paste0(label, "; n = ", nrow(xy))
  title(label, cex.main = 0.95, line = 0.3)
  edges <- graph$edges
  segments(xy[edges[, 1L], 1L], xy[edges[, 1L], 2L],
           xy[edges[, 2L], 1L], xy[edges[, 2L], 2L], col = "gray60", lwd = 0.65)
  edges <- edges[graph$boundary, , drop = FALSE]
  segments(xy[edges[, 1L], 1L], xy[edges[, 1L], 2L],
           xy[edges[, 2L], 1L], xy[edges[, 2L], 2L], col = "#1F5A94", lwd = 1.2)
  if (show_variance && method == 2L) mtext(sprintf("PCA retains %.1f%% of 4D variance",
                                 100 * result$retained.variance),
                          side = 1, line = 0.2, cex = 0.8)
}

mobius_draw_pair <- function(result) {
  op <- par(mfrow = c(1, 2), mar = c(0.5, 0.3, 2.2, 0.3))
  on.exit(par(op))
  for (method in 1:2) {
    mobius_panel(result, method, show_n = FALSE, show_variance = FALSE,
                 plot_limit = 1.4)
  }
}

mobius_draw_all <- function(results) {
  op <- par(mfrow = c(2, length(results)), mar = c(1.8, 0.3, 2.2, 0.3),
            oma = c(0.5, 0, 0.5, 0))
  on.exit(par(op))
  par(cex = 1)
  for (method in 1:2) {
    for (result in results) mobius_panel(result, method)
  }
}

mobius_generate <- function(out) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  sizes <- rbind(c(30L, 5L), c(50L, 6L), c(100L, 15L))
  results <- lapply(seq_len(nrow(sizes)), function(k) {
    mobius_compute(sizes[k, 1L], sizes[k, 2L])
  })
  saveRDS(list(results = results, settings = mobius_settings(),
              session = utils::sessionInfo()),
          file.path(out, "mobius-dimension-comparison.rds"))
  grDevices::pdf(file.path(out, "mobius-dimension-comparison.pdf"), width = 10, height = 6.3)
  mobius_draw_all(results)
  grDevices::dev.off()
  grDevices::png(file.path(out, "mobius-dimension-comparison.png"),
                width = 1800, height = 1134, res = 180)
  mobius_draw_all(results)
  grDevices::dev.off()
  print(data.frame(n = sizes[, 1L] * sizes[, 2L],
    retained.variance = vapply(results, function(x) x$retained.variance, numeric(1)),
    min.edge.ratio = vapply(results, function(x) min(x$projection.edge.ratio), numeric(1))))
  invisible(results)
}

if (sys.nframe() == 0L) {
  paper.lib <- Sys.getenv("GRIP_RJOURNAL_PACKAGE_LIBRARY")
  if (nzchar(paper.lib)) .libPaths(c(paper.lib, .libPaths()))
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) stop("Supply one output-directory argument.")
  mobius_generate(args[[1L]])
}
