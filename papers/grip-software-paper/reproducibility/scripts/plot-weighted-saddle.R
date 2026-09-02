# Common renderer for the manuscript and the 10x10 / 15x15 comparison figures.
# All coordinates are three-dimensional; only the display is projected to 2D.

weighted_saddle_panel_layouts <- function(saddle) {
  c(list("Target saddle" = saddle$target_coords), saddle$layouts[c(
    "Metric MDS", "Metric MDS + edge-KK", "Weighted KK (igraph)",
    "Weighted GRIP", "Weighted GRIP + edge-KK"
  )])
}

weighted_saddle_panel_scores <- function(saddle) {
  prepared <- saddle$prepared
  pairs <- prepared$pair_matrix
  g <- prepared$pair_graph_distance
  stopifnot(identical(prepared$tie_mode, "single"),
            nrow(pairs) == choose(saddle$n, 2L), all(is.finite(g)), all(g > 0))
  layouts <- weighted_saddle_panel_layouts(saddle)
  do.call(rbind, lapply(names(layouts), function(method) {
    z <- layouts[[method]]
    stopifnot(is.matrix(z), nrow(z) == saddle$n, ncol(z) == 3L,
              all(is.finite(z)))
    # Score raw 3D coordinates, not the display alignment or 2D projection.
    chord <- sqrt(rowSums((z[pairs[, 1L], , drop = FALSE] -
                            z[pairs[, 2L], , drop = FALSE])^2))
    chord_scale <- sum(chord * g) / sum(g^2)
    # Stress-1 has a configuration-distance denominator. The cached
    # metric.chord.stress instead uses scaled input distances in its denominator.
    stress1 <- sqrt(sum((chord - chord_scale * g)^2) / sum(chord^2))
    score <- grip::score.gmds(z, prepared = prepared, scale_mode = "profiled")
    stopifnot(is.finite(stress1), is.finite(score$gmds.stress))
    data.frame(method = method, mds.stress1 = stress1,
               gmds.path.rel.rmse = score$gmds.stress,
               chord.scale = chord_scale, path.scale = score$gmds.scale)
  }))
}

weighted_saddle_projection <- function(saddle) {
  target <- saddle$target_coords
  align <- function(z) {
    centered <- sweep(z, 2L, colMeans(z), "-")
    reference <- sweep(target, 2L, colMeans(target), "-")
    fit <- svd(crossprod(centered, reference))
    scale <- sum(fit$d) / sum(centered^2)
    sweep(scale * centered %*% fit$u %*% t(fit$v), 2L, colMeans(target), "+")
  }
  z <- weighted_saddle_panel_layouts(saddle)
  z[-1L] <- lapply(z[-1L], align)
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
  scores <- weighted_saddle_panel_scores(saddle)
  stopifnot(identical(names(projected), scores$method))
  score_label <- function(value, path = FALSE) {
    label <- if (path) quote("GMDS " * plain(RelRMSE)^plain(path)) else "MDS Stress-1"
    if (value < 1e-12) return(bquote(.(label) < 10^{-12}))
    number <- if (value < 0.001) {
      parts <- strsplit(formatC(value, format = "e", digits = 2L), "e")[[1L]]
      bquote(.(parts[[1L]]) %*% 10^{.(as.integer(parts[[2L]]))})
    } else formatC(value, format = "fg", digits = 3L, flag = "#")
    bquote(.(label) == .(number))
  }
  edges <- saddle$edges
  colors <- c("#666666", "#684A88", "#684A88", "#3A6B35", "#1F3B73", "#1F3B73")
  old <- par(mfrow = c(2, 3), mar = c(2.7, 0.3, 2.5, 0.3), cex = 1.1)
  on.exit(par(old))
  for (i in seq_along(projected)) {
    z <- projected[[i]]
    plot(z[, 1L], z[, 2L], type = "n", asp = 1, axes = FALSE,
         xlab = "", ylab = "", xlim = limits$x, ylim = limits$y,
         main = names(projected)[[i]])
    segments(z[edges[, 1L], 1L], z[edges[, 1L], 2L],
             z[edges[, 2L], 1L], z[edges[, 2L], 2L], col = "gray65", lwd = 1.2)
    points(z[, 1L], z[, 2L], pch = 16, cex = 0.35, col = colors[[i]])
    mtext(score_label(scores$mds.stress1[[i]]), side = 1, line = 0.15, cex = 1.05)
    mtext(score_label(scores$gmds.path.rel.rmse[[i]], path = TRUE),
          side = 1, line = 1.35, cex = 1.05)
  }
  invisible(limits)
}

# Re-export the saved comparisons without rerunning any layout or benchmark.
if (sys.nframe() == 0L) {
  paper_library <- Sys.getenv("GRIP_RJOURNAL_PACKAGE_LIBRARY")
  if (nzchar(paper_library)) .libPaths(c(paper_library, .libPaths()))
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) stop("Usage: plot-weighted-saddle.R BENCHMARK.rds OUTPUT_DIRECTORY")
  comparison <- readRDS(args[[1L]])$weighted_saddle_resolutions
  stopifnot(length(comparison$cases) > 0L)
  dir.create(args[[2L]], recursive = TRUE, showWarnings = FALSE)
  scores <- lapply(names(comparison$cases), function(name) {
    saddle <- comparison$cases[[name]]
    stem <- file.path(args[[2L]], paste0("weighted-saddle-", name))
    pdf(paste0(stem, ".pdf"), width = 11.2, height = 8.6, useDingbats = FALSE)
    plot_weighted_saddle(saddle, comparison$display_limits)
    dev.off()
    png(paste0(stem, ".png"), width = 11.2, height = 8.6, units = "in", res = 180)
    plot_weighted_saddle(saddle, comparison$display_limits)
    dev.off()
    cbind(mesh = name, weighted_saddle_panel_scores(saddle))
  })
  scores <- do.call(rbind, scores)
  write.csv(scores, file.path(args[[2L]], "weighted-saddle-panel-scores.csv"), row.names = FALSE)
  print(scores)
}
