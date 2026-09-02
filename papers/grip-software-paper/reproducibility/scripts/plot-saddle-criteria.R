#!/usr/bin/env Rscript
# Three-panel explanatory adaptation of the weighted-saddle comparison.
# Source plot-weighted-saddle.R before using these functions from the manuscript.
# The CLI loads that helper and reads saved coordinates; no layout is recomputed.

saddle_criteria_data <- function(saddle, limits = NULL) {
  if (is.null(limits)) limits <- weighted_saddle_limits(list(saddle))
  methods <- c("Target saddle", "Metric MDS", "Weighted GRIP + edge-KK")
  layouts <- weighted_saddle_panel_layouts(saddle)[methods]
  projected <- weighted_saddle_projection(saddle)[methods]
  scores <- weighted_saddle_panel_scores(saddle)
  scores <- scores[match(methods, scores$method), , drop = FALSE]
  target <- saddle$target_coords
  # Choose the grid cross-section nearest y = 0 on its nonnegative side.
  # This is a geometric choice, independent of the candidate scores.
  slice_y <- min(target[target[, 2L] >= 0, 2L])
  slice <- which(abs(target[, 2L] - slice_y) < 1e-12)
  endpoints <- sort(slice[c(which.min(target[slice, 1L]),
                            which.max(target[slice, 1L]))])
  prepared <- saddle$prepared
  pair_index <- which(prepared$pair_matrix[, 1L] == endpoints[[1L]] &
                        prepared$pair_matrix[, 2L] == endpoints[[2L]])
  stopifnot(length(pair_index) == 1L)
  path <- prepared$path_edges[[pair_index]]
  stopifnot(nrow(path) > 1L,
            identical(as.integer(c(path[1L, 1L], path[nrow(path), 2L])),
                      as.integer(endpoints)),
            all(path[-nrow(path), 2L] == path[-1L, 1L]),
            all(as.vector(path) %in% slice))
  edge_lengths <- function(z, e) sqrt(rowSums(
    (z[e[, 1L], , drop = FALSE] - z[e[, 2L], , drop = FALSE])^2))
  edge_scale <- median(edge_lengths(target, saddle$edges))
  g <- prepared$pair_graph_distance[[pair_index]]
  stopifnot(abs(sum(edge_lengths(target, path)) / edge_scale - g) < 1e-10)
  lengths <- do.call(rbind, lapply(methods, function(method) {
    z <- layouts[[method]]
    data.frame(method = method, vertex_i = endpoints[[1L]],
               vertex_j = endpoints[[2L]], graph_distance = g,
               embedded_path_length = sum(edge_lengths(z, path)),
               endpoint_chord = sqrt(sum((z[endpoints[[1L]], ] -
                                           z[endpoints[[2L]], ])^2)))
  }))
  list(projected = projected, scores = scores, lengths = lengths,
       path = path, endpoints = endpoints, slice_y = slice_y,
       edges = saddle$edges, pair_count = nrow(prepared$pair_matrix),
       limits = limits)
}

plot_saddle_criteria <- function(data) {
  score_label <- function(value, path = FALSE) {
    label <- if (path) quote("GMDS " * plain(RelRMSE)^plain(path)) else "MDS Stress-1"
    if (value < 1e-12) return(bquote(.(label) < 10^{-12}))
    number <- if (value < 0.001) {
      parts <- strsplit(formatC(value, format = "e", digits = 2L), "e")[[1L]]
      bquote(.(parts[[1L]]) %*% 10^{.(as.integer(parts[[2L]]))})
    } else formatC(value, format = "fg", digits = 3L, flag = "#")
    bquote(.(label) == .(number))
  }
  old <- par(no.readonly = TRUE)
  on.exit(par(old))
  layout(rbind(1:3, rep(4L, 3L)), heights = c(1, 0.16))
  par(mar = c(2.9, 0.4, 2.3, 0.4), cex = 1, font.main = 1)
  titles <- c("Generating saddle", "Metric MDS", "Weighted GRIP + edge-KK")
  path_color <- "#2166AC"
  chord_color <- "#B34A36"
  for (i in seq_along(data$projected)) {
    z <- data$projected[[i]]
    plot(z, type = "n", asp = 1, axes = FALSE, xlab = "", ylab = "",
         xlim = data$limits$x, ylim = data$limits$y,
         main = titles[[i]], cex.main = 1.05)
    bounds <- par("usr")
    stopifnot(all(z[, 1L] > bounds[[1L]] & z[, 1L] < bounds[[2L]]),
              all(z[, 2L] > bounds[[3L]] & z[, 2L] < bounds[[4L]]))
    e <- data$edges
    segments(z[e[, 1L], 1L], z[e[, 1L], 2L],
             z[e[, 2L], 1L], z[e[, 2L], 2L], col = "gray75", lwd = 0.8)
    points(z, pch = 16, cex = 0.25, col = "gray55")
    e <- data$path
    segments(z[e[, 1L], 1L], z[e[, 1L], 2L],
             z[e[, 2L], 1L], z[e[, 2L], 2L], col = path_color, lwd = 2.7)
    p <- data$endpoints
    segments(z[p[[1L]], 1L], z[p[[1L]], 2L], z[p[[2L]], 1L], z[p[[2L]], 2L],
             col = chord_color, lwd = 1.8, lty = 2)
    points(z[p, , drop = FALSE], pch = 21, bg = "white", col = "gray20",
           cex = 0.85, lwd = 1)
    mtext(score_label(data$scores$mds.stress1[[i]]),
          side = 1, line = 0.1, cex = 0.93)
    mtext(score_label(data$scores$gmds.path.rel.rmse[[i]], path = TRUE),
          side = 1, line = 1.4, cex = 0.93)
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  legend("top", legend = c("Same retained input-graph path", "Endpoint chord"),
         col = c(path_color, chord_color), lty = c(1, 2), lwd = c(2.7, 1.8),
         horiz = TRUE, bty = "n", cex = 0.92, x.intersp = 0.7)
  text(0.5, 0.15, sprintf("Scores use all %s vertex pairs and full 3D coordinates.",
                         format(data$pair_count, big.mark = ",", trim = TRUE)),
       cex = 0.85)
  invisible(data)
}

saddle_criteria_caption <- function() paste(
  "Chord and fixed-path fidelity answer different questions on the same weighted graph.",
  "The generating coordinates lie on z = 0.8(x^2 - y^2) at a regular 10 by 10 grid",
  "in [-1,1]^2; horizontal and vertical neighbors are joined, with Euclidean edge",
  "lengths normalized by their median. Each panel highlights the same retained",
  "input-graph shortest path (blue) and its endpoint chord (dashed red); white circles",
  "mark the same two vertices. The annotations summarize all 4,950 pairs, not just",
  "the highlighted path, using separately profiled scales. The generating saddle",
  "preserves fixed-path lengths but has nonzero chord stress. Metric MDS reduces",
  "chord stress while changing fixed-path lengths. The edge-refined weighted-GRIP",
  "configuration has almost zero fixed-path loss despite its flatter shape.",
  "The configurations use the alignment and orthographic view of the complete",
  "method comparison, with common display limits across panels. These are graph-path criteria,",
  "not measurements of continuous-surface geodesics or a test of recovered curvature.")

if (sys.nframe() == 0L) {
  package_library <- Sys.getenv("GRIP_RJOURNAL_PACKAGE_LIBRARY")
  if (nzchar(package_library)) .libPaths(c(package_library, .libPaths()))
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
  source(file.path(script_dir, "plot-weighted-saddle.R"))
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) stop("Usage: plot-saddle-criteria.R BENCHMARK.rds OUTPUT_DIRECTORY")
  benchmark <- readRDS(args[[1L]])
  saddle <- benchmark$weighted_saddle
  stopifnot(saddle$n == 100L, saddle$m == 180L)
  data <- saddle_criteria_data(saddle, benchmark$weighted_saddle_resolutions$display_limits)
  dir.create(args[[2L]], recursive = TRUE, showWarnings = FALSE)
  stem <- file.path(args[[2L]], "saddle-three-criteria")
  pdf(paste0(stem, ".pdf"), width = 10.8, height = 4.8, useDingbats = FALSE)
  plot_saddle_criteria(data)
  dev.off()
  png(paste0(stem, ".png"), width = 10.8, height = 4.8, units = "in", res = 180)
  plot_saddle_criteria(data)
  dev.off()
  write.csv(cbind(data$scores, data$lengths[, -1L]),
            paste0(stem, "-scores.csv"), row.names = FALSE)
  write.csv(data$path, paste0(stem, "-path.csv"), row.names = FALSE)
  writeLines(saddle_criteria_caption(), paste0(stem, "-caption.txt"))
  print(data$scores)
  cat("Retained path endpoints:", data$endpoints, "; grid cross-section y =",
      data$slice_y, "\n")
}
