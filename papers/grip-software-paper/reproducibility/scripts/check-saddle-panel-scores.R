#!/usr/bin/env Rscript
# Check the Figure 8 annotations without regenerating any candidate layout.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
paper_library <- Sys.getenv("GRIP_RJOURNAL_PACKAGE_LIBRARY")
if (nzchar(paper_library)) .libPaths(c(paper_library, .libPaths()))
source(file.path(script_dir, "plot-weighted-saddle.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: check-saddle-panel-scores.R BENCHMARK.rds")
cases <- readRDS(args[[1L]])$weighted_saddle_resolutions$cases
stopifnot(length(cases) > 0L)
tolerance <- 1e-10
for (name in names(cases)) {
  saddle <- cases[[name]]
  scores <- weighted_saddle_panel_scores(saddle)
  layouts <- weighted_saddle_panel_layouts(saddle)
  prepared <- saddle$prepared
  g <- prepared$pair_graph_distance
  pairs <- prepared$pair_matrix
  for (method in names(layouts)) {
    z <- layouts[[method]]
    # Independent distance-matrix lookup and QR least-squares fit.
    distances <- as.matrix(stats::dist(z))
    chord <- distances[pairs]
    fit <- stats::lm.fit(matrix(g, ncol = 1L), chord)
    stress1 <- sqrt(sum(fit$residuals^2) / sum(chord^2))
    # Match the package's documented numerical edge-length smoothing.
    path <- vapply(prepared$path_edges, function(edges) {
      sum(sqrt(distances[edges]^2 + 1e-16))
    }, numeric(1L))
    path_fit <- stats::lm.fit(matrix(g, ncol = 1L), path)
    path_rel <- sqrt(sum(path_fit$residuals^2) / sum(path_fit$fitted.values^2))
    actual <- scores[scores$method == method, ]
    stopifnot(abs(stress1 - actual$mds.stress1) < tolerance,
              abs(path_rel - actual$gmds.path.rel.rmse) < tolerance)
    if (method != "Target saddle") {
      cached <- saddle$scores[saddle$scores$method == method, ]
      # The two chord statistics are related, but are not interchangeable.
      chord_rel <- cached$metric.chord.stress
      stopifnot(abs(actual$gmds.path.rel.rmse - cached$gmds.stress) < tolerance,
                abs(stress1 - chord_rel / sqrt(1 + chord_rel^2)) < tolerance,
                abs(stress1 - chord_rel) > tolerance)
    }
  }
  stopifnot(scores$gmds.path.rel.rmse[[1L]] < 1e-12,
            scores$mds.stress1[[1L]] > 0.1)
  # A display-style similarity transformation must not alter either score.
  transform <- function(z) sweep(3.7 * z %*% diag(c(-1, 1, 1)),
                                 2L, c(4, -2, 7), "+")
  moved <- saddle
  moved$target_coords <- transform(saddle$target_coords)
  moved$layouts <- lapply(saddle$layouts, transform)
  moved_scores <- weighted_saddle_panel_scores(moved)
  fields <- c("mds.stress1", "gmds.path.rel.rmse")
  stopifnot(max(abs(as.matrix(scores[fields]) -
                       as.matrix(moved_scores[fields]))) < tolerance)
  cat(name, ": all six panel scores match independent calculations and cached",
      "path results; similarity invariance passed.\n")
}

# A two-edge path provides exact zero/nonzero controls for both criteria.
straight <- cbind(0:2, 0, 0)
bent <- rbind(c(0, 0, 0), c(1, 0, 0), c(1, 1, 0))
prepared <- grip::prepare.geodesic.kk(rbind(c(1L, 2L), c(2L, 3L)), n = 3L)
control <- cases[[1L]]
control$n <- 3L
control$prepared <- prepared
control$target_coords <- straight
control$layouts <- lapply(control$layouts, function(z) bent)
scores <- weighted_saddle_panel_scores(control)
stopifnot(scores$mds.stress1[[1L]] < tolerance,
          all(scores$mds.stress1[-1L] > 0.1),
          all(scores$gmds.path.rel.rmse < tolerance))
cat("Straight and bent path controls passed.\n")
