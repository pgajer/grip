#!/usr/bin/env Rscript

run_tag <- "lgkk-round4-sets-2026-03-30"
manual_root <- file.path("output", "gkk_lgkk_paper")
tmp_dir <- file.path(manual_root, "tmp", run_tag)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'devtools' or the 'grip' package to run this benchmark.")
}

if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Package 'igraph' is required to run this benchmark.")
}

helper_env <- new.env(parent = baseenv())
sys.source(file.path("tools", "benchmarks", "gkk_lgkk_paper", "benchmark-sierpinski-baseline.R"), envir = helper_env)

fmt_num <- function(x, digits = 4L) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}

carpet_spec <- function(level = 4L) {
  built <- helper_env$build_sierpinski_carpet(level)
  list(
    id = sprintf("carpet_level_%d", level),
    label = sprintf("Sierpinski carpet level %d", level),
    family = "sierpinski.carpet",
    edges = built$edges,
    edge_weights = NULL,
    canonical = built$coords
  )
}

triangle_spec <- function(level = 4L) {
  built <- helper_env$build_sierpinski_triangle(level)
  list(
    id = sprintf("triangle_level_%d", level),
    label = sprintf("Sierpinski triangle level %d", level),
    family = "sierpinski.triangle",
    edges = built$edges,
    edge_weights = NULL,
    canonical = built$coords
  )
}

mesh_spec <- function(h = 8L, w = 8L) {
  edges <- edges.mesh(h, w)
  coords <- matrix(0, nrow = h * w, ncol = 2L)
  for (i in seq_len(h)) {
    for (j in seq_len(w)) {
      id <- (i - 1L) * w + j
      coords[id, ] <- c(j, h - i)
    }
  }
  list(
    id = sprintf("mesh_%dx%d", h, w),
    label = sprintf("Mesh %dx%d", h, w),
    family = "mesh",
    edges = edges,
    edge_weights = NULL,
    canonical = coords
  )
}

weighted_mesh_spec <- function(h = 6L, w = 6L, horizontal = 1, vertical = 2) {
  edges <- edges.mesh(h, w)
  row_of <- function(v) ((v - 1L) %/% w) + 1L
  col_of <- function(v) ((v - 1L) %% w) + 1L
  edge_weights <- apply(edges, 1L, function(e) {
    r1 <- row_of(e[[1L]])
    r2 <- row_of(e[[2L]])
    c1 <- col_of(e[[1L]])
    c2 <- col_of(e[[2L]])
    if (r1 == r2 && abs(c1 - c2) == 1L) horizontal else vertical
  })
  coords <- matrix(0, nrow = h * w, ncol = 2L)
  for (i in seq_len(h)) {
    for (j in seq_len(w)) {
      id <- (i - 1L) * w + j
      coords[id, ] <- c(j, (h - i) * vertical)
    }
  }
  list(
    id = sprintf("weighted_mesh_%dx%d", h, w),
    label = sprintf("Weighted mesh %dx%d", h, w),
    family = "mesh",
    edges = edges,
    edge_weights = as.double(edge_weights),
    canonical = coords
  )
}

score_layout <- function(spec, coords, seed = 1L) {
  geom <- grip.geometry.diagnostics(
    coords = coords,
    target.coords = spec$canonical,
    edges = spec$edges,
    family = spec$family,
    sample.size.symmetry = 512L,
    sample.size.wedges = 2048L,
    rng.seed = seed
  )
  lgkk <- grip.score.landmark.geodesic.kk(
    coords = coords,
    edges = spec$edges,
    edge_weights = spec$edge_weights,
    n = nrow(spec$canonical),
    local_nbrs = 6L,
    landmark_count = 8L
  )
  edge_stats <- helper_env$edge_length_stats(coords, spec$edges)
  data.frame(
    rmse = geom$procrustes.rmse[[1L]],
    lgkk_rel_rmse = lgkk$lgkk.weighted.rel.rmse[[1L]],
    edge_axis_deviation = geom$edge.axis.deviation[[1L]],
    edge_length_cv = edge_stats$cv,
    central_hole_skew = geom$central.hole.skew[[1L]],
    central_hole_aspect_error = geom$central.hole.aspect.error[[1L]],
    central_hole_center_error = geom$central.hole.center.error[[1L]],
    stringsAsFactors = FALSE
  )
}

evaluate_grip <- function(spec, label, args, seed = 1L) {
  full_args <- c(
    list(edges = spec$edges, n = nrow(spec$canonical), dim = 2, seed = seed),
    if (!is.null(spec$edge_weights)) list(edge_weights = spec$edge_weights) else list(),
    args
  )
  timed <- system.time(coords <- do.call(grip.layout, full_args))
  metrics <- score_layout(spec, coords, seed = seed)
  metrics$graph <- spec$id
  metrics$graph_label <- spec$label
  metrics$label <- label
  metrics$elapsed_sec <- unname(timed[["elapsed"]])
  list(coords = coords, metrics = metrics)
}

evaluate_trace_variant <- function(spec, label, args, seed = 1L) {
  tr <- do.call(
    grip.layout.trace,
    c(
      list(edges = spec$edges, n = nrow(spec$canonical), dim = 2,
           trace = "round", trace.every = 1L, seed = seed),
      if (!is.null(spec$edge_weights)) list(edge_weights = spec$edge_weights) else list(),
      args
    )
  )
  first_full_idx <- which(tr$meta$active_vertices == nrow(spec$canonical))[1L]
  first_coords <- tr$frames[[first_full_idx]]
  final_coords <- tr$final

  first_metrics <- score_layout(spec, first_coords, seed = seed)
  first_metrics$graph <- spec$id
  first_metrics$graph_label <- spec$label
  first_metrics$label <- label
  first_metrics$stage <- "first_full"
  first_metrics$elapsed_sec <- NA_real_

  final_metrics <- score_layout(spec, final_coords, seed = seed)
  final_metrics$graph <- spec$id
  final_metrics$graph_label <- spec$label
  final_metrics$label <- label
  final_metrics$stage <- "final"
  final_metrics$elapsed_sec <- NA_real_

  list(
    trace = tr,
    first_full = first_coords,
    final = final_coords,
    metrics = rbind(first_metrics, final_metrics)
  )
}

evaluate_igraph_kk <- function(spec, seed = 1L) {
  graph <- igraph::graph_from_edgelist(as.matrix(spec$edges), directed = FALSE)
  timed <- system.time({
    coords <- if (is.null(spec$edge_weights)) {
      igraph::layout_with_kk(graph)
    } else {
      igraph::layout_with_kk(graph, weights = spec$edge_weights)
    }
  })
  metrics <- score_layout(spec, coords, seed = seed)
  metrics$graph <- spec$id
  metrics$graph_label <- spec$label
  metrics$label <- "igraph::KK"
  metrics$elapsed_sec <- unname(timed[["elapsed"]])
  list(coords = coords, metrics = metrics)
}

plot_panel_aligned <- function(spec, coords, title, subtitle) {
  aligned <- helper_env$align_to_target(coords, spec$canonical)$aligned
  helper_env$plot_layout_panel(aligned, spec$edges, title, subtitle)
}

rank_sum <- function(df, metric_cols) {
  if (nrow(df) == 0L) {
    return(df)
  }
  score <- rep(0, nrow(df))
  for (graph_id in unique(df$graph)) {
    rows <- df$graph == graph_id
    for (metric in metric_cols) {
      x <- df[[metric]][rows]
      score[rows] <- score[rows] + rank(x, ties.method = "min", na.last = "keep")
    }
  }
  score
}

choose_set_a_winner <- function(set_a_metrics) {
  cand <- subset(set_a_metrics, label %in% c("A2", "A3"))
  cand$score <- rank_sum(
    cand,
    c("rmse", "lgkk_rel_rmse", "edge_axis_deviation", "elapsed_sec", "central_hole_skew")
  )
  cand <- cand[order(cand$score, cand$label), ]
  as.character(cand$label[[1L]])
}

choose_set_b_winner <- function(set_b_metrics) {
  cand <- subset(set_b_metrics, label %in% c("B4", "B5") & stage == "first_full")
  cand <- cand[order(cand$central_hole_skew,
                     cand$rmse,
                     cand$lgkk_rel_rmse,
                     cand$edge_axis_deviation,
                     cand$label), ]
  as.character(cand$label[[1L]])
}

write_lines <- function(path, lines) {
  writeLines(lines, con = path)
}

carpet4 <- carpet_spec(4L)
weighted_mesh <- weighted_mesh_spec(6L, 6L, horizontal = 1, vertical = 2)
carpet3 <- carpet_spec(3L)
triangle4 <- triangle_spec(4L)
mesh8 <- mesh_spec(8L, 8L)

set_a_configs <- list(
  baseline = list(label = "Baseline", args = list()),
  shared_x4 = list(
    label = "Shared x4",
    args = list(
      lgkk_multiscale_rounds = 4L,
      lgkk_local_nbrs = 6L,
      lgkk_landmark_count = 8L,
      lgkk_multiscale_scope = "all",
      lgkk_active_limit = 4096L
    )
  ),
  A2 = list(
    label = "A2",
    args = list(
      lgkk_multiscale_rounds = 0L,
      lgkk_rounds_coarse = 0L,
      lgkk_rounds_pre_final = 2L,
      lgkk_rounds_final = 4L,
      lgkk_local_nbrs = 6L,
      lgkk_landmark_count = 8L,
      lgkk_multiscale_scope = "all",
      lgkk_active_limit = 4096L
    )
  ),
  A3 = list(
    label = "A3",
    args = list(
      lgkk_multiscale_rounds = 0L,
      lgkk_rounds_coarse = 1L,
      lgkk_rounds_pre_final = 2L,
      lgkk_rounds_final = 4L,
      lgkk_local_nbrs = 6L,
      lgkk_landmark_count = 8L,
      lgkk_multiscale_scope = "all",
      lgkk_active_limit = 4096L
    )
  ),
  post_x4 = list(
    label = "R polish x4",
    args = list(
      lgkk_polish_rounds = 4L,
      lgkk_local_nbrs = 6L,
      lgkk_landmark_count = 8L
    )
  )
)

set_b_configs <- list(
  baseline = list(label = "Baseline", args = list()),
  B4 = list(
    label = "B4",
    args = list(
      insertion_anchor_count = 6L,
      insertion_anchor_scope = "prev_misf",
      insertion_anchor_strategy = "spread_prev",
      level0_insertion_mode = "least_squares",
      level0_local_kk_steps = 0L
    )
  ),
  B5 = list(
    label = "B5",
    args = list(
      insertion_anchor_count = 6L,
      insertion_anchor_scope = "prev_misf",
      insertion_anchor_strategy = "spread_prev",
      level0_insertion_mode = "least_squares",
      level0_local_kk_steps = 1L
    )
  )
)

set_a_primary <- list(carpet4, weighted_mesh)
set_a_results <- lapply(set_a_primary, function(spec) {
  out <- lapply(set_a_configs, function(cfg) {
    evaluate_grip(spec, cfg$label, cfg$args, seed = 1L)
  })
  names(out) <- names(set_a_configs)
  out
})
names(set_a_results) <- vapply(set_a_primary, `[[`, character(1L), "id")
set_a_metrics <- do.call(
  rbind,
  unlist(lapply(set_a_results, function(x) lapply(x, `[[`, "metrics")), recursive = FALSE)
)
rownames(set_a_metrics) <- NULL
set_a_winner <- choose_set_a_winner(set_a_metrics)

set_b_result <- lapply(set_b_configs, function(cfg) {
  evaluate_trace_variant(carpet4, cfg$label, cfg$args, seed = 1L)
})
names(set_b_result) <- names(set_b_configs)
set_b_metrics <- do.call(rbind, lapply(set_b_result, `[[`, "metrics"))
rownames(set_b_metrics) <- NULL
set_b_winner <- choose_set_b_winner(set_b_metrics)

winner_a_args <- set_a_configs[[names(which(vapply(set_a_configs, function(x) identical(x$label, set_a_winner), logical(1L))))]]$args
winner_b_args <- set_b_configs[[names(which(vapply(set_b_configs, function(x) identical(x$label, set_b_winner), logical(1L))))]]$args

set_c_configs <- list(
  baseline = list(label = "Baseline", args = list()),
  winner_a = list(label = paste0(set_a_winner, " only"), args = winner_a_args),
  winner_b = list(label = paste0(set_b_winner, " only"), args = winner_b_args),
  combined = list(
    label = "Combined winner",
    args = c(winner_a_args, winner_b_args)
  ),
  post_x4 = list(
    label = "R polish x4",
    args = list(
      lgkk_polish_rounds = 4L,
      lgkk_local_nbrs = 6L,
      lgkk_landmark_count = 8L
    )
  )
)

set_c_specs <- list(carpet4, weighted_mesh, carpet3, triangle4, mesh8)
set_c_results <- lapply(set_c_specs, function(spec) {
  out <- lapply(set_c_configs, function(cfg) {
    evaluate_grip(spec, cfg$label, cfg$args, seed = 1L)
  })
  out$igraph_kk <- evaluate_igraph_kk(spec, seed = 1L)
  out
})
names(set_c_results) <- vapply(set_c_specs, `[[`, character(1L), "id")
set_c_metrics <- do.call(
  rbind,
  unlist(lapply(set_c_results, function(x) lapply(x, `[[`, "metrics")), recursive = FALSE)
)
rownames(set_c_metrics) <- NULL

write.csv(set_a_metrics, file.path(tmp_dir, "set-a-primary-metrics.csv"), row.names = FALSE)
write.csv(set_b_metrics, file.path(tmp_dir, "set-b-carpet-trace-metrics.csv"), row.names = FALSE)
write.csv(set_c_metrics, file.path(tmp_dir, "set-c-confirmation-metrics.csv"), row.names = FALSE)

set_a_summary <- file.path(tmp_dir, "set-a-summary.md")
set_b_summary <- file.path(tmp_dir, "set-b-summary.md")
set_c_summary <- file.path(tmp_dir, "set-c-summary.md")

write_lines(
  set_a_summary,
  c(
    "# LGKK Round 4 Set A Summary",
    "",
    "- Primary graphs: level-4 carpet, weighted mesh 6x6",
    sprintf("- Winner: `%s`", set_a_winner),
    "",
    "## Metrics",
    "",
    capture.output(print(set_a_metrics, row.names = FALSE))
  )
)

write_lines(
  set_b_summary,
  c(
    "# LGKK Round 4 Set B Summary",
    "",
    "- Primary graph: level-4 carpet",
    sprintf("- Winner: `%s`", set_b_winner),
    "",
    "## Metrics",
    "",
    capture.output(print(set_b_metrics, row.names = FALSE))
  )
)

write_lines(
  set_c_summary,
  c(
    "# LGKK Round 4 Set C Summary",
    "",
    sprintf("- Set A winner: `%s`", set_a_winner),
    sprintf("- Set B winner: `%s`", set_b_winner),
    "- Graphs: carpet level 4, weighted mesh 6x6, carpet level 3, triangle level 4, mesh 8x8",
    "",
    "## Metrics",
    "",
    capture.output(print(set_c_metrics, row.names = FALSE))
  )
)

grDevices::png(file.path(tmp_dir, "set-a-carpet-contact-sheet.png"),
               width = 2400, height = 1600, res = 180, bg = "#f7f3ea")
graphics::par(mfrow = c(2, 3), mar = c(0, 0, 3.2, 0))
plot_panel_aligned(carpet4, carpet4$canonical, "Canonical", "reference")
for (key in c("baseline", "shared_x4", "A2", "A3", "post_x4")) {
  row <- subset(set_a_metrics, graph == carpet4$id & label == set_a_configs[[key]]$label)
  plot_panel_aligned(
    carpet4,
    set_a_results[[carpet4$id]][[key]]$coords,
    set_a_configs[[key]]$label,
    sprintf("RMSE %s | skew %s | %.1fs",
            fmt_num(row$rmse[[1L]]),
            fmt_num(row$central_hole_skew[[1L]]),
            row$elapsed_sec[[1L]])
  )
}
grDevices::dev.off()

grDevices::png(file.path(tmp_dir, "set-a-weighted-mesh-contact-sheet.png"),
               width = 2400, height = 1600, res = 180, bg = "#f7f3ea")
graphics::par(mfrow = c(2, 3), mar = c(0, 0, 3.2, 0))
plot_panel_aligned(weighted_mesh, weighted_mesh$canonical, "Canonical", "reference")
for (key in c("baseline", "shared_x4", "A2", "A3", "post_x4")) {
  row <- subset(set_a_metrics, graph == weighted_mesh$id & label == set_a_configs[[key]]$label)
  plot_panel_aligned(
    weighted_mesh,
    set_a_results[[weighted_mesh$id]][[key]]$coords,
    set_a_configs[[key]]$label,
    sprintf("RMSE %s | LGKK %s | %.3fs",
            fmt_num(row$rmse[[1L]]),
            fmt_num(row$lgkk_rel_rmse[[1L]]),
            row$elapsed_sec[[1L]])
  )
}
grDevices::dev.off()

grDevices::png(file.path(tmp_dir, "set-b-carpet-contact-sheet.png"),
               width = 2800, height = 1600, res = 180, bg = "#f7f3ea")
graphics::par(mfrow = c(2, 4), mar = c(0, 0, 3.2, 0))
plot_panel_aligned(carpet4, carpet4$canonical, "Canonical", "reference")
for (key in c("baseline", "B4", "B5")) {
  first_row <- subset(set_b_metrics, label == set_b_configs[[key]]$label & stage == "first_full")
  final_row <- subset(set_b_metrics, label == set_b_configs[[key]]$label & stage == "final")
  plot_panel_aligned(
    carpet4,
    set_b_result[[key]]$first_full,
    paste(set_b_configs[[key]]$label, "| first full"),
    sprintf("RMSE %s | skew %s",
            fmt_num(first_row$rmse[[1L]]),
            fmt_num(first_row$central_hole_skew[[1L]]))
  )
  plot_panel_aligned(
    carpet4,
    set_b_result[[key]]$final,
    paste(set_b_configs[[key]]$label, "| final"),
    sprintf("RMSE %s | skew %s",
            fmt_num(final_row$rmse[[1L]]),
            fmt_num(final_row$central_hole_skew[[1L]]))
  )
}
igraph_carpet <- set_c_results[[carpet4$id]]$igraph_kk
igraph_row <- subset(set_c_metrics, graph == carpet4$id & label == "igraph::KK")
plot_panel_aligned(
  carpet4,
  igraph_carpet$coords,
  "igraph::KK",
  sprintf("RMSE %s | skew %s",
          fmt_num(igraph_row$rmse[[1L]]),
          fmt_num(igraph_row$central_hole_skew[[1L]]))
)
grDevices::dev.off()

grDevices::png(file.path(tmp_dir, "set-b-carpet-metric-plot.png"),
               width = 2200, height = 1200, res = 180, bg = "#f7f3ea")
graphics::par(mfrow = c(1, 2), mar = c(8, 4, 3, 1))
plot_df <- subset(set_b_metrics, stage %in% c("first_full", "final"))
plot_df$label <- factor(plot_df$label, levels = c("Baseline", "B4", "B5"))
stage_colors <- c(first_full = "#355070", final = "#b23a48")
for (metric_name in c("rmse", "central_hole_skew")) {
  ymax <- max(plot_df[[metric_name]], na.rm = TRUE) * 1.12
  graphics::plot(c(1, 3), c(0, ymax), type = "n",
                 axes = FALSE, xlab = "", ylab = metric_name,
                 main = if (metric_name == "rmse") "Set B RMSE" else "Set B central-hole skew")
  graphics::axis(2)
  graphics::axis(1, at = 1:3, labels = levels(plot_df$label), las = 2)
  for (i in seq_along(levels(plot_df$label))) {
    rows <- plot_df$label == levels(plot_df$label)[[i]]
    rf <- rows & plot_df$stage == "first_full"
    rl <- rows & plot_df$stage == "final"
    graphics::points(i - 0.1, plot_df[[metric_name]][rf], pch = 16, cex = 1.4, col = stage_colors[["first_full"]])
    graphics::points(i + 0.1, plot_df[[metric_name]][rl], pch = 17, cex = 1.4, col = stage_colors[["final"]])
    graphics::segments(i - 0.1, plot_df[[metric_name]][rf], i + 0.1, plot_df[[metric_name]][rl],
                       col = grDevices::adjustcolor("#444444", alpha.f = 0.7), lwd = 1.5)
  }
}
grDevices::dev.off()

make_set_c_sheet <- function(spec, path, include_igraph = TRUE) {
  panels <- list(
    list(coords = spec$canonical, title = "Canonical", subtitle = "reference")
  )
  for (key in c("baseline", "winner_a", "winner_b", "combined", "post_x4")) {
    row <- subset(set_c_metrics, graph == spec$id & label == set_c_configs[[key]]$label)
    panels[[length(panels) + 1L]] <- list(
      coords = set_c_results[[spec$id]][[key]]$coords,
      title = set_c_configs[[key]]$label,
      subtitle = sprintf("RMSE %s | LGKK %s",
                         fmt_num(row$rmse[[1L]]),
                         fmt_num(row$lgkk_rel_rmse[[1L]]))
    )
  }
  if (include_igraph) {
    row <- subset(set_c_metrics, graph == spec$id & label == "igraph::KK")
    panels[[length(panels) + 1L]] <- list(
      coords = set_c_results[[spec$id]]$igraph_kk$coords,
      title = "igraph::KK",
      subtitle = sprintf("RMSE %s | LGKK %s",
                         fmt_num(row$rmse[[1L]]),
                         fmt_num(row$lgkk_rel_rmse[[1L]]))
    )
  }

  n_panels <- length(panels)
  n_col <- if (n_panels <= 6L) 3L else 4L
  n_row <- ceiling(n_panels / n_col)
  grDevices::png(path, width = 2400, height = max(1200, 700 * n_row), res = 180, bg = "#f7f3ea")
  graphics::par(mfrow = c(n_row, n_col), mar = c(0, 0, 3.2, 0))
  for (panel in panels) {
    plot_panel_aligned(spec, panel$coords, panel$title, panel$subtitle)
  }
  for (i in seq_len(n_row * n_col - n_panels)) {
    graphics::plot.new()
  }
  grDevices::dev.off()
}

make_set_c_sheet(carpet4, file.path(tmp_dir, "set-c-carpet4-contact-sheet.png"), include_igraph = TRUE)
make_set_c_sheet(weighted_mesh, file.path(tmp_dir, "set-c-weighted-mesh-contact-sheet.png"), include_igraph = TRUE)

grDevices::png(file.path(tmp_dir, "set-c-confirmation-contact-sheet.png"),
               width = 3200, height = 2200, res = 180, bg = "#f7f3ea")
graphics::par(mfrow = c(3, 5), mar = c(0, 0, 3.0, 0))
for (spec in list(carpet3, triangle4, mesh8)) {
  plot_panel_aligned(spec, spec$canonical, paste(spec$label, "| canonical"), "reference")
  for (key in c("baseline", "winner_a", "winner_b", "combined")) {
    row <- subset(set_c_metrics, graph == spec$id & label == set_c_configs[[key]]$label)
    plot_panel_aligned(
      spec,
      set_c_results[[spec$id]][[key]]$coords,
      paste(spec$label, "|", set_c_configs[[key]]$label),
      sprintf("RMSE %s | LGKK %s",
              fmt_num(row$rmse[[1L]]),
              fmt_num(row$lgkk_rel_rmse[[1L]]))
    )
  }
}
grDevices::dev.off()

message("Set A winner: ", set_a_winner)
message("Set B winner: ", set_b_winner)
message("Wrote outputs under ", tmp_dir)
