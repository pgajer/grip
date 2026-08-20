#!/usr/bin/env Rscript

output_root <- file.path("output", "gkk_lgkk_paper", "tmp", "carpet-level4-finalstage-repulsion-grid")

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'devtools' or the 'grip' package to run the carpet final-stage grid postprocess.")
}

helper_env <- new.env(parent = baseenv())
sys.source(file.path("tools", "benchmarks", "gkk_lgkk_paper", "benchmark-sierpinski-baseline.R"), envir = helper_env)

level <- 4L
repulsion_grid <- c(0.00, 0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 2.50)
final_rounds_grid <- 1:32
rmse_target_threshold <- 0.0230

built <- helper_env$build_sierpinski_carpet(level)
edges <- built$edges
canonical <- built$coords
n <- nrow(canonical)

summary_csv_path <- file.path(output_root, "carpet-level4-finalstage-repulsion-grid-summary.csv")
summary_df <- utils::read.csv(summary_csv_path, stringsAsFactors = FALSE)

format_num <- function(x, digits = 4L) {
  ifelse(is.finite(x), format(round(x, digits), nsmall = digits, trim = TRUE), "NA")
}

build_metric_matrix <- function(summary_df, metric) {
  mat <- matrix(NA_real_,
                nrow = length(final_rounds_grid),
                ncol = length(repulsion_grid),
                dimnames = list(final_rounds_grid, format(repulsion_grid, nsmall = 2)))
  for (i in seq_len(nrow(summary_df))) {
    row <- summary_df[i, , drop = FALSE]
    mat[as.character(row$final_rounds[[1L]]),
        format(row$repulsion_factor[[1L]], nsmall = 2)] <- row[[metric]][[1L]]
  }
  mat
}

draw_heatmaps <- function(path, summary_df) {
  metrics <- list(
    procrustes_rmse_mean = "Mean RMSE",
    edge_length_cv_mean = "Mean edge CV",
    edge_axis_deviation_mean = "Mean axis deviation",
    boundary_waviness_mean = "Mean boundary waviness"
  )

  grDevices::png(path, width = 1900, height = 1500, res = 170, bg = "#f7f3ea")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(4, 5, 3, 1), oma = c(0, 0, 2, 0))

  for (metric in names(metrics)) {
    mat <- build_metric_matrix(summary_df, metric)
    graphics::image(
      x = final_rounds_grid,
      y = repulsion_grid,
      z = mat,
      col = grDevices::hcl.colors(48, "YlOrRd", rev = TRUE),
      xlab = "final_rounds",
      ylab = "repulsion_factor",
      main = metrics[[metric]]
    )
    graphics::contour(
      x = final_rounds_grid,
      y = repulsion_grid,
      z = mat,
      add = TRUE,
      drawlabels = FALSE,
      col = grDevices::adjustcolor("#16324f", alpha.f = 0.35)
    )
  }

  graphics::mtext("Level-4 carpet FR endgame grid", side = 3, outer = TRUE, cex = 1.2, font = 2)
}

draw_frontier <- function(path, summary_df, kk_row, ref_f16, ref_f32) {
  grDevices::png(path, width = 1600, height = 1200, res = 170, bg = "#f7f3ea")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(5, 5, 4, 1))

  cols <- grDevices::colorRampPalette(c("#bfe3c0", "#0f6a43"))(length(repulsion_grid))
  col_map <- setNames(cols, format(repulsion_grid, nsmall = 2))
  pt_cols <- col_map[format(summary_df$repulsion_factor, nsmall = 2)]

  graphics::plot(summary_df$procrustes_rmse_mean,
                 summary_df$edge_axis_deviation_mean,
                 pch = 19,
                 col = pt_cols,
                 cex = 0.9,
                 xlab = "Mean Procrustes RMSE",
                 ylab = "Mean axis deviation",
                 main = "RMSE vs rectilinearity across the FR endgame grid")
  graphics::abline(v = kk_row$procrustes_rmse_mean[[1L]],
                   col = "#355070", lty = 2, lwd = 2)
  graphics::abline(v = rmse_target_threshold,
                   col = "#b23a48", lty = 3, lwd = 2)
  graphics::points(ref_f16$procrustes_rmse_mean,
                   ref_f16$edge_axis_deviation_mean,
                   pch = 17, cex = 1.4, col = "#7c3aed")
  graphics::points(ref_f32$procrustes_rmse_mean,
                   ref_f32$edge_axis_deviation_mean,
                   pch = 15, cex = 1.4, col = "#f97316")
  graphics::legend("topright",
                   legend = c("repulsion grid", "KK mean RMSE", "RMSE threshold", "default f16", "default f32"),
                   col = c("#0f6a43", "#355070", "#b23a48", "#7c3aed", "#f97316"),
                   pch = c(19, NA, NA, 17, 15),
                   lty = c(NA, 2, 3, NA, NA),
                   bty = "n")
}

draw_panel_sheet <- function(path, panels) {
  n_cols <- 3L
  n_rows <- ceiling(length(panels) / n_cols)
  grDevices::png(path, width = 900L * n_cols, height = 900L * n_rows, res = 170, bg = "#f7f3ea")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(n_rows, n_cols), mar = c(0, 0, 3.3, 0), xaxs = "i", yaxs = "i")
  for (panel in panels) {
    helper_env$plot_layout_panel(
      panel$coords,
      edges,
      title_text = panel$title,
      subtitle_text = panel$subtitle
    )
  }
}

ref_f1 <- summary_df[summary_df$final_rounds == 1L & abs(summary_df$repulsion_factor - 2.5) < 1e-9, , drop = FALSE]
ref_f16 <- summary_df[summary_df$final_rounds == 16L & abs(summary_df$repulsion_factor - 2.5) < 1e-9, , drop = FALSE]
ref_f32 <- summary_df[summary_df$final_rounds == 32L & abs(summary_df$repulsion_factor - 2.5) < 1e-9, , drop = FALSE]

best_rmse <- summary_df[which.min(summary_df$procrustes_rmse_mean), , drop = FALSE]
best_quality <- summary_df[which.min(summary_df$quality_rank_sum), , drop = FALSE]
near_kk <- summary_df[summary_df$procrustes_rmse_mean <= rmse_target_threshold, , drop = FALSE]
best_near_kk_f16 <- near_kk[which.min(near_kk$distance_to_f16), , drop = FALSE]
best_near_kk_f32 <- near_kk[which.min(near_kk$distance_to_f32), , drop = FALSE]

kk_available <- requireNamespace("igraph", quietly = TRUE)
kk_row <- data.frame(
  method = "igraph::layout_with_kk()",
  procrustes_rmse_mean = NA_real_,
  edge_length_cv_mean = NA_real_,
  edge_axis_deviation_mean = NA_real_,
  boundary_waviness_mean = NA_real_,
  elapsed_sec = NA_real_,
  stringsAsFactors = FALSE
)
kk_coords <- NULL
if (kk_available) {
  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  started <- proc.time()[["elapsed"]]
  set.seed(1L)
  kk_coords <- igraph::layout_with_kk(g)
  kk_elapsed <- proc.time()[["elapsed"]] - started
  kk_geom <- grip.geometry.diagnostics(
    coords = kk_coords,
    target.coords = canonical,
    edges = edges,
    family = "sierpinski.carpet",
    rng.seed = 1L
  )
  kk_edge_stats <- helper_env$edge_length_stats(kk_coords, edges)
  kk_row <- data.frame(
    method = "igraph::layout_with_kk()",
    procrustes_rmse_mean = kk_geom$procrustes.rmse[[1L]],
    edge_length_cv_mean = kk_edge_stats$cv,
    edge_axis_deviation_mean = kk_geom$edge.axis.deviation[[1L]],
    boundary_waviness_mean = kk_geom$boundary.waviness[[1L]],
    elapsed_sec = kk_elapsed,
    stringsAsFactors = FALSE
  )
}

true_no_final_trace <- grip.layout.trace(
  edges = edges,
  n = n,
  dim = 2,
  seed = 1L,
  final_rounds = 0L,
  repulsion_factor = 2.5,
  trace = "round",
  trace.every = 1L
)
true_no_final_idx <- which(true_no_final_trace$meta$phase == "level_start" &
                             true_no_final_trace$meta$misf_level == 0L)
true_no_final_coords <- true_no_final_trace$frames[[true_no_final_idx[[1L]]]]
true_no_final_fit <- helper_env$align_to_target(true_no_final_coords, canonical)

heatmap_path <- file.path(output_root, "carpet-level4-finalstage-repulsion-grid-heatmaps.png")
frontier_path <- file.path(output_root, "carpet-level4-finalstage-repulsion-grid-frontier.png")
panel_png_path <- file.path(output_root, "carpet-level4-finalstage-repulsion-grid-selected-panels.png")
summary_md_path <- file.path(output_root, "carpet-level4-finalstage-repulsion-grid-summary.md")
panel_table_csv_path <- file.path(output_root, "carpet-level4-finalstage-repulsion-grid-panel-table.csv")

draw_heatmaps(heatmap_path, summary_df)
draw_frontier(frontier_path, summary_df, kk_row, ref_f16, ref_f32)

panel_rows <- list(
  list(label = "true_no_final",
       title = "true_no_final",
       subtitle = sprintf("seed1 pre-final | RMSE=%s", format_num(true_no_final_fit$rmse, 4L)),
       coords = true_no_final_fit$aligned,
       rmse = true_no_final_fit$rmse,
       final_rounds = 0L,
       repulsion_factor = 2.5),
  list(label = "igraph_kk",
       title = "igraph::KK",
       subtitle = if (kk_available) sprintf("RMSE=%s", format_num(kk_row$procrustes_rmse_mean, 4L)) else "not available",
       coords = if (kk_available) helper_env$align_to_target(kk_coords, canonical)$aligned else canonical,
       rmse = if (kk_available) kk_row$procrustes_rmse_mean else NA_real_,
       final_rounds = NA_integer_,
       repulsion_factor = NA_real_),
  list(label = "f1_r2.5",
       title = "f1 r=2.5",
       subtitle = sprintf("sampled_f0 equiv | mean RMSE=%s", format_num(ref_f1$procrustes_rmse_mean, 4L)),
       coords = helper_env$align_to_target(
         grip.layout(edges = edges, n = n, dim = 2, seed = 1L, final_rounds = 1L, repulsion_factor = 2.5),
         canonical
       )$aligned,
       rmse = ref_f1$procrustes_rmse_mean,
       final_rounds = 1L,
       repulsion_factor = 2.5),
  list(label = "f16_r2.5",
       title = "f16 r=2.5",
       subtitle = sprintf("default-style ref | mean RMSE=%s", format_num(ref_f16$procrustes_rmse_mean, 4L)),
       coords = helper_env$align_to_target(
         grip.layout(edges = edges, n = n, dim = 2, seed = 1L, final_rounds = 16L, repulsion_factor = 2.5),
         canonical
       )$aligned,
       rmse = ref_f16$procrustes_rmse_mean,
       final_rounds = 16L,
       repulsion_factor = 2.5),
  list(label = "f32_r2.5",
       title = "f32 r=2.5",
       subtitle = sprintf("default-style ref | mean RMSE=%s", format_num(ref_f32$procrustes_rmse_mean, 4L)),
       coords = helper_env$align_to_target(
         grip.layout(edges = edges, n = n, dim = 2, seed = 1L, final_rounds = 32L, repulsion_factor = 2.5),
         canonical
       )$aligned,
       rmse = ref_f32$procrustes_rmse_mean,
       final_rounds = 32L,
       repulsion_factor = 2.5),
  list(label = "best_rmse",
       title = sprintf("best RMSE f%d r=%.2f", best_rmse$final_rounds[[1L]], best_rmse$repulsion_factor[[1L]]),
       subtitle = sprintf("mean RMSE=%s", format_num(best_rmse$procrustes_rmse_mean, 4L)),
       coords = helper_env$align_to_target(
         grip.layout(edges = edges, n = n, dim = 2, seed = 1L,
                     final_rounds = best_rmse$final_rounds[[1L]],
                     repulsion_factor = best_rmse$repulsion_factor[[1L]]),
         canonical
       )$aligned,
       rmse = best_rmse$procrustes_rmse_mean,
       final_rounds = best_rmse$final_rounds[[1L]],
       repulsion_factor = best_rmse$repulsion_factor[[1L]]),
  list(label = "best_near_kk_f16",
       title = sprintf("near-KK close-to-f16 f%d r=%.2f", best_near_kk_f16$final_rounds[[1L]], best_near_kk_f16$repulsion_factor[[1L]]),
       subtitle = sprintf("mean RMSE=%s", format_num(best_near_kk_f16$procrustes_rmse_mean, 4L)),
       coords = helper_env$align_to_target(
         grip.layout(edges = edges, n = n, dim = 2, seed = 1L,
                     final_rounds = best_near_kk_f16$final_rounds[[1L]],
                     repulsion_factor = best_near_kk_f16$repulsion_factor[[1L]]),
         canonical
       )$aligned,
       rmse = best_near_kk_f16$procrustes_rmse_mean,
       final_rounds = best_near_kk_f16$final_rounds[[1L]],
       repulsion_factor = best_near_kk_f16$repulsion_factor[[1L]]),
  list(label = "best_near_kk_f32",
       title = sprintf("near-KK close-to-f32 f%d r=%.2f", best_near_kk_f32$final_rounds[[1L]], best_near_kk_f32$repulsion_factor[[1L]]),
       subtitle = sprintf("mean RMSE=%s", format_num(best_near_kk_f32$procrustes_rmse_mean, 4L)),
       coords = helper_env$align_to_target(
         grip.layout(edges = edges, n = n, dim = 2, seed = 1L,
                     final_rounds = best_near_kk_f32$final_rounds[[1L]],
                     repulsion_factor = best_near_kk_f32$repulsion_factor[[1L]]),
         canonical
       )$aligned,
       rmse = best_near_kk_f32$procrustes_rmse_mean,
       final_rounds = best_near_kk_f32$final_rounds[[1L]],
       repulsion_factor = best_near_kk_f32$repulsion_factor[[1L]]),
  list(label = "best_quality",
       title = sprintf("best quality f%d r=%.2f", best_quality$final_rounds[[1L]], best_quality$repulsion_factor[[1L]]),
       subtitle = sprintf("rank sum=%s", format_num(best_quality$quality_rank_sum, 1L)),
       coords = helper_env$align_to_target(
         grip.layout(edges = edges, n = n, dim = 2, seed = 1L,
                     final_rounds = best_quality$final_rounds[[1L]],
                     repulsion_factor = best_quality$repulsion_factor[[1L]]),
         canonical
       )$aligned,
       rmse = best_quality$procrustes_rmse_mean,
       final_rounds = best_quality$final_rounds[[1L]],
       repulsion_factor = best_quality$repulsion_factor[[1L]])
)

draw_panel_sheet(panel_png_path, lapply(panel_rows, function(x) x[c("coords", "title", "subtitle")]))

panel_table <- do.call(
  rbind,
  lapply(panel_rows, function(x) {
    data.frame(
      panel = x$label,
      final_rounds = x$final_rounds,
      repulsion_factor = x$repulsion_factor,
      rmse = x$rmse,
      stringsAsFactors = FALSE
    )
  })
)
utils::write.csv(panel_table, panel_table_csv_path, row.names = FALSE)

lines <- c(
  "# Carpet Level-4 Final-Stage Repulsion Grid",
  "",
  sprintf("- graph: level-%d Sierpinski carpet (`%d` vertices, `%d` edges)", level, n, nrow(edges)),
  sprintf("- final_rounds grid: `%d..%d`", min(final_rounds_grid), max(final_rounds_grid)),
  sprintf("- repulsion_factor grid: `%s`", paste(format(repulsion_grid, nsmall = 2), collapse = ", ")),
  sprintf("- near-KK RMSE threshold: `%.4f`", rmse_target_threshold),
  "",
  "## Headline results",
  "",
  sprintf("- igraph::KK mean RMSE: `%s`", format_num(kk_row$procrustes_rmse_mean, 4L)),
  sprintf("- best RMSE in FR grid: `%s` at `final_rounds=%d`, `repulsion_factor=%.2f`",
          format_num(best_rmse$procrustes_rmse_mean, 4L),
          best_rmse$final_rounds[[1L]],
          best_rmse$repulsion_factor[[1L]]),
  sprintf("- sampled_f0 equivalent (`f1 r=2.5`) mean RMSE: `%s`", format_num(ref_f1$procrustes_rmse_mean, 4L)),
  sprintf("- default-style `f16 r=2.5` mean RMSE: `%s`", format_num(ref_f16$procrustes_rmse_mean, 4L)),
  sprintf("- default-style `f32 r=2.5` mean RMSE: `%s`", format_num(ref_f32$procrustes_rmse_mean, 4L)),
  sprintf("- near-KK config closest to f16 smoothness: `f%d r=%.2f` with RMSE `%s`",
          best_near_kk_f16$final_rounds[[1L]],
          best_near_kk_f16$repulsion_factor[[1L]],
          format_num(best_near_kk_f16$procrustes_rmse_mean, 4L)),
  sprintf("- near-KK config closest to f32 smoothness: `f%d r=%.2f` with RMSE `%s`",
          best_near_kk_f32$final_rounds[[1L]],
          best_near_kk_f32$repulsion_factor[[1L]],
          format_num(best_near_kk_f32$procrustes_rmse_mean, 4L)),
  "",
  "## Outputs",
  "",
  sprintf("- summary CSV: `%s`", summary_csv_path),
  sprintf("- heatmaps: `%s`", heatmap_path),
  sprintf("- frontier plot: `%s`", frontier_path),
  sprintf("- selected panels: `%s`", panel_png_path),
  sprintf("- panel RMSE table: `%s`", panel_table_csv_path)
)
writeLines(lines, con = summary_md_path)

message(sprintf("Heatmaps written to %s", heatmap_path))
message(sprintf("Frontier plot written to %s", frontier_path))
message(sprintf("Selected-panel sheet written to %s", panel_png_path))
message(sprintf("Panel table written to %s", panel_table_csv_path))
message(sprintf("Markdown summary written to %s", summary_md_path))
