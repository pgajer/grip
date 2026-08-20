#!/usr/bin/env Rscript

output_root <- file.path("output", "gkk_lgkk_paper", "tmp", "carpet-level4-finalstage-structural-round3")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'devtools' or the 'grip' package to run the round-3 carpet structural study.")
}

helper_env <- new.env(parent = baseenv())
sys.source(file.path("tools", "benchmarks", "gkk_lgkk_paper", "benchmark-sierpinski-baseline.R"), envir = helper_env)

level <- 4L
seeds <- 1:6
final_rounds_grid <- c(8L, 16L, 32L)
anchor_grid <- c(0.00, 0.25, 0.50, 1.00, 2.00, 4.00)
move_scale_grid <- c(1.00, 0.75, 0.50, 0.25, 0.10)
repulsion_factor <- 2.5

built <- helper_env$build_sierpinski_carpet(level)
edges <- built$edges
canonical <- built$coords
n <- nrow(canonical)

format_num <- function(x, digits = 4L) {
  ifelse(is.finite(x), format(round(x, digits), nsmall = digits, trim = TRUE), "NA")
}

time_it <- function(expr) {
  started <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, elapsed = proc.time()[["elapsed"]] - started)
}

config_df <- expand.grid(
  final_rounds = final_rounds_grid,
  final_anchor_factor = anchor_grid,
  final_move_scale_after_first = move_scale_grid,
  stringsAsFactors = FALSE
)
config_df$config_id <- sprintf(
  "f%02d_a%s_m%s",
  config_df$final_rounds,
  gsub("\\.", "p", format(config_df$final_anchor_factor, nsmall = 2)),
  gsub("\\.", "p", format(config_df$final_move_scale_after_first, nsmall = 2))
)
config_df$is_structural <- with(
  config_df,
  final_anchor_factor > 0 | abs(final_move_scale_after_first - 1) > 1e-12
)
config_df <- config_df[, c(
  "config_id", "final_rounds", "final_anchor_factor",
  "final_move_scale_after_first", "is_structural"
)]

run_one <- function(config_row, seed) {
  timed <- time_it({
    grip.layout(
      edges = edges,
      n = n,
      dim = 2,
      seed = seed,
      final_rounds = config_row$final_rounds[[1L]],
      repulsion_factor = repulsion_factor,
      final_anchor_factor = config_row$final_anchor_factor[[1L]],
      final_move_scale_after_first = config_row$final_move_scale_after_first[[1L]]
    )
  })

  geom <- grip.geometry.diagnostics(
    coords = timed$value,
    target.coords = canonical,
    edges = edges,
    family = "sierpinski.carpet",
    rng.seed = 1L
  )
  edge_stats <- helper_env$edge_length_stats(timed$value, edges)

  data.frame(
    config_id = config_row$config_id[[1L]],
    final_rounds = config_row$final_rounds[[1L]],
    final_anchor_factor = config_row$final_anchor_factor[[1L]],
    final_move_scale_after_first = config_row$final_move_scale_after_first[[1L]],
    is_structural = config_row$is_structural[[1L]],
    seed = seed,
    procrustes_rmse = geom$procrustes.rmse[[1L]],
    global_symmetry_score = geom$global.symmetry.score[[1L]],
    local_angle_deviation = geom$local.angle.deviation[[1L]],
    edge_axis_deviation = geom$edge.axis.deviation[[1L]],
    boundary_waviness = geom$boundary.waviness[[1L]],
    corridor_waviness = geom$corridor.waviness[[1L]],
    hole_center_error = geom$hole.center.error[[1L]],
    edge_length_cv = edge_stats$cv,
    elapsed_sec = timed$elapsed,
    stringsAsFactors = FALSE
  )
}

summarize_configs <- function(raw_metrics) {
  summary_df <- do.call(
    rbind,
    lapply(split(raw_metrics, raw_metrics$config_id), function(df) {
      data.frame(
        config_id = df$config_id[[1L]],
        final_rounds = df$final_rounds[[1L]],
        final_anchor_factor = df$final_anchor_factor[[1L]],
        final_move_scale_after_first = df$final_move_scale_after_first[[1L]],
        is_structural = df$is_structural[[1L]],
        seeds = nrow(df),
        procrustes_rmse_mean = mean(df$procrustes_rmse),
        procrustes_rmse_sd = stats::sd(df$procrustes_rmse),
        edge_length_cv_mean = mean(df$edge_length_cv),
        edge_axis_deviation_mean = mean(df$edge_axis_deviation),
        boundary_waviness_mean = mean(df$boundary_waviness),
        corridor_waviness_mean = mean(df$corridor_waviness),
        hole_center_error_mean = mean(df$hole_center_error),
        elapsed_sec_mean = mean(df$elapsed_sec),
        stringsAsFactors = FALSE
      )
    })
  )

  rank_asc <- function(x) rank(x, ties.method = "average")
  summary_df$quality_rank_sum <- with(
    summary_df,
    rank_asc(procrustes_rmse_mean) +
      rank_asc(edge_length_cv_mean) +
      rank_asc(edge_axis_deviation_mean) +
      rank_asc(boundary_waviness_mean) +
      rank_asc(corridor_waviness_mean) +
      rank_asc(hole_center_error_mean)
  )

  summary_df[order(summary_df$quality_rank_sum, summary_df$procrustes_rmse_mean), , drop = FALSE]
}

build_metric_matrix <- function(summary_df, final_rounds_value, metric) {
  sub <- summary_df[summary_df$final_rounds == final_rounds_value, , drop = FALSE]
  mat <- matrix(
    NA_real_,
    nrow = length(move_scale_grid),
    ncol = length(anchor_grid),
    dimnames = list(format(move_scale_grid, nsmall = 2), format(anchor_grid, nsmall = 2))
  )
  for (i in seq_len(nrow(sub))) {
    row <- sub[i, , drop = FALSE]
    mat[format(row$final_move_scale_after_first[[1L]], nsmall = 2),
        format(row$final_anchor_factor[[1L]], nsmall = 2)] <- row[[metric]][[1L]]
  }
  mat
}

draw_heatmaps <- function(path, summary_df) {
  metrics <- list(
    procrustes_rmse_mean = "Mean RMSE",
    edge_length_cv_mean = "Mean edge CV",
    edge_axis_deviation_mean = "Mean axis deviation",
    hole_center_error_mean = "Mean hole-center error"
  )

  grDevices::png(path, width = 2200, height = 1800, res = 170, bg = "#f7f3ea")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(length(final_rounds_grid), length(metrics)),
                mar = c(4, 4, 2.7, 1), oma = c(0, 0, 3, 0))

  for (fr in final_rounds_grid) {
    for (metric in names(metrics)) {
      mat <- build_metric_matrix(summary_df, fr, metric)
      y_vals <- sort(move_scale_grid)
      mat_ord <- mat[format(y_vals, nsmall = 2), , drop = FALSE]
      graphics::image(
        x = anchor_grid,
        y = y_vals,
        z = t(mat_ord),
        col = grDevices::hcl.colors(48, "YlOrRd", rev = TRUE),
        xlab = "final_anchor_factor",
        ylab = "final_move_scale_after_first",
        main = sprintf("f%d | %s", fr, metrics[[metric]])
      )
      graphics::contour(
        x = anchor_grid,
        y = y_vals,
        z = t(mat_ord),
        add = TRUE,
        drawlabels = FALSE,
        col = grDevices::adjustcolor("#16324f", alpha.f = 0.35)
      )
    }
  }

  graphics::mtext("Round-3 structural grid on the level-4 carpet",
                  side = 3, outer = TRUE, cex = 1.2, font = 2)
}

draw_frontier <- function(path, summary_df, kk_row, reference_rows) {
  cols <- c("8" = "#355070", "16" = "#6d597a", "32" = "#b56576")
  pch_map <- c("1.00" = 1, "0.75" = 17, "0.50" = 15, "0.25" = 18, "0.10" = 19)

  grDevices::png(path, width = 1700, height = 1300, res = 170, bg = "#f7f3ea")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(5, 5, 4, 1))

  graphics::plot(summary_df$procrustes_rmse_mean,
                 summary_df$edge_axis_deviation_mean,
                 type = "n",
                 xlab = "Mean Procrustes RMSE",
                 ylab = "Mean axis deviation",
                 main = "Round-3 RMSE vs rectilinearity frontier")

  for (i in seq_len(nrow(summary_df))) {
    row <- summary_df[i, , drop = FALSE]
    graphics::points(
      row$procrustes_rmse_mean,
      row$edge_axis_deviation_mean,
      pch = pch_map[[format(row$final_move_scale_after_first[[1L]], nsmall = 2)]],
      col = cols[[as.character(row$final_rounds[[1L]])]],
      cex = 0.9
    )
  }

  graphics::points(kk_row$procrustes_rmse_mean,
                   kk_row$edge_axis_deviation_mean,
                   pch = 8, cex = 1.6, col = "#0a9396")

  for (nm in names(reference_rows)) {
    row <- reference_rows[[nm]]
    graphics::points(row$procrustes_rmse_mean,
                     row$edge_axis_deviation_mean,
                     pch = 4, cex = 1.4, col = "#ca6702", lwd = 1.5)
  }

  best_labels <- summary_df[order(summary_df$quality_rank_sum, summary_df$procrustes_rmse_mean), , drop = FALSE]
  best_labels <- best_labels[seq_len(min(6L, nrow(best_labels))), , drop = FALSE]
  graphics::text(best_labels$procrustes_rmse_mean,
                 best_labels$edge_axis_deviation_mean,
                 labels = best_labels$config_id,
                 pos = 4, cex = 0.75, offset = 0.5, col = "#16324f")

  graphics::legend(
    "topright",
    legend = c("f8", "f16", "f32", "move 1.00", "move 0.50", "move 0.10", "igraph::KK", "round-2 refs"),
    col = c(cols[["8"]], cols[["16"]], cols[["32"]], "black", "black", "black", "#0a9396", "#ca6702"),
    pch = c(19, 19, 19, pch_map[["1.00"]], pch_map[["0.50"]], pch_map[["0.10"]], 8, 4),
    bty = "n"
  )
}

draw_panel_sheet <- function(path, panels, n_cols = 4L) {
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

message("Running round-3 structural grid on level-4 carpet...")
raw_metrics <- do.call(
  rbind,
  lapply(seq_len(nrow(config_df)), function(i) {
    config_row <- config_df[i, , drop = FALSE]
    message(sprintf(
      "  %s (%d/%d)",
      config_row$config_id[[1L]],
      i,
      nrow(config_df)
    ))
    do.call(rbind, lapply(seeds, function(seed) run_one(config_row, seed)))
  })
)

message("Finished raw_metrics.")
raw_path <- file.path(output_root, "carpet-level4-finalstage-structural-round3-raw.csv")
utils::write.csv(raw_metrics, raw_path, row.names = FALSE)

summary_df <- summarize_configs(raw_metrics)
message("Finished summary_df.")

summary_path <- file.path(output_root, "carpet-level4-finalstage-structural-round3-summary.csv")
utils::write.csv(summary_df, summary_path, row.names = FALSE)
message("Wrote raw and summary CSV files.")

kk_row <- data.frame(
  method = "igraph::layout_with_kk()",
  procrustes_rmse_mean = NA_real_,
  edge_length_cv_mean = NA_real_,
  edge_axis_deviation_mean = NA_real_,
  boundary_waviness_mean = NA_real_,
  corridor_waviness_mean = NA_real_,
  hole_center_error_mean = NA_real_,
  elapsed_sec_mean = NA_real_,
  stringsAsFactors = FALSE
)

panel_seed <- 1L
if (requireNamespace("igraph", quietly = TRUE)) {
  message("Starting igraph KK reference.")
  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  timed_kk <- time_it(igraph::layout_with_kk(g))
  kk_coords <- as.matrix(timed_kk$value)
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
    corridor_waviness_mean = kk_geom$corridor.waviness[[1L]],
    hole_center_error_mean = kk_geom$hole.center.error[[1L]],
    elapsed_sec_mean = timed_kk$elapsed,
    stringsAsFactors = FALSE
  )
} else {
  kk_coords <- NULL
}
message("Finished igraph KK reference.")

reference_ids <- c(
  f8 = "f08_a0p00_m1p00",
  f16 = "f16_a0p00_m1p00",
  f32 = "f32_a0p00_m1p00"
)
reference_rows <- lapply(reference_ids, function(id) {
  summary_df[summary_df$config_id == id, , drop = FALSE]
})

f1_trace <- grip.layout.trace(
  edges = edges,
  n = n,
  dim = 2,
  seed = panel_seed,
  final_rounds = 0L,
  trace = "round",
  trace.every = 1L
)
true_no_final_idx <- which(
  f1_trace$meta$phase == "level_start" &
    f1_trace$meta$misf_level == 0L
)
if (length(true_no_final_idx) == 0L) {
  stop("Could not locate the pre-final full-graph snapshot.")
}
true_no_final_coords <- f1_trace$frames[[true_no_final_idx[[1L]]]]
true_no_final_fit <- helper_env$align_to_target(true_no_final_coords, canonical)
true_no_final_geom <- grip.geometry.diagnostics(
  coords = true_no_final_coords,
  target.coords = canonical,
  edges = edges,
  family = "sierpinski.carpet",
  rng.seed = 1L
)
message("Finished trace-derived references.")

best_overall_struct <- summary_df[summary_df$is_structural, , drop = FALSE]
best_overall_struct <- best_overall_struct[which.min(best_overall_struct$quality_rank_sum), , drop = FALSE]

best_by_budget <- lapply(final_rounds_grid, function(fr) {
  sub <- summary_df[summary_df$final_rounds == fr & summary_df$is_structural, , drop = FALSE]
  sub[which.min(sub$quality_rank_sum), , drop = FALSE]
})
names(best_by_budget) <- paste0("f", final_rounds_grid)

best_near_f1_struct <- summary_df[summary_df$is_structural, , drop = FALSE]
best_near_f1_struct <- best_near_f1_struct[order(best_near_f1_struct$procrustes_rmse_mean,
                                                 best_near_f1_struct$edge_axis_deviation_mean), , drop = FALSE]
best_near_f1_struct <- best_near_f1_struct[1L, , drop = FALSE]

panel_from_config <- function(config_row, seed = panel_seed, title = NULL) {
  coords <- grip.layout(
    edges = edges,
    n = n,
    dim = 2,
    seed = seed,
    final_rounds = config_row$final_rounds[[1L]],
    repulsion_factor = repulsion_factor,
    final_anchor_factor = config_row$final_anchor_factor[[1L]],
    final_move_scale_after_first = config_row$final_move_scale_after_first[[1L]]
  )
  fit <- helper_env$align_to_target(coords, canonical)
  list(
    coords = fit$aligned,
    title = if (is.null(title)) config_row$config_id[[1L]] else title,
    subtitle = sprintf(
      "RMSE=%s | edge CV=%s | axis=%s",
      format_num(config_row$procrustes_rmse_mean[[1L]], 4L),
      format_num(config_row$edge_length_cv_mean[[1L]], 4L),
      format_num(config_row$edge_axis_deviation_mean[[1L]], 4L)
    )
  )
}

panel_specs <- list(
  list(
    coords = canonical,
    title = "canonical",
    subtitle = "target embedding"
  ),
  list(
    coords = true_no_final_fit$aligned,
    title = "true_no_final",
    subtitle = sprintf("pre-final snapshot | RMSE=%s",
                       format_num(true_no_final_geom$procrustes.rmse[[1L]], 4L))
  )
)

if (!is.null(kk_coords)) {
  kk_fit <- helper_env$align_to_target(kk_coords, canonical)
  panel_specs[[length(panel_specs) + 1L]] <- list(
    coords = kk_fit$aligned,
    title = "igraph::KK",
    subtitle = sprintf("RMSE=%s | edge CV=%s | axis=%s",
                       format_num(kk_row$procrustes_rmse_mean[[1L]], 4L),
                       format_num(kk_row$edge_length_cv_mean[[1L]], 4L),
                       format_num(kk_row$edge_axis_deviation_mean[[1L]], 4L))
  )
}

f1_coords <- f1_trace$final
f1_geom <- grip.geometry.diagnostics(
  coords = f1_coords,
  target.coords = canonical,
  edges = edges,
  family = "sierpinski.carpet",
  rng.seed = 1L
)
f1_edge_stats <- helper_env$edge_length_stats(f1_coords, edges)
panel_specs[[length(panel_specs) + 1L]] <- list(
  coords = helper_env$align_to_target(f1_coords, canonical)$aligned,
  title = "f1 baseline",
  subtitle = sprintf("RMSE=%s | edge CV=%s | axis=%s",
                     format_num(f1_geom$procrustes.rmse[[1L]], 4L),
                     format_num(f1_edge_stats$cv, 4L),
                     format_num(f1_geom$edge.axis.deviation[[1L]], 4L))
)

panel_specs[[length(panel_specs) + 1L]] <- panel_from_config(reference_rows$f8, title = "f8 baseline")
panel_specs[[length(panel_specs) + 1L]] <- panel_from_config(best_by_budget$f8, title = "best f8 structural")
panel_specs[[length(panel_specs) + 1L]] <- panel_from_config(reference_rows$f16, title = "f16 baseline")
panel_specs[[length(panel_specs) + 1L]] <- panel_from_config(best_by_budget$f16, title = "best f16 structural")
panel_specs[[length(panel_specs) + 1L]] <- panel_from_config(reference_rows$f32, title = "f32 baseline")
panel_specs[[length(panel_specs) + 1L]] <- panel_from_config(best_by_budget$f32, title = "best f32 structural")
panel_specs[[length(panel_specs) + 1L]] <- panel_from_config(best_overall_struct, title = "best overall structural")
panel_specs[[length(panel_specs) + 1L]] <- panel_from_config(best_near_f1_struct, title = "closest structural to f1 RMSE")

panel_table <- do.call(
  rbind,
  lapply(panel_specs, function(panel) {
    coords <- panel$coords
    geom <- grip.geometry.diagnostics(
      coords = coords,
      target.coords = canonical,
      edges = edges,
      family = "sierpinski.carpet",
      rng.seed = 1L
    )
    edge_stats <- helper_env$edge_length_stats(coords, edges)
    data.frame(
      panel = panel$title,
      procrustes_rmse = geom$procrustes.rmse[[1L]],
      edge_length_cv = edge_stats$cv,
      edge_axis_deviation = geom$edge.axis.deviation[[1L]],
      stringsAsFactors = FALSE
    )
  })
)
message("Built panel table.")

heatmap_path <- file.path(output_root, "carpet-level4-finalstage-structural-round3-heatmaps.png")
frontier_path <- file.path(output_root, "carpet-level4-finalstage-structural-round3-frontier.png")
panels_path <- file.path(output_root, "carpet-level4-finalstage-structural-round3-selected-panels.png")
panel_table_path <- file.path(output_root, "carpet-level4-finalstage-structural-round3-panel-table.csv")
summary_md_path <- file.path(output_root, "carpet-level4-finalstage-structural-round3-summary.md")

draw_heatmaps(heatmap_path, summary_df)
draw_frontier(frontier_path, summary_df, kk_row, reference_rows)
draw_panel_sheet(panels_path, panel_specs, n_cols = 4L)
message("Rendered figures.")

utils::write.csv(panel_table, panel_table_path, row.names = FALSE)

lines <- c(
  "# Round-3 structural final-stage study",
  "",
  sprintf("- graph: level-%d Sierpinski carpet", level),
  sprintf("- seeds: %s", paste(seeds, collapse = ", ")),
  sprintf("- grid size: %d configurations, %d seeded runs", nrow(config_df), nrow(raw_metrics)),
  sprintf("- fixed repulsion_factor: %s", format_num(repulsion_factor, 2L)),
  "",
  "## Headline findings",
  "",
  sprintf("- best overall structural config: `%s` (RMSE `%s`, edge CV `%s`, axis `%s`, sec `%s`)",
          best_overall_struct$config_id,
          format_num(best_overall_struct$procrustes_rmse_mean),
          format_num(best_overall_struct$edge_length_cv_mean),
          format_num(best_overall_struct$edge_axis_deviation_mean),
          format_num(best_overall_struct$elapsed_sec_mean, 3L)),
  sprintf("- closest structural config to f1 RMSE: `%s` (RMSE `%s`, edge CV `%s`, axis `%s`)",
          best_near_f1_struct$config_id,
          format_num(best_near_f1_struct$procrustes_rmse_mean),
          format_num(best_near_f1_struct$edge_length_cv_mean),
          format_num(best_near_f1_struct$edge_axis_deviation_mean)),
  sprintf("- best f8 structural config: `%s` (RMSE `%s`, edge CV `%s`, axis `%s`)",
          best_by_budget$f8$config_id,
          format_num(best_by_budget$f8$procrustes_rmse_mean),
          format_num(best_by_budget$f8$edge_length_cv_mean),
          format_num(best_by_budget$f8$edge_axis_deviation_mean)),
  sprintf("- best f16 structural config: `%s` (RMSE `%s`, edge CV `%s`, axis `%s`)",
          best_by_budget$f16$config_id,
          format_num(best_by_budget$f16$procrustes_rmse_mean),
          format_num(best_by_budget$f16$edge_length_cv_mean),
          format_num(best_by_budget$f16$edge_axis_deviation_mean)),
  sprintf("- best f32 structural config: `%s` (RMSE `%s`, edge CV `%s`, axis `%s`)",
          best_by_budget$f32$config_id,
          format_num(best_by_budget$f32$procrustes_rmse_mean),
          format_num(best_by_budget$f32$edge_length_cv_mean),
          format_num(best_by_budget$f32$edge_axis_deviation_mean)),
  "",
  "## Files",
  "",
  sprintf("- raw metrics: `%s`", raw_path),
  sprintf("- summary: `%s`", summary_path),
  sprintf("- heatmaps: `%s`", heatmap_path),
  sprintf("- frontier: `%s`", frontier_path),
  sprintf("- selected panels: `%s`", panels_path),
  sprintf("- panel table: `%s`", panel_table_path)
)
writeLines(lines, summary_md_path)

message("Round-3 structural study complete.")
