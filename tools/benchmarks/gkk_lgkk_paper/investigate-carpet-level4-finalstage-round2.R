#!/usr/bin/env Rscript

output_root <- file.path("output", "gkk_lgkk_paper", "tmp", "carpet-level4-finalstage-round2")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'devtools' or the 'grip' package to run the carpet final-stage study.")
}

helper_env <- new.env(parent = baseenv())
sys.source(file.path("tools", "benchmarks", "gkk_lgkk_paper", "benchmark-sierpinski-baseline.R"), envir = helper_env)

level <- 4L
seeds <- 1:6
built <- helper_env$build_sierpinski_carpet(level)
edges <- built$edges
canonical <- built$coords
n <- nrow(canonical)

sampled_finals <- c(0L, 32L, 64L, 96L, 128L, 160L, 192L, 256L, 384L)
exact_finals <- c(0L, 64L, 96L, 128L, 160L, 192L)

config_df <- rbind(
  data.frame(
    config_id = paste0("sampled_f", sampled_finals),
    coarse_mode = "sampled",
    final_rounds = sampled_finals,
    coarse_repulsion_exact_below = 64L,
    stringsAsFactors = FALSE
  ),
  data.frame(
    config_id = paste0("exact_f", exact_finals),
    coarse_mode = "exact",
    final_rounds = exact_finals,
    coarse_repulsion_exact_below = n,
    stringsAsFactors = FALSE
  )
)

format_num <- function(x, digits = 4L) {
  ifelse(is.finite(x), format(round(x, digits), nsmall = digits, trim = TRUE), "NA")
}

edge_axis_deviation <- function(coords_aligned, edges) {
  diffs <- coords_aligned[edges[, 1L], , drop = FALSE] -
    coords_aligned[edges[, 2L], , drop = FALSE]
  lens <- sqrt(rowSums(diffs^2))
  keep <- is.finite(lens) & lens > 0
  if (!any(keep)) {
    return(NA_real_)
  }
  mean(pmin(abs(diffs[keep, 1L]), abs(diffs[keep, 2L])) / lens[keep])
}

time_it <- function(expr) {
  started <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, elapsed = proc.time()[["elapsed"]] - started)
}

run_one <- function(config_row, seed) {
  timed <- time_it({
    grip.layout(
      edges = edges,
      n = n,
      dim = 2,
      seed = seed,
      final_rounds = config_row$final_rounds[[1L]],
      coarse_repulsion_exact_below = config_row$coarse_repulsion_exact_below[[1L]]
    )
  })

  fit <- helper_env$align_to_target(timed$value, canonical)
  quality <- grip.score.layout(
    coords = timed$value,
    edges = edges,
    n = n,
    sample.size.stress = 2000L,
    sample.size.nonedge = 5000L,
    edge.crossings = "never"
  )

  data.frame(
    config_id = config_row$config_id[[1L]],
    coarse_mode = config_row$coarse_mode[[1L]],
    final_rounds = config_row$final_rounds[[1L]],
    coarse_repulsion_exact_below = config_row$coarse_repulsion_exact_below[[1L]],
    seed = seed,
    procrustes_rmse = fit$rmse,
    edge_length_cv = quality$edge.length.cv,
    sampled_nonedge_sep_ratio = quality$sampled.nonedge.sep.ratio,
    sampled_stress = quality$sampled.stress,
    axis_deviation = edge_axis_deviation(fit$aligned, edges),
    elapsed_sec = timed$elapsed,
    stringsAsFactors = FALSE
  )
}

draw_final_rounds_lines <- function(path, summary_df) {
  sampled <- summary_df[summary_df$coarse_mode == "sampled", , drop = FALSE]
  sampled <- sampled[order(sampled$final_rounds), , drop = FALSE]

  grDevices::png(path, width = 1800, height = 1200, res = 170, bg = "#f7f3ea")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

  plot_metric <- function(y, ylab, main, col) {
    graphics::plot(sampled$final_rounds, y,
                   type = "b",
                   pch = 19,
                   lwd = 2,
                   col = col,
                   xlab = "final_rounds",
                   ylab = ylab,
                   main = main)
  }

  plot_metric(sampled$procrustes_rmse_mean, "Mean RMSE", "Global symmetry vs final_rounds", "#b23a48")
  plot_metric(sampled$edge_length_cv_mean, "Mean edge CV", "Edge regularity vs final_rounds", "#355070")
  plot_metric(sampled$axis_deviation_mean, "Mean axis deviation", "Rectilinearity proxy vs final_rounds", "#3d7a57")
  plot_metric(sampled$elapsed_sec_mean, "Mean sec", "Runtime vs final_rounds", "#7c4d8b")

  graphics::mtext("Level-4 carpet, sampled coarse repulsion", side = 3, outer = TRUE, cex = 1.2, font = 2)
}

draw_tradeoff_scatter <- function(path, summary_df) {
  xx <- summary_df$procrustes_rmse_mean
  yy <- summary_df$axis_deviation_mean
  cols <- ifelse(summary_df$coarse_mode == "exact", "#355070", "#b23a48")
  pchs <- ifelse(summary_df$coarse_mode == "exact", 17, 19)

  grDevices::png(path, width = 1500, height = 1100, res = 170, bg = "#f7f3ea")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(5, 5, 4, 1))
  graphics::plot(xx, yy,
                 pch = pchs,
                 col = cols,
                 cex = 1.3,
                 xlab = "Mean Procrustes RMSE to canonical carpet",
                 ylab = "Mean axis deviation (lower is crisper)",
                 main = "Round 2 tradeoff map: symmetry vs rectilinearity")
  graphics::text(xx, yy,
                 labels = summary_df$short_label,
                 pos = 3,
                 cex = 0.85,
                 col = "#16324f")
  graphics::legend("topright",
                   legend = c("sampled coarse repulsion", "exact coarse repulsion"),
                   col = c("#b23a48", "#355070"),
                   pch = c(19, 17),
                   bty = "n")
}

draw_contact_sheet <- function(path, panels) {
  n_panels <- length(panels)
  n_cols <- ceiling(sqrt(n_panels))
  n_rows <- ceiling(n_panels / n_cols)

  grDevices::png(path,
                 width = 800L * n_cols,
                 height = 800L * n_rows,
                 res = 170,
                 bg = "#f7f3ea")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(n_rows, n_cols), mar = c(0, 0, 3, 0), xaxs = "i", yaxs = "i")

  for (panel in panels) {
    helper_env$plot_layout_panel(
      panel$coords,
      edges,
      title_text = panel$title,
      subtitle_text = panel$subtitle
    )
  }
}

message("Running level-4 carpet round-2 sweep...")
raw_metrics <- do.call(
  rbind,
  lapply(seq_len(nrow(config_df)), function(i) {
    config_row <- config_df[i, , drop = FALSE]
    do.call(
      rbind,
      lapply(seeds, function(seed) {
        message(sprintf(
          "  %s (%s, final_rounds=%d), seed %d/%d",
          config_row$config_id[[1L]],
          config_row$coarse_mode[[1L]],
          config_row$final_rounds[[1L]],
          seed,
          length(seeds)
        ))
        run_one(config_row, seed)
      })
    )
  })
)

summary_df <- do.call(
  rbind,
  lapply(split(raw_metrics, raw_metrics$config_id), function(df) {
    data.frame(
      config_id = df$config_id[[1L]],
      coarse_mode = df$coarse_mode[[1L]],
      final_rounds = df$final_rounds[[1L]],
      coarse_repulsion_exact_below = df$coarse_repulsion_exact_below[[1L]],
      seeds = nrow(df),
      procrustes_rmse_mean = mean(df$procrustes_rmse),
      procrustes_rmse_sd = stats::sd(df$procrustes_rmse),
      edge_length_cv_mean = mean(df$edge_length_cv),
      edge_length_cv_sd = stats::sd(df$edge_length_cv),
      sampled_nonedge_sep_ratio_mean = mean(df$sampled_nonedge_sep_ratio),
      sampled_stress_mean = mean(df$sampled_stress),
      axis_deviation_mean = mean(df$axis_deviation),
      axis_deviation_sd = stats::sd(df$axis_deviation),
      elapsed_sec_mean = mean(df$elapsed_sec),
      elapsed_sec_sd = stats::sd(df$elapsed_sec),
      stringsAsFactors = FALSE
    )
  })
)

rank_asc <- function(x) rank(x, ties.method = "average")
rank_desc <- function(x) rank(-x, ties.method = "average")

summary_df$quality_rank_sum <- with(summary_df,
  rank_asc(procrustes_rmse_mean) +
    rank_asc(edge_length_cv_mean) +
    rank_asc(axis_deviation_mean) +
    rank_desc(sampled_nonedge_sep_ratio_mean)
)
summary_df$short_label <- sprintf("%s-%d", substr(summary_df$coarse_mode, 1L, 1L), summary_df$final_rounds)
summary_df <- summary_df[order(summary_df$quality_rank_sum, summary_df$procrustes_rmse_mean), , drop = FALSE]
rownames(summary_df) <- NULL

raw_csv_path <- file.path(output_root, "carpet-level4-finalstage-round2-raw.csv")
summary_csv_path <- file.path(output_root, "carpet-level4-finalstage-round2-summary.csv")
summary_md_path <- file.path(output_root, "carpet-level4-finalstage-round2-summary.md")
line_plot_path <- file.path(output_root, "carpet-level4-finalstage-round2-final-rounds-lines.png")
tradeoff_plot_path <- file.path(output_root, "carpet-level4-finalstage-round2-tradeoff-scatter.png")
sheet_path <- file.path(output_root, "carpet-level4-finalstage-round2-contact-sheet.png")

utils::write.csv(raw_metrics, raw_csv_path, row.names = FALSE)
utils::write.csv(summary_df, summary_csv_path, row.names = FALSE)
draw_final_rounds_lines(line_plot_path, summary_df)
draw_tradeoff_scatter(tradeoff_plot_path, summary_df)

representative_ids <- c("sampled_f0", "sampled_f32", "sampled_f64", "sampled_f96",
                        "sampled_f128", "sampled_f192", "sampled_f384", "exact_f96")
representative_ids <- representative_ids[representative_ids %in% config_df$config_id]
panel_seed <- 1L

panel_specs <- list(
  list(
    coords = canonical,
    title = "Canonical carpet",
    subtitle = sprintf("level=%d | target", level)
  )
)

sampled_f0_trace <- grip.layout.trace(
  edges = edges,
  n = n,
  dim = 2,
  seed = panel_seed,
  final_rounds = 0L,
  coarse_repulsion_exact_below = 64L,
  trace = "round",
  trace.every = 1L
)
true_no_final_idx <- which(sampled_f0_trace$meta$phase == "level_start" &
                             sampled_f0_trace$meta$misf_level == 0L)
if (length(true_no_final_idx) == 0L) {
  stop("Could not locate the pre-final full-graph snapshot for sampled_f0.")
}
true_no_final_coords <- sampled_f0_trace$frames[[true_no_final_idx[[1L]]]]
true_no_final_fit <- helper_env$align_to_target(true_no_final_coords, canonical)
true_no_final_rmse <- true_no_final_fit$rmse

panel_specs[[length(panel_specs) + 1L]] <- list(
  coords = true_no_final_fit$aligned,
  title = "true_no_final",
  subtitle = sprintf("sampled schedule | pre-final snapshot | seed %d RMSE=%s",
                     panel_seed,
                     format_num(true_no_final_rmse, 4L))
)

for (config_id in representative_ids) {
  cfg <- config_df[config_df$config_id == config_id, , drop = FALSE]
  summary_row <- summary_df[summary_df$config_id == config_id, , drop = FALSE]
  if (identical(config_id, "sampled_f0")) {
    coords <- sampled_f0_trace$final
    panel_title <- "sampled_f0 (= 1 FR round)"
  } else {
    coords <- grip.layout(
      edges = edges,
      n = n,
      dim = 2,
      seed = panel_seed,
      final_rounds = cfg$final_rounds[[1L]],
      coarse_repulsion_exact_below = cfg$coarse_repulsion_exact_below[[1L]]
    )
    panel_title <- config_id
  }
  fit <- helper_env$align_to_target(coords, canonical)
  panel_specs[[length(panel_specs) + 1L]] <- list(
    coords = fit$aligned,
    title = panel_title,
    subtitle = sprintf("%s | mean RMSE=%s",
                       cfg$coarse_mode[[1L]],
                       format_num(summary_row$procrustes_rmse_mean, 4L))
  )
}

draw_contact_sheet(sheet_path, panel_specs)

best_idx <- which.min(summary_df$quality_rank_sum)
best_row <- summary_df[best_idx, , drop = FALSE]
default_row <- summary_df[summary_df$config_id == "sampled_f384", , drop = FALSE]
no_final_row <- summary_df[summary_df$config_id == "sampled_f0", , drop = FALSE]

lines <- c(
  "# Carpet Level-4 Round-2 Final-Stage Study",
  "",
  sprintf("- graph: level-%d Sierpinski carpet (`%d` vertices, `%d` edges)", level, n, nrow(edges)),
  sprintf("- seeds per configuration: `%d`", length(seeds)),
  "- Stage A: sampled coarse repulsion with `final_rounds = 0, 32, 64, 96, 128, 160, 192, 256, 384`",
  "- Stage B: exact coarse repulsion with `final_rounds = 0, 64, 96, 128, 160, 192`",
  "- rectilinearity proxy: mean axis deviation after Procrustes alignment; lower is better",
  "",
  "Top configurations by quality-rank sum:",
  "",
  "| Config | Coarse mode | final_rounds | Mean RMSE | Mean edge CV | Mean axis dev | Mean sep ratio | Mean sec | Quality rank sum |",
  "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"
)

for (i in seq_len(min(8L, nrow(summary_df)))) {
  row <- summary_df[i, , drop = FALSE]
  lines <- c(lines, sprintf(
    "| %s | %s | %d | %s | %s | %s | %s | %s | %s |",
    row$config_id,
    row$coarse_mode,
    row$final_rounds,
    format_num(row$procrustes_rmse_mean),
    format_num(row$edge_length_cv_mean),
    format_num(row$axis_deviation_mean),
    format_num(row$sampled_nonedge_sep_ratio_mean),
    format_num(row$elapsed_sec_mean, 3L),
    format_num(row$quality_rank_sum, 1L)
  ))
}

lines <- c(
  lines,
  "",
  "Reference rows:",
  "",
  sprintf("- default sampled_f384: RMSE `%s`, edge CV `%s`, axis deviation `%s`, sep ratio `%s`, sec `%s`",
          format_num(default_row$procrustes_rmse_mean),
          format_num(default_row$edge_length_cv_mean),
          format_num(default_row$axis_deviation_mean),
          format_num(default_row$sampled_nonedge_sep_ratio_mean),
          format_num(default_row$elapsed_sec_mean, 3L)),
  sprintf("- sampled_f0 (= one final FR round): RMSE `%s`, edge CV `%s`, axis deviation `%s`, sep ratio `%s`, sec `%s`",
          format_num(no_final_row$procrustes_rmse_mean),
          format_num(no_final_row$edge_length_cv_mean),
          format_num(no_final_row$axis_deviation_mean),
          format_num(no_final_row$sampled_nonedge_sep_ratio_mean),
          format_num(no_final_row$elapsed_sec_mean, 3L)),
  sprintf("- true_no_final pre-final snapshot (seed %d): RMSE `%s`",
          panel_seed,
          format_num(true_no_final_rmse, 4L)),
  sprintf("- best quality-rank config: `%s` (%s, final_rounds=%d)",
          best_row$config_id,
          best_row$coarse_mode,
          best_row$final_rounds),
  "",
  sprintf("- raw CSV: `%s`", raw_csv_path),
  sprintf("- summary CSV: `%s`", summary_csv_path),
  sprintf("- final-round line plot: `%s`", line_plot_path),
  sprintf("- tradeoff scatter: `%s`", tradeoff_plot_path),
  sprintf("- contact sheet: `%s`", sheet_path)
)

writeLines(lines, con = summary_md_path)

message(sprintf("Raw metrics written to %s", raw_csv_path))
message(sprintf("Summary metrics written to %s", summary_csv_path))
message(sprintf("Markdown summary written to %s", summary_md_path))
message(sprintf("Line plot written to %s", line_plot_path))
message(sprintf("Tradeoff plot written to %s", tradeoff_plot_path))
message(sprintf("Contact sheet written to %s", sheet_path))
