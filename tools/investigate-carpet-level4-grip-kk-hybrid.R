#!/usr/bin/env Rscript

output_root <- file.path("dev", "manual", "tmp", "carpet-level4-grip-kk-hybrid-study")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'devtools' or the 'grip' package to run the carpet hybrid study.")
}

if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Package 'igraph' is required to run the carpet hybrid study.")
}

helper_env <- new.env(parent = baseenv())
sys.source(file.path("tools", "benchmark-sierpinski-baseline.R"), envir = helper_env)

level <- 4L
seed <- 1L
trace_every <- 8L
kk_polish_iters <- c(5L, 10L, 20L, 50L, 100L)

built <- helper_env$build_sierpinski_carpet(level)
edges <- built$edges
canonical <- built$coords
n <- nrow(canonical)
graph <- igraph::graph_from_edgelist(as.matrix(edges), directed = FALSE)

time_it <- function(expr) {
  started <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, elapsed = proc.time()[["elapsed"]] - started)
}

score_candidate <- function(candidate_id, candidate_label, coords, elapsed_sec, kk_maxiter = NA_integer_) {
  fit <- helper_env$align_to_target(coords, canonical)
  graph_score <- grip.score.layout(
    coords = coords,
    edges = edges,
    n = n,
    sample.size.stress = 2000L,
    sample.size.nonedge = 5000L,
    edge.crossings = "never"
  )
  data.frame(
    candidate_id = candidate_id,
    candidate_label = candidate_label,
    kk_maxiter = kk_maxiter,
    procrustes_rmse = fit$rmse,
    edge_length_cv = graph_score$edge.length.cv,
    sampled_stress = graph_score$sampled.stress,
    sampled_nonedge_sep_ratio = graph_score$sampled.nonedge.sep.ratio,
    elapsed_sec = elapsed_sec,
    stringsAsFactors = FALSE
  )
}

format_num <- function(x, digits = 4L) {
  ifelse(is.finite(x), format(round(x, digits), nsmall = digits, trim = TRUE), "NA")
}

draw_trace_diagnostics <- function(path, diagnostics_df) {
  level_starts <- diagnostics_df$frame[diagnostics_df$phase %in% c("init", "level_start")]

  grDevices::png(path, width = 1800, height = 1200, res = 170, bg = "#f7f3ea")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

  plot_metric <- function(y, ylab, main, col) {
    graphics::plot(diagnostics_df$frame, y,
                   type = "l",
                   lwd = 2,
                   col = col,
                   xlab = "Trace frame",
                   ylab = ylab,
                   main = main)
    for (x in level_starts) {
      graphics::abline(v = x, col = grDevices::adjustcolor("#16324f", alpha.f = 0.12), lty = 3)
    }
  }

  plot_metric(diagnostics_df$procrustes.rmse,
              ylab = "RMSE",
              main = "Procrustes RMSE to canonical carpet",
              col = "#b23a48")
  plot_metric(diagnostics_df$edge.length.cv,
              ylab = "Edge CV",
              main = "Edge-length variation",
              col = "#355070")
  plot_metric(diagnostics_df$sampled.nonedge.sep.ratio,
              ylab = "Separation ratio",
              main = "Sampled non-edge separation",
              col = "#3d7a57")
  plot_metric(diagnostics_df$active_vertices,
              ylab = "Active vertices",
              main = "Active set size",
              col = "#7c4d8b")

  graphics::mtext("GRIP trace diagnostics, level-4 carpet", side = 3, outer = TRUE, cex = 1.2, font = 2)
}

draw_candidate_metrics <- function(path, metrics_df) {
  plot_df <- metrics_df
  plot_df$x <- seq_len(nrow(plot_df))
  labels <- plot_df$candidate_label

  grDevices::png(path, width = 1800, height = 1200, res = 170, bg = "#f7f3ea")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(7, 4, 3, 1), oma = c(0, 0, 2, 0))

  plot_series <- function(y, ylab, main, col) {
    graphics::plot(plot_df$x, y,
                   type = "b",
                   pch = 19,
                   lwd = 2,
                   col = col,
                   xaxt = "n",
                   xlab = "",
                   ylab = ylab,
                   main = main)
    graphics::axis(1, at = plot_df$x, labels = labels, las = 2, cex.axis = 0.9)
  }

  plot_series(plot_df$procrustes_rmse, "RMSE", "Final Procrustes RMSE", "#b23a48")
  plot_series(plot_df$edge_length_cv, "Edge CV", "Final edge-length variation", "#355070")
  plot_series(plot_df$sampled_nonedge_sep_ratio, "Separation ratio", "Final non-edge separation", "#3d7a57")
  plot_series(plot_df$elapsed_sec, "Seconds", "Elapsed time", "#7c4d8b")

  graphics::mtext("Level-4 carpet candidates: GRIP, GRIP->KK polishes, and pure KK", side = 3, outer = TRUE, cex = 1.2, font = 2)
}

draw_contact_sheet <- function(path, panel_specs) {
  grDevices::png(path, width = 2400, height = 1600, res = 170, bg = "#f7f3ea")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 4), mar = c(0, 0, 3, 0), xaxs = "i", yaxs = "i")

  for (spec in panel_specs) {
    helper_env$plot_layout_panel(
      spec$coords,
      edges,
      title_text = spec$title,
      subtitle_text = spec$subtitle
    )
  }
}

message("Tracing grip.layout() with per-frame diagnostics...")
grip_trace_timed <- time_it({
  grip.layout.trace(
    edges = edges,
    n = n,
    dim = 2,
    trace = "round",
    trace.every = trace_every,
    diagnostics = "light",
    target_coords = canonical,
    diagnostic_sample_size_nonedge = 500L,
    seed = seed
  )
})
grip_trace <- grip_trace_timed$value
grip_coords <- grip_trace$final

message("Timing plain grip.layout()...")
grip_plain <- time_it({
  grip.layout(edges = edges, n = n, dim = 2, seed = seed)
})

message("Running pure igraph KK...")
kk_default <- time_it({
  igraph::layout_with_kk(graph)
})

message("Running GRIP -> short KK warm-start polishes...")
hybrid_runs <- lapply(kk_polish_iters, function(maxiter) {
  message(sprintf("  KK polish maxiter=%d", maxiter))
  timed <- time_it({
    igraph::layout_with_kk(graph, coords = grip_coords, maxiter = maxiter)
  })
  list(maxiter = maxiter, coords = timed$value, elapsed = timed$elapsed)
})

candidate_metrics <- list()
candidate_panels <- list(
  list(
    coords = canonical,
    title = "Canonical carpet",
    subtitle = sprintf("level=%d", level)
  )
)

grip_fit <- helper_env$align_to_target(grip_coords, canonical)
candidate_metrics[[length(candidate_metrics) + 1L]] <- score_candidate(
  candidate_id = "grip",
  candidate_label = "grip",
  coords = grip_coords,
  elapsed_sec = grip_plain$elapsed,
  kk_maxiter = 0L
)
candidate_panels[[length(candidate_panels) + 1L]] <- list(
  coords = grip_fit$aligned,
  title = "grip.layout()",
  subtitle = sprintf("seed=%d | %.3fs", seed, grip_plain$elapsed)
)

for (run in hybrid_runs) {
  fit <- helper_env$align_to_target(run$coords, canonical)
  candidate_metrics[[length(candidate_metrics) + 1L]] <- score_candidate(
    candidate_id = sprintf("grip_to_kk_%d", run$maxiter),
    candidate_label = sprintf("grip->kk-%d", run$maxiter),
    coords = run$coords,
    elapsed_sec = run$elapsed,
    kk_maxiter = run$maxiter
  )
  candidate_panels[[length(candidate_panels) + 1L]] <- list(
    coords = fit$aligned,
    title = sprintf("GRIP -> KK %d", run$maxiter),
    subtitle = sprintf("seed=%d | %.3fs", seed, run$elapsed)
  )
}

kk_fit <- helper_env$align_to_target(kk_default$value, canonical)
candidate_metrics[[length(candidate_metrics) + 1L]] <- score_candidate(
  candidate_id = "kk_default",
  candidate_label = "kk-default",
  coords = kk_default$value,
  elapsed_sec = kk_default$elapsed,
  kk_maxiter = as.integer(50L * n)
)
candidate_panels[[length(candidate_panels) + 1L]] <- list(
  coords = kk_fit$aligned,
  title = "igraph::KK default",
  subtitle = sprintf("default maxiter=%d | %.3fs", 50L * n, kk_default$elapsed)
)

candidate_metrics <- do.call(rbind, candidate_metrics)
candidate_metrics <- candidate_metrics[
  match(c("grip", paste0("grip_to_kk_", kk_polish_iters), "kk_default"), candidate_metrics$candidate_id),
  ,
  drop = FALSE
]
rownames(candidate_metrics) <- NULL

trace_csv_path <- file.path(output_root, "carpet-level4-grip-trace-diagnostics.csv")
metrics_csv_path <- file.path(output_root, "carpet-level4-grip-kk-hybrid-metrics.csv")
trace_plot_path <- file.path(output_root, "carpet-level4-grip-trace-diagnostics.png")
metrics_plot_path <- file.path(output_root, "carpet-level4-grip-kk-hybrid-metrics.png")
sheet_path <- file.path(output_root, "carpet-level4-grip-kk-hybrid-contact-sheet.png")
summary_md_path <- file.path(output_root, "carpet-level4-grip-kk-hybrid-summary.md")

utils::write.csv(grip_trace$diagnostics, trace_csv_path, row.names = FALSE)
utils::write.csv(candidate_metrics, metrics_csv_path, row.names = FALSE)
draw_trace_diagnostics(trace_plot_path, grip_trace$diagnostics)
draw_candidate_metrics(metrics_plot_path, candidate_metrics)
draw_contact_sheet(sheet_path, candidate_panels)

best_idx <- which.min(candidate_metrics$procrustes_rmse)
lines <- c(
  "# Carpet Level-4 GRIP / KK Hybrid Study",
  "",
  sprintf("- graph: level-%d Sierpinski carpet (`%d` vertices, `%d` edges)", level, n, nrow(edges)),
  sprintf("- grip seed: `%d`", seed),
  sprintf("- grip trace: `trace = \"round\"`, `trace.every = %d`, diagnostics=`light`", trace_every),
  sprintf("- hybrid candidates: `%s`", paste(kk_polish_iters, collapse = ", ")),
  "",
  "Final candidate metrics:",
  "",
  "| Candidate | KK maxiter | Procrustes RMSE | Edge CV | Sampled stress | Non-edge sep ratio | Elapsed sec |",
  "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"
)

for (i in seq_len(nrow(candidate_metrics))) {
  row <- candidate_metrics[i, , drop = FALSE]
  lines <- c(lines, sprintf(
    "| %s | %s | %s | %s | %s | %s | %s |",
    row$candidate_label,
    if (is.na(row$kk_maxiter)) "NA" else row$kk_maxiter,
    format_num(row$procrustes_rmse),
    format_num(row$edge_length_cv),
    format_num(row$sampled_stress),
    format_num(row$sampled_nonedge_sep_ratio),
    format_num(row$elapsed_sec, 3L)
  ))
}

lines <- c(
  lines,
  "",
  sprintf("- best RMSE in this run: `%s`", candidate_metrics$candidate_label[[best_idx]]),
  sprintf("- trace diagnostics CSV: `%s`", trace_csv_path),
  sprintf("- hybrid metrics CSV: `%s`", metrics_csv_path),
  sprintf("- trace diagnostics plot: `%s`", trace_plot_path),
  sprintf("- hybrid metric plot: `%s`", metrics_plot_path),
  sprintf("- contact sheet: `%s`", sheet_path)
)
writeLines(lines, con = summary_md_path)

message(sprintf("Trace diagnostics written to %s", trace_csv_path))
message(sprintf("Hybrid metrics written to %s", metrics_csv_path))
message(sprintf("Trace diagnostics plot written to %s", trace_plot_path))
message(sprintf("Hybrid metrics plot written to %s", metrics_plot_path))
message(sprintf("Contact sheet written to %s", sheet_path))
message(sprintf("Markdown report written to %s", summary_md_path))
