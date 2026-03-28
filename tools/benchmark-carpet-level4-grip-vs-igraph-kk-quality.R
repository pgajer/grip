#!/usr/bin/env Rscript

output_root <- file.path(
  "dev", "manual", "tmp", "carpet-level4-grip-vs-igraph-kk-quality-20seeds"
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'devtools' or the 'grip' package to run the carpet quality benchmark.")
}

if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Package 'igraph' is required to run the carpet quality benchmark.")
}

helper_env <- new.env(parent = baseenv())
sys.source(file.path("tools", "benchmark-sierpinski-baseline.R"), envir = helper_env)

seeds <- 1:20
level <- 4L
built <- helper_env$build_sierpinski_carpet(level)
edges <- built$edges
canonical <- built$coords
n <- nrow(canonical)
graph <- igraph::graph_from_edgelist(as.matrix(edges), directed = FALSE)

run_one_layout <- function(method_id, seed) {
  if (identical(method_id, "grip.layout")) {
    started <- proc.time()[["elapsed"]]
    coords <- grip.layout(edges = edges, n = n, dim = 2, seed = seed)
    elapsed <- proc.time()[["elapsed"]] - started
  } else if (identical(method_id, "igraph_kk")) {
    set.seed(seed)
    started <- proc.time()[["elapsed"]]
    coords <- igraph::layout_with_kk(graph)
    elapsed <- proc.time()[["elapsed"]] - started
  } else {
    stop(sprintf("Unknown method_id: %s", method_id))
  }

  aligned <- helper_env$align_to_target(coords, canonical)
  data.frame(
    method_id = method_id,
    method_label = if (identical(method_id, "grip.layout")) "grip.layout()" else "igraph::layout_with_kk()",
    seed = seed,
    procrustes_rmse = aligned$rmse,
    elapsed_sec = elapsed,
    stringsAsFactors = FALSE
  )
}

format_num <- function(x, digits = 4L) {
  ifelse(is.finite(x), format(round(x, digits), nsmall = digits, trim = TRUE), "NA")
}

raw_metrics <- do.call(
  rbind,
  lapply(c("grip.layout", "igraph_kk"), function(method_id) {
    do.call(rbind, lapply(seeds, function(seed) {
      message(sprintf("Running %s on carpet level %d, seed %d/%d...",
                      method_id, level, seed, length(seeds)))
      run_one_layout(method_id, seed)
    }))
  })
)

raw_metrics <- raw_metrics[order(raw_metrics$method_label, raw_metrics$seed), , drop = FALSE]

summary_metrics <- do.call(
  rbind,
  lapply(split(raw_metrics, raw_metrics$method_label), function(df) {
    best_idx <- which.min(df$procrustes_rmse)
    worst_idx <- which.max(df$procrustes_rmse)
    data.frame(
      method_label = df$method_label[[1L]],
      seeds_tested = nrow(df),
      procrustes_rmse_mean = mean(df$procrustes_rmse),
      procrustes_rmse_sd = stats::sd(df$procrustes_rmse),
      procrustes_rmse_median = stats::median(df$procrustes_rmse),
      procrustes_rmse_min = min(df$procrustes_rmse),
      procrustes_rmse_max = max(df$procrustes_rmse),
      best_seed = df$seed[[best_idx]],
      worst_seed = df$seed[[worst_idx]],
      elapsed_sec_mean = mean(df$elapsed_sec),
      elapsed_sec_sd = stats::sd(df$elapsed_sec),
      stringsAsFactors = FALSE
    )
  })
)
summary_metrics <- summary_metrics[order(summary_metrics$procrustes_rmse_mean), , drop = FALSE]

plot_path <- file.path(output_root, "carpet-level4-grip-vs-igraph-kk-procrustes-boxplot.png")
grDevices::png(plot_path, width = 1600, height = 1000, res = 160, bg = "#f7f3ea")
graphics::par(mar = c(6, 5, 4, 1))
graphics::boxplot(
  procrustes_rmse ~ method_label,
  data = raw_metrics,
  border = "#16324f",
  col = c("#f05a28", "#0f3b5f"),
  outline = FALSE,
  ylab = "Procrustes RMSE to canonical carpet",
  xlab = "",
  main = "Level-4 Sierpinski Carpet: grip.layout() vs igraph::layout_with_kk()"
)
graphics::stripchart(
  procrustes_rmse ~ method_label,
  data = raw_metrics,
  method = "jitter",
  vertical = TRUE,
  add = TRUE,
  pch = 21,
  bg = grDevices::adjustcolor("#f7f3ea", alpha.f = 0.85),
  col = "#16324f"
)
graphics::mtext("20 seeds per method", side = 3, line = 0.6, col = "#466074")
grDevices::dev.off()

raw_csv_path <- file.path(output_root, "carpet-level4-grip-vs-igraph-kk-raw.csv")
summary_csv_path <- file.path(output_root, "carpet-level4-grip-vs-igraph-kk-summary.csv")
summary_md_path <- file.path(output_root, "carpet-level4-grip-vs-igraph-kk-summary.md")
utils::write.csv(raw_metrics, raw_csv_path, row.names = FALSE)
utils::write.csv(summary_metrics, summary_csv_path, row.names = FALSE)

lines <- c(
  "# Level-4 Carpet Quality: grip.layout vs igraph KK",
  "",
  sprintf("- graph: Sierpinski carpet level %d (`%d` vertices, `%d` edges)", level, n, nrow(edges)),
  sprintf("- seeds tested per method: `%d`", length(seeds)),
  "- quality metric: Procrustes RMSE to the canonical carpet embedding",
  "- igraph layout: `igraph::layout_with_kk()` with default settings and `set.seed(seed)` before each run",
  "- grip layout: `grip.layout()` with the current primary defaults",
  "",
  "Summary:",
  "",
  "| Method | Seeds | Mean RMSE | SD RMSE | Median RMSE | Best RMSE | Worst RMSE | Best seed | Worst seed | Mean elapsed sec |",
  "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"
)

for (i in seq_len(nrow(summary_metrics))) {
  row <- summary_metrics[i, , drop = FALSE]
  lines <- c(lines, sprintf(
    "| %s | %d | %s | %s | %s | %s | %s | %d | %d | %s |",
    row$method_label,
    row$seeds_tested,
    format_num(row$procrustes_rmse_mean),
    format_num(row$procrustes_rmse_sd),
    format_num(row$procrustes_rmse_median),
    format_num(row$procrustes_rmse_min),
    format_num(row$procrustes_rmse_max),
    row$best_seed,
    row$worst_seed,
    format_num(row$elapsed_sec_mean, 3L)
  ))
}

lines <- c(
  lines,
  "",
  sprintf("- distribution plot: `%s`", plot_path)
)
writeLines(lines, con = summary_md_path)

message(sprintf("Raw metrics written to %s", raw_csv_path))
message(sprintf("Summary metrics written to %s", summary_csv_path))
message(sprintf("Markdown report written to %s", summary_md_path))
message(sprintf("Distribution plot written to %s", plot_path))
