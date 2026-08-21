#!/usr/bin/env Rscript

output_root <- file.path("output", "gkk_lgkk_paper", "tmp", "landmark-geodesic-kk-carpet-level4")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'devtools' or the 'grip' package to run this evaluation.")
}

if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Package 'igraph' is required to run this evaluation.")
}

helper_env <- new.env(parent = baseenv())
sys.source(file.path("tools", "benchmarks", "gkk_lgkk_paper", "benchmark-sierpinski-baseline.R"), envir = helper_env)

level <- 4L
seed <- 1L
local_nbrs <- 20L
landmark_count <- 8L

built <- helper_env$build_sierpinski_carpet(level)
edges <- built$edges
canonical <- built$coords
n <- nrow(canonical)
graph <- igraph::graph_from_edgelist(as.matrix(edges), directed = FALSE)

prepared <- prepare.landmark.geodesic.kk(
  edges = edges,
  n = n,
  local_nbrs = local_nbrs,
  landmark_count = landmark_count
)

trace_fr <- trace.grip(
  edges = edges,
  n = n,
  dim = 2,
  seed = seed,
  final_mode = "fr",
  trace = "round",
  trace.every = 1L
)
pre_final_idx <- which(trace_fr$meta$phase == "level_start" &
                         trace_fr$meta$misf_level == 0L)
if (length(pre_final_idx) == 0L) {
  stop("Could not locate the pre-final full-graph snapshot in the trace.")
}

coords_true_no_final <- trace_fr$frames[[pre_final_idx[[1L]]]]
coords_f1 <- grip(
  edges = edges,
  n = n,
  dim = 2,
  seed = seed,
  final_mode = "fr",
  final_rounds = 1L
)
coords_default <- grip(
  edges = edges,
  n = n,
  dim = 2,
  seed = seed
)
coords_best_struct <- grip(
  edges = edges,
  n = n,
  dim = 2,
  seed = seed,
  final_mode = "fr",
  final_rounds = 8L,
  final_anchor_factor = 1.0,
  final_move_scale_after_first = 0.25
)
coords_kk <- igraph::layout_with_kk(graph)

fmt <- function(x, digits = 4L) {
  ifelse(is.finite(x), format(round(x, digits), nsmall = digits, trim = TRUE), "NA")
}

score_one <- function(candidate, title, coords) {
  lgkk <- score.landmark.geodesic.kk(coords, prepared = prepared)
  geom <- geometry.diagnostics(
    coords = coords,
    target.coords = canonical,
    edges = edges,
    family = "sierpinski.carpet",
    rng.seed = 1L
  )
  cbind(
    data.frame(
      candidate = candidate,
      title = title,
      stringsAsFactors = FALSE
    ),
    lgkk,
    geom,
    stringsAsFactors = FALSE
  )
}

scores <- do.call(rbind, list(
  score_one("canonical", "Canonical", canonical),
  score_one("true_no_final", "True No Final", coords_true_no_final),
  score_one("f1", "f1 (1 FR round)", coords_f1),
  score_one("best_struct_round3", "Best Round-3 Structural", coords_best_struct),
  score_one("grip_default", "Current grip()", coords_default),
  score_one("igraph_kk", "igraph::KK", coords_kk)
))
rownames(scores) <- NULL
scores <- scores[order(scores$lgkk.weighted.rel.rmse, scores$procrustes.rmse), , drop = FALSE]

panel_map <- list(
  canonical = canonical,
  true_no_final = coords_true_no_final,
  f1 = coords_f1,
  best_struct_round3 = coords_best_struct,
  grip_default = coords_default,
  igraph_kk = coords_kk
)

panel_specs <- lapply(seq_len(nrow(scores)), function(i) {
  row <- scores[i, , drop = FALSE]
  fit <- helper_env$align_to_target(panel_map[[row$candidate[[1L]]]], canonical)
  list(
    coords = fit$aligned,
    title = row$title[[1L]],
    subtitle = paste0(
      "relLGKK ", fmt(row$lgkk.weighted.rel.rmse[[1L]], 4L),
      " | RMSE ", fmt(row$procrustes.rmse[[1L]], 4L)
    )
  )
})

sheet_path <- file.path(output_root, "landmark-geodesic-kk-carpet-level4-contact-sheet.png")
rank_plot_path <- file.path(output_root, "landmark-geodesic-kk-carpet-level4-ranking.png")
csv_path <- file.path(output_root, "landmark-geodesic-kk-carpet-level4-scores.csv")
md_path <- file.path(output_root, "landmark-geodesic-kk-carpet-level4-summary.md")

grDevices::png(sheet_path, width = 2500, height = 1700, res = 180, bg = "#f7f3ea")
graphics::par(mfrow = c(2, 3), mar = c(0, 0, 3.3, 0), xaxs = "i", yaxs = "i")
for (panel in panel_specs) {
  helper_env$plot_layout_panel(
    panel$coords,
    edges,
    title_text = panel$title,
    subtitle_text = panel$subtitle
  )
}
grDevices::dev.off()

plot_df <- scores
plot_df$candidate_label <- plot_df$title
grDevices::png(rank_plot_path, width = 1800, height = 1200, res = 170, bg = "#f7f3ea")
graphics::par(mfrow = c(1, 2), mar = c(6.5, 4.5, 3, 1), oma = c(0, 0, 1.2, 0))
graphics::barplot(
  height = rev(plot_df$lgkk.weighted.rel.rmse),
  names.arg = rev(plot_df$candidate_label),
  horiz = TRUE,
  las = 1,
  col = "#355070",
  border = NA,
  main = "Relative landmark geodesic KK error",
  xlab = "smaller is better"
)
graphics::plot(
  plot_df$procrustes.rmse,
  plot_df$lgkk.weighted.rel.rmse,
  pch = 19,
  cex = 1.4,
  col = "#b23a48",
  xlab = "Procrustes RMSE to canonical",
  ylab = "Relative landmark geodesic KK error",
  main = "Relative LGKK error vs RMSE"
)
graphics::text(
  plot_df$procrustes.rmse,
  plot_df$lgkk.weighted.rel.rmse,
  labels = plot_df$candidate,
  pos = 4,
  cex = 0.9,
  col = "#16324f",
  offset = 0.5
)
graphics::mtext(
  sprintf("Level-%d carpet | local_nbrs=%d | landmark_count=%d", level, local_nbrs, landmark_count),
  side = 3,
  outer = TRUE,
  cex = 1.1,
  font = 2
)
grDevices::dev.off()

utils::write.csv(scores, csv_path, row.names = FALSE)

writeLines(
  c(
    "# Landmark Geodesic KK Evaluator: Level-4 Carpet",
    "",
    "## Setup",
    "",
    sprintf("- graph: `Sierpinski carpet level %d`", level),
    sprintf("- seed: `%d`", seed),
    sprintf("- local_nbrs: `%d`", local_nbrs),
    sprintf("- landmark_count: `%d`", landmark_count),
    "- candidate set: canonical, true pre-final GRIP snapshot, `f1`, best round-3 structural GRIP, current `grip()`, and `igraph::layout_with_kk()`.",
    "",
    "## Outputs",
    "",
    sprintf("- contact sheet: `%s`", sheet_path),
    sprintf("- ranking plot: `%s`", rank_plot_path),
    sprintf("- scores CSV: `%s`", csv_path),
    "",
    "## Ranked scores",
    "",
    paste(capture.output(print(scores[, c(
      "candidate", "lgkk.weighted.rel.rmse", "lgkk.mean.rel.path.error",
      "lgkk.energy", "lgkk.weighted.rmse", "scale.L0",
      "procrustes.rmse", "edge.axis.deviation", "boundary.waviness",
      "corridor.waviness", "hole.center.error"
    )])), collapse = "\n")
  ),
  con = md_path
)

message("Wrote:")
message("  ", sheet_path)
message("  ", rank_plot_path)
message("  ", csv_path)
message("  ", md_path)
