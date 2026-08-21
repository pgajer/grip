#!/usr/bin/env Rscript

output_root <- file.path("output", "gkk_lgkk_paper", "tmp", "carpet-level4-level0-insertion-ablation")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'devtools' or the 'grip' package to run this ablation study.")
}

if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Package 'igraph' is required to run this ablation study.")
}

helper_env <- new.env(parent = baseenv())
sys.source(file.path("tools", "benchmarks", "gkk_lgkk_paper", "benchmark-sierpinski-baseline.R"), envir = helper_env)

level <- 4L
seed <- 1L

built <- helper_env$build_sierpinski_carpet(level)
edges <- built$edges
canonical <- built$coords
n <- nrow(canonical)
graph <- igraph::graph_from_edgelist(as.matrix(edges), directed = FALSE)

prepared_lgkk <- prepare.landmark.geodesic.kk(
  edges = edges,
  n = n,
  local_nbrs = 20L,
  landmark_count = 8L
)

format_num <- function(x, digits = 4L) {
  ifelse(is.finite(x), format(round(x, digits), nsmall = digits, trim = TRUE), "NA")
}

score_layout <- function(coords, variant, stage) {
  geom <- geometry.diagnostics(
    coords = coords,
    target.coords = canonical,
    edges = edges,
    family = "sierpinski.carpet",
    sample.size.symmetry = 256L,
    sample.size.wedges = 1000L,
    rng.seed = 1L
  )
  lgkk <- score.landmark.geodesic.kk(
    coords = coords,
    prepared = prepared_lgkk
  )
  edge_stats <- helper_env$edge_length_stats(coords, edges)
  data.frame(
    variant = variant,
    stage = stage,
    procrustes.rmse = geom$procrustes.rmse[[1L]],
    global.symmetry.score = geom$global.symmetry.score[[1L]],
    local.angle.deviation = geom$local.angle.deviation[[1L]],
    edge.axis.deviation = geom$edge.axis.deviation[[1L]],
    edge.length.cv = edge_stats$cv,
    boundary.waviness = geom$boundary.waviness[[1L]],
    corridor.waviness = geom$corridor.waviness[[1L]],
    hole.center.error = geom$hole.center.error[[1L]],
    central.hole.skew = geom$central.hole.skew[[1L]],
    central.hole.aspect.error = geom$central.hole.aspect.error[[1L]],
    central.hole.center.error = geom$central.hole.center.error[[1L]],
    lgkk.weighted.rel.rmse = lgkk$lgkk.weighted.rel.rmse[[1L]],
    stringsAsFactors = FALSE
  )
}

run_variant <- function(name, args = list()) {
  tr <- do.call(
    trace.grip,
    c(
      list(
        edges = edges,
        n = n,
        dim = 2,
        placement = "barycenter",
        trace = "round",
        trace.every = 1L,
        seed = seed
      ),
      args
    )
  )

  full_idx <- which(tr$meta$active_vertices == n)
  if (length(full_idx) == 0L) {
    stop(sprintf("Variant '%s' never reached a fully active frame.", name))
  }
  first_full_idx <- full_idx[[1L]]
  final_idx <- length(tr$frames)

  first_full <- tr$frames[[first_full_idx]]
  final <- tr$final

  list(
    name = name,
    trace = tr,
    first_full_idx = first_full_idx,
    final_idx = final_idx,
    first_full = first_full,
    final = final,
    scores = rbind(
      score_layout(first_full, name, "first_full"),
      score_layout(final, name, "final")
    )
  )
}

variants <- list(
  baseline = list(
    label = "Baseline",
    args = list(
      level0_insertion_mode = "inherit",
      level0_anchor_count = 3L,
      level0_local_kk_steps = 3L
    )
  ),
  c1_barycenter_override = list(
    label = "C1: barycenter-only",
    args = list(
      level0_insertion_mode = "barycenter",
      level0_anchor_count = 3L,
      level0_local_kk_steps = 3L
    )
  ),
  c2_no_local_kk = list(
    label = "C2: no local KK",
    args = list(
      level0_insertion_mode = "inherit",
      level0_anchor_count = 3L,
      level0_local_kk_steps = 0L
    )
  ),
  c3_least_squares = list(
    label = "C3: LS + 6 anchors",
    args = list(
      level0_insertion_mode = "least_squares",
      level0_anchor_count = 6L,
      level0_local_kk_steps = 3L
    )
  )
)

results <- lapply(names(variants), function(name) {
  spec <- variants[[name]]
  out <- run_variant(name, spec$args)
  out$label <- spec$label
  out
})
names(results) <- names(variants)

kk_coords <- igraph::layout_with_kk(graph)
kk_scores <- score_layout(kk_coords, "igraph_kk", "final")

metrics <- do.call(rbind, lapply(results, `[[`, "scores"))
metrics <- rbind(
  score_layout(canonical, "canonical", "reference"),
  metrics,
  kk_scores
)
row.names(metrics) <- NULL

variant_labels <- c(
  canonical = "Canonical",
  baseline = variants$baseline$label,
  c1_barycenter_override = variants$c1_barycenter_override$label,
  c2_no_local_kk = variants$c2_no_local_kk$label,
  c3_least_squares = variants$c3_least_squares$label,
  igraph_kk = "igraph::KK"
)
metrics$variant_label <- unname(variant_labels[metrics$variant])

contact_sheet_path <- file.path(output_root, "carpet-level4-level0-insertion-contact-sheet.png")
comparison_plot_path <- file.path(output_root, "carpet-level4-level0-insertion-metrics.png")
csv_path <- file.path(output_root, "carpet-level4-level0-insertion-metrics.csv")
summary_path <- file.path(output_root, "carpet-level4-level0-insertion-summary.md")

make_panel <- function(coords, title, subtitle) {
  fit <- helper_env$align_to_target(coords, canonical)
  list(
    coords = fit$aligned,
    edges = edges,
    title = title,
    subtitle = subtitle
  )
}

panels <- list(
  make_panel(canonical, "Canonical", sprintf("level=%d carpet", level)),
  make_panel(kk_coords, "igraph::KK", sprintf(
    "RMSE %s | skew %s",
    format_num(kk_scores$procrustes.rmse[[1L]], 4L),
    format_num(kk_scores$central.hole.skew[[1L]], 4L)
  ))
)
for (name in names(results)) {
  res <- results[[name]]
  first_row <- subset(metrics, variant == name & stage == "first_full")
  final_row <- subset(metrics, variant == name & stage == "final")
  panels[[length(panels) + 1L]] <- make_panel(
    res$first_full,
    paste0(res$label, " | first full"),
    sprintf(
      "RMSE %s | skew %s",
      format_num(first_row$procrustes.rmse[[1L]], 4L),
      format_num(first_row$central.hole.skew[[1L]], 4L)
    )
  )
  panels[[length(panels) + 1L]] <- make_panel(
    res$final,
    paste0(res$label, " | final"),
    sprintf(
      "RMSE %s | skew %s",
      format_num(final_row$procrustes.rmse[[1L]], 4L),
      format_num(final_row$central.hole.skew[[1L]], 4L)
    )
  )
}

grDevices::png(contact_sheet_path, width = 3000, height = 1800, res = 180, bg = "#f7f3ea")
graphics::par(mfrow = c(2, 5), mar = c(0, 0, 3.3, 0), xaxs = "i", yaxs = "i")
for (panel in panels) {
  helper_env$plot_layout_panel(
    panel$coords,
    panel$edges,
    title_text = panel$title,
    subtitle_text = panel$subtitle
  )
}
grDevices::dev.off()

plot_df <- subset(metrics, variant %in% names(results) & stage %in% c("first_full", "final"))
plot_df$variant_label <- factor(
  unname(variant_labels[plot_df$variant]),
  levels = unname(variant_labels[names(results)])
)
stage_colors <- c(first_full = "#355070", final = "#b23a48")

grDevices::png(comparison_plot_path, width = 2200, height = 1200, res = 180, bg = "#f7f3ea")
graphics::par(mfrow = c(1, 2), mar = c(8, 4, 3, 1))
for (metric_name in c("procrustes.rmse", "central.hole.skew")) {
  y <- plot_df[[metric_name]]
  ymax <- max(y, na.rm = TRUE) * 1.12
  graphics::plot(
    x = c(1, length(levels(plot_df$variant_label))),
    y = c(0, ymax),
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = metric_name,
    main = if (metric_name == "procrustes.rmse") "RMSE by stage" else "Central-hole skew by stage"
  )
  graphics::axis(2)
  graphics::axis(1,
                 at = seq_along(levels(plot_df$variant_label)),
                 labels = levels(plot_df$variant_label),
                 las = 2)
  for (i in seq_along(levels(plot_df$variant_label))) {
    rows <- plot_df$variant_label == levels(plot_df$variant_label)[[i]]
    rows_first <- rows & plot_df$stage == "first_full"
    rows_final <- rows & plot_df$stage == "final"
    x_first <- i - 0.1
    x_final <- i + 0.1
    graphics::points(x_first, plot_df[[metric_name]][rows_first], pch = 16, cex = 1.4, col = stage_colors[["first_full"]])
    graphics::points(x_final, plot_df[[metric_name]][rows_final], pch = 17, cex = 1.4, col = stage_colors[["final"]])
    graphics::segments(
      x0 = x_first,
      y0 = plot_df[[metric_name]][rows_first],
      x1 = x_final,
      y1 = plot_df[[metric_name]][rows_final],
      col = grDevices::adjustcolor("#444444", alpha.f = 0.7),
      lwd = 1.5
    )
  }
  graphics::legend(
    "topright",
    legend = c("first_full", "final"),
    col = unname(stage_colors),
    pch = c(16, 17),
    bty = "n",
    cex = 0.9
  )
}
grDevices::dev.off()

utils::write.csv(metrics, csv_path, row.names = FALSE)

baseline_first <- subset(metrics, variant == "baseline" & stage == "first_full")
baseline_final <- subset(metrics, variant == "baseline" & stage == "final")
best_first_skew <- subset(
  plot_df[plot_df$stage == "first_full", ],
  central.hole.skew == min(plot_df$central.hole.skew[plot_df$stage == "first_full"])
)[1L, , drop = FALSE]
best_final_skew <- subset(
  plot_df[plot_df$stage == "final", ],
  central.hole.skew == min(plot_df$central.hole.skew[plot_df$stage == "final"])
)[1L, , drop = FALSE]
best_final_rmse <- subset(
  plot_df[plot_df$stage == "final", ],
  procrustes.rmse == min(plot_df$procrustes.rmse[plot_df$stage == "final"])
)[1L, , drop = FALSE]

writeLines(
  c(
    "# Level-0 Insertion Skewness Ablation on Sierpinski Carpet Level 4",
    "",
    sprintf("- Seed: `%d`", seed),
    "- Baseline uses the current `trace.grip()` default profile in 2D with explicit level-0 defaults (`inherit`, 3 anchors, 3 local KK steps).",
    "- Candidate 1 forces barycentric level-0 placement.",
    "- Candidate 2 removes the level-0 local KK micro-polish.",
    "- Candidate 3 uses level-0 least-squares placement with 6 anchors.",
    "- `igraph::KK` is included as a full-layout reference.",
    "",
    "## Highlights",
    "",
    sprintf(
      "- Baseline first fully active frame: RMSE `%s`, central-hole skew `%s`.",
      format_num(baseline_first$procrustes.rmse[[1L]], 4L),
      format_num(baseline_first$central.hole.skew[[1L]], 4L)
    ),
    sprintf(
      "- Baseline final layout: RMSE `%s`, central-hole skew `%s`.",
      format_num(baseline_final$procrustes.rmse[[1L]], 4L),
      format_num(baseline_final$central.hole.skew[[1L]], 4L)
    ),
    sprintf(
      "- Best first-full skew among GRIP variants: `%s` with skew `%s`.",
      best_first_skew$variant_label[[1L]],
      format_num(best_first_skew$central.hole.skew[[1L]], 4L)
    ),
    sprintf(
      "- Best final skew among GRIP variants: `%s` with skew `%s`.",
      best_final_skew$variant_label[[1L]],
      format_num(best_final_skew$central.hole.skew[[1L]], 4L)
    ),
    sprintf(
      "- Best final RMSE among GRIP variants: `%s` with RMSE `%s`.",
      best_final_rmse$variant_label[[1L]],
      format_num(best_final_rmse$procrustes.rmse[[1L]], 4L)
    ),
    sprintf(
      "- `igraph::KK` final reference: RMSE `%s`, central-hole skew `%s`.",
      format_num(kk_scores$procrustes.rmse[[1L]], 4L),
      format_num(kk_scores$central.hole.skew[[1L]], 4L)
    ),
    "",
    "## Files",
    "",
    sprintf("- Contact sheet: `%s`", contact_sheet_path),
    sprintf("- Paired metric plot: `%s`", comparison_plot_path),
    sprintf("- Metrics CSV: `%s`", csv_path)
  ),
  con = summary_path
)
