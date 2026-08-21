#!/usr/bin/env Rscript

output_root <- file.path("output", "gkk_lgkk_paper", "tmp", "carpet-level4-level0-ls-reduced-kk")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'devtools' or the 'grip' package to run this study.")
}

if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Package 'igraph' is required to run this study.")
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
  lgkk <- score.landmark.geodesic.kk(coords, prepared = prepared_lgkk)
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
  if (!length(full_idx)) {
    stop(sprintf("Variant '%s' never reached a fully active frame.", name))
  }
  first_full <- tr$frames[[full_idx[[1L]]]]
  final <- tr$final
  list(
    trace = tr,
    first_full = first_full,
    final = final,
    scores = rbind(
      score_layout(first_full, name, "first_full"),
      score_layout(final, name, "final")
    )
  )
}

reference_variants <- list(
  baseline = list(
    label = "Baseline",
    args = list(
      level0_insertion_mode = "inherit",
      level0_anchor_count = 3L,
      level0_local_kk_steps = 3L
    )
  ),
  no_local_kk = list(
    label = "No local KK",
    args = list(
      level0_insertion_mode = "inherit",
      level0_anchor_count = 3L,
      level0_local_kk_steps = 0L
    )
  ),
  ls_6_3 = list(
    label = "LS a6 k3",
    args = list(
      level0_insertion_mode = "least_squares",
      level0_anchor_count = 6L,
      level0_local_kk_steps = 3L
    )
  )
)

anchors <- c(6L, 8L, 10L)
kk_steps <- 0:3
grid_names <- c()
grid_specs <- list()
for (a in anchors) {
  for (k in kk_steps) {
    key <- sprintf("ls_a%d_k%d", a, k)
    grid_names <- c(grid_names, key)
    grid_specs[[key]] <- list(
      label = sprintf("LS a%d k%d", a, k),
      args = list(
        level0_insertion_mode = "least_squares",
        level0_anchor_count = a,
        level0_local_kk_steps = k
      )
    )
  }
}

all_specs <- c(reference_variants, grid_specs)
results <- lapply(names(all_specs), function(name) {
  out <- run_variant(name, all_specs[[name]]$args)
  out$label <- all_specs[[name]]$label
  out
})
names(results) <- names(all_specs)

kk_coords <- igraph::layout_with_kk(graph)
kk_scores <- score_layout(kk_coords, "igraph_kk", "final")

metrics <- do.call(rbind, lapply(results, `[[`, "scores"))
metrics <- rbind(
  score_layout(canonical, "canonical", "reference"),
  metrics,
  kk_scores
)
variant_labels <- c(
  canonical = "Canonical",
  igraph_kk = "igraph::KK",
  vapply(all_specs, `[[`, character(1L), "label")
)
metrics$variant_label <- unname(variant_labels[metrics$variant])
row.names(metrics) <- NULL

grid_metrics <- subset(metrics, grepl("^ls_a", variant) & stage %in% c("first_full", "final"))
grid_metrics$anchors <- as.integer(sub("^ls_a([0-9]+)_k([0-9]+)$", "\\1", grid_metrics$variant))
grid_metrics$kk_steps <- as.integer(sub("^ls_a([0-9]+)_k([0-9]+)$", "\\2", grid_metrics$variant))

best_grid_final_rmse <- grid_metrics[grid_metrics$stage == "final", ]
best_grid_final_rmse <- best_grid_final_rmse[which.min(best_grid_final_rmse$procrustes.rmse), , drop = FALSE]
best_grid_final_skew <- grid_metrics[grid_metrics$stage == "final", ]
best_grid_final_skew <- best_grid_final_skew[which.min(best_grid_final_skew$central.hole.skew), , drop = FALSE]
best_grid_first_rmse <- grid_metrics[grid_metrics$stage == "first_full", ]
best_grid_first_rmse <- best_grid_first_rmse[which.min(best_grid_first_rmse$procrustes.rmse), , drop = FALSE]
best_grid_final_lgkk <- grid_metrics[grid_metrics$stage == "final", ]
best_grid_final_lgkk <- best_grid_final_lgkk[which.min(best_grid_final_lgkk$lgkk.weighted.rel.rmse), , drop = FALSE]

csv_path <- file.path(output_root, "carpet-level4-level0-ls-reduced-kk-metrics.csv")
summary_path <- file.path(output_root, "carpet-level4-level0-ls-reduced-kk-summary.md")
heatmap_path <- file.path(output_root, "carpet-level4-level0-ls-reduced-kk-heatmaps.png")
contact_sheet_path <- file.path(output_root, "carpet-level4-level0-ls-reduced-kk-selected-panels.png")

utils::write.csv(metrics, csv_path, row.names = FALSE)

grid_matrix <- function(stage, metric_name) {
  z <- matrix(NA_real_, nrow = length(anchors), ncol = length(kk_steps),
              dimnames = list(paste0("a", anchors), paste0("k", kk_steps)))
  sub <- grid_metrics[grid_metrics$stage == stage, , drop = FALSE]
  for (i in seq_along(anchors)) {
    for (j in seq_along(kk_steps)) {
      row <- sub[sub$anchors == anchors[[i]] & sub$kk_steps == kk_steps[[j]], , drop = FALSE]
      if (nrow(row) == 1L) z[i, j] <- row[[metric_name]][[1L]]
    }
  }
  z
}

draw_heatmap <- function(z, title) {
  nr <- nrow(z)
  nc <- ncol(z)
  z_plot <- z[nr:1, , drop = FALSE]
  palette <- grDevices::colorRampPalette(c("#f7f3ea", "#d8c3a5", "#8c5c3a"))(100)
  finite_vals <- as.vector(z_plot[is.finite(z_plot)])
  rng <- range(finite_vals)
  graphics::plot(
    x = c(0.5, nc + 0.5),
    y = c(0.5, nr + 0.5),
    type = "n",
    axes = FALSE,
    xlab = "level0_local_kk_steps",
    ylab = "level0_anchor_count",
    main = title
  )
  if (length(finite_vals) == 0L) {
    graphics::axis(1, at = seq_len(nc), labels = kk_steps)
    graphics::axis(2, at = seq_len(nr), labels = rev(anchors))
    return(invisible(NULL))
  }
  if (diff(rng) <= 0) {
    rng[[2L]] <- rng[[1L]] + 1e-8
  }
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      val <- z_plot[i, j]
      fill <- "#f7f3ea"
      if (is.finite(val)) {
        idx <- 1L + floor((length(palette) - 1L) * (val - rng[[1L]]) / (rng[[2L]] - rng[[1L]]))
        idx <- max(1L, min(length(palette), idx))
        fill <- palette[[idx]]
      }
      graphics::rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5,
                     col = fill, border = "#d9d2c3")
    }
  }
  graphics::axis(1, at = seq_len(nc), labels = kk_steps)
  graphics::axis(2, at = seq_len(nr), labels = rev(anchors))
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      graphics::text(j, nr - i + 1, labels = format_num(z[i, j], 3L), cex = 0.8)
    }
  }
}

grDevices::png(heatmap_path, width = 2400, height = 1800, res = 180, bg = "#f7f3ea")
graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
draw_heatmap(grid_matrix("first_full", "procrustes.rmse"), "First-full RMSE")
draw_heatmap(grid_matrix("first_full", "central.hole.skew"), "First-full skew")
draw_heatmap(grid_matrix("final", "procrustes.rmse"), "Final RMSE")
draw_heatmap(grid_matrix("final", "central.hole.skew"), "Final skew")
grDevices::dev.off()

panel_variant_keys <- c(
  "baseline",
  "no_local_kk",
  "ls_6_3",
  best_grid_first_rmse$variant[[1L]],
  best_grid_final_rmse$variant[[1L]],
  best_grid_final_skew$variant[[1L]],
  best_grid_final_lgkk$variant[[1L]]
)
panel_variant_keys <- unique(panel_variant_keys)

make_panel <- function(coords, title, subtitle) {
  fit <- helper_env$align_to_target(coords, canonical)
  list(coords = fit$aligned, edges = edges, title = title, subtitle = subtitle)
}

panels <- list(
  make_panel(canonical, "Canonical", sprintf("level=%d carpet", level)),
  make_panel(kk_coords, "igraph::KK", sprintf(
    "RMSE %s | skew %s",
    format_num(kk_scores$procrustes.rmse[[1L]], 4L),
    format_num(kk_scores$central.hole.skew[[1L]], 4L)
  ))
)
for (key in panel_variant_keys) {
  res <- results[[key]]
  row_first <- subset(metrics, variant == key & stage == "first_full")
  row_final <- subset(metrics, variant == key & stage == "final")
  panels[[length(panels) + 1L]] <- make_panel(
    res$first_full,
    paste0(res$label, " | first full"),
    sprintf("RMSE %s | skew %s",
            format_num(row_first$procrustes.rmse[[1L]], 4L),
            format_num(row_first$central.hole.skew[[1L]], 4L))
  )
  panels[[length(panels) + 1L]] <- make_panel(
    res$final,
    paste0(res$label, " | final"),
    sprintf("RMSE %s | skew %s",
            format_num(row_final$procrustes.rmse[[1L]], 4L),
            format_num(row_final$central.hole.skew[[1L]], 4L))
  )
}

n_panels <- length(panels)
ncol_sheet <- 4L
nrow_sheet <- ceiling(n_panels / ncol_sheet)
grDevices::png(contact_sheet_path, width = 3200, height = 800L * nrow_sheet, res = 180, bg = "#f7f3ea")
graphics::par(mfrow = c(nrow_sheet, ncol_sheet), mar = c(0, 0, 3.2, 0), xaxs = "i", yaxs = "i")
for (panel in panels) {
  helper_env$plot_layout_panel(panel$coords, panel$edges,
                               title_text = panel$title,
                               subtitle_text = panel$subtitle)
}
if (n_panels < nrow_sheet * ncol_sheet) {
  for (i in seq_len(nrow_sheet * ncol_sheet - n_panels)) {
    graphics::plot.new()
  }
}
grDevices::dev.off()

writeLines(
  c(
    "# Level-0 Least-Squares + Reduced Local KK Follow-up",
    "",
    sprintf("- Seed: `%d`", seed),
    "- Focused grid: `level0_insertion_mode = \"least_squares\"`, anchors `{6, 8, 10}`, local KK steps `{0, 1, 2, 3}`.",
    "- References carried over from the first ablation: baseline, inherit/no-local-KK, LS with 6 anchors and 3 steps, and `igraph::KK`.",
    "",
    "## Best grid configurations",
    "",
    sprintf("- Best first-full RMSE: `%s` at `%s`.", format_num(best_grid_first_rmse$procrustes.rmse[[1L]], 4L), best_grid_first_rmse$variant_label[[1L]]),
    sprintf("- Best final RMSE: `%s` at `%s`.", format_num(best_grid_final_rmse$procrustes.rmse[[1L]], 4L), best_grid_final_rmse$variant_label[[1L]]),
    sprintf("- Best final skew: `%s` at `%s`.", format_num(best_grid_final_skew$central.hole.skew[[1L]], 4L), best_grid_final_skew$variant_label[[1L]]),
    sprintf("- Best final LGKK relative error: `%s` at `%s`.", format_num(best_grid_final_lgkk$lgkk.weighted.rel.rmse[[1L]], 4L), best_grid_final_lgkk$variant_label[[1L]]),
    "",
    "## References",
    "",
    sprintf("- Baseline final: RMSE `%s`, skew `%s`.", format_num(subset(metrics, variant == "baseline" & stage == "final")$procrustes.rmse[[1L]], 4L), format_num(subset(metrics, variant == "baseline" & stage == "final")$central.hole.skew[[1L]], 4L)),
    sprintf("- No local KK final: RMSE `%s`, skew `%s`.", format_num(subset(metrics, variant == "no_local_kk" & stage == "final")$procrustes.rmse[[1L]], 4L), format_num(subset(metrics, variant == "no_local_kk" & stage == "final")$central.hole.skew[[1L]], 4L)),
    sprintf("- LS a6 k3 final: RMSE `%s`, skew `%s`.", format_num(subset(metrics, variant == "ls_6_3" & stage == "final")$procrustes.rmse[[1L]], 4L), format_num(subset(metrics, variant == "ls_6_3" & stage == "final")$central.hole.skew[[1L]], 4L)),
    sprintf("- `igraph::KK` final: RMSE `%s`, skew `%s`.", format_num(kk_scores$procrustes.rmse[[1L]], 4L), format_num(kk_scores$central.hole.skew[[1L]], 4L)),
    "",
    "## Files",
    "",
    sprintf("- Metrics CSV: `%s`", csv_path),
    sprintf("- Heatmaps: `%s`", heatmap_path),
    sprintf("- Selected panels: `%s`", contact_sheet_path)
  ),
  con = summary_path
)
