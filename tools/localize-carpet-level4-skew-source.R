#!/usr/bin/env Rscript

output_root <- file.path("dev", "manual", "tmp", "carpet-level4-skew-localization")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'devtools' or the 'grip' package to run this localization study.")
}

if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Package 'igraph' is required to run this localization study.")
}

helper_env <- new.env(parent = baseenv())
sys.source(file.path("tools", "benchmark-sierpinski-baseline.R"), envir = helper_env)

level <- 4L
seed <- 1L

built <- helper_env$build_sierpinski_carpet(level)
edges <- built$edges
canonical <- built$coords
n <- nrow(canonical)
graph <- igraph::graph_from_edgelist(as.matrix(edges), directed = FALSE)
adj_list <- grip:::grip.build.adj.from.edges(edges, n)$adj_list

format_num <- function(x, digits = 4L) {
  ifelse(is.finite(x), format(round(x, digits), nsmall = digits, trim = TRUE), "NA")
}

trace_obj <- grip.layout.trace(
  edges = edges,
  n = n,
  dim = 2,
  trace = "round",
  trace.every = 1L,
  diagnostics = "none",
  seed = seed
)

score_frame <- function(frame.coords, meta.row) {
  active <- stats::complete.cases(frame.coords)
  active.count <- sum(active)
  row <- meta.row
  row$procrustes.rmse.active <- NA_real_
  row$global.symmetry.active <- NA_real_
  row$local.angle.deviation.active <- NA_real_
  row$edge.axis.deviation.active <- NA_real_
  row$boundary.waviness <- NA_real_
  row$corridor.waviness <- NA_real_
  row$hole.center.error <- NA_real_
  row$central.hole.skew <- NA_real_
  row$central.hole.aspect.error <- NA_real_
  row$central.hole.center.error <- NA_real_

  if (active.count >= 4L) {
    induced <- grip:::grip.induced.layout.inputs(adj_list, NULL, active)
    coords.active <- frame.coords[induced$vertex.ids, , drop = FALSE]
    target.active <- canonical[induced$vertex.ids, , drop = FALSE]
    geom.active <- grip.geometry.diagnostics(
      coords = coords.active,
      target.coords = target.active,
      edges = induced$edges,
      family = NULL,
      sample.size.symmetry = 256L,
      sample.size.wedges = 1000L,
      rng.seed = 1L
    )
    row$procrustes.rmse.active <- geom.active$procrustes.rmse[[1L]]
    row$global.symmetry.active <- geom.active$global.symmetry.score[[1L]]
    row$local.angle.deviation.active <- geom.active$local.angle.deviation[[1L]]
    row$edge.axis.deviation.active <- geom.active$edge.axis.deviation[[1L]]

    if (active.count == n) {
      geom.full <- grip.geometry.diagnostics(
        coords = coords.active,
        target.coords = canonical,
        edges = edges,
        family = "sierpinski.carpet",
        sample.size.symmetry = 256L,
        sample.size.wedges = 1000L,
        rng.seed = 1L
      )
      row$boundary.waviness <- geom.full$boundary.waviness[[1L]]
      row$corridor.waviness <- geom.full$corridor.waviness[[1L]]
      row$hole.center.error <- geom.full$hole.center.error[[1L]]
      row$central.hole.skew <- geom.full$central.hole.skew[[1L]]
      row$central.hole.aspect.error <- geom.full$central.hole.aspect.error[[1L]]
      row$central.hole.center.error <- geom.full$central.hole.center.error[[1L]]
    }
  }

  row
}

trace_df <- do.call(rbind, lapply(seq_along(trace_obj$frames), function(i) {
  score_frame(trace_obj$frames[[i]], trace_obj$meta[i, , drop = FALSE])
}))
rownames(trace_df) <- NULL

full_rows <- which(trace_df$active_vertices == n)
if (length(full_rows) == 0L) {
  stop("Trace never reaches a fully active frame.")
}

last_coarse_idx <- max(which(trace_df$phase == "level_start" & trace_df$misf_level > 0L))
first_full_idx <- full_rows[[1L]]
first_final_round_idx <- {
  idx <- which(trace_df$misf_level == 0L & trace_df$phase == "round")
  if (length(idx) > 0L) idx[[1L]] else first_full_idx
}
final_idx <- nrow(trace_df)
best_full_idx <- full_rows[[which.min(trace_df$central.hole.skew[full_rows])]]

delta_sym <- c(NA_real_, diff(trace_df$global.symmetry.active))
largest_sym_drop_idx <- which.min(delta_sym)

kk_coords <- igraph::layout_with_kk(graph)
kk_geom <- grip.geometry.diagnostics(
  coords = kk_coords,
  target.coords = canonical,
  edges = edges,
  family = "sierpinski.carpet",
  sample.size.symmetry = 256L,
  sample.size.wedges = 1000L,
  rng.seed = 1L
)

build_panel <- function(frame.idx, title, subtitle = NULL, active_only = FALSE) {
  coords <- trace_obj$frames[[frame.idx]]
  active <- stats::complete.cases(coords)
  if (active_only) {
    induced <- grip:::grip.induced.layout.inputs(adj_list, NULL, active)
    fit <- helper_env$align_to_target(coords[induced$vertex.ids, , drop = FALSE],
                                      canonical[induced$vertex.ids, , drop = FALSE])
    return(list(
      coords = fit$aligned,
      edges = induced$edges,
      title = title,
      subtitle = subtitle
    ))
  }
  fit <- helper_env$align_to_target(coords, canonical)
  list(
    coords = fit$aligned,
    edges = edges,
    title = title,
    subtitle = subtitle
  )
}

panels <- list(
  list(
    coords = canonical,
    edges = edges,
    title = "Canonical",
    subtitle = sprintf("level=%d carpet", level)
  ),
  build_panel(
    last_coarse_idx,
    "Last Coarse Level Start",
    subtitle = sprintf(
      "frame %d | active %d | sym %s",
      trace_df$frame[[last_coarse_idx]],
      trace_df$active_vertices[[last_coarse_idx]],
      format_num(trace_df$global.symmetry.active[[last_coarse_idx]], 3L)
    ),
    active_only = TRUE
  ),
  build_panel(
    first_full_idx,
    "First Full-Graph Frame",
    subtitle = sprintf(
      "skew %s | RMSE %s",
      format_num(trace_df$central.hole.skew[[first_full_idx]], 4L),
      format_num(trace_df$procrustes.rmse.active[[first_full_idx]], 4L)
    )
  ),
  build_panel(
    first_final_round_idx,
    "First Finest-Level Round",
    subtitle = sprintf(
      "skew %s | RMSE %s",
      format_num(trace_df$central.hole.skew[[first_final_round_idx]], 4L),
      format_num(trace_df$procrustes.rmse.active[[first_final_round_idx]], 4L)
    )
  ),
  build_panel(
    final_idx,
    "Final Layout",
    subtitle = sprintf(
      "skew %s | RMSE %s",
      format_num(trace_df$central.hole.skew[[final_idx]], 4L),
      format_num(trace_df$procrustes.rmse.active[[final_idx]], 4L)
    )
  ),
  list(
    coords = helper_env$align_to_target(kk_coords, canonical)$aligned,
    edges = edges,
    title = "igraph::KK",
    subtitle = sprintf(
      "skew %s | RMSE %s",
      format_num(kk_geom$central.hole.skew[[1L]], 4L),
      format_num(kk_geom$procrustes.rmse[[1L]], 4L)
    )
  )
)

trace_plot_path <- file.path(output_root, "carpet-level4-skew-localization-trace.png")
sheet_path <- file.path(output_root, "carpet-level4-skew-localization-contact-sheet.png")
csv_path <- file.path(output_root, "carpet-level4-skew-localization-trace.csv")
md_path <- file.path(output_root, "carpet-level4-skew-localization-summary.md")

grDevices::png(trace_plot_path, width = 2200, height = 1500, res = 180, bg = "#f7f3ea")
graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
level_starts <- trace_df$frame[trace_df$phase %in% c("init", "level_start")]
plot_metric <- function(y, ylab, main, col) {
  graphics::plot(trace_df$frame, y,
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
plot_metric(trace_df$procrustes.rmse.active,
            "RMSE", "Active-subgraph RMSE", "#b23a48")
plot_metric(trace_df$global.symmetry.active,
            "Symmetry score", "Active-subgraph symmetry", "#355070")
plot_metric(trace_df$central.hole.skew,
            "Central-hole skew", "Full-graph central-hole skew", "#3d7a57")
plot_metric(trace_df$central.hole.aspect.error,
            "Aspect error", "Full-graph central-hole aspect error", "#7c4d8b")
graphics::mtext("Level-4 carpet skew-source localization trace", side = 3, outer = TRUE, cex = 1.2, font = 2)
grDevices::dev.off()

grDevices::png(sheet_path, width = 2500, height = 1700, res = 180, bg = "#f7f3ea")
graphics::par(mfrow = c(2, 3), mar = c(0, 0, 3.2, 0), xaxs = "i", yaxs = "i")
for (panel in panels) {
  helper_env$plot_layout_panel(
    panel$coords,
    panel$edges,
    title_text = panel$title,
    subtitle_text = panel$subtitle
  )
}
grDevices::dev.off()

utils::write.csv(trace_df, csv_path, row.names = FALSE)

writeLines(
  c(
    "# Carpet Level-4 Skew-Source Localization",
    "",
    "## Setup",
    "",
    sprintf("- graph: `Sierpinski carpet level %d`", level),
    sprintf("- seed: `%d`", seed),
    "- method: trace the current `grip.layout()` pipeline round-by-round and score every frame against the canonical carpet.",
    "",
    "## Key localization points",
    "",
    sprintf(
      "- largest active-subgraph symmetry drop: frame `%d -> %d` (`%s -> %s`), delta `%s`",
      trace_df$frame[[largest_sym_drop_idx - 1L]],
      trace_df$frame[[largest_sym_drop_idx]],
      trace_df$phase[[largest_sym_drop_idx - 1L]],
      trace_df$phase[[largest_sym_drop_idx]],
      format_num(delta_sym[[largest_sym_drop_idx]], 4L)
    ),
    sprintf(
      "- first fully active frame: frame `%d`, central-hole skew `%s`, RMSE `%s`",
      trace_df$frame[[first_full_idx]],
      format_num(trace_df$central.hole.skew[[first_full_idx]], 4L),
      format_num(trace_df$procrustes.rmse.active[[first_full_idx]], 4L)
    ),
    sprintf(
      "- first finest-level round: frame `%d`, central-hole skew `%s`, RMSE `%s`",
      trace_df$frame[[first_final_round_idx]],
      format_num(trace_df$central.hole.skew[[first_final_round_idx]], 4L),
      format_num(trace_df$procrustes.rmse.active[[first_final_round_idx]], 4L)
    ),
    sprintf(
      "- final frame: frame `%d`, central-hole skew `%s`, RMSE `%s`",
      trace_df$frame[[final_idx]],
      format_num(trace_df$central.hole.skew[[final_idx]], 4L),
      format_num(trace_df$procrustes.rmse.active[[final_idx]], 4L)
    ),
    sprintf(
      "- best full-graph trace frame by central-hole skew: frame `%d`, skew `%s`, RMSE `%s`",
      trace_df$frame[[best_full_idx]],
      format_num(trace_df$central.hole.skew[[best_full_idx]], 4L),
      format_num(trace_df$procrustes.rmse.active[[best_full_idx]], 4L)
    ),
    sprintf(
      "- igraph::KK reference: central-hole skew `%s`, RMSE `%s`",
      format_num(kk_geom$central.hole.skew[[1L]], 4L),
      format_num(kk_geom$procrustes.rmse[[1L]], 4L)
    ),
    "",
    "## Outputs",
    "",
    sprintf("- trace plot: `%s`", trace_plot_path),
    sprintf("- contact sheet: `%s`", sheet_path),
    sprintf("- trace CSV: `%s`", csv_path)
  ),
  con = md_path
)

message("Wrote:")
message("  ", trace_plot_path)
message("  ", sheet_path)
message("  ", csv_path)
message("  ", md_path)
