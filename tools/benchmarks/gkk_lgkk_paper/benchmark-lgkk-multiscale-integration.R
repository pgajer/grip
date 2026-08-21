#!/usr/bin/env Rscript

run_tag <- "lgkk-multiscale-linesearch-2026-03-30"
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

helper_env <- new.env(parent = baseenv())
sys.source(file.path("tools", "benchmarks", "gkk_lgkk_paper", "benchmark-sierpinski-baseline.R"), envir = helper_env)

fmt_num <- function(x, digits = 4L) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
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

plot_layout_panel <- function(spec,
                              coords,
                              title,
                              subtitle,
                              point_cex = 0.18,
                              edge_col = "#9a9a9a55",
                              point_col = "#1f4e79") {
  fit <- helper_env$align_to_target(coords, spec$canonical)
  aligned <- fit$aligned
  xlim <- range(spec$canonical[, 1L], aligned[, 1L], finite = TRUE)
  ylim <- range(spec$canonical[, 2L], aligned[, 2L], finite = TRUE)
  pad.x <- diff(xlim) * 0.05
  pad.y <- diff(ylim) * 0.05
  plot(aligned[, 1L], aligned[, 2L],
       asp = 1,
       axes = FALSE,
       xlab = "",
       ylab = "",
       pch = 16,
       cex = point_cex,
       col = point_col,
       main = title,
       sub = subtitle,
       xlim = xlim + c(-pad.x, pad.x),
       ylim = ylim + c(-pad.y, pad.y))
  segments(aligned[spec$edges[, 1L], 1L],
           aligned[spec$edges[, 1L], 2L],
           aligned[spec$edges[, 2L], 1L],
           aligned[spec$edges[, 2L], 2L],
           col = edge_col)
  points(aligned[, 1L], aligned[, 2L],
         pch = 16,
         cex = point_cex,
         col = point_col)
}

evaluate_layout <- function(spec, label, args, seed = 1L) {
  full_args <- c(
    list(edges = spec$edges, n = nrow(spec$canonical), dim = 2, seed = seed),
    if (!is.null(spec$edge_weights)) list(edge_weights = spec$edge_weights) else list(),
    args
  )
  timed <- system.time(coords <- do.call(grip, full_args))
  diag <- geometry.diagnostics(
    coords = coords,
    target.coords = spec$canonical,
    edges = spec$edges,
    family = spec$family,
    sample.size.symmetry = 512L,
    sample.size.wedges = 2048L,
    rng.seed = seed
  )
  lgkk <- score.landmark.geodesic.kk(
    coords = coords,
    edges = spec$edges,
    edge_weights = spec$edge_weights,
    n = nrow(spec$canonical),
    local_nbrs = 6L,
    landmark_count = 8L
  )
  list(
    coords = coords,
    metrics = data.frame(
      graph = spec$id,
      label = label,
      elapsed_sec = unname(timed[["elapsed"]]),
      rmse = diag$procrustes.rmse[[1L]],
      edge_axis_deviation = diag$edge.axis.deviation[[1L]],
      central_hole_skew = diag$central.hole.skew[[1L]],
      lgkk_rel_rmse = lgkk$lgkk.weighted.rel.rmse[[1L]],
      stringsAsFactors = FALSE
    )
  )
}

carpet <- carpet_spec(4L)
weighted_mesh <- weighted_mesh_spec(6L, 6L, horizontal = 1, vertical = 2)

carpet_configs <- list(
  baseline = list(label = "Baseline", args = list()),
  compiled4 = list(
    label = "Compiled LGKK x4",
    args = list(
      lgkk_multiscale_rounds = 4L,
      lgkk_local_nbrs = 6L,
      lgkk_landmark_count = 8L,
      lgkk_multiscale_scope = "all",
      lgkk_active_limit = 4096L
    )
  ),
  compiled4_balanced = list(
    label = "Compiled LGKK x4 + balanced anchors",
    args = list(
      lgkk_multiscale_rounds = 4L,
      lgkk_local_nbrs = 6L,
      lgkk_landmark_count = 8L,
      lgkk_multiscale_scope = "all",
      lgkk_active_limit = 4096L,
      insertion_anchor_strategy = "balanced_band",
      insertion_anchor_count = 8L,
      insertion_anchor_scope = "prev_misf",
      level0_insertion_mode = "least_squares",
      level0_local_kk_steps = 1L
    )
  ),
  post4 = list(
    label = "R LGKK polish x4",
    args = list(
      lgkk_polish_rounds = 4L,
      lgkk_local_nbrs = 6L,
      lgkk_landmark_count = 8L
    )
  )
)

mesh_configs <- list(
  baseline = list(label = "Baseline", args = list()),
  compiled4 = list(
    label = "Compiled LGKK x4",
    args = list(
      lgkk_multiscale_rounds = 4L,
      lgkk_local_nbrs = 6L,
      lgkk_landmark_count = 8L,
      lgkk_multiscale_scope = "all",
      lgkk_active_limit = 512L
    )
  ),
  post4 = list(
    label = "R LGKK polish x4",
    args = list(
      lgkk_polish_rounds = 4L,
      lgkk_local_nbrs = 6L,
      lgkk_landmark_count = 8L
    )
  )
)

carpet_results <- lapply(carpet_configs, function(cfg) {
  evaluate_layout(carpet, cfg$label, cfg$args, seed = 1L)
})
mesh_results <- lapply(mesh_configs, function(cfg) {
  evaluate_layout(weighted_mesh, cfg$label, cfg$args, seed = 1L)
})

trace_compare <- function(args) {
  tr <- do.call(
    trace.grip,
    c(
      list(edges = carpet$edges, n = nrow(carpet$canonical), dim = 2, seed = 1L,
           trace = "round", trace.every = 1L),
      args
    )
  )
  first_idx <- which(tr$meta$active_vertices == nrow(carpet$canonical))[1L]
  first_diag <- geometry.diagnostics(
    coords = tr$frames[[first_idx]],
    target.coords = carpet$canonical,
    edges = carpet$edges,
    family = carpet$family,
    sample.size.symmetry = 512L,
    sample.size.wedges = 2048L,
    rng.seed = 1L
  )
  final_diag <- geometry.diagnostics(
    coords = tr$final,
    target.coords = carpet$canonical,
    edges = carpet$edges,
    family = carpet$family,
    sample.size.symmetry = 512L,
    sample.size.wedges = 2048L,
    rng.seed = 1L
  )
  data.frame(
    config = if (length(args) == 0L) "baseline" else "balanced_prev_ls",
    first_full_rmse = first_diag$procrustes.rmse[[1L]],
    first_full_skew = first_diag$central.hole.skew[[1L]],
    final_rmse = final_diag$procrustes.rmse[[1L]],
    final_skew = final_diag$central.hole.skew[[1L]],
    stringsAsFactors = FALSE
  )
}

trace_metrics <- rbind(
  trace_compare(list()),
  trace_compare(list(
    insertion_anchor_strategy = "balanced_band",
    insertion_anchor_count = 8L,
    insertion_anchor_scope = "prev_misf",
    level0_insertion_mode = "least_squares",
    level0_local_kk_steps = 1L
  ))
)

metrics <- do.call(
  rbind,
  c(
    lapply(carpet_results, `[[`, "metrics"),
    lapply(mesh_results, `[[`, "metrics")
  )
)

write.csv(metrics,
          file.path(tmp_dir, "lgkk-multiscale-integration-metrics.csv"),
          row.names = FALSE)
write.csv(trace_metrics,
          file.path(tmp_dir, "carpet-level4-insertion-trace-metrics.csv"),
          row.names = FALSE)

png(file.path(tmp_dir, "carpet-level4-lgkk-multiscale-contact-sheet.png"),
    width = 1800, height = 1600, res = 180)
par(mfrow = c(2, 3), mar = c(1.2, 1.2, 3.2, 1.2), oma = c(0, 0, 2, 0))
plot_layout_panel(
  carpet,
  carpet$canonical,
  "Canonical",
  "reference"
)
for (id in names(carpet_results)) {
  row <- carpet_results[[id]]$metrics[1L, ]
  plot_layout_panel(
    carpet,
    carpet_results[[id]]$coords,
    row$label[[1L]],
    sprintf("RMSE %s | skew %s | %.1fs",
            fmt_num(row$rmse[[1L]]),
      fmt_num(row$central_hole_skew[[1L]]),
      row$elapsed_sec[[1L]])
  )
}
plot.new()
mtext("Level-4 carpet: compiled LGKK line-search stage and symmetry-aware insertion",
      outer = TRUE, line = 0.5, cex = 1.2)
dev.off()

png(file.path(tmp_dir, "weighted-mesh-6x6-lgkk-multiscale-contact-sheet.png"),
    width = 1600, height = 1600, res = 180)
par(mfrow = c(2, 2), mar = c(1.2, 1.2, 3.2, 1.2), oma = c(0, 0, 2, 0))
plot_layout_panel(
  weighted_mesh,
  weighted_mesh$canonical,
  "Canonical",
  "reference"
)
for (id in names(mesh_results)) {
  row <- mesh_results[[id]]$metrics[1L, ]
  plot_layout_panel(
    weighted_mesh,
    mesh_results[[id]]$coords,
    row$label[[1L]],
    sprintf("RMSE %s | LGKK %s | %.3fs",
            fmt_num(row$rmse[[1L]]),
            fmt_num(row$lgkk_rel_rmse[[1L]]),
            row$elapsed_sec[[1L]])
  )
}
mtext("Weighted mesh 6x6: compiled LGKK line-search stage vs R post-polish",
      outer = TRUE, line = 0.5, cex = 1.2)
dev.off()

summary_path <- file.path(tmp_dir, "lgkk-multiscale-integration-summary.md")
lines <- c(
  "# LGKK Multiscale Integration Summary",
  "",
  "## Final Metrics",
  "",
  "| Graph | Layout | Elapsed sec | RMSE | Axis deviation | Central-hole skew | LGKK rel. RMSE |",
  "| --- | --- | ---: | ---: | ---: | ---: | ---: |"
)
for (i in seq_len(nrow(metrics))) {
  row <- metrics[i, , drop = FALSE]
  lines <- c(
    lines,
    sprintf(
      "| %s | %s | %s | %s | %s | %s | %s |",
      row$graph[[1L]],
      row$label[[1L]],
      fmt_num(row$elapsed_sec[[1L]], digits = 3L),
      fmt_num(row$rmse[[1L]]),
      fmt_num(row$edge_axis_deviation[[1L]]),
      fmt_num(row$central_hole_skew[[1L]]),
      fmt_num(row$lgkk_rel_rmse[[1L]])
    )
  )
}
lines <- c(
  lines,
  "",
  "## Carpet Insertion Trace Check",
  "",
  "| Config | First full RMSE | First full skew | Final RMSE | Final skew |",
  "| --- | ---: | ---: | ---: | ---: |"
)
for (i in seq_len(nrow(trace_metrics))) {
  row <- trace_metrics[i, , drop = FALSE]
  lines <- c(
    lines,
    sprintf(
      "| %s | %s | %s | %s | %s |",
      row$config[[1L]],
      fmt_num(row$first_full_rmse[[1L]]),
      fmt_num(row$first_full_skew[[1L]]),
      fmt_num(row$final_rmse[[1L]]),
      fmt_num(row$final_skew[[1L]])
    )
  )
}
lines <- c(
  lines,
  "",
  "## Figures",
  "",
  sprintf("- Carpet contact sheet: `%s`", file.path(tmp_dir, "carpet-level4-lgkk-multiscale-contact-sheet.png")),
  sprintf("- Weighted mesh contact sheet: `%s`", file.path(tmp_dir, "weighted-mesh-6x6-lgkk-multiscale-contact-sheet.png"))
)
writeLines(lines, summary_path)

message("Wrote benchmark outputs to ", tmp_dir)
