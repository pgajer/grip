#!/usr/bin/env Rscript

run_tag <- "lgkk-optimizer-suite-2026-03-30"
manual_root <- file.path("dev", "manual")
tmp_dir <- file.path(manual_root, "tmp", run_tag)
pdf_dir <- file.path(manual_root, "pdf", run_tag)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'devtools' or the 'grip' package to run the LGKK optimizer suite.")
}

if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("The 'igraph' package is required for the KK comparison in this benchmark.")
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
    if (r1 == r2 && abs(c1 - c2) == 1L) {
      horizontal
    } else {
      vertical
    }
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

unweighted_mesh_spec <- function(h = 8L, w = 8L) {
  edges <- edges.mesh(h, w)
  coords <- matrix(0, nrow = h * w, ncol = 2L)
  for (i in seq_len(h)) {
    for (j in seq_len(w)) {
      id <- (i - 1L) * w + j
      coords[id, ] <- c(j, h - i)
    }
  }
  list(
    id = sprintf("mesh_%dx%d", h, w),
    label = sprintf("Mesh %dx%d", h, w),
    family = "mesh",
    edges = edges,
    edge_weights = NULL,
    canonical = coords
  )
}

carpet_spec <- function(level) {
  built <- helper_env$build_sierpinski_carpet(level)
  list(
    id = sprintf("carpet_l%d", level),
    label = sprintf("Sierpinski carpet level %d", level),
    family = "sierpinski.carpet",
    edges = built$edges,
    edge_weights = NULL,
    canonical = built$coords
  )
}

triangle_spec <- function(level) {
  built <- helper_env$build_sierpinski_triangle(level)
  list(
    id = sprintf("triangle_l%d", level),
    label = sprintf("Sierpinski triangle level %d", level),
    family = "sierpinski.triangle",
    edges = built$edges,
    edge_weights = NULL,
    canonical = built$coords
  )
}

stage1_specs <- list(
  carpet_spec(3L),
  weighted_mesh_spec(6L, 6L, horizontal = 1, vertical = 2)
)

stage2_specs <- list(
  carpet_spec(3L),
  carpet_spec(4L),
  triangle_spec(4L),
  unweighted_mesh_spec(8L, 8L),
  weighted_mesh_spec(6L, 6L, horizontal = 1, vertical = 2)
)

candidate_grid <- expand.grid(
  lgkk_polish_rounds = c(2L, 4L, 8L),
  lgkk_local_nbrs = c(6L, 12L, 20L),
  lgkk_landmark_count = c(4L, 8L),
  stringsAsFactors = FALSE
)
candidate_grid$id <- sprintf(
  "lgkk_f%02d_q%02d_m%02d",
  candidate_grid$lgkk_polish_rounds,
  candidate_grid$lgkk_local_nbrs,
  candidate_grid$lgkk_landmark_count
)

baseline_method <- list(
  id = "baseline",
  label = "Baseline GRIP",
  kind = "grip",
  params = list()
)

igraph_method <- list(
  id = "igraph_kk",
  label = "igraph KK",
  kind = "igraph_kk",
  params = list()
)

lgkk_methods <- lapply(seq_len(nrow(candidate_grid)), function(i) {
  row <- candidate_grid[i, , drop = FALSE]
  list(
    id = row$id[[1L]],
    label = sprintf("LGKK f=%d q=%d M=%d",
                    row$lgkk_polish_rounds[[1L]],
                    row$lgkk_local_nbrs[[1L]],
                    row$lgkk_landmark_count[[1L]]),
    kind = "grip",
    params = list(
      lgkk_polish_rounds = row$lgkk_polish_rounds[[1L]],
      lgkk_local_nbrs = row$lgkk_local_nbrs[[1L]],
      lgkk_landmark_count = row$lgkk_landmark_count[[1L]]
    )
  )
})

stage1_methods <- c(list(baseline_method), lgkk_methods, list(igraph_method))

run_method <- function(spec, method, seed = 1L) {
  n <- nrow(spec$canonical)
  started <- proc.time()[["elapsed"]]
  coords <- switch(
    method$kind,
    grip = do.call(
      grip.layout,
      c(
        list(
          edges = spec$edges,
          n = n,
          dim = 2,
          edge_weights = spec$edge_weights,
          seed = seed
        ),
        method$params
      )
    ),
    igraph_kk = {
      graph <- igraph::graph_from_edgelist(spec$edges, directed = FALSE)
      if (!is.null(spec$edge_weights)) {
        igraph::E(graph)$weight <- spec$edge_weights
      }
      igraph::layout_with_kk(
        graph,
        weights = if (!is.null(spec$edge_weights)) igraph::E(graph)$weight else NULL
      )
    },
    stop(sprintf("Unsupported method kind: %s", method$kind))
  )
  elapsed <- proc.time()[["elapsed"]] - started
  list(coords = as.matrix(coords), elapsed = as.double(elapsed))
}

score_method <- function(spec, method, result) {
  n <- nrow(spec$canonical)
  lgkk_local_eval <- min(12L, n - 1L)
  lgkk_landmark_eval <- min(6L, n - 1L)
  geo <- grip.geometry.diagnostics(
    coords = result$coords,
    target.coords = spec$canonical,
    edges = spec$edges,
    family = spec$family,
    sample.size.symmetry = 256L,
    sample.size.wedges = 1000L,
    rng.seed = 1L
  )
  lgkk <- grip.score.landmark.geodesic.kk(
    coords = result$coords,
    edges = spec$edges,
    n = n,
    edge_weights = spec$edge_weights,
    local_nbrs = lgkk_local_eval,
    landmark_count = lgkk_landmark_eval
  )
  edge_stats <- helper_env$edge_length_stats(result$coords, spec$edges)
  data.frame(
    graph_id = spec$id,
    graph_label = spec$label,
    family = spec$family,
    method_id = method$id,
    method_label = method$label,
    method_kind = method$kind,
    lgkk_polish_rounds = if (!is.null(method$params$lgkk_polish_rounds)) method$params$lgkk_polish_rounds else 0L,
    lgkk_local_nbrs = if (!is.null(method$params$lgkk_local_nbrs)) method$params$lgkk_local_nbrs else 0L,
    lgkk_landmark_count = if (!is.null(method$params$lgkk_landmark_count)) method$params$lgkk_landmark_count else 0L,
    elapsed_sec = result$elapsed,
    procrustes_rmse = geo$procrustes.rmse[[1L]],
    edge_axis_deviation = geo$edge.axis.deviation[[1L]],
    edge_length_cv = edge_stats$cv,
    lgkk_weighted_rel_rmse = lgkk$lgkk.weighted.rel.rmse[[1L]],
    central_hole_skew = geo$central.hole.skew[[1L]],
    central_hole_aspect_error = geo$central.hole.aspect.error[[1L]],
    stringsAsFactors = FALSE
  )
}

run_suite <- function(specs, methods) {
  rows <- list()
  idx <- 0L
  layouts <- list()
  for (spec in specs) {
    layouts[[spec$id]] <- list()
    for (method in methods) {
      idx <- idx + 1L
      result <- run_method(spec, method, seed = 1L)
      layouts[[spec$id]][[method$id]] <- result$coords
      rows[[idx]] <- score_method(spec, method, result)
    }
  }
  list(
    metrics = do.call(rbind, rows),
    layouts = layouts
  )
}

stage1 <- run_suite(stage1_specs, stage1_methods)
stage1_metrics <- stage1$metrics
write.csv(stage1_metrics,
          file.path(tmp_dir, "stage1-raw-metrics.csv"),
          row.names = FALSE)

stage1_lgkk <- stage1_metrics[grepl("^lgkk_", stage1_metrics$method_id), , drop = FALSE]
stage1_summary <- aggregate(
  cbind(elapsed_sec,
        procrustes_rmse,
        edge_axis_deviation,
        edge_length_cv,
        lgkk_weighted_rel_rmse,
        central_hole_skew) ~ method_id + method_label +
    lgkk_polish_rounds + lgkk_local_nbrs + lgkk_landmark_count,
  data = stage1_lgkk,
  FUN = function(x) mean(x, na.rm = TRUE)
)
stage1_summary$quality_score <-
  rank(stage1_summary$procrustes_rmse, ties.method = "first") +
  rank(stage1_summary$lgkk_weighted_rel_rmse, ties.method = "first") +
  rank(stage1_summary$edge_axis_deviation, ties.method = "first")
stage1_summary$value_score <-
  stage1_summary$quality_score +
  rank(stage1_summary$elapsed_sec, ties.method = "first")
stage1_summary <- stage1_summary[order(stage1_summary$quality_score,
                                       stage1_summary$procrustes_rmse,
                                       stage1_summary$elapsed_sec), , drop = FALSE]
write.csv(stage1_summary,
          file.path(tmp_dir, "stage1-candidate-summary.csv"),
          row.names = FALSE)

best_quality <- stage1_summary[1L, , drop = FALSE]
best_value <- stage1_summary[order(stage1_summary$value_score,
                                   stage1_summary$procrustes_rmse,
                                   stage1_summary$elapsed_sec), , drop = FALSE]
if (best_value$method_id[[1L]] == best_quality$method_id[[1L]] && nrow(best_value) > 1L) {
  best_value <- best_value[2L, , drop = FALSE]
} else {
  best_value <- best_value[1L, , drop = FALSE]
}

make_method_from_row <- function(row, label_prefix) {
  list(
    id = row$method_id[[1L]],
    label = sprintf(
      "%s (f=%d, q=%d, M=%d)",
      label_prefix,
      row$lgkk_polish_rounds[[1L]],
      row$lgkk_local_nbrs[[1L]],
      row$lgkk_landmark_count[[1L]]
    ),
    kind = "grip",
    params = list(
      lgkk_polish_rounds = row$lgkk_polish_rounds[[1L]],
      lgkk_local_nbrs = row$lgkk_local_nbrs[[1L]],
      lgkk_landmark_count = row$lgkk_landmark_count[[1L]]
    )
  )
}

best_quality_method <- make_method_from_row(best_quality, "Best quality LGKK")
best_value_method <- make_method_from_row(best_value, "Best value LGKK")

stage2_methods <- list(
  baseline_method,
  best_quality_method,
  best_value_method,
  igraph_method
)
stage2 <- run_suite(stage2_specs, stage2_methods)
stage2_metrics <- stage2$metrics
write.csv(stage2_metrics,
          file.path(tmp_dir, "stage2-raw-metrics.csv"),
          row.names = FALSE)

stage2_summary <- stage2_metrics
stage2_summary$method_order <- match(
  stage2_summary$method_id,
  c("baseline", best_quality_method$id, best_value_method$id, "igraph_kk")
)
stage2_summary <- stage2_summary[order(stage2_summary$graph_id, stage2_summary$method_order), , drop = FALSE]
write.csv(stage2_summary,
          file.path(tmp_dir, "stage2-summary.csv"),
          row.names = FALSE)

plot_heatmap_metric <- function(df, metric_col, main_prefix, out_path) {
  png(out_path, width = 1600, height = 900, res = 150)
  old <- par(no.readonly = TRUE)
  on.exit({
    par(old)
    dev.off()
  }, add = TRUE)
  par(mfrow = c(2, 2), mar = c(4, 5, 4, 2))
  ms <- sort(unique(df$lgkk_landmark_count))
  for (m in ms) {
    sub <- df[df$lgkk_landmark_count == m, , drop = FALSE]
    vals <- xtabs(sub[[metric_col]] ~ sub$lgkk_local_nbrs + sub$lgkk_polish_rounds)
    image(
      x = seq_len(ncol(vals)),
      y = seq_len(nrow(vals)),
      z = t(vals[nrow(vals):1L, , drop = FALSE]),
      axes = FALSE,
      col = hcl.colors(20, "YlOrRd", rev = TRUE),
      main = sprintf("%s, M = %d", main_prefix, m),
      xlab = "LGKK polish rounds",
      ylab = "LGKK local neighbors"
    )
    axis(1, at = seq_len(ncol(vals)), labels = colnames(vals))
    axis(2, at = seq_len(nrow(vals)), labels = rev(rownames(vals)))
    for (i in seq_len(nrow(vals))) {
      for (j in seq_len(ncol(vals))) {
        text(j, nrow(vals) - i + 1L, labels = fmt_num(vals[i, j], 2L), cex = 0.9)
      }
    }
  }
}

plot_frontier <- function(df, out_path) {
  png(out_path, width = 1400, height = 900, res = 150)
  old <- par(no.readonly = TRUE)
  on.exit({
    par(old)
    dev.off()
  }, add = TRUE)
  par(mar = c(5, 5, 4, 2))
  cols <- c("steelblue", "tomato", "darkgreen")
  rounds_levels <- sort(unique(df$lgkk_polish_rounds))
  round_col <- cols[match(df$lgkk_polish_rounds, rounds_levels)]
  plot(df$elapsed_sec,
       df$procrustes_rmse,
       pch = 19,
       col = round_col,
       xlab = "Mean elapsed time (sec)",
       ylab = "Mean Procrustes RMSE",
       main = "LGKK optimizer frontier on the Stage-1 training set")
  text(df$elapsed_sec,
       df$procrustes_rmse,
       labels = sprintf("q%d/M%d", df$lgkk_local_nbrs, df$lgkk_landmark_count),
       pos = 4,
       cex = 0.7)
  legend("topright",
         legend = sprintf("rounds = %d", rounds_levels),
         col = cols[seq_along(rounds_levels)],
         pch = 19,
         bty = "n")
}

plot_heatmap_metric(stage1_summary, "quality_score", "Quality score", file.path(tmp_dir, "stage1-quality-heatmaps.png"))
plot_heatmap_metric(stage1_summary, "value_score", "Value score", file.path(tmp_dir, "stage1-value-heatmaps.png"))
plot_frontier(stage1_summary, file.path(tmp_dir, "stage1-frontier.png"))

render_contact_sheet <- function(spec, layouts, out_path) {
  methods <- list(
    list(id = "canonical", label = "Canonical", coords = spec$canonical, elapsed = NA_real_, rmse = 0, lgkk = 0),
    list(id = "baseline", label = "Baseline GRIP"),
    list(id = best_quality_method$id, label = "Best quality LGKK"),
    list(id = best_value_method$id, label = "Best value LGKK"),
    list(id = "igraph_kk", label = "igraph KK")
  )
  metrics_sub <- stage2_metrics[stage2_metrics$graph_id == spec$id, , drop = FALSE]
  panels <- lapply(methods, function(m) {
    if (m$id == "canonical") {
      list(label = m$label, coords = spec$canonical, subtitle = "reference")
    } else {
      row <- metrics_sub[metrics_sub$method_id == m$id, , drop = FALSE]
      fit <- helper_env$align_to_target(layouts[[m$id]], spec$canonical)
      list(
        label = m$label,
        coords = fit$aligned,
        subtitle = sprintf(
          "RMSE %s | relLGKK %s | time %ss",
          fmt_num(row$procrustes_rmse[[1L]], 4L),
          fmt_num(row$lgkk_weighted_rel_rmse[[1L]], 4L),
          fmt_num(row$elapsed_sec[[1L]], 2L)
        )
      )
    }
  })

  png(out_path, width = 2200, height = 550, res = 150)
  old <- par(no.readonly = TRUE)
  on.exit({
    par(old)
    dev.off()
  }, add = TRUE)
  layout(matrix(seq_len(5), nrow = 1L))
  for (panel in panels) {
    helper_env$plot_layout_panel(
      panel$coords,
      spec$edges,
      title_text = panel$label,
      subtitle_text = panel$subtitle
    )
  }
}

for (spec in stage2_specs) {
  render_contact_sheet(
    spec = spec,
    layouts = stage2$layouts[[spec$id]],
    out_path = file.path(tmp_dir, sprintf("%s-contact-sheet.png", spec$id))
  )
}

stage1_lines <- c(
  "# LGKK Optimizer Stage 1 Summary",
  "",
  sprintf("- Training graphs: %s", paste(vapply(stage1_specs, `[[`, character(1L), "label"), collapse = ", ")),
  sprintf("- Best quality candidate: `%s`", best_quality$method_id[[1L]]),
  sprintf("- Best value candidate: `%s`", best_value$method_id[[1L]]),
  "",
  "Top Stage-1 candidates:",
  ""
)
top_stage1 <- head(stage1_summary[, c("method_id", "procrustes_rmse", "lgkk_weighted_rel_rmse", "edge_axis_deviation", "elapsed_sec", "quality_score", "value_score")], 8L)
stage1_lines <- c(stage1_lines, apply(top_stage1, 1L, function(row) {
  sprintf(
    "- `%s`: mean RMSE `%s`, relLGKK `%s`, axis dev `%s`, time `%ss`, quality `%s`, value `%s`",
    row[["method_id"]],
    fmt_num(as.numeric(row[["procrustes_rmse"]]), 4L),
    fmt_num(as.numeric(row[["lgkk_weighted_rel_rmse"]]), 4L),
    fmt_num(as.numeric(row[["edge_axis_deviation"]]), 4L),
    fmt_num(as.numeric(row[["elapsed_sec"]]), 2L),
    fmt_num(as.numeric(row[["quality_score"]]), 1L),
    fmt_num(as.numeric(row[["value_score"]]), 1L)
  )
}))
writeLines(stage1_lines, con = file.path(tmp_dir, "stage1-summary.md"))

stage2_lines <- c(
  "# LGKK Optimizer Stage 2 Summary",
  "",
  sprintf("- Confirmation graphs: %s", paste(vapply(stage2_specs, `[[`, character(1L), "label"), collapse = ", ")),
  sprintf("- Best quality candidate carried forward: `%s`", best_quality$method_id[[1L]]),
  sprintf("- Best value candidate carried forward: `%s`", best_value$method_id[[1L]]),
  ""
)
for (spec in stage2_specs) {
  stage2_lines <- c(stage2_lines, sprintf("## %s", spec$label), "")
  sub <- stage2_summary[stage2_summary$graph_id == spec$id, , drop = FALSE]
  for (i in seq_len(nrow(sub))) {
    row <- sub[i, , drop = FALSE]
    stage2_lines <- c(stage2_lines, sprintf(
      "- `%s`: RMSE `%s`, relLGKK `%s`, axis dev `%s`, edge CV `%s`, skew `%s`, time `%ss`",
      row$method_id[[1L]],
      fmt_num(row$procrustes_rmse[[1L]], 4L),
      fmt_num(row$lgkk_weighted_rel_rmse[[1L]], 4L),
      fmt_num(row$edge_axis_deviation[[1L]], 4L),
      fmt_num(row$edge_length_cv[[1L]], 4L),
      fmt_num(row$central_hole_skew[[1L]], 4L),
      fmt_num(row$elapsed_sec[[1L]], 2L)
    ))
  }
  stage2_lines <- c(stage2_lines, "")
}
writeLines(stage2_lines, con = file.path(tmp_dir, "stage2-summary.md"))
