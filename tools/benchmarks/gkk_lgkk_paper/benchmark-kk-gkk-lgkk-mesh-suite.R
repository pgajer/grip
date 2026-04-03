#!/usr/bin/env Rscript

run_tag <- "kk-gkk-lgkk-mesh-suite-2026-03-31"
repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
output_root <- file.path(repo_root, "output", "benchmarks", run_tag)
figure_dir <- file.path(output_root, "figures")
table_dir <- file.path(output_root, "tables")
data_dir <- file.path(output_root, "data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'pkgload', 'devtools', or the 'grip' package to run this benchmark.")
}

if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Package 'igraph' is required to run this benchmark.")
}

grip_ns <- asNamespace("grip")
pkg_internal <- function(name) get(name, envir = grip_ns)

fmt_num <- function(x, digits = 4L) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}

normalize_coords <- function(coords) {
  coords <- as.matrix(coords)
  centered <- sweep(coords, 2L, colMeans(coords), "-", check.margin = FALSE)
  radius <- max(sqrt(rowSums(centered^2)))
  if (!is.finite(radius) || radius <= 0) {
    return(centered)
  }
  centered / radius
}

mean_sd_summary <- function(x) {
  ok <- is.finite(x)
  if (!any(ok)) {
    return(c(mean = NA_real_, sd = NA_real_))
  }
  vals <- x[ok]
  c(
    mean = mean(vals),
    sd = if (length(vals) > 1L) stats::sd(vals) else 0
  )
}

summarize_metrics <- function(df, group_cols, metric_cols) {
  split_key <- interaction(df[group_cols], drop = TRUE, lex.order = TRUE)
  pieces <- split(df, split_key)
  rows <- lapply(pieces, function(piece) {
    row <- piece[1L, group_cols, drop = FALSE]
    for (metric in metric_cols) {
      stats <- mean_sd_summary(piece[[metric]])
      row[[paste0(metric, "_mean")]] <- stats[["mean"]]
      row[[paste0(metric, "_sd")]] <- stats[["sd"]]
    }
    row
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

surface_configs <- list(
  list(
    id = "saddle_iso",
    label = "Saddle (isotropic)",
    surface = "saddle",
    amplitude = 0.75,
    x_scale = 1.00,
    y_scale = 1.00,
    connectivity = "orthogonal"
  ),
  list(
    id = "paraboloid_aniso",
    label = "Paraboloid (anisotropic)",
    surface = "paraboloid",
    amplitude = 0.75,
    x_scale = 1.00,
    y_scale = 1.25,
    connectivity = "orthogonal"
  ),
  list(
    id = "ripple",
    label = "Ripple",
    surface = "ripple",
    amplitude = 0.60,
    freq_u = 1.5,
    freq_v = 1.0,
    x_scale = 1.00,
    y_scale = 1.00,
    connectivity = "orthogonal"
  )
)

benchmark_cfg <- list(
  sizes = c(8L, 10L, 12L),
  dims = c(2L, 3L),
  seeds = 1:3,
  gkk_max_iter = 12L,
  lgkk_max_iter = 12L,
  lgkk_local_nbrs = 8L,
  lgkk_landmark_count = 10L,
  stress_sample = 2000L,
  jitter_xy = 0.08,
  jitter_z = 0.08
)

build_mesh_instance <- function(surface_cfg, side) {
  args <- c(
    list(
      h = side,
      w = side,
      normalize = "median",
      connectivity = surface_cfg$connectivity
    ),
    surface_cfg[setdiff(names(surface_cfg), c("id", "label", "connectivity"))]
  )
  spec <- do.call(mesh.surface.graph, args)
  spec$surface_id <- surface_cfg$id
  spec$surface_label <- surface_cfg$label
  spec$size <- side
  spec$size_label <- sprintf("%dx%d", side, side)
  spec$case_id <- sprintf("%s_%02dx%02d", surface_cfg$id, side, side)
  spec$case_label <- sprintf("%s %s", surface_cfg$label, spec$size_label)
  spec$graph <- igraph::graph_from_edgelist(as.matrix(spec$edges), directed = FALSE)
  spec$prepared_gkk <- grip.prepare.geodesic.kk(
    edges = spec$edges,
    n = spec$n,
    edge_weights = spec$edge_weights
  )
  spec$prepared_lgkk <- grip.prepare.landmark.geodesic.kk(
    edges = spec$edges,
    n = spec$n,
    edge_weights = spec$edge_weights,
    local_nbrs = benchmark_cfg$lgkk_local_nbrs,
    landmark_count = benchmark_cfg$lgkk_landmark_count
  )
  weight_stats <- mean_sd_summary(spec$edge_weights)
  spec$edge_weight_mean <- weight_stats[["mean"]]
  spec$edge_weight_sd <- weight_stats[["sd"]]
  spec$edge_weight_cv <- if (is.finite(spec$edge_weight_mean) && spec$edge_weight_mean > 0) {
    spec$edge_weight_sd / spec$edge_weight_mean
  } else {
    NA_real_
  }
  spec
}

build_initial_layout <- function(spec, dim, seed) {
  dim <- as.integer(dim)
  if (!(dim %in% c(2L, 3L))) {
    stop("dim must be 2 or 3")
  }
  base2d <- normalize_coords(spec$coords_param)
  draw_seed <- as.integer(
    100000L +
      1000L * seed +
      100L * dim +
      10L * spec$size +
      match(spec$surface_id, vapply(surface_configs, `[[`, "", "id"))
  )
  set.seed(draw_seed)
  if (dim == 2L) {
    coords <- base2d + matrix(
      stats::rnorm(spec$n * 2L, sd = benchmark_cfg$jitter_xy),
      ncol = 2L
    )
  } else {
    coords <- cbind(base2d, 0)
    jitter <- matrix(
      stats::rnorm(spec$n * 3L, sd = benchmark_cfg$jitter_z),
      ncol = 3L
    )
    jitter[, 1L:2L] <- 0.45 * jitter[, 1L:2L, drop = FALSE]
    coords <- coords + jitter
  }
  storage.mode(coords) <- "double"
  sweep(coords, 2L, colMeans(coords), "-", check.margin = FALSE)
}

score_classical_kk <- function(coords,
                               prepared,
                               stiffness = 1.0,
                               distance_floor = 1e-8) {
  g <- as.double(prepared$pair_graph_distance)
  pair_matrix <- prepared$pair_matrix
  if (length(g) == 0L) {
    return(list(
      scale.L0 = NA_real_,
      energy = NA_real_,
      weighted.rmse = NA_real_,
      weighted.rel.rmse = NA_real_,
      mean.abs.error = NA_real_,
      mean.rel.error = NA_real_
    ))
  }
  kk <- as.double(stiffness) / pmax(g, as.double(distance_floor))^2
  diffs <- coords[pair_matrix[, 1L], , drop = FALSE] -
    coords[pair_matrix[, 2L], , drop = FALSE]
  chord <- sqrt(rowSums(diffs^2))
  denom <- sum(kk * g * g)
  scale.L0 <- if (!is.finite(denom) || denom <= 0) {
    NA_real_
  } else {
    sum(kk * g * chord) / denom
  }
  if (!is.finite(scale.L0)) {
    return(list(
      scale.L0 = NA_real_,
      energy = NA_real_,
      weighted.rmse = NA_real_,
      weighted.rel.rmse = NA_real_,
      mean.abs.error = NA_real_,
      mean.rel.error = NA_real_
    ))
  }
  target <- scale.L0 * g
  resid <- chord - target
  rel.resid <- resid / pmax(target, distance_floor)
  weight_sum <- sum(kk)
  list(
    scale.L0 = scale.L0,
    energy = 0.5 * sum(kk * resid^2),
    weighted.rmse = sqrt(sum(kk * resid^2) / weight_sum),
    weighted.rel.rmse = sqrt(sum(kk * rel.resid^2) / weight_sum),
    mean.abs.error = mean(abs(resid)),
    mean.rel.error = mean(abs(rel.resid))
  )
}

score_layout <- function(spec, coords, dim, seed) {
  gkk <- grip.score.geodesic.kk(
    coords = coords,
    prepared = spec$prepared_gkk
  )
  lgkk <- grip.score.landmark.geodesic.kk(
    coords = coords,
    prepared = spec$prepared_lgkk
  )
  classic <- score_classical_kk(coords, prepared = spec$prepared_gkk)
  edge_stats <- pkg_internal("grip.edge.length.stats")(coords, spec$edges)
  sampled_stress <- pkg_internal("grip.sampled.stress")(
    coords = coords,
    adj.list = spec$prepared_gkk$adj_list,
    weight.list = spec$prepared_gkk$weight_list,
    sample.size = benchmark_cfg$stress_sample,
    rng.seed = 7000L + 100L * dim + 10L * spec$size + seed
  )
  target <- if (dim == 3L) spec$coords_surface else spec$coords_param
  target_space <- if (dim == 3L) "surface_3d" else "parameter_2d"
  procrustes_rmse <- pkg_internal("grip.align.to.target.nd")(coords, target)$rmse

  data.frame(
    case_id = spec$case_id,
    case_label = spec$case_label,
    surface_id = spec$surface_id,
    surface_label = spec$surface_label,
    surface = spec$surface,
    size = spec$size,
    size_label = spec$size_label,
    dim = dim,
    seed = seed,
    vertices = spec$n,
    edges = nrow(spec$edges),
    edge_weight_cv = spec$edge_weight_cv,
    target_space = target_space,
    kk_energy = classic$energy,
    kk_rel_rmse = classic$weighted.rel.rmse,
    kk_abs_error = classic$mean.abs.error,
    gkk_energy = gkk$gkk.energy[[1L]],
    gkk_rel_rmse = gkk$gkk.weighted.rel.rmse[[1L]],
    gkk_abs_error = gkk$gkk.mean.abs.path.error[[1L]],
    lgkk_sparse_energy = lgkk$lgkk.energy[[1L]],
    lgkk_sparse_rel_rmse = lgkk$lgkk.weighted.rel.rmse[[1L]],
    lgkk_sparse_abs_error = lgkk$lgkk.mean.abs.path.error[[1L]],
    procrustes_rmse = procrustes_rmse,
    edge_length_cv = edge_stats$cv,
    median_edge_length = edge_stats$median,
    sampled_stress = sampled_stress,
    stringsAsFactors = FALSE
  )
}

evaluate_layout <- function(spec, dim, seed) {
  initial <- build_initial_layout(spec, dim = dim, seed = seed)
  out <- list()

  row_init <- score_layout(spec, initial, dim = dim, seed = seed)
  row_init$method <- "start"
  row_init$method_label <- "Shared start"
  row_init$runtime_sec <- NA_real_
  row_init$accepted_steps <- NA_integer_
  out$start <- list(coords = initial, row = row_init)

  timed_kk <- system.time({
    coords_kk <- igraph::layout_with_kk(
      spec$graph,
      coords = initial,
      dim = dim,
      weights = spec$edge_weights
    )
  })
  row_kk <- score_layout(spec, coords_kk, dim = dim, seed = seed)
  row_kk$method <- "kk"
  row_kk$method_label <- "KK"
  row_kk$runtime_sec <- unname(timed_kk[["elapsed"]])
  row_kk$accepted_steps <- NA_integer_
  out$kk <- list(coords = coords_kk, row = row_kk)

  timed_gkk <- system.time({
    opt_gkk <- grip.optimize.geodesic.kk(
      coords = coords_kk,
      prepared = spec$prepared_gkk,
      max_iter = benchmark_cfg$gkk_max_iter,
      scale_mode = "profiled",
      return_trace = TRUE
    )
  })
  row_gkk <- score_layout(spec, opt_gkk$coords, dim = dim, seed = seed)
  row_gkk$method <- "gkk"
  row_gkk$method_label <- "KK->GKK"
  row_gkk$runtime_sec <- unname(timed_gkk[["elapsed"]])
  row_gkk$accepted_steps <- max(opt_gkk$trace$iteration[opt_gkk$trace$accepted], na.rm = TRUE)
  out$gkk <- list(coords = opt_gkk$coords, row = row_gkk, trace = opt_gkk$trace)

  timed_lgkk <- system.time({
    opt_lgkk <- grip.optimize.landmark.geodesic.kk(
      coords = coords_kk,
      prepared = spec$prepared_lgkk,
      max_iter = benchmark_cfg$lgkk_max_iter,
      return_trace = TRUE
    )
  })
  row_lgkk <- score_layout(spec, opt_lgkk$coords, dim = dim, seed = seed)
  row_lgkk$method <- "lgkk"
  row_lgkk$method_label <- "KK->LGKK"
  row_lgkk$runtime_sec <- unname(timed_lgkk[["elapsed"]])
  row_lgkk$accepted_steps <- max(opt_lgkk$trace$iteration[opt_lgkk$trace$accepted], na.rm = TRUE)
  out$lgkk <- list(coords = opt_lgkk$coords, row = row_lgkk, trace = opt_lgkk$trace)

  out
}

write_metric_panel_plot <- function(summary_df,
                                    dim_value,
                                    metric,
                                    ylab,
                                    path,
                                    log_scale = FALSE) {
  surface_ids <- vapply(surface_configs, `[[`, "", "id")
  surface_labels <- vapply(surface_configs, `[[`, "", "label")
  method_order <- c("start", "kk", "gkk", "lgkk")
  method_labels <- c(start = "Start", kk = "KK", gkk = "KK->GKK", lgkk = "KK->LGKK")
  method_cols <- c(start = "#bcb8b1", kk = "#355070", gkk = "#9c6644", lgkk = "#588157")

  panel_df <- summary_df[summary_df$dim == dim_value & summary_df$method %in% method_order, , drop = FALSE]
  if (nrow(panel_df) == 0L) {
    return(invisible(NULL))
  }

  upper <- panel_df[[paste0(metric, "_mean")]] + panel_df[[paste0(metric, "_sd")]]
  upper <- upper[is.finite(upper) & upper > 0]
  ymax <- if (length(upper) == 0L) 1 else max(upper)
  ymin <- if (isTRUE(log_scale)) {
    vals <- panel_df[[paste0(metric, "_mean")]]
    vals <- vals[is.finite(vals) & vals > 0]
    if (length(vals) == 0L) 1e-4 else min(vals) * 0.8
  } else {
    0
  }

  grDevices::png(path, width = 2400, height = 900, res = 180, bg = "#f7f3ea")
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(1, 3), mar = c(4.4, 4.8, 3.4, 1.0), oma = c(0, 0, 1.5, 0))

  for (i in seq_along(surface_ids)) {
    surf_df <- panel_df[panel_df$surface_id == surface_ids[[i]], , drop = FALSE]
    surf_df <- surf_df[order(surf_df$size, match(surf_df$method, method_order)), , drop = FALSE]
    graphics::plot(
      range(benchmark_cfg$sizes),
      c(ymin, if (isTRUE(log_scale)) ymax * 1.2 else ymax * 1.1),
      type = "n",
      log = if (isTRUE(log_scale)) "y" else "",
      xlab = "Mesh side length",
      ylab = ylab,
      xaxt = "n",
      main = surface_labels[[i]]
    )
    graphics::axis(1, at = benchmark_cfg$sizes, labels = benchmark_cfg$sizes)
    graphics::grid(col = "#dedbd2", lty = "dotted")
    for (method in method_order) {
      df_m <- surf_df[surf_df$method == method, , drop = FALSE]
      if (nrow(df_m) == 0L) {
        next
      }
      x <- df_m$size
      y <- df_m[[paste0(metric, "_mean")]]
      y_sd <- df_m[[paste0(metric, "_sd")]]
      graphics::lines(
        x,
        y,
        type = "b",
        lwd = if (identical(method, "start")) 2.0 else 2.5,
        lty = if (identical(method, "start")) "dashed" else "solid",
        pch = if (identical(method, "start")) 1 else 19,
        col = method_cols[[method]]
      )
      keep <- is.finite(y) & is.finite(y_sd)
      keep <- keep & is.finite(y_sd) & (y_sd > 0)
      if (any(keep)) {
        graphics::arrows(
          x0 = x[keep],
          y0 = pmax(ymin, y[keep] - y_sd[keep]),
          x1 = x[keep],
          y1 = y[keep] + y_sd[keep],
          angle = 90,
          code = 3,
          length = 0.04,
          col = method_cols[[method]]
        )
      }
    }
    if (i == 1L) {
      graphics::legend(
        "topleft",
        legend = unname(method_labels[method_order]),
        col = method_cols[method_order],
        lty = c("dashed", "solid", "solid", "solid"),
        lwd = c(2.0, 2.5, 2.5, 2.5),
        pch = c(1, 19, 19, 19),
        bty = "n",
        cex = 0.95
      )
    }
  }

  graphics::mtext(
    sprintf("Mesh-first KK / GKK / LGKK benchmark in d = %d", dim_value),
    outer = TRUE,
    cex = 1.2,
    font = 2
  )
  invisible(NULL)
}

write_gallery_png <- function(representatives, path) {
  methods <- c("target", "start", "kk", "gkk", "lgkk")
  method_titles <- c(
    target = "Target surface",
    start = "Shared start",
    kk = "KK",
    gkk = "KK->GKK",
    lgkk = "KK->LGKK"
  )
  surface_ids <- vapply(surface_configs, `[[`, "", "id")

  grDevices::png(path, width = 3000, height = 1800, res = 180, bg = "#f7f3ea")
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(length(surface_ids), length(methods)),
                mar = c(1.2, 1.2, 2.8, 0.4), oma = c(0, 0, 2.0, 0))

  for (surface_id in surface_ids) {
    rep_case <- representatives[[surface_id]]
    target <- rep_case$spec$coords_surface
    for (method in methods) {
      coords <- switch(
        method,
        target = target,
        start = rep_case$start,
        kk = rep_case$kk,
        gkk = rep_case$gkk,
        lgkk = rep_case$lgkk
      )
      if (!identical(method, "target")) {
        coords <- pkg_internal("grip.align.to.target.nd")(coords, target)$aligned
      }
      grip.plot(
        coords = coords,
        edges = rep_case$spec$edges,
        projection = "ortho",
        azimuth = 35,
        elevation = 24,
        vertex.col = if (identical(method, "target")) "#bc6c25" else "#3a5a40",
        edge.col = "#adb5bd",
        main = ""
      )
      graphics::title(
        main = sprintf("%s\n%s", rep_case$spec$surface_label, method_titles[[method]]),
        cex.main = 0.82
      )
    }
  }

  graphics::mtext(
    "Representative 3D mesh layouts (12x12, seed 1), shown as aligned orthographic projections",
    outer = TRUE,
    cex = 1.15,
    font = 2
  )
  invisible(NULL)
}

write_tabular <- function(df, column_labels, align, path, digits = 4L) {
  fmt_cell <- function(x) {
    if (is.numeric(x)) {
      if (length(x) == 1L && is.finite(x) && abs(x - round(x)) < 1e-10) {
        as.character(as.integer(round(x)))
      } else {
        fmt_num(x, digits = digits)
      }
    } else {
      as.character(x)
    }
  }

  lines <- c(
    sprintf("\\begin{tabular}{%s}", align),
    "\\toprule",
    paste(column_labels, collapse = " & "),
    "\\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(df))) {
    row <- vapply(df[i, , drop = FALSE], fmt_cell, character(1L))
    lines <- c(lines, paste(row, collapse = " & "), "\\\\")
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  writeLines(lines, con = path)
}

case_manifest_rows <- list()
raw_rows <- list()
representatives <- list()
row_used <- 0L

for (surface_cfg in surface_configs) {
  for (side in benchmark_cfg$sizes) {
    spec <- build_mesh_instance(surface_cfg, side)
    case_manifest_rows[[length(case_manifest_rows) + 1L]] <- data.frame(
      case_id = spec$case_id,
      case_label = spec$case_label,
      surface_id = spec$surface_id,
      surface_label = spec$surface_label,
      surface = spec$surface,
      size = spec$size,
      size_label = spec$size_label,
      connectivity = spec$connectivity,
      amplitude = if (!is.null(surface_cfg$amplitude)) surface_cfg$amplitude else NA_real_,
      freq_u = if (!is.null(surface_cfg$freq_u)) surface_cfg$freq_u else NA_real_,
      freq_v = if (!is.null(surface_cfg$freq_v)) surface_cfg$freq_v else NA_real_,
      x_scale = if (!is.null(surface_cfg$x_scale)) surface_cfg$x_scale else NA_real_,
      y_scale = if (!is.null(surface_cfg$y_scale)) surface_cfg$y_scale else NA_real_,
      vertices = spec$n,
      edges = nrow(spec$edges),
      edge_weight_cv = spec$edge_weight_cv,
      graph_diameter = spec$prepared_gkk$graph_diameter,
      stringsAsFactors = FALSE
    )

    for (dim_value in benchmark_cfg$dims) {
      for (seed in benchmark_cfg$seeds) {
        cat(sprintf(
          "[mesh-suite] case=%s dim=%d seed=%d\n",
          spec$case_id, dim_value, seed
        ))
        evaluated <- evaluate_layout(spec, dim = dim_value, seed = seed)
        for (piece in evaluated) {
          row_used <- row_used + 1L
          raw_rows[[row_used]] <- piece$row
        }
        if (side == max(benchmark_cfg$sizes) && dim_value == 3L && seed == 1L) {
          representatives[[spec$surface_id]] <- list(
            spec = spec,
            start = evaluated$start$coords,
            kk = evaluated$kk$coords,
            gkk = evaluated$gkk$coords,
            lgkk = evaluated$lgkk$coords
          )
        }
      }
    }
  }
}

case_manifest <- do.call(rbind, case_manifest_rows)
raw_results <- do.call(rbind, raw_rows)
raw_results <- raw_results[order(raw_results$dim,
                                 raw_results$surface_id,
                                 raw_results$size,
                                 raw_results$seed,
                                 match(raw_results$method, c("start", "kk", "gkk", "lgkk"))), ]
rownames(raw_results) <- NULL

main_results <- raw_results[raw_results$method != "start", , drop = FALSE]

case_summary <- summarize_metrics(
  raw_results,
  group_cols = c("case_id", "case_label", "surface_id", "surface_label", "size", "size_label", "dim", "method", "method_label"),
  metric_cols = c("kk_rel_rmse", "gkk_rel_rmse", "lgkk_sparse_rel_rmse",
                  "procrustes_rmse", "runtime_sec", "sampled_stress")
)
case_summary <- case_summary[order(case_summary$dim,
                                   case_summary$surface_id,
                                   case_summary$size,
                                   match(case_summary$method, c("start", "kk", "gkk", "lgkk"))), ]
rownames(case_summary) <- NULL

method_summary <- summarize_metrics(
  main_results,
  group_cols = c("dim", "method", "method_label"),
  metric_cols = c("kk_rel_rmse", "gkk_rel_rmse", "lgkk_sparse_rel_rmse",
                  "procrustes_rmse", "runtime_sec", "sampled_stress")
)
method_summary <- method_summary[order(method_summary$dim, match(method_summary$method, c("kk", "gkk", "lgkk"))), ]
rownames(method_summary) <- NULL

case_summary_dim3 <- case_summary[case_summary$dim == 3L & case_summary$method != "start", , drop = FALSE]
case_summary_dim2 <- case_summary[case_summary$dim == 2L & case_summary$method != "start", , drop = FALSE]

make_wide_case_table <- function(summary_df, metric = "gkk_rel_rmse", runtime_metric = "runtime_sec") {
  case_ids <- unique(summary_df$case_id)
  rows <- lapply(case_ids, function(case_id) {
    piece <- summary_df[summary_df$case_id == case_id, , drop = FALSE]
    piece <- piece[order(match(piece$method, c("kk", "gkk", "lgkk"))), , drop = FALSE]
    data.frame(
      geometry = piece$surface_label[[1L]],
      size = piece$size_label[[1L]],
      KK_rel = piece[[paste0(metric, "_mean")]][piece$method == "kk"],
      GKK_rel = piece[[paste0(metric, "_mean")]][piece$method == "gkk"],
      LGKK_rel = piece[[paste0(metric, "_mean")]][piece$method == "lgkk"],
      KK_sec = piece[[paste0(runtime_metric, "_mean")]][piece$method == "kk"],
      GKK_sec = piece[[paste0(runtime_metric, "_mean")]][piece$method == "gkk"],
      LGKK_sec = piece[[paste0(runtime_metric, "_mean")]][piece$method == "lgkk"],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

write.csv(case_manifest,
          file = file.path(data_dir, "case_manifest.csv"),
          row.names = FALSE)
write.csv(raw_results,
          file = file.path(data_dir, "raw_results.csv"),
          row.names = FALSE)
write.csv(case_summary,
          file = file.path(data_dir, "case_summary.csv"),
          row.names = FALSE)
write.csv(method_summary,
          file = file.path(data_dir, "method_summary.csv"),
          row.names = FALSE)

write_metric_panel_plot(
  summary_df = case_summary,
  dim_value = 3L,
  metric = "gkk_rel_rmse",
  ylab = "Full GKK relative RMSE",
  path = file.path(figure_dir, "mesh_gkk_rel_rmse_dim3.png"),
  log_scale = FALSE
)
write_metric_panel_plot(
  summary_df = case_summary,
  dim_value = 2L,
  metric = "gkk_rel_rmse",
  ylab = "Full GKK relative RMSE",
  path = file.path(figure_dir, "mesh_gkk_rel_rmse_dim2.png"),
  log_scale = FALSE
)
write_metric_panel_plot(
  summary_df = case_summary,
  dim_value = 3L,
  metric = "runtime_sec",
  ylab = "Elapsed time (sec)",
  path = file.path(figure_dir, "mesh_runtime_dim3.png"),
  log_scale = TRUE
)
write_metric_panel_plot(
  summary_df = case_summary,
  dim_value = 3L,
  metric = "kk_rel_rmse",
  ylab = "Classical KK relative RMSE",
  path = file.path(figure_dir, "mesh_classical_kk_rel_rmse_dim3.png"),
  log_scale = FALSE
)
write_gallery_png(
  representatives = representatives,
  path = file.path(figure_dir, "mesh_representative_3d_gallery.png")
)

write_tabular(
  df = data.frame(
    dimension = method_summary$dim,
    method = method_summary$method_label,
    gkk_rel_rmse = method_summary$gkk_rel_rmse_mean,
    kk_rel_rmse = method_summary$kk_rel_rmse_mean,
    lgkk_sparse_rel_rmse = method_summary$lgkk_sparse_rel_rmse_mean,
    procrustes_rmse = method_summary$procrustes_rmse_mean,
    runtime_sec = method_summary$runtime_sec_mean,
    stringsAsFactors = FALSE
  ),
  column_labels = c("Dim.", "Method", "Full GKK rel. RMSE", "Classical KK rel. RMSE",
                    "LGKK sparse rel. RMSE", "Target RMSE", "Runtime (s)"),
  align = "llrrrrr",
  path = file.path(table_dir, "mesh_method_summary.tex"),
  digits = 4L
)

write_tabular(
  df = make_wide_case_table(case_summary_dim3),
  column_labels = c("Geometry", "Size", "KK rel.", "GKK rel.", "LGKK rel.",
                    "KK (s)", "GKK (s)", "LGKK (s)"),
  align = "llrrrrrr",
  path = file.path(table_dir, "mesh_case_summary_dim3.tex"),
  digits = 4L
)

write_tabular(
  df = make_wide_case_table(case_summary_dim2),
  column_labels = c("Geometry", "Size", "KK rel.", "GKK rel.", "LGKK rel.",
                    "KK (s)", "GKK (s)", "LGKK (s)"),
  align = "llrrrrrr",
  path = file.path(table_dir, "mesh_case_summary_dim2.tex"),
  digits = 4L
)

findings_lines <- c(
  sprintf("Run tag: %s", run_tag),
  sprintf("Generated: %s", format(Sys.time(), tz = Sys.timezone(), usetz = TRUE)),
  "",
  "Method summary:",
  capture.output(print(method_summary, row.names = FALSE)),
  "",
  "Best method by case under full GKK relative RMSE:"
)

case_best <- do.call(rbind, lapply(split(main_results, list(main_results$dim, main_results$case_id), drop = TRUE), function(piece) {
  piece <- piece[order(piece$gkk_rel_rmse, piece$runtime_sec, piece$method_label), , drop = FALSE]
  piece[1L, c("dim", "case_label", "method_label", "gkk_rel_rmse", "runtime_sec")]
}))
rownames(case_best) <- NULL
findings_lines <- c(findings_lines, capture.output(print(case_best, row.names = FALSE)))

writeLines(findings_lines, con = file.path(output_root, "benchmark_notes.txt"))

cat(sprintf("Wrote benchmark outputs to %s\n", output_root))
