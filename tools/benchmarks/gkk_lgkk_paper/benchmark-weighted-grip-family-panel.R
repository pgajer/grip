#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% args

out_arg <- grep("^--out=", args, value = TRUE)
out_override <- if (length(out_arg) > 0L) sub("^--out=", "", out_arg[[1L]]) else NULL

run_tag <- if (smoke) {
  sprintf("weighted-grip-phase5-smoke-%s", format(Sys.Date(), "%Y-%m-%d"))
} else {
  sprintf("weighted-grip-phase5-family-panel-%s", format(Sys.Date(), "%Y-%m-%d"))
}

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
output_root <- if (!is.null(out_override) && nzchar(out_override)) {
  normalizePath(out_override, winslash = "/", mustWork = FALSE)
} else if (smoke) {
  file.path(tempdir(), run_tag)
} else {
  file.path(repo_root, "output", "gkk_lgkk_paper", "benchmarks", run_tag)
}
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

benchmark_cfg <- list(
  dims = c(2L, 3L),
  seeds = if (smoke) 1L else 1:2,
  gkk_max_iter = if (smoke) 4L else 10L,
  lgkk_max_iter = if (smoke) 4L else 10L,
  lgkk_local_nbrs = if (smoke) 8L else 10L,
  lgkk_landmark_count = if (smoke) 8L else 10L,
  weighted_core_lgkk_rounds = if (smoke) 1L else 3L,
  weighted_core_lgkk_scope = "all",
  weighted_core_lgkk_active_limit = if (smoke) 512L else 4096L,
  weighted_lgkk_polish_rounds = if (smoke) 4L else 8L,
  stress_sample = if (smoke) 500L else 1500L,
  jitter_xy = 0.08,
  jitter_z = 0.08
)

family_configs_all <- list(
  list(
    id = "mesh",
    label = "Mesh saddle",
    preset = "mesh",
    builder = function() mesh.surface.graph(
      if (smoke) 6 else 8,
      if (smoke) 6 else 8,
      surface = "saddle",
      amplitude = 0.75
    )
  ),
  list(
    id = "cylinder",
    label = "Cylinder hourglass",
    preset = "cylinder",
    builder = function() cylinder.surface.graph(
      if (smoke) 6 else 8,
      if (smoke) 8 else 12,
      surface = "hourglass",
      amplitude = 0.30
    )
  ),
  list(
    id = "torus",
    label = "Torus pinched",
    preset = "torus",
    builder = function() torus.surface.graph(
      if (smoke) 6 else 8,
      if (smoke) 6 else 8,
      surface = "pinched",
      amplitude = 0.22
    )
  ),
  list(
    id = "sphere",
    label = "Sphere wavy",
    preset = "sphere",
    builder = function() sphere.surface.graph(
      if (smoke) 6 else 8,
      if (smoke) 8 else 10,
      surface = "wavy",
      amplitude = 0.18
    )
  ),
  list(
    id = "sierpinski_carpet",
    label = "Sierpinski carpet ripple",
    preset = "carpet",
    builder = function() sierpinski.carpet.surface.graph(
      level = 2,
      surface = "ripple",
      amplitude = 0.70,
      freq_u = 2,
      freq_v = 1
    )
  ),
  list(
    id = "cube_channel_network",
    label = "Cube channel network twisted",
    preset = "irregular",
    builder = function() cube.channel.network.surface.graph(
      level = if (smoke) 1 else 2,
      surface = "twisted",
      amplitude = 0.16,
      twist = 0.55
    )
  ),
  list(
    id = "irregular_annulus",
    label = "Irregular annulus folded",
    preset = "irregular",
    builder = function() irregular.annulus.surface.graph(
      rings = if (smoke) 5 else 6,
      outer_count = if (smoke) 18 else 24,
      surface = "folded",
      amplitude = 0.45
    )
  ),
  list(
    id = "irregular_torus",
    label = "Irregular torus pinched",
    preset = "irregular",
    builder = function() irregular.torus.surface.graph(
      major_rings = if (smoke) 6 else 7,
      tube_count = if (smoke) 10 else 14,
      surface = "pinched",
      amplitude = 0.18
    )
  )
)

family_configs <- if (smoke) {
  keep <- c("mesh", "torus", "sierpinski_carpet", "cube_channel_network", "irregular_torus")
  family_configs_all[vapply(family_configs_all, function(cfg) cfg$id %in% keep, logical(1L))]
} else {
  family_configs_all
}

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
  c(mean = mean(vals), sd = if (length(vals) > 1L) stats::sd(vals) else 0)
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

prepare_family_spec <- function(cfg) {
  spec <- cfg$builder()
  spec$family_id <- cfg$id
  spec$family_label <- cfg$label
  spec$preset <- cfg$preset
  spec$graph <- igraph::graph_from_edgelist(as.matrix(spec$edges), directed = FALSE)
  spec$prepared_gkk <- prepare.geodesic.kk(
    edges = spec$edges,
    n = spec$n,
    edge_weights = spec$edge_weights
  )
  spec$prepared_lgkk <- prepare.landmark.geodesic.kk(
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

base_param_coords <- function(spec) {
  coords <- as.matrix(spec$coords_param)
  if (ncol(coords) < 2L) {
    coords <- cbind(coords, rep(0, nrow(coords)))
  }
  coords[, 1:2, drop = FALSE]
}

build_initial_layout <- function(spec, dim, seed) {
  dim <- as.integer(dim)
  base2d <- normalize_coords(base_param_coords(spec))
  fam_index <- match(spec$family_id, vapply(family_configs, `[[`, "", "id"))
  draw_seed <- as.integer(100000L + 1000L * seed + 100L * dim + 10L * fam_index)
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
    jitter[, 1:2] <- 0.45 * jitter[, 1:2, drop = FALSE]
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
  gkk <- score.geodesic.kk(coords = coords, prepared = spec$prepared_gkk)
  lgkk <- score.landmark.geodesic.kk(coords = coords, prepared = spec$prepared_lgkk)
  classic <- score_classical_kk(coords, prepared = spec$prepared_gkk)
  edge_stats <- pkg_internal("grip.edge.length.stats")(coords, spec$edges)
  sampled_stress <- pkg_internal("grip.sampled.stress")(
    coords = coords,
    adj.list = spec$prepared_gkk$adj_list,
    weight.list = spec$prepared_gkk$weight_list,
    sample.size = benchmark_cfg$stress_sample,
    rng.seed = 9000L + 100L * dim + 10L * seed + match(spec$family_id, vapply(family_configs, `[[`, "", "id"))
  )
  target <- if (dim == 3L) {
    spec$coords_surface
  } else {
    base_param_coords(spec)
  }
  target_space <- if (dim == 3L) "surface_3d" else "parameter_2d"
  procrustes_rmse <- pkg_internal("grip.align.to.target.nd")(coords, target)$rmse

  data.frame(
    family_id = spec$family_id,
    family_label = spec$family_label,
    family = spec$family,
    surface = if ("surface" %in% names(spec)) spec$surface else NA_character_,
    preset = spec$preset,
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

resolve_weighted_preset_args <- function(preset, dim) {
  resolver <- pkg_internal("grip.resolve.weighted.preset")
  resolver(
    preset = preset,
    dim = dim,
    placement = "barycenter",
    placement_missing = TRUE,
    rounds = 160L,
    rounds_missing = TRUE,
    final_rounds = 384L,
    final_rounds_missing = TRUE,
    num_init = 24L,
    num_init_missing = TRUE,
    num_nbrs = 20L,
    num_nbrs_missing = TRUE,
    r = 0.03,
    r_missing = TRUE,
    s = 7.5,
    s_missing = TRUE,
    repulsion_factor = 2.5,
    repulsion_factor_missing = TRUE
  )
}

evaluate_case <- function(spec, dim, seed) {
  tuning <- resolve_weighted_preset_args(spec$preset, dim = dim)
  initial <- build_initial_layout(spec, dim = dim, seed = seed)
  seed_method <- function(offset) as.integer(1000L * seed + 100L * dim + offset)
  out <- list()

  row_start <- score_layout(spec, initial, dim = dim, seed = seed)
  row_start$method <- "start"
  row_start$method_label <- "Shared start"
  row_start$runtime_sec <- NA_real_
  out$start <- row_start

  timed_grip <- system.time({
    coords_grip <- globalrep.grip(
      edges = spec$edges,
      edge_weights = spec$edge_weights,
      n = spec$n,
      dim = dim,
      placement = tuning$placement,
      rounds = tuning$rounds,
      final_rounds = tuning$final_rounds,
      num_init = tuning$num_init,
      num_nbrs = tuning$num_nbrs,
      r = tuning$r,
      s = tuning$s,
      repulsion_factor = tuning$repulsion_factor,
      seed = seed_method(1L)
    )
  })
  row_grip <- score_layout(spec, coords_grip, dim = dim, seed = seed)
  row_grip$method <- "grip"
  row_grip$method_label <- "GRIP"
  row_grip$runtime_sec <- unname(timed_grip[["elapsed"]])
  out$grip <- row_grip

  timed_wgrip <- system.time({
    coords_wgrip <- globalrep.weighted.grip(
      edges = spec$edges,
      edge_weights = spec$edge_weights,
      n = spec$n,
      dim = dim,
      preset = spec$preset,
      seed = seed_method(2L)
    )
  })
  row_wgrip <- score_layout(spec, coords_wgrip, dim = dim, seed = seed)
  row_wgrip$method <- "wgrip"
  row_wgrip$method_label <- "Weighted GRIP"
  row_wgrip$runtime_sec <- unname(timed_wgrip[["elapsed"]])
  out$wgrip <- row_wgrip

  timed_wgrip_core_lgkk <- system.time({
    coords_wgrip_core_lgkk <- globalrep.weighted.grip(
      edges = spec$edges,
      edge_weights = spec$edge_weights,
      n = spec$n,
      dim = dim,
      preset = spec$preset,
      lgkk_multiscale_rounds = benchmark_cfg$weighted_core_lgkk_rounds,
      lgkk_local_nbrs = benchmark_cfg$lgkk_local_nbrs,
      lgkk_landmark_count = benchmark_cfg$lgkk_landmark_count,
      lgkk_multiscale_scope = benchmark_cfg$weighted_core_lgkk_scope,
      lgkk_active_limit = benchmark_cfg$weighted_core_lgkk_active_limit,
      seed = seed_method(3L)
    )
  })
  row_wgrip_core_lgkk <- score_layout(spec, coords_wgrip_core_lgkk, dim = dim, seed = seed)
  row_wgrip_core_lgkk$method <- "wgrip_core_lgkk"
  row_wgrip_core_lgkk$method_label <- "Weighted GRIP + core LGKK"
  row_wgrip_core_lgkk$runtime_sec <- unname(timed_wgrip_core_lgkk[["elapsed"]])
  out$wgrip_core_lgkk <- row_wgrip_core_lgkk

  timed_wgrip_lgkk <- system.time({
    coords_wgrip_lgkk <- globalrep.weighted.grip(
      edges = spec$edges,
      edge_weights = spec$edge_weights,
      n = spec$n,
      dim = dim,
      preset = spec$preset,
      lgkk_polish_rounds = benchmark_cfg$weighted_lgkk_polish_rounds,
      lgkk_local_nbrs = benchmark_cfg$lgkk_local_nbrs,
      lgkk_landmark_count = benchmark_cfg$lgkk_landmark_count,
      seed = seed_method(4L)
    )
  })
  row_wgrip_lgkk <- score_layout(spec, coords_wgrip_lgkk, dim = dim, seed = seed)
  row_wgrip_lgkk$method <- "wgrip_polish_lgkk"
  row_wgrip_lgkk$method_label <- "Weighted GRIP + polish LGKK"
  row_wgrip_lgkk$runtime_sec <- unname(timed_wgrip_lgkk[["elapsed"]])
  out$wgrip_polish_lgkk <- row_wgrip_lgkk

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
  out$kk <- row_kk

  timed_gkk <- system.time({
    opt_gkk <- geodesic.kk(
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
  out$gkk <- row_gkk

  timed_lgkk <- system.time({
    opt_lgkk <- landmark.geodesic.kk(
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
  out$lgkk <- row_lgkk

  out
}

write_metric_plot <- function(summary_df, metric, ylab, path, log_scale = FALSE) {
  families <- vapply(family_configs, `[[`, "", "id")
  family_labels <- vapply(family_configs, `[[`, "", "label")
  method_order <- c("grip", "wgrip", "wgrip_core_lgkk", "wgrip_polish_lgkk", "kk", "gkk", "lgkk")
  method_labels <- c(
    grip = "GRIP",
    wgrip = "Weighted GRIP",
    wgrip_core_lgkk = "Weighted GRIP + core LGKK",
    wgrip_polish_lgkk = "Weighted GRIP + polish LGKK",
    kk = "KK",
    gkk = "KK->GKK",
    lgkk = "KK->LGKK"
  )
  method_cols <- c(
    grip = "#6c757d",
    wgrip = "#355070",
    wgrip_core_lgkk = "#588157",
    wgrip_polish_lgkk = "#2a9d8f",
    kk = "#9c6644",
    gkk = "#bc6c25",
    lgkk = "#6a4c93"
  )
  panel_df <- summary_df[summary_df$method %in% method_order, , drop = FALSE]
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
  offsets <- seq(-0.30, 0.30, length.out = length(method_order))

  grDevices::png(path, width = 2200, height = 1100, res = 180, bg = "#f7f3ea")
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(7.5, 4.8, 3.2, 1.0), oma = c(0, 0, 1.2, 0))

  for (dim_value in benchmark_cfg$dims) {
    df_dim <- panel_df[panel_df$dim == dim_value, , drop = FALSE]
    graphics::plot(
      c(0.5, length(families) + 0.5),
      c(ymin, if (isTRUE(log_scale)) ymax * 1.2 else ymax * 1.08),
      type = "n",
      log = if (isTRUE(log_scale)) "y" else "",
      xaxt = "n",
      xlab = "",
      ylab = ylab,
      main = sprintf("d = %d", dim_value)
    )
    graphics::axis(1, at = seq_along(families), labels = family_labels, las = 2, cex.axis = 0.85)
    graphics::grid(col = "#dedbd2", lty = "dotted")
    for (i in seq_along(method_order)) {
      method <- method_order[[i]]
      df_m <- df_dim[df_dim$method == method, , drop = FALSE]
      x <- match(df_m$family_id, families) + offsets[[i]]
      y <- df_m[[paste0(metric, "_mean")]]
      y_sd <- df_m[[paste0(metric, "_sd")]]
      if (isTRUE(log_scale)) {
        y <- pmax(y, ymin * 1.0001)
      }
      graphics::points(x, y, pch = 19, cex = 1.1, col = method_cols[[method]])
      keep <- is.finite(y) & is.finite(y_sd) & (y_sd > 0)
      if (any(keep)) {
        lower <- if (isTRUE(log_scale)) {
          pmax(ymin * 1.0001, y[keep] - y_sd[keep])
        } else {
          pmax(ymin, y[keep] - y_sd[keep])
        }
        graphics::arrows(
          x0 = x[keep],
          y0 = lower,
          x1 = x[keep],
          y1 = y[keep] + y_sd[keep],
          angle = 90,
          code = 3,
          length = 0.04,
          col = method_cols[[method]]
        )
      }
    }
    graphics::legend(
      "topright",
      legend = unname(method_labels[method_order]),
      col = method_cols[method_order],
      pch = 19,
      bty = "n",
      cex = 0.85
    )
  }

  graphics::mtext(
    "Phase 5 weighted GRIP family benchmark",
    outer = TRUE,
    cex = 1.15,
    font = 2
  )
  invisible(NULL)
}

write_readme <- function(path, specs, summary_df) {
  lines <- c(
    "# Weighted GRIP Phase 5 benchmark",
    "",
    sprintf("- Weighted core LGKK rounds: `%d`", benchmark_cfg$weighted_core_lgkk_rounds),
    sprintf("- Weighted core LGKK scope: `%s`", benchmark_cfg$weighted_core_lgkk_scope),
    sprintf("- Weighted core LGKK active limit: `%d`", benchmark_cfg$weighted_core_lgkk_active_limit),
    sprintf("- Weighted post LGKK polish rounds: `%d`", benchmark_cfg$weighted_lgkk_polish_rounds),
    "",
    sprintf("- Run tag: `%s`", run_tag),
    sprintf("- Smoke mode: `%s`", if (smoke) "yes" else "no"),
    sprintf("- Seeds: `%s`", paste(benchmark_cfg$seeds, collapse = ", ")),
    sprintf("- Dimensions: `%s`", paste(benchmark_cfg$dims, collapse = ", ")),
    "",
    "## Family panel",
    ""
  )
  for (spec in specs) {
    lines <- c(
        lines,
        sprintf("- `%s`: %s, preset `%s`, n = %d, m = %d", spec$family_id, spec$family_label, spec$preset, spec$n, nrow(spec$edges))
      )
  }
  lines <- c(
    lines,
    "",
    "## Output files",
    "",
    "- `data/raw_metrics.csv`: per-run metrics",
    "- `data/summary_metrics.csv`: mean/sd summary by family, dimension, and method",
    "- `tables/summary_gkk_rel_rmse.csv`: compact GKK-relative-RMSE table",
    "- `figures/gkk_rel_rmse.png`: method comparison by family and dimension",
    "- `figures/runtime_sec.png`: runtime comparison by family and dimension",
    ""
  )
  if (nrow(summary_df) > 0L) {
    best_rows <- do.call(rbind, lapply(split(summary_df, interaction(summary_df$family_id, summary_df$dim, drop = TRUE)), function(df) {
      metric <- df$gkk_rel_rmse_mean
      df[which.min(metric), , drop = FALSE]
    }))
    lines <- c(lines, "## Best method by family and dimension", "")
    for (i in seq_len(nrow(best_rows))) {
      row <- best_rows[i, , drop = FALSE]
      lines <- c(
        lines,
        sprintf(
          "- `%s`, d = %d: %s (mean GKK rel. RMSE = %s)",
          row$family_id[[1L]],
          row$dim[[1L]],
          row$method_label[[1L]],
          fmt_num(row$gkk_rel_rmse_mean[[1L]], digits = 4L)
        )
      )
    }
  }
  writeLines(lines, con = path)
}

family_specs <- lapply(family_configs, prepare_family_spec)

raw_rows <- list()
row_used <- 0L
for (spec in family_specs) {
  for (dim in benchmark_cfg$dims) {
    for (seed in benchmark_cfg$seeds) {
      result <- evaluate_case(spec, dim = dim, seed = seed)
      for (method in names(result)) {
        row_used <- row_used + 1L
        raw_rows[[row_used]] <- result[[method]]
      }
    }
  }
}
raw_df <- do.call(rbind, raw_rows)
raw_df <- raw_df[order(raw_df$family_id, raw_df$dim, raw_df$seed, raw_df$method), , drop = FALSE]

summary_df <- summarize_metrics(
  raw_df,
  group_cols = c("family_id", "family_label", "preset", "dim", "method", "method_label"),
  metric_cols = c(
    "runtime_sec",
    "gkk_rel_rmse",
    "lgkk_sparse_rel_rmse",
    "kk_rel_rmse",
    "sampled_stress",
    "procrustes_rmse",
    "edge_length_cv"
  )
)
summary_df <- summary_df[order(summary_df$family_id, summary_df$dim, summary_df$method), , drop = FALSE]

utils::write.csv(raw_df, file.path(data_dir, "raw_metrics.csv"), row.names = FALSE)
utils::write.csv(summary_df, file.path(data_dir, "summary_metrics.csv"), row.names = FALSE)

compact <- summary_df[, c("family_id", "dim", "method_label", "runtime_sec_mean", "gkk_rel_rmse_mean", "sampled_stress_mean", "procrustes_rmse_mean")]
utils::write.csv(compact, file.path(table_dir, "summary_gkk_rel_rmse.csv"), row.names = FALSE)

write_metric_plot(
  summary_df,
  metric = "gkk_rel_rmse",
  ylab = "Mean full-GKK relative RMSE",
  path = file.path(figure_dir, "gkk_rel_rmse.png"),
  log_scale = FALSE
)
write_metric_plot(
  summary_df,
  metric = "runtime_sec",
  ylab = "Mean runtime (sec)",
  path = file.path(figure_dir, "runtime_sec.png"),
  log_scale = TRUE
)
write_readme(file.path(output_root, "README.md"), specs = family_specs, summary_df = summary_df)

message(sprintf("Weighted GRIP Phase 5 benchmark completed: %s", output_root))
