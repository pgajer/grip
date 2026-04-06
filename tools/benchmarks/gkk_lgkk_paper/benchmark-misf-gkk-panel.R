#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% args

out_arg <- grep("^--out=", args, value = TRUE)
out_override <- if (length(out_arg) > 0L) sub("^--out=", "", out_arg[[1L]]) else NULL

run_tag <- if (smoke) {
  sprintf("misf-gkk-panel-smoke-%s", format(Sys.Date(), "%Y-%m-%d"))
} else {
  sprintf("misf-gkk-panel-%s", format(Sys.Date(), "%Y-%m-%d"))
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

options(rgl.useNULL = TRUE)

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

method_order <- c("start", "kk", "gkk", "lgkk", "misf_auto", "misf_full")
method_labels <- c(
  start = "Shared start",
  kk = "KK",
  gkk = "KK->GKK",
  lgkk = "KK->LGKK",
  misf_auto = "MISF-GKK-auto",
  misf_full = "MISF-GKK-full"
)
plot_methods <- c("kk", "gkk", "lgkk", "misf_auto", "misf_full")
plot_colors <- c(
  kk = "#355070",
  gkk = "#9c6644",
  lgkk = "#588157",
  misf_auto = "#b56576",
  misf_full = "#6d597a"
)
plot_pch <- c(
  kk = 16,
  gkk = 17,
  lgkk = 15,
  misf_auto = 18,
  misf_full = 8
)
regime_order <- c("overlap", "scale")
stage_order <- c("top_level", "insertion", "refinement", "final_polish")

benchmark_cfg <- list(
  dim = 3L,
  seeds = if (smoke) 1L else 1:3,
  gkk_max_iter = if (smoke) 4L else 12L,
  lgkk_max_iter = if (smoke) 4L else 12L,
  direct_lgkk_local_nbrs = 10L,
  direct_lgkk_landmark_count = 10L,
  stress_sample = if (smoke) 500L else 1500L,
  jitter_xy = 0.08,
  jitter_z = 0.08,
  misf_num_init = 24L,
  misf_num_nbrs = 20L,
  misf_top_level_full_limit = if (smoke) 64L else 256L,
  misf_top_level_local_nbrs = 12L,
  misf_top_level_landmark_count = 8L,
  misf_top_level_restarts = if (smoke) 2L else 8L,
  misf_top_level_max_iter = if (smoke) 4L else 16L,
  misf_insertion_mode = "weighted_kk",
  misf_insertion_max_iter = if (smoke) 8L else 32L,
  misf_refinement_full_limit = if (smoke) 64L else 192L,
  misf_refinement_local_nbrs = 8L,
  misf_refinement_landmark_count = 4L,
  misf_refinement_anchor_weight = 0.05,
  misf_refinement_anchor_weight_end = 0.05,
  misf_refinement_continuation = "constant",
  misf_refinement_max_iter = if (smoke) 4L else 8L,
  misf_final_full_limit = if (smoke) 64L else 256L,
  misf_final_local_nbrs = 12L,
  misf_final_landmark_count = 8L,
  misf_final_max_iter = if (smoke) 4L else 8L
)

family_configs_full <- list(
  list(
    id = "mesh",
    label = "Mesh saddle",
    tiers = list(
      list(size_id = "t1", size_rank = 1L, size_label = "12x12", regime = "overlap",
           builder = function() mesh.surface.graph(12, 12, surface = "saddle", amplitude = 0.75)),
      list(size_id = "t2", size_rank = 2L, size_label = "16x16", regime = "overlap",
           builder = function() mesh.surface.graph(16, 16, surface = "saddle", amplitude = 0.75)),
      list(size_id = "t3", size_rank = 3L, size_label = "20x20", regime = "scale",
           builder = function() mesh.surface.graph(20, 20, surface = "saddle", amplitude = 0.75))
    )
  ),
  list(
    id = "torus",
    label = "Torus pinched",
    tiers = list(
      list(size_id = "t1", size_rank = 1L, size_label = "12x12", regime = "overlap",
           builder = function() torus.surface.graph(12, 12, surface = "pinched", amplitude = 0.22)),
      list(size_id = "t2", size_rank = 2L, size_label = "16x16", regime = "overlap",
           builder = function() torus.surface.graph(16, 16, surface = "pinched", amplitude = 0.22)),
      list(size_id = "t3", size_rank = 3L, size_label = "20x20", regime = "scale",
           builder = function() torus.surface.graph(20, 20, surface = "pinched", amplitude = 0.22))
    )
  ),
  list(
    id = "irregular_annulus",
    label = "Irregular annulus folded",
    tiers = list(
      list(size_id = "t1", size_rank = 1L, size_label = "8x32", regime = "overlap",
           builder = function() irregular.annulus.surface.graph(
             rings = 8, outer_count = 32, surface = "folded", amplitude = 0.45
           )),
      list(size_id = "t2", size_rank = 2L, size_label = "10x40", regime = "overlap",
           builder = function() irregular.annulus.surface.graph(
             rings = 10, outer_count = 40, surface = "folded", amplitude = 0.45
           )),
      list(size_id = "t3", size_rank = 3L, size_label = "12x48", regime = "scale",
           builder = function() irregular.annulus.surface.graph(
             rings = 12, outer_count = 48, surface = "folded", amplitude = 0.45
           ))
    )
  ),
  list(
    id = "irregular_torus",
    label = "Irregular torus pinched",
    tiers = list(
      list(size_id = "t1", size_rank = 1L, size_label = "9x18", regime = "overlap",
           builder = function() irregular.torus.surface.graph(
             major_rings = 9, tube_count = 18, surface = "pinched", amplitude = 0.18
           )),
      list(size_id = "t2", size_rank = 2L, size_label = "12x24", regime = "overlap",
           builder = function() irregular.torus.surface.graph(
             major_rings = 12, tube_count = 24, surface = "pinched", amplitude = 0.18
           )),
      list(size_id = "t3", size_rank = 3L, size_label = "14x28", regime = "scale",
           builder = function() irregular.torus.surface.graph(
             major_rings = 14, tube_count = 28, surface = "pinched", amplitude = 0.18
           ))
    )
  )
)

family_configs <- if (smoke) {
  keep_ids <- c("mesh", "torus")
  lapply(
    family_configs_full[vapply(family_configs_full, function(cfg) cfg$id %in% keep_ids, logical(1L))],
    function(cfg) {
      cfg$tiers <- cfg$tiers[vapply(cfg$tiers, function(tier) tier$size_rank %in% c(1L, 3L), logical(1L))]
      cfg
    }
  )
} else {
  family_configs_full
}

fmt_num <- function(x, digits = 4L) {
  if (!is.finite(x)) {
    return("NA")
  }
  formatC(x, format = "f", digits = digits)
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

base_param_coords <- function(spec) {
  coords <- as.matrix(spec$coords_param)
  if (ncol(coords) < 2L) {
    coords <- cbind(coords, rep(0, nrow(coords)))
  }
  coords[, 1:2, drop = FALSE]
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
  if (!nrow(df)) {
    return(data.frame())
  }
  split_key <- interaction(df[group_cols], drop = TRUE, lex.order = TRUE)
  pieces <- split(df, split_key)
  rows <- lapply(pieces, function(piece) {
    row <- piece[1L, group_cols, drop = FALSE]
    row$n_cases <- nrow(piece)
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

write_tabular <- function(df, column_labels, align, path, digits = 4L) {
  fmt_cell <- function(x) {
    if (is.numeric(x)) {
      if (length(x) == 1L && is.finite(x) && abs(x - round(x)) < 1e-10) {
        return(as.character(as.integer(round(x))))
      }
      return(fmt_num(x, digits = digits))
    }
    as.character(x)
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

extract_accepted_steps <- function(trace) {
  if (is.null(trace) || !is.data.frame(trace) || !nrow(trace)) {
    return(NA_integer_)
  }
  if (!all(c("iteration", "accepted") %in% names(trace))) {
    return(NA_integer_)
  }
  keep <- isTRUE(trace$accepted) | (!is.na(trace$accepted) & trace$accepted)
  if (!any(keep)) {
    return(NA_integer_)
  }
  as.integer(max(trace$iteration[keep], na.rm = TRUE))
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

prepare_case_spec <- function(family_cfg, tier_cfg) {
  spec <- tier_cfg$builder()
  spec$family_id <- family_cfg$id
  spec$family_label <- family_cfg$label
  spec$size_id <- tier_cfg$size_id
  spec$size_rank <- tier_cfg$size_rank
  spec$size_label <- tier_cfg$size_label
  spec$regime <- tier_cfg$regime
  spec$case_id <- sprintf("%s_%s", spec$family_id, spec$size_id)
  spec$graph <- igraph::graph_from_edgelist(as.matrix(spec$edges), directed = FALSE)
  spec$prepared_gkk <- grip.prepare.geodesic.kk(
    edges = spec$edges,
    n = spec$n,
    edge_weights = spec$edge_weights,
    tie_mode = "average"
  )
  spec$prepared_lgkk <- grip.prepare.landmark.geodesic.kk(
    edges = spec$edges,
    n = spec$n,
    edge_weights = spec$edge_weights,
    local_nbrs = benchmark_cfg$direct_lgkk_local_nbrs,
    landmark_count = benchmark_cfg$direct_lgkk_landmark_count
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

build_initial_layout <- function(spec, seed) {
  base2d <- normalize_coords(base_param_coords(spec))
  fam_index <- match(spec$family_id, vapply(family_configs, `[[`, "", "id"))
  draw_seed <- as.integer(100000L + 1000L * seed + 100L * spec$size_rank + 10L * fam_index)
  set.seed(draw_seed)
  coords <- cbind(base2d, 0)
  jitter <- matrix(
    stats::rnorm(spec$n * 3L, sd = benchmark_cfg$jitter_z),
    ncol = 3L
  )
  jitter[, 1:2] <- 0.45 * jitter[, 1:2, drop = FALSE]
  coords <- coords + jitter
  storage.mode(coords) <- "double"
  sweep(coords, 2L, colMeans(coords), "-", check.margin = FALSE)
}

score_layout <- function(spec, coords, seed) {
  gkk <- grip.score.geodesic.kk(coords = coords, prepared = spec$prepared_gkk)
  lgkk <- grip.score.landmark.geodesic.kk(coords = coords, prepared = spec$prepared_lgkk)
  classic <- score_classical_kk(coords, prepared = spec$prepared_gkk)
  edge_stats <- pkg_internal("grip.edge.length.stats")(coords, spec$edges)
  sampled_stress <- pkg_internal("grip.sampled.stress")(
    coords = coords,
    adj.list = spec$prepared_gkk$adj_list,
    weight.list = spec$prepared_gkk$weight_list,
    sample.size = benchmark_cfg$stress_sample,
    rng.seed = 9000L + 100L * spec$size_rank + 10L * seed +
      match(spec$family_id, vapply(family_configs, `[[`, "", "id"))
  )
  target <- spec$coords_surface
  procrustes_rmse <- pkg_internal("grip.align.to.target.nd")(coords, target)$rmse

  data.frame(
    case_id = spec$case_id,
    family_id = spec$family_id,
    family_label = spec$family_label,
    regime = spec$regime,
    size_id = spec$size_id,
    size_rank = spec$size_rank,
    size_label = spec$size_label,
    dim = benchmark_cfg$dim,
    seed = seed,
    vertices = spec$n,
    edges = nrow(spec$edges),
    edge_weight_cv = spec$edge_weight_cv,
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

decorate_row <- function(row,
                         method,
                         method_label,
                         runtime_sec = NA_real_,
                         accepted_steps = NA_integer_) {
  row$method <- method
  row$method_label <- method_label
  row$runtime_sec <- as.double(runtime_sec)
  row$accepted_steps <- accepted_steps
  row$misf.requested.mode <- NA_character_
  row$misf.final.score.mode <- NA_character_
  row$misf.final.score.requested.mode <- NA_character_
  row$misf.final.score.pair.count <- NA_integer_
  row$misf.top.level <- NA_integer_
  row$misf.top.level.n <- NA_integer_
  row$misf.top.level.pair.mode <- NA_character_
  row$misf.refinement.level.count <- NA_integer_
  row$misf.inserted.vertex.count <- NA_integer_
  row$misf.elapsed.top.level <- NA_real_
  row$misf.elapsed.insertion <- NA_real_
  row$misf.elapsed.refinement <- NA_real_
  row$misf.elapsed.final.polish <- NA_real_
  row$misf.elapsed.total <- NA_real_
  row$misf.stage.rows <- NA_integer_
  row
}

augment_misf_row <- function(row,
                             fit,
                             method,
                             method_label,
                             requested_mode,
                             runtime_sec) {
  row <- decorate_row(
    row = row,
    method = method,
    method_label = method_label,
    runtime_sec = runtime_sec,
    accepted_steps = NA_integer_
  )
  score <- fit$score[1L, , drop = FALSE]
  row$misf.requested.mode <- requested_mode
  row$misf.final.score.mode <- as.character(score$final.score.mode[[1L]])
  row$misf.final.score.requested.mode <- as.character(score$final.score.requested.mode[[1L]])
  row$misf.final.score.pair.count <- as.integer(score$final.score.pair.count[[1L]])
  row$misf.top.level <- as.integer(score$top.level[[1L]])
  row$misf.top.level.n <- as.integer(score$top.level.n[[1L]])
  row$misf.top.level.pair.mode <- as.character(score$top.level.pair.mode[[1L]])
  row$misf.refinement.level.count <- as.integer(score$refinement.level.count[[1L]])
  row$misf.inserted.vertex.count <- as.integer(score$inserted.vertex.count[[1L]])
  row$misf.elapsed.top.level <- as.double(score$elapsed.top.level[[1L]])
  row$misf.elapsed.insertion <- as.double(score$elapsed.insertion[[1L]])
  row$misf.elapsed.refinement <- as.double(score$elapsed.refinement[[1L]])
  row$misf.elapsed.final.polish <- as.double(score$elapsed.final.polish[[1L]])
  row$misf.elapsed.total <- as.double(score$elapsed.total[[1L]])
  row$misf.stage.rows <- if (!is.null(fit$stage_trace)) nrow(fit$stage_trace) else NA_integer_
  row
}

append_stage_trace <- function(trace, spec, seed, method, method_label) {
  if (is.null(trace) || !is.data.frame(trace) || !nrow(trace)) {
    return(NULL)
  }
  cbind(
    data.frame(
      case_id = spec$case_id,
      family_id = spec$family_id,
      family_label = spec$family_label,
      regime = spec$regime,
      size_id = spec$size_id,
      size_rank = spec$size_rank,
      size_label = spec$size_label,
      seed = seed,
      method = method,
      method_label = method_label,
      stringsAsFactors = FALSE
    ),
    trace,
    stringsAsFactors = FALSE
  )
}

run_misf_method <- function(spec,
                            pair_mode,
                            seed) {
  grip.optimize.misf.geodesic.kk(
    edges = spec$edges,
    n = spec$n,
    edge_weights = spec$edge_weights,
    tie_mode = "average",
    num_init = benchmark_cfg$misf_num_init,
    num_nbrs = benchmark_cfg$misf_num_nbrs,
    dim = benchmark_cfg$dim,
    top_level_pair_mode = pair_mode,
    top_level_full_limit = benchmark_cfg$misf_top_level_full_limit,
    top_level_local_nbrs = benchmark_cfg$misf_top_level_local_nbrs,
    top_level_landmark_count = benchmark_cfg$misf_top_level_landmark_count,
    top_level_restarts = benchmark_cfg$misf_top_level_restarts,
    top_level_max_iter = benchmark_cfg$misf_top_level_max_iter,
    top_level_init = "cmdscale",
    insertion_mode = benchmark_cfg$misf_insertion_mode,
    insertion_max_iter = benchmark_cfg$misf_insertion_max_iter,
    refinement_pair_mode = pair_mode,
    refinement_full_limit = benchmark_cfg$misf_refinement_full_limit,
    refinement_local_nbrs = benchmark_cfg$misf_refinement_local_nbrs,
    refinement_landmark_count = benchmark_cfg$misf_refinement_landmark_count,
    refinement_anchor_weight = benchmark_cfg$misf_refinement_anchor_weight,
    refinement_anchor_weight_end = benchmark_cfg$misf_refinement_anchor_weight_end,
    refinement_continuation = benchmark_cfg$misf_refinement_continuation,
    refinement_max_iter = benchmark_cfg$misf_refinement_max_iter,
    final_pair_mode = pair_mode,
    final_full_limit = benchmark_cfg$misf_final_full_limit,
    final_local_nbrs = benchmark_cfg$misf_final_local_nbrs,
    final_landmark_count = benchmark_cfg$misf_final_landmark_count,
    final_max_iter = benchmark_cfg$misf_final_max_iter,
    return_trace = TRUE,
    seed = seed
  )
}

evaluate_case <- function(spec, seed) {
  initial <- build_initial_layout(spec, seed = seed)
  out <- list(rows = list(), stage_rows = list(), repr = list())

  row_start <- decorate_row(
    score_layout(spec, initial, seed = seed),
    method = "start",
    method_label = method_labels[["start"]],
    runtime_sec = NA_real_,
    accepted_steps = NA_integer_
  )
  out$rows[[length(out$rows) + 1L]] <- row_start
  out$repr$start <- initial

  timed_kk <- system.time({
    coords_kk <- igraph::layout_with_kk(
      spec$graph,
      coords = initial,
      dim = benchmark_cfg$dim,
      weights = spec$edge_weights
    )
  })
  row_kk <- decorate_row(
    score_layout(spec, coords_kk, seed = seed),
    method = "kk",
    method_label = method_labels[["kk"]],
    runtime_sec = timed_kk[["elapsed"]],
    accepted_steps = NA_integer_
  )
  out$rows[[length(out$rows) + 1L]] <- row_kk
  out$repr$kk <- coords_kk

  if (identical(spec$regime, "overlap")) {
    timed_gkk <- system.time({
      fit_gkk <- grip.optimize.geodesic.kk(
        coords = coords_kk,
        prepared = spec$prepared_gkk,
        max_iter = benchmark_cfg$gkk_max_iter,
        scale_mode = "profiled",
        return_trace = TRUE
      )
    })
    row_gkk <- decorate_row(
      score_layout(spec, fit_gkk$coords, seed = seed),
      method = "gkk",
      method_label = method_labels[["gkk"]],
      runtime_sec = timed_gkk[["elapsed"]],
      accepted_steps = extract_accepted_steps(fit_gkk$trace)
    )
    out$rows[[length(out$rows) + 1L]] <- row_gkk
    out$repr$gkk <- fit_gkk$coords
  }

  timed_lgkk <- system.time({
    fit_lgkk <- grip.optimize.landmark.geodesic.kk(
      coords = coords_kk,
      prepared = spec$prepared_lgkk,
      max_iter = benchmark_cfg$lgkk_max_iter,
      return_trace = TRUE
    )
  })
  row_lgkk <- decorate_row(
    score_layout(spec, fit_lgkk$coords, seed = seed),
    method = "lgkk",
    method_label = method_labels[["lgkk"]],
    runtime_sec = timed_lgkk[["elapsed"]],
    accepted_steps = extract_accepted_steps(fit_lgkk$trace)
  )
  out$rows[[length(out$rows) + 1L]] <- row_lgkk
  out$repr$lgkk <- fit_lgkk$coords

  timed_misf_auto <- system.time({
    fit_misf_auto <- run_misf_method(spec, pair_mode = "auto", seed = seed)
  })
  row_misf_auto <- augment_misf_row(
    row = score_layout(spec, fit_misf_auto$coords, seed = seed),
    fit = fit_misf_auto,
    method = "misf_auto",
    method_label = method_labels[["misf_auto"]],
    requested_mode = "auto",
    runtime_sec = timed_misf_auto[["elapsed"]]
  )
  out$rows[[length(out$rows) + 1L]] <- row_misf_auto
  out$stage_rows[[length(out$stage_rows) + 1L]] <- append_stage_trace(
    fit_misf_auto$stage_trace,
    spec = spec,
    seed = seed,
    method = "misf_auto",
    method_label = method_labels[["misf_auto"]]
  )
  out$repr$misf_auto <- fit_misf_auto$coords

  if (identical(spec$regime, "overlap")) {
    timed_misf_full <- system.time({
      fit_misf_full <- run_misf_method(spec, pair_mode = "full", seed = seed)
    })
    row_misf_full <- augment_misf_row(
      row = score_layout(spec, fit_misf_full$coords, seed = seed),
      fit = fit_misf_full,
      method = "misf_full",
      method_label = method_labels[["misf_full"]],
      requested_mode = "full",
      runtime_sec = timed_misf_full[["elapsed"]]
    )
    out$rows[[length(out$rows) + 1L]] <- row_misf_full
    out$stage_rows[[length(out$stage_rows) + 1L]] <- append_stage_trace(
      fit_misf_full$stage_trace,
      spec = spec,
      seed = seed,
      method = "misf_full",
      method_label = method_labels[["misf_full"]]
    )
    out$repr$misf_full <- fit_misf_full$coords
  }

  out
}

extract_value <- function(df, method, column) {
  keep <- df$method == method
  if (!any(keep)) {
    return(NA_real_)
  }
  df[[column]][keep][1L]
}

make_overlap_table <- function(summary_df) {
  case_ids <- unique(summary_df$case_id)
  rows <- lapply(case_ids, function(case_id) {
    piece <- summary_df[summary_df$case_id == case_id, , drop = FALSE]
    piece <- piece[order(match(piece$method, method_order)), , drop = FALSE]
    data.frame(
      Family = piece$family_label[[1L]],
      Size = piece$size_label[[1L]],
      n = piece$vertices[[1L]],
      GKK_rel = extract_value(piece, "gkk", "gkk_rel_rmse_mean"),
      LGKK_rel = extract_value(piece, "lgkk", "gkk_rel_rmse_mean"),
      MISF_auto_rel = extract_value(piece, "misf_auto", "gkk_rel_rmse_mean"),
      MISF_full_rel = extract_value(piece, "misf_full", "gkk_rel_rmse_mean"),
      GKK_sec = extract_value(piece, "gkk", "runtime_sec_mean"),
      LGKK_sec = extract_value(piece, "lgkk", "runtime_sec_mean"),
      MISF_auto_sec = extract_value(piece, "misf_auto", "runtime_sec_mean"),
      MISF_full_sec = extract_value(piece, "misf_full", "runtime_sec_mean"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

make_scale_table <- function(summary_df) {
  case_ids <- unique(summary_df$case_id)
  rows <- lapply(case_ids, function(case_id) {
    piece <- summary_df[summary_df$case_id == case_id, , drop = FALSE]
    piece <- piece[order(match(piece$method, method_order)), , drop = FALSE]
    data.frame(
      Family = piece$family_label[[1L]],
      Size = piece$size_label[[1L]],
      n = piece$vertices[[1L]],
      KK_rel = extract_value(piece, "kk", "gkk_rel_rmse_mean"),
      LGKK_rel = extract_value(piece, "lgkk", "gkk_rel_rmse_mean"),
      MISF_auto_rel = extract_value(piece, "misf_auto", "gkk_rel_rmse_mean"),
      KK_sec = extract_value(piece, "kk", "runtime_sec_mean"),
      LGKK_sec = extract_value(piece, "lgkk", "runtime_sec_mean"),
      MISF_auto_sec = extract_value(piece, "misf_auto", "runtime_sec_mean"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

representative_rank <- function(family_cfg) {
  overlap_ranks <- vapply(
    family_cfg$tiers[vapply(family_cfg$tiers, function(tier) identical(tier$regime, "overlap"), logical(1L))],
    function(tier) tier$size_rank,
    integer(1L)
  )
  max(overlap_ranks)
}

write_metric_panel_plot <- function(summary_df, metric, ylab, path, log_scale = FALSE) {
  n_families <- length(family_configs)
  n_row <- if (n_families <= 2L) 1L else 2L
  n_col <- ceiling(n_families / n_row)
  grDevices::png(path, width = 2600, height = if (n_row == 1L) 900 else 1600, res = 180, bg = "#f7f3ea")
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(n_row, n_col), mar = c(4.5, 4.8, 3.0, 1.0), oma = c(0, 0, 1.8, 0))

  for (i in seq_along(family_configs)) {
    cfg <- family_configs[[i]]
    fam_df <- summary_df[summary_df$family_id == cfg$id, , drop = FALSE]
    fam_df <- fam_df[order(fam_df$size_rank, match(fam_df$method, plot_methods)), , drop = FALSE]
    size_rows <- fam_df[!duplicated(fam_df$size_rank), c("size_rank", "size_label", "regime"), drop = FALSE]
    x <- seq_len(nrow(size_rows))
    upper <- fam_df[[paste0(metric, "_mean")]] + fam_df[[paste0(metric, "_sd")]]
    upper <- upper[is.finite(upper) & upper > 0]
    ymax <- if (length(upper)) max(upper) else 1
    if (isTRUE(log_scale)) {
      vals <- fam_df[[paste0(metric, "_mean")]]
      vals <- vals[is.finite(vals) & vals > 0]
      ymin <- if (length(vals)) min(vals) * 0.8 else 1e-4
      ytop <- ymax * 1.3
    } else {
      ymin <- 0
      ytop <- ymax * 1.1
    }
    graphics::plot(
      x,
      rep(NA_real_, length(x)),
      type = "n",
      xlab = "Size tier",
      ylab = ylab,
      xaxt = "n",
      yaxt = "s",
      ylim = c(ymin, ytop),
      log = if (isTRUE(log_scale)) "y" else "",
      main = cfg$label
    )
    graphics::axis(1, at = x, labels = size_rows$size_label)
    graphics::grid(col = "#dedbd2", lty = "dotted")
    overlap_n <- sum(size_rows$regime == "overlap")
    if (overlap_n > 0L && overlap_n < length(x)) {
      graphics::abline(v = overlap_n + 0.5, col = "#adb5bd", lty = "dashed")
    }
    for (method in plot_methods) {
      piece <- fam_df[fam_df$method == method, , drop = FALSE]
      if (!nrow(piece)) {
        next
      }
      xp <- match(piece$size_rank, size_rows$size_rank)
      yp <- piece[[paste0(metric, "_mean")]]
      graphics::lines(
        xp,
        yp,
        type = "b",
        lwd = 2.4,
        pch = plot_pch[[method]],
        col = plot_colors[[method]]
      )
    }
    if (i == 1L) {
      graphics::legend(
        "topleft",
        legend = unname(method_labels[plot_methods]),
        col = plot_colors[plot_methods],
        lwd = 2.4,
        pch = unname(plot_pch[plot_methods]),
        bty = "n",
        cex = 0.9
      )
    }
  }

  graphics::mtext(
    sprintf("MISF-GKK benchmark in d = %d", benchmark_cfg$dim),
    outer = TRUE,
    cex = 1.2,
    font = 2
  )
  invisible(NULL)
}

write_pareto_plot <- function(summary_df, path) {
  plot_df <- summary_df[summary_df$method %in% plot_methods, , drop = FALSE]
  if (!nrow(plot_df)) {
    return(invisible(NULL))
  }
  family_pch <- c(mesh = 16, torus = 15, irregular_annulus = 17, irregular_torus = 18)

  grDevices::png(path, width = 1800, height = 1200, res = 180, bg = "#f7f3ea")
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mar = c(4.8, 5.0, 3.0, 1.0))

  x <- plot_df$runtime_sec_mean
  y <- plot_df$gkk_rel_rmse_mean
  xmin <- min(x[x > 0], na.rm = TRUE) * 0.8
  xmax <- max(x, na.rm = TRUE) * 1.25
  ymax <- max(y, na.rm = TRUE) * 1.1

  graphics::plot(
    x,
    y,
    type = "n",
    log = "x",
    xlim = c(xmin, xmax),
    ylim = c(0, ymax),
    xlab = "Runtime (seconds, log scale)",
    ylab = "Exact full-GKK relative RMSE",
    main = "Runtime-quality tradeoff across the MISF-GKK panel"
  )
  graphics::grid(col = "#dedbd2", lty = "dotted")

  for (i in seq_len(nrow(plot_df))) {
    graphics::points(
      plot_df$runtime_sec_mean[[i]],
      plot_df$gkk_rel_rmse_mean[[i]],
      pch = family_pch[[plot_df$family_id[[i]]]],
      col = plot_colors[[plot_df$method[[i]]]],
      bg = plot_colors[[plot_df$method[[i]]]],
      cex = 1.3
    )
  }

  graphics::legend(
    "topright",
    legend = unname(method_labels[plot_methods]),
    col = plot_colors[plot_methods],
    pch = 16,
    bty = "n",
    cex = 0.9,
    title = "Method"
  )
  graphics::legend(
    "bottomleft",
    legend = vapply(family_configs, `[[`, "", "label"),
    pch = unname(family_pch[vapply(family_configs, `[[`, "", "id")]),
    col = "#3d405b",
    bty = "n",
    cex = 0.9,
    title = "Family"
  )
  invisible(NULL)
}

write_gallery_png <- function(representatives, path) {
  family_ids <- vapply(family_configs, `[[`, "", "id")
  family_ids <- family_ids[family_ids %in% names(representatives)]
  if (!length(family_ids)) {
    return(invisible(NULL))
  }
  methods <- c("target", "kk", "gkk", "lgkk", "misf_auto", "misf_full")
  method_titles <- c(
    target = "Target surface",
    kk = "KK",
    gkk = "KK->GKK",
    lgkk = "KK->LGKK",
    misf_auto = "MISF-GKK-auto",
    misf_full = "MISF-GKK-full"
  )

  grDevices::png(path, width = 3400, height = 2100, res = 180, bg = "#f7f3ea")
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(
    mfrow = c(length(family_ids), length(methods)),
    mar = c(1.2, 1.2, 2.8, 0.4),
    oma = c(0, 0, 2.0, 0)
  )

  for (family_id in family_ids) {
    rep_case <- representatives[[family_id]]
    target <- rep_case$spec$coords_surface
    for (method in methods) {
      coords <- switch(
        method,
        target = target,
        kk = rep_case$kk,
        gkk = rep_case$gkk,
        lgkk = rep_case$lgkk,
        misf_auto = rep_case$misf_auto,
        misf_full = rep_case$misf_full
      )
      if (is.null(coords)) {
        graphics::plot.new()
        graphics::title(
          main = sprintf("%s\n%s", rep_case$spec$family_label, method_titles[[method]]),
          cex.main = 0.82
        )
        next
      }
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
        main = sprintf("%s\n%s", rep_case$spec$family_label, method_titles[[method]]),
        cex.main = 0.82
      )
    }
  }

  graphics::mtext(
    "Representative 3D overlap-regime layouts, aligned to the target geometry",
    outer = TRUE,
    cex = 1.15,
    font = 2
  )
  invisible(NULL)
}

write_stage_breakdown_plot <- function(raw_df, path) {
  misf_df <- raw_df[raw_df$method %in% c("misf_auto", "misf_full"), , drop = FALSE]
  if (!nrow(misf_df)) {
    return(invisible(NULL))
  }
  timing_summary <- summarize_metrics(
    misf_df,
    group_cols = c("regime", "method", "method_label"),
    metric_cols = c(
      "misf.elapsed.top.level",
      "misf.elapsed.insertion",
      "misf.elapsed.refinement",
      "misf.elapsed.final.polish"
    )
  )
  timing_summary <- timing_summary[order(match(timing_summary$regime, regime_order), match(timing_summary$method, plot_methods)), , drop = FALSE]
  if (!nrow(timing_summary)) {
    return(invisible(NULL))
  }

  comp_labels <- c(
    "Top level" = "misf.elapsed.top.level_mean",
    "Insertion" = "misf.elapsed.insertion_mean",
    "Refinement" = "misf.elapsed.refinement_mean",
    "Final polish" = "misf.elapsed.final.polish_mean"
  )
  mat <- rbind(
    timing_summary[[comp_labels[[1L]]]],
    timing_summary[[comp_labels[[2L]]]],
    timing_summary[[comp_labels[[3L]]]],
    timing_summary[[comp_labels[[4L]]]]
  )
  rownames(mat) <- names(comp_labels)
  colnames(mat) <- sprintf(
    "%s\n%s",
    ifelse(timing_summary$regime == "overlap", "Overlap", "Scale"),
    timing_summary$method_label
  )

  grDevices::png(path, width = 1800, height = 1100, res = 180, bg = "#f7f3ea")
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mar = c(5.5, 5.0, 3.0, 1.0))
  cols <- c("#6d597a", "#b56576", "#355070", "#588157")
  graphics::barplot(
    mat,
    beside = FALSE,
    col = cols,
    border = NA,
    ylab = "Mean elapsed time (seconds)",
    main = "MISF-GKK timing breakdown"
  )
  graphics::legend(
    "topright",
    legend = rownames(mat),
    fill = cols,
    bty = "n"
  )
  invisible(NULL)
}

case_manifest_rows <- list()
raw_rows <- list()
stage_rows <- list()
representatives <- list()
row_used <- 0L
stage_used <- 0L

for (family_cfg in family_configs) {
  rep_rank <- representative_rank(family_cfg)
  for (tier_cfg in family_cfg$tiers) {
    spec <- prepare_case_spec(family_cfg, tier_cfg)
    case_manifest_rows[[length(case_manifest_rows) + 1L]] <- data.frame(
      case_id = spec$case_id,
      family_id = spec$family_id,
      family_label = spec$family_label,
      regime = spec$regime,
      size_id = spec$size_id,
      size_rank = spec$size_rank,
      size_label = spec$size_label,
      vertices = spec$n,
      edges = nrow(spec$edges),
      edge_weight_cv = spec$edge_weight_cv,
      graph_diameter = spec$prepared_gkk$graph_diameter,
      stringsAsFactors = FALSE
    )

    for (seed in benchmark_cfg$seeds) {
      cat(sprintf(
        "[misf-gkk-panel] case=%s regime=%s seed=%d\n",
        spec$case_id, spec$regime, seed
      ))
      evaluated <- evaluate_case(spec, seed = seed)
      for (piece in evaluated$rows) {
        row_used <- row_used + 1L
        raw_rows[[row_used]] <- piece
      }
      for (piece in evaluated$stage_rows) {
        if (!is.null(piece)) {
          stage_used <- stage_used + 1L
          stage_rows[[stage_used]] <- piece
        }
      }
      if (seed == 1L && identical(spec$regime, "overlap") && spec$size_rank == rep_rank) {
        representatives[[spec$family_id]] <- c(list(spec = spec), evaluated$repr)
      }
    }
  }
}

case_manifest <- do.call(rbind, case_manifest_rows)
raw_results <- do.call(rbind, raw_rows)
raw_results <- raw_results[order(
  match(raw_results$regime, regime_order),
  raw_results$family_id,
  raw_results$size_rank,
  raw_results$seed,
  match(raw_results$method, method_order)
), ]
rownames(raw_results) <- NULL

stage_trace_raw <- if (length(stage_rows)) {
  out <- do.call(rbind, stage_rows)
  out <- out[order(
    match(out$regime, regime_order),
    out$family_id,
    out$size_rank,
    out$seed,
    match(out$method, method_order),
    match(out$stage, stage_order),
    out$level
  ), ]
  rownames(out) <- NULL
  out
} else {
  data.frame()
}

main_results <- raw_results[raw_results$method != "start", , drop = FALSE]

family_size_summary <- summarize_metrics(
  main_results,
  group_cols = c("case_id", "family_id", "family_label", "regime", "size_id", "size_rank", "size_label",
                 "vertices", "edges", "method", "method_label"),
  metric_cols = c("kk_rel_rmse", "gkk_rel_rmse", "lgkk_sparse_rel_rmse",
                  "procrustes_rmse", "runtime_sec", "sampled_stress",
                  "misf.elapsed.top.level", "misf.elapsed.insertion",
                  "misf.elapsed.refinement", "misf.elapsed.final.polish")
)
family_size_summary <- family_size_summary[order(
  match(family_size_summary$regime, regime_order),
  family_size_summary$family_id,
  family_size_summary$size_rank,
  match(family_size_summary$method, method_order)
), ]
rownames(family_size_summary) <- NULL

method_summary <- summarize_metrics(
  main_results,
  group_cols = c("regime", "method", "method_label"),
  metric_cols = c("kk_rel_rmse", "gkk_rel_rmse", "lgkk_sparse_rel_rmse",
                  "procrustes_rmse", "runtime_sec", "sampled_stress")
)
method_summary <- method_summary[order(match(method_summary$regime, regime_order), match(method_summary$method, method_order)), ]
rownames(method_summary) <- NULL

misf_stage_summary <- if (nrow(stage_trace_raw)) {
  out <- summarize_metrics(
    stage_trace_raw,
    group_cols = c("family_id", "family_label", "regime", "size_id", "size_rank", "size_label",
                   "method", "method_label", "stage", "pair.mode"),
    metric_cols = c("active_n", "inserted_n", "pair_n", "energy",
                    "weighted.rel.rmse", "elapsed_sec", "trace_rows", "frame_count")
  )
  out[order(match(out$regime, regime_order), out$family_id, out$size_rank, match(out$method, method_order), match(out$stage, stage_order)), , drop = FALSE]
} else {
  data.frame()
}
rownames(misf_stage_summary) <- NULL

pair_mode_usage_summary <- if (nrow(stage_trace_raw)) {
  out <- summarize_metrics(
    stage_trace_raw[!is.na(stage_trace_raw$pair.mode), , drop = FALSE],
    group_cols = c("regime", "method", "method_label", "stage", "pair.mode"),
    metric_cols = c("active_n", "pair_n", "elapsed_sec")
  )
  out[order(match(out$regime, regime_order), match(out$method, method_order), match(out$stage, stage_order), out$pair.mode), , drop = FALSE]
} else {
  data.frame()
}
rownames(pair_mode_usage_summary) <- NULL

write.csv(case_manifest, file = file.path(data_dir, "case_manifest.csv"), row.names = FALSE)
write.csv(raw_results, file = file.path(data_dir, "raw_results.csv"), row.names = FALSE)
write.csv(family_size_summary, file = file.path(data_dir, "family_size_summary.csv"), row.names = FALSE)
write.csv(method_summary, file = file.path(data_dir, "method_summary.csv"), row.names = FALSE)
write.csv(stage_trace_raw, file = file.path(data_dir, "stage_trace_raw.csv"), row.names = FALSE)
write.csv(misf_stage_summary, file = file.path(data_dir, "misf_stage_summary.csv"), row.names = FALSE)
write.csv(pair_mode_usage_summary, file = file.path(data_dir, "pair_mode_usage_summary.csv"), row.names = FALSE)

write_metric_panel_plot(
  summary_df = family_size_summary,
  metric = "gkk_rel_rmse",
  ylab = "Exact full-GKK relative RMSE",
  path = file.path(figure_dir, "misf_gkk_quality_by_size_dim3.png"),
  log_scale = FALSE
)
write_metric_panel_plot(
  summary_df = family_size_summary,
  metric = "runtime_sec",
  ylab = "Elapsed time (sec)",
  path = file.path(figure_dir, "misf_gkk_runtime_by_size_dim3.png"),
  log_scale = TRUE
)
write_pareto_plot(
  summary_df = family_size_summary,
  path = file.path(figure_dir, "misf_gkk_quality_runtime_pareto_dim3.png")
)
write_gallery_png(
  representatives = representatives,
  path = file.path(figure_dir, "misf_gkk_representative_gallery_dim3.png")
)
write_stage_breakdown_plot(
  raw_df = raw_results,
  path = file.path(figure_dir, "misf_gkk_stage_breakdown_dim3.png")
)

method_table <- data.frame(
  Regime = ifelse(method_summary$regime == "overlap", "Overlap", "Scale"),
  Method = method_summary$method_label,
  GKK_rel_RMSE = method_summary$gkk_rel_rmse_mean,
  Runtime_sec = method_summary$runtime_sec_mean,
  Procrustes_RMSE = method_summary$procrustes_rmse_mean,
  stringsAsFactors = FALSE
)

overlap_table <- make_overlap_table(
  family_size_summary[family_size_summary$regime == "overlap", , drop = FALSE]
)
scale_table <- make_scale_table(
  family_size_summary[family_size_summary$regime == "scale", , drop = FALSE]
)

write_tabular(
  df = method_table,
  column_labels = c("Regime", "Method", "Exact full-GKK rel. RMSE", "Runtime (s)", "Target RMSE"),
  align = "llrrr",
  path = file.path(table_dir, "misf_gkk_method_summary_dim3.tex"),
  digits = 4L
)

write_tabular(
  df = overlap_table,
  column_labels = c("Family", "Size", "$n$", "GKK rel.", "LGKK rel.", "MISF-auto rel.",
                    "MISF-full rel.", "GKK (s)", "LGKK (s)", "MISF-auto (s)", "MISF-full (s)"),
  align = "llrrrrrrrrr",
  path = file.path(table_dir, "misf_gkk_overlap_summary_dim3.tex"),
  digits = 4L
)

write_tabular(
  df = scale_table,
  column_labels = c("Family", "Size", "$n$", "KK rel.", "LGKK rel.", "MISF-auto rel.",
                    "KK (s)", "LGKK (s)", "MISF-auto (s)"),
  align = "llrrrrrrr",
  path = file.path(table_dir, "misf_gkk_scale_summary_dim3.tex"),
  digits = 4L
)

notes_lines <- c(
  sprintf("Run tag: %s", run_tag),
  sprintf("Generated: %s", format(Sys.time(), tz = Sys.timezone(), usetz = TRUE)),
  "",
  "Method summary by regime:",
  capture.output(print(method_summary, row.names = FALSE)),
  ""
)
if (nrow(pair_mode_usage_summary)) {
  notes_lines <- c(
    notes_lines,
    "Pair-mode usage summary:",
    capture.output(print(pair_mode_usage_summary, row.names = FALSE)),
    ""
  )
}
writeLines(notes_lines, con = file.path(output_root, "benchmark_notes.txt"))

cat(sprintf("Wrote MISF-GKK benchmark outputs to %s\n", output_root))
