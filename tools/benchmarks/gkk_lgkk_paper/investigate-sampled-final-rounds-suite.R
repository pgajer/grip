#!/usr/bin/env Rscript

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Package 'pkgload' is required to run the sampled final-rounds suite.")
}
if (!requireNamespace("callr", quietly = TRUE)) {
  stop("Package 'callr' is required to run the sampled final-rounds suite.")
}

pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)
sys.source(file.path("tools", "benchmarks", "gkk_lgkk_paper", "benchmark-globalrep-fixed-candidate.R"), envir = environment())

output_root <- file.path("dev", "manual", "tmp", "sampled-final-rounds-suite-2026-03-28")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
comparison_root <- file.path(output_root, "png-comparisons")
smoothness_root <- file.path(output_root, "png-smoothness")
dir.create(comparison_root, recursive = TRUE, showWarnings = FALSE)
dir.create(smoothness_root, recursive = TRUE, showWarnings = FALSE)

package_internal <- function(name) get(name, envir = asNamespace("grip"))

all_suite_specs <- Filter(function(spec) {
  !identical(spec$family, "sphere") && max(spec$edges) <= 100000L
}, build_graph_specs())

sampled_finals <- 0:32

parse_named_args <- function(args) {
  out <- list()
  for (arg in args) {
    parts <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) {
      stop("Arguments must use key=value format")
    }
    out[[parts[[1L]]]] <- parts[[2L]]
  }
  out
}

parse_int_scalar <- function(x, name) {
  val <- suppressWarnings(as.integer(x))
  if (length(val) != 1L || is.na(val)) {
    stop(sprintf("%s must be a single integer", name))
  }
  val
}

parse_char_vector <- function(x) {
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  parts[nzchar(parts)]
}

seed_budget <- function(n_vertices) {
  if (n_vertices <= 500L) {
    return(1:3)
  }
  if (n_vertices <= 5000L) {
    return(1:2)
  }
  1L
}

time_limit_budget <- function(n_vertices) {
  if (n_vertices <= 1000L) {
    return(60L)
  }
  if (n_vertices <= 10000L) {
    return(180L)
  }
  180L
}

safe_slug <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "-", x)
}

sampled_cfg <- function(n, final_rounds) {
  cfg <- package_internal("grip.globalrep.base.defaults")(n)
  cfg$final_rounds <- as.integer(final_rounds)
  cfg$coarse_repulsion_exact_below <- 64L
  cfg
}

default_cfg <- function(n) {
  package_internal("grip.globalrep.base.defaults")(n)
}

run_layout_external <- function(spec, cfg, seed) {
  callr::r(
    function(edges, dim, cfg, seed) {
      pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)
      started <- proc.time()[["elapsed"]]
      coords <- grip.layout(
        edges = edges,
        n = max(edges),
        dim = dim,
        placement = cfg$placement,
        rounds = cfg$rounds,
        final_rounds = cfg$final_rounds,
        num_init = cfg$num_init,
        num_nbrs = cfg$num_nbrs,
        r = cfg$r,
        s = cfg$s,
        repulsion_factor = cfg$repulsion_factor,
        coarse_repulsion_factor = cfg$coarse_repulsion_factor,
        coarse_repulsion_sample = cfg$coarse_repulsion_sample,
        coarse_repulsion_exact_below = cfg$coarse_repulsion_exact_below,
        seed = seed
      )
      list(coords = coords, elapsed_sec = proc.time()[["elapsed"]] - started)
    },
    args = list(edges = spec$edges, dim = spec$dim, cfg = cfg, seed = seed),
    timeout = time_limit_budget(max(spec$edges)),
    stdout = "|",
    stderr = "|"
  )
}

run_layout_safe <- function(spec, cfg, seed) {
  n <- max(spec$edges)
  result <- tryCatch({
    if (n > 10000L) {
      run_layout_external(spec, cfg, seed)
    } else {
      started <- proc.time()[["elapsed"]]
      coords <- grip.layout(
        edges = spec$edges,
        n = n,
        dim = spec$dim,
        placement = cfg$placement,
        rounds = cfg$rounds,
        final_rounds = cfg$final_rounds,
        num_init = cfg$num_init,
        num_nbrs = cfg$num_nbrs,
        r = cfg$r,
        s = cfg$s,
        repulsion_factor = cfg$repulsion_factor,
        coarse_repulsion_factor = cfg$coarse_repulsion_factor,
        coarse_repulsion_sample = cfg$coarse_repulsion_sample,
        coarse_repulsion_exact_below = cfg$coarse_repulsion_exact_below,
        seed = seed
      )
      list(coords = coords, elapsed_sec = proc.time()[["elapsed"]] - started)
    }
  }, error = function(e) e)

  if (inherits(result, "error")) {
    status <- if (inherits(result, "callr_timeout_error") ||
      grepl("timeout|timed out|elapsed time limit", conditionMessage(result), ignore.case = TRUE)) {
      "timeout"
    } else {
      "error"
    }
    return(list(status = status, error_message = conditionMessage(result), coords = NULL, elapsed_sec = NA_real_))
  }

  list(status = "ok", error_message = "", coords = result$coords, elapsed_sec = result$elapsed_sec)
}

geometry_metrics <- function(coords, spec, seed) {
  package_internal("grip.geometry.diagnostics")(
    coords = coords,
    target.coords = spec$canonical,
    edges = spec$edges,
    family = spec$family,
    sample.size.symmetry = 512L,
    sample.size.wedges = 4000L,
    rng.seed = 5000L + seed
  )
}

score_metrics <- function(coords, spec, seed) {
  grip.score.layout(
    coords = coords,
    edges = spec$edges,
    n = max(spec$edges),
    sample.size.stress = sample_budget(max(spec$edges))$stress,
    sample.size.nonedge = sample_budget(max(spec$edges))$sep,
    edge.crossings = "never"
  )
}

run_one <- function(spec, final_rounds, seed) {
  cfg <- sampled_cfg(max(spec$edges), final_rounds)
  layout <- run_layout_safe(spec, cfg, seed)
  if (!identical(layout$status, "ok")) {
    return(data.frame(
      family = spec$family,
      graph_label = spec$graph_label,
      title = spec$title,
      dim = spec$dim,
      final_rounds = final_rounds,
      config_id = sprintf("sampled_f%d", final_rounds),
      seed = seed,
      vertices = max(spec$edges),
      edges = nrow(spec$edges),
      status = layout$status,
      error_message = layout$error_message,
      elapsed_sec = layout$elapsed_sec,
      procrustes_rmse = NA_real_,
      sampled_stress = NA_real_,
      edge_length_cv = NA_real_,
      sampled_nonedge_sep_ratio = NA_real_,
      global_symmetry_score = NA_real_,
      local_angle_deviation = NA_real_,
      edge_axis_concentration = NA_real_,
      edge_axis_deviation = NA_real_,
      boundary_waviness = NA_real_,
      corridor_waviness = NA_real_,
      hole_center_error = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  score <- score_metrics(layout$coords, spec, seed)
  geom <- geometry_metrics(layout$coords, spec, seed)
  data.frame(
    family = spec$family,
    graph_label = spec$graph_label,
    title = spec$title,
    dim = spec$dim,
    final_rounds = final_rounds,
    config_id = sprintf("sampled_f%d", final_rounds),
    seed = seed,
    vertices = max(spec$edges),
    edges = nrow(spec$edges),
    status = "ok",
    error_message = "",
    elapsed_sec = layout$elapsed_sec,
    procrustes_rmse = geom$procrustes.rmse[[1L]],
    sampled_stress = score$sampled.stress[[1L]],
    edge_length_cv = score$edge.length.cv[[1L]],
    sampled_nonedge_sep_ratio = score$sampled.nonedge.sep.ratio[[1L]],
    global_symmetry_score = geom$global.symmetry.score[[1L]],
    local_angle_deviation = geom$local.angle.deviation[[1L]],
    edge_axis_concentration = geom$edge.axis.concentration[[1L]],
    edge_axis_deviation = geom$edge.axis.deviation[[1L]],
    boundary_waviness = geom$boundary.waviness[[1L]],
    corridor_waviness = geom$corridor.waviness[[1L]],
    hole_center_error = geom$hole.center.error[[1L]],
    stringsAsFactors = FALSE
  )
}

rank01_na <- function(x, higher_better = FALSE) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) {
    return(out)
  }
  vals <- if (higher_better) -x[ok] else x[ok]
  if (length(vals) == 1L || all(abs(vals - vals[[1L]]) <= sqrt(.Machine$double.eps))) {
    out[ok] <- 0
    return(out)
  }
  ranks <- rank(vals, ties.method = "average")
  out[ok] <- (ranks - 1) / (length(ranks) - 1)
  out
}

add_rank_scores <- function(summary_df) {
  groups <- split(summary_df, paste(summary_df$family, summary_df$graph_label, sep = "::"))
  scored <- lapply(groups, function(df) {
    df <- df[order(df$final_rounds), , drop = FALSE]
    geom.cols <- list(
      global_symmetry_score = rank01_na(df$global_symmetry_score, higher_better = TRUE),
      local_angle_deviation = rank01_na(df$local_angle_deviation, higher_better = FALSE),
      edge_axis_deviation = rank01_na(df$edge_axis_deviation, higher_better = FALSE),
      boundary_waviness = rank01_na(df$boundary_waviness, higher_better = FALSE),
      corridor_waviness = rank01_na(df$corridor_waviness, higher_better = FALSE),
      hole_center_error = rank01_na(df$hole_center_error, higher_better = FALSE)
    )
    qual.cols <- list(
      procrustes_rmse = rank01_na(df$procrustes_rmse, higher_better = FALSE),
      edge_length_cv = rank01_na(df$edge_length_cv, higher_better = FALSE),
      sampled_stress = rank01_na(df$sampled_stress, higher_better = FALSE),
      sampled_nonedge_sep_ratio = rank01_na(df$sampled_nonedge_sep_ratio, higher_better = TRUE)
    )
    df$geometry_score <- vapply(seq_len(nrow(df)), function(i) {
      vals <- vapply(geom.cols, `[[`, numeric(1L), i)
      mean(vals[is.finite(vals)])
    }, numeric(1L))
    df$full_score <- vapply(seq_len(nrow(df)), function(i) {
      vals <- c(
        vapply(geom.cols, `[[`, numeric(1L), i),
        vapply(qual.cols, `[[`, numeric(1L), i)
      )
      mean(vals[is.finite(vals)])
    }, numeric(1L))
    df
  })
  do.call(rbind, scored)
}

summarize_runs <- function(raw_df) {
  ok <- raw_df[raw_df$status == "ok", , drop = FALSE]
  split_df <- split(ok, paste(ok$family, ok$graph_label, ok$config_id, sep = "::"))
  out <- do.call(rbind, lapply(split_df, function(df) {
    data.frame(
      family = df$family[[1L]],
      graph_label = df$graph_label[[1L]],
      title = df$title[[1L]],
      dim = df$dim[[1L]],
      final_rounds = df$final_rounds[[1L]],
      config_id = df$config_id[[1L]],
      vertices = df$vertices[[1L]],
      edges = df$edges[[1L]],
      seeds = nrow(df),
      elapsed_sec = mean(df$elapsed_sec),
      procrustes_rmse = mean(df$procrustes_rmse),
      sampled_stress = mean(df$sampled_stress),
      edge_length_cv = mean(df$edge_length_cv),
      sampled_nonedge_sep_ratio = mean(df$sampled_nonedge_sep_ratio),
      global_symmetry_score = mean(df$global_symmetry_score),
      local_angle_deviation = mean(df$local_angle_deviation),
      edge_axis_concentration = mean(df$edge_axis_concentration),
      edge_axis_deviation = mean(df$edge_axis_deviation),
      boundary_waviness = mean(df$boundary_waviness),
      corridor_waviness = mean(df$corridor_waviness),
      hole_center_error = mean(df$hole_center_error),
      stringsAsFactors = FALSE
    )
  }))
  add_rank_scores(out)
}

project_for_panel <- function(coords, spec) {
  if (spec$dim == 2L) {
    return(package_internal("grip.normalize.coords")(coords))
  }
  fit <- align_layout_to_spec(coords, spec)
  fit$aligned[, 1:2, drop = FALSE]
}

project_canonical_for_panel <- function(spec) {
  if (spec$dim == 2L) {
    return(package_internal("grip.normalize.coords")(spec$canonical))
  }
  fit <- align_layout_to_spec(spec$canonical, spec)
  fit$aligned[, 1:2, drop = FALSE]
}

draw_contact_sheet <- function(path, spec, panel_defs) {
  grDevices::png(path, width = 2400, height = 2400, res = 170, bg = "#f7f3ea")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(3, 3), mar = c(0, 0, 3, 0), xaxs = "i", yaxs = "i")
  for (panel in panel_defs) {
    plot_layout_panel(
      panel$coords,
      spec$edges,
      title_text = panel$title,
      subtitle_text = panel$subtitle
    )
  }
}

draw_metric_curve <- function(path, graph_df) {
  df <- graph_df[order(graph_df$final_rounds), , drop = FALSE]
  metric_names <- c("global_symmetry_score", "local_angle_deviation", "edge_axis_deviation",
                    "boundary_waviness", "corridor_waviness", "hole_center_error")
  labels <- c(
    global_symmetry_score = "symmetry",
    local_angle_deviation = "angle dev",
    edge_axis_deviation = "axis dev",
    boundary_waviness = "boundary",
    corridor_waviness = "corridor",
    hole_center_error = "hole center"
  )
  cols <- c(
    global_symmetry_score = "#355070",
    local_angle_deviation = "#b23a48",
    edge_axis_deviation = "#3d7a57",
    boundary_waviness = "#8c564b",
    corridor_waviness = "#7c4d8b",
    hole_center_error = "#d77a1f"
  )
  scaled <- lapply(metric_names, function(name) {
    x <- df[[name]]
    if (identical(name, "global_symmetry_score")) {
      rank01_na(x, higher_better = TRUE)
    } else {
      rank01_na(x, higher_better = FALSE)
    }
  })
  names(scaled) <- metric_names
  keep_metrics <- metric_names[vapply(scaled, function(x) any(is.finite(x)), logical(1L))]

  grDevices::png(path, width = 1800, height = 1200, res = 170, bg = "#f7f3ea")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 1), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  graphics::plot(df$final_rounds, df$geometry_score,
                 type = "b", pch = 19, lwd = 2, col = "#16324f",
                 xlab = "sampled final_rounds", ylab = "mean rank score",
                 main = "Geometry score vs sampled_f<index>")
  graphics::lines(df$final_rounds, df$full_score,
                  type = "b", pch = 17, lwd = 2, col = "#b23a48")
  graphics::legend("topright",
                   legend = c("geometry score", "full score"),
                   col = c("#16324f", "#b23a48"),
                   pch = c(19, 17),
                   bty = "n")

  graphics::plot(df$final_rounds, df$geometry_score,
                 type = "n",
                 ylim = c(0, 1),
                 xlab = "sampled final_rounds",
                 ylab = "scaled metric rank",
                 main = "Component smoothness diagnostics")
  for (name in keep_metrics) {
    graphics::lines(df$final_rounds, scaled[[name]],
                    type = "b", lwd = 2, pch = 19,
                    col = cols[[name]])
  }
  graphics::legend("topright",
                   legend = labels[keep_metrics],
                   col = cols[keep_metrics],
                   pch = rep(19, length(keep_metrics)),
                   bty = "n")
  graphics::mtext(df$title[[1L]], side = 3, outer = TRUE, cex = 1.2, font = 2)
}

select_panel_ids <- function(best_id) {
  base <- c("default_adaptive", "sampled_f0", "sampled_f4", "sampled_f8",
            "sampled_f16", "sampled_f24", "sampled_f32")
  if (!best_id %in% base) {
    return(c(base, best_id))
  }
  filler <- c("sampled_f12", "sampled_f20", "sampled_f28")
  extra <- filler[!filler %in% base][[1L]]
  c(base, extra)
}

run_panel_layout <- function(spec, panel_id, seed = 1L) {
  n <- max(spec$edges)
  cfg <- if (panel_id == "default_adaptive") {
    default_cfg(n)
  } else {
    final_rounds <- as.integer(sub("sampled_f", "", panel_id, fixed = TRUE))
    sampled_cfg(n, final_rounds)
  }
  layout <- run_layout_safe(spec, cfg, seed)
  if (!identical(layout$status, "ok")) {
    stop(sprintf("Failed to render %s for %s %s: %s",
                 panel_id, spec$family, spec$graph_label, layout$error_message))
  }
  layout$coords
}

write_appendix_tex <- function(path, graph_summary) {
  lines <- c("% Auto-generated appendix for the sampled final-rounds suite.", "")
  for (key in unique(paste(graph_summary$family, graph_summary$graph_label, sep = "::"))) {
    df <- graph_summary[paste(graph_summary$family, graph_summary$graph_label, sep = "::") == key, , drop = FALSE]
    family <- df$family[[1L]]
    label <- df$graph_label[[1L]]
    title <- df$title[[1L]]
    family_dir <- safe_slug(family)
    comp_rel <- sprintf("../tmp/sampled-final-rounds-suite-2026-03-28/png-comparisons/%s/%s-comparison.png",
                        family_dir, safe_slug(label))
    curve_rel <- sprintf("../tmp/sampled-final-rounds-suite-2026-03-28/png-smoothness/%s/%s-smoothness.png",
                         family_dir, safe_slug(label))
      lines <- c(
        lines,
        sprintf("\\subsection{%s}", title),
        "\\begin{figure}[H]",
        "  \\centering",
        sprintf("  \\includegraphics[width=0.92\\textwidth]{%s}", comp_rel),
        sprintf("  \\caption{3x3 comparison sheet for %s. Canonical, adaptive default, sampled\\_f0, sampled\\_f4, sampled\\_f8, sampled\\_f16, sampled\\_f24, sampled\\_f32, and the best sampled configuration by full score.}", title),
        "\\end{figure}",
        "\\begin{figure}[H]",
        "  \\centering",
        sprintf("  \\includegraphics[width=0.92\\textwidth]{%s}", curve_rel),
        sprintf("  \\caption{Smoothness-aware score traces for %s across sampled\\_f0 through sampled\\_f32.}", title),
        "\\end{figure}",
        ""
      )
  }
  writeLines(lines, con = path)
}

panel_subtitle_default <- function(coords, spec, seed = 1L) {
  geom <- geometry_metrics(coords, spec, seed)
  sprintf("adaptive default | rmse=%.4f | sym=%.3f | axis dev=%.3f",
          geom$procrustes.rmse[[1L]],
          geom$global.symmetry.score[[1L]],
          geom$edge.axis.deviation[[1L]])
}

run_sampled_final_rounds_suite <- function(max_graphs = NULL,
                                           family_filter = NULL,
                                           graph_filter = NULL) {
  suite_specs <- all_suite_specs
  if (!is.null(family_filter)) {
    suite_specs <- Filter(function(spec) spec$family %in% family_filter, suite_specs)
  }
  if (!is.null(graph_filter)) {
    suite_specs <- Filter(function(spec) spec$graph_label %in% graph_filter, suite_specs)
  }
  if (!is.null(max_graphs)) {
    max_graphs <- as.integer(max_graphs)
    suite_specs <- suite_specs[seq_len(min(max_graphs, length(suite_specs)))]
  }
  if (length(suite_specs) == 0L) {
    stop("No graphs selected for the sampled final-rounds suite")
  }

  message(sprintf("Running sampled_f0..32 suite on %d graph(s)...", length(suite_specs)))
  raw_rows <- list()
  idx <- 0L
  for (spec in suite_specs) {
    seeds <- seed_budget(max(spec$edges))
    message(sprintf("Graph: %s %s (%d vertices, seeds=%s)",
                    spec$family, spec$graph_label, max(spec$edges), paste(seeds, collapse = ",")))
    for (final_rounds in sampled_finals) {
      for (seed in seeds) {
        idx <- idx + 1L
        raw_rows[[idx]] <- run_one(spec, final_rounds, seed)
      }
    }
  }

  raw_metrics <- do.call(rbind, raw_rows)
  summary_metrics <- summarize_runs(raw_metrics)

  overall_final_summary <- do.call(
    rbind,
    lapply(split(summary_metrics, summary_metrics$final_rounds), function(df) {
      data.frame(
        final_rounds = df$final_rounds[[1L]],
        graphs = nrow(df),
        geometry_score_mean = mean(df$geometry_score),
        full_score_mean = mean(df$full_score),
        global_symmetry_score_mean = mean(df$global_symmetry_score, na.rm = TRUE),
        local_angle_deviation_mean = mean(df$local_angle_deviation, na.rm = TRUE),
        edge_axis_deviation_mean = mean(df$edge_axis_deviation, na.rm = TRUE),
        boundary_waviness_mean = mean(df$boundary_waviness, na.rm = TRUE),
        corridor_waviness_mean = mean(df$corridor_waviness, na.rm = TRUE),
        hole_center_error_mean = mean(df$hole_center_error, na.rm = TRUE),
        elapsed_sec_mean = mean(df$elapsed_sec, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )
  overall_final_summary <- overall_final_summary[order(overall_final_summary$final_rounds), , drop = FALSE]

  best_by_graph <- do.call(
    rbind,
    lapply(split(summary_metrics, paste(summary_metrics$family, summary_metrics$graph_label, sep = "::")), function(df) {
      df <- df[order(df$full_score, df$geometry_score, df$final_rounds), , drop = FALSE]
      best <- df[1L, , drop = FALSE]
      best$best_geometry_final_rounds <- df$final_rounds[[which.min(df$geometry_score)]]
      best
    })
  )

  best_by_family <- do.call(
    rbind,
    lapply(split(best_by_graph, best_by_graph$family), function(df) {
      tab <- sort(table(df$final_rounds), decreasing = TRUE)
      data.frame(
        family = df$family[[1L]],
        modal_best_final_rounds = as.integer(names(tab)[[1L]]),
        graphs = nrow(df),
        stringsAsFactors = FALSE
      )
    })
  )

  global_best_full <- overall_final_summary$final_rounds[[which.min(overall_final_summary$full_score_mean)]]
  global_best_geometry <- overall_final_summary$final_rounds[[which.min(overall_final_summary$geometry_score_mean)]]

  message("Rendering per-graph diagnostics...")
  for (spec in suite_specs) {
    graph_df <- summary_metrics[summary_metrics$family == spec$family &
                                  summary_metrics$graph_label == spec$graph_label, , drop = FALSE]
    graph_df <- graph_df[order(graph_df$final_rounds), , drop = FALSE]
    best_id <- graph_df$config_id[[which.min(graph_df$full_score)]]
    panel_ids <- select_panel_ids(best_id)
    family_dir <- file.path(comparison_root, safe_slug(spec$family))
    dir.create(family_dir, recursive = TRUE, showWarnings = FALSE)
    smooth_dir <- file.path(smoothness_root, safe_slug(spec$family))
    dir.create(smooth_dir, recursive = TRUE, showWarnings = FALSE)

    panel_defs <- list(list(
      coords = project_canonical_for_panel(spec),
      title = sprintf("%s - canonical", spec$title),
      subtitle = "target embedding"
    ))

    for (panel_id in panel_ids) {
      coords <- run_panel_layout(spec, panel_id, seed = 1L)
      subtitle <- if (panel_id == "default_adaptive") {
        panel_subtitle_default(coords, spec, seed = 1L)
      } else {
        summary_row <- graph_df[graph_df$config_id == panel_id, , drop = FALSE]
        sprintf("geom=%.3f | full=%.3f | rmse=%.4f",
                summary_row$geometry_score[[1L]],
                summary_row$full_score[[1L]],
                summary_row$procrustes_rmse[[1L]])
      }
      panel_defs[[length(panel_defs) + 1L]] <- list(
        coords = project_for_panel(coords, spec),
        title = sprintf("%s - %s", spec$title, panel_id),
        subtitle = subtitle
      )
    }

    draw_contact_sheet(
      path = file.path(family_dir, sprintf("%s-comparison.png", safe_slug(spec$graph_label))),
      spec = spec,
      panel_defs = panel_defs
    )
    draw_metric_curve(
      path = file.path(smooth_dir, sprintf("%s-smoothness.png", safe_slug(spec$graph_label))),
      graph_df = graph_df
    )
  }

  raw_csv <- file.path(output_root, "sampled-final-rounds-raw.csv")
  summary_csv <- file.path(output_root, "sampled-final-rounds-summary.csv")
  overall_csv <- file.path(output_root, "sampled-final-rounds-overall-summary.csv")
  best_graph_csv <- file.path(output_root, "sampled-final-rounds-best-by-graph.csv")
  best_family_csv <- file.path(output_root, "sampled-final-rounds-best-by-family.csv")
  appendix_tex <- file.path(output_root, "sampled-final-rounds-suite-appendix.tex")
  summary_md <- file.path(output_root, "sampled-final-rounds-suite-summary.md")

  utils::write.csv(raw_metrics, raw_csv, row.names = FALSE)
  utils::write.csv(summary_metrics, summary_csv, row.names = FALSE)
  utils::write.csv(overall_final_summary, overall_csv, row.names = FALSE)
  utils::write.csv(best_by_graph, best_graph_csv, row.names = FALSE)
  utils::write.csv(best_by_family, best_family_csv, row.names = FALSE)
  write_appendix_tex(appendix_tex, best_by_graph)

  summary_lines <- c(
    "# Sampled Final-Rounds Suite Summary",
    "",
    sprintf("- graphs benchmarked: `%d`", length(suite_specs)),
    "- configs benchmarked: `sampled_f0` through `sampled_f32`",
    sprintf("- global best by geometry score: `sampled_f%d`", global_best_geometry),
    sprintf("- global best by full score: `sampled_f%d`", global_best_full),
    "",
    "## Global summary",
    "",
    "| final_rounds | Mean geometry score | Mean full score | Mean symmetry | Mean angle dev | Mean axis dev | Mean sec |",
    "| ---: | ---: | ---: | ---: | ---: | ---: | ---: |"
  )
  for (i in seq_len(nrow(overall_final_summary))) {
    row <- overall_final_summary[i, , drop = FALSE]
    summary_lines <- c(summary_lines, sprintf(
      "| %d | %.4f | %.4f | %.4f | %.4f | %.4f | %.3f |",
      row$final_rounds,
      row$geometry_score_mean,
      row$full_score_mean,
      row$global_symmetry_score_mean,
      row$local_angle_deviation_mean,
      row$edge_axis_deviation_mean,
      row$elapsed_sec_mean
    ))
  }
  summary_lines <- c(summary_lines, "", "## Best sampled final_rounds by family", "")
  for (i in seq_len(nrow(best_by_family))) {
    row <- best_by_family[i, , drop = FALSE]
    summary_lines <- c(summary_lines, sprintf("- `%s`: modal best `sampled_f%d` across %d graphs",
                                              row$family, row$modal_best_final_rounds, row$graphs))
  }
  summary_lines <- c(summary_lines, "", "## Key outputs",
                     sprintf("- raw metrics: `%s`", raw_csv),
                     sprintf("- graph summary: `%s`", summary_csv),
                     sprintf("- overall summary: `%s`", overall_csv),
                     sprintf("- appendix fragment: `%s`", appendix_tex))
  writeLines(summary_lines, con = summary_md)

  message(sprintf("Summary written to %s", summary_md))
  invisible(list(
    suite_specs = suite_specs,
    raw_metrics = raw_metrics,
    summary_metrics = summary_metrics,
    overall_final_summary = overall_final_summary,
    best_by_graph = best_by_graph,
    best_by_family = best_by_family,
    output_root = output_root
  ))
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  parsed <- if (length(args) == 0L) list() else parse_named_args(args)
  max_graphs <- if (!is.null(parsed$max_graphs)) parse_int_scalar(parsed$max_graphs, "max_graphs") else NULL
  family_filter <- if (!is.null(parsed$family_filter)) parse_char_vector(parsed$family_filter) else NULL
  graph_filter <- if (!is.null(parsed$graph_filter)) parse_char_vector(parsed$graph_filter) else NULL
  run_sampled_final_rounds_suite(
    max_graphs = max_graphs,
    family_filter = family_filter,
    graph_filter = graph_filter
  )
}

if (sys.nframe() == 0L) {
  main()
}
