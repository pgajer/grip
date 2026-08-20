#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) == 0L) {
  stop("Unable to determine script path from commandArgs().")
}
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

run_tag <- "gmds-misf-paraboloid-2026-04-02"
manual_root <- file.path(repo_root, "output", "geodesic_mds_paper")
tmp_dir <- file.path(manual_root, "tmp", run_tag)
pdf_dir <- file.path(manual_root, "reports", run_tag)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

tex_path <- file.path(manual_root, "reports", "gmds_misf_paraboloid_report_2026-04-02.tex")
pdf_path <- file.path(manual_root, "reports", "gmds_misf_paraboloid_report_2026-04-02.pdf")
rds_path <- file.path(tmp_dir, "gmds_misf_paraboloid_results.rds")
metrics_csv <- file.path(tmp_dir, "gmds_misf_paraboloid_metrics.csv")
stage_trace_csv <- file.path(tmp_dir, "gmds_misf_paraboloid_stage_traces.csv")

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo_root, quiet = TRUE, helpers = FALSE)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(repo_root, quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'pkgload', 'devtools', or the 'grip' package to run this benchmark.")
}

ns <- asNamespace("grip")
align_to_target_nd <- get("grip.align.to.target.nd", envir = ns)
classical_mds_embedding <- get("grip.classical.mds.embedding", envir = ns)

fmt_num <- function(x, digits = 4L) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}

fmt_time <- function(x) {
  ifelse(
    is.finite(x),
    ifelse(x < 1, formatC(x, format = "f", digits = 3L), formatC(x, format = "f", digits = 2L)),
    "--"
  )
}

tex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x, perl = TRUE)
  x
}

grid_mesh_triangles <- function(h, w) {
  index <- matrix(seq_len(h * w), nrow = h, ncol = w, byrow = TRUE)
  triangles <- vector("list", 2L * (h - 1L) * (w - 1L))
  k <- 1L
  for (r in seq_len(h - 1L)) {
    for (c in seq_len(w - 1L)) {
      triangles[[k]] <- c(index[r, c], index[r + 1L, c], index[r, c + 1L])
      k <- k + 1L
      triangles[[k]] <- c(index[r + 1L, c], index[r + 1L, c + 1L], index[r, c + 1L])
      k <- k + 1L
    }
  }
  do.call(rbind, triangles)
}

occupied_mesh_triangles <- function(keep) {
  keep <- as.matrix(keep != 0)
  h <- nrow(keep)
  w <- ncol(keep)
  grid <- expand.grid(row = seq_len(h), col = seq_len(w))
  keep_idx <- keep[cbind(grid$row, grid$col)]
  cells <- data.frame(
    row = grid$row[keep_idx],
    col = grid$col[keep_idx]
  )
  id_map <- matrix(0L, nrow = h, ncol = w)
  if (nrow(cells) > 0L) {
    for (i in seq_len(nrow(cells))) {
      id_map[cells$row[[i]], cells$col[[i]]] <- i
    }
  }
  triangles <- list()
  k <- 1L
  if (h >= 2L && w >= 2L) {
    for (r in seq_len(h - 1L)) {
      for (c in seq_len(w - 1L)) {
        block_ids <- c(
          id_map[r, c],
          id_map[r + 1L, c],
          id_map[r, c + 1L],
          id_map[r + 1L, c + 1L]
        )
        if (all(block_ids > 0L)) {
          triangles[[k]] <- c(block_ids[[1L]], block_ids[[2L]], block_ids[[3L]])
          k <- k + 1L
          triangles[[k]] <- c(block_ids[[2L]], block_ids[[4L]], block_ids[[3L]])
          k <- k + 1L
        }
      }
    }
  }
  if (!length(triangles)) {
    return(matrix(integer(), ncol = 3L))
  }
  do.call(rbind, triangles)
}

triangle_areas <- function(coords, triangles) {
  coords <- as.matrix(coords)
  if (!nrow(triangles)) {
    return(numeric(0L))
  }
  if (ncol(coords) == 2L) {
    coords <- cbind(coords, 0)
  }
  v1 <- coords[triangles[, 2L], , drop = FALSE] - coords[triangles[, 1L], , drop = FALSE]
  v2 <- coords[triangles[, 3L], , drop = FALSE] - coords[triangles[, 1L], , drop = FALSE]
  cross <- cbind(
    v1[, 2L] * v2[, 3L] - v1[, 3L] * v2[, 2L],
    v1[, 3L] * v2[, 1L] - v1[, 1L] * v2[, 3L],
    v1[, 1L] * v2[, 2L] - v1[, 2L] * v2[, 1L]
  )
  0.5 * sqrt(rowSums(cross^2))
}

area_floor_ratio <- function(coords, triangles) {
  areas <- triangle_areas(coords, triangles)
  if (!length(areas)) {
    return(NA_real_)
  }
  med <- stats::median(areas)
  if (!is.finite(med) || med <= 0) {
    return(NA_real_)
  }
  as.double(stats::quantile(areas, probs = 0.05, names = FALSE)) / med
}

mesh_roughness <- function(coords, adj_list, edges) {
  centered <- sweep(coords, 2L, colMeans(coords), FUN = "-", check.margin = FALSE)
  median.edge <- stats::median(sqrt(rowSums(
    (centered[edges[, 1L], , drop = FALSE] - centered[edges[, 2L], , drop = FALSE])^2
  )))
  if (!is.finite(median.edge) || median.edge <= 0) {
    return(NA_real_)
  }
  residuals <- vapply(seq_len(nrow(centered)), function(i) {
    nbrs <- adj_list[[i]]
    if (length(nbrs) == 0L) {
      return(0)
    }
    delta <- centered[i, ] - colMeans(centered[nbrs, , drop = FALSE])
    sum(delta^2)
  }, numeric(1L))
  sqrt(mean(residuals)) / median.edge
}

case_family_label <- function(family) {
  if (identical(family, "regular")) "Regular orthogonal mesh" else "Irregular occupied mesh"
}

make_regular_case <- function(side, amplitude = 0.35) {
  bundle <- mesh.surface.graph(
    side,
    side,
    surface = "paraboloid",
    amplitude = amplitude,
    connectivity = "orthogonal",
    normalize = "median"
  )
  prepared <- grip.prepare.graph.geodesic.mds(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average"
  )
  started <- proc.time()[["elapsed"]]
  cmd <- classical_mds_embedding(prepared, dim = 3L, eig = TRUE)
  cmd_elapsed <- proc.time()[["elapsed"]] - started
  list(
    id = sprintf("paraboloid_regular_%dx%d", side, side),
    label = sprintf("Paraboloid regular mesh %dx%d", side, side),
    family = "regular",
    side = as.integer(side),
    n = bundle$n,
    edges = bundle$edges,
    keep = NULL,
    triangles = grid_mesh_triangles(side, side),
    truth = bundle$coords_surface,
    prepared = prepared,
    adj_list = prepared$adj_list,
    cmd = cmd,
    cmd_elapsed = as.double(cmd_elapsed)
  )
}

make_irregular_case <- function(side, amplitude = 0.35) {
  keep <- keep.asymmetric.notches(
    side,
    side,
    notch_depth = max(2L, as.integer(round(side / 4))),
    notch_width = max(1L, as.integer(round(side / 8)))
  )
  bundle <- occupied.mesh.surface.graph(
    keep = keep,
    surface = "paraboloid",
    amplitude = amplitude,
    connectivity = "orthogonal",
    normalize = "median"
  )
  prepared <- grip.prepare.graph.geodesic.mds(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average"
  )
  started <- proc.time()[["elapsed"]]
  cmd <- classical_mds_embedding(prepared, dim = 3L, eig = TRUE)
  cmd_elapsed <- proc.time()[["elapsed"]] - started
  list(
    id = sprintf("paraboloid_irregular_%dx%d", side, side),
    label = sprintf("Paraboloid irregular occupied mesh %dx%d", side, side),
    family = "irregular",
    side = as.integer(side),
    n = bundle$n,
    edges = bundle$edges,
    keep = keep,
    triangles = occupied_mesh_triangles(keep),
    truth = bundle$coords_surface,
    prepared = prepared,
    adj_list = prepared$adj_list,
    cmd = cmd,
    cmd_elapsed = as.double(cmd_elapsed)
  )
}

compute_metrics <- function(case,
                            spec,
                            coords,
                            elapsed_sec = NA_real_,
                            trace = NULL,
                            stage_trace = NULL) {
  score_df <- grip.score.geodesic.mds(coords = coords, prepared = case$prepared)
  aligned <- if (identical(spec$id, "reference")) {
    list(aligned = as.matrix(coords), rmse = 0)
  } else {
    align_to_target_nd(coords, case$truth, allow.reflection = TRUE)
  }
  data.frame(
    case_id = case$id,
    case_label = case$label,
    family = case$family,
    family_label = case_family_label(case$family),
    side = case$side,
    n = case$n,
    method_id = spec$id,
    method_label = spec$label,
    short_label = spec$short,
    elapsed_sec = as.double(elapsed_sec),
    gmds_energy = score_df$gmds.energy[[1L]],
    gmds_stress = score_df$gmds.stress[[1L]],
    gmds_raw_stress = score_df$gmds.raw_stress[[1L]],
    procrustes_rmse = aligned$rmse,
    roughness = mesh_roughness(coords, case$adj_list, case$edges),
    area_q05_ratio = area_floor_ratio(coords, case$triangles),
    iterations = if (is.null(trace) || nrow(trace) == 0L) 0L else max(trace$iteration),
    stage_rows = if (is.null(stage_trace) || nrow(stage_trace) == 0L) 0L else nrow(stage_trace),
    stringsAsFactors = FALSE
  )
}

method_specs <- list(
  list(id = "reference", label = "Reference surface", short = "Ref", kind = "reference"),
  list(id = "cmdscale", label = "cmdscale", short = "CMD", kind = "cmdscale"),
  list(
    id = "cmd_pure_gmds",
    label = "cmdscale -> pure GMDS",
    short = "CMD->GMDS",
    kind = "pure_gmds",
    args = list(engine = "cpp", max_iter = 8L, return_trace = TRUE, n_threads = 0L)
  ),
  list(
    id = "misf_gmds",
    label = "MISF-GMDS",
    short = "MISF",
    kind = "misf",
    args = list(
      dim = 3L,
      top_level_restarts = 2L,
      top_level_max_iter = 6L,
      top_level_engine = "cpp",
      insertion_anchor_policy = "prev_level_spread",
      insertion_max_iter = 20L,
      refinement_local_nbrs = 4L,
      refinement_landmark_count = 2L,
      refinement_pair_mode = "sparse",
      refinement_anchor_weight = 0.05,
      refinement_anchor_weight_end = 0.01,
      refinement_continuation = "linear",
      refinement_max_iter = 3L,
      refinement_engine = "cpp",
      final_polish_max_iter = 0L,
      final_polish_engine = "cpp",
      n_threads = 0L,
      return_trace = TRUE,
      return_frames = FALSE
    )
  ),
  list(
    id = "misf_gmds_polish",
    label = "MISF-GMDS + final polish",
    short = "MISF+P",
    kind = "misf",
    args = list(
      dim = 3L,
      top_level_restarts = 2L,
      top_level_max_iter = 6L,
      top_level_engine = "cpp",
      insertion_anchor_policy = "prev_level_spread",
      insertion_max_iter = 20L,
      refinement_local_nbrs = 4L,
      refinement_landmark_count = 2L,
      refinement_pair_mode = "sparse",
      refinement_anchor_weight = 0.05,
      refinement_anchor_weight_end = 0.01,
      refinement_continuation = "linear",
      refinement_max_iter = 3L,
      refinement_engine = "cpp",
      final_polish_max_iter = 4L,
      final_polish_engine = "cpp",
      n_threads = 0L,
      return_trace = TRUE,
      return_frames = FALSE
    )
  )
)

run_method <- function(case, spec) {
  if (identical(spec$kind, "reference")) {
    row <- compute_metrics(case = case, spec = spec, coords = case$truth)
    return(list(
      coords = case$truth,
      display_coords = case$truth,
      metrics = row,
      trace = NULL,
      stage_trace = NULL,
      fit = NULL
    ))
  }

  if (identical(spec$kind, "cmdscale")) {
    coords <- case$cmd$coords
    row <- compute_metrics(case = case, spec = spec, coords = coords, elapsed_sec = case$cmd_elapsed)
    return(list(
      coords = coords,
      display_coords = align_to_target_nd(coords, case$truth, allow.reflection = TRUE)$aligned,
      metrics = row,
      trace = NULL,
      stage_trace = NULL,
      fit = NULL
    ))
  }

  if (identical(spec$kind, "pure_gmds")) {
    started <- proc.time()[["elapsed"]]
    fit <- do.call(
      grip.optimize.geodesic.mds,
      c(list(coords = case$cmd$coords, prepared = case$prepared), spec$args)
    )
    elapsed <- proc.time()[["elapsed"]] - started
    row <- compute_metrics(
      case = case,
      spec = spec,
      coords = fit$coords,
      elapsed_sec = elapsed,
      trace = fit$trace,
      stage_trace = NULL
    )
    return(list(
      coords = fit$coords,
      display_coords = align_to_target_nd(fit$coords, case$truth, allow.reflection = TRUE)$aligned,
      metrics = row,
      trace = fit$trace,
      stage_trace = NULL,
      fit = fit
    ))
  }

  fit <- do.call(
    grip.optimize.misf.geodesic.mds,
    c(list(prepared = case$prepared, seed = case$side + if (identical(case$family, "irregular")) 100L else 0L), spec$args)
  )
  row <- compute_metrics(
    case = case,
    spec = spec,
    coords = fit$coords,
    elapsed_sec = fit$timing$total,
    trace = if (!is.null(fit$final_polish)) fit$final_polish$trace else NULL,
    stage_trace = fit$stage_trace
  )
  list(
    coords = fit$coords,
    display_coords = align_to_target_nd(fit$coords, case$truth, allow.reflection = TRUE)$aligned,
    metrics = row,
    trace = if (!is.null(fit$final_polish)) fit$final_polish$trace else NULL,
    stage_trace = fit$stage_trace,
    fit = fit
  )
}

build_panel_title <- function(method_result) {
  row <- method_result$metrics[1L, , drop = FALSE]
  if (identical(row$method_id[[1L]], "reference")) {
    return("Reference surface")
  }
  sprintf(
    "%s\nsigma %s, rho %s\nt %ss",
    row$method_label[[1L]],
    fmt_num(row$gmds_stress[[1L]], 4L),
    fmt_num(row$procrustes_rmse[[1L]], 4L),
    fmt_time(row$elapsed_sec[[1L]])
  )
}

save_case_panel_grid <- function(case_result, output_path) {
  methods <- case_result$methods
  n_panels <- length(methods)
  grDevices::png(output_path, width = 2400L, height = 1600L, res = 180, bg = "#ffffff")
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(2L, 3L), mar = c(1.2, 1.2, 2.9, 0.4), oma = c(0, 0, 1.2, 0))

  for (i in seq_len(6L)) {
    if (i > n_panels) {
      graphics::plot.new()
      next
    }
    method <- methods[[i]]
    grip.plot(
      coords = method$display_coords,
      edges = case_result$case$edges,
      projection = "ortho",
      azimuth = 35,
      elevation = 24,
      vertex.col = if (identical(method$metrics$method_id[[1L]], "reference")) "#bc6c25" else "#3a5a40",
      edge.col = "#adb5bd",
      main = ""
    )
    graphics::mtext(build_panel_title(method), side = 3L, line = 0.3, cex = 0.80)
  }
  graphics::mtext(case_result$case$label, side = 3L, outer = TRUE, line = -0.3, cex = 1.1, font = 2L)
}

stage_trace_rows <- function(case_result) {
  rows <- list()
  for (method in case_result$methods) {
    if (is.null(method$stage_trace) || !nrow(method$stage_trace)) {
      next
    }
    trace_df <- method$stage_trace
    trace_df <- trace_df[is.finite(trace_df$energy), , drop = FALSE]
    if (!nrow(trace_df)) {
      next
    }
    trace_df$case_id <- case_result$case$id
    trace_df$case_label <- case_result$case$label
    trace_df$family <- case_result$case$family
    trace_df$method_id <- method$metrics$method_id[[1L]]
    trace_df$method_label <- method$metrics$method_label[[1L]]
    trace_df$stage_index <- seq_len(nrow(trace_df))
    trace_df$stage_label <- ifelse(
      trace_df$stage == "refinement",
      sprintf("ref L%d", trace_df$level),
      ifelse(trace_df$stage == "top_level", "top", "final")
    )
    rows[[length(rows) + 1L]] <- trace_df
  }
  if (!length(rows)) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

save_stage_trace_grid <- function(stage_df, cases, output_path) {
  case_ids <- vapply(cases, `[[`, character(1L), "id")
  grDevices::png(output_path, width = 2400L, height = 1500L, res = 180, bg = "#ffffff")
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(2L, 3L), mar = c(3.5, 3.7, 2.2, 0.8), oma = c(0, 0, 1.0, 0))

  cols <- c("misf_gmds" = "#3a5a40", "misf_gmds_polish" = "#9c6644")
  for (case_id in case_ids) {
    df <- stage_df[stage_df$case_id == case_id, , drop = FALSE]
    if (!nrow(df)) {
      graphics::plot.new()
      next
    }
    yvals <- log10(pmax(df$energy, 1e-12))
    xlim <- c(1, max(df$stage_index))
    ylim <- range(yvals)
    graphics::plot(NA,
      xlim = xlim,
      ylim = ylim,
      xlab = "Stage step",
      ylab = "log10 energy",
      main = unique(df$case_label),
      xaxt = "n"
    )
    for (method_id in c("misf_gmds", "misf_gmds_polish")) {
      sub <- df[df$method_id == method_id, , drop = FALSE]
      if (!nrow(sub)) {
        next
      }
      graphics::lines(
        x = sub$stage_index,
        y = log10(pmax(sub$energy, 1e-12)),
        type = "b",
        pch = 16,
        lwd = 2,
        col = cols[[method_id]]
      )
      graphics::text(
        x = sub$stage_index,
        y = log10(pmax(sub$energy, 1e-12)),
        labels = sub$stage_label,
        pos = 3,
        cex = 0.65,
        col = cols[[method_id]]
      )
    }
    graphics::axis(1, at = sort(unique(df$stage_index)))
    graphics::legend(
      "topright",
      legend = c("MISF-GMDS", "MISF-GMDS + polish"),
      col = cols,
      lty = 1,
      lwd = 2,
      pch = 16,
      bty = "n",
      cex = 0.8
    )
  }
  graphics::mtext("Per-level MISF energy traces", side = 3L, outer = TRUE, line = -0.3, cex = 1.1, font = 2L)
}

write_case_table <- function(df) {
  rows <- vapply(seq_len(nrow(df)), function(i) {
    sprintf(
      "%s & %s & %s & %s & %s & %s & %s \\\\",
      tex_escape(df$method_label[[i]]),
      fmt_time(df$elapsed_sec[[i]]),
      fmt_num(df$gmds_stress[[i]], 4L),
      fmt_num(df$procrustes_rmse[[i]], 4L),
      fmt_num(df$roughness[[i]], 4L),
      fmt_num(df$area_q05_ratio[[i]], 4L),
      ifelse(df$stage_rows[[i]] > 0L, as.character(df$stage_rows[[i]]), "--")
    )
  }, character(1L))
  paste(
    "\\begin{table}[H]",
    "\\centering",
    "\\scriptsize",
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{lrrrrrr}",
    "\\toprule",
    "Method & $t$ (s) & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ & trace rows \\\\",
    "\\midrule",
    paste(rows, collapse = "\n"),
    "\\bottomrule",
    "\\end{tabular}}",
    "\\end{table}",
    sep = "\n"
  )
}

case_summary_paragraph <- function(df) {
  baseline <- df[df$method_id == "cmd_pure_gmds", , drop = FALSE]
  misf <- df[df$method_id == "misf_gmds", , drop = FALSE]
  polished <- df[df$method_id == "misf_gmds_polish", , drop = FALSE]
  sprintf(
    "On %s, the direct cmdscale initializer followed by pure GMDS lands at $\\sigma=%s$ and $\\rho=%s$. The multiscale MISF run without a final global polish reaches $\\sigma=%s$ and $\\rho=%s$, while the polished MISF variant finishes at $\\sigma=%s$ and $\\rho=%s$. This is the basin-selection comparison that motivated the MISF rewrite: the question is whether the coarse-to-fine path gets us into a smoother 3D basin than the global spectral start.",
    tex_escape(df$case_label[[1L]]),
    fmt_num(baseline$gmds_stress[[1L]], 4L),
    fmt_num(baseline$procrustes_rmse[[1L]], 4L),
    fmt_num(misf$gmds_stress[[1L]], 4L),
    fmt_num(misf$procrustes_rmse[[1L]], 4L),
    fmt_num(polished$gmds_stress[[1L]], 4L),
    fmt_num(polished$procrustes_rmse[[1L]], 4L)
  )
}

case_figure_rel <- function(case_id) {
  file.path(pdf_dir, sprintf("%s_grid.png", case_id))
}

stage_trace_figure_rel <- file.path(pdf_dir, "gmds_misf_paraboloid_stage_trace_grid.png")

cases <- c(
  lapply(c(12L, 15L, 20L), make_regular_case),
  lapply(c(12L, 15L, 20L), make_irregular_case)
)

case_results <- lapply(cases, function(case) {
  methods <- lapply(method_specs, run_method, case = case)
  metrics <- do.call(rbind, lapply(methods, `[[`, "metrics"))
  out <- list(case = case, methods = methods, metrics = metrics)
  save_case_panel_grid(out, file.path(pdf_dir, sprintf("%s_grid.png", case$id)))
  out
})

compact_case_result <- function(case_result) {
  list(
    case = list(
      id = case_result$case$id,
      label = case_result$case$label,
      family = case_result$case$family,
      side = case_result$case$side,
      n = case_result$case$n,
      edges = case_result$case$edges
    ),
    methods = lapply(case_result$methods, function(method) {
      list(
        display_coords = method$display_coords,
        metrics = method$metrics
      )
    }),
    metrics = case_result$metrics
  )
}

bundle_case_results <- lapply(case_results, compact_case_result)

metrics_df <- do.call(rbind, lapply(case_results, `[[`, "metrics"))
stage_trace_df <- do.call(
  rbind,
  Filter(function(x) nrow(x) > 0L, lapply(case_results, stage_trace_rows))
)

if (nrow(stage_trace_df) > 0L) {
  save_stage_trace_grid(stage_trace_df, cases = lapply(case_results, `[[`, "case"), output_path = stage_trace_figure_rel)
}

utils::write.csv(metrics_df, metrics_csv, row.names = FALSE)
utils::write.csv(stage_trace_df, stage_trace_csv, row.names = FALSE)
saveRDS(
  list(
    run_tag = run_tag,
    case_results = bundle_case_results,
    metrics = metrics_df,
    stage_traces = stage_trace_df
  ),
  rds_path
)

regular_summary <- subset(metrics_df, family == "regular" & method_id %in% c("cmd_pure_gmds", "misf_gmds", "misf_gmds_polish"))
irregular_summary <- subset(metrics_df, family == "irregular" & method_id %in% c("cmd_pure_gmds", "misf_gmds", "misf_gmds_polish"))

count_better <- function(df, method_a, method_b, field, smaller_is_better = TRUE) {
  wide <- reshape(
    df[, c("case_id", "method_id", field)],
    idvar = "case_id",
    timevar = "method_id",
    direction = "wide"
  )
  a <- wide[[sprintf("%s.%s", field, method_a)]]
  b <- wide[[sprintf("%s.%s", field, method_b)]]
  if (smaller_is_better) {
    sum(a < b, na.rm = TRUE)
  } else {
    sum(a > b, na.rm = TRUE)
  }
}

regular_better_rho <- count_better(regular_summary, "misf_gmds_polish", "cmd_pure_gmds", "procrustes_rmse", smaller_is_better = TRUE)
regular_better_sigma <- count_better(regular_summary, "misf_gmds_polish", "cmd_pure_gmds", "gmds_stress", smaller_is_better = TRUE)
irregular_better_rho <- count_better(irregular_summary, "misf_gmds_polish", "cmd_pure_gmds", "procrustes_rmse", smaller_is_better = TRUE)
irregular_better_sigma <- count_better(irregular_summary, "misf_gmds_polish", "cmd_pure_gmds", "gmds_stress", smaller_is_better = TRUE)

regular_sections <- vapply(case_results[vapply(case_results, function(x) identical(x$case$family, "regular"), logical(1L))], function(case_result) {
  df <- subset(case_result$metrics, method_id != "reference")
  paste(
    sprintf("\\subsection{%s}", tex_escape(case_result$case$label)),
    case_summary_paragraph(df),
    "",
    sprintf(
      "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{%s. The panel set keeps the comparison narrow: classical MDS, cmdscale-initialized GMDS, MISF-GMDS without final polish, and MISF-GMDS with the final full-graph polish.}\\label{fig:%s}\\end{figure}",
      case_figure_rel(case_result$case$id),
      tex_escape(case_result$case$label),
      case_result$case$id
    ),
    write_case_table(df),
    sep = "\n\n"
  )
}, character(1L))

irregular_sections <- vapply(case_results[vapply(case_results, function(x) identical(x$case$family, "irregular"), logical(1L))], function(case_result) {
  df <- subset(case_result$metrics, method_id != "reference")
  paste(
    sprintf("\\subsection{%s}", tex_escape(case_result$case$label)),
    case_summary_paragraph(df),
    "",
    sprintf(
      "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{%s. The same method family is rerun on an irregular occupied mesh to test whether MISF still helps once the global symmetry of the rectangular lattice is broken.}\\label{fig:%s}\\end{figure}",
      case_figure_rel(case_result$case$id),
      tex_escape(case_result$case$label),
      case_result$case$id
    ),
    write_case_table(df),
    sep = "\n\n"
  )
}, character(1L))

tex_lines <- c(
  "\\documentclass[11pt]{article}",
  "\\usepackage[margin=1in]{geometry}",
  "\\usepackage{graphicx}",
  "\\usepackage{booktabs}",
  "\\usepackage{float}",
  "\\usepackage{amsmath}",
  "\\usepackage[T1]{fontenc}",
  "\\usepackage[utf8]{inputenc}",
  "\\title{MISF-Based GMDS Paraboloid Validation Slice}",
  "\\author{Phase 6 of the MISF-GMDS implementation plan}",
  "\\date{2026-04-02}",
  "\\begin{document}",
  "\\maketitle",
  "\\section{Purpose}",
  "This report implements Phase 6 of the MISF-based GMDS plan: the first dedicated validation slice focused on the paraboloid families that originally exposed the basin problem. The central comparison is no longer between different regularizers, but between two initialization strategies for the same geodesic objective: a global classical-MDS start followed by pure GMDS, and a multiscale MISF start followed by sparse levelwise GMDS refinement.",
  "\\section{Benchmark design}",
  "The test set has regular orthogonal paraboloids and irregular occupied-mesh paraboloids at nominal side lengths $12$, $15$, and $20$. All graph weights come from the embedded paraboloid edge lengths after the helper-function normalization already used throughout the earlier experiments. The reported metrics are geodesic stress $\\sigma$, Procrustes RMSE $\\rho$, roughness $\\eta$, the area-floor ratio $\\alpha_{0.05}$, and wall-clock embedding time $t$ in seconds.",
  "\\section{Regular paraboloids}",
  sprintf(
    "Across the three regular cases, the polished MISF run improves $\\rho$ over the cmdscale baseline in %d of 3 cases and improves $\\sigma$ in %d of 3 cases. This is the direct basin test on the symmetric lattice that previously caused trouble.",
    regular_better_rho,
    regular_better_sigma
  ),
  paste(regular_sections, collapse = "\n\n"),
  "\\section{Irregular occupied paraboloids}",
  sprintf(
    "Across the three irregular cases, the polished MISF run improves $\\rho$ over the cmdscale baseline in %d of 3 cases and improves $\\sigma$ in %d of 3 cases. These runs test whether MISF still matters once the occupied mesh already breaks some of the lattice symmetry.",
    irregular_better_rho,
    irregular_better_sigma
  ),
  paste(irregular_sections, collapse = "\n\n"),
  "\\section{Per-level MISF energy traces}",
  "Figure~\\ref{fig:misf-stage-traces} condenses the multiscale traces into per-case energy ladders. Only stages with directly comparable GMDS energies are shown: the coarse top-level solve, each level-refinement endpoint, and the final polish. This makes it easier to see whether the multiscale pipeline descends smoothly into the full-graph objective rather than jumping immediately into the bad basin.",
  sprintf(
    "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{Per-case MISF energy traces for the no-polish and polished variants. The top-level energy and each level-refinement endpoint are plotted on a log scale.}\\label{fig:misf-stage-traces}\\end{figure}",
    stage_trace_figure_rel
  ),
  "\\section{Takeaways}",
  "This Phase 6 slice answers the implementation question the earlier experiments left open. The MISF-based initializer is now evaluated on exactly the families that motivated it, and it can be compared directly against the older cmdscale-initialized pure-GMDS path using the same geodesic score. The remaining question after this report is no longer whether MISF can be implemented cleanly in \\texttt{grip}; it is which graph families most clearly benefit from the new basin-selection strategy and how much final global polish should be kept by default.",
  "\\end{document}"
)

writeLines(tex_lines, tex_path)

message("Wrote benchmark bundle: ", rds_path)
message("Wrote metrics CSV: ", metrics_csv)
message("Wrote stage traces CSV: ", stage_trace_csv)
message("Wrote LaTeX report: ", tex_path)
message("Expected PDF path after latexmk: ", pdf_path)
