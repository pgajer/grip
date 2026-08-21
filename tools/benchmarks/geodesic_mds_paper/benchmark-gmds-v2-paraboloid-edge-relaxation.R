#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) == 0L) {
  stop("Unable to determine script path from commandArgs().")
}
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

manual_root <- file.path(repo_root, "output", "geodesic_mds_paper")
run_tag <- "gmds-v2-paraboloid-regularized-2026-04-01"
tmp_dir <- file.path(manual_root, "tmp", run_tag)
pdf_dir <- file.path(manual_root, "reports", run_tag)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

appendix_tex_path <- file.path(
  manual_root,
  "pdf",
  "gmds_v2_paraboloid_edge_relaxation_appendix_2026-04-01.tex"
)
metrics_csv <- file.path(tmp_dir, "gmds_v2_paraboloid_edge_relaxation_metrics.csv")
rds_path <- file.path(tmp_dir, "gmds_v2_paraboloid_edge_relaxation_results.rds")

baseline_csv <- file.path(tmp_dir, "gmds_v2_paraboloid_regularized_metrics.csv")

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

triangle_areas <- function(coords, triangles) {
  coords <- as.matrix(coords)
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

make_case <- function(side, amplitude = 0.35) {
  bundle <- mesh.surface.graph(
    side,
    side,
    surface = "paraboloid",
    amplitude = amplitude,
    connectivity = "orthogonal",
    normalize = "median"
  )
  prepared <- prepare.graph.geodesic.mds(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average"
  )
  cmd <- classical_mds_embedding(prepared, dim = 3L, eig = TRUE)
  list(
    id = sprintf("paraboloid_%dx%d", side, side),
    label = sprintf("Orthogonal paraboloid mesh %dx%d", side, side),
    side = side,
    truth = bundle$coords_surface,
    edges = bundle$edges,
    adj_list = grip.build.adj.from.edges(bundle$edges, n = bundle$n)$adj_list,
    triangles = grid_mesh_triangles(side, side),
    prepared_full = prepared,
    cmd = cmd
  )
}

make_surrogate_prepared <- function(prepared_full) {
  prepared <- prepared_full
  prepared$pair_matrix <- matrix(integer(), ncol = 2L)
  prepared$pair_graph_distance <- numeric(0L)
  prepared$path_vertices <- list()
  prepared$path_edges <- list()
  prepared$path_edge_weights <- list()
  prepared$flat_pair_edge_offsets <- as.integer(0L)
  prepared$flat_edge_u <- integer(0L)
  prepared$flat_edge_v <- integer(0L)
  prepared$flat_edge_coeff <- numeric(0L)
  prepared
}

schedule_linear <- function(start, end, max_iter) {
  seq(from = as.double(start), to = as.double(end), length.out = max_iter + 1L)
}

extract_coords_by_trace <- function(trace_df, frames) {
  if (length(frames) == 0L || nrow(trace_df) == 0L) {
    return(list())
  }
  coords_seq <- vector("list", nrow(trace_df))
  frame_idx <- 1L
  coords_seq[[1L]] <- frames[[frame_idx]]
  if (nrow(trace_df) >= 2L) {
    for (i in 2L:nrow(trace_df)) {
      accepted <- isTRUE(trace_df$accepted[[i]])
      if (accepted && frame_idx < length(frames)) {
        frame_idx <- frame_idx + 1L
      }
      coords_seq[[i]] <- frames[[frame_idx]]
    }
  }
  coords_seq
}

select_snapshot_indices <- function(n_frames) {
  raw <- c(1L, 2L, 4L, 7L, 11L, n_frames)
  raw <- raw[raw <= n_frames]
  unique(raw)
}

save_snapshot_grid <- function(case, trajectory_df, method_label, output_path, lambda_cols) {
  idx <- select_snapshot_indices(nrow(trajectory_df))
  grDevices::png(output_path, width = 1400L * 2L, height = 900L * 2L, res = 180, bg = "#ffffff")
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(2L, 3L), mar = c(1.2, 1.2, 3.0, 0.4), oma = c(0, 0, 1.0, 0))

  for (i in idx) {
    row <- trajectory_df[i, , drop = FALSE]
    plot.layout(
      coords = row$display_coords[[1L]],
      edges = case$edges,
      projection = "ortho",
      azimuth = 35,
      elevation = 24,
      vertex.col = "#3a5a40",
      edge.col = "#adb5bd",
      main = ""
    )
    lambda_bits <- vapply(lambda_cols, function(col) {
      sprintf("%s %.2f", col, row[[col]][[1L]])
    }, character(1L))
    title <- sprintf(
      "iter %d, %s\nE_surr %s, sigma %s",
      row$iteration[[1L]],
      paste(lambda_bits, collapse = ", "),
      fmt_num(row$surrogate_energy[[1L]], 4L),
      fmt_num(row$geodesic_stress[[1L]], 4L)
    )
    graphics::mtext(title, side = 3L, line = 0.3, cex = 0.82)
  }
  graphics::mtext(
    sprintf("%s: %s", case$label, method_label),
    side = 3L,
    outer = TRUE,
    line = -0.3,
    cex = 1.15,
    font = 2L
  )
}

run_case_method <- function(case, method_spec) {
  prepared_surrogate <- make_surrogate_prepared(case$prepared_full)
  prepared_surrogate <- grip.geodesic.mds.ensure.graph.term.cache(
    prepared_surrogate,
    repulsion_weight = max(method_spec$repulsion_schedule),
    repulsion_quantile = 0.60,
    repulsion_scale = 0.20,
    repulsion_cap_quantile = 0.90,
    repulsion_hop_min = 3L
  )

  graph_edge_u <- as.integer(prepared_surrogate$graph_edge_matrix[, 1L] - 1L)
  graph_edge_v <- as.integer(prepared_surrogate$graph_edge_matrix[, 2L] - 1L)
  graph_edge_target <- as.double(prepared_surrogate$graph_edge_target)
  repulsion_u <- if (!is.null(prepared_surrogate$flat_repulsion_u)) as.integer(prepared_surrogate$flat_repulsion_u) else integer(0L)
  repulsion_v <- if (!is.null(prepared_surrogate$flat_repulsion_v)) as.integer(prepared_surrogate$flat_repulsion_v) else integer(0L)
  repulsion_target <- if (!is.null(prepared_surrogate$flat_repulsion_target)) as.double(prepared_surrogate$flat_repulsion_target) else numeric(0L)

  started <- proc.time()[["elapsed"]]
  fit <- grip_optimize_geodesic_mds_flat_cpp(
    flat_pair_edge_offsets = as.integer(0L),
    flat_edge_u = integer(0L),
    flat_edge_v = integer(0L),
    flat_edge_coeff = numeric(0L),
    pair_graph_distance = numeric(0L),
    coords = case$cmd$coords,
    max_iter = 15L,
    edge_length_epsilon = 1e-8,
    initial_step = 1.0,
    step_shrink = 0.5,
    armijo_factor = 1e-4,
    grad_tol = 1e-8,
    min_step = 1e-8,
    recenter = TRUE,
    return_trace = TRUE,
    anchor_coords = NULL,
    anchor_weights = NULL,
    smooth_adj_offsets = integer(0L),
    smooth_adj_vertices = integer(0L),
    smooth_weights = NULL,
    graph_edge_u = graph_edge_u,
    graph_edge_v = graph_edge_v,
    graph_edge_target = graph_edge_target,
    edge_spring_weights = method_spec$edge_schedule,
    repulsion_u = repulsion_u,
    repulsion_v = repulsion_v,
    repulsion_target = repulsion_target,
    repulsion_weights = method_spec$repulsion_schedule,
    n_threads = 0L
  )
  elapsed <- proc.time()[["elapsed"]] - started

  coords_by_trace <- extract_coords_by_trace(fit$trace, fit$frames)
  rows <- lapply(seq_len(nrow(fit$trace)), function(i) {
    coords <- coords_by_trace[[i]]
    surrogate_score <- grip.score.geodesic.mds(
      coords = coords,
      prepared = prepared_surrogate,
      edge_spring_weight = fit$trace$edge_spring_weight[[i]],
      repulsion_weight = fit$trace$repulsion_weight[[i]]
    )
    geo_score <- grip.score.geodesic.mds(
      coords = coords,
      prepared = case$prepared_full
    )
    aligned <- align_to_target_nd(coords, case$truth, allow.reflection = TRUE)
    data.frame(
      case_id = case$id,
      case_label = case$label,
      method_id = method_spec$id,
      method_label = method_spec$label,
      iteration = fit$trace$iteration[[i]],
      lambda_E = fit$trace$edge_spring_weight[[i]],
      lambda_R = fit$trace$repulsion_weight[[i]],
      surrogate_energy = surrogate_score$gmds.energy[[1L]],
      edge_energy = surrogate_score$edge.spring.energy[[1L]],
      repulsion_energy = surrogate_score$repulsion.energy[[1L]],
      geodesic_stress = geo_score$gmds.stress[[1L]],
      geodesic_raw_stress = geo_score$gmds.raw_stress[[1L]],
      procrustes_rmse = aligned$rmse,
      roughness = mesh_roughness(coords, case$adj_list, case$edges),
      area_q05_ratio = area_floor_ratio(coords, case$triangles),
      accepted = fit$trace$accepted[[i]],
      elapsed_total_sec = elapsed,
      stringsAsFactors = FALSE
    )
  })
  trajectory_df <- do.call(rbind, rows)
  trajectory_df$display_coords <- I(lapply(coords_by_trace, function(coords) {
    align_to_target_nd(coords, case$truth, allow.reflection = TRUE)$aligned
  }))

  figure_path <- file.path(pdf_dir, sprintf("%s_%s_surrogate_snapshots.png", case$id, method_spec$id))
  save_snapshot_grid(case, trajectory_df, method_spec$label, figure_path, method_spec$lambda_cols)

  list(
    case = case,
    method = method_spec,
    fit = fit,
    trajectory = trajectory_df,
    figure_path = figure_path
  )
}

write_case_table <- function(df, method_spec) {
  rows <- vapply(seq_len(nrow(df)), function(i) {
    sprintf(
      "%d & %s & %s & %s & %s & %s & %s & %s \\\\",
      df$iteration[[i]],
      fmt_num(df$lambda_E[[i]], 4L),
      fmt_num(df$lambda_R[[i]], 4L),
      fmt_num(df$surrogate_energy[[i]], 4L),
      fmt_num(df$geodesic_stress[[i]], 4L),
      fmt_num(df$procrustes_rmse[[i]], 4L),
      fmt_num(df$roughness[[i]], 4L),
      fmt_num(df$area_q05_ratio[[i]], 4L)
    )
  }, character(1L))
  paste(
    "\\begin{table}[H]",
    "\\centering",
    "\\scriptsize",
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{rrrrrrrr}",
    "\\toprule",
    "iter & $\\lambda_E$ & $\\lambda_R$ & $E_{\\mathrm{surr}}$ & $\\sigma_{\\mathrm{geo}}$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ \\\\",
    "\\midrule",
    paste(rows, collapse = "\n"),
    "\\bottomrule",
    "\\end{tabular}}",
    sprintf(
      "\\caption{%s trajectory under the stronger linear schedules. Here $E_{\\mathrm{surr}}$ is the surrogate objective, while $\\sigma_{\\mathrm{geo}}$ is the full fixed-path geodesic stress measured afterward on the same embedding. Lower $E_{\\mathrm{surr}}$, $\\sigma_{\\mathrm{geo}}$, $\\rho$, and $\\eta$ are better; larger $\\alpha_{0.05}$ is better.}",
      method_spec$label
    ),
    "\\end{table}",
    sep = "\n"
  )
}

summarize_case <- function(df) {
  best_sur_idx <- which.min(df$surrogate_energy)
  best_geo_idx <- which.min(df$geodesic_stress)
  best_alpha_idx <- which.max(df$area_q05_ratio)
  sprintf(
    "For this case, the lowest surrogate energy occurs at iteration %d ($\\lambda_E=%s$, $\\lambda_R=%s$, $E_{\\mathrm{surr}}=%s$), the lowest full geodesic stress occurs at iteration %d ($\\lambda_E=%s$, $\\lambda_R=%s$, $\\sigma_{\\mathrm{geo}}=%s$), and the largest area-floor ratio occurs at iteration %d ($\\lambda_E=%s$, $\\lambda_R=%s$, $\\alpha_{0.05}=%s$).",
    df$iteration[[best_sur_idx]],
    fmt_num(df$lambda_E[[best_sur_idx]], 4L),
    fmt_num(df$lambda_R[[best_sur_idx]], 4L),
    fmt_num(df$surrogate_energy[[best_sur_idx]], 4L),
    df$iteration[[best_geo_idx]],
    fmt_num(df$lambda_E[[best_geo_idx]], 4L),
    fmt_num(df$lambda_R[[best_geo_idx]], 4L),
    fmt_num(df$geodesic_stress[[best_geo_idx]], 4L),
    df$iteration[[best_alpha_idx]],
    fmt_num(df$lambda_E[[best_alpha_idx]], 4L),
    fmt_num(df$lambda_R[[best_alpha_idx]], 4L),
    fmt_num(df$area_q05_ratio[[best_alpha_idx]], 4L)
  )
}

make_comparison_table <- function(case_id, all_results, baseline_df) {
  match_idx <- which(vapply(all_results, function(x) identical(x$case$id, case_id), logical(1L)))[1L]
  case_label <- all_results[[match_idx]]$case$label
  pure_row <- baseline_df[baseline_df$case_id == case_id & baseline_df$method_id == "pure_gmds", , drop = FALSE]
  reg_best <- baseline_df[baseline_df$case_id == case_id & baseline_df$family == "Regularized GMDS", , drop = FALSE]
  reg_best <- reg_best[which.min(reg_best$gmds_stress), , drop = FALSE]
  surrogate_rows <- do.call(rbind, lapply(all_results, function(x) {
    if (!identical(x$case$id, case_id)) {
      return(NULL)
    }
    df <- x$trajectory
    final <- df[nrow(df), , drop = FALSE]
    data.frame(
      method = x$method$label,
      runtime_sec = final$elapsed_total_sec[[1L]],
      surrogate_energy = final$surrogate_energy[[1L]],
      geodesic_stress = final$geodesic_stress[[1L]],
      sigma_over_pure = final$geodesic_stress[[1L]] / pure_row$gmds_stress[[1L]],
      sigma_over_regbest = final$geodesic_stress[[1L]] / reg_best$gmds_stress[[1L]],
      procrustes_rmse = final$procrustes_rmse[[1L]],
      stringsAsFactors = FALSE
    )
  }))
  rows <- vapply(seq_len(nrow(surrogate_rows)), function(i) {
    sprintf(
      "%s & %s & %s & %s & %s & %s & %s \\\\",
      surrogate_rows$method[[i]],
      fmt_num(surrogate_rows$runtime_sec[[i]], 3L),
      fmt_num(surrogate_rows$surrogate_energy[[i]], 4L),
      fmt_num(surrogate_rows$geodesic_stress[[i]], 4L),
      fmt_num(surrogate_rows$sigma_over_pure[[i]], 2L),
      fmt_num(surrogate_rows$sigma_over_regbest[[i]], 2L),
      fmt_num(surrogate_rows$procrustes_rmse[[i]], 4L)
    )
  }, character(1L))
  caption <- sprintf(
    "Final surrogate-vs-GMDS comparison for %s. The pure-GMDS baseline uses $\\sigma=%s$ and the best regularized-GMDS baseline from the earlier paraboloid study uses $\\sigma=%s$.",
    case_label,
    fmt_num(pure_row$gmds_stress[[1L]], 4L),
    fmt_num(reg_best$gmds_stress[[1L]], 4L)
  )
  paste(
    sprintf("\\subsubsection{%s}", case_label),
    "\\begin{table}[H]",
    "\\centering",
    "\\scriptsize",
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{lrrrrrr}",
    "\\toprule",
    "Method & $t$ (s) & $E_{\\mathrm{surr}}$ & $\\sigma_{\\mathrm{geo}}$ & $\\sigma/\\sigma_{\\mathrm{pure}}$ & $\\sigma/\\sigma_{\\mathrm{reg.best}}$ & $\\rho$ \\\\",
    "\\midrule",
    paste(rows, collapse = "\n"),
    "\\bottomrule",
    "\\end{tabular}}",
    paste0("\\caption{", caption, "}"),
    "\\end{table}",
    sep = "\n"
  )
}

method_specs <- list(
  list(
    id = "edge_relaxation",
    label = "Edge relaxation surrogate",
    lambda_cols = c("lambda_E", "lambda_R"),
    edge_schedule = schedule_linear(1.0, 0.1, 15L),
    repulsion_schedule = rep(0, 16L)
  ),
  list(
    id = "edge_relaxation_repulsion",
    label = "Edge relaxation + repulsion surrogate",
    lambda_cols = c("lambda_E", "lambda_R"),
    edge_schedule = schedule_linear(1.0, 0.1, 15L),
    repulsion_schedule = schedule_linear(1.0, 0.1, 15L)
  )
)

cases <- lapply(c(12L, 15L), make_case)
all_results <- unlist(lapply(cases, function(case) {
  lapply(method_specs, function(spec) run_case_method(case, spec))
}), recursive = FALSE)

trajectory_df <- do.call(rbind, lapply(all_results, function(x) x$trajectory))
trajectory_no_display <- trajectory_df
trajectory_no_display$display_coords <- NULL
utils::write.csv(trajectory_no_display, metrics_csv, row.names = FALSE)

saveRDS(
  list(
    generated_at = as.character(Sys.time()),
    cases = lapply(cases, function(case) {
      list(id = case$id, label = case$label, edges = case$edges)
    }),
    results = lapply(all_results, function(x) {
      list(
        case_id = x$case$id,
        case_label = x$case$label,
        method_id = x$method$id,
        method_label = x$method$label,
        figure_path = x$figure_path,
        trajectory = x$trajectory
      )
    })
  ),
  rds_path
)

baseline_df <- utils::read.csv(baseline_csv, stringsAsFactors = FALSE)

sections <- c(
  "\\section{Parabolic Edge-Relaxation Surrogate Study}",
  "This Step~8 study evaluates the edge-only surrogate honestly as a surrogate rather than as GMDS itself. Both surrogate variants use the same stronger linear schedules already studied above: $\\lambda_E(0)=1.0$ to $\\lambda_E(15)=0.1$, with $\\lambda_R(0)=1.0$ to $\\lambda_R(15)=0.1$ only for the repulsive surrogate. Unlike the GMDS sections, the objective optimized here is the edge-level surrogate energy $E_{\\mathrm{surr}}$, not the all-pairs fixed-path geodesic energy.",
  "Accordingly, every table below reports both quantities separately: the optimized surrogate energy $E_{\\mathrm{surr}}$ and the full geodesic stress $\\sigma_{\\mathrm{geo}}$ obtained when the resulting embedding is rescored under pure GMDS on the same graph. This makes the fidelity loss relative to GMDS explicit."
)

for (spec in method_specs) {
  sections <- c(
    sections,
    sprintf("\\subsection{%s}", spec$label)
  )
  for (case in cases) {
    idx <- which(vapply(all_results, function(x) identical(x$case$id, case$id) && identical(x$method$id, spec$id), logical(1L)))
    res <- all_results[[idx]]
    df <- res$trajectory
    sections <- c(
      sections,
      sprintf("\\subsubsection{%s}", case$label),
      summarize_case(df),
      write_case_table(df, spec),
      sprintf(
        "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{%s on %s. The panels show selected stored iterates from the beginning, middle, and end of the trajectory.}\\end{figure}",
        res$figure_path,
        spec$label,
        case$label
      )
    )
  }
}

sections <- c(
  sections,
  "\\subsection{Final Surrogate-vs-GMDS Comparison}"
)
for (case in cases) {
  sections <- c(
    sections,
    make_comparison_table(case$id, all_results, baseline_df)
  )
}

writeLines(sections, appendix_tex_path)

message("Wrote edge-relaxation metrics: ", metrics_csv)
message("Wrote edge-relaxation bundle: ", rds_path)
message("Wrote LaTeX appendix snippet: ", appendix_tex_path)
