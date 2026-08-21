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
base_run_tag <- "gmds-v2-paraboloid-regularized-2026-04-01"
tmp_dir <- file.path(manual_root, "tmp", base_run_tag)
pdf_dir <- file.path(manual_root, "reports", base_run_tag)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

appendix_tex_path <- file.path(manual_root, "reports", "gmds_v2_paraboloid_regularized_strong_anchor_appendix_2026-04-01.tex")
metrics_csv <- file.path(tmp_dir, "gmds_v2_paraboloid_regularized_strong_anchor_metrics.csv")
rds_path <- file.path(tmp_dir, "gmds_v2_paraboloid_regularized_strong_anchor_results.rds")

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
    prepared = prepared,
    cmd = cmd
  )
}

select_snapshot_indices <- function(n_frames) {
  raw <- c(1L, 2L, 4L, 7L, 11L, n_frames)
  raw <- raw[raw <= n_frames]
  unique(raw)
}

save_snapshot_grid <- function(case, trajectory_df, output_path) {
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
    title <- sprintf(
      "iter %d, lambda %.2f\nsigma %s, rho %s",
      row$iteration[[1L]],
      row$lambda_anchor[[1L]],
      fmt_num(row$gmds_stress[[1L]], 4L),
      fmt_num(row$procrustes_rmse[[1L]], 4L)
    )
    graphics::mtext(title, side = 3L, line = 0.3, cex = 0.82)
  }
  graphics::mtext(case$label, side = 3L, outer = TRUE, line = -0.3, cex = 1.15, font = 2L)
}

run_case <- function(case) {
  started <- proc.time()[["elapsed"]]
  fit <- grip.optimize.geodesic.mds(
    coords = case$cmd$coords,
    prepared = case$prepared,
    anchor_mode = "cmdscale",
    anchor_weight = 1.0,
    anchor_weight_end = 0.1,
    continuation = "linear",
    engine = "cpp",
    max_iter = 15L,
    return_trace = TRUE,
    n_threads = 0L
  )
  elapsed <- proc.time()[["elapsed"]] - started
  trace <- fit$trace
  frames <- fit$frames
  if (length(frames) != nrow(trace)) {
    stop(sprintf("Expected one stored frame per trace row for %s, got %d frames and %d trace rows", case$id, length(frames), nrow(trace)))
  }

  rows <- lapply(seq_len(nrow(trace)), function(i) {
    coords <- frames[[i]]
    score <- grip.score.geodesic.mds(coords = coords, prepared = case$prepared)
    aligned <- align_to_target_nd(coords, case$truth, allow.reflection = TRUE)
    data.frame(
      case_id = case$id,
      case_label = case$label,
      iteration = trace$iteration[[i]],
      lambda_anchor = trace$anchor_weight[[i]],
      gmds_stress = score$gmds.stress[[1L]],
      procrustes_rmse = aligned$rmse,
      roughness = mesh_roughness(coords, case$adj_list, case$edges),
      area_q05_ratio = area_floor_ratio(coords, case$triangles),
      accepted = trace$accepted[[i]],
      elapsed_total_sec = elapsed,
      stringsAsFactors = FALSE
    )
  })
  trajectory_df <- do.call(rbind, rows)
  trajectory_df$display_coords <- I(lapply(frames, function(coords) {
    align_to_target_nd(coords, case$truth, allow.reflection = TRUE)$aligned
  }))
  figure_path <- file.path(pdf_dir, sprintf("%s_strong_anchor_snapshots.png", case$id))
  save_snapshot_grid(case, trajectory_df, figure_path)

  list(
    case = case,
    fit = fit,
    trajectory = trajectory_df,
    figure_path = figure_path
  )
}

write_case_table <- function(df) {
  rows <- vapply(seq_len(nrow(df)), function(i) {
    sprintf(
      "%d & %s & %s & %s & %s & %s \\\\",
      df$iteration[[i]],
      fmt_num(df$lambda_anchor[[i]], 4L),
      fmt_num(df$gmds_stress[[i]], 4L),
      fmt_num(df$procrustes_rmse[[i]], 4L),
      fmt_num(df$roughness[[i]], 4L),
      fmt_num(df$area_q05_ratio[[i]], 4L)
    )
  }, character(1L))
  paste(
    "\\begin{table}[H]",
    "\\centering",
    "\\scriptsize",
    "\\begin{tabular}{rrrrrr}",
    "\\toprule",
    "iter & $\\lambda_A$ & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ \\\\",
    "\\midrule",
    paste(rows, collapse = "\n"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\caption{Anchor-only regularized GMDS trajectory under the stronger linear tether schedule. Each row corresponds to one stored iterate. Lower $\\sigma$, $\\rho$, and $\\eta$ are better; larger $\\alpha_{0.05}$ is better.}",
    "\\end{table}",
    sep = "\n"
  )
}

trajectory_results <- lapply(lapply(c(12L, 15L), make_case), run_case)
trajectory_df <- do.call(rbind, lapply(trajectory_results, function(x) {
  df <- x$trajectory
  df$display_coords <- NULL
  df
}))
utils::write.csv(trajectory_df, metrics_csv, row.names = FALSE)
saveRDS(trajectory_results, rds_path)

snippet_lines <- c(
  "\\section{Strong Anchor Schedule Experiment}",
  "We reran the anchor-only regularized GMDS model with a much stronger tether to the classical-MDS anchor:",
  "\\[\\lambda_A(0)=1.0,\\qquad \\lambda_A(15)=0.1,\\qquad \\text{linear continuation over }15\\text{ iterations}.\\]",
  "Here $\\alpha_{0.05}$ denotes the ratio of the 5th percentile triangle area to the median triangle area after triangulating each rectangular mesh cell into two deterministic triangles. Larger values therefore indicate fewer near-collapsed cells.",
  "The tables below report the whole trajectory, one row per stored iterate and one row per actual anchor weight. The snapshot figures show a representative subset of the trajectory in 3D after Procrustes alignment to the reference paraboloid.",
  ""
)

for (result in trajectory_results) {
  df <- result$trajectory[, c("iteration", "lambda_anchor", "gmds_stress", "procrustes_rmse", "roughness", "area_q05_ratio")]
  best_sigma_idx <- which.min(df$gmds_stress)
  best_rho_idx <- which.min(df$procrustes_rmse)
  best_alpha_idx <- which.max(df$area_q05_ratio)
  snippet_lines <- c(
    snippet_lines,
    sprintf("\\subsection{%s}", tex_escape(result$case$label)),
    sprintf(
      paste(
        "Under the stronger anchor schedule, the lowest geodesic stress for this case occurs at iteration %d ($\\lambda_A=%s$, $\\sigma=%s$),",
        "the best Procrustes fidelity occurs at iteration %d ($\\lambda_A=%s$, $\\rho=%s$),",
        "and the largest area-floor ratio occurs at iteration %d ($\\lambda_A=%s$, $\\alpha_{0.05}=%s$)."
      ),
      df$iteration[[best_sigma_idx]],
      fmt_num(df$lambda_anchor[[best_sigma_idx]], 4L),
      fmt_num(df$gmds_stress[[best_sigma_idx]], 4L),
      df$iteration[[best_rho_idx]],
      fmt_num(df$lambda_anchor[[best_rho_idx]], 4L),
      fmt_num(df$procrustes_rmse[[best_rho_idx]], 4L),
      df$iteration[[best_alpha_idx]],
      fmt_num(df$lambda_anchor[[best_alpha_idx]], 4L),
      fmt_num(df$area_q05_ratio[[best_alpha_idx]], 4L)
    ),
    write_case_table(df),
    sprintf(
      "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{%s under the stronger anchor-only schedule $\\lambda_A:1.0\\to0.1$. The panels show selected stored iterates from the beginning, middle, and end of the trajectory.}\\end{figure}",
      result$figure_path,
      tex_escape(result$case$label)
    ),
    ""
  )
}

writeLines(snippet_lines, appendix_tex_path)

message("Wrote anchor-trajectory metrics: ", metrics_csv)
message("Wrote anchor-trajectory bundle: ", rds_path)
message("Wrote LaTeX appendix snippet: ", appendix_tex_path)
