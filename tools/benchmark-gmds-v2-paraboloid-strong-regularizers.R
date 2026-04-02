#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) == 0L) {
  stop("Unable to determine script path from commandArgs().")
}
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

manual_root <- file.path(repo_root, "dev", "manual")
base_run_tag <- "gmds-v2-paraboloid-regularized-2026-04-01"
tmp_dir <- file.path(manual_root, "tmp", base_run_tag)
pdf_dir <- file.path(manual_root, "pdf", base_run_tag)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

appendix_tex_path <- file.path(
  manual_root,
  "pdf",
  "gmds_v2_paraboloid_regularized_strong_regularizers_appendix_2026-04-01.tex"
)
metrics_csv <- file.path(tmp_dir, "gmds_v2_paraboloid_regularized_strong_regularizers_metrics.csv")
rds_path <- file.path(tmp_dir, "gmds_v2_paraboloid_regularized_strong_regularizers_results.rds")

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
  prepared <- grip.prepare.graph.geodesic.mds(
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
    grip.plot(
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
      "iter %d, %s\nsigma %s, rho %s",
      row$iteration[[1L]],
      paste(lambda_bits, collapse = ", "),
      fmt_num(row$gmds_stress[[1L]], 4L),
      fmt_num(row$procrustes_rmse[[1L]], 4L)
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
  started <- proc.time()[["elapsed"]]
  fit <- do.call(
    grip.optimize.geodesic.mds,
    c(
      list(
        coords = case$cmd$coords,
        prepared = case$prepared
      ),
      method_spec$args
    )
  )
  elapsed <- proc.time()[["elapsed"]] - started
  trace <- fit$trace
  frames <- fit$frames
  if (length(frames) != nrow(trace)) {
    stop(sprintf(
      "Expected one stored frame per trace row for %s / %s, got %d frames and %d trace rows",
      case$id,
      method_spec$id,
      length(frames),
      nrow(trace)
    ))
  }

  rows <- lapply(seq_len(nrow(trace)), function(i) {
    coords <- frames[[i]]
    score <- grip.score.geodesic.mds(coords = coords, prepared = case$prepared)
    aligned <- align_to_target_nd(coords, case$truth, allow.reflection = TRUE)
    data.frame(
      case_id = case$id,
      case_label = case$label,
      method_id = method_spec$id,
      method_label = method_spec$label,
      iteration = trace$iteration[[i]],
      lambda_A = trace$anchor_weight[[i]],
      lambda_E = trace$edge_spring_weight[[i]],
      lambda_R = trace$repulsion_weight[[i]],
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
  figure_path <- file.path(pdf_dir, sprintf("%s_%s_strong_schedule_snapshots.png", case$id, method_spec$id))
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
  if (identical(method_spec$id, "anchor_repulsion")) {
    colspec <- "{rrrrrrr}"
    header <- "iter & $\\lambda_A$ & $\\lambda_R$ & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ \\\\"
    rows <- vapply(seq_len(nrow(df)), function(i) {
      sprintf(
        "%d & %s & %s & %s & %s & %s & %s \\\\",
        df$iteration[[i]],
        fmt_num(df$lambda_A[[i]], 4L),
        fmt_num(df$lambda_R[[i]], 4L),
        fmt_num(df$gmds_stress[[i]], 4L),
        fmt_num(df$procrustes_rmse[[i]], 4L),
        fmt_num(df$roughness[[i]], 4L),
        fmt_num(df$area_q05_ratio[[i]], 4L)
      )
    }, character(1L))
  } else {
    colspec <- "{rrrrrrr}"
    header <- "iter & $\\lambda_A$ & $\\lambda_E$ & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ \\\\"
    rows <- vapply(seq_len(nrow(df)), function(i) {
      sprintf(
        "%d & %s & %s & %s & %s & %s & %s \\\\",
        df$iteration[[i]],
        fmt_num(df$lambda_A[[i]], 4L),
        fmt_num(df$lambda_E[[i]], 4L),
        fmt_num(df$gmds_stress[[i]], 4L),
        fmt_num(df$procrustes_rmse[[i]], 4L),
        fmt_num(df$roughness[[i]], 4L),
        fmt_num(df$area_q05_ratio[[i]], 4L)
      )
    }, character(1L))
  }
  paste(
    "\\begin{table}[H]",
    "\\centering",
    "\\scriptsize",
    sprintf("\\begin{tabular}%s", colspec),
    "\\toprule",
    header,
    "\\midrule",
    paste(rows, collapse = "\n"),
    "\\bottomrule",
    "\\end{tabular}",
    sprintf(
      "\\caption{%s trajectory under the stronger linear schedules. Each row corresponds to one stored iterate. Lower $\\sigma$, $\\rho$, and $\\eta$ are better; larger $\\alpha_{0.05}$ is better.}",
      method_spec$label
    ),
    "\\end{table}",
    sep = "\n"
  )
}

summarize_case <- function(df, method_spec) {
  best_sigma_idx <- which.min(df$gmds_stress)
  best_rho_idx <- which.min(df$procrustes_rmse)
  best_alpha_idx <- which.max(df$area_q05_ratio)
  if (identical(method_spec$id, "anchor_repulsion")) {
    sprintf(
      "Under the stronger coupled schedules, the lowest geodesic stress for this case occurs at iteration %d ($\\lambda_A=%s$, $\\lambda_R=%s$, $\\sigma=%s$), the best Procrustes fidelity occurs at iteration %d ($\\lambda_A=%s$, $\\lambda_R=%s$, $\\rho=%s$), and the largest area-floor ratio occurs at iteration %d ($\\lambda_A=%s$, $\\lambda_R=%s$, $\\alpha_{0.05}=%s$).",
      df$iteration[[best_sigma_idx]],
      fmt_num(df$lambda_A[[best_sigma_idx]], 4L),
      fmt_num(df$lambda_R[[best_sigma_idx]], 4L),
      fmt_num(df$gmds_stress[[best_sigma_idx]], 4L),
      df$iteration[[best_rho_idx]],
      fmt_num(df$lambda_A[[best_rho_idx]], 4L),
      fmt_num(df$lambda_R[[best_rho_idx]], 4L),
      fmt_num(df$procrustes_rmse[[best_rho_idx]], 4L),
      df$iteration[[best_alpha_idx]],
      fmt_num(df$lambda_A[[best_alpha_idx]], 4L),
      fmt_num(df$lambda_R[[best_alpha_idx]], 4L),
      fmt_num(df$area_q05_ratio[[best_alpha_idx]], 4L)
    )
  } else {
    sprintf(
      "Under the stronger coupled schedules, the lowest geodesic stress for this case occurs at iteration %d ($\\lambda_A=%s$, $\\lambda_E=%s$, $\\sigma=%s$), the best Procrustes fidelity occurs at iteration %d ($\\lambda_A=%s$, $\\lambda_E=%s$, $\\rho=%s$), and the largest area-floor ratio occurs at iteration %d ($\\lambda_A=%s$, $\\lambda_E=%s$, $\\alpha_{0.05}=%s$).",
      df$iteration[[best_sigma_idx]],
      fmt_num(df$lambda_A[[best_sigma_idx]], 4L),
      fmt_num(df$lambda_E[[best_sigma_idx]], 4L),
      fmt_num(df$gmds_stress[[best_sigma_idx]], 4L),
      df$iteration[[best_rho_idx]],
      fmt_num(df$lambda_A[[best_rho_idx]], 4L),
      fmt_num(df$lambda_E[[best_rho_idx]], 4L),
      fmt_num(df$procrustes_rmse[[best_rho_idx]], 4L),
      df$iteration[[best_alpha_idx]],
      fmt_num(df$lambda_A[[best_alpha_idx]], 4L),
      fmt_num(df$lambda_E[[best_alpha_idx]], 4L),
      fmt_num(df$area_q05_ratio[[best_alpha_idx]], 4L)
    )
  }
}

method_specs <- list(
  list(
    id = "anchor_repulsion",
    label = "Anchor-plus-repulsion regularized GMDS",
    lambda_cols = c("lambda_A", "lambda_R"),
    args = list(
      anchor_mode = "cmdscale",
      anchor_weight = 1.0,
      anchor_weight_end = 0.1,
      repulsion_weight = 1.0,
      repulsion_weight_end = 0.1,
      continuation = "linear",
      repulsion_continuation = "linear",
      engine = "cpp",
      max_iter = 15L,
      return_trace = TRUE,
      n_threads = 0L
    )
  ),
  list(
    id = "anchor_edge_spring",
    label = "Anchor-plus-edge-spring regularized GMDS",
    lambda_cols = c("lambda_A", "lambda_E"),
    args = list(
      anchor_mode = "cmdscale",
      anchor_weight = 1.0,
      anchor_weight_end = 0.1,
      edge_spring_weight = 1.0,
      edge_spring_weight_end = 0.1,
      continuation = "linear",
      edge_spring_continuation = "linear",
      engine = "cpp",
      max_iter = 15L,
      return_trace = TRUE,
      n_threads = 0L
    )
  )
)

cases <- lapply(c(12L, 15L), make_case)
results <- unlist(lapply(cases, function(case) {
  lapply(method_specs, function(spec) run_case_method(case, spec))
}), recursive = FALSE)

metrics_df <- do.call(rbind, lapply(results, function(x) x$trajectory))
metrics_no_display <- metrics_df
metrics_no_display$display_coords <- NULL
utils::write.csv(metrics_no_display, metrics_csv, row.names = FALSE)

saveRDS(
  list(
    generated_at = as.character(Sys.time()),
    cases = lapply(cases, function(case) {
      list(id = case$id, label = case$label, edges = case$edges)
    }),
    results = lapply(results, function(x) {
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

sections <- c(
  "\\section{Strong Coupled Regularization Schedules}",
  "We reran the two mixed regularized-GMDS variants under a much stronger pair of linear continuation schedules. In both experiments the anchor term follows $\\lambda_A(0)=1.0$ and $\\lambda_A(15)=0.1$, while the companion regularizer follows the same schedule: $\\lambda_R(0)=1.0$ to $\\lambda_R(15)=0.1$ for the repulsion variant and $\\lambda_E(0)=1.0$ to $\\lambda_E(15)=0.1$ for the edge-spring variant.",
  "As above, $\\alpha_{0.05}$ denotes the ratio of the 5th percentile triangle area to the median triangle area after splitting each rectangular cell into two deterministic triangles. Larger values therefore indicate fewer near-collapsed cells."
)

for (spec in method_specs) {
  sections <- c(
    sections,
    sprintf("\\subsection{%s}", spec$label)
  )
  for (case in cases) {
    match_idx <- vapply(results, function(x) identical(x$case$id, case$id) && identical(x$method$id, spec$id), logical(1L))
    res <- results[[which(match_idx)]]
    df <- res$trajectory
    sections <- c(
      sections,
      sprintf("\\subsubsection{%s}", case$label),
      summarize_case(df, spec),
      write_case_table(df, spec),
      sprintf(
        "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{%s on %s under the stronger coupled schedules. The panels show selected stored iterates from the beginning, middle, and end of the trajectory.}\\end{figure}",
        res$figure_path,
        spec$label,
        case$label
      )
    )
  }
}

writeLines(sections, appendix_tex_path)

message("Wrote strong-regularizer metrics: ", metrics_csv)
message("Wrote strong-regularizer bundle: ", rds_path)
message("Wrote LaTeX appendix snippet: ", appendix_tex_path)
