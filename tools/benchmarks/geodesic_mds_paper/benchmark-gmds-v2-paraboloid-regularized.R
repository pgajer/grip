#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) == 0L) {
  stop("Unable to determine script path from commandArgs().")
}
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

run_tag <- "gmds-v2-paraboloid-regularized-2026-04-01"
manual_root <- file.path(repo_root, "dev", "manual")
tmp_dir <- file.path(manual_root, "tmp", run_tag)
pdf_dir <- file.path(manual_root, "pdf", run_tag)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

tex_path <- file.path(manual_root, "pdf", "gmds_v2_paraboloid_regularized_report_2026-04-01.tex")
pdf_path <- file.path(manual_root, "pdf", "gmds_v2_paraboloid_regularized_report_2026-04-01.pdf")
rds_path <- file.path(tmp_dir, "gmds_v2_paraboloid_regularized_results.rds")
metrics_csv <- file.path(tmp_dir, "gmds_v2_paraboloid_regularized_metrics.csv")

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

lambda_label <- function(anchor_start = 0,
                         anchor_end = anchor_start,
                         edge_start = 0,
                         edge_end = edge_start,
                         rep_start = 0,
                         rep_end = rep_start) {
  parts <- character(0L)
  if (anchor_start > 0 || anchor_end > 0) {
    parts <- c(parts, sprintf("A:%s->%s", fmt_num(anchor_start, 2L), fmt_num(anchor_end, 2L)))
  }
  if (edge_start > 0 || edge_end > 0) {
    parts <- c(parts, sprintf("E:%s->%s", fmt_num(edge_start, 2L), fmt_num(edge_end, 2L)))
  }
  if (rep_start > 0 || rep_end > 0) {
    parts <- c(parts, sprintf("R:%s->%s", fmt_num(rep_start, 2L), fmt_num(rep_end, 2L)))
  }
  if (length(parts) == 0L) {
    return("--")
  }
  paste(parts, collapse = ", ")
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
  started <- proc.time()[["elapsed"]]
  cmd <- classical_mds_embedding(prepared, dim = 3L, eig = TRUE)
  cmd_elapsed <- proc.time()[["elapsed"]] - started

  list(
    id = sprintf("paraboloid_%dx%d", side, side),
    label = sprintf("Orthogonal paraboloid mesh %dx%d", side, side),
    side = side,
    truth = bundle$coords_surface,
    param = bundle$coords_param,
    edges = bundle$edges,
    adj_list = grip.build.adj.from.edges(bundle$edges, n = bundle$n)$adj_list,
    triangles = grid_mesh_triangles(side, side),
    prepared = prepared,
    cmd = cmd,
    cmd_elapsed = cmd_elapsed
  )
}

compute_metrics <- function(case,
                            spec,
                            coords,
                            elapsed_sec = NA_real_,
                            fit = NULL) {
  aligned <- align_to_target_nd(coords, case$truth, allow.reflection = TRUE)
  score_df <- if (!is.null(fit) && !is.null(fit$score)) {
    fit$score
  } else {
    grip.score.geodesic.mds(coords = coords, prepared = case$prepared)
  }

  data.frame(
    case_id = case$id,
    case_label = case$label,
    side = case$side,
    n = nrow(case$truth),
    method_id = spec$id,
    method_label = spec$label,
    family = spec$family,
    lambda_label = spec$lambda_label,
    elapsed_sec = as.double(elapsed_sec),
    gmds_energy = score_df$gmds.energy[[1L]],
    gmds_stress = score_df$gmds.stress[[1L]],
    gmds_raw_stress = score_df$gmds.raw_stress[[1L]],
    procrustes_rmse = aligned$rmse,
    roughness = mesh_roughness(coords, case$adj_list, case$edges),
    area_q05_ratio = area_floor_ratio(coords, case$triangles),
    anchor_energy = score_df$anchor.energy[[1L]],
    edge_spring_energy = score_df$edge.spring.energy[[1L]],
    repulsion_energy = score_df$repulsion.energy[[1L]],
    iterations = if (is.null(fit$trace) || nrow(fit$trace) == 0L) 0L else max(fit$trace$iteration),
    stringsAsFactors = FALSE
  )
}

method_specs <- list(
  list(
    id = "reference",
    label = "Reference surface",
    family = "Reference",
    kind = "reference",
    lambda_label = "--"
  ),
  list(
    id = "cmdscale",
    label = "Classical MDS",
    family = "Baseline",
    kind = "cmdscale",
    lambda_label = "--"
  ),
  list(
    id = "pure_gmds",
    label = "Pure GMDS",
    family = "Pure GMDS",
    kind = "gmds",
    lambda_label = "--",
    args = list(engine = "cpp", max_iter = 15L, return_trace = TRUE, n_threads = 0L)
  ),
  list(
    id = "reg_anchor",
    label = "Regularized GMDS (anchor)",
    family = "Regularized GMDS",
    kind = "gmds",
    lambda_label = lambda_label(anchor_start = 0.10, anchor_end = 0.02),
    args = list(
      engine = "cpp",
      max_iter = 15L,
      return_trace = TRUE,
      n_threads = 0L,
      anchor_mode = "cmdscale",
      anchor_weight = 0.10,
      anchor_weight_end = 0.02,
      continuation = "linear"
    )
  ),
  list(
    id = "reg_anchor_rep",
    label = "Regularized GMDS (anchor + repulsion)",
    family = "Regularized GMDS",
    kind = "gmds",
    lambda_label = lambda_label(anchor_start = 0.10, anchor_end = 0.02, rep_start = 0.20, rep_end = 0.05),
    args = list(
      engine = "cpp",
      max_iter = 15L,
      return_trace = TRUE,
      n_threads = 0L,
      anchor_mode = "cmdscale",
      anchor_weight = 0.10,
      anchor_weight_end = 0.02,
      continuation = "linear",
      repulsion_weight = 0.20,
      repulsion_weight_end = 0.05,
      repulsion_continuation = "linear",
      repulsion_quantile = 0.40,
      repulsion_scale = 0.60,
      repulsion_cap_quantile = 1.00,
      repulsion_hop_min = 2L
    )
  ),
  list(
    id = "reg_anchor_edge",
    label = "Regularized GMDS (anchor + edge spring)",
    family = "Regularized GMDS",
    kind = "gmds",
    lambda_label = lambda_label(anchor_start = 0.10, anchor_end = 0.02, edge_start = 0.25, edge_end = 0.10),
    args = list(
      engine = "cpp",
      max_iter = 15L,
      return_trace = TRUE,
      n_threads = 0L,
      anchor_mode = "cmdscale",
      anchor_weight = 0.10,
      anchor_weight_end = 0.02,
      continuation = "linear",
      edge_spring_weight = 0.25,
      edge_spring_weight_end = 0.10,
      edge_spring_continuation = "linear"
    )
  )
)

run_method <- function(case, spec) {
  if (identical(spec$kind, "reference")) {
    row <- compute_metrics(case = case, spec = spec, coords = case$truth, elapsed_sec = NA_real_, fit = list(trace = NULL))
    return(list(coords = case$truth, display_coords = case$truth, metrics = row, fit = NULL))
  }

  if (identical(spec$kind, "cmdscale")) {
    coords <- case$cmd$coords
    row <- compute_metrics(case = case, spec = spec, coords = coords, elapsed_sec = case$cmd_elapsed, fit = list(trace = NULL))
    return(list(
      coords = coords,
      display_coords = align_to_target_nd(coords, case$truth, allow.reflection = TRUE)$aligned,
      metrics = row,
      fit = NULL
    ))
  }

  started <- proc.time()[["elapsed"]]
  fit <- do.call(
    grip.optimize.geodesic.mds,
    c(list(coords = case$cmd$coords, prepared = case$prepared), spec$args)
  )
  elapsed <- proc.time()[["elapsed"]] - started
  row <- compute_metrics(case = case, spec = spec, coords = fit$coords, elapsed_sec = elapsed, fit = fit)
  list(
    coords = fit$coords,
    display_coords = align_to_target_nd(fit$coords, case$truth, allow.reflection = TRUE)$aligned,
    metrics = row,
    fit = fit
  )
}

build_title <- function(method_result) {
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
  grDevices::png(
    output_path,
    width = 1280L * 2L,
    height = 820L * 2L,
    res = 180,
    bg = "#ffffff"
  )
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(2L, 3L), mar = c(1.2, 1.2, 3.0, 0.4), oma = c(0, 0, 1.2, 0))

  for (method in methods) {
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
    graphics::mtext(build_title(method), side = 3L, line = 0.3, cex = 0.78)
  }
  graphics::mtext(case_result$case$label, side = 3L, outer = TRUE, line = -0.3, cex = 1.15, font = 2L)
}

write_case_table <- function(df) {
  rows <- vapply(seq_len(nrow(df)), function(i) {
    sprintf(
      "%s & %s & %s & %s & %s & %s & %s \\\\",
      tex_escape(df$method_label[[i]]),
      tex_escape(df$lambda_label[[i]]),
      fmt_time(df$elapsed_sec[[i]]),
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
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{lllllll}",
    "\\toprule",
    "Method & $\\lambda$ & $t$ (s) & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ \\\\",
    "\\midrule",
    paste(rows, collapse = "\n"),
    "\\bottomrule",
    "\\end{tabular}}",
    "\\caption{Regularized-GMDS comparison on the selected paraboloid mesh. The $\\lambda$ column records the exact continuation schedule used for each embedding. Lower $\\sigma$, $\\rho$, and $\\eta$ are better; larger $\\alpha_{0.05}$ is better.}",
    "\\end{table}",
    sep = "\n"
  )
}

case_summary_paragraph <- function(df) {
  pure <- df[df$method_id == "pure_gmds", , drop = FALSE]
  regularized <- df[df$family == "Regularized GMDS", , drop = FALSE]
  best_sigma <- regularized[which.min(regularized$gmds_stress), , drop = FALSE]
  best_rho <- regularized[which.min(regularized$procrustes_rmse), , drop = FALSE]
  best_eta <- regularized[which.min(regularized$roughness), , drop = FALSE]

  sprintf(
    "Pure GMDS provides the unregularized baseline on %s: it reaches $\\sigma=%s$ with $\\rho=%s$ and roughness $\\eta=%s$ in %s seconds. Among the regularized variants, %s gives the lowest geodesic stress ($\\sigma=%s$, %s s), %s gives the best Procrustes fidelity ($\\rho=%s$), and %s gives the lowest roughness ($\\eta=%s$). This is the central Step~7 comparison: regularization is treated as its own family and judged against pure GMDS by explicit trade-offs rather than by folding it back into the base definition.",
    tex_escape(df$case_label[[1L]]),
    fmt_num(pure$gmds_stress[[1L]], 4L),
    fmt_num(pure$procrustes_rmse[[1L]], 4L),
    fmt_num(pure$roughness[[1L]], 4L),
    fmt_time(pure$elapsed_sec[[1L]]),
    tex_escape(best_sigma$method_label[[1L]]),
    fmt_num(best_sigma$gmds_stress[[1L]], 4L),
    fmt_time(best_sigma$elapsed_sec[[1L]]),
    tex_escape(best_rho$method_label[[1L]]),
    fmt_num(best_rho$procrustes_rmse[[1L]], 4L),
    tex_escape(best_eta$method_label[[1L]]),
    fmt_num(best_eta$roughness[[1L]], 4L)
  )
}

case_figure_path <- function(case_id) {
  file.path(pdf_dir, sprintf("%s_regularized_grid.png", case_id))
}

cases <- lapply(c(12L, 15L), make_case)
case_results <- lapply(cases, function(case) {
  methods <- lapply(method_specs, run_method, case = case)
  metrics <- do.call(rbind, lapply(methods, `[[`, "metrics"))
  out <- list(case = case, methods = methods, metrics = metrics)
  save_case_panel_grid(out, case_figure_path(case$id))
  out
})

metrics_df <- do.call(rbind, lapply(case_results, `[[`, "metrics"))
utils::write.csv(metrics_df, metrics_csv, row.names = FALSE)

saveRDS(
  list(
    run_tag = run_tag,
    case_results = case_results,
    metrics = metrics_df
  ),
  rds_path
)

case_sections <- vapply(case_results, function(case_result) {
  df <- subset(case_result$metrics, method_id != "reference")
  paste(
    sprintf("\\section{%s}", tex_escape(case_result$case$label)),
    case_summary_paragraph(df),
    "",
    sprintf(
      paste(
        "\\begin{figure}[p]\\centering",
        "\\includegraphics[width=0.98\\linewidth]{%s}",
        "\\caption{%s. Panels show the reference surface, classical MDS, pure GMDS, anchor-only regularized GMDS, anchor-plus-repulsion regularized GMDS, and anchor-plus-edge-spring regularized GMDS. The regularization schedules in this figure are explicit: A 0.10$\\\\to$0.02; A 0.10$\\\\to$0.02 with R 0.20$\\\\to$0.05; and A 0.10$\\\\to$0.02 with E 0.25$\\\\to$0.10.}",
        "\\label{fig:%s}",
        "\\end{figure}",
        sep = ""
      ),
      case_figure_path(case_result$case$id),
      tex_escape(case_result$case$label),
      case_result$case$id
    ),
    write_case_table(df),
    sep = "\n\n"
  )
}, character(1L))

takeaway_text <- paste(
  "\\section{Interpretation}",
  "The dedicated regularized-GMDS study separates the regularized family from both pure GMDS and the edge-relaxation surrogate. That separation is the main methodological goal of Step 7: regularization is no longer implicit or mixed into the base algorithm.",
  "Across the paraboloid meshes, the comparison should be read as a trade-off study. Pure GMDS remains the reference point for the fixed-path objective itself; anchor-based regularization can preserve more of the classical-MDS shape; repulsion and edge-spring terms can then be evaluated as distinct ways of controlling collapse or local distortion. The tables make those trade-offs explicit by keeping the continuation schedules and runtime visible for every embedding.",
  sep = "\n\n"
)

tex_lines <- c(
  "\\documentclass[11pt]{article}",
  "\\usepackage[margin=1in]{geometry}",
  "\\usepackage{graphicx}",
  "\\usepackage{booktabs}",
  "\\usepackage{float}",
  "\\usepackage{amsmath}",
  "\\usepackage[T1]{fontenc}",
  "\\usepackage[utf8]{inputenc}",
  "\\title{Paraboloid Regularized-GMDS Study}",
  "\\author{GMDS v2 cleanup: Step 7}",
  "\\date{2026-04-01}",
  "\\begin{document}",
  "\\maketitle",
  "\\section{Purpose}",
  "This report implements Step 7 of the GMDS v2 cleanup plan. The focus here is deliberately narrower than in the broader paraboloid suite: regularized GMDS is treated as its own method family, with pure GMDS retained only as the fixed-path baseline and classical MDS retained only as the initialization/reference baseline.",
  "\\section{Algorithms compared}",
  "Let $F_{\\mathrm{geo}}(Z)=\\tfrac12\\sum_{i<j}(h_{ij}(Z)-g_{ij})^2$ denote the pure fixed-path GMDS objective. The benchmark compares the following embedding procedures on the same prepared graph and from the same classical-MDS initialization:",
  "\\begin{itemize}",
  "\\item classical MDS on the graph-geodesic distance matrix,",
  "\\item pure GMDS minimizing $F_{\\mathrm{geo}}$,",
  "\\item anchor-only regularized GMDS with",
  "\\[F_A(Z)=F_{\\mathrm{geo}}(Z)+\\lambda_A(t)\\lVert Z-Z_{\\mathrm{cmd}}\\rVert_F^2,\\]",
  "\\item anchor-plus-repulsion regularized GMDS with",
  "\\[F_{AR}(Z)=F_{\\mathrm{geo}}(Z)+\\lambda_A(t)\\lVert Z-Z_{\\mathrm{cmd}}\\rVert_F^2+\\lambda_R(t)\\sum_{(u,v)}[r_{uv}-\\lVert z_u-z_v\\rVert]_+^2/2,\\]",
  "\\item anchor-plus-edge-spring regularized GMDS with",
  "\\[F_{AE}(Z)=F_{\\mathrm{geo}}(Z)+\\lambda_A(t)\\lVert Z-Z_{\\mathrm{cmd}}\\rVert_F^2+\\lambda_E(t)\\sum_{(u,v)\\in E}(\\lVert z_u-z_v\\rVert-w_{uv})^2/2.\\]",
  "\\end{itemize}",
  "All continuation schedules are linear, all runs use the compiled optimizer, and every reported $t$ is the wall-clock embedding time for that method only.",
  "\\section{Experiment design}",
  "The graph family is the orthogonal weighted paraboloid mesh, chosen because it remains simple enough to inspect visually while still exposing the pathologies that motivated the revised GMDS formulation. This Step 7 report concentrates on $12\\times12$ and $15\\times15$ meshes, so the study stays readable while still showing how the regularized family scales with mesh density.",
  paste(case_sections, collapse = "\n\n"),
  takeaway_text,
  "\\end{document}"
)

writeLines(tex_lines, tex_path)

message("Wrote benchmark bundle: ", rds_path)
message("Wrote metrics CSV: ", metrics_csv)
message("Wrote LaTeX report: ", tex_path)
message("Expected PDF path after latexmk: ", pdf_path)
