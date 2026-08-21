#!/usr/bin/env Rscript

run_tag <- "paraboloid-gmds-bending-2026-03-31"
manual_root <- file.path("output", "geodesic_mds_paper")
tmp_dir <- file.path(manual_root, "tmp", run_tag)
pdf_dir <- file.path(manual_root, "reports", run_tag)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'devtools' or the 'grip' package to run the paraboloid bending benchmark.")
}

fmt_num <- function(x, digits = 4L) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}

fmt_time <- function(x) {
  ifelse(
    is.finite(x),
    ifelse(
      x < 1,
      formatC(x, format = "f", digits = 3L),
      formatC(x, format = "f", digits = 2L)
    ),
    "NA"
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

triangle_areas <- function(coords, triangles) {
  v1 <- coords[triangles[, 2L], , drop = FALSE] - coords[triangles[, 1L], , drop = FALSE]
  v2 <- coords[triangles[, 3L], , drop = FALSE] - coords[triangles[, 1L], , drop = FALSE]
  cross <- cbind(
    v1[, 2L] * v2[, 3L] - v1[, 3L] * v2[, 2L],
    v1[, 3L] * v2[, 1L] - v1[, 1L] * v2[, 3L],
    v1[, 1L] * v2[, 2L] - v1[, 2L] * v2[, 1L]
  )
  0.5 * sqrt(rowSums(cross^2))
}

mesh_roughness <- function(coords, adj_list, edges) {
  centered <- sweep(coords, 2L, colMeans(coords), FUN = "-")
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

area_floor_ratio <- function(coords, triangles) {
  areas <- triangle_areas(coords, triangles)
  med <- stats::median(areas)
  if (!is.finite(med) || med <= 0) {
    return(NA_real_)
  }
  as.double(stats::quantile(areas, probs = 0.05, names = FALSE)) / med
}

new_case <- function(side = 12L, iter_budget = 25L) {
  bundle <- mesh.surface.graph(
    side,
    side,
    surface = "paraboloid",
    amplitude = 0.35,
    connectivity = "orthogonal",
    normalize = "median"
  )
  prepared <- prepare.geodesic.kk(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average"
  )
  started <- proc.time()[["elapsed"]]
  cmd <- grip:::grip.classical.mds.embedding(prepared, dim = 3L, eig = TRUE)
  cmd.elapsed <- proc.time()[["elapsed"]] - started

  list(
    side = side,
    iter_budget = as.integer(iter_budget),
    bundle = bundle,
    truth = bundle$coords_surface,
    edges = bundle$edges,
    adj_list = prepared$adj_list,
    triangles = grid_mesh_triangles(side, side),
    prepared = prepared,
    cmd = cmd,
    cmd_elapsed = cmd.elapsed,
    bending_stencils = grip:::grip.rectangular.grid.bending.stencils(side, side)
  )
}

evaluate_layout <- function(case,
                            coords,
                            setting_id,
                            setting_label,
                            elapsed_sec = NA_real_,
                            fit = NULL,
                            anchor_weight_override = NULL,
                            bending_weight_override = NULL) {
  anchor.coords <- if (is.null(fit)) NULL else fit$anchor_coords
  anchor.weight <- if (!is.null(anchor_weight_override)) {
    anchor_weight_override
  } else if (is.null(fit)) {
    0
  } else {
    fit$final_anchor_weight
  }
  bend.weight <- if (!is.null(bending_weight_override)) {
    bending_weight_override
  } else if (is.null(fit)) {
    0
  } else {
    fit$final_bending_weight
  }
  score <- grip.score.geodesic.mds(
    coords,
    prepared = case$prepared,
    anchor_coords = anchor.coords,
    anchor_weight = anchor.weight,
    bending_stencils = case$bending_stencils,
    bending_weight = bend.weight
  )
  aligned <- grip:::grip.align.to.target.nd(coords, case$truth, allow.reflection = TRUE)
  trace <- if (is.null(fit)) NULL else fit$trace
  steps <- if (is.null(trace) || nrow(trace) == 0L) numeric(0L) else trace$step[is.finite(trace$step)]
  data.frame(
    setting_id = setting_id,
    setting_label = setting_label,
    elapsed_sec = as.double(elapsed_sec),
    n_threads_used = if (is.null(fit)) 1L else fit$n_threads_used,
    gmds_stress = score$gmds.stress[[1L]],
    gmds_raw_stress = score$gmds.raw_stress[[1L]],
    procrustes_rmse = aligned$rmse,
    roughness = mesh_roughness(coords, case$adj_list, case$edges),
    area_q05_ratio = area_floor_ratio(coords, case$triangles),
    anchor_weight = anchor.weight,
    bend_weight = bend.weight,
    bend_raw_penalty = score$bend.raw.penalty[[1L]],
    bend_energy = score$bend.energy[[1L]],
    mean_step = if (length(steps) > 0L) mean(steps) else NA_real_,
    final_step = if (length(steps) > 0L) utils::tail(steps, 1L) else NA_real_,
    iterations = if (is.null(trace) || nrow(trace) == 0L) 0L else max(trace$iteration),
    stringsAsFactors = FALSE
  )
}

run_cmdscale <- function(case) {
  list(
    row = evaluate_layout(
      case = case,
      coords = case$cmd$coords,
      setting_id = "cmdscale",
      setting_label = "cmdscale",
      elapsed_sec = case$cmd_elapsed,
      fit = NULL
    ),
    coords = case$cmd$coords
  )
}

run_config <- function(case,
                       setting_id,
                       setting_label,
                       ...) {
  started <- proc.time()[["elapsed"]]
  fit <- grip.optimize.geodesic.mds(
    coords = case$cmd$coords,
    prepared = case$prepared,
    max_iter = case$iter_budget,
    engine = "cpp",
    return_trace = TRUE,
    n_threads = 1L,
    ...
  )
  elapsed <- proc.time()[["elapsed"]] - started
  cat(sprintf("completed %s in %ss\n", setting_label, fmt_time(elapsed)))
  list(
    row = evaluate_layout(
      case = case,
      coords = fit$coords,
      setting_id = setting_id,
      setting_label = setting_label,
      elapsed_sec = elapsed,
      fit = fit
    ),
    fit = fit
  )
}

run_beta_sweep <- function(case,
                           beta_values,
                           anchor_lambda = 0) {
  rows <- list()
  fits <- vector("list", length(beta_values))
  for (i in seq_along(beta_values)) {
    beta <- beta_values[[i]]
    label <- if (anchor_lambda > 0) {
      sprintf("lambda = %.3f, beta = %.3f", anchor_lambda, beta)
    } else {
      sprintf("beta = %.3f", beta)
    }
    setting_id <- if (anchor_lambda > 0) {
      sprintf("hybrid_lambda_%s_beta_%s", format(anchor_lambda, trim = TRUE), format(beta, trim = TRUE))
    } else {
      sprintf("bend_beta_%s", format(beta, trim = TRUE))
    }
    fit <- run_config(
      case = case,
      setting_id = setting_id,
      setting_label = label,
      bending_stencils = case$bending_stencils,
      bending_weight = beta,
      anchor_mode = if (anchor_lambda > 0) "cmdscale" else "none",
      anchor_weight = anchor_lambda,
      anchor_weight_end = anchor_lambda,
      continuation = "constant"
    )
    fit$row$beta <- beta
    fit$row$lambda <- anchor_lambda
    rows[[i]] <- fit$row
    fits[[i]] <- fit
  }
  list(
    metrics = do.call(rbind, rows),
    fits = fits
  )
}

pick_best_tradeoff <- function(df, cmd_sigma) {
  candidates <- df[df$bend_weight > 0 & df$gmds_stress < cmd_sigma, , drop = FALSE]
  if (nrow(candidates) == 0L) {
    candidates <- df[df$bend_weight > 0, , drop = FALSE]
  }
  candidates[order(
    candidates$roughness,
    -candidates$area_q05_ratio,
    candidates$procrustes_rmse,
    candidates$gmds_stress,
    candidates$elapsed_sec
  ), , drop = FALSE][1L, , drop = FALSE]
}

trajectory_metrics <- function(case,
                               fit,
                               setting_id_prefix,
                               setting_label_prefix) {
  do.call(rbind, lapply(seq_along(fit$frames), function(i) {
    iteration <- i - 1L
    row <- evaluate_layout(
      case = case,
      coords = fit$frames[[i]],
      setting_id = sprintf("%s_%02d", setting_id_prefix, iteration),
      setting_label = sprintf("%s iter %d", setting_label_prefix, iteration),
      elapsed_sec = NA_real_,
      fit = fit
    )
    row$iteration <- iteration
    row$accepted_step <- if (iteration == 0L) NA_real_ else fit$trace$step[[iteration + 1L]]
    row$gradient_norm_trace <- fit$trace$gradient_norm[[iteration + 1L]]
    row
  }))
}

save_surface_panel_grid <- function(case,
                                    coords_list,
                                    title_list,
                                    path,
                                    nrow,
                                    ncol,
                                    heading,
                                    truth_index = 1L) {
  grDevices::png(
    path,
    width = 1080L * ncol,
    height = 900L * nrow,
    res = 180,
    bg = "#ffffff"
  )
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(
    mfrow = c(nrow, ncol),
    mar = c(1.2, 1.2, 2.8, 0.4),
    oma = c(0, 0, 1.6, 0)
  )

  for (i in seq_len(nrow * ncol)) {
    if (i > length(coords_list)) {
      graphics::plot.new()
      next
    }
    coords <- coords_list[[i]]
    if (i != truth_index) {
      coords <- grip:::grip.align.to.target.nd(
        coords,
        case$truth,
        allow.reflection = TRUE
      )$aligned
    }
    plot.layout(
      coords = coords,
      edges = case$edges,
      projection = "ortho",
      azimuth = 35,
      elevation = 24,
      vertex.col = if (i == truth_index) "#bc6c25" else "#3a5a40",
      edge.col = "#adb5bd",
      main = ""
    )
    graphics::title(main = title_list[[i]], cex.main = 0.72)
  }

  graphics::mtext(heading, outer = TRUE, cex = 1.15, font = 2)
}

save_trace_compare_plot <- function(untethered_df,
                                    regularized_df,
                                    path) {
  grDevices::png(path, width = 2200, height = 1400, res = 180, bg = "#ffffff")
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(2L, 2L), mar = c(4, 4, 2.5, 1))

  plot_pair <- function(y1, y2, title, ylab) {
    graphics::plot(
      untethered_df$iteration,
      y1,
      type = "b",
      pch = 16,
      col = "#b23a48",
      xlab = "GMDS iteration",
      ylab = ylab,
      main = title
    )
    graphics::lines(
      regularized_df$iteration,
      y2,
      type = "b",
      pch = 17,
      col = "#2d6a4f"
    )
    graphics::legend(
      "topright",
      legend = c("untethered", "anchor + bend"),
      col = c("#b23a48", "#2d6a4f"),
      pch = c(16, 17),
      bty = "n",
      cex = 0.9
    )
  }

  plot_pair(untethered_df$gmds_stress, regularized_df$gmds_stress, expression(sigma ~ "over iterations"), expression(sigma))
  plot_pair(untethered_df$roughness, regularized_df$roughness, expression(eta ~ "over iterations"), expression(eta))
  plot_pair(untethered_df$area_q05_ratio, regularized_df$area_q05_ratio, expression(alpha[0.05] ~ "over iterations"), expression(alpha[0.05]))
  plot_pair(untethered_df$bend_raw_penalty, regularized_df$bend_raw_penalty, expression(B ~ "over iterations"), expression(B))
}

write_method_table <- function(df) {
  rows <- apply(df, 1L, function(row) {
    paste(
      tex_escape(row[["setting_label"]]),
      fmt_num(as.numeric(row[["anchor_weight"]]), digits = 3L),
      fmt_num(as.numeric(row[["bend_weight"]]), digits = 3L),
      fmt_time(as.numeric(row[["elapsed_sec"]])),
      fmt_num(as.numeric(row[["gmds_stress"]])),
      fmt_num(as.numeric(row[["procrustes_rmse"]])),
      fmt_num(as.numeric(row[["roughness"]])),
      fmt_num(as.numeric(row[["area_q05_ratio"]])),
      fmt_num(as.numeric(row[["bend_raw_penalty"]])),
      sep = " & "
    )
  })
  paste0(rows, " \\\\")
}

write_sweep_table <- function(df) {
  rows <- apply(df, 1L, function(row) {
    paste(
      fmt_num(as.numeric(row[["bend_weight"]]), digits = 3L),
      fmt_time(as.numeric(row[["elapsed_sec"]])),
      fmt_num(as.numeric(row[["gmds_stress"]])),
      fmt_num(as.numeric(row[["procrustes_rmse"]])),
      fmt_num(as.numeric(row[["roughness"]])),
      fmt_num(as.numeric(row[["area_q05_ratio"]])),
      fmt_num(as.numeric(row[["bend_raw_penalty"]])),
      fmt_time(as.numeric(row[["mean_step"]])),
      sep = " & "
    )
  })
  paste0(rows, " \\\\")
}

write_trajectory_table <- function(df) {
  rows <- apply(df, 1L, function(row) {
    paste(
      row[["iteration"]],
      fmt_num(as.numeric(row[["gmds_stress"]])),
      fmt_num(as.numeric(row[["procrustes_rmse"]])),
      fmt_num(as.numeric(row[["roughness"]])),
      fmt_num(as.numeric(row[["area_q05_ratio"]])),
      fmt_num(as.numeric(row[["bend_raw_penalty"]])),
      fmt_time(as.numeric(row[["accepted_step"]])),
      sep = " & "
    )
  })
  paste0(rows, " \\\\")
}

write_report <- function(case,
                         selected_df,
                         bending_df,
                         hybrid_df,
                         best_bending,
                         best_hybrid,
                         untethered_traj,
                         regularized_traj,
                         regularized_runtime,
                         tex_path) {
  selected_traj <- regularized_traj[
    regularized_traj$iteration %in% c(0L, 1L, 2L, 5L, 10L, 15L, 20L, 25L),
    ,
    drop = FALSE
  ]

  tex <- c(
    "\\documentclass[11pt]{article}",
    "\\usepackage[margin=1in]{geometry}",
    "\\usepackage{graphicx}",
    "\\usepackage{booktabs}",
    "\\usepackage{float}",
    "\\usepackage{amsmath}",
    "\\usepackage{pdflscape}",
    "\\usepackage{hyperref}",
    "\\title{Paraboloid GMDS Bending-Penalty Report}",
    "\\author{Focused 12x12 paraboloid study}",
    "\\date{March 31, 2026}",
    "",
    "\\begin{document}",
    "\\maketitle",
    "",
    "\\begin{abstract}",
    "This report tests the next regularization idea suggested by the earlier paraboloid pathology study: replace the first-order neighborhood-averaging penalty with a second-order rectangular-grid bending penalty. All runs use the same orthogonal $12\\times12$ paraboloid mesh, the same tie-averaged geodesic-MDS objective, and the same 25-step post-processing budget after the classical-MDS initialization. The central question is whether penalizing discrete curvature along the grid rows and columns suppresses the visually unsettling singularity more effectively than first-order smoothing did.",
    "\\end{abstract}",
    "",
    "\\section{Augmented Objective}",
    "All GMDS runs in this report start from the same classical-MDS embedding of the all-pairs graph-geodesic distance matrix and then run 25 extra GMDS correction steps with Armijo backtracking. The new objective is",
    "\\[",
    "E_{\\lambda,\\beta}(Z) = E_{\\mathrm{GMDS}}(Z) + \\lambda \\lVert Z - A \\rVert_F^2 + \\beta R_{\\mathrm{bend}}(Z),",
    "\\]",
    "where $A$ is the classical-MDS anchor and",
    "\\[",
    "R_{\\mathrm{bend}}(Z) = \\frac{1}{|\\mathcal{S}|} \\sum_{(a,b,c) \\in \\mathcal{S}} \\lVert z_a - 2 z_b + z_c \\rVert^2.",
    "\\]",
    "Here $\\mathcal{S}$ contains every horizontal and vertical three-vertex stencil on the rectangular mesh. In words, $R_{\\mathrm{bend}}$ penalizes concentrated second differences, so it directly targets sharp folds or creases instead of only discouraging first-order local variation.",
    "",
    "\\section{Selected Methods}",
    sprintf(
      "The selected comparison below shows the reference surface, classical MDS, untethered GMDS, anchor-only GMDS with $\\lambda=0.1$, the best bending-only run from the $\\beta$ sweep, and the best anchor-plus-bending run from the $(\\lambda,\\beta)$ sweep. In this study, the best bending-only run was %s and the best hybrid run was %s.",
      tex_escape(best_bending$setting_label[[1L]]),
      tex_escape(best_hybrid$setting_label[[1L]])
    ),
    sprintf(
      "Numerically, the best bending-only run reached $\\sigma=%s$, $\\rho=%s$, $\\eta=%s$, $\\alpha_{0.05}=%s$, and $B=%s$. The best hybrid run reached $\\sigma=%s$, $\\rho=%s$, $\\eta=%s$, $\\alpha_{0.05}=%s$, and $B=%s$.",
      fmt_num(best_bending$gmds_stress[[1L]]),
      fmt_num(best_bending$procrustes_rmse[[1L]]),
      fmt_num(best_bending$roughness[[1L]]),
      fmt_num(best_bending$area_q05_ratio[[1L]]),
      fmt_num(best_bending$bend_raw_penalty[[1L]]),
      fmt_num(best_hybrid$gmds_stress[[1L]]),
      fmt_num(best_hybrid$procrustes_rmse[[1L]]),
      fmt_num(best_hybrid$roughness[[1L]]),
      fmt_num(best_hybrid$area_q05_ratio[[1L]]),
      fmt_num(best_hybrid$bend_raw_penalty[[1L]])
    ),
    "Visually, however, the selected bending-only and hybrid layouts remain in the same folded basin as untethered GMDS. The second-order term changes the diagnostics only modestly and does not restore a smooth paraboloid.",
    "\\begin{figure}[H]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "paraboloid_bending_selected.png"))
    ),
    "\\caption{Reference surface, MDS, untethered GMDS, anchor-only GMDS, best bending-only GMDS, and best anchor-plus-bending GMDS.}",
    "\\end{figure}",
    "\\begin{table}[H]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\caption{Selected 12x12 paraboloid methods. Lower $\\sigma$, $\\rho$, $\\eta$, and $B$ are better; larger $\\alpha_{0.05}$ is better.}",
    "\\begin{tabular}{@{}lrrrrrrrr@{}}",
    "\\toprule",
    "Setting & $\\lambda$ & $\\beta$ & $t$ (s) & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ & $B$ \\\\",
    "\\midrule",
    write_method_table(selected_df),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    "\\section{Experiment 1: Bending Only}",
    "The first sweep tests the second-order penalty by itself. We keep $\\lambda = 0$ and vary the bending weight $\\beta$ across a moderate range. This isolates the effect of the curvature penalty from the MDS tether.",
    sprintf(
      "Across this sweep, the best bending-only tradeoff was %s, with $\\sigma=%s$, $\\rho=%s$, $\\eta=%s$, $\\alpha_{0.05}=%s$, and $B=%s$.",
      tex_escape(best_bending$setting_label[[1L]]),
      fmt_num(best_bending$gmds_stress[[1L]]),
      fmt_num(best_bending$procrustes_rmse[[1L]]),
      fmt_num(best_bending$roughness[[1L]]),
      fmt_num(best_bending$area_q05_ratio[[1L]]),
      fmt_num(best_bending$bend_raw_penalty[[1L]])
    ),
    "The key question here is whether second-order regularization alone can keep the correction inside a smooth paraboloid-like basin, or whether it only softens the same singular endpoint that untethered GMDS already prefers.",
    "In practice the bending-only sweep is almost flat numerically: from $\\beta=0$ to $\\beta=1$, the stress changes only from $0.0298$ to $0.0298$ at the displayed precision, while $\\alpha_{0.05}$ moves only from $0.3604$ to $0.3607$. The figure matches that story closely: the panels are visually almost indistinguishable, and the singularity remains.",
    "\\begin{figure}[H]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "paraboloid_bending_only_sweep.png"))
    ),
    "\\caption{Bending-only sweep with $\\lambda = 0$. The first two panels show the reference surface and the MDS initialization.}",
    "\\end{figure}",
    "\\begin{table}[H]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\caption{Bending-only sweep on the 12x12 paraboloid.}",
    "\\begin{tabular}{@{}rrrrrrrr@{}}",
    "\\toprule",
    "$\\beta$ & $t$ (s) & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ & $B$ & mean step \\\\",
    "\\midrule",
    write_sweep_table(bending_df),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    "\\section{Experiment 2: Anchor Plus Bending}",
    "The second sweep keeps the quadratic MDS tether active with fixed $\\lambda = 0.1$ and varies the bending weight $\\beta$. This is the most direct implementation of the working hypothesis from the previous report: preserve the smooth global MDS shape with the anchor while using the second-order penalty to oppose concentrated folds.",
    sprintf(
      "The best hybrid run in this sweep was %s, with $\\sigma=%s$, $\\rho=%s$, $\\eta=%s$, $\\alpha_{0.05}=%s$, and $B=%s$.",
      tex_escape(best_hybrid$setting_label[[1L]]),
      fmt_num(best_hybrid$gmds_stress[[1L]]),
      fmt_num(best_hybrid$procrustes_rmse[[1L]]),
      fmt_num(best_hybrid$roughness[[1L]]),
      fmt_num(best_hybrid$area_q05_ratio[[1L]]),
      fmt_num(best_hybrid$bend_raw_penalty[[1L]])
    ),
    "This is the decisive sweep for the present round, because it tests whether the second-order term is genuinely stronger than the earlier first-order smoothness term when both are given the same MDS anchor support.",
    "The hybrid sweep helps more than bending alone, but not in the decisive way we wanted. Relative to anchor-only GMDS, the strongest tested weight $\\beta=2$ improves $\\rho$ from $0.3224$ to $0.3218$ and $\\alpha_{0.05}$ from $0.5547$ to $0.5795$, but the fold remains visible in every panel.",
    "\\begin{figure}[H]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "paraboloid_bending_hybrid_sweep.png"))
    ),
    "\\caption{Anchor-plus-bending sweep with fixed $\\lambda = 0.1$. The second panel is anchor-only GMDS with $\\beta = 0$.}",
    "\\end{figure}",
    "\\begin{table}[H]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\caption{Hybrid sweep on the 12x12 paraboloid with fixed $\\lambda = 0.1$.}",
    "\\begin{tabular}{@{}rrrrrrrr@{}}",
    "\\toprule",
    "$\\beta$ & $t$ (s) & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ & $B$ & mean step \\\\",
    "\\midrule",
    write_sweep_table(hybrid_df),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    "\\section{Experiment 3: Trajectory Under the Best Hybrid}",
    sprintf(
      "To see whether the second-order penalty changes the optimization path itself, not only the endpoint, we compare the untethered trajectory against the best hybrid trajectory. The selected hybrid was %s and took %ss in total.",
      tex_escape(best_hybrid$setting_label[[1L]]),
      fmt_time(regularized_runtime)
    ),
    "The trace shows that the collapse happens almost immediately. After one GMDS correction step, the hybrid trajectory jumps from $\\rho=0.2077$ and $\\alpha_{0.05}=0.7754$ at the MDS start to $\\rho=0.4069$ and $\\alpha_{0.05}=0.1785$, and by iteration 2 the singular basin is already established. The bending penalty then only smooths that basin slightly rather than preventing it.",
    "\\begin{figure}[H]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=0.95\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "paraboloid_bending_trace_compare.png"))
    ),
    "\\caption{Untethered versus anchor-plus-bending trajectory metrics over 25 GMDS correction steps.}",
    "\\end{figure}",
    "\\begin{landscape}",
    "\\begin{figure}[p]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "paraboloid_bending_trajectory.png"))
    ),
    "\\caption{Reference surface plus selected checkpoints from the best anchor-plus-bending trajectory.}",
    "\\end{figure}",
    "\\end{landscape}",
    "\\begin{table}[H]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\caption{Selected checkpoints from the best hybrid trajectory.}",
    "\\begin{tabular}{@{}rrrrrrr@{}}",
    "\\toprule",
    "Iter & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ & $B$ & step \\\\",
    "\\midrule",
    write_trajectory_table(selected_traj),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    "\\section{Takeaways}",
    "The concrete conclusion from this round is negative but informative. The simple isotropic second-order bending penalty does not fix the paraboloid pathology on the rectangular mesh. Bending-only GMDS is almost indistinguishable from untethered GMDS, and anchor-plus-bending improves the diagnostics only modestly while leaving the same singular fold in place.",
    "That points to a stronger next step: use a regularizer or constraint that reacts to fold-over more directly, for example area-orientation barriers on mesh cells, signed-volume or self-intersection barriers, or a curvature penalty coupled to the reference MDS geometry instead of isotropic row-and-column second differences alone.",
    "",
    "\\end{document}"
  )

  writeLines(tex, tex_path)
}

case <- new_case(12L, 25L)
cmd_res <- run_cmdscale(case)
gmds_res <- run_config(case, "gmds_average", "GMDS avg")
anchor_res <- run_config(
  case,
  "gmds_anchor_0.1",
  "GMDS anchor-only (lambda = 0.1)",
  anchor_mode = "cmdscale",
  anchor_weight = 0.1,
  anchor_weight_end = 0.1,
  continuation = "constant"
)

bending_values <- c(0, 0.01, 0.05, 0.1, 0.2, 0.5, 1)
hybrid_values <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 2)
bending_sweep <- run_beta_sweep(case, beta_values = bending_values, anchor_lambda = 0)
hybrid_sweep <- run_beta_sweep(case, beta_values = hybrid_values, anchor_lambda = 0.1)

cmd_sigma <- cmd_res$row$gmds_stress[[1L]]
best_bending <- pick_best_tradeoff(bending_sweep$metrics, cmd_sigma = cmd_sigma)
best_hybrid <- pick_best_tradeoff(hybrid_sweep$metrics, cmd_sigma = cmd_sigma)

best_bending_fit <- bending_sweep$fits[[match(best_bending$setting_id, vapply(bending_sweep$fits, function(x) x$row$setting_id[[1L]], character(1L)))]]
best_hybrid_fit <- hybrid_sweep$fits[[match(best_hybrid$setting_id, vapply(hybrid_sweep$fits, function(x) x$row$setting_id[[1L]], character(1L)))]]

selected_target_cols <- unique(c(
  names(cmd_res$row),
  names(best_bending),
  names(best_hybrid)
))
selected_rows <- list(
  cmd_res$row,
  gmds_res$row,
  anchor_res$row,
  best_bending,
  best_hybrid
)
selected_df <- do.call(rbind, lapply(selected_rows, function(row) {
  missing_cols <- setdiff(selected_target_cols, names(row))
  for (nm in missing_cols) {
    row[[nm]] <- NA_real_
  }
  row[, selected_target_cols, drop = FALSE]
}))

utils::write.csv(selected_df, file.path(tmp_dir, "paraboloid_bending_selected.csv"), row.names = FALSE)
utils::write.csv(bending_sweep$metrics, file.path(tmp_dir, "paraboloid_bending_only_sweep.csv"), row.names = FALSE)
utils::write.csv(hybrid_sweep$metrics, file.path(tmp_dir, "paraboloid_bending_hybrid_sweep.csv"), row.names = FALSE)

untethered_traj <- trajectory_metrics(case, gmds_res$fit, "untethered", "untethered")
regularized_traj <- trajectory_metrics(case, best_hybrid_fit$fit, "regularized", "regularized")
utils::write.csv(untethered_traj, file.path(tmp_dir, "paraboloid_bending_untethered_trajectory.csv"), row.names = FALSE)
utils::write.csv(regularized_traj, file.path(tmp_dir, "paraboloid_bending_hybrid_trajectory.csv"), row.names = FALSE)

selected_titles <- c(
  "Reference surface",
  sprintf(
    "cmdscale\nsigma %s, rho %s\neta %s, alpha_0.05 %s\nt %ss",
    fmt_num(cmd_res$row$gmds_stress[[1L]]),
    fmt_num(cmd_res$row$procrustes_rmse[[1L]]),
    fmt_num(cmd_res$row$roughness[[1L]]),
    fmt_num(cmd_res$row$area_q05_ratio[[1L]]),
    fmt_time(cmd_res$row$elapsed_sec[[1L]])
  ),
  sprintf(
    "GMDS avg\nsigma %s, rho %s\neta %s, alpha_0.05 %s\nt %ss",
    fmt_num(gmds_res$row$gmds_stress[[1L]]),
    fmt_num(gmds_res$row$procrustes_rmse[[1L]]),
    fmt_num(gmds_res$row$roughness[[1L]]),
    fmt_num(gmds_res$row$area_q05_ratio[[1L]]),
    fmt_time(gmds_res$row$elapsed_sec[[1L]])
  ),
  sprintf(
    "anchor only\nsigma %s, rho %s\neta %s, alpha_0.05 %s\nt %ss",
    fmt_num(anchor_res$row$gmds_stress[[1L]]),
    fmt_num(anchor_res$row$procrustes_rmse[[1L]]),
    fmt_num(anchor_res$row$roughness[[1L]]),
    fmt_num(anchor_res$row$area_q05_ratio[[1L]]),
    fmt_time(anchor_res$row$elapsed_sec[[1L]])
  ),
  sprintf(
    "%s\nsigma %s, rho %s\neta %s, alpha_0.05 %s\nt %ss",
    best_bending$setting_label[[1L]],
    fmt_num(best_bending$gmds_stress[[1L]]),
    fmt_num(best_bending$procrustes_rmse[[1L]]),
    fmt_num(best_bending$roughness[[1L]]),
    fmt_num(best_bending$area_q05_ratio[[1L]]),
    fmt_time(best_bending$elapsed_sec[[1L]])
  ),
  sprintf(
    "%s\nsigma %s, rho %s\neta %s, alpha_0.05 %s\nt %ss",
    best_hybrid$setting_label[[1L]],
    fmt_num(best_hybrid$gmds_stress[[1L]]),
    fmt_num(best_hybrid$procrustes_rmse[[1L]]),
    fmt_num(best_hybrid$roughness[[1L]]),
    fmt_num(best_hybrid$area_q05_ratio[[1L]]),
    fmt_time(best_hybrid$elapsed_sec[[1L]])
  )
)
save_surface_panel_grid(
  case = case,
  coords_list = list(
    case$truth,
    cmd_res$coords,
    gmds_res$fit$coords,
    anchor_res$fit$coords,
    best_bending_fit$fit$coords,
    best_hybrid_fit$fit$coords
  ),
  title_list = selected_titles,
  path = file.path(pdf_dir, "paraboloid_bending_selected.png"),
  nrow = 2L,
  ncol = 3L,
  heading = "Paraboloid 12x12 selected bending-regularized methods",
  truth_index = 1L
)

bending_titles <- c(
  "Reference surface",
  sprintf(
    "cmdscale\nsigma %s, rho %s\neta %s, alpha_0.05 %s",
    fmt_num(cmd_res$row$gmds_stress[[1L]]),
    fmt_num(cmd_res$row$procrustes_rmse[[1L]]),
    fmt_num(cmd_res$row$roughness[[1L]]),
    fmt_num(cmd_res$row$area_q05_ratio[[1L]])
  ),
  vapply(seq_len(nrow(bending_sweep$metrics)), function(i) {
    sprintf(
      "beta = %s\nsigma %s, rho %s\neta %s, alpha_0.05 %s\nt %ss",
      fmt_num(bending_sweep$metrics$bend_weight[[i]], digits = 3L),
      fmt_num(bending_sweep$metrics$gmds_stress[[i]]),
      fmt_num(bending_sweep$metrics$procrustes_rmse[[i]]),
      fmt_num(bending_sweep$metrics$roughness[[i]]),
      fmt_num(bending_sweep$metrics$area_q05_ratio[[i]]),
      fmt_time(bending_sweep$metrics$elapsed_sec[[i]])
    )
  }, character(1L))
)
save_surface_panel_grid(
  case = case,
  coords_list = c(
    list(case$truth, cmd_res$coords),
    lapply(bending_sweep$fits, function(x) x$fit$coords)
  ),
  title_list = bending_titles,
  path = file.path(pdf_dir, "paraboloid_bending_only_sweep.png"),
  nrow = 3L,
  ncol = 3L,
  heading = "Paraboloid 12x12 bending-only sweep",
  truth_index = 1L
)

hybrid_titles <- c(
  "Reference surface",
  sprintf(
    "anchor only\nsigma %s, rho %s\neta %s, alpha_0.05 %s\nt %ss",
    fmt_num(anchor_res$row$gmds_stress[[1L]]),
    fmt_num(anchor_res$row$procrustes_rmse[[1L]]),
    fmt_num(anchor_res$row$roughness[[1L]]),
    fmt_num(anchor_res$row$area_q05_ratio[[1L]]),
    fmt_time(anchor_res$row$elapsed_sec[[1L]])
  ),
  vapply(seq_len(nrow(hybrid_sweep$metrics)), function(i) {
    sprintf(
      "beta = %s\nsigma %s, rho %s\neta %s, alpha_0.05 %s\nt %ss",
      fmt_num(hybrid_sweep$metrics$bend_weight[[i]], digits = 3L),
      fmt_num(hybrid_sweep$metrics$gmds_stress[[i]]),
      fmt_num(hybrid_sweep$metrics$procrustes_rmse[[i]]),
      fmt_num(hybrid_sweep$metrics$roughness[[i]]),
      fmt_num(hybrid_sweep$metrics$area_q05_ratio[[i]]),
      fmt_time(hybrid_sweep$metrics$elapsed_sec[[i]])
    )
  }, character(1L))
)
save_surface_panel_grid(
  case = case,
  coords_list = c(
    list(case$truth, anchor_res$fit$coords),
    lapply(hybrid_sweep$fits, function(x) x$fit$coords)
  ),
  title_list = hybrid_titles,
  path = file.path(pdf_dir, "paraboloid_bending_hybrid_sweep.png"),
  nrow = 3L,
  ncol = 3L,
  heading = "Paraboloid 12x12 anchor-plus-bending sweep",
  truth_index = 1L
)

save_trace_compare_plot(
  untethered_df = untethered_traj,
  regularized_df = regularized_traj,
  path = file.path(pdf_dir, "paraboloid_bending_trace_compare.png")
)

selected_iterations <- c(0L, 1L, 2L, 5L, 10L, 15L, 20L, 25L)
trajectory_titles <- c(
  "Reference surface",
  vapply(selected_iterations, function(iteration) {
    row <- regularized_traj[regularized_traj$iteration == iteration, , drop = FALSE]
    sprintf(
      "iter %d\nsigma %s, rho %s\neta %s, alpha_0.05 %s",
      iteration,
      fmt_num(row$gmds_stress[[1L]]),
      fmt_num(row$procrustes_rmse[[1L]]),
      fmt_num(row$roughness[[1L]]),
      fmt_num(row$area_q05_ratio[[1L]])
    )
  }, character(1L))
)
trajectory_coords <- c(
  list(case$truth),
  lapply(selected_iterations + 1L, function(i) best_hybrid_fit$fit$frames[[i]])
)
save_surface_panel_grid(
  case = case,
  coords_list = trajectory_coords,
  title_list = trajectory_titles,
  path = file.path(pdf_dir, "paraboloid_bending_trajectory.png"),
  nrow = 3L,
  ncol = 3L,
  heading = "Paraboloid 12x12 anchor-plus-bending trajectory checkpoints",
  truth_index = 1L
)

tex_path <- file.path(manual_root, "reports", "paraboloid_gmds_bending_report_2026-03-31.tex")
write_report(
  case = case,
  selected_df = selected_df,
  bending_df = bending_sweep$metrics,
  hybrid_df = hybrid_sweep$metrics,
  best_bending = best_bending,
  best_hybrid = best_hybrid,
  untethered_traj = untethered_traj,
  regularized_traj = regularized_traj,
  regularized_runtime = best_hybrid$elapsed_sec[[1L]],
  tex_path = tex_path
)

report_pdf_path <- sub("\\.tex$", ".pdf", tex_path)
cat(sprintf(
  paste(
    "Paraboloid GMDS bending benchmark complete.",
    "Selected methods: %s",
    "Bending-only sweep: %s",
    "Hybrid sweep: %s",
    "Untethered trajectory: %s",
    "Hybrid trajectory: %s",
    "Report (tex): %s",
    "Report (pdf target): %s",
    sep = "\n"
  ),
  file.path(tmp_dir, "paraboloid_bending_selected.csv"),
  file.path(tmp_dir, "paraboloid_bending_only_sweep.csv"),
  file.path(tmp_dir, "paraboloid_bending_hybrid_sweep.csv"),
  file.path(tmp_dir, "paraboloid_bending_untethered_trajectory.csv"),
  file.path(tmp_dir, "paraboloid_bending_hybrid_trajectory.csv"),
  tex_path,
  report_pdf_path
))
cat("\n")
if (nzchar(Sys.which("latexmk"))) {
  cat(sprintf(
    "Compile the PDF with: (cd %s && latexmk -pdf %s)\n",
    file.path(manual_root, "reports"),
    basename(tex_path)
  ))
}
