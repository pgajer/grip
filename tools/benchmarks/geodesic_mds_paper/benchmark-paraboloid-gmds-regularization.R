#!/usr/bin/env Rscript

run_tag <- "paraboloid-gmds-regularization-2026-03-31"
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
  stop("Install 'devtools' or the 'grip' package to run the paraboloid regularization benchmark.")
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
  prepared <- grip.prepare.geodesic.kk(
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
    cmd_elapsed = cmd.elapsed
  )
}

evaluate_layout <- function(case,
                            coords,
                            setting_id,
                            setting_label,
                            elapsed_sec = NA_real_,
                            fit = NULL) {
  anchor.coords <- if (is.null(fit)) NULL else fit$anchor_coords
  anchor.weight <- if (is.null(fit)) 0 else fit$final_anchor_weight
  smooth.weight <- if (is.null(fit)) 0 else fit$final_smoothness_weight
  score <- grip.score.geodesic.mds(
    coords,
    prepared = case$prepared,
    anchor_coords = anchor.coords,
    anchor_weight = anchor.weight,
    smoothness_weight = smooth.weight
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
    smoothness_weight = smooth.weight,
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

run_mu_sweep <- function(case,
                         mu_values,
                         anchor_lambda = 0) {
  rows <- list()
  fits <- vector("list", length(mu_values))
  for (i in seq_along(mu_values)) {
    mu <- mu_values[[i]]
    label <- if (anchor_lambda > 0) {
      sprintf("lambda = %.3f, mu = %.3f", anchor_lambda, mu)
    } else {
      sprintf("mu = %.3f", mu)
    }
    setting_id <- if (anchor_lambda > 0) {
      sprintf("hybrid_lambda_%s_mu_%s", format(anchor_lambda, trim = TRUE), format(mu, trim = TRUE))
    } else {
      sprintf("smooth_mu_%s", format(mu, trim = TRUE))
    }
    fit <- run_config(
      case = case,
      setting_id = setting_id,
      setting_label = label,
      smoothness_weight = mu,
      anchor_mode = if (anchor_lambda > 0) "cmdscale" else "none",
      anchor_weight = anchor_lambda,
      anchor_weight_end = anchor_lambda,
      continuation = "constant"
    )
    fit$row$mu <- mu
    fit$row$lambda <- anchor_lambda
    rows[[i]] <- fit$row
    fits[[i]] <- fit
  }
  list(
    metrics = do.call(rbind, rows),
    fits = fits
  )
}

pick_best_tradeoff <- function(df,
                               cmd_sigma) {
  candidates <- df[df$smoothness_weight > 0 & df$gmds_stress < cmd_sigma, , drop = FALSE]
  if (nrow(candidates) == 0L) {
    candidates <- df[df$smoothness_weight > 0, , drop = FALSE]
  }
  candidates[order(candidates$roughness,
                   -candidates$area_q05_ratio,
                   candidates$gmds_stress,
                   candidates$elapsed_sec), , drop = FALSE][1L, , drop = FALSE]
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
    grip.plot(
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
      legend = c("untethered", "anchor + smooth"),
      col = c("#b23a48", "#2d6a4f"),
      pch = c(16, 17),
      bty = "n",
      cex = 0.9
    )
  }

  plot_pair(untethered_df$gmds_stress, regularized_df$gmds_stress, expression(sigma ~ "over iterations"), expression(sigma))
  plot_pair(untethered_df$roughness, regularized_df$roughness, expression(eta ~ "over iterations"), expression(eta))
  plot_pair(untethered_df$area_q05_ratio, regularized_df$area_q05_ratio, expression(alpha[0.05] ~ "over iterations"), expression(alpha[0.05]))
  plot_pair(untethered_df$accepted_step, regularized_df$accepted_step, "Accepted step over iterations", "step")
}

write_method_table <- function(df) {
  rows <- apply(df, 1L, function(row) {
    paste(
      tex_escape(row[["setting_label"]]),
      fmt_num(as.numeric(row[["anchor_weight"]]), digits = 3L),
      fmt_num(as.numeric(row[["smoothness_weight"]]), digits = 3L),
      fmt_time(as.numeric(row[["elapsed_sec"]])),
      fmt_num(as.numeric(row[["gmds_stress"]])),
      fmt_num(as.numeric(row[["procrustes_rmse"]])),
      fmt_num(as.numeric(row[["roughness"]])),
      fmt_num(as.numeric(row[["area_q05_ratio"]])),
      sep = " & "
    )
  })
  paste0(rows, " \\\\")
}

write_sweep_table <- function(df) {
  rows <- apply(df, 1L, function(row) {
    paste(
      fmt_num(as.numeric(row[["smoothness_weight"]]), digits = 3L),
      fmt_time(as.numeric(row[["elapsed_sec"]])),
      fmt_num(as.numeric(row[["gmds_stress"]])),
      fmt_num(as.numeric(row[["procrustes_rmse"]])),
      fmt_num(as.numeric(row[["roughness"]])),
      fmt_num(as.numeric(row[["area_q05_ratio"]])),
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
      fmt_time(as.numeric(row[["accepted_step"]])),
      sep = " & "
    )
  })
  paste0(rows, " \\\\")
}

write_report <- function(case,
                         selected_df,
                         smooth_df,
                         hybrid_df,
                         best_smooth,
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
    "\\title{Paraboloid GMDS Regularization Report}",
    "\\author{Focused 12x12 paraboloid study}",
    "\\date{March 31, 2026}",
    "",
    "\\begin{document}",
    "\\maketitle",
    "",
    "\\begin{abstract}",
    "This report evaluates whether the unsettling singularities in the 3D paraboloid GMDS layouts are best treated as a pure step-size issue or as a missing regularization issue. We keep the same tie-averaged GMDS objective and add one new local-smoothness penalty, then compare untethered GMDS, anchor-only GMDS, smoothness-only GMDS, and the combined anchor-plus-smoothness model on the same $12\\times12$ orthogonal paraboloid mesh. The main outcome is negative but useful: this first-order smoothness term improves the diagnostics only modestly and does not remove the singularity.",
    "\\end{abstract}",
    "",
    "\\section{Augmented Objective}",
    "All GMDS runs in this report start from the same classical-MDS embedding of the graph-geodesic distance matrix and then run 25 extra GMDS iterations with Armijo backtracking. The new regularized objective is",
    "\\[",
    "E_{\\lambda,\\mu}(Z) = E_{\\mathrm{GMDS}}(Z) + \\lambda \\lVert Z - A \\rVert_F^2 + \\mu \\sum_{i=1}^n \\left\\lVert z_i - \\frac{1}{|N(i)|}\\sum_{j \\in N(i)} z_j \\right\\rVert^2,",
    "\\]",
    "where $A$ is the classical-MDS anchor and the new $\\mu$ term penalizes local departures from the neighborhood average. In words, $\\lambda$ keeps the embedding near the smooth global MDS shape, while $\\mu$ discourages curvature from concentrating into sharp local singularities.",
    "",
    "\\section{Selected Methods}",
    sprintf(
      "The selected comparison below shows the reference surface, classical MDS, untethered GMDS, anchor-only GMDS with $\\lambda=0.1$, the best smoothness-only run from the $\\mu$ sweep, and the best hybrid run from the $(\\lambda,\\mu)$ sweep. In this study, the best smoothness-only run was %s and the best hybrid run was %s.",
      tex_escape(best_smooth$setting_label[[1L]]),
      tex_escape(best_hybrid$setting_label[[1L]])
    ),
    "The selected figure already shows the main limitation of this regularizer: both the best smoothness-only and best hybrid layouts still sit in the same visually singular basin as the untethered GMDS solution.",
    "\\begin{figure}[H]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "paraboloid_regularized_selected.png"))
    ),
    "\\caption{Reference surface, MDS, untethered GMDS, anchor-only GMDS, best smoothness-only GMDS, and best anchor-plus-smoothness GMDS.}",
    "\\end{figure}",
    "\\begin{table}[H]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\caption{Selected 12x12 paraboloid methods. Lower $\\sigma$, $\\rho$, and $\\eta$ are better; larger $\\alpha_{0.05}$ is better.}",
    "\\begin{tabular}{@{}lrrrrrrr@{}}",
    "\\toprule",
    "Setting & $\\lambda$ & $\\mu$ & $t$ (s) & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ \\\\",
    "\\midrule",
    write_method_table(selected_df),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    "\\section{Experiment 1: Smoothness Only}",
    "The first sweep tests whether an explicit local smoothness penalty alone can suppress the singularity, without any MDS tether beyond the common initialization. We therefore keep $\\lambda = 0$ and vary $\\mu$ across a moderate range.",
    sprintf(
      "Across this sweep, the best smoothness-only tradeoff was %s, with $\\sigma=%s$, $\\rho=%s$, $\\eta=%s$, and $\\alpha_{0.05}=%s$.",
      tex_escape(best_smooth$setting_label[[1L]]),
      fmt_num(best_smooth$gmds_stress[[1L]]),
      fmt_num(best_smooth$procrustes_rmse[[1L]]),
      fmt_num(best_smooth$roughness[[1L]]),
      fmt_num(best_smooth$area_q05_ratio[[1L]])
    ),
    "Numerically, the roughness proxy changes only slightly across the sweep, and visually the singular fold remains present throughout. That means the first-order smoothness term is not strong enough, at least in this form, to pull the optimizer out of the bad basin.",
    "\\begin{figure}[H]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "paraboloid_regularized_smooth_sweep.png"))
    ),
    "\\caption{Smoothness-only sweep with $\\lambda = 0$. The first two panels show the reference surface and the MDS initialization.}",
    "\\end{figure}",
    "\\begin{table}[H]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\caption{Smoothness-only sweep on the 12x12 paraboloid.}",
    "\\begin{tabular}{@{}rrrrrrr@{}}",
    "\\toprule",
    "$\\mu$ & $t$ (s) & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ & mean step \\\\",
    "\\midrule",
    write_sweep_table(smooth_df),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    "\\section{Experiment 2: Anchor Plus Smoothness}",
    "The second sweep keeps the quadratic MDS tether active with $\\lambda = 0.1$ and varies the local smoothness weight $\\mu$. This directly tests the combined recommendation from the earlier pathology report: a global anchor to keep the overall shape smooth plus an explicit local regularizer to prevent concentrated folds.",
    sprintf(
      "The best hybrid run in this sweep was %s, with $\\sigma=%s$, $\\rho=%s$, $\\eta=%s$, and $\\alpha_{0.05}=%s$.",
      tex_escape(best_hybrid$setting_label[[1L]]),
      fmt_num(best_hybrid$gmds_stress[[1L]]),
      fmt_num(best_hybrid$procrustes_rmse[[1L]]),
      fmt_num(best_hybrid$roughness[[1L]]),
      fmt_num(best_hybrid$area_q05_ratio[[1L]])
    ),
    "This hybrid improves the area-floor diagnostic substantially relative to untethered GMDS and lowers roughness modestly, but the figure shows that it still does not restore a smooth paraboloid. So the MDS anchor and the present first-order smoothness term are not yet enough.",
    "\\begin{figure}[H]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "paraboloid_regularized_hybrid_sweep.png"))
    ),
    "\\caption{Anchor-plus-smoothness sweep with fixed $\\lambda = 0.1$. The second panel is anchor-only GMDS with $\\mu = 0$.}",
    "\\end{figure}",
    "\\begin{table}[H]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\caption{Hybrid sweep on the 12x12 paraboloid with fixed $\\lambda = 0.1$.}",
    "\\begin{tabular}{@{}rrrrrrr@{}}",
    "\\toprule",
    "$\\mu$ & $t$ (s) & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ & mean step \\\\",
    "\\midrule",
    write_sweep_table(hybrid_df),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    "\\section{Experiment 3: Regularized Trajectory}",
    sprintf(
    "To see whether the hybrid objective changes the optimization path itself, not only the final endpoint, we compare the untethered trajectory against the best hybrid trajectory. The regularized run used %s and took %ss in total.",
      tex_escape(best_hybrid$setting_label[[1L]]),
      fmt_time(regularized_runtime)
    ),
    "If the singularity is mostly objective-driven, the regularized trace should stay smoother throughout the 25-step correction, not only at the very end. That is exactly the behavior we want to check here.",
    "What the trace comparison actually shows is a partial but incomplete mitigation: the hybrid run has better area-floor behavior and slightly lower roughness than the untethered run, but it still converges to a singular-looking surface.",
    "\\begin{figure}[H]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=0.95\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "paraboloid_regularized_trace_compare.png"))
    ),
    "\\caption{Untethered versus regularized trajectory metrics over 25 GMDS correction steps.}",
    "\\end{figure}",
    "\\begin{landscape}",
    "\\begin{figure}[p]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "paraboloid_regularized_trajectory.png"))
    ),
    "\\caption{Reference surface plus selected checkpoints from the best hybrid trajectory.}",
    "\\end{figure}",
    "\\end{landscape}",
    "\\begin{table}[H]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{6pt}",
    "\\caption{Selected checkpoints from the best hybrid trajectory.}",
    "\\begin{tabular}{@{}rrrrrr@{}}",
    "\\toprule",
    "Iter & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ & step \\\\",
    "\\midrule",
    write_trajectory_table(selected_traj),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    "\\section{Takeaways}",
    "The new experiments separate three possibilities cleanly. If step size were the main problem, a local smoothness term would not systematically improve the trajectory. If the singularity were caused only by lack of global anchoring, then anchor-only GMDS would already be enough. Here we see something more specific: anchor-only GMDS helps a bit, and anchor plus first-order smoothness helps a bit more, but none of these variants removes the singularity.",
    "So the practical conclusion is that the pathology is not just a step-size issue, and the present first-order neighborhood-average penalty is not the right regularizer by itself. The natural next experiment is a second-order mesh bending penalty, for example a row-and-column discrete curvature term on the rectangular grid, because that targets sharp creases more directly than a first-order smoothing term does.",
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

mu_values <- c(0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2)
smooth_sweep <- run_mu_sweep(case, mu_values = mu_values, anchor_lambda = 0)
hybrid_sweep <- run_mu_sweep(case, mu_values = mu_values, anchor_lambda = 0.1)

cmd_sigma <- cmd_res$row$gmds_stress[[1L]]
best_smooth <- pick_best_tradeoff(smooth_sweep$metrics, cmd_sigma = cmd_sigma)
best_hybrid <- pick_best_tradeoff(hybrid_sweep$metrics, cmd_sigma = cmd_sigma)

best_smooth_fit <- smooth_sweep$fits[[match(best_smooth$setting_id, vapply(smooth_sweep$fits, function(x) x$row$setting_id[[1L]], character(1L)))]]
best_hybrid_fit <- hybrid_sweep$fits[[match(best_hybrid$setting_id, vapply(hybrid_sweep$fits, function(x) x$row$setting_id[[1L]], character(1L)))]]

selected_target_cols <- unique(c(
  names(cmd_res$row),
  names(best_smooth),
  names(best_hybrid)
))
selected_rows <- list(
  cmd_res$row,
  gmds_res$row,
  anchor_res$row,
  best_smooth,
  best_hybrid
)
selected_df <- do.call(rbind, lapply(selected_rows, function(row) {
  missing_cols <- setdiff(selected_target_cols, names(row))
  for (nm in missing_cols) {
    row[[nm]] <- NA_real_
  }
  row[, selected_target_cols, drop = FALSE]
}))
utils::write.csv(selected_df, file.path(tmp_dir, "paraboloid_regularized_selected.csv"), row.names = FALSE)
utils::write.csv(smooth_sweep$metrics, file.path(tmp_dir, "paraboloid_regularized_smooth_sweep.csv"), row.names = FALSE)
utils::write.csv(hybrid_sweep$metrics, file.path(tmp_dir, "paraboloid_regularized_hybrid_sweep.csv"), row.names = FALSE)

untethered_traj <- trajectory_metrics(case, gmds_res$fit, "untethered", "untethered")
regularized_traj <- trajectory_metrics(case, best_hybrid_fit$fit, "regularized", "regularized")
utils::write.csv(untethered_traj, file.path(tmp_dir, "paraboloid_regularized_untethered_trajectory.csv"), row.names = FALSE)
utils::write.csv(regularized_traj, file.path(tmp_dir, "paraboloid_regularized_hybrid_trajectory.csv"), row.names = FALSE)

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
    best_smooth$setting_label[[1L]],
    fmt_num(best_smooth$gmds_stress[[1L]]),
    fmt_num(best_smooth$procrustes_rmse[[1L]]),
    fmt_num(best_smooth$roughness[[1L]]),
    fmt_num(best_smooth$area_q05_ratio[[1L]]),
    fmt_time(best_smooth$elapsed_sec[[1L]])
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
    best_smooth_fit$fit$coords,
    best_hybrid_fit$fit$coords
  ),
  title_list = selected_titles,
  path = file.path(pdf_dir, "paraboloid_regularized_selected.png"),
  nrow = 2L,
  ncol = 3L,
  heading = "Paraboloid 12x12 selected regularized methods",
  truth_index = 1L
)

smooth_titles <- c(
  "Reference surface",
  sprintf(
    "cmdscale\nsigma %s, rho %s\neta %s, alpha_0.05 %s",
    fmt_num(cmd_res$row$gmds_stress[[1L]]),
    fmt_num(cmd_res$row$procrustes_rmse[[1L]]),
    fmt_num(cmd_res$row$roughness[[1L]]),
    fmt_num(cmd_res$row$area_q05_ratio[[1L]])
  ),
  vapply(seq_len(nrow(smooth_sweep$metrics)), function(i) {
    sprintf(
      "mu = %s\nsigma %s, rho %s\neta %s, alpha_0.05 %s\nt %ss",
      fmt_num(smooth_sweep$metrics$smoothness_weight[[i]], digits = 3L),
      fmt_num(smooth_sweep$metrics$gmds_stress[[i]]),
      fmt_num(smooth_sweep$metrics$procrustes_rmse[[i]]),
      fmt_num(smooth_sweep$metrics$roughness[[i]]),
      fmt_num(smooth_sweep$metrics$area_q05_ratio[[i]]),
      fmt_time(smooth_sweep$metrics$elapsed_sec[[i]])
    )
  }, character(1L))
)
save_surface_panel_grid(
  case = case,
  coords_list = c(
    list(case$truth, cmd_res$coords),
    lapply(smooth_sweep$fits, function(x) x$fit$coords)
  ),
  title_list = smooth_titles,
  path = file.path(pdf_dir, "paraboloid_regularized_smooth_sweep.png"),
  nrow = 3L,
  ncol = 3L,
  heading = "Paraboloid 12x12 smoothness-only sweep",
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
      "mu = %s\nsigma %s, rho %s\neta %s, alpha_0.05 %s\nt %ss",
      fmt_num(hybrid_sweep$metrics$smoothness_weight[[i]], digits = 3L),
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
  path = file.path(pdf_dir, "paraboloid_regularized_hybrid_sweep.png"),
  nrow = 3L,
  ncol = 3L,
  heading = "Paraboloid 12x12 anchor-plus-smoothness sweep",
  truth_index = 1L
)

save_trace_compare_plot(
  untethered_df = untethered_traj,
  regularized_df = regularized_traj,
  path = file.path(pdf_dir, "paraboloid_regularized_trace_compare.png")
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
  path = file.path(pdf_dir, "paraboloid_regularized_trajectory.png"),
  nrow = 3L,
  ncol = 3L,
  heading = "Paraboloid 12x12 regularized trajectory checkpoints",
  truth_index = 1L
)

tex_path <- file.path(manual_root, "reports", "paraboloid_gmds_regularization_report_2026-03-31.tex")
write_report(
  case = case,
  selected_df = selected_df,
  smooth_df = smooth_sweep$metrics,
  hybrid_df = hybrid_sweep$metrics,
  best_smooth = best_smooth,
  best_hybrid = best_hybrid,
  untethered_traj = untethered_traj,
  regularized_traj = regularized_traj,
  regularized_runtime = best_hybrid$elapsed_sec[[1L]],
  tex_path = tex_path
)

report_pdf_path <- sub("\\.tex$", ".pdf", tex_path)
cat(sprintf(
  paste(
    "Paraboloid GMDS regularization benchmark complete.",
    "Selected methods: %s",
    "Smoothness sweep: %s",
    "Hybrid sweep: %s",
    "Untethered trajectory: %s",
    "Regularized trajectory: %s",
    "Report (tex): %s",
    "Report (pdf target): %s",
    sep = "\n"
  ),
  file.path(tmp_dir, "paraboloid_regularized_selected.csv"),
  file.path(tmp_dir, "paraboloid_regularized_smooth_sweep.csv"),
  file.path(tmp_dir, "paraboloid_regularized_hybrid_sweep.csv"),
  file.path(tmp_dir, "paraboloid_regularized_untethered_trajectory.csv"),
  file.path(tmp_dir, "paraboloid_regularized_hybrid_trajectory.csv"),
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
