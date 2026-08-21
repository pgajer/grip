#!/usr/bin/env Rscript

run_tag <- "paraboloid-gmds-pathology-2026-03-31"
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
  stop("Install 'devtools' or the 'grip' package to run the paraboloid pathology benchmark.")
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
  centered <- coords
  centered <- sweep(centered, 2L, colMeans(centered), FUN = "-")
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

new_case <- function(side, iter_budget) {
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
    id = sprintf("paraboloid_%dx%d", side, side),
    label = sprintf("Paraboloid mesh %dx%d", side, side),
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
                            trace = NULL,
                            n_threads_used = NA_integer_) {
  score <- grip.score.geodesic.mds(coords, prepared = case$prepared)
  fit <- grip:::grip.align.to.target.nd(coords, case$truth, allow.reflection = TRUE)

  steps <- if (is.null(trace) || nrow(trace) == 0L) numeric(0L) else trace$step[is.finite(trace$step)]
  data.frame(
    case_id = case$id,
    case_label = case$label,
    side = case$side,
    n = case$bundle$n,
    setting_id = setting_id,
    setting_label = setting_label,
    elapsed_sec = as.double(elapsed_sec),
    n_threads_used = as.integer(n_threads_used),
    gmds_stress = score$gmds.stress[[1L]],
    gmds_raw_stress = score$gmds.raw_stress[[1L]],
    procrustes_rmse = fit$rmse,
    roughness = mesh_roughness(coords, case$adj_list, case$edges),
    area_q05_ratio = area_floor_ratio(coords, case$triangles),
    iterations = if (is.null(trace) || nrow(trace) == 0L) 0L else max(trace$iteration),
    mean_step = if (length(steps) > 0L) mean(steps) else NA_real_,
    final_step = if (length(steps) > 0L) utils::tail(steps, 1L) else NA_real_,
    gradient_norm = if (is.null(trace) || nrow(trace) == 0L) NA_real_ else utils::tail(trace$gradient_norm, 1L),
    stringsAsFactors = FALSE
  )
}

run_density_case <- function(case,
                             tether_weight = 0.05,
                             tether_weight_end = 0) {
  rows <- list()
  add_row <- function(row) {
    rows[[length(rows) + 1L]] <<- row
  }

  add_row(evaluate_layout(
    case = case,
    coords = case$truth,
    setting_id = "truth",
    setting_label = "Reference surface",
    elapsed_sec = NA_real_,
    trace = NULL,
    n_threads_used = 1L
  ))

  cmd.coords <- case$cmd$coords
  add_row(evaluate_layout(
    case = case,
    coords = cmd.coords,
    setting_id = "cmdscale",
    setting_label = "cmdscale (avg objective)",
    elapsed_sec = case$cmd_elapsed,
    trace = NULL,
    n_threads_used = 1L
  ))

  started <- proc.time()[["elapsed"]]
  avg.fit <- grip.optimize.geodesic.mds(
    coords = cmd.coords,
    prepared = case$prepared,
    max_iter = case$iter_budget,
    engine = "cpp",
    return_trace = TRUE
  )
  avg.elapsed <- proc.time()[["elapsed"]] - started
  add_row(evaluate_layout(
    case = case,
    coords = avg.fit$coords,
    setting_id = "gmds_average",
    setting_label = sprintf("GMDS avg (%d)", case$iter_budget),
    elapsed_sec = avg.elapsed,
    trace = avg.fit$trace,
    n_threads_used = avg.fit$n_threads_used
  ))

  started <- proc.time()[["elapsed"]]
  tether.fit <- grip.optimize.geodesic.mds(
    coords = cmd.coords,
    prepared = case$prepared,
    max_iter = case$iter_budget,
    engine = "cpp",
    return_trace = TRUE,
    anchor_mode = "cmdscale",
    anchor_weight = tether_weight,
    anchor_weight_end = tether_weight_end,
    continuation = "linear"
  )
  tether.elapsed <- proc.time()[["elapsed"]] - started
  add_row(evaluate_layout(
    case = case,
    coords = tether.fit$coords,
    setting_id = "gmds_linear_tether",
    setting_label = sprintf("GMDS avg + linear tether (%d)", case$iter_budget),
    elapsed_sec = tether.elapsed,
    trace = tether.fit$trace,
    n_threads_used = tether.fit$n_threads_used
  ))

  list(
    case = case,
    metrics = do.call(rbind, rows),
    layouts = list(
      truth = case$truth,
      cmdscale = cmd.coords,
      gmds_average = avg.fit$coords,
      gmds_linear_tether = tether.fit$coords
    )
  )
}

run_trajectory_case <- function(case,
                                max_iter = 25L,
                                initial_step = 1.0) {
  started <- proc.time()[["elapsed"]]
  fit <- grip.optimize.geodesic.mds(
    coords = case$cmd$coords,
    prepared = case$prepared,
    max_iter = max_iter,
    engine = "cpp",
    return_trace = TRUE,
    initial_step = initial_step
  )
  elapsed <- proc.time()[["elapsed"]] - started

  frame_metrics <- do.call(rbind, lapply(seq_along(fit$frames), function(i) {
    iteration <- i - 1L
    coords <- fit$frames[[i]]
    row <- evaluate_layout(
      case = case,
      coords = coords,
      setting_id = sprintf("trajectory_%02d", iteration),
      setting_label = sprintf("iter %d", iteration),
      elapsed_sec = if (iteration == max_iter) elapsed else NA_real_,
      trace = if (iteration == max_iter) fit$trace else NULL,
      n_threads_used = fit$n_threads_used
    )
    row$iteration <- iteration
    row$accepted_step <- if (iteration == 0L) NA_real_ else fit$trace$step[[iteration + 1L]]
    row$gradient_norm_trace <- fit$trace$gradient_norm[[iteration + 1L]]
    row
  }))

  list(
    case = case,
    fit = fit,
    metrics = frame_metrics
  )
}

run_step_sweep <- function(case,
                           initial_steps = c(0.05, 0.1, 0.25, 0.5, 1.0, 2.0),
                           max_iter = 25L) {
  rows <- list()
  layouts <- list()

  for (s0 in initial_steps) {
    started <- proc.time()[["elapsed"]]
    fit <- grip.optimize.geodesic.mds(
      coords = case$cmd$coords,
      prepared = case$prepared,
      max_iter = max_iter,
      engine = "cpp",
      return_trace = TRUE,
      initial_step = s0
    )
    elapsed <- proc.time()[["elapsed"]] - started
    row <- evaluate_layout(
      case = case,
      coords = fit$coords,
      setting_id = paste0("step_", format(s0, trim = TRUE)),
      setting_label = sprintf("s0 = %s", format(s0, trim = TRUE)),
      elapsed_sec = elapsed,
      trace = fit$trace,
      n_threads_used = fit$n_threads_used
    )
    row$initial_step <- s0
    rows[[length(rows) + 1L]] <- row
    layouts[[length(layouts) + 1L]] <- fit$coords
  }

  list(
    metrics = do.call(rbind, rows),
    layouts = layouts
  )
}

run_tether_sweep <- function(case,
                             lambdas = c(0, 0.01, 0.025, 0.05, 0.1, 0.2),
                             max_iter = 25L) {
  rows <- list()
  layouts <- list()

  for (lambda in lambdas) {
    started <- proc.time()[["elapsed"]]
    fit <- grip.optimize.geodesic.mds(
      coords = case$cmd$coords,
      prepared = case$prepared,
      max_iter = max_iter,
      engine = "cpp",
      return_trace = TRUE,
      anchor_mode = "cmdscale",
      anchor_weight = lambda,
      anchor_weight_end = lambda,
      continuation = "constant"
    )
    elapsed <- proc.time()[["elapsed"]] - started
    row <- evaluate_layout(
      case = case,
      coords = fit$coords,
      setting_id = paste0("lambda_", format(lambda, trim = TRUE)),
      setting_label = sprintf("lambda = %s", format(lambda, trim = TRUE)),
      elapsed_sec = elapsed,
      trace = fit$trace,
      n_threads_used = fit$n_threads_used
    )
    row$lambda_start <- lambda
    row$lambda_end <- lambda
    row$continuation <- "constant"
    rows[[length(rows) + 1L]] <- row
    layouts[[length(layouts) + 1L]] <- fit$coords
  }

  started <- proc.time()[["elapsed"]]
  linear.fit <- grip.optimize.geodesic.mds(
    coords = case$cmd$coords,
    prepared = case$prepared,
    max_iter = max_iter,
    engine = "cpp",
    return_trace = TRUE,
    anchor_mode = "cmdscale",
    anchor_weight = 0.05,
    anchor_weight_end = 0,
    continuation = "linear"
  )
  elapsed <- proc.time()[["elapsed"]] - started
  row <- evaluate_layout(
    case = case,
    coords = linear.fit$coords,
    setting_id = "lambda_linear_0.05_to_0",
    setting_label = "linear 0.05 -> 0",
    elapsed_sec = elapsed,
    trace = linear.fit$trace,
    n_threads_used = linear.fit$n_threads_used
  )
  row$lambda_start <- 0.05
  row$lambda_end <- 0
  row$continuation <- "linear"
  rows[[length(rows) + 1L]] <- row
  layouts[[length(layouts) + 1L]] <- linear.fit$coords

  list(
    metrics = do.call(rbind, rows),
    layouts = layouts
  )
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
    width = 1100L * ncol,
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
    mar = c(1.2, 1.2, 2.7, 0.4),
    oma = c(0, 0, 1.5, 0)
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
  invisible(NULL)
}

save_density_grid <- function(results, path) {
  grDevices::png(
    path,
    width = 4400,
    height = 950L * length(results),
    res = 180,
    bg = "#ffffff"
  )
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(
    mfrow = c(length(results), 4L),
    mar = c(1.2, 1.2, 2.7, 0.4),
    oma = c(0, 0, 1.5, 0)
  )

  for (result in results) {
    case <- result$case
    metrics <- result$metrics

    plot.layout(
      coords = case$truth,
      edges = case$edges,
      projection = "ortho",
      azimuth = 35,
      elevation = 24,
      vertex.col = "#bc6c25",
      edge.col = "#adb5bd",
      main = ""
    )
    graphics::title(main = sprintf("%s\nReference surface", case$label), cex.main = 0.74)

    for (setting_id in c("cmdscale", "gmds_average", "gmds_linear_tether")) {
      row <- metrics[metrics$setting_id == setting_id, , drop = FALSE]
      coords <- switch(
        setting_id,
        cmdscale = result$layouts$cmdscale,
        gmds_average = result$layouts$gmds_average,
        gmds_linear_tether = result$layouts$gmds_linear_tether
      )
      coords <- grip:::grip.align.to.target.nd(coords, case$truth, allow.reflection = TRUE)$aligned
      plot.layout(
        coords = coords,
        edges = case$edges,
        projection = "ortho",
        azimuth = 35,
        elevation = 24,
        vertex.col = "#3a5a40",
        edge.col = "#adb5bd",
        main = ""
      )
      label <- switch(
        setting_id,
        cmdscale = "cmdscale",
        gmds_average = "GMDS avg",
        gmds_linear_tether = "GMDS linear tether"
      )
      graphics::title(
        main = sprintf(
          "%s\nsigma %s, rho %s\neta %s, alpha_0.05 %s\nt %ss",
          label,
          fmt_num(row$gmds_stress[[1L]]),
          fmt_num(row$procrustes_rmse[[1L]]),
          fmt_num(row$roughness[[1L]]),
          fmt_num(row$area_q05_ratio[[1L]]),
          fmt_time(row$elapsed_sec[[1L]])
        ),
        cex.main = 0.72
      )
    }
  }

  graphics::mtext(
    "Paraboloid density context: reference, MDS, GMDS, and linear tether",
    outer = TRUE,
    cex = 1.15,
    font = 2
  )
}

save_trajectory_trace_plot <- function(traj_df, path) {
  grDevices::png(path, width = 2200, height = 1400, res = 180, bg = "#ffffff")
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(2L, 3L), mar = c(4, 4, 2.5, 1))

  panels <- list(
    list(y = traj_df$gmds_stress, title = expression(sigma ~ "over iterations"), ylab = expression(sigma)),
    list(y = traj_df$procrustes_rmse, title = expression(rho ~ "over iterations"), ylab = expression(rho)),
    list(y = traj_df$roughness, title = expression(eta ~ "over iterations"), ylab = expression(eta)),
    list(y = traj_df$area_q05_ratio, title = expression(alpha[0.05] ~ "over iterations"), ylab = expression(alpha[0.05])),
    list(y = traj_df$accepted_step, title = "Accepted step over iterations", ylab = "step"),
    list(y = traj_df$gradient_norm_trace, title = "Gradient norm over iterations", ylab = "||grad||")
  )

  for (panel in panels) {
    graphics::plot(
      traj_df$iteration,
      panel$y,
      type = "b",
      pch = 16,
      col = "#1d3557",
      xlab = "GMDS iteration",
      ylab = panel$ylab,
      main = panel$title
    )
  }
}

write_density_table <- function(df) {
  rows <- apply(df, 1L, function(row) {
    paste(
      row[["side"]],
      fmt_time(as.numeric(row[["elapsed_sec"]])),
      row[["n_threads_used"]],
      tex_escape(row[["setting_label"]]),
      fmt_num(as.numeric(row[["gmds_stress"]])),
      fmt_num(as.numeric(row[["procrustes_rmse"]])),
      fmt_num(as.numeric(row[["roughness"]])),
      fmt_num(as.numeric(row[["area_q05_ratio"]])),
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

write_step_table <- function(df) {
  rows <- apply(df, 1L, function(row) {
    paste(
      tex_escape(row[["setting_label"]]),
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

write_tether_table <- function(df) {
  rows <- apply(df, 1L, function(row) {
    paste(
      tex_escape(row[["setting_label"]]),
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

write_report <- function(density_df,
                         trajectory_df,
                         step_df,
                         tether_df,
                         tex_path) {
  selected_traj <- trajectory_df[
    trajectory_df$iteration %in% c(0L, 1L, 2L, 5L, 10L, 15L, 20L, 25L),
    ,
    drop = FALSE
  ]
  best_step <- step_df[which.min(step_df$roughness), , drop = FALSE]
  best_tether <- tether_df[which.min(tether_df$roughness), , drop = FALSE]
  lowest_stress_tether <- tether_df[which.min(tether_df$gmds_stress), , drop = FALSE]

  tex <- c(
    "\\documentclass[11pt]{article}",
    "\\usepackage[margin=1in]{geometry}",
    "\\usepackage{graphicx}",
    "\\usepackage{booktabs}",
    "\\usepackage{float}",
    "\\usepackage{amsmath}",
    "\\usepackage{pdflscape}",
    "\\usepackage{hyperref}",
    "\\title{Paraboloid GMDS Pathology Report}",
    "\\author{Focused paraboloid diagnostic study}",
    "\\date{March 31, 2026}",
    "",
    "\\begin{document}",
    "\\maketitle",
    "",
    "\\begin{abstract}",
    "This standalone report isolates the visually unsettling singularities that appear in the 3D paraboloid GMDS layouts. The goal is not only to show the pathology, but to determine whether it is driven mainly by the gradient-descent step schedule or by the stress objective itself, and to test whether stronger tethering to the classical-MDS shape can suppress it.",
    "\\end{abstract}",
    "",
    "\\section{Questions and Metrics}",
    "We study orthogonal weighted paraboloid meshes in 3D. The baseline observation is that classical MDS is smooth but globally distorted, whereas untethered GMDS lowers geodesic stress much more aggressively but can concentrate curvature into singular-looking regions.",
    "To make that pathology measurable we track four quantities: the GMDS stress $\\sigma$, Procrustes RMSE $\\rho$ against the reference paraboloid surface, a normalized roughness proxy",
    "\\[",
    "\\eta = \\frac{\\sqrt{\\frac{1}{n}\\sum_i \\lVert z_i - \\frac{1}{|N(i)|}\\sum_{j \\in N(i)} z_j \\rVert^2}}{\\operatorname{median}_{(u,v) \\in E} \\lVert z_u - z_v \\rVert},",
    "\\]",
    "and a cell-degeneracy proxy",
    "\\[",
    "\\alpha_{0.05} = \\frac{Q_{0.05}(\\text{triangle areas})}{\\operatorname{median}(\\text{triangle areas})}.",
    "\\]",
    "Lower $\\eta$ means smoother local geometry, and larger $\\alpha_{0.05}$ means fewer near-collapsed mesh cells.",
    "",
    "\\section{Density Context}",
    "Before the focused experiments, we reproduce the existing paraboloid density context for $12\\times12$, $15\\times15$, and $20\\times20$ meshes. The same pattern from the larger GMDS/MDS report remains visible here: untethered GMDS cuts $\\sigma$ dramatically, but can do so while increasing roughness and reducing the lower tail of local cell areas.",
    "\\begin{figure}[H]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "paraboloid_density_context.png"))
    ),
    "\\caption{Paraboloid density context. Each row shows the reference surface, classical MDS, untethered tie-averaged GMDS, and the current linear-tether variant.}",
    "\\end{figure}",
    "\\begin{table}[H]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\caption{Paraboloid density context. Lower $\\sigma$, $\\rho$, and $\\eta$ are better; larger $\\alpha_{0.05}$ is better.}",
    "\\begin{tabular}{@{}rrrlrrrr@{}}",
    "\\toprule",
    "Side & $t$ (s) & Threads & Setting & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ \\\\",
    "\\midrule",
    write_density_table(density_df),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    "\\section{Experiment 1: 12x12 MDS-to-GMDS Trajectory}",
    "The first question is whether the singularity emerges gradually as GMDS improves stress, or whether it is tied to a visibly unstable numerical jump. To answer that, we start from the 3D classical-MDS embedding and record every accepted GMDS frame for 25 extra iterations.",
    sprintf(
      "In this run, all 25 GMDS iterations were accepted. The accepted step size collapses quickly to a small value, yet $\\sigma$ continues to decrease from %s at iteration 0 to %s at iteration 25. Over the same interval, $\\rho$ moves from %s to %s, roughness $\\eta$ moves from %s to %s, and $\\alpha_{0.05}$ moves from %s to %s.",
      fmt_num(trajectory_df$gmds_stress[trajectory_df$iteration == 0L]),
      fmt_num(trajectory_df$gmds_stress[trajectory_df$iteration == 25L]),
      fmt_num(trajectory_df$procrustes_rmse[trajectory_df$iteration == 0L]),
      fmt_num(trajectory_df$procrustes_rmse[trajectory_df$iteration == 25L]),
      fmt_num(trajectory_df$roughness[trajectory_df$iteration == 0L]),
      fmt_num(trajectory_df$roughness[trajectory_df$iteration == 25L]),
      fmt_num(trajectory_df$area_q05_ratio[trajectory_df$iteration == 0L]),
      fmt_num(trajectory_df$area_q05_ratio[trajectory_df$iteration == 25L])
    ),
    "That pattern is important: if the accepted steps are already tiny while the surface continues to sharpen, the pathology is much more likely to be objective-driven than to be a simple overshoot artifact.",
    "\\begin{figure}[H]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=0.95\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "paraboloid_trajectory_metrics.png"))
    ),
    "\\caption{Metric traces for the 25-step GMDS correction of the $12\\times12$ paraboloid.}",
    "\\end{figure}",
    "\\begin{landscape}",
    "\\begin{figure}[p]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "paraboloid_trajectory_part1.png"))
    ),
    "\\caption{Trajectory snapshots, iterations 0 through 12. Iteration 0 is the classical-MDS initialization.}",
    "\\end{figure}",
    "\\end{landscape}",
    "\\begin{landscape}",
    "\\begin{figure}[p]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "paraboloid_trajectory_part2.png"))
    ),
    "\\caption{Trajectory snapshots, iterations 13 through 25.}",
    "\\end{figure}",
    "\\end{landscape}",
    "\\begin{table}[H]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{6pt}",
    "\\caption{Selected trajectory checkpoints on the $12\\times12$ paraboloid.}",
    "\\begin{tabular}{@{}rrrrrr@{}}",
    "\\toprule",
    "Iter & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ & step \\\\",
    "\\midrule",
    write_trajectory_table(selected_traj),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    "\\section{Experiment 2: Initial-Step Sweep}",
    "The second question is whether the singularity is primarily a gradient-descent tuning issue. Because the optimizer uses Armijo backtracking, the initial step $s_0$ only sets the trial scale; the accepted steps can still become much smaller. We therefore sweep $s_0 \\in \\{0.05, 0.1, 0.25, 0.5, 1, 2\\}$ on the same $12\\times12$ paraboloid for 25 GMDS iterations.",
    sprintf(
      "The smoothest final result in this sweep is %s, while the roughest is %s. The best roughness in the sweep is achieved by %s with $\\eta=%s$, $\\sigma=%s$, and $\\alpha_{0.05}=%s$.",
      tex_escape(best_step$setting_label[[1L]]),
      tex_escape(step_df$setting_label[[which.max(step_df$roughness)]]),
      tex_escape(best_step$setting_label[[1L]]),
      fmt_num(best_step$roughness[[1L]]),
      fmt_num(best_step$gmds_stress[[1L]]),
      fmt_num(best_step$area_q05_ratio[[1L]])
    ),
    "If the pathology were mostly due to overshooting, we would expect the smallest initial steps to look dramatically better. If instead the final surfaces remain qualitatively similar across the sweep, that points back to the GMDS objective rather than the outer step-size knob.",
    "\\begin{figure}[H]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "paraboloid_step_sweep.png"))
    ),
    "\\caption{Initial-step sweep on the $12\\times12$ paraboloid. The first two panels are the reference surface and the MDS initialization; the remaining panels show 25-step untethered GMDS runs for different initial step values.}",
    "\\end{figure}",
    "\\begin{table}[H]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\caption{Initial-step sweep on the $12\\times12$ paraboloid.}",
    "\\begin{tabular}{@{}lrrrrrr@{}}",
    "\\toprule",
    "Setting & $t$ (s) & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ & mean step \\\\",
    "\\midrule",
    write_step_table(step_df),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    "\\section{Experiment 3: Constant-Tether Sweep}",
    "The third question is whether a stronger shape prior from the classical-MDS initialization can suppress the singularity without destroying the GMDS correction entirely. We therefore hold the tether constant for all 25 GMDS iterations and test $\\lambda \\in \\{0, 0.01, 0.025, 0.05, 0.1, 0.2\\}$, plus the current linear schedule $0.05 \\to 0$ as a reference.",
    sprintf(
      "Among these runs, the lowest stress is still achieved by %s with $\\sigma=%s$, but the smoothest result is %s with $\\eta=%s$ and $\\alpha_{0.05}=%s$.",
      tex_escape(lowest_stress_tether$setting_label[[1L]]),
      fmt_num(lowest_stress_tether$gmds_stress[[1L]]),
      tex_escape(best_tether$setting_label[[1L]]),
      fmt_num(best_tether$roughness[[1L]]),
      fmt_num(best_tether$area_q05_ratio[[1L]])
    ),
    "This experiment is the most direct way to see the stress-versus-regularity tradeoff: larger $\\lambda$ keeps the surface closer to the smooth MDS shape, while smaller $\\lambda$ gives the GMDS objective more freedom to concentrate distortion.",
    "\\begin{figure}[H]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "paraboloid_tether_sweep.png"))
    ),
    "\\caption{Constant-tether sweep on the $12\\times12$ paraboloid. The first two panels are the reference surface and the MDS initialization.}",
    "\\end{figure}",
    "\\begin{table}[H]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\caption{Constant-tether sweep on the $12\\times12$ paraboloid.}",
    "\\begin{tabular}{@{}lrrrrrr@{}}",
    "\\toprule",
    "Setting & $t$ (s) & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ & mean step \\\\",
    "\\midrule",
    write_tether_table(tether_df),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    "\\section{Takeaways}",
    "The MDS-to-GMDS movie is the key diagnostic. If stress keeps falling while accepted steps are already tiny and the surface still sharpens into singular-looking regions, the pathology is not just a poorly tuned outer step length. The step-sweep section then tells us whether smaller or larger trial steps materially change the endpoint; if not, the dominant issue is the objective landscape. The constant-tether sweep is the most promising immediate mitigation because it directly trades stress reduction against geometric regularity.",
    "A practical next default for paraboloid-like families is therefore not to abandon GMDS, but to treat tether strength as a first-class regularization parameter and to tune it jointly with the iteration budget. If we later want a more intrinsic fix, the natural next step is to augment GMDS with an explicit local smoothness or bending penalty rather than relying only on the MDS anchor.",
    "",
    "\\end{document}"
  )

  writeLines(tex, tex_path)
}

density_cases <- list(
  new_case(12L, 24L),
  new_case(15L, 20L),
  new_case(20L, 16L)
)
density_results <- lapply(density_cases, run_density_case)
density_df <- do.call(rbind, lapply(density_results, `[[`, "metrics"))
density_df <- density_df[density_df$setting_id != "truth", , drop = FALSE]
density_path <- file.path(tmp_dir, "paraboloid_density_context.csv")
utils::write.csv(density_df, density_path, row.names = FALSE)

focus_case <- density_cases[[1L]]
trajectory <- run_trajectory_case(focus_case, max_iter = 25L, initial_step = 1.0)
trajectory_df <- trajectory$metrics
trajectory_path <- file.path(tmp_dir, "paraboloid_trajectory_metrics.csv")
utils::write.csv(trajectory_df, trajectory_path, row.names = FALSE)

step_sweep <- run_step_sweep(focus_case, max_iter = 25L)
step_df <- step_sweep$metrics
step_path <- file.path(tmp_dir, "paraboloid_step_sweep.csv")
utils::write.csv(step_df, step_path, row.names = FALSE)

tether_sweep <- run_tether_sweep(focus_case, max_iter = 25L)
tether_df <- tether_sweep$metrics
tether_path <- file.path(tmp_dir, "paraboloid_tether_sweep.csv")
utils::write.csv(tether_df, tether_path, row.names = FALSE)

save_density_grid(
  results = density_results,
  path = file.path(pdf_dir, "paraboloid_density_context.png")
)

save_trajectory_trace_plot(
  traj_df = trajectory_df,
  path = file.path(pdf_dir, "paraboloid_trajectory_metrics.png")
)

traj_titles <- vapply(seq_len(nrow(trajectory_df)), function(i) {
  sprintf(
    "iter %d\nsigma %s, rho %s\neta %s, alpha_0.05 %s",
    trajectory_df$iteration[[i]],
    fmt_num(trajectory_df$gmds_stress[[i]]),
    fmt_num(trajectory_df$procrustes_rmse[[i]]),
    fmt_num(trajectory_df$roughness[[i]]),
    fmt_num(trajectory_df$area_q05_ratio[[i]])
  )
}, character(1L))
traj_coords <- c(list(focus_case$truth), trajectory$fit$frames)
traj_titles_full <- c("Reference surface", traj_titles)

save_surface_panel_grid(
  case = focus_case,
  coords_list = traj_coords[1:14],
  title_list = traj_titles_full[1:14],
  path = file.path(pdf_dir, "paraboloid_trajectory_part1.png"),
  nrow = 4L,
  ncol = 4L,
  heading = "Paraboloid 12x12 trajectory: reference plus iterations 0 through 12",
  truth_index = 1L
)
save_surface_panel_grid(
  case = focus_case,
  coords_list = c(list(focus_case$truth), trajectory$fit$frames[14:26]),
  title_list = c("Reference surface", traj_titles[14:26]),
  path = file.path(pdf_dir, "paraboloid_trajectory_part2.png"),
  nrow = 4L,
  ncol = 4L,
  heading = "Paraboloid 12x12 trajectory: reference plus iterations 13 through 25",
  truth_index = 1L
)

step_titles <- c(
  "Reference surface",
  sprintf(
    "cmdscale\nsigma %s, rho %s\neta %s, alpha_0.05 %s",
    fmt_num(density_results[[1L]]$metrics$gmds_stress[density_results[[1L]]$metrics$setting_id == "cmdscale"]),
    fmt_num(density_results[[1L]]$metrics$procrustes_rmse[density_results[[1L]]$metrics$setting_id == "cmdscale"]),
    fmt_num(density_results[[1L]]$metrics$roughness[density_results[[1L]]$metrics$setting_id == "cmdscale"]),
    fmt_num(density_results[[1L]]$metrics$area_q05_ratio[density_results[[1L]]$metrics$setting_id == "cmdscale"])
  ),
  vapply(seq_len(nrow(step_df)), function(i) {
    sprintf(
      "%s\nsigma %s, rho %s\neta %s, alpha_0.05 %s",
      step_df$setting_label[[i]],
      fmt_num(step_df$gmds_stress[[i]]),
      fmt_num(step_df$procrustes_rmse[[i]]),
      fmt_num(step_df$roughness[[i]]),
      fmt_num(step_df$area_q05_ratio[[i]])
    )
  }, character(1L))
)
save_surface_panel_grid(
  case = focus_case,
  coords_list = c(list(focus_case$truth, focus_case$cmd$coords), step_sweep$layouts),
  title_list = step_titles,
  path = file.path(pdf_dir, "paraboloid_step_sweep.png"),
  nrow = 2L,
  ncol = 4L,
  heading = "Paraboloid 12x12 initial-step sweep",
  truth_index = 1L
)

tether_titles <- c(
  "Reference surface",
  sprintf(
    "cmdscale\nsigma %s, rho %s\neta %s, alpha_0.05 %s",
    fmt_num(density_results[[1L]]$metrics$gmds_stress[density_results[[1L]]$metrics$setting_id == "cmdscale"]),
    fmt_num(density_results[[1L]]$metrics$procrustes_rmse[density_results[[1L]]$metrics$setting_id == "cmdscale"]),
    fmt_num(density_results[[1L]]$metrics$roughness[density_results[[1L]]$metrics$setting_id == "cmdscale"]),
    fmt_num(density_results[[1L]]$metrics$area_q05_ratio[density_results[[1L]]$metrics$setting_id == "cmdscale"])
  ),
  vapply(seq_len(nrow(tether_df)), function(i) {
    sprintf(
      "%s\nsigma %s, rho %s\neta %s, alpha_0.05 %s",
      tether_df$setting_label[[i]],
      fmt_num(tether_df$gmds_stress[[i]]),
      fmt_num(tether_df$procrustes_rmse[[i]]),
      fmt_num(tether_df$roughness[[i]]),
      fmt_num(tether_df$area_q05_ratio[[i]])
    )
  }, character(1L))
)
save_surface_panel_grid(
  case = focus_case,
  coords_list = c(list(focus_case$truth, focus_case$cmd$coords), tether_sweep$layouts),
  title_list = tether_titles,
  path = file.path(pdf_dir, "paraboloid_tether_sweep.png"),
  nrow = 3L,
  ncol = 3L,
  heading = "Paraboloid 12x12 constant-tether sweep",
  truth_index = 1L
)

tex_path <- file.path(manual_root, "reports", "paraboloid_gmds_pathology_report_2026-03-31.tex")
write_report(
  density_df = density_df,
  trajectory_df = trajectory_df,
  step_df = step_df,
  tether_df = tether_df,
  tex_path = tex_path
)

report_pdf_path <- sub("\\.tex$", ".pdf", tex_path)
cat(sprintf(
  paste(
    "Paraboloid GMDS pathology benchmark complete.",
    "Density context: %s",
    "Trajectory metrics: %s",
    "Step sweep: %s",
    "Tether sweep: %s",
    "Report (tex): %s",
    "Report (pdf target): %s",
    sep = "\n"
  ),
  density_path,
  trajectory_path,
  step_path,
  tether_path,
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
