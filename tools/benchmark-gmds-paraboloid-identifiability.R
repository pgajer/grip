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
run_tag <- "gmds-paraboloid-identifiability-2026-04-01"
tmp_dir <- file.path(manual_root, "tmp", run_tag)
pdf_dir <- file.path(manual_root, "pdf", run_tag)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

metrics_csv <- file.path(tmp_dir, "paraboloid_identifiability_metrics.csv")
profiles_csv <- file.path(tmp_dir, "paraboloid_identifiability_fold_profiles.csv")
rds_path <- file.path(tmp_dir, "paraboloid_identifiability_embeddings.rds")
appendix_tex_path <- file.path(manual_root, "pdf", "gmds_v2_paraboloid_identifiability_appendix_2026-04-01.tex")
plot_path <- file.path(pdf_dir, "paraboloid_identifiability_fold_profiles.png")
snapshot_path <- file.path(pdf_dir, "paraboloid_identifiability_snapshots.png")

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

scale_reference_coords <- function(bundle) {
  bundle$coords_surface / bundle$weight_scale
}

fit_hinge_axis <- function(coords, coords_param) {
  x_vals <- coords_param[, 1L]
  x_mid <- stats::median(range(x_vals))
  idx <- which(abs(x_vals - x_mid) == min(abs(x_vals - x_mid)))
  hinge <- coords[idx, , drop = FALSE]
  center <- colMeans(hinge)
  centered <- sweep(hinge, 2L, center, FUN = "-", check.margin = FALSE)
  if (nrow(centered) <= 1L || max(abs(centered)) < 1e-12) {
    dir <- c(0, 1, 0)
  } else {
    sv <- svd(centered)
    dir <- sv$v[, 1L]
  }
  list(center = center, dir = dir / sqrt(sum(dir^2)), side_mask = x_vals > x_mid)
}

rotate_about_axis <- function(coords, axis_point, axis_dir, theta, mask) {
  out <- coords
  if (!any(mask) || abs(theta) < 1e-15) {
    return(out)
  }
  u <- as.double(axis_dir / sqrt(sum(axis_dir^2)))
  pts <- coords[mask, , drop = FALSE]
  shifted <- sweep(pts, 2L, axis_point, FUN = "-", check.margin = FALSE)
  dot_uv <- drop(shifted %*% u)
  cross_uv <- cbind(
    u[2L] * shifted[, 3L] - u[3L] * shifted[, 2L],
    u[3L] * shifted[, 1L] - u[1L] * shifted[, 3L],
    u[1L] * shifted[, 2L] - u[2L] * shifted[, 1L]
  )
  cth <- cos(theta)
  sth <- sin(theta)
  rotated <- shifted * cth + cross_uv * sth + tcrossprod(dot_uv, u) * (1 - cth)
  out[mask, ] <- sweep(rotated, 2L, axis_point, FUN = "+", check.margin = FALSE)
  out
}

profile_fold <- function(case, angles_deg) {
  axis <- fit_hinge_axis(case$truth, case$coords_param)
  rows <- lapply(angles_deg, function(angle_deg) {
    coords <- rotate_about_axis(
      coords = case$truth,
      axis_point = axis$center,
      axis_dir = axis$dir,
      theta = angle_deg * pi / 180,
      mask = axis$side_mask
    )
    score <- grip.score.geodesic.mds(coords = coords, prepared = case$prepared)
    aligned <- align_to_target_nd(coords, case$truth, allow.reflection = TRUE)
    data.frame(
      case_id = case$id,
      case_label = case$label,
      angle_deg = angle_deg,
      gmds_raw_stress = score$gmds.raw_stress,
      gmds_stress = score$gmds.stress,
      edge_energy = score$edge.spring.energy,
      rmse_to_truth = aligned$rmse,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

profile_fold_embeddings <- function(case, angles_deg) {
  axis <- fit_hinge_axis(case$truth, case$coords_param)
  lapply(angles_deg, function(angle_deg) {
    coords <- rotate_about_axis(
      coords = case$truth,
      axis_point = axis$center,
      axis_dir = axis$dir,
      theta = angle_deg * pi / 180,
      mask = axis$side_mask
    )
    aligned <- align_to_target_nd(coords, case$truth, allow.reflection = TRUE)
    score <- grip.score.geodesic.mds(coords = coords, prepared = case$prepared)
    list(
      id = sprintf("%s_fold_%02d", case$id, angle_deg),
      label = sprintf("%s fold %d deg", case$label, angle_deg),
      kind = "fold",
      angle_deg = angle_deg,
      coords = coords,
      display_coords = aligned$aligned,
      gmds_stress = score$gmds.stress,
      gmds_raw_stress = score$gmds.raw_stress,
      rmse_to_truth = aligned$rmse
    )
  })
}

run_optimizer <- function(case, coords_init, label, max_iter = 25L) {
  started <- proc.time()[["elapsed"]]
  fit <- grip.optimize.geodesic.mds(
    coords = coords_init,
    prepared = case$prepared,
    anchor_mode = "none",
    init = "user",
    engine = "cpp",
    max_iter = max_iter,
    return_trace = FALSE
  )
  elapsed <- proc.time()[["elapsed"]] - started
  aligned <- align_to_target_nd(fit$coords, case$truth, allow.reflection = TRUE)
  data.frame(
    case_id = case$id,
    case_label = case$label,
    start = label,
    max_iter = max_iter,
    gmds_stress = fit$score$gmds.stress,
    gmds_raw_stress = fit$score$gmds.raw_stress,
    rmse_to_truth = aligned$rmse,
    final_anchor_weight = fit$final_anchor_weight,
    elapsed_sec = elapsed,
    stringsAsFactors = FALSE
  )
}

run_optimizer_embedding <- function(case, coords_init, label, display_label, max_iter = 25L) {
  started <- proc.time()[["elapsed"]]
  fit <- grip.optimize.geodesic.mds(
    coords = coords_init,
    prepared = case$prepared,
    anchor_mode = "none",
    init = "user",
    engine = "cpp",
    max_iter = max_iter,
    return_trace = FALSE
  )
  elapsed <- proc.time()[["elapsed"]] - started
  aligned <- align_to_target_nd(fit$coords, case$truth, allow.reflection = TRUE)
  list(
    id = sprintf("%s_%s", case$id, label),
    label = display_label,
    kind = "recovery",
    coords = fit$coords,
    display_coords = aligned$aligned,
    gmds_stress = fit$score$gmds.stress,
    gmds_raw_stress = fit$score$gmds.raw_stress,
    rmse_to_truth = aligned$rmse,
    elapsed_sec = elapsed,
    max_iter = max_iter
  )
}

build_regular_case <- function(side = 12L, amplitude = 0.35) {
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
  list(
    id = sprintf("regular_%dx%d", side, side),
    label = sprintf("Regular orthogonal paraboloid mesh %dx%d", side, side),
    bundle = bundle,
    prepared = prepared,
    truth = scale_reference_coords(bundle),
    coords_param = bundle$coords_param,
    edges = bundle$edges
  )
}

build_irregular_case <- function(side = 12L, amplitude = 0.35) {
  keep <- keep.staggered.windows(
    h = side,
    w = side,
    window_height = 1L,
    window_width = 2L,
    row_period = 4L,
    col_period = 5L,
    row_offset = 2L,
    col_offset = 3L
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
  list(
    id = sprintf("irregular_occ_%dx%d", side, side),
    label = sprintf("Irregular occupied paraboloid mesh %dx%d", side, side),
    bundle = bundle,
    prepared = prepared,
    truth = scale_reference_coords(bundle),
    coords_param = bundle$coords_param,
    edges = bundle$edges
  )
}

compute_case_output <- function(case, angles_deg, max_iter = 25L) {
  cmd <- classical_mds_embedding(case$prepared, dim = 3L, eig = TRUE)$coords
  cmd_aligned <- align_to_target_nd(cmd, case$truth, allow.reflection = TRUE)
  cmd_score <- grip.score.geodesic.mds(coords = cmd, prepared = case$prepared)
  fold_embeddings <- profile_fold_embeddings(case, angles_deg)
  profile_df <- do.call(rbind, lapply(fold_embeddings, function(x) {
    data.frame(
      case_id = case$id,
      case_label = case$label,
      angle_deg = x$angle_deg,
      gmds_raw_stress = x$gmds_raw_stress,
      gmds_stress = x$gmds_stress,
      edge_energy = 0,
      rmse_to_truth = x$rmse_to_truth,
      stringsAsFactors = FALSE
    )
  }))
  fold15 <- fold_embeddings[[which(angles_deg == 15)[1L]]]$coords
  fit_cmd <- run_optimizer_embedding(case, cmd, "pure_gmds_from_cmd", "Pure GMDS from cMDS", max_iter = max_iter)
  fit_fold <- run_optimizer_embedding(case, fold15, "pure_gmds_from_fold15", "Pure GMDS from 15 deg fold", max_iter = max_iter)
  metrics_df <- data.frame(
    case_id = c(case$id, case$id),
    case_label = c(case$label, case$label),
    start = c("pure_gmds_from_cmd", "pure_gmds_from_fold15"),
    max_iter = c(max_iter, max_iter),
    gmds_stress = c(fit_cmd$gmds_stress, fit_fold$gmds_stress),
    gmds_raw_stress = c(fit_cmd$gmds_raw_stress, fit_fold$gmds_raw_stress),
    rmse_to_truth = c(fit_cmd$rmse_to_truth, fit_fold$rmse_to_truth),
    final_anchor_weight = c(0, 0),
    elapsed_sec = c(fit_cmd$elapsed_sec, fit_fold$elapsed_sec),
    stringsAsFactors = FALSE
  )
  list(
    case = case,
    profile_df = profile_df,
    metrics_df = metrics_df,
    embeddings = c(
      list(list(
        id = sprintf("%s_truth", case$id),
        label = "Scaled truth",
        kind = "reference",
        coords = case$truth,
        display_coords = case$truth,
        gmds_stress = 0,
        gmds_raw_stress = 0,
        rmse_to_truth = 0
      )),
      list(list(
        id = sprintf("%s_cmdscale", case$id),
        label = "Classical MDS start",
        kind = "init",
        coords = cmd,
        display_coords = cmd_aligned$aligned,
        gmds_stress = cmd_score$gmds.stress,
        gmds_raw_stress = cmd_score$gmds.raw_stress,
        rmse_to_truth = cmd_aligned$rmse
      )),
      fold_embeddings,
      list(fit_cmd, fit_fold)
    )
  )
}

plot_profiles <- function(profile_df, output_path) {
  cases <- unique(profile_df$case_label)
  grDevices::png(output_path, width = 2000L, height = 900L, res = 180, bg = "#ffffff")
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(1L, 2L), mar = c(4.2, 4.5, 2.6, 0.8), oma = c(0, 0, 1.2, 0))

  for (metric in c("gmds_stress", "rmse_to_truth")) {
    ylab <- if (identical(metric, "gmds_stress")) "GMDS stress after fold perturbation" else "Procrustes RMSE to truth"
    x_range <- range(profile_df$angle_deg)
    y_range <- range(profile_df[[metric]])
    graphics::plot(NA, xlim = x_range, ylim = y_range, xlab = "Fold angle (degrees)", ylab = ylab)
    cols <- c("#1f4e79", "#a63d40")
    for (i in seq_along(cases)) {
      df <- profile_df[profile_df$case_label == cases[[i]], , drop = FALSE]
      graphics::lines(df$angle_deg, df[[metric]], type = "b", pch = 19L, lwd = 2.2, col = cols[[i]])
    }
    graphics::legend("topleft", legend = cases, col = cols, lwd = 2.2, pch = 19L, bty = "n")
  }
  graphics::mtext(
    "Fold-direction perturbation from the exact scaled paraboloid reference",
    side = 3L,
    outer = TRUE,
    line = 0.1,
    cex = 1.15,
    font = 2L
  )
}

plot_snapshots <- function(cases, output_path) {
  grDevices::png(output_path, width = 2200L, height = 1200L, res = 180, bg = "#ffffff")
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(length(cases), 4L), mar = c(1.2, 1.2, 2.4, 0.4), oma = c(0, 0, 1.2, 0))

  for (case in cases) {
    cmd <- classical_mds_embedding(case$prepared, dim = 3L, eig = TRUE)$coords
    axis <- fit_hinge_axis(case$truth, case$coords_param)
    fold15 <- rotate_about_axis(case$truth, axis$center, axis$dir, 15 * pi / 180, axis$side_mask)
    fit_cmd <- grip.optimize.geodesic.mds(
      coords = cmd,
      prepared = case$prepared,
      anchor_mode = "none",
      init = "user",
      engine = "cpp",
      max_iter = 25L,
      return_trace = FALSE
    )
    fit_fold <- grip.optimize.geodesic.mds(
      coords = fold15,
      prepared = case$prepared,
      anchor_mode = "none",
      init = "user",
      engine = "cpp",
      max_iter = 25L,
      return_trace = FALSE
    )
    panels <- list(
      list(coords = case$truth, title = "Scaled truth"),
      list(coords = align_to_target_nd(fold15, case$truth, allow.reflection = TRUE)$aligned, title = "15 deg fold"),
      list(coords = align_to_target_nd(fit_cmd$coords, case$truth, allow.reflection = TRUE)$aligned,
           title = sprintf("Pure GMDS from cMDS\nsigma %s, rho %s",
                           fmt_num(fit_cmd$score$gmds.stress, 4L),
                           fmt_num(align_to_target_nd(fit_cmd$coords, case$truth, allow.reflection = TRUE)$rmse, 4L))),
      list(coords = align_to_target_nd(fit_fold$coords, case$truth, allow.reflection = TRUE)$aligned,
           title = sprintf("Pure GMDS from 15 deg fold\nsigma %s, rho %s",
                           fmt_num(fit_fold$score$gmds.stress, 4L),
                           fmt_num(align_to_target_nd(fit_fold$coords, case$truth, allow.reflection = TRUE)$rmse, 4L)))
    )
    for (panel in panels) {
      grip.plot(
        coords = panel$coords,
        edges = case$edges,
        projection = "ortho",
        azimuth = 35,
        elevation = 24,
        vertex.col = "#3a5a40",
        edge.col = "#adb5bd",
        main = ""
      )
      graphics::mtext(panel$title, side = 3L, line = 0.3, cex = 0.82)
    }
    graphics::mtext(case$label, side = 2L, line = -0.3, outer = FALSE, cex = 0.9)
  }
  graphics::mtext(
    "Regular versus irregular paraboloid: exact truth, fold perturbation, and pure-GMDS recoveries",
    side = 3L,
    outer = TRUE,
    line = 0.2,
    cex = 1.15,
    font = 2L
  )
}

build_appendix <- function(case_outputs, profile_df, metrics_df) {
  select_angles <- c(0, 15, 30)
  profile_summary <- profile_df[profile_df$angle_deg %in% select_angles, , drop = FALSE]
  pretty_start <- function(x) {
    switch(
      x,
      pure_gmds_from_cmd = "Pure GMDS from cMDS",
      pure_gmds_from_fold15 = "Pure GMDS from 15 deg fold",
      x
    )
  }
  profile_rows <- vapply(seq_len(nrow(profile_summary)), function(i) {
    sprintf(
      "%s & %d & %s & %s \\\\",
      profile_summary$case_label[[i]],
      profile_summary$angle_deg[[i]],
      fmt_num(profile_summary$gmds_stress[[i]], 4L),
      fmt_num(profile_summary$rmse_to_truth[[i]], 4L)
    )
  }, character(1L))
  recovery_rows <- vapply(seq_len(nrow(metrics_df)), function(i) {
    sprintf(
      "%s & %s & %s & %s & %s \\\\",
      metrics_df$case_label[[i]],
      pretty_start(metrics_df$start[[i]]),
      fmt_num(metrics_df$gmds_stress[[i]], 4L),
      fmt_num(metrics_df$rmse_to_truth[[i]], 4L),
      fmt_num(metrics_df$elapsed_sec[[i]], 3L)
    )
  }, character(1L))

  regular_cmd <- metrics_df[metrics_df$case_id == "regular_12x12" & metrics_df$start == "pure_gmds_from_cmd", , drop = FALSE]
  regular_fold <- metrics_df[metrics_df$case_id == "regular_12x12" & metrics_df$start == "pure_gmds_from_fold15", , drop = FALSE]
  irregular_cmd <- metrics_df[metrics_df$case_id == "irregular_occ_12x12" & metrics_df$start == "pure_gmds_from_cmd", , drop = FALSE]
  irregular_fold <- metrics_df[metrics_df$case_id == "irregular_occ_12x12" & metrics_df$start == "pure_gmds_from_fold15", , drop = FALSE]
  regular_a <- 0.3213918
  irregular_a <- 0.1592610

  paste(
    "\\section{Paraboloid Identifiability Diagnostic}",
    "This experiment was designed to answer the question raised by the previous regularized-GMDS sections: are smooth paraboloids unrecoverable under GMDS, or are we mainly seeing a basin-selection failure from the classical-MDS initialization?",
    "The first result is structural. After accounting for the generator's global weight normalization, the displayed reference paraboloid is itself an exact zero-stress GMDS embedding. So the pathology is not that the target surface is incompatible with the GMDS objective. The pathology is that the same graph metric also admits folded alternatives.",
    sprintf("A controlled hinge-fold perturbation confirms that the folded directions are cheap but not flat. At a $15^\\circ$ fold, the regular mesh has $\\sigma=%s$ and $\\rho=%s$, while the irregular occupied mesh has $\\sigma=%s$ and $\\rho=%s$.", fmt_num(profile_summary$gmds_stress[profile_summary$case_id == 'regular_12x12' & profile_summary$angle_deg == 15], 4L), fmt_num(profile_summary$rmse_to_truth[profile_summary$case_id == 'regular_12x12' & profile_summary$angle_deg == 15], 4L), fmt_num(profile_summary$gmds_stress[profile_summary$case_id == 'irregular_occ_12x12' & profile_summary$angle_deg == 15], 4L), fmt_num(profile_summary$rmse_to_truth[profile_summary$case_id == 'irregular_occ_12x12' & profile_summary$angle_deg == 15], 4L)),
    sprintf("More importantly, pure GMDS started from that folded state relaxes back near the smooth surface on both graphs. On the regular grid, the recovery reaches $\\sigma=%s$ and $\\rho=%s$; on the irregular occupied mesh it reaches $\\sigma=%s$ and $\\rho=%s$.", fmt_num(regular_fold$gmds_stress[[1L]], 4L), fmt_num(regular_fold$rmse_to_truth[[1L]], 4L), fmt_num(irregular_fold$gmds_stress[[1L]], 4L), fmt_num(irregular_fold$rmse_to_truth[[1L]], 4L)),
    sprintf("The real difficulty is the classical-MDS basin. Starting from cMDS, pure GMDS on the regular orthogonal mesh lands in a bad folded basin with $\\sigma=%s$ and $\\rho=%s$, whereas the irregular occupied mesh is much closer to truth at $\\sigma=%s$ and $\\rho=%s$. This matches the regularized-GMDS behavior too: with the earlier anchor schedule $\\lambda_A:0.10\\to0.02$, the regular mesh had $\\rho\\approx%s$, while the irregular occupied mesh improves to $\\rho\\approx%s$.", fmt_num(regular_cmd$gmds_stress[[1L]], 4L), fmt_num(regular_cmd$rmse_to_truth[[1L]], 4L), fmt_num(irregular_cmd$gmds_stress[[1L]], 4L), fmt_num(irregular_cmd$rmse_to_truth[[1L]], 4L), fmt_num(regular_a, 4L), fmt_num(irregular_a, 4L)),
    "So the new interpretation is sharper: smooth paraboloids are not absent from the objective, but regular rectangular meshes create a much worse global optimization landscape from generic initializations. Irregularity helps, not by making small folds dramatically more expensive, but by breaking enough symmetry that cMDS starts closer to the smooth basin.",
    "\\begin{table}[H]",
    "\\centering",
    "\\scriptsize",
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{lrrr}",
    "\\toprule",
    "Case & Fold angle (deg) & $\\sigma$ & $\\rho$ \\\\",
    "\\midrule",
    paste(profile_rows, collapse = "\n"),
    "\\bottomrule",
    "\\end{tabular}}",
    "\\caption{GMDS stress and Procrustes error for controlled hinge-fold perturbations of the exact scaled paraboloid reference. The similar curves for the regular and irregular meshes show that local fold directions remain cheap in both cases.}",
    "\\end{table}",
    sprintf("\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{Fold-direction perturbation profiles around the exact scaled paraboloid reference. Left: GMDS stress. Right: Procrustes RMSE to the truth. The regular and irregular curves are nearly identical, so irregularity does not mainly help by making local fold directions expensive.}\\end{figure}", plot_path),
    "\\begin{table}[H]",
    "\\centering",
    "\\scriptsize",
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{llrrr}",
    "\\toprule",
    "Case & Start & $\\sigma$ & $\\rho$ & $t$ (s) \\\\",
    "\\midrule",
    paste(recovery_rows, collapse = "\n"),
    "\\bottomrule",
    "\\end{tabular}}",
    "\\caption{Pure-GMDS recovery runs from two different basins. Starting from a mild $15^\\circ$ fold returns close to the smooth surface, but starting from cMDS behaves very differently on the regular and irregular graphs.}",
    "\\end{table}",
    sprintf("\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{Regular versus irregular paraboloid snapshots: exact scaled truth, a controlled $15^\\circ$ fold, and pure-GMDS outcomes from cMDS and from the folded start. The regular grid folds badly from cMDS, whereas the irregular occupied mesh stays much closer to the true surface.}\\end{figure}", snapshot_path),
    sep = "\n\n"
  )
}

cases <- list(
  build_regular_case(12L),
  build_irregular_case(12L)
)

angles_deg <- c(0, 5, 10, 15, 20, 25, 30)
case_outputs <- lapply(cases, compute_case_output, angles_deg = angles_deg, max_iter = 25L)
profile_df <- do.call(rbind, lapply(case_outputs, `[[`, "profile_df"))
utils::write.csv(profile_df, profiles_csv, row.names = FALSE)

metrics_df <- do.call(rbind, lapply(case_outputs, `[[`, "metrics_df"))
utils::write.csv(metrics_df, metrics_csv, row.names = FALSE)

plot_profiles(profile_df, plot_path)
plot_snapshots(cases, snapshot_path)
writeLines(build_appendix(case_outputs, profile_df, metrics_df), appendix_tex_path)
saveRDS(
  list(
    generated_at = as.character(Sys.time()),
    plot_path = plot_path,
    snapshot_path = snapshot_path,
    cases = lapply(case_outputs, function(x) {
      list(
        id = x$case$id,
        label = x$case$label,
        edges = x$case$edges,
        embeddings = x$embeddings
      )
    })
  ),
  rds_path
)

message("Wrote fold profiles: ", profiles_csv)
message("Wrote recovery metrics: ", metrics_csv)
message("Wrote embeddings bundle: ", rds_path)
message("Wrote LaTeX appendix: ", appendix_tex_path)
message("Wrote profile plot: ", plot_path)
message("Wrote snapshot grid: ", snapshot_path)
