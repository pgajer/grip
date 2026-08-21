#!/usr/bin/env Rscript

run_tag <- "gmds-mds-phase1-2026-03-31"
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
  stop("Install 'devtools' or the 'grip' package to run the GMDS/MDS Phase 1 benchmark.")
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

fmt_ratio <- function(x, digits = 2L) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}

tex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x, perl = TRUE)
  x
}

reference_geometry_label <- function(truth_label) {
  truth_label <- tolower(as.character(truth_label[[1L]]))
  switch(
    truth_label,
    parameter = "Reference geometry (u,v)",
    surface = "Reference geometry (surface)",
    sprintf("Reference geometry (%s)", truth_label)
  )
}

short_case_label <- function(label) {
  sub(" \\((2D|3D)\\)$", "", as.character(label))
}

pick_truth_coords <- function(bundle, dim) {
  if (dim == 2L) {
    if (!is.null(bundle$coords_param)) {
      return(as.matrix(bundle$coords_param))
    }
    if (!is.null(bundle$coords)) {
      return(as.matrix(bundle$coords))
    }
  }
  if (dim == 3L && !is.null(bundle$coords_surface)) {
    return(as.matrix(bundle$coords_surface))
  }
  NULL
}

new_case <- function(id,
                     label,
                     case_group,
                     bundle,
                     dim,
                     family,
                     truth_label,
                     comparison_key = id,
                     gmds_max_iter = 24L) {
  truth.coords <- pick_truth_coords(bundle, dim)
  if (is.null(truth.coords)) {
    stop(sprintf("No compatible truth coordinates found for case '%s'", id))
  }
  list(
    id = id,
    label = label,
    case_group = case_group,
    comparison_key = comparison_key,
    family = family,
    connectivity = if (!is.null(bundle$connectivity)) as.character(bundle$connectivity) else "orthogonal",
    dim = as.integer(dim),
    edges = bundle$edges,
    edge_weights = bundle$edge_weights,
    n = as.integer(bundle$n),
    truth_coords = truth.coords,
    truth_label = truth_label,
    gmds_max_iter = as.integer(gmds_max_iter)
  )
}

build_phase1_cases <- function() {
  list(
    new_case(
      id = "mesh_flat_10x10_2d",
      label = "Flat mesh 10x10 (2D)",
      case_group = "flat_regular",
      bundle = mesh.surface.graph(10, 10, surface = "saddle", amplitude = 0, normalize = "median"),
      dim = 2L,
      family = "mesh",
      truth_label = "parameter"
    ),
    new_case(
      id = "mesh_flat_16x16_2d",
      label = "Flat mesh 16x16 (2D)",
      case_group = "flat_regular",
      bundle = mesh.surface.graph(16, 16, surface = "saddle", amplitude = 0, normalize = "median"),
      dim = 2L,
      family = "mesh",
      truth_label = "parameter"
    ),
    new_case(
      id = "mesh_saddle_12x12_2d",
      label = "Saddle mesh 12x12 (2D)",
      case_group = "curved_regular",
      bundle = mesh.surface.graph(12, 12, surface = "saddle", amplitude = 0.35, normalize = "median"),
      dim = 2L,
      family = "mesh",
      truth_label = "parameter"
    ),
    new_case(
      id = "mesh_paraboloid_12x12_2d",
      label = "Paraboloid mesh 12x12 (2D)",
      case_group = "curved_regular",
      bundle = mesh.surface.graph(12, 12, surface = "paraboloid", amplitude = 0.35, normalize = "median"),
      dim = 2L,
      family = "mesh",
      truth_label = "parameter"
    ),
    new_case(
      id = "mesh_ripple_12x12_2d",
      label = "Ripple mesh 12x12 (2D)",
      case_group = "curved_regular",
      bundle = mesh.surface.graph(12, 12, surface = "ripple", amplitude = 0.50, freq_u = 2, freq_v = 2, normalize = "median"),
      dim = 2L,
      family = "mesh",
      truth_label = "parameter"
    ),
    new_case(
      id = "mesh_saddle_12x12_3d",
      label = "Saddle mesh 12x12 (3D)",
      case_group = "curved_regular_3d",
      bundle = mesh.surface.graph(12, 12, surface = "saddle", amplitude = 0.35, normalize = "median"),
      dim = 3L,
      family = "mesh",
      truth_label = "surface"
    ),
    new_case(
      id = "mesh_ripple_12x12_3d",
      label = "Ripple mesh 12x12 (3D)",
      case_group = "curved_regular_3d",
      bundle = mesh.surface.graph(12, 12, surface = "ripple", amplitude = 0.50, freq_u = 2, freq_v = 2, normalize = "median"),
      dim = 3L,
      family = "mesh",
      truth_label = "surface"
    ),
    new_case(
      id = "holes_flat_16x16_2d",
      label = "Periodic holes 16x16 (2D)",
      case_group = "flat_perforated",
      bundle = occupied.mesh.surface.graph(keep.periodic.holes(16, 16), surface = "saddle", amplitude = 0, normalize = "median"),
      dim = 2L,
      family = "occupied.mesh",
      truth_label = "parameter"
    ),
    new_case(
      id = "windows_flat_16x16_2d",
      label = "Staggered windows 16x16 (2D)",
      case_group = "flat_perforated",
      bundle = occupied.mesh.surface.graph(keep.staggered.windows(16, 16), surface = "saddle", amplitude = 0, normalize = "median"),
      dim = 2L,
      family = "occupied.mesh",
      truth_label = "parameter"
    ),
    new_case(
      id = "slit_flat_16x16_2d",
      label = "Slit channels 16x16 (2D)",
      case_group = "flat_perforated",
      bundle = occupied.mesh.surface.graph(
        keep.slit.channels(16, 16, orientation = "vertical", slit_period = 5, slit_width = 1, bridge_spacing = 4, bridge_size = 1, offset = 2),
        surface = "saddle",
        amplitude = 0,
        normalize = "median"
      ),
      dim = 2L,
      family = "occupied.mesh",
      truth_label = "parameter"
    ),
    new_case(
      id = "notches_flat_16x16_2d",
      label = "Asymmetric notches 16x16 (2D)",
      case_group = "flat_perforated",
      bundle = occupied.mesh.surface.graph(keep.asymmetric.notches(16, 16), surface = "saddle", amplitude = 0, normalize = "median"),
      dim = 2L,
      family = "occupied.mesh",
      truth_label = "parameter"
    ),
    new_case(
      id = "holes_ripple_14x14_2d",
      label = "Periodic holes ripple 14x14 (2D)",
      case_group = "curved_perforated",
      bundle = occupied.mesh.surface.graph(
        keep.periodic.holes(14, 14),
        surface = "ripple",
        amplitude = 0.50,
        freq_u = 2,
        freq_v = 2,
        normalize = "median"
      ),
      dim = 2L,
      family = "occupied.mesh",
      truth_label = "parameter"
    ),
    new_case(
      id = "slit_ripple_14x14_2d",
      label = "Slit channels ripple 14x14 (2D)",
      case_group = "curved_perforated",
      bundle = occupied.mesh.surface.graph(
        keep.slit.channels(14, 14, orientation = "vertical", slit_period = 5, slit_width = 1, bridge_spacing = 3, bridge_size = 1, offset = 2),
        surface = "ripple",
        amplitude = 0.50,
        freq_u = 2,
        freq_v = 2,
        normalize = "median"
      ),
      dim = 2L,
      family = "occupied.mesh",
      truth_label = "parameter"
    )
  )
}

build_connectivity_cases <- function() {
  list(
    new_case(
      id = "mesh_flat_10x10_orthogonal_2d",
      label = "Flat mesh 10x10 (2D)",
      case_group = "connectivity_refinement",
      comparison_key = "mesh_flat_10x10_2d",
      bundle = mesh.surface.graph(
        10, 10,
        surface = "saddle",
        amplitude = 0,
        connectivity = "orthogonal",
        normalize = "median"
      ),
      dim = 2L,
      family = "mesh",
      truth_label = "parameter"
    ),
    new_case(
      id = "mesh_flat_10x10_diagonal_2d",
      label = "Flat mesh 10x10 (2D)",
      case_group = "connectivity_refinement",
      comparison_key = "mesh_flat_10x10_2d",
      bundle = mesh.surface.graph(
        10, 10,
        surface = "saddle",
        amplitude = 0,
        connectivity = "diagonal",
        normalize = "median"
      ),
      dim = 2L,
      family = "mesh",
      truth_label = "parameter"
    ),
    new_case(
      id = "mesh_ripple_12x12_orthogonal_2d",
      label = "Ripple mesh 12x12 (2D)",
      case_group = "connectivity_refinement",
      comparison_key = "mesh_ripple_12x12_2d",
      bundle = mesh.surface.graph(
        12, 12,
        surface = "ripple",
        amplitude = 0.50,
        freq_u = 2,
        freq_v = 2,
        connectivity = "orthogonal",
        normalize = "median"
      ),
      dim = 2L,
      family = "mesh",
      truth_label = "parameter"
    ),
    new_case(
      id = "mesh_ripple_12x12_diagonal_2d",
      label = "Ripple mesh 12x12 (2D)",
      case_group = "connectivity_refinement",
      comparison_key = "mesh_ripple_12x12_2d",
      bundle = mesh.surface.graph(
        12, 12,
        surface = "ripple",
        amplitude = 0.50,
        freq_u = 2,
        freq_v = 2,
        connectivity = "diagonal",
        normalize = "median"
      ),
      dim = 2L,
      family = "mesh",
      truth_label = "parameter"
    ),
    new_case(
      id = "holes_flat_16x16_orthogonal_2d",
      label = "Periodic holes 16x16 (2D)",
      case_group = "connectivity_refinement",
      comparison_key = "holes_flat_16x16_2d",
      bundle = occupied.mesh.surface.graph(
        keep.periodic.holes(16, 16),
        surface = "saddle",
        amplitude = 0,
        connectivity = "orthogonal",
        normalize = "median"
      ),
      dim = 2L,
      family = "occupied.mesh",
      truth_label = "parameter"
    ),
    new_case(
      id = "holes_flat_16x16_diagonal_2d",
      label = "Periodic holes 16x16 (2D)",
      case_group = "connectivity_refinement",
      comparison_key = "holes_flat_16x16_2d",
      bundle = occupied.mesh.surface.graph(
        keep.periodic.holes(16, 16),
        surface = "saddle",
        amplitude = 0,
        connectivity = "diagonal",
        normalize = "median"
      ),
      dim = 2L,
      family = "occupied.mesh",
      truth_label = "parameter"
    ),
    new_case(
      id = "slit_flat_16x16_orthogonal_2d",
      label = "Slit channels 16x16 (2D)",
      case_group = "connectivity_refinement",
      comparison_key = "slit_flat_16x16_2d",
      bundle = occupied.mesh.surface.graph(
        keep.slit.channels(
          16, 16,
          orientation = "vertical",
          slit_period = 5,
          slit_width = 1,
          bridge_spacing = 4,
          bridge_size = 1,
          offset = 2
        ),
        surface = "saddle",
        amplitude = 0,
        connectivity = "orthogonal",
        normalize = "median"
      ),
      dim = 2L,
      family = "occupied.mesh",
      truth_label = "parameter"
    ),
    new_case(
      id = "slit_flat_16x16_diagonal_2d",
      label = "Slit channels 16x16 (2D)",
      case_group = "connectivity_refinement",
      comparison_key = "slit_flat_16x16_2d",
      bundle = occupied.mesh.surface.graph(
        keep.slit.channels(
          16, 16,
          orientation = "vertical",
          slit_period = 5,
          slit_width = 1,
          bridge_spacing = 4,
          bridge_size = 1,
          offset = 2
        ),
        surface = "saddle",
        amplitude = 0,
        connectivity = "diagonal",
        normalize = "median"
      ),
      dim = 2L,
      family = "occupied.mesh",
      truth_label = "parameter"
    )
  )
}

build_ripple_3d_density_cases <- function() {
  sizes <- c(12L, 15L, 20L)
  iter_budgets <- c(24L, 20L, 16L)
  lapply(seq_along(sizes), function(idx) {
    side <- sizes[[idx]]
    new_case(
      id = sprintf("mesh_ripple_%dx%d_3d_density", side, side),
      label = sprintf("Ripple mesh %dx%d (3D)", side, side),
      case_group = "ripple_surface_density_3d",
      bundle = mesh.surface.graph(
        side, side,
        surface = "ripple",
        amplitude = 0.50,
        freq_u = 2,
        freq_v = 2,
        normalize = "median"
      ),
      dim = 3L,
      family = "mesh",
      truth_label = "surface",
      gmds_max_iter = iter_budgets[[idx]]
    )
  })
}

build_paraboloid_3d_density_cases <- function() {
  sizes <- c(12L, 15L, 20L)
  iter_budgets <- c(24L, 20L, 16L)
  lapply(seq_along(sizes), function(idx) {
    side <- sizes[[idx]]
    new_case(
      id = sprintf("mesh_paraboloid_%dx%d_3d_density", side, side),
      label = sprintf("Paraboloid mesh %dx%d (3D)", side, side),
      case_group = "paraboloid_surface_density_3d",
      bundle = mesh.surface.graph(
        side, side,
        surface = "paraboloid",
        amplitude = 0.35,
        normalize = "median"
      ),
      dim = 3L,
      family = "mesh",
      truth_label = "surface",
      gmds_max_iter = iter_budgets[[idx]]
    )
  })
}

evaluate_method <- function(case,
                            prepared,
                            method_id,
                            method_label,
                            coords,
                            elapsed_sec,
                            n_threads_used = NA_integer_,
                            cmd_meta = NULL,
                            trace = NULL) {
  gmds <- grip.score.geodesic.mds(coords, prepared = prepared)
  classical <- grip:::grip.classical.mds.score.stats(coords, prepared)
  fit <- grip:::grip.align.to.target.nd(coords, case$truth_coords, allow.reflection = TRUE)

  data.frame(
    case_id = case$id,
    case_label = case$label,
    case_group = case$case_group,
    comparison_key = case$comparison_key,
    family = case$family,
    connectivity = case$connectivity,
    n = case$n,
    dim = case$dim,
    truth_label = case$truth_label,
    method_id = method_id,
    method_label = method_label,
    elapsed_sec = as.double(elapsed_sec),
    n_threads_used = as.integer(n_threads_used),
    gmds_raw_stress = gmds$gmds.raw_stress[[1L]],
    gmds_stress = gmds$gmds.stress[[1L]],
    gmds_mean_abs_path_error = gmds$gmds.mean.abs.path.error[[1L]],
    gmds_mean_rel_path_error = gmds$gmds.mean.rel.path.error[[1L]],
    euclidean_raw_stress = classical$raw_stress,
    euclidean_stress = classical$stress,
    euclidean_mean_abs_error = classical$mean.abs.error,
    euclidean_mean_rel_error = classical$mean.rel.error,
    graph_euclid_pearson = classical$pearson,
    graph_euclid_spearman = classical$spearman,
    procrustes_rmse = fit$rmse,
    positive_eigen_fraction = if (!is.null(cmd_meta)) cmd_meta$positive_eigen_fraction else NA_real_,
    negative_eigen_fraction = if (!is.null(cmd_meta)) cmd_meta$negative_eigen_fraction else NA_real_,
    additive_constant = if (!is.null(cmd_meta)) cmd_meta$additive_constant else NA_real_,
    gmds_iterations = if (!is.null(trace) && nrow(trace) > 0L) max(trace$iteration) else 0L,
    stringsAsFactors = FALSE
  )
}

run_case <- function(case) {
  prepared <- prepare.geodesic.kk(
    edges = case$edges,
    n = case$n,
    edge_weights = case$edge_weights
  )

  started <- proc.time()[["elapsed"]]
  cmd <- grip:::grip.classical.mds.embedding(prepared, dim = case$dim, eig = TRUE)
  cmd.elapsed <- proc.time()[["elapsed"]] - started

  started <- proc.time()[["elapsed"]]
  gmds <- grip.optimize.geodesic.mds(
    coords = cmd$coords,
    prepared = prepared,
    max_iter = case$gmds_max_iter,
    engine = "cpp",
    return_trace = TRUE
  )
  gmds.elapsed <- proc.time()[["elapsed"]] - started

  metrics <- rbind(
    evaluate_method(
      case = case,
      prepared = prepared,
      method_id = "cmdscale",
      method_label = "Classical MDS (cmdscale)",
      coords = cmd$coords,
      elapsed_sec = cmd.elapsed,
      n_threads_used = 1L,
      cmd_meta = cmd,
      trace = NULL
    ),
    evaluate_method(
      case = case,
      prepared = prepared,
      method_id = "gmds",
      method_label = "GMDS (cpp)",
      coords = gmds$coords,
      elapsed_sec = gmds.elapsed,
      n_threads_used = gmds$n_threads_used,
      cmd_meta = cmd,
      trace = gmds$trace
    )
  )

  list(
    case = case,
    metrics = metrics,
    layouts = list(
      truth = case$truth_coords,
      cmdscale = cmd$coords,
      gmds = gmds$coords
    )
  )
}

plot_layout_panel <- function(coords, edges, title) {
  coords <- as.matrix(coords)
  plot(
    coords,
    asp = 1,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = title,
    pch = 16,
    cex = 0.25
  )
  segments(
    coords[edges[, 1L], 1L],
    coords[edges[, 1L], 2L],
    coords[edges[, 2L], 1L],
    coords[edges[, 2L], 2L],
    col = "#00000033",
    lwd = 0.4
  )
  points(coords, pch = 16, cex = 0.25)
}

save_triptych <- function(case_result, summary_row, path) {
  truth <- case_result$layouts$truth
  if (ncol(truth) != 2L) {
    return(invisible(NULL))
  }
  truth.norm <- grip:::grip.normalize.coords(truth)
  cmd.aligned <- grip:::grip.align.to.target.nd(case_result$layouts$cmdscale, truth, allow.reflection = TRUE)$aligned
  gmds.aligned <- grip:::grip.align.to.target.nd(case_result$layouts$gmds, truth, allow.reflection = TRUE)$aligned

  grDevices::png(path, width = 1800, height = 600, res = 150)
  old.par <- graphics::par(mfrow = c(1L, 3L), mar = c(1.5, 1.5, 3, 1.5))
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)

  plot_layout_panel(
    truth.norm,
    case_result$case$edges,
    reference_geometry_label(case_result$case$truth_label)
  )
  plot_layout_panel(
    cmd.aligned,
    case_result$case$edges,
    sprintf(
      "cmdscale\nsigma %s, rho %s\nt %ss",
      fmt_num(summary_row$gmds_stress_cmdscale),
      fmt_num(summary_row$procrustes_rmse_cmdscale),
      fmt_time(summary_row$elapsed_sec_cmdscale)
    )
  )
  plot_layout_panel(
    gmds.aligned,
    case_result$case$edges,
    sprintf(
      "GMDS\nsigma %s, rho %s\nt %ss",
      fmt_num(summary_row$gmds_stress_gmds),
      fmt_num(summary_row$procrustes_rmse_gmds),
      fmt_time(summary_row$elapsed_sec_gmds)
    )
  )
}

save_connectivity_triptych <- function(orth_result,
                                       diag_result,
                                       comparison_row,
                                       path) {
  truth <- orth_result$layouts$truth
  if (ncol(truth) != 2L) {
    return(invisible(NULL))
  }
  truth.norm <- grip:::grip.normalize.coords(truth)
  orth.aligned <- grip:::grip.align.to.target.nd(
    orth_result$layouts$gmds,
    truth,
    allow.reflection = TRUE
  )$aligned
  diag.aligned <- grip:::grip.align.to.target.nd(
    diag_result$layouts$gmds,
    truth,
    allow.reflection = TRUE
  )$aligned

  grDevices::png(path, width = 1800, height = 600, res = 150)
  old.par <- graphics::par(mfrow = c(1L, 3L), mar = c(1.5, 1.5, 3, 1.5))
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)

  plot_layout_panel(
    truth.norm,
    orth_result$case$edges,
    reference_geometry_label(orth_result$case$truth_label)
  )
  plot_layout_panel(
    orth.aligned,
    orth_result$case$edges,
    sprintf(
      "Orthogonal GMDS\nsigma %s, rho %s\nt %ss",
      fmt_num(comparison_row$gmds_stress_orth[[1L]]),
      fmt_num(comparison_row$procrustes_rmse_orth[[1L]]),
      fmt_time(comparison_row$elapsed_sec_orth[[1L]])
    )
  )
  plot_layout_panel(
    diag.aligned,
    diag_result$case$edges,
    sprintf(
      "Diagonal GMDS\nsigma %s, rho %s\nt %ss",
      fmt_num(comparison_row$gmds_stress_diag[[1L]]),
      fmt_num(comparison_row$procrustes_rmse_diag[[1L]]),
      fmt_time(comparison_row$elapsed_sec_diag[[1L]])
    )
  )
}

save_named_triptych <- function(case,
                                truth,
                                middle_coords,
                                middle_title,
                                right_coords,
                                right_title,
                                path) {
  if (ncol(truth) != 2L) {
    return(invisible(NULL))
  }
  truth.norm <- grip:::grip.normalize.coords(truth)
  middle.aligned <- grip:::grip.align.to.target.nd(
    middle_coords,
    truth,
    allow.reflection = TRUE
  )$aligned
  right.aligned <- grip:::grip.align.to.target.nd(
    right_coords,
    truth,
    allow.reflection = TRUE
  )$aligned

  grDevices::png(path, width = 1800, height = 600, res = 150)
  old.par <- graphics::par(mfrow = c(1L, 3L), mar = c(1.5, 1.5, 3, 1.5))
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)

  plot_layout_panel(truth.norm, case$edges, reference_geometry_label(case$truth_label))
  plot_layout_panel(middle.aligned, case$edges, middle_title)
  plot_layout_panel(right.aligned, case$edges, right_title)
}

save_named_quadtych <- function(case,
                                truth,
                                coords_list,
                                title_list,
                                path) {
  if (ncol(truth) != 2L) {
    return(invisible(NULL))
  }
  if (length(coords_list) != 3L || length(title_list) != 3L) {
    stop("coords_list and title_list must both have length 3")
  }

  truth.norm <- grip:::grip.normalize.coords(truth)
  aligned.list <- lapply(coords_list, function(coords) {
    grip:::grip.align.to.target.nd(
      coords,
      truth,
      allow.reflection = TRUE
    )$aligned
  })

  grDevices::png(path, width = 2400, height = 600, res = 150)
  old.par <- graphics::par(mfrow = c(1L, 4L), mar = c(1.5, 1.5, 3, 1.5))
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)

  plot_layout_panel(truth.norm, case$edges, reference_geometry_label(case$truth_label))
  for (idx in seq_along(aligned.list)) {
    plot_layout_panel(aligned.list[[idx]], case$edges, title_list[[idx]])
  }
}

save_surface_3d_density_grid <- function(case_results,
                                         metrics.df,
                                         path,
                                         heading,
                                         azimuth = 35,
                                         elevation = 24) {
  if (length(case_results) == 0L) {
    return(invisible(NULL))
  }

  grDevices::png(
    path,
    width = 3200,
    height = max(1200L, 1000L * length(case_results)),
    res = 180,
    bg = "#ffffff"
  )
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(
    mfrow = c(length(case_results), 4L),
    mar = c(1.2, 1.2, 2.9, 0.5),
    oma = c(0, 0, 1.8, 0)
  )

  method_ids <- c("cmdscale_average", "gmds_average_50", "gmds_average_lin_50")

  for (case_result in case_results) {
    case <- case_result$case
    truth <- case$truth_coords
    short.label <- short_case_label(case$label)

    plot.layout(
      coords = truth,
      edges = case$edges,
      projection = "ortho",
      azimuth = azimuth,
      elevation = elevation,
      vertex.col = "#bc6c25",
      edge.col = "#adb5bd",
      main = ""
    )
    graphics::title(
      main = sprintf("%s\nReference surface", short.label),
      cex.main = 0.82
    )

    for (idx in seq_along(method_ids)) {
      row <- metrics.df[
        metrics.df$case_id == case$id &
          metrics.df$method_id == method_ids[[idx]],
        ,
        drop = FALSE
      ]
      coords <- switch(
        method_ids[[idx]],
        cmdscale_average = case_result$layouts$cmdscale,
        gmds_average_50 = case_result$layouts$gmds_average,
        gmds_average_lin_50 = case_result$layouts$gmds_average_linear_tether
      )
      coords <- grip:::grip.align.to.target.nd(
        coords,
        truth,
        allow.reflection = TRUE
      )$aligned
      plot.layout(
        coords = coords,
        edges = case$edges,
        projection = "ortho",
        azimuth = azimuth,
        elevation = elevation,
        vertex.col = "#3a5a40",
        edge.col = "#adb5bd",
        main = ""
      )
      graphics::title(
        main = sprintf(
          "%s\nsigma %s, rho %s\nt %ss",
          row$setting_label[[1L]],
          fmt_num(row$gmds_stress[[1L]]),
          fmt_num(row$procrustes_rmse[[1L]]),
          fmt_time(row$elapsed_sec[[1L]])
        ),
        cex.main = 0.76
      )
    }
  }

  graphics::mtext(
    heading,
    outer = TRUE,
    cex = 1.15,
    font = 2
  )
  invisible(NULL)
}

write_case_longtable <- function(summary.df) {
  rows <- apply(summary.df, 1L, function(row) {
    paste(
      tex_escape(row[["case_label"]]),
      row[["n"]],
      row[["dim"]],
      fmt_num(as.numeric(row[["gmds_stress_cmdscale"]])),
      fmt_num(as.numeric(row[["gmds_stress_gmds"]])),
      fmt_num(as.numeric(row[["gmds_improvement_pct"]]), digits = 1L),
      fmt_num(as.numeric(row[["procrustes_rmse_cmdscale"]])),
      fmt_num(as.numeric(row[["procrustes_rmse_gmds"]])),
      sep = " & "
    )
  })
  paste0(rows, " \\\\")
}

write_group_table <- function(group.df) {
  rows <- apply(group.df, 1L, function(row) {
    paste(
      tex_escape(row[["case_group"]]),
      row[["n_cases"]],
      fmt_num(as.numeric(row[["mean_gmds_improvement_pct"]]), digits = 1L),
      fmt_num(as.numeric(row[["mean_gmds_stress_cmdscale"]])),
      fmt_num(as.numeric(row[["mean_gmds_stress_gmds"]])),
      fmt_num(as.numeric(row[["mean_procrustes_rmse_cmdscale"]])),
      fmt_num(as.numeric(row[["mean_procrustes_rmse_gmds"]])),
      sep = " & "
    )
  })
  paste0(rows, " \\\\")
}

write_connectivity_table <- function(connectivity.df) {
  rows <- apply(connectivity.df, 1L, function(row) {
    paste(
      tex_escape(row[["case_label"]]),
      row[["n"]],
      fmt_num(as.numeric(row[["gmds_stress_orth"]])),
      fmt_num(as.numeric(row[["gmds_stress_diag"]])),
      fmt_num(as.numeric(row[["procrustes_rmse_orth"]])),
      fmt_num(as.numeric(row[["procrustes_rmse_diag"]])),
      fmt_num(as.numeric(row[["procrustes_improvement_pct"]]), digits = 1L),
      sep = " & "
    )
  })
  paste0(rows, " \\\\")
}

run_post_review_case <- function(case,
                                 single_iter = 50L,
                                 average_iter = 50L,
                                 tether_weight = 0.05,
                                 tether_weight_end = 0) {
  prepared.single <- prepare.geodesic.kk(
    edges = case$edges,
    n = case$n,
    edge_weights = case$edge_weights,
    tie_mode = "single"
  )
  prepared.average <- prepare.geodesic.kk(
    edges = case$edges,
    n = case$n,
    edge_weights = case$edge_weights,
    tie_mode = "average"
  )

  started <- proc.time()[["elapsed"]]
  cmd <- grip:::grip.classical.mds.embedding(prepared.single, dim = case$dim, eig = TRUE)
  cmd.elapsed <- proc.time()[["elapsed"]] - started

  rows <- list()
  layouts <- list(
    truth = case$truth_coords,
    cmdscale = cmd$coords
  )

  add_row <- function(row,
                      method_id,
                      setting_label,
                      tie_mode,
                      iter_budget,
                      anchor_start = 0,
                      anchor_end = 0,
                      continuation = "none") {
    row$method_id <- method_id
    row$setting_label <- setting_label
    row$tie_mode <- tie_mode
    row$iter_budget <- as.integer(iter_budget)
    row$anchor_weight_start <- as.double(anchor_start)
    row$anchor_weight_end <- as.double(anchor_end)
    row$continuation <- continuation
    row$case_short <- short_case_label(case$label)
    rows[[length(rows) + 1L]] <<- row
  }

  add_row(
    evaluate_method(
      case = case,
      prepared = prepared.single,
      method_id = "cmdscale_single",
      method_label = "cmdscale (single-path objective)",
      coords = cmd$coords,
      elapsed_sec = cmd.elapsed,
      n_threads_used = 1L,
      cmd_meta = cmd,
      trace = NULL
    ),
    method_id = "cmdscale_single",
    setting_label = "cmdscale (single)",
    tie_mode = "single",
    iter_budget = 0L
  )

  add_row(
    evaluate_method(
      case = case,
      prepared = prepared.average,
      method_id = "cmdscale_average",
      method_label = "cmdscale (tie-averaged objective)",
      coords = cmd$coords,
      elapsed_sec = cmd.elapsed,
      n_threads_used = 1L,
      cmd_meta = cmd,
      trace = NULL
    ),
    method_id = "cmdscale_average",
    setting_label = "cmdscale (avg)",
    tie_mode = "average",
    iter_budget = 0L
  )

  started <- proc.time()[["elapsed"]]
  single.fit <- grip.optimize.geodesic.mds(
    coords = cmd$coords,
    prepared = prepared.single,
    max_iter = single_iter,
    engine = "cpp",
    return_trace = TRUE
  )
  single.elapsed <- proc.time()[["elapsed"]] - started
  layouts$gmds_single <- single.fit$coords
  add_row(
    evaluate_method(
      case = case,
      prepared = prepared.single,
      method_id = "gmds_single_50",
      method_label = sprintf("GMDS single-path (%d iter)", single_iter),
      coords = single.fit$coords,
      elapsed_sec = single.elapsed,
      n_threads_used = single.fit$n_threads_used,
      cmd_meta = cmd,
      trace = single.fit$trace
    ),
    method_id = "gmds_single_50",
    setting_label = sprintf("single GMDS (%d)", single_iter),
    tie_mode = "single",
    iter_budget = single_iter
  )

  started <- proc.time()[["elapsed"]]
  average.fit <- grip.optimize.geodesic.mds(
    coords = cmd$coords,
    prepared = prepared.average,
    max_iter = average_iter,
    engine = "cpp",
    return_trace = TRUE
  )
  average.elapsed <- proc.time()[["elapsed"]] - started
  layouts$gmds_average <- average.fit$coords
  add_row(
    evaluate_method(
      case = case,
      prepared = prepared.average,
      method_id = "gmds_average_50",
      method_label = sprintf("GMDS tie-averaged (%d iter)", average_iter),
      coords = average.fit$coords,
      elapsed_sec = average.elapsed,
      n_threads_used = average.fit$n_threads_used,
      cmd_meta = cmd,
      trace = average.fit$trace
    ),
    method_id = "gmds_average_50",
    setting_label = sprintf("avg GMDS (%d)", average_iter),
    tie_mode = "average",
    iter_budget = average_iter
  )

  started <- proc.time()[["elapsed"]]
  tether.fit <- grip.optimize.geodesic.mds(
    coords = cmd$coords,
    prepared = prepared.average,
    max_iter = average_iter,
    engine = "cpp",
    return_trace = TRUE,
    anchor_mode = "cmdscale",
    anchor_weight = tether_weight,
    anchor_weight_end = tether_weight_end,
    continuation = "linear"
  )
  tether.elapsed <- proc.time()[["elapsed"]] - started
  layouts$gmds_average_linear_tether <- tether.fit$coords
  add_row(
    evaluate_method(
      case = case,
      prepared = prepared.average,
      method_id = "gmds_average_lin_50",
      method_label = sprintf(
        "GMDS tie-averaged + linear tether (%d iter)",
        average_iter
      ),
      coords = tether.fit$coords,
      elapsed_sec = tether.elapsed,
      n_threads_used = tether.fit$n_threads_used,
      cmd_meta = cmd,
      trace = tether.fit$trace
    ),
    method_id = "gmds_average_lin_50",
    setting_label = sprintf("avg+tether (%d)", average_iter),
    tie_mode = "average",
    iter_budget = average_iter,
    anchor_start = tether_weight,
    anchor_end = tether_weight_end,
    continuation = "linear"
  )

  if (identical(case$id, "mesh_flat_10x10_2d")) {
    started <- proc.time()[["elapsed"]]
    single.long.fit <- grip.optimize.geodesic.mds(
      coords = cmd$coords,
      prepared = prepared.single,
      max_iter = 200L,
      engine = "cpp",
      return_trace = TRUE
    )
    single.long.elapsed <- proc.time()[["elapsed"]] - started
    layouts$gmds_single_200 <- single.long.fit$coords
    add_row(
      evaluate_method(
        case = case,
        prepared = prepared.single,
        method_id = "gmds_single_200",
        method_label = "GMDS single-path (200 iter)",
        coords = single.long.fit$coords,
        elapsed_sec = single.long.elapsed,
        n_threads_used = single.long.fit$n_threads_used,
        cmd_meta = cmd,
        trace = single.long.fit$trace
      ),
      method_id = "gmds_single_200",
      setting_label = "single GMDS (200)",
      tie_mode = "single",
      iter_budget = 200L
    )
  }

  metrics <- do.call(rbind, rows)
  metrics <- metrics[order(metrics$case_id, metrics$tie_mode, metrics$iter_budget), , drop = FALSE]

  list(
    case = case,
    metrics = metrics,
    layouts = layouts
  )
}

run_density_3d_case <- function(case,
                                tether_weight = 0.05,
                                tether_weight_end = 0) {
  prepared.average <- prepare.geodesic.kk(
    edges = case$edges,
    n = case$n,
    edge_weights = case$edge_weights,
    tie_mode = "average"
  )

  started <- proc.time()[["elapsed"]]
  cmd <- grip:::grip.classical.mds.embedding(prepared.average, dim = case$dim, eig = TRUE)
  cmd.elapsed <- proc.time()[["elapsed"]] - started

  rows <- list()
  layouts <- list(
    truth = case$truth_coords,
    cmdscale = cmd$coords
  )

  add_row <- function(row,
                      method_id,
                      setting_label,
                      iter_budget,
                      anchor_start = 0,
                      anchor_end = 0,
                      continuation = "none") {
    row$method_id <- method_id
    row$setting_label <- setting_label
    row$tie_mode <- "average"
    row$iter_budget <- as.integer(iter_budget)
    row$anchor_weight_start <- as.double(anchor_start)
    row$anchor_weight_end <- as.double(anchor_end)
    row$continuation <- continuation
    row$case_short <- short_case_label(case$label)
    rows[[length(rows) + 1L]] <<- row
  }

  add_row(
    evaluate_method(
      case = case,
      prepared = prepared.average,
      method_id = "cmdscale_average",
      method_label = "cmdscale (tie-averaged objective)",
      coords = cmd$coords,
      elapsed_sec = cmd.elapsed,
      n_threads_used = 1L,
      cmd_meta = cmd,
      trace = NULL
    ),
    method_id = "cmdscale_average",
    setting_label = "cmdscale (avg)",
    iter_budget = 0L
  )

  started <- proc.time()[["elapsed"]]
  average.fit <- grip.optimize.geodesic.mds(
    coords = cmd$coords,
    prepared = prepared.average,
    max_iter = case$gmds_max_iter,
    engine = "cpp",
    return_trace = TRUE
  )
  average.elapsed <- proc.time()[["elapsed"]] - started
  layouts$gmds_average <- average.fit$coords
  add_row(
    evaluate_method(
      case = case,
      prepared = prepared.average,
      method_id = "gmds_average_50",
      method_label = sprintf("GMDS tie-averaged (%d iter)", case$gmds_max_iter),
      coords = average.fit$coords,
      elapsed_sec = average.elapsed,
      n_threads_used = average.fit$n_threads_used,
      cmd_meta = cmd,
      trace = average.fit$trace
    ),
    method_id = "gmds_average_50",
    setting_label = sprintf("avg GMDS (%d)", case$gmds_max_iter),
    iter_budget = case$gmds_max_iter
  )

  started <- proc.time()[["elapsed"]]
  tether.fit <- grip.optimize.geodesic.mds(
    coords = cmd$coords,
    prepared = prepared.average,
    max_iter = case$gmds_max_iter,
    engine = "cpp",
    return_trace = TRUE,
    anchor_mode = "cmdscale",
    anchor_weight = tether_weight,
    anchor_weight_end = tether_weight_end,
    continuation = "linear"
  )
  tether.elapsed <- proc.time()[["elapsed"]] - started
  layouts$gmds_average_linear_tether <- tether.fit$coords
  add_row(
    evaluate_method(
      case = case,
      prepared = prepared.average,
      method_id = "gmds_average_lin_50",
      method_label = sprintf(
        "GMDS tie-averaged + linear tether (%d iter)",
        case$gmds_max_iter
      ),
      coords = tether.fit$coords,
      elapsed_sec = tether.elapsed,
      n_threads_used = tether.fit$n_threads_used,
      cmd_meta = cmd,
      trace = tether.fit$trace
    ),
    method_id = "gmds_average_lin_50",
    setting_label = sprintf("avg+tether (%d)", case$gmds_max_iter),
    iter_budget = case$gmds_max_iter,
    anchor_start = tether_weight,
    anchor_end = tether_weight_end,
    continuation = "linear"
  )

  metrics <- do.call(rbind, rows)
  metrics <- metrics[order(metrics$case_id, metrics$iter_budget), , drop = FALSE]

  list(
    case = case,
    metrics = metrics,
    layouts = layouts
  )
}

write_post_review_table <- function(postreview.df) {
  rows <- apply(postreview.df, 1L, function(row) {
    paste(
      tex_escape(row[["case_short"]]),
      tex_escape(row[["setting_label"]]),
      fmt_time(as.numeric(row[["elapsed_sec"]])),
      fmt_num(as.numeric(row[["gmds_stress"]])),
      fmt_num(as.numeric(row[["procrustes_rmse"]])),
      sep = " & "
    )
  })
  paste0(rows, " \\\\")
}

write_cache_perf_table <- function(cache.df) {
  rows <- apply(cache.df, 1L, function(row) {
    paste(
      tex_escape(row[["case_short"]]),
      row[["n"]],
      row[["pair_count"]],
      row[["edge_terms"]],
      fmt_time(as.numeric(row[["elapsed_sec_r"]])),
      fmt_time(as.numeric(row[["elapsed_sec_cpp"]])),
      fmt_ratio(as.numeric(row[["speedup_cpp_vs_r"]])),
      sep = " & "
    )
  })
  paste0(rows, " \\\\")
}

write_optimizer_perf_table <- function(perf.df) {
  rows <- apply(perf.df, 1L, function(row) {
    paste(
      tex_escape(row[["case_short"]]),
      tex_escape(row[["setting_short"]]),
      fmt_time(as.numeric(row[["elapsed_sec_before"]])),
      fmt_time(as.numeric(row[["elapsed_sec_after"]])),
      fmt_ratio(as.numeric(row[["speedup_before_after"]])),
      row[["n_threads_used"]],
      sep = " & "
    )
  })
  paste0(rows, " \\\\")
}

benchmark_cache_build_case <- function(case) {
  base <- grip:::prepare.geodesic.kk.base(
    edges = case$edges,
    n = case$n,
    edge_weights = case$edge_weights,
    caller = "benchmark-gmds-mds-phase1"
  )
  pair.matrix <- grip:::grip.full.geodesic.kk.pair.matrix(base$n)

  started <- proc.time()[["elapsed"]]
  cache.r <- grip:::grip.build.tie.average.shortest.path.cache.r(
    pair.matrix = pair.matrix,
    adj.list = base$adj_list,
    weight.list = base$weight_list,
    dist.matrix = base$distance_matrix,
    parents = base$parents
  )
  elapsed.r <- proc.time()[["elapsed"]] - started

  started <- proc.time()[["elapsed"]]
  cache.cpp <- grip:::grip.build.tie.average.shortest.path.cache(
    pair.matrix = pair.matrix,
    adj.list = base$adj_list,
    weight.list = base$weight_list,
    dist.matrix = base$distance_matrix,
    parents = base$parents,
    cache_engine = "cpp"
  )
  elapsed.cpp <- proc.time()[["elapsed"]] - started

  stopifnot(
    isTRUE(all.equal(cache.cpp$pair_graph_distance, cache.r$pair_graph_distance, tolerance = 1e-10)),
    isTRUE(all.equal(cache.cpp$pair_path_count_log, cache.r$pair_path_count_log, tolerance = 1e-10))
  )

  data.frame(
    case_id = case$id,
    case_short = short_case_label(case$label),
    family = case$family,
    n = case$n,
    dim = case$dim,
    pair_count = nrow(pair.matrix),
    edge_terms = length(cache.cpp$flat_edge_u),
    elapsed_sec_r = as.double(elapsed.r),
    elapsed_sec_cpp = as.double(elapsed.cpp),
    speedup_cpp_vs_r = as.double(elapsed.r / pmax(elapsed.cpp, 1e-12)),
    stringsAsFactors = FALSE
  )
}

build_optimizer_perf_table <- function(current.df,
                                       baseline.paths) {
  baseline <- do.call(
    rbind,
    lapply(baseline.paths, utils::read.csv, stringsAsFactors = FALSE)
  )
  current <- current.df[
    current.df$method_id %in% c("gmds_average_50", "gmds_average_lin_50"),
    ,
    drop = FALSE
  ]
  baseline <- baseline[
    baseline$method_id %in% c("gmds_average_50", "gmds_average_lin_50"),
    ,
    drop = FALSE
  ]
  merged <- merge(
    baseline[, c("case_id", "method_id", "elapsed_sec")],
    current[, c("case_id", "method_id", "case_short", "setting_label", "elapsed_sec", "n_threads_used")],
    by = c("case_id", "method_id"),
    suffixes = c("_before", "_after"),
    sort = FALSE
  )
  merged$setting_short <- ifelse(
    merged$method_id == "gmds_average_50",
    "avg GMDS",
    "avg+tether"
  )
  merged$speedup_before_after <- merged$elapsed_sec_before / pmax(merged$elapsed_sec_after, 1e-12)
  case.order <- match(
    merged$case_id,
    c(
      "mesh_ripple_12x12_3d_density",
      "mesh_ripple_15x15_3d_density",
      "mesh_ripple_20x20_3d_density",
      "mesh_paraboloid_12x12_3d_density",
      "mesh_paraboloid_15x15_3d_density",
      "mesh_paraboloid_20x20_3d_density"
    )
  )
  method.order <- match(merged$method_id, c("gmds_average_50", "gmds_average_lin_50"))
  merged[order(case.order, method.order), , drop = FALSE]
}

write_report <- function(summary.df,
                         group.df,
                         connectivity.df,
                         postreview.df,
                         ripple3d.df,
                         paraboloid3d.df,
                         cacheperf.df,
                         optperf.df,
                         representative_ids,
                         connectivity_ids,
                         pdf_dir,
                         tex_path) {
  rep_sections <- unlist(lapply(representative_ids, function(case_id) {
    row <- summary.df[summary.df$case_id == case_id, , drop = FALSE]
    if (nrow(row) == 0L) {
      return(character(0L))
    }
    fig.name <- file.path(basename(pdf_dir), sprintf("%s-triptych.png", case_id))
    extra_text <- character(0L)
    if (identical(case_id, "mesh_ripple_12x12_2d")) {
      extra_text <- c(
        "This 2D ripple example should be interpreted only as a dimensionality-obstruction diagnostic. The mesh edge weights come from a lifted 3D ripple surface, so its graph-geodesic metric is not expected to admit a faithful planar realization. The parameter-plane reference is therefore a sampling reference, not a genuinely achievable isometric target.",
        "For geometric fidelity on ripple surfaces, the 3D ripple-density section below is the relevant comparison."
      )
    }
    c(
      sprintf("\\subsection{%s}", tex_escape(row$case_label[[1L]])),
      sprintf(
        "In the baseline single-path GMDS run, $\\sigma$ decreased from %s to %s, for $\\Delta\\sigma=%s\\%%$. The Procrustes error $\\rho$ moved from %s to %s.",
        fmt_num(row$gmds_stress_cmdscale[[1L]]),
        fmt_num(row$gmds_stress_gmds[[1L]]),
        fmt_num(row$gmds_improvement_pct[[1L]], digits = 1L),
        fmt_num(row$procrustes_rmse_cmdscale[[1L]]),
        fmt_num(row$procrustes_rmse_gmds[[1L]])
      ),
      extra_text,
      "\\begin{figure}[H]",
      "\\centering",
      sprintf("\\includegraphics[width=\\textwidth]{%s}", tex_escape(fig.name)),
      sprintf("\\caption{Reference geometry, classical MDS, and baseline single-path GMDS for %s.}", tex_escape(row$case_label[[1L]])),
      "\\end{figure}"
    )
  }))

  connectivity_sections <- unlist(lapply(connectivity_ids, function(case_id) {
    row <- connectivity.df[connectivity.df$comparison_key == case_id, , drop = FALSE]
    if (nrow(row) == 0L) {
      return(character(0L))
    }
    fig.name <- file.path(basename(pdf_dir), sprintf("%s-connectivity-triptych.png", case_id))
    c(
      sprintf("\\subsection{%s}", tex_escape(row$case_label[[1L]])),
      sprintf(
        "Using diagonal connectivity reduced $\\rho$ from %s to %s, a relative improvement of %s\\%%, while changing $\\sigma_G$ from %s to %s.",
        fmt_num(row$procrustes_rmse_orth[[1L]]),
        fmt_num(row$procrustes_rmse_diag[[1L]]),
        fmt_num(row$procrustes_improvement_pct[[1L]], digits = 1L),
        fmt_num(row$gmds_stress_orth[[1L]]),
        fmt_num(row$gmds_stress_diag[[1L]])
      ),
      "\\begin{figure}[H]",
      "\\centering",
      sprintf("\\includegraphics[width=\\textwidth]{%s}", tex_escape(fig.name)),
      sprintf("\\caption{Reference geometry, orthogonal-GMDS, and diagonal-GMDS for %s.}", tex_escape(row$case_label[[1L]])),
      "\\end{figure}"
    )
  }))

  post_lookup <- function(case_id, method_id) {
    postreview.df[
      postreview.df$case_id == case_id & postreview.df$method_id == method_id,
      ,
      drop = FALSE
    ]
  }
  post_num <- function(case_id, method_id, column, digits = 4L) {
    row <- post_lookup(case_id, method_id)
    if (nrow(row) == 0L) {
      return("NA")
    }
    fmt_num(row[[column]][[1L]], digits = digits)
  }
  surface_num <- function(df, case_id, method_id, column, digits = 4L) {
    row <- df[
      df$case_id == case_id & df$method_id == method_id,
      ,
      drop = FALSE
    ]
    if (nrow(row) == 0L) {
      return("NA")
    }
    fmt_num(row[[column]][[1L]], digits = digits)
  }
  surface_time <- function(df, case_id, method_id) {
    row <- df[
      df$case_id == case_id & df$method_id == method_id,
      ,
      drop = FALSE
    ]
    if (nrow(row) == 0L) {
      return("NA")
    }
    fmt_time(row$elapsed_sec[[1L]])
  }
  flat.base <- summary.df[summary.df$case_id == "mesh_flat_10x10_2d", , drop = FALSE]
  cache.best <- cacheperf.df[which.max(cacheperf.df$speedup_cpp_vs_r), , drop = FALSE]
  cache.worst <- cacheperf.df[which.min(cacheperf.df$speedup_cpp_vs_r), , drop = FALSE]
  opt.best <- optperf.df[which.max(optperf.df$speedup_before_after), , drop = FALSE]
  opt.worst <- optperf.df[which.min(optperf.df$speedup_before_after), , drop = FALSE]
  cache.median.speedup <- stats::median(cacheperf.df$speedup_cpp_vs_r)
  opt.median.speedup <- stats::median(optperf.df$speedup_before_after)

  tex <- c(
    "\\documentclass[11pt]{article}",
    "\\usepackage[margin=1in]{geometry}",
    "\\usepackage{graphicx}",
    "\\usepackage{booktabs}",
    "\\usepackage{longtable}",
    "\\usepackage{float}",
    "\\usepackage{amsmath}",
    "\\usepackage{hyperref}",
    "\\title{GMDS vs Classical MDS Comparison Report}",
    "\\author{GRIP Phase 1 benchmark}",
    "\\date{March 31, 2026}",
    "",
    "\\begin{document}",
    "\\maketitle",
    "",
    "\\begin{abstract}",
    "This report is the first populated version of the synthetic GMDS-versus-classical-MDS comparison suite proposed for the \\texttt{grip} package.",
    "At this stage only Phase~1 is implemented: regular meshes, curved meshes, and perforated meshes.",
    "This revision also records the main post-review modifications motivated by the Phase~1 figures: clearer reference-geometry terminology, an explicit single-path convergence diagnostic, tie-averaged GMDS on tied shortest-path families, an MDS-tethered continuation mode, and wall-clock runtime recording for each generated embedding.",
    "Future phases will extend the same report to wrapped surfaces, recursive families, triangulated manifolds, irregular samplings, porous 3D families, and intrinsic weighted trees.",
    "\\end{abstract}",
    "",
    "\\section{Protocol}",
    "Each Phase~1 baseline case was converted to a weighted graph using orthogonal mesh connectivity, all-pairs graph-geodesic distances were computed, classical MDS was fit with \\texttt{stats::cmdscale()}, and GMDS was run from the classical-MDS initialization with the compiled optimizer.",
    "Both embeddings were scored under the GMDS path objective and aligned to the synthetic truth coordinates for Procrustes RMSE.",
    "Wall-clock runtime $t$ was measured around each embedding call only, so the reported values exclude downstream scoring, Procrustes alignment, plotting, and LaTeX generation. These timings are therefore comparable within this report, but should still be interpreted as machine-dependent diagnostics rather than universal complexity claims.",
    "Throughout the report, $\\sigma$ denotes the fixed-path GMDS stress, $\\rho$ denotes Procrustes RMSE, subscripts $C$ and $G$ denote classical MDS and GMDS, and $\\Delta\\sigma$ is the relative stress reduction $100(\\sigma_C-\\sigma_G)/\\sigma_C$.",
    "",
    "\\section{Phase 1 Case Grid}",
    sprintf("Phase~1 contains %d cases across four groups: flat regular meshes, curved regular meshes, flat perforated meshes, and curved perforated meshes.", nrow(summary.df)),
    "",
    "\\begin{table}[H]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\caption{Phase 1 group summary.}",
    "\\begin{tabular}{@{}lrrrrrr@{}}",
    "\\toprule",
    "Group & $m$ & $\\overline{\\Delta\\sigma}$ (\\%) & $\\bar{\\sigma}_C$ & $\\bar{\\sigma}_G$ & $\\bar{\\rho}_C$ & $\\bar{\\rho}_G$ \\\\",
    "\\midrule",
    write_group_table(group.df),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    "\\begingroup",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\begin{longtable}{@{}p{0.26\\textwidth}rrrrrrr@{}}",
    "\\caption{Phase 1 case-level summary. Lower $\\sigma$ and $\\rho$ are better; $\\Delta\\sigma$ is computed relative to classical MDS under the GMDS path objective.}\\\\",
    "\\toprule",
    "Case & $n$ & $d$ & $\\sigma_C$ & $\\sigma_G$ & $\\Delta\\sigma$ (\\%) & $\\rho_C$ & $\\rho_G$ \\\\",
    "\\midrule",
    "\\endfirsthead",
    "\\toprule",
    "Case & $n$ & $d$ & $\\sigma_C$ & $\\sigma_G$ & $\\Delta\\sigma$ (\\%) & $\\rho_C$ & $\\rho_G$ \\\\",
    "\\midrule",
    "\\endhead",
    write_case_longtable(summary.df),
    "\\bottomrule",
    "\\end{longtable}",
    "\\endgroup",
    "",
    "\\begin{figure}[H]",
    "\\centering",
    sprintf("\\includegraphics[width=0.92\\textwidth]{%s}", tex_escape(file.path(basename(pdf_dir), "phase1-gmds-improvement-barplot.png"))),
    "\\caption{Relative stress reduction $\\Delta\\sigma$ over classical MDS under the GMDS path objective.}",
    "\\end{figure}",
    "",
    "\\begin{figure}[H]",
    "\\centering",
    sprintf("\\includegraphics[width=0.8\\textwidth]{%s}", tex_escape(file.path(basename(pdf_dir), "phase1-gmds-stress-scatter.png"))),
    "\\caption{Case-level $\\sigma_C$ versus $\\sigma_G$. Points below the diagonal favor GMDS.}",
    "\\end{figure}",
    "",
    "\\section{Representative Baseline Cases}",
    rep_sections,
    "",
    "\\section{Diagonal Connectivity Diagnostic}",
    "The orthogonal grid families use a 4-neighbor adjacency, so their graph-path metric is closer to an $L_1$ geometry than to the smooth parameter-domain truth. As an early diagnostic we repeated four representative cases with diagonal connectivity, which adds both diagonals of every available square and therefore makes the local graph metric more isotropic.",
    "This change is retained in the report only as a diagnosis of the artifact. It is not the adopted Phase~1 fix, because it changes the benchmark family itself. The post-review refinements in the next section keep the original rectangular 4-neighbor graphs unchanged.",
    "",
    "\\begin{table}[H]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\caption{Effect of diagonal connectivity on the optimized GMDS embeddings. Positive $\\Delta\\rho$ means the diagonal graph moves closer to truth.}",
    "\\begin{tabular}{@{}lrrrrrr@{}}",
    "\\toprule",
    "Case & $n$ & $\\sigma_G^{o}$ & $\\sigma_G^{d}$ & $\\rho_G^{o}$ & $\\rho_G^{d}$ & $\\Delta\\rho$ (\\%) \\\\",
    "\\midrule",
    write_connectivity_table(connectivity.df),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    connectivity_sections,
    "",
    "\\section{Post-Review GMDS Refinements}",
    "\\subsection{Reference-Panel Terminology}",
    "The left panel in the 2D mesh figures is an externally supplied benchmark reference, namely the synthetic parameter-plane coordinates $(u,v)$. It is not the unique intrinsic Euclidean geometry of the graph. For that reason the panel title was changed from \\texttt{Truth: parameter} to \\texttt{Reference geometry (u,v)}.",
    "",
    "\\subsection{Single-Path Diagnostic: More Iterations Are Not Enough}",
    sprintf(
      "In the flat $10 \\times 10$ orthogonal mesh, the original representative single-path GMDS run (24 iterations) reduced $\\sigma$ from %s to %s, but increased $\\rho$ from %s to %s. Extending the same single-path objective to 200 iterations lowered $\\sigma$ further to %s while only moving $\\rho$ to %s.",
      fmt_num(flat.base$gmds_stress_cmdscale[[1L]]),
      fmt_num(flat.base$gmds_stress_gmds[[1L]]),
      fmt_num(flat.base$procrustes_rmse_cmdscale[[1L]]),
      fmt_num(flat.base$procrustes_rmse_gmds[[1L]]),
      post_num("mesh_flat_10x10_2d", "gmds_single_200", "gmds_stress"),
      post_num("mesh_flat_10x10_2d", "gmds_single_200", "procrustes_rmse")
    ),
    "So the visually objectionable fold-over is not primarily an iteration-budget failure. The optimizer does keep lowering the single-path objective, but that fixed deterministic path choice can still reward asymmetric fold-overs on tied Manhattan shortest-path families.",
    "\\begin{figure}[H]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "mesh_flat_10x10_2d-triptych-gmds-200iter.png"))
    ),
    "\\caption{Flat orthogonal mesh after 200 iterations of single-path GMDS. The single-path stress keeps falling, but the final shape still drifts away from the reference square.}",
    "\\end{figure}",
    "",
    "\\subsection{Tie-Averaged Shortest Paths}",
    "The adopted rectangular-mesh fix is not to alter the graph family, but to change the path objective on tied shortest-path families. If $P_{ij}$ is the tied shortest-path family between vertices $i$ and $j$, we replace one deterministic chosen path by the exact uniform average",
    "\\[",
    "\\bar{h}_{ij}(Z) = \\frac{1}{|P_{ij}|} \\sum_{p \\in P_{ij}} h_p(Z) = \\sum_{e \\in E} \\pi_{ij}(e)\\,\\ell_e(Z).",
    "\\]",
    "Here $\\ell_e(Z)$ is the embedded Euclidean length of edge $e$, and $\\pi_{ij}(e)$ is the probability that $e$ belongs to a uniformly random shortest path from $i$ to $j$. The implementation computes these edge coefficients exactly from the shortest-path DAG, so no explicit path enumeration is required.",
    sprintf(
      "On the flat $10 \\times 10$ orthogonal mesh, single-path GMDS at 50 iterations reached $\\sigma=%s$ and $\\rho=%s$, whereas tie-averaged GMDS reached $\\sigma=%s$ and $\\rho=%s$ on the same 4-neighbor graph. On the ripple $12 \\times 12$ orthogonal mesh, single-path GMDS gave $\\sigma=%s$ and $\\rho=%s$, while tie-averaged GMDS gave $\\sigma=%s$ and $\\rho=%s$.",
      post_num("mesh_flat_10x10_2d", "gmds_single_50", "gmds_stress"),
      post_num("mesh_flat_10x10_2d", "gmds_single_50", "procrustes_rmse"),
      post_num("mesh_flat_10x10_2d", "gmds_average_50", "gmds_stress"),
      post_num("mesh_flat_10x10_2d", "gmds_average_50", "procrustes_rmse"),
      post_num("mesh_ripple_12x12_2d", "gmds_single_50", "gmds_stress"),
      post_num("mesh_ripple_12x12_2d", "gmds_single_50", "procrustes_rmse"),
      post_num("mesh_ripple_12x12_2d", "gmds_average_50", "gmds_stress"),
      post_num("mesh_ripple_12x12_2d", "gmds_average_50", "procrustes_rmse")
    ),
    "The flat case is the clearest indication that the earlier artifact came largely from deterministic single-path collapse rather than from an unavoidable limitation of orthogonal rectangular meshes.",
    "",
    "\\subsection{MDS Tether and Continuation}",
    "A second modification adds a quadratic tether to the classical-MDS solution $Z_{\\mathrm{cmd}}$:",
    "\\[",
    "E_t(Z) = E_{\\mathrm{avg}}(Z) + \\lambda_t \\lVert Z - Z_{\\mathrm{cmd}} \\rVert_F^2.",
    "\\]",
    "The continuation version lets the tether relax over the run. In the linear schedule used here,",
    "\\[",
    "\\lambda_t = \\left(1 - \\frac{t}{T}\\right)\\lambda_0 + \\frac{t}{T}\\lambda_T, \\qquad t = 0,1,\\dots,T.",
    "\\]",
    sprintf(
      "For the targeted orthogonal-mesh diagnostics below we used $(\\lambda_0, \\lambda_T) = (0.05, 0)$ and $T=50$. On the flat $10 \\times 10$ mesh, tie-averaged GMDS already resolves the pathology, so the tether changes $\\rho$ only from %s to %s. On the ripple $12 \\times 12$ mesh, the tether changes $\\sigma$ from %s to %s and $\\rho$ from %s to %s.",
      post_num("mesh_flat_10x10_2d", "gmds_average_50", "procrustes_rmse"),
      post_num("mesh_flat_10x10_2d", "gmds_average_lin_50", "procrustes_rmse"),
      post_num("mesh_ripple_12x12_2d", "gmds_average_50", "gmds_stress"),
      post_num("mesh_ripple_12x12_2d", "gmds_average_lin_50", "gmds_stress"),
      post_num("mesh_ripple_12x12_2d", "gmds_average_50", "procrustes_rmse"),
      post_num("mesh_ripple_12x12_2d", "gmds_average_lin_50", "procrustes_rmse")
    ),
    "So the tether should be interpreted as a shape-regularization control, not as a universal guarantee of lower Procrustes error. Its main role is to preserve the globally symmetric classical-MDS shape early in the optimization and then relax that bias gradually.",
    "The updated flat-mesh figure below corresponds to the same methods reported in Table~4, so the visual evidence now matches the post-review objective changes rather than the older baseline single-path figure.",
    "\\begin{figure}[H]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "mesh_flat_10x10_2d-postreview-quadtych.png"))
    ),
    "\\caption{Flat orthogonal mesh under the post-review methods. Panels show the reference geometry, classical MDS scored under the tie-averaged objective, tie-averaged GMDS, and tie-averaged GMDS with an MDS tether and linear continuation.}",
    "\\end{figure}",
    "",
    "\\begin{table}[H]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\caption{Targeted post-review diagnostics on orthogonal rectangular meshes. These rows keep the original 4-neighbor mesh graphs and vary only the path aggregation and hybrid objective. Lower $\\sigma$ and $\\rho$ are better.}",
    "\\begin{tabular}{@{}p{0.18\\textwidth}p{0.30\\textwidth}rrr@{}}",
    "\\toprule",
    "Case & Setting & $t$ (s) & $\\sigma$ & $\\rho$ \\\\",
    "\\midrule",
    write_post_review_table(postreview.df),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    "\\subsection{Ripple Surfaces Need 3D Embeddings}",
    "The 2D ripple example above is useful only for illustrating the limitation of planar embeddings. Its graph edge weights come from a genuinely curved 3D ripple surface, so there is no reason to expect the induced graph-geodesic metric to admit a faithful 2D Euclidean realization. The parameter plane is therefore only a sampling-domain reference, not an achievable intrinsic target.",
    "For that reason the ripple-family fidelity comparison should be read in 3D. We reran orthogonal ripple meshes in 3D for $12\\times12$, $15\\times15$, and $20\\times20$ grids, aligned the results to the reference surface, and compared classical MDS, tie-averaged GMDS, and tie-averaged GMDS with the same linear MDS tether.",
    sprintf(
      "At $12\\times12$, tie-averaged GMDS reduced $\\sigma$ from %s to %s in %ss while keeping $\\rho$ near %s; the tethered run took %ss and ended at $\\sigma=%s$ and $\\rho=%s$. The same pattern persists on denser meshes: at $20\\times20$, cmdscale took %ss for $\\sigma=%s$ and $\\rho=%s$, tie-averaged GMDS took %ss for $\\sigma=%s$ and $\\rho=%s$, and the tethered run took %ss for $\\sigma=%s$ and $\\rho=%s$.",
      surface_num(ripple3d.df, "mesh_ripple_12x12_3d_density", "cmdscale_average", "gmds_stress"),
      surface_num(ripple3d.df, "mesh_ripple_12x12_3d_density", "gmds_average_50", "gmds_stress"),
      surface_time(ripple3d.df, "mesh_ripple_12x12_3d_density", "gmds_average_50"),
      surface_num(ripple3d.df, "mesh_ripple_12x12_3d_density", "gmds_average_50", "procrustes_rmse"),
      surface_time(ripple3d.df, "mesh_ripple_12x12_3d_density", "gmds_average_lin_50"),
      surface_num(ripple3d.df, "mesh_ripple_12x12_3d_density", "gmds_average_lin_50", "gmds_stress"),
      surface_num(ripple3d.df, "mesh_ripple_12x12_3d_density", "gmds_average_lin_50", "procrustes_rmse"),
      surface_time(ripple3d.df, "mesh_ripple_20x20_3d_density", "cmdscale_average"),
      surface_num(ripple3d.df, "mesh_ripple_20x20_3d_density", "cmdscale_average", "gmds_stress"),
      surface_num(ripple3d.df, "mesh_ripple_20x20_3d_density", "cmdscale_average", "procrustes_rmse"),
      surface_time(ripple3d.df, "mesh_ripple_20x20_3d_density", "gmds_average_50"),
      surface_num(ripple3d.df, "mesh_ripple_20x20_3d_density", "gmds_average_50", "gmds_stress"),
      surface_num(ripple3d.df, "mesh_ripple_20x20_3d_density", "gmds_average_50", "procrustes_rmse"),
      surface_time(ripple3d.df, "mesh_ripple_20x20_3d_density", "gmds_average_lin_50"),
      surface_num(ripple3d.df, "mesh_ripple_20x20_3d_density", "gmds_average_lin_50", "gmds_stress"),
      surface_num(ripple3d.df, "mesh_ripple_20x20_3d_density", "gmds_average_lin_50", "procrustes_rmse")
    ),
    "\\begin{figure}[p]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "ripple_mesh_3d_density_postreview_grid.png"))
    ),
    "\\caption{Ripple meshes in 3D at increasing density, shown as aligned orthographic projections. Each row corresponds to one mesh size; columns show the reference surface, classical MDS scored under the tie-averaged objective, tie-averaged GMDS, and tie-averaged GMDS with an MDS tether and linear continuation.}",
    "\\end{figure}",
    "\\begin{table}[H]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\caption{3D ripple-surface diagnostics across mesh density. Lower $\\sigma$ and $\\rho$ are better.}",
    "\\begin{tabular}{@{}p{0.20\\textwidth}p{0.28\\textwidth}rrr@{}}",
    "\\toprule",
    "Case & Setting & $t$ (s) & $\\sigma$ & $\\rho$ \\\\",
    "\\midrule",
    write_post_review_table(ripple3d.df),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    "\\subsection{Paraboloid Surfaces Across Mesh Density}",
    "Paraboloid meshes give a second intrinsically curved family with a smoother, single-bowl geometry. As with ripple surfaces, the faithful comparison should be made in 3D rather than against the 2D sampling domain. We therefore repeated the same tie-averaged and tethered GMDS study on orthogonal paraboloid meshes at $12\\times12$, $15\\times15$, and $20\\times20$.",
    sprintf(
      "At $12\\times12$, tie-averaged GMDS reduced $\\sigma$ from %s to %s in %ss while moving $\\rho$ from %s to %s; the tethered run took %ss and finished at $\\sigma=%s$ and $\\rho=%s$. At $20\\times20$, cmdscale took %ss for $\\sigma=%s$ and $\\rho=%s$, tie-averaged GMDS took %ss for $\\sigma=%s$ and $\\rho=%s$, and the tethered run took %ss for $\\sigma=%s$ and $\\rho=%s$.",
      surface_num(paraboloid3d.df, "mesh_paraboloid_12x12_3d_density", "cmdscale_average", "gmds_stress"),
      surface_num(paraboloid3d.df, "mesh_paraboloid_12x12_3d_density", "gmds_average_50", "gmds_stress"),
      surface_time(paraboloid3d.df, "mesh_paraboloid_12x12_3d_density", "gmds_average_50"),
      surface_num(paraboloid3d.df, "mesh_paraboloid_12x12_3d_density", "cmdscale_average", "procrustes_rmse"),
      surface_num(paraboloid3d.df, "mesh_paraboloid_12x12_3d_density", "gmds_average_50", "procrustes_rmse"),
      surface_time(paraboloid3d.df, "mesh_paraboloid_12x12_3d_density", "gmds_average_lin_50"),
      surface_num(paraboloid3d.df, "mesh_paraboloid_12x12_3d_density", "gmds_average_lin_50", "gmds_stress"),
      surface_num(paraboloid3d.df, "mesh_paraboloid_12x12_3d_density", "gmds_average_lin_50", "procrustes_rmse"),
      surface_time(paraboloid3d.df, "mesh_paraboloid_20x20_3d_density", "cmdscale_average"),
      surface_num(paraboloid3d.df, "mesh_paraboloid_20x20_3d_density", "cmdscale_average", "gmds_stress"),
      surface_num(paraboloid3d.df, "mesh_paraboloid_20x20_3d_density", "cmdscale_average", "procrustes_rmse"),
      surface_time(paraboloid3d.df, "mesh_paraboloid_20x20_3d_density", "gmds_average_50"),
      surface_num(paraboloid3d.df, "mesh_paraboloid_20x20_3d_density", "gmds_average_50", "gmds_stress"),
      surface_num(paraboloid3d.df, "mesh_paraboloid_20x20_3d_density", "gmds_average_50", "procrustes_rmse"),
      surface_time(paraboloid3d.df, "mesh_paraboloid_20x20_3d_density", "gmds_average_lin_50"),
      surface_num(paraboloid3d.df, "mesh_paraboloid_20x20_3d_density", "gmds_average_lin_50", "gmds_stress"),
      surface_num(paraboloid3d.df, "mesh_paraboloid_20x20_3d_density", "gmds_average_lin_50", "procrustes_rmse")
    ),
    "\\begin{figure}[p]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "paraboloid_mesh_3d_density_postreview_grid.png"))
    ),
    "\\caption{Paraboloid meshes in 3D at increasing density, shown as aligned orthographic projections. Each row corresponds to one mesh size; columns show the reference surface, classical MDS scored under the tie-averaged objective, tie-averaged GMDS, and tie-averaged GMDS with an MDS tether and linear continuation.}",
    "\\end{figure}",
    "\\begin{table}[H]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\caption{3D paraboloid-surface diagnostics across mesh density. Lower $\\sigma$ and $\\rho$ are better.}",
    "\\begin{tabular}{@{}p{0.22\\textwidth}p{0.28\\textwidth}rrr@{}}",
    "\\toprule",
    "Case & Setting & $t$ (s) & $\\sigma$ & $\\rho$ \\\\",
    "\\midrule",
    write_post_review_table(paraboloid3d.df),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    "\\section{Performance Pass}",
    "This section addresses the runtime concern directly. The optimizer used throughout these GMDS experiments is compiled C++, not an R-level optimizer, but the original tie-averaged path still paid a heavy cost for R cache construction, list-based cache access, repeated inner-loop work, and a serial all-pairs evaluation.",
    "The performance pass therefore made four implementation changes: tie-average cache construction was moved into C++, a flattened contiguous cache was added alongside the legacy list cache, the flat optimizer path was wired into \\texttt{grip.optimize.geodesic.mds()}, and the flat all-pairs loop was parallelized with \\texttt{std::thread}.",
    "Even after these changes, GMDS should still remain much slower than \\texttt{cmdscale()}, because classical MDS solves one spectral problem while GMDS repeatedly evaluates and differentiates a nonconvex all-pairs path objective with Armijo backtracking.",
    "",
    "\\subsection{Tie-Average Cache Build: R vs C++}",
    sprintf(
      "On the dense orthogonal 3D mesh families used above, the new C++ tie-average cache builder is consistently faster than the legacy R builder. The median build speedup over these six cases is $%s\\times$. The slowest observed gain in this batch is $%s\\times$ on %s, while the largest gain is $%s\\times$ on %s.",
      fmt_ratio(cache.median.speedup),
      fmt_ratio(cache.worst$speedup_cpp_vs_r[[1L]]),
      tex_escape(cache.worst$case_short[[1L]]),
      fmt_ratio(cache.best$speedup_cpp_vs_r[[1L]]),
      tex_escape(cache.best$case_short[[1L]])
    ),
    "\\begin{figure}[H]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=0.88\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "phase1-cache-build-speedup.png"))
    ),
    "\\caption{Speedup of the C++ tie-average cache builder relative to the legacy R builder on the dense 3D mesh cases. Higher is better.}",
    "\\end{figure}",
    "\\begin{table}[H]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\caption{Tie-averaged shortest-path cache construction on dense orthogonal meshes. $P$ is the number of unordered vertex pairs and $|A|$ is the number of cached weighted edge terms in the flattened cache. $S=t_R/t_{C++}$ is the build speedup.}",
    "\\begin{tabular}{@{}p{0.24\\textwidth}rrrrrr@{}}",
    "\\toprule",
    "Case & $n$ & $P$ & $|A|$ & $t_R$ (s) & $t_{C++}$ (s) & $S$ \\\\",
    "\\midrule",
    write_cache_perf_table(cacheperf.df),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    "\\subsection{GMDS Optimizer: Before vs After}",
    "To isolate the optimizer change itself, we compare the current report timings against the saved pre-performance baseline CSVs from the same Phase~1 density runs. These comparisons therefore hold the benchmark families, iteration budgets, and reported methods fixed while changing only the implementation underneath.",
    sprintf(
      "Across the 3D ripple and paraboloid density runs, the median GMDS optimization speedup is $%s\\times$. The smallest speedup in this batch is $%s\\times$ on %s (%s), while the largest is $%s\\times$ on %s (%s).",
      fmt_ratio(opt.median.speedup),
      fmt_ratio(opt.worst$speedup_before_after[[1L]]),
      tex_escape(opt.worst$case_short[[1L]]),
      tex_escape(opt.worst$setting_short[[1L]]),
      fmt_ratio(opt.best$speedup_before_after[[1L]]),
      tex_escape(opt.best$case_short[[1L]]),
      tex_escape(opt.best$setting_short[[1L]])
    ),
    "\\begin{figure}[H]",
    "\\centering",
    sprintf(
      "\\includegraphics[width=0.92\\textwidth]{%s}",
      tex_escape(file.path(basename(pdf_dir), "phase1-optimizer-speedup.png"))
    ),
    "\\caption{Speedup of the optimized GMDS implementation relative to the saved pre-performance baseline timings. Higher is better.}",
    "\\end{figure}",
    "\\begin{table}[H]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\caption{GMDS optimizer timing before and after the flattened, threaded C++ path. $S=t_{old}/t_{new}$ is the speedup factor, and $q$ is the number of threads used by the new run.}",
    "\\begin{tabular}{@{}p{0.20\\textwidth}p{0.15\\textwidth}rrrr@{}}",
    "\\toprule",
    "Case & Setting & $t_{old}$ (s) & $t_{new}$ (s) & $S$ & $q$ \\\\",
    "\\midrule",
    write_optimizer_perf_table(optperf.df),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "The updated numbers reduce the implementation overhead substantially, but they do not erase the algorithmic gap to classical MDS. That remaining gap is expected: tie-averaged GMDS still solves a much harder optimization problem and still scales with repeated all-pairs path evaluations rather than a single eigendecomposition.",
    "",
    "\\clearpage",
    "",
    "\\section{Pending Phases}",
    "\\subsection{Phase 2: Wrapped Lattice Surfaces}",
    "Pending.",
    "\\subsection{Phase 3: Recursive Planar and Near-Planar Fractals}",
    "Pending.",
    "\\subsection{Phase 4: Triangulated Manifolds}",
    "Pending.",
    "\\subsection{Phase 5: Irregular Point-Sampled Manifolds}",
    "Pending.",
    "\\subsection{Phase 6: Recursive 3D Fractals and Porous Cubes}",
    "Pending.",
    "\\subsection{Phase 7: Irregular 3D Solids and Intrinsic Weighted Trees}",
    "Pending.",
    "",
    "\\end{document}"
  )

  writeLines(tex, tex_path)
}

cases <- build_phase1_cases()
results <- lapply(cases, run_case)
metrics <- do.call(rbind, lapply(results, `[[`, "metrics"))
metrics_path <- file.path(tmp_dir, "phase1-case-method-metrics.csv")
utils::write.csv(metrics, metrics_path, row.names = FALSE)

cmd_metrics <- metrics[metrics$method_id == "cmdscale", , drop = FALSE]
gmds_metrics <- metrics[metrics$method_id == "gmds", , drop = FALSE]
summary_df <- merge(
  cmd_metrics,
  gmds_metrics,
  by = c("case_id", "case_label", "case_group", "family", "n", "dim", "truth_label"),
  suffixes = c("_cmdscale", "_gmds"),
  sort = FALSE
)
summary_df$gmds_improvement_pct <- 100 * (
  summary_df$gmds_stress_cmdscale - summary_df$gmds_stress_gmds
) / pmax(summary_df$gmds_stress_cmdscale, 1e-12)

summary_df <- summary_df[order(summary_df$gmds_improvement_pct, decreasing = TRUE), , drop = FALSE]
summary_path <- file.path(tmp_dir, "phase1-case-summary.csv")
utils::write.csv(summary_df, summary_path, row.names = FALSE)

group_df <- aggregate(
  cbind(
    gmds_improvement_pct,
    gmds_stress_cmdscale,
    gmds_stress_gmds,
    procrustes_rmse_cmdscale,
    procrustes_rmse_gmds
  ) ~ case_group,
  data = summary_df,
  FUN = mean
)
group_df$n_cases <- as.integer(table(summary_df$case_group)[group_df$case_group])
names(group_df)[names(group_df) == "gmds_improvement_pct"] <- "mean_gmds_improvement_pct"
names(group_df)[names(group_df) == "gmds_stress_cmdscale"] <- "mean_gmds_stress_cmdscale"
names(group_df)[names(group_df) == "gmds_stress_gmds"] <- "mean_gmds_stress_gmds"
names(group_df)[names(group_df) == "procrustes_rmse_cmdscale"] <- "mean_procrustes_rmse_cmdscale"
names(group_df)[names(group_df) == "procrustes_rmse_gmds"] <- "mean_procrustes_rmse_gmds"
group_df <- group_df[, c(
  "case_group",
  "n_cases",
  "mean_gmds_improvement_pct",
  "mean_gmds_stress_cmdscale",
  "mean_gmds_stress_gmds",
  "mean_procrustes_rmse_cmdscale",
  "mean_procrustes_rmse_gmds"
)]
group_path <- file.path(tmp_dir, "phase1-group-summary.csv")
utils::write.csv(group_df, group_path, row.names = FALSE)

connectivity_cases <- build_connectivity_cases()
connectivity_results <- lapply(connectivity_cases, run_case)
connectivity_metrics <- do.call(rbind, lapply(connectivity_results, `[[`, "metrics"))
connectivity_metrics_path <- file.path(tmp_dir, "phase1-connectivity-case-method-metrics.csv")
utils::write.csv(connectivity_metrics, connectivity_metrics_path, row.names = FALSE)

connectivity_gmds <- connectivity_metrics[connectivity_metrics$method_id == "gmds", , drop = FALSE]
connectivity_orth <- connectivity_gmds[connectivity_gmds$connectivity == "orthogonal", , drop = FALSE]
connectivity_diag <- connectivity_gmds[connectivity_gmds$connectivity == "diagonal", , drop = FALSE]
connectivity_df <- merge(
  connectivity_orth,
  connectivity_diag,
  by = c("comparison_key", "case_label", "family", "n", "dim", "truth_label"),
  suffixes = c("_orth", "_diag"),
  sort = FALSE
)
connectivity_df$procrustes_improvement_pct <- 100 * (
  connectivity_df$procrustes_rmse_orth - connectivity_df$procrustes_rmse_diag
) / pmax(connectivity_df$procrustes_rmse_orth, 1e-12)
connectivity_df <- connectivity_df[order(connectivity_df$procrustes_improvement_pct, decreasing = TRUE), , drop = FALSE]
connectivity_path <- file.path(tmp_dir, "phase1-connectivity-summary.csv")
utils::write.csv(connectivity_df, connectivity_path, row.names = FALSE)

barplot_path <- file.path(pdf_dir, "phase1-gmds-improvement-barplot.png")
grDevices::png(barplot_path, width = 1800, height = 900, res = 150)
graphics::par(mar = c(10, 4, 3, 1))
graphics::barplot(
  summary_df$gmds_improvement_pct,
  names.arg = summary_df$case_id,
  las = 2L,
  col = "#3f7cac",
  border = NA,
  ylab = expression(Delta * sigma ~ "(%)"),
  main = expression("Phase 1: relative " * Delta * sigma)
)
graphics::abline(h = 0, lty = 2, col = "grey40")
invisible(grDevices::dev.off())

scatter_path <- file.path(pdf_dir, "phase1-gmds-stress-scatter.png")
grDevices::png(scatter_path, width = 1200, height = 900, res = 150)
graphics::plot(
  summary_df$gmds_stress_cmdscale,
  summary_df$gmds_stress_gmds,
  pch = 19,
  col = "#8f2d56",
  xlab = expression(sigma[C]),
  ylab = expression(sigma[G]),
  main = expression("Phase 1: " * sigma[C] * " vs " * sigma[G])
)
graphics::abline(0, 1, lty = 2, col = "grey40")
graphics::text(
  summary_df$gmds_stress_cmdscale,
  summary_df$gmds_stress_gmds,
  labels = summary_df$case_id,
  pos = 4,
  cex = 0.65
)
invisible(grDevices::dev.off())

representative_ids <- c(
  "mesh_flat_10x10_2d",
  "mesh_ripple_12x12_2d",
  "holes_flat_16x16_2d",
  "slit_flat_16x16_2d"
)
connectivity_ids <- c(
  "mesh_flat_10x10_2d",
  "mesh_ripple_12x12_2d",
  "holes_flat_16x16_2d",
  "slit_flat_16x16_2d"
)

result_map <- stats::setNames(results, vapply(results, function(x) x$case$id, character(1)))
postreview_case_ids <- c(
  "mesh_flat_10x10_2d",
  "mesh_ripple_12x12_2d"
)
postreview_results <- lapply(postreview_case_ids, function(case_id) {
  run_post_review_case(result_map[[case_id]]$case)
})
postreview_result_map <- stats::setNames(
  postreview_results,
  vapply(postreview_results, function(x) x$case$id, character(1))
)
postreview_df <- do.call(rbind, lapply(postreview_results, `[[`, "metrics"))
postreview_case_order <- match(postreview_df$case_id, postreview_case_ids)
postreview_method_order <- match(
  postreview_df$method_id,
  c(
    "cmdscale_single",
    "gmds_single_50",
    "gmds_single_200",
    "cmdscale_average",
    "gmds_average_50",
    "gmds_average_lin_50"
  )
)
postreview_df <- postreview_df[order(postreview_case_order, postreview_method_order), , drop = FALSE]
postreview_path <- file.path(tmp_dir, "phase1-postreview-diagnostics.csv")
utils::write.csv(postreview_df, postreview_path, row.names = FALSE)

ripple3d_cases <- build_ripple_3d_density_cases()
ripple3d_results <- lapply(ripple3d_cases, run_density_3d_case)
ripple3d_df <- do.call(rbind, lapply(ripple3d_results, `[[`, "metrics"))
ripple3d_case_order <- match(
  ripple3d_df$case_id,
  vapply(ripple3d_cases, function(x) x$id, character(1L))
)
ripple3d_method_order <- match(
  ripple3d_df$method_id,
  c("cmdscale_average", "gmds_average_50", "gmds_average_lin_50")
)
ripple3d_df <- ripple3d_df[order(ripple3d_case_order, ripple3d_method_order), , drop = FALSE]
ripple3d_path <- file.path(tmp_dir, "phase1-ripple-3d-density-diagnostics.csv")
utils::write.csv(ripple3d_df, ripple3d_path, row.names = FALSE)

paraboloid3d_cases <- build_paraboloid_3d_density_cases()
paraboloid3d_results <- lapply(paraboloid3d_cases, run_density_3d_case)
paraboloid3d_df <- do.call(rbind, lapply(paraboloid3d_results, `[[`, "metrics"))
paraboloid3d_case_order <- match(
  paraboloid3d_df$case_id,
  vapply(paraboloid3d_cases, function(x) x$id, character(1L))
)
paraboloid3d_method_order <- match(
  paraboloid3d_df$method_id,
  c("cmdscale_average", "gmds_average_50", "gmds_average_lin_50")
)
paraboloid3d_df <- paraboloid3d_df[
  order(paraboloid3d_case_order, paraboloid3d_method_order),
  ,
  drop = FALSE
]
paraboloid3d_path <- file.path(tmp_dir, "phase1-paraboloid-3d-density-diagnostics.csv")
utils::write.csv(paraboloid3d_df, paraboloid3d_path, row.names = FALSE)

cacheperf_cases <- c(ripple3d_cases, paraboloid3d_cases)
cacheperf_df <- do.call(rbind, lapply(cacheperf_cases, benchmark_cache_build_case))
cacheperf_case_order <- match(
  cacheperf_df$case_id,
  vapply(cacheperf_cases, function(x) x$id, character(1L))
)
cacheperf_df <- cacheperf_df[order(cacheperf_case_order), , drop = FALSE]
cacheperf_path <- file.path(tmp_dir, "phase1-cache-build-benchmark.csv")
utils::write.csv(cacheperf_df, cacheperf_path, row.names = FALSE)

optperf_df <- build_optimizer_perf_table(
  current.df = rbind(ripple3d_df, paraboloid3d_df),
  baseline.paths = c(
    file.path(tmp_dir, "phase1-ripple-3d-density-diagnostics-baseline-before-perf.csv"),
    file.path(tmp_dir, "phase1-paraboloid-3d-density-diagnostics-baseline-before-perf.csv")
  )
)
optperf_path <- file.path(tmp_dir, "phase1-optimizer-speedup-benchmark.csv")
utils::write.csv(optperf_df, optperf_path, row.names = FALSE)

cache_speedup_path <- file.path(pdf_dir, "phase1-cache-build-speedup.png")
grDevices::png(cache_speedup_path, width = 1800, height = 900, res = 150)
graphics::par(mar = c(11, 4, 3, 1))
graphics::barplot(
  cacheperf_df$speedup_cpp_vs_r,
  names.arg = cacheperf_df$case_short,
  las = 2L,
  col = "#588157",
  border = NA,
  ylab = "Speedup (x)",
  main = "Tie-average cache build: C++ vs R"
)
graphics::abline(h = 1, lty = 2, col = "grey40")
invisible(grDevices::dev.off())

opt_speedup_path <- file.path(pdf_dir, "phase1-optimizer-speedup.png")
opt_labels <- paste(optperf_df$case_short, optperf_df$setting_short, sep = "\n")
opt_cols <- ifelse(optperf_df$method_id == "gmds_average_50", "#6c8ebf", "#bc6c25")
grDevices::png(opt_speedup_path, width = 2200, height = 950, res = 150)
graphics::par(mar = c(12, 4, 3, 1))
graphics::barplot(
  optperf_df$speedup_before_after,
  names.arg = opt_labels,
  las = 2L,
  col = opt_cols,
  border = NA,
  ylab = "Speedup (x)",
  main = "GMDS optimizer speedup: before vs after performance pass"
)
graphics::abline(h = 1, lty = 2, col = "grey40")
graphics::legend(
  "topleft",
  legend = c("avg GMDS", "avg+tether"),
  fill = c("#6c8ebf", "#bc6c25"),
  bty = "n",
  cex = 0.9
)
invisible(grDevices::dev.off())

for (case_id in representative_ids) {
  case_result <- result_map[[case_id]]
  if (is.null(case_result)) {
    next
  }
  row <- summary_df[summary_df$case_id == case_id, , drop = FALSE]
  save_triptych(
    case_result = case_result,
    summary_row = row,
    path = file.path(pdf_dir, sprintf("%s-triptych.png", case_id))
  )
}

connectivity_result_map <- stats::setNames(
  connectivity_results,
  vapply(
    connectivity_results,
    function(x) paste(x$case$comparison_key, x$case$connectivity, sep = "::"),
    character(1)
  )
)
for (case_id in connectivity_ids) {
  orth_result <- connectivity_result_map[[paste(case_id, "orthogonal", sep = "::")]]
  diag_result <- connectivity_result_map[[paste(case_id, "diagonal", sep = "::")]]
  if (is.null(orth_result) || is.null(diag_result)) {
    next
  }
  row <- connectivity_df[connectivity_df$comparison_key == case_id, , drop = FALSE]
  save_connectivity_triptych(
    orth_result = orth_result,
    diag_result = diag_result,
    comparison_row = row,
    path = file.path(pdf_dir, sprintf("%s-connectivity-triptych.png", case_id))
  )
}

flat_post <- postreview_result_map[["mesh_flat_10x10_2d"]]
if (!is.null(flat_post) && !is.null(flat_post$layouts$gmds_single_200)) {
  cmd_row <- postreview_df[
    postreview_df$case_id == "mesh_flat_10x10_2d" &
      postreview_df$method_id == "cmdscale_single",
    ,
    drop = FALSE
  ]
  long_row <- postreview_df[
    postreview_df$case_id == "mesh_flat_10x10_2d" &
      postreview_df$method_id == "gmds_single_200",
    ,
    drop = FALSE
  ]
  save_named_triptych(
    case = flat_post$case,
    truth = flat_post$layouts$truth,
    middle_coords = flat_post$layouts$cmdscale,
    middle_title = sprintf(
      "cmdscale\nsigma %s, rho %s\nt %ss",
      fmt_num(cmd_row$gmds_stress[[1L]]),
      fmt_num(cmd_row$procrustes_rmse[[1L]]),
      fmt_time(cmd_row$elapsed_sec[[1L]])
    ),
    right_coords = flat_post$layouts$gmds_single_200,
    right_title = sprintf(
      "single GMDS (200)\nsigma %s, rho %s\nt %ss",
      fmt_num(long_row$gmds_stress[[1L]]),
      fmt_num(long_row$procrustes_rmse[[1L]]),
      fmt_time(long_row$elapsed_sec[[1L]])
    ),
    path = file.path(pdf_dir, "mesh_flat_10x10_2d-triptych-gmds-200iter.png")
  )
}

for (case_id in postreview_case_ids) {
  post_result <- postreview_result_map[[case_id]]
  if (is.null(post_result)) {
    next
  }
  cmd_row <- postreview_df[
    postreview_df$case_id == case_id &
      postreview_df$method_id == "cmdscale_average",
    ,
    drop = FALSE
  ]
  avg_row <- postreview_df[
    postreview_df$case_id == case_id &
      postreview_df$method_id == "gmds_average_50",
    ,
    drop = FALSE
  ]
  tether_row <- postreview_df[
    postreview_df$case_id == case_id &
      postreview_df$method_id == "gmds_average_lin_50",
    ,
    drop = FALSE
  ]
  if (nrow(cmd_row) == 0L || nrow(avg_row) == 0L || nrow(tether_row) == 0L) {
    next
  }
  save_named_quadtych(
    case = post_result$case,
    truth = post_result$layouts$truth,
    coords_list = list(
      post_result$layouts$cmdscale,
      post_result$layouts$gmds_average,
      post_result$layouts$gmds_average_linear_tether
    ),
    title_list = list(
      sprintf(
        "cmdscale (avg)\nsigma %s, rho %s\nt %ss",
        fmt_num(cmd_row$gmds_stress[[1L]]),
        fmt_num(cmd_row$procrustes_rmse[[1L]]),
        fmt_time(cmd_row$elapsed_sec[[1L]])
      ),
      sprintf(
        "avg GMDS (50)\nsigma %s, rho %s\nt %ss",
        fmt_num(avg_row$gmds_stress[[1L]]),
        fmt_num(avg_row$procrustes_rmse[[1L]]),
        fmt_time(avg_row$elapsed_sec[[1L]])
      ),
      sprintf(
        "avg+tether (50)\nsigma %s, rho %s\nt %ss",
        fmt_num(tether_row$gmds_stress[[1L]]),
        fmt_num(tether_row$procrustes_rmse[[1L]]),
        fmt_time(tether_row$elapsed_sec[[1L]])
      )
    ),
    path = file.path(pdf_dir, sprintf("%s-postreview-quadtych.png", case_id))
  )
}

save_surface_3d_density_grid(
  case_results = ripple3d_results,
  metrics.df = ripple3d_df,
  path = file.path(pdf_dir, "ripple_mesh_3d_density_postreview_grid.png"),
  heading = "Ripple surface meshes in 3D, shown as aligned orthographic projections"
)

save_surface_3d_density_grid(
  case_results = paraboloid3d_results,
  metrics.df = paraboloid3d_df,
  path = file.path(pdf_dir, "paraboloid_mesh_3d_density_postreview_grid.png"),
  heading = "Paraboloid surface meshes in 3D, shown as aligned orthographic projections"
)

tex_path <- file.path(manual_root, "reports", "gmds_mds_comparison_report_2026-03-31.tex")
write_report(
  summary.df = summary_df,
  group.df = group_df,
  connectivity.df = connectivity_df,
  postreview.df = postreview_df,
  ripple3d.df = ripple3d_df,
  paraboloid3d.df = paraboloid3d_df,
  cacheperf.df = cacheperf_df,
  optperf.df = optperf_df,
  representative_ids = representative_ids,
  connectivity_ids = connectivity_ids,
  pdf_dir = pdf_dir,
  tex_path = tex_path
)

report_pdf_path <- sub("\\.tex$", ".pdf", tex_path)

cat(sprintf("Phase 1 GMDS/MDS benchmark complete.\nMetrics: %s\nSummary: %s\nConnectivity summary: %s\nPost-review diagnostics: %s\nRipple 3D diagnostics: %s\nParaboloid 3D diagnostics: %s\nCache-build benchmark: %s\nOptimizer speedup benchmark: %s\nReport (tex): %s\nReport (pdf target): %s\n",
            metrics_path,
            summary_path,
            connectivity_path,
            postreview_path,
            ripple3d_path,
            paraboloid3d_path,
            cacheperf_path,
            optperf_path,
            tex_path,
            report_pdf_path))
if (nzchar(Sys.which("latexmk"))) {
  cat(sprintf("Compile the PDF with: (cd %s && latexmk -pdf %s)\n",
              file.path(manual_root, "reports"),
              basename(tex_path)))
}
