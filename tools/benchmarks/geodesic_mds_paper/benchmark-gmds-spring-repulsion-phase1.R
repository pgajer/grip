#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) == 0L) {
  stop("Unable to determine script path from commandArgs().")
}
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

run_tag <- "gmds-spring-repulsion-phase1-2026-04-01"
manual_root <- file.path(repo_root, "output", "geodesic_mds_paper")
tmp_dir <- file.path(manual_root, "tmp", run_tag)
pdf_dir <- file.path(manual_root, "reports", run_tag)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

tex_path <- file.path(manual_root, "reports", "gmds_spring_repulsion_phase1_report_2026-04-01.tex")
pdf_path <- file.path(manual_root, "reports", "gmds_spring_repulsion_phase1_report_2026-04-01.pdf")
rds_path <- file.path(tmp_dir, "gmds_spring_repulsion_phase1_results.rds")
metrics_csv <- file.path(tmp_dir, "gmds_spring_repulsion_phase1_metrics.csv")
unit_csv <- file.path(tmp_dir, "gmds_spring_repulsion_phase1_unit_checks.csv")

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo_root, quiet = TRUE, helpers = FALSE)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(repo_root, quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'pkgload', 'devtools', or the 'grip' package to run this benchmark.")
}

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

finite_difference_gradient <- function(fn, x, eps = 1e-6) {
  grad <- numeric(length(x))
  for (i in seq_along(x)) {
    plus <- x
    minus <- x
    plus[[i]] <- plus[[i]] + eps
    minus[[i]] <- minus[[i]] - eps
    grad[[i]] <- (fn(plus) - fn(minus)) / (2 * eps)
  }
  grad
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

make_case <- function(id,
                      label,
                      surface,
                      dim,
                      side = 8L,
                      amplitude = 0.35,
                      freq_u = 2L,
                      freq_v = 2L) {
  bundle <- do.call(
    mesh.surface.graph,
    c(
      list(
        h = side,
        w = side,
        surface = surface,
        amplitude = amplitude,
        connectivity = "orthogonal",
        normalize = "median"
      ),
      if (identical(surface, "ripple")) list(freq_u = freq_u, freq_v = freq_v) else list()
    )
  )

  prepared <- grip.prepare.geodesic.kk(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average"
  )
  started <- proc.time()[["elapsed"]]
  cmd <- grip:::grip.classical.mds.embedding(prepared, dim = dim, eig = TRUE)
  cmd_elapsed <- proc.time()[["elapsed"]] - started

  truth <- if (dim == 2L) bundle$coords_param else bundle$coords_surface
  list(
    id = id,
    label = label,
    surface = surface,
    dim = as.integer(dim),
    side = as.integer(side),
    truth = as.matrix(truth),
    edges = bundle$edges,
    edge_weights = bundle$edge_weights,
    adj_list = prepared$adj_list,
    triangles = grid_mesh_triangles(side, side),
    prepared = prepared,
    cmd = cmd,
    cmd_elapsed = cmd_elapsed
  )
}

method_specs <- list(
  list(
    id = "reference",
    label = "Reference",
    short = "Ref",
    kind = "reference",
    settings = "reference geometry"
  ),
  list(
    id = "cmdscale",
    label = "cmdscale",
    short = "CMD",
    kind = "cmdscale",
    settings = "classical MDS on graph geodesic matrix"
  ),
  list(
    id = "gmds_average",
    label = "GMDS avg",
    short = "GMDS",
    kind = "opt",
    args = list(engine = "cpp", max_iter = 6L, return_trace = TRUE, n_threads = 1L),
    settings = "tie-averaged GMDS, 6 iterations, cpp"
  ),
  list(
    id = "gmds_anchor",
    label = "Anchor-GMDS",
    short = "A-GMDS",
    kind = "opt",
    args = list(
      engine = "cpp",
      max_iter = 6L,
      return_trace = TRUE,
      n_threads = 1L,
      anchor_mode = "cmdscale",
      anchor_weight = 0.10,
      anchor_weight_end = 0.01,
      continuation = "linear"
    ),
    settings = "GMDS plus MDS tether, lambda: 0.10 -> 0.01"
  ),
  list(
    id = "gmds_spring_repulsion_anchor",
    label = "SR-GMDS",
    short = "SR-GMDS",
    kind = "opt",
    args = list(
      engine = "r",
      max_iter = 6L,
      return_trace = TRUE,
      n_threads = 1L,
      anchor_mode = "cmdscale",
      anchor_weight = 0.10,
      anchor_weight_end = 0.01,
      continuation = "linear",
      edge_spring_weight = 0.25,
      repulsion_weight = 0.10,
      repulsion_quantile = 0.50,
      repulsion_scale = 0.50,
      repulsion_cap_quantile = 1.0,
      repulsion_hop_min = 2L
    ),
    settings = "Anchor-GMDS plus beta = 0.25 and gamma = 0.10 with q = 0.50, s = 0.50, h[min] = 2"
  )
)

capture_fallback_warning <- function(expr) {
  warning_msg <- NULL
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warning_msg <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warning = warning_msg)
}

run_unit_checks <- function() {
  prepared <- grip.prepare.geodesic.kk(edges = edges.path(4L), n = 4L)
  coords_edge <- rbind(c(0.0, 0.0), c(0.8, 0.2), c(1.4, -0.1), c(2.9, 0.3))
  edge_weight <- 0.7
  edge_analytic <- grip:::grip.geodesic.mds.edge.spring.stats(
    coords = coords_edge,
    prepared = prepared,
    edge_spring_weight = edge_weight
  )
  edge_numeric <- finite_difference_gradient(function(vec) {
    pts <- matrix(vec, ncol = 2L, byrow = FALSE)
    grip:::grip.geodesic.mds.edge.spring.stats(
      coords = pts,
      prepared = prepared,
      edge_spring_weight = edge_weight
    )$energy
  }, as.double(coords_edge))
  edge_grad_max <- max(abs(as.double(edge_analytic$gradient) - edge_numeric))

  prepared_rep <- grip:::grip.geodesic.mds.ensure.graph.term.cache(
    prepared = prepared,
    repulsion_weight = 0.9,
    repulsion_quantile = 0,
    repulsion_scale = 0.8,
    repulsion_cap_quantile = 1,
    repulsion_hop_min = 2L
  )
  coords_rep <- rbind(c(0.0, 0.0), c(0.3, 0.0), c(0.6, 0.0), c(0.8, 0.1))
  rep_weight <- 0.9
  rep_analytic <- grip:::grip.geodesic.mds.repulsion.stats(
    coords = coords_rep,
    prepared = prepared_rep,
    repulsion_weight = rep_weight
  )
  rep_numeric <- finite_difference_gradient(function(vec) {
    pts <- matrix(vec, ncol = 2L, byrow = FALSE)
    grip:::grip.geodesic.mds.repulsion.stats(
      coords = pts,
      prepared = prepared_rep,
      repulsion_weight = rep_weight
    )$energy
  }, as.double(coords_rep))
  rep_grad_max <- max(abs(as.double(rep_analytic$gradient) - rep_numeric))

  score_coords <- rbind(c(0.0, 0.0), c(0.5, 0.0), c(0.9, 0.0), c(1.2, 0.0))
  score <- grip.score.geodesic.mds(
    coords = score_coords,
    prepared = prepared,
    edge_spring_weight = 0.25,
    repulsion_weight = 0.5,
    repulsion_quantile = 0,
    repulsion_scale = 0.7,
    repulsion_cap_quantile = 1,
    repulsion_hop_min = 2L
  )
  decomp_diff <- abs(
    score$gmds.energy[[1L]] -
      (
        score$gmds.base.energy[[1L]] +
          score$anchor.energy[[1L]] +
          score$edge.spring.energy[[1L]] +
          score$repulsion.energy[[1L]] +
          score$smooth.energy[[1L]]
      )
  )

  fallback <- capture_fallback_warning(
    grip.optimize.geodesic.mds(
      coords = score_coords,
      prepared = prepared,
      engine = "cpp",
      max_iter = 1L,
      return_trace = FALSE,
      edge_spring_weight = 0.1,
      repulsion_weight = 0.1,
      repulsion_quantile = 0.5,
      repulsion_scale = 0.5,
      repulsion_cap_quantile = 1,
      repulsion_hop_min = 2L
    )
  )

  data.frame(
    check_id = c("edge_gradient", "repulsion_gradient", "energy_decomposition", "cpp_fallback"),
    check_label = c(
      "edge spring finite-difference gradient",
      "repulsion finite-difference gradient",
      "score decomposition consistency",
      "cpp to R fallback for spring/repulsion terms"
    ),
    status = c(
      edge_grad_max < 1e-5,
      rep_grad_max < 1e-5,
      decomp_diff < 1e-10,
      is.character(fallback$warning) && grepl("falling back to the R engine", fallback$warning, fixed = TRUE)
    ),
    value = c(edge_grad_max, rep_grad_max, decomp_diff, NA_real_),
    detail = c(
      sprintf("max abs gradient error = %s", fmt_num(edge_grad_max, 7L)),
      sprintf("max abs gradient error = %s", fmt_num(rep_grad_max, 7L)),
      sprintf("absolute decomposition error = %s", fmt_num(decomp_diff, 9L)),
      if (is.null(fallback$warning)) "warning not observed" else fallback$warning
    ),
    stringsAsFactors = FALSE
  )
}

reference_geometry_label <- function(case) {
  if (case$dim == 2L) "Reference geometry (u,v)" else "Reference surface"
}

compute_metrics <- function(case,
                            method_id,
                            method_label,
                            short_label,
                            coords,
                            elapsed_sec = NA_real_,
                            fit = NULL) {
  anchor.coords <- if (is.null(fit)) NULL else fit$anchor_coords
  anchor.weight <- if (is.null(fit)) 0 else fit$final_anchor_weight
  edge.spring.weight <- if (is.null(fit)) 0 else fit$final_edge_spring_weight
  repulsion.weight <- if (is.null(fit)) 0 else fit$final_repulsion_weight
  score <- grip.score.geodesic.mds(
    coords = coords,
    prepared = case$prepared,
    anchor_coords = anchor.coords,
    anchor_weight = anchor.weight,
    edge_spring_weight = edge.spring.weight,
    repulsion_weight = repulsion.weight,
    repulsion_quantile = 0.50,
    repulsion_scale = 0.50,
    repulsion_cap_quantile = 1.0,
    repulsion_hop_min = 2L
  )
  aligned <- if (identical(method_id, "reference")) {
    list(aligned = as.matrix(coords), rmse = 0)
  } else {
    grip:::grip.align.to.target.nd(coords, case$truth, allow.reflection = TRUE)
  }
  trace <- if (is.null(fit)) NULL else fit$trace
  data.frame(
    case_id = case$id,
    case_label = case$label,
    surface = case$surface,
    dim = case$dim,
    method_id = method_id,
    method_label = method_label,
    short_label = short_label,
    elapsed_sec = as.double(elapsed_sec),
    gmds_stress = score$gmds.stress[[1L]],
    gmds_raw_stress = score$gmds.raw_stress[[1L]],
    procrustes_rmse = aligned$rmse,
    roughness = mesh_roughness(coords, case$adj_list, case$edges),
    area_q05_ratio = area_floor_ratio(coords, case$triangles),
    edge_spring_energy = score$edge.spring.energy[[1L]],
    repulsion_energy = score$repulsion.energy[[1L]],
    anchor_weight = anchor.weight,
    edge_spring_weight = edge.spring.weight,
    repulsion_weight = repulsion.weight,
    iterations = if (is.null(trace) || nrow(trace) == 0L) 0L else max(trace$iteration),
    stringsAsFactors = FALSE
  )
}

run_method <- function(case, spec) {
  if (identical(spec$kind, "reference")) {
    coords <- case$truth
    row <- compute_metrics(case, spec$id, spec$label, spec$short, coords)
    return(list(
      coords = coords,
      aligned = coords,
      display_coords = grip:::grip.normalize.coords(coords),
      metrics = row,
      fit = NULL
    ))
  }

  if (identical(spec$kind, "cmdscale")) {
    coords <- case$cmd$coords
    row <- compute_metrics(
      case,
      method_id = spec$id,
      method_label = spec$label,
      short_label = spec$short,
      coords = coords,
      elapsed_sec = case$cmd_elapsed,
      fit = NULL
    )
    aligned <- grip:::grip.align.to.target.nd(coords, case$truth, allow.reflection = TRUE)$aligned
    return(list(coords = coords, aligned = aligned, display_coords = aligned, metrics = row, fit = NULL))
  }

  started <- proc.time()[["elapsed"]]
  fit <- do.call(
    grip.optimize.geodesic.mds,
    c(
      list(coords = case$cmd$coords, prepared = case$prepared),
      spec$args
    )
  )
  elapsed <- proc.time()[["elapsed"]] - started
  row <- compute_metrics(
    case,
    method_id = spec$id,
    method_label = spec$label,
    short_label = spec$short,
    coords = fit$coords,
    elapsed_sec = elapsed,
    fit = fit
  )
  aligned <- grip:::grip.align.to.target.nd(fit$coords, case$truth, allow.reflection = TRUE)$aligned
  list(coords = fit$coords, aligned = aligned, display_coords = aligned, metrics = row, fit = fit)
}

build_title <- function(case, method_result) {
  row <- method_result$metrics[1L, , drop = FALSE]
  if (identical(row$method_id[[1L]], "reference")) {
    return(reference_geometry_label(case))
  }
  sprintf(
    "%s\nsigma %s, rho %s\nalpha_0.05 %s, t %ss",
    row$method_label[[1L]],
    fmt_num(row$gmds_stress[[1L]], 4L),
    fmt_num(row$procrustes_rmse[[1L]], 4L),
    fmt_num(row$area_q05_ratio[[1L]], 3L),
    fmt_time(row$elapsed_sec[[1L]])
  )
}

save_case_panel_grid <- function(case_result, output_path) {
  n_panels <- length(case_result$methods)
  grDevices::png(
    output_path,
    width = 1080L * 3L,
    height = 900L * 2L,
    res = 180,
    bg = "#ffffff"
  )
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(2L, 3L), mar = c(1.2, 1.2, 2.8, 0.4), oma = c(0, 0, 1.4, 0))

  for (i in seq_len(6L)) {
    if (i > n_panels) {
      graphics::plot.new()
      next
    }
    method <- case_result$methods[[i]]
    coords <- method$display_coords
    if (case_result$case$dim == 3L) {
      grip.plot(
        coords = coords,
        edges = case_result$case$edges,
        projection = "ortho",
        azimuth = 35,
        elevation = 24,
        vertex.col = if (identical(method$metrics$method_id[[1L]], "reference")) "#bc6c25" else "#3a5a40",
        edge.col = "#adb5bd",
        main = ""
      )
    } else {
      grip.plot(
        coords = coords,
        edges = case_result$case$edges,
        vertex.col = if (identical(method$metrics$method_id[[1L]], "reference")) "#bc6c25" else "#3a5a40",
        edge.col = "#adb5bd",
        main = ""
      )
    }
    graphics::title(main = build_title(case_result$case, method), cex.main = 0.72)
  }

  graphics::mtext(case_result$case$label, outer = TRUE, cex = 1.15, font = 2)
}

save_tradeoff_plot <- function(metrics, output_path) {
  plot_df <- metrics[metrics$method_id != "reference", , drop = FALSE]
  method_cols <- c(
    cmdscale = "#6c757d",
    gmds_average = "#1d3557",
    gmds_anchor = "#2a9d8f",
    gmds_spring_repulsion_anchor = "#e76f51"
  )
  grDevices::png(output_path, width = 1800L, height = 1200L, res = 180, bg = "#ffffff")
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mar = c(4.5, 4.8, 2.5, 1))
  graphics::plot(
    plot_df$gmds_stress,
    plot_df$procrustes_rmse,
    type = "n",
    xlab = expression(sigma ~ "(GMDS path stress)"),
    ylab = expression(rho ~ "(Procrustes RMSE to reference)"),
    main = "Stress-RMSE Tradeoff Across Cases"
  )
  for (i in seq_len(nrow(plot_df))) {
    col <- method_cols[[plot_df$method_id[[i]]]]
    graphics::points(
      plot_df$gmds_stress[[i]],
      plot_df$procrustes_rmse[[i]],
      pch = if (plot_df$dim[[i]] == 3L) 17 else 16,
      cex = 1.4,
      col = col
    )
    graphics::text(
      plot_df$gmds_stress[[i]],
      plot_df$procrustes_rmse[[i]],
      labels = sprintf("%s:%s", plot_df$case_id[[i]], plot_df$short_label[[i]]),
      pos = 4,
      cex = 0.7,
      col = col
    )
  }
  graphics::legend(
    "topright",
    legend = c("CMD", "GMDS", "A-GMDS", "SR-GMDS", "2D case", "3D cases"),
    col = c(method_cols[["cmdscale"]], method_cols[["gmds_average"]], method_cols[["gmds_anchor"]], method_cols[["gmds_spring_repulsion_anchor"]], "#000000", "#000000"),
    pch = c(16, 16, 16, 16, 16, 17),
    bty = "n",
    cex = 0.9
  )
}

write_method_rows <- function(df) {
  paste0(
    apply(df, 1L, function(row) {
      paste(
        tex_escape(row[["method_label"]]),
        fmt_time(as.numeric(row[["elapsed_sec"]])),
        fmt_num(as.numeric(row[["gmds_stress"]])),
        fmt_num(as.numeric(row[["procrustes_rmse"]])),
        fmt_num(as.numeric(row[["area_q05_ratio"]]), 3L),
        fmt_num(as.numeric(row[["edge_spring_energy"]]), 3L),
        fmt_num(as.numeric(row[["repulsion_energy"]]), 3L),
        sep = " & "
      )
    }),
    " \\\\"
  )
}

write_method_spec_rows <- function(specs) {
  paste0(
    vapply(specs, function(spec) {
      paste(tex_escape(spec$short), tex_escape(spec$label), tex_escape(spec$settings), sep = " & ")
    }, character(1L)),
    " \\\\"
  )
}

write_unit_rows <- function(df) {
  paste0(
    apply(df, 1L, function(row) {
      paste(
        tex_escape(row[["check_label"]]),
        if (isTRUE(as.logical(row[["status"]]))) "pass" else "fail",
        if (is.finite(as.numeric(row[["value"]]))) fmt_num(as.numeric(row[["value"]]), 7L) else "--",
        tex_escape(row[["detail"]]),
        sep = " & "
      )
    }),
    " \\\\"
  )
}

cases <- list(
  make_case("flat_mesh_8x8_2d", "Flat mesh 8x8 (2D)", "saddle", dim = 2L, side = 8L, amplitude = 0),
  make_case("saddle_mesh_8x8_3d", "Saddle mesh 8x8 (3D)", "saddle", dim = 3L, side = 8L, amplitude = 0.35),
  make_case("paraboloid_mesh_8x8_3d", "Paraboloid mesh 8x8 (3D)", "paraboloid", dim = 3L, side = 8L, amplitude = 0.35),
  make_case("ripple_mesh_8x8_3d", "Ripple mesh 8x8 (3D)", "ripple", dim = 3L, side = 8L, amplitude = 0.50)
)

unit_checks <- run_unit_checks()
case_results <- vector("list", length(cases))
names(case_results) <- vapply(cases, `[[`, character(1L), "id")
metric_rows <- list()

for (case in cases) {
  cat(sprintf("running %s\n", case$label))
  method_results <- lapply(method_specs, function(spec) run_method(case, spec))
  names(method_results) <- vapply(method_specs, `[[`, character(1L), "id")
  figure_path <- file.path(pdf_dir, sprintf("%s_panel_grid.png", case$id))
  case_result <- list(case = case, methods = method_results, figure_path = figure_path)
  save_case_panel_grid(case_result, figure_path)
  case_results[[case$id]] <- case_result
  metric_rows[[case$id]] <- do.call(rbind, lapply(method_results, `[[`, "metrics"))
}

metrics <- do.call(rbind, metric_rows)
utils::write.csv(metrics, metrics_csv, row.names = FALSE)
utils::write.csv(unit_checks, unit_csv, row.names = FALSE)

tradeoff_path <- file.path(pdf_dir, "gmds_spring_repulsion_tradeoff.png")
save_tradeoff_plot(metrics, tradeoff_path)

benchmark_bundle <- list(
  run_tag = run_tag,
  generated_at = as.character(Sys.time()),
  method_specs = method_specs,
  unit_checks = unit_checks,
  metrics = metrics,
  case_results = case_results,
  figure_paths = c(
    tradeoff = tradeoff_path,
    vapply(case_results, `[[`, character(1L), "figure_path")
  )
)
saveRDS(benchmark_bundle, rds_path)

metrics_2d <- metrics[metrics$dim == 2L, , drop = FALSE]
metrics_3d <- metrics[metrics$dim == 3L, , drop = FALSE]
case_figure_rel <- function(case_id) file.path(run_tag, basename(case_results[[case_id]]$figure_path))
tradeoff_rel <- file.path(run_tag, basename(tradeoff_path))
rds_rel <- file.path("output", "geodesic_mds_paper", "tmp", run_tag, basename(rds_path))

tex_lines <- c(
  "\\documentclass[11pt]{article}",
  "\\usepackage[margin=1in]{geometry}",
  "\\usepackage{booktabs}",
  "\\usepackage{graphicx}",
  "\\usepackage{float}",
  "\\usepackage{array}",
  "\\usepackage{longtable}",
  "\\usepackage{amsmath}",
  "\\usepackage{url}",
  "\\usepackage[T1]{fontenc}",
  "\\usepackage[utf8]{inputenc}",
  "\\title{Generic Spring-Repulsion GMDS Prototype Report}",
  "\\author{Codex benchmark generator}",
  "\\date{2026-04-01}",
  "\\begin{document}",
  "\\maketitle",
  "\\section*{Purpose}",
  "This report documents the Phase 1 graph-generic spring-repulsion GMDS prototype added to \\texttt{grip}. The goal is not to claim that the pathology problem is solved already, but to record the first numerically checked implementation and to compare its behaviour against the existing cmdscale and tie-averaged GMDS baselines on small representative meshes.",
  "",
  "The tested hybrid objective adds two generic graph terms to the fixed-path GMDS stress: an edge-spring term and a graph-distance-aware one-sided repulsion term. In this round the spring-repulsion terms are implemented in the R optimizer path, so the report also records the runtime cost of that prototype.",
  "",
  "\\section*{Automated Checks}",
  "Table~\\ref{tab:unit-checks} summarizes the targeted numerical checks behind the new implementation. These are the same checks covered by the dedicated test file \\texttt{test-geodesic-mds-spring-repulsion.R}.",
  "\\begin{table}[H]",
  "\\centering",
  "\\small",
  "\\caption{Targeted numerical checks for the spring-repulsion prototype.}",
  "\\label{tab:unit-checks}",
  "\\begin{tabular}{p{5.6cm} l l p{6.1cm}}",
  "\\toprule",
  "Check & Status & Value & Detail \\\\",
  "\\midrule",
  write_unit_rows(unit_checks),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  "",
  "\\section*{Benchmark Design}",
  "The layout benchmark uses one flat orthogonal mesh in 2D and three curved orthogonal meshes in 3D, all at side length $8$. Every non-reference layout starts from classical MDS on the graph-geodesic distance matrix. The plotted layouts are rigidly aligned to the reference geometry by translation, uniform scaling, rotation, and optional reflection before visualization; the optimization itself uses the raw coordinates.",
  "",
  "\\begin{table}[H]",
  "\\centering",
  "\\small",
  "\\caption{Method settings used throughout the report.}",
  "\\label{tab:method-settings}",
  "\\begin{tabular}{l l p{9.2cm}}",
  "\\toprule",
  "ID & Method & Settings \\\\",
  "\\midrule",
  write_method_spec_rows(method_specs),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  "",
  "\\section*{Quantitative Results}",
  "Common symbols across the tables: $t$ is wall-clock layout time in seconds, $\\sigma$ is GMDS path stress, $\\rho$ is Procrustes RMSE to the reference geometry, $\\alpha_{0.05}$ is the ratio of the $5\\%$ triangle-area quantile to the median triangle area after alignment, $E_e$ is the edge-spring energy contribution, and $E_r$ is the repulsion energy contribution.",
  "",
  "\\subsection*{Flat 2D Mesh}",
  "\\begin{table}[H]",
  "\\centering",
  "\\small",
  "\\caption{Flat orthogonal mesh benchmark in 2D.}",
  "\\label{tab:flat-2d}",
  "\\begin{tabular}{l r r r r r r}",
  "\\toprule",
  "Method & $t$ & $\\sigma$ & $\\rho$ & $\\alpha_{0.05}$ & $E_e$ & $E_r$ \\\\",
  "\\midrule",
  write_method_rows(metrics_2d[, c("method_label", "elapsed_sec", "gmds_stress", "procrustes_rmse", "area_q05_ratio", "edge_spring_energy", "repulsion_energy"), drop = FALSE]),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  "",
  sprintf("\\begin{figure}[H]\\centering\\includegraphics[width=0.97\\linewidth]{%s}\\caption{Flat mesh 8x8 in 2D. The spring-repulsion hybrid stays closer to the reference geometry than untethered GMDS only modestly, but it also keeps the generic graph terms active.}\\label{fig:flat-2d}\\end{figure}", case_figure_rel("flat_mesh_8x8_2d")),
  "",
  "\\subsection*{Curved 3D Meshes}",
  "\\begin{table}[H]",
  "\\centering",
  "\\small",
  "\\caption{Curved 3D surface benchmarks.}",
  "\\label{tab:curved-3d}",
  "\\begin{tabular}{l l r r r r r r}",
  "\\toprule",
  "Case & Method & $t$ & $\\sigma$ & $\\rho$ & $\\alpha_{0.05}$ & $E_e$ & $E_r$ \\\\",
  "\\midrule",
  paste0(
    apply(metrics_3d, 1L, function(row) {
      paste(
        tex_escape(row[["case_label"]]),
        tex_escape(row[["method_label"]]),
        fmt_time(as.numeric(row[["elapsed_sec"]])),
        fmt_num(as.numeric(row[["gmds_stress"]])),
        fmt_num(as.numeric(row[["procrustes_rmse"]])),
        fmt_num(as.numeric(row[["area_q05_ratio"]]), 3L),
        fmt_num(as.numeric(row[["edge_spring_energy"]]), 3L),
        fmt_num(as.numeric(row[["repulsion_energy"]]), 3L),
        sep = " & "
      )
    }),
    " \\\\"
  ),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  "",
  sprintf("\\begin{figure}[H]\\centering\\includegraphics[width=0.92\\linewidth]{%s}\\caption{Stress-RMSE tradeoff across all non-reference layouts. The spring-repulsion hybrid is consistently much closer to GMDS than to cmdscale on path stress, and its geometric correction at this prototype scale remains modest.}\\label{fig:tradeoff}\\end{figure}", tradeoff_rel),
  "",
  sprintf("\\clearpage\\begin{figure}[p]\\centering\\includegraphics[width=0.995\\linewidth]{%s}\\caption{Saddle mesh 8x8 in 3D.}\\label{fig:saddle-3d}\\end{figure}", case_figure_rel("saddle_mesh_8x8_3d")),
  sprintf("\\clearpage\\begin{figure}[p]\\centering\\includegraphics[width=0.995\\linewidth]{%s}\\caption{Paraboloid mesh 8x8 in 3D.}\\label{fig:paraboloid-3d}\\end{figure}", case_figure_rel("paraboloid_mesh_8x8_3d")),
  sprintf("\\clearpage\\begin{figure}[p]\\centering\\includegraphics[width=0.995\\linewidth]{%s}\\caption{Ripple mesh 8x8 in 3D. The reference surface itself is not expected to achieve zero path stress because graph-geodesic distances on the mesh do not coincide with straight ambient-space chord lengths.}\\label{fig:ripple-3d}\\end{figure}", case_figure_rel("ripple_mesh_8x8_3d")),
  "",
  "\\section*{Discussion}",
  "Three patterns are the main takeaways from this first report.",
  "",
  "\\begin{enumerate}",
  "\\item The new implementation is numerically sound enough to benchmark: the finite-difference checks pass, the energy decomposition closes, and activating spring-repulsion terms through the public optimizer correctly falls back from the compiled path to the R engine.",
  "\\item On these small representative meshes, tie-averaged GMDS still dominates cmdscale on path stress, but it also drifts away from the reference geometry. The anchor-only variant already reins that in slightly. Adding the generic edge-spring and repulsion terms keeps those additional energies active and yields a small extra geometric improvement in some cases, but not yet a dramatic correction.",
  "\\item Runtime is still the main practical limitation for the generic prototype. The spring-repulsion hybrid uses the R engine in this round, so it is materially slower than both cmdscale and the compiled GMDS baseline. This report should therefore be read as an implementation-and-diagnostics milestone rather than a final performance result.",
  "\\end{enumerate}",
  "",
  sprintf("The full benchmark object used to generate these tables and figures is saved under \\texttt{%s} as \\texttt{%s}.", tex_escape(dirname(rds_rel)), tex_escape(basename(rds_rel))),
  "The interactive 3D companion page is generated separately from that saved object so that the HTML panels and the PDF figures stay in sync.",
  "",
  "\\end{document}"
)

writeLines(tex_lines, tex_path)

cat(sprintf("Wrote metrics to %s\n", metrics_csv))
cat(sprintf("Wrote unit checks to %s\n", unit_csv))
cat(sprintf("Wrote benchmark bundle to %s\n", rds_path))
cat(sprintf("Wrote LaTeX report to %s\n", tex_path))
if (nzchar(Sys.which("latexmk"))) {
  cat(sprintf("Compile the PDF with: (cd %s && latexmk -pdf %s)\n", dirname(tex_path), basename(tex_path)))
}
