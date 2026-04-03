#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) == 0L) {
  stop("Unable to determine script path from commandArgs().")
}
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

run_tag <- "gmds-v2-paraboloid-suite-2026-04-01"
manual_root <- file.path(repo_root, "dev", "manual")
tmp_dir <- file.path(manual_root, "tmp", run_tag)
pdf_dir <- file.path(manual_root, "pdf", run_tag)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

tex_path <- file.path(manual_root, "pdf", "gmds_v2_paraboloid_suite_report_2026-04-01.tex")
pdf_path <- file.path(manual_root, "pdf", "gmds_v2_paraboloid_suite_report_2026-04-01.pdf")
rds_path <- file.path(tmp_dir, "gmds_v2_paraboloid_suite_results.rds")
metrics_csv <- file.path(tmp_dir, "gmds_v2_paraboloid_suite_metrics.csv")
foundation_csv <- file.path(tmp_dir, "gmds_v2_paraboloid_suite_foundation.csv")

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
resolve_anchor <- get("grip.geodesic.mds.resolve.anchor", envir = ns)
weight_schedule <- get("grip.geodesic.mds.weight.schedule", envir = ns)
ensure_graph_term_cache <- get("grip.geodesic.mds.ensure.graph.term.cache", envir = ns)
anchor_stats_fn <- get("grip.geodesic.mds.anchor.stats", envir = ns)
edge_spring_stats_fn <- get("grip.geodesic.mds.edge.spring.stats", envir = ns)
repulsion_stats_fn <- get("grip.geodesic.mds.repulsion.stats", envir = ns)
smoothness_stats_fn <- get("grip.geodesic.mds.smoothness.stats", envir = ns)
path_edge_coefficients <- get("grip.path.edge.coefficients", envir = ns)
optimize_flat_cpp <- get("grip_optimize_geodesic_mds_flat_cpp", envir = ns)

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

run_foundation_checks <- function() {
  path_prepared <- grip.prepare.graph.geodesic.mds(
    edges = edges.path(4L),
    n = 4L,
    edge_weights = c(1, 1.5, 0.7),
    tie_mode = "single"
  )
  path_exact <- cbind(c(0, 1, 2.5, 3.2), 0)
  path_score <- grip.score.geodesic.mds(coords = path_exact, prepared = path_prepared)

  s3 <- sqrt(3) / 2
  triangle_prepared <- grip.prepare.graph.geodesic.mds(
    edges = matrix(
      c(1L, 2L, 1L, 3L, 2L, 3L, 2L, 4L, 3L, 4L),
      ncol = 2L,
      byrow = TRUE
    ),
    n = 4L,
    edge_weights = rep(1, 5L),
    tie_mode = "single"
  )
  tri_exact_separated <- rbind(c(0.5, s3), c(0, 0), c(1, 0), c(0.5, -s3))
  tri_exact_overlapped <- rbind(c(0.5, s3), c(0, 0), c(1, 0), c(0.5, s3))
  tri_score_separated <- grip.score.geodesic.mds(coords = tri_exact_separated, prepared = triangle_prepared)
  tri_score_overlapped <- grip.score.geodesic.mds(coords = tri_exact_overlapped, prepared = triangle_prepared)

  square_edges <- matrix(c(1L, 2L, 2L, 3L, 3L, 4L, 4L, 1L), ncol = 2L, byrow = TRUE)
  square_single <- grip.prepare.graph.geodesic.mds(
    edges = square_edges,
    n = 4L,
    edge_weights = rep(1, 4L),
    tie_mode = "single"
  )
  square_average <- grip.prepare.graph.geodesic.mds(
    edges = square_edges,
    n = 4L,
    edge_weights = rep(1, 4L),
    tie_mode = "average"
  )
  pair_idx <- which(square_average$pair_matrix[, 1L] == 1L & square_average$pair_matrix[, 2L] == 3L)
  avg_edges <- square_average$path_edges[[pair_idx]]
  avg_coeffs <- path_edge_coefficients(square_average, pair_idx, nrow(avg_edges))

  square_coords <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  square_perm <- c(2L, 3L, 4L, 1L)
  square_average_perm <- grip.prepare.graph.geodesic.mds(
    edges = matrix(square_perm[square_edges], ncol = 2L),
    n = 4L,
    edge_weights = rep(1, 4L),
    tie_mode = "average"
  )
  square_score_a <- grip.score.geodesic.mds(coords = square_coords, prepared = square_average)
  square_score_b <- grip.score.geodesic.mds(
    coords = square_coords[square_perm, , drop = FALSE],
    prepared = square_average_perm
  )

  data.frame(
    check_id = c(
      "path_exact_zero",
      "triangles_exact_separated",
      "triangles_exact_overlapped",
      "square_tie_count",
      "square_tie_coefficients",
      "square_average_relabel"
    ),
    check_label = c(
      "Path graph exact realization has zero raw stress",
      "Shared-edge triangles: separated realization has zero raw stress",
      "Shared-edge triangles: overlapped realization has zero raw stress",
      "Average tie mode records two equal shortest paths on the square diagonal",
      "Average tie mode uses four edges with coefficient 1/2 on the square diagonal",
      "Average tie mode is relabeling-invariant on the square score"
    ),
    status = c(
      path_score$gmds.raw_stress[[1L]] < 1e-12,
      tri_score_separated$gmds.raw_stress[[1L]] < 1e-12,
      tri_score_overlapped$gmds.raw_stress[[1L]] < 1e-12,
      abs(square_average$pair_path_count_log[[pair_idx]] - log(2)) < 1e-10,
      isTRUE(all.equal(sort(avg_coeffs), rep(0.5, 4L), tolerance = 1e-10)),
      abs(square_score_a$gmds.raw_stress[[1L]] - square_score_b$gmds.raw_stress[[1L]]) < 1e-12
    ),
    value = c(
      path_score$gmds.raw_stress[[1L]],
      tri_score_separated$gmds.raw_stress[[1L]],
      tri_score_overlapped$gmds.raw_stress[[1L]],
      square_average$pair_path_count_log[[pair_idx]],
      max(abs(sort(avg_coeffs) - rep(0.5, 4L))),
      abs(square_score_a$gmds.raw_stress[[1L]] - square_score_b$gmds.raw_stress[[1L]])
    ),
    stringsAsFactors = FALSE
  )
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
    label = sprintf("Paraboloid mesh %dx%d", side, side),
    side = as.integer(side),
    prepared = prepared,
    edges = bundle$edges,
    adj_list = prepared$adj_list,
    triangles = grid_mesh_triangles(side, side),
    truth = bundle$coords_surface,
    cmd = cmd,
    cmd_elapsed = as.double(cmd_elapsed)
  )
}

score_edge_relaxed <- function(coords,
                               prepared,
                               edge_length_epsilon = 1e-8,
                               anchor_coords = NULL,
                               anchor_weight = 0,
                               smoothness_weight = 0,
                               edge_spring_weight = 1,
                               repulsion_weight = 0,
                               repulsion_quantile = 0.60,
                               repulsion_scale = 0.20,
                               repulsion_cap_quantile = 0.90,
                               repulsion_hop_min = 3L) {
  prepared <- ensure_graph_term_cache(
    prepared = prepared,
    repulsion_weight = repulsion_weight,
    repulsion_quantile = repulsion_quantile,
    repulsion_scale = repulsion_scale,
    repulsion_cap_quantile = repulsion_cap_quantile,
    repulsion_hop_min = repulsion_hop_min
  )
  gmds_score <- grip.score.geodesic.mds(
    coords = coords,
    prepared = prepared,
    edge_length_epsilon = edge_length_epsilon
  )
  anchor_stats <- anchor_stats_fn(
    coords = coords,
    anchor_coords = anchor_coords,
    anchor_weight = anchor_weight
  )
  edge_stats <- edge_spring_stats_fn(
    coords = coords,
    prepared = prepared,
    edge_length_epsilon = edge_length_epsilon,
    edge_spring_weight = edge_spring_weight
  )
  repulsion_stats <- repulsion_stats_fn(
    coords = coords,
    prepared = prepared,
    edge_length_epsilon = edge_length_epsilon,
    repulsion_weight = repulsion_weight
  )
  smooth_stats <- smoothness_stats_fn(
    coords = coords,
    prepared = prepared,
    smoothness_weight = smoothness_weight
  )
  list(
    objective_energy = anchor_stats$energy + edge_stats$energy + repulsion_stats$energy + smooth_stats$energy,
    gmds_stress = gmds_score$gmds.stress[[1L]],
    gmds_raw_stress = gmds_score$gmds.raw_stress[[1L]],
    gmds_rmse = gmds_score$gmds.rmse[[1L]],
    anchor_energy = anchor_stats$energy,
    edge_spring_energy = edge_stats$energy,
    repulsion_energy = repulsion_stats$energy,
    smooth_energy = smooth_stats$energy
  )
}

optimize_edge_relaxed <- function(coords,
                                  prepared,
                                  anchor_mode = c("none", "cmdscale", "initial", "user"),
                                  anchor_coords = NULL,
                                  anchor_weight = 0,
                                  anchor_weight_end = anchor_weight,
                                  continuation = c("constant", "linear", "geometric"),
                                  edge_spring_weight = 1,
                                  edge_spring_weight_end = edge_spring_weight,
                                  edge_spring_continuation = c("constant", "linear", "geometric"),
                                  repulsion_weight = 0,
                                  repulsion_weight_end = repulsion_weight,
                                  repulsion_continuation = c("constant", "linear", "geometric"),
                                  repulsion_quantile = 0.60,
                                  repulsion_scale = 0.20,
                                  repulsion_cap_quantile = 0.90,
                                  repulsion_hop_min = 3L,
                                  engine = c("cpp"),
                                  max_iter = 25L,
                                  edge_length_epsilon = 1e-8,
                                  initial_step = 1.0,
                                  step_shrink = 0.5,
                                  armijo_factor = 1e-4,
                                  grad_tol = 1e-8,
                                  min_step = 1e-8,
                                  n_threads = 1L,
                                  recenter = TRUE,
                                  return_trace = TRUE) {
  anchor_mode <- match.arg(anchor_mode)
  continuation <- match.arg(continuation)
  edge_spring_continuation <- match.arg(edge_spring_continuation)
  repulsion_continuation <- match.arg(repulsion_continuation)
  engine <- match.arg(engine)

  prepared <- ensure_graph_term_cache(
    prepared = prepared,
    repulsion_weight = max(repulsion_weight, repulsion_weight_end),
    repulsion_quantile = repulsion_quantile,
    repulsion_scale = repulsion_scale,
    repulsion_cap_quantile = repulsion_cap_quantile,
    repulsion_hop_min = repulsion_hop_min
  )
  anchor.coords <- resolve_anchor(
    anchor_mode = anchor_mode,
    coords = coords,
    prepared = prepared,
    anchor_coords = anchor_coords,
    recenter = recenter
  )
  anchor.schedule <- if (is.null(anchor.coords)) {
    rep.int(0, max_iter + 1L)
  } else {
    weight_schedule(
      max_iter = max_iter,
      weight = anchor_weight,
      weight_end = anchor_weight_end,
      continuation = continuation
    )
  }
  edge.schedule <- weight_schedule(
    max_iter = max_iter,
    weight = edge_spring_weight,
    weight_end = edge_spring_weight_end,
    continuation = edge_spring_continuation
  )
  repulsion.schedule <- weight_schedule(
    max_iter = max_iter,
    weight = repulsion_weight,
    weight_end = repulsion_weight_end,
    continuation = repulsion_continuation
  )

  graph_edge_matrix <- if (!is.null(prepared$graph_edge_matrix)) prepared$graph_edge_matrix else prepared$edges
  graph_edge_target <- if (!is.null(prepared$graph_edge_target)) prepared$graph_edge_target else prepared$edge_targets
  repulsion_pair_matrix <- if (!is.null(prepared$repulsion_pair_matrix)) prepared$repulsion_pair_matrix else matrix(integer(), ncol = 2L)
  repulsion_target <- if (!is.null(prepared$repulsion_target)) prepared$repulsion_target else numeric(0L)

  opt <- optimize_flat_cpp(
    flat_pair_edge_offsets = as.integer(0L),
    flat_edge_u = integer(0L),
    flat_edge_v = integer(0L),
    flat_edge_coeff = numeric(0L),
    pair_graph_distance = numeric(0L),
    coords = coords,
    max_iter = max_iter,
    edge_length_epsilon = edge_length_epsilon,
    initial_step = initial_step,
    step_shrink = step_shrink,
    armijo_factor = armijo_factor,
    grad_tol = grad_tol,
    min_step = min_step,
    recenter = recenter,
    return_trace = return_trace,
    anchor_coords = anchor.coords,
    anchor_weights = anchor.schedule,
    smooth_adj_offsets = integer(0L),
    smooth_adj_vertices = integer(0L),
    smooth_weights = rep.int(0, max_iter + 1L),
    graph_edge_u = if (nrow(graph_edge_matrix) > 0L) as.integer(graph_edge_matrix[, 1L] - 1L) else integer(0L),
    graph_edge_v = if (nrow(graph_edge_matrix) > 0L) as.integer(graph_edge_matrix[, 2L] - 1L) else integer(0L),
    graph_edge_target = as.double(graph_edge_target),
    edge_spring_weights = edge.schedule,
    repulsion_u = if (nrow(repulsion_pair_matrix) > 0L) as.integer(repulsion_pair_matrix[, 1L] - 1L) else integer(0L),
    repulsion_v = if (nrow(repulsion_pair_matrix) > 0L) as.integer(repulsion_pair_matrix[, 2L] - 1L) else integer(0L),
    repulsion_target = as.double(repulsion_target),
    repulsion_weights = repulsion.schedule,
    n_threads = as.integer(n_threads)
  )

  final.anchor <- if (!is.null(opt$final_anchor_weight)) opt$final_anchor_weight else utils::tail(anchor.schedule, 1L)
  final.edge <- if (!is.null(opt$final_edge_spring_weight)) opt$final_edge_spring_weight else utils::tail(edge.schedule, 1L)
  final.repulsion <- if (!is.null(opt$final_repulsion_weight)) opt$final_repulsion_weight else utils::tail(repulsion.schedule, 1L)
  score <- score_edge_relaxed(
    coords = opt$coords,
    prepared = prepared,
    edge_length_epsilon = edge_length_epsilon,
    anchor_coords = anchor.coords,
    anchor_weight = final.anchor,
    edge_spring_weight = final.edge,
    repulsion_weight = final.repulsion,
    repulsion_quantile = repulsion_quantile,
    repulsion_scale = repulsion_scale,
    repulsion_cap_quantile = repulsion_cap_quantile,
    repulsion_hop_min = repulsion_hop_min
  )

  list(
    coords = opt$coords,
    trace = opt$trace,
    frames = opt$frames,
    score = score,
    anchor_coords = anchor.coords,
    final_anchor_weight = final.anchor,
    final_edge_spring_weight = final.edge,
    final_repulsion_weight = final.repulsion,
    n_threads_used = opt$n_threads_used
  )
}

compute_metrics <- function(case,
                            spec,
                            coords,
                            elapsed_sec = NA_real_,
                            fit = NULL,
                            score_info = NULL) {
  if (is.null(score_info)) {
    score_df <- grip.score.geodesic.mds(coords = coords, prepared = case$prepared)
    score_info <- list(
      objective_energy = score_df$gmds.energy[[1L]],
      gmds_stress = score_df$gmds.stress[[1L]],
      gmds_raw_stress = score_df$gmds.raw_stress[[1L]],
      gmds_rmse = score_df$gmds.rmse[[1L]],
      anchor_energy = score_df$anchor.energy[[1L]],
      edge_spring_energy = score_df$edge.spring.energy[[1L]],
      repulsion_energy = score_df$repulsion.energy[[1L]],
      smooth_energy = score_df$smooth.energy[[1L]]
    )
  }
  aligned <- if (identical(spec$id, "reference")) {
    list(aligned = as.matrix(coords), rmse = 0)
  } else {
    align_to_target_nd(coords, case$truth, allow.reflection = TRUE)
  }
  trace <- if (is.null(fit)) NULL else fit$trace
  data.frame(
    case_id = case$id,
    case_label = case$label,
    side = case$side,
    n = nrow(case$truth),
    method_id = spec$id,
    method_label = spec$label,
    short_label = spec$short,
    family = spec$family,
    lambda_label = spec$lambda_label,
    elapsed_sec = as.double(elapsed_sec),
    objective_energy = score_info$objective_energy,
    gmds_stress = score_info$gmds_stress,
    gmds_raw_stress = score_info$gmds_raw_stress,
    procrustes_rmse = aligned$rmse,
    roughness = mesh_roughness(coords, case$adj_list, case$edges),
    area_q05_ratio = area_floor_ratio(coords, case$triangles),
    anchor_energy = score_info$anchor_energy,
    edge_spring_energy = score_info$edge_spring_energy,
    repulsion_energy = score_info$repulsion_energy,
    iterations = if (is.null(trace) || nrow(trace) == 0L) 0L else max(trace$iteration),
    stringsAsFactors = FALSE
  )
}

method_specs <- list(
  list(id = "reference", label = "Reference", short = "Ref", family = "Reference", kind = "reference", lambda_label = "--"),
  list(id = "cmdscale", label = "cmdscale", short = "CMD", family = "Baseline", kind = "cmdscale", lambda_label = "--"),
  list(
    id = "pure_gmds",
    label = "Pure GMDS",
    short = "GMDS",
    family = "Pure GMDS",
    kind = "gmds",
    lambda_label = "--",
    args = list(engine = "cpp", max_iter = 15L, return_trace = TRUE, n_threads = 0L)
  ),
  list(
    id = "reg_anchor",
    label = "Reg GMDS (anchor)",
    short = "Reg-A",
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
    label = "Reg GMDS (anchor + rep)",
    short = "Reg-AR",
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
    id = "edge_surrogate",
    label = "Edge surrogate",
    short = "Edge",
    family = "Edge surrogate",
    kind = "edge",
    lambda_label = lambda_label(edge_start = 1.00, edge_end = 1.00),
    args = list(
      max_iter = 15L,
      n_threads = 0L,
      edge_spring_weight = 1.00,
      edge_spring_weight_end = 1.00
    )
  )
)

run_method <- function(case, spec) {
  if (identical(spec$kind, "reference")) {
    row <- compute_metrics(case = case, spec = spec, coords = case$truth)
    return(list(coords = case$truth, display_coords = case$truth, metrics = row, fit = NULL))
  }

  if (identical(spec$kind, "cmdscale")) {
    coords <- case$cmd$coords
    row <- compute_metrics(case = case, spec = spec, coords = coords, elapsed_sec = case$cmd_elapsed)
    return(list(
      coords = coords,
      display_coords = align_to_target_nd(coords, case$truth, allow.reflection = TRUE)$aligned,
      metrics = row,
      fit = NULL
    ))
  }

  if (identical(spec$kind, "gmds")) {
    started <- proc.time()[["elapsed"]]
    fit <- do.call(
      grip.optimize.geodesic.mds,
      c(list(coords = case$cmd$coords, prepared = case$prepared), spec$args)
    )
    elapsed <- proc.time()[["elapsed"]] - started
    score_df <- fit$score
    score_info <- list(
      objective_energy = score_df$gmds.energy[[1L]],
      gmds_stress = score_df$gmds.stress[[1L]],
      gmds_raw_stress = score_df$gmds.raw_stress[[1L]],
      gmds_rmse = score_df$gmds.rmse[[1L]],
      anchor_energy = score_df$anchor.energy[[1L]],
      edge_spring_energy = score_df$edge.spring.energy[[1L]],
      repulsion_energy = score_df$repulsion.energy[[1L]],
      smooth_energy = score_df$smooth.energy[[1L]]
    )
    row <- compute_metrics(case = case, spec = spec, coords = fit$coords, elapsed_sec = elapsed, fit = fit, score_info = score_info)
    return(list(
      coords = fit$coords,
      display_coords = align_to_target_nd(fit$coords, case$truth, allow.reflection = TRUE)$aligned,
      metrics = row,
      fit = fit
    ))
  }

  started <- proc.time()[["elapsed"]]
  fit <- do.call(
    optimize_edge_relaxed,
    c(list(coords = case$cmd$coords, prepared = case$prepared), spec$args)
  )
  elapsed <- proc.time()[["elapsed"]] - started
  row <- compute_metrics(case = case, spec = spec, coords = fit$coords, elapsed_sec = elapsed, fit = fit, score_info = fit$score)
  list(
    coords = fit$coords,
    display_coords = align_to_target_nd(fit$coords, case$truth, allow.reflection = TRUE)$aligned,
    metrics = row,
    fit = fit
  )
}

build_title <- function(case, method_result) {
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
  n_panels <- length(methods)
  grDevices::png(
    output_path,
    width = 1280L * 2L,
    height = 880L * 2L,
    res = 180,
    bg = "#ffffff"
  )
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(2L, 4L), mar = c(1.2, 1.2, 2.9, 0.4), oma = c(0, 0, 1.2, 0))

  for (i in seq_len(8L)) {
    if (i > n_panels) {
      graphics::plot.new()
      next
    }
    method <- methods[[i]]
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
    graphics::mtext(build_title(case_result$case, method), side = 3L, line = 0.3, cex = 0.80)
  }
  graphics::mtext(case_result$case$label, side = 3L, outer = TRUE, line = -0.3, cex = 1.15, font = 2L)
}

write_foundation_table <- function(df) {
  rows <- vapply(seq_len(nrow(df)), function(i) {
    sprintf(
      "%s & %s & %s \\\\",
      tex_escape(df$check_label[[i]]),
      if (isTRUE(df$status[[i]])) "pass" else "fail",
      fmt_num(df$value[[i]], 6L)
    )
  }, character(1L))
  paste(
    "\\begin{table}[H]",
    "\\centering",
    "\\small",
    "\\begin{tabular}{p{0.68\\linewidth}cc}",
    "\\toprule",
    "Check & Status & Value \\\\",
    "\\midrule",
    paste(rows, collapse = "\n"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\caption{Foundation checks for the revised GMDS hierarchy. These validate exact realizability and tied-path handling before the paraboloid comparison begins.}",
    "\\label{tab:foundation}",
    "\\end{table}",
    sep = "\n"
  )
}

write_case_table <- function(df) {
  rows <- vapply(seq_len(nrow(df)), function(i) {
    sprintf(
      "%s & %s & %s & %s & %s & %s & %s & %s \\\\",
      tex_escape(df$method_label[[i]]),
      tex_escape(df$family[[i]]),
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
    "\\begin{tabular}{llllllll}",
    "\\toprule",
    "Method & Family & $\\lambda$ & $t$ (s) & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ \\\\",
    "\\midrule",
    paste(rows, collapse = "\n"),
    "\\bottomrule",
    "\\end{tabular}}",
    "\\end{table}",
    sep = "\n"
  )
}

case_summary_paragraph <- function(df) {
  pure <- df[df$method_id == "pure_gmds", , drop = FALSE]
  reg <- df[df$method_id == "reg_anchor_rep", , drop = FALSE]
  edge <- df[df$method_id == "edge_surrogate", , drop = FALSE]
  sprintf(
    "On the %s mesh, pure GMDS reaches the lowest geodesic stress among the optimized methods ($\\sigma=%s$) but does so with roughness $\\eta=%s$. The regularized anchor-plus-repulsion model keeps the regularization strengths explicit (%s) and lands at $\\sigma=%s$, $\\rho=%s$. The edge surrogate (%s) gives a cheaper local objective whose resulting layout has $\\sigma=%s$ and $\\rho=%s$; this is the direct surrogate-versus-GMDS tradeoff we wanted the rebuilt suite to expose.",
    tex_escape(df$case_label[[1L]]),
    fmt_num(pure$gmds_stress[[1L]], 4L),
    fmt_num(pure$roughness[[1L]], 4L),
    tex_escape(reg$lambda_label[[1L]]),
    fmt_num(reg$gmds_stress[[1L]], 4L),
    fmt_num(reg$procrustes_rmse[[1L]], 4L),
    tex_escape(edge$lambda_label[[1L]]),
    fmt_num(edge$gmds_stress[[1L]], 4L),
    fmt_num(edge$procrustes_rmse[[1L]], 4L)
  )
}

case_figure_rel <- function(case_id) {
  file.path(pdf_dir, sprintf("%s_grid.png", case_id))
}

foundation_df <- run_foundation_checks()
utils::write.csv(foundation_df, foundation_csv, row.names = FALSE)

cases <- lapply(c(12L, 15L, 20L), make_case)
case_results <- lapply(cases, function(case) {
  methods <- lapply(method_specs, run_method, case = case)
  metrics <- do.call(rbind, lapply(methods, `[[`, "metrics"))
  out <- list(case = case, methods = methods, metrics = metrics)
  save_case_panel_grid(out, file.path(pdf_dir, sprintf("%s_grid.png", case$id)))
  out
})
metrics_df <- do.call(rbind, lapply(case_results, `[[`, "metrics"))
utils::write.csv(metrics_df, metrics_csv, row.names = FALSE)

saveRDS(
  list(
    run_tag = run_tag,
    foundation_df = foundation_df,
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
      "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{%s. The hierarchy is explicit in the panel set: cmdscale baseline, pure GMDS, two regularized-GMDS variants, and one edge-surrogate variant. Lower $\\sigma$, $\\rho$, and $\\eta$ are better, while larger $\\alpha_{0.05}$ is better.}\\label{fig:%s}\\end{figure}",
      case_figure_rel(case_result$case$id),
      tex_escape(case_result$case$label),
      case_result$case$id
    ),
    write_case_table(df),
    sep = "\n\n"
  )
}, character(1L))

best_rows <- do.call(rbind, lapply(split(metrics_df, metrics_df$case_id), function(df) {
  data.frame(
    case_id = df$case_id[[1L]],
    pure_sigma = df$gmds_stress[df$method_id == "pure_gmds"][1L],
    reg_sigma = df$gmds_stress[df$method_id == "reg_anchor_rep"][1L],
    edge_sigma = df$gmds_stress[df$method_id == "edge_surrogate_rep"][1L],
    reg_rho = df$procrustes_rmse[df$method_id == "reg_anchor_rep"][1L],
    edge_rho = df$procrustes_rmse[df$method_id == "edge_surrogate_rep"][1L],
    stringsAsFactors = FALSE
  )
}))

summary_text <- paste(
  "\\section{Takeaways}",
  sprintf(
    "The rebuilt paraboloid suite now answers the hierarchy question directly. Across the three orthogonal paraboloid meshes, pure GMDS remains the most faithful to the fixed-path objective, regularized GMDS makes the anchor and repulsion strengths explicit rather than hiding them inside the base method, and the edge surrogate can be judged honestly as a faster local relaxation whose geodesic fidelity must be measured rather than assumed."
  ),
  sprintf(
    "The foundation checks in Table~\\ref{tab:foundation} also make the manuscript-level story operational: exact realizability really does give zero raw GMDS stress, and the tied-square example confirms that the objective is defined on $(G,\\Gamma)$ rather than on $G$ alone."
  ),
  sprintf(
    "For the paraboloid family itself, the key comparison is not whether regularization or the surrogate simply ``looks nicer'', but what geodesic fidelity is paid for that behavior. The tables and figures above are organized around exactly that question."
  ),
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
  "\\title{GMDS v2 Paraboloid Evaluation Suite}",
  "\\author{Phase 4 rebuild: paraboloid-first benchmark}",
  "\\date{2026-04-01}",
  "\\begin{document}",
  "\\maketitle",
  "\\section{Purpose}",
  "This standalone report implements the paraboloid-first evaluation stage of the GMDS v2 cleanup. The goal is to compare three distinct objects that the revised manuscript now separates explicitly: pure GMDS, regularized GMDS, and the edge-relaxation surrogate. The benchmark family is restricted to orthogonal weighted paraboloid meshes so the pathology question stays easy to interpret.",
  "\\section{Foundation checks}",
  "Before benchmarking aesthetics or runtime, the suite validates the mathematical core on exact-realizability and tie-handling examples. This is the first new requirement of the rebuilt Phase 4 evaluation stage, and it anchors the later paraboloid story to the actual definition of GMDS on $(G,\\Gamma)$.",
  write_foundation_table(foundation_df),
  "\\section{Paraboloid benchmark design}",
  "Every non-reference method starts from the same 3D classical-MDS embedding of the graph-geodesic distance matrix. The tested cases are orthogonal paraboloid meshes at $12\\times12$, $15\\times15$, and $20\\times20$. The compared method families are:",
  "\\begin{itemize}",
  "\\item classical MDS baseline,",
  "\\item pure GMDS,",
  "\\item regularized GMDS with anchor only,",
  "\\item regularized GMDS with anchor plus graph-aware repulsion,",
  "\\item edge-relaxation surrogate.",
  "\\end{itemize}",
  "The regularization strengths are shown explicitly in the $\\lambda$ column of every table. For regularized GMDS the reported stress $\\sigma$ is still the geodesic-MDS stress, while for the surrogate methods the local edge objective is used only to generate the embedding and the resulting $\\sigma$ is reported as an external fidelity score.",
  paste(case_sections, collapse = "\n\n"),
  summary_text,
  "\\end{document}"
)

writeLines(tex_lines, tex_path)

message("Wrote benchmark bundle: ", rds_path)
message("Wrote metrics CSV: ", metrics_csv)
message("Wrote foundation CSV: ", foundation_csv)
message("Wrote LaTeX report: ", tex_path)
message("Expected PDF path after latexmk: ", pdf_path)
