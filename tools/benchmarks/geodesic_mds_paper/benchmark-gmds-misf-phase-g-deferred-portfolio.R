#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) == 0L) {
  stop("Unable to determine script path from commandArgs().")
}
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

args <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% args

run_tag <- if (smoke) {
  sprintf("gmds-misf-phase-g-deferred-portfolio-smoke-%s", format(Sys.Date(), "%Y-%m-%d"))
} else {
  "gmds-misf-phase-g-deferred-portfolio-2026-04-05"
}

design_root <- file.path(repo_root, "dev", "design")
paper_design_root <- file.path(design_root, "geodesic_mds")
report_name <- "gmds_misf_phase_g_deferred_portfolio_report_2026-04-05"
report_source_dir <- file.path(paper_design_root, "reports", report_name)
report_output_dir <- file.path(repo_root, "output", "geodesic_mds_paper", "reports", report_name)
tmp_dir <- file.path(design_root, "tmp", run_tag)
pdf_dir <- file.path(report_source_dir, run_tag)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

tex_path <- file.path(report_source_dir, paste0(report_name, ".tex"))
pdf_path <- file.path(report_output_dir, paste0(report_name, ".pdf"))
rds_path <- file.path(tmp_dir, "gmds_misf_phase_g_deferred_portfolio_results.rds")
selection_csv <- file.path(tmp_dir, "gmds_misf_phase_g_deferred_portfolio_selection_metrics.csv")
final_csv <- file.path(tmp_dir, "gmds_misf_phase_g_deferred_portfolio_final_metrics.csv")
summary_csv <- file.path(tmp_dir, "gmds_misf_phase_g_deferred_portfolio_summary_metrics.csv")
shortlist_csv <- file.path(tmp_dir, "gmds_misf_phase_g_deferred_portfolio_shortlist_metrics.csv")
portfolio_csv <- file.path(tmp_dir, "gmds_misf_phase_g_deferred_portfolio_portfolio_metrics.csv")

for (pkg in c("igraph", "interp")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required to run this benchmark.", pkg))
  }
}

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
active_level_vertices <- get("grip.geodesic.misf.active.level.vertices", envir = ns)
induced_active_graph <- get("grip.geodesic.misf.induced_active_graph", envir = ns)
place_level_with_layout <- get("grip.geodesic.misf.place.level.with.layout", envir = ns)
prepare_active_level <- get("grip.geodesic.misf.prepare.active.level", envir = ns)
classical_mds_embedding <- get("grip.classical.mds.embedding", envir = ns)

cfg <- list(
  phase_seed = 20260402L,
  amplitude = 0.35,
  dim = 3L,
  candidate_ids = c("cmdscale", "kk", "grip", "weighted_grip", "weighted_grip_polish_lgkk"),
  correction_anchor_weight = 0.10,
  correction_anchor_weight_end = 0.02,
  correction_continuation = "linear",
  correction_regularized_max_iter = if (smoke) 4L else 8L,
  correction_pure_max_iter = if (smoke) 4L else 8L,
  insertion_layout_k = 6L,
  lookahead_refine_local_nbrs = if (smoke) 3L else 4L,
  lookahead_refine_landmarks = 2L,
  lookahead_refine_max_iter = if (smoke) 2L else 3L,
  full_refine_local_nbrs = if (smoke) 3L else 4L,
  full_refine_landmarks = 2L,
  full_refine_max_iter = if (smoke) 2L else 3L,
  final_polish_max_iter = if (smoke) 2L else 4L,
  trace_frames = if (smoke) 3L else 4L,
  fr_niter = if (smoke) 240L else 640L,
  weighted_polish_lgkk_rounds = if (smoke) 3L else 8L,
  n_threads = 0L
)

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

rank_smallest <- function(x) {
  x <- as.double(x)
  finite <- is.finite(x)
  if (!any(finite)) {
    return(rep.int(1L, length(x)))
  }
  fill <- max(x[finite], na.rm = TRUE) + max(1, stats::sd(x[finite]), na.rm = TRUE) + 1
  x[!finite] <- fill
  rank(x, ties.method = "min")
}

rank_largest <- function(x) {
  rank_smallest(-as.double(x))
}

tex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x, perl = TRUE)
  x
}

normalize_coords <- function(coords) {
  coords <- as.matrix(coords)
  centered <- sweep(coords, 2L, colMeans(coords), FUN = "-", check.margin = FALSE)
  radius <- max(sqrt(rowSums(centered^2)))
  if (!is.finite(radius) || radius <= 0) {
    return(centered)
  }
  centered / radius
}

center_coords <- function(coords) {
  coords <- as.matrix(coords)
  sweep(coords, 2L, colMeans(coords), "-", check.margin = FALSE)
}

build_shared_start <- function(n, dim, seed) {
  set.seed(as.integer(seed))
  coords <- matrix(stats::rnorm(as.integer(n) * as.integer(dim)), ncol = as.integer(dim))
  storage.mode(coords) <- "double"
  normalize_coords(coords)
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

triangles_from_param <- function(param_coords) {
  param_coords <- as.matrix(param_coords)
  if (nrow(param_coords) < 3L) {
    return(matrix(integer(), ncol = 3L))
  }
  tri <- tryCatch({
    mesh <- interp::tri.mesh(param_coords[, 1L], param_coords[, 2L])
    out <- interp::triangles(mesh)[, 1:3, drop = FALSE]
    matrix(as.integer(out), ncol = 3L)
  }, error = function(e) {
    matrix(integer(), ncol = 3L)
  })
  if (!nrow(tri)) {
    return(tri)
  }
  tri
}

edges_from_triangles <- function(triangles) {
  triangles <- as.matrix(triangles)
  if (!nrow(triangles)) {
    return(matrix(integer(), ncol = 2L))
  }
  edges <- rbind(
    triangles[, c(1L, 2L), drop = FALSE],
    triangles[, c(2L, 3L), drop = FALSE],
    triangles[, c(1L, 3L), drop = FALSE]
  )
  edges <- t(apply(edges, 1L, sort))
  edges <- unique(edges)
  matrix(as.integer(edges), ncol = 2L)
}

build_display_mesh <- function(param_coords) {
  param_coords <- as.matrix(param_coords)
  if (ncol(param_coords) < 2L) {
    param_coords <- cbind(param_coords, rep(0, nrow(param_coords)))
  }
  triangles <- triangles_from_param(param_coords[, 1:2, drop = FALSE])
  edges <- edges_from_triangles(triangles)
  if (!nrow(edges) && nrow(param_coords) >= 2L) {
    dmat <- as.matrix(stats::dist(param_coords[, 1:2, drop = FALSE]))
    diag(dmat) <- Inf
    edges <- do.call(rbind, lapply(seq_len(nrow(param_coords)), function(i) {
      j <- order(dmat[i, ])[seq_len(min(4L, nrow(param_coords) - 1L))]
      cbind(i, j)
    }))
    edges <- t(apply(edges, 1L, sort))
    edges <- unique(edges)
    edges <- matrix(as.integer(edges), ncol = 2L)
  }
  adj <- grip.build.adj.from.edges(edges, n = nrow(param_coords))$adj_list
  list(
    triangles = triangles,
    edges = edges,
    adj_list = adj
  )
}

build_sparse_layout_graph <- function(distance_matrix, k = 6L) {
  distance_matrix <- as.matrix(distance_matrix)
  n <- nrow(distance_matrix)
  if (n <= 1L) {
    return(list(edges = matrix(integer(), ncol = 2L), edge_weights = numeric(0L)))
  }
  k <- max(1L, min(as.integer(k), n - 1L))
  diag(distance_matrix) <- Inf

  knn_edges <- do.call(rbind, lapply(seq_len(n), function(i) {
    nbrs <- order(distance_matrix[i, ])[seq_len(k)]
    cbind(i, nbrs)
  }))
  knn_edges <- t(apply(knn_edges, 1L, sort))
  knn_edges <- unique(knn_edges)

  complete_edges <- t(utils::combn(n, 2L))
  complete_graph <- igraph::graph_from_edgelist(complete_edges, directed = FALSE)
  igraph::E(complete_graph)$weight <- distance_matrix[complete_edges]
  mst_graph <- igraph::mst(complete_graph, weights = igraph::E(complete_graph)$weight)
  mst_edges <- igraph::ends(mst_graph, igraph::E(mst_graph), names = FALSE)

  edges <- unique(rbind(knn_edges, mst_edges))
  edges <- edges[order(edges[, 1L], edges[, 2L]), , drop = FALSE]
  list(
    edges = matrix(as.integer(edges), ncol = 2L),
    edge_weights = as.double(distance_matrix[edges])
  )
}

triangle_areas <- function(coords, triangles) {
  coords <- as.matrix(coords)
  triangles <- as.matrix(triangles)
  if (!nrow(triangles)) {
    return(numeric(0L))
  }
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
  if (!length(areas)) {
    return(NA_real_)
  }
  med <- stats::median(areas)
  if (!is.finite(med) || med <= 0) {
    return(NA_real_)
  }
  as.double(stats::quantile(areas, probs = 0.05, names = FALSE)) / med
}

edge_floor_ratio <- function(coords, edges) {
  coords <- as.matrix(coords)
  edges <- as.matrix(edges)
  if (!nrow(edges)) {
    return(NA_real_)
  }
  lengths <- sqrt(rowSums(
    (coords[edges[, 1L], , drop = FALSE] - coords[edges[, 2L], , drop = FALSE])^2
  ))
  med <- stats::median(lengths)
  if (!is.finite(med) || med <= 0) {
    return(NA_real_)
  }
  as.double(stats::quantile(lengths, probs = 0.05, names = FALSE)) / med
}

mesh_roughness <- function(coords, adj_list, edges) {
  centered <- sweep(coords, 2L, colMeans(coords), FUN = "-", check.margin = FALSE)
  if (!nrow(edges)) {
    return(NA_real_)
  }
  median.edge <- stats::median(sqrt(rowSums(
    (centered[edges[, 1L], , drop = FALSE] - centered[edges[, 2L], , drop = FALSE])^2
  )))
  if (!is.finite(median.edge) || median.edge <= 0) {
    return(NA_real_)
  }
  residuals <- vapply(seq_len(nrow(centered)), function(i) {
    nbrs <- adj_list[[i]]
    if (!length(nbrs)) {
      return(0)
    }
    delta <- centered[i, ] - colMeans(centered[nbrs, , drop = FALSE])
    sum(delta^2)
  }, numeric(1L))
  sqrt(mean(residuals)) / median.edge
}

align_partial_to_truth <- function(coords, target) {
  coords <- as.matrix(coords)
  keep <- stats::complete.cases(coords)
  out <- matrix(NA_real_, nrow = nrow(coords), ncol = ncol(coords))
  if (sum(keep) >= 3L) {
    aligned <- align_to_target_nd(coords[keep, , drop = FALSE], target[keep, , drop = FALSE], allow.reflection = TRUE)
    out[keep, ] <- aligned$aligned
    return(list(aligned = out, rmse = aligned$rmse))
  }
  out[keep, ] <- coords[keep, , drop = FALSE]
  list(aligned = out, rmse = NA_real_)
}

complete_partial_coords <- function(coords) {
  coords <- as.matrix(coords)
  keep <- stats::complete.cases(coords)
  if (all(keep)) {
    return(list(coords = coords, missing = which(!keep)))
  }
  if (!any(keep)) {
    stop("Cannot complete a partial coordinate matrix with no finite rows.")
  }
  fill.center <- colMeans(coords[keep, , drop = FALSE])
  out <- coords
  out[!keep, ] <- matrix(fill.center, nrow = sum(!keep), ncol = ncol(coords), byrow = TRUE)
  list(coords = out, missing = which(!keep))
}

plot_partial_layout <- function(coords,
                                edges,
                                projection = "ortho",
                                azimuth = 35,
                                elevation = 24,
                                vertex.col = "#355070",
                                edge.col = "#adb5bd",
                                main = "") {
  coords <- as.matrix(coords)
  keep <- stats::complete.cases(coords)
  if (!any(keep)) {
    graphics::plot.new()
    return(invisible(NULL))
  }
  index.map <- integer(nrow(coords))
  index.map[keep] <- seq_len(sum(keep))
  display.edges <- as.matrix(edges)
  if (nrow(display.edges)) {
    good.edges <- keep[display.edges[, 1L]] & keep[display.edges[, 2L]]
    display.edges <- display.edges[good.edges, , drop = FALSE]
    if (nrow(display.edges)) {
      display.edges <- cbind(
        as.integer(index.map[display.edges[, 1L]]),
        as.integer(index.map[display.edges[, 2L]])
      )
    }
  }
  grip.plot(
    coords = coords[keep, , drop = FALSE],
    edges = display.edges,
    projection = projection,
    azimuth = azimuth,
    elevation = elevation,
    vertex.col = vertex.col,
    edge.col = edge.col,
    main = main
  )
  invisible(NULL)
}

select_trace_indices <- function(n_frames, max_frames = 4L) {
  if (n_frames <= 0L) {
    return(integer(0L))
  }
  max_frames <- max(1L, as.integer(max_frames))
  idx <- unique(round(seq(1, n_frames, length.out = min(max_frames, n_frames))))
  as.integer(idx)
}

extract_gmds_trace_frames <- function(fit, target, max_frames = 4L) {
  if (is.null(fit$frames) || !length(fit$frames)) {
    return(list())
  }
  idx <- select_trace_indices(length(fit$frames), max_frames = max_frames)
  trace_df <- if (is.data.frame(fit$trace)) fit$trace else NULL
  lapply(idx, function(i) {
    frame.coords <- as.matrix(fit$frames[[i]])
    aligned <- align_to_target_nd(frame.coords, target, allow.reflection = TRUE)
    trace_row <- if (!is.null(trace_df) && nrow(trace_df) >= i) trace_df[i, , drop = FALSE] else NULL
    list(
      frame_index = as.integer(i),
      total_frames = as.integer(length(fit$frames)),
      iteration = if (is.null(trace_row) || !"iteration" %in% names(trace_row)) NA_integer_ else as.integer(trace_row$iteration[[1L]]),
      energy = if (is.null(trace_row) || !"energy" %in% names(trace_row)) NA_real_ else as.double(trace_row$energy[[1L]]),
      gmds_energy = if (is.null(trace_row) || !"gmds_energy" %in% names(trace_row)) NA_real_ else as.double(trace_row$gmds_energy[[1L]]),
      anchor_weight = if (is.null(trace_row) || !"anchor_weight" %in% names(trace_row)) NA_real_ else as.double(trace_row$anchor_weight[[1L]]),
      display_coords = aligned$aligned
    )
  })
}

resolve_unweighted_grip_args <- function(family) {
  if (identical(family, "irregular_rectangle")) {
    return(list(
      placement = "barycenter",
      rounds = 96L,
      final_rounds = 128L,
      num_init = 18L,
      num_nbrs = 24L,
      r = 0.05,
      s = 6.5,
      repulsion_factor = 1.10
    ))
  }
  list(
    placement = "barycenter",
    rounds = 72L,
    final_rounds = 96L,
    num_init = 12L,
    num_nbrs = 18L,
    r = 0.08,
    s = 5.0,
    repulsion_factor = 1.25
  )
}

resolve_weighted_grip_args <- function(top_n) {
  list(
    rounds = 96L,
    final_rounds = 128L,
    num_init = min(18L, max(4L, top_n - 1L)),
    num_nbrs = min(20L, max(4L, top_n - 1L))
  )
}

surface_title <- function(surface) {
  switch(
    surface,
    flat = "flat",
    saddle = "saddle",
    paraboloid = "paraboloid",
    ripple = "ripple",
    surface
  )
}

make_case <- function(side, surface, family = c("regular", "irregular_rectangle")) {
  family <- match.arg(family)
  bundle <- if (identical(family, "regular")) {
    if (identical(surface, "flat")) {
      irregular.rectangle.surface.graph(
        side,
        side,
        surface = "flat",
        amplitude = cfg$amplitude,
        row_irregularity = 0,
        col_irregularity = 0,
        interior_warp = 0,
        shear = 0,
        min_step_ratio = 1,
        connectivity = "orthogonal",
        normalize = "median"
      )
    } else {
      mesh.surface.graph(
        side,
        side,
        surface = surface,
        amplitude = cfg$amplitude,
        connectivity = "orthogonal",
        normalize = "median"
      )
    }
  } else {
    irregular.rectangle.surface.graph(
      side,
      side,
      surface = surface,
      amplitude = cfg$amplitude,
      connectivity = "orthogonal",
      normalize = "median"
    )
  }
  surface.label <- surface_title(surface)
  misf.seed <- if (identical(family, "regular")) {
    cfg$phase_seed + side + match(surface, c("flat", "saddle", "paraboloid", "ripple")) * 100L
  } else {
    cfg$phase_seed + side + 1000L + match(surface, c("flat", "saddle", "paraboloid", "ripple")) * 100L
  }
  prepared <- grip.prepare.misf.geodesic.mds(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average",
    dim = cfg$dim,
    top_level_mode = "skip",
    seed = misf.seed
  )
  top_graph <- prepared$top_level_graph
  layout_graph <- build_sparse_layout_graph(top_graph$distance_matrix, k = cfg$insertion_layout_k)
  top_param <- as.matrix(bundle$coords_param[top_graph$vertex_ids, , drop = FALSE])
  top_truth <- as.matrix(bundle$coords_surface[top_graph$vertex_ids, , drop = FALSE])
  top_display <- build_display_mesh(top_param)
  started <- proc.time()[["elapsed"]]
  cmd <- classical_mds_embedding(prepared, dim = cfg$dim, eig = TRUE)
  cmd.elapsed <- proc.time()[["elapsed"]] - started
  list(
    id = sprintf("%s_%s_%dx%d", family, surface, side, side),
    label = sprintf(
      "%s %s mesh %dx%d",
      if (identical(family, "regular")) "Regular" else "Irregular rectangle",
      surface.label,
      side,
      side
    ),
    family = family,
    surface = surface,
    side = as.integer(side),
    n = bundle$n,
    edges = bundle$edges,
    triangles = grid_mesh_triangles(side, side),
    truth = bundle$coords_surface,
    prepared = prepared,
    adj_list = prepared$adj_list,
    cmd = cmd,
    cmd_elapsed = as.double(cmd.elapsed),
    misf_seed = misf.seed,
    top_level = prepared$top_level_level,
    top_n = top_graph$n,
    top_vertex_ids = top_graph$vertex_ids,
    top_graph = top_graph,
    top_prepared = prepared$top_level_prepared,
    layout_graph = layout_graph,
    top_display_edges = top_display$edges,
    top_display_triangles = top_display$triangles,
    top_display_adj = top_display$adj_list,
    weighted_preset = if (identical(family, "regular")) "mesh" else "irregular",
    unweighted_grip_args = resolve_unweighted_grip_args(family),
    weighted_grip_args = resolve_weighted_grip_args(top_graph$n),
    shared_start = build_shared_start(top_graph$n, cfg$dim, misf.seed + 77L)
  )
}

method_specs <- list(
  list(id = "cmdscale", label = "cMDS"),
  list(id = "kk", label = "KK"),
  list(id = "grip", label = "GRIP"),
  list(id = "weighted_grip", label = "Weighted GRIP"),
  list(id = "weighted_grip_polish_lgkk", label = "Weighted GRIP + polish LGKK")
)

run_top_initializer <- function(case, method_spec) {
  method_id <- method_spec$id
  method_label <- method_spec$label
  graph_obj <- igraph::graph_from_edgelist(case$layout_graph$edges, directed = FALSE)
  note <- NULL

  if (identical(method_id, "cmdscale")) {
    started <- proc.time()[["elapsed"]]
    fit <- classical_mds_embedding(case$top_prepared, dim = cfg$dim, eig = TRUE)
    elapsed <- proc.time()[["elapsed"]] - started
    coords <- fit$coords
  } else if (identical(method_id, "kk")) {
    started <- proc.time()[["elapsed"]]
    coords <- igraph::layout_with_kk(
      graph_obj,
      coords = case$shared_start,
      dim = cfg$dim,
      weights = case$layout_graph$edge_weights
    )
    elapsed <- proc.time()[["elapsed"]] - started
    note <- "Shared generic random-normal start"
  } else if (identical(method_id, "grip")) {
    args <- case$unweighted_grip_args
    started <- proc.time()[["elapsed"]]
    tr <- grip.layout.trace(
      edges = case$layout_graph$edges,
      n = case$top_n,
      dim = cfg$dim,
      placement = args$placement,
      rounds = args$rounds,
      final_rounds = args$final_rounds,
      num_init = min(args$num_init, max(2L, case$top_n - 1L)),
      num_nbrs = min(args$num_nbrs, max(2L, case$top_n - 1L)),
      r = args$r,
      s = args$s,
      repulsion_factor = args$repulsion_factor,
      seed = case$misf_seed + 301L
    )
    elapsed <- proc.time()[["elapsed"]] - started
    coords <- tr$final
    note <- sprintf("Trace snapshots: %d selected of %d total frames", length(tr$frames), length(tr$frames))
  } else if (identical(method_id, "weighted_grip")) {
    started <- proc.time()[["elapsed"]]
    tr <- grip.layout.trace.weighted(
      edges = case$layout_graph$edges,
      edge_weights = case$layout_graph$edge_weights,
      n = case$top_n,
      dim = cfg$dim,
      preset = case$weighted_preset,
      seed = case$misf_seed + 401L
    )
    elapsed <- proc.time()[["elapsed"]] - started
    coords <- tr$final
    note <- sprintf("Trace snapshots: %d selected of %d total frames", length(tr$frames), length(tr$frames))
  } else if (identical(method_id, "weighted_grip_polish_lgkk")) {
    started <- proc.time()[["elapsed"]]
    tr <- grip.layout.trace.weighted(
      edges = case$layout_graph$edges,
      edge_weights = case$layout_graph$edge_weights,
      n = case$top_n,
      dim = cfg$dim,
      preset = case$weighted_preset,
      lgkk_polish_rounds = cfg$weighted_polish_lgkk_rounds,
      seed = case$misf_seed + 601L
    )
    elapsed <- proc.time()[["elapsed"]] - started
    coords <- tr$final
    note <- sprintf("Trace snapshots: %d selected of %d total frames", length(tr$frames), length(tr$frames))
  } else {
    stop("Unknown method id: ", method_id)
  }

  top.score <- grip.score.geodesic.mds(coords = coords, prepared = case$top_prepared)
  top.aligned <- align_to_target_nd(coords, case$truth[case$top_vertex_ids, , drop = FALSE], allow.reflection = TRUE)
  metrics <- data.frame(
    case_id = case$id,
    case_label = case$label,
    family = case$family,
    surface = case$surface,
    side = case$side,
    method_id = method_id,
    method_label = method_label,
    elapsed_sec = as.double(elapsed),
    top_sigma = as.double(top.score$gmds.stress[[1L]]),
    top_rho = as.double(top.aligned$rmse),
    note = if (is.null(note)) "" else as.character(note),
    stringsAsFactors = FALSE
  )

  list(
    method_id = method_id,
    method_label = method_label,
    coords = as.matrix(coords),
    metrics = metrics
  )
}

run_top_correction <- function(case, initializer_result) {
  initializer.coords <- center_coords(initializer_result$coords)

  started <- proc.time()[["elapsed"]]
  anchor.fit <- grip.optimize.geodesic.mds(
    coords = initializer.coords,
    prepared = case$top_prepared,
    init = "user",
    anchor_mode = "user",
    anchor_coords = initializer.coords,
    anchor_weight = cfg$correction_anchor_weight,
    anchor_weight_end = cfg$correction_anchor_weight_end,
    continuation = cfg$correction_continuation,
    engine = "cpp",
    n_threads = cfg$n_threads,
    max_iter = cfg$correction_regularized_max_iter,
    return_trace = TRUE,
    recenter = TRUE
  )
  anchor.elapsed <- proc.time()[["elapsed"]] - started

  started <- proc.time()[["elapsed"]]
  pure.fit <- grip.optimize.geodesic.mds(
    coords = anchor.fit$coords,
    prepared = case$top_prepared,
    init = "user",
    anchor_mode = "none",
    engine = "cpp",
    n_threads = cfg$n_threads,
    max_iter = cfg$correction_pure_max_iter,
    return_trace = TRUE,
    recenter = TRUE
  )
  pure.elapsed <- proc.time()[["elapsed"]] - started

  pure.display <- align_partial_to_truth(
    coords = grip.geodesic.misf.partial.coords(pure.fit$coords, case$top_vertex_ids, case$n),
    target = case$truth
  )$aligned

  list(
    method_id = initializer_result$method_id,
    method_label = initializer_result$method_label,
    initializer = initializer_result,
    anchor = list(
      fit = anchor.fit,
      elapsed_sec = as.double(anchor.elapsed),
      trace_selected = extract_gmds_trace_frames(
        anchor.fit,
        case$truth[case$top_vertex_ids, , drop = FALSE],
        max_frames = cfg$trace_frames
      )
    ),
    pure = list(
      fit = pure.fit,
      elapsed_sec = as.double(pure.elapsed),
      trace_selected = extract_gmds_trace_frames(
        pure.fit,
        case$truth[case$top_vertex_ids, , drop = FALSE],
        max_frames = cfg$trace_frames
      ),
      display = pure.display
    )
  )
}

build_injected_top_fit <- function(case, corrected_method) {
  top.coords <- as.matrix(corrected_method$pure$fit$coords)
  top.prepared <- case$top_prepared
  top.score <- grip.score.geodesic.mds(coords = top.coords, prepared = top.prepared)
  seed.elapsed <- sum(c(
    corrected_method$initializer$metrics$elapsed_sec[[1L]],
    corrected_method$anchor$elapsed_sec,
    corrected_method$pure$elapsed_sec
  ), na.rm = TRUE)
  list(
    fit = list(
      coords = top.coords,
      trace = data.frame(iter = 0L, stringsAsFactors = FALSE),
      frames = list(top.coords, top.coords),
      prepared = top.prepared,
      score = top.score,
      restart_summary = data.frame(),
      best_restart = 1L,
      best_restart_row = data.frame(),
      vertex_ids = case$top_vertex_ids,
      coords_full = grip.geodesic.misf.partial.coords(top.coords, case$top_vertex_ids, case$n),
      injected = TRUE,
      injected_seed_label = corrected_method$method_label,
      injected_seed_elapsed_sec = seed.elapsed
    ),
    seed_label = corrected_method$method_label,
    seed_elapsed_sec = as.double(seed.elapsed),
    method_id = corrected_method$method_id
  )
}

compute_final_metrics <- function(case, method_id, method_label, coords, elapsed_sec, note = "") {
  score_df <- grip.score.geodesic.mds(coords = coords, prepared = case$prepared)
  aligned <- align_to_target_nd(coords, case$truth, allow.reflection = TRUE)
  data.frame(
    case_id = case$id,
    case_label = case$label,
    family = case$family,
    surface = case$surface,
    side = case$side,
    n = case$n,
    method_id = method_id,
    method_label = method_label,
    elapsed_sec = as.double(elapsed_sec),
    gmds_stress = score_df$gmds.stress[[1L]],
    gmds_energy = score_df$gmds.energy[[1L]],
    procrustes_rmse = aligned$rmse,
    roughness = mesh_roughness(coords, case$prepared$adj_list, case$edges),
    edge_floor = edge_floor_ratio(coords, case$edges),
    area_q05_ratio = area_floor_ratio(coords, case$triangles),
    note = note,
    stringsAsFactors = FALSE
  )
}

compute_lookahead_metrics <- function(case, prepared, coords, level) {
  active.vertices <- active_level_vertices(prepared, level)
  active.full.coords <- as.matrix(coords)
  active.coords <- active.full.coords[active.vertices, , drop = FALSE]
  active.prepared <- prepare_active_level(
    prepared = prepared,
    active_vertices = active.vertices,
    local_nbrs = cfg$lookahead_refine_local_nbrs,
    landmark_count = cfg$lookahead_refine_landmarks,
    pair_mode = "full"
  )
  active.graph <- induced_active_graph(prepared, active.vertices)
  active.score <- grip.score.geodesic.mds(active.full.coords, prepared = active.prepared)
  active.align <- align_to_target_nd(active.coords, case$truth[active.vertices, , drop = FALSE], allow.reflection = TRUE)
  data.frame(
    active_level = level,
    active_n = length(active.vertices),
    lookahead_sigma = active.score$gmds.stress[[1L]],
    lookahead_energy = active.score$gmds.energy[[1L]],
    lookahead_rho = active.align$rmse,
    lookahead_roughness = mesh_roughness(active.coords, active.graph$adj_list, active.graph$edges),
    lookahead_edge_floor = edge_floor_ratio(active.coords, active.graph$edges),
    stringsAsFactors = FALSE
  )
}

run_additional_lookahead <- function(case, prepared, coords, current_level) {
  next.level <- as.integer(current_level) - 1L
  if (!is.finite(next.level) || next.level < 0L) {
    return(NULL)
  }

  placement <- place_level_with_layout(
    prepared = prepared,
    coords = coords,
    level = next.level,
    method = "weighted_kk",
    layout_k = cfg$insertion_layout_k,
    seed = cfg$phase_seed + 5000L + case$side + next.level
  )
  placed.partial <- placement$coords
  placed.completed <- complete_partial_coords(placed.partial)
  refined <- grip:::grip.geodesic.misf.refine.level(
    prepared = prepared,
    coords = placed.completed$coords,
    level = next.level,
    local_nbrs = cfg$lookahead_refine_local_nbrs,
    landmark_count = cfg$lookahead_refine_landmarks,
    pair_mode = "sparse",
    anchor_weight = 0.05,
    anchor_weight_end = 0.01,
    continuation = "linear",
    max_iter = cfg$lookahead_refine_max_iter,
    engine = "cpp",
    n_threads = cfg$n_threads,
    recenter = FALSE,
    return_trace = FALSE
  )
  display.coords <- refined$coords
  if (length(placed.completed$missing)) {
    display.coords[placed.completed$missing, ] <- NA_real_
  }
  list(
    level = next.level,
    coords_full = refined$coords,
    display = align_partial_to_truth(display.coords, case$truth)$aligned,
    metrics = compute_lookahead_metrics(case, prepared, refined$coords, next.level)
  )
}

run_candidate <- function(case, corrected_method) {
  top.seed <- build_injected_top_fit(case, corrected_method)
  prepared <- case$prepared
  prepared$top_level_fit <- top.seed$fit
  top.display <- corrected_method$pure$display

  next.level <- prepared$top_level_level - 1L
  lookahead.coords <- top.seed$fit$coords_full
  lookahead.coords.full <- top.seed$fit$coords_full
  if (next.level >= 0L) {
    placement <- place_level_with_layout(
      prepared = prepared,
      coords = lookahead.coords,
      level = next.level,
      method = "weighted_kk",
      layout_k = cfg$insertion_layout_k,
      seed = cfg$phase_seed + case$side + match(corrected_method$method_id, cfg$candidate_ids)
    )
    lookahead.partial <- placement$coords
    lookahead.completed <- complete_partial_coords(lookahead.partial)
    refined <- grip:::grip.geodesic.misf.refine.level(
      prepared = prepared,
      coords = lookahead.completed$coords,
      level = next.level,
      local_nbrs = cfg$lookahead_refine_local_nbrs,
      landmark_count = cfg$lookahead_refine_landmarks,
      pair_mode = "sparse",
      anchor_weight = 0.05,
      anchor_weight_end = 0.01,
      continuation = "linear",
      max_iter = cfg$lookahead_refine_max_iter,
      engine = "cpp",
      n_threads = cfg$n_threads,
      recenter = FALSE,
      return_trace = FALSE
    )
    lookahead.metrics <- compute_lookahead_metrics(case, prepared, refined$coords, next.level)
    lookahead.coords.full <- refined$coords
    lookahead.coords <- refined$coords
    if (length(lookahead.completed$missing)) {
      lookahead.coords[lookahead.completed$missing, ] <- NA_real_
    }
  } else {
    lookahead.metrics <- data.frame(
      active_level = prepared$top_level_level,
      active_n = length(prepared$top_level_vertices),
      lookahead_sigma = top.seed$fit$score$gmds.stress[[1L]],
      lookahead_energy = top.seed$fit$score$gmds.energy[[1L]],
      lookahead_rho = align_to_target_nd(top.seed$fit$coords, case$truth[prepared$top_level_vertices, , drop = FALSE], allow.reflection = TRUE)$rmse,
      lookahead_roughness = NA_real_,
      lookahead_edge_floor = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  lookahead.display <- align_partial_to_truth(lookahead.coords, case$truth)$aligned

  full.start <- proc.time()[["elapsed"]]
  full.fit <- grip.optimize.misf.geodesic.mds(
    prepared = prepared,
    insertion_mode = "weighted_kk",
    insertion_layout_k = cfg$insertion_layout_k,
    refinement_local_nbrs = cfg$full_refine_local_nbrs,
    refinement_landmark_count = cfg$full_refine_landmarks,
    refinement_pair_mode = "sparse",
    refinement_anchor_weight = 0.05,
    refinement_anchor_weight_end = 0.01,
    refinement_continuation = "linear",
    refinement_max_iter = cfg$full_refine_max_iter,
    refinement_engine = "cpp",
    final_polish_max_iter = cfg$final_polish_max_iter,
    final_polish_engine = "cpp",
    n_threads = cfg$n_threads,
    return_trace = FALSE,
    return_frames = FALSE,
    seed = cfg$phase_seed + case$side
  )
  full.elapsed <- proc.time()[["elapsed"]] - full.start
  final.metrics <- compute_final_metrics(
    case = case,
    method_id = corrected_method$method_id,
    method_label = corrected_method$method_label,
    coords = full.fit$coords,
    elapsed_sec = full.elapsed,
    note = "Phase G benchmark full weighted-KK lower-level placement pipeline"
  )

  list(
    method_id = corrected_method$method_id,
    method_label = corrected_method$method_label,
    top_seed = top.seed,
    top_display = top.display,
    lookahead_display = lookahead.display,
    lookahead_coords_full = lookahead.coords.full,
    lookahead_level = next.level,
    final_display = align_to_target_nd(full.fit$coords, case$truth, allow.reflection = TRUE)$aligned,
    lookahead_metrics = lookahead.metrics,
    final_metrics = final.metrics
  )
}

run_candidate_with_frames <- function(case, corrected_method) {
  top.seed <- build_injected_top_fit(case, corrected_method)
  prepared <- case$prepared
  prepared$top_level_fit <- top.seed$fit
  fit <- grip.optimize.misf.geodesic.mds(
    prepared = prepared,
    insertion_mode = "weighted_kk",
    insertion_layout_k = cfg$insertion_layout_k,
    refinement_local_nbrs = cfg$full_refine_local_nbrs,
    refinement_landmark_count = cfg$full_refine_landmarks,
    refinement_pair_mode = "sparse",
    refinement_anchor_weight = 0.05,
    refinement_anchor_weight_end = 0.01,
    refinement_continuation = "linear",
    refinement_max_iter = cfg$full_refine_max_iter,
    refinement_engine = "cpp",
    final_polish_max_iter = cfg$final_polish_max_iter,
    final_polish_engine = "cpp",
    n_threads = cfg$n_threads,
    return_trace = FALSE,
    return_frames = TRUE,
    seed = cfg$phase_seed + case$side
  )
  list(
    top_level = align_partial_to_truth(fit$frames$after_top_level, case$truth)$aligned,
    after_insertion = align_partial_to_truth(fit$frames$after_insertion, case$truth)$aligned,
    after_refinement = align_partial_to_truth(fit$frames$after_refinement, case$truth)$aligned,
    final = align_partial_to_truth(fit$frames$final, case$truth)$aligned
  )
}

run_direct_control <- function(case) {
  started <- proc.time()[["elapsed"]]
  fit <- grip.optimize.geodesic.mds(
    coords = case$cmd$coords,
    prepared = case$prepared,
    init = "user",
    anchor_mode = "none",
    engine = "cpp",
    max_iter = cfg$final_polish_max_iter,
    n_threads = cfg$n_threads,
    return_trace = FALSE,
    recenter = TRUE
  )
  elapsed <- proc.time()[["elapsed"]] - started
  list(
    coords = fit$coords,
    display_coords = align_to_target_nd(fit$coords, case$truth, allow.reflection = TRUE)$aligned,
    metrics = compute_final_metrics(
      case = case,
      method_id = "cmd_pure_gmds",
      method_label = "cmdscale -> pure GMDS",
      coords = fit$coords,
      elapsed_sec = elapsed + case$cmd_elapsed,
      note = "Full-graph cMDS start followed by pure GMDS"
    )
  )
}

fixed_control_seed_id <- function(case) {
  if (identical(case$family, "irregular_rectangle")) {
    return("cmdscale")
  }
  if (identical(case$family, "regular") && case$side <= 12L) {
    return("grip")
  }
  "weighted_grip"
}

save_candidate_grid <- function(case_result, stage = c("top", "lookahead"), output_path) {
  stage <- match.arg(stage)
  entries <- c(
    list(list(
      display = case_result$case$truth,
      title = "Reference surface",
      subtitle = case_result$case$label
    )),
    lapply(case_result$candidates, function(candidate) {
      look <- candidate$lookahead_metrics[1L, , drop = FALSE]
      list(
        display = if (identical(stage, "top")) candidate$top_display else candidate$lookahead_display,
        title = candidate$method_label,
        subtitle = if (identical(stage, "top")) {
          "corrected top seed"
        } else {
          sprintf("look sigma %s, rho %s", fmt_num(look$lookahead_sigma[[1L]], 4L), fmt_num(look$lookahead_rho[[1L]], 4L))
        }
      )
    })
  )
  n <- length(entries)
  ncol.panel <- 3L
  nrow.panel <- ceiling(n / ncol.panel)
  grDevices::png(output_path, width = 2400L, height = 820L * nrow.panel, res = 180, bg = "#ffffff", type = "cairo")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(nrow.panel, ncol.panel), mar = c(1.1, 1.1, 3.0, 0.4), oma = c(0, 0, 1.0, 0))
  for (i in seq_len(nrow.panel * ncol.panel)) {
    if (i > length(entries)) {
      graphics::plot.new()
      next
    }
    entry <- entries[[i]]
    plot_partial_layout(
      coords = entry$display,
      edges = case_result$case$edges,
      projection = "ortho",
      azimuth = 35,
      elevation = 24,
      vertex.col = if (i == 1L) "#bc6c25" else "#3a5a40",
      edge.col = "#adb5bd",
      main = ""
    )
    graphics::mtext(sprintf("%s\n%s", entry$title, entry$subtitle), side = 3L, line = 0.3, cex = 0.82)
  }
  ttl <- if (identical(stage, "top")) "Phase G corrected top-level seeds" else "Phase E one-level weighted-KK lookahead"
  graphics::mtext(sprintf("%s: %s", case_result$case$label, ttl), side = 3L, outer = TRUE, line = -0.3, cex = 1.1, font = 2L)
}

save_shortlist_grid <- function(case_result, output_path) {
  shortlist <- case_result$shortlist_results
  if (!length(shortlist)) {
    return(invisible(NULL))
  }
  entries <- c(
    list(list(
      display = case_result$case$truth,
      title = "Reference surface",
      subtitle = "target geometry"
    )),
    lapply(shortlist, function(entry) {
      second.metrics <- entry$second_lookahead$metrics[1L, , drop = FALSE]
      phase_f_rank <- case_result$shortlist_metrics$phase_f_rank_sum[[match(entry$method_id, case_result$shortlist_metrics$method_id)]]
      list(
        display = entry$second_lookahead$display,
        title = sprintf("%s\nPhase F shortlist rank %d", entry$method_label, as.integer(phase_f_rank)),
        subtitle = sprintf(
          "2-step sigma %s, rough %s, edge floor %s",
          fmt_num(second.metrics$lookahead_sigma[[1L]], 4L),
          fmt_num(second.metrics$lookahead_roughness[[1L]], 4L),
          fmt_num(second.metrics$lookahead_edge_floor[[1L]], 4L)
        )
      )
    })
  )
  ncol.panel <- 3L
  nrow.panel <- ceiling(length(entries) / ncol.panel)
  grDevices::png(output_path, width = 2400L, height = 820L * nrow.panel, res = 180, bg = "#ffffff", type = "cairo")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(nrow.panel, ncol.panel), mar = c(1.1, 1.1, 3.0, 0.4), oma = c(0, 0, 1.0, 0))
  for (i in seq_len(nrow.panel * ncol.panel)) {
    if (i > length(entries)) {
      graphics::plot.new()
      next
    }
    entry <- entries[[i]]
    plot_partial_layout(
      coords = entry$display,
      edges = case_result$case$edges,
      projection = "ortho",
      azimuth = 35,
      elevation = 24,
      vertex.col = if (i == 1L) "#bc6c25" else "#355070",
      edge.col = "#adb5bd",
      main = ""
    )
    graphics::mtext(sprintf("%s\n%s", entry$title, entry$subtitle), side = 3L, line = 0.3, cex = 0.8)
  }
  graphics::mtext(sprintf("%s: Phase F shortlist second lookahead", case_result$case$label), side = 3L, outer = TRUE, line = -0.3, cex = 1.1, font = 2L)
}

save_portfolio_grid <- function(case_result, output_path) {
  entries <- c(
    list(list(
      display = case_result$case$truth,
      title = "Reference surface",
      subtitle = "target geometry"
    )),
    lapply(case_result$portfolio_results, function(entry) {
      list(
        display = entry$display,
        title = sprintf("%s\nPhase G rank %d", entry$method_label, as.integer(entry$phase_g_rank_sum)),
        subtitle = sprintf(
          "sigma %s, rough %s, edge floor %s, alpha %s",
          fmt_num(entry$final_sigma, 4L),
          fmt_num(entry$final_roughness, 4L),
          fmt_num(entry$final_edge_floor, 4L),
          fmt_num(entry$final_area_q05_ratio, 4L)
        )
      )
    })
  )
  ncol.panel <- 3L
  nrow.panel <- ceiling(length(entries) / ncol.panel)
  grDevices::png(output_path, width = 2400L, height = 820L * nrow.panel, res = 180, bg = "#ffffff", type = "cairo")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(nrow.panel, ncol.panel), mar = c(1.1, 1.1, 3.0, 0.4), oma = c(0, 0, 1.0, 0))
  for (i in seq_len(nrow.panel * ncol.panel)) {
    if (i > length(entries)) {
      graphics::plot.new()
      next
    }
    entry <- entries[[i]]
    plot_partial_layout(
      coords = entry$display,
      edges = case_result$case$edges,
      projection = "ortho",
      azimuth = 35,
      elevation = 24,
      vertex.col = if (i == 1L) "#bc6c25" else "#4a4e69",
      edge.col = "#adb5bd",
      main = ""
    )
    graphics::mtext(sprintf("%s\n%s", entry$title, entry$subtitle), side = 3L, line = 0.3, cex = 0.78)
  }
  graphics::mtext(sprintf("%s: Phase G deferred full-portfolio ranking", case_result$case$label), side = 3L, outer = TRUE, line = -0.3, cex = 1.1, font = 2L)
}

save_final_comparison_grid <- function(case_result, output_path) {
  entries <- case_result$comparison_entries
  n <- length(entries)
  ncol.panel <- 3L
  nrow.panel <- ceiling(n / ncol.panel)
  grDevices::png(output_path, width = 2400L, height = 820L * nrow.panel, res = 180, bg = "#ffffff", type = "cairo")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(nrow.panel, ncol.panel), mar = c(1.1, 1.1, 3.0, 0.4), oma = c(0, 0, 1.0, 0))
  for (i in seq_len(nrow.panel * ncol.panel)) {
    if (i > length(entries)) {
      graphics::plot.new()
      next
    }
    entry <- entries[[i]]
    plot_partial_layout(
      coords = entry$display_coords,
      edges = case_result$case$edges,
      projection = "ortho",
      azimuth = 35,
      elevation = 24,
      vertex.col = if (entry$method_id == "reference") "#bc6c25" else "#355070",
      edge.col = "#adb5bd",
      main = ""
    )
    graphics::mtext(
      sprintf("%s\nsigma %s, rho %s", entry$method_label, fmt_num(entry$gmds_stress, 4L), fmt_num(entry$procrustes_rmse, 4L)),
      side = 3L,
      line = 0.3,
      cex = 0.78
    )
  }
  graphics::mtext(sprintf("%s: final pipeline comparison", case_result$case$label), side = 3L, outer = TRUE, line = -0.3, cex = 1.1, font = 2L)
}

save_selected_stage_grid <- function(case_result, output_path) {
  methods <- case_result$selected_stage_methods
  stage.names <- c("top_level", "after_insertion", "after_refinement", "final")
  grDevices::png(output_path, width = 420L * length(stage.names), height = 320L * length(methods), res = 180, bg = "#ffffff", type = "cairo")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(length(methods), length(stage.names)), mar = c(0.8, 0.8, 2.7, 0.3), oma = c(0, 0, 1.0, 0))
  for (method in methods) {
    for (stage.name in stage.names) {
      coords <- method$stage_display[[stage.name]]
      plot_partial_layout(
        coords = coords,
        edges = case_result$case$edges,
        projection = "ortho",
        azimuth = 35,
        elevation = 24,
        vertex.col = "#355070",
        edge.col = "#c6d1db",
        main = ""
      )
      label <- switch(
        stage.name,
        top_level = "top",
        after_insertion = "after insertion",
        after_refinement = "after refinement",
        final = "final",
        stage.name
      )
      graphics::mtext(sprintf("%s\n%s", method$method_label, label), side = 3L, line = 0.3, cex = 0.72)
    }
  }
  graphics::mtext(sprintf("%s: selected pipeline stage layouts", case_result$case$label), side = 3L, outer = TRUE, line = -0.3, cex = 1.05, font = 2L)
}

write_candidate_table <- function(df) {
  rows <- vapply(seq_len(nrow(df)), function(i) {
    sprintf(
      "%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
      tex_escape(df$method_label[[i]]),
      fmt_num(df$lookahead_sigma[[i]], 4L),
      fmt_num(df$lookahead_rho[[i]], 4L),
      fmt_num(df$lookahead_roughness[[i]], 4L),
      fmt_num(df$lookahead_edge_floor[[i]], 4L),
      as.integer(df$proxy_rank_sum[[i]]),
      fmt_num(df$final_sigma[[i]], 4L),
      fmt_num(df$final_rho[[i]], 4L),
      fmt_time(df$elapsed_sec[[i]])
    )
  }, character(1L))
  paste(
    "\\begin{table}[H]",
    "\\centering",
    "\\scriptsize",
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{lrrrrrrrr}",
    "\\toprule",
    "Candidate & $\\sigma_{look}$ & $\\rho_{look}$ & $\\eta_{look}$ & $q_{0.05}/q_{0.50}$ & proxy rank & $\\sigma$ & $\\rho$ & $t$ (s) \\\\",
    "\\midrule",
    paste(rows, collapse = "\n"),
    "\\bottomrule",
    "\\end{tabular}}",
    "\\end{table}",
    sep = "\n"
  )
}

write_shortlist_table <- function(df) {
  rows <- vapply(seq_len(nrow(df)), function(i) {
    sprintf(
      "%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
      tex_escape(df$method_label[[i]]),
      as.integer(df$phase_e_proxy_rank[[i]]),
      fmt_num(df$second_sigma[[i]], 4L),
      fmt_num(df$second_roughness[[i]], 4L),
      fmt_num(df$second_edge_floor[[i]], 4L),
      fmt_num(df$energy_stability[[i]], 4L),
      as.integer(df$phase_f_rank_sum[[i]]),
      fmt_num(df$final_sigma[[i]], 4L),
      fmt_num(df$final_rho[[i]], 4L)
    )
  }, character(1L))
  paste(
    "\\begin{table}[H]",
    "\\centering",
    "\\scriptsize",
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{lrrrrrrrr}",
    "\\toprule",
    "Shortlist candidate & E-rank & $\\sigma_{2}$ & $\\eta_{2}$ & $q_{0.05}/q_{0.50,2}$ & stab. & F-rank & $\\sigma$ & $\\rho$ \\\\",
    "\\midrule",
    paste(rows, collapse = "\n"),
    "\\bottomrule",
    "\\end{tabular}}",
    "\\end{table}",
    sep = "\n"
  )
}

write_portfolio_table <- function(df) {
  rows <- vapply(seq_len(nrow(df)), function(i) {
    sprintf(
      "%s & %s & %s & %s & %s & %s & %s \\\\",
      tex_escape(df$method_label[[i]]),
      fmt_num(df$final_sigma[[i]], 4L),
      fmt_num(df$final_roughness[[i]], 4L),
      fmt_num(df$final_edge_floor[[i]], 4L),
      fmt_num(df$final_area_q05_ratio[[i]], 4L),
      as.integer(df$phase_g_rank_sum[[i]]),
      fmt_num(df$final_rho[[i]], 4L)
    )
  }, character(1L))
  paste(
    "\\begin{table}[H]",
    "\\centering",
    "\\scriptsize",
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{lrrrrrr}",
    "\\toprule",
    "Portfolio candidate & $\\sigma$ & $\\eta$ & edge floor & $\\alpha_{0.05}$ & G-rank & $\\rho$ \\\\",
    "\\midrule",
    paste(rows, collapse = "\n"),
    "\\bottomrule",
    "\\end{tabular}}",
    "\\end{table}",
    sep = "\n"
  )
}

cases <- if (smoke) {
  list(
    make_case(12L, "flat", "regular"),
    make_case(15L, "paraboloid", "irregular_rectangle")
  )
} else {
  c(
    lapply(c("flat", "saddle", "paraboloid", "ripple"), function(surface) make_case(12L, surface, "regular")),
    lapply(c("flat", "saddle", "paraboloid", "ripple"), function(surface) make_case(15L, surface, "regular")),
    lapply(c("flat", "saddle", "paraboloid", "ripple"), function(surface) make_case(15L, surface, "irregular_rectangle"))
  )
}

case_results <- lapply(cases, function(case) {
  message("Running Phase G deferred portfolio selector study for: ", case$label)

  initializers <- lapply(method_specs, function(spec) {
    message("  initializer: ", spec$id)
    run_top_initializer(case, spec)
  })
  names(initializers) <- vapply(initializers, `[[`, character(1L), "method_id")

  corrected <- lapply(initializers, function(init) run_top_correction(case, init))
  prepared <- case$prepared
  candidates <- lapply(cfg$candidate_ids, function(method_id) {
    message("  candidate: ", method_id)
    run_candidate(case, corrected[[method_id]])
  })
  names(candidates) <- cfg$candidate_ids

  selection_df <- do.call(rbind, lapply(candidates, function(candidate) {
    look <- candidate$lookahead_metrics[1L, , drop = FALSE]
    fin <- candidate$final_metrics[1L, , drop = FALSE]
    data.frame(
      case_id = case$id,
      case_label = case$label,
      family = case$family,
      surface = case$surface,
      side = case$side,
      method_id = candidate$method_id,
      method_label = candidate$method_label,
      lookahead_sigma = look$lookahead_sigma[[1L]],
      lookahead_rho = look$lookahead_rho[[1L]],
      lookahead_roughness = look$lookahead_roughness[[1L]],
      lookahead_edge_floor = look$lookahead_edge_floor[[1L]],
      final_sigma = fin$gmds_stress[[1L]],
      final_rho = fin$procrustes_rmse[[1L]],
      final_roughness = fin$roughness[[1L]],
      final_edge_floor = fin$edge_floor[[1L]],
      final_area_q05_ratio = fin$area_q05_ratio[[1L]],
      elapsed_sec = fin$elapsed_sec[[1L]],
      stringsAsFactors = FALSE
    )
  }))
  selection_df$rank_sigma <- rank(selection_df$lookahead_sigma, ties.method = "min")
  selection_df$rank_roughness <- rank(selection_df$lookahead_roughness, ties.method = "min")
  selection_df$rank_edge_floor <- rank(-selection_df$lookahead_edge_floor, ties.method = "min")
  selection_df$proxy_rank_sum <- selection_df$rank_sigma + selection_df$rank_roughness + selection_df$rank_edge_floor
  selection_df <- selection_df[order(selection_df$proxy_rank_sum, selection_df$lookahead_sigma, selection_df$method_label), , drop = FALSE]
  rownames(selection_df) <- NULL
  selection_df$phase_e_proxy_rank <- seq_len(nrow(selection_df))

  phase_e.method.id <- selection_df$method_id[[1L]]
  oracle.sigma.id <- selection_df$method_id[[which.min(selection_df$final_sigma)]]
  oracle.rho.id <- selection_df$method_id[[which.min(selection_df$final_rho)]]
  fixed.method.id <- fixed_control_seed_id(case)
  shortlist.ids <- selection_df$method_id[seq_len(min(2L, nrow(selection_df)))]

  shortlist.results <- lapply(shortlist.ids, function(method_id) {
    second <- run_additional_lookahead(
      case = case,
      prepared = prepared,
      coords = candidates[[method_id]]$lookahead_coords_full,
      current_level = candidates[[method_id]]$lookahead_level
    )
    first.look <- candidates[[method_id]]$lookahead_metrics[1L, , drop = FALSE]
    if (is.null(second)) {
      second <- list(
        level = candidates[[method_id]]$lookahead_level,
        coords_full = candidates[[method_id]]$lookahead_coords_full,
        display = candidates[[method_id]]$lookahead_display,
        metrics = first.look
      )
    }
    energy.stability <- abs(log((second$metrics$lookahead_energy[[1L]] + 1e-8) / (first.look$lookahead_energy[[1L]] + 1e-8)))
    list(
      method_id = method_id,
      method_label = candidates[[method_id]]$method_label,
      phase_e_proxy_rank = selection_df$phase_e_proxy_rank[[match(method_id, selection_df$method_id)]],
      first_lookahead = first.look,
      second_lookahead = second,
      energy_stability = energy.stability,
      final_metrics = candidates[[method_id]]$final_metrics
    )
  })

  shortlist_df <- do.call(rbind, lapply(shortlist.results, function(entry) {
    second.metrics <- entry$second_lookahead$metrics[1L, , drop = FALSE]
    data.frame(
      case_id = case$id,
      case_label = case$label,
      family = case$family,
      surface = case$surface,
      side = case$side,
      method_id = entry$method_id,
      method_label = entry$method_label,
      phase_e_proxy_rank = entry$phase_e_proxy_rank,
      second_level = entry$second_lookahead$level,
      second_sigma = second.metrics$lookahead_sigma[[1L]],
      second_energy = second.metrics$lookahead_energy[[1L]],
      second_roughness = second.metrics$lookahead_roughness[[1L]],
      second_edge_floor = second.metrics$lookahead_edge_floor[[1L]],
      energy_stability = entry$energy_stability,
      final_sigma = entry$final_metrics$gmds_stress[[1L]],
      final_rho = entry$final_metrics$procrustes_rmse[[1L]],
      stringsAsFactors = FALSE
    )
  }))
  shortlist_df$rank_phase_e <- rank_smallest(shortlist_df$phase_e_proxy_rank)
  shortlist_df$rank_second_sigma <- rank_smallest(shortlist_df$second_sigma)
  shortlist_df$rank_second_roughness <- rank_smallest(shortlist_df$second_roughness)
  shortlist_df$rank_second_edge_floor <- rank_largest(shortlist_df$second_edge_floor)
  shortlist_df$rank_energy_stability <- rank_smallest(shortlist_df$energy_stability)
  shortlist_df$phase_f_rank_sum <- shortlist_df$rank_phase_e +
    shortlist_df$rank_second_sigma +
    shortlist_df$rank_second_roughness +
    shortlist_df$rank_second_edge_floor +
    shortlist_df$rank_energy_stability
  shortlist_df <- shortlist_df[order(shortlist_df$phase_f_rank_sum, shortlist_df$phase_e_proxy_rank, shortlist_df$method_label), , drop = FALSE]
  rownames(shortlist_df) <- NULL
  phase_f.method.id <- shortlist_df$method_id[[1L]]

  portfolio_df <- selection_df[, c(
    "case_id",
    "case_label",
    "family",
    "surface",
    "side",
    "method_id",
    "method_label",
    "final_sigma",
    "final_rho",
    "final_roughness",
    "final_edge_floor",
    "final_area_q05_ratio",
    "elapsed_sec"
  ), drop = FALSE]
  portfolio_df$rank_final_sigma <- rank_smallest(portfolio_df$final_sigma)
  portfolio_df$rank_final_roughness <- rank_smallest(portfolio_df$final_roughness)
  portfolio_df$rank_final_edge_floor <- rank_largest(portfolio_df$final_edge_floor)
  portfolio_df$rank_final_area_q05_ratio <- rank_largest(portfolio_df$final_area_q05_ratio)
  portfolio_df$phase_g_rank_sum <- portfolio_df$rank_final_sigma +
    portfolio_df$rank_final_roughness +
    portfolio_df$rank_final_edge_floor +
    portfolio_df$rank_final_area_q05_ratio
  portfolio_df <- portfolio_df[order(portfolio_df$phase_g_rank_sum, portfolio_df$final_sigma, portfolio_df$method_label), , drop = FALSE]
  rownames(portfolio_df) <- NULL
  phase_g.method.id <- portfolio_df$method_id[[1L]]

  direct.control <- run_direct_control(case)

  selected.ids <- unique(c(phase_f.method.id, phase_g.method.id, fixed.method.id, oracle.rho.id))
  selected.stage.methods <- lapply(selected.ids, function(method_id) {
    stage.display <- run_candidate_with_frames(case, corrected[[method_id]])
    list(method_id = method_id, method_label = corrected[[method_id]]$method_label, stage_display = stage.display)
  })

  candidate_grid_top <- file.path(pdf_dir, sprintf("%s_top_seed_grid.png", case$id))
  candidate_grid_lookahead <- file.path(pdf_dir, sprintf("%s_lookahead_grid.png", case$id))
  shortlist_grid <- file.path(pdf_dir, sprintf("%s_shortlist_grid.png", case$id))
  portfolio_grid <- file.path(pdf_dir, sprintf("%s_portfolio_grid.png", case$id))
  final_grid <- file.path(pdf_dir, sprintf("%s_final_comparison_grid.png", case$id))
  stage_grid <- file.path(pdf_dir, sprintf("%s_selected_stage_grid.png", case$id))

  portfolio_results <- lapply(seq_len(nrow(portfolio_df)), function(i) {
    method_id <- portfolio_df$method_id[[i]]
    list(
      method_id = method_id,
      method_label = portfolio_df$method_label[[i]],
      phase_g_rank_sum = portfolio_df$phase_g_rank_sum[[i]],
      final_sigma = portfolio_df$final_sigma[[i]],
      final_rho = portfolio_df$final_rho[[i]],
      final_roughness = portfolio_df$final_roughness[[i]],
      final_edge_floor = portfolio_df$final_edge_floor[[i]],
      final_area_q05_ratio = portfolio_df$final_area_q05_ratio[[i]],
      display = candidates[[method_id]]$final_display
    )
  })

  comparison_entries <- list(
    list(
      method_id = "reference",
      method_label = "Reference surface",
      display_coords = case$truth,
      gmds_stress = 0,
      procrustes_rmse = 0
    ),
    list(
      method_id = "cmd_pure_gmds",
      method_label = direct.control$metrics$method_label[[1L]],
      display_coords = direct.control$display_coords,
      gmds_stress = direct.control$metrics$gmds_stress[[1L]],
      procrustes_rmse = direct.control$metrics$procrustes_rmse[[1L]]
    ),
    list(
      method_id = fixed.method.id,
      method_label = paste0("Fixed-top W-KK (", corrected[[fixed.method.id]]$method_label, ")"),
      display_coords = candidates[[fixed.method.id]]$final_display,
      gmds_stress = candidates[[fixed.method.id]]$final_metrics$gmds_stress[[1L]],
      procrustes_rmse = candidates[[fixed.method.id]]$final_metrics$procrustes_rmse[[1L]]
    ),
    list(
      method_id = phase_f.method.id,
      method_label = paste0("Phase F shortlist winner (", corrected[[phase_f.method.id]]$method_label, ")"),
      display_coords = candidates[[phase_f.method.id]]$final_display,
      gmds_stress = candidates[[phase_f.method.id]]$final_metrics$gmds_stress[[1L]],
      procrustes_rmse = candidates[[phase_f.method.id]]$final_metrics$procrustes_rmse[[1L]]
    ),
    list(
      method_id = phase_g.method.id,
      method_label = paste0("Phase G deferred portfolio winner (", corrected[[phase_g.method.id]]$method_label, ")"),
      display_coords = candidates[[phase_g.method.id]]$final_display,
      gmds_stress = candidates[[phase_g.method.id]]$final_metrics$gmds_stress[[1L]],
      procrustes_rmse = candidates[[phase_g.method.id]]$final_metrics$procrustes_rmse[[1L]]
    )
  )
  if (!oracle.rho.id %in% vapply(comparison_entries, `[[`, character(1L), "method_id")) {
    comparison_entries[[length(comparison_entries) + 1L]] <- list(
      method_id = oracle.rho.id,
      method_label = paste0("Oracle best final rho (", corrected[[oracle.rho.id]]$method_label, ")"),
      display_coords = candidates[[oracle.rho.id]]$final_display,
      gmds_stress = candidates[[oracle.rho.id]]$final_metrics$gmds_stress[[1L]],
      procrustes_rmse = candidates[[oracle.rho.id]]$final_metrics$procrustes_rmse[[1L]]
    )
  }

  case.out <- list(
    case = case,
    corrected = corrected,
    candidates = candidates,
    selection_metrics = selection_df,
    phase_e_proxy_method_id = phase_e.method.id,
    phase_f_method_id = phase_f.method.id,
    phase_g_method_id = phase_g.method.id,
    shortlist_metrics = shortlist_df,
    shortlist_results = shortlist.results,
    portfolio_metrics = portfolio_df,
    portfolio_results = portfolio_results,
    oracle_sigma_method_id = oracle.sigma.id,
    oracle_rho_method_id = oracle.rho.id,
    fixed_method_id = fixed.method.id,
    direct_control = direct.control,
    comparison_entries = comparison_entries,
    selected_stage_methods = selected.stage.methods
  )

  save_candidate_grid(case.out, stage = "top", output_path = candidate_grid_top)
  save_candidate_grid(case.out, stage = "lookahead", output_path = candidate_grid_lookahead)
  save_shortlist_grid(case.out, output_path = shortlist_grid)
  save_portfolio_grid(case.out, output_path = portfolio_grid)
  save_final_comparison_grid(case.out, output_path = final_grid)
  save_selected_stage_grid(case.out, output_path = stage_grid)

  case.out
})

selection_df <- do.call(rbind, lapply(case_results, `[[`, "selection_metrics"))
shortlist_df_all <- do.call(rbind, lapply(case_results, `[[`, "shortlist_metrics"))
portfolio_df_all <- do.call(rbind, lapply(case_results, `[[`, "portfolio_metrics"))

final_rows <- do.call(rbind, lapply(case_results, function(case_result) {
  data.frame(
    case_id = case_result$case$id,
    case_label = case_result$case$label,
    family = case_result$case$family,
    surface = case_result$case$surface,
    side = case_result$case$side,
    baseline = c("direct_cmd_pure_gmds", "fixed_top_weighted_kk", "phase_e_proxy", "phase_f_shortlist", "phase_g_deferred_portfolio", "oracle_best_rho"),
    method_id = c(
      "cmd_pure_gmds",
      case_result$fixed_method_id,
      case_result$phase_e_proxy_method_id,
      case_result$phase_f_method_id,
      case_result$phase_g_method_id,
      case_result$oracle_rho_method_id
    ),
    method_label = c(
      case_result$direct_control$metrics$method_label[[1L]],
      paste0("Fixed-top W-KK (", case_result$fixed_method_id, ")"),
      paste0("Phase E proxy winner (", case_result$phase_e_proxy_method_id, ")"),
      paste0("Phase F shortlist winner (", case_result$phase_f_method_id, ")"),
      paste0("Phase G deferred portfolio winner (", case_result$phase_g_method_id, ")"),
      paste0("Oracle best final rho (", case_result$oracle_rho_method_id, ")")
    ),
    gmds_stress = c(
      case_result$direct_control$metrics$gmds_stress[[1L]],
      case_result$candidates[[case_result$fixed_method_id]]$final_metrics$gmds_stress[[1L]],
      case_result$candidates[[case_result$phase_e_proxy_method_id]]$final_metrics$gmds_stress[[1L]],
      case_result$candidates[[case_result$phase_f_method_id]]$final_metrics$gmds_stress[[1L]],
      case_result$candidates[[case_result$phase_g_method_id]]$final_metrics$gmds_stress[[1L]],
      case_result$candidates[[case_result$oracle_rho_method_id]]$final_metrics$gmds_stress[[1L]]
    ),
    procrustes_rmse = c(
      case_result$direct_control$metrics$procrustes_rmse[[1L]],
      case_result$candidates[[case_result$fixed_method_id]]$final_metrics$procrustes_rmse[[1L]],
      case_result$candidates[[case_result$phase_e_proxy_method_id]]$final_metrics$procrustes_rmse[[1L]],
      case_result$candidates[[case_result$phase_f_method_id]]$final_metrics$procrustes_rmse[[1L]],
      case_result$candidates[[case_result$phase_g_method_id]]$final_metrics$procrustes_rmse[[1L]],
      case_result$candidates[[case_result$oracle_rho_method_id]]$final_metrics$procrustes_rmse[[1L]]
    ),
    stringsAsFactors = FALSE
  )
}))

summary_rows <- do.call(rbind, lapply(case_results, function(case_result) {
  data.frame(
    case_id = case_result$case$id,
    case_label = case_result$case$label,
    family = case_result$case$family,
    surface = case_result$case$surface,
    side = case_result$case$side,
    phase_e_proxy_method_id = case_result$phase_e_proxy_method_id,
    phase_f_method_id = case_result$phase_f_method_id,
    phase_g_method_id = case_result$phase_g_method_id,
    fixed_method_id = case_result$fixed_method_id,
    oracle_sigma_method_id = case_result$oracle_sigma_method_id,
    oracle_rho_method_id = case_result$oracle_rho_method_id,
    phase_e_matches_oracle_rho = identical(case_result$phase_e_proxy_method_id, case_result$oracle_rho_method_id),
    phase_f_matches_oracle_rho = identical(case_result$phase_f_method_id, case_result$oracle_rho_method_id),
    phase_f_matches_oracle_sigma = identical(case_result$phase_f_method_id, case_result$oracle_sigma_method_id),
    phase_g_matches_oracle_rho = identical(case_result$phase_g_method_id, case_result$oracle_rho_method_id),
    phase_g_matches_oracle_sigma = identical(case_result$phase_g_method_id, case_result$oracle_sigma_method_id),
    phase_f_beats_phase_e_rho =
      case_result$candidates[[case_result$phase_f_method_id]]$final_metrics$procrustes_rmse[[1L]] <
      case_result$candidates[[case_result$phase_e_proxy_method_id]]$final_metrics$procrustes_rmse[[1L]],
    phase_f_beats_phase_e_sigma =
      case_result$candidates[[case_result$phase_f_method_id]]$final_metrics$gmds_stress[[1L]] <
      case_result$candidates[[case_result$phase_e_proxy_method_id]]$final_metrics$gmds_stress[[1L]],
    phase_g_beats_phase_e_rho =
      case_result$candidates[[case_result$phase_g_method_id]]$final_metrics$procrustes_rmse[[1L]] <
      case_result$candidates[[case_result$phase_e_proxy_method_id]]$final_metrics$procrustes_rmse[[1L]],
    phase_g_beats_phase_e_sigma =
      case_result$candidates[[case_result$phase_g_method_id]]$final_metrics$gmds_stress[[1L]] <
      case_result$candidates[[case_result$phase_e_proxy_method_id]]$final_metrics$gmds_stress[[1L]],
    phase_g_beats_phase_f_rho =
      case_result$candidates[[case_result$phase_g_method_id]]$final_metrics$procrustes_rmse[[1L]] <
      case_result$candidates[[case_result$phase_f_method_id]]$final_metrics$procrustes_rmse[[1L]],
    phase_g_beats_phase_f_sigma =
      case_result$candidates[[case_result$phase_g_method_id]]$final_metrics$gmds_stress[[1L]] <
      case_result$candidates[[case_result$phase_f_method_id]]$final_metrics$gmds_stress[[1L]],
    phase_f_beats_fixed_rho =
      case_result$candidates[[case_result$phase_f_method_id]]$final_metrics$procrustes_rmse[[1L]] <
      case_result$candidates[[case_result$fixed_method_id]]$final_metrics$procrustes_rmse[[1L]],
    phase_f_beats_fixed_sigma =
      case_result$candidates[[case_result$phase_f_method_id]]$final_metrics$gmds_stress[[1L]] <
      case_result$candidates[[case_result$fixed_method_id]]$final_metrics$gmds_stress[[1L]],
    phase_g_beats_fixed_rho =
      case_result$candidates[[case_result$phase_g_method_id]]$final_metrics$procrustes_rmse[[1L]] <
      case_result$candidates[[case_result$fixed_method_id]]$final_metrics$procrustes_rmse[[1L]],
    phase_g_beats_fixed_sigma =
      case_result$candidates[[case_result$phase_g_method_id]]$final_metrics$gmds_stress[[1L]] <
      case_result$candidates[[case_result$fixed_method_id]]$final_metrics$gmds_stress[[1L]],
    stringsAsFactors = FALSE
  )
}))

utils::write.csv(selection_df, selection_csv, row.names = FALSE)
utils::write.csv(final_rows, final_csv, row.names = FALSE)
utils::write.csv(summary_rows, summary_csv, row.names = FALSE)
utils::write.csv(shortlist_df_all, shortlist_csv, row.names = FALSE)
utils::write.csv(portfolio_df_all, portfolio_csv, row.names = FALSE)

compact_case_result <- function(case_result) {
  list(
    case = list(
      id = case_result$case$id,
      label = case_result$case$label,
      family = case_result$case$family,
      surface = case_result$case$surface,
      side = case_result$case$side,
      n = case_result$case$n,
      edges = case_result$case$edges,
      truth = case_result$case$truth
    ),
    candidates = lapply(case_result$candidates, function(candidate) {
      list(
        method_id = candidate$method_id,
        method_label = candidate$method_label,
        top_display = candidate$top_display,
        lookahead_display = candidate$lookahead_display,
        lookahead_level = candidate$lookahead_level,
        final_display = candidate$final_display,
        lookahead_metrics = candidate$lookahead_metrics,
        final_metrics = candidate$final_metrics
      )
    }),
    selection_metrics = case_result$selection_metrics,
    phase_e_proxy_method_id = case_result$phase_e_proxy_method_id,
    phase_f_method_id = case_result$phase_f_method_id,
    phase_g_method_id = case_result$phase_g_method_id,
    shortlist_metrics = case_result$shortlist_metrics,
    shortlist_results = lapply(case_result$shortlist_results, function(entry) {
      list(
        method_id = entry$method_id,
        method_label = entry$method_label,
        phase_e_proxy_rank = entry$phase_e_proxy_rank,
        energy_stability = entry$energy_stability,
        second_lookahead = list(
          level = entry$second_lookahead$level,
          display = entry$second_lookahead$display,
          metrics = entry$second_lookahead$metrics
        )
      )
    }),
    portfolio_metrics = case_result$portfolio_metrics,
    portfolio_results = case_result$portfolio_results,
    oracle_sigma_method_id = case_result$oracle_sigma_method_id,
    oracle_rho_method_id = case_result$oracle_rho_method_id,
    fixed_method_id = case_result$fixed_method_id,
    direct_control = list(
      display_coords = case_result$direct_control$display_coords,
      metrics = case_result$direct_control$metrics
    ),
    comparison_entries = case_result$comparison_entries,
    selected_stage_methods = case_result$selected_stage_methods
  )
}

saveRDS(
  list(
    run_tag = run_tag,
    cfg = cfg,
    case_results = lapply(case_results, compact_case_result),
    selection_metrics = selection_df,
    shortlist_metrics = shortlist_df_all,
    portfolio_metrics = portfolio_df_all,
    final_metrics = final_rows,
    summary_metrics = summary_rows
  ),
  rds_path
)

summary_by_surface <- do.call(rbind, lapply(split(summary_rows, summary_rows$surface), function(df) {
  data.frame(
    surface = df$surface[[1L]],
    n_cases = nrow(df),
    phase_e_matches_oracle_rho = sum(df$phase_e_matches_oracle_rho),
    phase_f_matches_oracle_rho = sum(df$phase_f_matches_oracle_rho),
    phase_g_matches_oracle_rho = sum(df$phase_g_matches_oracle_rho),
    phase_g_beats_phase_f_rho = sum(df$phase_g_beats_phase_f_rho),
    phase_g_beats_fixed_rho = sum(df$phase_g_beats_fixed_rho),
    stringsAsFactors = FALSE
  )
}))

summary_lines <- paste(
  "\\begin{table}[H]",
  "\\centering",
  "\\small",
  "\\begin{tabular}{lrrrrr}",
  "\\toprule",
  "Surface & E = oracle $\\rho$ & F = oracle $\\rho$ & G = oracle $\\rho$ & G beats F $\\rho$ & G beats fixed $\\rho$ \\\\",
  "\\midrule",
  paste(
    vapply(seq_len(nrow(summary_by_surface)), function(i) {
      sprintf(
        "%s & %d/%d & %d/%d & %d/%d & %d/%d & %d/%d \\\\",
        tex_escape(surface_title(summary_by_surface$surface[[i]])),
        summary_by_surface$phase_e_matches_oracle_rho[[i]],
        summary_by_surface$n_cases[[i]],
        summary_by_surface$phase_f_matches_oracle_rho[[i]],
        summary_by_surface$n_cases[[i]],
        summary_by_surface$phase_g_matches_oracle_rho[[i]],
        summary_by_surface$n_cases[[i]],
        summary_by_surface$phase_g_beats_phase_f_rho[[i]],
        summary_by_surface$n_cases[[i]],
        summary_by_surface$phase_g_beats_fixed_rho[[i]],
        summary_by_surface$n_cases[[i]]
      )
    }, character(1L)),
    collapse = "\n"
  ),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  sep = "\n"
)

case_sections <- vapply(case_results, function(case_result) {
  candidate.table <- case_result$selection_metrics[order(case_result$selection_metrics$proxy_rank_sum, case_result$selection_metrics$final_rho), , drop = FALSE]
  portfolio.table <- case_result$portfolio_metrics[order(case_result$portfolio_metrics$phase_g_rank_sum, case_result$portfolio_metrics$final_rho), , drop = FALSE]
  phase_e.row <- candidate.table[candidate.table$method_id == case_result$phase_e_proxy_method_id, , drop = FALSE]
  phase_f.row <- candidate.table[candidate.table$method_id == case_result$phase_f_method_id, , drop = FALSE]
  phase_g.row <- candidate.table[candidate.table$method_id == case_result$phase_g_method_id, , drop = FALSE]
  fixed.row <- candidate.table[candidate.table$method_id == case_result$fixed_method_id, , drop = FALSE]
  oracle.row <- candidate.table[candidate.table$method_id == case_result$oracle_rho_method_id, , drop = FALSE]
  shortlist.table <- case_result$shortlist_metrics[order(case_result$shortlist_metrics$phase_f_rank_sum, case_result$shortlist_metrics$final_rho), , drop = FALSE]
  direct.row <- case_result$direct_control$metrics
  paste(
    sprintf("\\section*{%s}", tex_escape(case_result$case$label)),
    sprintf(
      paste(
        "This Phase~G case keeps the same seeded MISF-GMDS pipeline and compares three selector layers:",
        "the original Phase~E one-level proxy, the Phase~F shortlist second lookahead, and a deferred full-portfolio Phase~G rank.",
        "The fixed-top weighted-KK control is \\textbf{%s};",
        "the Phase~E proxy winner is \\textbf{%s};",
        "the Phase~F shortlist winner is \\textbf{%s};",
        "the Phase~G deferred portfolio winner is \\textbf{%s};",
        "the oracle best final $\\\\rho$ seed is \\textbf{%s}.",
        "The direct full-graph cMDS control finishes at $\\\\sigma=%s$, $\\\\rho=%s$.",
        "The fixed-top control finishes at $\\\\sigma=%s$, $\\\\rho=%s$.",
        "The Phase~E winner finishes at $\\\\sigma=%s$, $\\\\rho=%s$.",
        "The Phase~F winner finishes at $\\\\sigma=%s$, $\\\\rho=%s$.",
        "The Phase~G winner finishes at $\\\\sigma=%s$, $\\\\rho=%s$."
      ),
      tex_escape(fixed.row$method_label[[1L]]),
      tex_escape(phase_e.row$method_label[[1L]]),
      tex_escape(phase_f.row$method_label[[1L]]),
      tex_escape(phase_g.row$method_label[[1L]]),
      tex_escape(oracle.row$method_label[[1L]]),
      fmt_num(direct.row$gmds_stress[[1L]], 4L),
      fmt_num(direct.row$procrustes_rmse[[1L]], 4L),
      fmt_num(fixed.row$final_sigma[[1L]], 4L),
      fmt_num(fixed.row$final_rho[[1L]], 4L),
      fmt_num(phase_e.row$final_sigma[[1L]], 4L),
      fmt_num(phase_e.row$final_rho[[1L]], 4L),
      fmt_num(phase_f.row$final_sigma[[1L]], 4L),
      fmt_num(phase_f.row$final_rho[[1L]], 4L),
      fmt_num(phase_g.row$final_sigma[[1L]], 4L),
      fmt_num(phase_g.row$final_rho[[1L]], 4L)
    ),
    sprintf(
      "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{%s. Corrected top-level seed candidates.}\\end{figure}",
      file.path(pdf_dir, sprintf("%s_top_seed_grid.png", case_result$case$id)),
      tex_escape(case_result$case$label)
    ),
    sprintf(
      "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{%s. One-level weighted-KK lookahead used for the original Phase~E proxy ranking.}\\end{figure}",
      file.path(pdf_dir, sprintf("%s_lookahead_grid.png", case_result$case$id)),
      tex_escape(case_result$case$label)
    ),
    write_candidate_table(candidate.table),
    sprintf(
      "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{%s. Phase~F shortlist second lookahead for the top-two Phase~E candidates.}\\end{figure}",
      file.path(pdf_dir, sprintf("%s_shortlist_grid.png", case_result$case$id)),
      tex_escape(case_result$case$label)
    ),
    write_shortlist_table(shortlist.table),
    sprintf(
      "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{%s. Deferred full-portfolio ranking over the five corrected top-seed candidates.}\\end{figure}",
      file.path(pdf_dir, sprintf("%s_portfolio_grid.png", case_result$case$id)),
      tex_escape(case_result$case$label)
    ),
    write_portfolio_table(portfolio.table),
    sprintf(
      "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{%s. Final comparison: direct cmdscale control, fixed-top weighted-KK control, Phase~F shortlist winner, Phase~G deferred portfolio winner, and oracle best final $\\rho$.}\\end{figure}",
      file.path(pdf_dir, sprintf("%s_final_comparison_grid.png", case_result$case$id)),
      tex_escape(case_result$case$label)
    ),
    sprintf(
      "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{%s. Selected stage layouts for the fixed, Phase~F, Phase~G, and oracle seeded pipelines.}\\end{figure}",
      file.path(pdf_dir, sprintf("%s_selected_stage_grid.png", case_result$case$id)),
      tex_escape(case_result$case$label)
    ),
    sep = "\n\n"
  )
}, character(1L))

tex_lines <- c(
  "\\documentclass[11pt]{article}",
  "\\usepackage[margin=1in]{geometry}",
  "\\usepackage{graphicx}",
  "\\usepackage{booktabs}",
  "\\usepackage{float}",
  "\\usepackage{amsmath}",
  "\\usepackage[T1]{fontenc}",
  "\\usepackage[utf8]{inputenc}",
  "\\title{Phase G: Deferred Portfolio Selection for Seeded MISF-GMDS}",
  "\\author{Selector-improvement study on the mesh-geometry panel}",
  "\\date{2026-04-05}",
  "\\begin{document}",
  "\\maketitle",
  "\\section{Purpose}",
  "Phase~G asks whether seeded MISF-GMDS should stop choosing the top seed early. The same rectangular mesh-geometry panel is reused: regular $12\\times 12$, regular $15\\times 15$, and irregular rectangle $15\\times 15$, each with flat, saddle, paraboloid, and ripple geometries.",
  "\\section{Protocol}",
  "For each case we compare the same five corrected top-level seeds (cMDS, KK, GRIP, Weighted GRIP, and Weighted GRIP + polish LGKK). Each seed receives the same short top-level anchor-regularized GMDS correction, the same weighted-KK lower-level placement, the same sparse GMDS refinement, and the same final pure-GMDS polish. Phase~G keeps the earlier Phase~E and Phase~F selectors as explicit baselines, but defers the new decision until the full portfolio has finished. The deferred portfolio is then ranked by final full-graph GMDS stress, final roughness, final edge-floor, and the final $\\alpha_{0.05}$ area-floor ratio.",
  "\\section{Block Summary}",
  "The table below aggregates how often the original Phase~E selector, the Phase~F shortlist selector, and the new Phase~G deferred portfolio selector match the oracle best final $\\rho$ seed, and how often the Phase~G selector improves over the Phase~F selector and the current fixed-top weighted-KK control on final $\\rho$ within each geometry class.",
  summary_lines,
  sprintf(
    "The companion HTML gallery is generated by \\\\texttt{%s}.",
    tex_escape("tools/reports/geodesic_mds_paper/render-gmds-misf-phase-g-deferred-portfolio-html.R")
  ),
  unname(case_sections),
  "\\end{document}"
)

writeLines(tex_lines, tex_path)

message("Wrote Phase G selection metrics: ", selection_csv)
message("Wrote Phase G shortlist metrics: ", shortlist_csv)
message("Wrote Phase G deferred portfolio metrics: ", portfolio_csv)
message("Wrote Phase G final metrics: ", final_csv)
message("Wrote Phase G summary metrics: ", summary_csv)
message("Wrote Phase G bundle: ", rds_path)
message("Wrote Phase G report TeX: ", tex_path)
