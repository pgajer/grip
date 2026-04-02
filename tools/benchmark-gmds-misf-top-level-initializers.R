#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) == 0L) {
  stop("Unable to determine script path from commandArgs().")
}
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

args <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% args

run_tag <- if (smoke) {
  sprintf("gmds-misf-top-level-initializers-smoke-%s", format(Sys.Date(), "%Y-%m-%d"))
} else {
  "gmds-misf-top-level-initializers-2026-04-02"
}

manual_root <- file.path(repo_root, "dev", "manual")
tmp_dir <- file.path(manual_root, "tmp", run_tag)
pdf_dir <- file.path(manual_root, "pdf", run_tag)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

tex_path <- file.path(manual_root, "pdf", "gmds_misf_top_level_initializer_report_2026-04-02.tex")
pdf_path <- file.path(manual_root, "pdf", "gmds_misf_top_level_initializer_report_2026-04-02.pdf")
rds_path <- file.path(tmp_dir, "gmds_misf_top_level_initializer_results.rds")
metrics_csv <- file.path(tmp_dir, "gmds_misf_top_level_initializer_metrics.csv")

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
classical_mds_embedding <- get("grip.classical.mds.embedding", envir = ns)

cfg <- list(
  dim = 3L,
  seed = 20260402L,
  trace_frames = if (smoke) 4L else 6L,
  fr_niter = if (smoke) 300L else 800L,
  edge_relax_max_iter = if (smoke) 12L else 24L,
  weighted_core_lgkk_rounds = if (smoke) 1L else 3L,
  weighted_polish_lgkk_rounds = if (smoke) 3L else 8L,
  lgkk_local_nbrs = if (smoke) 6L else 10L,
  lgkk_landmark_count = if (smoke) 6L else 10L,
  lgkk_active_limit = if (smoke) 256L else 4096L
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

build_shared_start <- function(n, dim, seed) {
  set.seed(as.integer(seed))
  coords <- matrix(stats::rnorm(as.integer(n) * as.integer(dim)), ncol = as.integer(dim))
  storage.mode(coords) <- "double"
  normalize_coords(coords)
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

sample_roughness <- function(coords, adj_list, edges) {
  coords <- as.matrix(coords)
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
    if (length(nbrs) == 0L) {
      return(0)
    }
    delta <- centered[i, ] - colMeans(centered[nbrs, , drop = FALSE])
    sum(delta^2)
  }, numeric(1L))
  sqrt(mean(residuals)) / median.edge
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

make_surrogate_prepared <- function(prepared_full) {
  prepared <- prepared_full
  prepared$pair_matrix <- matrix(integer(), ncol = 2L)
  prepared$pair_graph_distance <- numeric(0L)
  prepared$path_vertices <- list()
  prepared$path_edges <- list()
  prepared$path_edge_weights <- list()
  prepared$pair_path_count_log <- numeric(0L)
  prepared$flat_pair_edge_offsets <- as.integer(0L)
  prepared$flat_edge_u <- integer(0L)
  prepared$flat_edge_v <- integer(0L)
  prepared$flat_edge_coeff <- numeric(0L)
  prepared
}

resolve_unweighted_grip_args <- function(family) {
  if (identical(family, "irregular_rectangle")) {
    return(list(
      placement = "barycenter",
      rounds = 192L,
      final_rounds = 256L,
      num_init = 18L,
      num_nbrs = 24L,
      r = 0.05,
      s = 6.5,
      repulsion_factor = 1.10
    ))
  }
  list(
    placement = "barycenter",
    rounds = 128L,
    final_rounds = 128L,
    num_init = 12L,
    num_nbrs = 20L,
    r = 0.10,
    s = 4.5,
    repulsion_factor = 1.50
  )
}

select_trace_indices <- function(meta, n_frames, max_frames = 6L) {
  if (n_frames <= 0L) {
    return(integer(0L))
  }
  max_frames <- max(1L, as.integer(max_frames))
  base <- unique(round(seq(1, n_frames, length.out = min(max_frames, n_frames))))
  if (!is.null(meta) && nrow(meta) >= 1L && "phase" %in% names(meta)) {
    phase <- as.character(meta$phase)
    phase_breaks <- c(1L, which(phase[-1L] != phase[-length(phase)]) + 1L, n_frames)
    base <- sort(unique(c(base, phase_breaks)))
  }
  if (length(base) <= max_frames) {
    return(as.integer(base))
  }
  keep <- unique(round(seq(1, length(base), length.out = max_frames)))
  as.integer(base[keep])
}

extract_trace_frames <- function(trace_obj, target, max_frames = 6L) {
  if (is.null(trace_obj$frames) || !length(trace_obj$frames)) {
    return(list())
  }
  idx <- select_trace_indices(trace_obj$meta, length(trace_obj$frames), max_frames = max_frames)
  lapply(idx, function(i) {
    frame.coords <- as.matrix(trace_obj$frames[[i]])
    keep <- stats::complete.cases(frame.coords)
    if (sum(keep) >= 3L) {
      aligned.partial <- align_to_target_nd(frame.coords[keep, , drop = FALSE], target[keep, , drop = FALSE], allow.reflection = TRUE)
      display.coords <- matrix(NA_real_, nrow = nrow(frame.coords), ncol = ncol(frame.coords))
      display.coords[keep, ] <- aligned.partial$aligned
    } else {
      display.coords <- frame.coords
    }
    row <- if (!is.null(trace_obj$meta) && nrow(trace_obj$meta) >= i) {
      trace_obj$meta[i, , drop = FALSE]
    } else {
      NULL
    }
    list(
      frame_index = as.integer(i),
      total_frames = as.integer(length(trace_obj$frames)),
      phase = if (is.null(row) || !"phase" %in% names(row)) NA_character_ else as.character(row$phase[[1L]]),
      round_in_level = if (is.null(row) || !"round_in_level" %in% names(row)) NA_integer_ else as.integer(row$round_in_level[[1L]]),
      level_index = if (is.null(row) || !"level_index" %in% names(row)) NA_integer_ else as.integer(row$level_index[[1L]]),
      display_coords = display.coords
    )
  })
}

compute_metrics <- function(case,
                            method_id,
                            method_label,
                            coords,
                            elapsed_sec = NA_real_,
                            note = NULL) {
  if (is.null(coords) || !is.matrix(coords) || any(!is.finite(coords))) {
    return(data.frame(
      case_id = case$id,
      case_label = case$label,
      family = case$family,
      top_level = case$top_level,
      top_n = case$top_n,
      method_id = method_id,
      method_label = method_label,
      elapsed_sec = as.double(elapsed_sec),
      gmds_energy = NA_real_,
      sigma_geo = NA_real_,
      rho = NA_real_,
      eta = NA_real_,
      alpha_0_05 = NA_real_,
      note = if (is.null(note)) "non-finite coordinates" else as.character(note),
      stringsAsFactors = FALSE
    ))
  }
  score <- grip.score.geodesic.mds(coords = coords, prepared = case$top_prepared)
  aligned <- if (identical(method_id, "reference")) {
    list(aligned = as.matrix(coords), rmse = 0)
  } else {
    align_to_target_nd(coords, case$truth, allow.reflection = TRUE)
  }
  display.coords <- aligned$aligned
  data.frame(
    case_id = case$id,
    case_label = case$label,
    family = case$family,
    top_level = case$top_level,
    top_n = case$top_n,
    method_id = method_id,
    method_label = method_label,
    elapsed_sec = as.double(elapsed_sec),
    gmds_energy = as.double(score$gmds.energy[[1L]]),
    sigma_geo = as.double(score$gmds.stress[[1L]]),
    rho = as.double(aligned$rmse),
    eta = as.double(sample_roughness(display.coords, case$display_adj, case$display_edges)),
    alpha_0_05 = as.double(area_floor_ratio(display.coords, case$display_triangles)),
    note = if (is.null(note)) "" else as.character(note),
    stringsAsFactors = FALSE
  )
}

make_case_regular <- function(side, amplitude = 0.35) {
  bundle <- mesh.surface.graph(
    side,
    side,
    surface = "paraboloid",
    amplitude = amplitude,
    connectivity = "orthogonal",
    normalize = "median"
  )
  prepared.misf <- grip.prepare.misf.geodesic.mds(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average",
    dim = cfg$dim,
    top_level_mode = "skip",
    seed = cfg$seed + side
  )
  top <- grip.geodesic.misf.induced_level_graph(prepared.misf)
  top.prepared <- grip.prepare.graph.geodesic.mds(
    edges = top$edges,
    n = top$n,
    edge_weights = top$edge_weights,
    tie_mode = "average"
  )
  layout.graph <- build_sparse_layout_graph(top$distance_matrix, k = 6L)
  layout.prepared <- grip.prepare.graph.geodesic.mds(
    edges = layout.graph$edges,
    n = top$n,
    edge_weights = layout.graph$edge_weights,
    tie_mode = "average"
  )
  param <- as.matrix(bundle$coords_param[top$vertex_ids, , drop = FALSE])
  truth <- as.matrix(bundle$coords_surface[top$vertex_ids, , drop = FALSE])
  display <- build_display_mesh(param)
  list(
    id = sprintf("paraboloid_regular_%dx%d_top", side, side),
    label = sprintf("Regular paraboloid mesh %dx%d: top MISF level", side, side),
    family = "regular",
    side = as.integer(side),
    weighted_preset = "mesh",
    unweighted_args = resolve_unweighted_grip_args("regular"),
    top_level = as.integer(top$level),
    top_n = as.integer(top$n),
    top_vertex_ids = as.integer(top$vertex_ids),
    top_graph = top,
    top_prepared = top.prepared,
    layout_graph = layout.graph,
    layout_prepared = layout.prepared,
    truth = truth,
    param = param,
    display_triangles = display$triangles,
    display_edges = display$edges,
    display_adj = display$adj_list,
    shared_start = build_shared_start(top$n, cfg$dim, cfg$seed + side * 11L)
  )
}

make_case_irregular_rectangle <- function(side, amplitude = 0.35) {
  bundle <- irregular.rectangle.surface.graph(
    side,
    side,
    surface = "paraboloid",
    amplitude = amplitude,
    connectivity = "orthogonal",
    normalize = "median"
  )
  prepared.misf <- grip.prepare.misf.geodesic.mds(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average",
    dim = cfg$dim,
    top_level_mode = "skip",
    seed = cfg$seed + side + 100L
  )
  top <- grip.geodesic.misf.induced_level_graph(prepared.misf)
  top.prepared <- grip.prepare.graph.geodesic.mds(
    edges = top$edges,
    n = top$n,
    edge_weights = top$edge_weights,
    tie_mode = "average"
  )
  layout.graph <- build_sparse_layout_graph(top$distance_matrix, k = 6L)
  layout.prepared <- grip.prepare.graph.geodesic.mds(
    edges = layout.graph$edges,
    n = top$n,
    edge_weights = layout.graph$edge_weights,
    tie_mode = "average"
  )
  param <- as.matrix(bundle$coords_param[top$vertex_ids, , drop = FALSE])
  truth <- as.matrix(bundle$coords_surface[top$vertex_ids, , drop = FALSE])
  display <- build_display_mesh(param)
  list(
    id = sprintf("paraboloid_irregular_rectangle_%dx%d_top", side, side),
    label = sprintf("Irregular rectangle paraboloid %dx%d: top MISF level", side, side),
    family = "irregular_rectangle",
    side = as.integer(side),
    weighted_preset = "irregular",
    unweighted_args = resolve_unweighted_grip_args("irregular_rectangle"),
    top_level = as.integer(top$level),
    top_n = as.integer(top$n),
    top_vertex_ids = as.integer(top$vertex_ids),
    top_graph = top,
    top_prepared = top.prepared,
    layout_graph = layout.graph,
    layout_prepared = layout.prepared,
    truth = truth,
    param = param,
    display_triangles = display$triangles,
    display_edges = display$edges,
    display_adj = display$adj_list,
    shared_start = build_shared_start(top$n, cfg$dim, cfg$seed + side * 17L)
  )
}

method_specs <- list(
  list(id = "cmdscale", label = "cMDS"),
  list(id = "kk", label = "KK"),
  list(id = "fr", label = "FR"),
  list(id = "grip", label = "GRIP"),
  list(id = "weighted_grip", label = "Weighted GRIP"),
  list(id = "weighted_grip_core_lgkk", label = "Weighted GRIP + core LGKK"),
  list(id = "weighted_grip_polish_lgkk", label = "Weighted GRIP + polish LGKK"),
  list(id = "edge_relaxation", label = "Edge relaxation surrogate")
)

run_method <- function(case, method_spec) {
  method_id <- method_spec$id
  method_label <- method_spec$label
  top_graph <- case$top_graph
  top_prepared <- case$top_prepared
  layout_graph <- case$layout_graph
  graph_obj <- igraph::graph_from_edgelist(layout_graph$edges, directed = FALSE)
  trace_selected <- list()
  note <- NULL

  if (identical(method_id, "cmdscale")) {
    started <- proc.time()[["elapsed"]]
    fit <- classical_mds_embedding(top_prepared, dim = cfg$dim, eig = TRUE)
    elapsed <- proc.time()[["elapsed"]] - started
    coords <- fit$coords
  } else if (identical(method_id, "kk")) {
    started <- proc.time()[["elapsed"]]
    coords <- igraph::layout_with_kk(
      graph_obj,
      coords = case$shared_start,
      dim = cfg$dim,
      weights = layout_graph$edge_weights
    )
    elapsed <- proc.time()[["elapsed"]] - started
    note <- "Shared generic random-normal start"
  } else if (identical(method_id, "fr")) {
    started <- proc.time()[["elapsed"]]
    coords <- igraph::layout_with_fr(
      graph_obj,
      coords = case$shared_start,
      dim = cfg$dim,
      weights = layout_graph$edge_weights,
      niter = cfg$fr_niter
    )
    elapsed <- proc.time()[["elapsed"]] - started
    note <- "Shared generic random-normal start"
  } else if (identical(method_id, "grip")) {
    args <- case$unweighted_args
    started <- proc.time()[["elapsed"]]
    tr <- grip.layout.trace(
      edges = layout_graph$edges,
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
      seed = cfg$seed + 301L + case$side
    )
    elapsed <- proc.time()[["elapsed"]] - started
    coords <- tr$final
    trace_selected <- extract_trace_frames(tr, target = case$truth, max_frames = cfg$trace_frames)
    note <- sprintf("Trace snapshots: %d selected of %d total frames", length(trace_selected), length(tr$frames))
  } else if (identical(method_id, "weighted_grip")) {
    started <- proc.time()[["elapsed"]]
    tr <- grip.layout.trace.weighted(
      edges = layout_graph$edges,
      edge_weights = layout_graph$edge_weights,
      n = case$top_n,
      dim = cfg$dim,
      preset = case$weighted_preset,
      seed = cfg$seed + 401L + case$side
    )
    elapsed <- proc.time()[["elapsed"]] - started
    coords <- tr$final
    trace_selected <- extract_trace_frames(tr, target = case$truth, max_frames = cfg$trace_frames)
    note <- sprintf("Trace snapshots: %d selected of %d total frames", length(trace_selected), length(tr$frames))
  } else if (identical(method_id, "weighted_grip_core_lgkk")) {
    started <- proc.time()[["elapsed"]]
    tr <- grip.layout.trace.weighted(
      edges = layout_graph$edges,
      edge_weights = layout_graph$edge_weights,
      n = case$top_n,
      dim = cfg$dim,
      preset = case$weighted_preset,
      lgkk_multiscale_rounds = cfg$weighted_core_lgkk_rounds,
      lgkk_local_nbrs = min(cfg$lgkk_local_nbrs, max(2L, top_graph$n - 1L)),
      lgkk_landmark_count = min(cfg$lgkk_landmark_count, max(2L, top_graph$n - 1L)),
      lgkk_multiscale_scope = "all",
      lgkk_active_limit = cfg$lgkk_active_limit,
      seed = cfg$seed + 501L + case$side
    )
    elapsed <- proc.time()[["elapsed"]] - started
    coords <- tr$final
    trace_selected <- extract_trace_frames(tr, target = case$truth, max_frames = cfg$trace_frames)
    note <- sprintf("Trace snapshots: %d selected of %d total frames", length(trace_selected), length(tr$frames))
  } else if (identical(method_id, "weighted_grip_polish_lgkk")) {
    started <- proc.time()[["elapsed"]]
    tr <- grip.layout.trace.weighted(
      edges = layout_graph$edges,
      edge_weights = layout_graph$edge_weights,
      n = case$top_n,
      dim = cfg$dim,
      preset = case$weighted_preset,
      lgkk_polish_rounds = cfg$weighted_polish_lgkk_rounds,
      lgkk_local_nbrs = min(cfg$lgkk_local_nbrs, max(2L, top_graph$n - 1L)),
      lgkk_landmark_count = min(cfg$lgkk_landmark_count, max(2L, top_graph$n - 1L)),
      seed = cfg$seed + 601L + case$side
    )
    elapsed <- proc.time()[["elapsed"]] - started
    coords <- tr$final
    trace_selected <- extract_trace_frames(tr, target = case$truth, max_frames = cfg$trace_frames)
    note <- sprintf("Trace snapshots: %d selected of %d total frames", length(trace_selected), length(tr$frames))
  } else if (identical(method_id, "edge_relaxation")) {
    prepared_surrogate <- make_surrogate_prepared(case$layout_prepared)
    prepared_surrogate <- grip.geodesic.mds.ensure.graph.term.cache(prepared_surrogate)
    started <- proc.time()[["elapsed"]]
    fit <- grip.optimize.geodesic.mds(
      coords = case$shared_start,
      prepared = prepared_surrogate,
      init = "user",
      anchor_mode = "none",
      engine = "cpp",
      max_iter = cfg$edge_relax_max_iter,
      edge_spring_weight = 1.0,
      return_trace = TRUE,
      recenter = TRUE
    )
    elapsed <- proc.time()[["elapsed"]] - started
    coords <- fit$coords
    note <- "Shared generic random-normal start"
  } else {
    stop("Unknown method id: ", method_id)
  }

  metrics <- compute_metrics(
    case = case,
    method_id = method_id,
    method_label = method_label,
    coords = coords,
    elapsed_sec = elapsed,
    note = note
  )
  aligned <- if (identical(method_id, "reference")) {
    list(aligned = coords, rmse = 0)
  } else {
    align_to_target_nd(coords, case$truth, allow.reflection = TRUE)
  }

  list(
    method_id = method_id,
    method_label = method_label,
    coords = coords,
    display_coords = aligned$aligned,
    metrics = metrics,
    trace_selected = trace_selected
  )
}

safe_run_method <- function(case, method_spec) {
  message("  initializer: ", method_spec$label)
  tryCatch(
    run_method(case, method_spec),
    error = function(e) {
      msg <- conditionMessage(e)
      coords <- matrix(NA_real_, nrow = case$top_n, ncol = cfg$dim)
      list(
        method_id = method_spec$id,
        method_label = method_spec$label,
        coords = coords,
        display_coords = coords,
        metrics = compute_metrics(
          case = case,
          method_id = method_spec$id,
          method_label = method_spec$label,
          coords = coords,
          elapsed_sec = NA_real_,
          note = paste("failed:", msg)
        ),
        trace_selected = list()
      )
    }
  )
}

save_final_panel <- function(case_result, file_path) {
  entries <- c(
    list(list(
      method_id = "reference",
      method_label = "Reference sample",
      display_coords = case_result$case$truth,
      metrics = compute_metrics(
        case_result$case,
        method_id = "reference",
        method_label = "Reference sample",
        coords = case_result$case$truth,
        elapsed_sec = NA_real_,
        note = "Reference coarse sample"
      )
    )),
    case_result$methods
  )

  grDevices::png(file_path, width = 2200, height = 2200, res = 220, bg = "#ffffff")
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(3L, 3L), mar = c(1.1, 1.1, 3.2, 0.4), oma = c(0, 0, 1.0, 0))

  for (entry in entries) {
    row <- entry$metrics[1L, , drop = FALSE]
    if (all(is.finite(entry$display_coords))) {
      grip.plot(
        coords = entry$display_coords,
        edges = case_result$case$display_edges,
        projection = "ortho",
        azimuth = 35,
        elevation = 24,
        vertex.col = if (identical(entry$method_id, "reference")) "#bc6c25" else "#3a5a40",
        edge.col = "#adb5bd",
        main = ""
      )
    } else {
      graphics::plot.new()
      graphics::text(0.5, 0.55, labels = row$method_label[[1L]], cex = 1.0, font = 2L)
      graphics::text(0.5, 0.40, labels = "initializer failed", cex = 0.95, col = "#8d0801")
    }
    ttl <- sprintf(
      "%s\nsigma %s, rho %s",
      row$method_label[[1L]],
      fmt_num(row$sigma_geo[[1L]], 4L),
      fmt_num(row$rho[[1L]], 4L)
    )
    graphics::mtext(ttl, side = 3L, line = 0.3, cex = 0.82)
  }
  graphics::mtext(case_result$case$label, side = 3L, outer = TRUE, line = -0.3, cex = 1.15, font = 2L)
}

save_trace_panel <- function(case_result, file_path) {
  trace.methods <- Filter(function(x) length(x$trace_selected) > 0L, case_result$methods)
  if (!length(trace.methods)) {
    return(invisible(NULL))
  }
  ncol.panel <- max(vapply(trace.methods, function(x) length(x$trace_selected), integer(1L)))
  nrow.panel <- length(trace.methods)
  grDevices::png(file_path, width = 360L * ncol.panel, height = 290L * nrow.panel, res = 190, bg = "#ffffff")
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(nrow.panel, ncol.panel), mar = c(0.8, 0.8, 2.6, 0.3), oma = c(0, 0, 1.0, 0))

  for (method in trace.methods) {
    for (j in seq_len(ncol.panel)) {
      if (j <= length(method$trace_selected)) {
        tr <- method$trace_selected[[j]]
        if (all(is.finite(tr$display_coords))) {
          grip.plot(
            coords = tr$display_coords,
            edges = case_result$case$display_edges,
            projection = "ortho",
            azimuth = 35,
            elevation = 24,
            vertex.col = "#355070",
            edge.col = "#c6d1db",
            main = ""
          )
        } else {
          graphics::plot.new()
        }
        ttl <- sprintf(
          "%s\nf %d/%d%s",
          method$method_label,
          tr$frame_index,
          tr$total_frames,
          if (is.na(tr$phase) || !nzchar(tr$phase)) "" else paste0(", ", tr$phase)
        )
        graphics::mtext(ttl, side = 3L, line = 0.25, cex = 0.72)
      } else {
        graphics::plot.new()
      }
    }
  }
  graphics::mtext(
    sprintf("%s: selected multiscale trace snapshots", case_result$case$label),
    side = 3L,
    outer = TRUE,
    line = -0.3,
    cex = 1.05,
    font = 2L
  )
}

make_case_summary_table <- function(df) {
  rows <- vapply(seq_len(nrow(df)), function(i) {
    paste(
      tex_escape(df$method_label[[i]]),
      fmt_time(df$elapsed_sec[[i]]),
      fmt_num(df$sigma_geo[[i]], 4L),
      fmt_num(df$rho[[i]], 4L),
      fmt_num(df$eta[[i]], 4L),
      fmt_num(df$alpha_0_05[[i]], 4L),
      sep = " & "
    )
  }, character(1L))
  paste(rows, collapse = " \\\\\n")
}

cases <- if (smoke) {
  list(
    make_case_regular(12),
    make_case_irregular_rectangle(15)
  )
} else {
  list(
    make_case_regular(12),
    make_case_regular(15),
    make_case_irregular_rectangle(15)
  )
}

case_results <- lapply(cases, function(case) {
  message("Running Phase A top-level panel for: ", case$label)
  method_results <- lapply(method_specs, function(spec) safe_run_method(case, spec))
  metrics <- do.call(rbind, lapply(method_results, `[[`, "metrics"))
  final_figure <- file.path(pdf_dir, sprintf("%s_top_level_initializer_grid.png", case$id))
  trace_figure <- file.path(pdf_dir, sprintf("%s_top_level_trace_grid.png", case$id))
  case_result <- list(
    case = case,
    methods = method_results,
    metrics = metrics,
    final_figure = final_figure,
    trace_figure = trace_figure
  )
  save_final_panel(case_result, final_figure)
  save_trace_panel(case_result, trace_figure)
  case_result
})

metrics_df <- do.call(rbind, lapply(case_results, `[[`, "metrics"))
utils::write.csv(metrics_df, metrics_csv, row.names = FALSE)

summary_df <- do.call(rbind, lapply(case_results, function(case_result) {
  df <- case_result$metrics
  best.sigma <- df$method_label[[which.min(df$sigma_geo)]]
  best.rho <- df$method_label[[which.min(df$rho)]]
  data.frame(
    case_id = case_result$case$id,
    case_label = case_result$case$label,
    top_level = case_result$case$top_level,
    top_n = case_result$case$top_n,
    best_sigma_method = best.sigma,
    best_rho_method = best.rho,
    stringsAsFactors = FALSE
  )
}))

bundle <- list(
  run_tag = run_tag,
  generated_at = as.character(Sys.time()),
  cfg = cfg,
  case_results = case_results,
  metrics = metrics_df,
  summary = summary_df,
  output = list(
    tex = tex_path,
    pdf = pdf_path,
    rds = rds_path,
    metrics_csv = metrics_csv
  )
)
saveRDS(bundle, rds_path)

overall_lines <- c(
  "\\documentclass[11pt]{article}",
  "\\usepackage[margin=1in]{geometry}",
  "\\usepackage{graphicx}",
  "\\usepackage{booktabs}",
  "\\usepackage{longtable}",
  "\\usepackage{xcolor}",
  "\\usepackage{hyperref}",
  "\\usepackage{float}",
  "\\hypersetup{colorlinks=true,linkcolor=blue,urlcolor=blue,citecolor=blue}",
  "\\setlength{\\parskip}{0.6em}",
  "\\setlength{\\parindent}{0pt}",
  "\\begin{document}",
  "\\begin{center}",
  "{\\LARGE Phase A: MISF-GMDS Top-Level Initializer Benchmark}\\\\[0.5em]",
  "{\\large 2026-04-02}",
  "\\end{center}",
  paste(
    "This report implements Group~A of the MISF-GMDS top-level initializer plan.",
    "Each case is reduced to its coarsest MISF active set, the induced complete weighted coarse graph is embedded by a panel of initializer methods, and every initializer is scored under the same coarse pure-GMDS objective.",
    "The displayed lines are \\emph{not} the complete coarse graph edges; for readability, every 3D snapshot is rendered with a fixed Delaunay-style triangulation of the active-set reference sample in parameter space."
  ),
  paste(
    "The compared initializers are cMDS, weighted KK, weighted FR, unweighted GRIP, Weighted GRIP, Weighted GRIP + core LGKK, Weighted GRIP + polish LGKK, and the edge-relaxation surrogate.",
    "For KK, FR, and the edge surrogate, the start is the same centered random-normal seed for each case."
  ),
  "\\section*{Overall Summary}",
  "\\small",
  "\\begin{longtable}{p{5.0cm}rrp{3.3cm}p{3.3cm}}",
  "\\toprule",
  "Case & $L$ & $|V_L|$ & Best $\\sigma_{\\mathrm{geo}}$ & Best $\\rho$ \\\\",
  "\\midrule",
  "\\endhead"
)

overall_rows <- vapply(seq_len(nrow(summary_df)), function(i) {
  paste(
    tex_escape(summary_df$case_label[[i]]),
    summary_df$top_level[[i]],
    summary_df$top_n[[i]],
    tex_escape(summary_df$best_sigma_method[[i]]),
    tex_escape(summary_df$best_rho_method[[i]]),
    sep = " & "
  )
}, character(1L))

overall_lines <- c(
  overall_lines,
  paste(overall_rows, collapse = " \\\\\n"),
  " \\\\",
  "\\bottomrule",
  "\\end{longtable}",
  "\\normalsize"
)

for (case_result in case_results) {
  case <- case_result$case
  df <- case_result$metrics[order(case_result$metrics$rho, case_result$metrics$sigma_geo), , drop = FALSE]
  best.sigma.idx <- which.min(df$sigma_geo)
  best.rho.idx <- which.min(df$rho)
  trace.methods <- Filter(function(x) length(x$trace_selected) > 0L, case_result$methods)
  note.trace <- if (length(trace.methods)) {
    paste(
      "Selected multiscale trace snapshots are shown for",
      paste(vapply(trace.methods, `[[`, "", "method_label"), collapse = ", "),
      "in the second figure."
    )
  } else {
    "No multiscale trace figure is included for this case."
  }

  overall_lines <- c(
    overall_lines,
    sprintf("\\section*{%s}", tex_escape(case$label)),
    sprintf(
      "Top-level MISF level $L=%d$ with $|V_L|=%d$ active vertices. The best coarse $\\sigma_{\\mathrm{geo}}$ in this panel is achieved by \\textbf{%s}, while the best Procrustes error $\\rho$ is achieved by \\textbf{%s}. %s",
      case$top_level,
      case$top_n,
      tex_escape(df$method_label[[best.sigma.idx]]),
      tex_escape(df$method_label[[best.rho.idx]]),
      note.trace
    ),
    "\\begin{table}[H]",
    "\\centering",
    "\\begin{tabular}{lrrrrr}",
    "\\toprule",
    "Method & $t$ (s) & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ \\\\",
    "\\midrule",
    make_case_summary_table(df),
    " \\\\",
    "\\bottomrule",
    "\\end{tabular}",
    sprintf("\\caption{Top-level-only metrics for %s. Lower $\\sigma$, $\\rho$, and $\\eta$ are better; larger $\\alpha_{0.05}$ is better.}", tex_escape(case$label)),
    "\\end{table}",
    "\\begin{figure}[H]",
    "\\centering",
    sprintf("\\includegraphics[width=0.94\\textwidth]{%s}", normalizePath(case_result$final_figure, winslash = "/", mustWork = TRUE)),
    sprintf("\\caption{Reference coarse sample and all Phase~A top-level initializers for %s. Every non-reference embedding is Procrustes-aligned to the reference sample before plotting.}", tex_escape(case$label)),
    "\\end{figure}"
  )
  if (file.exists(case_result$trace_figure)) {
    overall_lines <- c(
      overall_lines,
      "\\begin{figure}[H]",
      "\\centering",
      sprintf("\\includegraphics[width=0.98\\textwidth]{%s}", normalizePath(case_result$trace_figure, winslash = "/", mustWork = TRUE)),
      sprintf("\\caption{Selected trace snapshots for the multistage GRIP-family initializers on %s.}", tex_escape(case$label)),
      "\\end{figure}"
    )
  }
}

overall_lines <- c(
  overall_lines,
  "\\section*{Interactive Companion}",
  paste(
    "The companion HTML gallery is generated by",
    "\\texttt{tools/render-gmds-misf-top-level-initializers-html.R}.",
    "It displays the same final coarse layouts plus selected GRIP-family trace snapshots as interactive \\texttt{rglwidget} 3D panels."
  ),
  "\\end{document}"
)

writeLines(overall_lines, tex_path)

message("Wrote Phase A metrics: ", metrics_csv)
message("Wrote Phase A bundle: ", rds_path)
message("Wrote Phase A report TeX: ", tex_path)
