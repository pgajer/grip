# Baseline-comparison helpers for the MISF edge-KK reference prototype.
#
# This file is development-only. It compares the reference wrapper against the
# current package workflows used for stretched isometric graph embeddings.

mek_all_candidate_pairs <- function(n) {
  pairs <- t(utils::combn(seq_len(n), 2L))
  storage.mode(pairs) <- "integer"
  pairs
}

mek_cmdscale_insertion_coords <- function(n, edges, edge_lengths, dim = 2L) {
  dist <- mek_graph_distance_matrix(n, edges, edge_lengths)
  coords <- stats::cmdscale(stats::as.dist(dist), k = dim, eig = FALSE)
  coords <- as.matrix(coords)
  if (ncol(coords) < dim) {
    coords <- cbind(coords, matrix(0, nrow = n, ncol = dim - ncol(coords)))
  }
  coords <- coords[, seq_len(dim), drop = FALSE]
  sweep(coords, 2L, colMeans(coords), "-", check.margin = FALSE)
}

mek_grid_edges <- function(nx, ny) {
  id <- function(x, y) (y - 1L) * nx + x
  edges <- list()
  for (y in seq_len(ny)) {
    for (x in seq_len(nx)) {
      if (x < nx) {
        edges[[length(edges) + 1L]] <- c(id(x, y), id(x + 1L, y))
      }
      if (y < ny) {
        edges[[length(edges) + 1L]] <- c(id(x, y), id(x, y + 1L))
      }
    }
  }
  edges <- do.call(rbind, edges)
  storage.mode(edges) <- "integer"
  edges
}

mek_fixture <- function(name, n, edges, edge_lengths, misf_order, level_sizes,
                        dim = 2L) {
  list(
    name = name,
    n = as.integer(n),
    edges = mek_canonical_pairs(edges),
    edge_lengths = as.double(edge_lengths),
    misf_order = as.integer(misf_order),
    level_sizes = as.integer(level_sizes),
    candidate_metric_pairs = mek_all_candidate_pairs(n),
    coords_post_insertion = mek_cmdscale_insertion_coords(
      n, edges, edge_lengths, dim = dim
    ),
    dim = as.integer(dim)
  )
}

mek_baseline_fixtures <- function(dim = 2L) {
  path_edges <- cbind(1:4, 2:5)
  cycle_edges <- rbind(cbind(1:5, 2:6), c(6L, 1L))
  grid_edges <- mek_grid_edges(3L, 3L)
  grid_lengths <- rep(1, nrow(grid_edges))
  weighted_grid_lengths <- vapply(seq_len(nrow(grid_edges)), function(i) {
    u <- grid_edges[i, 1L]
    v <- grid_edges[i, 2L]
    horizontal <- abs(u - v) == 1L
    if (horizontal) {
      0.85 + 0.1 * ((min(u, v) - 1L) %% 3L)
    } else {
      1.15 + 0.08 * floor((min(u, v) - 1L) / 3L)
    }
  }, numeric(1L))

  list(
    mek_fixture(
      name = "path5",
      n = 5L,
      edges = path_edges,
      edge_lengths = c(1, 1.2, 0.9, 1.1),
      misf_order = c(1L, 3L, 5L, 2L, 4L),
      level_sizes = c(3L, 4L, 5L),
      dim = dim
    ),
    mek_fixture(
      name = "cycle6",
      n = 6L,
      edges = cycle_edges,
      edge_lengths = rep(1, 6L),
      misf_order = c(1L, 3L, 5L, 2L, 4L, 6L),
      level_sizes = c(3L, 5L, 6L),
      dim = dim
    ),
    mek_fixture(
      name = "grid3x3",
      n = 9L,
      edges = grid_edges,
      edge_lengths = grid_lengths,
      misf_order = c(1L, 3L, 7L, 9L, 5L, 2L, 4L, 6L, 8L),
      level_sizes = c(4L, 7L, 9L),
      dim = dim
    ),
    mek_fixture(
      name = "weighted_grid3x3",
      n = 9L,
      edges = grid_edges,
      edge_lengths = weighted_grid_lengths,
      misf_order = c(1L, 3L, 7L, 9L, 5L, 2L, 4L, 6L, 8L),
      level_sizes = c(4L, 7L, 9L),
      dim = dim
    )
  )
}

mek_procrustes_rmse <- function(coords, reference) {
  coords <- as.matrix(coords)
  reference <- as.matrix(reference)
  if (!identical(dim(coords), dim(reference)) ||
      any(!is.finite(coords)) || any(!is.finite(reference))) {
    return(NA_real_)
  }
  x <- sweep(coords, 2L, colMeans(coords), "-", check.margin = FALSE)
  y <- sweep(reference, 2L, colMeans(reference), "-", check.margin = FALSE)
  xnorm <- sqrt(sum(x^2))
  ynorm <- sqrt(sum(y^2))
  if (!is.finite(xnorm) || !is.finite(ynorm) || xnorm <= 0 || ynorm <= 0) {
    return(NA_real_)
  }
  x <- x / xnorm
  y <- y / ynorm
  sv <- svd(t(x) %*% y)
  rotation <- sv$u %*% t(sv$v)
  aligned <- x %*% rotation
  sqrt(mean(rowSums((aligned - y)^2)))
}

mek_min_pair_distance <- function(coords) {
  coords <- as.matrix(coords)
  if (nrow(coords) < 2L || any(!is.finite(coords))) {
    return(NA_real_)
  }
  min(stats::dist(coords))
}

mek_edge_energy_from_score <- function(coords, prepared, diagnostics) {
  edge_lengths <- prepared$edge_targets
  edges <- prepared$edges
  scale <- diagnostics$edge.scale[[1L]]
  state <- mek_edge_stress_state(
    coords = coords,
    edges = edges,
    edge_lengths = edge_lengths,
    stiffness = rep(1, nrow(edges)),
    scale = scale
  )
  state$energy
}

mek_score_baseline_coords <- function(coords, prepared, reference_coords = NULL) {
  if (is.null(coords) || any(!is.finite(coords))) {
    stop("coordinates are missing or non-finite", call. = FALSE)
  }
  diagnostics <- grip.score.gmds.layout(
    coords = coords,
    prepared = prepared,
    scale_mode = "profiled"
  )
  data.frame(
    finite = all(is.finite(coords)),
    edge_rel_rmse = diagnostics$edge.rel.rmse[[1L]],
    edge_rmse = diagnostics$edge.rmse[[1L]],
    gmds_stress = diagnostics$gmds.stress[[1L]],
    edge_energy = mek_edge_energy_from_score(coords, prepared, diagnostics),
    shape_rmse_vs_mds_edge_kk = if (is.null(reference_coords)) {
      NA_real_
    } else {
      mek_procrustes_rmse(coords, reference_coords)
    },
    stringsAsFactors = FALSE
  )
}

mek_run_reference_method <- function(fixture) {
  level_count <- length(fixture$level_sizes)
  mek_layout_misf_edge_kk_reference(
    n = fixture$n,
    edges = fixture$edges,
    edge_lengths = fixture$edge_lengths,
    coords_post_insertion = fixture$coords_post_insertion,
    misf_order = fixture$misf_order,
    level_sizes = fixture$level_sizes,
    candidate_metric_pairs = fixture$candidate_metric_pairs,
    rho_schedule = seq(0.6, 0, length.out = level_count),
    lambda_schedule = seq(0.1, 0, length.out = level_count),
    anchor_weight_schedule = seq(0.25, 0, length.out = level_count),
    exact_repulsion_below = 16L,
    repulsion_sample_count = 32L,
    repulsion_seed = 11L,
    level_max_iter = 12L,
    level_initial_step = 0.1,
    polish_max_iter = 30L,
    polish_initial_step = 0.1,
    polish_scale_mode = "profiled"
  )
}

mek_run_mds_edge_kk_method <- function(fixture, prepared) {
  mds <- grip.metric.mds.layout(
    prepared = prepared,
    dim = fixture$dim,
    diagnostics = FALSE
  )
  grip.optimize.edge.kk.layout(
    coords = mds$coords,
    prepared = prepared,
    dim = fixture$dim,
    stiffness_method = "uniform",
    density_mix_schedule = 1,
    scale_mode = "profiled",
    max_iter = 30L,
    initial_step = 0.2,
    return_trace = TRUE,
    diagnostics = TRUE,
    engine = "cpp"
  )
}

mek_run_weighted_edge_kk_method <- function(fixture, prepared) {
  start <- grip.layout.weighted(
    edges = fixture$edges,
    n = fixture$n,
    edge_weights = fixture$edge_lengths,
    dim = fixture$dim,
    rounds = 30L,
    final_rounds = 50L,
    num_init = min(4L, fixture$n),
    num_nbrs = min(6L, fixture$n - 1L),
    lgkk_polish_rounds = 0L,
    seed = 17L
  )
  grip.optimize.edge.kk.layout(
    coords = start,
    prepared = prepared,
    dim = fixture$dim,
    stiffness_method = "uniform",
    density_mix_schedule = 1,
    scale_mode = "profiled",
    max_iter = 30L,
    initial_step = 0.2,
    return_trace = TRUE,
    diagnostics = TRUE,
    engine = "cpp"
  )
}

mek_timed <- function(expr) {
  elapsed <- system.time(value <- force(expr))[["elapsed"]]
  list(value = value, elapsed = elapsed)
}
