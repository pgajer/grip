grip_weighted_nd_trace_parity_align <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  x_center <- sweep(x, 2L, colMeans(x), "-")
  y_center <- sweep(y, 2L, colMeans(y), "-")
  denom <- sum(y_center^2)
  if (!is.finite(denom) || denom <= 0) {
    return(list(
      x_center = x_center,
      y_center = y_center,
      y_aligned = y_center,
      scale = NA_real_,
      ok = FALSE
    ))
  }
  decomp <- svd(t(y_center) %*% x_center)
  rotation <- decomp$u %*% t(decomp$v)
  y_rot <- y_center %*% rotation
  scale <- sum(x_center * y_rot) / sum(y_rot^2)
  list(
    x_center = x_center,
    y_center = y_center,
    y_aligned = y_rot * scale,
    scale = scale,
    ok = is.finite(scale)
  )
}

grip_weighted_nd_trace_parity_tuning <- function(dim = 2L) {
  list(
    dim = as.integer(dim),
    placement = "barycenter",
    rounds = 3L,
    final_rounds = 2L,
    num_init = if (dim == 2L) 5L else 6L,
    num_nbrs = 6L,
    r = 0.03,
    s = 6.0,
    repulsion_factor = 1.5,
    length_normalization = "median",
    tinit_factor = 6L,
    seed = 81L,
    trace.every = 1L
  )
}

grip_weighted_nd_trace_parity_fixtures <- function() {
  list(
    list(
      fixture_id = "mesh_3x3_uniform",
      n = 9L,
      edges = edges.mesh(3, 3),
      edge_weights = rep(1, nrow(edges.mesh(3, 3)))
    ),
    list(
      fixture_id = "cycle_8_weighted",
      n = 8L,
      edges = edges.cycle(8),
      edge_weights = seq_len(nrow(edges.cycle(8))) / nrow(edges.cycle(8)) + 1
    )
  )
}

grip_weighted_nd_trace_parity_compare_one <- function(fixture,
                                                      dim = 2L,
                                                      tuning = NULL) {
  if (is.null(tuning)) {
    tuning <- grip_weighted_nd_trace_parity_tuning(dim = dim)
  }
  common <- c(
    list(
      edges = fixture$edges,
      edge_weights = fixture$edge_weights,
      n = fixture$n
    ),
    tuning
  )
  final_anchor_factor <- if (is.null(tuning$final_anchor_factor)) {
    0
  } else {
    tuning$final_anchor_factor
  }
  final_move_scale_after_first <- if (is.null(tuning$final_move_scale_after_first)) {
    1
  } else {
    tuning$final_move_scale_after_first
  }
  legacy_common <- common[
    !names(common) %in% c("final_anchor_factor", "final_move_scale_after_first")
  ]
  legacy_args <- c(
    legacy_common,
    list(
      trace = "round",
      coarse_repulsion_factor = tuning$repulsion_factor,
      coarse_repulsion_sample = 100000L,
      coarse_repulsion_exact_below = 100000L,
      final_anchor_factor = final_anchor_factor,
      final_move_scale_after_first = final_move_scale_after_first,
      final_mode = "fr",
      lgkk_polish_rounds = 0L,
      lgkk_multiscale_rounds = 0L,
      metric_neighbor_cap = NULL,
      diagnostics = "none"
    )
  )
  nd_args <- common

  legacy <- do.call(grip.layout.trace.weighted, legacy_args)
  nd_trace_fn <- get("grip.layout.weighted.nd.trace", asNamespace("grip"))
  nd <- do.call(nd_trace_fn, nd_args)

  frame_count <- min(length(legacy$frames), length(nd$frames))
  rows <- vector("list", frame_count)
  for (idx in seq_len(frame_count)) {
    legacy_frame <- as.matrix(legacy$frames[[idx]])
    nd_frame <- as.matrix(nd$frames[[idx]])
    active <- stats::complete.cases(legacy_frame) &
      stats::complete.cases(nd_frame)
    legacy_active <- legacy_frame[active, , drop = FALSE]
    nd_active <- nd_frame[active, , drop = FALSE]
    direct_delta <- legacy_active - nd_active
    centered_rmse <- NA_real_
    centered_max_abs <- NA_real_
    procrustes_rmse <- NA_real_
    procrustes_max_abs <- NA_real_
    procrustes_scale <- NA_real_
    if (nrow(legacy_active) >= 2L && ncol(legacy_active) >= 1L) {
      aligned <- grip_weighted_nd_trace_parity_align(legacy_active, nd_active)
      centered_delta <- aligned$x_center - aligned$y_center
      centered_rmse <- sqrt(mean(centered_delta^2))
      centered_max_abs <- max(abs(centered_delta))
      if (isTRUE(aligned$ok)) {
        aligned_delta <- aligned$x_center - aligned$y_aligned
        procrustes_rmse <- sqrt(mean(aligned_delta^2))
        procrustes_max_abs <- max(abs(aligned_delta))
        procrustes_scale <- aligned$scale
      }
    }

    rows[[idx]] <- data.frame(
      fixture_id = fixture$fixture_id,
      dim = as.integer(dim),
      frame = idx,
      phase = as.character(legacy$meta$phase[[idx]]),
      nd_phase = as.character(nd$meta$phase[[idx]]),
      level_index = as.integer(legacy$meta$level_index[[idx]]),
      nd_level_index = as.integer(nd$meta$level_index[[idx]]),
      misf_level = as.integer(legacy$meta$misf_level[[idx]]),
      nd_misf_level = as.integer(nd$meta$misf_level[[idx]]),
      round_in_level = as.integer(legacy$meta$round_in_level[[idx]]),
      nd_round_in_level = as.integer(nd$meta$round_in_level[[idx]]),
      active_vertices = sum(active),
      direct_rmse = sqrt(mean(direct_delta^2)),
      direct_max_abs = max(abs(direct_delta)),
      centered_rmse = centered_rmse,
      centered_max_abs = centered_max_abs,
      procrustes_rmse = procrustes_rmse,
      procrustes_max_abs = procrustes_max_abs,
      procrustes_scale = procrustes_scale,
      metadata_match = identical(legacy$meta[idx, c("phase", "level_index", "misf_level", "round_in_level", "active_vertices")],
                                 nd$meta[idx, c("phase", "level_index", "misf_level", "round_in_level", "active_vertices")]),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  out$frame_count_match <- length(legacy$frames) == length(nd$frames)
  rownames(out) <- NULL
  out
}

grip_weighted_nd_trace_parity_run <- function(dims = c(2L, 3L),
                                              fixtures = grip_weighted_nd_trace_parity_fixtures(),
                                              tuning_fun = grip_weighted_nd_trace_parity_tuning) {
  rows <- list()
  i <- 1L
  for (fixture in fixtures) {
    for (dim in dims) {
      rows[[i]] <- grip_weighted_nd_trace_parity_compare_one(
        fixture,
        dim = dim,
        tuning = tuning_fun(dim)
      )
      i <- i + 1L
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

grip_weighted_nd_trace_parity_first_divergence <- function(metrics,
                                                           tolerance = 1e-8) {
  split_metrics <- split(metrics, list(metrics$fixture_id, metrics$dim), drop = TRUE)
  rows <- lapply(split_metrics, function(df) {
    finite <- is.finite(df$procrustes_rmse)
    divergent <- finite & df$procrustes_rmse > tolerance
    if (!any(divergent)) {
      idx <- nrow(df)
      status <- "within_tolerance"
    } else {
      idx <- which(divergent)[[1L]]
      status <- "diverged"
    }
    data.frame(
      fixture_id = df$fixture_id[[idx]],
      dim = df$dim[[idx]],
      first_divergent_frame = if (identical(status, "diverged")) df$frame[[idx]] else NA_integer_,
      first_divergent_phase = if (identical(status, "diverged")) df$phase[[idx]] else NA_character_,
      first_divergent_round = if (identical(status, "diverged")) df$round_in_level[[idx]] else NA_integer_,
      first_divergent_active_vertices = if (identical(status, "diverged")) df$active_vertices[[idx]] else NA_integer_,
      first_divergent_procrustes_rmse = if (identical(status, "diverged")) df$procrustes_rmse[[idx]] else NA_real_,
      final_procrustes_rmse = utils::tail(df$procrustes_rmse, 1L),
      metadata_all_match = all(df$metadata_match) && all(df$frame_count_match),
      status = status,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$fixture_id, out$dim), , drop = FALSE]
}

grip_weighted_nd_trace_parity_weighted_adj <- function(edges, edge_weights, n) {
  adj <- vector("list", n)
  weights <- vector("list", n)
  for (i in seq_len(nrow(edges))) {
    u <- edges[i, 1L]
    v <- edges[i, 2L]
    w <- edge_weights[[i]]
    adj[[u]] <- c(adj[[u]], v)
    adj[[v]] <- c(adj[[v]], u)
    weights[[u]] <- c(weights[[u]], w)
    weights[[v]] <- c(weights[[v]], w)
  }
  list(adj = adj, weights = weights)
}

grip_weighted_nd_trace_parity_dijkstra <- function(adj, weights, root) {
  n <- length(adj)
  dist <- rep(Inf, n)
  settled <- rep(FALSE, n)
  order <- integer()
  dist[[root]] <- 0
  repeat {
    unsettled <- which(!settled & is.finite(dist))
    if (!length(unsettled)) {
      break
    }
    best <- unsettled[order(dist[unsettled], unsettled)[[1L]]]
    settled[[best]] <- TRUE
    if (best != root) {
      order <- c(order, best)
    }
    for (idx in seq_along(adj[[best]])) {
      nbr <- adj[[best]][[idx]]
      alt <- dist[[best]] + weights[[best]][[idx]]
      scale <- max(1, abs(alt), if (is.finite(dist[[nbr]])) abs(dist[[nbr]]) else 0)
      if (!is.finite(dist[[nbr]]) || alt + 1e-10 * scale < dist[[nbr]]) {
        dist[[nbr]] <- alt
      }
    }
  }
  list(order = order, dist = dist)
}

grip_weighted_nd_trace_parity_neighbor_table <- function(misf,
                                                         adj,
                                                         weights,
                                                         level,
                                                         active_vertices = NULL) {
  if (is.null(active_vertices)) {
    active_vertices <- misf$mish_order[seq_len(misf$misf_size[[level + 1L]])]
  }
  target_count <- misf$num_nbrs_schedule[[level + 1L]]
  rows <- list()
  row_id <- 1L
  for (root in active_vertices) {
    d <- grip_weighted_nd_trace_parity_dijkstra(adj, weights, root)
    eligible <- d$order[misf$vertex_depth[d$order] >= level]
    eligible <- eligible[seq_len(min(target_count, length(eligible)))]
    if (!length(eligible)) {
      rows[[row_id]] <- data.frame(
        root = root,
        level = level,
        rank = NA_integer_,
        neighbor = NA_integer_,
        distance = NA_real_,
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
    } else {
      for (rank in seq_along(eligible)) {
        rows[[row_id]] <- data.frame(
          root = root,
          level = level,
          rank = rank,
          neighbor = eligible[[rank]],
          distance = d$dist[[eligible[[rank]]]],
          stringsAsFactors = FALSE
        )
        row_id <- row_id + 1L
      }
    }
  }
  do.call(rbind, rows)
}

grip_weighted_nd_trace_parity_neighbor_diagnostic <- function(fixture,
                                                             dim = 2L) {
  tuning <- grip_weighted_nd_trace_parity_tuning(dim = dim)
  graph <- grip_weighted_nd_trace_parity_weighted_adj(
    fixture$edges,
    fixture$edge_weights,
    fixture$n
  )
  legacy <- grip.build.misf.weighted(
    edges = fixture$edges,
    edge_weights = fixture$edge_weights,
    n = fixture$n,
    num_init = tuning$num_init,
    num_nbrs = tuning$num_nbrs,
    length_normalization = tuning$length_normalization,
    seed = tuning$seed
  )
  nd <- grip:::grip.build.misf.weighted.nd(
    edges = fixture$edges,
    edge_weights = fixture$edge_weights,
    n = fixture$n,
    num_init = tuning$num_init,
    num_nbrs = tuning$num_nbrs,
    length_normalization = tuning$length_normalization,
    seed = tuning$seed
  )
  level <- max(legacy$misf_height, 0L)
  active <- legacy$mish_order[seq_len(legacy$misf_size[[level + 1L]])]
  legacy_neighbors <- grip_weighted_nd_trace_parity_neighbor_table(
    legacy,
    graph$adj,
    graph$weights,
    level = level,
    active_vertices = active
  )
  nd_neighbors <- grip_weighted_nd_trace_parity_neighbor_table(
    nd,
    graph$adj,
    graph$weights,
    level = level,
    active_vertices = active
  )
  comparison <- merge(
    legacy_neighbors,
    nd_neighbors,
    by = c("root", "level", "rank"),
    suffixes = c("_legacy", "_nd"),
    all = TRUE
  )
  comparison$neighbor_match <- comparison$neighbor_legacy == comparison$neighbor_nd
  comparison$distance_delta <- comparison$distance_legacy - comparison$distance_nd
  comparison$distance_match <- abs(comparison$distance_delta) <= 1e-10
  comparison$fixture_id <- fixture$fixture_id
  comparison$dim <- as.integer(dim)
  comparison[, c(
    "fixture_id",
    "dim",
    "root",
    "level",
    "rank",
    "neighbor_legacy",
    "neighbor_nd",
    "distance_legacy",
    "distance_nd",
    "distance_delta",
    "neighbor_match",
    "distance_match"
  )]
}

grip_weighted_nd_trace_parity_neighbor_run <- function(
    dims = c(2L, 3L),
    fixtures = grip_weighted_nd_trace_parity_fixtures()) {
  rows <- list()
  i <- 1L
  for (fixture in fixtures) {
    for (dim in dims) {
      rows[[i]] <- grip_weighted_nd_trace_parity_neighbor_diagnostic(fixture, dim = dim)
      i <- i + 1L
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
