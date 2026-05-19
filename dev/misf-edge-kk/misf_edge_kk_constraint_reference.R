# Reference MISF and constraint helpers for grip.layout.misf.edge.kk().
#
# This file is development-only. It validates the contract that active-level
# optimization sees exactly the intended vertices and constraint sets.

mek_canonical_pairs <- function(pairs) {
  if (!is.matrix(pairs) || ncol(pairs) != 2L) {
    stop("pairs must be a two-column matrix", call. = FALSE)
  }
  out <- cbind(pmin(pairs[, 1L], pairs[, 2L]),
               pmax(pairs[, 1L], pairs[, 2L]))
  storage.mode(out) <- "integer"
  out
}

mek_pair_key <- function(pairs) {
  pairs <- mek_canonical_pairs(pairs)
  paste(pairs[, 1L], pairs[, 2L], sep = "-")
}

mek_unique_pairs <- function(pairs) {
  pairs <- mek_canonical_pairs(pairs)
  if (nrow(pairs) == 0L) {
    return(matrix(integer(0L), ncol = 2L))
  }
  pairs[!duplicated(mek_pair_key(pairs)), , drop = FALSE]
}

mek_misf_active_set <- function(misf_order, level_sizes, level_index) {
  if (level_index < 1L || level_index > length(level_sizes)) {
    stop("level_index is out of range", call. = FALSE)
  }
  misf_order <- as.integer(misf_order)
  level_sizes <- as.integer(level_sizes)
  size <- level_sizes[[level_index]]
  if (size < 1L || size > length(misf_order)) {
    stop("level size must be within misf_order length", call. = FALSE)
  }
  misf_order[seq_len(size)]
}

mek_misf_new_vertices <- function(misf_order, level_sizes, level_index) {
  active <- mek_misf_active_set(misf_order, level_sizes, level_index)
  if (level_index == 1L) {
    return(active)
  }
  previous <- mek_misf_active_set(misf_order, level_sizes, level_index - 1L)
  active[!active %in% previous]
}

mek_active_edges <- function(edges, active) {
  edges <- mek_canonical_pairs(edges)
  keep <- edges[, 1L] %in% active & edges[, 2L] %in% active
  edges[keep, , drop = FALSE]
}

mek_graph_distance_matrix <- function(n, edges, edge_lengths) {
  edges <- mek_canonical_pairs(edges)
  edge_lengths <- as.double(edge_lengths)
  if (length(edge_lengths) != nrow(edges) ||
      any(!is.finite(edge_lengths)) || any(edge_lengths <= 0)) {
    stop("edge_lengths must contain finite values > 0 parallel to edges",
         call. = FALSE)
  }
  dist <- matrix(Inf, n, n)
  diag(dist) <- 0
  for (e in seq_len(nrow(edges))) {
    u <- edges[e, 1L]
    v <- edges[e, 2L]
    w <- edge_lengths[[e]]
    if (w < dist[u, v]) {
      dist[u, v] <- w
      dist[v, u] <- w
    }
  }
  for (k in seq_len(n)) {
    for (i in seq_len(n)) {
      dik <- dist[i, k]
      if (!is.finite(dik)) next
      alt <- dik + dist[k, ]
      improve <- alt < dist[i, ]
      dist[i, improve] <- alt[improve]
    }
  }
  dist
}

mek_metric_constraints <- function(n,
                                   edges,
                                   edge_lengths,
                                   active,
                                   candidate_pairs) {
  candidate_pairs <- mek_unique_pairs(candidate_pairs)
  active <- as.integer(active)
  keep_active <- candidate_pairs[, 1L] %in% active &
    candidate_pairs[, 2L] %in% active
  candidate_pairs <- candidate_pairs[keep_active, , drop = FALSE]
  if (nrow(candidate_pairs) == 0L) {
    return(data.frame(i = integer(), j = integer(), distance = numeric()))
  }

  active_edge_keys <- mek_pair_key(mek_active_edges(edges, active))
  keep_metric <- !(mek_pair_key(candidate_pairs) %in% active_edge_keys)
  candidate_pairs <- candidate_pairs[keep_metric, , drop = FALSE]
  if (nrow(candidate_pairs) == 0L) {
    return(data.frame(i = integer(), j = integer(), distance = numeric()))
  }

  dist <- mek_graph_distance_matrix(n, edges, edge_lengths)
  metric_distance <- vapply(seq_len(nrow(candidate_pairs)), function(idx) {
    dist[candidate_pairs[idx, 1L], candidate_pairs[idx, 2L]]
  }, numeric(1L))
  finite <- is.finite(metric_distance)
  data.frame(
    i = candidate_pairs[finite, 1L],
    j = candidate_pairs[finite, 2L],
    distance = metric_distance[finite],
    stringsAsFactors = FALSE
  )
}

mek_exact_repulsion_pairs <- function(active) {
  active <- sort(as.integer(active))
  if (length(active) < 2L) {
    return(matrix(integer(0L), ncol = 2L))
  }
  pairs <- t(utils::combn(active, 2L))
  storage.mode(pairs) <- "integer"
  pairs
}

mek_sampled_repulsion_pairs <- function(active,
                                        sample_count,
                                        seed = 1L,
                                        level_index = 1L) {
  all_pairs <- mek_exact_repulsion_pairs(active)
  if (nrow(all_pairs) <= sample_count) {
    return(all_pairs)
  }
  u <- all_pairs[, 1L]
  v <- all_pairs[, 2L]
  score <- (u * 1103515245 + v * 12345 +
              as.integer(seed) * 2654435761 +
              as.integer(level_index) * 97531) %% 2147483647
  all_pairs[order(score, u, v)[seq_len(sample_count)], , drop = FALSE]
}

mek_repulsion_pairs <- function(active,
                                exact_below = 64L,
                                sample_count = 128L,
                                seed = 1L,
                                level_index = 1L) {
  if (length(active) <= exact_below) {
    mek_exact_repulsion_pairs(active)
  } else {
    mek_sampled_repulsion_pairs(active, sample_count, seed, level_index)
  }
}

mek_capture_anchors <- function(coords, active, weights = NULL) {
  coords <- mek_validate_coords(coords)
  active <- as.integer(active)
  if (any(active < 1L | active > nrow(coords))) {
    stop("active contains vertices outside coords", call. = FALSE)
  }
  if (is.null(weights)) {
    weights <- rep(1, length(active))
  }
  weights <- as.double(weights)
  if (length(weights) != length(active) ||
      any(!is.finite(weights)) || any(weights < 0)) {
    stop("weights must contain finite values >= 0 parallel to active",
         call. = FALSE)
  }
  list(
    vertices = active,
    weights = weights,
    coords = coords[active, , drop = FALSE]
  )
}

mek_contract_fixture <- function() {
  edges <- matrix(c(
    1L, 2L,
    1L, 3L,
    2L, 4L,
    3L, 5L,
    4L, 6L,
    5L, 7L,
    6L, 7L
  ), ncol = 2L, byrow = TRUE)
  list(
    n = 7L,
    edges = edges,
    edge_lengths = c(1.0, 1.2, 0.9, 1.4, 1.1, 1.3, 0.8),
    misf_order = c(4L, 2L, 6L, 1L, 3L, 5L, 7L),
    level_sizes = c(3L, 5L, 7L),
    candidate_metric_pairs = matrix(c(
      4L, 1L,
      1L, 4L,
      2L, 6L,
      1L, 3L,
      3L, 7L,
      5L, 6L,
      6L, 7L,
      7L, 1L
    ), ncol = 2L, byrow = TRUE),
    coords_post_insertion = matrix(seq(0.1, 2.1, length.out = 14L),
                                   nrow = 7L, ncol = 2L)
  )
}
