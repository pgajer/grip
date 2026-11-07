# Geodesic-MDS preparation, base energies, and base optimizer.
# Bending extensions call these explicitly; no load-time function capture.

grip.geodesic.mds.resolve.anchor <- function(anchor.mode,
                                             coords,
                                             prepared,
                                             anchor.coords = NULL,
                                             recenter = TRUE) {
  anchor.mode <- match.arg(anchor.mode, c("none", "cmdscale", "initial", "user"))
  if (identical(anchor.mode, "none")) {
    return(NULL)
  }
  out <- switch(
    anchor.mode,
    cmdscale = grip.geodesic.mds.cmdscale.init(prepared, ncol(coords)),
    initial = coords,
    user = anchor.coords
  )
  if (is.null(out)) {
    stop("anchor.coords must be supplied when anchor.mode = 'user'")
  }
  out <- grip.validate.coords(out)
  if (!identical(dim(out), dim(coords))) {
    stop("anchor.coords must have the same dimensions as coords")
  }
  if (isTRUE(recenter)) {
    out <- sweep(out, 2L, colMeans(out), "-", check.margin = FALSE)
  }
  out
}

grip.geodesic.mds.weight.schedule <- function(max.iter,
                                              weight,
                                              weight.end = weight,
                                              continuation = c("constant", "linear", "geometric")) {
  continuation <- match.arg(continuation)
  weight <- as.double(weight)
  weight.end <- as.double(weight.end)
  grip.validate.scalar(weight, "weight", lower = 0)
  grip.validate.scalar(weight.end, "weight.end", lower = 0)
  if (max.iter <= 0L) {
    return(weight)
  }
  s <- seq.int(0, max.iter) / max.iter
  if (identical(continuation, "constant")) {
    return(rep.int(weight, max.iter + 1L))
  }
  if (identical(continuation, "linear")) {
    return((1 - s) * weight + s * weight.end)
  }
  if (weight <= 0 || weight.end <= 0) {
    stop("geometric continuation requires weight and weight.end to both be > 0")
  }
  weight * (weight.end / weight)^s
}

grip.geodesic.mds.anchor.schedule <- function(max.iter,
                                              anchor.weight,
                                              anchor.weight.end = anchor.weight,
                                              continuation = c("constant", "linear", "geometric")) {
  grip.geodesic.mds.weight.schedule(
    max.iter = max.iter,
    weight = anchor.weight,
    weight.end = anchor.weight.end,
    continuation = continuation
  )
}

grip.build.graph.repulsion.cache <- function(prepared,
                                             repulsion.quantile = 0.60,
                                             repulsion.scale = 0.20,
                                             repulsion.cap.quantile = 0.90,
                                             repulsion.hop.min = 3L) {
  grip.validate.scalar(repulsion.quantile, "repulsion.quantile", lower = 0, upper = 1)
  grip.validate.scalar(repulsion.scale, "repulsion.scale", lower = 0)
  grip.validate.scalar(repulsion.cap.quantile, "repulsion.cap.quantile", lower = 0, upper = 1)
  repulsion.hop.min <- grip.validate.count(repulsion.hop.min, "repulsion.hop.min")
  if (repulsion.hop.min < 2L) {
    stop("repulsion.hop.min must be at least 2")
  }

  pair.matrix <- prepared$pair_matrix
  graph.dist <- as.double(prepared$pair_graph_distance)
  if (nrow(pair.matrix) == 0L) {
    return(list(
      repulsion_pair_matrix = matrix(integer(), ncol = 2L),
      flat_repulsion_u = integer(0L),
      flat_repulsion_v = integer(0L),
      flat_repulsion_target = numeric(0L),
      repulsion_target = numeric(0L),
      repulsion_source_distance = numeric(0L),
      repulsion_threshold = NA_real_,
      repulsion_cap = NA_real_,
      repulsion_hop_min = repulsion.hop.min
    ))
  }

  hop.matrix <- prepared$hop_distance_matrix
  if (is.null(hop.matrix)) {
    hop.matrix <- grip.graph.hop.distance.matrix(prepared$adj_list)
  }
  hop.dist <- hop.matrix[cbind(pair.matrix[, 1L], pair.matrix[, 2L])]
  eligible <- is.finite(graph.dist) & graph.dist > 0 & is.finite(hop.dist) & hop.dist >= repulsion.hop.min
  if (!any(eligible)) {
    return(list(
      repulsion_pair_matrix = matrix(integer(), ncol = 2L),
      flat_repulsion_u = integer(0L),
      flat_repulsion_v = integer(0L),
      flat_repulsion_target = numeric(0L),
      repulsion_target = numeric(0L),
      repulsion_source_distance = numeric(0L),
      repulsion_threshold = NA_real_,
      repulsion_cap = NA_real_,
      repulsion_hop_min = repulsion.hop.min
    ))
  }

  eligible.dist <- graph.dist[eligible]
  threshold <- as.double(stats::quantile(eligible.dist, probs = repulsion.quantile, names = FALSE))
  cap <- as.double(stats::quantile(eligible.dist, probs = repulsion.cap.quantile, names = FALSE))
  keep <- eligible & graph.dist >= threshold
  target <- repulsion.scale * pmin(graph.dist[keep], cap)

  list(
    repulsion_pair_matrix = pair.matrix[keep, , drop = FALSE],
    flat_repulsion_u = as.integer(pair.matrix[keep, 1L] - 1L),
    flat_repulsion_v = as.integer(pair.matrix[keep, 2L] - 1L),
    flat_repulsion_target = as.double(target),
    repulsion_target = as.double(target),
    repulsion_source_distance = graph.dist[keep],
    repulsion_threshold = threshold,
    repulsion_cap = cap,
    repulsion_hop_min = repulsion.hop.min,
    hop_distance_matrix = hop.matrix
  )
}

grip.geodesic.mds.ensure.graph.term.cache <- function(prepared,
                                                      repulsion.weight = 0,
                                                      repulsion.quantile = 0.60,
                                                      repulsion.scale = 0.20,
                                                      repulsion.cap.quantile = 0.90,
                                                      repulsion.hop.min = 3L) {
  if (is.null(prepared$graph_edge_matrix)) {
    prepared$graph_edge_matrix <- prepared$edges
  }
  if (is.null(prepared$graph_edge_target)) {
    prepared$graph_edge_target <- if (!is.null(prepared$edge_targets)) {
      as.double(prepared$edge_targets)
    } else {
      grip.edge.weights.from.adj.list(prepared$adj_list, prepared$weight_list)
    }
  }

  if (is.finite(repulsion.weight) && repulsion.weight > 0) {
    settings <- list(
      repulsion_quantile = as.double(repulsion.quantile),
      repulsion_scale = as.double(repulsion.scale),
      repulsion_cap_quantile = as.double(repulsion.cap.quantile),
      repulsion_hop_min = as.integer(repulsion.hop.min)
    )
    needs.cache <- is.null(prepared$repulsion_pair_matrix) ||
      is.null(prepared$repulsion_target) ||
      !isTRUE(identical(prepared$repulsion_settings, settings))
    if (needs.cache) {
      cache <- grip.build.graph.repulsion.cache(
        prepared = prepared,
        repulsion.quantile = repulsion.quantile,
        repulsion.scale = repulsion.scale,
        repulsion.cap.quantile = repulsion.cap.quantile,
        repulsion.hop.min = repulsion.hop.min
      )
      prepared$repulsion_pair_matrix <- cache$repulsion_pair_matrix
      prepared$flat_repulsion_u <- cache$flat_repulsion_u
      prepared$flat_repulsion_v <- cache$flat_repulsion_v
      prepared$flat_repulsion_target <- cache$flat_repulsion_target
      prepared$repulsion_target <- cache$repulsion_target
      prepared$repulsion_source_distance <- cache$repulsion_source_distance
      prepared$repulsion_threshold <- cache$repulsion_threshold
      prepared$repulsion_cap <- cache$repulsion_cap
      prepared$repulsion_settings <- settings
      prepared$hop_distance_matrix <- cache$hop_distance_matrix
    }
  }

  prepared
}

grip.geodesic.mds.resolve.anchor.vertex.weight <- function(anchor.vertex.weight,
                                                           coords) {
  if (is.null(anchor.vertex.weight)) {
    return(NULL)
  }
  weights <- as.double(anchor.vertex.weight)
  if (length(weights) != nrow(coords)) {
    stop("anchor.vertex.weight must have length nrow(coords)")
  }
  if (any(!is.finite(weights) | weights < 0)) {
    stop("anchor.vertex.weight must contain finite values >= 0")
  }
  weights
}

grip.geodesic.mds.anchor.stats <- function(coords,
                                           anchor.coords = NULL,
                                           anchor.weight = 0,
                                           anchor.vertex.weight = NULL) {
  anchor.vertex.weight <- grip.geodesic.mds.resolve.anchor.vertex.weight(
    anchor.vertex.weight = anchor.vertex.weight,
    coords = coords
  )
  if (is.null(anchor.coords) || !is.finite(anchor.weight) || anchor.weight <= 0) {
    return(list(
      anchor_weight = as.double(anchor.weight),
      raw_penalty = 0,
      energy = 0,
      gradient = matrix(0, nrow = nrow(coords), ncol = ncol(coords))
    ))
  }
  diff <- coords - anchor.coords
  if (!is.null(anchor.vertex.weight)) {
    diff <- diff * matrix(anchor.vertex.weight, nrow = nrow(coords), ncol = ncol(coords))
    raw.penalty <- sum((coords - anchor.coords)^2 * anchor.vertex.weight)
  } else {
    raw.penalty <- sum(diff^2)
  }
  list(
    anchor_weight = as.double(anchor.weight),
    raw_penalty = raw.penalty,
    energy = as.double(anchor.weight) * raw.penalty,
    gradient = 2 * as.double(anchor.weight) * diff
  )
}

grip.geodesic.mds.edge.spring.stats <- function(coords,
                                                prepared,
                                                edge.length.epsilon = 1e-8,
                                                edge.spring.weight = 0,
                                                graph.edge.matrix = NULL,
                                                graph.edge.target = NULL) {
  grip.validate.scalar(edge.spring.weight, "edge.spring.weight", lower = 0)
  grad <- matrix(0, nrow = nrow(coords), ncol = ncol(coords))
  if (!is.finite(edge.spring.weight) || edge.spring.weight <= 0) {
    return(list(
      edge_spring_weight = as.double(edge.spring.weight),
      raw_penalty = 0,
      energy = 0,
      gradient = grad,
      edge_count = 0L
    ))
  }

  edge.matrix <- graph.edge.matrix
  if (is.null(edge.matrix)) {
    edge.matrix <- prepared$graph_edge_matrix
  }
  if (is.null(edge.matrix)) {
    edge.matrix <- prepared$edges
  }
  edge.matrix <- as.matrix(edge.matrix)
  if (nrow(edge.matrix) == 0L) {
    return(list(
      edge_spring_weight = as.double(edge.spring.weight),
      raw_penalty = 0,
      energy = 0,
      gradient = grad,
      edge_count = 0L
    ))
  }

  edge.target <- graph.edge.target
  if (is.null(edge.target)) {
    edge.target <- prepared$graph_edge_target
  }
  if (is.null(edge.target)) {
    edge.target <- if (!is.null(prepared$edge_targets)) {
      as.double(prepared$edge_targets)
    } else {
      grip.edge.weights.from.adj.list(prepared$adj_list, prepared$weight_list)
    }
  }
  edge.target <- as.double(edge.target)
  if (length(edge.target) != nrow(edge.matrix)) {
    stop("graph.edge.target must be parallel to graph.edge.matrix")
  }

  diffs <- coords[edge.matrix[, 1L], , drop = FALSE] - coords[edge.matrix[, 2L], , drop = FALSE]
  edge.lengths <- sqrt(rowSums(diffs^2) + edge.length.epsilon^2)
  resid <- edge.lengths - edge.target
  unit.vecs <- diffs / edge.lengths
  for (j in seq_len(nrow(edge.matrix))) {
    u <- edge.matrix[j, 1L]
    v <- edge.matrix[j, 2L]
    step <- as.double(edge.spring.weight) * resid[[j]] * unit.vecs[j, ]
    grad[u, ] <- grad[u, ] + step
    grad[v, ] <- grad[v, ] - step
  }

  raw.penalty <- sum(resid^2)
  list(
    edge_spring_weight = as.double(edge.spring.weight),
    raw_penalty = raw.penalty,
    energy = 0.5 * as.double(edge.spring.weight) * raw.penalty,
    gradient = grad,
    edge_count = nrow(edge.matrix)
  )
}

grip.geodesic.mds.repulsion.stats <- function(coords,
                                              prepared,
                                              edge.length.epsilon = 1e-8,
                                              repulsion.weight = 0,
                                              repulsion.pair.matrix = NULL,
                                              repulsion.target = NULL) {
  grip.validate.scalar(repulsion.weight, "repulsion.weight", lower = 0)
  grad <- matrix(0, nrow = nrow(coords), ncol = ncol(coords))
  if (!is.finite(repulsion.weight) || repulsion.weight <= 0) {
    return(list(
      repulsion_weight = as.double(repulsion.weight),
      raw_penalty = 0,
      energy = 0,
      gradient = grad,
      pair_count = 0L,
      active_pair_count = 0L
    ))
  }

  pair.matrix <- repulsion.pair.matrix
  if (is.null(pair.matrix)) {
    pair.matrix <- prepared$repulsion_pair_matrix
  }
  target <- repulsion.target
  if (is.null(target)) {
    target <- prepared$repulsion_target
  }
  if (is.null(pair.matrix) || is.null(target)) {
    stop("repulsion pair cache is missing; call grip.geodesic.mds.ensure.graph.term.cache() first")
  }

  pair.matrix <- as.matrix(pair.matrix)
  target <- as.double(target)
  if (nrow(pair.matrix) == 0L) {
    return(list(
      repulsion_weight = as.double(repulsion.weight),
      raw_penalty = 0,
      energy = 0,
      gradient = grad,
      pair_count = 0L,
      active_pair_count = 0L
    ))
  }
  if (length(target) != nrow(pair.matrix)) {
    stop("repulsion.target must be parallel to repulsion.pair.matrix")
  }

  diffs <- coords[pair.matrix[, 1L], , drop = FALSE] - coords[pair.matrix[, 2L], , drop = FALSE]
  pair.lengths <- sqrt(rowSums(diffs^2) + edge.length.epsilon^2)
  resid <- pmax(target - pair.lengths, 0)
  active <- which(resid > 0)
  if (length(active) > 0L) {
    unit.vecs <- diffs[active, , drop = FALSE] / pair.lengths[active]
    for (idx in seq_along(active)) {
      row <- active[[idx]]
      u <- pair.matrix[row, 1L]
      v <- pair.matrix[row, 2L]
      step <- -as.double(repulsion.weight) * resid[[row]] * unit.vecs[idx, ]
      grad[u, ] <- grad[u, ] + step
      grad[v, ] <- grad[v, ] - step
    }
  }

  raw.penalty <- sum(resid^2)
  list(
    repulsion_weight = as.double(repulsion.weight),
    raw_penalty = raw.penalty,
    energy = 0.5 * as.double(repulsion.weight) * raw.penalty,
    gradient = grad,
    pair_count = nrow(pair.matrix),
    active_pair_count = length(active)
  )
}

grip.flatten.adj.list.zero.based <- function(adj.list) {
  n <- length(adj.list)
  deg <- lengths(adj.list)
  offsets <- integer(n + 1L)
  if (n > 0L) {
    offsets[-1L] <- cumsum(as.integer(deg))
  }
  vertices <- integer(sum(as.integer(deg)))
  cursor <- 0L
  for (i in seq_len(n)) {
    nbrs <- as.integer(adj.list[[i]])
    if (length(nbrs) == 0L) {
      next
    }
    idx <- seq.int(cursor + 1L, cursor + length(nbrs))
    vertices[idx] <- nbrs - 1L
    cursor <- cursor + length(nbrs)
  }
  list(
    flat_adj_offsets = offsets,
    flat_adj_vertices = vertices
  )
}

grip.geodesic.mds.smoothness.stats <- function(coords,
                                               prepared,
                                               smoothness.weight = 0) {
  if (!is.finite(smoothness.weight) || smoothness.weight <= 0) {
    return(list(
      smoothness_weight = as.double(smoothness.weight),
      raw_penalty = 0,
      energy = 0,
      gradient = matrix(0, nrow = nrow(coords), ncol = ncol(coords))
    ))
  }

  adj.list <- prepared$adj_list
  if (is.null(adj.list) || length(adj.list) != nrow(coords)) {
    stop("prepared must contain an adjacency list parallel to coords for smoothness regularization")
  }

  residual <- matrix(0, nrow = nrow(coords), ncol = ncol(coords))
  grad <- matrix(0, nrow = nrow(coords), ncol = ncol(coords))
  raw.penalty <- 0

  for (i in seq_len(nrow(coords))) {
    nbrs <- adj.list[[i]]
    if (length(nbrs) == 0L) {
      next
    }
    local.resid <- coords[i, ] - colMeans(coords[nbrs, , drop = FALSE])
    residual[i, ] <- local.resid
    raw.penalty <- raw.penalty + sum(local.resid^2)
  }

  grad <- 2 * as.double(smoothness.weight) * residual
  for (i in seq_len(nrow(coords))) {
    nbrs <- adj.list[[i]]
    deg <- length(nbrs)
    if (deg == 0L) {
      next
    }
    grad[nbrs, ] <- grad[nbrs, , drop = FALSE] -
      matrix(
        (2 * as.double(smoothness.weight) / deg) * residual[i, ],
        nrow = deg,
        ncol = ncol(coords),
        byrow = TRUE
      )
  }

  list(
    smoothness_weight = as.double(smoothness.weight),
    raw_penalty = raw.penalty,
    energy = as.double(smoothness.weight) * raw.penalty,
    gradient = grad
  )
}

grip.geodesic.mds.energy.gradient.base <- function(coords,
                                              prepared,
                                              edge.length.epsilon = 1e-8,
                                              anchor.coords = NULL,
                                              anchor.weight = 0,
                                              anchor.vertex.weight = NULL,
                                              smoothness.weight = 0,
                                              edge.spring.weight = 0,
                                              repulsion.weight = 0,
                                              repulsion.quantile = 0.60,
                                              repulsion.scale = 0.20,
                                              repulsion.cap.quantile = 0.90,
                                              repulsion.hop.min = 3L) {
  prepared <- grip.geodesic.mds.ensure.graph.term.cache(
    prepared = prepared,
    repulsion.weight = repulsion.weight,
    repulsion.quantile = repulsion.quantile,
    repulsion.scale = repulsion.scale,
    repulsion.cap.quantile = repulsion.cap.quantile,
    repulsion.hop.min = repulsion.hop.min
  )
  g <- as.double(prepared$pair_graph_distance)
  n.pairs <- length(g)
  grad <- matrix(0, nrow = nrow(coords), ncol = ncol(coords))
  gmds.energy <- 0
  path.lengths <- numeric(n.pairs)

  if (n.pairs == 0L) {
    anchor.stats <- grip.geodesic.mds.anchor.stats(
      coords = coords,
      anchor.coords = anchor.coords,
      anchor.weight = anchor.weight,
      anchor.vertex.weight = anchor.vertex.weight
    )
    edge.spring.stats <- grip.geodesic.mds.edge.spring.stats(
      coords = coords,
      prepared = prepared,
      edge.length.epsilon = edge.length.epsilon,
      edge.spring.weight = edge.spring.weight
    )
    repulsion.stats <- grip.geodesic.mds.repulsion.stats(
      coords = coords,
      prepared = prepared,
      edge.length.epsilon = edge.length.epsilon,
      repulsion.weight = repulsion.weight
    )
    smooth.stats <- grip.geodesic.mds.smoothness.stats(
      coords = coords,
      prepared = prepared,
      smoothness.weight = smoothness.weight
    )
    total.grad <- anchor.stats$gradient + edge.spring.stats$gradient +
      repulsion.stats$gradient + smooth.stats$gradient
    return(list(
      energy = anchor.stats$energy + edge.spring.stats$energy +
        repulsion.stats$energy + smooth.stats$energy,
      gmds_energy = 0,
      anchor_energy = anchor.stats$energy,
      edge_spring_energy = edge.spring.stats$energy,
      repulsion_energy = repulsion.stats$energy,
      smooth_energy = smooth.stats$energy,
      anchor_raw_penalty = anchor.stats$raw_penalty,
      edge_spring_raw_penalty = edge.spring.stats$raw_penalty,
      edge_spring_edge_count = edge.spring.stats$edge_count,
      repulsion_raw_penalty = repulsion.stats$raw_penalty,
      repulsion_pair_count = repulsion.stats$pair_count,
      repulsion_active_pair_count = repulsion.stats$active_pair_count,
      smooth_raw_penalty = smooth.stats$raw_penalty,
      gradient = total.grad,
      gradient_norm = sqrt(sum(total.grad^2)),
      path_lengths = path.lengths,
      target = g
    ))
  }

  for (i in seq_len(n.pairs)) {
    edges <- prepared$path_edges[[i]]
    if (nrow(edges) == 0L) {
      next
    }
    coeffs <- grip.path.edge.coefficients(prepared, i, nrow(edges))
    diffs <- coords[edges[, 1L], , drop = FALSE] - coords[edges[, 2L], , drop = FALSE]
    edge.lengths <- sqrt(rowSums(diffs^2) + edge.length.epsilon^2)
    h <- sum(coeffs * edge.lengths)
    path.lengths[i] <- h
    resid <- h - g[i]
    gmds.energy <- gmds.energy + 0.5 * resid^2
    unit.vecs <- diffs / edge.lengths
    for (j in seq_len(nrow(edges))) {
      u <- edges[j, 1L]
      v <- edges[j, 2L]
      grad[u, ] <- grad[u, ] + resid * coeffs[j] * unit.vecs[j, ]
      grad[v, ] <- grad[v, ] - resid * coeffs[j] * unit.vecs[j, ]
    }
  }

  anchor.stats <- grip.geodesic.mds.anchor.stats(
    coords = coords,
    anchor.coords = anchor.coords,
    anchor.weight = anchor.weight,
    anchor.vertex.weight = anchor.vertex.weight
  )
  edge.spring.stats <- grip.geodesic.mds.edge.spring.stats(
    coords = coords,
    prepared = prepared,
    edge.length.epsilon = edge.length.epsilon,
    edge.spring.weight = edge.spring.weight
  )
  repulsion.stats <- grip.geodesic.mds.repulsion.stats(
    coords = coords,
    prepared = prepared,
    edge.length.epsilon = edge.length.epsilon,
    repulsion.weight = repulsion.weight
  )
  smooth.stats <- grip.geodesic.mds.smoothness.stats(
    coords = coords,
    prepared = prepared,
    smoothness.weight = smoothness.weight
  )
  grad <- grad + anchor.stats$gradient + edge.spring.stats$gradient +
    repulsion.stats$gradient + smooth.stats$gradient

  list(
    energy = gmds.energy + anchor.stats$energy + edge.spring.stats$energy +
      repulsion.stats$energy + smooth.stats$energy,
    gmds_energy = gmds.energy,
    anchor_energy = anchor.stats$energy,
    edge_spring_energy = edge.spring.stats$energy,
    repulsion_energy = repulsion.stats$energy,
    smooth_energy = smooth.stats$energy,
    anchor_raw_penalty = anchor.stats$raw_penalty,
    edge_spring_raw_penalty = edge.spring.stats$raw_penalty,
    edge_spring_edge_count = edge.spring.stats$edge_count,
    repulsion_raw_penalty = repulsion.stats$raw_penalty,
    repulsion_pair_count = repulsion.stats$pair_count,
    repulsion_active_pair_count = repulsion.stats$active_pair_count,
    smooth_raw_penalty = smooth.stats$raw_penalty,
    gradient = grad,
    gradient_norm = sqrt(sum(grad^2)),
    path_lengths = path.lengths,
    target = g
  )
}

grip.geodesic.mds.score.stats.base <- function(coords,
                                          prepared,
                                          edge.length.epsilon = 1e-8,
                                          anchor.coords = NULL,
                                          anchor.weight = 0,
                                          anchor.vertex.weight = NULL,
                                          smoothness.weight = 0,
                                          edge.spring.weight = 0,
                                          repulsion.weight = 0,
                                          repulsion.quantile = 0.60,
                                          repulsion.scale = 0.20,
                                          repulsion.cap.quantile = 0.90,
                                          repulsion.hop.min = 3L) {
  grip.validate.scalar(edge.length.epsilon, "edge.length.epsilon", lower = 0)
  prepared <- grip.geodesic.mds.ensure.graph.term.cache(
    prepared = prepared,
    repulsion.weight = repulsion.weight,
    repulsion.quantile = repulsion.quantile,
    repulsion.scale = repulsion.scale,
    repulsion.cap.quantile = repulsion.cap.quantile,
    repulsion.hop.min = repulsion.hop.min
  )
  g <- as.double(prepared$pair_graph_distance)
  anchor.stats <- grip.geodesic.mds.anchor.stats(
    coords = coords,
    anchor.coords = anchor.coords,
    anchor.weight = anchor.weight,
    anchor.vertex.weight = anchor.vertex.weight
  )
  edge.spring.stats <- grip.geodesic.mds.edge.spring.stats(
    coords = coords,
    prepared = prepared,
    edge.length.epsilon = edge.length.epsilon,
    edge.spring.weight = edge.spring.weight
  )
  repulsion.stats <- grip.geodesic.mds.repulsion.stats(
    coords = coords,
    prepared = prepared,
    edge.length.epsilon = edge.length.epsilon,
    repulsion.weight = repulsion.weight
  )
  smooth.stats <- grip.geodesic.mds.smoothness.stats(
    coords = coords,
    prepared = prepared,
    smoothness.weight = smoothness.weight
  )

  if (length(g) == 0L) {
    return(list(
      n.pairs = 0L,
      energy = anchor.stats$energy + edge.spring.stats$energy +
        repulsion.stats$energy + smooth.stats$energy,
      gmds.energy = 0,
      anchor.weight = anchor.stats$anchor_weight,
      anchor.raw.penalty = anchor.stats$raw_penalty,
      anchor.energy = anchor.stats$energy,
      edge.spring.weight = edge.spring.stats$edge_spring_weight,
      edge.spring.raw.penalty = edge.spring.stats$raw_penalty,
      edge.spring.energy = edge.spring.stats$energy,
      edge.spring.edge.count = edge.spring.stats$edge_count,
      repulsion.weight = repulsion.stats$repulsion_weight,
      repulsion.raw.penalty = repulsion.stats$raw_penalty,
      repulsion.energy = repulsion.stats$energy,
      repulsion.pair.count = repulsion.stats$pair_count,
      repulsion.active.pair.count = repulsion.stats$active_pair_count,
      smooth.weight = smooth.stats$smoothness_weight,
      smooth.raw.penalty = smooth.stats$raw_penalty,
      smooth.energy = smooth.stats$energy,
      raw_stress = NA_real_,
      stress = NA_real_,
      rmse = NA_real_,
      mean.abs.path.error = NA_real_,
      mean.rel.path.error = NA_real_,
      path.lengths = numeric(0L),
      target = numeric(0L),
      residual = numeric(0L),
      relative.residual = numeric(0L)
    ))
  }

  h <- grip.geodesic.kk.path.lengths(
    coords = coords,
    prepared = prepared,
    edge.length.epsilon = edge.length.epsilon
  )
  resid <- h - g
  rel.resid <- resid / pmax(g, edge.length.epsilon)
  raw.stress <- sum(resid^2)
  denom <- sum(g^2)

  list(
    n.pairs = length(g),
    energy = 0.5 * raw.stress + anchor.stats$energy + edge.spring.stats$energy +
      repulsion.stats$energy + smooth.stats$energy,
    gmds.energy = 0.5 * raw.stress,
    anchor.weight = anchor.stats$anchor_weight,
    anchor.raw.penalty = anchor.stats$raw_penalty,
    anchor.energy = anchor.stats$energy,
    edge.spring.weight = edge.spring.stats$edge_spring_weight,
    edge.spring.raw.penalty = edge.spring.stats$raw_penalty,
    edge.spring.energy = edge.spring.stats$energy,
    edge.spring.edge.count = edge.spring.stats$edge_count,
    repulsion.weight = repulsion.stats$repulsion_weight,
    repulsion.raw.penalty = repulsion.stats$raw_penalty,
    repulsion.energy = repulsion.stats$energy,
    repulsion.pair.count = repulsion.stats$pair_count,
    repulsion.active.pair.count = repulsion.stats$active_pair_count,
    smooth.weight = smooth.stats$smoothness_weight,
    smooth.raw.penalty = smooth.stats$raw_penalty,
    smooth.energy = smooth.stats$energy,
    raw_stress = raw.stress,
    stress = if (is.finite(denom) && denom > 0) sqrt(raw.stress / denom) else NA_real_,
    rmse = sqrt(mean(resid^2)),
    mean.abs.path.error = mean(abs(resid)),
    mean.rel.path.error = mean(abs(rel.resid)),
    path.lengths = h,
    target = g,
    residual = resid,
    relative.residual = rel.resid
  )
}

grip.geodesic.mds.pair.details <- function(prepared, stats) {
  out <- data.frame(
    i = prepared$pair_matrix[, 1L],
    j = prepared$pair_matrix[, 2L],
    graph.distance = as.double(prepared$pair_graph_distance),
    embedded.path.length = stats$path.lengths,
    target.length = stats$target,
    residual = stats$residual,
    relative.residual = stats$relative.residual,
    stringsAsFactors = FALSE
  )
  if (!is.null(prepared$tie_mode)) {
    out$tie.mode <- as.character(prepared$tie_mode)
  }
  if (!is.null(prepared$pair_path_count_log)) {
    out$shortest.path.count.log <- as.double(prepared$pair_path_count_log)
    out$shortest.path.count <- grip.safe.exp(out$shortest.path.count.log)
  }
  out
}

grip.geodesic.mds.evaluate.state.base <- function(coords,
                                             prepared,
                                             edge.length.epsilon = 1e-8,
                                             anchor.coords = NULL,
                                             anchor.weight = 0,
                                             anchor.vertex.weight = NULL,
                                             smoothness.weight = 0,
                                             edge.spring.weight = 0,
                                             repulsion.weight = 0,
                                             repulsion.quantile = 0.60,
                                             repulsion.scale = 0.20,
                                             repulsion.cap.quantile = 0.90,
                                             repulsion.hop.min = 3L) {
  grip.geodesic.mds.energy.gradient(
    coords = coords,
    prepared = prepared,
    edge.length.epsilon = edge.length.epsilon,
    anchor.coords = anchor.coords,
    anchor.weight = anchor.weight,
    anchor.vertex.weight = anchor.vertex.weight,
    smoothness.weight = smoothness.weight,
    edge.spring.weight = edge.spring.weight,
    repulsion.weight = repulsion.weight,
    repulsion.quantile = repulsion.quantile,
    repulsion.scale = repulsion.scale,
    repulsion.cap.quantile = repulsion.cap.quantile,
    repulsion.hop.min = repulsion.hop.min
  )
}

grip.geodesic.mds.cmdscale.init <- function(prepared, dim) {
  dim <- grip.validate.count(dim, "dim")
  if (!(dim %in% c(2L, 3L))) {
    stop("dim must be 2 or 3")
  }
  fit <- stats::cmdscale(stats::as.dist(prepared$distance_matrix), k = dim)
  fit <- as.matrix(fit)
  if (ncol(fit) < dim) {
    fit <- cbind(fit, matrix(0, nrow = nrow(fit), ncol = dim - ncol(fit)))
  }
  storage.mode(fit) <- "double"
  fit
}

grip.canonical.edge.bundle.from.adj.list <- function(adj.list,
                                                     weight.list = NULL,
                                                     caller = "prepare.edge.kk") {
  edges <- list()
  weight.chunks <- vector("list", length(adj.list))
  for (u in seq_along(adj.list)) {
    nb <- as.integer(adj.list[[u]])
    if (length(nb) == 0L) next
    keep <- nb > u
    if (!any(keep)) next
    chosen <- nb[keep]
    edges[[length(edges) + 1L]] <- cbind(rep.int(u, length(chosen)), chosen)
    if (is.null(weight.list)) {
      weight.chunks[[u]] <- rep.int(1, length(chosen))
    } else {
      weight.chunks[[u]] <- as.double(weight.list[[u]][keep])
    }
  }

  edge.weights <- as.double(unlist(weight.chunks, use.names = FALSE))
  edges <- .bind.edges(edges)
  if (nrow(edges) == 0L) {
    return(list(edges = .empty.edge.matrix(), edge_targets = numeric(0L)))
  }
  edges <- cbind(pmin(edges[, 1L], edges[, 2L]), pmax(edges[, 1L], edges[, 2L]))
  keep <- edges[, 1L] != edges[, 2L]
  edges <- edges[keep, , drop = FALSE]
  edge.weights <- edge.weights[keep]
  if (nrow(edges) == 0L) {
    return(list(edges = .empty.edge.matrix(), edge_targets = numeric(0L)))
  }
  keys <- paste(edges[, 1L], edges[, 2L], sep = "-")
  if (any(duplicated(keys))) {
    stop(
      sprintf(
        "%s() does not support duplicate undirected edges; collapse duplicate edges before preparing edge-KK",
        caller
      ),
      call. = FALSE
    )
  }
  storage.mode(edges) <- "integer"
  list(edges = edges, edge_targets = as.double(edge.weights))
}

#' Prepare an edge-only graph for edge-KK repair
#'
#' \code{prepare.edge.kk()} validates an undirected weighted graph and
#' returns the lightweight prepared object used by
#' \code{\link{edge.kk}()} when only graph-edge targets are
#' needed. Unlike \code{\link{prepare.graph.geodesic.mds}()}, this helper
#' does not compute all-pairs shortest paths, path caches, or a dense graph
#' distance matrix.
#'
#' Use this helper when a starting layout is already available, for example from
#' \code{\link{grip}(..., metric = "edge_length")}, and the next step is
#' scalable edge-KK
#' local repair. Edge-only objects support edge-fidelity diagnostics through
#' \code{\link{score.gmds}()}; all-pairs GMDS path and chord
#' diagnostics are unavailable and are reported as \code{NA}.
#'
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#' @param n Number of vertices. If omitted with \code{adj.list}, defaults to
#'   \code{length(adj.list)}. If omitted with \code{edges}, defaults to
#'   \code{max(edges)}.
#' @param adj.list Adjacency list (1-based) for an undirected graph.
#' @param weight.list Optional parallel list of positive edge weights.
#' @param edge.weights Optional positive edge-weight vector parallel to
#'   \code{edges}.
#'
#' @return A lightweight prepared object of class
#'   \code{"grip_edge_kk_prepared"} layered on the common
#'   \code{"grip_gmds_prepared"} graph class. It contains canonical graph
#'   edges and edge targets but no all-pairs geodesic cache.
#' @export
prepare.edge.kk <- function(edges = NULL,
                            n = NULL,
                            adj.list = NULL,
                            weight.list = NULL,
                            edge.weights = NULL) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights)
  n <- grip.resolve.graph.n(n, edges, adj.list)

  validated <- grip.validate.layout.inputs(
    edges = edges,
    n = n,
    adj.list = adj.list,
    weight.list = weight.list,
    edge.weights = edge.weights,
    dim = 2L,
    placement = "barycenter",
    seed = 1L
  )
  sorted <- grip.sort.adj.with.weights(validated$adj_list, validated$weight_list)
  bundle <- grip.canonical.edge.bundle.from.adj.list(
    adj.list = sorted$adj_list,
    weight.list = sorted$weight_list,
    caller = "prepare.edge.kk"
  )
  comp <- grip.connected.components(sorted$adj_list, validated$n)

  prepared <- list(
    n = validated$n,
    edges = bundle$edges,
    edge_targets = bundle$edge_targets,
    adj_list = sorted$adj_list,
    weight_list = sorted$weight_list,
    pair_matrix = matrix(integer(), ncol = 2L),
    pair_graph_distance = numeric(0L),
    path_vertices = list(),
    path_edges = list(),
    path_edge_weights = list(),
    pair_path_count_log = numeric(0L),
    graph_diameter = NA_real_,
    distance_matrix = NULL,
    pair_mode = "edge_only",
    graph_build_mode = "edge_input",
    n_components = length(unique(comp)),
    component_id = comp
  )
  class(prepared) <- c(
    "grip_edge_kk_prepared",
    "grip_gmds_prepared",
    "grip_gkk_prepared",
    "grip_geodesic_kk_prepared",
    "list"
  )
  prepared
}

#' Prepare a graph-first geodesic-MDS path cache
#'
#' \code{prepare.graph.geodesic.mds()} prepares the full all-pairs chosen
#' geodesic cache for an arbitrary connected weighted graph. This is the
#' graph-first entry point corresponding to the manuscript's definition of GMDS
#' on a connected weighted graph together with a chosen geodesic family
#' \eqn{(G, \Gamma)}.
#'
#' The graph can be supplied either as an edge list plus parallel weights or as
#' an adjacency-list representation. The returned object stores the all-pairs
#' graph distances, the chosen shortest-path family, and the flattened edge-path
#' cache reused by the GMDS scorer and optimizer.
#'
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#' @param n Number of vertices. If omitted with \code{adj.list}, defaults to
#'   \code{length(adj.list)}. If omitted with \code{edges}, defaults to
#'   \code{max(edges)}.
#' @param adj.list Adjacency list (1-based) for an undirected graph.
#' @param weight.list Optional parallel list of positive edge weights.
#' @param edge.weights Optional positive edge-weight vector parallel to
#'   \code{edges}.
#' @param tie.mode Shortest-path aggregation mode. \code{"single"} uses one
#'   deterministic chosen shortest path per pair. \code{"average"} replaces
#'   each tied shortest-path family by the exact uniform average over all
#'   shortest paths between the pair.
#'
#' @details The all-pairs `distance.matrix` contains symmetric strict graph
#'   distances. Retained-route lengths in `pair_graph_distance` follow the
#'   deterministic near-tie path convention and can differ slightly from these
#'   targets. Route selection and path diagnostics are not changed by this
#'   distinction.
#'
#' @return A prepared object with class \code{"grip_gmds_prepared"} layered on
#'   top of the existing full geodesic path-cache structure.
#' @export
#' @md
prepare.graph.geodesic.mds <- function(edges = NULL,
                                            n = NULL,
                                            adj.list = NULL,
                                            weight.list = NULL,
                                            edge.weights = NULL,
                                            tie.mode = c("single", "average")) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights)
  tie.mode <- match.arg(tie.mode)
  base <- grip.prepare.geodesic.kk.base(
    edges = edges,
    n = n,
    adj.list = adj.list,
    weight.list = weight.list,
    edge.weights = edge.weights,
    caller = "prepare.graph.geodesic.mds"
  )
  pair.matrix <- grip.full.geodesic.kk.pair.matrix(base$n)
  cache <- grip.build.geodesic.mds.path.cache(
    pair.matrix = pair.matrix,
    adj.list = base$adj_list,
    weight.list = base$weight_list,
    dist.matrix = base$distance_matrix,
    parents = base$parents,
    tie.mode = tie.mode
  )
  flat.cache <- if (!is.null(cache$flat_pair_edge_offsets)) {
    cache[c("flat_pair_edge_offsets", "flat_edge_u", "flat_edge_v", "flat_edge_coeff")]
  } else {
    grip.flatten.geodesic.path.cache(
      path.edges = cache$path_edges,
      path.edge.weights = cache$path_edge_weights
    )
  }
  prepared <- list(
    n = base$n,
    edges = base$edges,
    edge_targets = base$edge_targets,
    adj_list = base$adj_list,
    weight_list = base$weight_list,
    pair_matrix = pair.matrix,
    pair_graph_distance = cache$pair_graph_distance,
    path_vertices = cache$path_vertices,
    path_edges = cache$path_edges,
    path_edge_weights = cache$path_edge_weights,
    pair_path_count_log = cache$pair_path_count_log,
    flat_pair_edge_offsets = flat.cache$flat_pair_edge_offsets,
    flat_edge_u = flat.cache$flat_edge_u,
    flat_edge_v = flat.cache$flat_edge_v,
    flat_edge_coeff = flat.cache$flat_edge_coeff,
    graph_diameter = base$graph_diameter,
    distance_matrix = base$mds_distance_matrix,
    pair_mode = "all_pairs",
    graph_build_mode = "graph_input",
    tie_mode = tie.mode
  )
  class(prepared) <- c("grip_gmds_prepared", "grip_gkk_prepared", "grip_geodesic_kk_prepared", "list")
  prepared
}

#' Prepare a geodesic-MDS graph and fixed path family from data
#'
#' \code{grip.prepare.geodesic.mds()} builds a deterministic symmetric
#' \eqn{k}-nearest-neighbor graph from an input data matrix, augments it to
#' connectedness when requested, and then delegates to
#' \code{\link{prepare.graph.geodesic.mds}()} to prepare the full all-pairs
#' chosen geodesic cache used by the geodesic-MDS scorer and optimizer.
#'
#' The current implementation uses Euclidean distances in the input space to
#' weight graph edges. If the symmetric \eqn{k}-NN graph is disconnected and
#' \code{connect = "mst"}, the Euclidean minimum spanning tree is unioned with
#' the \eqn{k}-NN graph before the full path cache is built. This function is
#' the data-native convenience wrapper; the graph-first API is
#' \code{\link{prepare.graph.geodesic.mds}()}.
#'
#' @param data Numeric matrix whose rows are observations.
#' @param k Symmetric \eqn{k}-NN neighborhood size.
#' @param connect Connectivity policy. \code{"mst"} augments a disconnected
#'   \eqn{k}-NN graph with Euclidean MST edges; \code{"error"} stops instead.
#' @param tie.mode Shortest-path aggregation mode. \code{"single"} uses one
#'   deterministic chosen shortest path per pair. \code{"average"} replaces
#'   each tied shortest-path family by the exact uniform average over all
#'   shortest paths between the pair.
#'
#' @return A prepared object with class \code{"grip_gmds_prepared"} layered on
#'   top of the existing full geodesic path-cache structure.
#' @noRd
grip.prepare.geodesic.mds <- function(data,
                                      k,
                                      connect = c("mst", "error"),
                                      tie.mode = c("single", "average")) {
  tie.mode <- match.arg(tie.mode)
  built <- grip.prepare.geodesic.mds.graph(
    data = data,
    k = k,
    connect = connect
  )
  prepared <- prepare.graph.geodesic.mds(
    edges = built$edges,
    n = nrow(built$data),
    edge.weights = built$edge_weights,
    tie.mode = tie.mode
  )
  prepared$input_data <- built$data
  prepared$k <- built$k
  prepared$connect <- built$connect
  prepared$knn_distance_matrix <- built$distance_matrix
  prepared$knn_edges <- built$knn_edges
  prepared$knn_edge_weights <- built$knn_edge_weights
  prepared$mst_added_edges <- built$mst_added_edges
  prepared$mst_added_edge_weights <- built$mst_added_edge_weights
  prepared$graph_build_mode <- "symmetric_knn"
  prepared
}

#' Score a layout under the geodesic-MDS objective
#'
#' \code{grip.score.geodesic.mds()} evaluates an embedding using the fixed-path
#' geodesic-MDS criterion from the manuscript: the target for each unordered
#' vertex pair is the corresponding graph geodesic itself, with no fitted scale
#' factor and no KK-style inverse-distance weighting. The scorer also allows a
#' graph-generic edge-spring term and a graph-distance-aware
#' one-sided repulsion term to be included in the reported total energy.
#'
#' @param coords Numeric coordinate matrix with 2 or 3 columns.
#' @param prepared Optional prepared geodesic object from
#'   \code{\link{prepare.graph.geodesic.mds}()},
#'   \code{\link{grip.prepare.geodesic.mds}()}, or
#'   \code{\link{prepare.geodesic.kk}()}.
#' @param data Optional data matrix used when \code{prepared} is omitted.
#' @param k Optional \eqn{k}-NN neighborhood size used when \code{prepared} is
#'   omitted.
#' @param connect Connectivity policy used when \code{prepared} is omitted.
#' @param tie.mode Shortest-path aggregation mode used when \code{prepared} is
#'   omitted.
#' @param edge.length.epsilon Small non-negative stabilizer added inside each
#'   embedded edge length.
#' @param anchor.coords Optional anchor embedding used to add the quadratic
#'   tether term \eqn{\lambda \|Z - A\|_F^2}.
#' @param anchor.weight Non-negative anchor weight \eqn{\lambda}.
#' @param anchor.vertex.weight Optional non-negative per-vertex weights for the
#'   anchor term. When supplied, only vertices with positive weights are
#'   tethered, and each tether uses the corresponding multiplier.
#' @param smoothness.weight Non-negative local smoothness weight \eqn{\mu}
#'   applied to \eqn{\sum_i \|z_i - |N(i)|^{-1}\sum_{j \in N(i)} z_j\|^2}.
#' @param edge.spring.weight Non-negative coefficient for the graph-edge spring
#'   term \eqn{\frac{\beta}{2}\sum_{(u,v)\in E}(\|z_u-z_v\|-w_{uv})^2}.
#' @param repulsion.weight Non-negative coefficient for the graph-distance-aware
#'   repulsion term applied to graph-distant vertex pairs.
#' @param repulsion.quantile Graph-distance quantile used to select the
#'   repulsion pair family from all eligible nonlocal pairs.
#' @param repulsion.scale Positive scale factor converting graph distances into
#'   one-sided Euclidean separation targets.
#' @param repulsion.cap.quantile Upper graph-distance quantile used to cap the
#'   repulsion targets before scaling.
#' @param repulsion.hop.min Minimum graph hop distance required for a pair to be
#'   eligible for the repulsion family.
#' @param bending.stencils Optional list describing discrete bending stencils for
#'   the layout. Each stencil should identify the three vertices forming the
#'   local bending constraint and may optionally supply per-stencil weights.
#' @param bending.weight Non-negative coefficient for the bending penalty term.
#' @param return.pair.details If \code{TRUE}, attach per-pair residual details.
#'
#' @return A one-row data frame summarizing the geodesic-MDS fit.
#' @noRd
grip.score.geodesic.mds.base <- function(coords,
                                    prepared = NULL,
                                    data = NULL,
                                    k = NULL,
                                    connect = c("mst", "error"),
                                    tie.mode = c("single", "average"),
                                    edge.length.epsilon = 1e-8,
                                    anchor.coords = NULL,
                                    anchor.weight = 0,
                                    anchor.vertex.weight = NULL,
                                    smoothness.weight = 0,
                                    edge.spring.weight = 0,
                                    repulsion.weight = 0,
                                    repulsion.quantile = 0.60,
                                    repulsion.scale = 0.20,
                                    repulsion.cap.quantile = 0.90,
                                    repulsion.hop.min = 3L,
                                    return.pair.details = FALSE) {
  coords <- grip.validate.coords(coords)
  tie.mode <- match.arg(tie.mode)
  if (is.null(prepared)) {
    prepared <- grip.prepare.geodesic.mds(
      data = data,
      k = k,
      connect = connect,
      tie.mode = tie.mode
    )
  }
  prepared <- grip.validate.geodesic.mds.prepared(prepared, coords = coords)
  prepared <- grip.geodesic.mds.ensure.graph.term.cache(
    prepared = prepared,
    repulsion.weight = repulsion.weight,
    repulsion.quantile = repulsion.quantile,
    repulsion.scale = repulsion.scale,
    repulsion.cap.quantile = repulsion.cap.quantile,
    repulsion.hop.min = repulsion.hop.min
  )
  anchor.coords <- if (is.null(anchor.coords)) NULL else {
    grip.geodesic.mds.resolve.anchor(
      anchor.mode = "user",
      coords = coords,
      prepared = prepared,
      anchor.coords = anchor.coords,
      recenter = FALSE
    )
  }
  stats <- grip.geodesic.mds.score.stats(
    coords = coords,
    prepared = prepared,
    edge.length.epsilon = edge.length.epsilon,
    anchor.coords = anchor.coords,
    anchor.weight = anchor.weight,
    anchor.vertex.weight = anchor.vertex.weight,
    smoothness.weight = smoothness.weight,
    edge.spring.weight = edge.spring.weight,
    repulsion.weight = repulsion.weight,
    repulsion.quantile = repulsion.quantile,
    repulsion.scale = repulsion.scale,
    repulsion.cap.quantile = repulsion.cap.quantile,
    repulsion.hop.min = repulsion.hop.min
  )

  out <- data.frame(
    n = prepared$n,
    n.pairs = stats$n.pairs,
    gmds.energy = stats$energy,
    gmds.base.energy = stats$gmds.energy,
    gmds.raw_stress = stats$raw_stress,
    gmds.stress = stats$stress,
    gmds.rmse = stats$rmse,
    gmds.mean.abs.path.error = stats$mean.abs.path.error,
    gmds.mean.rel.path.error = stats$mean.rel.path.error,
    anchor.weight = stats$anchor.weight,
    anchor.raw.penalty = stats$anchor.raw.penalty,
    anchor.energy = stats$anchor.energy,
    edge.spring.weight = stats$edge.spring.weight,
    edge.spring.raw.penalty = stats$edge.spring.raw.penalty,
    edge.spring.energy = stats$edge.spring.energy,
    edge.spring.edge.count = stats$edge.spring.edge.count,
    repulsion.weight = stats$repulsion.weight,
    repulsion.raw.penalty = stats$repulsion.raw.penalty,
    repulsion.energy = stats$repulsion.energy,
    repulsion.pair.count = stats$repulsion.pair.count,
    repulsion.active.pair.count = stats$repulsion.active.pair.count,
    smooth.weight = stats$smooth.weight,
    smooth.raw.penalty = stats$smooth.raw.penalty,
    smooth.energy = stats$smooth.energy,
    tie.mode = if (!is.null(prepared$tie_mode)) prepared$tie_mode else "single",
    stringsAsFactors = FALSE
  )

  if (isTRUE(return.pair.details)) {
    out$pair.details <- list(grip.geodesic.mds.pair.details(prepared, stats))
  }

  out
}

#' Optimize a layout under the geodesic-MDS objective
#'
#' \code{grip.optimize.geodesic.mds()} optimizes a fixed-path geodesic-MDS
#' objective using deterministic gradient descent with Armijo backtracking. By
#' default it initializes from classical MDS on the graph geodesic distance
#' matrix and then runs a compiled optimizer. When the flattened all-pairs path
#' cache is available, the compiled path also supports the optional
#' graph-edge-spring, graph-distance-aware repulsion, anchor, and smoothness
#' continuation terms.
#'
#' @param coords Optional numeric coordinate matrix with 2 or 3 columns. If
#'   omitted, coordinates are initialized according to \code{init}.
#' @param prepared Optional prepared object from
#'   \code{\link{prepare.graph.geodesic.mds}()},
#'   \code{\link{grip.prepare.geodesic.mds}()}, or
#'   \code{\link{prepare.geodesic.kk}()}.
#' @param data Optional input data matrix used when \code{prepared} is omitted.
#' @param k Optional \eqn{k}-NN neighborhood size used when \code{prepared} is
#'   omitted.
#' @param dim Target output dimension used when \code{coords} is omitted.
#' @param connect Connectivity policy used when \code{prepared} is omitted.
#' @param tie.mode Shortest-path aggregation mode used when \code{prepared} is
#'   omitted.
#' @param init Initialization mode: \code{"cmdscale"}, \code{"random"}, or
#'   \code{"user"}.
#' @param anchor.mode Anchor family used by the quadratic tether term.
#'   \code{"cmdscale"} tethers to the classical-MDS embedding of the graph
#'   distance matrix, \code{"initial"} tethers to the starting layout, and
#'   \code{"user"} uses \code{anchor.coords}.
#' @param anchor.coords Optional anchor coordinates used when
#'   \code{anchor.mode = "user"}.
#' @param anchor.weight Initial non-negative tether weight.
#' @param anchor.weight.end Final non-negative tether weight used at the end of
#'   the continuation schedule.
#' @param anchor.vertex.weight Optional non-negative per-vertex weights for the
#'   anchor term. This can be used to pin only selected vertices during
#'   refinement while leaving the remaining vertices free.
#' @param continuation Continuation schedule for the tether weight. With
#'   \code{"constant"}, the tether weight stays fixed; \code{"linear"} and
#'   \code{"geometric"} gradually relax the tether from
#'   \code{anchor.weight} to \code{anchor.weight.end}.
#' @param smoothness.weight Initial non-negative local smoothness weight.
#' @param smoothness.weight.end Final non-negative smoothness weight used at the
#'   end of the continuation schedule.
#' @param smoothness.continuation Continuation schedule for the smoothness
#'   weight. The same options and semantics as \code{continuation}.
#' @param edge.spring.weight Initial non-negative graph-edge spring weight.
#' @param edge.spring.weight.end Final non-negative edge-spring weight used at
#'   the end of the continuation schedule.
#' @param edge.spring.continuation Continuation schedule for the edge-spring
#'   weight. The same options and semantics as \code{continuation}.
#' @param repulsion.weight Initial non-negative graph-distance-aware repulsion
#'   weight.
#' @param repulsion.weight.end Final non-negative repulsion weight used at the
#'   end of the continuation schedule.
#' @param repulsion.continuation Continuation schedule for the repulsion
#'   weight. The same options and semantics as \code{continuation}.
#' @param repulsion.quantile Graph-distance quantile used to define the
#'   repulsion pair family.
#' @param repulsion.scale Positive scale factor converting selected graph
#'   distances into one-sided Euclidean separation targets.
#' @param repulsion.cap.quantile Upper graph-distance quantile used to cap the
#'   repulsion targets before scaling.
#' @param repulsion.hop.min Minimum graph hop distance required for a pair to be
#'   eligible for repulsion.
#' @param bending.stencils Optional list describing discrete bending stencils for
#'   the layout. Each stencil should identify the three vertices forming the
#'   local bending constraint and may optionally supply per-stencil weights.
#' @param bending.weight Initial non-negative bending penalty weight.
#' @param bending.weight.end Final non-negative bending penalty weight used at
#'   the end of the continuation schedule.
#' @param bending.continuation Continuation schedule for the bending penalty
#'   weight. The same options and semantics as \code{continuation}.
#' @param engine Optimization engine: \code{"cpp"} or \code{"r"}.
#' @param max.iter Maximum number of gradient-descent iterations.
#' @param edge.length.epsilon Small non-negative stabilizer added inside each
#'   embedded edge length.
#' @param initial.step Initial line-search step size.
#' @param step.shrink Multiplicative backtracking shrink factor in `(0, 1)`.
#' @param armijo.factor Non-negative Armijo decrease constant.
#' @param grad.tol Non-negative stopping tolerance on the gradient norm.
#' @param min.step Positive minimum accepted line-search step.
#' @param n.threads Number of CPU threads used by the flattened compiled
#'   optimizer. \code{0} picks an automatic value, \code{1} forces serial
#'   evaluation, and all values are capped at two threads.
#'   With \code{n.threads = 0}, a positive integer in the environment variable
#'   \env{GRIP_GMDS_THREADS} takes precedence over hardware-based selection.
#'   An explicit positive \code{n.threads} overrides the environment variable.
#'   Set the variable to \code{1} or \code{2}; non-positive or nonnumeric values
#'   fall back to hardware-based selection, and larger values are capped at two.
#'   R-engine optimization is serial; the setting applies to the compiled engine.
#' @param recenter If \code{TRUE}, recenter accepted proposals to zero mean.
#' @param return.trace If \code{TRUE}, include per-iteration diagnostics and
#'   accepted frames.
#' @param seed Optional integer seed used only for random initialization.
#'
#' @return A list with \code{coords}, \code{trace}, \code{frames},
#'   \code{prepared}, and \code{score}.
#' @noRd
grip.optimize.geodesic.mds.base <- function(coords = NULL,
                                       prepared = NULL,
                                       data = NULL,
                                       k = NULL,
                                       dim = 2L,
                                       connect = c("mst", "error"),
                                       tie.mode = c("single", "average"),
                                       init = c("cmdscale", "random", "user"),
                                       anchor.mode = c("none", "cmdscale", "initial", "user"),
                                       anchor.coords = NULL,
                                       anchor.weight = 0,
                                       anchor.weight.end = anchor.weight,
                                       anchor.vertex.weight = NULL,
                                       continuation = c("constant", "linear", "geometric"),
                                       smoothness.weight = 0,
                                       smoothness.weight.end = smoothness.weight,
                                       smoothness.continuation = c("constant", "linear", "geometric"),
                                       edge.spring.weight = 0,
                                       edge.spring.weight.end = edge.spring.weight,
                                       edge.spring.continuation = c("constant", "linear", "geometric"),
                                       repulsion.weight = 0,
                                       repulsion.weight.end = repulsion.weight,
                                       repulsion.continuation = c("constant", "linear", "geometric"),
                                       repulsion.quantile = 0.60,
                                       repulsion.scale = 0.20,
                                       repulsion.cap.quantile = 0.90,
                                       repulsion.hop.min = 3L,
                                       engine = c("cpp", "r"),
                                       max.iter = 16L,
                                       edge.length.epsilon = 1e-8,
                                       initial.step = 1.0,
                                       step.shrink = 0.5,
                                       armijo.factor = 1e-4,
                                       grad.tol = 1e-8,
                                       min.step = 1e-8,
                                       n.threads = 0L,
                                       recenter = TRUE,
                                       return.trace = FALSE,
                                       seed = NULL) {
  init <- match.arg(init)
  tie.mode <- match.arg(tie.mode)
  anchor.mode <- match.arg(anchor.mode)
  continuation <- match.arg(continuation)
  smoothness.continuation <- match.arg(smoothness.continuation)
  edge.spring.continuation <- match.arg(edge.spring.continuation)
  repulsion.continuation <- match.arg(repulsion.continuation)
  engine <- match.arg(engine)
  grip.validate.scalar(max.iter, "max.iter", lower = 0)
  grip.validate.scalar(edge.length.epsilon, "edge.length.epsilon", lower = 0)
  grip.validate.scalar(initial.step, "initial.step", lower = 0, open.lower = TRUE)
  grip.validate.scalar(step.shrink, "step.shrink", lower = 0, upper = 1, open.lower = TRUE, open.upper = TRUE)
  grip.validate.scalar(armijo.factor, "armijo.factor", lower = 0)
  grip.validate.scalar(grad.tol, "grad.tol", lower = 0)
  grip.validate.scalar(min.step, "min.step", lower = 0, open.lower = TRUE)
  grip.validate.scalar(n.threads, "n.threads", lower = 0)
  grip.validate.scalar(anchor.weight, "anchor.weight", lower = 0)
  grip.validate.scalar(anchor.weight.end, "anchor.weight.end", lower = 0)
  grip.validate.scalar(smoothness.weight, "smoothness.weight", lower = 0)
  grip.validate.scalar(smoothness.weight.end, "smoothness.weight.end", lower = 0)
  grip.validate.scalar(edge.spring.weight, "edge.spring.weight", lower = 0)
  grip.validate.scalar(edge.spring.weight.end, "edge.spring.weight.end", lower = 0)
  grip.validate.scalar(repulsion.weight, "repulsion.weight", lower = 0)
  grip.validate.scalar(repulsion.weight.end, "repulsion.weight.end", lower = 0)
  grip.validate.scalar(repulsion.quantile, "repulsion.quantile", lower = 0, upper = 1)
  grip.validate.scalar(repulsion.scale, "repulsion.scale", lower = 0)
  grip.validate.scalar(repulsion.cap.quantile, "repulsion.cap.quantile", lower = 0, upper = 1)
  repulsion.hop.min <- grip.validate.count(repulsion.hop.min, "repulsion.hop.min")
  if (repulsion.hop.min < 2L) {
    stop("repulsion.hop.min must be at least 2")
  }
  if (identical(anchor.mode, "none") && (anchor.weight > 0 || anchor.weight.end > 0)) {
    stop("anchor.mode must not be 'none' when anchor.weight or anchor.weight.end is positive")
  }
  max.iter <- as.integer(round(max.iter))
  n.threads <- as.integer(round(n.threads))
  if (is.null(prepared)) {
    prepared <- grip.prepare.geodesic.mds(
      data = data,
      k = k,
      connect = connect,
      tie.mode = tie.mode
    )
  }

  if (is.null(coords)) {
    dim <- grip.validate.count(dim, "dim")
    if (!(dim %in% c(2L, 3L))) {
      stop("dim must be 2 or 3")
    }
    if (identical(init, "user")) {
      stop("coords must be supplied when init = 'user'")
    }
    if (identical(init, "cmdscale")) {
      coords <- grip.geodesic.mds.cmdscale.init(prepared, dim)
    } else {
      if (!is.null(seed)) {
        set.seed(as.integer(seed))
      }
      coords <- matrix(stats::rnorm(prepared$n * dim), ncol = dim)
    }
  } else {
    coords <- grip.validate.coords(coords)
  }

  prepared <- grip.validate.geodesic.mds.prepared(prepared, coords = coords)
  anchor.vertex.weight <- grip.geodesic.mds.resolve.anchor.vertex.weight(
    anchor.vertex.weight = anchor.vertex.weight,
    coords = coords
  )
  anchor.coords <- grip.geodesic.mds.resolve.anchor(
    anchor.mode = anchor.mode,
    coords = coords,
    prepared = prepared,
    anchor.coords = anchor.coords,
    recenter = recenter
  )
  anchor.schedule <- if (is.null(anchor.coords)) {
    rep.int(0, max.iter + 1L)
  } else {
    grip.geodesic.mds.weight.schedule(
      max.iter = max.iter,
      weight = anchor.weight,
      weight.end = anchor.weight.end,
      continuation = continuation
    )
  }
  smoothness.schedule <- grip.geodesic.mds.weight.schedule(
    max.iter = max.iter,
    weight = smoothness.weight,
    weight.end = smoothness.weight.end,
    continuation = smoothness.continuation
  )
  edge.spring.schedule <- grip.geodesic.mds.weight.schedule(
    max.iter = max.iter,
    weight = edge.spring.weight,
    weight.end = edge.spring.weight.end,
    continuation = edge.spring.continuation
  )
  repulsion.schedule <- grip.geodesic.mds.weight.schedule(
    max.iter = max.iter,
    weight = repulsion.weight,
    weight.end = repulsion.weight.end,
    continuation = repulsion.continuation
  )
  prepared <- grip.geodesic.mds.ensure.graph.term.cache(
    prepared = prepared,
    repulsion.weight = max(repulsion.schedule),
    repulsion.quantile = repulsion.quantile,
    repulsion.scale = repulsion.scale,
    repulsion.cap.quantile = repulsion.cap.quantile,
    repulsion.hop.min = repulsion.hop.min
  )

  if (nrow(coords) <= 1L || length(prepared$pair_graph_distance) == 0L || max.iter == 0L) {
    score <- grip.score.geodesic.mds(
      coords = coords,
      prepared = prepared,
      edge.length.epsilon = edge.length.epsilon,
      anchor.coords = anchor.coords,
      anchor.weight = anchor.schedule[[1L]],
      anchor.vertex.weight = anchor.vertex.weight,
      smoothness.weight = smoothness.schedule[[1L]],
      edge.spring.weight = edge.spring.schedule[[1L]],
      repulsion.weight = repulsion.schedule[[1L]],
      repulsion.quantile = repulsion.quantile,
      repulsion.scale = repulsion.scale,
      repulsion.cap.quantile = repulsion.cap.quantile,
      repulsion.hop.min = repulsion.hop.min
    )
    return(list(
      coords = coords,
      trace = data.frame(),
      frames = list(coords),
      prepared = prepared,
      score = score,
      anchor_coords = anchor.coords,
      anchor_schedule = anchor.schedule,
      anchor_vertex_weight = anchor.vertex.weight,
      smoothness_schedule = smoothness.schedule,
      edge_spring_schedule = edge.spring.schedule,
      repulsion_schedule = repulsion.schedule,
      final_anchor_weight = anchor.schedule[[1L]],
      final_smoothness_weight = smoothness.schedule[[1L]],
      final_edge_spring_weight = edge.spring.schedule[[1L]],
      final_repulsion_weight = repulsion.schedule[[1L]],
      n_threads_used = 1L
    ))
  }

  if (identical(engine, "cpp")) {
    if ((!is.null(anchor.vertex.weight) ||
         any(smoothness.schedule > 0) ||
         any(edge.spring.schedule > 0) ||
         any(repulsion.schedule > 0)) &&
        (is.null(prepared$flat_pair_edge_offsets) ||
         is.null(prepared$flat_edge_u) ||
         is.null(prepared$flat_edge_v) ||
         is.null(prepared$flat_edge_coeff))) {
      warning("graph regularization in the compiled engine requires the flattened cache; falling back to the R engine")
      engine <- "r"
    }
  }

  if (identical(engine, "cpp")) {
    smooth.flat <- if (any(smoothness.schedule > 0)) {
      grip.flatten.adj.list.zero.based(prepared$adj_list)
    } else {
      list(flat_adj_offsets = integer(), flat_adj_vertices = integer())
    }
    graph.edge.u <- if (!is.null(prepared$graph_edge_matrix) && nrow(prepared$graph_edge_matrix) > 0L) {
      as.integer(prepared$graph_edge_matrix[, 1L] - 1L)
    } else {
      integer(0L)
    }
    graph.edge.v <- if (!is.null(prepared$graph_edge_matrix) && nrow(prepared$graph_edge_matrix) > 0L) {
      as.integer(prepared$graph_edge_matrix[, 2L] - 1L)
    } else {
      integer(0L)
    }
    graph.edge.target <- if (!is.null(prepared$graph_edge_target)) {
      as.double(prepared$graph_edge_target)
    } else {
      numeric(0L)
    }
    repulsion.u <- if (!is.null(prepared$flat_repulsion_u)) {
      as.integer(prepared$flat_repulsion_u)
    } else if (!is.null(prepared$repulsion_pair_matrix) && nrow(prepared$repulsion_pair_matrix) > 0L) {
      as.integer(prepared$repulsion_pair_matrix[, 1L] - 1L)
    } else {
      integer(0L)
    }
    repulsion.v <- if (!is.null(prepared$flat_repulsion_v)) {
      as.integer(prepared$flat_repulsion_v)
    } else if (!is.null(prepared$repulsion_pair_matrix) && nrow(prepared$repulsion_pair_matrix) > 0L) {
      as.integer(prepared$repulsion_pair_matrix[, 2L] - 1L)
    } else {
      integer(0L)
    }
    repulsion.target <- if (!is.null(prepared$flat_repulsion_target)) {
      as.double(prepared$flat_repulsion_target)
    } else if (!is.null(prepared$repulsion_target)) {
      as.double(prepared$repulsion_target)
    } else {
      numeric(0L)
    }
    opt <- if (!is.null(prepared$flat_pair_edge_offsets) &&
               !is.null(prepared$flat_edge_u) &&
               !is.null(prepared$flat_edge_v) &&
               !is.null(prepared$flat_edge_coeff)) {
      grip_optimize_geodesic_mds_flat_cpp(
        flat_pair_edge_offsets = prepared$flat_pair_edge_offsets,
        flat_edge_u = prepared$flat_edge_u,
        flat_edge_v = prepared$flat_edge_v,
        flat_edge_coeff = prepared$flat_edge_coeff,
        pair_graph_distance = prepared$pair_graph_distance,
        coords = coords,
        max_iter = max.iter,
        edge_length_epsilon = edge.length.epsilon,
        initial_step = initial.step,
        step_shrink = step.shrink,
        armijo_factor = armijo.factor,
        grad_tol = grad.tol,
        min_step = min.step,
        recenter = recenter,
        return_trace = return.trace,
        anchor_coords = anchor.coords,
        anchor_weights = anchor.schedule,
        anchor_vertex_weight = anchor.vertex.weight,
        smooth_adj_offsets = smooth.flat$flat_adj_offsets,
        smooth_adj_vertices = smooth.flat$flat_adj_vertices,
        smooth_weights = smoothness.schedule,
        graph_edge_u = graph.edge.u,
        graph_edge_v = graph.edge.v,
        graph_edge_target = graph.edge.target,
        edge_spring_weights = edge.spring.schedule,
        repulsion_u = repulsion.u,
        repulsion_v = repulsion.v,
        repulsion_target = repulsion.target,
        repulsion_weights = repulsion.schedule,
        n_threads = n.threads
      )
    } else {
      out <- grip_optimize_geodesic_mds_cache_cpp(
        path_edges = prepared$path_edges,
        path_edge_weights = prepared$path_edge_weights,
        pair_graph_distance = prepared$pair_graph_distance,
        coords = coords,
        max_iter = max.iter,
        edge_length_epsilon = edge.length.epsilon,
        initial_step = initial.step,
        step_shrink = step.shrink,
        armijo_factor = armijo.factor,
        grad_tol = grad.tol,
        min_step = min.step,
        recenter = recenter,
        return_trace = return.trace,
        anchor_coords = anchor.coords,
        anchor_weights = anchor.schedule
      )
      out$n_threads_used <- 1L
      out
    }
    final.anchor.weight <- opt$final_anchor_weight
    final.edge.spring.weight <- if (!is.null(opt$final_edge_spring_weight)) opt$final_edge_spring_weight else 0
    final.repulsion.weight <- if (!is.null(opt$final_repulsion_weight)) opt$final_repulsion_weight else 0
    score <- grip.score.geodesic.mds(
      coords = opt$coords,
      prepared = prepared,
      edge.length.epsilon = edge.length.epsilon,
      anchor.coords = anchor.coords,
      anchor.weight = final.anchor.weight,
      anchor.vertex.weight = anchor.vertex.weight,
      smoothness.weight = opt$final_smoothness_weight,
      edge.spring.weight = final.edge.spring.weight,
      repulsion.weight = final.repulsion.weight,
      repulsion.quantile = repulsion.quantile,
      repulsion.scale = repulsion.scale,
      repulsion.cap.quantile = repulsion.cap.quantile,
      repulsion.hop.min = repulsion.hop.min
    )
    return(list(
      coords = opt$coords,
      trace = opt$trace,
      frames = opt$frames,
      prepared = prepared,
      score = score,
      anchor_coords = anchor.coords,
      anchor_schedule = anchor.schedule,
      anchor_vertex_weight = anchor.vertex.weight,
      smoothness_schedule = smoothness.schedule,
      edge_spring_schedule = edge.spring.schedule,
      repulsion_schedule = repulsion.schedule,
      final_anchor_weight = final.anchor.weight,
      final_smoothness_weight = opt$final_smoothness_weight,
      final_edge_spring_weight = final.edge.spring.weight,
      final_repulsion_weight = final.repulsion.weight,
      n_threads_used = opt$n_threads_used
    ))
  }

  current <- coords
  trace.rows <- vector("list", max.iter + 1L)
  accepted.frames <- list(current)
  state <- grip.geodesic.mds.evaluate.state(
    coords = current,
    prepared = prepared,
    edge.length.epsilon = edge.length.epsilon,
    anchor.coords = anchor.coords,
    anchor.weight = anchor.schedule[[1L]],
    anchor.vertex.weight = anchor.vertex.weight,
    smoothness.weight = smoothness.schedule[[1L]],
    edge.spring.weight = edge.spring.schedule[[1L]],
    repulsion.weight = repulsion.schedule[[1L]],
    repulsion.quantile = repulsion.quantile,
    repulsion.scale = repulsion.scale,
    repulsion.cap.quantile = repulsion.cap.quantile,
    repulsion.hop.min = repulsion.hop.min
  )
  trace.rows[[1L]] <- data.frame(
    iteration = 0L,
    energy = state$energy,
    gmds_energy = state$gmds_energy,
    anchor_energy = state$anchor_energy,
    edge_spring_energy = state$edge_spring_energy,
    repulsion_energy = state$repulsion_energy,
    repulsion_pair_count = state$repulsion_pair_count,
    repulsion_active_pair_count = state$repulsion_active_pair_count,
    smooth_energy = state$smooth_energy,
    gradient_norm = state$gradient_norm,
    step = NA_real_,
    accepted = TRUE,
    anchor_weight = anchor.schedule[[1L]],
    edge_spring_weight = edge.spring.schedule[[1L]],
    repulsion_weight = repulsion.schedule[[1L]],
    smooth_weight = smoothness.schedule[[1L]],
    stringsAsFactors = FALSE
  )
  used <- 1L

  for (iter in seq_len(max.iter)) {
    iter.anchor.weight <- anchor.schedule[[iter + 1L]]
    iter.smooth.weight <- smoothness.schedule[[iter + 1L]]
    iter.edge.spring.weight <- edge.spring.schedule[[iter + 1L]]
    iter.repulsion.weight <- repulsion.schedule[[iter + 1L]]
    state <- grip.geodesic.mds.evaluate.state(
      coords = current,
      prepared = prepared,
      edge.length.epsilon = edge.length.epsilon,
      anchor.coords = anchor.coords,
      anchor.weight = iter.anchor.weight,
      anchor.vertex.weight = anchor.vertex.weight,
      smoothness.weight = iter.smooth.weight,
      edge.spring.weight = iter.edge.spring.weight,
      repulsion.weight = iter.repulsion.weight,
      repulsion.quantile = repulsion.quantile,
      repulsion.scale = repulsion.scale,
      repulsion.cap.quantile = repulsion.cap.quantile,
      repulsion.hop.min = repulsion.hop.min
    )
    if (!is.finite(state$gradient_norm) || state$gradient_norm <= grad.tol) {
      break
    }
    step <- as.double(initial.step)
    accepted <- FALSE
    candidate <- current
    candidate.state <- state

    while (is.finite(step) && step >= min.step) {
      proposal <- current - step * state$gradient
      if (isTRUE(recenter)) {
        proposal <- sweep(proposal, 2L, colMeans(proposal), "-", check.margin = FALSE)
      }
      proposal.state <- grip.geodesic.mds.evaluate.state(
        coords = proposal,
        prepared = prepared,
        edge.length.epsilon = edge.length.epsilon,
      anchor.coords = anchor.coords,
      anchor.weight = iter.anchor.weight,
      anchor.vertex.weight = anchor.vertex.weight,
      smoothness.weight = iter.smooth.weight,
        edge.spring.weight = iter.edge.spring.weight,
        repulsion.weight = iter.repulsion.weight,
        repulsion.quantile = repulsion.quantile,
        repulsion.scale = repulsion.scale,
        repulsion.cap.quantile = repulsion.cap.quantile,
        repulsion.hop.min = repulsion.hop.min
      )
      target.energy <- state$energy - armijo.factor * step * state$gradient_norm^2
      if (is.finite(proposal.state$energy) && proposal.state$energy <= target.energy) {
        candidate <- proposal
        candidate.state <- proposal.state
        accepted <- TRUE
        break
      }
      step <- step * step.shrink
    }

    used <- used + 1L
    trace.rows[[used]] <- data.frame(
      iteration = iter,
      energy = if (accepted) candidate.state$energy else state$energy,
      gmds_energy = if (accepted) candidate.state$gmds_energy else state$gmds_energy,
      anchor_energy = if (accepted) candidate.state$anchor_energy else state$anchor_energy,
      edge_spring_energy = if (accepted) candidate.state$edge_spring_energy else state$edge_spring_energy,
      repulsion_energy = if (accepted) candidate.state$repulsion_energy else state$repulsion_energy,
      repulsion_pair_count = if (accepted) candidate.state$repulsion_pair_count else state$repulsion_pair_count,
      repulsion_active_pair_count = if (accepted) candidate.state$repulsion_active_pair_count else state$repulsion_active_pair_count,
      smooth_energy = if (accepted) candidate.state$smooth_energy else state$smooth_energy,
      gradient_norm = if (accepted) candidate.state$gradient_norm else state$gradient_norm,
      step = if (accepted) step else NA_real_,
      accepted = accepted,
      anchor_weight = iter.anchor.weight,
      edge_spring_weight = iter.edge.spring.weight,
      repulsion_weight = iter.repulsion.weight,
      smooth_weight = iter.smooth.weight,
      stringsAsFactors = FALSE
    )

    if (!accepted) {
      break
    }

    current <- candidate
    state <- candidate.state
    accepted.frames[[length(accepted.frames) + 1L]] <- current
  }

  trace.df <- do.call(rbind, trace.rows[seq_len(used)])
  if (!isTRUE(return.trace)) {
    trace.df <- trace.df[, c("iteration", "energy", "gmds_energy", "anchor_energy", "edge_spring_energy", "repulsion_energy", "repulsion_pair_count", "repulsion_active_pair_count", "smooth_energy", "gradient_norm", "step", "accepted", "anchor_weight", "edge_spring_weight", "repulsion_weight", "smooth_weight"), drop = FALSE]
    accepted.frames <- list(current)
  }
  final.anchor.weight <- trace.df$anchor_weight[[nrow(trace.df)]]
  final.edge.spring.weight <- trace.df$edge_spring_weight[[nrow(trace.df)]]
  final.repulsion.weight <- trace.df$repulsion_weight[[nrow(trace.df)]]
  final.smooth.weight <- trace.df$smooth_weight[[nrow(trace.df)]]
  score <- grip.score.geodesic.mds(
    coords = current,
    prepared = prepared,
    edge.length.epsilon = edge.length.epsilon,
    anchor.coords = anchor.coords,
    anchor.weight = final.anchor.weight,
    anchor.vertex.weight = anchor.vertex.weight,
    smoothness.weight = final.smooth.weight,
    edge.spring.weight = final.edge.spring.weight,
    repulsion.weight = final.repulsion.weight,
    repulsion.quantile = repulsion.quantile,
    repulsion.scale = repulsion.scale,
    repulsion.cap.quantile = repulsion.cap.quantile,
    repulsion.hop.min = repulsion.hop.min
  )

  list(
    coords = current,
    trace = trace.df,
    frames = accepted.frames,
    prepared = prepared,
    score = score,
    anchor_coords = anchor.coords,
    anchor_schedule = anchor.schedule,
    anchor_vertex_weight = anchor.vertex.weight,
    smoothness_schedule = smoothness.schedule,
    edge_spring_schedule = edge.spring.schedule,
    repulsion_schedule = repulsion.schedule,
    final_anchor_weight = final.anchor.weight,
    final_smoothness_weight = final.smooth.weight,
    final_edge_spring_weight = final.edge.spring.weight,
    final_repulsion_weight = final.repulsion.weight,
    n_threads_used = 1L
  )
}

