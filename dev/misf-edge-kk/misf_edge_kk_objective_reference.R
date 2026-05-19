# Reference objective pieces for the proposed grip.layout.misf.edge.kk().
#
# This file is intentionally kept under dev/ rather than R/. It is a small,
# auditable R prototype used to validate the objective before any package API or
# compiled backend is added.

mek_validate_coords <- function(coords) {
  if (!is.matrix(coords) || !is.numeric(coords)) {
    stop("coords must be a numeric matrix", call. = FALSE)
  }
  if (nrow(coords) < 2L || ncol(coords) < 2L) {
    stop("coords must have at least two rows and two columns", call. = FALSE)
  }
  if (any(!is.finite(coords))) {
    stop("coords must contain only finite values", call. = FALSE)
  }
  coords
}

mek_validate_pair_matrix <- function(pairs, n, name = "pairs") {
  if (!is.matrix(pairs) || ncol(pairs) != 2L) {
    stop(sprintf("%s must be a two-column matrix", name), call. = FALSE)
  }
  storage.mode(pairs) <- "integer"
  if (any(!is.finite(pairs)) || any(pairs < 1L) || any(pairs > n)) {
    stop(sprintf("%s must contain finite 1-based vertex ids", name), call. = FALSE)
  }
  if (any(pairs[, 1L] == pairs[, 2L])) {
    stop(sprintf("%s must not contain self pairs", name), call. = FALSE)
  }
  pairs
}

mek_validate_positive_vector <- function(x, len, name) {
  x <- as.double(x)
  if (length(x) != len || any(!is.finite(x)) || any(x <= 0)) {
    stop(sprintf("%s must contain finite values > 0 parallel to pairs", name),
         call. = FALSE)
  }
  x
}

mek_pair_stress_state <- function(coords,
                                  pairs,
                                  targets,
                                  stiffness = rep(1, nrow(pairs)),
                                  scale = 1,
                                  eps = 1e-8) {
  coords <- mek_validate_coords(coords)
  pairs <- mek_validate_pair_matrix(pairs, nrow(coords), "pairs")
  targets <- mek_validate_positive_vector(targets, nrow(pairs), "targets")
  stiffness <- mek_validate_positive_vector(stiffness, nrow(pairs), "stiffness")
  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    stop("scale must be a single finite value > 0", call. = FALSE)
  }
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps < 0) {
    stop("eps must be a single finite value >= 0", call. = FALSE)
  }

  grad <- matrix(0, nrow(coords), ncol(coords))
  lengths <- numeric(nrow(pairs))
  residuals <- numeric(nrow(pairs))
  energy <- 0
  denom <- 0

  for (e in seq_len(nrow(pairs))) {
    u <- pairs[e, 1L]
    v <- pairs[e, 2L]
    delta <- coords[u, ] - coords[v, ]
    len <- sqrt(sum(delta^2) + eps^2)
    target <- scale * targets[[e]]
    residual <- len - target
    coeff <- stiffness[[e]] * residual / len
    contribution <- coeff * delta

    grad[u, ] <- grad[u, ] + contribution
    grad[v, ] <- grad[v, ] - contribution
    lengths[[e]] <- len
    residuals[[e]] <- residual
    energy <- energy + 0.5 * stiffness[[e]] * residual^2
    denom <- denom + stiffness[[e]] * target^2
  }

  rel_rmse <- if (is.finite(denom) && denom > 0) {
    sqrt(sum(stiffness * residuals^2) / denom)
  } else {
    NA_real_
  }
  list(
    energy = energy,
    gradient = grad,
    gradient_norm = sqrt(sum(grad^2)),
    lengths = lengths,
    residuals = residuals,
    rel_rmse = rel_rmse,
    scale = scale
  )
}

mek_edge_stress_state <- function(coords, edges, edge_lengths, stiffness,
                                  scale = 1, eps = 1e-8) {
  mek_pair_stress_state(
    coords = coords,
    pairs = edges,
    targets = edge_lengths,
    stiffness = stiffness,
    scale = scale,
    eps = eps
  )
}

mek_metric_stress_state <- function(coords, metric_pairs, metric_distances,
                                    stiffness, scale = 1, eps = 1e-8) {
  mek_pair_stress_state(
    coords = coords,
    pairs = metric_pairs,
    targets = metric_distances,
    stiffness = stiffness,
    scale = scale,
    eps = eps
  )
}

mek_log_repulsion_state <- function(coords,
                                    pairs,
                                    pair_weights,
                                    eps = 1e-8) {
  coords <- mek_validate_coords(coords)
  pairs <- mek_validate_pair_matrix(pairs, nrow(coords), "pairs")
  pair_weights <- as.double(pair_weights)
  if (length(pair_weights) != nrow(pairs) ||
      any(!is.finite(pair_weights)) || any(pair_weights < 0)) {
    stop("pair_weights must contain finite values >= 0 parallel to pairs",
         call. = FALSE)
  }
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("eps must be a single finite value > 0", call. = FALSE)
  }

  grad <- matrix(0, nrow(coords), ncol(coords))
  lengths <- numeric(nrow(pairs))
  energy <- 0
  eps2 <- eps^2

  for (e in seq_len(nrow(pairs))) {
    u <- pairs[e, 1L]
    v <- pairs[e, 2L]
    delta <- coords[u, ] - coords[v, ]
    denom <- sum(delta^2) + eps2
    len <- sqrt(denom)
    w <- pair_weights[[e]]
    contribution <- -w * delta / denom

    grad[u, ] <- grad[u, ] + contribution
    grad[v, ] <- grad[v, ] - contribution
    lengths[[e]] <- len
    energy <- energy - w * log(len)
  }

  list(
    energy = energy,
    gradient = grad,
    gradient_norm = sqrt(sum(grad^2)),
    lengths = lengths
  )
}

mek_anchor_state <- function(coords,
                             anchor_targets,
                             anchor_weights) {
  coords <- mek_validate_coords(coords)
  anchor_targets <- mek_validate_coords(anchor_targets)
  if (!identical(dim(coords), dim(anchor_targets))) {
    stop("anchor_targets must have the same dimensions as coords", call. = FALSE)
  }
  anchor_weights <- as.double(anchor_weights)
  if (length(anchor_weights) != nrow(coords) ||
      any(!is.finite(anchor_weights)) || any(anchor_weights < 0)) {
    stop("anchor_weights must contain finite values >= 0 parallel to coords rows",
         call. = FALSE)
  }

  delta <- coords - anchor_targets
  grad <- delta * anchor_weights
  energy <- 0.5 * sum(anchor_weights * rowSums(delta^2))
  list(
    energy = energy,
    gradient = grad,
    gradient_norm = sqrt(sum(grad^2))
  )
}

mek_total_state <- function(coords,
                            edges,
                            edge_lengths,
                            edge_stiffness,
                            metric_pairs,
                            metric_distances,
                            metric_stiffness,
                            repulsion_pairs,
                            repulsion_weights,
                            anchor_targets,
                            anchor_weights,
                            rho = 0,
                            lambda = 0,
                            edge_scale = 1,
                            metric_scale = 1,
                            edge_eps = 1e-8,
                            repulsion_eps = 1e-8) {
  coords <- mek_validate_coords(coords)
  zero_state <- list(
    energy = 0,
    gradient = matrix(0, nrow(coords), ncol(coords)),
    gradient_norm = 0,
    rel_rmse = NA_real_
  )

  edge <- if (nrow(edges)) {
    mek_edge_stress_state(coords, edges, edge_lengths, edge_stiffness,
                          scale = edge_scale, eps = edge_eps)
  } else {
    zero_state
  }
  metric <- if (nrow(metric_pairs)) {
    mek_metric_stress_state(coords, metric_pairs, metric_distances,
                            metric_stiffness, scale = metric_scale,
                            eps = edge_eps)
  } else {
    zero_state
  }
  repulsion <- if (nrow(repulsion_pairs)) {
    mek_log_repulsion_state(coords, repulsion_pairs, repulsion_weights,
                            eps = repulsion_eps)
  } else {
    zero_state
  }
  anchor <- mek_anchor_state(coords, anchor_targets, anchor_weights)

  gradient <- edge$gradient +
    rho * metric$gradient +
    lambda * repulsion$gradient +
    anchor$gradient
  total_energy <- edge$energy +
    rho * metric$energy +
    lambda * repulsion$energy +
    anchor$energy

  list(
    energy = total_energy,
    gradient = gradient,
    gradient_norm = sqrt(sum(gradient^2)),
    edge = edge,
    metric = metric,
    repulsion = repulsion,
    anchor = anchor,
    rho = rho,
    lambda = lambda
  )
}

mek_finite_difference_gradient <- function(coords, fn, eps = 1e-6) {
  coords <- mek_validate_coords(coords)
  grad <- matrix(0, nrow(coords), ncol(coords))
  for (i in seq_len(nrow(coords))) {
    for (j in seq_len(ncol(coords))) {
      plus <- coords
      minus <- coords
      plus[i, j] <- plus[i, j] + eps
      minus[i, j] <- minus[i, j] - eps
      grad[i, j] <- (fn(plus) - fn(minus)) / (2 * eps)
    }
  }
  grad
}

mek_gradient_check <- function(name,
                               coords,
                               state_fn,
                               eps = 1e-6,
                               tolerance = 1e-5) {
  analytic <- state_fn(coords)
  fd <- mek_finite_difference_gradient(coords, function(x) state_fn(x)$energy,
                                       eps = eps)
  abs_error <- max(abs(analytic$gradient - fd))
  rel_error <- sqrt(sum((analytic$gradient - fd)^2)) /
    max(1, sqrt(sum(fd^2)))
  data.frame(
    gate = name,
    passed = is.finite(rel_error) && rel_error <= tolerance,
    abs_error = abs_error,
    rel_error = rel_error,
    tolerance = tolerance,
    stringsAsFactors = FALSE
  )
}

mek_normalize_repulsion_weights <- function(pairs, active_count) {
  if (nrow(pairs) == 0L) {
    return(numeric(0L))
  }
  rep(active_count / nrow(pairs), nrow(pairs))
}

mek_validate_misf_scale_policy <- function(scale_mode, rho, lambda) {
  scale_mode <- match.arg(scale_mode, c("identity", "global_fixed",
                                        "level_fixed", "user", "profiled"))
  if (identical(scale_mode, "profiled") && (rho > 0 || lambda > 0)) {
    stop("profiled scale is not allowed when rho > 0 or lambda > 0",
         call. = FALSE)
  }
  TRUE
}

mek_objective_fixture <- function(dim = 3L, repeated = FALSE) {
  if (dim < 2L) {
    stop("dim must be >= 2", call. = FALSE)
  }
  set.seed(1103 + dim + if (repeated) 100L else 0L)
  coords <- matrix(rnorm(5L * dim), nrow = 5L, ncol = dim)
  if (repeated) {
    coords[2L, ] <- coords[1L, ]
  }
  edges <- matrix(c(
    1L, 2L,
    2L, 3L,
    3L, 4L,
    4L, 5L
  ), ncol = 2L, byrow = TRUE)
  metric_pairs <- matrix(c(
    1L, 3L,
    1L, 4L,
    2L, 5L
  ), ncol = 2L, byrow = TRUE)
  repulsion_pairs <- utils::combn(5L, 2L)
  repulsion_pairs <- t(repulsion_pairs)
  anchor_targets <- coords + matrix(seq_len(length(coords)), nrow(coords), dim) * 1e-3

  list(
    coords = coords,
    edges = edges,
    edge_lengths = c(1.0, 1.4, 0.8, 1.2),
    edge_stiffness = c(1.0, 1.5, 0.75, 1.25),
    metric_pairs = metric_pairs,
    metric_distances = c(2.4, 3.1, 2.7),
    metric_stiffness = c(0.5, 0.75, 0.6),
    repulsion_pairs = repulsion_pairs,
    repulsion_weights = mek_normalize_repulsion_weights(repulsion_pairs, 5L),
    anchor_targets = anchor_targets,
    anchor_weights = c(0.05, 0.1, 0, 0.2, 0.15)
  )
}
