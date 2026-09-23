# Display-only registration of triangulated surfaces. No optimizer is rerun.
# Deterministic, area-weighted surface samples approximate the two surfaces.
comparison.surface.samples <- function(coords, triangles, subdivisions = 3L) {
  stopifnot(is.matrix(coords), ncol(coords) == 3L, all(is.finite(coords)),
    is.matrix(triangles), ncol(triangles) == 3L, nrow(triangles) > 0L,
    all(triangles == as.integer(triangles)), all(triangles >= 1L),
    all(triangles <= nrow(coords)), length(subdivisions) == 1L,
    subdivisions >= 1L, subdivisions == as.integer(subdivisions))
  m <- as.integer(subdivisions)
  bary <- list()
  for (i in 0:(m - 1L)) for (j in 0:(m - 1L - i)) {
    bary[[length(bary) + 1L]] <- c(i + 1/3, j + 1/3) / m
    if (i + j < m - 1L)
      bary[[length(bary) + 1L]] <- c(i + 2/3, j + 2/3) / m
  }
  uv <- do.call(rbind, bary)
  bary <- cbind(1 - rowSums(uv), uv)
  a <- coords[triangles[, 2L], , drop = FALSE] - coords[triangles[, 1L], , drop = FALSE]
  b <- coords[triangles[, 3L], , drop = FALSE] - coords[triangles[, 1L], , drop = FALSE]
  cross <- cbind(a[,2]*b[,3]-a[,3]*b[,2], a[,3]*b[,1]-a[,1]*b[,3],
                 a[,1]*b[,2]-a[,2]*b[,1])
  area <- sqrt(rowSums(cross^2)) / 2
  stopifnot(sum(area) > 0)
  keep <- which(area > 0)
  list(points = do.call(rbind, lapply(keep, function(i)
      bary %*% coords[triangles[i, ], , drop = FALSE])),
    weights = rep(area[keep] / sum(area) / nrow(bary), each = nrow(bary)))
}

comparison.surface.matches <- function(source, target) {
  forward <- FNN::get.knnx(target$points, source$points, k = 1L)
  backward <- FNN::get.knnx(source$points, target$points, k = 1L)
  list(forward = forward$nn.index[,1], backward = backward$nn.index[,1],
    mse = (sum(source$weights * forward$nn.dist[,1]^2) +
      sum(target$weights * backward$nn.dist[,1]^2)) / 2)
}

comparison.surface.align <- function(coords, reference, triangles,
                                     subdivisions = 3L, max.iter = 80L,
                                     tolerance = 1e-9,
                                     angles = seq(0, 345, by = 15) * pi / 180) {
  if (!requireNamespace('FNN', quietly = TRUE)) stop('Surface alignment requires FNN.')
  stopifnot(max.iter >= 1L, tolerance > 0, length(angles) > 0L, all(is.finite(angles)))
  # Vertex Procrustes supplies the initial coordinate frame, not the matches
  # optimized below. Its reflection is permitted, as in the original overlays.
  initial <- comparison_align(coords, reference)
  source <- comparison.surface.samples(initial, triangles, subdivisions)
  target <- comparison.surface.samples(reference, triangles, subdivisions)
  center <- colSums(target$points * target$weights)
  initial.mse <- comparison.surface.matches(source, target)$mse
  best <- list(coords = initial, mse = initial.mse, iterations = 0L,
    converged = NA, start = 0L)
  runs <- list()
  weights <- c(source$weights, target$weights) / 2
  # The reference saddle has vertical axis z. Search both vertical reflections
  # and a full turn about that axis; each ICP update can rotate about any axis.
  for (flip in c(1, -1)) for (angle in angles) {
    rotation <- matrix(c(cos(angle), sin(angle), 0, -sin(angle), cos(angle), 0,
                         0, 0, flip), 3L, 3L)
    transform <- function(x) sweep(sweep(x, 2L, center) %*% rotation, 2L, center, '+')
    moved <- source; moved$points <- transform(source$points)
    vertices <- transform(initial)
    match <- comparison.surface.matches(moved, target)
    converged <- FALSE
    for (iteration in seq_len(max.iter)) {
      A <- rbind(moved$points, moved$points[match$backward, , drop = FALSE])
      B <- rbind(target$points[match$forward, , drop = FALSE], target$points)
      ca <- colSums(A * weights); cb <- colSums(B * weights)
      A <- sweep(A, 2L, ca); B <- sweep(B, 2L, cb)
      s <- svd(crossprod(A * weights, B))
      Q <- s$u %*% t(s$v)
      move <- function(x) sweep(sweep(x, 2L, ca) %*% Q, 2L, cb, '+')
      proposal <- moved; proposal$points <- move(moved$points)
      next.match <- comparison.surface.matches(proposal, target)
      if (next.match$mse > match$mse) {
        # Keep the previous iterate; a roundoff-sized increase is stationarity.
        converged <- next.match$mse - match$mse <=
          tolerance * max(initial.mse, .Machine$double.eps)
        break
      }
      gain <- match$mse - next.match$mse
      moved <- proposal; vertices <- move(vertices); match <- next.match
      if (gain <= tolerance * max(initial.mse, .Machine$double.eps)) {
        converged <- TRUE
        break
      }
    }
    index <- length(runs) + 1L
    runs[[index]] <- data.frame(angle = angle, reflection = flip,
      iterations = iteration, converged = converged, rmse = sqrt(match$mse))
    if (match$mse < best$mse) best <- list(coords = vertices, mse = match$mse,
      iterations = iteration, converged = converged, start = index)
  }
  best$initial.rmse <- sqrt(initial.mse)
  best$rmse <- sqrt(best$mse)
  best$subdivisions <- subdivisions
  best$runs <- do.call(rbind, runs)
  best
}

comparison.surface.example <- function(case, fits) {
  fits <- fits[vapply(fits, function(f) !is.null(f$coords), logical(1))]
  registrations <- lapply(fits, function(f)
    comparison.surface.align(f$coords, case$X, case$triangles))
  audit <- do.call(rbind, lapply(names(fits), function(method) {
    data.frame(backend = toupper(method), comparison.surface.audit(
      fits[[method]]$coords, case$X, case$triangles, registrations[[method]]))
  }))
  # Both displays use these same limits, including all original and moved fits.
  points <- rbind(case$X, do.call(rbind, lapply(fits, function(f)
    comparison_align(f$coords, case$X))),
    do.call(rbind, lapply(registrations, `[[`, 'coords')))
  bounds <- t(apply(points, 2L, range))
  padding <- pmax(bounds[,2] - bounds[,1], 1e-8) * .04
  list(registrations = registrations, audit = audit,
    limits = bounds + cbind(-padding, padding))
}

# Assess both alignments on denser samples than those used for registration.
# This is a sampled symmetric surface distance, not a vertex-error or MDS stress.
comparison.surface.audit <- function(coords, reference, triangles, registered,
                                     subdivisions = c(6L, 9L)) {
  do.call(rbind, lapply(subdivisions, function(m) {
    target <- comparison.surface.samples(reference, triangles, m)
    score <- function(x) sqrt(comparison.surface.matches(
      comparison.surface.samples(x, triangles, m), target)$mse)
    data.frame(subdivisions = m, vertex.alignment = score(comparison_align(coords, reference)),
      surface.alignment = score(registered$coords))
  }))
}
