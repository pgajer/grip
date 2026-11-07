grip.geodesic.mds.bending.stats <- function(coords,
                                            bending.stencils = NULL,
                                            bending.weight = 0) {
  if (!is.finite(bending.weight) || bending.weight <= 0) {
    return(list(
      bending_weight = as.double(bending.weight),
      raw_penalty = 0,
      energy = 0,
      gradient = matrix(0, nrow = nrow(coords), ncol = ncol(coords))
    ))
  }
  stencils <- grip.validate.bending.stencils(
    bending.stencils,
    n = nrow(coords)
  )
  if (is.null(stencils) || nrow(stencils) == 0L) {
    stop("bending.stencils must be provided when bending.weight is positive")
  }

  residual <- coords[stencils[, 1L], , drop = FALSE] -
    2 * coords[stencils[, 2L], , drop = FALSE] +
    coords[stencils[, 3L], , drop = FALSE]
  raw.penalty <- mean(rowSums(residual^2))
  scale <- 2 * as.double(bending.weight) / nrow(stencils)
  grad <- matrix(0, nrow = nrow(coords), ncol = ncol(coords))
  for (i in seq_len(nrow(stencils))) {
    a <- stencils[i, 1L]
    b <- stencils[i, 2L]
    c <- stencils[i, 3L]
    grad[a, ] <- grad[a, ] + scale * residual[i, ]
    grad[b, ] <- grad[b, ] - 2 * scale * residual[i, ]
    grad[c, ] <- grad[c, ] + scale * residual[i, ]
  }

  list(
    bending_weight = as.double(bending.weight),
    raw_penalty = raw.penalty,
    energy = as.double(bending.weight) * raw.penalty,
    gradient = grad
  )
}
