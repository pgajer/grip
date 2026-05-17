grip.safe.correlation <- function(x, y, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  x <- as.double(x)
  y <- as.double(y)
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2L) {
    return(NA_real_)
  }
  x <- x[keep]
  y <- y[keep]
  if (stats::sd(x) <= 0 || stats::sd(y) <= 0) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x, y, method = method))
}

grip.classical.mds.embedding <- function(prepared,
                                         dim = 2L,
                                         add = FALSE,
                                         eig = TRUE) {
  prepared <- grip.validate.geodesic.mds.prepared(prepared)
  dim <- grip.validate.count(dim, "dim")
  if (dim < 2L) {
    stop("dim must be at least 2")
  }

  fit <- stats::cmdscale(
    stats::as.dist(prepared$distance_matrix),
    k = dim,
    eig = eig,
    add = add
  )

  if (is.list(fit)) {
    coords <- as.matrix(fit$points)
    eigvals <- as.double(fit$eig)
    additive.constant <- if (!is.null(fit$ac)) as.double(fit$ac) else 0
    gof <- if (!is.null(fit$GOF)) as.double(fit$GOF) else c(NA_real_, NA_real_)
  } else {
    coords <- as.matrix(fit)
    eigvals <- numeric(0L)
    additive.constant <- 0
    gof <- c(NA_real_, NA_real_)
  }

  if (ncol(coords) < dim) {
    coords <- cbind(coords, matrix(0, nrow = nrow(coords), ncol = dim - ncol(coords)))
  }
  storage.mode(coords) <- "double"

  pos.mass <- sum(pmax(eigvals, 0))
  neg.mass <- sum(abs(pmin(eigvals, 0)))
  total.mass <- pos.mass + neg.mass

  list(
    coords = coords,
    eig = eigvals,
    additive_constant = additive.constant,
    gof = gof,
    positive_eigen_fraction = if (total.mass > 0) pos.mass / total.mass else NA_real_,
    negative_eigen_fraction = if (total.mass > 0) neg.mass / total.mass else NA_real_
  )
}

grip.classical.mds.score.stats <- function(coords,
                                           prepared,
                                           distance_floor = 1e-8) {
  coords <- grip.validate.coords.nd(coords)
  prepared <- grip.validate.geodesic.mds.prepared(prepared, coords = coords)
  grip.validate.scalar(distance_floor, "distance_floor", lower = 0, open.lower = TRUE)

  g <- as.double(prepared$pair_graph_distance)
  if (length(g) == 0L) {
    return(list(
      n.pairs = 0L,
      raw_stress = NA_real_,
      stress = NA_real_,
      rmse = NA_real_,
      mean.abs.error = NA_real_,
      mean.rel.error = NA_real_,
      pearson = NA_real_,
      spearman = NA_real_,
      euclidean_distance = numeric(0L),
      target = numeric(0L),
      residual = numeric(0L),
      relative.residual = numeric(0L)
    ))
  }

  pair.matrix <- prepared$pair_matrix
  diffs <- coords[pair.matrix[, 1L], , drop = FALSE] - coords[pair.matrix[, 2L], , drop = FALSE]
  d <- sqrt(rowSums(diffs^2))
  resid <- d - g
  rel.resid <- resid / pmax(g, distance_floor)
  raw.stress <- sum(resid^2)
  denom <- sum(g^2)

  list(
    n.pairs = length(g),
    raw_stress = raw.stress,
    stress = if (is.finite(denom) && denom > 0) sqrt(raw.stress / denom) else NA_real_,
    rmse = sqrt(mean(resid^2)),
    mean.abs.error = mean(abs(resid)),
    mean.rel.error = mean(abs(rel.resid)),
    pearson = grip.safe.correlation(g, d, method = "pearson"),
    spearman = grip.safe.correlation(g, d, method = "spearman"),
    euclidean_distance = d,
    target = g,
    residual = resid,
    relative.residual = rel.resid
  )
}
