.reference.coords <- function(x, name, dim = NULL) {
  x <- as.matrix(x)
  if (!is.numeric(x) || is.complex(x) || nrow(x) < 1L || ncol(x) < 1L ||
      any(!is.finite(x)) || (!is.null(dim) && ncol(x) != dim))
    stop(name, " must be a nonempty finite numeric coordinate matrix")
  x
}

#' Score coordinates against corresponding reference observations
#'
#' Measures recovery of known coordinates, not preservation of graph distances.
#' Existing layout stability and Figure 7 diagnostics are not changed.
#' @param coords,reference Numeric matrices of identical dimensions, with rows
#'   in the same observation order. Row names do not trigger reordering.
#' @param alignment One of `"rigid"` (translation and orthogonal transformation),
#'   `"none"`, or `"similarity"` (also fit a nonnegative uniform scale).
#' @param allow_reflection Allow a reflection in the fitted orthogonal matrix.
#' @return A list with `rmse`, `relative_rmse`, `coords` (aligned coordinates),
#'   `rotation`, `translation`, `scale`, `alignment`, and `reference_radius`.
#'   RMSE is the square root of the mean squared Euclidean vertex displacement,
#'   not the mean over individual coordinate entries. Relative RMSE divides it
#'   by the reference RMS radius about its centroid; it is NA for zero radius.
#'   The transformation is `scale * coords %*% rotation + translation`.
#' @details A rigid fit preserves all distances. A similarity fit removes
#'   uniform scale as well. No vertex correspondence is estimated. A collapsed
#'   source is allowed for rigid alignment but has no identifiable similarity
#'   scale and is rejected for similarity alignment. Nonunique rotations in
#'   rank-deficient configurations can have the same minimum error.
#' @export
#' @examples
#' x <- rbind(c(0, 0), c(1, 0), c(0, 1))
#' score.coordinates(x + 2, x)$rmse
#' score.coordinates(2 * x, x, alignment = "similarity")$scale
score.coordinates <- function(coords, reference,
                              alignment = c("rigid", "none", "similarity"),
                              allow_reflection = TRUE) {
  coords <- .reference.coords(coords, "coords")
  reference <- .reference.coords(reference, "reference")
  if (!identical(dim(coords), dim(reference))) stop("Coordinate dimensions must agree")
  if (!is.logical(allow_reflection) || length(allow_reflection) != 1L || is.na(allow_reflection))
    stop("allow_reflection must be TRUE or FALSE")
  alignment <- match.arg(alignment)
  d <- ncol(coords); rotation <- diag(d); translation <- rep(0, d); scale <- 1
  center <- colMeans(reference)
  target <- sweep(reference, 2, center, "-")
  if (alignment != "none") {
    source.center <- colMeans(coords)
    source <- sweep(coords, 2, source.center, "-")
    # Scale the SVD inputs to reduce overflow and underflow risk.
    a <- max(abs(source)); b <- max(abs(target))
    if (alignment == "similarity" && a == 0) stop("A collapsed source has no identifiable similarity scale")
    if (a > 0 && b > 0) {
      sv <- svd(crossprod(source / a, target / b))
      signs <- rep(1, d)
      if (!allow_reflection && det(sv$u %*% t(sv$v)) < 0) signs[d] <- -1
      rotation <- sv$u %*% diag(signs, d) %*% t(sv$v)
      if (alignment == "similarity")
        scale <- max(0, sum(sv$d * signs)) / sum((source / a)^2) * (b / a)
    } else if (alignment == "similarity") scale <- 0
    translation <- center - as.vector(scale * source.center %*% rotation)
  }
  aligned <- sweep(scale * coords %*% rotation, 2, translation, "+")
  rmse <- sqrt(mean(rowSums((aligned - reference)^2)))
  radius <- sqrt(mean(rowSums(target^2)))
  if (!is.finite(rmse) || !is.finite(radius)) stop("Coordinates exceed the supported numerical range")
  list(rmse = rmse, relative_rmse = if (radius > 0) rmse / radius else NA_real_,
       coords = aligned, rotation = rotation, translation = translation,
       scale = scale, alignment = alignment, reference_radius = radius)
}

.reference.mesh <- function(coords, triangles, name) {
  coords <- .reference.coords(coords, name, 3L)
  if (!is.matrix(triangles) || !is.numeric(triangles) || is.complex(triangles) ||
      ncol(triangles) != 3L || nrow(triangles) < 1L || any(!is.finite(triangles)) ||
      any(triangles != floor(triangles) | triangles < 1 | triangles > nrow(coords)))
    stop(name, " triangles must contain three valid vertex indices per row")
  if (any(triangles[, 1] == triangles[, 2] | triangles[, 2] == triangles[, 3] |
          triangles[, 1] == triangles[, 3])) stop("Triangle indices must be distinct")
  if (anyDuplicated(t(apply(triangles, 1, sort)))) stop("Duplicate triangles are not supported")
  storage.mode(triangles) <- "integer"
  a <- coords[triangles[, 2], , drop = FALSE] - coords[triangles[, 1], , drop = FALSE]
  b <- coords[triangles[, 3], , drop = FALSE] - coords[triangles[, 1], , drop = FALSE]
  cross <- cbind(a[, 2]*b[, 3]-a[, 3]*b[, 2], a[, 3]*b[, 1]-a[, 1]*b[, 3],
                 a[, 1]*b[, 2]-a[, 2]*b[, 1])
  areas <- sqrt(rowSums(cross^2)) / 2
  if (any(!is.finite(areas)) || !is.finite(sum(areas)) || sum(areas) <= 0)
    stop(name, " must have finite positive total triangle area")
  if (any(areas == 0)) warning(name, " contains zero-area faces; excluded from area sampling but retained as distance targets")
  list(coords = coords, triangles = triangles, areas = areas)
}

.reference.sample <- function(mesh, n) {
  # Interleaved draws make a larger sample share the smaller sample's prefix.
  u <- matrix(stats::runif(3L * n), ncol = 3L, byrow = TRUE)
  keep <- which(mesh$areas > 0)
  cumulative <- c(0, cumsum(mesh$areas[keep] / sum(mesh$areas)))
  cumulative[length(cumulative)] <- 1
  index <- keep[pmin(length(keep), findInterval(u[, 1], cumulative))]
  t <- mesh$triangles[index, , drop = FALSE]
  r <- sqrt(u[, 2]); s <- u[, 3]
  mesh$coords[t[, 1], , drop = FALSE] * (1-r) +
    mesh$coords[t[, 2], , drop = FALSE] * (r*(1-s)) + mesh$coords[t[, 3], , drop = FALSE] * (r*s)
}

#' Estimate symmetric area-weighted distance between triangular surfaces
#'
#' Compares already aligned surfaces using closest points on triangles, not
#' nearest vertices. No registration, rescaling, or triangulation is performed.
#' @param coords,reference_coords Numeric three-column vertex matrices.
#' @param triangles,reference_triangles Three-column matrices of one-based
#'   triangle indices into the corresponding vertex matrix. Open surfaces and
#'   different triangulations are supported. Duplicate faces are rejected.
#' @param sample_size Number of independent area-uniform samples per surface,
#'   at least two. Increase this to assess Monte Carlo convergence.
#' @param seed Nonnegative integer seed, at most 2147483646. Separate fixed
#'   streams are used for the two directions; the caller's RNG state is restored.
#' @return A list with `rms`, `forward_rms`, `reverse_rms`, `forward_mean`,
#'   `reverse_mean`, `rms_mc_se`, surface areas, zero-area face counts,
#'   `sample_size`, and `seed`. Forward means coords to reference.
#' @details For surfaces A and B, squared symmetric RMS is one half of the sum
#'   of the area-normalized integrals of squared closest-point distance in the
#'   two directions. Each direction has equal weight regardless of total area.
#'   The returned score has coordinate-distance units. The Monte Carlo standard
#'   error uses a delta-method approximation; it does not quantify mesh
#'   discretization error, alignment uncertainty, or between-cloud variability.
#'   Zero observed error gives a standard error of zero, not a proof of identity.
#'   Samples are uniform over triangle area, so overlapping faces count with
#'   multiplicity. Self-intersections are not repaired or detected. Zero-area
#'   faces have no sampling mass but remain distance targets. A surface with
#'   no positive-area faces is rejected. No Hausdorff maximum is estimated.
#' @export
#' @examples
#' x <- rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0))
#' f <- matrix(1:3, nrow = 1)
#' score.surface(x, f, x, f, sample_size = 100)$rms
score.surface <- function(coords, triangles, reference_coords, reference_triangles,
                          sample_size = 5000L, seed = 1L) {
  grip.validate.scalar(sample_size, "sample_size", lower = 2, upper = 1e7)
  grip.validate.scalar(seed, "seed", lower = 0, upper = 2147483646)
  if (sample_size != floor(sample_size) || seed != floor(seed)) stop("sample_size and seed must be integers")
  source <- .reference.mesh(coords, triangles, "coords")
  target <- .reference.mesh(reference_coords, reference_triangles, "reference_coords")
  old.kind <- RNGkind()
  has.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has.seed) old.seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old.kind))
    if (has.seed) assign(".Random.seed", old.seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  })
  RNGkind("Mersenne-Twister", "Inversion")
  set.seed(seed)
  p <- .reference.sample(source, sample_size)
  set.seed(seed + 1L)
  q <- .reference.sample(target, sample_size)
  forward <- grip_surface_distances_cpp(p, target$coords, target$triangles)
  reverse <- grip_surface_distances_cpp(q, source$coords, source$triangles)
  rms <- sqrt((mean(forward^2) + mean(reverse^2)) / 2)
  if (!is.finite(rms)) stop("Surface distances exceed the supported numerical range")
  se2 <- sqrt((stats::var(forward^2) + stats::var(reverse^2)) / (4 * sample_size))
  list(rms = rms, forward_rms = sqrt(mean(forward^2)), reverse_rms = sqrt(mean(reverse^2)),
       forward_mean = mean(forward), reverse_mean = mean(reverse),
       rms_mc_se = if (rms > 0) se2 / (2*rms) else 0,
       area = sum(source$areas), reference_area = sum(target$areas),
       zero_area_faces = sum(source$areas == 0), reference_zero_area_faces = sum(target$areas == 0),
       sample_size = sample_size, seed = seed)
}
