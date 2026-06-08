#' Project 3D coordinates to 2D for static plotting
#'
#' Applies a rotation defined by azimuth and elevation angles and returns the
#' first two columns of the rotated matrix.
#'
#' This is the projection used internally by \code{\link{plot.layout}} when
#' \code{projection = "ortho"} and is exported so that users can pre-project
#' coordinates for custom plotting.
#'
#' @param coords Numeric matrix with at least 3 columns (only the first 3 are
#'   used).
#' @param azimuth Rotation around the vertical axis in degrees. Default 35.
#' @param elevation Rotation around the horizontal axis in degrees. Default 22.
#' @return A two-column numeric matrix of projected \code{(x, y)} coordinates.
#' @examples
#' edges <- edges.torus(6, 10)
#' coords3d <- grip(edges, n = max(edges), dim = 3,
#'                         preset = "torus", seed = 1)
#' xy <- project.3d(coords3d)
#' plot(xy, asp = 1, pch = 16, cex = 0.5)
#' @export
project.3d <- function(coords, azimuth = 35, elevation = 22) {
  if (!is.matrix(coords) || ncol(coords) < 3) {
    stop("coords must be a numeric matrix with at least 3 columns")
  }
  az <- azimuth * pi / 180
  el <- elevation * pi / 180
  rz <- matrix(c(
    cos(az), -sin(az), 0,
    sin(az),  cos(az), 0,
    0,        0,       1
  ), nrow = 3, byrow = TRUE)
  rx <- matrix(c(
    1, 0,       0,
    0, cos(el), -sin(el),
    0, sin(el),  cos(el)
  ), nrow = 3, byrow = TRUE)
  rotated <- coords[, 1:3, drop = FALSE] %*% rz %*% rx
  rotated[, 1:2, drop = FALSE]
}
