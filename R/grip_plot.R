#' Quick plot for GRIP layouts
#'
#' @param coords Numeric matrix from \code{grip.layout()}.
#' @param edges Optional two-column matrix of edges (1-based vertex ids).
#' @param ... Additional parameters passed to the plotting function.
#' @return NULL (called for side effects).
#' @examples
#' edges <- cbind(1:5, 2:6)
#' coords <- grip.layout(edges, n = 6, dim = 2,
#'                       engine = "mish_v5",
#'                       placement = "barycenter",
#'                       rounds = 5, final_rounds = 5,
#'                       num_init = 3, num_nbrs = 4,
#'                       seed = 1)
#' grip.plot(coords, edges, main = "Path graph", pch = 16, cex = 0.8)
#' @importFrom graphics segments
#' @export
grip.plot <- function(coords, edges = NULL, ...) {
  if (!is.matrix(coords) || ncol(coords) < 2) {
    stop("coords must be a numeric matrix with at least 2 columns")
  }
  d <- ncol(coords)

  if (d == 2) {
    plot(coords[, 1], coords[, 2], asp = 1, ...)
    if (!is.null(edges) && nrow(edges) > 0) {
      apply(edges, 1, function(e) {
        graphics::segments(coords[e[1], 1], coords[e[1], 2],
                           coords[e[2], 1], coords[e[2], 2],
                           col = "gray70")
      })
    }
    return(invisible(NULL))
  }

  if (d >= 3 && requireNamespace("rgl", quietly = TRUE)) {
    rgl::plot3d(coords[, 1], coords[, 2], coords[, 3], ...)
    if (!is.null(edges) && nrow(edges) > 0) {
      apply(edges, 1, function(e) {
        rgl::segments3d(coords[e, 1], coords[e, 2], coords[e, 3],
                        col = "gray70")
      })
    }
    return(invisible(NULL))
  }

  warning("3D plotting requires the 'rgl' package; plotting first two dimensions instead")
  plot(coords[, 1], coords[, 2], asp = 1, ...)
  invisible(NULL)
}
