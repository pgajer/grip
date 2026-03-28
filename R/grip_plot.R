#' Quick plot for GRIP layouts
#'
#' For 2D layouts, draws vertices and edges using base graphics.
#' For 3D layouts, the behaviour depends on the \code{projection} argument:
#' \describe{
#'   \item{\code{"rgl"} (default)}{Opens an interactive \pkg{rgl} scene when
#'     that package is installed. Falls back to \code{"ortho"} with a warning
#'     if \pkg{rgl} is not available.}
#'   \item{\code{"ortho"}}{Projects the 3D coordinates to 2D with
#'     \code{\link{grip.project.3d}} and draws them with base graphics, giving
#'     a static figure suitable for vignettes and non-interactive reports.}
#' }
#'
#' @param coords Numeric matrix from \code{grip.layout()}.
#' @param edges Optional two-column matrix of edges (1-based vertex ids).
#' @param projection Character string controlling how 3D layouts are displayed:
#'   \code{"rgl"} (default) for an interactive \pkg{rgl} widget, or
#'   \code{"ortho"} for a static orthographic projection via base graphics.
#'   Ignored for 2D layouts.
#' @param azimuth Rotation around the vertical axis in degrees. Only used when
#'   \code{projection = "ortho"}. Default 35.
#' @param elevation Rotation around the horizontal axis in degrees. Only used
#'   when \code{projection = "ortho"}. Default 22.
#' @param vertex.col Colour(s) for the vertices. Recycled to match the number
#'   of vertices. Default \code{"black"}.
#' @param edge.col Colour for the edges. Default \code{"gray70"}.
#' @param ... Additional parameters passed to the underlying \code{plot()} call
#'   for 2D layouts or the \pkg{rgl} 3D plotting call.
#' @return NULL (called for side effects).
#' @examples
#' edges <- cbind(1:5, 2:6)
#' coords <- grip.layout(edges, n = 6, dim = 2,
#'                       placement = "barycenter",
#'                       rounds = 5, final_rounds = 5,
#'                       num_init = 3, num_nbrs = 4,
#'                       seed = 1)
#' grip.plot(coords, edges, main = "Path graph", pch = 16, cex = 0.8)
#' @importFrom graphics segments points
#' @export
grip.plot <- function(coords,
                      edges = NULL,
                      projection = c("rgl", "ortho"),
                      azimuth = 35,
                      elevation = 22,
                      vertex.col = "black",
                      edge.col = "gray70",
                      ...) {
  if (!is.matrix(coords) || ncol(coords) < 2) {
    stop("coords must be a numeric matrix with at least 2 columns")
  }
  projection <- match.arg(projection)
  d <- ncol(coords)

  ## ---- 2D path ----
  if (d == 2) {
    plot(coords[, 1], coords[, 2], asp = 1,
         col = vertex.col, ...)
    if (!is.null(edges) && nrow(edges) > 0) {
      apply(edges, 1, function(e) {
        graphics::segments(coords[e[1], 1], coords[e[1], 2],
                           coords[e[2], 1], coords[e[2], 2],
                           col = edge.col)
      })
    }
    return(invisible(NULL))
  }

  ## ---- 3D: orthographic static projection ----
  if (d >= 3 && projection == "ortho") {
    xy <- grip.project.3d(coords, azimuth = azimuth, elevation = elevation)
    xlim <- range(xy[, 1])
    ylim <- range(xy[, 2])
    xpad <- 0.08 * diff(xlim)
    ypad <- 0.08 * diff(ylim)
    if (!is.finite(xpad) || xpad == 0) xpad <- 0.2
    if (!is.finite(ypad) || ypad == 0) ypad <- 0.2

    plot(xy[, 1], xy[, 2],
         type = "n", asp = 1, axes = FALSE,
         xlab = "", ylab = "",
         xlim = xlim + c(-xpad, xpad),
         ylim = ylim + c(-ypad, ypad),
         ...)
    if (!is.null(edges) && nrow(edges) > 0) {
      apply(edges, 1, function(e) {
        graphics::segments(xy[e[1], 1], xy[e[1], 2],
                           xy[e[2], 1], xy[e[2], 2],
                           col = edge.col)
      })
    }
    graphics::points(xy[, 1], xy[, 2], pch = 16, cex = 0.55,
                     col = vertex.col)
    return(invisible(NULL))
  }

  ## ---- 3D: interactive rgl ----
  if (d >= 3 && requireNamespace("rgl", quietly = TRUE)) {
    rgl::plot3d(coords[, 1], coords[, 2], coords[, 3],
                col = vertex.col, ...)
    if (!is.null(edges) && nrow(edges) > 0) {
      apply(edges, 1, function(e) {
        rgl::segments3d(coords[e, 1], coords[e, 2], coords[e, 3],
                        col = edge.col)
      })
    }
    return(invisible(NULL))
  }

  warning("3D plotting requires the 'rgl' package; using orthographic ",
          "projection instead")
  grip.plot(coords, edges, projection = "ortho",
            azimuth = azimuth, elevation = elevation,
            vertex.col = vertex.col, edge.col = edge.col, ...)
  invisible(NULL)
}
