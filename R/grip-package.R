#' grip: Graph dRawing with Intelligent Placement
#'
#' Fast multiscale graph layouts in 2D and 3D.
#'
#' The recommended workflows are \code{\link{grip}()} for unweighted
#' graphs and \code{\link{weighted.grip}()} for weighted graphs. The
#' package also
#' provides \code{\link{score.layout}()} and
#' \code{\link{compare.layouts}()} for real-data layout selection,
#' \code{\link{trace.grip}()} and
#' \code{\link{trace.weighted.grip}()} for multiscale diagnostics,
#' and advanced public experimental geodesic-KK helpers such as
#' \code{\link{prepare.edge.kk}()},
#' \code{\link{prepare.geodesic.kk}()},
#' \code{\link{score.geodesic.kk}()},
#' \code{\link{prepare.landmark.geodesic.kk}()}, and
#' \code{\link{score.landmark.geodesic.kk}()} for weighted-layout
#' evaluation and refinement. Convenience graph generators such as
#' \code{\link{edges.path}()}, \code{\link{edges.sierpinski.triangle}()},
#' and \code{\link{edges.sierpinski.tetrahedron}()}, quick plotting via
#' \code{\link{plot.layout}()}, and optional Shiny explorers round out the
#' package.
#'
#' @references
#' Gajer, P. and Kobourov, S.G. (2002). GRIP: Graph dRawing with Intelligent
#' Placement. \emph{Journal of Graph Algorithms and Applications}, 6(3),
#' 203--224. doi:10.7155/jgaa.00052.
#'
#' Gajer, P., Goodrich, M.T. and Kobourov, S.G. (2004). A multi-dimensional
#' approach to force-directed layouts of large graphs.
#' \emph{Computational Geometry}, 29(1), 3--18.
#' doi:10.1016/j.comgeo.2004.03.014.
#'
#' @useDynLib grip, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
