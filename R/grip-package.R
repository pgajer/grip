#' grip: Graph dRawing with Intelligent Placement
#'
#' Fast multiscale graph layouts in 2D and 3D.
#'
#' The package provides \code{\link{grip.layout}()} to compute graph layouts
#' from edge lists or adjacency lists, \code{\link{grip.plot}()} for quick
#' visual inspection of the resulting coordinates, and convenience graph
#' generators such as \code{\link{edges.path}()},
#' \code{\link{edges.sierpinski.triangle}()}, and
#' \code{\link{edges.sierpinski.tetrahedron}()}.
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
