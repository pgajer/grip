#' grip: Graph dRawing with Intelligent Placement
#'
#' Fast multiscale graph layouts in 2D and 3D.
#'
#' The recommended workflow is \code{\link{grip}()} with
#' \code{metric = "hop"} for topology-first layouts or
#' \code{metric = "edge_length"} when positive edge lengths define the graph
#' metric. The package also
#' provides \code{\link{score.layout}()} and
#' \code{\link{compare.layouts}()} for real-data layout selection,
#' \code{\link{trace.grip}()} for multiscale diagnostics under either metric,
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
#' @section GMDS thread control:
#' The internal compiled GMDS optimizer, also used by the GMDS exploration
#' apps, uses at most two CPU threads. With automatic selection
#' (\code{n_threads = 0}), \env{GRIP_GMDS_THREADS} can be set to \code{1} or
#' \code{2} to choose the thread count before hardware concurrency is consulted.
#' An explicit positive \code{n_threads} overrides the environment variable;
#' values above two are capped at two. Non-positive or nonnumeric environment
#' values fall back to hardware-based selection. The R optimization engine is
#' serial. This setting is not a package-wide control for GRIP layouts or for
#' third-party numerical libraries.
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
