#' grip: Graph Embedding by Metric Multidimensional Scaling and GRIP
#'
#' Embeds graphs in Euclidean space and helps assess how well the resulting
#' coordinates preserve graph distances and local structure.
#'
#' The package offers two main ways to compute coordinates.
#' \code{\link{metric.mds}()} fits the coordinates to shortest-path distances
#' by metric multidimensional scaling, by default with stochastic gradient
#' descent; its sparse mode avoids a full matrix of pairwise distances on
#' larger graphs. This is the recommended starting point. As an alternative,
#' \code{\link{grip}()} computes a multiscale GRIP (Graph dRawing with
#' Intelligent Placement) layout in two or three dimensions. In both cases,
#' use \code{\link{score.layout}()} to assess the result.
#'
#' @section Start here:
#' Open \code{vignette("grip-examples", package = "grip")} for a small
#' complete workflow, or \code{vignette("function-guide", package = "grip")}
#' to choose functions by task. For a first embedding of your own graph, call
#' \code{metric.mds(edges = edges, n = n, dim = 3, backend = "sgd")}, and add
#' \code{approximation = "sparse"} when the graph is too large for a full
#' matrix of pairwise distances. The example below uses the alternative
#' \code{grip()} to draw a 20-vertex mesh and prints its sampled
#' graph-distance stress. A stress of
#' zero means perfect agreement, after fitting a common scale, on the pairs
#' of vertices that were sampled; it says nothing about the other pairs.
#'
#' In \code{grip()}, use \code{metric = "hop"} to lay out a graph by its
#' connectivity alone, or \code{metric = "edge_length"} when positive edge
#' lengths define the distances. See
#' \code{vignette("weighted-grip-intro", package = "grip")}.
#'
#' @section Find the next workflow:
#' \describe{
#'   \item{Choose settings and seeds}{\code{vignette("grip-real-data", package = "grip")}.}
#'   \item{Trace the construction}{\code{vignette("grip-trace-and-diagnostics", package = "grip")}.}
#'   \item{Construct graph examples}{\code{vignette("synthetic-graph-families", package = "grip")}.}
#' }
#' List available guides with \code{vignette(package = "grip")}.
#' Development installations from GitHub include these guides only when built
#' with \code{build_vignettes = TRUE}. The
#' \href{https://pgajer.github.io/grip/}{website} also provides rendered guides,
#' the graph gallery, comparisons, case study, and explorer walkthrough from
#' the same installed vignette sources. Check
#' \code{packageVersion("grip")} when using development documentation.
#'
#' @section Advanced refinement:
#' The edge-KK and geodesic-KK functions evaluate and refine layouts of
#' graphs with edge lengths. They are exported but experimental. Try the main
#' layout, scoring, and tracing functions first, and consult the function
#' guide before using them.
#'
#' @section Graph input validation:
#' Vertex indices and vertex counts must be finite integers within the range
#' of R integers. In an adjacency list, every edge must be listed from both
#' of its ends, with the same length and the same number of copies.
#' Self-loops are not accepted; remove them before calling a function. Give
#' a graph either as \code{edges} with \code{edge.weights} or as
#' \code{adj.list} with \code{weight.list}, not both. A function that takes a
#' \code{prepared} graph does not accept these inputs in addition, and an
#' explicit \code{n} must equal the number of vertices in the prepared graph.
#' If the edges, the order of the vertices, or the edge lengths change,
#' prepare the graph again.
#'
#' @section GMDS thread control:
#' The internal compiled GMDS optimizer, also used by the GMDS exploration
#' apps, uses at most two CPU threads. With automatic selection
#' (\code{n.threads = 0}), \env{GRIP_GMDS_THREADS} can be set to \code{1} or
#' \code{2} to choose the thread count before hardware concurrency is consulted.
#' An explicit positive \code{n.threads} overrides the environment variable;
#' values above two are capped at two. Non-positive or nonnumeric environment
#' values fall back to hardware-based selection. The R optimization engine is
#' serial. This setting is not a package-wide control for GRIP layouts or for
#' third-party numerical libraries.
#'
#' @examples
#' edges <- edges.mesh(4, 5)
#' coords <- grip(edges, n = 20, dim = 2, preset = "mesh", seed = 11)
#' plot.layout(coords, edges = edges, pch = 16, main = "A 4 by 5 mesh")
#' quality <- score.layout(coords, edges = edges, n = 20,
#'                         sample.size.stress = 200, stress.seed = 11,
#'                         edge.crossings = "never")
#' quality[, c("sampled.stress", "edge.length.cv")]
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
