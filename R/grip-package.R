#' grip: Graph dRawing with Intelligent Placement
#'
#' Fast multiscale graph layouts in 2D and 3D.
#'
#' Use \code{\link{grip}()} to draw a graph, then
#' \code{\link{score.layout}()} to assess the result.
#'
#' @section Start here:
#' Open \code{vignette("grip-examples", package = "grip")} for a small
#' complete workflow, or \code{vignette("function-guide", package = "grip")}
#' to choose functions by task. The example below draws a 20-vertex mesh and
#' prints sampled graph-distance stress: zero means perfect agreement on the
#' sampled pairs after fitting a common scale, not a certificate for every pair.
#'
#' Use \code{metric = "hop"} for connectivity or
#' \code{metric = "edge_length"} when positive lengths define the graph
#' metric. See \code{vignette("weighted-grip-intro", package = "grip")}.
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
#' The public experimental edge-KK and geodesic-KK helpers provide
#' weighted-layout evaluation and refinement. Choose these through the
#' function guide after trying the main layout, scoring, and tracing workflows.
#'
#' @section Graph input validation:
#' Vertex indices and counts must be finite integers within the R integer
#' range. Adjacency must be reciprocal with matching lengths and multiplicity.
#' Self-loops are rejected; remove them explicitly. Supply either \code{edges}/\code{edge.weights} or
#' \code{adj.list}/\code{weight.list}, not both. Functions accepting a
#' \code{prepared} graph reject additional raw graph inputs; an explicit
#' \code{n} must match the stored vertex count. Rebuild the preparation when
#' changing topology, vertex order, or edge lengths.
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
