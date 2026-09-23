#' Six Graph Examples Used by Zheng, Pawar and Goodman
#'
#' Undirected, unit-length graphs derived from the nonzero off-diagonal
#' patterns of six SuiteSparse matrices: `dwt_66`, `lesmis`, `dwt_307`,
#' `494_bus`, `dwt_1005`, and `1138_bus`. All source vertices are retained.
#' Matrix values are not used as traversal lengths. In particular, the
#' Les Miserables coappearance counts are strengths, not distances.
#'
#' @docType data
#' @format A named list of six graphs. Each graph contains:
#' \describe{
#'   \item{n}{Number of vertices, including any isolated vertices.}
#'   \item{edges}{Integer matrix with two columns; unique undirected edges,
#'     with 1-based endpoints, sorted lexicographically.}
#'   \item{edge_weights}{Unit traversal lengths parallel to the edge rows.}
#'   \item{vertex_data}{Vertex index, original matrix row index, and label.
#'     Character names are supplied for Les Miserables; other labels are indices.}
#'   \item{graph_info}{Source URL, license, original Matrix Market header,
#'     file checksum, and a description of the conversion. Derived graph names
#'     end in `_unweighted_graph` to distinguish them from the original matrices.}
#' }
#' @source SuiteSparse Matrix Collection, \url{https://sparse.tamu.edu/}.
#'   Matrices and their derived graph data are distributed under CC BY 4.0;
#'   see `extdata/zheng-graphs/PROVENANCE.md` for individual credits and sources.
#'   Original matrices, including their values and metadata, are archived there.
#' @references Zheng, J. X., Pawar, S. and Goodman, D. F. M. (2018).
#'   Graph Drawing by Stochastic Gradient Descent.
#'   \doi{10.1109/TVCG.2018.2859997}.
#'
#'   Davis, T. A. and Hu, Y. (2011). The University of Florida Sparse Matrix
#'   Collection. ACM Transactions on Mathematical Software, 38(1), Article 1.
#'   \doi{10.1145/2049662.2049663}.
#' @keywords datasets
#' @examples
#' data(zheng.graphs)
#' g <- zheng.graphs[["lesmis"]]
#' c(vertices = g$n, edges = nrow(g$edges))
#' head(g$vertex_data)
#' @md
zheng.graphs <- NULL
