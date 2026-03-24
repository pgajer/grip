#' HMP/U01 Coarsened Giant-Component Graph
#'
#' A bundled coarsened graph derived from the giant connected component of the
#' HMP+U01 16S amplicon iKNN graph used in the HMP/U01 layout-selection study.
#' The upstream graph came from the `>=1%` relative-abundance plus PCA branch
#' with `k = 3`, after restricting to the biologically meaningful major
#' component. The graph was then greedily coarsened in two rounds from `6474`
#' vertices to `1828` supernodes.
#'
#' The object is intended for examples and vignettes on real-world 3D layout
#' comparison. It is small enough to ship with the package while preserving the
#' arm-like graph structure that motivated the tuned preset analysis.
#'
#' @docType data
#' @format A named list with components:
#' \describe{
#'   \item{adj_list}{Adjacency list of length `1828`, with 1-based integer
#'   neighbor indices.}
#'   \item{weight_list}{Parallel list of positive edge weights.}
#'   \item{vertex_data}{Data frame with one row per coarse vertex and columns
#'   `vertex_id`, `size`, `cst`, `subcst`, `dcst.depth1.absorb`,
#'   `dcst.depth1.rare`, `dcst.depth2.absorb`, `dcst.depth2.rare`, `ph`, and
#'   `log10_reads`.}
#'   \item{graph_info}{Named list summarizing provenance and graph-level
#'   metadata, including the original and coarsened vertex counts, selected
#'   `k`, edge count, and coarsening rounds.}
#' }
#'
#' @source Derived from the publicly available HMP+U01 16S amplicon analysis
#'   workflow used in the `chm_paper` project. Raw bundled artifacts and a
#'   provenance note are available under
#'   `inst/extdata/hmp_u01_gc_coarse/`.
#' @keywords datasets
#'
#' @examples
#' data(hmp.u01.gc.coarse)
#' length(hmp.u01.gc.coarse$adj_list)
#' head(hmp.u01.gc.coarse$vertex_data[, c("vertex_id", "size", "cst", "subcst")])
hmp.u01.gc.coarse <- NULL
