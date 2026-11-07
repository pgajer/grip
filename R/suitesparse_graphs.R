#' Bundled SuiteSparse Graphs
#'
#' Load a real graph offline, or browse the frozen package catalogue. No network,
#' Python, or external viewer is needed. The catalogue includes the curated
#' comparison gallery and all six examples used by Zheng, Pawar and Goodman.
#' Gallery membership and original-six membership are separate columns.
#'
#' @param name Full graph ID (for example `"HB/494_bus"`) or a unique short
#'   name (`"494_bus"`). Names are case-sensitive. Component-derived IDs include
#'   the `__LCC_` suffix and never resolve as the original full graph.
#' @return `suitesparse.graphs()` returns a data frame with `id`, `name`, `group`,
#'   `vertices`, `edges`, `conversion_type`, `components`, `original_vertices`,
#'   `retained_fraction`, `gallery`, and `original_six`.
#'
#'   `suitesparse.graph()` returns a list with `n` (including isolated vertices),
#'   `edges` (a two-column integer matrix of unique, lexicographically sorted,
#'   one-based undirected edges), `edge_weights` (unit traversal lengths),
#'   `vertex_data`, and `provenance`. Vertex data contain consecutive `vertex`,
#'   one-based `source_index`, `label`, `source_partition` (`row` or `column`),
#'   `original_vertex` (one-based index before component extraction), and
#'   one-based `component`, `vertex_id` (saved graph identity), and
#'   `source_vertex_id` (identity before extraction). In bipartite graphs, source indices refer separately
#'   to matrix rows and columns. Les Miserables retains character names.
#'
#'   Provenance includes source and download URLs, source matrix headers with
#'   original credits, CC BY 4.0 licensing, SHA-256 checksums, conversion recipe,
#'   original graph identity, and component retention information. The original
#'   six also retain their earlier metadata under `original_example`.
#' @details Square matrices use the union of off-diagonal numerical support
#'   with its transpose, after summing duplicates and removing numerical zeros.
#'   Rectangular matrices use distinct row and column vertices. Matrix signs,
#'   magnitudes and physical units are not traversal lengths: every edge has
#'   length one. Some full graphs are disconnected. Largest-component graphs
#'   are explicitly identified and retain their original vertex mapping.
#'
#'   The installed selection is fixed; changing viewer favorites does not change
#'   package data. See `extdata/suitesparse/selection-30_sept_2026.json` and
#'   `extdata/suitesparse/PROVENANCE.html` in the installed package.
#' @source SuiteSparse Matrix Collection, \url{https://sparse.tamu.edu/}.
#' @references Davis, T. A. and Hu, Y. (2011). The University of Florida Sparse
#'   Matrix Collection. \doi{10.1145/2049662.2049663}.
#'
#'   Zheng, J. X., Pawar, S. and Goodman, D. F. M. (2018). Graph Drawing by
#'   Stochastic Gradient Descent. \doi{10.1109/TVCG.2018.2859997}.
#' @examples
#' catalogue <- suitesparse.graphs()
#' head(catalogue)
#' g <- suitesparse.graph("494_bus")
#' identical(g, suitesparse.graph("HB/494_bus"))
#' c(vertices = g$n, edges = nrow(g$edges))
#' head(suitesparse.graph("lesmis")$vertex_data)
#' @export
suitesparse.graph <- function(name) {
  bundle <- .suitesparse_bundle()
  id <- .suitesparse_resolve(name, bundle$catalogue)
  bundle$graphs[[id]]
}

#' @rdname suitesparse.graph
#' @export
suitesparse.graphs <- function() {
  .suitesparse_bundle()$catalogue
}

.suitesparse_bundle <- local({
  bundle <- NULL
  function() {
    if (is.null(bundle)) {
      path <- system.file('extdata', 'suitesparse', 'graphs.rds', package = 'grip')
      if (!nzchar(path)) stop('Bundled SuiteSparse data are missing; reinstall grip.', call. = FALSE)
      bundle <<- readRDS(path)
    }
    bundle
  }
})

.suitesparse_resolve <- function(name, catalogue) {
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop('name must be one non-empty graph ID or short name.', call. = FALSE)
  }
  if (name %in% catalogue$id) return(name)
  candidates <- catalogue$id[catalogue$name == name]
  if (length(candidates) == 1L) return(candidates)
  if (length(candidates) > 1L) {
    stop('Ambiguous SuiteSparse graph name: ', name, '. Use a full ID: ',
         paste(candidates, collapse = ', '), call. = FALSE)
  }
  nearest <- order(utils::adist(name, catalogue$name))[seq_len(min(5L, nrow(catalogue)))]
  stop('Unknown bundled SuiteSparse graph: ', name, '. Candidate IDs: ',
       paste(catalogue$id[nearest], collapse = ', '),
       '. Use suitesparse.graphs() to browse the catalogue.', call. = FALSE)
}
