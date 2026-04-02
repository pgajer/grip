grip.validate.misf.count <- function(x, name, lower = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop(sprintf("%s must be a single finite numeric value", name))
  }
  x <- as.integer(round(x))
  if (x < lower) {
    stop(sprintf("%s must be >= %d", name, lower))
  }
  x
}

grip.complete.vertex_depth.from.levels <- function(levels, n) {
  depth <- integer(n)
  if (!length(levels)) {
    return(depth)
  }
  for (level in seq_along(levels)) {
    verts <- levels[[level]]
    if (length(verts)) {
      depth[verts] <- level - 1L
    }
  }
  depth
}

#' Build a maximal independent set filtration (MISF) for a graph
#'
#' `grip.build.misf()` exposes the maximal independent set filtration already
#' constructed internally by the GRIP layout engine. The result is graph-first:
#' you can supply either an edge list plus `n`, or an adjacency/weight-list
#' pair.
#'
#' The current implementation returns the same MISF structure used by GRIP's
#' multiscale layout core. Edge weights are accepted and passed through the
#' underlying graph object, but the current MISF construction itself is driven
#' by the graph topology and hop-distance thresholds rather than weighted
#' shortest-path distances.
#'
#' @param edges Two-column integer matrix of undirected edges (1-based vertex
#'   ids).
#' @param n Number of vertices.
#' @param adj_list Optional adjacency list (1-based integer vectors).
#' @param weight_list Optional positive edge-weight list parallel to `adj_list`.
#' @param edge_weights Optional positive vector parallel to `edges`.
#' @param num_init Target top-level active-set size used by the current GRIP
#'   MISF builder. The returned highest MISF level has size at most
#'   `min(num_init, n)`.
#' @param num_nbrs Retained local neighborhood budget per MISF level, matching
#'   the current GRIP refinement schedule metadata.
#' @param seed Optional integer seed passed to the current GRIP graph RNG.
#'
#' @return An object of class `"grip_misf"` containing:
#' \describe{
#'   \item{levels}{Named list `V0, V1, ...` of nested vertex sets (1-based).}
#'   \item{vertex_depth}{Integer vector giving the highest MISF level containing
#'   each vertex.}
#'   \item{mish_order}{The internal GRIP MISF order (1-based vertex ids).}
#'   \item{misf_size}{Sizes of the nested levels.}
#'   \item{num_nbrs_schedule}{Per-level retained neighborhood counts used by the
#'   current GRIP engine.}
#'   \item{misf_height}{Highest MISF level index.}
#'   \item{top_level_size}{Size of the highest MISF level.}
#' }
#'
#' @examples
#' edges <- edges.mesh(4, 4)
#' misf <- grip.build.misf(edges = edges, n = 16, num_init = 6, seed = 1)
#' misf$misf_size
#' misf$levels[[1L]]
#' @export
grip.build.misf <- function(edges = NULL,
                            n = NULL,
                            adj_list = NULL,
                            weight_list = NULL,
                            edge_weights = NULL,
                            num_init = 24L,
                            num_nbrs = 20L,
                            seed = 6L) {
  validated <- grip.validate.layout.inputs(
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights,
    dim = 2L,
    placement = "barycenter",
    seed = seed
  )
  num_init <- grip.validate.misf.count(num_init, "num_init", lower = 1L)
  num_nbrs <- grip.validate.misf.count(num_nbrs, "num_nbrs", lower = 1L)
  if (!is.null(validated$seed)) {
    seed <- as.integer(validated$seed)
  }

  raw <- grip_build_misf_adj_cpp(
    adj_list = validated$adj_list,
    weight_list = validated$weight_list,
    n = validated$n,
    num_init = num_init,
    num_nbrs = num_nbrs,
    seed = seed
  )

  levels <- lapply(raw$levels, as.integer)
  names(levels) <- names(raw$levels)
  vertex.depth <- as.integer(raw$vertex_depth)
  completed.depth <- grip.complete.vertex_depth.from.levels(levels, validated$n)
  if (length(vertex.depth) == length(completed.depth)) {
    vertex.depth <- pmax(vertex.depth, completed.depth)
  } else {
    vertex.depth <- completed.depth
  }

  out <- list(
    levels = levels,
    vertex_depth = vertex.depth,
    mish_order = as.integer(raw$mish_order),
    misf_size = as.integer(raw$misf_size),
    num_nbrs_schedule = as.integer(raw$num_nbrs_schedule),
    misf_height = as.integer(raw$misf_height),
    top_level_size = as.integer(raw$top_level_size),
    n = validated$n,
    num_init = as.integer(min(num_init, validated$n)),
    num_nbrs = num_nbrs,
    seed = if (is.null(seed)) NULL else as.integer(seed)
  )
  class(out) <- "grip_misf"
  out
}

print.grip_misf <- function(x, ...) {
  cat("<grip_misf>\n")
  cat("  n:", x$n, "\n")
  cat("  levels:", length(x$levels), "\n")
  cat("  sizes:", paste(x$misf_size, collapse = " -> "), "\n")
  cat("  top level size:", x$top_level_size, "\n")
  invisible(x)
}
