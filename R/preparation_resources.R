#' Estimate preparation storage before allocating graph distances
#'
#' Reports a lower bound for dense distance storage, not peak memory or runtime.
#' This arithmetic-only helper does not construct a graph or allocate its cache.
#' @param n Number of vertices.
#' @param n.edges Number of undirected edges, or NA if unknown.
#' @param pair_mode Preparation mode: full paths, landmark paths, distances only,
#'   or edges only.
#' @param local_nbrs,landmark_count Counts used by landmark preparation.
#' @return A one-row data frame containing graph counts, pair mode, a pair-count
#'   upper bound, and the dense-distance storage lower bound in bytes and GiB.
#' @details One dense double matrix needs `8 * n^2` bytes. Full and landmark
#'   preparations currently compute dense distances; landmark mode reduces
#'   retained pairs, not that matrix. Trees, routes, copies, R object overhead,
#'   and solver workspace require additional memory and are not estimated.
#'   Edge-only preparation has no dense distance matrix, but still stores the
#'   graph. Use [prepare.edge.kk()] when only edge constraints are needed, or
#'   ordinary [grip()] when a full geodesic cache is unnecessary.
#'
#'   Dense preparation warns before its searches/allocations when this bound
#'   exceeds `getOption("grip.preparation.warn.bytes", 512 * 1024^2)`.
#'   Set that option to another positive byte count, or `Inf` to acknowledge and
#'   disable the advisory warning. It is not a hard memory limit. Small calls are
#'   quiet. Native refinement checks for cancellation between rounds or optimizer
#'   evaluations; worker threads finish an evaluation before interruption is
#'   raised on the main thread. Cancellation does not return a partial fit.
#' @examples
#' estimate.preparation(10000, n.edges = 20000)
#' estimate.preparation(10000, n.edges = 20000, pair_mode = "edge_only")
#' @seealso [prepare.graph.geodesic.mds()], [prepare.landmark.geodesic.kk()]
#' @md
#' @export
estimate.preparation <- function(n, n.edges = NA_integer_,
                                 pair_mode = c("all_pairs", "landmark_sparse",
                                               "distance_matrix_only", "edge_only"),
                                 local_nbrs = 20L, landmark_count = 8L) {
  n <- grip.validate.vertex.count(n)
  if (!(length(n.edges) == 1L && is.na(n.edges))) {
    n.edges <- grip.validate.resource.count(n.edges, "n.edges")
  }
  pair_mode <- match.arg(pair_mode)
  local_nbrs <- grip.validate.resource.count(local_nbrs, "local_nbrs")
  landmark_count <- grip.validate.resource.count(landmark_count, "landmark_count")
  all.pairs <- as.double(n) * (n - 1) / 2
  pairs <- switch(pair_mode, all_pairs = all.pairs,
                  landmark_sparse = min(all.pairs, as.double(n) * (as.double(local_nbrs) + landmark_count)),
                  distance_matrix_only = 0, edge_only = 0)
  bytes <- if (pair_mode == "edge_only") 0 else 8 * as.double(n)^2
  data.frame(n.vertices = n, n.edges = n.edges, pair.mode = pair_mode,
             pair.count.upper.bound = pairs, dense.distance.bytes.lower.bound = bytes,
             dense.distance.GiB.lower.bound = bytes / 1024^3)
}

grip.preflight.preparation <- function(n, n.edges = NA_integer_, pair_mode = "all_pairs") {
  estimate <- estimate.preparation(n, n.edges, pair_mode)
  budget <- getOption("grip.preparation.warn.bytes", 512 * 1024^2)
  if (!is.numeric(budget) || length(budget) != 1L || is.na(budget) || budget <= 0) {
    stop("option grip.preparation.warn.bytes must be a positive byte count or Inf")
  }
  if (estimate$dense.distance.bytes.lower.bound > budget) {
    warning(sprintf(paste0("%s preparation for %d vertices (%s edges) needs at least %.3f GiB ",
      "for one dense distance matrix, before routes, copies, and workspace. ",
      "Consider prepare.edge.kk() for edge-only constraints or grip() for ordinary layouts. ",
      "Landmark preparation still uses dense distances. Inspect estimate.preparation(); ",
      "set options(grip.preparation.warn.bytes = Inf) to acknowledge this allocation."),
      pair_mode, n, as.character(n.edges), estimate$dense.distance.GiB.lower.bound), call. = FALSE)
  }
  invisible(estimate)
}

grip.validate.resource.count <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 ||
      x > .Machine$integer.max || x != trunc(x)) {
    stop(name, " must be a single finite nonnegative integer within the R integer range")
  }
  as.integer(x)
}
