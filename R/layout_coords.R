#' Extract coordinates from a layout result
#'
#' Use one extraction step for ordinary layouts, traces, MDS fits, and supported
#' refinement results. Existing return types and direct field access are unchanged.
#' @param x A numeric coordinate matrix, a `grip_layout_trace`, a
#'   `grip_gmds_layout`, or a refinement list with the documented
#'   `coords/trace/frames/prepared/score` or `coords/state/trace/frames` structure.
#' @return The coordinate matrix, unchanged, including vertex order, dimensions,
#'   dimnames, and attributes. No projection, alignment, scaling, or scoring is
#'   performed. Malformed, ambiguous, or unsupported objects produce an error.
#' @details For traces the final frame is selected. For MDS and refinement fits
#'   the fitted `coords` are selected. A plain list containing only `coords` is
#'   not a recognized result; pass that matrix explicitly instead.
#'
#'   Keep graph vertex order when reusing coordinates. MDS result fields `method`,
#'   `coords`, `diagnostics`, and `metadata` remain directly available. Stopping
#'   information in `metadata` is method-dependent. Trace `$final` and refinement
#'   `$coords` remain supported. Prepared caches and internal state fields are
#'   version-dependent: rebuild them when the graph or package contract changes.
#' @examples
#' edges <- edges.path(6)
#' fits <- list(grip(edges, n = 6, dim = 2, seed = 1),
#'              trace.grip(edges, n = 6, dim = 2, seed = 1),
#'              classical.mds(edges = edges, n = 6))
#' first <- layout.coords(fits[[1]])
#' coords <- lapply(fits, layout.coords)
#' plot.layout(coords[[1]], edges = edges)
#' vapply(coords, function(z) score.layout(z, edges = edges, n = 6,
#'     sample.size.stress = 15, stress.seed = 1,
#'     edge.crossings = "never")$sampled.stress, numeric(1))
#' @seealso [grip()], [trace.grip()], [classical.mds()], [geodesic.kk()]
#' @md
#' @export
layout.coords <- function(x) {
  if (is.matrix(x)) {
    coords <- x
  } else if (is.list(x)) {
    if (anyDuplicated(names(x)) || all(c("coords", "final") %in% names(x))) {
      stop("Ambiguous layout result: duplicate names or both coords and final are present")
    }
    if (inherits(x, "grip_layout_trace")) {
      coords <- x[["final", exact = TRUE]]
    } else if (inherits(x, "grip_gmds_layout") ||
               all(c("coords", "trace", "frames", "prepared", "score") %in% names(x)) ||
               all(c("coords", "state", "trace", "frames") %in% names(x))) {
      coords <- x[["coords", exact = TRUE]]
    } else {
      stop("Unsupported layout result; supply a coordinate matrix or a documented grip result")
    }
    if (!is.null(x$prepared$n) && is.matrix(coords) && nrow(coords) != x$prepared$n) {
      stop("Malformed layout result: coordinate rows do not match prepared$n")
    }
  } else {
    stop("Unsupported layout result; supply a coordinate matrix or a documented grip result")
  }
  if (!is.matrix(coords) || !is.numeric(coords) || ncol(coords) < 2L ||
      any(!is.finite(coords))) {
    stop("Malformed layout result: coordinates must be a finite numeric matrix with at least two columns")
  }
  coords
}
