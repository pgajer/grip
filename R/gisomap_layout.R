#' GISOMAP layout by graph initialization followed by edge-KK polish
#'
#' `grip.layout.gisomap()` is a small orchestration wrapper for
#' geodesic-isometric graph embeddings. It initializes a weighted graph layout
#' with either metric MDS or weighted-GRIP, then polishes the coordinates with
#' [grip.optimize.edge.kk.layout()].
#'
#' The name GISOMAP is used here for "geodesic-isometric mapping" rather than
#' classical Isomap alone. Use `init = "metric_mds"` for the Isomap-like path
#' and `init = "weighted_grip"` for the multiscale weighted-GRIP path.
#'
#' @inheritParams grip.score.gmds.layout
#' @param dim Target embedding dimension. Must be at least 2.
#' @param init Initializer used before edge-KK polish.
#' @param init_args Named list of extra arguments passed to the initializer.
#'   For `init = "metric_mds"`, these are passed to
#'   [grip.metric.mds.layout()]. For `init = "weighted_grip"`, these are
#'   passed to [grip.layout.weighted()].
#' @param edge_kk_args Named list of extra arguments passed to
#'   [grip.optimize.edge.kk.layout()].
#' @param diagnostics If `TRUE`, attach the common GMDS diagnostic panel to the
#'   returned edge-KK-polished result and to the stored initializer metadata
#'   when available.
#' @param seed Random seed passed to the weighted-GRIP initializer unless
#'   overridden in `init_args`. Metric MDS does not use this seed.
#'
#' @return A `"grip_gisomap_layout"` object, also inheriting from
#'   `"grip_gmds_layout"`. The returned `coords` are the edge-KK-polished
#'   coordinates. The `metadata` field records the initializer, dimension,
#'   edge-KK method, initializer diagnostics, and edge-KK metadata.
#' @export
grip.layout.gisomap <- function(prepared = NULL,
                                edges = NULL,
                                n = NULL,
                                adj_list = NULL,
                                weight_list = NULL,
                                edge_weights = NULL,
                                dim = 2L,
                                init = c("metric_mds", "weighted_grip"),
                                init_args = list(),
                                edge_kk_args = list(),
                                diagnostics = TRUE,
                                seed = 1L) {
  init <- match.arg(init)
  dim <- grip.validate.nd.scalar.integer(dim, "dim", min = 2L)
  grip.gisomap.validate.named.list(init_args, "init_args")
  grip.gisomap.validate.named.list(edge_kk_args, "edge_kk_args")
  if (!is.logical(diagnostics) || length(diagnostics) != 1L || is.na(diagnostics)) {
    stop("diagnostics must be TRUE or FALSE")
  }

  graph_arg_names <- c(
    "prepared", "edges", "n", "adj_list", "weight_list", "edge_weights", "dim"
  )
  bad_init_args <- intersect(names(init_args), graph_arg_names)
  if (length(bad_init_args) > 0L) {
    stop(
      "init_args must not override graph inputs or dim: ",
      paste(bad_init_args, collapse = ", ")
    )
  }
  bad_edge_args <- intersect(names(edge_kk_args), c("coords", graph_arg_names))
  if (length(bad_edge_args) > 0L) {
    stop(
      "edge_kk_args must not override coords, graph inputs, or dim: ",
      paste(bad_edge_args, collapse = ", ")
    )
  }

  prepared <- grip.gmds.require.prepared(
    prepared = prepared,
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights
  )
  if (is.null(prepared$edges) || is.null(prepared$edge_targets)) {
    stop("prepared object must contain edges and edge_targets")
  }

  if (identical(init, "metric_mds")) {
    init_call <- utils::modifyList(
      list(prepared = prepared, dim = dim, diagnostics = diagnostics),
      init_args,
      keep.null = TRUE
    )
    init_fit <- do.call(grip.metric.mds.layout, init_call)
    init_coords <- init_fit$coords
    init_diagnostics <- init_fit$diagnostics
    init_metadata <- init_fit$metadata
  } else {
    weighted_call <- utils::modifyList(
      list(
        edges = prepared$edges,
        n = prepared$n,
        edge_weights = prepared$edge_targets,
        dim = dim,
        seed = seed
      ),
      init_args,
      keep.null = TRUE
    )
    init_coords <- do.call(grip.layout.weighted, weighted_call)
    init_coords <- grip.validate.coords.nd(init_coords)
    init_diagnostics <- if (isTRUE(diagnostics)) {
      grip.score.gmds.layout(coords = init_coords, prepared = prepared)
    } else {
      NULL
    }
    init_metadata <- list(engine = "weighted_grip")
  }

  edge_call <- utils::modifyList(
    list(
      coords = init_coords,
      prepared = prepared,
      dim = dim,
      diagnostics = diagnostics
    ),
    edge_kk_args,
    keep.null = TRUE
  )
  edge_fit <- do.call(grip.optimize.edge.kk.layout, edge_call)

  edge_metadata <- edge_fit$metadata
  edge_method <- edge_fit$method
  edge_fit$method <- "gisomap"
  edge_fit$metadata <- c(
    list(
      initializer = init,
      polish = "edge_kk",
      dim = dim,
      edge_kk_method = edge_method,
      initial_diagnostics = init_diagnostics,
      initial_metadata = init_metadata,
      edge_kk_metadata = edge_metadata
    ),
    edge_metadata
  )
  class(edge_fit) <- unique(c("grip_gisomap_layout", class(edge_fit)))
  edge_fit
}

grip.gisomap.validate.named.list <- function(x, name) {
  if (!is.list(x)) {
    stop(name, " must be a named list")
  }
  if (length(x) > 0L && (is.null(names(x)) || any(!nzchar(names(x))))) {
    stop(name, " must be a named list")
  }
  invisible(x)
}
