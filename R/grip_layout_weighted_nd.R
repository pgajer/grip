grip.validate.nd.scalar.integer <- function(x,
                                            name,
                                            min = NULL,
                                            allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) return(NULL)
    stop(sprintf("%s must be supplied", name))
  }
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop(sprintf("%s must be a single finite numeric value", name))
  }
  out <- as.integer(x)
  if (is.na(out)) {
    stop(sprintf("%s must be an integer", name))
  }
  if (abs(as.double(x) - as.double(out)) > sqrt(.Machine$double.eps)) {
    stop(sprintf("%s must be an integer", name))
  }
  if (!is.null(min) && out < min) {
    stop(sprintf("%s must be >= %d", name, as.integer(min)))
  }
  out
}

grip.validate.weighted.nd.graph <- function(edges = NULL,
                                            n = NULL,
                                            adj_list = NULL,
                                            weight_list = NULL,
                                            edge_weights = NULL,
                                            caller = "grip.layout.weighted.nd") {
  if (!is.null(edges) && !is.null(adj_list)) {
    stop(sprintf("%s() accepts either edges or adj_list, not both", caller))
  }

  if (is.null(adj_list)) {
    if (is.null(edges)) {
      stop(sprintf("%s() requires edges or adj_list", caller))
    }
    edges <- as.matrix(edges)
    if (is.null(n)) {
      if (length(edges) == 0L) {
        stop("n is required when edges is empty")
      }
      n <- max(edges)
    }
    n <- grip.validate.nd.scalar.integer(n, "n", min = 1L)
    built <- grip.build.adj.from.edges(
      edges = edges,
      n = n,
      edge_weights = edge_weights
    )
    adj_list <- built$adj_list
    weight_list <- built$weight_list
  } else {
    if (!is.list(adj_list)) {
      stop("adj_list must be a list")
    }
    if (is.null(n)) {
      n <- length(adj_list)
    }
    n <- grip.validate.nd.scalar.integer(n, "n", min = 1L)
    if (length(adj_list) != n) {
      stop("adj_list length must equal n")
    }
    if (!is.null(edge_weights)) {
      stop("edge_weights is only used with edges; provide weight_list with adj_list")
    }
  }

  if (is.null(weight_list)) {
    stop(sprintf("%s() requires edge weights", caller))
  }
  if (!is.list(weight_list)) {
    stop("weight_list must be a list")
  }
  if (length(weight_list) != n) {
    stop("weight_list length must equal n")
  }

  adj_list <- lapply(seq_len(n), function(i) {
    nb <- adj_list[[i]]
    if (length(nb) == 0L) return(integer(0))
    if (!is.numeric(nb)) stop("adj_list entries must be numeric vertex ids")
    nb <- as.integer(nb)
    if (any(!is.finite(nb)) || any(nb < 1L | nb > n)) {
      stop("adj_list entries must be finite 1-based vertex ids within [1, n]")
    }
    nb
  })

  weight_list <- lapply(seq_len(n), function(i) {
    w <- weight_list[[i]]
    if (length(w) == 0L) return(numeric(0))
    if (!is.numeric(w)) stop("weight_list entries must be numeric")
    w <- as.double(w)
    if (length(w) != length(adj_list[[i]])) {
      stop("weight_list entries must be parallel to adj_list entries")
    }
    if (any(!is.finite(w) | w <= 0)) {
      stop("weight_list entries must be finite values > 0")
    }
    w
  })

  for (i in seq_len(n)) {
    if (length(adj_list[[i]]) == 0L) next
    keep <- adj_list[[i]] != i
    adj_list[[i]] <- adj_list[[i]][keep]
    weight_list[[i]] <- weight_list[[i]][keep]
  }

  list(adj_list = adj_list, weight_list = weight_list, n = n)
}

grip.validate.weighted.nd.layout.inputs <- function(edges = NULL,
                                                    n = NULL,
                                                    adj_list = NULL,
                                                    weight_list = NULL,
                                                    edge_weights = NULL,
                                                    dim = 3,
                                                    seed = 6,
                                                    length_normalization = c("median", "mean", "none"),
                                                    caller = "grip.layout.weighted.nd") {
  length_normalization <- match.arg(length_normalization)
  dim <- grip.validate.nd.scalar.integer(dim, "dim", min = 2L)
  graph <- grip.validate.weighted.nd.graph(
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights,
    caller = caller
  )

  normalized <- grip.normalize.weight_list(
    weight_list = graph$weight_list,
    mode = length_normalization
  )

  if (!is.null(seed)) {
    seed <- grip.validate.nd.scalar.integer(seed, "seed")
  }

  list(
    adj_list = graph$adj_list,
    weight_list = normalized$weight_list,
    n = graph$n,
    dim = dim,
    seed = seed,
    weight_scale = normalized$scale,
    length_normalization = length_normalization
  )
}

#' Weighted GRIP layout in arbitrary dimensions
#'
#' \code{grip.layout.weighted.nd()} is an opt-in weighted GRIP layout backend
#' for embeddings in dimensions \code{dim >= 2}. It is implemented beside the
#' legacy 2D/3D GRIP code path so that existing weighted-GRIP entry points and
#' their dimensionality checks are unchanged.
#'
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#' @param n Number of vertices.
#' @param adj_list Adjacency list (1-based) for an undirected graph.
#' @param weight_list Parallel list of strictly positive edge lengths.
#' @param edge_weights Optional vector of edge lengths for \code{edges}.
#' @param dim Embedding dimension. Must be at least 2.
#' @param preset Optional weighted layout preset: \code{"carpet"},
#'   \code{"mesh"}, \code{"cylinder"}, \code{"torus"}, \code{"sphere"},
#'   \code{"irregular"}, or \code{"tree"}.
#' @param placement Reserved for API symmetry with legacy weighted GRIP.
#'   Currently only \code{"barycenter"} is accepted.
#' @param rounds Number of weighted refinement rounds before the final phase.
#' @param final_rounds Number of final weighted refinement rounds.
#' @param num_init Coarsest-level size control. Defaults to at least
#'   \code{dim + 1}.
#' @param num_nbrs Local-neighborhood control used by the ND backend.
#' @param r Movement-rate parameter in \code{[0, 1]}.
#' @param s Repulsion scale parameter.
#' @param repulsion_factor Non-edge repulsion multiplier.
#' @param tinit_factor Initial spread multiplier.
#' @param length_normalization Global edge-length normalization:
#'   \code{"median"} (default), \code{"mean"}, or \code{"none"}.
#' @param disconnected How to handle disconnected graphs:
#'   \code{"components"} lays out each component separately and packs them;
#'   \code{"error"} rejects disconnected graphs.
#' @param seed Optional RNG seed for reproducibility. If NULL, uses current time.
#' @return Numeric matrix with one row per vertex and \code{dim} columns.
#' @export
grip.layout.weighted.nd <- function(edges = NULL,
                                    n = NULL,
                                    adj_list = NULL,
                                    weight_list = NULL,
                                    edge_weights = NULL,
                                    dim = 3,
                                    preset = NULL,
                                    placement = "barycenter",
                                    rounds = 160,
                                    final_rounds = 256,
                                    num_init = NULL,
                                    num_nbrs = 24,
                                    r = 0.03,
                                    s = 6.0,
                                    repulsion_factor = 1.5,
                                    tinit_factor = 2,
                                    length_normalization = c("median", "mean", "none"),
                                    disconnected = c("components", "error"),
                                    seed = 6) {
  disconnected <- match.arg(disconnected)
  length_normalization <- match.arg(length_normalization)
  dim <- grip.validate.nd.scalar.integer(dim, "dim", min = 2L)
  preset <- grip.normalize.weighted.preset(
    preset,
    fn = "grip.layout.weighted.nd"
  )
  resolved <- grip.resolve.weighted.preset(
    preset = preset,
    dim = dim,
    placement = placement,
    placement_missing = missing(placement),
    rounds = rounds,
    rounds_missing = missing(rounds),
    final_rounds = final_rounds,
    final_rounds_missing = missing(final_rounds),
    num_init = if (is.null(num_init)) max(24L, as.integer(dim) + 1L) else num_init,
    num_init_missing = is.null(num_init),
    num_nbrs = num_nbrs,
    num_nbrs_missing = missing(num_nbrs),
    r = r,
    r_missing = missing(r),
    s = s,
    s_missing = missing(s),
    repulsion_factor = repulsion_factor,
    repulsion_factor_missing = missing(repulsion_factor)
  )

  if (identical(preset, "tree") && missing(placement) &&
      identical(resolved$placement, "circle")) {
    resolved$placement <- "barycenter"
  }
  if (!identical(resolved$placement, "barycenter")) {
    stop("grip.layout.weighted.nd() currently supports placement = 'barycenter'")
  }

  validated <- grip.validate.weighted.nd.layout.inputs(
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights,
    dim = dim,
    seed = seed,
    length_normalization = length_normalization,
    caller = "grip.layout.weighted.nd"
  )

  rounds <- grip.validate.nd.scalar.integer(resolved$rounds, "rounds", min = 0L)
  final_rounds <- grip.validate.nd.scalar.integer(resolved$final_rounds, "final_rounds", min = 0L)
  num_init <- grip.validate.nd.scalar.integer(resolved$num_init, "num_init", min = validated$dim + 1L)
  num_nbrs <- grip.validate.nd.scalar.integer(resolved$num_nbrs, "num_nbrs", min = 1L)
  tinit_factor <- grip.validate.nd.scalar.integer(tinit_factor, "tinit_factor", min = 1L)
  if (!is.numeric(r) || length(r) != 1L || !is.finite(r) || r < 0 || r > 1) {
    stop("r must be finite and in [0, 1]")
  }
  if (!is.numeric(s) || length(s) != 1L || !is.finite(s) || s < 0) {
    stop("s must be finite and >= 0")
  }
  if (!is.numeric(repulsion_factor) || length(repulsion_factor) != 1L ||
      !is.finite(repulsion_factor) || repulsion_factor < 0) {
    stop("repulsion_factor must be finite and >= 0")
  }

  comp <- grip.connected.components(
    adj_list = validated$adj_list,
    n = validated$n
  )
  if (length(unique(comp)) > 1L && identical(disconnected, "error")) {
    stop("graph must be connected; set disconnected = 'components' to layout components separately")
  }

  layout_one <- function(adj, weights, nn) {
    z <- grip_layout_weighted_nd_adj_cpp(
      adj,
      weights,
      nn,
      validated$dim,
      rounds,
      final_rounds,
      num_init,
      num_nbrs,
      as.double(r),
      as.double(s),
      as.double(repulsion_factor),
      tinit_factor,
      validated$seed
    )
    colnames(z) <- paste0("Dim", seq_len(validated$dim))
    z
  }

  if (length(unique(comp)) == 1L) {
    return(layout_one(
      validated$adj_list,
      validated$weight_list,
      validated$n
    ))
  }

  comp.ids <- sort(unique(comp))
  layouts <- vector("list", length(comp.ids))
  for (i in seq_along(comp.ids)) {
    vertices <- which(comp == comp.ids[[i]])
    sub <- grip.induce.subgraph(
      adj_list = validated$adj_list,
      weight_list = validated$weight_list,
      vertices = vertices,
      n = validated$n
    )
    layouts[[i]] <- layout_one(
      sub$adj_list,
      sub$weight_list,
      length(vertices)
    )
  }
  out <- grip.pack.component.layouts(
    layouts = layouts,
    comp = comp,
    n = validated$n,
    dim = validated$dim
  )
  colnames(out) <- paste0("Dim", seq_len(validated$dim))
  out
}
