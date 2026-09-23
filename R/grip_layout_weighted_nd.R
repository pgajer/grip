grip.validate.nd.scalar.integer <- function(x,
                                            name,
                                            min = NULL,
                                            allow.null = FALSE) {
  if (is.null(x)) {
    if (allow.null) return(NULL)
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

grip.validate.weighted.nd.graph <- function(edges = NULL, n = NULL,
                                            adj.list = NULL, weight.list = NULL,
                                            edge.weights = NULL,
                                            caller = "weighted.grip.nd") {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights)
  n <- grip.resolve.graph.n(n, edges, adj.list)
  graph <- grip.validate.layout.inputs(edges, n, adj.list, weight.list,
                                       edge.weights, dim = 2L, seed = NULL)
  if (is.null(graph$weight_list)) stop(sprintf("%s() requires edge weights", caller))
  graph[c("adj_list", "weight_list", "n")]
}

grip.validate.weighted.nd.layout.inputs <- function(edges = NULL,
                                                    n = NULL,
                                                    adj.list = NULL,
                                                    weight.list = NULL,
                                                    edge.weights = NULL,
                                                    dim = 3,
                                                    seed = 6,
                                                    length.normalization = c("median", "mean", "none"),
                                                    caller = "weighted.grip.nd") {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights)
  length.normalization <- match.arg(length.normalization)
  dim <- grip.validate.nd.scalar.integer(dim, "dim", min = 2L)
  graph <- grip.validate.weighted.nd.graph(
    edges = edges,
    n = n,
    adj.list = adj.list,
    weight.list = weight.list,
    edge.weights = edge.weights,
    caller = caller
  )

  normalized <- grip.normalize.weight.list(
    weight.list = graph$weight_list,
    mode = length.normalization
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
    length_normalization = length.normalization
  )
}

#' Weighted GRIP layout in arbitrary dimensions
#'
#' \code{weighted.grip.nd()} is an opt-in weighted GRIP layout backend
#' for embeddings in dimensions \code{dim >= 2}. It is implemented beside the
#' legacy 2D/3D GRIP code path so that existing weighted-GRIP entry points and
#' their dimensionality checks are unchanged.
#'
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#' @param n Number of vertices.
#' @param adj.list Adjacency list (1-based) for an undirected graph.
#' @param weight.list Parallel list of strictly positive edge lengths.
#' @param edge.weights Optional vector of edge lengths for \code{edges}.
#' @param dim Embedding dimension. Must be at least 2.
#' @param preset Optional weighted layout preset: \code{"carpet"},
#'   \code{"mesh"}, \code{"cylinder"}, \code{"torus"}, \code{"sphere"},
#'   \code{"irregular"}, or \code{"tree"}.
#' @param placement Initial insertion placement strategy. \code{"circle"} is
#'   only available for \code{dim = 2}; higher dimensions use
#'   \code{"barycenter"}.
#' @param rounds Number of weighted refinement rounds before the final phase.
#' @param final.rounds Number of final weighted refinement rounds.
#' @param num.init Coarsest-level size control. Defaults to at least
#'   \code{dim + 1}.
#' @param num.nbrs Local-neighborhood control used by the ND backend.
#' @param r Movement-rate parameter in \code{[0, 1]}.
#' @param s Repulsion scale parameter.
#' @param repulsion.factor Non-edge repulsion multiplier.
#' @param tinit.factor Initial spread multiplier.
#' @param final.move.scale.after.first Final-stage FR displacement multiplier
#'   applied after the first final round. Must be in \code{[0, 1]}.
#' @param final.mode Final refinement mode: \code{"fr"} or
#'   \code{"kk_repulse"}.
#' @param metric.neighbor.cap Optional cap on the number of settled Dijkstra
#'   vertices used when building weighted neighborhood caches for inserted
#'   vertices. \code{NULL} keeps the exact weighted neighborhood search.
#' @param insertion.anchor.count Number of already placed anchor vertices used
#'   when inserting non-final MISF levels.
#' @param insertion.anchor.scope Anchor eligibility rule:
#'   \code{"any_higher"} or \code{"prev_misf"}.
#' @param insertion.anchor.strategy Anchor selection rule:
#'   \code{"first"}, \code{"distance_band"}, \code{"balanced_band"}, or
#'   \code{"spread_prev"}.
#' @param level0.insertion.mode Level-0 insertion placement override:
#'   \code{"inherit"}, \code{"barycenter"}, or \code{"least_squares"}.
#' @param level0.anchor.count Number of anchors used during level-0 insertion.
#' @param level0.local.kk.steps Number of local weighted-KK polish steps used
#'   immediately after level-0 insertion.
#' @param length.normalization Global edge-length normalization:
#'   \code{"median"} (default), \code{"mean"}, or \code{"none"}.
#' @param disconnected How to handle disconnected graphs:
#'   \code{"components"} lays out each component separately and packs them;
#'   \code{"error"} rejects disconnected graphs.
#' @param seed Optional RNG seed for reproducibility. If NULL, uses current time.
#' @return Numeric matrix with one row per vertex and \code{dim} columns.
#' @export
weighted.grip.nd <- function(edges = NULL,
                                    n = NULL,
                                    adj.list = NULL,
                                    weight.list = NULL,
                                    edge.weights = NULL,
                                    dim = 3,
                                    preset = NULL,
                                    placement = c("barycenter", "circle"),
                                    rounds = 160,
                                    final.rounds = 256,
                                    num.init = NULL,
                                    num.nbrs = 24,
                                    r = 0.03,
                                    s = 6.0,
                                    repulsion.factor = 1.5,
                                    tinit.factor = 2,
                                    final.move.scale.after.first = 1,
                                    final.mode = c("fr", "kk_repulse"),
                                    metric.neighbor.cap = NULL,
                                    insertion.anchor.count = 3,
                                    insertion.anchor.scope = c("any_higher", "prev_misf"),
                                    insertion.anchor.strategy = c("first", "distance_band", "balanced_band", "spread_prev"),
                                    level0.insertion.mode = c("inherit", "barycenter", "least_squares"),
                                    level0.anchor.count = insertion.anchor.count,
                                    level0.local.kk.steps = 3,
                                    length.normalization = c("median", "mean", "none"),
                                    disconnected = c("components", "error"),
                                    seed = 6) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights)
  disconnected <- match.arg(disconnected)
  length.normalization <- match.arg(length.normalization)
  dim <- grip.validate.nd.scalar.integer(dim, "dim", min = 2L)
  preset <- grip.normalize.weighted.preset(
    preset,
    fn = "weighted.grip.nd"
  )
  resolved <- grip.resolve.weighted.preset(
    preset = preset,
    dim = dim,
    placement = placement,
    placement.missing = missing(placement),
    rounds = rounds,
    rounds.missing = missing(rounds),
    final.rounds = final.rounds,
    final.rounds.missing = missing(final.rounds),
    num.init = if (is.null(num.init)) max(24L, as.integer(dim) + 1L) else num.init,
    num.init.missing = is.null(num.init),
    num.nbrs = num.nbrs,
    num.nbrs.missing = missing(num.nbrs),
    r = r,
    r.missing = missing(r),
    s = s,
    s.missing = missing(s),
    repulsion.factor = repulsion.factor,
    repulsion.factor.missing = missing(repulsion.factor)
  )

  if (identical(preset, "tree") && missing(placement) &&
      identical(resolved$placement, "circle") && dim != 2L) {
    resolved$placement <- "barycenter"
  }
  placement <- match.arg(resolved$placement, c("barycenter", "circle"))
  if (identical(placement, "circle") && dim != 2L) {
    stop("placement = 'circle' requires dim = 2")
  }
  insertion.anchor.scope <- match.arg(insertion.anchor.scope)
  insertion.anchor.strategy <- match.arg(insertion.anchor.strategy)
  level0.insertion.mode <- match.arg(level0.insertion.mode)
  final.mode <- match.arg(final.mode)

  validated <- grip.validate.weighted.nd.layout.inputs(
    edges = edges,
    n = n,
    adj.list = adj.list,
    weight.list = weight.list,
    edge.weights = edge.weights,
    dim = dim,
    seed = seed,
    length.normalization = length.normalization,
    caller = "weighted.grip.nd"
  )

  rounds <- grip.validate.nd.scalar.integer(resolved$rounds, "rounds", min = 0L)
  final.rounds <- grip.validate.nd.scalar.integer(resolved$final_rounds, "final.rounds", min = 0L)
  num.init <- grip.validate.nd.scalar.integer(resolved$num_init, "num.init", min = validated$dim + 1L)
  num.nbrs <- grip.validate.nd.scalar.integer(resolved$num_nbrs, "num.nbrs", min = 1L)
  tinit.factor <- grip.validate.nd.scalar.integer(tinit.factor, "tinit.factor", min = 1L)
  metric.neighbor.cap <- grip.validate.weighted.metric.search.inputs(
    metric.neighbor.cap = metric.neighbor.cap,
    caller = "weighted.grip.nd"
  )
  insertion.anchor.count <- grip.validate.nd.scalar.integer(
    insertion.anchor.count, "insertion.anchor.count", min = 1L
  )
  level0.anchor.count <- grip.validate.nd.scalar.integer(
    level0.anchor.count, "level0.anchor.count", min = 1L
  )
  level0.local.kk.steps <- grip.validate.nd.scalar.integer(
    level0.local.kk.steps, "level0.local.kk.steps", min = 0L
  )
  if (!is.numeric(r) || length(r) != 1L || !is.finite(r) || r < 0 || r > 1) {
    stop("r must be finite and in [0, 1]")
  }
  if (!is.numeric(s) || length(s) != 1L || !is.finite(s) || s < 0) {
    stop("s must be finite and >= 0")
  }
  if (!is.numeric(repulsion.factor) || length(repulsion.factor) != 1L ||
      !is.finite(repulsion.factor) || repulsion.factor < 0) {
    stop("repulsion.factor must be finite and >= 0")
  }
  if (!is.numeric(final.move.scale.after.first) ||
      length(final.move.scale.after.first) != 1L ||
      !is.finite(final.move.scale.after.first) ||
      final.move.scale.after.first < 0 ||
      final.move.scale.after.first > 1) {
    stop("final.move.scale.after.first must be in [0, 1]")
  }

  comp <- grip.connected.components(
    adj.list = validated$adj_list,
    n = validated$n
  )
  if (length(unique(comp)) > 1L && identical(disconnected, "error")) {
    stop("graph must be connected; set disconnected = 'components' to layout components separately")
  }

  layout.one <- function(adj, weights, nn) {
    z <- grip_layout_weighted_nd_adj_cpp(
      adj,
      weights,
      nn,
      validated$dim,
      placement,
      rounds,
      final.rounds,
      num.init,
      num.nbrs,
      as.double(r),
      as.double(s),
      as.double(repulsion.factor),
      tinit.factor,
      as.double(final.move.scale.after.first),
      final.mode,
      if (is.null(metric.neighbor.cap)) 0L else metric.neighbor.cap,
      insertion.anchor.count,
      insertion.anchor.scope,
      insertion.anchor.strategy,
      level0.insertion.mode,
      level0.anchor.count,
      level0.local.kk.steps,
      validated$seed
    )
    colnames(z) <- paste0("Dim", seq_len(validated$dim))
    z
  }

  if (length(unique(comp)) == 1L) {
    return(layout.one(
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
      adj.list = validated$adj_list,
      weight.list = validated$weight_list,
      vertices = vertices,
      n = validated$n
    )
    layouts[[i]] <- layout.one(
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

grip.layout.weighted.nd.trace <- function(edges = NULL,
                                          n = NULL,
                                          adj.list = NULL,
                                          weight.list = NULL,
                                          edge.weights = NULL,
                                          dim = 3,
                                          preset = NULL,
                                          placement = c("barycenter", "circle"),
                                          rounds = 160,
                                          final.rounds = 256,
                                          num.init = NULL,
                                          num.nbrs = 24,
                                          r = 0.03,
                                          s = 6.0,
                                          repulsion.factor = 1.5,
                                          tinit.factor = 2,
                                          final.move.scale.after.first = 1,
                                          final.mode = c("fr", "kk_repulse"),
                                          metric.neighbor.cap = NULL,
                                          insertion.anchor.count = 3,
                                          insertion.anchor.scope = c("any_higher", "prev_misf"),
                                          insertion.anchor.strategy = c("first", "distance_band", "balanced_band", "spread_prev"),
                                          level0.insertion.mode = c("inherit", "barycenter", "least_squares"),
                                          level0.anchor.count = insertion.anchor.count,
                                          level0.local.kk.steps = 3,
                                          length.normalization = c("median", "mean", "none"),
                                          seed = 6,
                                          trace.every = 1) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights)
  length.normalization <- match.arg(length.normalization)
  dim <- grip.validate.nd.scalar.integer(dim, "dim", min = 2L)
  trace.every <- grip.validate.nd.scalar.integer(trace.every, "trace.every", min = 1L)
  preset <- grip.normalize.weighted.preset(
    preset,
    fn = "grip.layout.weighted.nd.trace"
  )
  resolved <- grip.resolve.weighted.preset(
    preset = preset,
    dim = dim,
    placement = placement,
    placement.missing = missing(placement),
    rounds = rounds,
    rounds.missing = missing(rounds),
    final.rounds = final.rounds,
    final.rounds.missing = missing(final.rounds),
    num.init = if (is.null(num.init)) max(24L, as.integer(dim) + 1L) else num.init,
    num.init.missing = is.null(num.init),
    num.nbrs = num.nbrs,
    num.nbrs.missing = missing(num.nbrs),
    r = r,
    r.missing = missing(r),
    s = s,
    s.missing = missing(s),
    repulsion.factor = repulsion.factor,
    repulsion.factor.missing = missing(repulsion.factor)
  )

  if (identical(preset, "tree") && missing(placement) &&
      identical(resolved$placement, "circle") && dim != 2L) {
    resolved$placement <- "barycenter"
  }
  placement <- match.arg(resolved$placement, c("barycenter", "circle"))
  if (identical(placement, "circle") && dim != 2L) {
    stop("placement = 'circle' requires dim = 2")
  }
  insertion.anchor.scope <- match.arg(insertion.anchor.scope)
  insertion.anchor.strategy <- match.arg(insertion.anchor.strategy)
  level0.insertion.mode <- match.arg(level0.insertion.mode)
  final.mode <- match.arg(final.mode)

  validated <- grip.validate.weighted.nd.layout.inputs(
    edges = edges,
    n = n,
    adj.list = adj.list,
    weight.list = weight.list,
    edge.weights = edge.weights,
    dim = dim,
    seed = seed,
    length.normalization = length.normalization,
    caller = "grip.layout.weighted.nd.trace"
  )

  rounds <- grip.validate.nd.scalar.integer(resolved$rounds, "rounds", min = 0L)
  final.rounds <- grip.validate.nd.scalar.integer(resolved$final_rounds, "final.rounds", min = 0L)
  num.init <- grip.validate.nd.scalar.integer(resolved$num_init, "num.init", min = validated$dim + 1L)
  num.nbrs <- grip.validate.nd.scalar.integer(resolved$num_nbrs, "num.nbrs", min = 1L)
  tinit.factor <- grip.validate.nd.scalar.integer(tinit.factor, "tinit.factor", min = 1L)
  metric.neighbor.cap <- grip.validate.weighted.metric.search.inputs(
    metric.neighbor.cap = metric.neighbor.cap,
    caller = "grip.layout.weighted.nd.trace"
  )
  insertion.anchor.count <- grip.validate.nd.scalar.integer(
    insertion.anchor.count, "insertion.anchor.count", min = 1L
  )
  level0.anchor.count <- grip.validate.nd.scalar.integer(
    level0.anchor.count, "level0.anchor.count", min = 1L
  )
  level0.local.kk.steps <- grip.validate.nd.scalar.integer(
    level0.local.kk.steps, "level0.local.kk.steps", min = 0L
  )
  if (!is.numeric(r) || length(r) != 1L || !is.finite(r) || r < 0 || r > 1) {
    stop("r must be finite and in [0, 1]")
  }
  if (!is.numeric(s) || length(s) != 1L || !is.finite(s) || s < 0) {
    stop("s must be finite and >= 0")
  }
  if (!is.numeric(repulsion.factor) || length(repulsion.factor) != 1L ||
      !is.finite(repulsion.factor) || repulsion.factor < 0) {
    stop("repulsion.factor must be finite and >= 0")
  }
  if (!is.numeric(final.move.scale.after.first) ||
      length(final.move.scale.after.first) != 1L ||
      !is.finite(final.move.scale.after.first) ||
      final.move.scale.after.first < 0 ||
      final.move.scale.after.first > 1) {
    stop("final.move.scale.after.first must be in [0, 1]")
  }

  comp <- grip.connected.components(
    adj.list = validated$adj_list,
    n = validated$n
  )
  if (length(unique(comp)) > 1L) {
    stop("grip.layout.weighted.nd.trace() currently supports only connected graphs")
  }

  out <- grip_layout_weighted_nd_trace_adj_cpp(
    validated$adj_list,
    validated$weight_list,
    validated$n,
    validated$dim,
    placement,
    rounds,
    final.rounds,
    num.init,
    num.nbrs,
    as.double(r),
    as.double(s),
    as.double(repulsion.factor),
    tinit.factor,
    as.double(final.move.scale.after.first),
    final.mode,
    if (is.null(metric.neighbor.cap)) 0L else metric.neighbor.cap,
    insertion.anchor.count,
    insertion.anchor.scope,
    insertion.anchor.strategy,
    level0.insertion.mode,
    level0.anchor.count,
    level0.local.kk.steps,
    validated$seed,
    trace.every,
    FALSE,
    -1L,
    -1L,
    -1L,
    -1L
  )
  out$final <- as.matrix(out$final)
  colnames(out$final) <- paste0("Dim", seq_len(validated$dim))
  out$frames <- lapply(out$frames, function(frame) {
    frame <- as.matrix(frame)
    colnames(frame) <- paste0("Dim", seq_len(validated$dim))
    frame
  })
  out$meta <- data.frame(
    frame = seq_along(out$frames),
    phase = as.character(out$phase),
    level_index = as.integer(out$level_index),
    misf_level = as.integer(out$misf_level),
    round_in_level = as.integer(out$round),
    active_vertices = as.integer(out$active_vertices),
    stringsAsFactors = FALSE
  )
  out$trace.every <- trace.every
  out$phase <- NULL
  out$round <- NULL
  out$level_index <- NULL
  out$misf_level <- NULL
  out$active_vertices <- NULL
  out$refinement_step_trace <- NULL
  class(out) <- c("grip_layout_weighted_nd_trace", class(out))
  out
}

grip.build.misf.weighted.nd <- function(edges = NULL,
                                        n = NULL,
                                        adj.list = NULL,
                                        weight.list = NULL,
                                        edge.weights = NULL,
                                        num.init = 24,
                                        num.nbrs = 20,
                                        length.normalization = c("median", "mean", "none"),
                                        seed = 6) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights)
  validated <- grip.validate.weighted.nd.layout.inputs(
    edges = edges,
    n = n,
    adj.list = adj.list,
    weight.list = weight.list,
    edge.weights = edge.weights,
    dim = 2,
    seed = seed,
    length.normalization = length.normalization,
    caller = "grip.build.misf.weighted.nd"
  )
  num.init <- grip.validate.nd.scalar.integer(num.init, "num.init", min = 1L)
  num.nbrs <- grip.validate.nd.scalar.integer(num.nbrs, "num.nbrs", min = 1L)

  out <- grip_build_weighted_misf_nd_adj_cpp(
    adj_list = validated$adj_list,
    weight_list = validated$weight_list,
    n = validated$n,
    num_init = as.integer(num.init),
    num_nbrs = as.integer(num.nbrs),
    seed = validated$seed
  )
  out$levels <- lapply(out$levels, as.integer)
  names(out$levels) <- names(out$levels)
  out$vertex_depth <- as.integer(out$vertex_depth)
  out$mish_order <- as.integer(out$mish_order)
  out$misf_size <- as.integer(out$misf_size)
  out$num_nbrs_schedule <- as.integer(out$num_nbrs_schedule)
  out$misf_height <- as.integer(out$misf_height)
  out$top_level_size <- as.integer(out$top_level_size)
  out$weight_scale <- validated$weight_scale
  out$length_normalization <- validated$length_normalization
  out
}
