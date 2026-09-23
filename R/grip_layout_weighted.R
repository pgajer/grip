grip.normalize.weight.list <- function(weight.list,
                                       mode = c("median", "mean", "none")) {
  mode <- match.arg(mode)
  if (is.null(weight.list)) {
    stop("weight_list is required")
  }

  values <- unlist(weight.list, use.names = FALSE)
  if (length(values) == 0L) {
    return(list(weight_list = weight.list, scale = 1))
  }

  scale <- switch(
    mode,
    median = stats::median(values),
    mean = mean(values),
    none = 1
  )
  if (!is.finite(scale) || scale <= 0) {
    stop("could not derive a positive finite normalization scale from weight_list")
  }
  if (identical(mode, "none")) {
    return(list(weight_list = weight.list, scale = 1))
  }

  list(
    weight_list = lapply(weight.list, function(w) as.double(w / scale)),
    scale = as.double(scale)
  )
}

grip.normalize.weighted.preset <- function(preset, fn = "grip") {
  if (is.null(preset)) {
    return(NULL)
  }
  allowed <- c("carpet", "mesh", "cylinder", "torus", "sphere", "irregular", "tree")
  if (!is.character(preset) || length(preset) != 1L || is.na(preset)) {
    stop(sprintf(
      "preset for %s must be NULL, %s",
      fn,
      paste(sprintf("'%s'", allowed), collapse = ", ")
    ))
  }
  if (preset %in% allowed) {
    return(preset)
  }
  stop(sprintf(
    "preset for %s must be NULL, %s",
    fn,
    paste(sprintf("'%s'", allowed), collapse = ", ")
  ))
}

grip.weighted.carpet.preset.defaults <- function() {
  grip.carpet.preset.defaults()
}

grip.weighted.mesh.preset.defaults <- function() {
  grip.mesh.preset.defaults()
}

grip.weighted.cylinder.preset.defaults <- function() {
  list(
    placement = "barycenter",
    rounds = 160L,
    final_rounds = 224L,
    num_init = 14L,
    num_nbrs = 22L,
    r = 0.08,
    s = 5.8,
    repulsion_factor = 1.10
  )
}

grip.weighted.torus.preset.defaults <- function() {
  grip.torus.preset.defaults()
}

grip.weighted.sphere.preset.defaults <- function() {
  list(
    placement = "barycenter",
    rounds = 176L,
    final_rounds = 240L,
    num_init = 14L,
    num_nbrs = 24L,
    r = 0.06,
    s = 6.5,
    repulsion_factor = 0.90
  )
}

grip.weighted.irregular.preset.defaults <- function() {
  list(
    placement = "barycenter",
    rounds = 192L,
    final_rounds = 256L,
    num_init = 18L,
    num_nbrs = 24L,
    r = 0.05,
    s = 6.5,
    repulsion_factor = 1.10
  )
}

grip.weighted.tree.preset.defaults <- function(dim = 2L) {
  grip.tree.preset.defaults(dim = dim)
}

grip.resolve.weighted.preset <- function(preset,
                                         dim = 2L,
                                         placement,
                                         placement.missing,
                                         rounds,
                                         rounds.missing,
                                         final.rounds,
                                         final.rounds.missing,
                                         num.init,
                                         num.init.missing,
                                         num.nbrs,
                                         num.nbrs.missing,
                                         r,
                                         r.missing,
                                         s,
                                         s.missing,
                                         repulsion.factor,
                                         repulsion.factor.missing) {
  if (is.null(preset)) {
    return(list(
      placement = placement,
      rounds = rounds,
      final_rounds = final.rounds,
      num_init = num.init,
      num_nbrs = num.nbrs,
      r = r,
      s = s,
      repulsion_factor = repulsion.factor
    ))
  }

  defaults <- switch(
    preset,
    carpet = grip.weighted.carpet.preset.defaults(),
    mesh = grip.weighted.mesh.preset.defaults(),
    cylinder = grip.weighted.cylinder.preset.defaults(),
    torus = grip.weighted.torus.preset.defaults(),
    sphere = grip.weighted.sphere.preset.defaults(),
    irregular = grip.weighted.irregular.preset.defaults(),
    tree = grip.weighted.tree.preset.defaults(dim = dim),
    stop("unknown weighted preset")
  )

  if (placement.missing) placement <- defaults$placement
  if (rounds.missing) rounds <- defaults$rounds
  if (final.rounds.missing) final.rounds <- defaults$final_rounds
  if (num.init.missing) num.init <- defaults$num_init
  if (num.nbrs.missing) num.nbrs <- defaults$num_nbrs
  if (r.missing) r <- defaults$r
  if (s.missing) s <- defaults$s
  if (repulsion.factor.missing) repulsion.factor <- defaults$repulsion_factor

  list(
    placement = placement,
    rounds = rounds,
    final_rounds = final.rounds,
    num_init = num.init,
    num_nbrs = num.nbrs,
    r = r,
    s = s,
    repulsion_factor = repulsion.factor
  )
}

grip.validate.weighted.layout.inputs <- function(edges = NULL,
                                                 n = NULL,
                                                 adj.list = NULL,
                                                 weight.list = NULL,
                                                 edge.weights = NULL,
                                                 dim = 3,
                                                 placement = "barycenter",
                                                 seed = 6,
                                                 length.normalization = c("median", "mean", "none"),
                                                 caller = "globalrep.weighted.grip") {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights)
  length.normalization <- match.arg(length.normalization)
  validated <- grip.validate.layout.inputs(
    edges = edges,
    n = n,
    adj.list = adj.list,
    weight.list = weight.list,
    edge.weights = edge.weights,
    dim = dim,
    placement = placement,
    seed = seed
  )
  if (is.null(validated$weight_list)) {
    stop(sprintf("%s() requires edge weights", caller))
  }

  normalized <- grip.normalize.weight.list(
    weight.list = validated$weight_list,
    mode = length.normalization
  )
  validated$weight_list <- normalized$weight_list
  validated$weight_scale <- normalized$scale
  validated$length_normalization <- length.normalization
  validated
}

grip.validate.weighted.metric.search.inputs <- function(metric.neighbor.cap = NULL,
                                                        caller = "grip") {
  if (is.null(metric.neighbor.cap)) {
    return(0L)
  }
  if (!is.numeric(metric.neighbor.cap) || length(metric.neighbor.cap) != 1L || !is.finite(metric.neighbor.cap)) {
    stop(sprintf("%s() metric.neighbor.cap must be NULL or a single finite numeric value", caller))
  }
  metric.neighbor.cap <- as.integer(metric.neighbor.cap)
  if (is.na(metric.neighbor.cap) || metric.neighbor.cap <= 0L) {
    stop(sprintf("%s() metric.neighbor.cap must be a positive integer when supplied", caller))
  }
  metric.neighbor.cap
}

#' Build a weighted MISF hierarchy
#'
#' \code{build.weighted.misf()} builds the weighted max-independent-set
#' filtration used by the weighted GRIP layout core. It is primarily a
#' diagnostic and benchmarking helper for comparing weighted and combinatorial
#' hierarchies on the same graph.
#'
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#' @param n Number of vertices.
#' @param adj.list Adjacency list (1-based) for undirected graphs.
#' @param weight.list Parallel list of strictly positive edge lengths.
#' @param edge.weights Optional vector of edge lengths for \code{edges}.
#' @param num.init Number of initial vertices in the coarsest level.
#' @param num.nbrs Maximum number of retained local neighbors per level.
#' @param length.normalization Global edge-length normalization:
#'   \code{"median"} (default), \code{"mean"}, or \code{"none"}.
#' @param seed Optional RNG seed for reproducibility. If NULL, uses current time.
#' @return A list with weighted MISF levels, \code{vertex_depth},
#'   \code{mish_order}, \code{misf_size}, \code{num_nbrs_schedule},
#'   \code{misf_height}, \code{top_level_size}, \code{weight_scale}, and
#'   \code{length_normalization}.
#' @export
build.weighted.misf <- function(edges = NULL,
                                     n = NULL,
                                     adj.list = NULL,
                                     weight.list = NULL,
                                     edge.weights = NULL,
                                     num.init = 24,
                                     num.nbrs = 20,
                                     length.normalization = c("median", "mean", "none"),
                                     seed = 6) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights)
  validated <- grip.validate.weighted.layout.inputs(
    edges = edges,
    n = n,
    adj.list = adj.list,
    weight.list = weight.list,
    edge.weights = edge.weights,
    dim = 2,
    placement = "barycenter",
    seed = seed,
    length.normalization = length.normalization,
    caller = "build.weighted.misf"
  )

  out <- grip_build_weighted_misf_adj_cpp(
    adj_list = validated$adj_list,
    weight_list = validated$weight_list,
    n = validated$n,
    num_init = as.integer(num.init),
    num_nbrs = as.integer(num.nbrs),
    seed = validated$seed
  )
  out$weight_scale <- validated$weight_scale
  out$length_normalization <- validated$length_normalization
  out
}

#' Compute a weighted geometry-aware GRIP layout
#'
#' \code{globalrep.weighted.grip()} is the compatibility backend equivalent to
#' \code{\link{grip}(..., metric = "edge_length")}. It runs the edge-length
#' multiscale core, which uses edge lengths in the filtration, local
#' neighborhoods, insertion, refinement forces, and optional in-core
#' multiscale LGKK refinement.
#'
#' @inheritParams globalrep.grip
#' @param weight.list Parallel list of positive edge lengths for
#'   \code{adj.list}; required when adjacency-list input is used.
#' @param edge.weights Positive edge lengths for \code{edges}, in row order;
#'   required when edge-list input is used.
#' @param preset Optional weighted tuning preset. \code{NULL} uses the
#'   quality-first defaults for the weighted core. \code{"mesh"} targets
#'   rectangular and near-mesh weighted surfaces, \code{"cylinder"} targets
#'   cylindrical grids, \code{"torus"} targets wrapped surface grids,
#'   \code{"sphere"} targets closed near-spherical meshes,
#'   \code{"irregular"} targets irregular manifold-like weighted families,
#'   \code{"tree"} targets intrinsic weighted trees, and \code{"carpet"} keeps
#'   a high-neighborhood profile for carpet-like recursive lattices. Explicit
#'   tuning arguments override the preset field by field.
#' @param metric.neighbor.cap Optional cap on the number of settled Dijkstra
#'   vertices used when building weighted neighborhood caches for inserted
#'   vertices. \code{NULL} (default) keeps the exact weighted neighborhood
#'   search, but now stops as soon as the required weighted neighbors and
#'   anchors are filled. Supplying a positive integer enables an approximate
#'   weighted neighborhood mode for larger graphs.
#' @param length.normalization Global edge-length normalization:
#'   \code{"median"} (default), \code{"mean"}, or \code{"none"}.
#' @return A numeric matrix with \code{n} rows and \code{dim} columns.
#' @export
globalrep.weighted.grip <- function(edges = NULL,
                                           n = NULL,
                                           adj.list = NULL,
                                           weight.list = NULL,
                                           edge.weights = NULL,
                                           dim = 3,
                                           placement = c("barycenter", "circle"),
                                           preset = NULL,
                                           rounds = 160,
                                           final.rounds = 384,
                                           num.init = 24,
                                           num.nbrs = 20,
                                           r = 0.03,
                                           s = 7.5,
                                           repulsion.factor = 2.5,
                                           coarse.repulsion.factor = 1.5,
                                           coarse.repulsion.sample = 16,
                                           coarse.repulsion.exact.below = 64,
                                           final.anchor.factor = 0,
                                           final.move.scale.after.first = 1,
                                           final.mode = c("fr", "kk_repulse"),
                                           insertion.anchor.count = 3,
                                           insertion.anchor.scope = c("any_higher", "prev_misf"),
                                           insertion.anchor.strategy = c("first", "distance_band", "balanced_band", "spread_prev"),
                                           level0.insertion.mode = c("inherit", "barycenter", "least_squares"),
                                           level0.anchor.count = insertion.anchor.count,
                                           level0.local.kk.steps = 3,
                                           lgkk.polish.rounds = 0L,
                                           lgkk.multiscale.rounds = 0L,
                                           lgkk.rounds.coarse = NULL,
                                           lgkk.rounds.pre.final = NULL,
                                           lgkk.rounds.final = NULL,
                                           lgkk.local.nbrs = 20L,
                                           lgkk.landmark.count = 8L,
                                           lgkk.multiscale.scope = c("all", "coarse"),
                                           lgkk.active.limit = 4096L,
                                           metric.neighbor.cap = NULL,
                                           length.normalization = c("median", "mean", "none"),
                                           tinit.factor = 6,
                                           seed = 6,
                                           disconnected = c("components", "error")) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights)
  placement_missing <- missing(placement)
  rounds_missing <- missing(rounds)
  final_rounds_missing <- missing(final.rounds)
  num_init_missing <- missing(num.init)
  num_nbrs_missing <- missing(num.nbrs)
  r_missing <- missing(r)
  s_missing <- missing(s)
  repulsion_factor_missing <- missing(repulsion.factor)

  preset <- grip.normalize.weighted.preset(
    preset,
    fn = "globalrep.weighted.grip"
  )

  resolved <- grip.resolve.weighted.preset(
    preset = preset,
    dim = dim,
    placement = placement,
    placement.missing = placement_missing,
    rounds = rounds,
    rounds.missing = rounds_missing,
    final.rounds = final.rounds,
    final.rounds.missing = final_rounds_missing,
    num.init = num.init,
    num.init.missing = num_init_missing,
    num.nbrs = num.nbrs,
    num.nbrs.missing = num_nbrs_missing,
    r = r,
    r.missing = r_missing,
    s = s,
    s.missing = s_missing,
    repulsion.factor = repulsion.factor,
    repulsion.factor.missing = repulsion_factor_missing
  )
  placement <- resolved$placement
  rounds <- resolved$rounds
  final.rounds <- resolved$final_rounds
  num.init <- resolved$num_init
  num.nbrs <- resolved$num_nbrs
  r <- resolved$r
  s <- resolved$s
  repulsion.factor <- resolved$repulsion_factor
  placement <- match.arg(placement)
  final.mode <- match.arg(final.mode)
  insertion.anchor.scope <- match.arg(insertion.anchor.scope)
  insertion.anchor.strategy <- match.arg(insertion.anchor.strategy)
  level0.insertion.mode <- match.arg(level0.insertion.mode)
  lgkk.multiscale.scope <- match.arg(lgkk.multiscale.scope)
  disconnected <- match.arg(disconnected)

  validated <- grip.validate.weighted.layout.inputs(
    edges = edges,
    n = n,
    adj.list = adj.list,
    weight.list = weight.list,
    edge.weights = edge.weights,
    dim = dim,
    placement = placement,
    seed = seed,
    length.normalization = length.normalization,
    caller = "globalrep.weighted.grip"
  )
  adj.list <- validated$adj_list
  weight.list <- validated$weight_list
  n <- validated$n
  dim <- validated$dim
  seed <- validated$seed
  metric.neighbor.cap <- grip.validate.weighted.metric.search.inputs(
    metric.neighbor.cap = metric.neighbor.cap,
    caller = "globalrep.weighted.grip"
  )

  if (is.null(preset) && final_rounds_missing) {
    final.rounds <- grip.globalrep.default.final.rounds(n)
  }

  tuning <- grip.validate.globalrep.tuning.inputs(
    num.nbrs = num.nbrs,
    r = r,
    s = s,
    repulsion.factor = repulsion.factor,
    coarse.repulsion.factor = coarse.repulsion.factor,
    coarse.repulsion.sample = coarse.repulsion.sample,
    coarse.repulsion.exact.below = coarse.repulsion.exact.below,
    final.anchor.factor = final.anchor.factor,
    final.move.scale.after.first = final.move.scale.after.first,
    insertion.anchor.count = insertion.anchor.count,
    insertion.anchor.scope = insertion.anchor.scope,
    insertion.anchor.strategy = insertion.anchor.strategy,
    level0.insertion.mode = level0.insertion.mode,
    level0.anchor.count = level0.anchor.count,
    level0.local.kk.steps = level0.local.kk.steps
  )
  num.nbrs <- tuning$num_nbrs
  r <- tuning$r
  s <- tuning$s
  repulsion.factor <- tuning$repulsion_factor
  coarse.repulsion.factor <- tuning$coarse_repulsion_factor
  coarse.repulsion.sample <- tuning$coarse_repulsion_sample
  coarse.repulsion.exact.below <- tuning$coarse_repulsion_exact_below
  final.anchor.factor <- tuning$final_anchor_factor
  final.move.scale.after.first <- tuning$final_move_scale_after_first
  insertion.anchor.count <- tuning$insertion_anchor_count
  insertion.anchor.scope <- tuning$insertion_anchor_scope
  insertion.anchor.strategy <- tuning$insertion_anchor_strategy
  level0.insertion.mode <- tuning$level0_insertion_mode
  level0.anchor.count <- tuning$level0_anchor_count
  level0.local.kk.steps <- tuning$level0_local_kk_steps

  lgkk <- grip.validate.lgkk.polish.inputs(
    lgkk.polish.rounds = lgkk.polish.rounds,
    lgkk.multiscale.rounds = lgkk.multiscale.rounds,
    lgkk.rounds.coarse = lgkk.rounds.coarse,
    lgkk.rounds.pre.final = lgkk.rounds.pre.final,
    lgkk.rounds.final = lgkk.rounds.final,
    lgkk.local.nbrs = lgkk.local.nbrs,
    lgkk.landmark.count = lgkk.landmark.count,
    lgkk.multiscale.scope = lgkk.multiscale.scope,
    lgkk.active.limit = lgkk.active.limit
  )
  lgkk.polish.rounds <- lgkk$lgkk_polish_rounds
  lgkk.multiscale.rounds <- lgkk$lgkk_multiscale_rounds
  lgkk.rounds.coarse <- lgkk$lgkk_rounds_coarse
  lgkk.rounds.pre.final <- lgkk$lgkk_rounds_pre_final
  lgkk.rounds.final <- lgkk$lgkk_rounds_final
  lgkk.local.nbrs <- lgkk$lgkk_local_nbrs
  lgkk.landmark.count <- lgkk$lgkk_landmark_count
  lgkk.multiscale.scope <- lgkk$lgkk_multiscale_scope
  lgkk.active.limit <- lgkk$lgkk_active_limit

  layout.adj <- function(adj.list, weight.list, n) {
    coords <- grip_layout_globalrep_weighted_adj_cpp(
      adj_list = adj.list,
      weight_list = weight.list,
      n = n,
      dim = dim,
      placement = placement,
      rounds = as.integer(rounds),
      final_rounds = as.integer(final.rounds),
      num_init = as.integer(num.init),
      num_nbrs = num.nbrs,
      r = r,
      s = s,
      repulsion_factor = repulsion.factor,
      coarse_repulsion_factor = coarse.repulsion.factor,
      coarse_repulsion_sample = coarse.repulsion.sample,
      coarse_repulsion_exact_below = coarse.repulsion.exact.below,
      final_anchor_factor = final.anchor.factor,
      final_move_scale_after_first = final.move.scale.after.first,
      insertion_anchor_count = insertion.anchor.count,
      insertion_anchor_scope = insertion.anchor.scope,
      insertion_anchor_strategy = insertion.anchor.strategy,
      level0_insertion_mode = level0.insertion.mode,
      level0_anchor_count = level0.anchor.count,
      level0_local_kk_steps = level0.local.kk.steps,
      lgkk_multiscale_rounds = lgkk.multiscale.rounds,
      lgkk_rounds_coarse = lgkk.rounds.coarse,
      lgkk_rounds_pre_final = lgkk.rounds.pre.final,
      lgkk_rounds_final = lgkk.rounds.final,
      lgkk_local_nbrs = lgkk.local.nbrs,
      lgkk_landmark_count = lgkk.landmark.count,
      lgkk_multiscale_scope = lgkk.multiscale.scope,
      lgkk_active_limit = lgkk.active.limit,
      final_mode = final.mode,
      tinit_factor = as.integer(tinit.factor),
      seed = seed,
      metric_neighbor_cap = metric.neighbor.cap
    )
    polished <- grip.apply.lgkk.polish(
      coords = coords,
      adj.list = adj.list,
      weight.list = weight.list,
      rounds = lgkk.polish.rounds,
      lgkk.local.nbrs = lgkk.local.nbrs,
      lgkk.landmark.count = lgkk.landmark.count,
      return.trace = FALSE
    )
    polished$coords
  }

  comp <- grip.connected.components(adj.list = adj.list, n = n)
  n.comp <- length(unique(comp))

  if (n.comp == 1L) {
    return(layout.adj(adj.list = adj.list, weight.list = weight.list, n = n))
  }

  if (identical(disconnected, "error")) {
    stop(sprintf(
      "Input graph has %d connected components; the weighted GRIP layout core assumes connected graphs. Use disconnected = 'components' to lay out each component safely.",
      n.comp
    ))
  }

  warning(
    sprintf(
      "Input graph has %d connected components; laying out components separately to avoid disconnected-graph instability.",
      n.comp
    ),
    call. = FALSE
  )

  comp.ids <- sort(unique(comp))
  layouts <- vector("list", length(comp.ids))
  for (k in seq_along(comp.ids)) {
    rows <- which(comp == comp.ids[[k]])
    sub <- grip.induce.subgraph(
      adj.list = adj.list,
      weight.list = weight.list,
      vertices = rows,
      n = n
    )
    layouts[[k]] <- layout.adj(
      adj.list = sub$adj_list,
      weight.list = sub$weight_list,
      n = length(rows)
    )
  }

  grip.pack.component.layouts(layouts = layouts, comp = comp, n = n, dim = dim)
}

#' Internal edge-length-metric GRIP trace engine
#' @noRd
grip.trace.edge.length <- function(edges = NULL,
                                       n = NULL,
                                       adj.list = NULL,
                                       weight.list = NULL,
                                       edge.weights = NULL,
                                       dim = 3,
                                       placement = c("barycenter", "circle"),
                                       preset = NULL,
                                       rounds = 160,
                                       final.rounds = 384,
                                       num.init = 24,
                                       num.nbrs = 20,
                                       r = 0.03,
                                       s = 7.5,
                                       repulsion.factor = 2.5,
                                       coarse.repulsion.factor = 1.5,
                                       coarse.repulsion.sample = 16,
                                       coarse.repulsion.exact.below = 64,
                                       final.anchor.factor = 0,
                                       final.move.scale.after.first = 1,
                                       final.mode = c("fr", "kk_repulse"),
                                       insertion.anchor.count = 3,
                                       insertion.anchor.scope = c("any_higher", "prev_misf"),
                                       insertion.anchor.strategy = c("first", "distance_band", "balanced_band", "spread_prev"),
                                       level0.insertion.mode = c("inherit", "barycenter", "least_squares"),
                                       level0.anchor.count = insertion.anchor.count,
                                       level0.local.kk.steps = 3,
                                       lgkk.polish.rounds = 0L,
                                       lgkk.multiscale.rounds = 0L,
                                       lgkk.rounds.coarse = NULL,
                                       lgkk.rounds.pre.final = NULL,
                                       lgkk.rounds.final = NULL,
                                       lgkk.local.nbrs = 20L,
                                       lgkk.landmark.count = 8L,
                                       lgkk.multiscale.scope = c("all", "coarse"),
                                       lgkk.active.limit = 4096L,
                                       metric.neighbor.cap = NULL,
                                       length.normalization = c("median", "mean", "none"),
                                       tinit.factor = 6,
                                       seed = 6,
                                       trace = c("round", "level"),
                                       trace.every = 1,
                                       diagnostics = c("none", "light", "full"),
                                       target.coords = NULL,
                                       diagnostic.sample.size.nonedge = 1000L,
                                       diagnostic.sample.size.stress = 500L,
                                       diagnostic.nonedge.seed = 1L,
                                       diagnostic.stress.seed = 1L) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights)
  placement_missing <- missing(placement)
  rounds_missing <- missing(rounds)
  final_rounds_missing <- missing(final.rounds)
  num_init_missing <- missing(num.init)
  num_nbrs_missing <- missing(num.nbrs)
  r_missing <- missing(r)
  s_missing <- missing(s)
  repulsion_factor_missing <- missing(repulsion.factor)

  preset <- grip.normalize.weighted.preset(
    preset,
    fn = "trace.grip"
  )

  resolved <- grip.resolve.weighted.preset(
    preset = preset,
    dim = dim,
    placement = placement,
    placement.missing = placement_missing,
    rounds = rounds,
    rounds.missing = rounds_missing,
    final.rounds = final.rounds,
    final.rounds.missing = final_rounds_missing,
    num.init = num.init,
    num.init.missing = num_init_missing,
    num.nbrs = num.nbrs,
    num.nbrs.missing = num_nbrs_missing,
    r = r,
    r.missing = r_missing,
    s = s,
    s.missing = s_missing,
    repulsion.factor = repulsion.factor,
    repulsion.factor.missing = repulsion_factor_missing
  )
  placement <- resolved$placement
  rounds <- resolved$rounds
  final.rounds <- resolved$final_rounds
  num.init <- resolved$num_init
  num.nbrs <- resolved$num_nbrs
  r <- resolved$r
  s <- resolved$s
  repulsion.factor <- resolved$repulsion_factor
  placement <- match.arg(placement)
  final.mode <- match.arg(final.mode)
  insertion.anchor.scope <- match.arg(insertion.anchor.scope)
  insertion.anchor.strategy <- match.arg(insertion.anchor.strategy)
  level0.insertion.mode <- match.arg(level0.insertion.mode)
  trace <- match.arg(trace)
  diagnostics <- match.arg(diagnostics)
  lgkk.multiscale.scope <- match.arg(lgkk.multiscale.scope)

  if (!is.numeric(trace.every) || length(trace.every) != 1L || !is.finite(trace.every)) {
    stop("trace.every must be a single finite numeric value")
  }
  trace.every <- as.integer(trace.every)
  if (is.na(trace.every) || trace.every <= 0L) {
    stop("trace.every must be a positive integer")
  }

  validated <- grip.validate.weighted.layout.inputs(
    edges = edges,
    n = n,
    adj.list = adj.list,
    weight.list = weight.list,
    edge.weights = edge.weights,
    dim = dim,
    placement = placement,
    seed = seed,
    length.normalization = length.normalization,
    caller = "trace.grip"
  )
  adj.list <- validated$adj_list
  weight.list <- validated$weight_list
  n <- validated$n
  dim <- validated$dim
  seed <- validated$seed
  trace.edges <- grip.edges.from.adj.list(adj.list)
  metric.neighbor.cap <- grip.validate.weighted.metric.search.inputs(
    metric.neighbor.cap = metric.neighbor.cap,
    caller = "trace.grip"
  )

  if (is.null(preset) && final_rounds_missing) {
    final.rounds <- grip.globalrep.default.final.rounds(n)
  }

  tuning <- grip.validate.globalrep.tuning.inputs(
    num.nbrs = num.nbrs,
    r = r,
    s = s,
    repulsion.factor = repulsion.factor,
    coarse.repulsion.factor = coarse.repulsion.factor,
    coarse.repulsion.sample = coarse.repulsion.sample,
    coarse.repulsion.exact.below = coarse.repulsion.exact.below,
    final.anchor.factor = final.anchor.factor,
    final.move.scale.after.first = final.move.scale.after.first,
    insertion.anchor.count = insertion.anchor.count,
    insertion.anchor.scope = insertion.anchor.scope,
    insertion.anchor.strategy = insertion.anchor.strategy,
    level0.insertion.mode = level0.insertion.mode,
    level0.anchor.count = level0.anchor.count,
    level0.local.kk.steps = level0.local.kk.steps
  )
  num.nbrs <- tuning$num_nbrs
  r <- tuning$r
  s <- tuning$s
  repulsion.factor <- tuning$repulsion_factor
  coarse.repulsion.factor <- tuning$coarse_repulsion_factor
  coarse.repulsion.sample <- tuning$coarse_repulsion_sample
  coarse.repulsion.exact.below <- tuning$coarse_repulsion_exact_below
  final.anchor.factor <- tuning$final_anchor_factor
  final.move.scale.after.first <- tuning$final_move_scale_after_first
  insertion.anchor.count <- tuning$insertion_anchor_count
  insertion.anchor.scope <- tuning$insertion_anchor_scope
  insertion.anchor.strategy <- tuning$insertion_anchor_strategy
  level0.insertion.mode <- tuning$level0_insertion_mode
  level0.anchor.count <- tuning$level0_anchor_count
  level0.local.kk.steps <- tuning$level0_local_kk_steps

  lgkk <- grip.validate.lgkk.polish.inputs(
    lgkk.polish.rounds = lgkk.polish.rounds,
    lgkk.multiscale.rounds = lgkk.multiscale.rounds,
    lgkk.rounds.coarse = lgkk.rounds.coarse,
    lgkk.rounds.pre.final = lgkk.rounds.pre.final,
    lgkk.rounds.final = lgkk.rounds.final,
    lgkk.local.nbrs = lgkk.local.nbrs,
    lgkk.landmark.count = lgkk.landmark.count,
    lgkk.multiscale.scope = lgkk.multiscale.scope,
    lgkk.active.limit = lgkk.active.limit
  )
  lgkk.polish.rounds <- lgkk$lgkk_polish_rounds
  lgkk.multiscale.rounds <- lgkk$lgkk_multiscale_rounds
  lgkk.rounds.coarse <- lgkk$lgkk_rounds_coarse
  lgkk.rounds.pre.final <- lgkk$lgkk_rounds_pre_final
  lgkk.rounds.final <- lgkk$lgkk_rounds_final
  lgkk.local.nbrs <- lgkk$lgkk_local_nbrs
  lgkk.landmark.count <- lgkk$lgkk_landmark_count
  lgkk.multiscale.scope <- lgkk$lgkk_multiscale_scope
  lgkk.active.limit <- lgkk$lgkk_active_limit

  comp <- grip.connected.components(adj.list = adj.list, n = n)
  n.comp <- length(unique(comp))
  if (n.comp != 1L) {
    stop(sprintf(
      "trace.grip() currently supports only connected graphs; input graph has %d connected components.",
      n.comp
    ))
  }

  out <- grip_layout_globalrep_weighted_trace_adj_cpp(
    adj_list = adj.list,
    weight_list = weight.list,
    n = n,
    dim = dim,
    placement = placement,
    rounds = as.integer(rounds),
    final_rounds = as.integer(final.rounds),
    num_init = as.integer(num.init),
    num_nbrs = num.nbrs,
    r = r,
    s = s,
    repulsion_factor = repulsion.factor,
    coarse_repulsion_factor = coarse.repulsion.factor,
    coarse_repulsion_sample = coarse.repulsion.sample,
    coarse_repulsion_exact_below = coarse.repulsion.exact.below,
    final_anchor_factor = final.anchor.factor,
    final_move_scale_after_first = final.move.scale.after.first,
    insertion_anchor_count = insertion.anchor.count,
    insertion_anchor_scope = insertion.anchor.scope,
    insertion_anchor_strategy = insertion.anchor.strategy,
    level0_insertion_mode = level0.insertion.mode,
    level0_anchor_count = level0.anchor.count,
    level0_local_kk_steps = level0.local.kk.steps,
    lgkk_multiscale_rounds = lgkk.multiscale.rounds,
    lgkk_rounds_coarse = lgkk.rounds.coarse,
    lgkk_rounds_pre_final = lgkk.rounds.pre.final,
    lgkk_rounds_final = lgkk.rounds.final,
    lgkk_local_nbrs = lgkk.local.nbrs,
    lgkk_landmark_count = lgkk.landmark.count,
    lgkk_multiscale_scope = lgkk.multiscale.scope,
    lgkk_active_limit = lgkk.active.limit,
    final_mode = final.mode,
    tinit_factor = as.integer(tinit.factor),
    seed = seed,
    trace = trace,
    trace_every = trace.every,
    metric_neighbor_cap = metric.neighbor.cap,
    refinement_step_trace = FALSE,
    refinement_step_level_index = -1L,
    refinement_step_misf_level = -1L,
    refinement_step_round_start = -1L,
    refinement_step_round_end = -1L
  )
  out$refinement_step_trace <- NULL

  if (lgkk.polish.rounds > 0L) {
    polished <- grip.apply.lgkk.polish(
      coords = out$final,
      adj.list = adj.list,
      weight.list = weight.list,
      rounds = lgkk.polish.rounds,
      lgkk.local.nbrs = lgkk.local.nbrs,
      lgkk.landmark.count = lgkk.landmark.count,
      return.trace = TRUE
    )
    if (length(polished$frames) > 1L) {
      add.frames <- polished$frames[-1L]
      add.meta <- data.frame(
        frame = seq.int(nrow(out$meta) + 1L, nrow(out$meta) + length(add.frames)),
        phase = rep("lgkk", length(add.frames)),
        level_index = rep(utils::tail(out$meta$level_index, 1L), length(add.frames)),
        misf_level = rep(utils::tail(out$meta$misf_level, 1L), length(add.frames)),
        round_in_level = seq_len(length(add.frames)),
        active_vertices = rep(n, length(add.frames)),
        stringsAsFactors = FALSE
      )
      out$frames <- c(out$frames, add.frames)
      out$meta <- rbind(out$meta, add.meta)
    } else {
      out$frames[[length(out$frames)]] <- polished$coords
    }
    out$final <- polished$coords
    out$lgkk.polish <- polished$trace
  } else {
    out$lgkk.polish <- data.frame()
  }

  out$trace <- trace
  out$trace.every <- trace.every
  out$diagnostics <- grip.trace.compute.diagnostics(
    frames = out$frames,
    meta = out$meta,
    adj.list = adj.list,
    weight.list = weight.list,
    diagnostics = diagnostics,
    target.coords = target.coords,
    sample.size.nonedge = diagnostic.sample.size.nonedge,
    sample.size.stress = diagnostic.sample.size.stress,
    nonedge.seed = diagnostic.nonedge.seed,
    stress.seed = diagnostic.stress.seed
  )
  stage.bundle <- grip.layout.trace.as.stage.bundle(
    trace = out,
    edges = trace.edges
  )
  out$stage_trace <- stage.bundle$stage_trace
  out$stage_data <- stage.bundle$stage_data
  class(out) <- c("grip_layout_trace", class(out))
  out
}
