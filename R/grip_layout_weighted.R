grip.normalize.weight_list <- function(weight_list,
                                       mode = c("median", "mean", "none")) {
  mode <- match.arg(mode)
  if (is.null(weight_list)) {
    stop("weight_list is required")
  }

  values <- unlist(weight_list, use.names = FALSE)
  if (length(values) == 0L) {
    return(list(weight_list = weight_list, scale = 1))
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
    return(list(weight_list = weight_list, scale = 1))
  }

  list(
    weight_list = lapply(weight_list, function(w) as.double(w / scale)),
    scale = as.double(scale)
  )
}

grip.normalize.weighted.preset <- function(preset, fn = "grip.layout.weighted") {
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
                                         placement_missing,
                                         rounds,
                                         rounds_missing,
                                         final_rounds,
                                         final_rounds_missing,
                                         num_init,
                                         num_init_missing,
                                         num_nbrs,
                                         num_nbrs_missing,
                                         r,
                                         r_missing,
                                         s,
                                         s_missing,
                                         repulsion_factor,
                                         repulsion_factor_missing) {
  if (is.null(preset)) {
    return(list(
      placement = placement,
      rounds = rounds,
      final_rounds = final_rounds,
      num_init = num_init,
      num_nbrs = num_nbrs,
      r = r,
      s = s,
      repulsion_factor = repulsion_factor
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

  if (placement_missing) placement <- defaults$placement
  if (rounds_missing) rounds <- defaults$rounds
  if (final_rounds_missing) final_rounds <- defaults$final_rounds
  if (num_init_missing) num_init <- defaults$num_init
  if (num_nbrs_missing) num_nbrs <- defaults$num_nbrs
  if (r_missing) r <- defaults$r
  if (s_missing) s <- defaults$s
  if (repulsion_factor_missing) repulsion_factor <- defaults$repulsion_factor

  list(
    placement = placement,
    rounds = rounds,
    final_rounds = final_rounds,
    num_init = num_init,
    num_nbrs = num_nbrs,
    r = r,
    s = s,
    repulsion_factor = repulsion_factor
  )
}

grip.validate.weighted.layout.inputs <- function(edges = NULL,
                                                 n = NULL,
                                                 adj_list = NULL,
                                                 weight_list = NULL,
                                                 edge_weights = NULL,
                                                 dim = 3,
                                                 placement = "barycenter",
                                                 seed = 6,
                                                 length_normalization = c("median", "mean", "none"),
                                                 caller = "grip.layout.globalrep.weighted") {
  length_normalization <- match.arg(length_normalization)
  validated <- grip.validate.layout.inputs(
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights,
    dim = dim,
    placement = placement,
    seed = seed
  )
  if (is.null(validated$weight_list)) {
    stop(sprintf("%s() requires edge weights", caller))
  }

  normalized <- grip.normalize.weight_list(
    weight_list = validated$weight_list,
    mode = length_normalization
  )
  validated$weight_list <- normalized$weight_list
  validated$weight_scale <- normalized$scale
  validated$length_normalization <- length_normalization
  validated
}

grip.validate.weighted.metric.search.inputs <- function(metric_neighbor_cap = NULL,
                                                        caller = "grip.layout.weighted") {
  if (is.null(metric_neighbor_cap)) {
    return(0L)
  }
  if (!is.numeric(metric_neighbor_cap) || length(metric_neighbor_cap) != 1L || !is.finite(metric_neighbor_cap)) {
    stop(sprintf("%s() metric_neighbor_cap must be NULL or a single finite numeric value", caller))
  }
  metric_neighbor_cap <- as.integer(metric_neighbor_cap)
  if (is.na(metric_neighbor_cap) || metric_neighbor_cap <= 0L) {
    stop(sprintf("%s() metric_neighbor_cap must be a positive integer when supplied", caller))
  }
  metric_neighbor_cap
}

#' Build a weighted MISF hierarchy
#'
#' \code{grip.build.misf.weighted()} builds the weighted max-independent-set
#' filtration used by the weighted GRIP layout core. It is primarily a
#' developer and benchmarking helper for comparing weighted and combinatorial
#' hierarchies on the same graph.
#'
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#' @param n Number of vertices.
#' @param adj_list Adjacency list (1-based) for undirected graphs.
#' @param weight_list Parallel list of strictly positive edge lengths.
#' @param edge_weights Optional vector of edge lengths for \code{edges}.
#' @param num_init Number of initial vertices in the coarsest level.
#' @param num_nbrs Maximum number of retained local neighbors per level.
#' @param length_normalization Global edge-length normalization:
#'   \code{"median"} (default), \code{"mean"}, or \code{"none"}.
#' @param seed Optional RNG seed for reproducibility. If NULL, uses current time.
#' @return A list with weighted MISF levels, \code{vertex_depth},
#'   \code{mish_order}, \code{misf_size}, \code{num_nbrs_schedule},
#'   \code{misf_height}, \code{top_level_size}, \code{weight_scale}, and
#'   \code{length_normalization}.
#' @export
grip.build.misf.weighted <- function(edges = NULL,
                                     n = NULL,
                                     adj_list = NULL,
                                     weight_list = NULL,
                                     edge_weights = NULL,
                                     num_init = 24,
                                     num_nbrs = 20,
                                     length_normalization = c("median", "mean", "none"),
                                     seed = 6) {
  validated <- grip.validate.weighted.layout.inputs(
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights,
    dim = 2,
    placement = "barycenter",
    seed = seed,
    length_normalization = length_normalization,
    caller = "grip.build.misf.weighted"
  )

  out <- grip_build_weighted_misf_adj_cpp(
    adj_list = validated$adj_list,
    weight_list = validated$weight_list,
    n = validated$n,
    num_init = as.integer(num_init),
    num_nbrs = as.integer(num_nbrs),
    seed = validated$seed
  )
  out$weight_scale <- validated$weight_scale
  out$length_normalization <- validated$length_normalization
  out
}

#' Compute a weighted geometry-aware GRIP layout
#'
#' \code{grip.layout.globalrep.weighted()} is the weighted counterpart to
#' \code{\link{grip.layout.globalrep}()}. It keeps the current combinatorial
#' GRIP entry points intact and runs a separate weighted multiscale core that
#' uses edge lengths in the filtration, local neighborhoods, insertion,
#' refinement forces, and optional in-core multiscale LGKK refinement.
#'
#' @inheritParams grip.layout.globalrep
#' @param preset Optional weighted tuning preset. \code{NULL} uses the
#'   quality-first defaults for the weighted core. \code{"mesh"} targets
#'   rectangular and near-mesh weighted surfaces, \code{"cylinder"} targets
#'   cylindrical grids, \code{"torus"} targets wrapped surface grids,
#'   \code{"sphere"} targets closed near-spherical meshes,
#'   \code{"irregular"} targets irregular manifold-like weighted families,
#'   \code{"tree"} targets intrinsic weighted trees, and \code{"carpet"} keeps
#'   a high-neighborhood profile for carpet-like recursive lattices. Explicit
#'   tuning arguments override the preset field by field.
#' @param metric_neighbor_cap Optional cap on the number of settled Dijkstra
#'   vertices used when building weighted neighborhood caches for inserted
#'   vertices. \code{NULL} (default) keeps the exact weighted neighborhood
#'   search, but now stops as soon as the required weighted neighbors and
#'   anchors are filled. Supplying a positive integer enables an approximate
#'   weighted neighborhood mode for larger graphs.
#' @param length_normalization Global edge-length normalization:
#'   \code{"median"} (default), \code{"mean"}, or \code{"none"}.
#' @return A numeric matrix with \code{n} rows and \code{dim} columns.
#' @export
grip.layout.globalrep.weighted <- function(edges = NULL,
                                           n = NULL,
                                           adj_list = NULL,
                                           weight_list = NULL,
                                           edge_weights = NULL,
                                           dim = 3,
                                           placement = c("barycenter", "circle"),
                                           preset = NULL,
                                           rounds = 160,
                                           final_rounds = 384,
                                           num_init = 24,
                                           num_nbrs = 20,
                                           r = 0.03,
                                           s = 7.5,
                                           repulsion_factor = 2.5,
                                           coarse_repulsion_factor = 1.5,
                                           coarse_repulsion_sample = 16,
                                           coarse_repulsion_exact_below = 64,
                                           final_anchor_factor = 0,
                                           final_move_scale_after_first = 1,
                                           final_mode = c("fr", "kk_repulse"),
                                           insertion_anchor_count = 3,
                                           insertion_anchor_scope = c("any_higher", "prev_misf"),
                                           insertion_anchor_strategy = c("first", "distance_band", "balanced_band", "spread_prev"),
                                           level0_insertion_mode = c("inherit", "barycenter", "least_squares"),
                                           level0_anchor_count = insertion_anchor_count,
                                           level0_local_kk_steps = 3,
                                           lgkk_polish_rounds = 0L,
                                           lgkk_multiscale_rounds = 0L,
                                           lgkk_rounds_coarse = NULL,
                                           lgkk_rounds_pre_final = NULL,
                                           lgkk_rounds_final = NULL,
                                           lgkk_local_nbrs = 20L,
                                           lgkk_landmark_count = 8L,
                                           lgkk_multiscale_scope = c("all", "coarse"),
                                           lgkk_active_limit = 4096L,
                                           metric_neighbor_cap = NULL,
                                           length_normalization = c("median", "mean", "none"),
                                           tinit_factor = 6,
                                           seed = 6,
                                           disconnected = c("components", "error")) {
  placement_missing <- missing(placement)
  rounds_missing <- missing(rounds)
  final_rounds_missing <- missing(final_rounds)
  num_init_missing <- missing(num_init)
  num_nbrs_missing <- missing(num_nbrs)
  r_missing <- missing(r)
  s_missing <- missing(s)
  repulsion_factor_missing <- missing(repulsion_factor)

  preset <- grip.normalize.weighted.preset(
    preset,
    fn = "grip.layout.globalrep.weighted"
  )

  resolved <- grip.resolve.weighted.preset(
    preset = preset,
    dim = dim,
    placement = placement,
    placement_missing = placement_missing,
    rounds = rounds,
    rounds_missing = rounds_missing,
    final_rounds = final_rounds,
    final_rounds_missing = final_rounds_missing,
    num_init = num_init,
    num_init_missing = num_init_missing,
    num_nbrs = num_nbrs,
    num_nbrs_missing = num_nbrs_missing,
    r = r,
    r_missing = r_missing,
    s = s,
    s_missing = s_missing,
    repulsion_factor = repulsion_factor,
    repulsion_factor_missing = repulsion_factor_missing
  )
  placement <- resolved$placement
  rounds <- resolved$rounds
  final_rounds <- resolved$final_rounds
  num_init <- resolved$num_init
  num_nbrs <- resolved$num_nbrs
  r <- resolved$r
  s <- resolved$s
  repulsion_factor <- resolved$repulsion_factor
  placement <- match.arg(placement)
  final_mode <- match.arg(final_mode)
  insertion_anchor_scope <- match.arg(insertion_anchor_scope)
  insertion_anchor_strategy <- match.arg(insertion_anchor_strategy)
  level0_insertion_mode <- match.arg(level0_insertion_mode)
  lgkk_multiscale_scope <- match.arg(lgkk_multiscale_scope)
  disconnected <- match.arg(disconnected)

  validated <- grip.validate.weighted.layout.inputs(
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights,
    dim = dim,
    placement = placement,
    seed = seed,
    length_normalization = length_normalization,
    caller = "grip.layout.globalrep.weighted"
  )
  adj_list <- validated$adj_list
  weight_list <- validated$weight_list
  n <- validated$n
  dim <- validated$dim
  seed <- validated$seed
  metric_neighbor_cap <- grip.validate.weighted.metric.search.inputs(
    metric_neighbor_cap = metric_neighbor_cap,
    caller = "grip.layout.globalrep.weighted"
  )

  if (is.null(preset) && final_rounds_missing) {
    final_rounds <- grip.globalrep.default.final_rounds(n)
  }

  tuning <- grip.validate.globalrep.tuning.inputs(
    num_nbrs = num_nbrs,
    r = r,
    s = s,
    repulsion_factor = repulsion_factor,
    coarse_repulsion_factor = coarse_repulsion_factor,
    coarse_repulsion_sample = coarse_repulsion_sample,
    coarse_repulsion_exact_below = coarse_repulsion_exact_below,
    final_anchor_factor = final_anchor_factor,
    final_move_scale_after_first = final_move_scale_after_first,
    insertion_anchor_count = insertion_anchor_count,
    insertion_anchor_scope = insertion_anchor_scope,
    insertion_anchor_strategy = insertion_anchor_strategy,
    level0_insertion_mode = level0_insertion_mode,
    level0_anchor_count = level0_anchor_count,
    level0_local_kk_steps = level0_local_kk_steps
  )
  num_nbrs <- tuning$num_nbrs
  r <- tuning$r
  s <- tuning$s
  repulsion_factor <- tuning$repulsion_factor
  coarse_repulsion_factor <- tuning$coarse_repulsion_factor
  coarse_repulsion_sample <- tuning$coarse_repulsion_sample
  coarse_repulsion_exact_below <- tuning$coarse_repulsion_exact_below
  final_anchor_factor <- tuning$final_anchor_factor
  final_move_scale_after_first <- tuning$final_move_scale_after_first
  insertion_anchor_count <- tuning$insertion_anchor_count
  insertion_anchor_scope <- tuning$insertion_anchor_scope
  insertion_anchor_strategy <- tuning$insertion_anchor_strategy
  level0_insertion_mode <- tuning$level0_insertion_mode
  level0_anchor_count <- tuning$level0_anchor_count
  level0_local_kk_steps <- tuning$level0_local_kk_steps

  lgkk <- grip.validate.lgkk.polish.inputs(
    lgkk_polish_rounds = lgkk_polish_rounds,
    lgkk_multiscale_rounds = lgkk_multiscale_rounds,
    lgkk_rounds_coarse = lgkk_rounds_coarse,
    lgkk_rounds_pre_final = lgkk_rounds_pre_final,
    lgkk_rounds_final = lgkk_rounds_final,
    lgkk_local_nbrs = lgkk_local_nbrs,
    lgkk_landmark_count = lgkk_landmark_count,
    lgkk_multiscale_scope = lgkk_multiscale_scope,
    lgkk_active_limit = lgkk_active_limit
  )
  lgkk_polish_rounds <- lgkk$lgkk_polish_rounds
  lgkk_multiscale_rounds <- lgkk$lgkk_multiscale_rounds
  lgkk_rounds_coarse <- lgkk$lgkk_rounds_coarse
  lgkk_rounds_pre_final <- lgkk$lgkk_rounds_pre_final
  lgkk_rounds_final <- lgkk$lgkk_rounds_final
  lgkk_local_nbrs <- lgkk$lgkk_local_nbrs
  lgkk_landmark_count <- lgkk$lgkk_landmark_count
  lgkk_multiscale_scope <- lgkk$lgkk_multiscale_scope
  lgkk_active_limit <- lgkk$lgkk_active_limit

  layout.adj <- function(adj_list, weight_list, n) {
    coords <- grip_layout_globalrep_weighted_adj_cpp(
      adj_list = adj_list,
      weight_list = weight_list,
      n = n,
      dim = dim,
      placement = placement,
      rounds = as.integer(rounds),
      final_rounds = as.integer(final_rounds),
      num_init = as.integer(num_init),
      num_nbrs = num_nbrs,
      r = r,
      s = s,
      repulsion_factor = repulsion_factor,
      coarse_repulsion_factor = coarse_repulsion_factor,
      coarse_repulsion_sample = coarse_repulsion_sample,
      coarse_repulsion_exact_below = coarse_repulsion_exact_below,
      final_anchor_factor = final_anchor_factor,
      final_move_scale_after_first = final_move_scale_after_first,
      insertion_anchor_count = insertion_anchor_count,
      insertion_anchor_scope = insertion_anchor_scope,
      insertion_anchor_strategy = insertion_anchor_strategy,
      level0_insertion_mode = level0_insertion_mode,
      level0_anchor_count = level0_anchor_count,
      level0_local_kk_steps = level0_local_kk_steps,
      lgkk_multiscale_rounds = lgkk_multiscale_rounds,
      lgkk_rounds_coarse = lgkk_rounds_coarse,
      lgkk_rounds_pre_final = lgkk_rounds_pre_final,
      lgkk_rounds_final = lgkk_rounds_final,
      lgkk_local_nbrs = lgkk_local_nbrs,
      lgkk_landmark_count = lgkk_landmark_count,
      lgkk_multiscale_scope = lgkk_multiscale_scope,
      lgkk_active_limit = lgkk_active_limit,
      final_mode = final_mode,
      tinit_factor = as.integer(tinit_factor),
      seed = seed,
      metric_neighbor_cap = metric_neighbor_cap
    )
    polished <- grip.apply.lgkk.polish(
      coords = coords,
      adj_list = adj_list,
      weight_list = weight_list,
      rounds = lgkk_polish_rounds,
      lgkk_local_nbrs = lgkk_local_nbrs,
      lgkk_landmark_count = lgkk_landmark_count,
      return_trace = FALSE
    )
    polished$coords
  }

  comp <- grip.connected.components(adj_list = adj_list, n = n)
  n.comp <- length(unique(comp))

  if (n.comp == 1L) {
    return(layout.adj(adj_list = adj_list, weight_list = weight_list, n = n))
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
      adj_list = adj_list,
      weight_list = weight_list,
      vertices = rows,
      n = n
    )
    layouts[[k]] <- layout.adj(
      adj_list = sub$adj_list,
      weight_list = sub$weight_list,
      n = length(rows)
    )
  }

  grip.pack.component.layouts(layouts = layouts, comp = comp, n = n, dim = dim)
}

#' Trace a weighted geometry-aware GRIP layout
#'
#' \code{grip.layout.trace.weighted()} records the weighted GRIP layout
#' trajectory. It mirrors \code{\link{grip.layout.trace}()} but uses the
#' weighted GRIP sister core, including optional in-core multiscale LGKK
#' refinement.
#'
#' @inheritParams grip.layout.globalrep.weighted
#' @inheritParams grip.layout.trace
#' @return A list with \code{final}, \code{frames}, \code{meta}, \code{trace},
#'   \code{trace.every}, and \code{diagnostics}. \code{final} is the final
#'   coordinate matrix. \code{frames} is a list of coordinate matrices with
#'   \code{NA} rows for vertices that have not yet been introduced by GRIP.
#'   \code{meta} is a data frame describing each frame with columns
#'   \code{frame}, \code{phase}, \code{level_index}, \code{misf_level},
#'   \code{round_in_level}, and \code{active_vertices}. When diagnostics are
#'   requested, \code{diagnostics} is a data frame parallel to \code{meta} that
#'   appends per-frame quality metrics such as \code{edge.length.cv},
#'   \code{sampled.nonedge.sep.ratio}, and optional \code{procrustes.rmse}.
#' @export
grip.layout.trace.weighted <- function(edges = NULL,
                                       n = NULL,
                                       adj_list = NULL,
                                       weight_list = NULL,
                                       edge_weights = NULL,
                                       dim = 3,
                                       placement = c("barycenter", "circle"),
                                       preset = NULL,
                                       rounds = 160,
                                       final_rounds = 384,
                                       num_init = 24,
                                       num_nbrs = 20,
                                       r = 0.03,
                                       s = 7.5,
                                       repulsion_factor = 2.5,
                                       coarse_repulsion_factor = 1.5,
                                       coarse_repulsion_sample = 16,
                                       coarse_repulsion_exact_below = 64,
                                       final_anchor_factor = 0,
                                       final_move_scale_after_first = 1,
                                       final_mode = c("fr", "kk_repulse"),
                                       insertion_anchor_count = 3,
                                       insertion_anchor_scope = c("any_higher", "prev_misf"),
                                       insertion_anchor_strategy = c("first", "distance_band", "balanced_band", "spread_prev"),
                                       level0_insertion_mode = c("inherit", "barycenter", "least_squares"),
                                       level0_anchor_count = insertion_anchor_count,
                                       level0_local_kk_steps = 3,
                                       lgkk_polish_rounds = 0L,
                                       lgkk_multiscale_rounds = 0L,
                                       lgkk_rounds_coarse = NULL,
                                       lgkk_rounds_pre_final = NULL,
                                       lgkk_rounds_final = NULL,
                                       lgkk_local_nbrs = 20L,
                                       lgkk_landmark_count = 8L,
                                       lgkk_multiscale_scope = c("all", "coarse"),
                                       lgkk_active_limit = 4096L,
                                       metric_neighbor_cap = NULL,
                                       length_normalization = c("median", "mean", "none"),
                                       tinit_factor = 6,
                                       seed = 6,
                                       trace = c("round", "level"),
                                       trace.every = 1,
                                       diagnostics = c("none", "light", "full"),
                                       target_coords = NULL,
                                       diagnostic_sample_size_nonedge = 1000L,
                                       diagnostic_sample_size_stress = 500L,
                                       diagnostic_nonedge_seed = 1L,
                                       diagnostic_stress_seed = 1L) {
  placement_missing <- missing(placement)
  rounds_missing <- missing(rounds)
  final_rounds_missing <- missing(final_rounds)
  num_init_missing <- missing(num_init)
  num_nbrs_missing <- missing(num_nbrs)
  r_missing <- missing(r)
  s_missing <- missing(s)
  repulsion_factor_missing <- missing(repulsion_factor)

  preset <- grip.normalize.weighted.preset(
    preset,
    fn = "grip.layout.trace.weighted"
  )

  resolved <- grip.resolve.weighted.preset(
    preset = preset,
    dim = dim,
    placement = placement,
    placement_missing = placement_missing,
    rounds = rounds,
    rounds_missing = rounds_missing,
    final_rounds = final_rounds,
    final_rounds_missing = final_rounds_missing,
    num_init = num_init,
    num_init_missing = num_init_missing,
    num_nbrs = num_nbrs,
    num_nbrs_missing = num_nbrs_missing,
    r = r,
    r_missing = r_missing,
    s = s,
    s_missing = s_missing,
    repulsion_factor = repulsion_factor,
    repulsion_factor_missing = repulsion_factor_missing
  )
  placement <- resolved$placement
  rounds <- resolved$rounds
  final_rounds <- resolved$final_rounds
  num_init <- resolved$num_init
  num_nbrs <- resolved$num_nbrs
  r <- resolved$r
  s <- resolved$s
  repulsion_factor <- resolved$repulsion_factor
  placement <- match.arg(placement)
  final_mode <- match.arg(final_mode)
  insertion_anchor_scope <- match.arg(insertion_anchor_scope)
  insertion_anchor_strategy <- match.arg(insertion_anchor_strategy)
  level0_insertion_mode <- match.arg(level0_insertion_mode)
  trace <- match.arg(trace)
  diagnostics <- match.arg(diagnostics)
  lgkk_multiscale_scope <- match.arg(lgkk_multiscale_scope)

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
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights,
    dim = dim,
    placement = placement,
    seed = seed,
    length_normalization = length_normalization,
    caller = "grip.layout.trace.weighted"
  )
  adj_list <- validated$adj_list
  weight_list <- validated$weight_list
  n <- validated$n
  dim <- validated$dim
  seed <- validated$seed
  metric_neighbor_cap <- grip.validate.weighted.metric.search.inputs(
    metric_neighbor_cap = metric_neighbor_cap,
    caller = "grip.layout.trace.weighted"
  )

  if (is.null(preset) && final_rounds_missing) {
    final_rounds <- grip.globalrep.default.final_rounds(n)
  }

  tuning <- grip.validate.globalrep.tuning.inputs(
    num_nbrs = num_nbrs,
    r = r,
    s = s,
    repulsion_factor = repulsion_factor,
    coarse_repulsion_factor = coarse_repulsion_factor,
    coarse_repulsion_sample = coarse_repulsion_sample,
    coarse_repulsion_exact_below = coarse_repulsion_exact_below,
    final_anchor_factor = final_anchor_factor,
    final_move_scale_after_first = final_move_scale_after_first,
    insertion_anchor_count = insertion_anchor_count,
    insertion_anchor_scope = insertion_anchor_scope,
    insertion_anchor_strategy = insertion_anchor_strategy,
    level0_insertion_mode = level0_insertion_mode,
    level0_anchor_count = level0_anchor_count,
    level0_local_kk_steps = level0_local_kk_steps
  )
  num_nbrs <- tuning$num_nbrs
  r <- tuning$r
  s <- tuning$s
  repulsion_factor <- tuning$repulsion_factor
  coarse_repulsion_factor <- tuning$coarse_repulsion_factor
  coarse_repulsion_sample <- tuning$coarse_repulsion_sample
  coarse_repulsion_exact_below <- tuning$coarse_repulsion_exact_below
  final_anchor_factor <- tuning$final_anchor_factor
  final_move_scale_after_first <- tuning$final_move_scale_after_first
  insertion_anchor_count <- tuning$insertion_anchor_count
  insertion_anchor_scope <- tuning$insertion_anchor_scope
  insertion_anchor_strategy <- tuning$insertion_anchor_strategy
  level0_insertion_mode <- tuning$level0_insertion_mode
  level0_anchor_count <- tuning$level0_anchor_count
  level0_local_kk_steps <- tuning$level0_local_kk_steps

  lgkk <- grip.validate.lgkk.polish.inputs(
    lgkk_polish_rounds = lgkk_polish_rounds,
    lgkk_multiscale_rounds = lgkk_multiscale_rounds,
    lgkk_rounds_coarse = lgkk_rounds_coarse,
    lgkk_rounds_pre_final = lgkk_rounds_pre_final,
    lgkk_rounds_final = lgkk_rounds_final,
    lgkk_local_nbrs = lgkk_local_nbrs,
    lgkk_landmark_count = lgkk_landmark_count,
    lgkk_multiscale_scope = lgkk_multiscale_scope,
    lgkk_active_limit = lgkk_active_limit
  )
  lgkk_polish_rounds <- lgkk$lgkk_polish_rounds
  lgkk_multiscale_rounds <- lgkk$lgkk_multiscale_rounds
  lgkk_rounds_coarse <- lgkk$lgkk_rounds_coarse
  lgkk_rounds_pre_final <- lgkk$lgkk_rounds_pre_final
  lgkk_rounds_final <- lgkk$lgkk_rounds_final
  lgkk_local_nbrs <- lgkk$lgkk_local_nbrs
  lgkk_landmark_count <- lgkk$lgkk_landmark_count
  lgkk_multiscale_scope <- lgkk$lgkk_multiscale_scope
  lgkk_active_limit <- lgkk$lgkk_active_limit

  comp <- grip.connected.components(adj_list = adj_list, n = n)
  n.comp <- length(unique(comp))
  if (n.comp != 1L) {
    stop(sprintf(
      "grip.layout.trace.weighted() currently supports only connected graphs; input graph has %d connected components.",
      n.comp
    ))
  }

  out <- grip_layout_globalrep_weighted_trace_adj_cpp(
    adj_list = adj_list,
    weight_list = weight_list,
    n = n,
    dim = dim,
    placement = placement,
    rounds = as.integer(rounds),
    final_rounds = as.integer(final_rounds),
    num_init = as.integer(num_init),
    num_nbrs = num_nbrs,
    r = r,
    s = s,
    repulsion_factor = repulsion_factor,
    coarse_repulsion_factor = coarse_repulsion_factor,
    coarse_repulsion_sample = coarse_repulsion_sample,
    coarse_repulsion_exact_below = coarse_repulsion_exact_below,
    final_anchor_factor = final_anchor_factor,
    final_move_scale_after_first = final_move_scale_after_first,
    insertion_anchor_count = insertion_anchor_count,
    insertion_anchor_scope = insertion_anchor_scope,
    insertion_anchor_strategy = insertion_anchor_strategy,
    level0_insertion_mode = level0_insertion_mode,
    level0_anchor_count = level0_anchor_count,
    level0_local_kk_steps = level0_local_kk_steps,
    lgkk_multiscale_rounds = lgkk_multiscale_rounds,
    lgkk_rounds_coarse = lgkk_rounds_coarse,
    lgkk_rounds_pre_final = lgkk_rounds_pre_final,
    lgkk_rounds_final = lgkk_rounds_final,
    lgkk_local_nbrs = lgkk_local_nbrs,
    lgkk_landmark_count = lgkk_landmark_count,
    lgkk_multiscale_scope = lgkk_multiscale_scope,
    lgkk_active_limit = lgkk_active_limit,
    final_mode = final_mode,
    tinit_factor = as.integer(tinit_factor),
    seed = seed,
    trace = trace,
    trace_every = trace.every,
    metric_neighbor_cap = metric_neighbor_cap
  )

  if (lgkk_polish_rounds > 0L) {
    polished <- grip.apply.lgkk.polish(
      coords = out$final,
      adj_list = adj_list,
      weight_list = weight_list,
      rounds = lgkk_polish_rounds,
      lgkk_local_nbrs = lgkk_local_nbrs,
      lgkk_landmark_count = lgkk_landmark_count,
      return_trace = TRUE
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
    adj.list = adj_list,
    weight.list = weight_list,
    diagnostics = diagnostics,
    target.coords = target_coords,
    sample.size.nonedge = diagnostic_sample_size_nonedge,
    sample.size.stress = diagnostic_sample_size_stress,
    nonedge.seed = diagnostic_nonedge_seed,
    stress.seed = diagnostic_stress_seed
  )
  class(out) <- c("grip_layout_trace", class(out))
  out
}

#' Alias of \code{grip.layout.globalrep.weighted()}
#'
#' @inheritParams grip.layout.globalrep.weighted
#' @return A numeric matrix with \code{n} rows and \code{dim} columns.
#' @export
grip.layout.weighted <- function(edges = NULL,
                                 n = NULL,
                                 adj_list = NULL,
                                 weight_list = NULL,
                                 edge_weights = NULL,
                                 dim = 3,
                                 placement = c("barycenter", "circle"),
                                 preset = NULL,
                                 rounds = 160,
                                 final_rounds = 384,
                                 num_init = 24,
                                 num_nbrs = 20,
                                 r = 0.03,
                                 s = 7.5,
                                 repulsion_factor = 2.5,
                                 coarse_repulsion_factor = 1.5,
                                 coarse_repulsion_sample = 16,
                                 coarse_repulsion_exact_below = 64,
                                 final_anchor_factor = 0,
                                 final_move_scale_after_first = 1,
                                 final_mode = c("fr", "kk_repulse"),
                                 insertion_anchor_count = 3,
                                 insertion_anchor_scope = c("any_higher", "prev_misf"),
                                 insertion_anchor_strategy = c("first", "distance_band", "balanced_band", "spread_prev"),
                                 level0_insertion_mode = c("inherit", "barycenter", "least_squares"),
                                 level0_anchor_count = insertion_anchor_count,
                                 level0_local_kk_steps = 3,
                                 lgkk_polish_rounds = 0L,
                                 lgkk_multiscale_rounds = 0L,
                                 lgkk_rounds_coarse = NULL,
                                 lgkk_rounds_pre_final = NULL,
                                 lgkk_rounds_final = NULL,
                                 lgkk_local_nbrs = 20L,
                                 lgkk_landmark_count = 8L,
                                 lgkk_multiscale_scope = c("all", "coarse"),
                                 lgkk_active_limit = 4096L,
                                 metric_neighbor_cap = NULL,
                                 length_normalization = c("median", "mean", "none"),
                                 tinit_factor = 6,
                                 seed = 6,
                                 disconnected = c("components", "error")) {
  grip.forward_call(grip.layout.globalrep.weighted, match.call(expand.dots = FALSE), env = parent.frame())
}
