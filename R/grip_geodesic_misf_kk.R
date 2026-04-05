grip.new.misf.geodesic.kk.prepared <- function(x) {
  class(x) <- unique(c("grip_misf_gkk_prepared", class(x), "list"))
  x
}

grip.new.misf.geodesic.kk.fit <- function(x = list()) {
  class(x) <- unique(c("grip_misf_gkk_fit", class(x), "list"))
  x
}

grip.validate.misf.geodesic.kk.prepared <- function(prepared, coords = NULL) {
  if (!inherits(prepared, "grip_misf_gkk_prepared")) {
    stop("prepared must be an object from grip.prepare.misf.geodesic.kk()")
  }
  if (!is.null(coords) && nrow(coords) != prepared$n) {
    stop("nrow(coords) must match the graph size stored in prepared")
  }
  prepared
}

grip.validate.misf.geodesic.kk.fit <- function(fit) {
  if (!inherits(fit, "grip_misf_gkk_fit")) {
    stop("fit must be an object from grip.optimize.misf.geodesic.kk()")
  }
  fit
}

grip.geodesic.misf.kk.resolve.pair.mode <- function(pair_mode = c("auto", "full", "landmark"),
                                                     active_n,
                                                     full_limit = 512L) {
  pair_mode <- match.arg(pair_mode)
  active_n <- grip.validate.count(active_n, "active_n")
  full_limit <- grip.validate.misf.count(full_limit, "full_limit", lower = 1L)
  effective <- if (identical(pair_mode, "auto")) {
    if (active_n <= full_limit) "full" else "landmark"
  } else {
    pair_mode
  }
  list(
    requested = pair_mode,
    effective = effective,
    full_limit = full_limit,
    active_n = active_n
  )
}

grip.resolve.misf.geodesic.kk.prepared <- function(prepared = NULL,
                                                   edges = NULL,
                                                   n = NULL,
                                                   adj_list = NULL,
                                                   weight_list = NULL,
                                                   edge_weights = NULL,
                                                   tie_mode = NULL,
                                                   num_init = 24L,
                                                   num_nbrs = 20L,
                                                   dim = NULL,
                                                   top_level_pair_mode = c("auto", "full", "landmark"),
                                                   top_level_full_limit = 512L,
                                                   top_level_local_nbrs = 20L,
                                                   top_level_landmark_count = 8L,
                                                   top_level_restarts = 8L,
                                                   top_level_max_iter = 16L,
                                                   top_level_init = c("cmdscale", "random"),
                                                   seed = 6L) {
  top_level_pair_mode <- match.arg(top_level_pair_mode)
  top_level_init <- match.arg(top_level_init)
  if (is.null(prepared)) {
    resolved.dim <- if (is.null(dim)) 2L else grip.validate.count(dim, "dim")
    resolved.tie.mode <- if (is.null(tie_mode)) "average" else {
      match.arg(tie_mode, c("single", "average"))
    }
    return(grip.prepare.misf.geodesic.kk(
      edges = edges,
      n = n,
      adj_list = adj_list,
      weight_list = weight_list,
      edge_weights = edge_weights,
      tie_mode = resolved.tie.mode,
      num_init = num_init,
      num_nbrs = num_nbrs,
      dim = resolved.dim,
      top_level_mode = "skip",
      top_level_pair_mode = top_level_pair_mode,
      top_level_full_limit = top_level_full_limit,
      top_level_local_nbrs = top_level_local_nbrs,
      top_level_landmark_count = top_level_landmark_count,
      top_level_restarts = top_level_restarts,
      top_level_max_iter = top_level_max_iter,
      top_level_init = top_level_init,
      seed = seed
    ))
  }
  if (inherits(prepared, "grip_misf_gkk_prepared")) {
    return(grip.validate.misf.geodesic.kk.prepared(prepared))
  }
  if (!inherits(prepared, "grip_gkk_prepared")) {
    stop(
      "prepared must be NULL, an object from grip.prepare.geodesic.kk(), ",
      "or an object from grip.prepare.misf.geodesic.kk()"
    )
  }
  resolved.dim <- if (is.null(dim)) 2L else grip.validate.count(dim, "dim")
  resolved.tie.mode <- if (is.null(tie_mode)) {
    if (!is.null(prepared$tie_mode)) prepared$tie_mode else "average"
  } else {
    match.arg(tie_mode, c("single", "average"))
  }
  grip.prepare.misf.geodesic.kk(
    n = prepared$n,
    adj_list = prepared$adj_list,
    weight_list = prepared$weight_list,
    tie_mode = resolved.tie.mode,
    num_init = num_init,
    num_nbrs = num_nbrs,
    dim = resolved.dim,
    top_level_mode = "skip",
    top_level_pair_mode = top_level_pair_mode,
    top_level_full_limit = top_level_full_limit,
    top_level_local_nbrs = top_level_local_nbrs,
    top_level_landmark_count = top_level_landmark_count,
    top_level_restarts = top_level_restarts,
    top_level_max_iter = top_level_max_iter,
    top_level_init = top_level_init,
    seed = seed
  )
}

grip.geodesic.misf.kk.extract.score.metric <- function(score, names) {
  for (name in names) {
    if (!is.null(score[[name]]) && length(score[[name]]) >= 1L) {
      return(as.double(score[[name]][[1L]]))
    }
  }
  NA_real_
}

#' Prepare a MISF-based multiscale geodesic-KK object
#'
#' `grip.prepare.misf.geodesic.kk()` scaffolds a MISF-based geodesic-KK
#' prepared object by layering the maximal independent set filtration (MISF) on
#' top of the existing full geodesic-KK cache. The current scaffold builds the
#' graph metadata, the top-level coarse graph, and the exact and sparse
#' top-level prepared objects, but it does not yet run the multiscale solve.
#'
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#' @param n Number of vertices. If omitted with `adj_list`, defaults to
#'   `length(adj_list)`. If omitted with `edges`, defaults to `max(edges)`.
#' @param adj_list Adjacency list (1-based) for an undirected graph.
#' @param weight_list Optional parallel list of positive edge weights.
#' @param edge_weights Optional positive edge-weight vector parallel to `edges`.
#' @param tie_mode Shortest-path aggregation mode inherited from
#'   [grip.prepare.geodesic.kk()].
#' @param num_init Target top-level active-set size passed to `grip.build.misf`.
#' @param num_nbrs Per-level local-neighborhood schedule metadata passed to
#'   `grip.build.misf`.
#' @param dim Target embedding dimension (`2` or `3`) for the future multiscale
#'   solve.
#' @param top_level_mode Either `"solve"` or `"skip"`. The current scaffold
#'   records the requested mode but does not yet execute the top-level solve.
#' @param top_level_pair_mode Pair policy for the top MISF level: `"auto"`,
#'   `"full"`, or `"landmark"`.
#' @param top_level_full_limit Active-set size threshold used when
#'   `top_level_pair_mode = "auto"`.
#' @param top_level_local_nbrs Number of nearest graph-metric neighbors retained
#'   per vertex in the sparse top-level LGKK cache.
#' @param top_level_landmark_count Number of farthest-point landmarks retained
#'   per vertex in the sparse top-level LGKK cache.
#' @param top_level_restarts Planned number of top-level restarts for the future
#'   optimizer.
#' @param top_level_max_iter Planned top-level iteration budget for the future
#'   optimizer.
#' @param top_level_init Planned top-level initializer (`"cmdscale"` or
#'   `"random"`).
#' @param seed Optional integer seed reused for MISF extraction.
#'
#' @return An object of class `"grip_misf_gkk_prepared"` layered on top of the
#'   full geodesic-KK prepared structure, with added MISF metadata, top-level
#'   coarse-graph caches, and stored default controls for the future MISF-GKK
#'   optimizer.
#'
#' @examples
#' edges <- edges.mesh(4, 4)
#' prepared <- grip.prepare.misf.geodesic.kk(
#'   edges = edges,
#'   n = 16,
#'   tie_mode = "average",
#'   num_init = 4,
#'   top_level_mode = "skip",
#'   seed = 1
#' )
#' prepared$top_level_vertices
#' @export
grip.prepare.misf.geodesic.kk <- function(edges = NULL,
                                          n = NULL,
                                          adj_list = NULL,
                                          weight_list = NULL,
                                          edge_weights = NULL,
                                          tie_mode = c("single", "average"),
                                          num_init = 24L,
                                          num_nbrs = 20L,
                                          dim = 2L,
                                          top_level_mode = c("solve", "skip"),
                                          top_level_pair_mode = c("auto", "full", "landmark"),
                                          top_level_full_limit = 512L,
                                          top_level_local_nbrs = 20L,
                                          top_level_landmark_count = 8L,
                                          top_level_restarts = 8L,
                                          top_level_max_iter = 16L,
                                          top_level_init = c("cmdscale", "random"),
                                          seed = 6L) {
  tie_mode <- match.arg(tie_mode)
  top_level_mode <- match.arg(top_level_mode)
  top_level_pair_mode <- match.arg(top_level_pair_mode)
  top_level_init <- match.arg(top_level_init)
  dim <- grip.validate.count(dim, "dim")
  if (!(dim %in% c(2L, 3L))) {
    stop("dim must be 2 or 3")
  }
  top_level_full_limit <- grip.validate.misf.count(
    top_level_full_limit,
    "top_level_full_limit",
    lower = 1L
  )
  top_level_local_nbrs <- grip.validate.count(top_level_local_nbrs, "top_level_local_nbrs")
  top_level_landmark_count <- grip.validate.count(top_level_landmark_count, "top_level_landmark_count")
  top_level_restarts <- grip.validate.misf.count(top_level_restarts, "top_level_restarts", lower = 1L)
  grip.validate.scalar(top_level_max_iter, "top_level_max_iter", lower = 0)
  top_level_max_iter <- as.integer(round(top_level_max_iter))
  if (!is.null(seed)) {
    seed <- grip.validate.count(seed, "seed")
  }

  prepared <- grip.prepare.geodesic.kk(
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights,
    tie_mode = tie_mode
  )
  misf <- grip.build.misf(
    edges = prepared$edges,
    n = prepared$n,
    adj_list = prepared$adj_list,
    weight_list = prepared$weight_list,
    num_init = num_init,
    num_nbrs = num_nbrs,
    seed = seed
  )
  top_level_index <- length(misf$levels)
  top_level_id <- as.integer(top_level_index - 1L)
  top_level_vertices <- as.integer(misf$levels[[top_level_index]])
  top_level_graph <- grip.geodesic.misf.induced_level_graph(
    prepared = prepared,
    vertex_ids = top_level_vertices,
    level = top_level_id
  )
  top_level_resolution <- grip.geodesic.misf.kk.resolve.pair.mode(
    pair_mode = top_level_pair_mode,
    active_n = top_level_graph$n,
    full_limit = top_level_full_limit
  )
  top_level_prepared_full <- grip.prepare.geodesic.kk(
    edges = top_level_graph$edges,
    n = top_level_graph$n,
    edge_weights = top_level_graph$edge_weights,
    tie_mode = tie_mode
  )
  top_level_prepared_sparse <- grip.prepare.landmark.geodesic.kk(
    edges = top_level_graph$edges,
    n = top_level_graph$n,
    edge_weights = top_level_graph$edge_weights,
    local_nbrs = top_level_local_nbrs,
    landmark_count = top_level_landmark_count
  )

  prepared$misf <- misf
  prepared$level_vertices <- misf$levels
  prepared$active_levels <- misf$levels
  prepared$insertion_order <- misf$mish_order
  prepared$top_level_index <- top_level_index
  prepared$top_level_level <- top_level_id
  prepared$top_level_vertices <- top_level_vertices
  prepared$top_level_graph <- top_level_graph
  prepared$top_level_pair_mode <- top_level_resolution$requested
  prepared$top_level_effective_pair_mode <- top_level_resolution$effective
  prepared$top_level_full_limit <- top_level_resolution$full_limit
  prepared$top_level_local_nbrs <- top_level_local_nbrs
  prepared$top_level_landmark_count <- top_level_landmark_count
  prepared$top_level_prepared_full <- top_level_prepared_full
  prepared$top_level_prepared_sparse <- top_level_prepared_sparse
  prepared$top_level_prepared <- if (identical(top_level_resolution$effective, "full")) {
    top_level_prepared_full
  } else {
    top_level_prepared_sparse
  }
  prepared$top_level_dim <- as.integer(dim)
  prepared$top_level_mode <- top_level_mode
  prepared$top_level_restarts <- as.integer(top_level_restarts)
  prepared$top_level_max_iter <- as.integer(top_level_max_iter)
  prepared$top_level_init <- top_level_init
  prepared$multiscale_mode <- "misf"
  prepared$insertion_anchor_policy <- "prev_level_spread"
  prepared$insertion_anchor_count <- grip.geodesic.misf.default.anchor.count(dim)
  prepared$insertion_anchor_weight_mode <- "inverse_graph_distance_sq"
  prepared$insertion_max_iter <- 32L
  prepared$insertion_mode <- "weighted_kk"
  prepared$insertion_layout_k <- 6L
  prepared$insertion_weighted_preset <- NULL
  prepared$insertion_grip_args <- list()
  prepared$insertion_weighted_args <- list()
  prepared$insertion_fr_niter <- 800L
  prepared$refinement_pair_mode <- "auto"
  prepared$refinement_full_limit <- 256L
  prepared$refinement_local_nbrs <- 8L
  prepared$refinement_landmark_count <- 4L
  prepared$refinement_max_iter <- 8L
  prepared$final_pair_mode <- "auto"
  prepared$final_full_limit <- 1024L
  prepared$final_local_nbrs <- 20L
  prepared$final_landmark_count <- 8L
  prepared$final_max_iter <- 8L
  prepared$stiffness <- 1.0
  prepared$distance_floor <- 1e-8
  prepared$edge_length_epsilon <- 1e-8
  prepared$initial_step <- 1.0
  prepared$step_shrink <- 0.5
  prepared$armijo_factor <- 1e-4
  prepared$grad_tol <- 1e-8
  prepared$min_step <- 1e-8
  prepared$recenter <- TRUE
  prepared$misf_seed <- seed
  prepared$top_level_fit <- NULL
  prepared <- grip.new.misf.geodesic.kk.prepared(prepared)

  if (identical(top_level_mode, "solve")) {
    warning(
      "top_level_mode = 'solve' is scaffolded but not implemented yet; ",
      "returning a prepared object with top_level_fit = NULL"
    )
  }
  prepared
}

#' Optimize an embedding with the MISF-based geodesic-KK pipeline
#'
#' `grip.optimize.misf.geodesic.kk()` is currently scaffolded but not yet
#' implemented. The public signature is in place so the future multiscale
#' geodesic-KK solver can be added without another API rename.
#'
#' @inheritParams grip.prepare.misf.geodesic.kk
#' @param prepared Optional prepared object. This can be either a full
#'   geodesic-KK prepared object or a MISF-GKK prepared object.
#' @param top_level_pair_mode Optional override for the top-level pair policy.
#' @param top_level_full_limit Optional override for the top-level exact/sparse
#'   threshold.
#' @param top_level_local_nbrs Optional override for the sparse top-level local
#'   neighborhood size.
#' @param top_level_landmark_count Optional override for the sparse top-level
#'   landmark count.
#' @param top_level_init Optional override for the top-level initializer.
#' @param insertion_anchor_policy Optional insertion anchor policy.
#' @param insertion_anchor_count Optional insertion anchor count.
#' @param insertion_anchor_weight_mode Optional insertion anchor weighting mode.
#' @param insertion_max_iter Optional per-vertex insertion iteration budget.
#' @param insertion_mode Optional insertion warm-start mode.
#' @param insertion_layout_k Optional active-level layout graph size used by
#'   layout-based insertion modes.
#' @param insertion_weighted_preset Optional weighted preset forwarded by future
#'   weighted insertion modes.
#' @param insertion_grip_args Optional named list of extra arguments for future
#'   combinatorial GRIP insertion modes.
#' @param insertion_weighted_args Optional named list of extra arguments for
#'   future weighted GRIP insertion modes.
#' @param insertion_fr_niter Optional FR iteration budget for future FR-based
#'   insertion modes.
#' @param refinement_pair_mode Optional active-level pair policy.
#' @param refinement_full_limit Optional active-level exact/sparse threshold.
#' @param refinement_local_nbrs Optional active-level sparse local neighborhood
#'   size.
#' @param refinement_landmark_count Optional active-level sparse landmark count.
#' @param refinement_max_iter Optional active-level refinement iteration budget.
#' @param final_pair_mode Optional final full-graph pair policy.
#' @param final_full_limit Optional final full-graph exact/sparse threshold.
#' @param final_local_nbrs Optional final full-graph sparse local neighborhood
#'   size.
#' @param final_landmark_count Optional final full-graph sparse landmark count.
#' @param final_max_iter Optional final polish iteration budget.
#' @param stiffness Global stiffness constant \(K\).
#' @param distance_floor Small positive floor used in
#'   \code{k_ij = K / max(g_ij, distance_floor)^2}.
#' @param edge_length_epsilon Small positive stabilizer added inside each
#'   embedded edge length.
#' @param initial_step Initial line-search step size.
#' @param step_shrink Multiplicative shrink factor in `(0, 1)` for backtracking.
#' @param armijo_factor Non-negative Armijo decrease constant.
#' @param grad_tol Non-negative stopping tolerance on the gradient norm.
#' @param min_step Positive minimum accepted line-search step before giving up.
#' @param recenter If `TRUE`, recenter accepted proposals to zero mean after
#'   each accepted step.
#' @param return_trace If `TRUE`, retain detailed stage traces when the solver is
#'   implemented.
#' @param return_frames If `TRUE`, retain intermediate coordinate frames when
#'   the solver is implemented.
#'
#' @return This scaffold currently stops with a not-yet-implemented error after
#'   validating or building the corresponding prepared object.
#'
#' @export
grip.optimize.misf.geodesic.kk <- function(prepared = NULL,
                                           edges = NULL,
                                           n = NULL,
                                           adj_list = NULL,
                                           weight_list = NULL,
                                           edge_weights = NULL,
                                           tie_mode = NULL,
                                           num_init = 24L,
                                           num_nbrs = 20L,
                                           dim = NULL,
                                           top_level_pair_mode = NULL,
                                           top_level_full_limit = NULL,
                                           top_level_local_nbrs = NULL,
                                           top_level_landmark_count = NULL,
                                           top_level_restarts = NULL,
                                           top_level_max_iter = NULL,
                                           top_level_init = NULL,
                                           insertion_anchor_policy = NULL,
                                           insertion_anchor_count = NULL,
                                           insertion_anchor_weight_mode = NULL,
                                           insertion_max_iter = NULL,
                                           insertion_mode = NULL,
                                           insertion_layout_k = NULL,
                                           insertion_weighted_preset = NULL,
                                           insertion_grip_args = NULL,
                                           insertion_weighted_args = NULL,
                                           insertion_fr_niter = NULL,
                                           refinement_pair_mode = NULL,
                                           refinement_full_limit = NULL,
                                           refinement_local_nbrs = NULL,
                                           refinement_landmark_count = NULL,
                                           refinement_max_iter = NULL,
                                           final_pair_mode = NULL,
                                           final_full_limit = NULL,
                                           final_local_nbrs = NULL,
                                           final_landmark_count = NULL,
                                           final_max_iter = NULL,
                                           stiffness = 1.0,
                                           distance_floor = 1e-8,
                                           edge_length_epsilon = 1e-8,
                                           initial_step = 1.0,
                                           step_shrink = 0.5,
                                           armijo_factor = 1e-4,
                                           grad_tol = 1e-8,
                                           min_step = 1e-8,
                                           recenter = TRUE,
                                           return_trace = FALSE,
                                           return_frames = FALSE,
                                           seed = 6L) {
  resolved.pair.mode <- if (is.null(top_level_pair_mode)) "auto" else {
    match.arg(top_level_pair_mode, c("auto", "full", "landmark"))
  }
  resolved.top.level.init <- if (is.null(top_level_init)) "cmdscale" else {
    match.arg(top_level_init, c("cmdscale", "random"))
  }
  prepared <- grip.resolve.misf.geodesic.kk.prepared(
    prepared = prepared,
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights,
    tie_mode = tie_mode,
    num_init = num_init,
    num_nbrs = num_nbrs,
    dim = dim,
    top_level_pair_mode = resolved.pair.mode,
    top_level_full_limit = if (is.null(top_level_full_limit)) 512L else top_level_full_limit,
    top_level_local_nbrs = if (is.null(top_level_local_nbrs)) 20L else top_level_local_nbrs,
    top_level_landmark_count = if (is.null(top_level_landmark_count)) 8L else top_level_landmark_count,
    top_level_restarts = if (is.null(top_level_restarts)) 8L else top_level_restarts,
    top_level_max_iter = if (is.null(top_level_max_iter)) 16L else top_level_max_iter,
    top_level_init = resolved.top.level.init,
    seed = seed
  )
  stop(
    "grip.optimize.misf.geodesic.kk() is scaffolded but not implemented yet. ",
    "Use grip.prepare.misf.geodesic.kk() to build a prepared object and ",
    "grip.score.misf.geodesic.kk() to rescore external coordinates for now."
  )
}

#' Score a MISF-based geodesic-KK scaffold
#'
#' `grip.score.misf.geodesic.kk()` summarizes external coordinates against a
#' MISF-GKK prepared object using the current exact geodesic-KK scorer. The
#' multiscale optimizer itself is not implemented yet, so this scaffold focuses
#' on the final full-graph score and the prepared MISF metadata.
#'
#' @param fit Optional fit from [grip.optimize.misf.geodesic.kk()].
#' @param coords Optional coordinate matrix used when `fit` is omitted.
#' @param prepared Optional MISF-GKK prepared object used when `fit` is omitted.
#' @param stiffness Global stiffness constant \(K\).
#' @param distance_floor Small positive floor used in
#'   \code{k_ij = K / max(g_ij, distance_floor)^2}.
#' @param edge_length_epsilon Small positive stabilizer added inside each
#'   embedded edge length.
#' @param score_pair_mode Requested scoring pair policy. The current scaffold
#'   always falls back to exact full geodesic-KK scoring.
#' @param score_full_limit Active-set size threshold used when
#'   `score_pair_mode = "auto"`.
#' @param score_local_nbrs Planned sparse local-neighborhood size for future
#'   LGKK-based scoring.
#' @param score_landmark_count Planned sparse landmark count for future
#'   LGKK-based scoring.
#' @param return_trace If `TRUE`, attach any trace tables stored inside `fit`.
#'
#' @return A one-row data frame with final full-graph geodesic-KK metrics and
#'   MISF-stage summary fields. If `return_trace = TRUE`, trace list-columns are
#'   attached when available.
#'
#' @export
grip.score.misf.geodesic.kk <- function(fit = NULL,
                                        coords = NULL,
                                        prepared = NULL,
                                        stiffness = 1.0,
                                        distance_floor = 1e-8,
                                        edge_length_epsilon = 1e-8,
                                        score_pair_mode = c("full", "landmark", "auto"),
                                        score_full_limit = 2048L,
                                        score_local_nbrs = 20L,
                                        score_landmark_count = 8L,
                                        return_trace = FALSE) {
  score_pair_mode <- match.arg(score_pair_mode)
  score_full_limit <- grip.validate.misf.count(score_full_limit, "score_full_limit", lower = 1L)
  score_local_nbrs <- grip.validate.count(score_local_nbrs, "score_local_nbrs")
  score_landmark_count <- grip.validate.count(score_landmark_count, "score_landmark_count")
  grip.validate.scalar(stiffness, "stiffness", lower = 0, open.lower = TRUE)
  grip.validate.scalar(distance_floor, "distance_floor", lower = 0, open.lower = TRUE)
  grip.validate.scalar(edge_length_epsilon, "edge_length_epsilon", lower = 0)
  if (!is.logical(return_trace) || length(return_trace) != 1L || is.na(return_trace)) {
    stop("return_trace must be TRUE or FALSE")
  }

  top.level.fit <- NULL
  stage.trace <- data.frame()
  timing <- list(
    top_level = NA_real_,
    insertion = NA_real_,
    refinement = NA_real_,
    final_polish = NA_real_,
    total = NA_real_
  )

  if (!is.null(fit)) {
    fit <- grip.validate.misf.geodesic.kk.fit(fit)
    if (is.null(prepared)) {
      prepared <- fit$prepared
    }
    if (is.null(coords)) {
      coords <- fit$coords
    }
    if (!is.null(fit$top_level_fit)) {
      top.level.fit <- fit$top_level_fit
    }
    if (!is.null(fit$stage_trace)) {
      stage.trace <- fit$stage_trace
    }
    if (!is.null(fit$timing)) {
      timing <- fit$timing
    }
  }

  coords <- grip.validate.coords(coords)
  prepared <- grip.validate.misf.geodesic.kk.prepared(prepared, coords = coords)
  score.resolution <- grip.geodesic.misf.kk.resolve.pair.mode(
    pair_mode = score_pair_mode,
    active_n = prepared$n,
    full_limit = score_full_limit
  )
  if (!identical(score.resolution$effective, "full")) {
    warning(
      "score_pair_mode = '", score.resolution$requested,
      "' is scaffolded but not implemented yet; using exact geodesic-KK scoring"
    )
  }

  final.score <- grip.score.geodesic.kk(
    coords = coords,
    prepared = prepared,
    stiffness = stiffness,
    distance_floor = distance_floor,
    edge_length_epsilon = edge_length_epsilon
  )
  names(final.score) <- paste0("final.", names(final.score))

  top.level.energy <- NA_real_
  top.level.weighted.rmse <- NA_real_
  top.level.weighted.rel.rmse <- NA_real_
  if (!is.null(top.level.fit) && !is.null(top.level.fit$score)) {
    top.level.energy <- grip.geodesic.misf.kk.extract.score.metric(
      top.level.fit$score,
      c("gkk.energy", "lgkk.energy")
    )
    top.level.weighted.rmse <- grip.geodesic.misf.kk.extract.score.metric(
      top.level.fit$score,
      c("gkk.weighted.rmse", "lgkk.weighted.rmse")
    )
    top.level.weighted.rel.rmse <- grip.geodesic.misf.kk.extract.score.metric(
      top.level.fit$score,
      c("gkk.weighted.rel.rmse", "lgkk.weighted.rel.rmse")
    )
  }

  out <- data.frame(
    multiscale.mode = "misf",
    n = prepared$n,
    dim = ncol(coords),
    top.level = prepared$top_level_level,
    top.level.n = length(prepared$top_level_vertices),
    top.level.pair.mode = prepared$top_level_effective_pair_mode,
    top.level.energy = top.level.energy,
    top.level.weighted.rmse = top.level.weighted.rmse,
    top.level.weighted.rel.rmse = top.level.weighted.rel.rmse,
    insertion.level.count = as.integer(max(0L, prepared$top_level_level)),
    inserted.vertex.count = as.integer(max(0L, prepared$n - length(prepared$top_level_vertices))),
    refinement.level.count = if (nrow(stage.trace)) sum(stage.trace$stage == "refinement") else NA_integer_,
    final.polish.trace.rows = NA_integer_,
    elapsed.top.level = as.double(timing$top_level),
    elapsed.insertion = as.double(timing$insertion),
    elapsed.refinement = as.double(timing$refinement),
    elapsed.final.polish = as.double(timing$final_polish),
    elapsed.total = as.double(timing$total),
    stringsAsFactors = FALSE
  )
  out <- cbind(out, final.score)

  if (isTRUE(return_trace)) {
    out$stage.trace <- list(stage.trace)
    out$top.restart.summary <- list(data.frame())
    out$insertion.level.trace <- list(data.frame())
    out$insertion.vertex.trace <- list(data.frame())
    out$refinement.level.trace <- list(data.frame())
    out$final.polish.trace <- list(data.frame())
  }
  out
}
