grip.validate.misf.geodesic.prepared <- function(prepared, coords = NULL) {
  if (!inherits(prepared, "grip_misf_gmds_prepared")) {
    stop("prepared must be an object from grip.prepare.misf.geodesic.mds()")
  }
  grip.validate.geodesic.mds.prepared(prepared, coords = coords)
}

grip.geodesic.misf.level.to.index <- function(misf, level = NULL) {
  if (is.null(level)) {
    return(length(misf$levels))
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level)) {
    stop("level must be a single finite numeric value")
  }
  level <- as.integer(round(level))
  if (level < 0L) {
    stop("level must be >= 0")
  }
  if (level > misf$misf_height) {
    stop(sprintf("level must be in [0, %d]", misf$misf_height))
  }
  as.integer(level + 1L)
}

grip.geodesic.misf.complete.edge.matrix <- function(n) {
  n <- as.integer(n)
  if (is.na(n) || n < 2L) {
    return(matrix(integer(), ncol = 2L))
  }
  out <- t(utils::combn(n, 2L))
  matrix(as.integer(out), ncol = 2L)
}

grip.geodesic.misf.partial.coords <- function(coords, vertex_ids, n) {
  out <- matrix(NA_real_, nrow = n, ncol = ncol(coords))
  out[vertex_ids, ] <- coords
  out
}

grip.geodesic.misf.default.anchor.count <- function(dim) {
  as.integer(max(grip.validate.count(dim, "dim") + 1L, 6L))
}

grip.geodesic.misf.partial.coords.from.top.level <- function(prepared) {
  prepared <- grip.validate.misf.geodesic.prepared(prepared)
  if (is.null(prepared$top_level_fit) || is.null(prepared$top_level_fit$coords_full)) {
    stop("prepared must contain a top-level fit before MISF insertion can start")
  }
  coords <- prepared$top_level_fit$coords_full
  storage.mode(coords) <- "double"
  coords
}

grip.geodesic.misf.validate.partial.coords <- function(coords,
                                                       prepared,
                                                       dim = NULL) {
  prepared <- grip.validate.misf.geodesic.prepared(prepared)
  if (is.null(coords)) {
    coords <- grip.geodesic.misf.partial.coords.from.top.level(prepared)
  }
  if (!is.matrix(coords) || !is.numeric(coords)) {
    stop("coords must be NULL or a numeric matrix")
  }
  if (nrow(coords) != prepared$n) {
    stop("nrow(coords) must match prepared$n")
  }
  if (!is.null(dim) && ncol(coords) != as.integer(dim)) {
    stop("ncol(coords) must match the requested dim")
  }
  finite.mask <- is.finite(coords)
  partial.rows <- rowSums(finite.mask) > 0L & rowSums(finite.mask) < ncol(coords)
  if (any(partial.rows)) {
    stop("each coords row must be either fully finite or fully NA")
  }
  storage.mode(coords) <- "double"
  coords
}

grip.geodesic.misf.level.insert.vertices <- function(prepared, level) {
  prepared <- grip.validate.misf.geodesic.prepared(prepared)
  level.index <- grip.geodesic.misf.level.to.index(prepared$misf, level)
  level.id <- as.integer(level.index - 1L)
  if (level.id >= prepared$misf$misf_height) {
    return(integer(0L))
  }
  current <- prepared$misf$levels[[level.index]]
  next.level <- prepared$misf$levels[[level.index + 1L]]
  vertices <- setdiff(current, next.level)
  as.integer(prepared$insertion_order[prepared$insertion_order %in% vertices])
}

grip.geodesic.misf.previous.level.vertices <- function(prepared, level) {
  prepared <- grip.validate.misf.geodesic.prepared(prepared)
  level.index <- grip.geodesic.misf.level.to.index(prepared$misf, level)
  level.id <- as.integer(level.index - 1L)
  if (level.id >= prepared$misf$misf_height) {
    stop("the top MISF level has no coarser anchor level")
  }
  as.integer(prepared$misf$levels[[level.index + 1L]])
}

grip.geodesic.misf.anchor.weights <- function(anchor_distances,
                                              mode = c("inverse_graph_distance_sq", "uniform")) {
  mode <- match.arg(mode)
  anchor_distances <- as.double(anchor_distances)
  if (!length(anchor_distances)) {
    return(numeric(0L))
  }
  if (identical(mode, "uniform")) {
    return(rep.int(1, length(anchor_distances)))
  }
  scale <- pmax(anchor_distances, sqrt(.Machine$double.eps))
  1 / (scale * scale)
}

grip.geodesic.misf.distance.band.order <- function(count, n) {
  if (n <= 0L || count <= 0L) {
    return(integer(0L))
  }
  base <- unique(pmax(1L, pmin(n, as.integer(round(seq(1, n, length.out = min(count, n)))))))
  if (length(base) < min(count, n)) {
    base <- unique(c(base, seq_len(n)))
  }
  as.integer(base[seq_len(min(count, length(base)))])
}

grip.geodesic.misf.spread.order <- function(candidate_ids,
                                            candidate_coords,
                                            candidate_distances,
                                            count) {
  if (!length(candidate_ids) || count <= 0L) {
    return(integer(0L))
  }
  ord.near <- order(candidate_distances, candidate_ids)
  sorted.ids <- candidate_ids[ord.near]
  sorted.coords <- candidate_coords[ord.near, , drop = FALSE]
  sorted.dist <- candidate_distances[ord.near]
  selected <- sorted.ids[[1L]]
  selected.idx <- 1L

  while (length(selected) < min(count, length(sorted.ids))) {
    remaining <- setdiff(seq_along(sorted.ids), selected.idx)
    if (!length(remaining)) {
      break
    }
    selected.coords <- sorted.coords[selected.idx, , drop = FALSE]
    scores <- vapply(remaining, function(idx) {
      diffs <- sweep(selected.coords, 2L, sorted.coords[idx, ], FUN = "-")
      min(sqrt(rowSums(diffs^2)))
    }, numeric(1L))
    ord <- order(-scores, sorted.dist[remaining], sorted.ids[remaining])
    choice <- remaining[[ord[[1L]]]]
    selected.idx <- c(selected.idx, choice)
    selected <- c(selected, sorted.ids[[choice]])
  }

  as.integer(selected)
}

grip.geodesic.misf.build.restart.row <- function(restart,
                                                  seed,
                                                  init.score,
                                                  fit) {
  data.frame(
    restart = as.integer(restart),
    seed = if (is.null(seed)) NA_integer_ else as.integer(seed),
    initial.energy = as.double(init.score$gmds.energy[[1L]]),
    initial.stress = as.double(init.score$gmds.stress[[1L]]),
    final.energy = as.double(fit$score$gmds.energy[[1L]]),
    final.stress = as.double(fit$score$gmds.stress[[1L]]),
    improved = as.logical(fit$score$gmds.energy[[1L]] <= init.score$gmds.energy[[1L]] + 1e-12),
    trace.rows = as.integer(if (is.null(fit$trace)) 0L else nrow(fit$trace)),
    stringsAsFactors = FALSE
  )
}

#' Build the MISF-induced coarse graph for a GMDS level
#'
#' `grip.geodesic.misf.induced_level_graph()` extracts a level of the maximal
#' independent set filtration and turns it into the weighted complete graph
#' whose edge weights are the original full-graph geodesic distances restricted
#' to that level. This is the coarse graph used by the Phase 2 MISF-GMDS
#' initializer.
#'
#' @param prepared A graph-first GMDS prepared object or a MISF-GMDS prepared
#'   object.
#' @param level Optional MISF level index in the filtration numbering `V_0,
#'   V_1, ...`. If omitted, the coarsest level is used.
#' @param vertex_ids Optional explicit vertex ids from the original graph. When
#'   supplied, `level` is ignored.
#'
#' @return A list describing the coarse graph with local vertex numbering,
#'   original vertex ids, complete weighted edges, and the restricted graph
#'   distance matrix.
grip.geodesic.misf.induced_level_graph <- function(prepared,
                                                   level = NULL,
                                                   vertex_ids = NULL) {
  prepared <- grip.validate.geodesic.mds.prepared(prepared)

  if (is.null(vertex_ids)) {
    if (!inherits(prepared, "grip_misf_gmds_prepared")) {
      stop("vertex_ids must be supplied when prepared is not a MISF-GMDS prepared object")
    }
    level.index <- grip.geodesic.misf.level.to.index(prepared$misf, level)
    vertex_ids <- prepared$misf$levels[[level.index]]
    level.id <- as.integer(level.index - 1L)
  } else {
    vertex_ids <- as.integer(vertex_ids)
    if (length(vertex_ids) == 0L) {
      stop("vertex_ids must contain at least one vertex id")
    }
    if (any(!is.finite(vertex_ids)) || any(vertex_ids < 1L) || any(vertex_ids > prepared$n)) {
      stop("vertex_ids must be valid 1-based vertex ids from the prepared graph")
    }
    if (anyDuplicated(vertex_ids)) {
      stop("vertex_ids must be unique")
    }
    level.id <- if (!is.null(level)) as.integer(level) else NA_integer_
  }

  sub.dist <- prepared$distance_matrix[vertex_ids, vertex_ids, drop = FALSE]
  edge.matrix <- grip.geodesic.misf.complete.edge.matrix(length(vertex_ids))
  edge.weights <- if (nrow(edge.matrix) == 0L) {
    numeric(0L)
  } else {
    as.double(sub.dist[cbind(edge.matrix[, 1L], edge.matrix[, 2L])])
  }
  global.edge.matrix <- if (nrow(edge.matrix) == 0L) {
    matrix(integer(), ncol = 2L)
  } else {
    cbind(vertex_ids[edge.matrix[, 1L]], vertex_ids[edge.matrix[, 2L]])
  }

  list(
    level = level.id,
    n = length(vertex_ids),
    vertex_ids = as.integer(vertex_ids),
    distance_matrix = sub.dist,
    edges = edge.matrix,
    edge_weights = edge.weights,
    global_edge_matrix = matrix(as.integer(global.edge.matrix), ncol = 2L)
  )
}

#' Solve the coarsest MISF level with restartable pure GMDS
#'
#' `grip.geodesic.misf.solve.top.level()` runs a pure-GMDS solve on the coarse
#' top MISF graph using multiple random restarts. It is the Phase 2 replacement
#' for ad hoc random or spectral top-level initialization.
#'
#' @param prepared A MISF-GMDS prepared object from
#'   `grip.prepare.misf.geodesic.mds()`, or directly a graph-first GMDS
#'   prepared object representing a coarse level.
#' @param dim Target embedding dimension (`2` or `3`).
#' @param n_restarts Number of random restarts.
#' @param max_iter Maximum number of pure-GMDS iterations per restart.
#' @param engine Optimization engine passed through to
#'   `grip.optimize.geodesic.mds()`.
#' @param edge_length_epsilon Small non-negative edge-length stabilizer.
#' @param initial_step Initial Armijo line-search step.
#' @param step_shrink Backtracking shrink factor.
#' @param armijo_factor Armijo decrease constant.
#' @param grad_tol Gradient-norm stopping tolerance.
#' @param min_step Minimum accepted line-search step.
#' @param n_threads Number of compiled-engine threads.
#' @param recenter Whether to recenter accepted proposals to zero mean.
#' @param return_trace Whether to retain per-iteration traces/frames for the
#'   best restart.
#' @param seed Optional base seed; restart `r` uses `seed + r - 1`.
#'
#' @return The best restart fit, together with a restart summary and the coarse
#'   vertex ids.
grip.geodesic.misf.solve.top.level <- function(prepared,
                                               dim = 2L,
                                               n_restarts = 8L,
                                               max_iter = 16L,
                                               engine = c("cpp", "r"),
                                               edge_length_epsilon = 1e-8,
                                               initial_step = 1.0,
                                               step_shrink = 0.5,
                                               armijo_factor = 1e-4,
                                               grad_tol = 1e-8,
                                               min_step = 1e-8,
                                               n_threads = 0L,
                                               recenter = TRUE,
                                               return_trace = FALSE,
                                               seed = 6L) {
  engine <- match.arg(engine)
  dim <- grip.validate.count(dim, "dim")
  if (!(dim %in% c(2L, 3L))) {
    stop("dim must be 2 or 3")
  }
  n_restarts <- grip.validate.misf.count(n_restarts, "n_restarts", lower = 1L)
  grip.validate.scalar(max_iter, "max_iter", lower = 0)
  max_iter <- as.integer(round(max_iter))
  grip.validate.scalar(edge_length_epsilon, "edge_length_epsilon", lower = 0)
  grip.validate.scalar(initial_step, "initial_step", lower = 0, open.lower = TRUE)
  grip.validate.scalar(step_shrink, "step_shrink", lower = 0, upper = 1, open.lower = TRUE, open.upper = TRUE)
  grip.validate.scalar(armijo_factor, "armijo_factor", lower = 0)
  grip.validate.scalar(grad_tol, "grad_tol", lower = 0)
  grip.validate.scalar(min_step, "min_step", lower = 0, open.lower = TRUE)
  grip.validate.scalar(n_threads, "n_threads", lower = 0)
  n_threads <- as.integer(round(n_threads))
  if (!is.null(seed)) {
    seed <- grip.validate.count(seed, "seed")
  }

  is.misf.prepared <- inherits(prepared, "grip_misf_gmds_prepared") ||
    (!is.null(prepared$top_level_prepared) && !is.null(prepared$top_level_vertices))
  coarse.prepared <- if (is.misf.prepared) {
    prepared$top_level_prepared
  } else {
    prepared
  }
  coarse.prepared <- grip.validate.geodesic.mds.prepared(coarse.prepared)
  vertex.ids <- if (is.misf.prepared) {
    prepared$top_level_vertices
  } else {
    seq_len(coarse.prepared$n)
  }
  full.n <- if (is.misf.prepared) prepared$n else coarse.prepared$n

  if (coarse.prepared$n <= 1L) {
    coords <- matrix(0, nrow = coarse.prepared$n, ncol = dim)
    score <- grip.score.geodesic.mds(
      coords = coords,
      prepared = coarse.prepared,
      edge_length_epsilon = edge_length_epsilon
    )
    return(list(
      coords = coords,
      trace = data.frame(),
      frames = list(coords),
      prepared = coarse.prepared,
      score = score,
      restart_summary = data.frame(
        restart = 1L,
        seed = if (is.null(seed)) NA_integer_ else seed,
        initial.energy = score$gmds.energy[[1L]],
        initial.stress = score$gmds.stress[[1L]],
        final.energy = score$gmds.energy[[1L]],
        final.stress = score$gmds.stress[[1L]],
        improved = TRUE,
        trace.rows = 0L
      ),
      best_restart = 1L,
      vertex_ids = as.integer(vertex.ids),
      coords_full = grip.geodesic.misf.partial.coords(coords, vertex.ids, full.n)
    ))
  }

  restart.rows <- vector("list", n_restarts)
  best.fit <- NULL
  best.row <- NULL
  best.restart <- 1L
  best.energy <- Inf

  for (restart in seq_len(n_restarts)) {
    restart.seed <- if (is.null(seed)) NULL else as.integer(seed + restart - 1L)
    if (!is.null(restart.seed)) {
      set.seed(restart.seed)
    }
    init.coords <- matrix(stats::rnorm(coarse.prepared$n * dim), ncol = dim)
    storage.mode(init.coords) <- "double"
    if (isTRUE(recenter)) {
      init.coords <- sweep(init.coords, 2L, colMeans(init.coords), FUN = "-")
    }
    init.score <- grip.score.geodesic.mds(
      coords = init.coords,
      prepared = coarse.prepared,
      edge_length_epsilon = edge_length_epsilon
    )
    fit <- grip.optimize.geodesic.mds(
      coords = init.coords,
      prepared = coarse.prepared,
      init = "user",
      anchor_mode = "none",
      engine = engine,
      max_iter = max_iter,
      edge_length_epsilon = edge_length_epsilon,
      initial_step = initial_step,
      step_shrink = step_shrink,
      armijo_factor = armijo_factor,
      grad_tol = grad_tol,
      min_step = min_step,
      n_threads = n_threads,
      recenter = recenter,
      return_trace = return_trace
    )
    restart.rows[[restart]] <- grip.geodesic.misf.build.restart.row(
      restart = restart,
      seed = restart.seed,
      init.score = init.score,
      fit = fit
    )
    final.energy <- fit$score$gmds.energy[[1L]]
    if (!is.finite(best.energy) || final.energy < best.energy - 1e-12) {
      best.fit <- fit
      best.row <- restart.rows[[restart]]
      best.restart <- restart
      best.energy <- final.energy
    }
  }

  best.fit$restart_summary <- do.call(rbind, restart.rows)
  best.fit$best_restart <- as.integer(best.restart)
  best.fit$best_restart_row <- best.row
  best.fit$vertex_ids <- as.integer(vertex.ids)
  best.fit$coords_full <- grip.geodesic.misf.partial.coords(best.fit$coords, vertex.ids, full.n)
  best.fit
}

grip.geodesic.misf.select.anchors <- function(prepared,
                                              coords,
                                              vertex,
                                              level = NULL,
                                              anchor_policy = c(
                                                "prev_level_first",
                                                "prev_level_distance_band",
                                                "prev_level_spread"
                                              ),
                                              anchor_count = NULL) {
  prepared <- grip.validate.misf.geodesic.prepared(prepared)
  coords <- grip.geodesic.misf.validate.partial.coords(coords, prepared)
  anchor_policy <- match.arg(anchor_policy)
  vertex <- grip.validate.count(vertex, "vertex")
  if (vertex > prepared$n) {
    stop("vertex must be a valid 1-based vertex id from prepared")
  }

  if (is.null(level)) {
    level <- prepared$misf$vertex_depth[[vertex]]
  }
  level.index <- grip.geodesic.misf.level.to.index(prepared$misf, level)
  level.id <- as.integer(level.index - 1L)
  if (level.id >= prepared$misf$misf_height) {
    stop("cannot insert a top-level MISF vertex; it is already handled by the top-level solve")
  }

  candidate.ids <- grip.geodesic.misf.previous.level.vertices(prepared, level.id)
  placed <- rowSums(is.finite(coords[candidate.ids, , drop = FALSE])) == ncol(coords)
  candidate.ids <- candidate.ids[placed]
  candidate.ids <- candidate.ids[candidate.ids != vertex]
  if (!length(candidate.ids)) {
    stop("no placed previous-level anchors are available for this vertex")
  }

  if (is.null(anchor_count)) {
    anchor_count <- grip.geodesic.misf.default.anchor.count(ncol(coords))
  } else {
    anchor_count <- grip.validate.misf.count(anchor_count, "anchor_count", lower = 1L)
  }
  anchor_count <- min(anchor_count, length(candidate.ids))

  candidate.dist <- as.double(prepared$distance_matrix[vertex, candidate.ids])
  ord.near <- order(candidate.dist, candidate.ids)
  selected <- switch(
    anchor_policy,
    prev_level_first = candidate.ids[ord.near[seq_len(anchor_count)]],
    prev_level_distance_band = {
      band.order <- grip.geodesic.misf.distance.band.order(anchor_count, length(candidate.ids))
      candidate.ids[ord.near[band.order]]
    },
    prev_level_spread = grip.geodesic.misf.spread.order(
      candidate_ids = candidate.ids,
      candidate_coords = coords[candidate.ids, , drop = FALSE],
      candidate_distances = candidate.dist,
      count = anchor_count
    )
  )
  selected <- as.integer(selected)
  selected.dist <- as.double(prepared$distance_matrix[vertex, selected])

  list(
    vertex = vertex,
    level = level.id,
    anchor_policy = anchor_policy,
    anchor_ids = selected,
    anchor_coords = coords[selected, , drop = FALSE],
    anchor_distances = selected.dist,
    candidate_ids = as.integer(candidate.ids)
  )
}

grip.geodesic.misf.insert.vertex <- function(prepared,
                                             coords = NULL,
                                             vertex,
                                             level = NULL,
                                             anchor_policy = c(
                                               "prev_level_first",
                                               "prev_level_distance_band",
                                               "prev_level_spread"
                                             ),
                                             anchor_count = NULL,
                                             anchor_weight_mode = c(
                                               "inverse_graph_distance_sq",
                                               "uniform"
                                             ),
                                             max_iter = 64L,
                                             initial_step = 1.0,
                                             step_shrink = 0.5,
                                             armijo_factor = 1e-4,
                                             grad_tol = 1e-8,
                                             min_step = 1e-8) {
  prepared <- grip.validate.misf.geodesic.prepared(prepared)
  coords <- grip.geodesic.misf.validate.partial.coords(coords, prepared)
  anchor_policy <- match.arg(anchor_policy)
  anchor_weight_mode <- match.arg(anchor_weight_mode)
  vertex <- grip.validate.count(vertex, "vertex")
  if (vertex > prepared$n) {
    stop("vertex must be a valid 1-based vertex id from prepared")
  }
  if (is.null(level)) {
    level <- prepared$misf$vertex_depth[[vertex]]
  }
  selection <- grip.geodesic.misf.select.anchors(
    prepared = prepared,
    coords = coords,
    vertex = vertex,
    level = level,
    anchor_policy = anchor_policy,
    anchor_count = anchor_count
  )
  anchor.weights <- grip.geodesic.misf.anchor.weights(
    selection$anchor_distances,
    mode = anchor_weight_mode
  )
  init.coord <- if (all(is.finite(coords[vertex, ]))) {
    as.double(coords[vertex, ])
  } else {
    NULL
  }
  fit <- grip_geodesic_misf_insert_vertex_cpp(
    anchor_coords = selection$anchor_coords,
    anchor_distance = selection$anchor_distances,
    anchor_weights = anchor.weights,
    init_coord = init.coord,
    max_iter = as.integer(max_iter),
    initial_step = initial_step,
    step_shrink = step_shrink,
    armijo_factor = armijo_factor,
    grad_tol = grad_tol,
    min_step = min_step
  )

  coords[vertex, ] <- as.double(fit$coord)
  list(
    vertex = vertex,
    level = as.integer(selection$level),
    coords = coords,
    coord = as.double(fit$coord),
    initial_coord = as.double(fit$initial_coord),
    anchor_ids = selection$anchor_ids,
    anchor_coords = selection$anchor_coords,
    anchor_distances = selection$anchor_distances,
    anchor_weights = anchor.weights,
    anchor_policy = anchor_policy,
    anchor_weight_mode = anchor_weight_mode,
    objective = as.double(fit$objective),
    initial_objective = as.double(fit$initial_objective),
    grad_norm = as.double(fit$grad_norm),
    iterations = as.integer(fit$iterations),
    converged = isTRUE(fit$converged),
    accepted_step = as.double(fit$accepted_step)
  )
}

grip.geodesic.misf.insert.level <- function(prepared,
                                            coords = NULL,
                                            level = NULL,
                                            anchor_policy = c(
                                              "prev_level_first",
                                              "prev_level_distance_band",
                                              "prev_level_spread"
                                            ),
                                            anchor_count = NULL,
                                            anchor_weight_mode = c(
                                              "inverse_graph_distance_sq",
                                              "uniform"
                                            ),
                                            max_iter = 64L,
                                            initial_step = 1.0,
                                            step_shrink = 0.5,
                                            armijo_factor = 1e-4,
                                            grad_tol = 1e-8,
                                            min_step = 1e-8) {
  prepared <- grip.validate.misf.geodesic.prepared(prepared)
  coords <- grip.geodesic.misf.validate.partial.coords(coords, prepared)
  anchor_policy <- match.arg(anchor_policy)
  anchor_weight_mode <- match.arg(anchor_weight_mode)

  if (is.null(level)) {
    level <- prepared$top_level_level - 1L
  }
  level.index <- grip.geodesic.misf.level.to.index(prepared$misf, level)
  level.id <- as.integer(level.index - 1L)
  vertices <- grip.geodesic.misf.level.insert.vertices(prepared, level.id)
  if (!length(vertices)) {
    return(list(
      coords = coords,
      level = level.id,
      inserted_vertices = integer(0L),
      vertex_results = list(),
      vertex_trace = data.frame()
    ))
  }

  vertex.results <- vector("list", length(vertices))
  for (i in seq_along(vertices)) {
    vertex.results[[i]] <- grip.geodesic.misf.insert.vertex(
      prepared = prepared,
      coords = coords,
      vertex = vertices[[i]],
      level = level.id,
      anchor_policy = anchor_policy,
      anchor_count = anchor_count,
      anchor_weight_mode = anchor_weight_mode,
      max_iter = max_iter,
      initial_step = initial_step,
      step_shrink = step_shrink,
      armijo_factor = armijo_factor,
      grad_tol = grad_tol,
      min_step = min_step
    )
    coords <- vertex.results[[i]]$coords
  }

  vertex.trace <- do.call(rbind, lapply(vertex.results, function(result) {
    data.frame(
      level = result$level,
      vertex = result$vertex,
      anchor_count = length(result$anchor_ids),
      objective = result$objective,
      initial_objective = result$initial_objective,
      grad_norm = result$grad_norm,
      iterations = result$iterations,
      converged = result$converged,
      stringsAsFactors = FALSE
    )
  }))

  list(
    coords = coords,
    level = level.id,
    inserted_vertices = as.integer(vertices),
    vertex_results = vertex.results,
    vertex_trace = vertex.trace
  )
}

grip.geodesic.misf.insert.all.levels <- function(prepared,
                                                 coords = NULL,
                                                 anchor_policy = c(
                                                   "prev_level_first",
                                                   "prev_level_distance_band",
                                                   "prev_level_spread"
                                                 ),
                                                 anchor_count = NULL,
                                                 anchor_weight_mode = c(
                                                   "inverse_graph_distance_sq",
                                                   "uniform"
                                                 ),
                                                 max_iter = 64L,
                                                 initial_step = 1.0,
                                                 step_shrink = 0.5,
                                                 armijo_factor = 1e-4,
                                                 grad_tol = 1e-8,
                                                 min_step = 1e-8) {
  prepared <- grip.validate.misf.geodesic.prepared(prepared)
  coords <- grip.geodesic.misf.validate.partial.coords(coords, prepared)
  anchor_policy <- match.arg(anchor_policy)
  anchor_weight_mode <- match.arg(anchor_weight_mode)

  if (prepared$top_level_level <= 0L) {
    return(list(
      coords = coords,
      level_results = list(),
      level_trace = data.frame()
    ))
  }

  level.ids <- seq.int(from = prepared$top_level_level - 1L, to = 0L, by = -1L)
  level.results <- vector("list", length(level.ids))
  for (i in seq_along(level.ids)) {
    level.results[[i]] <- grip.geodesic.misf.insert.level(
      prepared = prepared,
      coords = coords,
      level = level.ids[[i]],
      anchor_policy = anchor_policy,
      anchor_count = anchor_count,
      anchor_weight_mode = anchor_weight_mode,
      max_iter = max_iter,
      initial_step = initial_step,
      step_shrink = step_shrink,
      armijo_factor = armijo_factor,
      grad_tol = grad_tol,
      min_step = min_step
    )
    coords <- level.results[[i]]$coords
  }

  level.trace <- do.call(rbind, lapply(level.results, function(result) {
    if (!nrow(result$vertex_trace)) {
      return(NULL)
    }
    data.frame(
      level = result$level,
      inserted = nrow(result$vertex_trace),
      mean_objective = mean(result$vertex_trace$objective),
      max_grad_norm = max(result$vertex_trace$grad_norm),
      all_converged = all(result$vertex_trace$converged),
      stringsAsFactors = FALSE
    )
  }))

  list(
    coords = coords,
    level_results = level.results,
    level_trace = if (is.null(level.trace)) data.frame() else level.trace
  )
}

grip.geodesic.misf.active.level.vertices <- function(prepared, level) {
  prepared <- grip.validate.misf.geodesic.prepared(prepared)
  level.index <- grip.geodesic.misf.level.to.index(prepared$misf, level)
  as.integer(prepared$misf$levels[[level.index]])
}

grip.geodesic.misf.induced_active_graph <- function(prepared, active_vertices) {
  prepared <- grip.validate.misf.geodesic.prepared(prepared)
  active.vertices <- as.integer(active_vertices)
  if (!length(active.vertices)) {
    stop("active_vertices must contain at least one vertex id")
  }
  if (any(!is.finite(active.vertices)) || any(active.vertices < 1L) || any(active.vertices > prepared$n)) {
    stop("active_vertices must be valid 1-based vertex ids from prepared")
  }
  if (anyDuplicated(active.vertices)) {
    stop("active_vertices must be unique")
  }

  local.count <- length(active.vertices)
  global.to.local <- integer(prepared$n)
  global.to.local[active.vertices] <- seq_len(local.count)
  active.set <- logical(prepared$n)
  active.set[active.vertices] <- TRUE

  adj.list <- vector("list", local.count)
  weight.list <- vector("list", local.count)
  for (local.i in seq_along(active.vertices)) {
    global.i <- active.vertices[[local.i]]
    nbrs <- as.integer(prepared$adj_list[[global.i]])
    weights <- as.double(prepared$weight_list[[global.i]])
    keep <- active.set[nbrs]
    adj.list[[local.i]] <- if (any(keep)) as.integer(global.to.local[nbrs[keep]]) else integer(0L)
    weight.list[[local.i]] <- if (any(keep)) as.double(weights[keep]) else numeric(0L)
  }

  sorted <- grip.sort.adj.with.weights(adj.list, weight.list)
  adj.list <- sorted$adj_list
  weight.list <- sorted$weight_list
  edges <- grip.edges.from.adj.list(adj.list)
  edge.targets <- grip.edge.weights.from.adj.list(adj.list, weight.list)

  list(
    n = local.count,
    vertex_ids = active.vertices,
    global_to_local = global.to.local,
    local_to_global = active.vertices,
    adj_list = adj.list,
    weight_list = weight.list,
    edges = edges,
    edge_targets = edge.targets
  )
}

grip.geodesic.misf.prepare.active.level <- function(prepared,
                                                    active_vertices,
                                                    local_nbrs = 8L,
                                                    landmark_count = 4L,
                                                    pair_mode = c("sparse", "full")) {
  prepared <- grip.validate.misf.geodesic.prepared(prepared)
  pair_mode <- match.arg(pair_mode)
  local_nbrs <- grip.validate.count(local_nbrs, "local_nbrs")
  landmark_count <- grip.validate.count(landmark_count, "landmark_count")
  active.vertices <- as.integer(active_vertices)

  active.graph <- grip.geodesic.misf.induced_active_graph(prepared, active.vertices)
  active.comp <- grip.connected.components(active.graph$adj_list, active.graph$n)
  active.distance <- prepared$distance_matrix[active.vertices, active.vertices, drop = FALSE]
  pair.matrix.local <- if (identical(pair_mode, "full")) {
    grip.full.geodesic.kk.pair.matrix(length(active.vertices))
  } else {
    grip.landmark.geodesic.kk.pair.matrix(
      dist.matrix = active.distance,
      local_nbrs = local_nbrs,
      landmark_count = landmark_count
    )
  }
  pair.matrix <- if (nrow(pair.matrix.local) == 0L) {
    matrix(integer(), ncol = 2L)
  } else {
    cbind(
      as.integer(active.vertices[pair.matrix.local[, 1L]]),
      as.integer(active.vertices[pair.matrix.local[, 2L]])
    )
  }
  cache <- grip.build.geodesic.mds.path.cache(
    pair.matrix = pair.matrix,
    adj.list = prepared$adj_list,
    weight.list = prepared$weight_list,
    dist.matrix = prepared$distance_matrix,
    parents = prepared$parents,
    tie_mode = prepared$tie_mode,
    cache_engine = "r"
  )
  flat.cache <- grip.flatten.geodesic.path.cache(
    path.edges = cache$path_edges,
    path.edge.weights = cache$path_edge_weights
  )

  out <- list(
    n = prepared$n,
    edges = prepared$edges,
    edge_targets = prepared$edge_targets,
    adj_list = prepared$adj_list,
    weight_list = prepared$weight_list,
    pair_matrix = pair.matrix,
    pair_matrix_local = pair.matrix.local,
    pair_graph_distance = cache$pair_graph_distance,
    path_vertices = cache$path_vertices,
    path_edges = cache$path_edges,
    path_edge_weights = cache$path_edge_weights,
    pair_path_count_log = cache$pair_path_count_log,
    flat_pair_edge_offsets = flat.cache$flat_pair_edge_offsets,
    flat_edge_u = flat.cache$flat_edge_u,
    flat_edge_v = flat.cache$flat_edge_v,
    flat_edge_coeff = flat.cache$flat_edge_coeff,
    graph_diameter = prepared$graph_diameter,
    distance_matrix = prepared$distance_matrix,
    pair_mode = if (identical(pair_mode, "full")) "all_pairs" else "misf_sparse",
    graph_build_mode = "misf_active_level",
    tie_mode = prepared$tie_mode,
    active_vertex_ids = active.graph$vertex_ids,
    global_to_local = active.graph$global_to_local,
    local_to_global = active.graph$local_to_global,
    local_nbrs = local_nbrs,
    landmark_count = landmark_count,
    active_component = as.integer(active.comp)
  )
  class(out) <- c("grip_gmds_prepared", "grip_gkk_prepared", "grip_geodesic_kk_prepared", "list")
  out
}

grip.geodesic.misf.build.level.pairs <- function(prepared,
                                                 level = NULL,
                                                 local_nbrs = 8L,
                                                 landmark_count = 4L,
                                                 pair_mode = c("sparse", "full")) {
  prepared <- grip.validate.misf.geodesic.prepared(prepared)
  pair_mode <- match.arg(pair_mode)
  if (is.null(level)) {
    level <- prepared$top_level_level
  }
  active.vertices <- grip.geodesic.misf.active.level.vertices(prepared, level)
  active.prepared <- grip.geodesic.misf.prepare.active.level(
    prepared = prepared,
    active_vertices = active.vertices,
    local_nbrs = local_nbrs,
    landmark_count = landmark_count,
    pair_mode = pair_mode
  )
  global.pair.matrix <- if (nrow(active.prepared$pair_matrix) == 0L) {
    matrix(integer(), ncol = 2L)
  } else if (!is.null(active.prepared$pair_matrix_local)) {
    active.prepared$pair_matrix
  } else {
    active.prepared$pair_matrix
  }
  list(
    level = as.integer(level),
    active_vertices = active.vertices,
    active_prepared = active.prepared,
    pair_matrix = if (!is.null(active.prepared$pair_matrix_local)) active.prepared$pair_matrix_local else active.prepared$pair_matrix,
    global_pair_matrix = matrix(as.integer(global.pair.matrix), ncol = 2L)
  )
}

grip.geodesic.misf.refine.level <- function(prepared,
                                            coords,
                                            level = NULL,
                                            local_nbrs = 8L,
                                            landmark_count = 4L,
                                            pair_mode = c("sparse", "full"),
                                            anchor_weight = 0.05,
                                            anchor_weight_end = anchor_weight,
                                            continuation = c("constant", "linear", "geometric"),
                                            max_iter = 8L,
                                            engine = c("cpp", "r"),
                                            edge_length_epsilon = 1e-8,
                                            initial_step = 1.0,
                                            step_shrink = 0.5,
                                            armijo_factor = 1e-4,
                                            grad_tol = 1e-8,
                                            min_step = 1e-8,
                                            n_threads = 0L,
                                            recenter = TRUE,
                                            return_trace = TRUE) {
  prepared <- grip.validate.misf.geodesic.prepared(prepared)
  coords <- grip.validate.coords(coords)
  if (nrow(coords) != prepared$n) {
    stop("nrow(coords) must match prepared$n")
  }
  pair_mode <- match.arg(pair_mode)
  continuation <- match.arg(continuation)
  engine <- match.arg(engine)
  if (is.null(level)) {
    level <- prepared$top_level_level
  }
  level <- as.integer(level)

  built <- grip.geodesic.misf.build.level.pairs(
    prepared = prepared,
    level = level,
    local_nbrs = local_nbrs,
    landmark_count = landmark_count,
    pair_mode = pair_mode
  )
  active.vertices <- built$active_vertices

  anchor.coords <- coords
  anchor.vertex.weight <- numeric(nrow(coords))
  inactive.vertices <- setdiff(seq_len(prepared$n), active.vertices)
  if (length(inactive.vertices) > 0L) {
    anchor.vertex.weight[inactive.vertices] <- 1
  }
  if (level < prepared$top_level_level) {
    pinned.global <- prepared$misf$levels[[level + 2L]]
    anchor.vertex.weight[pinned.global] <- 1
  }
  use.anchor <- any(anchor.vertex.weight > 0)

  before <- grip.score.geodesic.mds(
    coords = coords,
    prepared = built$active_prepared
  )
  fit <- grip.optimize.geodesic.mds(
    coords = coords,
    prepared = built$active_prepared,
    init = "user",
    anchor_mode = if (use.anchor) "user" else "none",
    anchor_coords = if (use.anchor) anchor.coords else NULL,
    anchor_weight = if (use.anchor) anchor_weight else 0,
    anchor_weight_end = if (use.anchor) anchor_weight_end else 0,
    anchor_vertex_weight = if (use.anchor) anchor.vertex.weight else NULL,
    continuation = continuation,
    engine = engine,
    max_iter = max_iter,
    edge_length_epsilon = edge_length_epsilon,
    initial_step = initial_step,
    step_shrink = step_shrink,
    armijo_factor = armijo_factor,
    grad_tol = grad_tol,
    min_step = min_step,
    n_threads = n_threads,
    recenter = recenter,
    return_trace = return_trace
  )
  after <- grip.score.geodesic.mds(
    coords = fit$coords,
    prepared = built$active_prepared
  )
  coords <- fit$coords

  list(
    level = level,
    active_vertices = active.vertices,
    pinned_vertices = if (use.anchor) which(anchor.vertex.weight > 0) else integer(0L),
    active_prepared = built$active_prepared,
    pair_matrix = built$pair_matrix,
    global_pair_matrix = built$global_pair_matrix,
    coords = coords,
    fit = fit,
    before = before,
    after = after,
    anchor_vertex_weight = anchor.vertex.weight
  )
}

grip.geodesic.misf.refine.all.levels <- function(prepared,
                                                 coords,
                                                 local_nbrs = 8L,
                                                 landmark_count = 4L,
                                                 pair_mode = c("sparse", "full"),
                                                 anchor_weight = 0.05,
                                                 anchor_weight_end = anchor_weight,
                                                 continuation = c("constant", "linear", "geometric"),
                                                 max_iter = 8L,
                                                 engine = c("cpp", "r"),
                                                 edge_length_epsilon = 1e-8,
                                                 initial_step = 1.0,
                                                 step_shrink = 0.5,
                                                 armijo_factor = 1e-4,
                                                 grad_tol = 1e-8,
                                                 min_step = 1e-8,
                                                 n_threads = 0L,
                                                 recenter = TRUE,
                                                 return_trace = TRUE) {
  prepared <- grip.validate.misf.geodesic.prepared(prepared)
  coords <- grip.validate.coords(coords)
  pair_mode <- match.arg(pair_mode)
  continuation <- match.arg(continuation)
  engine <- match.arg(engine)

  level.ids <- seq.int(from = prepared$top_level_level, to = 0L, by = -1L)
  level.results <- vector("list", length(level.ids))
  for (i in seq_along(level.ids)) {
    level.results[[i]] <- grip.geodesic.misf.refine.level(
      prepared = prepared,
      coords = coords,
      level = level.ids[[i]],
      local_nbrs = local_nbrs,
      landmark_count = landmark_count,
      pair_mode = pair_mode,
      anchor_weight = anchor_weight,
      anchor_weight_end = anchor_weight_end,
      continuation = continuation,
      max_iter = max_iter,
      engine = engine,
      edge_length_epsilon = edge_length_epsilon,
      initial_step = initial_step,
      step_shrink = step_shrink,
      armijo_factor = armijo_factor,
      grad_tol = grad_tol,
      min_step = min_step,
      n_threads = n_threads,
      recenter = recenter,
      return_trace = return_trace
    )
    coords <- level.results[[i]]$coords
  }

  level.trace <- do.call(rbind, lapply(level.results, function(result) {
    data.frame(
      level = result$level,
      active_n = length(result$active_vertices),
      pair_n = nrow(result$pair_matrix),
      before_energy = result$before$gmds.energy[[1L]],
      after_energy = result$after$gmds.energy[[1L]],
      after_stress = result$after$gmds.stress[[1L]],
      stringsAsFactors = FALSE
    )
  }))

  list(
    coords = coords,
    level_results = level.results,
    level_trace = level.trace
  )
}

grip.geodesic.misf.final.polish <- function(prepared,
                                            coords,
                                            max_iter = 8L,
                                            engine = c("cpp", "r"),
                                            edge_length_epsilon = 1e-8,
                                            initial_step = 1.0,
                                            step_shrink = 0.5,
                                            armijo_factor = 1e-4,
                                            grad_tol = 1e-8,
                                            min_step = 1e-8,
                                            n_threads = 0L,
                                            recenter = TRUE,
                                            return_trace = TRUE) {
  prepared <- grip.validate.misf.geodesic.prepared(prepared)
  coords <- grip.validate.coords(coords)
  if (nrow(coords) != prepared$n) {
    stop("nrow(coords) must match prepared$n")
  }
  engine <- match.arg(engine)
  fit <- grip.optimize.geodesic.mds(
    coords = coords,
    prepared = prepared,
    init = "user",
    anchor_mode = "none",
    engine = engine,
    max_iter = max_iter,
    edge_length_epsilon = edge_length_epsilon,
    initial_step = initial_step,
    step_shrink = step_shrink,
    armijo_factor = armijo_factor,
    grad_tol = grad_tol,
    min_step = min_step,
    n_threads = n_threads,
    recenter = recenter,
    return_trace = return_trace
  )
  fit
}

#' Prepare a MISF-based multiscale GMDS object
#'
#' `grip.prepare.misf.geodesic.mds()` augments the graph-first GMDS prepared
#' object with the maximal independent set filtration extracted from GRIP and a
#' restartable pure-GMDS solve on the coarsest MISF level. This is the Phase 2
#' preparation layer for the new MISF-based GMDS initializer.
#'
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#' @param n Number of vertices. If omitted with `adj_list`, defaults to
#'   `length(adj_list)`. If omitted with `edges`, defaults to `max(edges)`.
#' @param adj_list Adjacency list (1-based) for an undirected graph.
#' @param weight_list Optional parallel list of positive edge weights.
#' @param edge_weights Optional positive edge-weight vector parallel to
#'   `edges`.
#' @param tie_mode Shortest-path aggregation mode inherited by both the full
#'   graph cache and the top-level coarse graph.
#' @param num_init Target top-level active-set size passed to `grip.build.misf`.
#' @param num_nbrs Per-level local-neighborhood schedule metadata passed to
#'   `grip.build.misf`.
#' @param dim Target coarse-level embedding dimension used by the top-level
#'   pure-GMDS solve.
#' @param top_level_mode Either `"solve"` to run the coarse pure-GMDS solve
#'   immediately, or `"skip"` to prepare the MISF object without solving it.
#' @param top_level_restarts Number of random restarts used by the coarse solve.
#' @param top_level_max_iter Maximum number of pure-GMDS iterations per restart
#'   on the coarse graph.
#' @param top_level_engine Optimization engine used by the coarse solve.
#' @param seed Optional integer seed reused for both the MISF extraction and the
#'   top-level restart family.
#'
#' @return An object of class `"grip_misf_gmds_prepared"` layered on top of the
#'   graph-first GMDS prepared structure, with added MISF metadata and the
#'   optional top-level fit.
#' @examples
#' edges <- edges.mesh(4, 4)
#' prepared <- grip.prepare.misf.geodesic.mds(
#'   edges = edges,
#'   n = 16,
#'   tie_mode = "average",
#'   num_init = 4,
#'   top_level_restarts = 2,
#'   top_level_max_iter = 2,
#'   seed = 1
#' )
#' prepared$top_level_vertices
#' @export
grip.prepare.misf.geodesic.mds <- function(edges = NULL,
                                           n = NULL,
                                           adj_list = NULL,
                                           weight_list = NULL,
                                           edge_weights = NULL,
                                           tie_mode = c("single", "average"),
                                           num_init = 24L,
                                           num_nbrs = 20L,
                                           dim = 2L,
                                           top_level_mode = c("solve", "skip"),
                                           top_level_restarts = 8L,
                                           top_level_max_iter = 16L,
                                           top_level_engine = c("cpp", "r"),
                                           seed = 6L) {
  tie_mode <- match.arg(tie_mode)
  top_level_mode <- match.arg(top_level_mode)
  top_level_engine <- match.arg(top_level_engine)
  dim <- grip.validate.count(dim, "dim")
  if (!(dim %in% c(2L, 3L))) {
    stop("dim must be 2 or 3")
  }
  top_level_restarts <- grip.validate.misf.count(top_level_restarts, "top_level_restarts", lower = 1L)
  grip.validate.scalar(top_level_max_iter, "top_level_max_iter", lower = 0)
  top_level_max_iter <- as.integer(round(top_level_max_iter))
  if (!is.null(seed)) {
    seed <- grip.validate.count(seed, "seed")
  }

  prepared <- grip.prepare.graph.geodesic.mds(
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
  top.level.index <- length(misf$levels)
  top.level.id <- as.integer(top.level.index - 1L)
  top.level.vertices <- as.integer(misf$levels[[top.level.index]])
  top.level.graph <- grip.geodesic.misf.induced_level_graph(
    prepared = prepared,
    vertex_ids = top.level.vertices,
    level = top.level.id
  )
  top.level.prepared <- grip.prepare.graph.geodesic.mds(
    edges = top.level.graph$edges,
    n = top.level.graph$n,
    edge_weights = top.level.graph$edge_weights,
    tie_mode = tie_mode
  )

  prepared$misf <- misf
  prepared$level_vertices <- misf$levels
  prepared$active_levels <- misf$levels
  prepared$insertion_order <- misf$mish_order
  prepared$top_level_index <- top.level.index
  prepared$top_level_level <- top.level.id
  prepared$top_level_vertices <- top.level.vertices
  prepared$top_level_graph <- top.level.graph
  prepared$top_level_prepared <- top.level.prepared
  prepared$top_level_dim <- as.integer(dim)
  prepared$top_level_mode <- top_level_mode
  prepared$top_level_restarts <- as.integer(top_level_restarts)
  prepared$top_level_max_iter <- as.integer(top_level_max_iter)
  prepared$top_level_engine <- top_level_engine
  prepared$multiscale_mode <- "misf"
  prepared$insertion_anchor_policy <- "prev_level_spread"
  prepared$insertion_anchor_count <- grip.geodesic.misf.default.anchor.count(dim)
  prepared$insertion_anchor_weight_mode <- "inverse_graph_distance_sq"
  prepared$insertion_max_iter <- 64L
  prepared$top_level_fit <- NULL
  class(prepared) <- c("grip_misf_gmds_prepared", class(prepared))

  if (identical(top_level_mode, "solve")) {
    prepared$top_level_fit <- grip.geodesic.misf.solve.top.level(
      prepared = prepared,
      dim = dim,
      n_restarts = top_level_restarts,
      max_iter = top_level_max_iter,
      engine = top_level_engine,
      seed = seed
    )
  }
  prepared
}
