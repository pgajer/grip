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
