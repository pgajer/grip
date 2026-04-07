grip.validate.misf.geodesic.prepared <- function(prepared, coords = NULL) {
  if (!(inherits(prepared, "grip_misf_gmds_prepared") ||
        inherits(prepared, "grip_misf_gkk_prepared"))) {
    stop(
      "prepared must be an object from grip.prepare.misf.geodesic.mds() ",
      "or grip.prepare.misf.geodesic.kk()"
    )
  }
  grip.validate.geodesic.mds.prepared(prepared, coords = coords)
}

grip.validate.misf.geodesic.fit <- function(fit) {
  if (!inherits(fit, "grip_misf_gmds_fit")) {
    stop("fit must be an object from grip.optimize.misf.geodesic.mds()")
  }
  fit
}

grip.resolve.misf.geodesic.prepared <- function(prepared = NULL,
                                                edges = NULL,
                                                n = NULL,
                                                adj_list = NULL,
                                                weight_list = NULL,
                                                edge_weights = NULL,
                                                tie_mode = NULL,
                                                num_init = 24L,
                                                num_nbrs = 20L,
                                                dim = NULL,
                                                top_level_init = c("geometric", "random"),
                                                top_level_restarts = 8L,
                                                top_level_max_iter = 16L,
                                                top_level_engine = c("cpp", "r"),
                                                seed = 6L) {
  top_level_init <- match.arg(top_level_init)
  top_level_engine <- match.arg(top_level_engine)
  if (is.null(prepared)) {
    resolved.dim <- if (is.null(dim)) 2L else grip.validate.count(dim, "dim")
    resolved.tie.mode <- if (is.null(tie_mode)) "average" else {
      match.arg(tie_mode, c("single", "average"))
    }
    return(grip.prepare.misf.geodesic.mds(
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
      top_level_init = top_level_init,
      top_level_restarts = top_level_restarts,
      top_level_max_iter = top_level_max_iter,
      top_level_engine = top_level_engine,
      seed = seed
    ))
  }
  if (inherits(prepared, "grip_misf_gmds_prepared")) {
    return(grip.validate.misf.geodesic.prepared(prepared))
  }
  prepared <- grip.validate.geodesic.mds.prepared(prepared)
  resolved.dim <- if (is.null(dim)) 2L else grip.validate.count(dim, "dim")
  resolved.tie.mode <- if (is.null(tie_mode)) {
    if (!is.null(prepared$tie_mode)) prepared$tie_mode else "average"
  } else {
    match.arg(tie_mode, c("single", "average"))
  }
  grip.prepare.misf.geodesic.mds(
    n = prepared$n,
    adj_list = prepared$adj_list,
    weight_list = prepared$weight_list,
    tie_mode = resolved.tie.mode,
    num_init = num_init,
    num_nbrs = num_nbrs,
    dim = resolved.dim,
    top_level_mode = "skip",
    top_level_init = top_level_init,
    top_level_restarts = top_level_restarts,
    top_level_max_iter = top_level_max_iter,
    top_level_engine = top_level_engine,
    seed = seed
  )
}

grip.geodesic.misf.expand.top.level.frames <- function(top_level_fit, n) {
  if (is.null(top_level_fit$frames) || !length(top_level_fit$frames)) {
    return(list())
  }
  lapply(top_level_fit$frames, function(frame) {
    grip.geodesic.misf.partial.coords(
      coords = frame,
      vertex_ids = top_level_fit$vertex_ids,
      n = n
    )
  })
}

grip.geodesic.misf.frame.count <- function(frames) {
  if (is.null(frames)) {
    return(0L)
  }
  as.integer(length(frames))
}

grip.geodesic.misf.collect.level.frames <- function(level_results) {
  if (!length(level_results)) {
    return(list())
  }
  stats::setNames(
    lapply(level_results, `[[`, "coords"),
    paste0("level_", vapply(level_results, function(result) result$level, integer(1L)))
  )
}

grip.geodesic.misf.collect.insertion.vertex.trace <- function(level_results) {
  if (!length(level_results)) {
    return(data.frame())
  }
  rows <- lapply(level_results, `[[`, "vertex_trace")
  rows <- Filter(function(row) !is.null(row) && nrow(row) > 0L, rows)
  if (!length(rows)) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

grip.geodesic.misf.build.stage.trace <- function(prepared,
                                                 top_level_fit,
                                                 top_level_elapsed,
                                                 top_level_frames,
                                                 insertion,
                                                 insertion_elapsed,
                                                 insertion_frames,
                                                 refinement,
                                                 refinement_elapsed,
                                                 refinement_frames,
                                                 final_polish,
                                                 final_polish_elapsed) {
  rows <- list()

  rows[[length(rows) + 1L]] <- data.frame(
    stage = "top_level",
    level = prepared$top_level_level,
    active_n = length(prepared$top_level_vertices),
    inserted_n = 0L,
    pair_n = length(top_level_fit$prepared$pair_graph_distance),
    energy = top_level_fit$score$gmds.energy[[1L]],
    stress = top_level_fit$score$gmds.stress[[1L]],
    mean_objective = NA_real_,
    max_grad_norm = NA_real_,
    all_converged = NA,
    elapsed_sec = as.double(top_level_elapsed),
    trace_rows = if (is.null(top_level_fit$trace)) 0L else nrow(top_level_fit$trace),
    frame_count = grip.geodesic.misf.frame.count(top_level_frames),
    stringsAsFactors = FALSE
  )

  if (!is.null(insertion$level_trace) && nrow(insertion$level_trace) > 0L) {
    insert.rows <- data.frame(
      stage = "insertion",
      level = insertion$level_trace$level,
      active_n = NA_integer_,
      inserted_n = insertion$level_trace$inserted,
      pair_n = NA_integer_,
      energy = insertion$level_trace$mean_objective,
      stress = NA_real_,
      mean_objective = insertion$level_trace$mean_objective,
      max_grad_norm = insertion$level_trace$max_grad_norm,
      all_converged = insertion$level_trace$all_converged,
      elapsed_sec = NA_real_,
      trace_rows = insertion$level_trace$inserted,
      frame_count = 1L,
      stringsAsFactors = FALSE
    )
    if (nrow(insert.rows) > 0L) {
      insert.rows$elapsed_sec[[nrow(insert.rows)]] <- as.double(insertion_elapsed)
      rows[[length(rows) + 1L]] <- insert.rows
    }
  }

  if (!is.null(refinement$level_trace) && nrow(refinement$level_trace) > 0L) {
    refine.trace.rows <- vapply(refinement$level_results, function(result) {
      if (is.null(result$fit$trace)) 0L else nrow(result$fit$trace)
    }, integer(1L))
    refine.rows <- data.frame(
      stage = "refinement",
      level = refinement$level_trace$level,
      active_n = refinement$level_trace$active_n,
      inserted_n = 0L,
      pair_n = refinement$level_trace$pair_n,
      energy = refinement$level_trace$after_energy,
      stress = refinement$level_trace$after_stress,
      mean_objective = NA_real_,
      max_grad_norm = NA_real_,
      all_converged = NA,
      elapsed_sec = NA_real_,
      trace_rows = refine.trace.rows,
      frame_count = 1L,
      stringsAsFactors = FALSE
    )
    if (nrow(refine.rows) > 0L) {
      refine.rows$elapsed_sec[[nrow(refine.rows)]] <- as.double(refinement_elapsed)
      rows[[length(rows) + 1L]] <- refine.rows
    }
  }

  rows[[length(rows) + 1L]] <- data.frame(
    stage = "final_polish",
    level = 0L,
    active_n = prepared$n,
    inserted_n = 0L,
    pair_n = length(prepared$pair_graph_distance),
    energy = final_polish$score$gmds.energy[[1L]],
    stress = final_polish$score$gmds.stress[[1L]],
    mean_objective = NA_real_,
    max_grad_norm = NA_real_,
    all_converged = NA,
    elapsed_sec = as.double(final_polish_elapsed),
    trace_rows = if (is.null(final_polish$trace)) 0L else nrow(final_polish$trace),
    frame_count = grip.geodesic.misf.frame.count(final_polish$frames),
    stringsAsFactors = FALSE
  )

  do.call(rbind, rows)
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

grip.geodesic.misf.build.layout.graph <- function(distance_matrix, k = 6L) {
  distance_matrix <- as.matrix(distance_matrix)
  n <- nrow(distance_matrix)
  if (!is.numeric(distance_matrix) || n != ncol(distance_matrix)) {
    stop("distance_matrix must be a square numeric matrix")
  }
  if (n <= 1L) {
    return(list(
      n = as.integer(n),
      edges = matrix(integer(), ncol = 2L),
      edge_weights = numeric(0L),
      adj_list = vector("list", n),
      weight_list = vector("list", n)
    ))
  }

  k <- grip.validate.count(k, "k")
  k <- min(k, n - 1L)
  knn.edges <- grip.knn.edge.matrix(distance_matrix, k = k)
  mst.edges <- grip.minimum.spanning.tree.edges(distance_matrix)
  final.edges <- grip.union.edge.matrix(knn.edges, mst.edges)
  edge.weights <- grip.edge.weights.from.distance.matrix(final.edges, distance_matrix)
  built <- grip.build.adj.from.edges(final.edges, n = n, edge_weights = edge.weights)

  list(
    n = as.integer(n),
    edges = final.edges,
    edge_weights = edge.weights,
    adj_list = built$adj_list,
    weight_list = built$weight_list
  )
}

grip.geodesic.misf.local.start.coords <- function(coords,
                                                  active_vertices,
                                                  anchor_vertices,
                                                  seed = NULL) {
  coords <- as.matrix(coords)
  active.vertices <- as.integer(active_vertices)
  anchor.vertices <- as.integer(anchor_vertices)
  dim <- ncol(coords)
  out <- matrix(NA_real_, nrow = length(active.vertices), ncol = dim)
  anchor.local <- match(anchor.vertices, active.vertices)
  out[anchor.local, ] <- coords[anchor.vertices, , drop = FALSE]

  pending.local <- which(!stats::complete.cases(out))
  if (!length(pending.local)) {
    return(out)
  }

  anchor.coords <- out[anchor.local, , drop = FALSE]
  center <- if (nrow(anchor.coords)) {
    colMeans(anchor.coords)
  } else {
    rep(0, dim)
  }
  scale <- if (nrow(anchor.coords) >= 2L) {
    stats::median(stats::dist(anchor.coords))
  } else {
    1
  }
  if (!is.finite(scale) || scale <= 0) {
    scale <- 1
  }
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }
  noise <- matrix(stats::rnorm(length(pending.local) * dim, sd = 0.05 * scale), ncol = dim)
  out[pending.local, ] <- matrix(center, nrow = length(pending.local), ncol = dim, byrow = TRUE) + noise
  storage.mode(out) <- "double"
  out
}

grip.geodesic.misf.align.active.layout.to.anchors <- function(active_coords,
                                                              anchor_local,
                                                              target_anchor_coords,
                                                              allow.reflection = TRUE) {
  active.coords <- as.matrix(active_coords)
  anchor.local <- as.integer(anchor_local)
  target.anchor.coords <- as.matrix(target_anchor_coords)
  if (!length(anchor.local)) {
    return(active.coords)
  }
  if (length(anchor.local) != nrow(target.anchor.coords)) {
    stop("anchor_local and target_anchor_coords must have matching sizes")
  }

  source.anchor <- active.coords[anchor.local, , drop = FALSE]
  src.meta <- grip.normalize.coords.with.meta(source.anchor)
  dst.meta <- grip.normalize.coords.with.meta(target.anchor.coords)
  cross <- t(src.meta$normalized) %*% dst.meta$normalized
  sv <- svd(cross)
  rot <- sv$u %*% t(sv$v)
  if (!allow.reflection && ncol(rot) > 1L && det(rot) < 0) {
    fix <- diag(ncol(rot))
    fix[ncol(fix), ncol(fix)] <- -1
    rot <- sv$u %*% fix %*% t(sv$v)
  }

  src.all <- sweep(active.coords, 2L, src.meta$center, check.margin = FALSE) / src.meta$radius
  aligned <- src.all %*% rot
  sweep(aligned * dst.meta$radius, 2L, dst.meta$center, FUN = "+", check.margin = FALSE)
}

grip.geodesic.misf.place.level.with.layout <- function(prepared,
                                                       coords = NULL,
                                                       level = NULL,
                                                       method = c("kk", "weighted_kk", "fr", "grip", "weighted_grip"),
                                                       layout_k = 6L,
                                                       weighted_preset = NULL,
                                                       grip_args = list(),
                                                       weighted_args = list(),
                                                       fr_niter = 800L,
                                                       seed = NULL) {
  prepared <- grip.validate.misf.geodesic.prepared(prepared)
  coords <- grip.geodesic.misf.validate.partial.coords(coords, prepared)
  method <- match.arg(method)
  layout_k <- grip.validate.count(layout_k, "layout_k")
  fr_niter <- grip.validate.count(fr_niter, "fr_niter")
  if (!is.null(seed)) {
    seed <- grip.validate.count(seed, "seed")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for layout-based MISF placement")
  }

  if (is.null(level)) {
    level <- prepared$top_level_level - 1L
  }
  level.index <- grip.geodesic.misf.level.to.index(prepared$misf, level)
  level.id <- as.integer(level.index - 1L)
  active.vertices <- grip.geodesic.misf.active.level.vertices(prepared, level.id)
  anchor.vertices <- grip.geodesic.misf.previous.level.vertices(prepared, level.id)
  anchor.vertices <- anchor.vertices[stats::complete.cases(coords[anchor.vertices, , drop = FALSE])]
  placed.vertices <- setdiff(active.vertices, anchor.vertices)
  if (!length(placed.vertices)) {
    return(list(
      level = level.id,
      method = method,
      coords = coords,
      active_vertices = active.vertices,
      anchor_vertices = anchor.vertices,
      placed_vertices = integer(0L),
      layout_graph = grip.geodesic.misf.build.layout.graph(
        prepared$distance_matrix[active.vertices, active.vertices, drop = FALSE],
        k = layout_k
      ),
      local_coords = coords[active.vertices, , drop = FALSE],
      aligned_active_coords = coords[active.vertices, , drop = FALSE]
    ))
  }

  active.distance <- prepared$distance_matrix[active.vertices, active.vertices, drop = FALSE]
  layout.graph <- grip.geodesic.misf.build.layout.graph(active.distance, k = layout_k)
  local.start <- grip.geodesic.misf.local.start.coords(
    coords = coords,
    active_vertices = active.vertices,
    anchor_vertices = anchor.vertices,
    seed = seed
  )
  graph.obj <- igraph::graph_from_edgelist(layout.graph$edges, directed = FALSE)

  local.coords <- switch(
    method,
    kk = {
      igraph::layout_with_kk(
        graph.obj,
        coords = local.start,
        dim = ncol(coords)
      )
    },
    weighted_kk = {
      igraph::layout_with_kk(
        graph.obj,
        coords = local.start,
        dim = ncol(coords),
        weights = layout.graph$edge_weights
      )
    },
    fr = {
      igraph::layout_with_fr(
        graph.obj,
        coords = local.start,
        dim = ncol(coords),
        niter = as.integer(fr_niter)
      )
    },
    grip = {
      args <- c(
        list(
          edges = layout.graph$edges,
          n = layout.graph$n,
          dim = ncol(coords),
          seed = seed
        ),
        grip_args
      )
      do.call(grip.layout.globalrep, args)
    },
    weighted_grip = {
      args <- c(
        list(
          edges = layout.graph$edges,
          edge_weights = layout.graph$edge_weights,
          n = layout.graph$n,
          dim = ncol(coords),
          seed = seed
        ),
        if (!is.null(weighted_preset)) list(preset = weighted_preset) else list(),
        weighted_args
      )
      do.call(grip.layout.globalrep.weighted, args)
    }
  )
  local.coords <- as.matrix(local.coords)
  aligned.active <- grip.geodesic.misf.align.active.layout.to.anchors(
    active_coords = local.coords,
    anchor_local = match(anchor.vertices, active.vertices),
    target_anchor_coords = coords[anchor.vertices, , drop = FALSE]
  )

  placed.local <- match(placed.vertices, active.vertices)
  coords[placed.vertices, ] <- aligned.active[placed.local, , drop = FALSE]

  list(
    level = level.id,
    method = method,
    coords = coords,
    active_vertices = active.vertices,
    anchor_vertices = anchor.vertices,
    placed_vertices = as.integer(placed.vertices),
    layout_graph = layout.graph,
    local_coords = local.coords,
    aligned_active_coords = aligned.active
  )
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

grip.geodesic.misf.recenter.coords <- function(coords) {
  coords <- as.matrix(coords)
  if (!nrow(coords) || !ncol(coords)) {
    return(coords)
  }
  sweep(coords, 2L, colMeans(coords), FUN = "-", check.margin = FALSE)
}

grip.geodesic.misf.classical.mds.stats <- function(distance_matrix,
                                                   dim,
                                                   tol = 1e-8) {
  distance_matrix <- as.matrix(distance_matrix)
  dim <- grip.validate.count(dim, "dim")
  if (!(dim %in% c(2L, 3L))) {
    stop("dim must be 2 or 3")
  }
  n <- nrow(distance_matrix)
  if (!n || ncol(distance_matrix) != n) {
    stop("distance_matrix must be a non-empty square matrix")
  }

  d2 <- distance_matrix^2
  row.mean <- rowMeans(d2)
  gram <- -0.5 * (
    d2 -
      outer(row.mean, rep.int(1, n)) -
      outer(rep.int(1, n), row.mean) +
      mean(d2)
  )
  eig <- eigen(gram, symmetric = TRUE)
  eig.values <- as.double(eig$values)
  pos.idx <- which(eig.values > tol)
  use.count <- min(length(pos.idx), dim)
  coords <- matrix(0, nrow = n, ncol = dim)
  if (use.count > 0L) {
    vecs <- eig$vectors[, pos.idx[seq_len(use.count)], drop = FALSE]
    coords[, seq_len(use.count)] <- sweep(
      vecs,
      2L,
      sqrt(eig.values[pos.idx[seq_len(use.count)]]),
      FUN = "*",
      check.margin = FALSE
    )
  }
  coords <- grip.geodesic.misf.recenter.coords(coords)
  storage.mode(coords) <- "double"
  list(
    coords = coords,
    eigenvalues = eig.values,
    positive_rank = length(pos.idx),
    target_rank = min(dim, max(n - 1L, 0L)),
    realized_rank = use.count
  )
}

grip.geodesic.misf.score.seed.metric <- function(distance_matrix,
                                                 dim,
                                                 tol = 1e-8) {
  stats <- grip.geodesic.misf.classical.mds.stats(
    distance_matrix = distance_matrix,
    dim = dim,
    tol = tol
  )
  target.rank <- stats$target_rank
  eig.values <- stats$eigenvalues
  target.eig <- if (target.rank > 0L && length(eig.values) >= target.rank) {
    as.double(eig.values[[target.rank]])
  } else {
    -Inf
  }
  pair.mean <- if (nrow(distance_matrix) > 1L) {
    mean(distance_matrix[upper.tri(distance_matrix)])
  } else {
    0
  }
  list(
    positive_rank = stats$positive_rank,
    target_rank = target.rank,
    target_eig = target.eig,
    pair_mean = as.double(pair.mean),
    coords = stats$coords
  )
}

grip.geodesic.misf.select.seed.vertices <- function(distance_matrix,
                                                    count,
                                                    dim,
                                                    max_combinations = 50000L,
                                                    tol = 1e-8) {
  distance_matrix <- as.matrix(distance_matrix)
  count <- grip.validate.misf.count(count, "count", lower = 1L)
  dim <- grip.validate.count(dim, "dim")
  if (!(dim %in% c(2L, 3L))) {
    stop("dim must be 2 or 3")
  }
  n <- nrow(distance_matrix)
  if (!n || ncol(distance_matrix) != n) {
    stop("distance_matrix must be a non-empty square matrix")
  }
  count <- min(as.integer(count), n)
  if (count <= 1L) {
    return(as.integer(1L))
  }

  fallback <- grip.geodesic.misf.select.spread.seed.vertices(
    distance_matrix = distance_matrix,
    count = count
  )
  combo.count <- choose(n, count)
  if (!is.finite(combo.count) || combo.count <= 1L || combo.count > max_combinations) {
    return(as.integer(fallback))
  }

  combos <- combn(n, count)
  best.idx <- NULL
  best.score <- NULL
  fallback.key <- paste(sort(fallback), collapse = ",")

  for (j in seq_len(ncol(combos))) {
    idx <- as.integer(combos[, j])
    score <- grip.geodesic.misf.score.seed.metric(
      distance_matrix = distance_matrix[idx, idx, drop = FALSE],
      dim = dim,
      tol = tol
    )
    key <- paste(sort(idx), collapse = ",")
    better <- is.null(best.score) ||
      score$positive_rank > best.score$positive_rank ||
      (score$positive_rank == best.score$positive_rank && score$target_eig > best.score$target_eig + tol) ||
      (score$positive_rank == best.score$positive_rank &&
         abs(score$target_eig - best.score$target_eig) <= tol &&
         score$pair_mean > best.score$pair_mean + tol) ||
      (score$positive_rank == best.score$positive_rank &&
         abs(score$target_eig - best.score$target_eig) <= tol &&
         abs(score$pair_mean - best.score$pair_mean) <= tol &&
         identical(key, fallback.key))
    if (better) {
      best.idx <- idx
      best.score <- score
    }
  }

  if (is.null(best.idx)) {
    return(as.integer(fallback))
  }
  as.integer(best.idx)
}

grip.geodesic.misf.select.spread.seed.vertices <- function(distance_matrix,
                                                           count) {
  distance_matrix <- as.matrix(distance_matrix)
  count <- grip.validate.misf.count(count, "count", lower = 1L)
  n <- nrow(distance_matrix)
  if (!n || ncol(distance_matrix) != n) {
    stop("distance_matrix must be a non-empty square matrix")
  }
  count <- min(as.integer(count), n)
  if (count == 1L) {
    return(as.integer(1L))
  }

  work <- distance_matrix
  diag(work) <- -Inf
  pair.idx <- which(work == max(work, na.rm = TRUE), arr.ind = TRUE)[1L, ]
  selected <- unique(as.integer(pair.idx))
  if (length(selected) < 2L) {
    selected <- c(selected, setdiff(seq_len(n), selected)[1L])
  }

  while (length(selected) < count) {
    remaining <- setdiff(seq_len(n), selected)
    min.dist <- vapply(remaining, function(idx) {
      min(distance_matrix[idx, selected, drop = TRUE])
    }, numeric(1L))
    mean.dist <- vapply(remaining, function(idx) {
      mean(distance_matrix[idx, selected, drop = TRUE])
    }, numeric(1L))
    choice <- remaining[[order(-min.dist, -mean.dist, remaining)[1L]]]
    selected <- c(selected, choice)
  }

  as.integer(selected[seq_len(count)])
}

grip.geodesic.misf.embed.small.metric <- function(distance_matrix, dim) {
  distance_matrix <- as.matrix(distance_matrix)
  dim <- grip.validate.count(dim, "dim")
  if (!(dim %in% c(2L, 3L))) {
    stop("dim must be 2 or 3")
  }
  n <- nrow(distance_matrix)
  if (!n || ncol(distance_matrix) != n) {
    stop("distance_matrix must be a non-empty square matrix")
  }

  coords <- matrix(0, nrow = n, ncol = dim)
  if (n == 1L) {
    storage.mode(coords) <- "double"
    return(coords)
  }
  if (n == 2L) {
    span <- as.double(distance_matrix[1L, 2L]) / 2
    coords[1L, 1L] <- -span
    coords[2L, 1L] <- span
    storage.mode(coords) <- "double"
    return(coords)
  }
  if (n == 3L && dim >= 2L && is.finite(distance_matrix[1L, 2L]) &&
      distance_matrix[1L, 2L] > sqrt(.Machine$double.eps)) {
    a <- as.double(distance_matrix[1L, 2L])
    b <- as.double(distance_matrix[1L, 3L])
    c <- as.double(distance_matrix[2L, 3L])
    x3 <- (b * b + a * a - c * c) / (2 * a)
    y3.sq <- max(b * b - x3 * x3, 0)
    coords[2L, 1L] <- a
    coords[3L, 1L] <- x3
    coords[3L, 2L] <- sqrt(y3.sq)
    coords <- grip.geodesic.misf.recenter.coords(coords)
    storage.mode(coords) <- "double"
    return(coords)
  }

  stats <- grip.geodesic.misf.classical.mds.stats(
    distance_matrix = distance_matrix,
    dim = dim
  )
  stats$coords
}

grip.geodesic.misf.jitter.coords <- function(coords,
                                             restart = 1L,
                                             seed = NULL,
                                             scale = 0.05) {
  coords <- as.matrix(coords)
  storage.mode(coords) <- "double"
  restart <- grip.validate.misf.count(restart, "restart", lower = 1L)
  if (restart <= 1L || !length(coords)) {
    return(grip.geodesic.misf.recenter.coords(coords))
  }
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }
  spread <- stats::sd(as.double(coords))
  if (!is.finite(spread) || spread <= 0) {
    spread <- 1.0
  }
  jitter <- matrix(
    stats::rnorm(length(coords), sd = scale * spread),
    nrow = nrow(coords),
    ncol = ncol(coords)
  )
  coords <- coords + jitter
  coords <- grip.geodesic.misf.recenter.coords(coords)
  storage.mode(coords) <- "double"
  coords
}

grip.geodesic.misf.build.geometric.seed.coords <- function(distance_matrix,
                                                           dim,
                                                           vertex_ids = NULL,
                                                           insertion_order = NULL,
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
  distance_matrix <- as.matrix(distance_matrix)
  dim <- grip.validate.count(dim, "dim")
  if (!(dim %in% c(2L, 3L))) {
    stop("dim must be 2 or 3")
  }
  if (!nrow(distance_matrix) || ncol(distance_matrix) != nrow(distance_matrix)) {
    stop("distance_matrix must be a non-empty square matrix")
  }
  anchor_weight_mode <- match.arg(anchor_weight_mode)
  n <- nrow(distance_matrix)
  if (is.null(vertex_ids)) {
    vertex_ids <- seq_len(n)
  } else {
    vertex_ids <- as.integer(vertex_ids)
    if (length(vertex_ids) != n) {
      stop("length(vertex_ids) must match nrow(distance_matrix)")
    }
  }
  if (is.null(anchor_count)) {
    anchor_count <- grip.geodesic.misf.default.anchor.count(dim)
  } else {
    anchor_count <- grip.validate.misf.count(anchor_count, "anchor_count", lower = 1L)
  }

  seed_size <- min(n, dim + 1L)
  seed_local <- grip.geodesic.misf.select.seed.vertices(
    distance_matrix = distance_matrix,
    count = seed_size,
    dim = dim
  )
  coords <- matrix(NA_real_, nrow = n, ncol = dim)
  coords[seed_local, ] <- grip.geodesic.misf.embed.small.metric(
    distance_matrix[seed_local, seed_local, drop = FALSE],
    dim = dim
  )

  if (is.null(insertion_order)) {
    remaining_local <- setdiff(seq_len(n), seed_local)
    if (length(remaining_local)) {
      seed_min <- vapply(remaining_local, function(idx) {
        min(distance_matrix[idx, seed_local, drop = TRUE])
      }, numeric(1L))
      seed_mean <- vapply(remaining_local, function(idx) {
        mean(distance_matrix[idx, seed_local, drop = TRUE])
      }, numeric(1L))
      remaining_local <- remaining_local[order(seed_min, -seed_mean, vertex_ids[remaining_local])]
    }
  } else {
    insertion_order <- as.integer(insertion_order)
    order_global <- insertion_order[insertion_order %in% vertex_ids]
    order_local <- match(order_global, vertex_ids)
    order_local <- order_local[!is.na(order_local)]
    remaining_local <- setdiff(order_local, seed_local)
    if (length(remaining_local) != n - seed_size) {
      remaining_local <- setdiff(seq_len(n), seed_local)
    }
  }

  vertex_rows <- vector("list", length(remaining_local))
  placed_local <- seed_local
  placement_order_local <- seed_local

  for (idx in seq_along(remaining_local)) {
    vertex_local <- remaining_local[[idx]]
    candidate_local <- as.integer(placed_local)
    candidate_dist <- as.double(distance_matrix[vertex_local, candidate_local, drop = TRUE])
    selected_local <- if (length(candidate_local) <= anchor_count) {
      candidate_local
    } else {
      grip.geodesic.misf.spread.order(
        candidate_ids = candidate_local,
        candidate_coords = coords[candidate_local, , drop = FALSE],
        candidate_distances = candidate_dist,
        count = anchor_count
      )
    }
    selected_local <- as.integer(selected_local)
    selected_dist <- as.double(distance_matrix[vertex_local, selected_local, drop = TRUE])
    anchor_weights <- grip.geodesic.misf.anchor.weights(
      selected_dist,
      mode = anchor_weight_mode
    )
    fit <- grip_geodesic_misf_insert_vertex_cpp(
      anchor_coords = coords[selected_local, , drop = FALSE],
      anchor_distance = selected_dist,
      anchor_weights = anchor_weights,
      max_iter = as.integer(max_iter),
      initial_step = initial_step,
      step_shrink = step_shrink,
      armijo_factor = armijo_factor,
      grad_tol = grad_tol,
      min_step = min_step
    )
    coords[vertex_local, ] <- as.double(fit$coord)
    placed_local <- c(placed_local, vertex_local)
    placement_order_local <- c(placement_order_local, vertex_local)
    vertex_rows[[idx]] <- data.frame(
      local_vertex = as.integer(vertex_local),
      vertex = as.integer(vertex_ids[[vertex_local]]),
      placement_step = as.integer(length(placement_order_local)),
      anchor_count = as.integer(length(selected_local)),
      objective = as.double(fit$objective),
      initial_objective = as.double(fit$initial_objective),
      grad_norm = as.double(fit$grad_norm),
      iterations = as.integer(fit$iterations),
      converged = isTRUE(fit$converged),
      stringsAsFactors = FALSE
    )
  }

  coords <- grip.geodesic.misf.recenter.coords(coords)
  storage.mode(coords) <- "double"
  list(
    coords = coords,
    vertex_ids = as.integer(vertex_ids),
    seed_local = as.integer(seed_local),
    seed_vertices = as.integer(vertex_ids[seed_local]),
    placement_order_local = as.integer(placement_order_local),
    placement_order_vertices = as.integer(vertex_ids[placement_order_local]),
    vertex_trace = if (length(vertex_rows)) do.call(rbind, vertex_rows) else data.frame()
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
#' @param init Top-level initialization mode. `"geometric"` seeds the coarse
#'   level from a spread `d+1`-vertex geometric placement and inserts the
#'   remaining coarse vertices before refinement. `"random"` keeps the previous
#'   random-restart behavior.
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
                                               init = c("geometric", "random"),
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
  init <- match.arg(init)
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
  geometric.init <- NULL
  if (identical(init, "geometric")) {
    insertion.order <- if (is.misf.prepared && !is.null(prepared$insertion_order)) {
      prepared$insertion_order[prepared$insertion_order %in% vertex.ids]
    } else {
      vertex.ids
    }
    anchor.count <- if (is.misf.prepared && !is.null(prepared$insertion_anchor_count)) {
      prepared$insertion_anchor_count
    } else {
      grip.geodesic.misf.default.anchor.count(dim)
    }
    anchor.mode <- if (is.misf.prepared && !is.null(prepared$insertion_anchor_weight_mode)) {
      prepared$insertion_anchor_weight_mode
    } else {
      "inverse_graph_distance_sq"
    }
    geometric.init <- grip.geodesic.misf.build.geometric.seed.coords(
      distance_matrix = coarse.prepared$distance_matrix,
      dim = dim,
      vertex_ids = vertex.ids,
      insertion_order = insertion.order,
      anchor_count = anchor.count,
      anchor_weight_mode = anchor.mode
    )
  }

  for (restart in seq_len(n_restarts)) {
    restart.seed <- if (is.null(seed)) NULL else as.integer(seed + restart - 1L)
    init.coords <- if (identical(init, "geometric")) {
      grip.geodesic.misf.jitter.coords(
        geometric.init$coords,
        restart = restart,
        seed = restart.seed
      )
    } else {
      if (!is.null(restart.seed)) {
        set.seed(restart.seed)
      }
      init.coords <- matrix(stats::rnorm(coarse.prepared$n * dim), ncol = dim)
      storage.mode(init.coords) <- "double"
      if (isTRUE(recenter)) {
        init.coords <- sweep(init.coords, 2L, colMeans(init.coords), FUN = "-", check.margin = FALSE)
      }
      init.coords
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
  best.fit$top_level_init <- init
  best.fit$initial_placement <- geometric.init
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

grip.geodesic.misf.insert.all.levels.with.layout <- function(prepared,
                                                             coords = NULL,
                                                             method = c("kk", "weighted_kk", "fr", "grip", "weighted_grip"),
                                                             layout_k = 6L,
                                                             weighted_preset = NULL,
                                                             grip_args = list(),
                                                             weighted_args = list(),
                                                             fr_niter = 800L,
                                                             seed = NULL) {
  prepared <- grip.validate.misf.geodesic.prepared(prepared)
  coords <- grip.geodesic.misf.validate.partial.coords(coords, prepared)
  method <- match.arg(method)
  layout_k <- grip.validate.count(layout_k, "layout_k")
  fr_niter <- grip.validate.count(fr_niter, "fr_niter")
  if (!is.null(seed)) {
    seed <- grip.validate.count(seed, "seed")
  }

  if (prepared$top_level_level <= 0L) {
    return(list(
      coords = coords,
      level_results = list(),
      level_trace = data.frame()
    ))
  }

  level.ids <- seq.int(from = prepared$top_level_level - 1L, to = 0L, by = -1L)
  level.results <- vector("list", length(level.ids))
  level.trace.rows <- vector("list", length(level.ids))
  for (i in seq_along(level.ids)) {
    level.id <- level.ids[[i]]
    level.seed <- if (is.null(seed)) NULL else as.integer(seed + i - 1L)
    level.results[[i]] <- grip.geodesic.misf.place.level.with.layout(
      prepared = prepared,
      coords = coords,
      level = level.id,
      method = method,
      layout_k = layout_k,
      weighted_preset = weighted_preset,
      grip_args = grip_args,
      weighted_args = weighted_args,
      fr_niter = fr_niter,
      seed = level.seed
    )
    coords <- level.results[[i]]$coords
    level.trace.rows[[i]] <- data.frame(
      level = level.results[[i]]$level,
      inserted = length(level.results[[i]]$placed_vertices),
      mean_objective = NA_real_,
      max_grad_norm = NA_real_,
      all_converged = NA,
      stringsAsFactors = FALSE
    )
  }

  list(
    coords = coords,
    level_results = level.results,
    level_trace = do.call(rbind, level.trace.rows)
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
#' @param top_level_init Coarse-level initializer. `"geometric"` builds a
#'   spread `d+1`-vertex seed and inserts the remaining coarse vertices before
#'   pure-GMDS refinement. `"random"` keeps the legacy random restart family.
#' @param top_level_restarts Number of coarse-level restarts used by the
#'   pure-GMDS refinement stage.
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
                                           top_level_init = c("geometric", "random"),
                                           top_level_restarts = 8L,
                                           top_level_max_iter = 16L,
                                           top_level_engine = c("cpp", "r"),
                                           seed = 6L) {
  tie_mode <- match.arg(tie_mode)
  top_level_mode <- match.arg(top_level_mode)
  top_level_init <- match.arg(top_level_init)
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
  prepared$top_level_init <- top_level_init
  prepared$top_level_restarts <- as.integer(top_level_restarts)
  prepared$top_level_max_iter <- as.integer(top_level_max_iter)
  prepared$top_level_engine <- top_level_engine
  prepared$multiscale_mode <- "misf"
  prepared$insertion_anchor_policy <- "prev_level_spread"
  prepared$insertion_anchor_count <- grip.geodesic.misf.default.anchor.count(dim)
  prepared$insertion_anchor_weight_mode <- "inverse_graph_distance_sq"
  prepared$insertion_max_iter <- 64L
  prepared$insertion_mode <- "geodesic"
  prepared$insertion_layout_k <- 6L
  prepared$insertion_weighted_preset <- NULL
  prepared$insertion_grip_args <- list()
  prepared$insertion_weighted_args <- list()
  prepared$insertion_fr_niter <- 800L
  prepared$refinement_local_nbrs <- 8L
  prepared$refinement_landmark_count <- 4L
  prepared$refinement_pair_mode <- "sparse"
  prepared$refinement_anchor_weight <- 0.05
  prepared$refinement_anchor_weight_end <- 0.05
  prepared$refinement_continuation <- "constant"
  prepared$refinement_max_iter <- 8L
  prepared$refinement_engine <- top_level_engine
  prepared$final_polish_max_iter <- 8L
  prepared$final_polish_engine <- top_level_engine
  prepared$misf_seed <- seed
  prepared$top_level_fit <- NULL
  class(prepared) <- c("grip_misf_gmds_prepared", class(prepared))

  if (identical(top_level_mode, "solve")) {
    prepared$top_level_fit <- grip.geodesic.misf.solve.top.level(
      prepared = prepared,
      dim = dim,
      n_restarts = top_level_restarts,
      max_iter = top_level_max_iter,
      init = top_level_init,
      engine = top_level_engine,
      seed = seed
    )
  }
  prepared
}

#' Optimize an embedding with the MISF-based GMDS pipeline
#'
#' `grip.optimize.misf.geodesic.mds()` runs the experimental multiscale GMDS
#' pipeline built on top of GRIP's maximal independent set filtration (MISF).
#' It solves the coarsest MISF level with pure GMDS, inserts successive levels
#' either by geodesic anchor trilateration or by layout-based active-level warm
#' starts, refines each active level with sparse pure GMDS, and finishes with a
#' short full-graph pure-GMDS polish.
#'
#' The function accepts either a graph-first GMDS prepared object from
#' [grip.prepare.graph.geodesic.mds()], a MISF-GMDS prepared object from
#' [grip.prepare.misf.geodesic.mds()], or raw graph inputs.
#'
#' @param prepared Optional prepared object. This can be either a graph-first
#'   GMDS prepared object or a MISF-GMDS prepared object.
#' @param edges Two-column integer matrix of edges (1-based vertex ids) used
#'   when `prepared` is omitted.
#' @param n Number of vertices used when `prepared` is omitted.
#' @param adj_list Optional adjacency list used when `prepared` is omitted.
#' @param weight_list Optional edge-weight list parallel to `adj_list`.
#' @param edge_weights Optional positive edge-weight vector parallel to `edges`.
#' @param tie_mode Optional shortest-path aggregation mode used when a new MISF
#'   prepared object must be built.
#' @param num_init MISF top-level target size used when a new MISF prepared
#'   object must be built.
#' @param num_nbrs MISF neighborhood schedule parameter used when a new MISF
#'   prepared object must be built.
#' @param dim Optional target embedding dimension. If omitted, reuse the
#'   dimension stored in `prepared` when available.
#' @param top_level_init Initializer used for the coarse-level placement before
#'   pure-GMDS refinement.
#' @param top_level_restarts Number of coarse-level restarts used by the
#'   top-level pure-GMDS refinement.
#' @param top_level_max_iter Maximum number of top-level pure-GMDS iterations
#'   per restart.
#' @param top_level_engine Engine used by the top-level pure-GMDS solve.
#' @param insertion_anchor_policy Anchor-selection policy used during the
#'   insertion stage.
#' @param insertion_anchor_count Number of anchors used for each inserted
#'   vertex. If omitted, reuse the stored MISF-GMDS default.
#' @param insertion_anchor_weight_mode Anchor-weight schedule used during the
#'   insertion stage.
#' @param insertion_max_iter Maximum number of trilateration iterations per
#'   inserted vertex.
#' @param insertion_mode Lower-level placement mode. `"geodesic"` keeps the
#'   original anchor-trilateration insertion, while the other modes use a
#'   layout-based warm start on the active level before GMDS refinement.
#' @param insertion_layout_k Neighborhood size used to build the sparse active
#'   graph for layout-based insertion.
#' @param insertion_weighted_preset Optional weighted GRIP preset forwarded to
#'   `grip.layout.globalrep.weighted()` when `insertion_mode =
#'   "weighted_grip"`.
#' @param insertion_grip_args Optional named list of extra arguments forwarded
#'   to `grip.layout.globalrep()` when `insertion_mode = "grip"`.
#' @param insertion_weighted_args Optional named list of extra arguments
#'   forwarded to `grip.layout.globalrep.weighted()` when `insertion_mode =
#'   "weighted_grip"`.
#' @param insertion_fr_niter Number of FR iterations used when `insertion_mode
#'   = "fr"`.
#' @param refinement_local_nbrs Number of sparse local graph-neighbor pairs used
#'   during active-level refinement.
#' @param refinement_landmark_count Number of sparse landmark pairs used during
#'   active-level refinement.
#' @param refinement_pair_mode Sparse pair policy used during level refinement.
#' @param refinement_anchor_weight Initial anchor weight used to pin previously
#'   placed vertices during active-level refinement.
#' @param refinement_anchor_weight_end Final anchor weight used at the end of
#'   the level-refinement continuation schedule.
#' @param refinement_continuation Continuation schedule used during active-level
#'   refinement.
#' @param refinement_max_iter Maximum number of sparse pure-GMDS iterations per
#'   active level.
#' @param refinement_engine Engine used during active-level refinement.
#' @param final_polish_max_iter Maximum number of full-graph pure-GMDS polish
#'   iterations.
#' @param final_polish_engine Engine used during the final full-graph polish.
#' @param edge_length_epsilon Small non-negative stabilizer added inside
#'   embedded edge lengths during the refinement and final-polish stages.
#' @param n_threads Number of compiled-engine threads used by the refinement and
#'   final-polish stages.
#' @param return_trace If `TRUE`, include detailed per-stage traces.
#' @param return_frames If `TRUE`, retain intermediate coordinate frames for the
#'   multiscale stages.
#' @param seed Optional integer seed reused when a top-level restart solve must
#'   be computed.
#'
#' @return A list of class `"grip_misf_gmds_fit"` containing the final
#'   coordinates, the prepared MISF-GMDS object, stage summaries, optional
#'   detailed traces, and optional intermediate frames.
#' @export
grip.optimize.misf.geodesic.mds <- function(prepared = NULL,
                                            edges = NULL,
                                            n = NULL,
                                            adj_list = NULL,
                                            weight_list = NULL,
                                            edge_weights = NULL,
                                            tie_mode = NULL,
                                            num_init = 24L,
                                            num_nbrs = 20L,
                                            dim = NULL,
                                            top_level_init = NULL,
                                            top_level_restarts = NULL,
                                            top_level_max_iter = NULL,
                                            top_level_engine = NULL,
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
                                            refinement_local_nbrs = NULL,
                                            refinement_landmark_count = NULL,
                                            refinement_pair_mode = NULL,
                                            refinement_anchor_weight = NULL,
                                            refinement_anchor_weight_end = NULL,
                                            refinement_continuation = NULL,
                                            refinement_max_iter = NULL,
                                            refinement_engine = NULL,
                                            final_polish_max_iter = NULL,
                                            final_polish_engine = NULL,
                                            edge_length_epsilon = 1e-8,
                                            n_threads = 0L,
                                            return_trace = FALSE,
                                            return_frames = FALSE,
                                            seed = 6L) {
  if (!is.null(seed)) {
    seed <- grip.validate.count(seed, "seed")
  }

  dim.resolved <- if (is.null(dim) && !is.null(prepared) && inherits(prepared, "grip_misf_gmds_prepared")) {
    prepared$top_level_dim
  } else if (is.null(dim)) {
    2L
  } else {
    grip.validate.count(dim, "dim")
  }
  top.level.init <- if (is.null(top_level_init) && !is.null(prepared) && inherits(prepared, "grip_misf_gmds_prepared")) {
    if (!is.null(prepared$top_level_init)) prepared$top_level_init else "geometric"
  } else if (is.null(top_level_init)) {
    "geometric"
  } else {
    match.arg(top_level_init, c("geometric", "random"))
  }
  if (!(dim.resolved %in% c(2L, 3L))) {
    stop("dim must be 2 or 3")
  }

  top.level.restarts <- if (is.null(top_level_restarts) && !is.null(prepared) && inherits(prepared, "grip_misf_gmds_prepared")) {
    prepared$top_level_restarts
  } else if (is.null(top_level_restarts)) {
    8L
  } else {
    grip.validate.misf.count(top_level_restarts, "top_level_restarts", lower = 1L)
  }
  top.level.max.iter <- if (is.null(top_level_max_iter) && !is.null(prepared) && inherits(prepared, "grip_misf_gmds_prepared")) {
    prepared$top_level_max_iter
  } else if (is.null(top_level_max_iter)) {
    16L
  } else {
    grip.validate.scalar(top_level_max_iter, "top_level_max_iter", lower = 0)
    as.integer(round(top_level_max_iter))
  }
  top.level.engine <- if (is.null(top_level_engine) && !is.null(prepared) && inherits(prepared, "grip_misf_gmds_prepared")) {
    prepared$top_level_engine
  } else if (is.null(top_level_engine)) {
    "cpp"
  } else {
    match.arg(top_level_engine, c("cpp", "r"))
  }

  prepared <- grip.resolve.misf.geodesic.prepared(
    prepared = prepared,
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights,
    tie_mode = tie_mode,
    num_init = num_init,
    num_nbrs = num_nbrs,
    dim = dim.resolved,
    top_level_init = top.level.init,
    top_level_restarts = top.level.restarts,
    top_level_max_iter = top.level.max.iter,
    top_level_engine = top.level.engine,
    seed = seed
  )

  insertion.anchor.policy <- if (is.null(insertion_anchor_policy)) {
    prepared$insertion_anchor_policy
  } else {
    match.arg(insertion_anchor_policy, c("prev_level_first", "prev_level_distance_band", "prev_level_spread"))
  }
  insertion.anchor.count <- if (is.null(insertion_anchor_count)) {
    prepared$insertion_anchor_count
  } else {
    grip.validate.misf.count(insertion_anchor_count, "insertion_anchor_count", lower = 1L)
  }
  insertion.anchor.weight.mode <- if (is.null(insertion_anchor_weight_mode)) {
    prepared$insertion_anchor_weight_mode
  } else {
    match.arg(insertion_anchor_weight_mode, c("inverse_graph_distance_sq", "uniform"))
  }
  insertion.max.iter <- if (is.null(insertion_max_iter)) {
    prepared$insertion_max_iter
  } else {
    grip.validate.scalar(insertion_max_iter, "insertion_max_iter", lower = 0)
    as.integer(round(insertion_max_iter))
  }
  insertion.mode <- if (is.null(insertion_mode)) {
    prepared$insertion_mode
  } else {
    match.arg(insertion_mode, c("geodesic", "kk", "weighted_kk", "fr", "grip", "weighted_grip"))
  }
  insertion.layout.k <- if (is.null(insertion_layout_k)) {
    prepared$insertion_layout_k
  } else {
    grip.validate.count(insertion_layout_k, "insertion_layout_k")
  }
  insertion.weighted.preset <- if (is.null(insertion_weighted_preset)) {
    prepared$insertion_weighted_preset
  } else {
    as.character(insertion_weighted_preset)
  }
  insertion.grip.args <- if (is.null(insertion_grip_args)) {
    prepared$insertion_grip_args
  } else {
    if (!is.list(insertion_grip_args)) {
      stop("insertion_grip_args must be NULL or a list")
    }
    insertion_grip_args
  }
  insertion.weighted.args <- if (is.null(insertion_weighted_args)) {
    prepared$insertion_weighted_args
  } else {
    if (!is.list(insertion_weighted_args)) {
      stop("insertion_weighted_args must be NULL or a list")
    }
    insertion_weighted_args
  }
  insertion.fr.niter <- if (is.null(insertion_fr_niter)) {
    prepared$insertion_fr_niter
  } else {
    grip.validate.count(insertion_fr_niter, "insertion_fr_niter")
  }
  prepared$insertion_anchor_policy <- insertion.anchor.policy
  prepared$insertion_anchor_count <- insertion.anchor.count
  prepared$insertion_anchor_weight_mode <- insertion.anchor.weight.mode
  prepared$insertion_max_iter <- insertion.max.iter
  prepared$insertion_mode <- insertion.mode
  prepared$insertion_layout_k <- insertion.layout.k
  prepared$insertion_weighted_preset <- insertion.weighted.preset
  prepared$insertion_grip_args <- insertion.grip.args
  prepared$insertion_weighted_args <- insertion.weighted.args
  prepared$insertion_fr_niter <- insertion.fr.niter

  refinement.local.nbrs <- if (is.null(refinement_local_nbrs)) {
    prepared$refinement_local_nbrs
  } else {
    grip.validate.count(refinement_local_nbrs, "refinement_local_nbrs")
  }
  refinement.landmark.count <- if (is.null(refinement_landmark_count)) {
    prepared$refinement_landmark_count
  } else {
    grip.validate.count(refinement_landmark_count, "refinement_landmark_count")
  }
  refinement.pair.mode <- if (is.null(refinement_pair_mode)) {
    prepared$refinement_pair_mode
  } else {
    match.arg(refinement_pair_mode, c("sparse", "full"))
  }
  refinement.anchor.weight <- if (is.null(refinement_anchor_weight)) {
    prepared$refinement_anchor_weight
  } else {
    grip.validate.scalar(refinement_anchor_weight, "refinement_anchor_weight", lower = 0)
    as.double(refinement_anchor_weight)
  }
  refinement.anchor.weight.end <- if (is.null(refinement_anchor_weight_end)) {
    if (!is.null(prepared$refinement_anchor_weight_end)) prepared$refinement_anchor_weight_end else refinement.anchor.weight
  } else {
    grip.validate.scalar(refinement_anchor_weight_end, "refinement_anchor_weight_end", lower = 0)
    as.double(refinement_anchor_weight_end)
  }
  refinement.schedule <- if (is.null(refinement_continuation)) {
    prepared$refinement_continuation
  } else {
    match.arg(refinement_continuation, c("constant", "linear", "geometric"))
  }
  refinement.max.iter <- if (is.null(refinement_max_iter)) {
    prepared$refinement_max_iter
  } else {
    grip.validate.scalar(refinement_max_iter, "refinement_max_iter", lower = 0)
    as.integer(round(refinement_max_iter))
  }
  refinement.engine.resolved <- if (is.null(refinement_engine)) {
    prepared$refinement_engine
  } else {
    match.arg(refinement_engine, c("cpp", "r"))
  }

  final.polish.max.iter <- if (is.null(final_polish_max_iter)) {
    prepared$final_polish_max_iter
  } else {
    grip.validate.scalar(final_polish_max_iter, "final_polish_max_iter", lower = 0)
    as.integer(round(final_polish_max_iter))
  }
  final.polish.engine <- if (is.null(final_polish_engine)) {
    prepared$final_polish_engine
  } else {
    match.arg(final_polish_engine, c("cpp", "r"))
  }
  grip.validate.scalar(edge_length_epsilon, "edge_length_epsilon", lower = 0)
  grip.validate.scalar(n_threads, "n_threads", lower = 0)
  n_threads <- as.integer(round(n_threads))

  need.top.trace <- isTRUE(return_trace) || isTRUE(return_frames)
  need.top.solve <- is.null(prepared$top_level_fit)
  if (!need.top.solve) {
    need.top.solve <- ncol(prepared$top_level_fit$coords) != dim.resolved
    if (!need.top.solve && isTRUE(return_trace)) {
      need.top.solve <- is.null(prepared$top_level_fit$trace) || nrow(prepared$top_level_fit$trace) == 0L
    }
    if (!need.top.solve && isTRUE(return_frames)) {
      need.top.solve <- is.null(prepared$top_level_fit$frames) || length(prepared$top_level_fit$frames) <= 1L
    }
  }

  top.level.elapsed <- 0
  if (need.top.solve) {
    top.level.start <- proc.time()[["elapsed"]]
    prepared$top_level_fit <- grip.geodesic.misf.solve.top.level(
      prepared = prepared,
      dim = dim.resolved,
      n_restarts = top.level.restarts,
      max_iter = top.level.max.iter,
      init = top.level.init,
      engine = top.level.engine,
      edge_length_epsilon = edge_length_epsilon,
      n_threads = n_threads,
      return_trace = need.top.trace,
      seed = seed
    )
    top.level.elapsed <- proc.time()[["elapsed"]] - top.level.start
  }
  top.level.fit <- prepared$top_level_fit
  coords <- top.level.fit$coords_full

  insertion.elapsed <- 0
  if (any(!is.finite(coords))) {
    insertion.start <- proc.time()[["elapsed"]]
    insertion <- if (identical(insertion.mode, "geodesic")) {
      grip.geodesic.misf.insert.all.levels(
        prepared = prepared,
        coords = coords,
        anchor_policy = insertion.anchor.policy,
        anchor_count = insertion.anchor.count,
        anchor_weight_mode = insertion.anchor.weight.mode,
        max_iter = insertion.max.iter
      )
    } else {
      grip.geodesic.misf.insert.all.levels.with.layout(
        prepared = prepared,
        coords = coords,
        method = insertion.mode,
        layout_k = insertion.layout.k,
        weighted_preset = insertion.weighted.preset,
        grip_args = insertion.grip.args,
        weighted_args = insertion.weighted.args,
        fr_niter = insertion.fr.niter,
        seed = seed
      )
    }
    insertion.elapsed <- proc.time()[["elapsed"]] - insertion.start
  } else {
    insertion <- list(
      coords = coords,
      level_results = list(),
      level_trace = data.frame()
    )
  }

  refinement.start <- proc.time()[["elapsed"]]
  refinement <- grip.geodesic.misf.refine.all.levels(
    prepared = prepared,
    coords = insertion$coords,
    local_nbrs = refinement.local.nbrs,
    landmark_count = refinement.landmark.count,
    pair_mode = refinement.pair.mode,
    anchor_weight = refinement.anchor.weight,
    anchor_weight_end = refinement.anchor.weight.end,
    continuation = refinement.schedule,
    max_iter = refinement.max.iter,
    engine = refinement.engine.resolved,
    edge_length_epsilon = edge_length_epsilon,
    n_threads = n_threads,
    return_trace = isTRUE(return_trace) || isTRUE(return_frames)
  )
  refinement.elapsed <- proc.time()[["elapsed"]] - refinement.start

  final.polish.start <- proc.time()[["elapsed"]]
  final.polish <- grip.geodesic.misf.final.polish(
    prepared = prepared,
    coords = refinement$coords,
    max_iter = final.polish.max.iter,
    engine = final.polish.engine,
    edge_length_epsilon = edge_length_epsilon,
    n_threads = n_threads,
    return_trace = isTRUE(return_trace) || isTRUE(return_frames)
  )
  final.polish.elapsed <- proc.time()[["elapsed"]] - final.polish.start

  top.level.frames <- if (isTRUE(return_frames)) {
    grip.geodesic.misf.expand.top.level.frames(top.level.fit, prepared$n)
  } else {
    NULL
  }
  insertion.level.frames <- if (isTRUE(return_frames)) {
    grip.geodesic.misf.collect.level.frames(insertion$level_results)
  } else {
    NULL
  }
  refinement.level.frames <- if (isTRUE(return_frames)) {
    grip.geodesic.misf.collect.level.frames(refinement$level_results)
  } else {
    NULL
  }
  stage.trace <- grip.geodesic.misf.build.stage.trace(
    prepared = prepared,
    top_level_fit = top.level.fit,
    top_level_elapsed = top.level.elapsed,
    top_level_frames = top.level.frames,
    insertion = insertion,
    insertion_elapsed = insertion.elapsed,
    insertion_frames = insertion.level.frames,
    refinement = refinement,
    refinement_elapsed = refinement.elapsed,
    refinement_frames = refinement.level.frames,
    final_polish = final.polish,
    final_polish_elapsed = final.polish.elapsed
  )

  trace.detail <- if (isTRUE(return_trace)) {
    list(
      top_level_trace = top.level.fit$trace,
      top_restart_summary = top.level.fit$restart_summary,
      insertion_level_trace = insertion$level_trace,
      insertion_vertex_trace = grip.geodesic.misf.collect.insertion.vertex.trace(insertion$level_results),
      refinement_level_trace = refinement$level_trace,
      final_polish_trace = final.polish$trace
    )
  } else {
    NULL
  }
  frames <- if (isTRUE(return_frames)) {
    list(
      top_level = top.level.frames,
      after_top_level = top.level.fit$coords_full,
      insertion_levels = insertion.level.frames,
      after_insertion = insertion$coords,
      refinement_levels = refinement.level.frames,
      after_refinement = refinement$coords,
      final_polish = final.polish$frames,
      final = final.polish$coords
    )
  } else {
    NULL
  }

  fit <- list(
    coords = final.polish$coords,
    prepared = prepared,
    top_level_fit = top.level.fit,
    insertion = insertion,
    refinement = refinement,
    final_polish = final.polish,
    stage_trace = stage.trace,
    trace = trace.detail,
    frames = frames,
    timing = list(
      top_level = as.double(top.level.elapsed),
      insertion = as.double(insertion.elapsed),
      refinement = as.double(refinement.elapsed),
      final_polish = as.double(final.polish.elapsed),
      total = as.double(top.level.elapsed + insertion.elapsed + refinement.elapsed + final.polish.elapsed)
    )
  )
  class(fit) <- c("grip_misf_gmds_fit", "list")
  fit$score <- grip.score.misf.geodesic.mds(fit = fit, return_trace = isTRUE(return_trace))
  fit
}

#' Score a MISF-based GMDS fit
#'
#' `grip.score.misf.geodesic.mds()` summarizes the final full-graph GMDS fit
#' together with multiscale stage metadata from
#' [grip.optimize.misf.geodesic.mds()].
#'
#' @param fit Optional fit from [grip.optimize.misf.geodesic.mds()].
#' @param coords Optional final coordinate matrix used when `fit` is omitted.
#' @param prepared Optional MISF-GMDS prepared object used when `fit` is
#'   omitted.
#' @param edge_length_epsilon Small non-negative stabilizer added inside each
#'   embedded edge length when rescoring `coords`.
#' @param return_trace If `TRUE`, attach the multiscale trace tables as list
#'   columns.
#'
#' @return A one-row data frame with final GMDS metrics and MISF-stage summary
#'   fields. When `return_trace = TRUE`, the trace tables are attached as list
#'   columns.
#' @export
grip.score.misf.geodesic.mds <- function(fit = NULL,
                                         coords = NULL,
                                         prepared = NULL,
                                         edge_length_epsilon = 1e-8,
                                         return_trace = FALSE) {
  top.level.fit <- NULL
  insertion <- NULL
  refinement <- NULL
  final.polish <- NULL
  stage.trace <- data.frame()
  timing <- list(
    top_level = NA_real_,
    insertion = NA_real_,
    refinement = NA_real_,
    final_polish = NA_real_,
    total = NA_real_
  )

  if (!is.null(fit)) {
    fit <- grip.validate.misf.geodesic.fit(fit)
    if (is.null(prepared)) {
      prepared <- fit$prepared
    }
    if (is.null(coords)) {
      coords <- fit$coords
    }
    top.level.fit <- fit$top_level_fit
    insertion <- fit$insertion
    refinement <- fit$refinement
    final.polish <- fit$final_polish
    if (!is.null(fit$stage_trace)) {
      stage.trace <- fit$stage_trace
    }
    if (!is.null(fit$timing)) {
      timing <- fit$timing
    }
  }

  coords <- grip.validate.coords(coords)
  prepared <- grip.validate.misf.geodesic.prepared(prepared, coords = coords)
  grip.validate.scalar(edge_length_epsilon, "edge_length_epsilon", lower = 0)

  final.score <- grip.score.geodesic.mds(
    coords = coords,
    prepared = prepared,
    edge_length_epsilon = edge_length_epsilon
  )
  names(final.score) <- paste0("final.", names(final.score))

  insertion.level.count <- if (!is.null(insertion) && !is.null(insertion$level_trace)) {
    nrow(insertion$level_trace)
  } else if (prepared$top_level_level > 0L) {
    prepared$top_level_level
  } else {
    0L
  }
  inserted.vertex.count <- if (!is.null(insertion) && !is.null(insertion$level_trace) && nrow(insertion$level_trace) > 0L) {
    sum(insertion$level_trace$inserted)
  } else {
    max(0L, prepared$n - length(prepared$top_level_vertices))
  }
  refinement.level.count <- if (!is.null(refinement) && !is.null(refinement$level_trace)) {
    nrow(refinement$level_trace)
  } else {
    0L
  }

  out <- data.frame(
    multiscale.mode = "misf",
    n = prepared$n,
    dim = ncol(coords),
    top.level = prepared$top_level_level,
    top.level.n = length(prepared$top_level_vertices),
    top.level.energy = if (!is.null(top.level.fit)) top.level.fit$score$gmds.energy[[1L]] else NA_real_,
    top.level.stress = if (!is.null(top.level.fit)) top.level.fit$score$gmds.stress[[1L]] else NA_real_,
    insertion.level.count = as.integer(insertion.level.count),
    inserted.vertex.count = as.integer(inserted.vertex.count),
    refinement.level.count = as.integer(refinement.level.count),
    final.polish.trace.rows = if (!is.null(final.polish) && !is.null(final.polish$trace)) nrow(final.polish$trace) else NA_integer_,
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
    out$top.restart.summary <- list(if (!is.null(top.level.fit)) top.level.fit$restart_summary else data.frame())
    out$insertion.level.trace <- list(if (!is.null(insertion)) insertion$level_trace else data.frame())
    out$insertion.vertex.trace <- list(
      if (!is.null(insertion)) grip.geodesic.misf.collect.insertion.vertex.trace(insertion$level_results) else data.frame()
    )
    out$refinement.level.trace <- list(if (!is.null(refinement)) refinement$level_trace else data.frame())
    out$final.polish.trace <- list(if (!is.null(final.polish)) final.polish$trace else data.frame())
  }

  out
}
