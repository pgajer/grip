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
    stop("prepared must be an object from prepare.misf.geodesic.kk()")
  }
  if (!is.null(coords) && nrow(coords) != prepared$n) {
    stop("nrow(coords) must match the graph size stored in prepared")
  }
  prepared
}

grip.validate.misf.geodesic.kk.fit <- function(fit) {
  if (!inherits(fit, "grip_misf_gkk_fit")) {
    stop("fit must be an object from misf.geodesic.kk()")
  }
  fit
}

grip.geodesic.misf.kk.resolve.pair.mode <- function(pair.mode = c("auto", "full", "landmark"),
                                                     active.n,
                                                     full.limit = 512L) {
  pair.mode <- match.arg(pair.mode)
  active.n <- grip.validate.count(active.n, "active.n")
  full.limit <- grip.validate.misf.count(full.limit, "full.limit", lower = 1L)
  effective <- if (identical(pair.mode, "auto")) {
    if (active.n <= full.limit) "full" else "landmark"
  } else {
    pair.mode
  }
  list(
    requested = pair.mode,
    effective = effective,
    full_limit = full.limit,
    active_n = active.n
  )
}

grip.resolve.misf.geodesic.kk.prepared <- function(prepared = NULL,
                                                   edges = NULL,
                                                   n = NULL,
                                                   adj.list = NULL,
                                                   weight.list = NULL,
                                                   edge.weights = NULL,
                                                   tie.mode = NULL,
                                                   num.init = 24L,
                                                   num.nbrs = 20L,
                                                   dim = NULL,
                                                   top.level.pair.mode = c("auto", "full", "landmark"),
                                                   top.level.full.limit = 512L,
                                                   top.level.local.nbrs = 20L,
                                                   top.level.landmark.count = 8L,
                                                   top.level.restarts = 8L,
                                                   top.level.max.iter = 16L,
                                                   top.level.init = c("geometric", "cmdscale", "random"),
                                                   seed = 6L) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights, prepared)
  top.level.pair.mode <- match.arg(top.level.pair.mode)
  top.level.init <- match.arg(top.level.init)
  if (is.null(prepared)) {
    resolved.dim <- if (is.null(dim)) 2L else grip.validate.count(dim, "dim")
    resolved.tie.mode <- if (is.null(tie.mode)) "average" else {
      match.arg(tie.mode, c("single", "average"))
    }
    return(prepare.misf.geodesic.kk(
      edges = edges,
      n = n,
      adj.list = adj.list,
      weight.list = weight.list,
      edge.weights = edge.weights,
      tie.mode = resolved.tie.mode,
      num.init = num.init,
      num.nbrs = num.nbrs,
      dim = resolved.dim,
      top.level.mode = "skip",
      top.level.pair.mode = top.level.pair.mode,
      top.level.full.limit = top.level.full.limit,
      top.level.local.nbrs = top.level.local.nbrs,
      top.level.landmark.count = top.level.landmark.count,
      top.level.restarts = top.level.restarts,
      top.level.max.iter = top.level.max.iter,
      top.level.init = top.level.init,
      seed = seed
    ))
  }
  if (inherits(prepared, "grip_misf_gkk_prepared")) {
    return(grip.validate.misf.geodesic.kk.prepared(prepared))
  }
  if (!inherits(prepared, "grip_gkk_prepared")) {
    stop(
      "prepared must be NULL, an object from prepare.geodesic.kk(), ",
      "or an object from prepare.misf.geodesic.kk()"
    )
  }
  resolved.dim <- if (is.null(dim)) 2L else grip.validate.count(dim, "dim")
  resolved.tie.mode <- if (is.null(tie.mode)) {
    if (!is.null(prepared$tie_mode)) prepared$tie_mode else "average"
  } else {
    match.arg(tie.mode, c("single", "average"))
  }
  prepare.misf.geodesic.kk(
    n = prepared$n,
    adj.list = prepared$adj_list,
    weight.list = prepared$weight_list,
    tie.mode = resolved.tie.mode,
    num.init = num.init,
    num.nbrs = num.nbrs,
    dim = resolved.dim,
    top.level.mode = "skip",
    top.level.pair.mode = top.level.pair.mode,
    top.level.full.limit = top.level.full.limit,
    top.level.local.nbrs = top.level.local.nbrs,
    top.level.landmark.count = top.level.landmark.count,
    top.level.restarts = top.level.restarts,
    top.level.max.iter = top.level.max.iter,
    top.level.init = top.level.init,
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

grip.geodesic.misf.kk.score.kind <- function(prepared) {
  pair.mode <- if (!is.null(prepared$pair_mode)) prepared$pair_mode else "all_pairs"
  if (pair.mode %in% c("landmark_sparse", "misf_sparse")) "landmark" else "full"
}

grip.geodesic.misf.kk.add.prepared.class <- function(prepared,
                                                     effective.pair.mode = NULL) {
  if (is.null(effective.pair.mode)) {
    effective.pair.mode <- grip.geodesic.misf.kk.score.kind(prepared)
  }
  extra.class <- if (identical(effective.pair.mode, "landmark")) {
    "grip_lgkk_prepared"
  } else {
    "grip_gkk_prepared"
  }
  class(prepared) <- unique(c(extra.class, class(prepared), "grip_geodesic_kk_prepared", "list"))
  prepared
}

grip.geodesic.misf.kk.score.prepared <- function(coords,
                                                 prepared,
                                                 stiffness = 1.0,
                                                 distance.floor = 1e-8,
                                                 edge.length.epsilon = 1e-8) {
  effective.pair.mode <- grip.geodesic.misf.kk.score.kind(prepared)
  prepared <- grip.geodesic.misf.kk.add.prepared.class(prepared, effective.pair.mode)
  if (identical(effective.pair.mode, "landmark")) {
    score.landmark.geodesic.kk(
      coords = coords,
      prepared = prepared,
      stiffness = stiffness,
      distance.floor = distance.floor,
      edge.length.epsilon = edge.length.epsilon
    )
  } else {
    score.geodesic.kk(
      coords = coords,
      prepared = prepared,
      stiffness = stiffness,
      distance.floor = distance.floor,
      edge.length.epsilon = edge.length.epsilon,
      scale.mode = "profiled"
    )
  }
}

grip.geodesic.misf.kk.resolve.score.prepared <- function(prepared,
                                                         pair.mode = c("full", "landmark", "auto"),
                                                         full.limit = 2048L,
                                                         local.nbrs = 20L,
                                                         landmark.count = 8L) {
  prepared <- grip.validate.misf.geodesic.kk.prepared(prepared)
  pair.mode <- match.arg(pair.mode)
  resolution <- grip.geodesic.misf.kk.resolve.pair.mode(
    pair.mode = pair.mode,
    active.n = prepared$n,
    full.limit = full.limit
  )
  score.prepared <- if (identical(resolution$effective, "full")) {
    grip.geodesic.misf.kk.add.prepared.class(prepared, "full")
  } else {
    grip.geodesic.misf.kk.add.prepared.class(
      prepare.landmark.geodesic.kk(
        n = prepared$n,
        adj.list = prepared$adj_list,
        weight.list = prepared$weight_list,
        local.nbrs = local.nbrs,
        landmark.count = landmark.count
      ),
      "landmark"
    )
  }
  list(
    prepared = score.prepared,
    pair_resolution = resolution
  )
}

grip.geodesic.misf.kk.evaluate.state <- function(coords,
                                                 prepared,
                                                 stiffness = 1.0,
                                                 distance.floor = 1e-8,
                                                 edge.length.epsilon = 1e-8,
                                                 scale.mode = c("profiled", "fixed"),
                                                 scale.L0 = NULL,
                                                 anchor.coords = NULL,
                                                 anchor.weight = 0,
                                                 anchor.vertex.weight = NULL) {
  scale.mode <- match.arg(scale.mode)
  anchor.coords <- if (is.null(anchor.coords)) {
    NULL
  } else {
    grip.geodesic.mds.resolve.anchor(
      anchor.mode = "user",
      coords = coords,
      prepared = prepared,
      anchor.coords = anchor.coords,
      recenter = FALSE
    )
  }
  kk.state <- grip.geodesic.kk.evaluate.state(
    coords = coords,
    prepared = prepared,
    stiffness = stiffness,
    distance.floor = distance.floor,
    edge.length.epsilon = edge.length.epsilon,
    scale.mode = scale.mode,
    scale.L0 = scale.L0
  )
  anchor.stats <- grip.geodesic.mds.anchor.stats(
    coords = coords,
    anchor.coords = anchor.coords,
    anchor.weight = anchor.weight,
    anchor.vertex.weight = anchor.vertex.weight
  )
  kk.state$kk_energy <- kk.state$energy
  kk.state$anchor_energy <- anchor.stats$energy
  kk.state$anchor_raw_penalty <- anchor.stats$raw_penalty
  kk.state$anchor_weight <- anchor.stats$anchor_weight
  kk.state$energy <- kk.state$kk_energy + anchor.stats$energy
  kk.state$gradient <- kk.state$gradient + anchor.stats$gradient
  kk.state$gradient_norm <- sqrt(sum(kk.state$gradient^2))
  kk.state
}

grip.geodesic.misf.kk.optimize.prepared <- function(coords,
                                                    prepared,
                                                    max.iter = 8L,
                                                    stiffness = 1.0,
                                                    distance.floor = 1e-8,
                                                    edge.length.epsilon = 1e-8,
                                                    initial.step = 1.0,
                                                    step.shrink = 0.5,
                                                    armijo.factor = 1e-4,
                                                    grad.tol = 1e-8,
                                                    min.step = 1e-8,
                                                    anchor.coords = NULL,
                                                    anchor.weight = 0,
                                                    anchor.weight.end = anchor.weight,
                                                    anchor.vertex.weight = NULL,
                                                    continuation = c("constant", "linear", "geometric"),
                                                    recenter = TRUE,
                                                    return.trace = FALSE) {
  continuation <- match.arg(continuation)
  effective.pair.mode <- grip.geodesic.misf.kk.score.kind(prepared)
  prepared <- grip.geodesic.misf.kk.add.prepared.class(prepared, effective.pair.mode)
  grip.validate.scalar(anchor.weight, "anchor.weight", lower = 0)
  grip.validate.scalar(anchor.weight.end, "anchor.weight.end", lower = 0)
  use.anchor <- !is.null(anchor.coords) &&
    is.finite(anchor.weight) &&
    is.finite(anchor.weight.end) &&
    max(anchor.weight, anchor.weight.end) > 0
  if (!use.anchor) {
    anchor.weight <- 0
    anchor.weight.end <- 0
  }

  if (!use.anchor) {
    fit <- if (identical(effective.pair.mode, "landmark")) {
      landmark.geodesic.kk(
        coords = coords,
        prepared = prepared,
        max.iter = max.iter,
        stiffness = stiffness,
        distance.floor = distance.floor,
        edge.length.epsilon = edge.length.epsilon,
        initial.step = initial.step,
        step.shrink = step.shrink,
        armijo.factor = armijo.factor,
        grad.tol = grad.tol,
        min.step = min.step,
        recenter = recenter,
        return.trace = return.trace
      )
    } else {
      geodesic.kk(
        coords = coords,
        prepared = prepared,
        max.iter = max.iter,
        stiffness = stiffness,
        distance.floor = distance.floor,
        edge.length.epsilon = edge.length.epsilon,
        initial.step = initial.step,
        step.shrink = step.shrink,
        armijo.factor = armijo.factor,
        grad.tol = grad.tol,
        min.step = min.step,
        recenter = recenter,
        return.trace = return.trace,
        scale.mode = "profiled"
      )
    }
    fit$effective_pair_mode <- effective.pair.mode
    return(fit)
  }

  anchor.coords <- grip.geodesic.mds.resolve.anchor(
    anchor.mode = "user",
    coords = coords,
    prepared = prepared,
    anchor.coords = anchor.coords,
    recenter = FALSE
  )
  anchor.vertex.weight <- grip.geodesic.mds.resolve.anchor.vertex.weight(
    anchor.vertex.weight = anchor.vertex.weight,
    coords = coords
  )
  anchor.schedule <- grip.geodesic.mds.anchor.schedule(
    max.iter = max.iter,
    anchor.weight = anchor.weight,
    anchor.weight.end = anchor.weight.end,
    continuation = continuation
  )
  fixed.scale.L0 <- NULL
  scale.mode <- if (identical(effective.pair.mode, "landmark")) "fixed" else "profiled"
  if (identical(scale.mode, "fixed") && length(prepared$pair_graph_distance) > 0L) {
    fixed.scale.L0 <- grip.geodesic.kk.fit.scale(
      path.lengths = grip.geodesic.kk.path.lengths(
        coords = coords,
        prepared = prepared,
        edge.length.epsilon = edge.length.epsilon
      ),
      graph.distances = prepared$pair_graph_distance,
      stiffness = stiffness,
      distance.floor = distance.floor
    )
    if (!is.finite(fixed.scale.L0)) {
      stop("failed to fit an initial LGKK scale")
    }
  }

  if (nrow(coords) <= 1L || length(prepared$pair_graph_distance) == 0L || max.iter == 0L) {
    score <- grip.geodesic.misf.kk.score.prepared(
      coords = coords,
      prepared = prepared,
      stiffness = stiffness,
      distance.floor = distance.floor,
      edge.length.epsilon = edge.length.epsilon
    )
    out <- list(
      coords = coords,
      trace = data.frame(),
      frames = list(coords),
      prepared = prepared,
      score = score,
      anchor_coords = anchor.coords,
      anchor_schedule = anchor.schedule,
      anchor_vertex_weight = anchor.vertex.weight,
      final_anchor_weight = anchor.schedule[[1L]]
    )
    out$effective_pair_mode <- effective.pair.mode
    return(out)
  }

  current <- coords
  trace.rows <- vector("list", max.iter + 1L)
  accepted.frames <- list(current)
  state <- grip.geodesic.misf.kk.evaluate.state(
    coords = current,
    prepared = prepared,
    stiffness = stiffness,
    distance.floor = distance.floor,
    edge.length.epsilon = edge.length.epsilon,
    scale.mode = scale.mode,
    scale.L0 = fixed.scale.L0,
    anchor.coords = anchor.coords,
    anchor.weight = anchor.schedule[[1L]],
    anchor.vertex.weight = anchor.vertex.weight
  )
  trace.rows[[1L]] <- data.frame(
    iteration = 0L,
    energy = state$energy,
    kk_energy = state$kk_energy,
    anchor_energy = state$anchor_energy,
    gradient_norm = state$gradient_norm,
    step = NA_real_,
    accepted = TRUE,
    scale.L0 = state$scale.L0,
    anchor_weight = anchor.schedule[[1L]],
    stringsAsFactors = FALSE
  )
  used <- 1L

  for (iter in seq_len(max.iter)) {
    iter.anchor.weight <- anchor.schedule[[iter + 1L]]
    state <- grip.geodesic.misf.kk.evaluate.state(
      coords = current,
      prepared = prepared,
      stiffness = stiffness,
      distance.floor = distance.floor,
      edge.length.epsilon = edge.length.epsilon,
      scale.mode = scale.mode,
      scale.L0 = fixed.scale.L0,
      anchor.coords = anchor.coords,
      anchor.weight = iter.anchor.weight,
      anchor.vertex.weight = anchor.vertex.weight
    )
    if (!is.finite(state$gradient_norm) || state$gradient_norm <= grad.tol) {
      break
    }
    step <- as.double(initial.step)
    accepted <- FALSE
    candidate <- current
    candidate.state <- state

    while (is.finite(step) && step >= min.step) {
      proposal <- current - step * state$gradient
      if (isTRUE(recenter)) {
        proposal <- sweep(proposal, 2L, colMeans(proposal), "-", check.margin = FALSE)
      }
      proposal.state <- grip.geodesic.misf.kk.evaluate.state(
        coords = proposal,
        prepared = prepared,
        stiffness = stiffness,
        distance.floor = distance.floor,
        edge.length.epsilon = edge.length.epsilon,
        scale.mode = scale.mode,
        scale.L0 = fixed.scale.L0,
        anchor.coords = anchor.coords,
        anchor.weight = iter.anchor.weight,
        anchor.vertex.weight = anchor.vertex.weight
      )
      target.energy <- state$energy - armijo.factor * step * state$gradient_norm^2
      if (is.finite(proposal.state$energy) && proposal.state$energy <= target.energy) {
        candidate <- proposal
        candidate.state <- proposal.state
        accepted <- TRUE
        break
      }
      step <- step * step.shrink
    }

    used <- used + 1L
    trace.rows[[used]] <- data.frame(
      iteration = iter,
      energy = if (accepted) candidate.state$energy else state$energy,
      kk_energy = if (accepted) candidate.state$kk_energy else state$kk_energy,
      anchor_energy = if (accepted) candidate.state$anchor_energy else state$anchor_energy,
      gradient_norm = if (accepted) candidate.state$gradient_norm else state$gradient_norm,
      step = if (accepted) step else NA_real_,
      accepted = accepted,
      scale.L0 = if (accepted) candidate.state$scale.L0 else state$scale.L0,
      anchor_weight = iter.anchor.weight,
      stringsAsFactors = FALSE
    )

    if (!accepted) {
      break
    }

    current <- candidate
    state <- candidate.state
    accepted.frames[[length(accepted.frames) + 1L]] <- current
  }

  trace.df <- do.call(rbind, trace.rows[seq_len(used)])
  if (!isTRUE(return.trace)) {
    trace.df <- trace.df[, c("iteration", "energy", "kk_energy", "anchor_energy", "gradient_norm", "step", "accepted", "scale.L0", "anchor_weight"), drop = FALSE]
    accepted.frames <- list(current)
  }
  score <- grip.geodesic.misf.kk.score.prepared(
    coords = current,
    prepared = prepared,
    stiffness = stiffness,
    distance.floor = distance.floor,
    edge.length.epsilon = edge.length.epsilon
  )
  fit <- list()
  fit$coords <- current
  fit$trace <- trace.df
  fit$frames <- accepted.frames
  fit$prepared <- prepared
  fit$score <- score
  fit$anchor_coords <- anchor.coords
  fit$anchor_schedule <- anchor.schedule
  fit$anchor_vertex_weight <- anchor.vertex.weight
  fit$final_anchor_weight <- trace.df$anchor_weight[[nrow(trace.df)]]
  fit$effective_pair_mode <- effective.pair.mode
  fit
}

grip.geodesic.misf.kk.init.coords <- function(prepared,
                                              dim = 2L,
                                              init = c("cmdscale", "random"),
                                              restart = 1L,
                                              seed = NULL) {
  prepared <- grip.validate.geodesic.mds.prepared(prepared)
  init <- match.arg(init)
  restart <- grip.validate.misf.count(restart, "restart", lower = 1L)
  dim <- grip.validate.count(dim, "dim")
  if (!(dim %in% c(2L, 3L))) {
    stop("dim must be 2 or 3")
  }
  if (!is.null(seed)) {
    set.seed(as.integer(seed + restart - 1L))
  }

  coords <- if (identical(init, "cmdscale")) {
    grip.geodesic.mds.cmdscale.init(prepared, dim)
  } else {
    matrix(stats::rnorm(prepared$n * dim), nrow = prepared$n, ncol = dim)
  }
  coords <- as.matrix(coords)
  storage.mode(coords) <- "double"
  if (identical(init, "cmdscale") && restart > 1L) {
    spread <- stats::sd(as.double(coords))
    if (!is.finite(spread) || spread <= 0) {
      spread <- 1.0
    }
    coords <- coords + matrix(
      stats::rnorm(length(coords), sd = 0.05 * spread),
      nrow = nrow(coords),
      ncol = ncol(coords)
    )
  }
  if (nrow(coords) > 0L) {
    coords <- sweep(coords, 2L, colMeans(coords), FUN = "-", check.margin = FALSE)
  }
  coords
}

grip.geodesic.misf.kk.build.restart.row <- function(restart,
                                                    seed,
                                                    init.score,
                                                    fit) {
  data.frame(
    restart = as.integer(restart),
    seed = if (is.null(seed)) NA_integer_ else as.integer(seed),
    pair.mode = if (!is.null(fit$effective_pair_mode)) fit$effective_pair_mode else NA_character_,
    initial.energy = grip.geodesic.misf.kk.extract.score.metric(
      init.score,
      c("gkk.energy", "lgkk.energy")
    ),
    initial.weighted.rel.rmse = grip.geodesic.misf.kk.extract.score.metric(
      init.score,
      c("gkk.weighted.rel.rmse", "lgkk.weighted.rel.rmse")
    ),
    final.energy = grip.geodesic.misf.kk.extract.score.metric(
      fit$score,
      c("gkk.energy", "lgkk.energy")
    ),
    final.weighted.rel.rmse = grip.geodesic.misf.kk.extract.score.metric(
      fit$score,
      c("gkk.weighted.rel.rmse", "lgkk.weighted.rel.rmse")
    ),
    improved = grip.geodesic.misf.kk.extract.score.metric(
      fit$score,
      c("gkk.energy", "lgkk.energy")
    ) <= grip.geodesic.misf.kk.extract.score.metric(
      init.score,
      c("gkk.energy", "lgkk.energy")
    ) + 1e-12,
    trace.rows = if (is.null(fit$trace)) 0L else nrow(fit$trace),
    stringsAsFactors = FALSE
  )
}

grip.geodesic.misf.kk.resolve.top.level.prepared <- function(prepared,
                                                             pair.mode = c("auto", "full", "landmark"),
                                                             full.limit = 512L,
                                                             local.nbrs = 20L,
                                                             landmark.count = 8L) {
  prepared <- grip.validate.misf.geodesic.kk.prepared(prepared)
  pair.mode <- match.arg(pair.mode)
  resolution <- grip.geodesic.misf.kk.resolve.pair.mode(
    pair.mode = pair.mode,
    active.n = prepared$top_level_graph$n,
    full.limit = full.limit
  )
  coarse.prepared <- if (identical(resolution$effective, "full")) {
    if (identical(local.nbrs, prepared$top_level_local_nbrs) &&
        identical(landmark.count, prepared$top_level_landmark_count) &&
        !is.null(prepared$top_level_prepared_full)) {
      prepared$top_level_prepared_full
    } else {
      prepare.geodesic.kk(
        edges = prepared$top_level_graph$edges,
        n = prepared$top_level_graph$n,
        edge.weights = prepared$top_level_graph$edge_weights,
        tie.mode = prepared$tie_mode
      )
    }
  } else {
    if (identical(local.nbrs, prepared$top_level_local_nbrs) &&
        identical(landmark.count, prepared$top_level_landmark_count) &&
        !is.null(prepared$top_level_prepared_sparse)) {
      prepared$top_level_prepared_sparse
    } else {
      prepare.landmark.geodesic.kk(
        edges = prepared$top_level_graph$edges,
        n = prepared$top_level_graph$n,
        edge.weights = prepared$top_level_graph$edge_weights,
        local.nbrs = local.nbrs,
        landmark.count = landmark.count
      )
    }
  }
  coarse.prepared <- grip.geodesic.misf.kk.add.prepared.class(coarse.prepared, resolution$effective)
  list(
    prepared = coarse.prepared,
    pair_resolution = resolution,
    vertex_ids = prepared$top_level_vertices
  )
}

grip.geodesic.misf.kk.solve.top.level <- function(prepared,
                                                  dim = 2L,
                                                  pair.mode = c("auto", "full", "landmark"),
                                                  full.limit = 512L,
                                                  local.nbrs = 20L,
                                                  landmark.count = 8L,
                                                  n.restarts = 8L,
                                                  max.iter = 16L,
                                                  init = c("geometric", "cmdscale", "random"),
                                                  stiffness = 1.0,
                                                  distance.floor = 1e-8,
                                                  edge.length.epsilon = 1e-8,
                                                  initial.step = 1.0,
                                                  step.shrink = 0.5,
                                                  armijo.factor = 1e-4,
                                                  grad.tol = 1e-8,
                                                  min.step = 1e-8,
                                                  recenter = TRUE,
                                                  return.trace = FALSE,
                                                  seed = 6L) {
  prepared <- grip.validate.misf.geodesic.kk.prepared(prepared)
  pair.mode <- match.arg(pair.mode)
  init <- match.arg(init)
  dim <- grip.validate.count(dim, "dim")
  if (!(dim %in% c(2L, 3L))) {
    stop("dim must be 2 or 3")
  }
  n.restarts <- grip.validate.misf.count(n.restarts, "n.restarts", lower = 1L)
  grip.validate.scalar(max.iter, "max.iter", lower = 0)
  max.iter <- as.integer(round(max.iter))

  resolved <- grip.geodesic.misf.kk.resolve.top.level.prepared(
    prepared = prepared,
    pair.mode = pair.mode,
    full.limit = full.limit,
    local.nbrs = local.nbrs,
    landmark.count = landmark.count
  )
  coarse.prepared <- resolved$prepared
  vertex.ids <- resolved$vertex_ids

  if (coarse.prepared$n <= 1L) {
    coords <- matrix(0, nrow = coarse.prepared$n, ncol = dim)
    score <- grip.geodesic.misf.kk.score.prepared(
      coords = coords,
      prepared = coarse.prepared,
      stiffness = stiffness,
      distance.floor = distance.floor,
      edge.length.epsilon = edge.length.epsilon
    )
    return(list(
      coords = coords,
      trace = data.frame(),
      frames = list(coords),
      prepared = coarse.prepared,
      score = score,
      restart_summary = data.frame(
        restart = 1L,
        seed = if (is.null(seed)) NA_integer_ else as.integer(seed),
        pair.mode = resolved$pair_resolution$effective,
        initial.energy = grip.geodesic.misf.kk.extract.score.metric(score, c("gkk.energy", "lgkk.energy")),
        initial.weighted.rel.rmse = grip.geodesic.misf.kk.extract.score.metric(score, c("gkk.weighted.rel.rmse", "lgkk.weighted.rel.rmse")),
        final.energy = grip.geodesic.misf.kk.extract.score.metric(score, c("gkk.energy", "lgkk.energy")),
        final.weighted.rel.rmse = grip.geodesic.misf.kk.extract.score.metric(score, c("gkk.weighted.rel.rmse", "lgkk.weighted.rel.rmse")),
        improved = TRUE,
        trace.rows = 0L,
        stringsAsFactors = FALSE
      ),
      best_restart = 1L,
      vertex_ids = as.integer(vertex.ids),
      coords_full = grip.geodesic.misf.partial.coords(coords, vertex.ids, prepared$n),
      effective_pair_mode = resolved$pair_resolution$effective,
      requested_pair_mode = resolved$pair_resolution$requested
    ))
  }

  restart.rows <- vector("list", n.restarts)
  best.fit <- NULL
  best.row <- NULL
  best.restart <- 1L
  best.energy <- Inf
  geometric.init <- NULL
  if (identical(init, "geometric")) {
    insertion.order <- if (!is.null(prepared$insertion_order)) {
      prepared$insertion_order[prepared$insertion_order %in% vertex.ids]
    } else {
      vertex.ids
    }
    anchor.count <- if (!is.null(prepared$insertion_anchor_count)) {
      prepared$insertion_anchor_count
    } else {
      grip.geodesic.misf.default.anchor.count(dim)
    }
    anchor.mode <- if (!is.null(prepared$insertion_anchor_weight_mode)) {
      prepared$insertion_anchor_weight_mode
    } else {
      "inverse_graph_distance_sq"
    }
    geometric.init <- grip.geodesic.misf.build.geometric.seed.coords(
      distance.matrix = prepared$top_level_graph$distance_matrix,
      dim = dim,
      vertex.ids = vertex.ids,
      insertion.order = insertion.order,
      anchor.count = anchor.count,
      anchor.weight.mode = anchor.mode
    )
  }

  for (restart in seq_len(n.restarts)) {
    restart.seed <- if (is.null(seed)) NULL else as.integer(seed + restart - 1L)
    init.coords <- if (identical(init, "geometric")) {
      grip.geodesic.misf.jitter.coords(
        geometric.init$coords,
        restart = restart,
        seed = restart.seed
      )
    } else {
      restart.init <- if (identical(init, "cmdscale") && restart == 1L) "cmdscale" else "random"
      grip.geodesic.misf.kk.init.coords(
        prepared = coarse.prepared,
        dim = dim,
        init = restart.init,
        restart = restart,
        seed = restart.seed
      )
    }
    init.score <- grip.geodesic.misf.kk.score.prepared(
      coords = init.coords,
      prepared = coarse.prepared,
      stiffness = stiffness,
      distance.floor = distance.floor,
      edge.length.epsilon = edge.length.epsilon
    )
    fit <- grip.geodesic.misf.kk.optimize.prepared(
      coords = init.coords,
      prepared = coarse.prepared,
      max.iter = max.iter,
      stiffness = stiffness,
      distance.floor = distance.floor,
      edge.length.epsilon = edge.length.epsilon,
      initial.step = initial.step,
      step.shrink = step.shrink,
      armijo.factor = armijo.factor,
      grad.tol = grad.tol,
      min.step = min.step,
      recenter = recenter,
      return.trace = return.trace
    )
    restart.rows[[restart]] <- grip.geodesic.misf.kk.build.restart.row(
      restart = restart,
      seed = restart.seed,
      init.score = init.score,
      fit = fit
    )
    final.energy <- grip.geodesic.misf.kk.extract.score.metric(
      fit$score,
      c("gkk.energy", "lgkk.energy")
    )
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
  best.fit$coords_full <- grip.geodesic.misf.partial.coords(best.fit$coords, vertex.ids, prepared$n)
  best.fit$effective_pair_mode <- resolved$pair_resolution$effective
  best.fit$requested_pair_mode <- resolved$pair_resolution$requested
  best.fit$pair_resolution <- resolved$pair_resolution
  best.fit$top_level_init <- init
  best.fit$initial_placement <- geometric.init
  best.fit
}

grip.geodesic.misf.kk.prepare.active.level <- function(prepared,
                                                       level = NULL,
                                                       local.nbrs = 8L,
                                                       landmark.count = 4L,
                                                       pair.mode = c("auto", "full", "landmark"),
                                                       full.limit = 256L) {
  prepared <- grip.validate.misf.geodesic.kk.prepared(prepared)
  pair.mode <- match.arg(pair.mode)
  if (is.null(level)) {
    level <- prepared$top_level_level
  }
  level <- as.integer(level)
  active.vertices <- grip.geodesic.misf.active.level.vertices(prepared, level)
  resolution <- grip.geodesic.misf.kk.resolve.pair.mode(
    pair.mode = pair.mode,
    active.n = length(active.vertices),
    full.limit = full.limit
  )
  active.prepared <- grip.geodesic.misf.prepare.active.level(
    prepared = prepared,
    active.vertices = active.vertices,
    local.nbrs = local.nbrs,
    landmark.count = landmark.count,
    pair.mode = if (identical(resolution$effective, "full")) "full" else "sparse"
  )
  active.prepared <- grip.geodesic.misf.kk.add.prepared.class(
    active.prepared,
    resolution$effective
  )
  list(
    level = level,
    active_vertices = active.vertices,
    active_prepared = active.prepared,
    pair_resolution = resolution,
    pair_matrix = if (!is.null(active.prepared$pair_matrix_local)) active.prepared$pair_matrix_local else active.prepared$pair_matrix,
    global_pair_matrix = matrix(as.integer(active.prepared$pair_matrix), ncol = 2L)
  )
}

grip.geodesic.misf.kk.refine.level <- function(prepared,
                                               coords,
                                               level = NULL,
                                               local.nbrs = 8L,
                                               landmark.count = 4L,
                                               pair.mode = c("auto", "full", "landmark"),
                                               full.limit = 256L,
                                               max.iter = 8L,
                                               anchor.weight = 0.05,
                                               anchor.weight.end = anchor.weight,
                                               continuation = c("constant", "linear", "geometric"),
                                               stiffness = 1.0,
                                               distance.floor = 1e-8,
                                               edge.length.epsilon = 1e-8,
                                               initial.step = 1.0,
                                               step.shrink = 0.5,
                                               armijo.factor = 1e-4,
                                               grad.tol = 1e-8,
                                               min.step = 1e-8,
                                               recenter = FALSE,
                                               return.trace = TRUE) {
  prepared <- grip.validate.misf.geodesic.kk.prepared(prepared)
  coords <- grip.validate.coords(coords)
  if (nrow(coords) != prepared$n) {
    stop("nrow(coords) must match prepared$n")
  }
  continuation <- match.arg(continuation)
  built <- grip.geodesic.misf.kk.prepare.active.level(
    prepared = prepared,
    level = level,
    local.nbrs = local.nbrs,
    landmark.count = landmark.count,
    pair.mode = pair.mode,
    full.limit = full.limit
  )
  before <- grip.geodesic.misf.kk.score.prepared(
    coords = coords,
    prepared = built$active_prepared,
    stiffness = stiffness,
    distance.floor = distance.floor,
    edge.length.epsilon = edge.length.epsilon
  )
  anchor.coords <- coords
  anchor.vertex.weight <- numeric(nrow(coords))
  inactive.vertices <- setdiff(seq_len(prepared$n), built$active_vertices)
  if (length(inactive.vertices) > 0L) {
    anchor.vertex.weight[inactive.vertices] <- 1
  }
  if (built$level < prepared$top_level_level) {
    pinned.global <- prepared$misf$levels[[built$level + 2L]]
    anchor.vertex.weight[pinned.global] <- 1
  }
  use.anchor <- any(anchor.vertex.weight > 0) && is.finite(anchor.weight.end) && anchor.weight.end > 0
  fit <- grip.geodesic.misf.kk.optimize.prepared(
    coords = coords,
    prepared = built$active_prepared,
    max.iter = max.iter,
    anchor.coords = if (use.anchor) anchor.coords else NULL,
    anchor.weight = if (use.anchor) anchor.weight else 0,
    anchor.weight.end = if (use.anchor) anchor.weight.end else 0,
    anchor.vertex.weight = if (use.anchor) anchor.vertex.weight else NULL,
    continuation = continuation,
    stiffness = stiffness,
    distance.floor = distance.floor,
    edge.length.epsilon = edge.length.epsilon,
    initial.step = initial.step,
    step.shrink = step.shrink,
    armijo.factor = armijo.factor,
    grad.tol = grad.tol,
    min.step = min.step,
    recenter = recenter,
    return.trace = return.trace
  )
  after <- grip.geodesic.misf.kk.score.prepared(
    coords = fit$coords,
    prepared = built$active_prepared,
    stiffness = stiffness,
    distance.floor = distance.floor,
    edge.length.epsilon = edge.length.epsilon
  )
  list(
    level = built$level,
    active_vertices = built$active_vertices,
    active_prepared = built$active_prepared,
    pair_matrix = built$pair_matrix,
    global_pair_matrix = built$global_pair_matrix,
    effective_pair_mode = built$pair_resolution$effective,
    requested_pair_mode = built$pair_resolution$requested,
    pinned_vertices = if (use.anchor) which(anchor.vertex.weight > 0) else integer(0L),
    anchor_vertex_weight = anchor.vertex.weight,
    anchor_weight = if (use.anchor) as.double(anchor.weight) else 0,
    anchor_weight_end = if (use.anchor) as.double(anchor.weight.end) else 0,
    anchor_continuation = continuation,
    coords = fit$coords,
    fit = fit,
    before = before,
    after = after
  )
}

grip.geodesic.misf.kk.refine.all.levels <- function(prepared,
                                                    coords,
                                                    local.nbrs = 8L,
                                                    landmark.count = 4L,
                                                    pair.mode = c("auto", "full", "landmark"),
                                                    full.limit = 256L,
                                                    max.iter = 8L,
                                                    anchor.weight = 0.05,
                                                    anchor.weight.end = anchor.weight,
                                                    continuation = c("constant", "linear", "geometric"),
                                                    stiffness = 1.0,
                                                    distance.floor = 1e-8,
                                                    edge.length.epsilon = 1e-8,
                                                    initial.step = 1.0,
                                                    step.shrink = 0.5,
                                                    armijo.factor = 1e-4,
                                                    grad.tol = 1e-8,
                                                    min.step = 1e-8,
                                                    recenter = FALSE,
                                                    return.trace = TRUE) {
  prepared <- grip.validate.misf.geodesic.kk.prepared(prepared)
  coords <- grip.validate.coords(coords)
  pair.mode <- match.arg(pair.mode)
  continuation <- match.arg(continuation)

  level.ids <- seq.int(from = prepared$top_level_level, to = 0L, by = -1L)
  level.results <- vector("list", length(level.ids))
  for (i in seq_along(level.ids)) {
    level.results[[i]] <- grip.geodesic.misf.kk.refine.level(
      prepared = prepared,
      coords = coords,
      level = level.ids[[i]],
      local.nbrs = local.nbrs,
      landmark.count = landmark.count,
      pair.mode = pair.mode,
      full.limit = full.limit,
      max.iter = max.iter,
      anchor.weight = anchor.weight,
      anchor.weight.end = anchor.weight.end,
      continuation = continuation,
      stiffness = stiffness,
      distance.floor = distance.floor,
      edge.length.epsilon = edge.length.epsilon,
      initial.step = initial.step,
      step.shrink = step.shrink,
      armijo.factor = armijo.factor,
      grad.tol = grad.tol,
      min.step = min.step,
      recenter = recenter,
      return.trace = return.trace
    )
    coords <- level.results[[i]]$coords
  }

  level.trace <- do.call(rbind, lapply(level.results, function(result) {
    data.frame(
      level = result$level,
      active_n = length(result$active_vertices),
      pair_n = nrow(result$pair_matrix),
      pair.mode = result$effective_pair_mode,
      before.energy = grip.geodesic.misf.kk.extract.score.metric(
        result$before,
        c("gkk.energy", "lgkk.energy")
      ),
      after.energy = grip.geodesic.misf.kk.extract.score.metric(
        result$after,
        c("gkk.energy", "lgkk.energy")
      ),
      after.weighted.rel.rmse = grip.geodesic.misf.kk.extract.score.metric(
        result$after,
        c("gkk.weighted.rel.rmse", "lgkk.weighted.rel.rmse")
      ),
      stringsAsFactors = FALSE
    )
  }))

  list(
    coords = coords,
    level_results = level.results,
    level_trace = level.trace
  )
}

grip.geodesic.misf.kk.resolve.final.prepared <- function(prepared,
                                                         pair.mode = c("auto", "full", "landmark"),
                                                         full.limit = 1024L,
                                                         local.nbrs = 20L,
                                                         landmark.count = 8L) {
  prepared <- grip.validate.misf.geodesic.kk.prepared(prepared)
  pair.mode <- match.arg(pair.mode)
  resolution <- grip.geodesic.misf.kk.resolve.pair.mode(
    pair.mode = pair.mode,
    active.n = prepared$n,
    full.limit = full.limit
  )
  final.prepared <- if (identical(resolution$effective, "full")) {
    grip.geodesic.misf.kk.add.prepared.class(prepared, "full")
  } else {
    grip.geodesic.misf.kk.add.prepared.class(
      prepare.landmark.geodesic.kk(
        n = prepared$n,
        adj.list = prepared$adj_list,
        weight.list = prepared$weight_list,
        local.nbrs = local.nbrs,
        landmark.count = landmark.count
      ),
      "landmark"
    )
  }
  list(
    prepared = final.prepared,
    pair_resolution = resolution
  )
}

grip.geodesic.misf.kk.final.polish <- function(prepared,
                                               coords,
                                               pair.mode = c("auto", "full", "landmark"),
                                               full.limit = 1024L,
                                               local.nbrs = 20L,
                                               landmark.count = 8L,
                                               max.iter = 8L,
                                               stiffness = 1.0,
                                               distance.floor = 1e-8,
                                               edge.length.epsilon = 1e-8,
                                               initial.step = 1.0,
                                               step.shrink = 0.5,
                                               armijo.factor = 1e-4,
                                               grad.tol = 1e-8,
                                               min.step = 1e-8,
                                               recenter = TRUE,
                                               return.trace = TRUE) {
  prepared <- grip.validate.misf.geodesic.kk.prepared(prepared)
  coords <- grip.validate.coords(coords)
  if (nrow(coords) != prepared$n) {
    stop("nrow(coords) must match prepared$n")
  }
  resolved <- grip.geodesic.misf.kk.resolve.final.prepared(
    prepared = prepared,
    pair.mode = pair.mode,
    full.limit = full.limit,
    local.nbrs = local.nbrs,
    landmark.count = landmark.count
  )
  fit <- grip.geodesic.misf.kk.optimize.prepared(
    coords = coords,
    prepared = resolved$prepared,
    max.iter = max.iter,
    stiffness = stiffness,
    distance.floor = distance.floor,
    edge.length.epsilon = edge.length.epsilon,
    initial.step = initial.step,
    step.shrink = step.shrink,
    armijo.factor = armijo.factor,
    grad.tol = grad.tol,
    min.step = min.step,
    recenter = recenter,
    return.trace = return.trace
  )
  fit$requested_pair_mode <- resolved$pair_resolution$requested
  fit$effective_pair_mode <- resolved$pair_resolution$effective
  fit
}

grip.geodesic.misf.kk.build.stage.bundle <- function(prepared,
                                                     top.level.fit,
                                                     top.level.elapsed,
                                                     top.level.frames,
                                                     insertion,
                                                     insertion.elapsed,
                                                     insertion.frames,
                                                     refinement,
                                                     refinement.elapsed,
                                                     refinement.frames,
                                                     final.polish,
                                                     final.polish.elapsed) {
  records <- list()

  top.initial <- top.level.fit$initial_placement
  if (!is.null(top.initial) && !is.null(top.initial$coords) && !is.null(top.initial$vertex_ids)) {
    top.vertex.ids <- as.integer(top.initial$vertex_ids)
    top.coords.full <- grip.geodesic.misf.partial.coords(
      coords = top.initial$coords,
      vertex.ids = top.vertex.ids,
      n = prepared$n
    )
    if (!is.null(top.initial$seed_vertices) && length(top.initial$seed_vertices)) {
      seed.vertices <- as.integer(top.initial$seed_vertices)
      seed.local <- match(seed.vertices, top.vertex.ids)
      seed.coords.full <- grip.geodesic.misf.partial.coords(
        coords = top.initial$coords[seed.local, , drop = FALSE],
        vertex.ids = seed.vertices,
        n = prepared$n
      )
      records[[length(records) + 1L]] <- grip.geodesic.misf.new.stage.record(
        stage = "seed",
        level = prepared$top_level_level,
        method.family = "geodesic_kk",
        coords.full = seed.coords.full,
        active.vertices = seed.vertices,
        frames = list(seed.coords.full),
        summary = list(
          active_n = length(seed.vertices),
          pair_n = if (length(seed.vertices) >= 2L) choose(length(seed.vertices), 2L) else 0L,
          frame_count = 1L
        )
      )
    }
    init.trace <- if (!is.null(top.initial$vertex_trace)) top.initial$vertex_trace else data.frame()
    records[[length(records) + 1L]] <- grip.geodesic.misf.new.stage.record(
      stage = "initial_placement",
      level = prepared$top_level_level,
      method.family = "geodesic_kk",
      coords.full = top.coords.full,
      active.vertices = top.vertex.ids,
      inserted.vertices = if (!is.null(top.initial$seed_vertices)) setdiff(top.vertex.ids, top.initial$seed_vertices) else integer(0L),
      trace = init.trace,
      frames = list(top.coords.full),
      summary = list(
        pair_n = NA_integer_,
        energy = if (nrow(init.trace)) mean(init.trace$objective) else NA_real_,
        mean_objective = if (nrow(init.trace)) mean(init.trace$objective) else NA_real_,
        max_grad_norm = if (nrow(init.trace)) max(init.trace$grad_norm) else NA_real_,
        all_converged = if (nrow(init.trace)) all(init.trace$converged) else NA,
        trace_rows = if (nrow(init.trace)) nrow(init.trace) else 0L,
        frame_count = 1L
      )
    )
  }

  records[[length(records) + 1L]] <- grip.geodesic.misf.new.stage.record(
    stage = "top_level",
    level = prepared$top_level_level,
    method.family = "geodesic_kk",
    coords.full = top.level.fit$coords_full,
    active.vertices = prepared$top_level_vertices,
    pair.mode = if (!is.null(top.level.fit$effective_pair_mode)) top.level.fit$effective_pair_mode else NA_character_,
    trace = top.level.fit$trace,
    frames = top.level.frames,
    summary = list(
      pair_n = length(top.level.fit$prepared$pair_graph_distance),
      energy = grip.geodesic.misf.kk.extract.score.metric(
        top.level.fit$score,
        c("gkk.energy", "lgkk.energy")
      ),
      weighted_rel_rmse = grip.geodesic.misf.kk.extract.score.metric(
        top.level.fit$score,
        c("gkk.weighted.rel.rmse", "lgkk.weighted.rel.rmse")
      ),
      elapsed_sec = as.double(top.level.elapsed)
    )
  )

  if (!is.null(insertion$level_results) && length(insertion$level_results)) {
    for (i in seq_along(insertion$level_results)) {
      result <- insertion$level_results[[i]]
      vertex.trace <- if (!is.null(result$vertex_trace)) result$vertex_trace else data.frame()
      active.vertices <- grip.geodesic.misf.active.level.vertices(prepared, result$level)
      records[[length(records) + 1L]] <- grip.geodesic.misf.new.stage.record(
        stage = "insertion",
        level = result$level,
        method.family = "geodesic_kk",
        coords.full = result$coords,
        active.vertices = active.vertices,
        inserted.vertices = result$inserted_vertices,
        trace = vertex.trace,
        frames = if (!is.null(insertion.frames[[paste0("level_", result$level)]])) {
          insertion.frames[[paste0("level_", result$level)]]
        } else {
          list(result$coords)
        },
        summary = list(
          active_n = length(active.vertices),
          inserted_n = length(result$inserted_vertices),
          energy = if (nrow(vertex.trace)) mean(vertex.trace$objective) else NA_real_,
          mean_objective = if (nrow(vertex.trace)) mean(vertex.trace$objective) else NA_real_,
          max_grad_norm = if (nrow(vertex.trace)) max(vertex.trace$grad_norm) else NA_real_,
          all_converged = if (nrow(vertex.trace)) all(vertex.trace$converged) else NA,
          elapsed_sec = if (i == length(insertion$level_results)) as.double(insertion.elapsed) else NA_real_,
          frame_count = 1L
        )
      )
    }
  }

  if (!is.null(refinement$level_results) && length(refinement$level_results)) {
    for (i in seq_along(refinement$level_results)) {
      result <- refinement$level_results[[i]]
      stage.trace <- if (!is.null(result$fit$trace)) result$fit$trace else data.frame()
      records[[length(records) + 1L]] <- grip.geodesic.misf.new.stage.record(
        stage = "refinement",
        level = result$level,
        method.family = "geodesic_kk",
        coords.full = result$coords,
        active.vertices = result$active_vertices,
        pair.mode = result$effective_pair_mode,
        trace = stage.trace,
        frames = if (!is.null(refinement.frames[[paste0("level_", result$level)]])) {
          refinement.frames[[paste0("level_", result$level)]]
        } else {
          list(result$coords)
        },
        summary = list(
          pair_n = nrow(result$pair_matrix),
          energy = grip.geodesic.misf.kk.extract.score.metric(
            result$after,
            c("gkk.energy", "lgkk.energy")
          ),
          weighted_rel_rmse = grip.geodesic.misf.kk.extract.score.metric(
            result$after,
            c("gkk.weighted.rel.rmse", "lgkk.weighted.rel.rmse")
          ),
          elapsed_sec = if (i == length(refinement$level_results)) as.double(refinement.elapsed) else NA_real_,
          frame_count = 1L
        )
      )
    }
  }

  records[[length(records) + 1L]] <- grip.geodesic.misf.new.stage.record(
    stage = "final_polish",
    level = 0L,
    method.family = "geodesic_kk",
    coords.full = final.polish$coords,
    active.vertices = seq_len(prepared$n),
    pair.mode = if (!is.null(final.polish$effective_pair_mode)) final.polish$effective_pair_mode else NA_character_,
    trace = final.polish$trace,
    frames = final.polish$frames,
    summary = list(
      pair_n = length(final.polish$prepared$pair_graph_distance),
      energy = grip.geodesic.misf.kk.extract.score.metric(
        final.polish$score,
        c("gkk.energy", "lgkk.energy")
      ),
      weighted_rel_rmse = grip.geodesic.misf.kk.extract.score.metric(
        final.polish$score,
        c("gkk.weighted.rel.rmse", "lgkk.weighted.rel.rmse")
      ),
      elapsed_sec = as.double(final.polish.elapsed)
    )
  )

  names(records) <- vapply(records, `[[`, character(1L), "stage_key")
  list(
    stage_trace = grip.geodesic.misf.stage.trace.from.records(records),
    stage_data = records
  )
}

grip.geodesic.misf.kk.build.stage.trace <- function(prepared,
                                                    top.level.fit,
                                                    top.level.elapsed,
                                                    top.level.frames,
                                                    insertion,
                                                    insertion.elapsed,
                                                    insertion.frames,
                                                    refinement,
                                                    refinement.elapsed,
                                                    refinement.frames,
                                                    final.polish,
                                                    final.polish.elapsed) {
  grip.geodesic.misf.kk.build.stage.bundle(
    prepared = prepared,
    top.level.fit = top.level.fit,
    top.level.elapsed = top.level.elapsed,
    top.level.frames = top.level.frames,
    insertion = insertion,
    insertion.elapsed = insertion.elapsed,
    insertion.frames = insertion.frames,
    refinement = refinement,
    refinement.elapsed = refinement.elapsed,
    refinement.frames = refinement.frames,
    final.polish = final.polish,
    final.polish.elapsed = final.polish.elapsed
  )$stage_trace
}

#' Prepare a MISF-based multiscale geodesic-KK object
#'
#' `prepare.misf.geodesic.kk()` builds the prepared state for a MISF-based
#' geodesic-KK pipeline by layering the maximal independent set filtration
#' (MISF) on top of the existing full geodesic-KK cache. The prepared object
#' stores the graph metadata, the coarsest admissible top-level graph, the
#' exact and sparse top-level prepared caches, and the default controls used by
#' the multiscale MISF-GKK optimizer.
#'
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#' @param n Number of vertices. If omitted with `adj.list`, defaults to
#'   `length(adj.list)`. If omitted with `edges`, defaults to `max(edges)`.
#' @param adj.list Adjacency list (1-based) for an undirected graph.
#' @param weight.list Optional parallel list of positive edge weights.
#' @param edge.weights Optional positive edge-weight vector parallel to `edges`.
#' @param tie.mode Shortest-path aggregation mode inherited from
#'   [prepare.geodesic.kk()].
#' @param num.init Target top-level active-set size passed to `build.misf`.
#' @param num.nbrs Per-level local-neighborhood schedule metadata passed to
#'   `build.misf`.
#' @param dim Target embedding dimension (`2` or `3`) for the multiscale solve.
#' @param top.level.mode Either `"solve"` or `"skip"`. When set to `"solve"`,
#'   the coarsest admissible MISF level is optimized immediately and stored in
#'   `prepared$top_level_fit`.
#' @param top.level.pair.mode Pair policy for the top MISF level: `"auto"`,
#'   `"full"`, or `"landmark"`.
#' @param top.level.full.limit Active-set size threshold used when
#'   `top.level.pair.mode = "auto"`.
#' @param top.level.local.nbrs Number of nearest graph-metric neighbors retained
#'   per vertex in the sparse top-level LGKK cache.
#' @param top.level.landmark.count Number of farthest-point landmarks retained
#'   per vertex in the sparse top-level LGKK cache.
#' @param top.level.restarts Number of top-level restarts used by the coarse
#'   MISF-GKK solve.
#' @param top.level.max.iter Top-level iteration budget.
#' @param top.level.init Top-level initializer (`"geometric"`, `"cmdscale"`, or
#'   `"random"`).
#' @param seed Optional integer seed reused for MISF extraction.
#'
#' @return An object of class `"grip_misf_gkk_prepared"` layered on top of the
#'   full geodesic-KK prepared structure, with added MISF metadata, top-level
#'   coarse-graph caches, and stored default controls for the MISF-GKK
#'   optimizer.
#'
#' @examples
#' edges <- edges.mesh(4, 4)
#' prepared <- prepare.misf.geodesic.kk(
#'   edges = edges,
#'   n = 16,
#'   tie.mode = "average",
#'   num.init = 4,
#'   top.level.mode = "skip",
#'   seed = 1
#' )
#' prepared$top_level_vertices
#' @export
#' @md
prepare.misf.geodesic.kk <- function(edges = NULL,
                                          n = NULL,
                                          adj.list = NULL,
                                          weight.list = NULL,
                                          edge.weights = NULL,
                                          tie.mode = c("single", "average"),
                                          num.init = 24L,
                                          num.nbrs = 20L,
                                          dim = 2L,
                                          top.level.mode = c("solve", "skip"),
                                          top.level.pair.mode = c("auto", "full", "landmark"),
                                          top.level.full.limit = 512L,
                                          top.level.local.nbrs = 20L,
                                          top.level.landmark.count = 8L,
                                          top.level.restarts = 8L,
                                          top.level.max.iter = 16L,
                                          top.level.init = c("geometric", "cmdscale", "random"),
                                          seed = 6L) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights)
  tie.mode <- match.arg(tie.mode)
  top.level.mode <- match.arg(top.level.mode)
  top.level.pair.mode <- match.arg(top.level.pair.mode)
  top.level.init <- match.arg(top.level.init)
  dim <- grip.validate.count(dim, "dim")
  if (!(dim %in% c(2L, 3L))) {
    stop("dim must be 2 or 3")
  }
  top.level.full.limit <- grip.validate.misf.count(
    top.level.full.limit,
    "top.level.full.limit",
    lower = 1L
  )
  top.level.local.nbrs <- grip.validate.count(top.level.local.nbrs, "top.level.local.nbrs")
  top.level.landmark.count <- grip.validate.count(top.level.landmark.count, "top.level.landmark.count")
  top.level.restarts <- grip.validate.misf.count(top.level.restarts, "top.level.restarts", lower = 1L)
  grip.validate.scalar(top.level.max.iter, "top.level.max.iter", lower = 0)
  top.level.max.iter <- as.integer(round(top.level.max.iter))
  if (!is.null(seed)) {
    seed <- grip.validate.count(seed, "seed")
  }

  prepared <- prepare.geodesic.kk(
    edges = edges,
    n = n,
    adj.list = adj.list,
    weight.list = weight.list,
    edge.weights = edge.weights,
    tie.mode = tie.mode
  )
  misf <- build.misf(
    n = prepared$n,
    adj.list = prepared$adj_list,
    weight.list = prepared$weight_list,
    num.init = num.init,
    num.nbrs = num.nbrs,
    seed = seed
  )
  top.level.selection <- grip.geodesic.misf.resolve.top.level(
    prepared = prepared,
    misf = misf,
    dim = dim,
    tie.mode = tie.mode
  )
  top_level_index <- top.level.selection$level_index
  top_level_id <- top.level.selection$level
  top_level_vertices <- top.level.selection$vertices
  top_level_graph <- top.level.selection$graph
  top_level_resolution <- grip.geodesic.misf.kk.resolve.pair.mode(
    pair.mode = top.level.pair.mode,
    active.n = top_level_graph$n,
    full.limit = top.level.full.limit
  )
  top_level_prepared_full <- prepare.geodesic.kk(
    edges = top_level_graph$edges,
    n = top_level_graph$n,
    edge.weights = top_level_graph$edge_weights,
    tie.mode = tie.mode
  )
  top_level_prepared_sparse <- prepare.landmark.geodesic.kk(
    edges = top_level_graph$edges,
    n = top_level_graph$n,
    edge.weights = top_level_graph$edge_weights,
    local.nbrs = top.level.local.nbrs,
    landmark.count = top.level.landmark.count
  )

  prepared$misf <- misf
  prepared$level_vertices <- misf$levels
  prepared$active_levels <- misf$levels
  prepared$insertion_order <- misf$mish_order
  prepared$coarsest_level_index <- length(misf$levels)
  prepared$coarsest_level_level <- as.integer(length(misf$levels) - 1L)
  prepared$top_level_index <- top_level_index
  prepared$top_level_level <- top_level_id
  prepared$top_level_vertices <- top_level_vertices
  prepared$top_level_graph <- top_level_graph
  prepared$top_level_min_required_size <- top.level.selection$min_required_size
  prepared$top_level_selection_reason <- top.level.selection$selection_reason
  prepared$top_level_seed_vertices <- top.level.selection$seed_vertices
  prepared$top_level_seed_positive_rank <- top.level.selection$seed_positive_rank
  prepared$top_level_seed_required_rank <- top.level.selection$seed_required_rank
  prepared$top_level_pair_mode <- top_level_resolution$requested
  prepared$top_level_effective_pair_mode <- top_level_resolution$effective
  prepared$top_level_full_limit <- top_level_resolution$full_limit
  prepared$top_level_local_nbrs <- top.level.local.nbrs
  prepared$top_level_landmark_count <- top.level.landmark.count
  prepared$top_level_prepared_full <- top_level_prepared_full
  prepared$top_level_prepared_sparse <- top_level_prepared_sparse
  prepared$top_level_prepared <- if (identical(top_level_resolution$effective, "full")) {
    top_level_prepared_full
  } else {
    top_level_prepared_sparse
  }
  prepared$top_level_dim <- as.integer(dim)
  prepared$top_level_mode <- top.level.mode
  prepared$top_level_restarts <- as.integer(top.level.restarts)
  prepared$top_level_max_iter <- as.integer(top.level.max.iter)
  prepared$top_level_init <- top.level.init
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
  prepared$refinement_anchor_weight <- 0.05
  prepared$refinement_anchor_weight_end <- 0.05
  prepared$refinement_continuation <- "constant"
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

  if (identical(top.level.mode, "solve")) {
    prepared$top_level_fit <- grip.geodesic.misf.kk.solve.top.level(
      prepared = prepared,
      dim = dim,
      pair.mode = top.level.pair.mode,
      full.limit = top.level.full.limit,
      local.nbrs = top.level.local.nbrs,
      landmark.count = top.level.landmark.count,
      n.restarts = top.level.restarts,
      max.iter = top.level.max.iter,
      init = top.level.init,
      stiffness = prepared$stiffness,
      distance.floor = prepared$distance_floor,
      edge.length.epsilon = prepared$edge_length_epsilon,
      initial.step = prepared$initial_step,
      step.shrink = prepared$step_shrink,
      armijo.factor = prepared$armijo_factor,
      grad.tol = prepared$grad_tol,
      min.step = prepared$min_step,
      recenter = prepared$recenter,
      return.trace = FALSE,
      seed = seed
    )
  }
  prepared
}

#' Optimize an embedding with the MISF-based geodesic-KK pipeline
#'
#' `misf.geodesic.kk()` runs a multiscale MISF-based geodesic-KK
#' pipeline. It solves the top MISF graph with either full GKK or sparse LGKK,
#' inserts lower-level vertices with the existing MISF placement helpers, refines
#' each active MISF level under a GKK/LGKK objective, and finishes with a final
#' full-graph polish.
#'
#' @inheritParams prepare.misf.geodesic.kk
#' @param prepared Optional prepared object. This can be either a full
#'   geodesic-KK prepared object or a MISF-GKK prepared object.
#' @param top.level.pair.mode Optional override for the top-level pair policy.
#' @param top.level.full.limit Optional override for the top-level exact/sparse
#'   threshold.
#' @param top.level.local.nbrs Optional override for the sparse top-level local
#'   neighborhood size.
#' @param top.level.landmark.count Optional override for the sparse top-level
#'   landmark count.
#' @param top.level.init Optional override for the top-level initializer.
#' @param insertion.anchor.policy Optional insertion anchor policy.
#' @param insertion.anchor.count Optional insertion anchor count.
#' @param insertion.anchor.weight.mode Optional insertion anchor weighting mode.
#' @param insertion.max.iter Optional per-vertex insertion iteration budget.
#' @param insertion.mode Optional insertion warm-start mode.
#' @param insertion.layout.k Optional active-level layout graph size used by
#'   layout-based insertion modes.
#' @param insertion.weighted.preset Optional preset passed to weighted GRIP
#'   insertion.
#' @param insertion.grip.args Optional named list of extra arguments passed to
#'   combinatorial GRIP insertion.
#' @param insertion.weighted.args Optional named list of extra arguments for
#'   weighted GRIP insertion.
#' @param insertion.fr.niter Optional FR iteration budget for FR-based
#'   insertion.
#' @param refinement.pair.mode Optional active-level pair policy.
#' @param refinement.full.limit Optional active-level exact/sparse threshold.
#' @param refinement.local.nbrs Optional active-level sparse local neighborhood
#'   size.
#' @param refinement.landmark.count Optional active-level sparse landmark count.
#' @param refinement.anchor.weight Optional initial anchor weight used to pin
#'   inactive and coarser-level vertices during KK refinement.
#' @param refinement.anchor.weight.end Optional final anchor weight used at the
#'   end of the active-level continuation schedule.
#' @param refinement.continuation Optional continuation schedule used for the
#'   active-level anchor penalty.
#' @param refinement.max.iter Optional active-level refinement iteration budget.
#' @param final.pair.mode Optional final full-graph pair policy.
#' @param final.full.limit Optional final full-graph exact/sparse threshold.
#' @param final.local.nbrs Optional final full-graph sparse local neighborhood
#'   size.
#' @param final.landmark.count Optional final full-graph sparse landmark count.
#' @param final.max.iter Optional final polish iteration budget.
#' @param stiffness Global stiffness constant \(K\).
#' @param distance.floor Small positive floor used in
#'   \code{k_ij = K / max(g_ij, distance.floor)^2}.
#' @param edge.length.epsilon Small positive stabilizer added inside each
#'   embedded edge length.
#' @param initial.step Initial line-search step size.
#' @param step.shrink Multiplicative shrink factor in `(0, 1)` for backtracking.
#' @param armijo.factor Non-negative Armijo decrease constant.
#' @param grad.tol Non-negative stopping tolerance on the gradient norm.
#' @param min.step Positive minimum accepted line-search step before giving up.
#' @param recenter If `TRUE`, recenter accepted proposals to zero mean after
#'   each accepted step.
#' @param return.trace If `TRUE`, retain detailed stage traces.
#' @param return.frames If `TRUE`, retain intermediate coordinate frames.
#'
#' @return A list with `coords`, `prepared`, the per-stage multiscale results,
#'   the stage trace, optional trace/frame details, timing, and the final MISF
#'   score summary.
#'
#' @export
#' @md
misf.geodesic.kk <- function(prepared = NULL,
                                           edges = NULL,
                                           n = NULL,
                                           adj.list = NULL,
                                           weight.list = NULL,
                                           edge.weights = NULL,
                                           tie.mode = NULL,
                                           num.init = 24L,
                                           num.nbrs = 20L,
                                           dim = NULL,
                                           top.level.pair.mode = NULL,
                                           top.level.full.limit = NULL,
                                           top.level.local.nbrs = NULL,
                                           top.level.landmark.count = NULL,
                                           top.level.restarts = NULL,
                                           top.level.max.iter = NULL,
                                           top.level.init = NULL,
                                           insertion.anchor.policy = NULL,
                                           insertion.anchor.count = NULL,
                                           insertion.anchor.weight.mode = NULL,
                                           insertion.max.iter = NULL,
                                           insertion.mode = NULL,
                                           insertion.layout.k = NULL,
                                           insertion.weighted.preset = NULL,
                                           insertion.grip.args = NULL,
                                           insertion.weighted.args = NULL,
                                           insertion.fr.niter = NULL,
                                           refinement.pair.mode = NULL,
                                           refinement.full.limit = NULL,
                                           refinement.local.nbrs = NULL,
                                           refinement.landmark.count = NULL,
                                           refinement.anchor.weight = NULL,
                                           refinement.anchor.weight.end = NULL,
                                           refinement.continuation = NULL,
                                           refinement.max.iter = NULL,
                                           final.pair.mode = NULL,
                                           final.full.limit = NULL,
                                           final.local.nbrs = NULL,
                                           final.landmark.count = NULL,
                                           final.max.iter = NULL,
                                           stiffness = 1.0,
                                           distance.floor = 1e-8,
                                           edge.length.epsilon = 1e-8,
                                           initial.step = 1.0,
                                           step.shrink = 0.5,
                                           armijo.factor = 1e-4,
                                           grad.tol = 1e-8,
                                           min.step = 1e-8,
                                           recenter = TRUE,
                                           return.trace = FALSE,
                                           return.frames = FALSE,
                                           seed = 6L) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights, prepared)
  if (!is.null(seed)) {
    seed <- grip.validate.count(seed, "seed")
  }

  dim.resolved <- if (is.null(dim) && !is.null(prepared) && inherits(prepared, "grip_misf_gkk_prepared")) {
    prepared$top_level_dim
  } else if (is.null(dim)) {
    2L
  } else {
    grip.validate.count(dim, "dim")
  }
  if (!(dim.resolved %in% c(2L, 3L))) {
    stop("dim must be 2 or 3")
  }

  top.level.pair.mode <- if (is.null(top.level.pair.mode) && !is.null(prepared) && inherits(prepared, "grip_misf_gkk_prepared")) {
    prepared$top_level_pair_mode
  } else if (is.null(top.level.pair.mode)) {
    "auto"
  } else {
    match.arg(top.level.pair.mode, c("auto", "full", "landmark"))
  }
  top.level.full.limit <- if (is.null(top.level.full.limit) && !is.null(prepared) && inherits(prepared, "grip_misf_gkk_prepared")) {
    prepared$top_level_full_limit
  } else if (is.null(top.level.full.limit)) {
    512L
  } else {
    grip.validate.misf.count(top.level.full.limit, "top.level.full.limit", lower = 1L)
  }
  top.level.local.nbrs <- if (is.null(top.level.local.nbrs) && !is.null(prepared) && inherits(prepared, "grip_misf_gkk_prepared")) {
    prepared$top_level_local_nbrs
  } else if (is.null(top.level.local.nbrs)) {
    20L
  } else {
    grip.validate.count(top.level.local.nbrs, "top.level.local.nbrs")
  }
  top.level.landmark.count <- if (is.null(top.level.landmark.count) && !is.null(prepared) && inherits(prepared, "grip_misf_gkk_prepared")) {
    prepared$top_level_landmark_count
  } else if (is.null(top.level.landmark.count)) {
    8L
  } else {
    grip.validate.count(top.level.landmark.count, "top.level.landmark.count")
  }
  top.level.restarts <- if (is.null(top.level.restarts) && !is.null(prepared) && inherits(prepared, "grip_misf_gkk_prepared")) {
    prepared$top_level_restarts
  } else if (is.null(top.level.restarts)) {
    8L
  } else {
    grip.validate.misf.count(top.level.restarts, "top.level.restarts", lower = 1L)
  }
  top.level.max.iter <- if (is.null(top.level.max.iter) && !is.null(prepared) && inherits(prepared, "grip_misf_gkk_prepared")) {
    prepared$top_level_max_iter
  } else if (is.null(top.level.max.iter)) {
    16L
  } else {
    grip.validate.scalar(top.level.max.iter, "top.level.max.iter", lower = 0)
    as.integer(round(top.level.max.iter))
  }
  top.level.init <- if (is.null(top.level.init) && !is.null(prepared) && inherits(prepared, "grip_misf_gkk_prepared")) {
    if (!is.null(prepared$top_level_init)) prepared$top_level_init else "geometric"
  } else if (is.null(top.level.init)) {
    "geometric"
  } else {
    match.arg(top.level.init, c("geometric", "cmdscale", "random"))
  }

  prepared <- grip.resolve.misf.geodesic.kk.prepared(
    prepared = prepared,
    edges = edges,
    n = n,
    adj.list = adj.list,
    weight.list = weight.list,
    edge.weights = edge.weights,
    tie.mode = tie.mode,
    num.init = num.init,
    num.nbrs = num.nbrs,
    dim = dim.resolved,
    top.level.pair.mode = top.level.pair.mode,
    top.level.full.limit = top.level.full.limit,
    top.level.local.nbrs = top.level.local.nbrs,
    top.level.landmark.count = top.level.landmark.count,
    top.level.restarts = top.level.restarts,
    top.level.max.iter = top.level.max.iter,
    top.level.init = top.level.init,
    seed = seed
  )

  insertion.anchor.policy <- if (is.null(insertion.anchor.policy)) {
    prepared$insertion_anchor_policy
  } else {
    match.arg(insertion.anchor.policy, c("prev_level_first", "prev_level_distance_band", "prev_level_spread"))
  }
  insertion.anchor.count <- if (is.null(insertion.anchor.count)) {
    prepared$insertion_anchor_count
  } else {
    grip.validate.misf.count(insertion.anchor.count, "insertion.anchor.count", lower = 1L)
  }
  insertion.anchor.weight.mode <- if (is.null(insertion.anchor.weight.mode)) {
    prepared$insertion_anchor_weight_mode
  } else {
    match.arg(insertion.anchor.weight.mode, c("inverse_graph_distance_sq", "uniform"))
  }
  insertion.max.iter <- if (is.null(insertion.max.iter)) {
    prepared$insertion_max_iter
  } else {
    grip.validate.scalar(insertion.max.iter, "insertion.max.iter", lower = 0)
    as.integer(round(insertion.max.iter))
  }
  insertion.mode.resolved <- if (is.null(insertion.mode)) {
    prepared$insertion_mode
  } else {
    match.arg(insertion.mode, c("geodesic", "kk", "weighted_kk", "fr", "grip", "weighted_grip"))
  }
  insertion.layout.k <- if (is.null(insertion.layout.k)) {
    prepared$insertion_layout_k
  } else {
    grip.validate.count(insertion.layout.k, "insertion.layout.k")
  }
  insertion.weighted.preset <- if (is.null(insertion.weighted.preset)) {
    prepared$insertion_weighted_preset
  } else {
    as.character(insertion.weighted.preset)
  }
  insertion.grip.args <- if (is.null(insertion.grip.args)) {
    prepared$insertion_grip_args
  } else {
    if (!is.list(insertion.grip.args)) {
      stop("insertion.grip.args must be NULL or a list")
    }
    insertion.grip.args
  }
  insertion.weighted.args <- if (is.null(insertion.weighted.args)) {
    prepared$insertion_weighted_args
  } else {
    if (!is.list(insertion.weighted.args)) {
      stop("insertion.weighted.args must be NULL or a list")
    }
    insertion.weighted.args
  }
  insertion.fr.niter <- if (is.null(insertion.fr.niter)) {
    prepared$insertion_fr_niter
  } else {
    grip.validate.count(insertion.fr.niter, "insertion.fr.niter")
  }

  refinement.pair.mode <- if (is.null(refinement.pair.mode)) {
    prepared$refinement_pair_mode
  } else {
    match.arg(refinement.pair.mode, c("auto", "full", "landmark"))
  }
  refinement.full.limit <- if (is.null(refinement.full.limit)) {
    prepared$refinement_full_limit
  } else {
    grip.validate.misf.count(refinement.full.limit, "refinement.full.limit", lower = 1L)
  }
  refinement.local.nbrs <- if (is.null(refinement.local.nbrs)) {
    prepared$refinement_local_nbrs
  } else {
    grip.validate.count(refinement.local.nbrs, "refinement.local.nbrs")
  }
  refinement.landmark.count <- if (is.null(refinement.landmark.count)) {
    prepared$refinement_landmark_count
  } else {
    grip.validate.count(refinement.landmark.count, "refinement.landmark.count")
  }
  refinement.anchor.weight <- if (is.null(refinement.anchor.weight)) {
    prepared$refinement_anchor_weight
  } else {
    grip.validate.scalar(refinement.anchor.weight, "refinement.anchor.weight", lower = 0)
    as.double(refinement.anchor.weight)
  }
  refinement.anchor.weight.end <- if (is.null(refinement.anchor.weight.end)) {
    if (!is.null(prepared$refinement_anchor_weight_end)) prepared$refinement_anchor_weight_end else refinement.anchor.weight
  } else {
    grip.validate.scalar(refinement.anchor.weight.end, "refinement.anchor.weight.end", lower = 0)
    as.double(refinement.anchor.weight.end)
  }
  refinement.continuation <- if (is.null(refinement.continuation)) {
    prepared$refinement_continuation
  } else {
    match.arg(refinement.continuation, c("constant", "linear", "geometric"))
  }
  refinement.max.iter <- if (is.null(refinement.max.iter)) {
    prepared$refinement_max_iter
  } else {
    grip.validate.scalar(refinement.max.iter, "refinement.max.iter", lower = 0)
    as.integer(round(refinement.max.iter))
  }

  final.pair.mode <- if (is.null(final.pair.mode)) {
    prepared$final_pair_mode
  } else {
    match.arg(final.pair.mode, c("auto", "full", "landmark"))
  }
  final.full.limit <- if (is.null(final.full.limit)) {
    prepared$final_full_limit
  } else {
    grip.validate.misf.count(final.full.limit, "final.full.limit", lower = 1L)
  }
  final.local.nbrs <- if (is.null(final.local.nbrs)) {
    prepared$final_local_nbrs
  } else {
    grip.validate.count(final.local.nbrs, "final.local.nbrs")
  }
  final.landmark.count <- if (is.null(final.landmark.count)) {
    prepared$final_landmark_count
  } else {
    grip.validate.count(final.landmark.count, "final.landmark.count")
  }
  final.max.iter <- if (is.null(final.max.iter)) {
    prepared$final_max_iter
  } else {
    grip.validate.scalar(final.max.iter, "final.max.iter", lower = 0)
    as.integer(round(final.max.iter))
  }

  grip.validate.scalar(stiffness, "stiffness", lower = 0, open.lower = TRUE)
  grip.validate.scalar(distance.floor, "distance.floor", lower = 0, open.lower = TRUE)
  grip.validate.scalar(edge.length.epsilon, "edge.length.epsilon", lower = 0)
  grip.validate.scalar(initial.step, "initial.step", lower = 0, open.lower = TRUE)
  grip.validate.scalar(step.shrink, "step.shrink", lower = 0, upper = 1, open.lower = TRUE, open.upper = TRUE)
  grip.validate.scalar(armijo.factor, "armijo.factor", lower = 0)
  grip.validate.scalar(grad.tol, "grad.tol", lower = 0)
  grip.validate.scalar(min.step, "min.step", lower = 0, open.lower = TRUE)
  if (!is.logical(recenter) || length(recenter) != 1L || is.na(recenter)) {
    stop("recenter must be TRUE or FALSE")
  }
  if (!is.logical(return.trace) || length(return.trace) != 1L || is.na(return.trace)) {
    stop("return.trace must be TRUE or FALSE")
  }
  if (!is.logical(return.frames) || length(return.frames) != 1L || is.na(return.frames)) {
    stop("return.frames must be TRUE or FALSE")
  }

  prepared$top_level_pair_mode <- top.level.pair.mode
  prepared$top_level_full_limit <- top.level.full.limit
  prepared$top_level_local_nbrs <- top.level.local.nbrs
  prepared$top_level_landmark_count <- top.level.landmark.count
  prepared$top_level_restarts <- top.level.restarts
  prepared$top_level_max_iter <- top.level.max.iter
  prepared$top_level_init <- top.level.init
  prepared$top_level_dim <- dim.resolved
  prepared$insertion_anchor_policy <- insertion.anchor.policy
  prepared$insertion_anchor_count <- insertion.anchor.count
  prepared$insertion_anchor_weight_mode <- insertion.anchor.weight.mode
  prepared$insertion_max_iter <- insertion.max.iter
  prepared$insertion_mode <- insertion.mode.resolved
  prepared$insertion_layout_k <- insertion.layout.k
  prepared$insertion_weighted_preset <- insertion.weighted.preset
  prepared$insertion_grip_args <- insertion.grip.args
  prepared$insertion_weighted_args <- insertion.weighted.args
  prepared$insertion_fr_niter <- insertion.fr.niter
  prepared$refinement_pair_mode <- refinement.pair.mode
  prepared$refinement_full_limit <- refinement.full.limit
  prepared$refinement_local_nbrs <- refinement.local.nbrs
  prepared$refinement_landmark_count <- refinement.landmark.count
  prepared$refinement_anchor_weight <- refinement.anchor.weight
  prepared$refinement_anchor_weight_end <- refinement.anchor.weight.end
  prepared$refinement_continuation <- refinement.continuation
  prepared$refinement_max_iter <- refinement.max.iter
  prepared$final_pair_mode <- final.pair.mode
  prepared$final_full_limit <- final.full.limit
  prepared$final_local_nbrs <- final.local.nbrs
  prepared$final_landmark_count <- final.landmark.count
  prepared$final_max_iter <- final.max.iter
  prepared$stiffness <- stiffness
  prepared$distance_floor <- distance.floor
  prepared$edge_length_epsilon <- edge.length.epsilon
  prepared$initial_step <- initial.step
  prepared$step_shrink <- step.shrink
  prepared$armijo_factor <- armijo.factor
  prepared$grad_tol <- grad.tol
  prepared$min_step <- min.step
  prepared$recenter <- recenter

  need.top.trace <- isTRUE(return.trace) || isTRUE(return.frames)
  top.override.requested <- any(!vapply(
    list(
      top.level.pair.mode,
      top.level.full.limit,
      top.level.local.nbrs,
      top.level.landmark.count,
      top.level.restarts,
      top.level.max.iter,
      top.level.init,
      dim
    ),
    is.null,
    logical(1L)
  ))
  need.top.solve <- is.null(prepared$top_level_fit) || isTRUE(top.override.requested)
  if (!need.top.solve) {
    need.top.solve <- ncol(prepared$top_level_fit$coords) != dim.resolved
    if (!need.top.solve && isTRUE(return.trace)) {
      need.top.solve <- is.null(prepared$top_level_fit$trace) || nrow(prepared$top_level_fit$trace) == 0L
    }
    if (!need.top.solve && isTRUE(return.frames)) {
      need.top.solve <- is.null(prepared$top_level_fit$frames) || length(prepared$top_level_fit$frames) <= 1L
    }
  }

  top.level.elapsed <- 0
  if (need.top.solve) {
    top.level.start <- proc.time()[["elapsed"]]
    prepared$top_level_fit <- grip.geodesic.misf.kk.solve.top.level(
      prepared = prepared,
      dim = dim.resolved,
      pair.mode = top.level.pair.mode,
      full.limit = top.level.full.limit,
      local.nbrs = top.level.local.nbrs,
      landmark.count = top.level.landmark.count,
      n.restarts = top.level.restarts,
      max.iter = top.level.max.iter,
      init = top.level.init,
      stiffness = stiffness,
      distance.floor = distance.floor,
      edge.length.epsilon = edge.length.epsilon,
      initial.step = initial.step,
      step.shrink = step.shrink,
      armijo.factor = armijo.factor,
      grad.tol = grad.tol,
      min.step = min.step,
      recenter = recenter,
      return.trace = need.top.trace,
      seed = seed
    )
    top.level.elapsed <- proc.time()[["elapsed"]] - top.level.start
  }
  top.level.fit <- prepared$top_level_fit
  coords <- top.level.fit$coords_full

  insertion.elapsed <- 0
  if (any(!is.finite(coords))) {
    insertion.start <- proc.time()[["elapsed"]]
    insertion <- if (identical(insertion.mode.resolved, "geodesic")) {
      grip.geodesic.misf.insert.all.levels(
        prepared = prepared,
        coords = coords,
        anchor.policy = insertion.anchor.policy,
        anchor.count = insertion.anchor.count,
        anchor.weight.mode = insertion.anchor.weight.mode,
        max.iter = insertion.max.iter,
        initial.step = initial.step,
        step.shrink = step.shrink,
        armijo.factor = armijo.factor,
        grad.tol = grad.tol,
        min.step = min.step
      )
    } else {
      grip.geodesic.misf.insert.all.levels.with.layout(
        prepared = prepared,
        coords = coords,
        method = insertion.mode.resolved,
        layout.k = insertion.layout.k,
        weighted.preset = insertion.weighted.preset,
        grip.args = insertion.grip.args,
        weighted.args = insertion.weighted.args,
        fr.niter = insertion.fr.niter,
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
  refinement <- grip.geodesic.misf.kk.refine.all.levels(
    prepared = prepared,
    coords = insertion$coords,
    local.nbrs = refinement.local.nbrs,
    landmark.count = refinement.landmark.count,
    pair.mode = refinement.pair.mode,
    full.limit = refinement.full.limit,
    max.iter = refinement.max.iter,
    anchor.weight = refinement.anchor.weight,
    anchor.weight.end = refinement.anchor.weight.end,
    continuation = refinement.continuation,
    stiffness = stiffness,
    distance.floor = distance.floor,
    edge.length.epsilon = edge.length.epsilon,
    initial.step = initial.step,
    step.shrink = step.shrink,
    armijo.factor = armijo.factor,
    grad.tol = grad.tol,
    min.step = min.step,
    recenter = FALSE,
    return.trace = isTRUE(return.trace) || isTRUE(return.frames)
  )
  refinement.elapsed <- proc.time()[["elapsed"]] - refinement.start

  final.polish.start <- proc.time()[["elapsed"]]
  final.polish <- grip.geodesic.misf.kk.final.polish(
    prepared = prepared,
    coords = refinement$coords,
    pair.mode = final.pair.mode,
    full.limit = final.full.limit,
    local.nbrs = final.local.nbrs,
    landmark.count = final.landmark.count,
    max.iter = final.max.iter,
    stiffness = stiffness,
    distance.floor = distance.floor,
    edge.length.epsilon = edge.length.epsilon,
    initial.step = initial.step,
    step.shrink = step.shrink,
    armijo.factor = armijo.factor,
    grad.tol = grad.tol,
    min.step = min.step,
    recenter = recenter,
    return.trace = isTRUE(return.trace) || isTRUE(return.frames)
  )
  final.polish.elapsed <- proc.time()[["elapsed"]] - final.polish.start

  top.level.frames <- if (isTRUE(return.frames)) {
    grip.geodesic.misf.expand.top.level.frames(top.level.fit, prepared$n)
  } else {
    NULL
  }
  insertion.level.frames <- if (isTRUE(return.frames)) {
    grip.geodesic.misf.collect.level.frames(insertion$level_results)
  } else {
    NULL
  }
  refinement.level.frames <- if (isTRUE(return.frames)) {
    grip.geodesic.misf.collect.level.frames(refinement$level_results)
  } else {
    NULL
  }
  stage.bundle <- grip.geodesic.misf.kk.build.stage.bundle(
    prepared = prepared,
    top.level.fit = top.level.fit,
    top.level.elapsed = top.level.elapsed,
    top.level.frames = top.level.frames,
    insertion = insertion,
    insertion.elapsed = insertion.elapsed,
    insertion.frames = insertion.level.frames,
    refinement = refinement,
    refinement.elapsed = refinement.elapsed,
    refinement.frames = refinement.level.frames,
    final.polish = final.polish,
    final.polish.elapsed = final.polish.elapsed
  )
  stage.trace <- stage.bundle$stage_trace

  trace.detail <- if (isTRUE(return.trace)) {
    list(
      trace_schema_version = grip.geodesic.misf.stage.trace.schema.version(),
      stage_data = stage.bundle$stage_data,
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
  frames <- if (isTRUE(return.frames)) {
    list(
      stage_data = stage.bundle$stage_data,
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

  fit <- grip.new.misf.geodesic.kk.fit(list(
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
  ))
  fit$score <- score.misf.geodesic.kk(
    fit = fit,
    return.trace = isTRUE(return.trace)
  )
  fit
}

#' Score a MISF-based geodesic-KK fit
#'
#' `score.misf.geodesic.kk()` summarizes external coordinates against a
#' MISF-GKK prepared object using either the exact full geodesic-KK scorer or
#' the sparse landmark geodesic-KK scorer, depending on `score.pair.mode`. When
#' a full MISF-GKK fit is supplied, the function also carries over the
#' multiscale stage metadata and any trace tables retained by the optimizer.
#'
#' @param fit Optional fit from [misf.geodesic.kk()].
#' @param coords Optional coordinate matrix used when `fit` is omitted.
#' @param prepared Optional MISF-GKK prepared object used when `fit` is omitted.
#' @param stiffness Global stiffness constant \(K\).
#' @param distance.floor Small positive floor used in
#'   \code{k_ij = K / max(g_ij, distance.floor)^2}.
#' @param edge.length.epsilon Small positive stabilizer added inside each
#'   embedded edge length.
#' @param score.pair.mode Requested scoring pair policy: `"full"`,
#'   `"landmark"`, or `"auto"`.
#' @param score.full.limit Active-set size threshold used when
#'   `score.pair.mode = "auto"`.
#' @param score.local.nbrs Sparse local-neighborhood size used for
#'   landmark-based scoring.
#' @param score.landmark.count Sparse landmark count used for landmark-based
#'   scoring.
#' @param return.trace If `TRUE`, attach any trace tables stored inside `fit`.
#'
#' @return A one-row data frame with final MISF-GKK summary fields and either
#'   full-GKK or landmark-GKK score columns, depending on the resolved scoring
#'   policy. If `return.trace = TRUE`, trace list-columns are attached when
#'   available.
#'
#' @export
#' @md
score.misf.geodesic.kk <- function(fit = NULL,
                                        coords = NULL,
                                        prepared = NULL,
                                        stiffness = 1.0,
                                        distance.floor = 1e-8,
                                        edge.length.epsilon = 1e-8,
                                        score.pair.mode = c("full", "landmark", "auto"),
                                        score.full.limit = 2048L,
                                        score.local.nbrs = 20L,
                                        score.landmark.count = 8L,
                                        return.trace = FALSE) {
  score.pair.mode <- match.arg(score.pair.mode)
  score.full.limit <- grip.validate.misf.count(score.full.limit, "score.full.limit", lower = 1L)
  score.local.nbrs <- grip.validate.count(score.local.nbrs, "score.local.nbrs")
  score.landmark.count <- grip.validate.count(score.landmark.count, "score.landmark.count")
  grip.validate.scalar(stiffness, "stiffness", lower = 0, open.lower = TRUE)
  grip.validate.scalar(distance.floor, "distance.floor", lower = 0, open.lower = TRUE)
  grip.validate.scalar(edge.length.epsilon, "edge.length.epsilon", lower = 0)
  if (!is.logical(return.trace) || length(return.trace) != 1L || is.na(return.trace)) {
    stop("return.trace must be TRUE or FALSE")
  }

  top.level.fit <- NULL
  stage.trace <- data.frame()
  trace.detail <- list()
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
    if (!is.null(fit$trace)) {
      trace.detail <- fit$trace
    }
    if (!is.null(fit$timing)) {
      timing <- fit$timing
    }
  }

  coords <- grip.validate.coords(coords)
  prepared <- grip.validate.misf.geodesic.kk.prepared(prepared, coords = coords)
  score.resolved <- grip.geodesic.misf.kk.resolve.score.prepared(
    prepared = prepared,
    pair.mode = score.pair.mode,
    full.limit = score.full.limit,
    local.nbrs = score.local.nbrs,
    landmark.count = score.landmark.count
  )
  score.prepared <- score.resolved$prepared
  score.mode <- score.resolved$pair_resolution$effective
  final.score <- if (identical(score.mode, "landmark")) {
    score.landmark.geodesic.kk(
      coords = coords,
      prepared = score.prepared,
      stiffness = stiffness,
      distance.floor = distance.floor,
      edge.length.epsilon = edge.length.epsilon
    )
  } else {
    score.geodesic.kk(
      coords = coords,
      prepared = score.prepared,
      stiffness = stiffness,
      distance.floor = distance.floor,
      edge.length.epsilon = edge.length.epsilon
    )
  }
  final.energy <- grip.geodesic.misf.kk.extract.score.metric(
    final.score,
    c("gkk.energy", "lgkk.energy")
  )
  final.weighted.rmse <- grip.geodesic.misf.kk.extract.score.metric(
    final.score,
    c("gkk.weighted.rmse", "lgkk.weighted.rmse")
  )
  final.weighted.rel.rmse <- grip.geodesic.misf.kk.extract.score.metric(
    final.score,
    c("gkk.weighted.rel.rmse", "lgkk.weighted.rel.rmse")
  )
  final.mean.abs.path.error <- grip.geodesic.misf.kk.extract.score.metric(
    final.score,
    c("gkk.mean.abs.path.error", "lgkk.mean.abs.path.error")
  )
  final.mean.rel.path.error <- grip.geodesic.misf.kk.extract.score.metric(
    final.score,
    c("gkk.mean.rel.path.error", "lgkk.mean.rel.path.error")
  )
  final.score.mode <- score.mode
  final.score.requested.mode <- score.resolved$pair_resolution$requested
  final.score.pair.count <- length(score.prepared$pair_graph_distance)
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
    final.score.mode = final.score.mode,
    final.score.requested.mode = final.score.requested.mode,
    final.score.pair.count = as.integer(final.score.pair.count),
    final.energy = final.energy,
    final.weighted.rmse = final.weighted.rmse,
    final.weighted.rel.rmse = final.weighted.rel.rmse,
    final.mean.abs.path.error = final.mean.abs.path.error,
    final.mean.rel.path.error = final.mean.rel.path.error,
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

  if (isTRUE(return.trace)) {
    out$stage.trace <- list(stage.trace)
    out$top.restart.summary <- list(
      if (!is.null(trace.detail$top_restart_summary)) trace.detail$top_restart_summary else data.frame()
    )
    out$insertion.level.trace <- list(
      if (!is.null(trace.detail$insertion_level_trace)) trace.detail$insertion_level_trace else data.frame()
    )
    out$insertion.vertex.trace <- list(
      if (!is.null(trace.detail$insertion_vertex_trace)) trace.detail$insertion_vertex_trace else data.frame()
    )
    out$refinement.level.trace <- list(
      if (!is.null(trace.detail$refinement_level_trace)) trace.detail$refinement_level_trace else data.frame()
    )
    out$final.polish.trace <- list(
      if (!is.null(trace.detail$final_polish_trace)) trace.detail$final_polish_trace else data.frame()
    )
  }
  out
}
