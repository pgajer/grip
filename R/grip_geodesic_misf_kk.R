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

grip.geodesic.misf.kk.score.kind <- function(prepared) {
  pair.mode <- if (!is.null(prepared$pair_mode)) prepared$pair_mode else "all_pairs"
  if (pair.mode %in% c("landmark_sparse", "misf_sparse")) "landmark" else "full"
}

grip.geodesic.misf.kk.add.prepared.class <- function(prepared,
                                                     effective_pair_mode = NULL) {
  if (is.null(effective_pair_mode)) {
    effective_pair_mode <- grip.geodesic.misf.kk.score.kind(prepared)
  }
  extra.class <- if (identical(effective_pair_mode, "landmark")) {
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
                                                 distance_floor = 1e-8,
                                                 edge_length_epsilon = 1e-8) {
  effective.pair.mode <- grip.geodesic.misf.kk.score.kind(prepared)
  prepared <- grip.geodesic.misf.kk.add.prepared.class(prepared, effective.pair.mode)
  if (identical(effective.pair.mode, "landmark")) {
    grip.score.landmark.geodesic.kk(
      coords = coords,
      prepared = prepared,
      stiffness = stiffness,
      distance_floor = distance_floor,
      edge_length_epsilon = edge_length_epsilon
    )
  } else {
    grip.score.geodesic.kk(
      coords = coords,
      prepared = prepared,
      stiffness = stiffness,
      distance_floor = distance_floor,
      edge_length_epsilon = edge_length_epsilon,
      scale_mode = "profiled"
    )
  }
}

grip.geodesic.misf.kk.resolve.score.prepared <- function(prepared,
                                                         pair_mode = c("full", "landmark", "auto"),
                                                         full_limit = 2048L,
                                                         local_nbrs = 20L,
                                                         landmark_count = 8L) {
  prepared <- grip.validate.misf.geodesic.kk.prepared(prepared)
  pair_mode <- match.arg(pair_mode)
  resolution <- grip.geodesic.misf.kk.resolve.pair.mode(
    pair_mode = pair_mode,
    active_n = prepared$n,
    full_limit = full_limit
  )
  score.prepared <- if (identical(resolution$effective, "full")) {
    grip.geodesic.misf.kk.add.prepared.class(prepared, "full")
  } else {
    grip.geodesic.misf.kk.add.prepared.class(
      grip.prepare.landmark.geodesic.kk(
        n = prepared$n,
        adj_list = prepared$adj_list,
        weight_list = prepared$weight_list,
        local_nbrs = local_nbrs,
        landmark_count = landmark_count
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
                                                 distance_floor = 1e-8,
                                                 edge_length_epsilon = 1e-8,
                                                 scale_mode = c("profiled", "fixed"),
                                                 scale.L0 = NULL,
                                                 anchor_coords = NULL,
                                                 anchor_weight = 0,
                                                 anchor_vertex_weight = NULL) {
  scale_mode <- match.arg(scale_mode)
  anchor.coords <- if (is.null(anchor_coords)) {
    NULL
  } else {
    grip.geodesic.mds.resolve.anchor(
      anchor_mode = "user",
      coords = coords,
      prepared = prepared,
      anchor_coords = anchor_coords,
      recenter = FALSE
    )
  }
  kk.state <- grip.geodesic.kk.evaluate.state(
    coords = coords,
    prepared = prepared,
    stiffness = stiffness,
    distance_floor = distance_floor,
    edge_length_epsilon = edge_length_epsilon,
    scale_mode = scale_mode,
    scale.L0 = scale.L0
  )
  anchor.stats <- grip.geodesic.mds.anchor.stats(
    coords = coords,
    anchor_coords = anchor.coords,
    anchor_weight = anchor_weight,
    anchor_vertex_weight = anchor_vertex_weight
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
                                                    max_iter = 8L,
                                                    stiffness = 1.0,
                                                    distance_floor = 1e-8,
                                                    edge_length_epsilon = 1e-8,
                                                    initial_step = 1.0,
                                                    step_shrink = 0.5,
                                                    armijo_factor = 1e-4,
                                                    grad_tol = 1e-8,
                                                    min_step = 1e-8,
                                                    anchor_coords = NULL,
                                                    anchor_weight = 0,
                                                    anchor_weight_end = anchor_weight,
                                                    anchor_vertex_weight = NULL,
                                                    continuation = c("constant", "linear", "geometric"),
                                                    recenter = TRUE,
                                                    return_trace = FALSE) {
  continuation <- match.arg(continuation)
  effective.pair.mode <- grip.geodesic.misf.kk.score.kind(prepared)
  prepared <- grip.geodesic.misf.kk.add.prepared.class(prepared, effective.pair.mode)
  grip.validate.scalar(anchor_weight, "anchor_weight", lower = 0)
  grip.validate.scalar(anchor_weight_end, "anchor_weight_end", lower = 0)
  use.anchor <- !is.null(anchor_coords) &&
    is.finite(anchor_weight) &&
    is.finite(anchor_weight_end) &&
    max(anchor_weight, anchor_weight_end) > 0
  if (!use.anchor) {
    anchor_weight <- 0
    anchor_weight_end <- 0
  }

  if (!use.anchor) {
    fit <- if (identical(effective.pair.mode, "landmark")) {
      grip.optimize.landmark.geodesic.kk(
        coords = coords,
        prepared = prepared,
        max_iter = max_iter,
        stiffness = stiffness,
        distance_floor = distance_floor,
        edge_length_epsilon = edge_length_epsilon,
        initial_step = initial_step,
        step_shrink = step_shrink,
        armijo_factor = armijo_factor,
        grad_tol = grad_tol,
        min_step = min_step,
        recenter = recenter,
        return_trace = return_trace
      )
    } else {
      grip.optimize.geodesic.kk(
        coords = coords,
        prepared = prepared,
        max_iter = max_iter,
        stiffness = stiffness,
        distance_floor = distance_floor,
        edge_length_epsilon = edge_length_epsilon,
        initial_step = initial_step,
        step_shrink = step_shrink,
        armijo_factor = armijo_factor,
        grad_tol = grad_tol,
        min_step = min_step,
        recenter = recenter,
        return_trace = return_trace,
        scale_mode = "profiled"
      )
    }
    fit$effective_pair_mode <- effective.pair.mode
    return(fit)
  }

  anchor.coords <- grip.geodesic.mds.resolve.anchor(
    anchor_mode = "user",
    coords = coords,
    prepared = prepared,
    anchor_coords = anchor_coords,
    recenter = FALSE
  )
  anchor.vertex.weight <- grip.geodesic.mds.resolve.anchor.vertex.weight(
    anchor_vertex_weight = anchor_vertex_weight,
    coords = coords
  )
  anchor.schedule <- grip.geodesic.mds.anchor.schedule(
    max_iter = max_iter,
    anchor_weight = anchor_weight,
    anchor_weight_end = anchor_weight_end,
    continuation = continuation
  )
  fixed.scale.L0 <- NULL
  scale.mode <- if (identical(effective.pair.mode, "landmark")) "fixed" else "profiled"
  if (identical(scale.mode, "fixed") && length(prepared$pair_graph_distance) > 0L) {
    fixed.scale.L0 <- grip.geodesic.kk.fit.scale(
      path.lengths = grip.geodesic.kk.path.lengths(
        coords = coords,
        prepared = prepared,
        edge_length_epsilon = edge_length_epsilon
      ),
      graph.distances = prepared$pair_graph_distance,
      stiffness = stiffness,
      distance_floor = distance_floor
    )
    if (!is.finite(fixed.scale.L0)) {
      stop("failed to fit an initial LGKK scale")
    }
  }

  if (nrow(coords) <= 1L || length(prepared$pair_graph_distance) == 0L || max_iter == 0L) {
    score <- grip.geodesic.misf.kk.score.prepared(
      coords = coords,
      prepared = prepared,
      stiffness = stiffness,
      distance_floor = distance_floor,
      edge_length_epsilon = edge_length_epsilon
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
  trace.rows <- vector("list", max_iter + 1L)
  accepted.frames <- list(current)
  state <- grip.geodesic.misf.kk.evaluate.state(
    coords = current,
    prepared = prepared,
    stiffness = stiffness,
    distance_floor = distance_floor,
    edge_length_epsilon = edge_length_epsilon,
    scale_mode = scale.mode,
    scale.L0 = fixed.scale.L0,
    anchor_coords = anchor.coords,
    anchor_weight = anchor.schedule[[1L]],
    anchor_vertex_weight = anchor.vertex.weight
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

  for (iter in seq_len(max_iter)) {
    iter.anchor.weight <- anchor.schedule[[iter + 1L]]
    state <- grip.geodesic.misf.kk.evaluate.state(
      coords = current,
      prepared = prepared,
      stiffness = stiffness,
      distance_floor = distance_floor,
      edge_length_epsilon = edge_length_epsilon,
      scale_mode = scale.mode,
      scale.L0 = fixed.scale.L0,
      anchor_coords = anchor.coords,
      anchor_weight = iter.anchor.weight,
      anchor_vertex_weight = anchor.vertex.weight
    )
    if (!is.finite(state$gradient_norm) || state$gradient_norm <= grad_tol) {
      break
    }
    step <- as.double(initial_step)
    accepted <- FALSE
    candidate <- current
    candidate.state <- state

    while (is.finite(step) && step >= min_step) {
      proposal <- current - step * state$gradient
      if (isTRUE(recenter)) {
        proposal <- sweep(proposal, 2L, colMeans(proposal), "-", check.margin = FALSE)
      }
      proposal.state <- grip.geodesic.misf.kk.evaluate.state(
        coords = proposal,
        prepared = prepared,
        stiffness = stiffness,
        distance_floor = distance_floor,
        edge_length_epsilon = edge_length_epsilon,
        scale_mode = scale.mode,
        scale.L0 = fixed.scale.L0,
        anchor_coords = anchor.coords,
        anchor_weight = iter.anchor.weight,
        anchor_vertex_weight = anchor.vertex.weight
      )
      target.energy <- state$energy - armijo_factor * step * state$gradient_norm^2
      if (is.finite(proposal.state$energy) && proposal.state$energy <= target.energy) {
        candidate <- proposal
        candidate.state <- proposal.state
        accepted <- TRUE
        break
      }
      step <- step * step_shrink
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
  if (!isTRUE(return_trace)) {
    trace.df <- trace.df[, c("iteration", "energy", "kk_energy", "anchor_energy", "gradient_norm", "step", "accepted", "scale.L0", "anchor_weight"), drop = FALSE]
    accepted.frames <- list(current)
  }
  score <- grip.geodesic.misf.kk.score.prepared(
    coords = current,
    prepared = prepared,
    stiffness = stiffness,
    distance_floor = distance_floor,
    edge_length_epsilon = edge_length_epsilon
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
                                                             pair_mode = c("auto", "full", "landmark"),
                                                             full_limit = 512L,
                                                             local_nbrs = 20L,
                                                             landmark_count = 8L) {
  prepared <- grip.validate.misf.geodesic.kk.prepared(prepared)
  pair_mode <- match.arg(pair_mode)
  resolution <- grip.geodesic.misf.kk.resolve.pair.mode(
    pair_mode = pair_mode,
    active_n = prepared$top_level_graph$n,
    full_limit = full_limit
  )
  coarse.prepared <- if (identical(resolution$effective, "full")) {
    if (identical(local_nbrs, prepared$top_level_local_nbrs) &&
        identical(landmark_count, prepared$top_level_landmark_count) &&
        !is.null(prepared$top_level_prepared_full)) {
      prepared$top_level_prepared_full
    } else {
      grip.prepare.geodesic.kk(
        edges = prepared$top_level_graph$edges,
        n = prepared$top_level_graph$n,
        edge_weights = prepared$top_level_graph$edge_weights,
        tie_mode = prepared$tie_mode
      )
    }
  } else {
    if (identical(local_nbrs, prepared$top_level_local_nbrs) &&
        identical(landmark_count, prepared$top_level_landmark_count) &&
        !is.null(prepared$top_level_prepared_sparse)) {
      prepared$top_level_prepared_sparse
    } else {
      grip.prepare.landmark.geodesic.kk(
        edges = prepared$top_level_graph$edges,
        n = prepared$top_level_graph$n,
        edge_weights = prepared$top_level_graph$edge_weights,
        local_nbrs = local_nbrs,
        landmark_count = landmark_count
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
                                                  pair_mode = c("auto", "full", "landmark"),
                                                  full_limit = 512L,
                                                  local_nbrs = 20L,
                                                  landmark_count = 8L,
                                                  n_restarts = 8L,
                                                  max_iter = 16L,
                                                  init = c("cmdscale", "random"),
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
                                                  seed = 6L) {
  prepared <- grip.validate.misf.geodesic.kk.prepared(prepared)
  pair_mode <- match.arg(pair_mode)
  init <- match.arg(init)
  dim <- grip.validate.count(dim, "dim")
  if (!(dim %in% c(2L, 3L))) {
    stop("dim must be 2 or 3")
  }
  n_restarts <- grip.validate.misf.count(n_restarts, "n_restarts", lower = 1L)
  grip.validate.scalar(max_iter, "max_iter", lower = 0)
  max_iter <- as.integer(round(max_iter))

  resolved <- grip.geodesic.misf.kk.resolve.top.level.prepared(
    prepared = prepared,
    pair_mode = pair_mode,
    full_limit = full_limit,
    local_nbrs = local_nbrs,
    landmark_count = landmark_count
  )
  coarse.prepared <- resolved$prepared
  vertex.ids <- resolved$vertex_ids

  if (coarse.prepared$n <= 1L) {
    coords <- matrix(0, nrow = coarse.prepared$n, ncol = dim)
    score <- grip.geodesic.misf.kk.score.prepared(
      coords = coords,
      prepared = coarse.prepared,
      stiffness = stiffness,
      distance_floor = distance_floor,
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

  restart.rows <- vector("list", n_restarts)
  best.fit <- NULL
  best.row <- NULL
  best.restart <- 1L
  best.energy <- Inf

  for (restart in seq_len(n_restarts)) {
    restart.seed <- if (is.null(seed)) NULL else as.integer(seed + restart - 1L)
    restart.init <- if (identical(init, "cmdscale") && restart == 1L) "cmdscale" else "random"
    init.coords <- grip.geodesic.misf.kk.init.coords(
      prepared = coarse.prepared,
      dim = dim,
      init = restart.init,
      restart = restart,
      seed = restart.seed
    )
    init.score <- grip.geodesic.misf.kk.score.prepared(
      coords = init.coords,
      prepared = coarse.prepared,
      stiffness = stiffness,
      distance_floor = distance_floor,
      edge_length_epsilon = edge_length_epsilon
    )
    fit <- grip.geodesic.misf.kk.optimize.prepared(
      coords = init.coords,
      prepared = coarse.prepared,
      max_iter = max_iter,
      stiffness = stiffness,
      distance_floor = distance_floor,
      edge_length_epsilon = edge_length_epsilon,
      initial_step = initial_step,
      step_shrink = step_shrink,
      armijo_factor = armijo_factor,
      grad_tol = grad_tol,
      min_step = min_step,
      recenter = recenter,
      return_trace = return_trace
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
  best.fit
}

grip.geodesic.misf.kk.prepare.active.level <- function(prepared,
                                                       level = NULL,
                                                       local_nbrs = 8L,
                                                       landmark_count = 4L,
                                                       pair_mode = c("auto", "full", "landmark"),
                                                       full_limit = 256L) {
  prepared <- grip.validate.misf.geodesic.kk.prepared(prepared)
  pair_mode <- match.arg(pair_mode)
  if (is.null(level)) {
    level <- prepared$top_level_level
  }
  level <- as.integer(level)
  active.vertices <- grip.geodesic.misf.active.level.vertices(prepared, level)
  resolution <- grip.geodesic.misf.kk.resolve.pair.mode(
    pair_mode = pair_mode,
    active_n = length(active.vertices),
    full_limit = full_limit
  )
  active.prepared <- grip.geodesic.misf.prepare.active.level(
    prepared = prepared,
    active_vertices = active.vertices,
    local_nbrs = local_nbrs,
    landmark_count = landmark_count,
    pair_mode = if (identical(resolution$effective, "full")) "full" else "sparse"
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
                                               local_nbrs = 8L,
                                               landmark_count = 4L,
                                               pair_mode = c("auto", "full", "landmark"),
                                               full_limit = 256L,
                                               max_iter = 8L,
                                               anchor_weight = 0.05,
                                               anchor_weight_end = anchor_weight,
                                               continuation = c("constant", "linear", "geometric"),
                                               stiffness = 1.0,
                                               distance_floor = 1e-8,
                                               edge_length_epsilon = 1e-8,
                                               initial_step = 1.0,
                                               step_shrink = 0.5,
                                               armijo_factor = 1e-4,
                                               grad_tol = 1e-8,
                                               min_step = 1e-8,
                                               recenter = FALSE,
                                               return_trace = TRUE) {
  prepared <- grip.validate.misf.geodesic.kk.prepared(prepared)
  coords <- grip.validate.coords(coords)
  if (nrow(coords) != prepared$n) {
    stop("nrow(coords) must match prepared$n")
  }
  continuation <- match.arg(continuation)
  built <- grip.geodesic.misf.kk.prepare.active.level(
    prepared = prepared,
    level = level,
    local_nbrs = local_nbrs,
    landmark_count = landmark_count,
    pair_mode = pair_mode,
    full_limit = full_limit
  )
  before <- grip.geodesic.misf.kk.score.prepared(
    coords = coords,
    prepared = built$active_prepared,
    stiffness = stiffness,
    distance_floor = distance_floor,
    edge_length_epsilon = edge_length_epsilon
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
  use.anchor <- any(anchor.vertex.weight > 0) && is.finite(anchor_weight_end) && anchor_weight_end > 0
  fit <- grip.geodesic.misf.kk.optimize.prepared(
    coords = coords,
    prepared = built$active_prepared,
    max_iter = max_iter,
    anchor_coords = if (use.anchor) anchor.coords else NULL,
    anchor_weight = if (use.anchor) anchor_weight else 0,
    anchor_weight_end = if (use.anchor) anchor_weight_end else 0,
    anchor_vertex_weight = if (use.anchor) anchor.vertex.weight else NULL,
    continuation = continuation,
    stiffness = stiffness,
    distance_floor = distance_floor,
    edge_length_epsilon = edge_length_epsilon,
    initial_step = initial_step,
    step_shrink = step_shrink,
    armijo_factor = armijo_factor,
    grad_tol = grad_tol,
    min_step = min_step,
    recenter = recenter,
    return_trace = return_trace
  )
  after <- grip.geodesic.misf.kk.score.prepared(
    coords = fit$coords,
    prepared = built$active_prepared,
    stiffness = stiffness,
    distance_floor = distance_floor,
    edge_length_epsilon = edge_length_epsilon
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
    anchor_weight = if (use.anchor) as.double(anchor_weight) else 0,
    anchor_weight_end = if (use.anchor) as.double(anchor_weight_end) else 0,
    anchor_continuation = continuation,
    coords = fit$coords,
    fit = fit,
    before = before,
    after = after
  )
}

grip.geodesic.misf.kk.refine.all.levels <- function(prepared,
                                                    coords,
                                                    local_nbrs = 8L,
                                                    landmark_count = 4L,
                                                    pair_mode = c("auto", "full", "landmark"),
                                                    full_limit = 256L,
                                                    max_iter = 8L,
                                                    anchor_weight = 0.05,
                                                    anchor_weight_end = anchor_weight,
                                                    continuation = c("constant", "linear", "geometric"),
                                                    stiffness = 1.0,
                                                    distance_floor = 1e-8,
                                                    edge_length_epsilon = 1e-8,
                                                    initial_step = 1.0,
                                                    step_shrink = 0.5,
                                                    armijo_factor = 1e-4,
                                                    grad_tol = 1e-8,
                                                    min_step = 1e-8,
                                                    recenter = FALSE,
                                                    return_trace = TRUE) {
  prepared <- grip.validate.misf.geodesic.kk.prepared(prepared)
  coords <- grip.validate.coords(coords)
  pair_mode <- match.arg(pair_mode)
  continuation <- match.arg(continuation)

  level.ids <- seq.int(from = prepared$top_level_level, to = 0L, by = -1L)
  level.results <- vector("list", length(level.ids))
  for (i in seq_along(level.ids)) {
    level.results[[i]] <- grip.geodesic.misf.kk.refine.level(
      prepared = prepared,
      coords = coords,
      level = level.ids[[i]],
      local_nbrs = local_nbrs,
      landmark_count = landmark_count,
      pair_mode = pair_mode,
      full_limit = full_limit,
      max_iter = max_iter,
      anchor_weight = anchor_weight,
      anchor_weight_end = anchor_weight_end,
      continuation = continuation,
      stiffness = stiffness,
      distance_floor = distance_floor,
      edge_length_epsilon = edge_length_epsilon,
      initial_step = initial_step,
      step_shrink = step_shrink,
      armijo_factor = armijo_factor,
      grad_tol = grad_tol,
      min_step = min_step,
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
                                                         pair_mode = c("auto", "full", "landmark"),
                                                         full_limit = 1024L,
                                                         local_nbrs = 20L,
                                                         landmark_count = 8L) {
  prepared <- grip.validate.misf.geodesic.kk.prepared(prepared)
  pair_mode <- match.arg(pair_mode)
  resolution <- grip.geodesic.misf.kk.resolve.pair.mode(
    pair_mode = pair_mode,
    active_n = prepared$n,
    full_limit = full_limit
  )
  final.prepared <- if (identical(resolution$effective, "full")) {
    grip.geodesic.misf.kk.add.prepared.class(prepared, "full")
  } else {
    grip.geodesic.misf.kk.add.prepared.class(
      grip.prepare.landmark.geodesic.kk(
        n = prepared$n,
        adj_list = prepared$adj_list,
        weight_list = prepared$weight_list,
        local_nbrs = local_nbrs,
        landmark_count = landmark_count
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
                                               pair_mode = c("auto", "full", "landmark"),
                                               full_limit = 1024L,
                                               local_nbrs = 20L,
                                               landmark_count = 8L,
                                               max_iter = 8L,
                                               stiffness = 1.0,
                                               distance_floor = 1e-8,
                                               edge_length_epsilon = 1e-8,
                                               initial_step = 1.0,
                                               step_shrink = 0.5,
                                               armijo_factor = 1e-4,
                                               grad_tol = 1e-8,
                                               min_step = 1e-8,
                                               recenter = TRUE,
                                               return_trace = TRUE) {
  prepared <- grip.validate.misf.geodesic.kk.prepared(prepared)
  coords <- grip.validate.coords(coords)
  if (nrow(coords) != prepared$n) {
    stop("nrow(coords) must match prepared$n")
  }
  resolved <- grip.geodesic.misf.kk.resolve.final.prepared(
    prepared = prepared,
    pair_mode = pair_mode,
    full_limit = full_limit,
    local_nbrs = local_nbrs,
    landmark_count = landmark_count
  )
  fit <- grip.geodesic.misf.kk.optimize.prepared(
    coords = coords,
    prepared = resolved$prepared,
    max_iter = max_iter,
    stiffness = stiffness,
    distance_floor = distance_floor,
    edge_length_epsilon = edge_length_epsilon,
    initial_step = initial_step,
    step_shrink = step_shrink,
    armijo_factor = armijo_factor,
    grad_tol = grad_tol,
    min_step = min_step,
    recenter = recenter,
    return_trace = return_trace
  )
  fit$requested_pair_mode <- resolved$pair_resolution$requested
  fit$effective_pair_mode <- resolved$pair_resolution$effective
  fit
}

grip.geodesic.misf.kk.build.stage.trace <- function(prepared,
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
    pair.mode = if (!is.null(top_level_fit$effective_pair_mode)) top_level_fit$effective_pair_mode else NA_character_,
    energy = grip.geodesic.misf.kk.extract.score.metric(
      top_level_fit$score,
      c("gkk.energy", "lgkk.energy")
    ),
    weighted.rel.rmse = grip.geodesic.misf.kk.extract.score.metric(
      top_level_fit$score,
      c("gkk.weighted.rel.rmse", "lgkk.weighted.rel.rmse")
    ),
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
      pair.mode = NA_character_,
      energy = insertion$level_trace$mean_objective,
      weighted.rel.rmse = NA_real_,
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
      pair.mode = refinement$level_trace$pair.mode,
      energy = refinement$level_trace$after.energy,
      weighted.rel.rmse = refinement$level_trace$after.weighted.rel.rmse,
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
    pair_n = length(final_polish$prepared$pair_graph_distance),
    pair.mode = if (!is.null(final_polish$effective_pair_mode)) final_polish$effective_pair_mode else NA_character_,
    energy = grip.geodesic.misf.kk.extract.score.metric(
      final_polish$score,
      c("gkk.energy", "lgkk.energy")
    ),
    weighted.rel.rmse = grip.geodesic.misf.kk.extract.score.metric(
      final_polish$score,
      c("gkk.weighted.rel.rmse", "lgkk.weighted.rel.rmse")
    ),
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

#' Prepare a MISF-based multiscale geodesic-KK object
#'
#' `grip.prepare.misf.geodesic.kk()` builds the prepared state for a MISF-based
#' geodesic-KK pipeline by layering the maximal independent set filtration
#' (MISF) on top of the existing full geodesic-KK cache. The prepared object
#' stores the graph metadata, the top-level coarse graph, the exact and sparse
#' top-level prepared caches, and the default controls used by the multiscale
#' MISF-GKK optimizer.
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
#' @param dim Target embedding dimension (`2` or `3`) for the multiscale solve.
#' @param top_level_mode Either `"solve"` or `"skip"`. When set to `"solve"`,
#'   the top MISF level is optimized immediately and stored in
#'   `prepared$top_level_fit`.
#' @param top_level_pair_mode Pair policy for the top MISF level: `"auto"`,
#'   `"full"`, or `"landmark"`.
#' @param top_level_full_limit Active-set size threshold used when
#'   `top_level_pair_mode = "auto"`.
#' @param top_level_local_nbrs Number of nearest graph-metric neighbors retained
#'   per vertex in the sparse top-level LGKK cache.
#' @param top_level_landmark_count Number of farthest-point landmarks retained
#'   per vertex in the sparse top-level LGKK cache.
#' @param top_level_restarts Number of top-level restarts used by the coarse
#'   MISF-GKK solve.
#' @param top_level_max_iter Top-level iteration budget.
#' @param top_level_init Top-level initializer (`"cmdscale"` or `"random"`).
#' @param seed Optional integer seed reused for MISF extraction.
#'
#' @return An object of class `"grip_misf_gkk_prepared"` layered on top of the
#'   full geodesic-KK prepared structure, with added MISF metadata, top-level
#'   coarse-graph caches, and stored default controls for the MISF-GKK
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

  if (identical(top_level_mode, "solve")) {
    prepared$top_level_fit <- grip.geodesic.misf.kk.solve.top.level(
      prepared = prepared,
      dim = dim,
      pair_mode = top_level_pair_mode,
      full_limit = top_level_full_limit,
      local_nbrs = top_level_local_nbrs,
      landmark_count = top_level_landmark_count,
      n_restarts = top_level_restarts,
      max_iter = top_level_max_iter,
      init = top_level_init,
      stiffness = prepared$stiffness,
      distance_floor = prepared$distance_floor,
      edge_length_epsilon = prepared$edge_length_epsilon,
      initial_step = prepared$initial_step,
      step_shrink = prepared$step_shrink,
      armijo_factor = prepared$armijo_factor,
      grad_tol = prepared$grad_tol,
      min_step = prepared$min_step,
      recenter = prepared$recenter,
      return_trace = FALSE,
      seed = seed
    )
  }
  prepared
}

#' Optimize an embedding with the MISF-based geodesic-KK pipeline
#'
#' `grip.optimize.misf.geodesic.kk()` runs a multiscale MISF-based geodesic-KK
#' pipeline. It solves the top MISF graph with either full GKK or sparse LGKK,
#' inserts lower-level vertices with the existing MISF placement helpers, refines
#' each active MISF level under a GKK/LGKK objective, and finishes with a final
#' full-graph polish.
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
#' @param refinement_anchor_weight Optional initial anchor weight used to pin
#'   inactive and coarser-level vertices during KK refinement.
#' @param refinement_anchor_weight_end Optional final anchor weight used at the
#'   end of the active-level continuation schedule.
#' @param refinement_continuation Optional continuation schedule used for the
#'   active-level anchor penalty.
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
#' @param return_trace If `TRUE`, retain detailed stage traces.
#' @param return_frames If `TRUE`, retain intermediate coordinate frames.
#'
#' @return A list with `coords`, `prepared`, the per-stage multiscale results,
#'   the stage trace, optional trace/frame details, timing, and the final MISF
#'   score summary.
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
                                           refinement_anchor_weight = NULL,
                                           refinement_anchor_weight_end = NULL,
                                           refinement_continuation = NULL,
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

  top.level.pair.mode <- if (is.null(top_level_pair_mode) && !is.null(prepared) && inherits(prepared, "grip_misf_gkk_prepared")) {
    prepared$top_level_pair_mode
  } else if (is.null(top_level_pair_mode)) {
    "auto"
  } else {
    match.arg(top_level_pair_mode, c("auto", "full", "landmark"))
  }
  top.level.full.limit <- if (is.null(top_level_full_limit) && !is.null(prepared) && inherits(prepared, "grip_misf_gkk_prepared")) {
    prepared$top_level_full_limit
  } else if (is.null(top_level_full_limit)) {
    512L
  } else {
    grip.validate.misf.count(top_level_full_limit, "top_level_full_limit", lower = 1L)
  }
  top.level.local.nbrs <- if (is.null(top_level_local_nbrs) && !is.null(prepared) && inherits(prepared, "grip_misf_gkk_prepared")) {
    prepared$top_level_local_nbrs
  } else if (is.null(top_level_local_nbrs)) {
    20L
  } else {
    grip.validate.count(top_level_local_nbrs, "top_level_local_nbrs")
  }
  top.level.landmark.count <- if (is.null(top_level_landmark_count) && !is.null(prepared) && inherits(prepared, "grip_misf_gkk_prepared")) {
    prepared$top_level_landmark_count
  } else if (is.null(top_level_landmark_count)) {
    8L
  } else {
    grip.validate.count(top_level_landmark_count, "top_level_landmark_count")
  }
  top.level.restarts <- if (is.null(top_level_restarts) && !is.null(prepared) && inherits(prepared, "grip_misf_gkk_prepared")) {
    prepared$top_level_restarts
  } else if (is.null(top_level_restarts)) {
    8L
  } else {
    grip.validate.misf.count(top_level_restarts, "top_level_restarts", lower = 1L)
  }
  top.level.max.iter <- if (is.null(top_level_max_iter) && !is.null(prepared) && inherits(prepared, "grip_misf_gkk_prepared")) {
    prepared$top_level_max_iter
  } else if (is.null(top_level_max_iter)) {
    16L
  } else {
    grip.validate.scalar(top_level_max_iter, "top_level_max_iter", lower = 0)
    as.integer(round(top_level_max_iter))
  }
  top.level.init <- if (is.null(top_level_init) && !is.null(prepared) && inherits(prepared, "grip_misf_gkk_prepared")) {
    prepared$top_level_init
  } else if (is.null(top_level_init)) {
    "cmdscale"
  } else {
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
    dim = dim.resolved,
    top_level_pair_mode = top.level.pair.mode,
    top_level_full_limit = top.level.full.limit,
    top_level_local_nbrs = top.level.local.nbrs,
    top_level_landmark_count = top.level.landmark.count,
    top_level_restarts = top.level.restarts,
    top_level_max_iter = top.level.max.iter,
    top_level_init = top.level.init,
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
  insertion.mode.resolved <- if (is.null(insertion_mode)) {
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

  refinement.pair.mode <- if (is.null(refinement_pair_mode)) {
    prepared$refinement_pair_mode
  } else {
    match.arg(refinement_pair_mode, c("auto", "full", "landmark"))
  }
  refinement.full.limit <- if (is.null(refinement_full_limit)) {
    prepared$refinement_full_limit
  } else {
    grip.validate.misf.count(refinement_full_limit, "refinement_full_limit", lower = 1L)
  }
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
  refinement.continuation <- if (is.null(refinement_continuation)) {
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

  final.pair.mode <- if (is.null(final_pair_mode)) {
    prepared$final_pair_mode
  } else {
    match.arg(final_pair_mode, c("auto", "full", "landmark"))
  }
  final.full.limit <- if (is.null(final_full_limit)) {
    prepared$final_full_limit
  } else {
    grip.validate.misf.count(final_full_limit, "final_full_limit", lower = 1L)
  }
  final.local.nbrs <- if (is.null(final_local_nbrs)) {
    prepared$final_local_nbrs
  } else {
    grip.validate.count(final_local_nbrs, "final_local_nbrs")
  }
  final.landmark.count <- if (is.null(final_landmark_count)) {
    prepared$final_landmark_count
  } else {
    grip.validate.count(final_landmark_count, "final_landmark_count")
  }
  final.max.iter <- if (is.null(final_max_iter)) {
    prepared$final_max_iter
  } else {
    grip.validate.scalar(final_max_iter, "final_max_iter", lower = 0)
    as.integer(round(final_max_iter))
  }

  grip.validate.scalar(stiffness, "stiffness", lower = 0, open.lower = TRUE)
  grip.validate.scalar(distance_floor, "distance_floor", lower = 0, open.lower = TRUE)
  grip.validate.scalar(edge_length_epsilon, "edge_length_epsilon", lower = 0)
  grip.validate.scalar(initial_step, "initial_step", lower = 0, open.lower = TRUE)
  grip.validate.scalar(step_shrink, "step_shrink", lower = 0, upper = 1, open.lower = TRUE, open.upper = TRUE)
  grip.validate.scalar(armijo_factor, "armijo_factor", lower = 0)
  grip.validate.scalar(grad_tol, "grad_tol", lower = 0)
  grip.validate.scalar(min_step, "min_step", lower = 0, open.lower = TRUE)
  if (!is.logical(recenter) || length(recenter) != 1L || is.na(recenter)) {
    stop("recenter must be TRUE or FALSE")
  }
  if (!is.logical(return_trace) || length(return_trace) != 1L || is.na(return_trace)) {
    stop("return_trace must be TRUE or FALSE")
  }
  if (!is.logical(return_frames) || length(return_frames) != 1L || is.na(return_frames)) {
    stop("return_frames must be TRUE or FALSE")
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
  prepared$distance_floor <- distance_floor
  prepared$edge_length_epsilon <- edge_length_epsilon
  prepared$initial_step <- initial_step
  prepared$step_shrink <- step_shrink
  prepared$armijo_factor <- armijo_factor
  prepared$grad_tol <- grad_tol
  prepared$min_step <- min_step
  prepared$recenter <- recenter

  need.top.trace <- isTRUE(return_trace) || isTRUE(return_frames)
  top.override.requested <- any(!vapply(
    list(
      top_level_pair_mode,
      top_level_full_limit,
      top_level_local_nbrs,
      top_level_landmark_count,
      top_level_restarts,
      top_level_max_iter,
      top_level_init,
      dim
    ),
    is.null,
    logical(1L)
  ))
  need.top.solve <- is.null(prepared$top_level_fit) || isTRUE(top.override.requested)
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
    prepared$top_level_fit <- grip.geodesic.misf.kk.solve.top.level(
      prepared = prepared,
      dim = dim.resolved,
      pair_mode = top.level.pair.mode,
      full_limit = top.level.full.limit,
      local_nbrs = top.level.local.nbrs,
      landmark_count = top.level.landmark.count,
      n_restarts = top.level.restarts,
      max_iter = top.level.max.iter,
      init = top.level.init,
      stiffness = stiffness,
      distance_floor = distance_floor,
      edge_length_epsilon = edge_length_epsilon,
      initial_step = initial_step,
      step_shrink = step_shrink,
      armijo_factor = armijo_factor,
      grad_tol = grad_tol,
      min_step = min_step,
      recenter = recenter,
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
    insertion <- if (identical(insertion.mode.resolved, "geodesic")) {
      grip.geodesic.misf.insert.all.levels(
        prepared = prepared,
        coords = coords,
        anchor_policy = insertion.anchor.policy,
        anchor_count = insertion.anchor.count,
        anchor_weight_mode = insertion.anchor.weight.mode,
        max_iter = insertion.max.iter,
        initial_step = initial_step,
        step_shrink = step_shrink,
        armijo_factor = armijo_factor,
        grad_tol = grad_tol,
        min_step = min_step
      )
    } else {
      grip.geodesic.misf.insert.all.levels.with.layout(
        prepared = prepared,
        coords = coords,
        method = insertion.mode.resolved,
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
  refinement <- grip.geodesic.misf.kk.refine.all.levels(
    prepared = prepared,
    coords = insertion$coords,
    local_nbrs = refinement.local.nbrs,
    landmark_count = refinement.landmark.count,
    pair_mode = refinement.pair.mode,
    full_limit = refinement.full.limit,
    max_iter = refinement.max.iter,
    anchor_weight = refinement.anchor.weight,
    anchor_weight_end = refinement.anchor.weight.end,
    continuation = refinement.continuation,
    stiffness = stiffness,
    distance_floor = distance_floor,
    edge_length_epsilon = edge_length_epsilon,
    initial_step = initial_step,
    step_shrink = step_shrink,
    armijo_factor = armijo_factor,
    grad_tol = grad_tol,
    min_step = min_step,
    recenter = FALSE,
    return_trace = isTRUE(return_trace) || isTRUE(return_frames)
  )
  refinement.elapsed <- proc.time()[["elapsed"]] - refinement.start

  final.polish.start <- proc.time()[["elapsed"]]
  final.polish <- grip.geodesic.misf.kk.final.polish(
    prepared = prepared,
    coords = refinement$coords,
    pair_mode = final.pair.mode,
    full_limit = final.full.limit,
    local_nbrs = final.local.nbrs,
    landmark_count = final.landmark.count,
    max_iter = final.max.iter,
    stiffness = stiffness,
    distance_floor = distance_floor,
    edge_length_epsilon = edge_length_epsilon,
    initial_step = initial_step,
    step_shrink = step_shrink,
    armijo_factor = armijo_factor,
    grad_tol = grad_tol,
    min_step = min_step,
    recenter = recenter,
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
  stage.trace <- grip.geodesic.misf.kk.build.stage.trace(
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
  fit$score <- grip.score.misf.geodesic.kk(
    fit = fit,
    return_trace = isTRUE(return_trace)
  )
  fit
}

#' Score a MISF-based geodesic-KK fit
#'
#' `grip.score.misf.geodesic.kk()` summarizes external coordinates against a
#' MISF-GKK prepared object using either the exact full geodesic-KK scorer or
#' the sparse landmark geodesic-KK scorer, depending on `score_pair_mode`. When
#' a full MISF-GKK fit is supplied, the function also carries over the
#' multiscale stage metadata and any trace tables retained by the optimizer.
#'
#' @param fit Optional fit from [grip.optimize.misf.geodesic.kk()].
#' @param coords Optional coordinate matrix used when `fit` is omitted.
#' @param prepared Optional MISF-GKK prepared object used when `fit` is omitted.
#' @param stiffness Global stiffness constant \(K\).
#' @param distance_floor Small positive floor used in
#'   \code{k_ij = K / max(g_ij, distance_floor)^2}.
#' @param edge_length_epsilon Small positive stabilizer added inside each
#'   embedded edge length.
#' @param score_pair_mode Requested scoring pair policy: `"full"`,
#'   `"landmark"`, or `"auto"`.
#' @param score_full_limit Active-set size threshold used when
#'   `score_pair_mode = "auto"`.
#' @param score_local_nbrs Planned sparse local-neighborhood size for future
#'   LGKK-based scoring.
#' @param score_landmark_count Planned sparse landmark count for future
#'   LGKK-based scoring.
#' @param return_trace If `TRUE`, attach any trace tables stored inside `fit`.
#'
#' @return A one-row data frame with final MISF-GKK summary fields and either
#'   full-GKK or landmark-GKK score columns, depending on the resolved scoring
#'   policy. If `return_trace = TRUE`, trace list-columns are attached when
#'   available.
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
    pair_mode = score_pair_mode,
    full_limit = score_full_limit,
    local_nbrs = score_local_nbrs,
    landmark_count = score_landmark_count
  )
  score.prepared <- score.resolved$prepared
  score.mode <- score.resolved$pair_resolution$effective
  final.score <- if (identical(score.mode, "landmark")) {
    grip.score.landmark.geodesic.kk(
      coords = coords,
      prepared = score.prepared,
      stiffness = stiffness,
      distance_floor = distance_floor,
      edge_length_epsilon = edge_length_epsilon
    )
  } else {
    grip.score.geodesic.kk(
      coords = coords,
      prepared = score.prepared,
      stiffness = stiffness,
      distance_floor = distance_floor,
      edge_length_epsilon = edge_length_epsilon
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

  if (isTRUE(return_trace)) {
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
