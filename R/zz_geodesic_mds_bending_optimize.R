grip.optimize.geodesic.mds <- function(coords = NULL,
                                       prepared = NULL,
                                       data = NULL,
                                       k = NULL,
                                       dim = 2L,
                                       connect = c("mst", "error"),
                                       tie_mode = c("single", "average"),
                                       init = c("cmdscale", "random", "user"),
                                       anchor_mode = c("none", "cmdscale", "initial", "user"),
                                       anchor_coords = NULL,
                                       anchor_weight = 0,
                                       anchor_weight_end = anchor_weight,
                                       continuation = c("constant", "linear", "geometric"),
                                       smoothness_weight = 0,
                                       smoothness_weight_end = smoothness_weight,
                                       smoothness_continuation = c("constant", "linear", "geometric"),
                                       bending_stencils = NULL,
                                       bending_weight = 0,
                                       bending_weight_end = bending_weight,
                                       bending_continuation = c("constant", "linear", "geometric"),
                                       engine = c("cpp", "r"),
                                       max_iter = 16L,
                                       edge_length_epsilon = 1e-8,
                                       initial_step = 1.0,
                                       step_shrink = 0.5,
                                       armijo_factor = 1e-4,
                                       grad_tol = 1e-8,
                                       min_step = 1e-8,
                                       n_threads = 0L,
                                       recenter = TRUE,
                                       return_trace = FALSE,
                                       seed = NULL) {
  if (bending_weight <= 0 && bending_weight_end <= 0 && is.null(bending_stencils)) {
    out <- grip.optimize.geodesic.mds.base(
      coords = coords, prepared = prepared, data = data, k = k, dim = dim,
      connect = connect, tie_mode = tie_mode, init = init, anchor_mode = anchor_mode,
      anchor_coords = anchor_coords, anchor_weight = anchor_weight,
      anchor_weight_end = anchor_weight_end, continuation = continuation,
      smoothness_weight = smoothness_weight, smoothness_weight_end = smoothness_weight_end,
      smoothness_continuation = smoothness_continuation, engine = engine,
      max_iter = max_iter, edge_length_epsilon = edge_length_epsilon,
      initial_step = initial_step, step_shrink = step_shrink,
      armijo_factor = armijo_factor, grad_tol = grad_tol, min_step = min_step,
      n_threads = n_threads, recenter = recenter, return_trace = return_trace, seed = seed
    )
    out$bending_schedule <- rep.int(0, max_iter + 1L)
    out$bending_stencils <- NULL
    out$final_bending_weight <- 0
    return(out)
  }

  if (smoothness_weight > 0 || smoothness_weight_end > 0) {
    stop("combined smoothness and bending regularization is not implemented in this round")
  }

  init <- match.arg(init)
  tie_mode <- match.arg(tie_mode)
  anchor_mode <- match.arg(anchor_mode)
  continuation <- match.arg(continuation)
  bending_continuation <- match.arg(bending_continuation)
  engine <- match.arg(engine)
  grip.validate.scalar(max_iter, "max_iter", lower = 0)
  grip.validate.scalar(edge_length_epsilon, "edge_length_epsilon", lower = 0)
  grip.validate.scalar(initial_step, "initial_step", lower = 0, open.lower = TRUE)
  grip.validate.scalar(step_shrink, "step_shrink", lower = 0, upper = 1, open.lower = TRUE, open.upper = TRUE)
  grip.validate.scalar(armijo_factor, "armijo_factor", lower = 0)
  grip.validate.scalar(grad_tol, "grad_tol", lower = 0)
  grip.validate.scalar(min_step, "min_step", lower = 0, open.lower = TRUE)
  grip.validate.scalar(anchor_weight, "anchor_weight", lower = 0)
  grip.validate.scalar(anchor_weight_end, "anchor_weight_end", lower = 0)
  grip.validate.scalar(bending_weight, "bending_weight", lower = 0)
  grip.validate.scalar(bending_weight_end, "bending_weight_end", lower = 0)
  if (identical(anchor_mode, "none") && (anchor_weight > 0 || anchor_weight_end > 0)) {
    stop("anchor_mode must not be 'none' when anchor_weight or anchor_weight_end is positive")
  }
  if (is.null(prepared)) {
    prepared <- grip.prepare.geodesic.mds(data = data, k = k, connect = connect, tie_mode = tie_mode)
  }
  if (is.null(coords)) {
    dim <- grip.validate.count(dim, "dim")
    if (identical(init, "user")) stop("coords must be supplied when init = 'user'")
    coords <- if (identical(init, "cmdscale")) {
      grip.geodesic.mds.cmdscale.init(prepared, dim)
    } else {
      if (!is.null(seed)) set.seed(as.integer(seed))
      matrix(stats::rnorm(prepared$n * dim), ncol = dim)
    }
  } else {
    coords <- grip.validate.coords(coords)
  }
  prepared <- grip.validate.geodesic.mds.prepared(prepared, coords = coords)
  bend.stencils <- grip.validate.bending.stencils(bending_stencils, n = nrow(coords))
  if (is.null(bend.stencils) || nrow(bend.stencils) == 0L) {
    stop("bending_stencils must be provided when bending regularization is used")
  }

  anchor.coords <- grip.geodesic.mds.resolve.anchor(
    anchor_mode = anchor_mode, coords = coords, prepared = prepared,
    anchor_coords = anchor_coords, recenter = recenter
  )
  anchor.schedule <- if (is.null(anchor.coords)) rep.int(0, max_iter + 1L) else
    grip.geodesic.mds.weight.schedule(max_iter, anchor_weight, anchor_weight_end, continuation)
  bending.schedule <- grip.geodesic.mds.weight.schedule(
    max_iter = max_iter, weight = bending_weight,
    weight_end = bending_weight_end, continuation = bending_continuation
  )

  if (identical(engine, "cpp") &&
      !is.null(prepared$flat_pair_edge_offsets) &&
      !is.null(prepared$flat_edge_u) &&
      !is.null(prepared$flat_edge_v) &&
      !is.null(prepared$flat_edge_coeff)) {
    bend.flat <- grip.flatten.bending.stencils.zero.based(bend.stencils)
    opt <- grip_optimize_geodesic_mds_flat_bending_cpp(
      flat_pair_edge_offsets = prepared$flat_pair_edge_offsets,
      flat_edge_u = prepared$flat_edge_u,
      flat_edge_v = prepared$flat_edge_v,
      flat_edge_coeff = prepared$flat_edge_coeff,
      pair_graph_distance = prepared$pair_graph_distance,
      coords = coords, max_iter = max_iter,
      edge_length_epsilon = edge_length_epsilon,
      initial_step = initial_step, step_shrink = step_shrink,
      armijo_factor = armijo_factor, grad_tol = grad_tol, min_step = min_step,
      recenter = recenter, return_trace = return_trace,
      anchor_coords = anchor.coords, anchor_weights = anchor.schedule,
      bend_a = bend.flat$flat_bend_a, bend_b = bend.flat$flat_bend_b,
      bend_c = bend.flat$flat_bend_c, bend_weights = bending.schedule
    )
  } else {
    opt <- grip.optimize.geodesic.mds.bending.r(
      coords = coords, prepared = prepared, anchor.coords = anchor.coords,
      anchor.schedule = anchor.schedule, bending.stencils = bend.stencils,
      bending.schedule = bending.schedule, edge_length_epsilon = edge_length_epsilon,
      max_iter = max_iter, initial_step = initial_step, step_shrink = step_shrink,
      armijo_factor = armijo_factor, grad_tol = grad_tol, min_step = min_step,
      recenter = recenter, return_trace = return_trace
    )
    opt$final_anchor_weight <- utils::tail(anchor.schedule, 1L)
    opt$final_bending_weight <- utils::tail(bending.schedule, 1L)
    opt$n_threads_used <- 1L
  }

  score <- grip.score.geodesic.mds(
    coords = opt$coords, prepared = prepared,
    edge_length_epsilon = edge_length_epsilon,
    anchor_coords = anchor.coords, anchor_weight = opt$final_anchor_weight,
    smoothness_weight = 0, bending_stencils = bend.stencils,
    bending_weight = opt$final_bending_weight
  )
  list(
    coords = opt$coords, trace = opt$trace, frames = opt$frames,
    prepared = prepared, score = score, anchor_coords = anchor.coords,
    anchor_schedule = anchor.schedule, smoothness_schedule = rep.int(0, max_iter + 1L),
    bending_schedule = bending.schedule, bending_stencils = bend.stencils,
    final_anchor_weight = opt$final_anchor_weight,
    final_smoothness_weight = 0, final_bending_weight = opt$final_bending_weight,
    n_threads_used = if (!is.null(opt$n_threads_used)) opt$n_threads_used else 1L
  )
}
