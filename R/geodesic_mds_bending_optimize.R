grip.optimize.geodesic.mds <- function(coords = NULL,
                                       prepared = NULL,
                                       data = NULL,
                                       k = NULL,
                                       dim = 2L,
                                       connect = c("mst", "error"),
                                       tie.mode = c("single", "average"),
                                       init = c("cmdscale", "random", "user"),
                                       anchor.mode = c("none", "cmdscale", "initial", "user"),
                                       anchor.coords = NULL,
                                       anchor.weight = 0,
                                       anchor.weight.end = anchor.weight,
                                       anchor.vertex.weight = NULL,
                                       continuation = c("constant", "linear", "geometric"),
                                       smoothness.weight = 0,
                                       smoothness.weight.end = smoothness.weight,
                                       smoothness.continuation = c("constant", "linear", "geometric"),
                                       edge.spring.weight = 0,
                                       edge.spring.weight.end = edge.spring.weight,
                                       edge.spring.continuation = c("constant", "linear", "geometric"),
                                       repulsion.weight = 0,
                                       repulsion.weight.end = repulsion.weight,
                                       repulsion.continuation = c("constant", "linear", "geometric"),
                                       repulsion.quantile = 0.60,
                                       repulsion.scale = 0.20,
                                       repulsion.cap.quantile = 0.90,
                                       repulsion.hop.min = 3L,
                                       bending.stencils = NULL,
                                       bending.weight = 0,
                                       bending.weight.end = bending.weight,
                                       bending.continuation = c("constant", "linear", "geometric"),
                                       engine = c("cpp", "r"),
                                       max.iter = 16L,
                                       edge.length.epsilon = 1e-8,
                                       initial.step = 1.0,
                                       step.shrink = 0.5,
                                       armijo.factor = 1e-4,
                                       grad.tol = 1e-8,
                                       min.step = 1e-8,
                                       n.threads = 0L,
                                       recenter = TRUE,
                                       return.trace = FALSE,
                                       seed = NULL) {
  if (bending.weight <= 0 && bending.weight.end <= 0 && is.null(bending.stencils)) {
    out <- grip.optimize.geodesic.mds.base(
      coords = coords, prepared = prepared, data = data, k = k, dim = dim,
      connect = connect, tie.mode = tie.mode, init = init, anchor.mode = anchor.mode,
      anchor.coords = anchor.coords, anchor.weight = anchor.weight,
      anchor.weight.end = anchor.weight.end, anchor.vertex.weight = anchor.vertex.weight,
      continuation = continuation,
      smoothness.weight = smoothness.weight, smoothness.weight.end = smoothness.weight.end,
      smoothness.continuation = smoothness.continuation,
      edge.spring.weight = edge.spring.weight,
      edge.spring.weight.end = edge.spring.weight.end,
      edge.spring.continuation = edge.spring.continuation,
      repulsion.weight = repulsion.weight,
      repulsion.weight.end = repulsion.weight.end,
      repulsion.continuation = repulsion.continuation,
      repulsion.quantile = repulsion.quantile,
      repulsion.scale = repulsion.scale,
      repulsion.cap.quantile = repulsion.cap.quantile,
      repulsion.hop.min = repulsion.hop.min,
      engine = engine,
      max.iter = max.iter, edge.length.epsilon = edge.length.epsilon,
      initial.step = initial.step, step.shrink = step.shrink,
      armijo.factor = armijo.factor, grad.tol = grad.tol, min.step = min.step,
      n.threads = n.threads, recenter = recenter, return.trace = return.trace, seed = seed
    )
    out$bending_schedule <- rep.int(0, max.iter + 1L)
    out$bending_stencils <- NULL
    out$final_bending_weight <- 0
    return(out)
  }

  if (smoothness.weight > 0 || smoothness.weight.end > 0) {
    stop("combined smoothness and bending regularization is not implemented in this round")
  }

  init <- match.arg(init)
  tie.mode <- match.arg(tie.mode)
  anchor.mode <- match.arg(anchor.mode)
  continuation <- match.arg(continuation)
  smoothness.continuation <- match.arg(smoothness.continuation)
  edge.spring.continuation <- match.arg(edge.spring.continuation)
  repulsion.continuation <- match.arg(repulsion.continuation)
  bending.continuation <- match.arg(bending.continuation)
  engine <- match.arg(engine)
  grip.validate.scalar(max.iter, "max.iter", lower = 0)
  grip.validate.scalar(edge.length.epsilon, "edge.length.epsilon", lower = 0)
  grip.validate.scalar(initial.step, "initial.step", lower = 0, open.lower = TRUE)
  grip.validate.scalar(step.shrink, "step.shrink", lower = 0, upper = 1, open.lower = TRUE, open.upper = TRUE)
  grip.validate.scalar(armijo.factor, "armijo.factor", lower = 0)
  grip.validate.scalar(grad.tol, "grad.tol", lower = 0)
  grip.validate.scalar(min.step, "min.step", lower = 0, open.lower = TRUE)
  grip.validate.scalar(anchor.weight, "anchor.weight", lower = 0)
  grip.validate.scalar(anchor.weight.end, "anchor.weight.end", lower = 0)
  grip.validate.scalar(smoothness.weight, "smoothness.weight", lower = 0)
  grip.validate.scalar(smoothness.weight.end, "smoothness.weight.end", lower = 0)
  grip.validate.scalar(edge.spring.weight, "edge.spring.weight", lower = 0)
  grip.validate.scalar(edge.spring.weight.end, "edge.spring.weight.end", lower = 0)
  grip.validate.scalar(repulsion.weight, "repulsion.weight", lower = 0)
  grip.validate.scalar(repulsion.weight.end, "repulsion.weight.end", lower = 0)
  grip.validate.scalar(repulsion.quantile, "repulsion.quantile", lower = 0, upper = 1)
  grip.validate.scalar(repulsion.scale, "repulsion.scale", lower = 0)
  grip.validate.scalar(repulsion.cap.quantile, "repulsion.cap.quantile", lower = 0, upper = 1)
  repulsion.hop.min <- grip.validate.count(repulsion.hop.min, "repulsion.hop.min")
  if (repulsion.hop.min < 2L) {
    stop("repulsion.hop.min must be at least 2")
  }
  grip.validate.scalar(bending.weight, "bending.weight", lower = 0)
  grip.validate.scalar(bending.weight.end, "bending.weight.end", lower = 0)
  if (identical(anchor.mode, "none") && (anchor.weight > 0 || anchor.weight.end > 0)) {
    stop("anchor.mode must not be 'none' when anchor.weight or anchor.weight.end is positive")
  }
  if (is.null(prepared)) {
    prepared <- grip.prepare.geodesic.mds(data = data, k = k, connect = connect, tie.mode = tie.mode)
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
  bend.stencils <- grip.validate.bending.stencils(bending.stencils, n = nrow(coords))
  if (is.null(bend.stencils) || nrow(bend.stencils) == 0L) {
    stop("bending.stencils must be provided when bending regularization is used")
  }

  anchor.coords <- grip.geodesic.mds.resolve.anchor(
    anchor.mode = anchor.mode, coords = coords, prepared = prepared,
    anchor.coords = anchor.coords, recenter = recenter
  )
  anchor.schedule <- if (is.null(anchor.coords)) rep.int(0, max.iter + 1L) else
    grip.geodesic.mds.weight.schedule(max.iter, anchor.weight, anchor.weight.end, continuation)
  smoothness.schedule <- grip.geodesic.mds.weight.schedule(
    max.iter = max.iter, weight = smoothness.weight,
    weight.end = smoothness.weight.end, continuation = smoothness.continuation
  )
  edge.spring.schedule <- grip.geodesic.mds.weight.schedule(
    max.iter = max.iter, weight = edge.spring.weight,
    weight.end = edge.spring.weight.end, continuation = edge.spring.continuation
  )
  repulsion.schedule <- grip.geodesic.mds.weight.schedule(
    max.iter = max.iter, weight = repulsion.weight,
    weight.end = repulsion.weight.end, continuation = repulsion.continuation
  )
  bending.schedule <- grip.geodesic.mds.weight.schedule(
    max.iter = max.iter, weight = bending.weight,
    weight.end = bending.weight.end, continuation = bending.continuation
  )
  prepared <- grip.geodesic.mds.ensure.graph.term.cache(
    prepared = prepared,
    repulsion.weight = max(repulsion.schedule),
    repulsion.quantile = repulsion.quantile,
    repulsion.scale = repulsion.scale,
    repulsion.cap.quantile = repulsion.cap.quantile,
    repulsion.hop.min = repulsion.hop.min
  )

  if (identical(engine, "cpp") &&
      (any(edge.spring.schedule > 0) || any(repulsion.schedule > 0))) {
    warning("edge_spring_weight/repulsion_weight are currently implemented only in the R engine; falling back to the R engine")
    engine <- "r"
  }

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
      coords = coords, max_iter = max.iter,
      edge_length_epsilon = edge.length.epsilon,
      initial_step = initial.step, step_shrink = step.shrink,
      armijo_factor = armijo.factor, grad_tol = grad.tol, min_step = min.step,
      recenter = recenter, return_trace = return.trace,
      anchor_coords = anchor.coords, anchor_weights = anchor.schedule,
      bend_a = bend.flat$flat_bend_a, bend_b = bend.flat$flat_bend_b,
      bend_c = bend.flat$flat_bend_c, bend_weights = bending.schedule
    )
    opt$final_smoothness_weight <- 0
    opt$final_edge_spring_weight <- 0
    opt$final_repulsion_weight <- 0
  } else {
    opt <- grip.optimize.geodesic.mds.bending.r(
      coords = coords, prepared = prepared, anchor.coords = anchor.coords,
      anchor.schedule = anchor.schedule, bending.stencils = bend.stencils,
      bending.schedule = bending.schedule, edge.length.epsilon = edge.length.epsilon,
      smoothness.schedule = smoothness.schedule,
      edge.spring.schedule = edge.spring.schedule,
      repulsion.schedule = repulsion.schedule,
      repulsion.quantile = repulsion.quantile,
      repulsion.scale = repulsion.scale,
      repulsion.cap.quantile = repulsion.cap.quantile,
      repulsion.hop.min = repulsion.hop.min,
      max.iter = max.iter, initial.step = initial.step, step.shrink = step.shrink,
      armijo.factor = armijo.factor, grad.tol = grad.tol, min.step = min.step,
      recenter = recenter, return.trace = return.trace
    )
    opt$final_anchor_weight <- utils::tail(anchor.schedule, 1L)
    opt$final_smoothness_weight <- utils::tail(smoothness.schedule, 1L)
    opt$final_edge_spring_weight <- utils::tail(edge.spring.schedule, 1L)
    opt$final_repulsion_weight <- utils::tail(repulsion.schedule, 1L)
    opt$final_bending_weight <- utils::tail(bending.schedule, 1L)
    opt$n_threads_used <- 1L
  }

  score <- grip.score.geodesic.mds(
    coords = opt$coords, prepared = prepared,
    edge.length.epsilon = edge.length.epsilon,
    anchor.coords = anchor.coords, anchor.weight = opt$final_anchor_weight,
    anchor.vertex.weight = anchor.vertex.weight,
    smoothness.weight = opt$final_smoothness_weight,
    edge.spring.weight = opt$final_edge_spring_weight,
    repulsion.weight = opt$final_repulsion_weight,
    repulsion.quantile = repulsion.quantile,
    repulsion.scale = repulsion.scale,
    repulsion.cap.quantile = repulsion.cap.quantile,
    repulsion.hop.min = repulsion.hop.min,
    bending.stencils = bend.stencils,
    bending.weight = opt$final_bending_weight
  )
  list(
    coords = opt$coords, trace = opt$trace, frames = opt$frames,
    prepared = prepared, score = score, anchor_coords = anchor.coords,
    anchor_schedule = anchor.schedule, smoothness_schedule = smoothness.schedule,
    edge_spring_schedule = edge.spring.schedule,
    repulsion_schedule = repulsion.schedule,
    bending_schedule = bending.schedule, bending_stencils = bend.stencils,
    final_anchor_weight = opt$final_anchor_weight,
    final_smoothness_weight = opt$final_smoothness_weight,
    final_edge_spring_weight = opt$final_edge_spring_weight,
    final_repulsion_weight = opt$final_repulsion_weight,
    final_bending_weight = opt$final_bending_weight,
    n_threads_used = if (!is.null(opt$n_threads_used)) opt$n_threads_used else 1L
  )
}
