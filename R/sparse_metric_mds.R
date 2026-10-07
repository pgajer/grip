# Private sparse implementation; public dispatch and documentation: metric.mds().
.sparse.metric.mds <- function(edges = NULL, n = NULL, adj.list = NULL,
                              weight.list = NULL, edge.weights = NULL,
                              dim = 2L, n.pivots = 200L, pivots = NULL,
                              init = "random", max.iter = 30L, seed = 1L,
                              sgd.control = list(), constraints = NULL) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights)
  n <- grip.resolve.graph.n(n, edges, adj.list)
  if (n < 2L) stop("Sparse SGD needs at least two vertices", call. = FALSE)
  if (!is.numeric(dim) || length(dim) != 1L || !is.finite(dim) || dim < 2 || dim != trunc(dim) || dim > .Machine$integer.max)
    stop("dim must be an integer at least 2", call. = FALSE)
  dim <- as.integer(dim)
  supplied.n_pivots <- !missing(n.pivots)
  n.pivots <- min(n, grip.validate.vertex.count(n.pivots))
  if (!is.null(pivots)) {
    grip.validate.vertex.ids(pivots, "pivots", n)
    if (!length(pivots) || anyDuplicated(pivots) || !is.null(base::dim(pivots)))
      stop("pivots must be a nonempty vector of distinct vertex ids", call. = FALSE)
    if (supplied.n_pivots && n.pivots != length(pivots))
      stop("n.pivots must match the supplied pivots", call. = FALSE)
    n.pivots <- length(pivots)
  }
  max.iter <- grip.validate.vertex.count(max.iter)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
       seed != trunc(seed) || abs(seed) > .Machine$integer.max))
    stop("seed must be an integer or NULL", call. = FALSE)
  if (is.list(sgd.control) && is.null(sgd.control[["max.workspace.bytes"]]))
    sgd.control$max.workspace.bytes <- 512 * 1024^2
  control <- grip.mds.sgd.control(sgd.control, max.iter)
  if (!is.null(seed)) {
    had.seed <- exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
    old.seed <- if (had.seed) get(".Random.seed", envir=.GlobalEnv) else NULL
    on.exit({
      if (had.seed) assign(".Random.seed", old.seed, envir=.GlobalEnv)
      else if (exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))
        rm(".Random.seed", envir=.GlobalEnv)
    }, add=TRUE)
    set.seed(seed)
  }
  preparation.seed <- sample.int(.Machine$integer.max,1L)-1L
  fitting.seed <- sample.int(.Machine$integer.max,1L)-1L
  began <- proc.time()[["elapsed"]]
  if (is.null(constraints)) {
    prepared <- prepare.edge.kk(edges, n, adj.list, weight.list, edge.weights)
    if (prepared$n_components != 1L) stop("Sparse SGD requires a connected graph", call. = FALSE)
    sparse <- grip_sparse_prepare_cpp(n, prepared$edges, prepared$edge_targets,
      n.pivots, as.integer(pivots), preparation.seed, control$max.workspace.bytes)
  } else {
    sparse <- .metric.mds.validate.constraints(constraints, n, control$max.workspace.bytes)
    prepared <- grip_constraint_graph_cpp(n,sparse$pairs,sparse$targets)
    if (prepared$n_components != 1L) stop("Sparse constraints must form a connected network", call. = FALSE)
    prepared <- c(prepared,list(n=n,edges=sparse$pairs,edge_targets=sparse$targets,
      pair_matrix=matrix(integer(),ncol=2L),pair_graph_distance=numeric(),
      path_vertices=list(),path_edges=list(),path_edge_weights=list(),pair_path_count_log=numeric(),
      graph_diameter=NA_real_,distance_matrix=NULL,pair_mode="edge_only",graph_build_mode="distance_constraints"))
    class(prepared)<-c("grip_edge_kk_prepared","grip_gmds_prepared","grip_gkk_prepared","grip_geodesic_kk_prepared","list")
  }
  preparation.seconds <- proc.time()[["elapsed"]] - began
  began <- proc.time()[["elapsed"]]
  target <- sparse$targets
  maximum <- max(target)
  rms <- maximum * sqrt(mean((target / maximum)^2))
  targets <- target / rms
  wi <- sparse$count_i * (1 / targets)^2
  wj <- sparse$count_j * (1 / targets)^2
  if (any(!is.finite(wi)) || any(!is.finite(wj)) ||
      any(wi[sparse$count_i > 0] <= 0) || any(wj[sparse$count_j > 0] <= 0))
    stop("Sparse weights are outside the representable numeric range", call. = FALSE)
  rates <- grip.mds.sgd.rates(control,max.iter,c(wi,wj))
  if (is.matrix(init)) {
    if (!is.numeric(init) || !identical(base::dim(init), c(n,dim)) || any(!is.finite(init)))
      stop("init must be a finite n by dim matrix", call. = FALSE)
    start <- sweep(init, 2L, colMeans(init), "-") / rms
    initialization <- "supplied"
  } else {
    if (!identical(init,"random")) stop("init must be 'random' or a coordinate matrix", call. = FALSE)
    start <- matrix(stats::runif(n*dim),n,dim)
    initialization <- "random"
  }
  fit <- grip_sgd_mds_cpp(start, targets, rates, fitting.seed,
    control$checkpoint.every, control$max.workspace.bytes, TRUE,
    wi, sparse$pairs, wj, FALSE)
  coords <- sweep(fit$terminal_conf,2L,colMeans(fit$terminal_conf),"-") * rms
  if (any(!is.finite(coords))) stop("Sparse coordinates exceed the numeric range", call. = FALSE)
  proxy <- utils::tail(fit$trace$raw_stress,1L)
  energy <- sum((sparse$count_i + sparse$count_j)/2)
  trace <- fit$trace[, c("epoch","raw_stress","pair_updates","learning_rate")]
  names(trace)[2L] <- "sparse_proxy_stress"
  gmds.result(coords, "metric_mds", prepared=prepared, trace=trace,
    metadata=list(engine="sgd", approximation="sparse", backend_version="grip-sparse-sgd-v1",
      constraint_source=if (is.null(constraints)) "graph_pivots" else "explicit_distances",
      objective=if (all(sparse$count_i == sparse$count_j)) "symmetric sparse inverse-squared stress" else
        "asymmetric sparse constraints (endpoint-average proxy reported separately)",
      pair_weights="inverse_squared", sparse_proxy_stress=unname(proxy),
      sparse_proxy_normalized_rmse=sqrt(unname(proxy)/energy),
      input_rms_distance=rms, sparse=sparse, pair_updates=fit$pair_updates,
      preparation_seconds=preparation.seconds,
      fitting_seconds=proc.time()[["elapsed"]]-began,
      fitting_workspace_bytes=fit$workspace_bytes,
      preparation_seed=preparation.seed, fitting_seed=fitting.seed,
      termination="iteration_limit", converged=FALSE,
      settings=list(dim=dim,n.pivots=if (is.null(constraints)) n.pivots else NULL,init=initialization,
                    max_iter=max.iter,seed=seed,sgd_control=control)))
}

# Canonicalize explicit constraints without interpreting them as graph lengths.
.metric.mds.validate.constraints <- function(x, n, budget) {
  allowed <- c("pairs", "targets", "count_i", "count_j")
  if (!is.list(x) || is.null(names(x)) || anyNA(names(x)) ||
      anyDuplicated(names(x)) || any(!names(x) %in% allowed) ||
      !all(c("pairs", "targets") %in% names(x)))
    stop("constraints must contain pairs, targets and optional count_i/count_j", call. = FALSE)
  pairs <- x$pairs
  if (!is.matrix(pairs) || !is.numeric(pairs) || ncol(pairs) != 2L || !nrow(pairs))
    stop("constraints pairs must be a nonempty two-column numeric matrix", call. = FALSE)
  grip.validate.vertex.ids(pairs, "constraints pairs", n)
  m <- nrow(pairs)
  estimate <- 256 * as.double(m) + 128 * as.double(n)
  if (estimate > budget)
    stop("Explicit constraint preparation exceeds max.workspace.bytes", call. = FALSE)
  if (any(pairs[,1] == pairs[,2])) stop("constraints cannot contain self-pairs", call. = FALSE)
  targets <- x$targets
  if (!is.numeric(targets) || !is.null(dim(targets)) || length(targets) != m ||
      any(!is.finite(targets)) || any(targets <= 0))
    stop("constraints targets must be positive finite distances, one per pair", call. = FALSE)
  counts <- lapply(c("count_i", "count_j"), function(key) {
    v <- x[[key]]
    if (is.null(v)) return(rep(1, m))
    if (!is.numeric(v) || !is.null(dim(v)) || length(v) != m ||
        any(!is.finite(v)) || any(v < 0))
      stop("constraints endpoint multiplicities must be finite and nonnegative", call. = FALSE)
    v
  })
  a <- counts[[1]]; b <- counts[[2]]
  if (any(!is.finite(a+b)) || any(a+b <= 0) || !is.finite(sum(a+b)))
    stop("constraints require positive finite total multiplicity per pair", call. = FALSE)
  reverse <- pairs[,1] > pairs[,2]
  pairs[reverse,] <- pairs[reverse,2:1,drop=FALSE]
  tmp <- a[reverse]; a[reverse] <- b[reverse]; b[reverse] <- tmp
  order <- order(pairs[,1],pairs[,2])
  a.ids <- pairs[order,1]; b.ids <- pairs[order,2]
  if(m>1L && any(a.ids[-1L]==a.ids[-m] & b.ids[-1L]==b.ids[-m]))
    stop("constraints pairs must be unique", call. = FALSE)
  storage.mode(pairs) <- "integer"
  list(pairs=pairs[order,,drop=FALSE], targets=targets[order],
       count_i=a[order], count_j=b[order], pivots=integer(), region=integer(),
       workspace_estimate_bytes=estimate)
}
