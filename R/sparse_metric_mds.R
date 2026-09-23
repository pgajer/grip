# Private sparse implementation; public dispatch and documentation: metric.mds().
.sparse.metric.mds <- function(edges = NULL, n = NULL, adj_list = NULL,
                              weight_list = NULL, edge_weights = NULL,
                              dim = 2L, n_pivots = 200L, pivots = NULL,
                              init = "random", max_iter = 30L, seed = 1L,
                              sgd_control = list()) {
  grip.validate.graph.arguments(edges, n, adj_list, weight_list, edge_weights)
  n <- grip.resolve.graph.n(n, edges, adj_list)
  if (n < 2L) stop("Sparse SGD needs at least two vertices", call. = FALSE)
  if (!is.numeric(dim) || length(dim) != 1L || !is.finite(dim) || !dim %in% c(2,3))
    stop("dim must be 2 or 3", call. = FALSE)
  dim <- as.integer(dim)
  supplied.n_pivots <- !missing(n_pivots)
  n_pivots <- min(n, grip.validate.vertex.count(n_pivots))
  if (!is.null(pivots)) {
    grip.validate.vertex.ids(pivots, "pivots", n)
    if (!length(pivots) || anyDuplicated(pivots) || !is.null(base::dim(pivots)))
      stop("pivots must be a nonempty vector of distinct vertex ids", call. = FALSE)
    if (supplied.n_pivots && n_pivots != length(pivots))
      stop("n_pivots must match the supplied pivots", call. = FALSE)
    n_pivots <- length(pivots)
  }
  max_iter <- grip.validate.vertex.count(max_iter)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
       seed != trunc(seed) || abs(seed) > .Machine$integer.max))
    stop("seed must be an integer or NULL", call. = FALSE)
  if (is.list(sgd_control) && is.null(sgd_control[["max_workspace_bytes"]]))
    sgd_control$max_workspace_bytes <- 512 * 1024^2
  control <- grip.mds.sgd.control(sgd_control, max_iter)
  rates <- grip.mds.sgd.rates(control, max_iter)
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
  prepared <- prepare.edge.kk(edges, n, adj_list, weight_list, edge_weights)
  if (prepared$n_components != 1L) stop("Sparse SGD requires a connected graph", call. = FALSE)
  sparse <- grip_sparse_prepare_cpp(n, prepared$edges, prepared$edge_targets,
    n_pivots, as.integer(pivots), preparation.seed, control$max_workspace_bytes)
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
    control$checkpoint_every, control$max_workspace_bytes, TRUE,
    wi, sparse$pairs, wj, FALSE)
  coords <- sweep(fit$terminal_conf,2L,colMeans(fit$terminal_conf),"-") * rms
  if (any(!is.finite(coords))) stop("Sparse coordinates exceed the numeric range", call. = FALSE)
  proxy <- utils::tail(fit$trace$raw_stress,1L)
  energy <- sum((sparse$count_i + sparse$count_j)/2)
  trace <- fit$trace[, c("epoch","raw_stress","pair_updates","learning_rate")]
  names(trace)[2L] <- "sparse_proxy_stress"
  gmds.result(coords, "metric_mds", prepared=prepared, trace=trace,
    metadata=list(engine="sgd", approximation="sparse", backend_version="grip-sparse-sgd-v1",
      objective="asymmetric sparse constraints (endpoint-average proxy reported separately)",
      pair_weights="inverse_squared", sparse_proxy_stress=unname(proxy),
      sparse_proxy_normalized_rmse=sqrt(unname(proxy)/energy),
      input_rms_distance=rms, sparse=sparse, pair_updates=fit$pair_updates,
      preparation_seconds=preparation.seconds,
      fitting_seconds=proc.time()[["elapsed"]]-began,
      fitting_workspace_bytes=fit$workspace_bytes,
      preparation_seed=preparation.seed, fitting_seed=fitting.seed,
      termination="iteration_limit", converged=FALSE,
      settings=list(dim=dim,n_pivots=n_pivots,init=initialization,
                    max_iter=max_iter,seed=seed,sgd_control=control)))
}
