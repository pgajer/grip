#' Sparse stochastic-gradient graph layout
#'
#' `sparse.metric.mds()` uses graph edges and a small set of representative
#' vertices (pivots) to approximate long-range stress interactions. It implements
#' the asymmetric sparse SGD updates of Zheng et al. (2018), based on the sparse
#' model of Ortmann et al. (2017), in 2D or 3D. Distances are prepared only
#' from pivots, without a separate all-pairs matrix.
#'
#' @param edges Two-column matrix of undirected edges, using 1-based vertex ids.
#' @param n Number of vertices; inferred from edges or adjacency when omitted.
#' @param adj_list Reciprocal adjacency list, instead of `edges`.
#' @param weight_list Positive lengths parallel to `adj_list`.
#' @param edge_weights Positive finite lengths parallel to `edges`; omitted
#'   lengths are one.
#' @param dim Output dimension, 2 or 3.
#' @param n_pivots Number of pivots, capped at `n`. Default 200 is a provisional
#'   calibration choice. Increasing it increases memory and work per epoch.
#' @param pivots Optional distinct vertex ids in selection order. Overrides
#'   the default `n_pivots`; an explicitly supplied inconsistent count is an error.
#' @param init `"random"` for a uniform unit cube in normalized-distance units,
#'   or a finite `n` by `dim` coordinate matrix in input-distance units.
#' @param max_iter Number of epochs (passes through the sparse constraints).
#'   Default 30 is a provisional calibration choice, not a convergence criterion.
#' @param seed Optional integer seed, preserving the caller's R RNG state.
#' @param sgd_control Controls as for [metric.mds()]: scheduler, learning_rate,
#'   final_rate, switch_ratio, checkpoint_every, and max_workspace_bytes.
#'   The sparse default native allowance is 512 MiB. Separate preparation and
#'   fitting estimates are checked; this is not a cap on total R process memory.
#'
#' @details The graph must be connected. Disconnected graphs are rejected without
#' repair. Pivots are chosen with probability proportional to distance to the
#' nearest selected pivot, starting from a uniform random vertex. Graph lengths
#' are normalized by their maximum before shortest-path calculations. Exact distance
#' ties in region assignment go to the earliest selected pivot. Near floating-point
#' boundaries, general weighted graphs can still be sensitive to arithmetic. Each nonadjacent
#' vertex receives an influence from a pivot weighted by the number of vertices
#' in that pivot's region within half their separation, divided by squared
#' separation. The reverse influence can be zero or different. Edges have
#' inverse-squared weights at both endpoints.
#'
#' Edge targets are the supplied edge lengths, as in the upstream sparse
#' implementation. For a weighted edge with a shorter alternative route, that
#' target differs from the shortest-path target of [metric.mds()]. Exact
#' full-pivot update parity requires every edge to be a shortest path, as with
#' unit lengths or Euclidean chord lengths, and matched starts, rates and order.
#'
#' Targets are normalized by RMS over the retained unordered pairs; endpoint
#' weights use the same units. The default hybrid schedule is grip's existing
#' schedule, not the paper's weight-dependent annealing. The terminal iterate is
#' returned, centered and converted back to input units, with no fitted scale.
#' A checkpoint proxy averages the two endpoint weights for each pair. Asymmetric
#' updates are not claimed to follow the gradient of that proxy, or to decrease
#' it monotonically. This diagnostic is not full stress or equation (19) of
#' Zheng et al. Use independent distance evaluation to compare sparse and full
#' layouts. No dense diagnostics or classical initialization run implicitly.
#'
#' @return A `grip_gmds_layout` with `method = "sparse_metric_mds"`.
#'   `coords` is the terminal layout. `prepared` contains an edge-only graph,
#'   without a distance matrix. `metadata$sparse` records retained pairs,
#'   original-unit targets, endpoint multiplicities, selected pivots and region
#'   owners. Metadata also records seeds, controls, preparation/fitting seconds,
#'   native workspace estimates, `sparse_proxy_stress`, its normalized root
#'   error, and `termination = "iteration_limit"`, `converged = FALSE`.
#'   `trace` reports checkpoint proxy stress, epochs, rates and pair updates.
#' @references
#' Ortmann, M., Klimenta, M., and Brandes, U. (2017). A Sparse Stress Model.
#' Journal of Graph Algorithms and Applications, 21(5), 791--821.
#' \doi{10.7155/jgaa.00440}.
#'
#' Zheng, J. X., Pawar, S., and Goodman, D. F. M. (2018).
#' Graph Drawing by Stochastic Gradient Descent.
#' \doi{10.1109/TVCG.2018.2859997}, Algorithm 2.
#' @seealso [metric.mds()], [grip()], [edge.kk()], [layout.coords()]
#' @examples
#' fit <- sparse.metric.mds(edges.path(30), n = 30, n_pivots = 5,
#'                          dim = 3, max_iter = 20, seed = 1)
#' plot.layout(layout.coords(fit), edges = edges.path(30))
#' fit$metadata$sparse_proxy_stress
#' @export
sparse.metric.mds <- function(edges = NULL, n = NULL, adj_list = NULL,
                              weight_list = NULL, edge_weights = NULL,
                              dim = 3L, n_pivots = 200L, pivots = NULL,
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
  gmds.result(coords, "sparse_metric_mds", prepared=prepared, trace=trace,
    metadata=list(engine="sgd", backend_version="grip-sparse-sgd-v1",
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
