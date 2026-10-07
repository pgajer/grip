#' Prepare sparse metric MDS directly from observations
#'
#' Selects pivots and computes fitting weights internally, without constructing
#' a full distance matrix. Pass the result to [metric.mds()] with
#' `approximation = "sparse"`. Preparation can be reused across dimensions and
#' fit seeds. The default fit schedule for this preparation is `"zheng"`.
#'
#' @param X Finite numeric observation matrix with distinct rows.
#' @param method `"euclidean"` uses regional Euclidean sparse stress;
#'   `"uniform"` uses uniformly selected landmarks without regional multipliers
#'   or local edges; `"geodesic"` uses graph shortest-path targets and regions.
#' @param n.pivots Positive pivot count, capped at `nrow(X) - 1` to retain a
#'   strictly sparse preparation. All-pivot reduction is tested internally.
#' @param pivot.selection Distance-proportional `"randomized"` selection or
#'   `"farthest"` (Euclidean regional method only). Uniform mode selects
#'   uniformly without replacement regardless of the omitted default.
#' @param q Number of non-self neighbors for the symmetric union graph.
#' @param neighbor.variance.target Fraction of PCA variance retained for neighbor
#'   search, default 0.9; `NULL` searches original coordinates. Scores are not
#'   whitened. Only search changes: edge lengths and bridges use original `X`.
#'   A target of one retains all numerically resolved components.
#' @param seed Nonnegative integer preparation seed. The caller's RNG is untouched.
#' @param max.workspace.bytes Conservative preparation allocation allowance;
#'   this is not a total process memory cap.
#' @return A `grip_sparse_mds_preparation` list containing constraints, pivot
#'   distance rows, pivots, regions, graph (when used), PCA information and
#'   metadata. Input coordinates and settings are retained for cache validation.
#' @details The geodesic method uses exact component-MST repair in `dgraphs`.
#'   Its streaming implementation can evaluate quadratically many distances but
#'   never retains their full matrix. A recent local `dgraphs` implementation
#'   supporting `search.data` and streaming component repair is required.
#'   Uniform and Euclidean modes compute one Euclidean distance row per pivot.
#'   All modes use inverse-squared target weights. Regional counts are asymmetric
#'   endpoint multipliers, not ordinary symmetric stress weights.
#' @seealso [metric.mds()], [landmark.mds.constraints()]
#' @examples
#' X <- cbind(seq_len(12), sin(seq_len(12)))
#' prep <- prepare.sparse.mds(X, method = "uniform", n.pivots = 3)
#' fit <- metric.mds(prepared = prep, approximation = "sparse", max.iter = 5)
#' fit$metadata$data_preparation$method
#' @export
prepare.sparse.mds <- function(X, method = c("euclidean", "uniform", "geodesic"),
                               n.pivots = 200L,
                               pivot.selection = c("randomized", "farthest"),
                               q = 10L, neighbor.variance.target = 0.9,
                               seed = 1L, max.workspace.bytes = 512 * 1024^2) {
  method <- match.arg(method)
  pivot.selection <- match.arg(pivot.selection)
  if (method != "euclidean" && pivot.selection != "randomized")
    stop("farthest selection is supported only by the euclidean method", call. = FALSE)
  if (!(is.matrix(X) || is.data.frame(X))) stop("X must be a numeric matrix", call. = FALSE)
  X <- as.matrix(X)
  if (!is.numeric(X) || nrow(X)<3L || ncol(X)<1L || any(!is.finite(X)))
    stop("X must have at least three finite numeric rows and one column", call. = FALSE)
  storage.mode(X) <- "double"
  if (anyDuplicated(as.data.frame(X))) stop("X contains duplicate rows; consolidate duplicates", call. = FALSE)
  n <- nrow(X)
  n.pivots <- min(grip.validate.vertex.count(n.pivots), n-1L)
  q <- grip.validate.vertex.count(q)
  if (method != "uniform" && q >= n) stop("q must be smaller than nrow(X)", call. = FALSE)
  if (!is.numeric(seed) || length(seed)!=1L || !is.finite(seed) || seed<0 ||
      seed!=floor(seed) || seed>.Machine$integer.max)
    stop("seed must be a nonnegative integer", call. = FALSE)
  grip.validate.scalar(max.workspace.bytes,"max.workspace.bytes",lower=0,open.lower=TRUE)
  if (!is.null(neighbor.variance.target))
    grip.validate.scalar(neighbor.variance.target,"neighbor.variance.target",lower=0,upper=1,open.lower=TRUE)
  # Include input/search copies, PCA work, saved rows, native terms and output.
  estimate <- 64 * as.double(n) * ncol(X) +
    64 * min(n,ncol(X))^2 + 512 * as.double(n) * (n.pivots + if(method=="uniform") 0 else q)
  if (!is.finite(estimate) || estimate > max.workspace.bytes)
    stop("Sparse preparation exceeds max.workspace.bytes before allocation", call. = FALSE)
  began <- proc.time()[["elapsed"]]
  graph <- NULL; pca <- NULL
  edges <- matrix(integer(), ncol=2L); lengths <- numeric()
  if (method != "uniform") {
    if (!requireNamespace("dgraphs",quietly=TRUE)) stop("dgraphs is required for local graph construction",call.=FALSE)
    if (!"search.data" %in% names(formals(dgraphs::create.sknn.graph)))
      stop("Update dgraphs: search.data and streaming component repair are required",call.=FALSE)
    U <- X
    if (!is.null(neighbor.variance.target)) {
      fit <- stats::prcomp(X,center=TRUE,scale.=FALSE)
      variance <- fit$sdev^2
      total <- sum(variance)
      if (!is.finite(total) || total<=0) stop("PCA variance is not positive and representable",call.=FALSE)
      cumulative <- cumsum(variance)/total
      resolved <- which(variance > max(dim(X)) * .Machine$double.eps * variance[1])
      r <- if (neighbor.variance.target == 1) max(resolved) else
        which(cumulative >= neighbor.variance.target)[1L]
      if (is.na(r)) stop("PCA retained-variance target was not reached",call.=FALSE)
      U <- fit$x[,seq_len(r),drop=FALSE]
      pca <- list(center=fit$center,loadings=fit$rotation[,seq_len(r),drop=FALSE],
        dimension=r,achieved=cumulative[r],target=neighbor.variance.target,
        variance=variance)
    }
    graph <- dgraphs::create.sknn.graph(X,k=q,search.data=U,neighbor.method="ann",
      ann.eps=0,connect.components=method=="geodesic",connect.method="component.mst",
      prune.method="none",prune.edges=FALSE,
      graph.detail=if(method=="geodesic") "full" else "minimal")
    if (!identical(graph$metadata$dense_distance_matrix_allocated,FALSE))
      stop("dgraphs did not certify matrix-free graph construction; update dgraphs",call.=FALSE)
    ed <- dgraphs::graph.edges(graph)
    edges <- as.matrix(ed[,c("from","to")]); storage.mode(edges)<-"integer"
    lengths <- ed$length
    if (any(!is.finite(lengths)) || any(lengths<=0)) stop("Invalid original-space edge lengths",call.=FALSE)
  }
  selection <- if (method=="uniform") "uniform" else pivot.selection
  sparse <- grip_sparse_prepare_cpp(n,edges,lengths,n.pivots,integer(),as.integer(seed),
    max.workspace.bytes, observations=if(method=="geodesic") NULL else X,
    selection=selection,region_weighting=method!="uniform",save_distances=TRUE)
  constraints <- sparse[c("pairs","targets","count_i","count_j")]
  metadata <- list(method=method,target_metric=if(method=="geodesic") "graph_geodesic" else "euclidean",
    pivot_selection=selection,n_pivots=n.pivots,q=if(method=="uniform") NULL else q,
    neighbor_variance_target=if(method=="uniform") NULL else neighbor.variance.target,
    seed=as.integer(seed),workspace_estimate_bytes=max(estimate,sparse$workspace_estimate_bytes),
    max_workspace_bytes=max.workspace.bytes,dense_distance_matrix_allocated=FALSE,
    preparation_seconds=proc.time()[["elapsed"]]-began,algorithm_version="observation-sparse-v1")
  structure(list(n=n,input=X,constraints=constraints,pivots=sparse$pivots,
    region=sparse$region,distance_rows=sparse$distance_rows,graph=graph,pca=pca,
    metadata=metadata),class="grip_sparse_mds_preparation")
}
