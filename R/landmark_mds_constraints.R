#' Construct sparse MDS targets from landmark distances
#'
#' Converts supplied distances into soft pairwise fitting targets for
#' [metric.mds()]. No distances are recomputed and no graph paths are inferred.
#' With `weighting = "uniform"`, every retained unordered pair has endpoint
#' multiplicities one, giving inverse-squared distance stress on those pairs.
#' With `weighting = "region"`, samples belong to their nearest landmark (ties
#' use landmark row order). The update at a non-landmark endpoint represents
#' the number of samples in the landmark's region at most half the target
#' distance from that landmark. Opposite landmark directions are combined.
#' Optional local pairs replace landmark contributions for that pair with
#' multiplicities one at each endpoint. Actual fitting weights divide these
#' multiplicities by squared target distances. Region weighting can give
#' asymmetric endpoint updates; it is not ordinary symmetric weighted stress.
#'
#' @param distance.matrix Numeric landmark-by-sample distance matrix, finite
#'   and nonnegative. Each row has zero at its landmark and strictly positive
#'   values elsewhere. Consolidate zero-distance duplicates before use.
#' @param landmarks Distinct integer sample indices identifying matrix rows,
#'   in row order; sample indices correspond to matrix columns.
#' @param local.pairs Optional two-column integer matrix of sample indices.
#' @param local.distances Positive finite distances for `local.pairs`, from
#'   the same metric as `distance.matrix`. Repeated pairs and overlapping
#'   landmark targets must agree to relative tolerance 1e-10.
#' @param weighting Endpoint weighting scheme: `"uniform"` or `"region"`.
#' @return List with `pairs`, `targets`, `count_i`, and `count_j`, suitable
#'   for `metric.mds(n = ncol(distance.matrix), constraints = ...,
#'   approximation = "sparse")`. Distances are retained in their supplied units.
#' @export
#' @examples
#' X <- matrix(sin(seq_len(24)), 8, 3)
#' D <- as.matrix(dist(X))
#' C <- landmark.mds.constraints(D[c(1, 5), ], c(1, 5))
#' fit <- metric.mds(n = 8, constraints = C, approximation = "sparse",
#'                   dim = 3, max.iter = 5)
landmark.mds.constraints <- function(distance.matrix, landmarks,
                                     local.pairs = NULL, local.distances = NULL,
                                     weighting = c("uniform", "region")) {
  weighting <- match.arg(weighting)
  D <- distance.matrix
  if (!is.matrix(D) || !is.numeric(D) || !nrow(D) || ncol(D) < 2L ||
      any(!is.finite(D)) || any(D < 0))
    stop("distance.matrix must be a finite nonnegative landmark-by-sample matrix", call. = FALSE)
  n <- ncol(D)
  if (!is.null(dim(landmarks)))
    stop("landmarks must be an index vector", call. = FALSE)
  grip.validate.vertex.ids(landmarks, "landmarks", n)
  if (length(landmarks) != nrow(D) || anyDuplicated(landmarks))
    stop("landmarks must be distinct with one index per matrix row", call. = FALSE)
  self <- cbind(seq_len(nrow(D)), landmarks)
  if (any(D[self] != 0) || any(rowSums(D == 0) != 1L))
    stop("Each landmark row must have only its self-distance zero; consolidate duplicates", call. = FALSE)
  if (xor(is.null(local.pairs), is.null(local.distances)))
    stop("Supply both local.pairs and local.distances", call. = FALSE)
  close <- function(a, b) abs(a-b) <= 1e-10 * pmax(abs(a), abs(b))
  if (any(!close(D[, landmarks, drop = FALSE], t(D[, landmarks, drop = FALSE]))))
    stop("Opposite landmark distances must agree", call. = FALSE)
  # Avoid materializing and sorting a six-column double table for tens of
  # millions of landmark targets. Emit each unordered pair once, in the same
  # lexicographic order as the general adapter below.
  if(weighting == "uniform" && is.null(local.pairs)) {
    h <- length(landmarks); size <- as.double(h)*n-h*(h+1)/2
    if(size>.Machine$integer.max)stop("Too many landmark pairs for an R matrix",call.=FALSE)
    pairs<-matrix(0L,size,2L);targets<-numeric(size)
    landmark.row<-match(seq_len(n),landmarks);sorted<-sort(landmarks);cursor<-0L
    for(i in seq_len(n-1L)) {
      js<-if(!is.na(landmark.row[i]))seq.int(i+1L,n) else sorted[sorted>i]
      if(!length(js))next
      rows<-seq.int(cursor+1L,cursor+length(js));pairs[rows,1]<-i;pairs[rows,2]<-js
      targets[rows]<-if(!is.na(landmark.row[i]))D[landmark.row[i],js] else D[landmark.row[js],i]
      cursor<-cursor+length(js)
    }
    return(list(pairs=pairs,targets=targets,count_i=rep(1,size),count_j=rep(1,size)))
  }
  owner <- if (weighting == "region") max.col(-t(D), ties.method = "first") else NULL
  rows <- lapply(seq_along(landmarks), function(j) {
    pivot <- landmarks[j]
    ids <- seq_len(n)[seq_len(n) != pivot]
    counts <- if (weighting == "region")
      findInterval(D[j, ids] / 2, sort(D[j, owner == j])) else rep(1, length(ids))
    a <- if (weighting == "region") ifelse(ids < pivot, counts, 0) else counts
    b <- if (weighting == "region") ifelse(ids > pivot, counts, 0) else counts
    cbind(pmin(ids, pivot), pmax(ids, pivot), D[j, ids], a, b, 0)
  })
  if (!is.null(local.pairs)) {
    if (!is.matrix(local.pairs) || !is.numeric(local.pairs) || ncol(local.pairs) != 2L)
      stop("local.pairs must have two numeric columns", call. = FALSE)
    grip.validate.vertex.ids(local.pairs, "local.pairs", n)
    if (any(local.pairs[, 1] == local.pairs[, 2]) ||
        !is.numeric(local.distances) || !is.null(dim(local.distances)) ||
        length(local.distances) != nrow(local.pairs) ||
        any(!is.finite(local.distances)) || any(local.distances <= 0))
      stop("Local targets require distinct endpoints and positive finite distances", call. = FALSE)
    rows[[length(rows) + 1L]] <- cbind(pmin(local.pairs[, 1], local.pairs[, 2]),
      pmax(local.pairs[, 1], local.pairs[, 2]), local.distances, 1, 1, 1)
  }
  rows <- do.call(rbind, rows)
  rows <- rows[order(rows[, 1], rows[, 2]), , drop = FALSE]
  first <- c(TRUE, rows[-1, 1] != rows[-nrow(rows), 1] |
                    rows[-1, 2] != rows[-nrow(rows), 2])
  group <- cumsum(first)
  if (any(!close(rows[, 3], rows[first, 3][group])))
    stop("Repeated pair distances must agree", call. = FALSE)
  counts <- rowsum(rows[, 4:6, drop = FALSE], group, reorder = FALSE)
  if (weighting == "uniform") counts[, 1:2] <- 1
  counts[counts[, 3] > 0, 1:2] <- 1
  pairs <- unname(rows[first, 1:2, drop = FALSE])
  storage.mode(pairs) <- "integer"
  list(pairs = pairs, targets = unname(rows[first, 3]),
       count_i = unname(counts[, 1]), count_j = unname(counts[, 2]))
}
