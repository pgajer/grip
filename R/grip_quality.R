# Layout-quality metrics and comparison helpers for real-world graphs.

grip.validate.coords <- function(coords) {
  coords <- as.matrix(coords)
  if (!is.numeric(coords) || !(ncol(coords) %in% c(2L, 3L))) {
    stop("coords must be a numeric matrix with 2 or 3 columns")
  }
  if (nrow(coords) < 2L) {
    stop("coords must have at least 2 rows")
  }
  if (any(!is.finite(coords))) {
    stop("coords must contain only finite values")
  }
  coords
}

grip.normalize.coords <- function(coords) {
  centered <- scale(coords, center = TRUE, scale = FALSE)
  radius <- max(sqrt(rowSums(centered^2)))
  if (!is.finite(radius) || radius <= 0) {
    return(centered)
  }
  centered / radius
}

grip.align.to.target.nd <- function(source, target, allow.reflection = TRUE) {
  src <- grip.normalize.coords(source)
  dst <- grip.normalize.coords(target)
  cross <- t(src) %*% dst
  sv <- svd(cross)
  rot <- sv$u %*% t(sv$v)
  if (!allow.reflection && ncol(rot) > 1L && det(rot) < 0) {
    fix <- diag(ncol(rot))
    fix[ncol(fix), ncol(fix)] <- -1
    rot <- sv$u %*% fix %*% t(sv$v)
  }
  aligned <- src %*% rot
  list(
    aligned = aligned,
    target = dst,
    rmse = sqrt(mean(rowSums((aligned - dst)^2)))
  )
}

grip.edges.from.adj.list <- function(adj.list) {
  edges <- list()
  for (u in seq_along(adj.list)) {
    nb <- adj.list[[u]]
    if (length(nb) == 0L) next
    keep <- nb > u
    if (!any(keep)) next
    for (v in nb[keep]) {
      edges[[length(edges) + 1L]] <- c(u, v)
    }
  }
  .normalize_undirected_edges(.bind_edges(edges))
}

grip.make.adj.list <- function(edges, n) {
  adj <- vector("list", n)
  for (i in seq_len(nrow(edges))) {
    u <- edges[i, 1L]
    v <- edges[i, 2L]
    adj[[u]] <- c(adj[[u]], v)
    adj[[v]] <- c(adj[[v]], u)
  }
  lapply(adj, as.integer)
}

grip.bfs.distances <- function(adj.list, source) {
  n <- length(adj.list)
  dist <- rep.int(-1L, n)
  q <- integer(n)
  head <- 1L
  tail <- 1L
  q[[tail]] <- source
  dist[[source]] <- 0L
  while (head <= tail) {
    v <- q[[head]]
    head <- head + 1L
    for (u in adj.list[[v]]) {
      if (dist[[u]] == -1L) {
        tail <- tail + 1L
        q[[tail]] <- u
        dist[[u]] <- dist[[v]] + 1L
      }
    }
  }
  as.double(dist)
}

grip.dijkstra.distances <- function(adj.list, weight.list, source) {
  n <- length(adj.list)
  dist <- rep.int(Inf, n)
  visited <- rep.int(FALSE, n)
  dist[[source]] <- 0
  for (iter in seq_len(n)) {
    remaining <- which(!visited)
    if (length(remaining) == 0L) break
    v <- remaining[[which.min(dist[remaining])]]
    if (!is.finite(dist[[v]])) break
    visited[[v]] <- TRUE
    nb <- adj.list[[v]]
    if (length(nb) == 0L) next
    ww <- weight.list[[v]]
    alt <- dist[[v]] + ww
    improve <- alt < dist[nb]
    if (any(improve)) {
      dist[nb[improve]] <- alt[improve]
    }
  }
  dist
}

grip.graph.distances <- function(adj.list, weight.list, source) {
  if (is.null(weight.list)) {
    return(grip.bfs.distances(adj.list, source))
  }
  grip.dijkstra.distances(adj.list, weight.list, source)
}

grip.sample.vertex.pairs <- function(n, sample.size, rng.seed) {
  sample.size <- as.integer(sample.size)
  if (is.na(sample.size) || sample.size <= 0L || n < 2L) {
    return(matrix(integer(), ncol = 2L))
  }
  set.seed(rng.seed)
  out <- matrix(0L, nrow = sample.size, ncol = 2L)
  for (i in seq_len(sample.size)) {
    pair <- sample.int(n, 2L, replace = FALSE)
    out[i, ] <- sort(pair)
  }
  unique(out)
}

grip.edge.length.stats <- function(coords, edges) {
  if (nrow(edges) == 0L) {
    return(list(lengths = numeric(0L), median = NA_real_, cv = NA_real_))
  }
  diffs <- coords[edges[, 1L], , drop = FALSE] - coords[edges[, 2L], , drop = FALSE]
  lengths <- sqrt(rowSums(diffs^2))
  mean.length <- mean(lengths)
  cv <- if (!is.finite(mean.length) || mean.length <= 0) NA_real_ else stats::sd(lengths) / mean.length
  list(
    lengths = lengths,
    median = stats::median(lengths),
    cv = cv
  )
}

grip.sampled.stress <- function(coords,
                                adj.list,
                                weight.list = NULL,
                                sample.size = 2000L,
                                rng.seed = 1L) {
  pairs <- grip.sample.vertex.pairs(nrow(coords), sample.size, rng.seed)
  if (nrow(pairs) == 0L) {
    return(NA_real_)
  }

  src.ids <- unique(pairs[, 1L])
  dist.map <- vector("list", length(src.ids))
  names(dist.map) <- as.character(src.ids)
  for (src in src.ids) {
    dist.map[[as.character(src)]] <- grip.graph.distances(adj.list, weight.list, src)
  }

  gd <- numeric(nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    gd[[i]] <- dist.map[[as.character(pairs[i, 1L])]][[pairs[i, 2L]]]
  }
  keep <- is.finite(gd) & gd > 0
  if (!any(keep)) {
    return(NA_real_)
  }

  pairs <- pairs[keep, , drop = FALSE]
  gd <- as.double(gd[keep])
  ed <- sqrt(rowSums((coords[pairs[, 1L], , drop = FALSE] -
                      coords[pairs[, 2L], , drop = FALSE])^2))
  denom <- sum(gd * gd)
  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }
  scale.factor <- sum(ed * gd) / denom
  sqrt(mean((ed - scale.factor * gd)^2))
}

grip.sampled.nonedge.sep.ratio <- function(coords,
                                           edges,
                                           sample.size = 5000L,
                                           rng.seed = 1L) {
  n <- nrow(coords)
  edge.stats <- grip.edge.length.stats(coords, edges)
  if (!is.finite(edge.stats$median) || edge.stats$median <= 0) {
    return(NA_real_)
  }

  edge.keys <- paste(edges[, 1L], edges[, 2L], sep = "-")
  edge.set <- unique(edge.keys)
  set.seed(rng.seed)
  samples <- numeric(0L)
  attempts <- 0L
  max.attempts <- as.integer(sample.size) * 20L
  while (length(samples) < sample.size && attempts < max.attempts) {
    pair <- sort(sample.int(n, 2L, replace = FALSE))
    key <- paste(pair[[1L]], pair[[2L]], sep = "-")
    attempts <- attempts + 1L
    if (key %in% edge.set) next
    diff <- coords[pair[[1L]], ] - coords[pair[[2L]], ]
    samples <- c(samples, sqrt(sum(diff * diff)))
  }

  if (length(samples) == 0L) {
    return(NA_real_)
  }
  min(samples) / edge.stats$median
}

grip.segment.orientation <- function(a, b, c) {
  (b[[1L]] - a[[1L]]) * (c[[2L]] - a[[2L]]) -
    (b[[2L]] - a[[2L]]) * (c[[1L]] - a[[1L]])
}

grip.segments.intersect.strict <- function(a, b, c, d, tol = 1e-12) {
  o1 <- grip.segment.orientation(a, b, c)
  o2 <- grip.segment.orientation(a, b, d)
  o3 <- grip.segment.orientation(c, d, a)
  o4 <- grip.segment.orientation(c, d, b)
  if (abs(o1) <= tol || abs(o2) <= tol || abs(o3) <= tol || abs(o4) <= tol) {
    return(FALSE)
  }
  ((o1 > 0) != (o2 > 0)) && ((o3 > 0) != (o4 > 0))
}

grip.count.edge.crossings <- function(coords, edges) {
  if (ncol(coords) != 2L) {
    return(NA_integer_)
  }
  edges <- .normalize_undirected_edges(edges)
  m <- nrow(edges)
  if (m < 2L) {
    return(0L)
  }

  crossings <- 0L
  for (i in seq_len(m - 1L)) {
    e1 <- edges[i, ]
    for (j in (i + 1L):m) {
      e2 <- edges[j, ]
      if (length(intersect(e1, e2)) > 0L) next
      if (grip.segments.intersect.strict(coords[e1[1L], ], coords[e1[2L], ],
                                         coords[e2[1L], ], coords[e2[2L], ])) {
        crossings <- crossings + 1L
      }
    }
  }
  crossings
}

grip.cluster.separation <- function(coords, clusters) {
  if (is.null(clusters)) {
    return(NA_real_)
  }
  if (length(clusters) != nrow(coords)) {
    stop("clusters must have length equal to nrow(coords)")
  }
  keep <- !is.na(clusters)
  coords <- coords[keep, , drop = FALSE]
  clusters <- as.factor(clusters[keep])
  if (nlevels(clusters) < 2L || nrow(coords) < 2L) {
    return(NA_real_)
  }

  split.idx <- split(seq_len(nrow(coords)), clusters)
  centers <- do.call(rbind, lapply(split.idx, function(idx) colMeans(coords[idx, , drop = FALSE])))
  within <- mean(vapply(split.idx, function(idx) {
    grp <- coords[idx, , drop = FALSE]
    ctr <- colMeans(grp)
    mean(sqrt(rowSums((grp - matrix(ctr, nrow(grp), ncol(grp), byrow = TRUE))^2)))
  }, numeric(1L)))

  if (nrow(centers) < 2L) {
    return(NA_real_)
  }
  pair.idx <- utils::combn(seq_len(nrow(centers)), 2L)
  between <- mean(apply(pair.idx, 2L, function(ii) {
    sqrt(sum((centers[ii[1L], ] - centers[ii[2L], ])^2))
  }))

  if (!is.finite(within) || within < 0) {
    return(NA_real_)
  }
  if (within == 0) {
    return(if (between > 0) Inf else NA_real_)
  }
  between / within
}

grip.rank01 <- function(x, higher.better = FALSE) {
  n <- length(x)
  out <- rep(1, n)
  ok <- is.finite(x)
  if (!any(ok)) {
    return(out)
  }
  vals <- if (higher.better) -x[ok] else x[ok]
  if (length(vals) == 1L || all(abs(vals - vals[[1L]]) <= sqrt(.Machine$double.eps))) {
    out[ok] <- 0
    return(out)
  }
  ranks <- rank(vals, ties.method = "average")
  out[ok] <- (ranks - 1) / (length(ranks) - 1)
  out
}

grip.default.compare.score.weights <- function() {
  c(
    sampled.stress = 0.35,
    edge.length.cv = 0.20,
    sampled.nonedge.sep.ratio = 0.20,
    stability.procrustes.mean = 0.15,
    edge.crossings = 0.05,
    cluster.separation = 0.05
  )
}

grip.stability.procrustes <- function(layouts) {
  if (length(layouts) < 2L) {
    return(c(stability.procrustes.mean = NA_real_, stability.procrustes.max = NA_real_))
  }
  pair.idx <- utils::combn(seq_along(layouts), 2L)
  values <- apply(pair.idx, 2L, function(ii) {
    grip.align.to.target.nd(layouts[[ii[1L]]], layouts[[ii[2L]]], allow.reflection = TRUE)$rmse
  })
  c(
    stability.procrustes.mean = mean(values),
    stability.procrustes.max = max(values)
  )
}

grip.normalize.compare.candidates <- function(candidates) {
  if (is.null(candidates)) {
    return(list(default = list()))
  }

  if (is.character(candidates)) {
    vals <- as.character(candidates)
    nm <- names(candidates)
    out <- vector("list", length(vals))
    for (i in seq_along(vals)) {
      value <- vals[[i]]
      out[[i]] <- if (identical(value, "default")) list() else list(preset = value)
      if (is.null(nm) || !nzchar(nm[[i]])) {
        nm_i <- value
      } else {
        nm_i <- nm[[i]]
      }
      names(out)[[i]] <- nm_i
    }
    if (anyDuplicated(names(out))) {
      stop("candidate names must be unique")
    }
    return(out)
  }

  if (!is.list(candidates) || length(candidates) == 0L) {
    stop("candidates must be NULL, a character vector, or a named list")
  }

  allowed <- c("placement", "preset", "rounds", "final_rounds", "num_init",
               "num_nbrs", "r", "s", "repulsion_factor", "tinit_factor")
  nm <- names(candidates)
  out <- vector("list", length(candidates))
  for (i in seq_along(candidates)) {
    cand <- candidates[[i]]
    if (is.null(cand)) {
      args <- list()
    } else if (is.character(cand) && length(cand) == 1L && !is.na(cand)) {
      args <- if (identical(cand, "default")) list() else list(preset = cand)
    } else if (is.list(cand)) {
      args <- cand
    } else {
      stop("each candidate must be NULL, a length-1 character value, or a named list")
    }

    bad <- setdiff(names(args), allowed)
    if (length(bad) > 0L) {
      stop(sprintf(
        "candidate '%s' contains unsupported layout arguments: %s",
        if (!is.null(nm) && nzchar(nm[[i]])) nm[[i]] else paste0("candidate.", i),
        paste(bad, collapse = ", ")
      ))
    }

    out[[i]] <- args
    names(out)[[i]] <- if (!is.null(nm) && nzchar(nm[[i]])) nm[[i]] else paste0("candidate.", i)
  }

  if (anyDuplicated(names(out))) {
    stop("candidate names must be unique")
  }
  out
}

grip.resolve.compare.candidate <- function(candidate) {
  defaults <- list(
    placement = "barycenter",
    rounds = 20,
    final_rounds = 25,
    num_init = 36,
    num_nbrs = 10,
    r = 0.15,
    s = 3.0,
    repulsion_factor = 1.0,
    tinit_factor = 6
  )

  if (is.null(candidate)) {
    candidate <- list()
  }
  placement.missing <- !("placement" %in% names(candidate))
  rounds.missing <- !("rounds" %in% names(candidate))
  final.rounds.missing <- !("final_rounds" %in% names(candidate))
  num.init.missing <- !("num_init" %in% names(candidate))
  num.nbrs.missing <- !("num_nbrs" %in% names(candidate))
  r.missing <- !("r" %in% names(candidate))
  s.missing <- !("s" %in% names(candidate))
  repulsion.missing <- !("repulsion_factor" %in% names(candidate))

  preset <- if ("preset" %in% names(candidate)) candidate$preset else NULL
  preset <- grip.normalize.preset(preset, fn = "grip.compare.layouts")

  resolved <- grip.resolve.preset(
    preset = preset,
    placement = if (placement.missing) defaults$placement else candidate$placement,
    placement_missing = placement.missing,
    rounds = if (rounds.missing) defaults$rounds else candidate$rounds,
    rounds_missing = rounds.missing,
    final_rounds = if (final.rounds.missing) defaults$final_rounds else candidate$final_rounds,
    final_rounds_missing = final.rounds.missing,
    num_init = if (num.init.missing) defaults$num_init else candidate$num_init,
    num_init_missing = num.init.missing,
    num_nbrs = if (num.nbrs.missing) defaults$num_nbrs else candidate$num_nbrs,
    num_nbrs_missing = num.nbrs.missing,
    r = if (r.missing) defaults$r else candidate$r,
    r_missing = r.missing,
    s = if (s.missing) defaults$s else candidate$s,
    s_missing = s.missing,
    repulsion_factor = if (repulsion.missing) defaults$repulsion_factor else candidate$repulsion_factor,
    repulsion_factor_missing = repulsion.missing
  )

  tuning <- grip.validate.tuning.inputs(
    num_nbrs = resolved$num_nbrs,
    r = resolved$r,
    s = resolved$s,
    repulsion_factor = resolved$repulsion_factor
  )

  tinit.factor <- if ("tinit_factor" %in% names(candidate)) candidate$tinit_factor else defaults$tinit_factor
  if (!is.numeric(tinit.factor) || length(tinit.factor) != 1L || !is.finite(tinit.factor) || tinit.factor <= 0) {
    stop("tinit_factor must be a single finite numeric value > 0")
  }

  list(
    preset = preset,
    placement = match.arg(resolved$placement, c("barycenter", "circle")),
    rounds = as.integer(resolved$rounds),
    final_rounds = as.integer(resolved$final_rounds),
    num_init = as.integer(resolved$num_init),
    num_nbrs = tuning$num_nbrs,
    r = tuning$r,
    s = tuning$s,
    repulsion_factor = tuning$repulsion_factor,
    tinit_factor = as.double(tinit.factor)
  )
}

grip.compare.summary <- function(runs, layouts.by.candidate, score.weights) {
  groups <- split(runs, runs$candidate)
  summary <- do.call(rbind, lapply(groups, function(df) {
    ok <- df$status == "ok"
    good <- df[ok, , drop = FALSE]
    stability <- grip.stability.procrustes(layouts.by.candidate[[df$candidate[[1L]]]])
    out <- data.frame(
      candidate = df$candidate[[1L]],
      preset = if (all(is.na(df$preset) | df$preset == "")) NA_character_ else df$preset[[which.max(nzchar(df$preset))]],
      n.runs = nrow(df),
      n.ok = sum(ok),
      n.fail = sum(!ok),
      sampled.stress.mean = if (nrow(good) > 0L) mean(good$sampled.stress) else NA_real_,
      sampled.stress.sd = if (nrow(good) > 1L) stats::sd(good$sampled.stress) else 0,
      edge.length.cv.mean = if (nrow(good) > 0L) mean(good$edge.length.cv) else NA_real_,
      edge.length.cv.sd = if (nrow(good) > 1L) stats::sd(good$edge.length.cv) else 0,
      sampled.nonedge.sep.ratio.mean = if (nrow(good) > 0L) mean(good$sampled.nonedge.sep.ratio) else NA_real_,
      sampled.nonedge.sep.ratio.sd = if (nrow(good) > 1L) stats::sd(good$sampled.nonedge.sep.ratio) else 0,
      edge.crossings.mean = if (nrow(good) > 0L) mean(good$edge.crossings, na.rm = TRUE) else NA_real_,
      edge.crossings.sd = if (sum(is.finite(good$edge.crossings)) > 1L) stats::sd(good$edge.crossings, na.rm = TRUE) else 0,
      cluster.separation.mean = if (nrow(good) > 0L) mean(good$cluster.separation, na.rm = TRUE) else NA_real_,
      cluster.separation.sd = if (sum(is.finite(good$cluster.separation)) > 1L) stats::sd(good$cluster.separation, na.rm = TRUE) else 0,
      stability.procrustes.mean = stability[["stability.procrustes.mean"]],
      stability.procrustes.max = stability[["stability.procrustes.max"]],
      elapsed.sec.mean = if (nrow(good) > 0L) mean(good$elapsed.sec) else NA_real_,
      stringsAsFactors = FALSE
    )
    out
  }))

  if (!is.null(score.weights)) {
    metric.map <- list(
      sampled.stress = list(column = "sampled.stress.mean", higher.better = FALSE),
      edge.length.cv = list(column = "edge.length.cv.mean", higher.better = FALSE),
      sampled.nonedge.sep.ratio = list(column = "sampled.nonedge.sep.ratio.mean", higher.better = TRUE),
      stability.procrustes.mean = list(column = "stability.procrustes.mean", higher.better = FALSE),
      edge.crossings = list(column = "edge.crossings.mean", higher.better = FALSE),
      cluster.separation = list(column = "cluster.separation.mean", higher.better = TRUE)
    )
    score.weights <- score.weights[names(score.weights) %in% names(metric.map)]
    present <- vapply(names(score.weights), function(metric) {
      col <- metric.map[[metric]]$column
      any(is.finite(summary[[col]]))
    }, logical(1L))
    if (any(present)) {
      score.weights <- score.weights[present]
      score.weights <- score.weights / sum(score.weights)
      score <- rep(0, nrow(summary))
      for (metric in names(score.weights)) {
        col <- metric.map[[metric]]$column
        rank.col <- grip.rank01(summary[[col]], higher.better = metric.map[[metric]]$higher.better)
        score <- score + score.weights[[metric]] * rank.col
      }
      summary$score.composite <- score
      summary <- summary[order(summary$score.composite, summary$sampled.stress.mean), , drop = FALSE]
    } else {
      summary$score.composite <- NA_real_
    }
  }

  rownames(summary) <- NULL
  summary
}

#' Score a single layout using graph-aware quality heuristics
#'
#' \code{grip.score.layout()} evaluates a realized layout without assuming a
#' canonical embedding. It is intended for real-world graphs where visual
#' quality is judged by graph-distance faithfulness, edge-length consistency,
#' separation of non-neighbors, and optionally edge crossings or cluster
#' separation.
#'
#' @param coords Numeric coordinate matrix with 2 or 3 columns.
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#' @param n Number of vertices. If omitted with \code{adj_list}, defaults to
#'   \code{length(adj_list)}. If omitted with \code{edges}, defaults to
#'   \code{nrow(coords)}.
#' @param adj_list Adjacency list (1-based) for undirected graphs.
#' @param weight_list Optional parallel list of positive edge weights.
#' @param edge_weights Optional positive edge-weight vector parallel to
#'   \code{edges}.
#' @param clusters Optional cluster or community labels of length
#'   \code{nrow(coords)}. When supplied, \code{cluster.separation} is reported.
#' @param sample.size.stress Number of vertex pairs sampled for
#'   \code{sampled.stress}.
#' @param sample.size.nonedge Number of non-edge pairs sampled for
#'   \code{sampled.nonedge.sep.ratio}.
#' @param stress.seed RNG seed used for the stress sample.
#' @param nonedge.seed RNG seed used for the non-edge sample.
#' @param edge.crossings How to compute \code{edge.crossings} for 2D layouts:
#'   \code{"auto"} computes exact crossings only when the graph is small enough,
#'   \code{"always"} always computes them, and \code{"never"} skips them.
#' @param edge.crossings.max.edges Edge-count threshold used by
#'   \code{edge.crossings = "auto"}.
#'
#' @return A one-row data frame with dot-delimited metric names.
#' @examples
#' edges <- edges.mesh(5, 5)
#' coords <- grip.layout(edges, n = 25, dim = 2, preset = "mesh", seed = 1)
#' grip.score.layout(coords, edges = edges, n = 25)
#' @export
grip.score.layout <- function(coords,
                              edges = NULL,
                              n = NULL,
                              adj_list = NULL,
                              weight_list = NULL,
                              edge_weights = NULL,
                              clusters = NULL,
                              sample.size.stress = 2000L,
                              sample.size.nonedge = 5000L,
                              stress.seed = 1L,
                              nonedge.seed = 1L,
                              edge.crossings = c("auto", "always", "never"),
                              edge.crossings.max.edges = 1000L) {
  coords <- grip.validate.coords(coords)
  edge.crossings <- match.arg(edge.crossings)
  if (is.null(n)) {
    n <- nrow(coords)
  }
  validated <- grip.validate.layout.inputs(
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights,
    dim = ncol(coords),
    placement = "barycenter",
    seed = 1
  )
  if (validated$n != nrow(coords)) {
    stop("nrow(coords) must match the graph size")
  }

  edges <- grip.edges.from.adj.list(validated$adj_list)
  edge.stats <- grip.edge.length.stats(coords, edges)
  crossings <- NA_integer_
  if (ncol(coords) == 2L) {
    if (edge.crossings == "always" ||
        (edge.crossings == "auto" && nrow(edges) <= as.integer(edge.crossings.max.edges))) {
      crossings <- grip.count.edge.crossings(coords, edges)
    }
  }

  data.frame(
    n.vertices = nrow(coords),
    n.edges = nrow(edges),
    dim = ncol(coords),
    sampled.stress = grip.sampled.stress(
      coords = coords,
      adj.list = validated$adj_list,
      weight.list = validated$weight_list,
      sample.size = sample.size.stress,
      rng.seed = stress.seed
    ),
    edge.length.cv = edge.stats$cv,
    median.edge.length = edge.stats$median,
    sampled.nonedge.sep.ratio = grip.sampled.nonedge.sep.ratio(
      coords = coords,
      edges = edges,
      sample.size = sample.size.nonedge,
      rng.seed = nonedge.seed
    ),
    edge.crossings = crossings,
    cluster.separation = grip.cluster.separation(coords, clusters),
    stringsAsFactors = FALSE
  )
}

#' Compare multiple layout candidates across seeds
#'
#' \code{grip.compare.layouts()} computes layouts for several candidate presets
#' or parameter lists, scores each run with \code{\link{grip.score.layout}()},
#' and summarizes both quality metrics and seed-to-seed stability. This is
#' useful for real-world datasets where no canonical embedding is known.
#'
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#' @param n Number of vertices.
#' @param adj_list Adjacency list (1-based) for undirected graphs.
#' @param weight_list Optional parallel list of positive edge weights.
#' @param edge_weights Optional positive edge-weight vector parallel to
#'   \code{edges}.
#' @param dim Layout dimension (2 or 3).
#' @param candidates Either a character vector such as
#'   \code{c("default", "mesh", "tree")} or a named list of candidate layout
#'   specifications. Each list element may be \code{NULL} (use defaults), a
#'   single preset name, or a named list of \code{\link{grip.layout}()}
#'   tuning arguments such as \code{preset}, \code{placement},
#'   \code{rounds}, or \code{repulsion_factor}.
#' @param clusters Optional cluster or community labels used to compute
#'   \code{cluster.separation}.
#' @param seeds Integer seeds used for repeated runs.
#' @param sample.size.stress Number of sampled pairs used for
#'   \code{sampled.stress}.
#' @param sample.size.nonedge Number of sampled non-edge pairs used for
#'   \code{sampled.nonedge.sep.ratio}.
#' @param edge.crossings How to compute \code{edge.crossings} for 2D layouts.
#' @param edge.crossings.max.edges Edge-count threshold for
#'   \code{edge.crossings = "auto"}.
#' @param score.weights Optional named numeric vector used to compute
#'   \code{score.composite}. Set to \code{NULL} to omit the composite score.
#' @param return.layouts If \code{TRUE}, include the realized coordinate
#'   matrices in the return value.
#' @param disconnected Passed through to \code{\link{grip.layout}()}.
#'
#' @return A list with \code{runs} and \code{summary} data frames and,
#'   optionally, realized \code{layouts}.
#' @examples
#' edges <- edges.mesh(6, 6)
#' cmp <- grip.compare.layouts(
#'   edges = edges,
#'   n = 36,
#'   dim = 2,
#'   candidates = c("default", "mesh"),
#'   seeds = 1:2
#' )
#' cmp$summary
#' @export
grip.compare.layouts <- function(edges = NULL,
                                 n = NULL,
                                 adj_list = NULL,
                                 weight_list = NULL,
                                 edge_weights = NULL,
                                 dim = 2,
                                 candidates = c("default"),
                                 clusters = NULL,
                                 seeds = 1:3,
                                 sample.size.stress = 2000L,
                                 sample.size.nonedge = 5000L,
                                 edge.crossings = c("auto", "always", "never"),
                                 edge.crossings.max.edges = 1000L,
                                 score.weights = grip.default.compare.score.weights(),
                                 return.layouts = FALSE,
                                 disconnected = c("components", "error")) {
  edge.crossings <- match.arg(edge.crossings)
  disconnected <- match.arg(disconnected)
  if (is.null(n) && !is.null(edges)) {
    n <- max(as.matrix(edges))
  }
  validated <- grip.validate.layout.inputs(
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights,
    dim = dim,
    placement = "barycenter",
    seed = 1
  )
  seeds <- as.integer(seeds)
  if (length(seeds) == 0L || any(is.na(seeds))) {
    stop("seeds must be a non-empty integer vector")
  }
  candidate.list <- grip.normalize.compare.candidates(candidates)
  graph.args <- list(
    adj_list = validated$adj_list,
    weight_list = validated$weight_list,
    n = validated$n
  )

  layouts.by.candidate <- vector("list", length(candidate.list))
  names(layouts.by.candidate) <- names(candidate.list)
  run.rows <- list()
  idx <- 0L

  for (candidate.name in names(candidate.list)) {
    resolved <- grip.resolve.compare.candidate(candidate.list[[candidate.name]])
    layouts.by.candidate[[candidate.name]] <- list()
    for (seed in seeds) {
      idx <- idx + 1L
      started <- proc.time()[["elapsed"]]
      row <- tryCatch({
        coords <- do.call(
          grip.layout,
          c(
            graph.args,
            list(
              dim = dim,
              placement = resolved$placement,
              preset = resolved$preset,
              rounds = resolved$rounds,
              final_rounds = resolved$final_rounds,
              num_init = resolved$num_init,
              num_nbrs = resolved$num_nbrs,
              r = resolved$r,
              s = resolved$s,
              repulsion_factor = resolved$repulsion_factor,
              tinit_factor = resolved$tinit_factor,
              seed = seed,
              disconnected = disconnected
            )
          )
        )
        layouts.by.candidate[[candidate.name]][[as.character(seed)]] <- coords
        score <- grip.score.layout(
          coords = coords,
          adj_list = validated$adj_list,
          weight_list = validated$weight_list,
          n = validated$n,
          clusters = clusters,
          sample.size.stress = sample.size.stress,
          sample.size.nonedge = sample.size.nonedge,
          stress.seed = 1000L + seed,
          nonedge.seed = 2000L + seed,
          edge.crossings = edge.crossings,
          edge.crossings.max.edges = edge.crossings.max.edges
        )
        cbind(
          data.frame(
            candidate = candidate.name,
            seed = seed,
            status = "ok",
            error.message = "",
            preset = if (is.null(resolved$preset)) "" else resolved$preset,
            placement = resolved$placement,
            rounds = resolved$rounds,
            final.rounds = resolved$final_rounds,
            num.init = resolved$num_init,
            num.nbrs = resolved$num_nbrs,
            r = resolved$r,
            s = resolved$s,
            repulsion.factor = resolved$repulsion_factor,
            tinit.factor = resolved$tinit_factor,
            elapsed.sec = proc.time()[["elapsed"]] - started,
            stringsAsFactors = FALSE
          ),
          score
        )
      }, error = function(e) {
        data.frame(
          candidate = candidate.name,
          seed = seed,
          status = "error",
          error.message = conditionMessage(e),
          preset = if (is.null(resolved$preset)) "" else resolved$preset,
          placement = resolved$placement,
          rounds = resolved$rounds,
          final.rounds = resolved$final_rounds,
          num.init = resolved$num_init,
          num.nbrs = resolved$num_nbrs,
          r = resolved$r,
          s = resolved$s,
          repulsion.factor = resolved$repulsion_factor,
          tinit.factor = resolved$tinit_factor,
          elapsed.sec = proc.time()[["elapsed"]] - started,
          n.vertices = validated$n,
          n.edges = nrow(grip.edges.from.adj.list(validated$adj_list)),
          dim = dim,
          sampled.stress = NA_real_,
          edge.length.cv = NA_real_,
          median.edge.length = NA_real_,
          sampled.nonedge.sep.ratio = NA_real_,
          edge.crossings = NA_real_,
          cluster.separation = NA_real_,
          stringsAsFactors = FALSE
        )
      })
      run.rows[[idx]] <- row
    }
  }

  runs <- do.call(rbind, run.rows)
  rownames(runs) <- NULL
  summary <- grip.compare.summary(runs, layouts.by.candidate, score.weights = score.weights)

  out <- list(
    runs = runs,
    summary = summary
  )
  if (isTRUE(return.layouts)) {
    out$layouts <- layouts.by.candidate
  }
  out
}
