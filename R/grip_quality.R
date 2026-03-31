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

grip.normalize.coords.with.meta <- function(coords) {
  coords <- grip.validate.coords(coords)
  center <- colMeans(coords)
  centered <- sweep(coords, 2L, center, check.margin = FALSE)
  radius <- max(sqrt(rowSums(centered^2)))
  if (!is.finite(radius) || radius <= 0) {
    radius <- 1
  }
  list(
    center = center,
    radius = radius,
    normalized = centered / radius
  )
}

grip.align.to.target.nd <- function(source, target, allow.reflection = TRUE) {
  src.meta <- grip.normalize.coords.with.meta(source)
  dst.meta <- grip.normalize.coords.with.meta(target)
  src <- src.meta$normalized
  dst <- dst.meta$normalized
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

grip.axis.permutations <- function(dim) {
  if (dim == 2L) {
    return(list(c(1L, 2L), c(2L, 1L)))
  }
  if (dim == 3L) {
    return(list(
      c(1L, 2L, 3L),
      c(1L, 3L, 2L),
      c(2L, 1L, 3L),
      c(2L, 3L, 1L),
      c(3L, 1L, 2L),
      c(3L, 2L, 1L)
    ))
  }
  stop("dim must be 2 or 3")
}

grip.sample.symmetry.points <- function(coords, sample.size, rng.seed) {
  n <- nrow(coords)
  sample.size <- max(1L, min(as.integer(sample.size), n))
  if (sample.size >= n) {
    return(coords)
  }
  set.seed(rng.seed)
  coords[sample.int(n, sample.size, replace = FALSE), , drop = FALSE]
}

grip.min.match.distance <- function(source, target) {
  src2 <- rowSums(source^2)
  dst2 <- rowSums(target^2)
  d2 <- outer(src2, dst2, "+") - 2 * (source %*% t(target))
  mean(sqrt(pmax(0, apply(d2, 1L, min))))
}

grip.global.symmetry.score <- function(coords,
                                       sample.size = 512L,
                                       rng.seed = 1L) {
  coords <- grip.validate.coords(coords)
  centered <- grip.normalize.coords(coords)
  sample.coords <- grip.sample.symmetry.points(centered, sample.size, rng.seed)
  dim <- ncol(sample.coords)
  perms <- grip.axis.permutations(dim)
  sign.grid <- expand.grid(rep(list(c(-1, 1)), dim), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  best.error <- Inf

  for (perm in perms) {
    for (i in seq_len(nrow(sign.grid))) {
      signs <- as.numeric(sign.grid[i, , drop = TRUE])
      if (identical(perm, seq_len(dim)) && all(signs == 1)) {
        next
      }
      transformed <- sweep(sample.coords[, perm, drop = FALSE], 2L, signs, `*`)
      err <- grip.min.match.distance(transformed, sample.coords)
      if (is.finite(err) && err < best.error) {
        best.error <- err
      }
    }
  }

  if (!is.finite(best.error)) {
    return(NA_real_)
  }
  max(0, 1 - best.error / 2)
}

grip.sample.wedges <- function(adj.list, sample.size, rng.seed) {
  eligible <- which(lengths(adj.list) >= 2L)
  if (length(eligible) == 0L) {
    return(matrix(integer(), ncol = 3L))
  }
  sample.size <- as.integer(sample.size)
  if (is.na(sample.size) || sample.size <= 0L) {
    return(matrix(integer(), ncol = 3L))
  }
  set.seed(rng.seed)
  out <- matrix(0L, nrow = sample.size, ncol = 3L)
  for (i in seq_len(sample.size)) {
    v <- sample(eligible, 1L)
    pair <- sample(adj.list[[v]], 2L, replace = FALSE)
    out[i, ] <- c(pair[[1L]], v, pair[[2L]])
  }
  unique(out)
}

grip.angle.between <- function(a, b) {
  na <- sqrt(sum(a * a))
  nb <- sqrt(sum(b * b))
  if (!is.finite(na) || !is.finite(nb) || na <= 0 || nb <= 0) {
    return(NA_real_)
  }
  cos.theta <- sum(a * b) / (na * nb)
  acos(max(-1, min(1, cos.theta)))
}

grip.local.angle.deviation <- function(coords,
                                       target.coords,
                                       edges,
                                       sample.size = 4000L,
                                       rng.seed = 1L) {
  coords <- grip.validate.coords(coords)
  target.coords <- grip.validate.coords(target.coords)
  if (!identical(dim(coords), dim(target.coords))) {
    stop("coords and target.coords must have the same dimensions")
  }
  fit <- grip.align.to.target.nd(coords, target.coords, allow.reflection = TRUE)
  adj.list <- grip.make.adj.list(edges, nrow(coords))
  wedges <- grip.sample.wedges(adj.list, sample.size, rng.seed)
  if (nrow(wedges) == 0L) {
    return(NA_real_)
  }

  diffs <- numeric(nrow(wedges))
  keep <- logical(nrow(wedges))
  for (i in seq_len(nrow(wedges))) {
    u <- wedges[i, 1L]
    v <- wedges[i, 2L]
    w <- wedges[i, 3L]
    layout.angle <- grip.angle.between(fit$aligned[u, ] - fit$aligned[v, ],
                                       fit$aligned[w, ] - fit$aligned[v, ])
    target.angle <- grip.angle.between(fit$target[u, ] - fit$target[v, ],
                                       fit$target[w, ] - fit$target[v, ])
    if (is.finite(layout.angle) && is.finite(target.angle)) {
      diffs[[i]] <- abs(layout.angle - target.angle) / pi
      keep[[i]] <- TRUE
    }
  }
  if (!any(keep)) {
    return(NA_real_)
  }
  mean(diffs[keep])
}

grip.edge.axis.concentration <- function(coords, edges) {
  coords <- grip.validate.coords(coords)
  edges <- as.matrix(edges)
  if (nrow(edges) == 0L) {
    return(NA_real_)
  }
  diffs <- coords[edges[, 1L], , drop = FALSE] - coords[edges[, 2L], , drop = FALSE]
  lens <- sqrt(rowSums(diffs^2))
  keep <- is.finite(lens) & lens > 0
  if (!any(keep)) {
    return(NA_real_)
  }
  unit <- diffs[keep, , drop = FALSE] / lens[keep]
  mean(apply(abs(unit), 1L, max))
}

grip.carpet.metadata <- function(target.coords) {
  target.coords <- grip.validate.coords(target.coords)
  if (ncol(target.coords) != 2L) {
    stop("carpet diagnostics require 2D target coordinates")
  }

  side <- as.integer(round(max(target.coords[, 1L]) + 0.5))
  cell.x <- as.integer(round(target.coords[, 1L] - 0.5))
  cell.y <- as.integer(round((side - 0.5) - target.coords[, 2L]))
  if (any(cell.x < 0L | cell.x >= side | cell.y < 0L | cell.y >= side)) {
    stop("target.coords do not look like Sierpinski carpet cell centers")
  }

  id.map <- matrix(0L, nrow = side, ncol = side)
  for (i in seq_len(nrow(target.coords))) {
    id.map[cell.x[[i]] + 1L, cell.y[[i]] + 1L] <- i
  }
  missing <- id.map == 0L
  visited <- matrix(FALSE, nrow = side, ncol = side)
  holes <- list()

  norm.meta <- grip.normalize.coords.with.meta(target.coords)
  normalize.point <- function(pt) {
    (pt - norm.meta$center) / norm.meta$radius
  }

  for (x in seq_len(side)) {
    for (y in seq_len(side)) {
      if (!missing[x, y] || visited[x, y]) {
        next
      }
      queue <- matrix(c(x, y), ncol = 2L)
      visited[x, y] <- TRUE
      head <- 1L
      component <- queue
      adjacent <- integer(0L)
      corridor.refs <- list()

      while (head <= nrow(queue)) {
        cur.x <- queue[head, 1L]
        cur.y <- queue[head, 2L]
        head <- head + 1L
        nbrs <- rbind(
          c(cur.x - 1L, cur.y),
          c(cur.x + 1L, cur.y),
          c(cur.x, cur.y - 1L),
          c(cur.x, cur.y + 1L)
        )
        for (k in seq_len(nrow(nbrs))) {
          nx <- nbrs[k, 1L]
          ny <- nbrs[k, 2L]
          if (nx < 1L || nx > side || ny < 1L || ny > side) {
            next
          }
          if (missing[nx, ny]) {
            if (!visited[nx, ny]) {
              visited[nx, ny] <- TRUE
              queue <- rbind(queue, c(nx, ny))
              component <- rbind(component, c(nx, ny))
            }
          } else {
            vid <- id.map[nx, ny]
            adjacent <- c(adjacent, vid)
            axis <- if (nx != cur.x) "x" else "y"
            ref.val <- if (axis == "x") target.coords[vid, 1L] else target.coords[vid, 2L]
            corridor.refs[[length(corridor.refs) + 1L]] <- list(vertex = vid, axis = axis, ref = ref.val)
          }
        }
      }

      hole.center <- c(
        mean(component[, 1L] - 0.5),
        mean((side + 0.5) - component[, 2L])
      )
      holes[[length(holes) + 1L]] <- list(
        size = nrow(component),
        center = normalize.point(hole.center),
        boundary.vertices = unique(adjacent),
        corridor.refs = corridor.refs
      )
    }
  }

  list(
    target.normalized = norm.meta$normalized,
    outer.boundary = list(
      left = which(cell.x == 0L),
      right = which(cell.x == side - 1L),
      top = which(cell.y == 0L),
      bottom = which(cell.y == side - 1L)
    ),
    holes = holes
  )
}

grip.carpet.central.hole.metrics <- function(aligned, target.norm, holes) {
  empty <- list(
    central.hole.skew = NA_real_,
    central.hole.aspect.error = NA_real_,
    central.hole.center.error = NA_real_
  )
  if (length(holes) == 0L) {
    return(empty)
  }

  hole.sizes <- vapply(holes, function(hole) hole$size, numeric(1L))
  hole <- holes[[which.max(hole.sizes)]]
  ids <- hole$boundary.vertices
  if (length(ids) == 0L) {
    return(empty)
  }

  target.boundary <- target.norm[ids, , drop = FALSE]
  x.min <- min(target.boundary[, 1L])
  x.max <- max(target.boundary[, 1L])
  y.min <- min(target.boundary[, 2L])
  y.max <- max(target.boundary[, 2L])
  width.target <- x.max - x.min
  height.target <- y.max - y.min
  scale <- mean(c(width.target, height.target))
  if (!is.finite(scale) || scale <= 0) {
    scale <- max(width.target, height.target, 1)
  }

  tol.x <- max(width.target * 0.05, 1e-8)
  tol.y <- max(height.target * 0.05, 1e-8)
  left.ids <- ids[abs(target.boundary[, 1L] - x.min) <= tol.x]
  right.ids <- ids[abs(target.boundary[, 1L] - x.max) <= tol.x]
  bottom.ids <- ids[abs(target.boundary[, 2L] - y.min) <= tol.y]
  top.ids <- ids[abs(target.boundary[, 2L] - y.max) <= tol.y]
  if (length(left.ids) == 0L || length(right.ids) == 0L ||
      length(bottom.ids) == 0L || length(top.ids) == 0L) {
    return(empty)
  }

  top.x <- mean(aligned[top.ids, 1L])
  bottom.x <- mean(aligned[bottom.ids, 1L])
  left.y <- mean(aligned[left.ids, 2L])
  right.y <- mean(aligned[right.ids, 2L])
  width.est <- mean(aligned[right.ids, 1L]) - mean(aligned[left.ids, 1L])
  height.est <- mean(aligned[top.ids, 2L]) - mean(aligned[bottom.ids, 2L])
  center.est <- colMeans(aligned[ids, , drop = FALSE])
  target.center <- colMeans(target.boundary)

  list(
    central.hole.skew = sqrt((top.x - bottom.x)^2 + (left.y - right.y)^2) / scale,
    central.hole.aspect.error = abs(abs(width.est) - abs(height.est)) / scale,
    central.hole.center.error = sqrt(sum((center.est - target.center)^2))
  )
}

grip.carpet.diagnostics <- function(coords, target.coords) {
  fit <- grip.align.to.target.nd(coords, target.coords, allow.reflection = TRUE)
  meta <- grip.carpet.metadata(target.coords)
  target.norm <- meta$target.normalized
  aligned <- fit$aligned

  side.dev <- function(ids, axis) {
    if (length(ids) == 0L) {
      return(numeric(0L))
    }
    abs(aligned[ids, axis] - target.norm[ids, axis])
  }

  boundary.devs <- c(
    side.dev(meta$outer.boundary$left, 1L),
    side.dev(meta$outer.boundary$right, 1L),
    side.dev(meta$outer.boundary$top, 2L),
    side.dev(meta$outer.boundary$bottom, 2L)
  )

  corridor.devs <- numeric(0L)
  hole.center.err <- numeric(0L)
  for (hole in meta$holes) {
    refs <- hole$corridor.refs
    if (length(refs) > 0L) {
      for (ref in refs) {
        axis <- if (identical(ref$axis, "x")) 1L else 2L
        norm.ref <- if (axis == 1L) {
          target.norm[ref$vertex, 1L]
        } else {
          target.norm[ref$vertex, 2L]
        }
        corridor.devs <- c(corridor.devs, abs(aligned[ref$vertex, axis] - norm.ref))
      }
    }
    vids <- hole$boundary.vertices
    if (length(vids) > 0L) {
      center.est <- colMeans(aligned[vids, , drop = FALSE])
      hole.center.err <- c(hole.center.err, sqrt(sum((center.est - hole$center)^2)))
    }
  }

  central <- grip.carpet.central.hole.metrics(aligned, target.norm, meta$holes)

  list(
    procrustes.rmse = fit$rmse,
    boundary.waviness = if (length(boundary.devs) > 0L) mean(boundary.devs) else NA_real_,
    corridor.waviness = if (length(corridor.devs) > 0L) mean(corridor.devs) else NA_real_,
    hole.center.error = if (length(hole.center.err) > 0L) mean(hole.center.err) else NA_real_,
    central.hole.skew = central$central.hole.skew,
    central.hole.aspect.error = central$central.hole.aspect.error,
    central.hole.center.error = central$central.hole.center.error,
    aligned = fit$aligned,
    target = fit$target
  )
}

#' Geometry-aware diagnostics against a canonical target
#'
#' \code{grip.geometry.diagnostics()} augments the graph-aware quality measures
#' in \code{\link{grip.score.layout}()} with target-aware geometric diagnostics.
#' The function aligns \code{coords} to \code{target.coords} using an orthogonal
#' Procrustes fit, then reports global symmetry, local angle preservation,
#' edge-axis concentration, and, for Sierpinski carpet layouts, boundary,
#' corridor, and hole-center diagnostics.
#'
#' Metrics are reported so that larger \code{global.symmetry.score} and
#' \code{edge.axis.concentration} are better, while smaller
#' \code{procrustes.rmse}, \code{local.angle.deviation},
#' \code{boundary.waviness}, \code{corridor.waviness},
#' \code{hole.center.error}, \code{central.hole.skew},
#' \code{central.hole.aspect.error}, and \code{central.hole.center.error} are
#' better.
#'
#' @param coords Numeric layout matrix with 2 or 3 columns.
#' @param target.coords Canonical target coordinates with the same shape as
#'   \code{coords}.
#' @param edges Two-column integer edge matrix.
#' @param family Optional graph-family label. Use \code{"sierpinski.carpet"} to
#'   enable carpet-specific diagnostics.
#' @param sample.size.symmetry Number of vertices sampled when evaluating global
#'   symmetry.
#' @param sample.size.wedges Number of wedges sampled when evaluating local
#'   angle deviation.
#' @param rng.seed Integer seed used for the symmetry and wedge sampling.
#'
#' @return A one-row data frame of geometric diagnostics.
#' @export
grip.geometry.diagnostics <- function(coords,
                                      target.coords,
                                      edges,
                                      family = NULL,
                                      sample.size.symmetry = 512L,
                                      sample.size.wedges = 4000L,
                                      rng.seed = 1L) {
  coords <- grip.validate.coords(coords)
  target.coords <- grip.validate.coords(target.coords)
  if (!identical(dim(coords), dim(target.coords))) {
    stop("coords and target.coords must have the same dimensions")
  }
  edges <- as.matrix(edges)
  if (ncol(edges) != 2L) {
    stop("edges must be a two-column matrix")
  }

  fit <- grip.align.to.target.nd(coords, target.coords, allow.reflection = TRUE)
  axis.conc <- grip.edge.axis.concentration(fit$aligned, edges)
  out <- data.frame(
    procrustes.rmse = fit$rmse,
    global.symmetry.score = grip.global.symmetry.score(fit$aligned,
                                                       sample.size = sample.size.symmetry,
                                                       rng.seed = rng.seed),
    local.angle.deviation = grip.local.angle.deviation(coords,
                                                       target.coords,
                                                       edges,
                                                       sample.size = sample.size.wedges,
                                                       rng.seed = rng.seed + 1000L),
    edge.axis.concentration = axis.conc,
    edge.axis.deviation = if (is.finite(axis.conc)) 1 - axis.conc else NA_real_,
    boundary.waviness = NA_real_,
    corridor.waviness = NA_real_,
    hole.center.error = NA_real_,
    central.hole.skew = NA_real_,
    central.hole.aspect.error = NA_real_,
    central.hole.center.error = NA_real_,
    stringsAsFactors = FALSE
  )

  if (identical(family, "sierpinski.carpet")) {
    carpet <- grip.carpet.diagnostics(coords, target.coords)
    out$boundary.waviness[[1L]] <- carpet$boundary.waviness
    out$corridor.waviness[[1L]] <- carpet$corridor.waviness
    out$hole.center.error[[1L]] <- carpet$hole.center.error
    out$central.hole.skew[[1L]] <- carpet$central.hole.skew
    out$central.hole.aspect.error[[1L]] <- carpet$central.hole.aspect.error
    out$central.hole.center.error[[1L]] <- carpet$central.hole.center.error
  }

  out
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

grip.sort.adj.with.weights <- function(adj.list, weight.list = NULL) {
  out.adj <- vector("list", length(adj.list))
  out.w <- if (is.null(weight.list)) NULL else vector("list", length(adj.list))
  for (i in seq_along(adj.list)) {
    nb <- as.integer(adj.list[[i]])
    if (length(nb) == 0L) {
      out.adj[[i]] <- integer(0L)
      if (!is.null(out.w)) out.w[[i]] <- numeric(0L)
      next
    }
    ord <- order(nb)
    out.adj[[i]] <- nb[ord]
    if (!is.null(out.w)) {
      out.w[[i]] <- as.double(weight.list[[i]][ord])
    }
  }
  list(adj_list = out.adj, weight_list = out.w)
}

grip.bfs.tree <- function(adj.list, source) {
  n <- length(adj.list)
  dist <- rep.int(-1L, n)
  parent <- integer(n)
  q <- integer(n)
  head <- 1L
  tail <- 1L
  q[[tail]] <- source
  dist[[source]] <- 0L
  while (head <= tail) {
    v <- q[[head]]
    head <- head + 1L
    nb <- adj.list[[v]]
    if (length(nb) == 0L) next
    for (u in nb) {
      cand.dist <- dist[[v]] + 1L
      if (dist[[u]] == -1L) {
        dist[[u]] <- cand.dist
        parent[[u]] <- v
        tail <- tail + 1L
        q[[tail]] <- u
      } else if (dist[[u]] == cand.dist && (parent[[u]] == 0L || v < parent[[u]])) {
        parent[[u]] <- v
      }
    }
  }
  list(dist = as.double(dist), parent = as.integer(parent))
}

grip.dijkstra.tree <- function(adj.list, weight.list, source) {
  n <- length(adj.list)
  dist <- rep.int(Inf, n)
  parent <- integer(n)
  visited <- rep.int(FALSE, n)
  dist[[source]] <- 0
  tol <- sqrt(.Machine$double.eps)

  for (iter in seq_len(n)) {
    remaining <- which(!visited)
    if (length(remaining) == 0L) break
    current.dist <- dist[remaining]
    v <- remaining[[which.min(current.dist)]]
    if (!is.finite(dist[[v]])) break
    visited[[v]] <- TRUE
    nb <- adj.list[[v]]
    if (length(nb) == 0L) next
    ww <- weight.list[[v]]
    for (k in seq_along(nb)) {
      u <- nb[[k]]
      alt <- dist[[v]] + ww[[k]]
      best <- dist[[u]]
      equal.dist <- is.finite(best) &&
        abs(alt - best) <= tol * max(1, abs(alt), abs(best))
      if (alt + tol < best || (equal.dist && (parent[[u]] == 0L || v < parent[[u]]))) {
        dist[[u]] <- alt
        parent[[u]] <- v
      }
    }
  }

  list(dist = as.double(dist), parent = as.integer(parent))
}

grip.shortest.path.tree <- function(adj.list, weight.list, source) {
  if (is.null(weight.list)) {
    return(grip.bfs.tree(adj.list, source))
  }
  grip.dijkstra.tree(adj.list, weight.list, source)
}

grip.closest.active.vertices <- function(dist.row, source, count) {
  count <- as.integer(count)
  if (is.na(count) || count <= 0L) {
    return(integer(0L))
  }
  ids <- seq_along(dist.row)
  keep <- ids != source & is.finite(dist.row)
  if (!any(keep)) {
    return(integer(0L))
  }
  ids <- ids[keep]
  vals <- dist.row[keep]
  ord <- order(vals, ids)
  ids[ord][seq_len(min(count, length(ids)))]
}

grip.farthest.landmarks <- function(source, dist.matrix, count) {
  count <- as.integer(count)
  if (is.na(count) || count <= 0L) {
    return(integer(0L))
  }
  n <- nrow(dist.matrix)
  ids <- seq_len(n)
  keep <- ids != source & is.finite(dist.matrix[source, ])
  candidates <- ids[keep]
  if (length(candidates) == 0L) {
    return(integer(0L))
  }

  selected <- integer(0L)
  coverage <- dist.matrix[source, candidates]
  for (step in seq_len(min(count, length(candidates)))) {
    if (step == 1L) {
      scores <- dist.matrix[source, candidates]
    } else {
      scores <- coverage
    }
    ord <- order(-scores, candidates)
    choice <- candidates[ord[[1L]]]
    selected <- c(selected, choice)
    keep.next <- candidates != choice
    candidates <- candidates[keep.next]
    if (length(candidates) == 0L) {
      break
    }
    coverage <- pmin(coverage[keep.next], dist.matrix[choice, candidates])
  }
  as.integer(selected)
}

grip.reconstruct.path.vertices <- function(parent, source, target) {
  if (source == target) {
    return(as.integer(source))
  }
  path <- integer(0L)
  cur <- as.integer(target)
  seen <- logical(length(parent))
  while (cur != 0L && cur != source) {
    if (seen[[cur]]) {
      stop("detected a cycle while reconstructing a chosen shortest path")
    }
    seen[[cur]] <- TRUE
    path <- c(cur, path)
    cur <- parent[[cur]]
  }
  if (cur != source) {
    stop(sprintf("could not reconstruct a shortest path from %d to %d", source, target))
  }
  as.integer(c(source, path))
}

grip.path.euclidean.length <- function(coords, vertices, edge_length_epsilon = 1e-8) {
  if (length(vertices) <= 1L) {
    return(0)
  }
  diffs <- coords[vertices[-1L], , drop = FALSE] - coords[vertices[-length(vertices)], , drop = FALSE]
  sum(sqrt(rowSums(diffs^2) + edge_length_epsilon^2))
}

grip.path.vertices.to.edges <- function(vertices) {
  if (length(vertices) <= 1L) {
    return(matrix(integer(), ncol = 2L))
  }
  cbind(
    vertices[-length(vertices)],
    vertices[-1L]
  )
}

grip.validate.scalar <- function(x,
                                 name,
                                 lower = -Inf,
                                 upper = Inf,
                                 open.lower = FALSE,
                                 open.upper = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop(sprintf("%s must be a single finite numeric value", name))
  }
  if ((open.lower && x <= lower) || (!open.lower && x < lower) ||
      (open.upper && x >= upper) || (!open.upper && x > upper)) {
    lower.txt <- if (is.finite(lower)) {
      if (open.lower) sprintf("(%s", format(lower, digits = 16)) else sprintf("[%s", format(lower, digits = 16))
    } else {
      "(-Inf"
    }
    upper.txt <- if (is.finite(upper)) {
      if (open.upper) sprintf("%s)", format(upper, digits = 16)) else sprintf("%s]", format(upper, digits = 16))
    } else {
      "Inf)"
    }
    stop(sprintf("%s must lie in %s, %s", name, lower.txt, upper.txt))
  }
}

grip.validate.count <- function(x, name) {
  grip.validate.scalar(x, name, lower = 0)
  x <- as.integer(round(x))
  if (is.na(x) || x < 0L) {
    stop(sprintf("%s must be >= 0", name))
  }
  x
}

grip.validate.prepared.object <- function(prepared,
                                          class_name,
                                          prepare_fun_name,
                                          coords = NULL) {
  if (!inherits(prepared, class_name)) {
    stop(sprintf("prepared must be NULL or an object from %s()", prepare_fun_name))
  }
  if (!is.null(coords) && nrow(coords) != prepared$n) {
    stop("nrow(coords) must match the graph size stored in prepared")
  }
  prepared
}

grip.prepare.geodesic.kk.base <- function(edges = NULL,
                                          n = NULL,
                                          adj_list = NULL,
                                          weight_list = NULL,
                                          edge_weights = NULL,
                                          caller = "grip.prepare.geodesic.kk") {
  if (is.null(n) && is.null(adj_list) && !is.null(edges)) {
    n <- max(as.integer(edges), na.rm = TRUE)
  }
  if (is.null(n) && !is.null(adj_list)) {
    n <- length(adj_list)
  }
  if (is.null(n) || !is.finite(n) || n <= 0L) {
    stop("n must be provided or inferable from edges/adj_list")
  }
  n <- as.integer(n)

  validated <- grip.validate.layout.inputs(
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights,
    dim = 2L,
    placement = "barycenter",
    seed = 1L
  )

  sorted <- grip.sort.adj.with.weights(validated$adj_list, validated$weight_list)
  adj.list <- sorted$adj_list
  weight.list <- sorted$weight_list
  edges <- grip.edges.from.adj.list(adj.list)

  trees <- lapply(seq_len(validated$n), function(source) {
    grip.shortest.path.tree(adj.list, weight.list, source)
  })
  dist.matrix <- do.call(rbind, lapply(trees, `[[`, "dist"))
  finite.mask <- row(dist.matrix) != col(dist.matrix)
  if (any(!is.finite(dist.matrix[finite.mask]))) {
    stop(sprintf("%s() currently requires a connected graph", caller))
  }

  list(
    n = validated$n,
    edges = edges,
    adj_list = adj.list,
    weight_list = weight.list,
    trees = trees,
    parents = lapply(trees, `[[`, "parent"),
    distance_matrix = dist.matrix,
    graph_diameter = max(dist.matrix)
  )
}

grip.landmark.geodesic.kk.pair.matrix <- function(dist.matrix,
                                                  local_nbrs,
                                                  landmark_count) {
  n <- nrow(dist.matrix)
  pair.keys <- character(0L)

  for (source in seq_len(n)) {
    local <- grip.closest.active.vertices(dist.matrix[source, ], source, local_nbrs)
    landmarks <- grip.farthest.landmarks(source, dist.matrix, landmark_count)
    chosen <- unique(c(local, landmarks))
    chosen <- chosen[chosen != source]
    if (length(chosen) == 0L) {
      next
    }
    pairs <- cbind(pmin(source, chosen), pmax(source, chosen))
    pairs <- unique(pairs)
    pair.keys <- c(pair.keys, paste(pairs[, 1L], pairs[, 2L], sep = "-"))
  }

  if (length(pair.keys) == 0L) {
    return(matrix(integer(), ncol = 2L))
  }

  pair.matrix <- do.call(rbind, strsplit(unique(pair.keys), "-", fixed = TRUE))
  pair.matrix <- matrix(as.integer(pair.matrix), ncol = 2L)
  pair.matrix[order(pair.matrix[, 1L], pair.matrix[, 2L]), , drop = FALSE]
}

grip.full.geodesic.kk.pair.matrix <- function(n) {
  n <- as.integer(n)
  if (is.na(n) || n < 2L) {
    return(matrix(integer(), ncol = 2L))
  }
  pair.matrix <- t(utils::combn(seq_len(n), 2L))
  matrix(as.integer(pair.matrix), ncol = 2L)
}

grip.build.geodesic.kk.path.cache <- function(pair.matrix,
                                              parents,
                                              dist.matrix) {
  n.pairs <- nrow(pair.matrix)
  path.vertices <- vector("list", n.pairs)
  path.edges <- vector("list", n.pairs)
  pair.graph.distance <- numeric(n.pairs)

  if (n.pairs == 0L) {
    return(list(
      path_vertices = path.vertices,
      path_edges = path.edges,
      pair_graph_distance = pair.graph.distance
    ))
  }

  for (i in seq_len(n.pairs)) {
    src <- pair.matrix[i, 1L]
    dst <- pair.matrix[i, 2L]
    path.vertices[[i]] <- grip.reconstruct.path.vertices(parents[[src]], src, dst)
    path.edges[[i]] <- grip.path.vertices.to.edges(path.vertices[[i]])
    pair.graph.distance[[i]] <- dist.matrix[src, dst]
  }

  list(
    path_vertices = path.vertices,
    path_edges = path.edges,
    pair_graph_distance = as.double(pair.graph.distance)
  )
}

grip.geodesic.kk.path.lengths <- function(coords,
                                          prepared,
                                          edge_length_epsilon = 1e-8) {
  vapply(prepared$path_vertices, function(path) {
    grip.path.euclidean.length(
      coords,
      path,
      edge_length_epsilon = edge_length_epsilon
    )
  }, numeric(1L))
}

grip.geodesic.kk.fit.scale <- function(path.lengths,
                                       graph.distances,
                                       stiffness = 1.0,
                                       distance_floor = 1e-8) {
  kk <- as.double(stiffness) / pmax(as.double(graph.distances), as.double(distance_floor))^2
  denom <- sum(kk * graph.distances * graph.distances)
  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }
  sum(kk * graph.distances * path.lengths) / denom
}

grip.geodesic.kk.energy.gradient <- function(coords,
                                             prepared,
                                             scale.L0,
                                             stiffness = 1.0,
                                             distance_floor = 1e-8,
                                             edge_length_epsilon = 1e-8) {
  g <- as.double(prepared$pair_graph_distance)
  n.pairs <- length(g)
  kk <- as.double(stiffness) / pmax(g, as.double(distance_floor))^2
  target <- as.double(scale.L0) * g
  grad <- matrix(0, nrow = nrow(coords), ncol = ncol(coords))
  energy <- 0
  path.lengths <- numeric(n.pairs)

  if (n.pairs == 0L) {
    return(list(
      energy = energy,
      gradient = grad,
      gradient_norm = 0,
      path_lengths = path.lengths,
      target = target
    ))
  }

  for (i in seq_len(n.pairs)) {
    edges <- prepared$path_edges[[i]]
    if (nrow(edges) == 0L) {
      next
    }
    diffs <- coords[edges[, 1L], , drop = FALSE] - coords[edges[, 2L], , drop = FALSE]
    edge.lengths <- sqrt(rowSums(diffs^2) + edge_length_epsilon^2)
    h <- sum(edge.lengths)
    path.lengths[[i]] <- h
    resid <- h - target[[i]]
    coeff <- kk[[i]] * resid
    energy <- energy + 0.5 * kk[[i]] * resid^2
    unit.vecs <- diffs / edge.lengths
    for (j in seq_len(nrow(edges))) {
      u <- edges[j, 1L]
      v <- edges[j, 2L]
      grad[u, ] <- grad[u, ] + coeff * unit.vecs[j, ]
      grad[v, ] <- grad[v, ] - coeff * unit.vecs[j, ]
    }
  }

  list(
    energy = energy,
    gradient = grad,
    gradient_norm = sqrt(sum(grad^2)),
    path_lengths = path.lengths,
    target = target
  )
}

grip.lgkk.path.lengths <- function(coords,
                                   prepared,
                                   edge_length_epsilon = 1e-8) {
  grip.geodesic.kk.path.lengths(
    coords = coords,
    prepared = prepared,
    edge_length_epsilon = edge_length_epsilon
  )
}

grip.lgkk.fit.scale <- function(path.lengths,
                                graph.distances,
                                stiffness = 1.0,
                                distance_floor = 1e-8) {
  grip.geodesic.kk.fit.scale(
    path.lengths = path.lengths,
    graph.distances = graph.distances,
    stiffness = stiffness,
    distance_floor = distance_floor
  )
}

grip.lgkk.energy.gradient <- function(coords,
                                      prepared,
                                      scale.L0,
                                      stiffness = 1.0,
                                      distance_floor = 1e-8,
                                      edge_length_epsilon = 1e-8) {
  grip.geodesic.kk.energy.gradient(
    coords = coords,
    prepared = prepared,
    scale.L0 = scale.L0,
    stiffness = stiffness,
    distance_floor = distance_floor,
    edge_length_epsilon = edge_length_epsilon
  )
}

grip.geodesic.kk.score.stats <- function(coords,
                                         prepared,
                                         stiffness = 1.0,
                                         distance_floor = 1e-8,
                                         edge_length_epsilon = 1e-8,
                                         scale_mode = c("profiled", "user"),
                                         scale.L0 = NULL) {
  scale_mode <- match.arg(scale_mode)
  grip.validate.scalar(stiffness, "stiffness", lower = 0, open.lower = TRUE)
  grip.validate.scalar(distance_floor, "distance_floor", lower = 0, open.lower = TRUE)
  grip.validate.scalar(edge_length_epsilon, "edge_length_epsilon", lower = 0)
  if (identical(scale_mode, "user")) {
    grip.validate.scalar(scale.L0, "scale.L0", lower = 0, open.lower = TRUE)
  }

  g <- as.double(prepared$pair_graph_distance)
  kk <- as.double(stiffness) / pmax(g, as.double(distance_floor))^2

  if (length(g) == 0L) {
    return(list(
      n.pairs = 0L,
      scale_mode = scale_mode,
      scale.L0 = if (identical(scale_mode, "user")) as.double(scale.L0) else NA_real_,
      scale.L = NA_real_,
      energy = NA_real_,
      weighted.rmse = NA_real_,
      weighted.rel.rmse = NA_real_,
      mean.abs.path.error = NA_real_,
      mean.rel.path.error = NA_real_,
      path.lengths = numeric(0L),
      target = numeric(0L),
      residual = numeric(0L),
      relative.residual = numeric(0L),
      stiffnesses = kk
    ))
  }

  h <- grip.geodesic.kk.path.lengths(
    coords = coords,
    prepared = prepared,
    edge_length_epsilon = edge_length_epsilon
  )
  if (identical(scale_mode, "profiled")) {
    scale.L0 <- grip.geodesic.kk.fit.scale(
      path.lengths = h,
      graph.distances = g,
      stiffness = stiffness,
      distance_floor = distance_floor
    )
  }
  if (!is.finite(scale.L0)) {
    stop("failed to fit a geodesic KK scale")
  }

  target <- scale.L0 * g
  resid <- h - target
  rel.resid <- resid / pmax(target, distance_floor)

  list(
    n.pairs = length(g),
    scale_mode = scale_mode,
    scale.L0 = scale.L0,
    scale.L = scale.L0 * prepared$graph_diameter,
    energy = 0.5 * sum(kk * resid^2),
    weighted.rmse = sqrt(sum(kk * resid^2) / sum(kk)),
    weighted.rel.rmse = sqrt(sum(kk * rel.resid^2) / sum(kk)),
    mean.abs.path.error = mean(abs(resid)),
    mean.rel.path.error = mean(abs(rel.resid)),
    path.lengths = h,
    target = target,
    residual = resid,
    relative.residual = rel.resid,
    stiffnesses = kk
  )
}

grip.geodesic.kk.pair.details <- function(prepared, stats) {
  data.frame(
    i = prepared$pair_matrix[, 1L],
    j = prepared$pair_matrix[, 2L],
    graph.distance = as.double(prepared$pair_graph_distance),
    embedded.path.length = stats$path.lengths,
    target.length = stats$target,
    residual = stats$residual,
    relative.residual = stats$relative.residual,
    stiffness = stats$stiffnesses,
    stringsAsFactors = FALSE
  )
}

grip.geodesic.kk.evaluate.state <- function(coords,
                                            prepared,
                                            stiffness = 1.0,
                                            distance_floor = 1e-8,
                                            edge_length_epsilon = 1e-8,
                                            scale_mode = c("fixed", "profiled"),
                                            scale.L0 = NULL) {
  scale_mode <- match.arg(scale_mode)
  if (identical(scale_mode, "profiled")) {
    scale.L0 <- grip.geodesic.kk.fit.scale(
      path.lengths = grip.geodesic.kk.path.lengths(
        coords = coords,
        prepared = prepared,
        edge_length_epsilon = edge_length_epsilon
      ),
      graph.distances = prepared$pair_graph_distance,
      stiffness = stiffness,
      distance_floor = distance_floor
    )
    if (!is.finite(scale.L0)) {
      stop("failed to fit a geodesic KK scale")
    }
  } else {
    grip.validate.scalar(scale.L0, "scale.L0", lower = 0, open.lower = TRUE)
  }

  state <- grip.geodesic.kk.energy.gradient(
    coords = coords,
    prepared = prepared,
    scale.L0 = scale.L0,
    stiffness = stiffness,
    distance_floor = distance_floor,
    edge_length_epsilon = edge_length_epsilon
  )
  state$scale.L0 <- scale.L0
  state
}

#' Prepare sparse landmark-geodesic KK data for repeated layout evaluation
#'
#' \code{grip.prepare.landmark.geodesic.kk()} builds the deterministic shortest
#' path trees, graph-distance cache, sparse local-plus-landmark pair set, and
#' chosen path realizations needed to evaluate the landmark geodesic KK energy
#' repeatedly on the same graph. This is intended as a reusable preparation step
#' for experiments where many layouts of the same graph are compared.
#'
#' The sparse set follows the implementation choice recorded in
#' \code{landmark\_geodesic\_kk\_spec\_2026-03-30.tex}: each vertex contributes
#' its \code{local_nbrs} nearest active vertices in graph distance and its
#' \code{landmark_count} farthest-point landmarks, both selected
#' deterministically.
#'
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#' @param n Number of vertices. If omitted with \code{adj_list}, defaults to
#'   \code{length(adj_list)}. If omitted with \code{edges}, defaults to
#'   \code{max(edges)}.
#' @param adj_list Adjacency list (1-based) for an undirected graph.
#' @param weight_list Optional parallel list of positive edge weights.
#' @param edge_weights Optional positive edge-weight vector parallel to
#'   \code{edges}.
#' @param local_nbrs Number of nearest graph-metric neighbors retained per
#'   vertex.
#' @param landmark_count Number of farthest-point landmarks retained per vertex.
#'
#' @return A list with the sparse pair set, graph distances, chosen paths, and
#'   other cached data. The object has class \code{"grip_lgkk_prepared"}.
#' @export
grip.prepare.landmark.geodesic.kk <- function(edges = NULL,
                                              n = NULL,
                                              adj_list = NULL,
                                              weight_list = NULL,
                                              edge_weights = NULL,
                                              local_nbrs = 20L,
                                              landmark_count = 8L) {
  local_nbrs <- grip.validate.count(local_nbrs, "local_nbrs")
  landmark_count <- grip.validate.count(landmark_count, "landmark_count")

  base <- grip.prepare.geodesic.kk.base(
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights,
    caller = "grip.prepare.landmark.geodesic.kk"
  )
  pair.matrix <- grip.landmark.geodesic.kk.pair.matrix(
    dist.matrix = base$distance_matrix,
    local_nbrs = local_nbrs,
    landmark_count = landmark_count
  )
  cache <- grip.build.geodesic.kk.path.cache(
    pair.matrix = pair.matrix,
    parents = base$parents,
    dist.matrix = base$distance_matrix
  )

  out <- list(
    n = base$n,
    edges = base$edges,
    adj_list = base$adj_list,
    weight_list = base$weight_list,
    local_nbrs = local_nbrs,
    landmark_count = landmark_count,
    pair_matrix = pair.matrix,
    pair_graph_distance = cache$pair_graph_distance,
    path_vertices = cache$path_vertices,
    path_edges = cache$path_edges,
    graph_diameter = base$graph_diameter,
    distance_matrix = base$distance_matrix,
    pair_mode = "landmark_sparse"
  )
  class(out) <- c("grip_lgkk_prepared", "grip_geodesic_kk_prepared")
  out
}

#' Prepare full geodesic KK data for repeated layout evaluation
#'
#' \code{grip.prepare.geodesic.kk()} builds the deterministic all-pairs shortest
#' path cache needed to evaluate or optimize the full geodesic Kamada--Kawai
#' objective repeatedly on the same connected graph.
#'
#' For each unordered vertex pair, the prepared object stores one deterministic
#' chosen graph shortest path, its graph distance, and the corresponding path
#' edge sequence. This is the full all-pairs analogue of the sparse landmark
#' cache used by \code{\link{grip.prepare.landmark.geodesic.kk}()}.
#'
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#' @param n Number of vertices. If omitted with \code{adj_list}, defaults to
#'   \code{length(adj_list)}. If omitted with \code{edges}, defaults to
#'   \code{max(edges)}.
#' @param adj_list Adjacency list (1-based) for an undirected graph.
#' @param weight_list Optional parallel list of positive edge weights.
#' @param edge_weights Optional positive edge-weight vector parallel to
#'   \code{edges}.
#'
#' @return A list with the all-pairs graph distances, chosen paths, and cached
#'   path-edge realizations. The object has class \code{"grip_gkk_prepared"}.
#' @export
grip.prepare.geodesic.kk <- function(edges = NULL,
                                     n = NULL,
                                     adj_list = NULL,
                                     weight_list = NULL,
                                     edge_weights = NULL) {
  base <- grip.prepare.geodesic.kk.base(
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights,
    caller = "grip.prepare.geodesic.kk"
  )
  pair.matrix <- grip.full.geodesic.kk.pair.matrix(base$n)
  cache <- grip.build.geodesic.kk.path.cache(
    pair.matrix = pair.matrix,
    parents = base$parents,
    dist.matrix = base$distance_matrix
  )

  out <- list(
    n = base$n,
    edges = base$edges,
    adj_list = base$adj_list,
    weight_list = base$weight_list,
    pair_matrix = pair.matrix,
    pair_graph_distance = cache$pair_graph_distance,
    path_vertices = cache$path_vertices,
    path_edges = cache$path_edges,
    graph_diameter = base$graph_diameter,
    distance_matrix = base$distance_matrix,
    pair_mode = "all_pairs"
  )
  class(out) <- c("grip_gkk_prepared", "grip_geodesic_kk_prepared")
  out
}

#' Score a layout under the landmark geodesic KK energy
#'
#' \code{grip.score.landmark.geodesic.kk()} evaluates a layout using the
#' landmark geodesic KK objective described in
#' \code{landmark\_geodesic\_kk\_spec\_2026-03-30.tex}. Distances are measured
#' along fixed chosen graph shortest paths, not by straight-line chord length.
#' The target scale factor \code{L0} is fit analytically for the supplied layout
#' so that rankings are not dominated by an arbitrary global drawing scale.
#'
#' This is a scoring and comparison helper, not an optimizer.
#'
#' @param coords Numeric coordinate matrix with 2 or 3 columns.
#' @param prepared Optional object returned by
#'   \code{\link{grip.prepare.landmark.geodesic.kk}()}.
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#' @param n Number of vertices.
#' @param adj_list Adjacency list (1-based) for an undirected graph.
#' @param weight_list Optional parallel list of positive edge weights.
#' @param edge_weights Optional positive edge-weight vector parallel to
#'   \code{edges}.
#' @param local_nbrs Number of nearest graph-metric neighbors retained per
#'   vertex when \code{prepared} is not supplied.
#' @param landmark_count Number of farthest-point landmarks retained per vertex
#'   when \code{prepared} is not supplied.
#' @param stiffness Global stiffness constant \(K\).
#' @param distance_floor Small positive floor used in
#'   \code{k_ij = K / max(g_ij, distance_floor)^2}.
#' @param edge_length_epsilon Small positive stabilizer added inside each
#'   embedded edge length.
#' @param return_pair_details If \code{TRUE}, include pairwise path lengths and
#'   residuals in a list column.
#'
#' @return A one-row data frame with the fitted scale factor and landmark
#'   geodesic KK energy summary.
#' @export
grip.score.landmark.geodesic.kk <- function(coords,
                                            prepared = NULL,
                                            edges = NULL,
                                            n = NULL,
                                            adj_list = NULL,
                                            weight_list = NULL,
                                            edge_weights = NULL,
                                            local_nbrs = 20L,
                                            landmark_count = 8L,
                                            stiffness = 1.0,
                                            distance_floor = 1e-8,
                                            edge_length_epsilon = 1e-8,
                                            return_pair_details = FALSE) {
  coords <- grip.validate.coords(coords)
  if (is.null(prepared)) {
    prepared <- grip.prepare.landmark.geodesic.kk(
      edges = edges,
      n = if (is.null(n)) nrow(coords) else n,
      adj_list = adj_list,
      weight_list = weight_list,
      edge_weights = edge_weights,
      local_nbrs = local_nbrs,
      landmark_count = landmark_count
    )
  }
  prepared <- grip.validate.prepared.object(
    prepared = prepared,
    class_name = "grip_lgkk_prepared",
    prepare_fun_name = "grip.prepare.landmark.geodesic.kk",
    coords = coords
  )

  stats <- grip.geodesic.kk.score.stats(
    coords = coords,
    prepared = prepared,
    stiffness = stiffness,
    distance_floor = distance_floor,
    edge_length_epsilon = edge_length_epsilon,
    scale_mode = "profiled"
  )
  if (stats$n.pairs == 0L) {
    out <- data.frame(
      n.vertices = prepared$n,
      n.pairs = 0L,
      local.nbrs = prepared$local_nbrs,
      landmark.count = prepared$landmark_count,
      scale.L0 = NA_real_,
      scale.L = NA_real_,
      lgkk.energy = NA_real_,
      lgkk.weighted.rmse = NA_real_,
      lgkk.weighted.rel.rmse = NA_real_,
      lgkk.mean.abs.path.error = NA_real_,
      lgkk.mean.rel.path.error = NA_real_,
      stringsAsFactors = FALSE
    )
    if (isTRUE(return_pair_details)) {
      out$pair.details <- list(data.frame())
    }
    return(out)
  }

  out <- data.frame(
    n.vertices = prepared$n,
    n.pairs = stats$n.pairs,
    local.nbrs = prepared$local_nbrs,
    landmark.count = prepared$landmark_count,
    scale.L0 = stats$scale.L0,
    scale.L = stats$scale.L,
    lgkk.energy = stats$energy,
    lgkk.weighted.rmse = stats$weighted.rmse,
    lgkk.weighted.rel.rmse = stats$weighted.rel.rmse,
    lgkk.mean.abs.path.error = stats$mean.abs.path.error,
    lgkk.mean.rel.path.error = stats$mean.rel.path.error,
    stringsAsFactors = FALSE
  )

  if (isTRUE(return_pair_details)) {
    out$pair.details <- list(grip.geodesic.kk.pair.details(prepared, stats))
  }
  out
}

#' Score a layout under the full geodesic KK energy
#'
#' \code{grip.score.geodesic.kk()} evaluates a layout using the full all-pairs
#' geodesic Kamada--Kawai objective. Distances are measured along fixed chosen
#' graph shortest paths, not by straight-line chord length.
#'
#' By default the global target scale \code{L0} is fit analytically for the
#' supplied layout. Alternatively, the score can be evaluated at a user-supplied
#' scale.
#'
#' @param coords Numeric coordinate matrix with 2 or 3 columns.
#' @param prepared Optional object returned by
#'   \code{\link{grip.prepare.geodesic.kk}()}.
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#' @param n Number of vertices.
#' @param adj_list Adjacency list (1-based) for an undirected graph.
#' @param weight_list Optional parallel list of positive edge weights.
#' @param edge_weights Optional positive edge-weight vector parallel to
#'   \code{edges}.
#' @param stiffness Global stiffness constant \(K\).
#' @param distance_floor Small positive floor used in
#'   \code{k_ij = K / max(g_ij, distance_floor)^2}.
#' @param edge_length_epsilon Small positive stabilizer added inside each
#'   embedded edge length.
#' @param scale_mode Either \code{"profiled"} to fit \code{L0} analytically
#'   for the supplied layout or \code{"user"} to use \code{scale.L0}.
#' @param scale.L0 Optional user-supplied geodesic KK scale, required when
#'   \code{scale_mode = "user"}.
#' @param return_pair_details If \code{TRUE}, include pairwise path lengths and
#'   residuals in a list column.
#'
#' @return A one-row data frame with the fitted or user-supplied scale factor
#'   and full geodesic KK energy summary.
#' @export
grip.score.geodesic.kk <- function(coords,
                                   prepared = NULL,
                                   edges = NULL,
                                   n = NULL,
                                   adj_list = NULL,
                                   weight_list = NULL,
                                   edge_weights = NULL,
                                   stiffness = 1.0,
                                   distance_floor = 1e-8,
                                   edge_length_epsilon = 1e-8,
                                   scale_mode = c("profiled", "user"),
                                   scale.L0 = NULL,
                                   return_pair_details = FALSE) {
  coords <- grip.validate.coords(coords)
  scale_mode <- match.arg(scale_mode)
  if (is.null(prepared)) {
    prepared <- grip.prepare.geodesic.kk(
      edges = edges,
      n = if (is.null(n)) nrow(coords) else n,
      adj_list = adj_list,
      weight_list = weight_list,
      edge_weights = edge_weights
    )
  }
  prepared <- grip.validate.prepared.object(
    prepared = prepared,
    class_name = "grip_gkk_prepared",
    prepare_fun_name = "grip.prepare.geodesic.kk",
    coords = coords
  )

  stats <- grip.geodesic.kk.score.stats(
    coords = coords,
    prepared = prepared,
    stiffness = stiffness,
    distance_floor = distance_floor,
    edge_length_epsilon = edge_length_epsilon,
    scale_mode = scale_mode,
    scale.L0 = scale.L0
  )

  out <- data.frame(
    n.vertices = prepared$n,
    n.pairs = stats$n.pairs,
    pair.mode = if (!is.null(prepared$pair_mode)) prepared$pair_mode else "all_pairs",
    scale.mode = stats$scale_mode,
    scale.L0 = stats$scale.L0,
    scale.L = stats$scale.L,
    gkk.energy = stats$energy,
    gkk.weighted.rmse = stats$weighted.rmse,
    gkk.weighted.rel.rmse = stats$weighted.rel.rmse,
    gkk.mean.abs.path.error = stats$mean.abs.path.error,
    gkk.mean.rel.path.error = stats$mean.rel.path.error,
    stringsAsFactors = FALSE
  )

  if (isTRUE(return_pair_details)) {
    out$pair.details <- list(grip.geodesic.kk.pair.details(prepared, stats))
  }
  out
}

#' Optimize a layout under the landmark geodesic KK energy
#'
#' \code{grip.optimize.landmark.geodesic.kk()} applies a deterministic
#' warm-started gradient-descent polish under the sparse landmark geodesic KK
#' energy. This is the first experimental optimizer prototype: it starts from an
#' existing layout and refines it, rather than replacing the full multiscale
#' GRIP refinement pipeline.
#'
#' The current prototype fits the LGKK target scale \code{L0} once from the
#' starting layout and then optimizes against those fixed target path lengths.
#' That keeps the gradient simple and makes the resulting line search robust for
#' an initial implementation.
#'
#' @param coords Numeric coordinate matrix with 2 or 3 columns.
#' @param prepared Optional object returned by
#'   \code{\link{grip.prepare.landmark.geodesic.kk}()}.
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#' @param n Number of vertices.
#' @param adj_list Adjacency list (1-based) for an undirected graph.
#' @param weight_list Optional parallel list of positive edge weights.
#' @param edge_weights Optional positive edge-weight vector parallel to
#'   \code{edges}.
#' @param local_nbrs Number of nearest graph-metric neighbors retained per
#'   vertex when \code{prepared} is not supplied.
#' @param landmark_count Number of farthest-point landmarks retained per vertex
#'   when \code{prepared} is not supplied.
#' @param max_iter Maximum number of gradient-descent iterations.
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
#' @param recenter If \code{TRUE}, recenter the layout to zero mean after each
#'   accepted step.
#' @param return_trace If \code{TRUE}, include per-iteration diagnostics and the
#'   accepted intermediate coordinate frames.
#'
#' @return A list with \code{coords}, \code{trace}, \code{frames},
#'   \code{prepared}, and \code{score}.
#' @export
grip.optimize.landmark.geodesic.kk <- function(coords,
                                               prepared = NULL,
                                               edges = NULL,
                                               n = NULL,
                                               adj_list = NULL,
                                               weight_list = NULL,
                                               edge_weights = NULL,
                                               local_nbrs = 20L,
                                               landmark_count = 8L,
                                               max_iter = 16L,
                                               stiffness = 1.0,
                                               distance_floor = 1e-8,
                                               edge_length_epsilon = 1e-8,
                                               initial_step = 1.0,
                                               step_shrink = 0.5,
                                               armijo_factor = 1e-4,
                                               grad_tol = 1e-8,
                                               min_step = 1e-8,
                                               recenter = TRUE,
                                               return_trace = FALSE) {
  coords <- grip.validate.coords(coords)
  if (is.null(prepared)) {
    prepared <- grip.prepare.landmark.geodesic.kk(
      edges = edges,
      n = if (is.null(n)) nrow(coords) else n,
      adj_list = adj_list,
      weight_list = weight_list,
      edge_weights = edge_weights,
      local_nbrs = local_nbrs,
      landmark_count = landmark_count
    )
  }
  prepared <- grip.validate.prepared.object(
    prepared = prepared,
    class_name = "grip_lgkk_prepared",
    prepare_fun_name = "grip.prepare.landmark.geodesic.kk",
    coords = coords
  )

  grip.validate.scalar(max_iter, "max_iter", lower = 0)
  grip.validate.scalar(stiffness, "stiffness", lower = 0, open.lower = TRUE)
  grip.validate.scalar(distance_floor, "distance_floor", lower = 0, open.lower = TRUE)
  grip.validate.scalar(edge_length_epsilon, "edge_length_epsilon", lower = 0)
  grip.validate.scalar(initial_step, "initial_step", lower = 0, open.lower = TRUE)
  grip.validate.scalar(step_shrink, "step_shrink", lower = 0, upper = 1, open.lower = TRUE, open.upper = TRUE)
  grip.validate.scalar(armijo_factor, "armijo_factor", lower = 0)
  grip.validate.scalar(grad_tol, "grad_tol", lower = 0)
  grip.validate.scalar(min_step, "min_step", lower = 0, open.lower = TRUE)
  max_iter <- as.integer(round(max_iter))
  if (is.na(max_iter) || max_iter < 0L) {
    stop("max_iter must be a non-negative integer")
  }
  if (!is.logical(recenter) || length(recenter) != 1L || is.na(recenter)) {
    stop("recenter must be TRUE or FALSE")
  }
  if (!is.logical(return_trace) || length(return_trace) != 1L || is.na(return_trace)) {
    stop("return_trace must be TRUE or FALSE")
  }

  if (nrow(coords) <= 1L || length(prepared$pair_graph_distance) == 0L || max_iter == 0L) {
    score <- grip.score.landmark.geodesic.kk(
      coords = coords,
      prepared = prepared,
      stiffness = stiffness,
      distance_floor = distance_floor,
      edge_length_epsilon = edge_length_epsilon
    )
    return(list(
      coords = coords,
      trace = data.frame(),
      frames = list(coords),
      prepared = prepared,
      score = score
    ))
  }

  current <- coords
  initial.path.lengths <- grip.geodesic.kk.path.lengths(
    current,
    prepared,
    edge_length_epsilon = edge_length_epsilon
  )
  scale.L0 <- grip.geodesic.kk.fit.scale(
    path.lengths = initial.path.lengths,
    graph.distances = prepared$pair_graph_distance,
    stiffness = stiffness,
    distance_floor = distance_floor
  )
  if (!is.finite(scale.L0)) {
    stop("failed to fit an initial LGKK scale")
  }

  trace.rows <- vector("list", max_iter + 1L)
  accepted.frames <- list(current)
  state <- grip.geodesic.kk.energy.gradient(
    current,
    prepared = prepared,
    scale.L0 = scale.L0,
    stiffness = stiffness,
    distance_floor = distance_floor,
    edge_length_epsilon = edge_length_epsilon
  )
  trace.rows[[1L]] <- data.frame(
    iteration = 0L,
    energy = state$energy,
    gradient_norm = state$gradient_norm,
    step = NA_real_,
    accepted = TRUE,
    scale.L0 = scale.L0,
    stringsAsFactors = FALSE
  )
  used <- 1L

  for (iter in seq_len(max_iter)) {
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
      proposal.state <- grip.geodesic.kk.energy.gradient(
        proposal,
        prepared = prepared,
        scale.L0 = scale.L0,
        stiffness = stiffness,
        distance_floor = distance_floor,
        edge_length_epsilon = edge_length_epsilon
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
      gradient_norm = if (accepted) candidate.state$gradient_norm else state$gradient_norm,
      step = if (accepted) step else NA_real_,
      accepted = accepted,
      scale.L0 = scale.L0,
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
  score <- grip.score.landmark.geodesic.kk(
    coords = current,
    prepared = prepared,
    stiffness = stiffness,
    distance_floor = distance_floor,
    edge_length_epsilon = edge_length_epsilon
  )
  if (!isTRUE(return_trace)) {
    trace.df <- trace.df[, c("iteration", "energy", "gradient_norm", "step", "accepted"), drop = FALSE]
    accepted.frames <- list(current)
  }

  list(
    coords = current,
    trace = trace.df,
    frames = accepted.frames,
    prepared = prepared,
    score = score
  )
}

#' Optimize a layout under the full geodesic KK energy
#'
#' \code{grip.optimize.geodesic.kk()} applies a deterministic gradient-descent
#' polish under the full all-pairs geodesic Kamada--Kawai objective.
#'
#' The optimizer supports three scale policies. With
#' \code{scale_mode = "fixed_initial"} the target scale is fit once from the
#' starting layout and then held fixed during optimization, matching the current
#' LGKK prototype behavior. With \code{scale_mode = "profiled"} the scale is
#' re-fit analytically at each evaluation. With \code{scale_mode = "user"}, a
#' fixed user-supplied \code{scale.L0} is used throughout.
#'
#' @param coords Numeric coordinate matrix with 2 or 3 columns.
#' @param prepared Optional object returned by
#'   \code{\link{grip.prepare.geodesic.kk}()}.
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#' @param n Number of vertices.
#' @param adj_list Adjacency list (1-based) for an undirected graph.
#' @param weight_list Optional parallel list of positive edge weights.
#' @param edge_weights Optional positive edge-weight vector parallel to
#'   \code{edges}.
#' @param max_iter Maximum number of gradient-descent iterations.
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
#' @param recenter If \code{TRUE}, recenter the layout to zero mean after each
#'   accepted step.
#' @param return_trace If \code{TRUE}, include per-iteration diagnostics and the
#'   accepted intermediate coordinate frames.
#' @param scale_mode One of \code{"fixed_initial"}, \code{"profiled"}, or
#'   \code{"user"}.
#' @param scale.L0 Optional user-supplied fixed scale, required when
#'   \code{scale_mode = "user"}.
#'
#' @return A list with \code{coords}, \code{trace}, \code{frames},
#'   \code{prepared}, and \code{score}.
#' @export
grip.optimize.geodesic.kk <- function(coords,
                                      prepared = NULL,
                                      edges = NULL,
                                      n = NULL,
                                      adj_list = NULL,
                                      weight_list = NULL,
                                      edge_weights = NULL,
                                      max_iter = 16L,
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
                                      scale_mode = c("fixed_initial", "profiled", "user"),
                                      scale.L0 = NULL) {
  coords <- grip.validate.coords(coords)
  scale_mode <- match.arg(scale_mode)
  if (is.null(prepared)) {
    prepared <- grip.prepare.geodesic.kk(
      edges = edges,
      n = if (is.null(n)) nrow(coords) else n,
      adj_list = adj_list,
      weight_list = weight_list,
      edge_weights = edge_weights
    )
  }
  prepared <- grip.validate.prepared.object(
    prepared = prepared,
    class_name = "grip_gkk_prepared",
    prepare_fun_name = "grip.prepare.geodesic.kk",
    coords = coords
  )

  grip.validate.scalar(max_iter, "max_iter", lower = 0)
  grip.validate.scalar(stiffness, "stiffness", lower = 0, open.lower = TRUE)
  grip.validate.scalar(distance_floor, "distance_floor", lower = 0, open.lower = TRUE)
  grip.validate.scalar(edge_length_epsilon, "edge_length_epsilon", lower = 0)
  grip.validate.scalar(initial_step, "initial_step", lower = 0, open.lower = TRUE)
  grip.validate.scalar(step_shrink, "step_shrink", lower = 0, upper = 1, open.lower = TRUE, open.upper = TRUE)
  grip.validate.scalar(armijo_factor, "armijo_factor", lower = 0)
  grip.validate.scalar(grad_tol, "grad_tol", lower = 0)
  grip.validate.scalar(min_step, "min_step", lower = 0, open.lower = TRUE)
  max_iter <- as.integer(round(max_iter))
  if (is.na(max_iter) || max_iter < 0L) {
    stop("max_iter must be a non-negative integer")
  }
  if (!is.logical(recenter) || length(recenter) != 1L || is.na(recenter)) {
    stop("recenter must be TRUE or FALSE")
  }
  if (!is.logical(return_trace) || length(return_trace) != 1L || is.na(return_trace)) {
    stop("return_trace must be TRUE or FALSE")
  }
  if (identical(scale_mode, "user")) {
    grip.validate.scalar(scale.L0, "scale.L0", lower = 0, open.lower = TRUE)
  }

  current <- coords
  fixed.scale.L0 <- NULL
  if (length(prepared$pair_graph_distance) > 0L) {
    if (identical(scale_mode, "fixed_initial")) {
      fixed.scale.L0 <- grip.geodesic.kk.fit.scale(
        path.lengths = grip.geodesic.kk.path.lengths(
          coords = current,
          prepared = prepared,
          edge_length_epsilon = edge_length_epsilon
        ),
        graph.distances = prepared$pair_graph_distance,
        stiffness = stiffness,
        distance_floor = distance_floor
      )
      if (!is.finite(fixed.scale.L0)) {
        stop("failed to fit an initial GKK scale")
      }
    } else if (identical(scale_mode, "user")) {
      fixed.scale.L0 <- as.double(scale.L0)
    }
  } else if (identical(scale_mode, "user")) {
    fixed.scale.L0 <- as.double(scale.L0)
  }

  score.mode <- if (identical(scale_mode, "profiled")) "profiled" else "user"
  score.scale <- if (identical(score.mode, "user")) fixed.scale.L0 else NULL

  if (nrow(coords) <= 1L || length(prepared$pair_graph_distance) == 0L || max_iter == 0L) {
    score <- grip.score.geodesic.kk(
      coords = coords,
      prepared = prepared,
      stiffness = stiffness,
      distance_floor = distance_floor,
      edge_length_epsilon = edge_length_epsilon,
      scale_mode = score.mode,
      scale.L0 = score.scale
    )
    return(list(
      coords = coords,
      trace = data.frame(),
      frames = list(coords),
      prepared = prepared,
      score = score
    ))
  }

  trace.rows <- vector("list", max_iter + 1L)
  accepted.frames <- list(current)
  state <- if (identical(scale_mode, "profiled")) {
    grip.geodesic.kk.evaluate.state(
      coords = current,
      prepared = prepared,
      stiffness = stiffness,
      distance_floor = distance_floor,
      edge_length_epsilon = edge_length_epsilon,
      scale_mode = "profiled"
    )
  } else {
    grip.geodesic.kk.evaluate.state(
      coords = current,
      prepared = prepared,
      stiffness = stiffness,
      distance_floor = distance_floor,
      edge_length_epsilon = edge_length_epsilon,
      scale_mode = "fixed",
      scale.L0 = fixed.scale.L0
    )
  }
  trace.rows[[1L]] <- data.frame(
    iteration = 0L,
    energy = state$energy,
    gradient_norm = state$gradient_norm,
    step = NA_real_,
    accepted = TRUE,
    scale.L0 = state$scale.L0,
    stringsAsFactors = FALSE
  )
  used <- 1L

  for (iter in seq_len(max_iter)) {
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
      proposal.state <- if (identical(scale_mode, "profiled")) {
        grip.geodesic.kk.evaluate.state(
          coords = proposal,
          prepared = prepared,
          stiffness = stiffness,
          distance_floor = distance_floor,
          edge_length_epsilon = edge_length_epsilon,
          scale_mode = "profiled"
        )
      } else {
        grip.geodesic.kk.evaluate.state(
          coords = proposal,
          prepared = prepared,
          stiffness = stiffness,
          distance_floor = distance_floor,
          edge_length_epsilon = edge_length_epsilon,
          scale_mode = "fixed",
          scale.L0 = fixed.scale.L0
        )
      }
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
      gradient_norm = if (accepted) candidate.state$gradient_norm else state$gradient_norm,
      step = if (accepted) step else NA_real_,
      accepted = accepted,
      scale.L0 = if (accepted) candidate.state$scale.L0 else state$scale.L0,
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
  score <- grip.score.geodesic.kk(
    coords = current,
    prepared = prepared,
    stiffness = stiffness,
    distance_floor = distance_floor,
    edge_length_epsilon = edge_length_epsilon,
    scale_mode = score.mode,
    scale.L0 = score.scale
  )
  if (!isTRUE(return_trace)) {
    trace.df <- trace.df[, c("iteration", "energy", "gradient_norm", "step", "accepted"), drop = FALSE]
    accepted.frames <- list(current)
  }

  list(
    coords = current,
    trace = trace.df,
    frames = accepted.frames,
    prepared = prepared,
    score = score
  )
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

grip.induced.layout.inputs <- function(adj.list, weight.list = NULL, active) {
  active <- as.logical(active)
  if (length(active) != length(adj.list)) {
    stop("active must be parallel to adj.list")
  }
  vertex.ids <- which(active)
  map <- integer(length(active))
  map[vertex.ids] <- seq_along(vertex.ids)

  sub.adj <- vector("list", length(vertex.ids))
  sub.weights <- if (is.null(weight.list)) NULL else vector("list", length(vertex.ids))

  for (i in seq_along(vertex.ids)) {
    old.id <- vertex.ids[[i]]
    nb <- adj.list[[old.id]]
    keep <- active[nb]
    sub.adj[[i]] <- if (any(keep)) {
      as.integer(map[nb[keep]])
    } else {
      integer(0L)
    }
    if (!is.null(weight.list)) {
      ww <- weight.list[[old.id]]
      sub.weights[[i]] <- if (any(keep)) as.double(ww[keep]) else numeric(0L)
    }
  }

  list(
    vertex.ids = vertex.ids,
    adj.list = sub.adj,
    weight.list = sub.weights,
    edges = grip.edges.from.adj.list(sub.adj)
  )
}

grip.validate.trace.diagnostics <- function(diagnostics = c("none", "light", "full"),
                                            target.coords = NULL,
                                            n = NULL,
                                            dim = NULL,
                                            sample.size.nonedge = 1000L,
                                            sample.size.stress = 500L,
                                            nonedge.seed = 1L,
                                            stress.seed = 1L) {
  diagnostics <- match.arg(diagnostics)

  if (!is.null(target.coords)) {
    target.coords <- as.matrix(target.coords)
    if (!is.numeric(target.coords) || any(!is.finite(target.coords))) {
      stop("target_coords must be a finite numeric matrix when supplied")
    }
    if (!is.null(n) && nrow(target.coords) != n) {
      stop("nrow(target_coords) must match the graph size")
    }
    if (!is.null(dim) && ncol(target.coords) != dim) {
      stop("ncol(target_coords) must match dim")
    }
  }

  sample.size.nonedge <- as.integer(sample.size.nonedge)
  if (is.na(sample.size.nonedge) || sample.size.nonedge <= 0L) {
    stop("diagnostic_sample_size_nonedge must be a positive integer")
  }

  sample.size.stress <- as.integer(sample.size.stress)
  if (is.na(sample.size.stress) || sample.size.stress <= 0L) {
    stop("diagnostic_sample_size_stress must be a positive integer")
  }

  nonedge.seed <- as.integer(nonedge.seed)
  stress.seed <- as.integer(stress.seed)
  if (is.na(nonedge.seed) || is.na(stress.seed)) {
    stop("diagnostic seeds must be finite integers")
  }

  list(
    diagnostics = diagnostics,
    target.coords = target.coords,
    sample.size.nonedge = sample.size.nonedge,
    sample.size.stress = sample.size.stress,
    nonedge.seed = nonedge.seed,
    stress.seed = stress.seed
  )
}

grip.trace.compute.diagnostics <- function(frames,
                                           meta,
                                           adj.list,
                                           weight.list = NULL,
                                           diagnostics = c("none", "light", "full"),
                                           target.coords = NULL,
                                           sample.size.nonedge = 1000L,
                                           sample.size.stress = 500L,
                                           nonedge.seed = 1L,
                                           stress.seed = 1L) {
  validated <- grip.validate.trace.diagnostics(
    diagnostics = diagnostics,
    target.coords = target.coords,
    n = length(adj.list),
    dim = ncol(frames[[1L]]),
    sample.size.nonedge = sample.size.nonedge,
    sample.size.stress = sample.size.stress,
    nonedge.seed = nonedge.seed,
    stress.seed = stress.seed
  )
  diagnostics <- validated$diagnostics
  if (identical(diagnostics, "none")) {
    return(NULL)
  }

  rows <- vector("list", length(frames))
  for (i in seq_along(frames)) {
    frame.coords <- frames[[i]]
    active <- stats::complete.cases(frame.coords)
    active.count <- sum(active)

    row <- meta[i, , drop = FALSE]
    row$active.edges <- NA_integer_
    row$edge.length.cv <- NA_real_
    row$median.edge.length <- NA_real_
    row$sampled.nonedge.sep.ratio <- NA_real_
    row$sampled.stress <- NA_real_
    row$procrustes.rmse <- NA_real_

    if (active.count >= 2L) {
      induced <- grip.induced.layout.inputs(adj.list, weight.list, active)
      coords.active <- frame.coords[induced$vertex.ids, , drop = FALSE]
      edge.stats <- grip.edge.length.stats(coords.active, induced$edges)

      row$active.edges <- nrow(induced$edges)
      row$edge.length.cv <- edge.stats$cv
      row$median.edge.length <- edge.stats$median
      row$sampled.nonedge.sep.ratio <- grip.sampled.nonedge.sep.ratio(
        coords.active,
        induced$edges,
        sample.size = validated$sample.size.nonedge,
        rng.seed = validated$nonedge.seed + i - 1L
      )

      if (identical(diagnostics, "full")) {
        row$sampled.stress <- grip.sampled.stress(
          coords.active,
          adj.list = induced$adj.list,
          weight.list = induced$weight.list,
          sample.size = validated$sample.size.stress,
          rng.seed = validated$stress.seed + i - 1L
        )
      }

      if (!is.null(validated$target.coords)) {
        target.active <- validated$target.coords[induced$vertex.ids, , drop = FALSE]
        row$procrustes.rmse <- grip.align.to.target.nd(coords.active, target.active)$rmse
      }
    }

    rows[[i]] <- row
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
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

grip.compare.allowed.args <- function() {
  c(
    "placement", "preset", "rounds", "final_rounds", "num_init",
    "num_nbrs", "r", "s", "repulsion_factor", "tinit_factor"
  )
}

grip.expand.compare.search <- function(search) {
  if (is.null(search)) {
    return(list())
  }
  if (!is.list(search) || length(search) == 0L || is.null(names(search))) {
    stop("search must be a named list")
  }

  search <- search[!vapply(search, is.null, logical(1L))]
  if (length(search) == 0L) {
    stop("search must contain at least one non-NULL field")
  }

  candidate.prefix <- if ("candidate.prefix" %in% names(search)) {
    prefix <- search[["candidate.prefix"]]
    if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix) || !nzchar(prefix)) {
      stop("search$candidate.prefix must be a single non-empty character value")
    }
    prefix
  } else {
    "search"
  }

  include.base <- if ("include.base" %in% names(search)) {
    include <- search[["include.base"]]
    if (!is.logical(include) || length(include) != 1L || is.na(include)) {
      stop("search$include.base must be TRUE or FALSE")
    }
    include
  } else {
    FALSE
  }

  tuning.names <- setdiff(names(search), c("candidate.prefix", "include.base"))
  bad <- setdiff(tuning.names, grip.compare.allowed.args())
  if (length(bad) > 0L) {
    stop("search contains unsupported layout arguments: ", paste(bad, collapse = ", "))
  }
  if (length(tuning.names) == 0L) {
    stop("search must contain at least one layout argument to vary or fix")
  }

  values <- search[tuning.names]
  values <- lapply(values, function(x) {
    if (length(x) == 0L) {
      stop("search parameter vectors must be non-empty")
    }
    x
  })

  grid <- do.call(expand.grid, c(values, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE))
  if (include.base) {
    base.row <- as.data.frame(lapply(values, function(x) x[[1L]]), stringsAsFactors = FALSE)
    names(base.row) <- names(values)
    grid <- unique(rbind(base.row, grid))
  }

  varying.names <- names(values)[vapply(values, function(x) length(unique(as.character(x))) > 1L, logical(1L))]
  out <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    args <- as.list(grid[i, , drop = FALSE])
    args <- args[!vapply(args, function(x) length(x) == 1L && is.character(x) && is.na(x), logical(1L))]
    out[[i]] <- args
    if (length(varying.names) == 0L) {
      names(out)[[i]] <- candidate.prefix
    } else {
      label.bits <- vapply(varying.names, function(name) {
        value <- args[[name]]
        value.txt <- gsub("[^[:alnum:]]+", "", as.character(value))
        paste0(name, ".", value.txt)
      }, character(1L))
      names(out)[[i]] <- paste(c(candidate.prefix, label.bits), collapse = ".")
    }
  }

  if (anyDuplicated(names(out))) {
    seq.names <- sprintf("%s.%03d", candidate.prefix, seq_along(out))
    names(out) <- seq.names
  }

  out
}

grip.resolve.compare.candidate <- function(candidate, dim = 2L) {
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
    dim = dim,
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
      placement = df$placement[[1L]],
      rounds = df$rounds[[1L]],
      final.rounds = df$final.rounds[[1L]],
      num.init = df$num.init[[1L]],
      num.nbrs = df$num.nbrs[[1L]],
      r = df$r[[1L]],
      s = df$s[[1L]],
      repulsion.factor = df$repulsion.factor[[1L]],
      tinit.factor = df$tinit.factor[[1L]],
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
#' canonical embedding. It is the low-level scoring helper behind
#' \code{\link{grip.compare.layouts}()} and is most useful when you already have
#' one realized layout in hand, for example from a cached run or another graph
#' drawing tool. For real-world graphs, quality is judged by graph-distance
#' faithfulness, edge-length consistency, separation of non-neighbors, and
#' optionally edge crossings or cluster separation.
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
#' or parameter lists, optionally expands a local parameter search, scores each
#' run with \code{\link{grip.score.layout}()}, and summarizes both quality
#' metrics and seed-to-seed stability. This is the main real-data workflow for
#' graphs where no canonical embedding is known.
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
#' @param search Optional named list describing a grid search over layout
#'   settings. Any of \code{preset}, \code{placement}, \code{rounds},
#'   \code{final_rounds}, \code{num_init}, \code{num_nbrs}, \code{r},
#'   \code{s}, \code{repulsion_factor}, and \code{tinit_factor} may be supplied
#'   as vectors. All combinations are expanded into candidates. Special fields
#'   \code{candidate.prefix} and \code{include.base} control candidate naming
#'   and whether the all-first-values setting is guaranteed to appear.
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
#' edges <- edges.path(5)
#' cmp <- grip.compare.layouts(
#'   edges = edges,
#'   n = 5,
#'   dim = 2,
#'   candidates = c("default", "tree"),
#'   seeds = 1L
#' )
#' cmp$summary[, c("candidate", "rounds", "final.rounds")]
#'
#' search.cmp <- grip.compare.layouts(
#'   edges = edges,
#'   n = 5,
#'   dim = 2,
#'   search = list(
#'     candidate.prefix = "path.search",
#'     rounds = c(4L, 6L),
#'     final_rounds = c(4L, 6L)
#'   ),
#'   seeds = 1L
#' )
#' search.cmp$summary[, c("candidate", "rounds", "final.rounds")]
#' @export
grip.compare.layouts <- function(edges = NULL,
                                 n = NULL,
                                 adj_list = NULL,
                                 weight_list = NULL,
                                 edge_weights = NULL,
                                 dim = 2,
                                 candidates = c("default"),
                                 search = NULL,
                                 clusters = NULL,
                                 seeds = 1:3,
                                 sample.size.stress = 2000L,
                                 sample.size.nonedge = 5000L,
                                 edge.crossings = c("auto", "always", "never"),
                                 edge.crossings.max.edges = 1000L,
                                 score.weights = grip.default.compare.score.weights(),
                                 return.layouts = FALSE,
                                 disconnected = c("components", "error")) {
  candidates.missing <- missing(candidates)
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
  explicit.candidates <- if (is.null(search) && candidates.missing) {
    grip.normalize.compare.candidates(candidates)
  } else if (!candidates.missing && !is.null(candidates)) {
    grip.normalize.compare.candidates(candidates)
  } else {
    list()
  }
  search.candidates <- grip.expand.compare.search(search)
  candidate.list <- c(explicit.candidates, search.candidates)
  if (length(candidate.list) == 0L) {
    stop("no layout candidates were generated")
  }
  if (anyDuplicated(names(candidate.list))) {
    stop("candidate names must be unique after combining candidates and search")
  }
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
    resolved <- grip.resolve.compare.candidate(candidate.list[[candidate.name]], dim = dim)
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
