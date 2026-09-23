grip.null.coalesce <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

grip.validate.star.structure <- function(star, n = NULL) {
  if (!inherits(star, "grip_riemannian_star")) {
    stop("star must be an object returned by graph.riemannian.star.structure()")
  }
  required <- c(
    "center", "neighbor_1", "neighbor_2",
    "w_1", "w_2", "cos_angle", "reliability", "angle_weight"
  )
  missing <- setdiff(required, names(star$pairs))
  if (length(missing) > 0L) {
    stop("star$pairs is missing required column(s): ", paste(missing, collapse = ", "))
  }
  if (!is.null(n) && !identical(as.integer(star$n), as.integer(n))) {
    stop("star structure vertex count does not match the graph")
  }
  star
}

grip.graph.adjacency.from.edges <- function(edges, edge.weights, n) {
  adj <- vector("list", n)
  weights <- vector("list", n)
  for (i in seq_len(n)) {
    adj[[i]] <- integer(0L)
    weights[[i]] <- numeric(0L)
  }
  if (nrow(edges) > 0L) {
    for (e in seq_len(nrow(edges))) {
      u <- edges[e, 1L]
      v <- edges[e, 2L]
      w <- edge.weights[[e]]
      adj[[u]] <- c(adj[[u]], v)
      weights[[u]] <- c(weights[[u]], w)
      adj[[v]] <- c(adj[[v]], u)
      weights[[v]] <- c(weights[[v]], w)
    }
  }
  list(adj_list = adj, weight_list = weights)
}

#' Build a local Riemannian star structure for Gram-gKK layouts
#'
#' `graph.riemannian.star.structure()` builds the local star-pair table used by
#' kernel Gram-gKK. For each center vertex `u`, it considers unordered neighbor
#' pairs `(v, v')` in the graph star and stores the target edge lengths, target
#' cosine angle, and kernel weight
#' \deqn{
#'   r_u(v,v')\left(\frac{1-\cos\alpha_u(v,v')}{2}\right)^q.
#' }
#' In the first implementation, target angles are estimated from the ambient
#' coordinates `X`. The optimizer is intentionally agnostic to where those
#' cosines came from, so later graph constructors can attach a different
#' Riemannian source without changing the layout backend.
#'
#' @param graph Optional graph/prepared object containing `adj.list` and
#'   `weight.list`, or `edges`/`edge_targets` when it is a prepared GMDS graph.
#' @param X Numeric matrix used to estimate local target angles.
#' @param prepared Optional prepared GMDS object. Used when `graph` is omitted.
#' @param edges,n,edge.weights Optional edge representation used when `graph`
#'   and `prepared` are omitted.
#' @param adj.list,weight.list Optional adjacency-list representation.
#' @param angle.power Non-negative exponent `q` in the antipodal kernel.
#' @param reliability Reliability weighting rule. `"length.balance"` multiplies
#'   by `min(w_1, w_2) / max(w_1, w_2)`; `"none"` uses one.
#' @param min.angle.weight Pairs with final angle weight at or below this value
#'   are omitted from the returned table.
#' @param star.quantile Optional quantile in `[0, 1)`. When positive, retain
#'   only star pairs whose angle weight is at least this empirical quantile
#'   after applying `min.angle.weight`.
#'
#' @return A list of class `"grip_riemannian_star"` with the star-pair table and
#'   construction metadata.
#' @export
#' @md
graph.riemannian.star.structure <- function(graph = NULL,
                                            X,
                                            prepared = NULL,
                                            edges = NULL,
                                            n = NULL,
                                            adj.list = NULL,
                                            weight.list = NULL,
                                            edge.weights = NULL,
                                            angle.power = 4,
                                            reliability = c("length.balance", "none"),
                                            min.angle.weight = 0,
                                            star.quantile = 0) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights, prepared)
  X <- grip.validate.coords(X)
  reliability <- match.arg(reliability)
  grip.validate.scalar(angle.power, "angle.power", lower = 0)
  grip.validate.scalar(min.angle.weight, "min.angle.weight", lower = 0)
  grip.validate.scalar(star.quantile, "star.quantile", lower = 0, upper = 1, open.upper = TRUE)

  if (!is.null(graph)) {
    if (any(!vapply(list(prepared, edges, adj.list, weight.list, edge.weights), is.null, logical(1)))) {
      stop("graph cannot be combined with prepared or raw graph inputs; supply one graph representation",
           call. = FALSE)
    }
    if (!is.null(n) && !is.null(graph$n) &&
        !identical(grip.validate.vertex.count(n), grip.validate.vertex.count(graph$n))) {
      stop("n must match the graph size stored in graph", call. = FALSE)
    }
    if (inherits(graph, "grip_geodesic_kk_prepared")) {
      prepared <- graph
    } else {
      adj.list <- grip.null.coalesce(graph$adj_list, graph$adj.list)
      weight.list <- grip.null.coalesce(graph$weight_list, graph$weight.list)
      edges <- grip.null.coalesce(graph$edges, edges)
      edge.weights <- grip.null.coalesce(graph$edge_targets, grip.null.coalesce(graph$edge_weights, edge.weights))
      n <- grip.null.coalesce(graph$n, n)
      # Bundles may store both representations, but they must describe the
      # same weighted graph. Preserve the adjacency ordering after checking.
      if (!is.null(edges) && !is.null(adj.list)) {
        from.edges <- prepare.edge.kk(edges = edges, n = n, edge.weights = edge.weights)
        from.adj <- prepare.edge.kk(adj.list = adj.list, n = n, weight.list = weight.list)
        if (!identical(from.edges$edges, from.adj$edges) ||
            !identical(from.edges$edge_targets, from.adj$edge_targets)) {
          stop("graph contains contradictory edge and adjacency representations", call. = FALSE)
        }
        edges <- edge.weights <- NULL
      }
    }
  }
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights, prepared)
  if (!is.null(prepared)) {
    prepared <- grip.validate.geodesic.mds.prepared(prepared)
    n <- prepared$n
    adj.list <- prepared$adj_list
    weight.list <- prepared$weight_list
    if (is.null(adj.list) || is.null(weight.list)) {
      aw <- grip.graph.adjacency.from.edges(prepared$edges, prepared$edge_targets, prepared$n)
      adj.list <- aw$adj_list
      weight.list <- aw$weight_list
    }
  } else if (is.null(adj.list) || is.null(weight.list)) {
    if (is.null(edges) || is.null(edge.weights)) {
      stop("provide graph/prepared, adj_list/weight_list, or edges/edge_weights")
    }
    if (is.null(n)) {
      n <- grip.resolve.graph.n(n, edges, adj.list)
    }
    aw <- grip.graph.adjacency.from.edges(matrix(as.integer(edges), ncol = 2L), edge.weights, as.integer(n))
    adj.list <- aw$adj_list
    weight.list <- aw$weight_list
  }

  n <- grip.validate.vertex.count(grip.null.coalesce(n, length(adj.list)))
  if (nrow(X) != n) {
    stop("nrow(X) must match the graph vertex count")
  }
  if (!is.list(adj.list) || !is.list(weight.list) || length(adj.list) != n || length(weight.list) != n) {
    stop("adj.list and weight.list must be lists parallel to the graph vertices")
  }

  rows <- vector("list", n)
  row.idx <- 0L
  for (u in seq_len(n)) {
    nbrs <- as.integer(adj.list[[u]])
    weights <- as.double(weight.list[[u]])
    if (length(nbrs) != length(weights)) {
      stop("weight.list entries must be parallel to adj.list entries")
    }
    if (length(nbrs) < 2L) {
      next
    }
    pair.idx <- utils::combn(seq_along(nbrs), 2L)
    center <- X[u, , drop = FALSE]
    for (p in seq_len(ncol(pair.idx))) {
      i1 <- pair.idx[1L, p]
      i2 <- pair.idx[2L, p]
      v1 <- nbrs[[i1]]
      v2 <- nbrs[[i2]]
      if (v1 < 1L || v1 > n || v2 < 1L || v2 > n || v1 == u || v2 == u || v1 == v2) {
        stop("adj_list contains invalid vertex indices")
      }
      a <- X[v1, , drop = FALSE] - center
      b <- X[v2, , drop = FALSE] - center
      na <- sqrt(sum(a^2))
      nb <- sqrt(sum(b^2))
      if (!is.finite(na) || !is.finite(nb) || na <= 0 || nb <= 0) {
        next
      }
      cos.angle <- as.double(sum(a * b) / (na * nb))
      cos.angle <- max(-1, min(1, cos.angle))
      w1 <- weights[[i1]]
      w2 <- weights[[i2]]
      if (!is.finite(w1) || !is.finite(w2) || w1 <= 0 || w2 <= 0) {
        stop("weight.list must contain finite positive edge lengths")
      }
      rel <- if (identical(reliability, "length.balance")) {
        min(w1, w2) / max(w1, w2)
      } else {
        1
      }
      angle.weight <- rel * ((1 - cos.angle) / 2)^angle.power
      if (!is.finite(angle.weight) || angle.weight <= min.angle.weight) {
        next
      }
      row.idx <- row.idx + 1L
      rows[[row.idx]] <- data.frame(
        center = u,
        neighbor_1 = v1,
        neighbor_2 = v2,
        w_1 = w1,
        w_2 = w2,
        cos_angle = cos.angle,
        reliability = rel,
        angle_weight = angle.weight,
        stringsAsFactors = FALSE
      )
    }
  }
  pairs <- if (row.idx == 0L) {
    data.frame(
      center = integer(0L),
      neighbor_1 = integer(0L),
      neighbor_2 = integer(0L),
      w_1 = numeric(0L),
      w_2 = numeric(0L),
      cos_angle = numeric(0L),
      reliability = numeric(0L),
      angle_weight = numeric(0L)
    )
  } else {
    do.call(rbind, rows[seq_len(row.idx)])
  }
  if (nrow(pairs) > 0L && star.quantile > 0) {
    cutoff <- as.double(stats::quantile(pairs$angle_weight, probs = star.quantile, names = FALSE))
    pairs <- pairs[pairs$angle_weight >= cutoff, , drop = FALSE]
  }
  out <- list(
    n = n,
    pairs = pairs,
    angle.power = as.double(angle.power),
    reliability = reliability,
    min.angle.weight = as.double(min.angle.weight),
    star.quantile = as.double(star.quantile),
    diagnostics = data.frame(
      n.vertices = n,
      n.star.pairs = nrow(pairs),
      angle.weight.mean = if (nrow(pairs)) mean(pairs$angle_weight) else NA_real_,
      angle.weight.max = if (nrow(pairs)) max(pairs$angle_weight) else NA_real_,
      cos.angle.min = if (nrow(pairs)) min(pairs$cos_angle) else NA_real_,
      cos.angle.q10 = if (nrow(pairs)) as.double(stats::quantile(pairs$cos_angle, 0.1, names = FALSE)) else NA_real_,
      stringsAsFactors = FALSE
    )
  )
  class(out) <- c("grip_riemannian_star", "list")
  out
}

grip.kernel.gram.energy.gradient <- function(coords,
                                             edges,
                                             edge.weights,
                                             edge.stiffness,
                                             star,
                                             edge.scale = 1,
                                             lambda.edge = 1,
                                             lambda.gram = 1,
                                             edge.length.epsilon = 1e-8) {
  coords <- grip.validate.coords(coords)
  star <- grip.validate.star.structure(star, nrow(coords))
  grip.validate.scalar(edge.scale, "edge.scale", lower = 0, open.lower = TRUE)
  grip.validate.scalar(lambda.edge, "lambda.edge", lower = 0)
  grip.validate.scalar(lambda.gram, "lambda.gram", lower = 0)
  edge.state <- grip.edge.isometric.energy.gradient(
    coords = coords,
    edges = edges,
    edge.weights = edge.weights,
    stiffness = lambda.edge * edge.stiffness,
    scale = edge.scale,
    edge.length.epsilon = edge.length.epsilon
  )
  gradient <- edge.state$gradient
  gram.energy <- 0
  gram.resid2 <- 0
  gram.denom <- 0
  pairs <- star$pairs
  if (nrow(pairs) > 0L && lambda.gram > 0) {
    for (r in seq_len(nrow(pairs))) {
      u <- pairs$center[[r]]
      v <- pairs$neighbor_1[[r]]
      w <- pairs$neighbor_2[[r]]
      a <- coords[v, ] - coords[u, ]
      b <- coords[w, ] - coords[u, ]
      observed <- sum(a * b)
      target <- edge.scale^2 * pairs$w_1[[r]] * pairs$w_2[[r]] * pairs$cos_angle[[r]]
      residual <- observed - target
      k <- lambda.gram * pairs$angle_weight[[r]]
      gram.energy <- gram.energy + 0.5 * k * residual^2
      gram.resid2 <- gram.resid2 + k * residual^2
      gram.denom <- gram.denom + k * target^2
      coeff <- k * residual
      gradient[v, ] <- gradient[v, ] + coeff * b
      gradient[w, ] <- gradient[w, ] + coeff * a
      gradient[u, ] <- gradient[u, ] - coeff * (a + b)
    }
  }
  list(
    energy = edge.state$energy + gram.energy,
    edge_energy = edge.state$energy,
    gram_energy = gram.energy,
    gradient = gradient,
    gradient_norm = sqrt(sum(gradient^2)),
    edge_lengths = edge.state$edge_lengths,
    edge_residuals = edge.state$residuals,
    edge_rel_rmse = grip.gmds.residual.stats(
      observed = edge.state$edge_lengths,
      target = edge.scale * edge.weights,
      weights = edge.stiffness
    )$stress,
    gram_rel_rmse = if (is.finite(gram.denom) && gram.denom > 0) sqrt(gram.resid2 / gram.denom) else NA_real_,
    edge_scale = edge.scale
  )
}

grip.kernel.gram.score <- function(coords,
                                   star,
                                   edge.scale = 1,
                                   lambda.gram = 1,
                                   distance.floor = 1e-8) {
  coords <- grip.validate.coords(coords)
  star <- grip.validate.star.structure(star, nrow(coords))
  grip.validate.scalar(edge.scale, "edge.scale", lower = 0, open.lower = TRUE)
  grip.validate.scalar(lambda.gram, "lambda.gram", lower = 0)
  pairs <- star$pairs
  if (nrow(pairs) == 0L) {
    return(data.frame(
      n.star.pairs = 0L,
      gram.rel.rmse = NA_real_,
      gram.rmse = NA_real_,
      gram.mean.rel.abs.error = NA_real_,
      gram.signed.bias = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  observed <- numeric(nrow(pairs))
  target <- numeric(nrow(pairs))
  for (r in seq_len(nrow(pairs))) {
    u <- pairs$center[[r]]
    v <- pairs$neighbor_1[[r]]
    w <- pairs$neighbor_2[[r]]
    observed[[r]] <- sum((coords[v, ] - coords[u, ]) * (coords[w, ] - coords[u, ]))
    target[[r]] <- edge.scale^2 * pairs$w_1[[r]] * pairs$w_2[[r]] * pairs$cos_angle[[r]]
  }
  weights <- lambda.gram * pairs$angle_weight
  denom <- sum(weights * target^2)
  residual <- observed - target
  rel <- residual / pmax(abs(target), distance.floor)
  data.frame(
    n.star.pairs = nrow(pairs),
    gram.rel.rmse = if (is.finite(denom) && denom > 0) sqrt(sum(weights * residual^2) / denom) else NA_real_,
    gram.rmse = sqrt(sum(weights * residual^2) / sum(weights)),
    gram.mean.rel.abs.error = sum(weights * abs(rel)) / sum(weights),
    gram.signed.bias = sum(weights * rel) / sum(weights),
    stringsAsFactors = FALSE
  )
}

#' Optimize a kernel Gram-gKK layout
#'
#' `kernel.gram.gkk()` extends edge-only gKK with a local
#' Riemannian star penalty. The off-diagonal Gram term preserves target inner
#' products between pairs of incident edge directions:
#' \deqn{
#'   \frac{\lambda_{\mathrm{gram}}}{2}\sum_{u}\sum_{v<v'\in N(u)}
#'   \omega_u(v,v')\left(
#'   \langle z_v-z_u,z_{v'}-z_u\rangle
#'   - s^2 w_{uv}w_{uv'}\cos\alpha_u(v,v')\right)^2.
#' }
#' The star weights are built by [graph.riemannian.star.structure()], usually
#' with an antipodal kernel controlled by `angle.power`.
#'
#' @inheritParams edge.kk
#' @param init Starting layout used when `coords` is omitted. `"classical_mds"`
#'   uses classical scaling (the default), `"metric_mds"` uses stress MDS
#'   through the default native SGD backend, both from an all-pairs prepared
#'   object; `"random"`
#'   uses centered Gaussian coordinates.
#' @param return.trace If `TRUE`, keep per-iteration trace rows and coordinate
#'   frames.
#' @param seed Random seed used only for `init = "random"`.
#' @param X Optional ambient/source coordinates used to build `star` when `star`
#'   is omitted.
#' @param star Optional object from [graph.riemannian.star.structure()].
#' @param angle.power,reliability,min.angle.weight Passed to
#'   [graph.riemannian.star.structure()] when `star` is omitted.
#' @param star.quantile Optional quantile filter passed to
#'   [graph.riemannian.star.structure()] when `star` is omitted.
#' @param lambda.edge,lambda.gram Non-negative weights for the diagonal
#'   edge-length and off-diagonal Gram penalties.
#' @param stiffness.method,stiffness.transform,density.mix,bandwidth,density.n
#'   Parameters passed to [edge.length.density.stiffness()] to construct
#'   edge-length stiffnesses for the diagonal edge term.
#' @param distance.power,stiffness.floor,stiffness.ceiling Additional stiffness
#'   constructor parameters.
#'
#' @return A `"grip_gmds_layout"` object with method `"kernel_gram_gkk"`.
#' @export
#' @md
kernel.gram.gkk <- function(coords = NULL,
                                                 prepared = NULL,
                                                 edges = NULL,
                                                 n = NULL,
                                                 adj.list = NULL,
                                                 weight.list = NULL,
                                                 edge.weights = NULL,
                                                 X = NULL,
                                                 star = NULL,
                                                 dim = 2L,
                                                 init = c("classical_mds", "metric_mds", "random"),
                                                 angle.power = 4,
                                                 reliability = c("length.balance", "none"),
                                                 min.angle.weight = 0,
                                                 star.quantile = 0,
                                                 lambda.edge = 1,
                                                 lambda.gram = 1,
                                                 stiffness.method = c("density", "uniform", "distance_power"),
                                                 stiffness.transform = c("identity", "sqrt", "log"),
                                                 density.mix = 1,
                                                 bandwidth = NULL,
                                                 density.n = 512L,
                                                 distance.power = 0,
                                                 stiffness.floor = 0,
                                                 stiffness.ceiling = Inf,
                                                 scale.mode = c("profiled", "identity", "user"),
                                                 scale = NULL,
                                                 max.iter = 50L,
                                                 initial.step = 0.1,
                                                 step.shrink = 0.5,
                                                 armijo.factor = 1e-4,
                                                 grad.tol = 1e-8,
                                                 min.step = 1e-8,
                                                 edge.length.epsilon = 1e-8,
                                                 distance.floor = 1e-8,
                                                 recenter = TRUE,
                                                 return.trace = TRUE,
                                                 diagnostics = TRUE,
                                                 seed = 1L,
                                                 engine = c("cpp", "R")) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights, prepared)
  init <- match.arg(init)
  reliability <- match.arg(reliability)
  stiffness.method <- match.arg(stiffness.method)
  stiffness.transform <- match.arg(stiffness.transform)
  scale.mode <- match.arg(scale.mode)
  engine <- match.arg(engine)
  grip.validate.scalar(lambda.edge, "lambda.edge", lower = 0)
  grip.validate.scalar(lambda.gram, "lambda.gram", lower = 0)
  grip.validate.scalar(dim, "dim", lower = 2, upper = 3)
  dim <- as.integer(round(dim))
  grip.validate.scalar(max.iter, "max.iter", lower = 0)
  max.iter <- as.integer(round(max.iter))
  grip.validate.scalar(initial.step, "initial.step", lower = 0, open.lower = TRUE)
  grip.validate.scalar(step.shrink, "step.shrink", lower = 0, upper = 1, open.lower = TRUE, open.upper = TRUE)
  grip.validate.scalar(armijo.factor, "armijo.factor", lower = 0)
  grip.validate.scalar(grad.tol, "grad.tol", lower = 0)
  grip.validate.scalar(min.step, "min.step", lower = 0, open.lower = TRUE)
  grip.validate.scalar(edge.length.epsilon, "edge.length.epsilon", lower = 0)
  grip.validate.scalar(distance.floor, "distance.floor", lower = 0, open.lower = TRUE)
  if (!is.logical(recenter) || length(recenter) != 1L || is.na(recenter)) {
    stop("recenter must be TRUE or FALSE")
  }
  if (!is.logical(return.trace) || length(return.trace) != 1L || is.na(return.trace)) {
    stop("return.trace must be TRUE or FALSE")
  }
  if (!is.logical(diagnostics) || length(diagnostics) != 1L || is.na(diagnostics)) {
    stop("diagnostics must be TRUE or FALSE")
  }
  if (identical(scale.mode, "user")) {
    grip.validate.scalar(scale, "scale", lower = 0, open.lower = TRUE)
  }

  prepared <- grip.gmds.require.prepared(
    prepared = prepared,
    coords = coords,
    edges = edges,
    n = n,
    adj.list = adj.list,
    weight.list = weight.list,
    edge.weights = edge.weights
  )
  if (is.null(prepared$edges) || is.null(prepared$edge_targets)) {
    stop("prepared object must contain edges and edge_targets")
  }
  if (is.null(star)) {
    if (is.null(X)) {
      stop("X is required when star is omitted")
    }
    star <- graph.riemannian.star.structure(
      prepared = prepared,
      X = X,
      angle.power = angle.power,
      reliability = reliability,
      min.angle.weight = min.angle.weight,
      star.quantile = star.quantile
    )
  } else {
    star <- grip.validate.star.structure(star, prepared$n)
  }

  current <- grip.edge.isometric.initial.coords(
    prepared = prepared,
    coords = coords,
    init = init,
    dim = dim,
    seed = seed
  )
  if (isTRUE(recenter)) {
    current <- sweep(current, 2L, colMeans(current), "-", check.margin = FALSE)
  }
  edges <- prepared$edges
  edge.targets <- as.double(prepared$edge_targets)
  stiff <- edge.length.density.stiffness(
    edge.weights = edge.targets,
    method = stiffness.method,
    mix = density.mix,
    bandwidth = bandwidth,
    density.n = density.n,
    transform = stiffness.transform,
    distance.power = distance.power,
    stiffness.floor = stiffness.floor,
    stiffness.ceiling = stiffness.ceiling
  )
  edge.stiffness <- stiff$stiffness

  state.scale.mode <- switch(
    scale.mode,
    profiled = "profiled",
    identity = "identity",
    user = "user"
  )
  state.scale <- if (identical(scale.mode, "user")) as.double(scale) else NA_real_

  if (identical(engine, "cpp")) {
    fit <- grip_optimize_kernel_gram_gkk_layout_cpp(
      edges = matrix(as.integer(edges), ncol = 2L),
      edge_weights = edge.targets,
      edge_stiffness = edge.stiffness,
      star_center = as.integer(star$pairs$center),
      star_v1 = as.integer(star$pairs$neighbor_1),
      star_v2 = as.integer(star$pairs$neighbor_2),
      star_w1 = as.double(star$pairs$w_1),
      star_w2 = as.double(star$pairs$w_2),
      star_cos = as.double(star$pairs$cos_angle),
      star_weight = as.double(star$pairs$angle_weight),
      coords = current,
      max_iter = max.iter,
      scale_mode = state.scale.mode,
      scale = state.scale,
      lambda_edge = lambda.edge,
      lambda_gram = lambda.gram,
      edge_length_epsilon = edge.length.epsilon,
      initial_step = initial.step,
      step_shrink = step.shrink,
      armijo_factor = armijo.factor,
      grad_tol = grad.tol,
      min_step = min.step,
      distance_floor = distance.floor,
      recenter = recenter,
      return_trace = return.trace
    )
    current <- fit$coords
    trace.df <- fit$trace
    edge.scale <- if (nrow(trace.df)) trace.df$edge.scale[[nrow(trace.df)]] else 1
    gram.diag <- grip.kernel.gram.score(
      coords = current,
      star = star,
      edge.scale = edge.scale,
      lambda.gram = lambda.gram,
      distance.floor = distance.floor
    )
    diag <- if (isTRUE(diagnostics)) {
      cbind(
        score.gmds(
          coords = current,
          prepared = prepared,
          scale.mode = if (identical(scale.mode, "identity")) "identity" else "profiled",
          distance.floor = distance.floor,
          edge.length.epsilon = edge.length.epsilon
        ),
        gram.diag
      )
    } else {
      NULL
    }
    return(gmds.result(
      coords = current,
      method = "kernel_gram_gkk",
      prepared = prepared,
      trace = trace.df,
      diagnostics = diag,
      metadata = list(
        engine = "cpp_gradient_descent_armijo",
        initialization = if (is.null(coords)) init else "supplied",
        star = star,
        angle.power = star$angle.power,
        reliability = star$reliability,
        star.quantile = star$star.quantile,
        lambda.edge = lambda.edge,
        lambda.gram = lambda.gram,
        stiffness = stiff,
        frames = fit$frames,
        scale_mode = scale.mode
      )
    ))
  }

  trace.rows <- list()
  frames <- list(current)
  for (iter in 0:max.iter) {
    edge.lengths <- if (nrow(edges) == 0L) {
      numeric(0L)
    } else {
      diffs <- current[edges[, 1L], , drop = FALSE] - current[edges[, 2L], , drop = FALSE]
      sqrt(rowSums(diffs^2) + edge.length.epsilon^2)
    }
    edge.scale <- switch(
      scale.mode,
      profiled = grip.edge.isometric.fit.scale(edge.lengths, edge.targets, edge.stiffness, distance.floor),
      identity = 1.0,
      user = as.double(scale)
    )
    if (!is.finite(edge.scale)) {
      edge.scale <- 1.0
    }
    state <- grip.kernel.gram.energy.gradient(
      coords = current,
      edges = edges,
      edge.weights = edge.targets,
      edge.stiffness = edge.stiffness,
      star = star,
      edge.scale = edge.scale,
      lambda.edge = lambda.edge,
      lambda.gram = lambda.gram,
      edge.length.epsilon = edge.length.epsilon
    )
    if (iter == 0L) {
      trace.rows[[length(trace.rows) + 1L]] <- data.frame(
        iteration = iter,
        energy = state$energy,
        edge.energy = state$edge_energy,
        gram.energy = state$gram_energy,
        gradient_norm = state$gradient_norm,
        step = NA_real_,
        accepted = TRUE,
        edge.scale = state$edge_scale,
        edge.rel.rmse = state$edge_rel_rmse,
        gram.rel.rmse = state$gram_rel_rmse,
        stringsAsFactors = FALSE
      )
      next
    }
    if (!is.finite(state$gradient_norm) || state$gradient_norm <= grad.tol) {
      break
    }
    step <- as.double(initial.step)
    accepted <- FALSE
    candidate <- current
    candidate.state <- state
    while (is.finite(step) && step >= min.step) {
      proposal <- current - step * state$gradient
      if (isTRUE(recenter)) {
        proposal <- sweep(proposal, 2L, colMeans(proposal), "-", check.margin = FALSE)
      }
      edge.lengths.p <- if (nrow(edges) == 0L) {
        numeric(0L)
      } else {
        diffs <- proposal[edges[, 1L], , drop = FALSE] - proposal[edges[, 2L], , drop = FALSE]
        sqrt(rowSums(diffs^2) + edge.length.epsilon^2)
      }
      edge.scale.p <- switch(
        scale.mode,
        profiled = grip.edge.isometric.fit.scale(edge.lengths.p, edge.targets, edge.stiffness, distance.floor),
        identity = 1.0,
        user = as.double(scale)
      )
      if (!is.finite(edge.scale.p)) {
        edge.scale.p <- 1.0
      }
      proposal.state <- grip.kernel.gram.energy.gradient(
        coords = proposal,
        edges = edges,
        edge.weights = edge.targets,
        edge.stiffness = edge.stiffness,
        star = star,
        edge.scale = edge.scale.p,
        lambda.edge = lambda.edge,
        lambda.gram = lambda.gram,
        edge.length.epsilon = edge.length.epsilon
      )
      target.energy <- state$energy - armijo.factor * step * state$gradient_norm^2
      if (is.finite(proposal.state$energy) && proposal.state$energy <= target.energy) {
        candidate <- proposal
        candidate.state <- proposal.state
        accepted <- TRUE
        break
      }
      step <- step * step.shrink
    }
    trace.rows[[length(trace.rows) + 1L]] <- data.frame(
      iteration = iter,
      energy = if (accepted) candidate.state$energy else state$energy,
      edge.energy = if (accepted) candidate.state$edge_energy else state$edge_energy,
      gram.energy = if (accepted) candidate.state$gram_energy else state$gram_energy,
      gradient_norm = if (accepted) candidate.state$gradient_norm else state$gradient_norm,
      step = if (accepted) step else NA_real_,
      accepted = accepted,
      edge.scale = if (accepted) candidate.state$edge_scale else state$edge_scale,
      edge.rel.rmse = if (accepted) candidate.state$edge_rel_rmse else state$edge_rel_rmse,
      gram.rel.rmse = if (accepted) candidate.state$gram_rel_rmse else state$gram_rel_rmse,
      stringsAsFactors = FALSE
    )
    if (!accepted) {
      break
    }
    current <- candidate
    if (isTRUE(return.trace)) {
      frames[[length(frames) + 1L]] <- current
    }
  }
  trace.df <- do.call(rbind, trace.rows)
  if (!isTRUE(return.trace)) {
    frames <- list(current)
  }
  edge.scale <- if (nrow(trace.df)) trace.df$edge.scale[[nrow(trace.df)]] else 1
  gram.diag <- grip.kernel.gram.score(current, star, edge.scale = edge.scale, lambda.gram = lambda.gram)
  diag <- if (isTRUE(diagnostics)) {
    cbind(
      score.gmds(
        coords = current,
        prepared = prepared,
        scale.mode = if (identical(scale.mode, "identity")) "identity" else "profiled",
        distance.floor = distance.floor,
        edge.length.epsilon = edge.length.epsilon
      ),
      gram.diag
    )
  } else {
    NULL
  }
  gmds.result(
    coords = current,
    method = "kernel_gram_gkk",
    prepared = prepared,
    trace = trace.df,
    diagnostics = diag,
    metadata = list(
      engine = "r_gradient_descent_armijo",
      initialization = if (is.null(coords)) init else "supplied",
      star = star,
      angle.power = star$angle.power,
      reliability = star$reliability,
      star.quantile = star$star.quantile,
      lambda.edge = lambda.edge,
      lambda.gram = lambda.gram,
      stiffness = stiff,
      frames = frames,
      scale_mode = scale.mode
    )
  )
}
