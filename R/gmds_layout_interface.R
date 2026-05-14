grip.gmds.require.prepared <- function(prepared = NULL,
                                       coords = NULL,
                                       edges = NULL,
                                       n = NULL,
                                       adj_list = NULL,
                                       weight_list = NULL,
                                       edge_weights = NULL) {
  if (!is.null(prepared)) {
    return(grip.validate.geodesic.mds.prepared(prepared, coords = coords))
  }
  grip.prepare.graph.geodesic.mds(
    edges = edges,
    n = if (is.null(n) && !is.null(coords)) nrow(coords) else n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights
  )
}

grip.gmds.fit.scale.unweighted <- function(observed,
                                           target,
                                           distance_floor = 1e-8) {
  observed <- as.double(observed)
  target <- as.double(target)
  keep <- is.finite(observed) & is.finite(target) & target > distance_floor
  if (!any(keep)) {
    return(NA_real_)
  }
  denom <- sum(target[keep]^2)
  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }
  sum(observed[keep] * target[keep]) / denom
}

grip.gmds.resolve.scale <- function(observed,
                                    target,
                                    scale_mode = c("profiled", "identity", "user"),
                                    scale = NULL,
                                    distance_floor = 1e-8) {
  scale_mode <- match.arg(scale_mode)
  if (identical(scale_mode, "identity")) {
    return(1.0)
  }
  if (identical(scale_mode, "user")) {
    grip.validate.scalar(scale, "scale", lower = 0, open.lower = TRUE)
    return(as.double(scale))
  }
  grip.gmds.fit.scale.unweighted(
    observed = observed,
    target = target,
    distance_floor = distance_floor
  )
}

grip.gmds.residual.stats <- function(observed,
                                     target,
                                     weights = NULL,
                                     distance_floor = 1e-8) {
  observed <- as.double(observed)
  target <- as.double(target)
  if (is.null(weights)) {
    weights <- rep(1, length(target))
  } else {
    weights <- as.double(weights)
  }
  keep <- is.finite(observed) & is.finite(target) & is.finite(weights) &
    weights > 0 & target > distance_floor
  if (!any(keep)) {
    return(list(
      n = 0L,
      stress = NA_real_,
      rmse = NA_real_,
      mean_abs_error = NA_real_,
      mean_rel_abs_error = NA_real_,
      signed_bias = NA_real_,
      shortcut_fraction = NA_real_
    ))
  }
  resid <- observed[keep] - target[keep]
  rel <- resid / pmax(target[keep], distance_floor)
  w <- weights[keep]
  denom <- sum(w * target[keep]^2)
  list(
    n = sum(keep),
    stress = if (is.finite(denom) && denom > 0) sqrt(sum(w * resid^2) / denom) else NA_real_,
    rmse = sqrt(sum(w * resid^2) / sum(w)),
    mean_abs_error = sum(w * abs(resid)) / sum(w),
    mean_rel_abs_error = sum(w * abs(rel)) / sum(w),
    signed_bias = sum(w * rel) / sum(w),
    shortcut_fraction = mean(rel < -0.05)
  )
}

grip.gmds.edge.lengths <- function(coords, prepared) {
  if (is.null(prepared$edges) || nrow(prepared$edges) == 0L) {
    return(numeric(0L))
  }
  diffs <- coords[prepared$edges[, 1L], , drop = FALSE] -
    coords[prepared$edges[, 2L], , drop = FALSE]
  sqrt(rowSums(diffs^2))
}

grip.gmds.pair.chord.lengths <- function(coords, prepared) {
  if (is.null(prepared$pair_matrix) || nrow(prepared$pair_matrix) == 0L) {
    return(numeric(0L))
  }
  diffs <- coords[prepared$pair_matrix[, 1L], , drop = FALSE] -
    coords[prepared$pair_matrix[, 2L], , drop = FALSE]
  sqrt(rowSums(diffs^2))
}

grip.gmds.band.weights <- function(graph.distances,
                                   band = c("all", "short", "mid", "long"),
                                   quantiles = c(1 / 3, 2 / 3)) {
  band <- match.arg(band)
  g <- as.double(graph.distances)
  if (identical(band, "all")) {
    return(rep(1, length(g)))
  }
  if (length(g) == 0L || all(!is.finite(g))) {
    return(numeric(length(g)))
  }
  q <- as.double(stats::quantile(g[is.finite(g)], probs = quantiles, names = FALSE))
  if (identical(band, "short")) {
    as.double(g <= q[[1L]])
  } else if (identical(band, "mid")) {
    as.double(g > q[[1L]] & g <= q[[2L]])
  } else {
    as.double(g > q[[2L]])
  }
}

#' Score a layout with the common GMDS diagnostic panel
#'
#' `grip.score.gmds.layout()` is the Phase 0 common scorer for GMDS-oriented
#' layout experiments. It summarizes edge-length fidelity, fixed-path geodesic
#' stress, short/mid/long geodesic bands, ordinary metric-MDS chord stress, and
#' simple spread/shortcut diagnostics for any coordinate matrix on a prepared
#' graph.
#'
#' @param coords Numeric coordinate matrix.
#' @param prepared Optional object returned by
#'   [grip.prepare.graph.geodesic.mds()] or [grip.prepare.geodesic.kk()].
#' @param edges Two-column edge matrix used when `prepared` is omitted.
#' @param n Number of vertices used when `prepared` is omitted.
#' @param adj_list Optional adjacency list used when `prepared` is omitted.
#' @param weight_list Optional edge-weight list parallel to `adj_list`.
#' @param edge_weights Optional positive edge weights parallel to `edges`.
#' @param scale_mode Scale policy for edge, path, and chord targets:
#'   `"profiled"` fits one scalar for each diagnostic family, `"identity"`
#'   uses scale one, and `"user"` uses the corresponding supplied scale.
#' @param edge_scale,path_scale,chord_scale Optional user scales used when
#'   `scale_mode = "user"`.
#' @param distance_floor Positive floor for relative residuals.
#' @param edge_length_epsilon Small stabilizer for fixed-path embedded lengths.
#' @param band_quantiles Two quantiles splitting graph distances into short,
#'   mid, and long bands.
#'
#' @return A one-row data frame with the common diagnostic columns.
#' @export
grip.score.gmds.layout <- function(coords,
                                   prepared = NULL,
                                   edges = NULL,
                                   n = NULL,
                                   adj_list = NULL,
                                   weight_list = NULL,
                                   edge_weights = NULL,
                                   scale_mode = c("profiled", "identity", "user"),
                                   edge_scale = NULL,
                                   path_scale = NULL,
                                   chord_scale = NULL,
                                   distance_floor = 1e-8,
                                   edge_length_epsilon = 1e-8,
                                   band_quantiles = c(1 / 3, 2 / 3)) {
  coords <- grip.validate.coords(coords)
  scale_mode <- match.arg(scale_mode)
  grip.validate.scalar(distance_floor, "distance_floor", lower = 0, open.lower = TRUE)
  grip.validate.scalar(edge_length_epsilon, "edge_length_epsilon", lower = 0)
  if (!is.numeric(band_quantiles) || length(band_quantiles) != 2L ||
      any(!is.finite(band_quantiles)) ||
      any(band_quantiles <= 0 | band_quantiles >= 1) ||
      band_quantiles[[1L]] >= band_quantiles[[2L]]) {
    stop("band_quantiles must be two increasing probabilities in (0, 1)")
  }

  prepared <- grip.gmds.require.prepared(
    prepared = prepared,
    coords = coords,
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights
  )

  edge.observed <- grip.gmds.edge.lengths(coords, prepared)
  edge.target.raw <- if (!is.null(prepared$edge_targets)) {
    as.double(prepared$edge_targets)
  } else {
    numeric(0L)
  }
  edge.scale <- grip.gmds.resolve.scale(
    observed = edge.observed,
    target = edge.target.raw,
    scale_mode = scale_mode,
    scale = edge_scale,
    distance_floor = distance_floor
  )
  edge.stats <- grip.gmds.residual.stats(
    observed = edge.observed,
    target = edge.scale * edge.target.raw,
    distance_floor = distance_floor
  )

  path.observed <- grip.geodesic.kk.path.lengths(
    coords = coords,
    prepared = prepared,
    edge_length_epsilon = edge_length_epsilon
  )
  path.target.raw <- as.double(prepared$pair_graph_distance)
  path.scale <- grip.gmds.resolve.scale(
    observed = path.observed,
    target = path.target.raw,
    scale_mode = scale_mode,
    scale = path_scale,
    distance_floor = distance_floor
  )
  path.target <- path.scale * path.target.raw
  path.stats <- grip.gmds.residual.stats(
    observed = path.observed,
    target = path.target,
    distance_floor = distance_floor
  )
  short.stats <- grip.gmds.residual.stats(
    observed = path.observed,
    target = path.target,
    weights = grip.gmds.band.weights(path.target.raw, "short", band_quantiles),
    distance_floor = distance_floor
  )
  mid.stats <- grip.gmds.residual.stats(
    observed = path.observed,
    target = path.target,
    weights = grip.gmds.band.weights(path.target.raw, "mid", band_quantiles),
    distance_floor = distance_floor
  )
  long.stats <- grip.gmds.residual.stats(
    observed = path.observed,
    target = path.target,
    weights = grip.gmds.band.weights(path.target.raw, "long", band_quantiles),
    distance_floor = distance_floor
  )

  chord.observed <- grip.gmds.pair.chord.lengths(coords, prepared)
  chord.scale <- grip.gmds.resolve.scale(
    observed = chord.observed,
    target = path.target.raw,
    scale_mode = scale_mode,
    scale = chord_scale,
    distance_floor = distance_floor
  )
  chord.stats <- grip.gmds.residual.stats(
    observed = chord.observed,
    target = chord.scale * path.target.raw,
    distance_floor = distance_floor
  )

  spread.score <- if (length(chord.observed) > 0L && length(path.target.raw) > 0L) {
    chord.rms <- sqrt(mean(chord.observed^2))
    graph.rms <- sqrt(mean(path.target.raw^2))
    if (is.finite(graph.rms) && graph.rms > 0) chord.rms / graph.rms else NA_real_
  } else {
    NA_real_
  }

  data.frame(
    n.vertices = prepared$n,
    n.edges = length(edge.observed),
    n.pairs = length(path.target.raw),
    scale.mode = scale_mode,
    edge.scale = edge.scale,
    edge.rel.rmse = edge.stats$stress,
    edge.rmse = edge.stats$rmse,
    edge.mean.abs.error = edge.stats$mean_abs_error,
    edge.mean.rel.abs.error = edge.stats$mean_rel_abs_error,
    edge.signed.bias = edge.stats$signed_bias,
    gmds.scale = path.scale,
    gmds.stress = path.stats$stress,
    gmds.rmse = path.stats$rmse,
    gmds.mean.abs.error = path.stats$mean_abs_error,
    gmds.mean.rel.abs.error = path.stats$mean_rel_abs_error,
    gmds.signed.bias = path.stats$signed_bias,
    gmds.short.stress = short.stats$stress,
    gmds.mid.stress = mid.stats$stress,
    gmds.long.stress = long.stats$stress,
    gmds.short.signed.bias = short.stats$signed_bias,
    gmds.mid.signed.bias = mid.stats$signed_bias,
    gmds.long.signed.bias = long.stats$signed_bias,
    shortcut.fraction = path.stats$shortcut_fraction,
    metric.chord.scale = chord.scale,
    metric.chord.stress = chord.stats$stress,
    metric.chord.rmse = chord.stats$rmse,
    metric.chord.mean.rel.abs.error = chord.stats$mean_rel_abs_error,
    spread.score = spread.score,
    stringsAsFactors = FALSE
  )
}

#' Construct a common GMDS layout result
#'
#' `grip.gmds.layout.result()` wraps coordinates, method metadata, optional
#' traces, and optional diagnostics in a common result shape used by the
#' experimental GMDS layout program.
#'
#' @param coords Numeric coordinate matrix.
#' @param method Short method identifier.
#' @param prepared Optional prepared graph object.
#' @param trace Optional trace table or list.
#' @param diagnostics Optional one-row diagnostic data frame.
#' @param metadata Optional named list with method-specific metadata.
#'
#' @return A list of class `"grip_gmds_layout"` with fields `coords`, `method`,
#'   `prepared`, `trace`, `diagnostics`, and `metadata`.
#' @export
grip.gmds.layout.result <- function(coords,
                                    method,
                                    prepared = NULL,
                                    trace = NULL,
                                    diagnostics = NULL,
                                    metadata = list()) {
  coords <- grip.validate.coords(coords)
  if (!is.character(method) || length(method) != 1L || is.na(method) || !nzchar(method)) {
    stop("method must be a non-empty character scalar")
  }
  if (!is.null(prepared)) {
    prepared <- grip.validate.geodesic.mds.prepared(prepared, coords = coords)
  }
  if (!is.null(diagnostics) && !is.data.frame(diagnostics)) {
    stop("diagnostics must be NULL or a data frame")
  }
  if (!is.list(metadata)) {
    stop("metadata must be a list")
  }
  out <- list(
    coords = coords,
    method = method,
    prepared = prepared,
    trace = trace,
    diagnostics = diagnostics,
    metadata = metadata
  )
  class(out) <- c("grip_gmds_layout", "list")
  out
}

#' @export
print.grip_gmds_layout <- function(x, ...) {
  cat("<grip_gmds_layout>\n")
  cat("  method:", x$method, "\n")
  cat("  vertices:", nrow(x$coords), "\n")
  cat("  dimension:", ncol(x$coords), "\n")
  if (!is.null(x$diagnostics) && nrow(x$diagnostics) > 0L) {
    if ("edge.rel.rmse" %in% names(x$diagnostics)) {
      cat("  edge.rel.rmse:", format(x$diagnostics$edge.rel.rmse[[1L]], digits = 4), "\n")
    }
    if ("gmds.stress" %in% names(x$diagnostics)) {
      cat("  gmds.stress:", format(x$diagnostics$gmds.stress[[1L]], digits = 4), "\n")
    }
  }
  invisible(x)
}

#' Metric-MDS baseline layout for GMDS experiments
#'
#' `grip.metric.mds.layout()` computes the current metric-MDS baseline using
#' `stats::cmdscale()` on the graph shortest-path distance matrix and returns a
#' common `"grip_gmds_layout"` result.
#'
#' @inheritParams grip.score.gmds.layout
#' @param dim Target embedding dimension, currently 2 or 3.
#' @param add,eig Passed to `stats::cmdscale()`.
#' @param diagnostics If `TRUE`, attach the common GMDS diagnostic panel.
#'
#' @return A `"grip_gmds_layout"` object.
#' @export
grip.metric.mds.layout <- function(prepared = NULL,
                                   edges = NULL,
                                   n = NULL,
                                   adj_list = NULL,
                                   weight_list = NULL,
                                   edge_weights = NULL,
                                   dim = 2L,
                                   add = FALSE,
                                   eig = TRUE,
                                   diagnostics = TRUE,
                                   scale_mode = c("profiled", "identity", "user"),
                                   distance_floor = 1e-8,
                                   edge_length_epsilon = 1e-8,
                                   band_quantiles = c(1 / 3, 2 / 3)) {
  scale_mode <- match.arg(scale_mode)
  prepared <- grip.gmds.require.prepared(
    prepared = prepared,
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights
  )
  fit <- grip.classical.mds.embedding(
    prepared = prepared,
    dim = dim,
    add = add,
    eig = eig
  )
  diag <- if (isTRUE(diagnostics)) {
    grip.score.gmds.layout(
      coords = fit$coords,
      prepared = prepared,
      scale_mode = scale_mode,
      distance_floor = distance_floor,
      edge_length_epsilon = edge_length_epsilon,
      band_quantiles = band_quantiles
    )
  } else {
    NULL
  }
  grip.gmds.layout.result(
    coords = fit$coords,
    method = "metric_mds",
    prepared = prepared,
    trace = NULL,
    diagnostics = diag,
    metadata = list(
      engine = "cmdscale",
      eig = fit$eig,
      additive_constant = fit$additive_constant,
      gof = fit$gof,
      positive_eigen_fraction = fit$positive_eigen_fraction,
      negative_eigen_fraction = fit$negative_eigen_fraction
    )
  )
}

grip.validate.infinite.ceiling <- function(x, name, lower = 0) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= lower) {
    stop(sprintf("%s must be a numeric scalar greater than %s", name, format(lower, digits = 16)))
  }
  if (!is.finite(x) && !is.infinite(x)) {
    stop(sprintf("%s must be finite or Inf", name))
  }
}

grip.normalize.mean.one <- function(x, name) {
  x <- as.double(x)
  if (length(x) == 0L) {
    return(x)
  }
  m <- mean(x)
  if (!is.finite(m) || m <= 0) {
    stop(sprintf("%s could not be normalized to mean one", name))
  }
  x / m
}

grip.apply.stiffness.transform <- function(x, transform) {
  x <- pmax(as.double(x), 0)
  switch(
    transform,
    identity = x,
    sqrt = sqrt(x),
    log = log1p(x)
  )
}

#' Construct edge-length stiffnesses for edge-isometric GMDS layouts
#'
#' `grip.edge.length.density.stiffness()` is the Phase 1 constructor for
#' edge-only geodesic-MDS experiments. It turns positive edge lengths into
#' spring stiffnesses normalized to mean one. The `"density"` method emphasizes
#' edge lengths near the empirical edge-length density mode, and `mix` provides
#' a continuation path from that density-weighted signal to uniform stiffness.
#'
#' @param edge_weights Positive numeric edge lengths.
#' @param method Stiffness rule. `"density"` estimates an empirical density,
#'   `"uniform"` returns equal stiffnesses, and `"distance_power"` uses
#'   `(w / median(w))^distance_power`.
#' @param mix Continuation parameter in `[0, 1]`. `0` uses the selected method;
#'   `1` returns uniform stiffness.
#' @param bandwidth Optional bandwidth passed to `stats::density()`.
#' @param density_n Number of evaluation points for `stats::density()`.
#' @param transform Optional transformation of the raw density/power signal
#'   before mixing and normalization.
#' @param distance_power Exponent for `method = "distance_power"`.
#' @param stiffness_floor,stiffness_ceiling Optional clipping bounds applied
#'   before the final mean-one normalization.
#'
#' @return A list with `stiffness`, raw signal diagnostics, estimated mode, and
#'   clipping/normalization metadata.
#' @export
grip.edge.length.density.stiffness <- function(edge_weights,
                                               method = c("density", "uniform", "distance_power"),
                                               mix = 0,
                                               bandwidth = NULL,
                                               density_n = 512L,
                                               transform = c("identity", "sqrt", "log"),
                                               distance_power = 0,
                                               stiffness_floor = 0,
                                               stiffness_ceiling = Inf) {
  method <- match.arg(method)
  transform <- match.arg(transform)
  if (!is.numeric(edge_weights)) {
    stop("edge_weights must be a numeric vector")
  }
  edge_weights <- as.double(edge_weights)
  bad <- which(!is.finite(edge_weights) | edge_weights <= 0)
  if (length(bad) > 0L) {
    stop(sprintf("edge_weights must contain finite values > 0; first invalid at edge_weights[%d]", bad[[1L]]))
  }
  grip.validate.scalar(mix, "mix", lower = 0, upper = 1)
  if (!is.null(bandwidth)) {
    grip.validate.scalar(bandwidth, "bandwidth", lower = 0, open.lower = TRUE)
  }
  grip.validate.scalar(density_n, "density_n", lower = 16)
  density_n <- as.integer(round(density_n))
  grip.validate.scalar(distance_power, "distance_power")
  grip.validate.scalar(stiffness_floor, "stiffness_floor", lower = 0)
  grip.validate.infinite.ceiling(stiffness_ceiling, "stiffness_ceiling", lower = 0)
  if (is.finite(stiffness_ceiling) && stiffness_ceiling < stiffness_floor) {
    stop("stiffness_ceiling must be greater than or equal to stiffness_floor")
  }

  n.edge <- length(edge_weights)
  if (n.edge == 0L) {
    return(list(
      stiffness = numeric(0L),
      raw_signal = numeric(0L),
      mixed_signal = numeric(0L),
      edge_weights = edge_weights,
      mode = NA_real_,
      method = method,
      mix = mix,
      transform = transform,
      bandwidth = bandwidth,
      diagnostics = data.frame(
        n.edges = 0L,
        raw.min = NA_real_,
        raw.max = NA_real_,
        stiffness.min = NA_real_,
        stiffness.max = NA_real_,
        stiffness.mean = NA_real_,
        stringsAsFactors = FALSE
      )
    ))
  }

  mode <- stats::median(edge_weights)
  raw <- rep(1, n.edge)
  density.grid <- NULL
  if (identical(method, "density")) {
    if (n.edge >= 2L && diff(range(edge_weights)) > 0) {
      dens <- stats::density(
        edge_weights,
        bw = if (is.null(bandwidth)) "nrd0" else bandwidth,
        n = density_n,
        from = min(edge_weights),
        to = max(edge_weights)
      )
      raw <- stats::approx(dens$x, dens$y, xout = edge_weights, rule = 2)$y
      raw[!is.finite(raw) | raw < 0] <- 0
      mode <- dens$x[which.max(dens$y)]
      density.grid <- data.frame(x = dens$x, y = dens$y)
      if (is.null(bandwidth)) {
        bandwidth <- dens$bw
      }
    }
  } else if (identical(method, "distance_power")) {
    center <- stats::median(edge_weights)
    if (!is.finite(center) || center <= 0) {
      center <- mean(edge_weights)
    }
    raw <- (edge_weights / center)^distance_power
    raw[!is.finite(raw) | raw < 0] <- 0
  }

  raw <- grip.apply.stiffness.transform(raw, transform)
  raw <- grip.normalize.mean.one(raw, "raw stiffness signal")
  mixed <- (1 - mix) * raw + mix
  if (stiffness_floor > 0) {
    mixed <- pmax(mixed, stiffness_floor)
  }
  if (is.finite(stiffness_ceiling)) {
    mixed <- pmin(mixed, stiffness_ceiling)
  }
  stiffness <- grip.normalize.mean.one(mixed, "stiffness")

  out <- list(
    stiffness = stiffness,
    raw_signal = raw,
    mixed_signal = mixed,
    edge_weights = edge_weights,
    mode = as.double(mode),
    method = method,
    mix = mix,
    transform = transform,
    bandwidth = bandwidth,
    distance_power = distance_power,
    stiffness_floor = stiffness_floor,
    stiffness_ceiling = stiffness_ceiling,
    density = density.grid,
    diagnostics = data.frame(
      n.edges = n.edge,
      raw.min = min(raw),
      raw.max = max(raw),
      stiffness.min = min(stiffness),
      stiffness.max = max(stiffness),
      stiffness.mean = mean(stiffness),
      stringsAsFactors = FALSE
    )
  )
  class(out) <- c("grip_edge_length_stiffness", "list")
  out
}

grip.edge.isometric.fit.scale <- function(observed,
                                          target,
                                          stiffness = NULL,
                                          distance_floor = 1e-8) {
  observed <- as.double(observed)
  target <- as.double(target)
  if (is.null(stiffness)) {
    stiffness <- rep(1, length(target))
  } else {
    stiffness <- as.double(stiffness)
  }
  keep <- is.finite(observed) & is.finite(target) & is.finite(stiffness) &
    stiffness > 0 & target > distance_floor
  if (!any(keep)) {
    return(NA_real_)
  }
  denom <- sum(stiffness[keep] * target[keep]^2)
  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }
  sum(stiffness[keep] * observed[keep] * target[keep]) / denom
}

grip.edge.isometric.energy.gradient <- function(coords,
                                                edges,
                                                edge_weights,
                                                stiffness,
                                                scale = 1,
                                                edge_length_epsilon = 1e-8) {
  coords <- grip.validate.coords(coords)
  edge_weights <- as.double(edge_weights)
  stiffness <- as.double(stiffness)
  scale <- as.double(scale)
  grip.validate.scalar(scale, "scale", lower = 0, open.lower = TRUE)
  grip.validate.scalar(edge_length_epsilon, "edge_length_epsilon", lower = 0)
  if (!is.matrix(edges) || ncol(edges) != 2L) {
    stop("edges must be a two-column matrix")
  }
  if (length(edge_weights) != nrow(edges) || length(stiffness) != nrow(edges)) {
    stop("edge_weights and stiffness must be parallel to edges")
  }
  if (nrow(edges) == 0L) {
    return(list(
      energy = 0,
      gradient = matrix(0, nrow(coords), ncol(coords)),
      gradient_norm = 0,
      edge_lengths = numeric(0L),
      residuals = numeric(0L),
      scale = scale
    ))
  }

  diffs <- coords[edges[, 1L], , drop = FALSE] - coords[edges[, 2L], , drop = FALSE]
  edge.lengths <- sqrt(rowSums(diffs^2) + edge_length_epsilon^2)
  residuals <- edge.lengths - scale * edge_weights
  coeff <- stiffness * residuals / edge.lengths
  gradient <- matrix(0, nrow(coords), ncol(coords))
  for (e in seq_len(nrow(edges))) {
    contribution <- coeff[[e]] * diffs[e, ]
    gradient[edges[e, 1L], ] <- gradient[edges[e, 1L], ] + contribution
    gradient[edges[e, 2L], ] <- gradient[edges[e, 2L], ] - contribution
  }
  energy <- 0.5 * sum(stiffness * residuals^2)
  list(
    energy = energy,
    gradient = gradient,
    gradient_norm = sqrt(sum(gradient^2)),
    edge_lengths = edge.lengths,
    residuals = residuals,
    scale = scale
  )
}

grip.edge.isometric.evaluate.state <- function(coords,
                                               edges,
                                               edge_weights,
                                               stiffness,
                                               scale_mode = c("profiled", "identity", "fixed", "user"),
                                               scale = NULL,
                                               edge_length_epsilon = 1e-8,
                                               distance_floor = 1e-8) {
  scale_mode <- match.arg(scale_mode)
  edge.lengths <- if (nrow(edges) == 0L) {
    numeric(0L)
  } else {
    diffs <- coords[edges[, 1L], , drop = FALSE] - coords[edges[, 2L], , drop = FALSE]
    sqrt(rowSums(diffs^2) + edge_length_epsilon^2)
  }
  edge.scale <- switch(
    scale_mode,
    profiled = grip.edge.isometric.fit.scale(
      observed = edge.lengths,
      target = edge_weights,
      stiffness = stiffness,
      distance_floor = distance_floor
    ),
    identity = 1.0,
    fixed = {
      grip.validate.scalar(scale, "scale", lower = 0, open.lower = TRUE)
      as.double(scale)
    },
    user = {
      grip.validate.scalar(scale, "scale", lower = 0, open.lower = TRUE)
      as.double(scale)
    }
  )
  if (!is.finite(edge.scale)) {
    edge.scale <- 1.0
  }
  state <- grip.edge.isometric.energy.gradient(
    coords = coords,
    edges = edges,
    edge_weights = edge_weights,
    stiffness = stiffness,
    scale = edge.scale,
    edge_length_epsilon = edge_length_epsilon
  )
  target <- edge.scale * edge_weights
  edge.stats <- grip.gmds.residual.stats(
    observed = state$edge_lengths,
    target = target,
    weights = stiffness,
    distance_floor = distance_floor
  )
  state$edge_scale <- edge.scale
  state$edge_rel_rmse <- edge.stats$stress
  state$edge_signed_bias <- edge.stats$signed_bias
  state
}

grip.edge.isometric.initial.coords <- function(prepared,
                                               coords = NULL,
                                               init = c("metric_mds", "random"),
                                               dim = 2L,
                                               seed = 1L) {
  grip.validate.scalar(dim, "dim", lower = 2, upper = 3)
  dim <- as.integer(round(dim))
  if (!is.null(coords)) {
    coords <- grip.validate.coords(coords)
    if (nrow(coords) != prepared$n) {
      stop("nrow(coords) must match the graph size stored in prepared")
    }
    if (ncol(coords) != dim) {
      stop("ncol(coords) must match dim")
    }
    return(coords)
  }
  init <- match.arg(init)
  if (identical(init, "metric_mds")) {
    return(grip.metric.mds.layout(
      prepared = prepared,
      dim = dim,
      diagnostics = FALSE
    )$coords)
  }
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }
  out <- matrix(stats::rnorm(prepared$n * dim), nrow = prepared$n, ncol = dim)
  sweep(out, 2L, colMeans(out), "-", check.margin = FALSE)
}

#' Optimize an edge-isometric GMDS layout
#'
#' `grip.optimize.edge.isometric.layout()` is the edge-KK local repair operator
#' for the experimental GMDS layout program. Earlier notes called this operator
#' edge-gKK; edge-KK is the preferred name because the objective is restricted
#' to graph edges rather than all graph-geodesic pairs. It minimizes
#' weighted edge-length stress
#' \deqn{
#'   \frac{1}{2}\sum_{(i,j)\in E} k_{ij}
#'   \left(\|z_i-z_j\|_2 - s w_{ij}\right)^2
#' }
#' using deterministic gradient descent with Armijo backtracking. The default
#' `density_mix_schedule` runs a continuation from density-weighted stiffnesses
#' toward uniform stiffnesses.
#'
#' @inheritParams grip.score.gmds.layout
#' @param coords Optional starting coordinates. If omitted, `init` is used.
#' @param dim Target embedding dimension.
#' @param init Starting layout used when `coords` is omitted.
#' @param stiffness_method,stiffness_transform,density_mix_schedule,bandwidth,density_n
#'   Parameters passed to [grip.edge.length.density.stiffness()].
#' @param distance_power,stiffness_floor,stiffness_ceiling Additional stiffness
#'   constructor parameters.
#' @param scale_mode Scale policy for edge targets. `"profiled"` analytically
#'   refits `s` at every state evaluation, `"identity"` fixes `s = 1`,
#'   `"fixed_initial"` fits `s` once at the first continuation stage, and
#'   `"user"` uses `scale`.
#' @param scale User scale for `scale_mode = "user"`.
#' @param max_iter Maximum iterations per continuation stage.
#' @param initial_step,step_shrink,armijo_factor,grad_tol,min_step Line-search
#'   controls.
#' @param recenter If `TRUE`, recenter the layout after accepted steps.
#' @param return_trace If `TRUE`, keep per-iteration trace rows.
#' @param diagnostics If `TRUE`, attach the common GMDS diagnostic panel.
#' @param seed Random seed used only for `init = "random"`.
#' @param engine Optimizer engine. `"cpp"` uses the Rcpp backend for the hot
#'   edge-stress loop; `"R"` keeps the reference prototype.
#'
#' @return A `"grip_gmds_layout"` object with method `"edge_isometric_gkk"`
#'   for compatibility with existing experiment assets.
#' @export
grip.optimize.edge.isometric.layout <- function(coords = NULL,
                                                prepared = NULL,
                                                edges = NULL,
                                                n = NULL,
                                                adj_list = NULL,
                                                weight_list = NULL,
                                                edge_weights = NULL,
                                                dim = 2L,
                                                init = c("metric_mds", "random"),
                                                stiffness_method = c("density", "uniform", "distance_power"),
                                                stiffness_transform = c("identity", "sqrt", "log"),
                                                density_mix_schedule = c(0, 0.25, 0.5, 0.75, 1),
                                                bandwidth = NULL,
                                                density_n = 512L,
                                                distance_power = 0,
                                                stiffness_floor = 0,
                                                stiffness_ceiling = Inf,
                                                scale_mode = c("profiled", "identity", "fixed_initial", "user"),
                                                scale = NULL,
                                                max_iter = 50L,
                                                initial_step = 1.0,
                                                step_shrink = 0.5,
                                                armijo_factor = 1e-4,
                                                grad_tol = 1e-8,
                                                min_step = 1e-8,
                                                edge_length_epsilon = 1e-8,
                                                distance_floor = 1e-8,
                                                recenter = TRUE,
                                                return_trace = TRUE,
                                                diagnostics = TRUE,
                                                seed = 1L,
                                                engine = c("cpp", "R")) {
  init <- match.arg(init)
  stiffness_method <- match.arg(stiffness_method)
  stiffness_transform <- match.arg(stiffness_transform)
  scale_mode <- match.arg(scale_mode)
  engine <- match.arg(engine)
  prepared <- grip.gmds.require.prepared(
    prepared = prepared,
    coords = coords,
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights
  )
  if (is.null(prepared$edges) || is.null(prepared$edge_targets)) {
    stop("prepared object must contain edges and edge_targets")
  }
  grip.validate.scalar(dim, "dim", lower = 2, upper = 3)
  dim <- as.integer(round(dim))
  if (!is.numeric(density_mix_schedule) || length(density_mix_schedule) < 1L ||
      any(!is.finite(density_mix_schedule)) ||
      any(density_mix_schedule < 0 | density_mix_schedule > 1)) {
    stop("density_mix_schedule must contain values in [0, 1]")
  }
  density_mix_schedule <- as.double(density_mix_schedule)
  grip.validate.scalar(max_iter, "max_iter", lower = 0)
  max_iter <- as.integer(round(max_iter))
  grip.validate.scalar(initial_step, "initial_step", lower = 0, open.lower = TRUE)
  grip.validate.scalar(step_shrink, "step_shrink", lower = 0, upper = 1, open.lower = TRUE, open.upper = TRUE)
  grip.validate.scalar(armijo_factor, "armijo_factor", lower = 0)
  grip.validate.scalar(grad_tol, "grad_tol", lower = 0)
  grip.validate.scalar(min_step, "min_step", lower = 0, open.lower = TRUE)
  grip.validate.scalar(edge_length_epsilon, "edge_length_epsilon", lower = 0)
  grip.validate.scalar(distance_floor, "distance_floor", lower = 0, open.lower = TRUE)
  if (!is.logical(recenter) || length(recenter) != 1L || is.na(recenter)) {
    stop("recenter must be TRUE or FALSE")
  }
  if (!is.logical(return_trace) || length(return_trace) != 1L || is.na(return_trace)) {
    stop("return_trace must be TRUE or FALSE")
  }
  if (!is.logical(diagnostics) || length(diagnostics) != 1L || is.na(diagnostics)) {
    stop("diagnostics must be TRUE or FALSE")
  }
  if (identical(scale_mode, "user")) {
    grip.validate.scalar(scale, "scale", lower = 0, open.lower = TRUE)
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
  stiffness.objects <- lapply(density_mix_schedule, function(mix) {
    grip.edge.length.density.stiffness(
      edge_weights = edge.targets,
      method = stiffness_method,
      mix = mix,
      bandwidth = bandwidth,
      density_n = density_n,
      transform = stiffness_transform,
      distance_power = distance_power,
      stiffness_floor = stiffness_floor,
      stiffness_ceiling = stiffness_ceiling
    )
  })
  stiffness.matrix <- do.call(cbind, lapply(stiffness.objects, `[[`, "stiffness"))

  if (identical(engine, "cpp")) {
    cpp.scale.mode <- switch(
      scale_mode,
      profiled = "profiled",
      identity = "identity",
      fixed_initial = "fixed",
      user = "user"
    )
    cpp.scale <- NA_real_
    if (identical(scale_mode, "user")) {
      cpp.scale <- as.double(scale)
    } else if (identical(scale_mode, "fixed_initial")) {
      initial.lengths <- if (nrow(edges) == 0L) {
        numeric(0L)
      } else {
        diffs <- current[edges[, 1L], , drop = FALSE] - current[edges[, 2L], , drop = FALSE]
        sqrt(rowSums(diffs^2) + edge_length_epsilon^2)
      }
      cpp.scale <- grip.edge.isometric.fit.scale(
        observed = initial.lengths,
        target = edge.targets,
        stiffness = stiffness.matrix[, 1L],
        distance_floor = distance_floor
      )
      if (!is.finite(cpp.scale)) {
        cpp.scale <- 1.0
      }
    }

    fit <- grip_optimize_edge_isometric_layout_cpp(
      edges = matrix(as.integer(edges), ncol = 2L),
      edge_weights = edge.targets,
      stiffness_matrix = stiffness.matrix,
      mix_schedule = density_mix_schedule,
      coords = current,
      max_iter = max_iter,
      scale_mode = cpp.scale.mode,
      scale = cpp.scale,
      edge_length_epsilon = edge_length_epsilon,
      initial_step = initial_step,
      step_shrink = step_shrink,
      armijo_factor = armijo_factor,
      grad_tol = grad_tol,
      min_step = min_step,
      distance_floor = distance_floor,
      recenter = recenter,
      return_trace = return_trace
    )
    current <- fit$coords
    trace.df <- fit$trace
    stage.summaries <- lapply(seq_along(stiffness.objects), function(stage.idx) {
      stiff <- stiffness.objects[[stage.idx]]
      rows <- trace.df[trace.df$stage == stage.idx, , drop = FALSE]
      final <- rows[nrow(rows), , drop = FALSE]
      data.frame(
        stage = stage.idx,
        mix = density_mix_schedule[[stage.idx]],
        method = stiff$method,
        transform = stiff$transform,
        mode = stiff$mode,
        stiffness.min = min(stiff$stiffness),
        stiffness.max = max(stiff$stiffness),
        final.energy = final$energy[[1L]],
        final.edge.scale = final$edge.scale[[1L]],
        final.edge.rel.rmse = final$edge.rel.rmse[[1L]],
        stringsAsFactors = FALSE
      )
    })
    diag <- if (isTRUE(diagnostics)) {
      grip.score.gmds.layout(
        coords = current,
        prepared = prepared,
        scale_mode = if (identical(scale_mode, "identity")) "identity" else "profiled",
        distance_floor = distance_floor,
        edge_length_epsilon = edge_length_epsilon
      )
    } else {
      NULL
    }
    return(grip.gmds.layout.result(
      coords = current,
      method = "edge_isometric_gkk",
      prepared = prepared,
      trace = trace.df,
      diagnostics = diag,
      metadata = list(
        engine = "cpp_gradient_descent_armijo",
        stiffness_method = stiffness_method,
        stiffness_transform = stiffness_transform,
        density_mix_schedule = density_mix_schedule,
        stage_summaries = do.call(rbind, stage.summaries),
        frames = fit$frames,
        scale_mode = scale_mode
      )
    ))
  }

  fixed.scale <- NULL
  trace.rows <- list()
  stage.summaries <- vector("list", length(density_mix_schedule))
  frames <- list(current)

  for (stage.idx in seq_along(density_mix_schedule)) {
    mix <- density_mix_schedule[[stage.idx]]
    stiff <- stiffness.objects[[stage.idx]]
    stiffness <- stiff$stiffness
    state.scale.mode <- switch(
      scale_mode,
      profiled = "profiled",
      identity = "identity",
      fixed_initial = "fixed",
      user = "user"
    )
    if (identical(scale_mode, "fixed_initial") && is.null(fixed.scale)) {
      initial.lengths <- if (nrow(edges) == 0L) {
        numeric(0L)
      } else {
        diffs <- current[edges[, 1L], , drop = FALSE] - current[edges[, 2L], , drop = FALSE]
        sqrt(rowSums(diffs^2) + edge_length_epsilon^2)
      }
      fixed.scale <- grip.edge.isometric.fit.scale(
        observed = initial.lengths,
        target = edge.targets,
        stiffness = stiffness,
        distance_floor = distance_floor
      )
      if (!is.finite(fixed.scale)) {
        fixed.scale <- 1.0
      }
    }
    stage.scale <- if (identical(scale_mode, "user")) {
      as.double(scale)
    } else if (identical(scale_mode, "fixed_initial")) {
      fixed.scale
    } else {
      NULL
    }
    state <- grip.edge.isometric.evaluate.state(
      coords = current,
      edges = edges,
      edge_weights = edge.targets,
      stiffness = stiffness,
      scale_mode = state.scale.mode,
      scale = stage.scale,
      edge_length_epsilon = edge_length_epsilon,
      distance_floor = distance_floor
    )
    trace.rows[[length(trace.rows) + 1L]] <- data.frame(
      stage = stage.idx,
      mix = mix,
      iteration = 0L,
      energy = state$energy,
      gradient_norm = state$gradient_norm,
      step = NA_real_,
      accepted = TRUE,
      edge.scale = state$edge_scale,
      edge.rel.rmse = state$edge_rel_rmse,
      stringsAsFactors = FALSE
    )

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
        proposal.state <- grip.edge.isometric.evaluate.state(
          coords = proposal,
          edges = edges,
          edge_weights = edge.targets,
          stiffness = stiffness,
          scale_mode = state.scale.mode,
          scale = stage.scale,
          edge_length_epsilon = edge_length_epsilon,
          distance_floor = distance_floor
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
      trace.rows[[length(trace.rows) + 1L]] <- data.frame(
        stage = stage.idx,
        mix = mix,
        iteration = iter,
        energy = if (accepted) candidate.state$energy else state$energy,
        gradient_norm = if (accepted) candidate.state$gradient_norm else state$gradient_norm,
        step = if (accepted) step else NA_real_,
        accepted = accepted,
        edge.scale = if (accepted) candidate.state$edge_scale else state$edge_scale,
        edge.rel.rmse = if (accepted) candidate.state$edge_rel_rmse else state$edge_rel_rmse,
        stringsAsFactors = FALSE
      )
      if (!accepted) {
        break
      }
      current <- candidate
      state <- candidate.state
      if (isTRUE(return_trace)) {
        frames[[length(frames) + 1L]] <- current
      }
    }
    stage.summaries[[stage.idx]] <- data.frame(
      stage = stage.idx,
      mix = mix,
      method = stiff$method,
      transform = stiff$transform,
      mode = stiff$mode,
      stiffness.min = min(stiffness),
      stiffness.max = max(stiffness),
      final.energy = state$energy,
      final.edge.scale = state$edge_scale,
      final.edge.rel.rmse = state$edge_rel_rmse,
      stringsAsFactors = FALSE
    )
  }

  trace.df <- do.call(rbind, trace.rows)
  if (!isTRUE(return_trace)) {
    frames <- list(current)
  }
  diag <- if (isTRUE(diagnostics)) {
    grip.score.gmds.layout(
      coords = current,
      prepared = prepared,
      scale_mode = if (identical(scale_mode, "identity")) "identity" else "profiled",
      distance_floor = distance_floor,
      edge_length_epsilon = edge_length_epsilon
    )
  } else {
    NULL
  }
  grip.gmds.layout.result(
    coords = current,
    method = "edge_isometric_gkk",
    prepared = prepared,
    trace = trace.df,
    diagnostics = diag,
    metadata = list(
      engine = "r_gradient_descent_armijo",
      stiffness_method = stiffness_method,
      stiffness_transform = stiffness_transform,
      density_mix_schedule = density_mix_schedule,
      stage_summaries = do.call(rbind, stage.summaries),
      frames = frames,
      scale_mode = scale_mode
    )
  )
}

#' Optimize an edge-KK local repair layout
#'
#' `grip.optimize.edge.kk.layout()` is the preferred public name for the
#' edge-restricted Kamada--Kawai local repair operator implemented by
#' [grip.optimize.edge.isometric.layout()]. It is a thin compatibility wrapper:
#' the optimization path is identical, but the returned `method` field is
#' normalized to `"edge_kk"` for new reports and experiment manifests.
#'
#' @inheritParams grip.optimize.edge.isometric.layout
#' @return A `"grip_gmds_layout"` object with method `"edge_kk"`.
#' @export
grip.optimize.edge.kk.layout <- function(...) {
  out <- grip.optimize.edge.isometric.layout(...)
  out$metadata$legacy_method <- out$method
  out$method <- "edge_kk"
  out
}

#' Legacy alias for edge-KK local repair
#'
#' `grip.optimize.edge.gkk.layout()` is kept as a quiet compatibility alias for
#' older notes that used the name edge-gKK. New code should use
#' [grip.optimize.edge.kk.layout()].
#'
#' @inheritParams grip.optimize.edge.kk.layout
#' @return A `"grip_gmds_layout"` object with method `"edge_kk"`.
#' @export
grip.optimize.edge.gkk.layout <- function(...) {
  grip.optimize.edge.kk.layout(...)
}
