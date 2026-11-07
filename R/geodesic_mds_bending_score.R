grip.score.geodesic.mds <- function(coords,
                                    prepared = NULL,
                                    data = NULL,
                                    k = NULL,
                                    connect = c("mst", "error"),
                                    tie.mode = c("single", "average"),
                                    edge.length.epsilon = 1e-8,
                                    anchor.coords = NULL,
                                    anchor.weight = 0,
                                    anchor.vertex.weight = NULL,
                                    smoothness.weight = 0,
                                    edge.spring.weight = 0,
                                    repulsion.weight = 0,
                                    repulsion.quantile = 0.60,
                                    repulsion.scale = 0.20,
                                    repulsion.cap.quantile = 0.90,
                                    repulsion.hop.min = 3L,
                                    bending.stencils = NULL,
                                    bending.weight = 0,
                                    return.pair.details = FALSE) {
  coords <- grip.validate.coords(coords)
  tie.mode <- match.arg(tie.mode)
  if (is.null(prepared)) {
    prepared <- grip.prepare.geodesic.mds(
      data = data,
      k = k,
      connect = connect,
      tie.mode = tie.mode
    )
  }
  prepared <- grip.validate.geodesic.mds.prepared(prepared, coords = coords)
  bend.stencils <- grip.validate.bending.stencils(
    bending.stencils,
    n = nrow(coords)
  )
  anchor.coords <- if (is.null(anchor.coords)) NULL else {
    grip.geodesic.mds.resolve.anchor(
      anchor.mode = "user",
      coords = coords,
      prepared = prepared,
      anchor.coords = anchor.coords,
      recenter = FALSE
    )
  }
  stats <- grip.geodesic.mds.score.stats(
    coords = coords,
    prepared = prepared,
    edge.length.epsilon = edge.length.epsilon,
    anchor.coords = anchor.coords,
    anchor.weight = anchor.weight,
    anchor.vertex.weight = anchor.vertex.weight,
    smoothness.weight = smoothness.weight,
    edge.spring.weight = edge.spring.weight,
    repulsion.weight = repulsion.weight,
    repulsion.quantile = repulsion.quantile,
    repulsion.scale = repulsion.scale,
    repulsion.cap.quantile = repulsion.cap.quantile,
    repulsion.hop.min = repulsion.hop.min,
    bending.stencils = bend.stencils,
    bending.weight = bending.weight
  )
  out <- data.frame(
    n = prepared$n, n.pairs = stats$n.pairs, gmds.energy = stats$energy,
    gmds.base.energy = stats$gmds.energy, gmds.raw_stress = stats$raw_stress,
    gmds.stress = stats$stress, gmds.rmse = stats$rmse,
    gmds.mean.abs.path.error = stats$mean.abs.path.error,
    gmds.mean.rel.path.error = stats$mean.rel.path.error,
    anchor.weight = stats$anchor.weight, anchor.raw.penalty = stats$anchor.raw.penalty,
    anchor.energy = stats$anchor.energy,
    edge.spring.weight = stats$edge.spring.weight,
    edge.spring.raw.penalty = stats$edge.spring.raw.penalty,
    edge.spring.energy = stats$edge.spring.energy,
    edge.spring.edge.count = stats$edge.spring.edge.count,
    repulsion.weight = stats$repulsion.weight,
    repulsion.raw.penalty = stats$repulsion.raw.penalty,
    repulsion.energy = stats$repulsion.energy,
    repulsion.pair.count = stats$repulsion.pair.count,
    repulsion.active.pair.count = stats$repulsion.active.pair.count,
    smooth.weight = stats$smooth.weight,
    smooth.raw.penalty = stats$smooth.raw.penalty, smooth.energy = stats$smooth.energy,
    bend.weight = stats$bend.weight, bend.raw.penalty = stats$bend.raw.penalty,
    bend.energy = stats$bend.energy,
    tie.mode = if (!is.null(prepared$tie_mode)) prepared$tie_mode else "single",
    stringsAsFactors = FALSE
  )
  if (isTRUE(return.pair.details)) {
    out$pair.details <- list(grip.geodesic.mds.pair.details(prepared, stats))
  }
  out
}
