grip.score.geodesic.mds <- function(coords,
                                    prepared = NULL,
                                    data = NULL,
                                    k = NULL,
                                    connect = c("mst", "error"),
                                    tie_mode = c("single", "average"),
                                    edge_length_epsilon = 1e-8,
                                    anchor_coords = NULL,
                                    anchor_weight = 0,
                                    smoothness_weight = 0,
                                    edge_spring_weight = 0,
                                    repulsion_weight = 0,
                                    repulsion_quantile = 0.60,
                                    repulsion_scale = 0.20,
                                    repulsion_cap_quantile = 0.90,
                                    repulsion_hop_min = 3L,
                                    bending_stencils = NULL,
                                    bending_weight = 0,
                                    return_pair_details = FALSE) {
  coords <- grip.validate.coords(coords)
  tie_mode <- match.arg(tie_mode)
  if (is.null(prepared)) {
    prepared <- grip.prepare.geodesic.mds(
      data = data,
      k = k,
      connect = connect,
      tie_mode = tie_mode
    )
  }
  prepared <- grip.validate.geodesic.mds.prepared(prepared, coords = coords)
  bend.stencils <- grip.validate.bending.stencils(
    bending_stencils,
    n = nrow(coords)
  )
  anchor.coords <- if (is.null(anchor_coords)) NULL else {
    grip.geodesic.mds.resolve.anchor(
      anchor_mode = "user",
      coords = coords,
      prepared = prepared,
      anchor_coords = anchor_coords,
      recenter = FALSE
    )
  }
  stats <- grip.geodesic.mds.score.stats(
    coords = coords,
    prepared = prepared,
    edge_length_epsilon = edge_length_epsilon,
    anchor_coords = anchor.coords,
    anchor_weight = anchor_weight,
    smoothness_weight = smoothness_weight,
    edge_spring_weight = edge_spring_weight,
    repulsion_weight = repulsion_weight,
    repulsion_quantile = repulsion_quantile,
    repulsion_scale = repulsion_scale,
    repulsion_cap_quantile = repulsion_cap_quantile,
    repulsion_hop_min = repulsion_hop_min,
    bending_stencils = bend.stencils,
    bending_weight = bending_weight
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
  if (isTRUE(return_pair_details)) {
    out$pair.details <- list(grip.geodesic.mds.pair.details(prepared, stats))
  }
  out
}
