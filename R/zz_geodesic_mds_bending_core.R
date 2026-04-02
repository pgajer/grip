grip.geodesic.mds.energy.gradient.base <- grip.geodesic.mds.energy.gradient
grip.geodesic.mds.score.stats.base <- grip.geodesic.mds.score.stats
grip.geodesic.mds.evaluate.state.base <- grip.geodesic.mds.evaluate.state
grip.score.geodesic.mds.base <- grip.score.geodesic.mds
grip.optimize.geodesic.mds.base <- grip.optimize.geodesic.mds

grip.geodesic.mds.energy.gradient <- function(coords,
                                              prepared,
                                              edge_length_epsilon = 1e-8,
                                              anchor_coords = NULL,
                                              anchor_weight = 0,
                                              anchor_vertex_weight = NULL,
                                              smoothness_weight = 0,
                                              edge_spring_weight = 0,
                                              repulsion_weight = 0,
                                              repulsion_quantile = 0.60,
                                              repulsion_scale = 0.20,
                                              repulsion_cap_quantile = 0.90,
                                              repulsion_hop_min = 3L,
                                              bending_stencils = NULL,
                                              bending_weight = 0) {
  base <- grip.geodesic.mds.energy.gradient.base(
    coords = coords,
    prepared = prepared,
    edge_length_epsilon = edge_length_epsilon,
    anchor_coords = anchor_coords,
    anchor_weight = anchor_weight,
    anchor_vertex_weight = anchor_vertex_weight,
    smoothness_weight = smoothness_weight,
    edge_spring_weight = edge_spring_weight,
    repulsion_weight = repulsion_weight,
    repulsion_quantile = repulsion_quantile,
    repulsion_scale = repulsion_scale,
    repulsion_cap_quantile = repulsion_cap_quantile,
    repulsion_hop_min = repulsion_hop_min
  )
  bend.stats <- grip.geodesic.mds.bending.stats(
    coords = coords,
    bending_stencils = bending_stencils,
    bending_weight = bending_weight
  )
  total.grad <- base$gradient + bend.stats$gradient
  base$energy <- base$energy + bend.stats$energy
  base$bend_energy <- bend.stats$energy
  base$bend_raw_penalty <- bend.stats$raw_penalty
  base$gradient <- total.grad
  base$gradient_norm <- sqrt(sum(total.grad^2))
  base
}

grip.geodesic.mds.score.stats <- function(coords,
                                          prepared,
                                          edge_length_epsilon = 1e-8,
                                          anchor_coords = NULL,
                                          anchor_weight = 0,
                                          anchor_vertex_weight = NULL,
                                          smoothness_weight = 0,
                                          edge_spring_weight = 0,
                                          repulsion_weight = 0,
                                          repulsion_quantile = 0.60,
                                          repulsion_scale = 0.20,
                                          repulsion_cap_quantile = 0.90,
                                          repulsion_hop_min = 3L,
                                          bending_stencils = NULL,
                                          bending_weight = 0) {
  base <- grip.geodesic.mds.score.stats.base(
    coords = coords,
    prepared = prepared,
    edge_length_epsilon = edge_length_epsilon,
    anchor_coords = anchor_coords,
    anchor_weight = anchor_weight,
    anchor_vertex_weight = anchor_vertex_weight,
    smoothness_weight = smoothness_weight,
    edge_spring_weight = edge_spring_weight,
    repulsion_weight = repulsion_weight,
    repulsion_quantile = repulsion_quantile,
    repulsion_scale = repulsion_scale,
    repulsion_cap_quantile = repulsion_cap_quantile,
    repulsion_hop_min = repulsion_hop_min
  )
  bend.stats <- grip.geodesic.mds.bending.stats(
    coords = coords,
    bending_stencils = bending_stencils,
    bending_weight = bending_weight
  )
  base$energy <- base$energy + bend.stats$energy
  base$bend.weight <- bend.stats$bending_weight
  base$bend.raw.penalty <- bend.stats$raw_penalty
  base$bend.energy <- bend.stats$energy
  base
}

grip.geodesic.mds.evaluate.state <- function(coords,
                                             prepared,
                                             edge_length_epsilon = 1e-8,
                                             anchor_coords = NULL,
                                             anchor_weight = 0,
                                             anchor_vertex_weight = NULL,
                                             smoothness_weight = 0,
                                             edge_spring_weight = 0,
                                             repulsion_weight = 0,
                                             repulsion_quantile = 0.60,
                                             repulsion_scale = 0.20,
                                             repulsion_cap_quantile = 0.90,
                                             repulsion_hop_min = 3L,
                                             bending_stencils = NULL,
                                             bending_weight = 0) {
  grip.geodesic.mds.energy.gradient(
    coords = coords,
    prepared = prepared,
    edge_length_epsilon = edge_length_epsilon,
    anchor_coords = anchor_coords,
    anchor_weight = anchor_weight,
    anchor_vertex_weight = anchor_vertex_weight,
    smoothness_weight = smoothness_weight,
    edge_spring_weight = edge_spring_weight,
    repulsion_weight = repulsion_weight,
    repulsion_quantile = repulsion_quantile,
    repulsion_scale = repulsion_scale,
    repulsion_cap_quantile = repulsion_cap_quantile,
    repulsion_hop_min = repulsion_hop_min,
    bending_stencils = bending_stencils,
    bending_weight = bending_weight
  )
}
