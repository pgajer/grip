grip.geodesic.mds.energy.gradient <- function(coords,
                                              prepared,
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
                                              bending.weight = 0) {
  base <- grip.geodesic.mds.energy.gradient.base(
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
    repulsion.hop.min = repulsion.hop.min
  )
  bend.stats <- grip.geodesic.mds.bending.stats(
    coords = coords,
    bending.stencils = bending.stencils,
    bending.weight = bending.weight
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
                                          bending.weight = 0) {
  base <- grip.geodesic.mds.score.stats.base(
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
    repulsion.hop.min = repulsion.hop.min
  )
  bend.stats <- grip.geodesic.mds.bending.stats(
    coords = coords,
    bending.stencils = bending.stencils,
    bending.weight = bending.weight
  )
  base$energy <- base$energy + bend.stats$energy
  base$bend.weight <- bend.stats$bending_weight
  base$bend.raw.penalty <- bend.stats$raw_penalty
  base$bend.energy <- bend.stats$energy
  base
}

grip.geodesic.mds.evaluate.state <- function(coords,
                                             prepared,
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
                                             bending.weight = 0) {
  grip.geodesic.mds.energy.gradient(
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
    bending.stencils = bending.stencils,
    bending.weight = bending.weight
  )
}
