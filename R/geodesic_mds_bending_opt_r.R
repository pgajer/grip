grip.optimize.geodesic.mds.bending.r <- function(coords,
                                                 prepared,
                                                 anchor.coords,
                                                 anchor.schedule,
                                                 smoothness.schedule,
                                                 edge.spring.schedule,
                                                 repulsion.schedule,
                                                 bending.stencils,
                                                 bending.schedule,
                                                 repulsion.quantile,
                                                 repulsion.scale,
                                                 repulsion.cap.quantile,
                                                 repulsion.hop.min,
                                                 edge.length.epsilon,
                                                 max.iter,
                                                 initial.step,
                                                 step.shrink,
                                                 armijo.factor,
                                                 grad.tol,
                                                 min.step,
                                                 recenter,
                                                 return.trace) {
  current <- coords
  trace.rows <- vector("list", max.iter + 1L)
  accepted.frames <- list(current)
  state <- grip.geodesic.mds.evaluate.state(
    coords = current, prepared = prepared,
    edge.length.epsilon = edge.length.epsilon,
    anchor.coords = anchor.coords,
    anchor.weight = anchor.schedule[[1L]],
    smoothness.weight = smoothness.schedule[[1L]],
    edge.spring.weight = edge.spring.schedule[[1L]],
    repulsion.weight = repulsion.schedule[[1L]],
    repulsion.quantile = repulsion.quantile,
    repulsion.scale = repulsion.scale,
    repulsion.cap.quantile = repulsion.cap.quantile,
    repulsion.hop.min = repulsion.hop.min,
    bending.stencils = bending.stencils,
    bending.weight = bending.schedule[[1L]]
  )
  trace.rows[[1L]] <- data.frame(
    iteration = 0L, energy = state$energy, gmds_energy = state$gmds_energy,
    anchor_energy = state$anchor_energy, edge_spring_energy = state$edge_spring_energy,
    repulsion_energy = state$repulsion_energy, smooth_energy = state$smooth_energy,
    bend_energy = state$bend_energy,
    gradient_norm = state$gradient_norm, step = NA_real_, accepted = TRUE,
    anchor_weight = anchor.schedule[[1L]], edge_spring_weight = edge.spring.schedule[[1L]],
    repulsion_weight = repulsion.schedule[[1L]], smooth_weight = smoothness.schedule[[1L]],
    bend_weight = bending.schedule[[1L]], stringsAsFactors = FALSE
  )
  used <- 1L

  for (iter in seq_len(max.iter)) {
    iter.anchor.weight <- anchor.schedule[[iter + 1L]]
    iter.smooth.weight <- smoothness.schedule[[iter + 1L]]
    iter.edge.spring.weight <- edge.spring.schedule[[iter + 1L]]
    iter.repulsion.weight <- repulsion.schedule[[iter + 1L]]
    iter.bend.weight <- bending.schedule[[iter + 1L]]
    state <- grip.geodesic.mds.evaluate.state(
      coords = current, prepared = prepared,
      edge.length.epsilon = edge.length.epsilon,
      anchor.coords = anchor.coords,
      anchor.weight = iter.anchor.weight,
      smoothness.weight = iter.smooth.weight,
      edge.spring.weight = iter.edge.spring.weight,
      repulsion.weight = iter.repulsion.weight,
      repulsion.quantile = repulsion.quantile,
      repulsion.scale = repulsion.scale,
      repulsion.cap.quantile = repulsion.cap.quantile,
      repulsion.hop.min = repulsion.hop.min,
      bending.stencils = bending.stencils,
      bending.weight = iter.bend.weight
    )
    if (!is.finite(state$gradient_norm) || state$gradient_norm <= grad.tol) break
    step <- as.double(initial.step)
    accepted <- FALSE
    candidate <- current
    candidate.state <- state
    while (is.finite(step) && step >= min.step) {
      proposal <- current - step * state$gradient
      if (isTRUE(recenter)) {
        proposal <- sweep(proposal, 2L, colMeans(proposal), "-", check.margin = FALSE)
      }
      proposal.state <- grip.geodesic.mds.evaluate.state(
        coords = proposal, prepared = prepared,
        edge.length.epsilon = edge.length.epsilon,
        anchor.coords = anchor.coords,
        anchor.weight = iter.anchor.weight,
        smoothness.weight = iter.smooth.weight,
        edge.spring.weight = iter.edge.spring.weight,
        repulsion.weight = iter.repulsion.weight,
        repulsion.quantile = repulsion.quantile,
        repulsion.scale = repulsion.scale,
        repulsion.cap.quantile = repulsion.cap.quantile,
        repulsion.hop.min = repulsion.hop.min,
        bending.stencils = bending.stencils,
        bending.weight = iter.bend.weight
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
    used <- used + 1L
    trace.rows[[used]] <- data.frame(
      iteration = iter,
      energy = if (accepted) candidate.state$energy else state$energy,
      gmds_energy = if (accepted) candidate.state$gmds_energy else state$gmds_energy,
      anchor_energy = if (accepted) candidate.state$anchor_energy else state$anchor_energy,
      edge_spring_energy = if (accepted) candidate.state$edge_spring_energy else state$edge_spring_energy,
      repulsion_energy = if (accepted) candidate.state$repulsion_energy else state$repulsion_energy,
      smooth_energy = if (accepted) candidate.state$smooth_energy else state$smooth_energy,
      bend_energy = if (accepted) candidate.state$bend_energy else state$bend_energy,
      gradient_norm = if (accepted) candidate.state$gradient_norm else state$gradient_norm,
      step = if (accepted) step else NA_real_, accepted = accepted,
      anchor_weight = iter.anchor.weight,
      edge_spring_weight = iter.edge.spring.weight,
      repulsion_weight = iter.repulsion.weight,
      smooth_weight = iter.smooth.weight,
      bend_weight = iter.bend.weight, stringsAsFactors = FALSE
    )
    if (!accepted) break
    current <- candidate
    state <- candidate.state
    accepted.frames[[length(accepted.frames) + 1L]] <- current
  }

  trace.df <- do.call(rbind, trace.rows[seq_len(used)])
  if (!isTRUE(return.trace)) {
    trace.df <- trace.df[, c("iteration", "energy", "gmds_energy", "anchor_energy", "edge_spring_energy", "repulsion_energy", "smooth_energy", "bend_energy", "gradient_norm", "step", "accepted", "anchor_weight", "edge_spring_weight", "repulsion_weight", "smooth_weight", "bend_weight"), drop = FALSE]
    accepted.frames <- list(current)
  }
  list(coords = current, trace = trace.df, frames = accepted.frames)
}
