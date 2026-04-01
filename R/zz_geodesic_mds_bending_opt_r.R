grip.optimize.geodesic.mds.bending.r <- function(coords,
                                                 prepared,
                                                 anchor.coords,
                                                 anchor.schedule,
                                                 smoothness.schedule,
                                                 edge.spring.schedule,
                                                 repulsion.schedule,
                                                 bending.stencils,
                                                 bending.schedule,
                                                 repulsion_quantile,
                                                 repulsion_scale,
                                                 repulsion_cap_quantile,
                                                 repulsion_hop_min,
                                                 edge_length_epsilon,
                                                 max_iter,
                                                 initial_step,
                                                 step_shrink,
                                                 armijo_factor,
                                                 grad_tol,
                                                 min_step,
                                                 recenter,
                                                 return_trace) {
  current <- coords
  trace.rows <- vector("list", max_iter + 1L)
  accepted.frames <- list(current)
  state <- grip.geodesic.mds.evaluate.state(
    coords = current, prepared = prepared,
    edge_length_epsilon = edge_length_epsilon,
    anchor_coords = anchor.coords,
    anchor_weight = anchor.schedule[[1L]],
    smoothness_weight = smoothness.schedule[[1L]],
    edge_spring_weight = edge.spring.schedule[[1L]],
    repulsion_weight = repulsion.schedule[[1L]],
    repulsion_quantile = repulsion_quantile,
    repulsion_scale = repulsion_scale,
    repulsion_cap_quantile = repulsion_cap_quantile,
    repulsion_hop_min = repulsion_hop_min,
    bending_stencils = bending.stencils,
    bending_weight = bending.schedule[[1L]]
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

  for (iter in seq_len(max_iter)) {
    iter.anchor.weight <- anchor.schedule[[iter + 1L]]
    iter.smooth.weight <- smoothness.schedule[[iter + 1L]]
    iter.edge.spring.weight <- edge.spring.schedule[[iter + 1L]]
    iter.repulsion.weight <- repulsion.schedule[[iter + 1L]]
    iter.bend.weight <- bending.schedule[[iter + 1L]]
    state <- grip.geodesic.mds.evaluate.state(
      coords = current, prepared = prepared,
      edge_length_epsilon = edge_length_epsilon,
      anchor_coords = anchor.coords,
      anchor_weight = iter.anchor.weight,
      smoothness_weight = iter.smooth.weight,
      edge_spring_weight = iter.edge.spring.weight,
      repulsion_weight = iter.repulsion.weight,
      repulsion_quantile = repulsion_quantile,
      repulsion_scale = repulsion_scale,
      repulsion_cap_quantile = repulsion_cap_quantile,
      repulsion_hop_min = repulsion_hop_min,
      bending_stencils = bending.stencils,
      bending_weight = iter.bend.weight
    )
    if (!is.finite(state$gradient_norm) || state$gradient_norm <= grad_tol) break
    step <- as.double(initial_step)
    accepted <- FALSE
    candidate <- current
    candidate.state <- state
    while (is.finite(step) && step >= min_step) {
      proposal <- current - step * state$gradient
      if (isTRUE(recenter)) {
        proposal <- sweep(proposal, 2L, colMeans(proposal), "-", check.margin = FALSE)
      }
      proposal.state <- grip.geodesic.mds.evaluate.state(
        coords = proposal, prepared = prepared,
        edge_length_epsilon = edge_length_epsilon,
        anchor_coords = anchor.coords,
        anchor_weight = iter.anchor.weight,
        smoothness_weight = iter.smooth.weight,
        edge_spring_weight = iter.edge.spring.weight,
        repulsion_weight = iter.repulsion.weight,
        repulsion_quantile = repulsion_quantile,
        repulsion_scale = repulsion_scale,
        repulsion_cap_quantile = repulsion_cap_quantile,
        repulsion_hop_min = repulsion_hop_min,
        bending_stencils = bending.stencils,
        bending_weight = iter.bend.weight
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
  if (!isTRUE(return_trace)) {
    trace.df <- trace.df[, c("iteration", "energy", "gmds_energy", "anchor_energy", "edge_spring_energy", "repulsion_energy", "smooth_energy", "bend_energy", "gradient_norm", "step", "accepted", "anchor_weight", "edge_spring_weight", "repulsion_weight", "smooth_weight", "bend_weight"), drop = FALSE]
    accepted.frames <- list(current)
  }
  list(coords = current, trace = trace.df, frames = accepted.frames)
}
