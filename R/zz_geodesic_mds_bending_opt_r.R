grip.optimize.geodesic.mds.bending.r <- function(coords,
                                                 prepared,
                                                 anchor.coords,
                                                 anchor.schedule,
                                                 bending.stencils,
                                                 bending.schedule,
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
    smoothness_weight = 0,
    bending_stencils = bending.stencils,
    bending_weight = bending.schedule[[1L]]
  )
  trace.rows[[1L]] <- data.frame(
    iteration = 0L, energy = state$energy, gmds_energy = state$gmds_energy,
    anchor_energy = state$anchor_energy, smooth_energy = 0, bend_energy = state$bend_energy,
    gradient_norm = state$gradient_norm, step = NA_real_, accepted = TRUE,
    anchor_weight = anchor.schedule[[1L]], smooth_weight = 0,
    bend_weight = bending.schedule[[1L]], stringsAsFactors = FALSE
  )
  used <- 1L

  for (iter in seq_len(max_iter)) {
    iter.anchor.weight <- anchor.schedule[[iter + 1L]]
    iter.bend.weight <- bending.schedule[[iter + 1L]]
    state <- grip.geodesic.mds.evaluate.state(
      coords = current, prepared = prepared,
      edge_length_epsilon = edge_length_epsilon,
      anchor_coords = anchor.coords,
      anchor_weight = iter.anchor.weight,
      smoothness_weight = 0,
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
        smoothness_weight = 0,
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
      smooth_energy = 0,
      bend_energy = if (accepted) candidate.state$bend_energy else state$bend_energy,
      gradient_norm = if (accepted) candidate.state$gradient_norm else state$gradient_norm,
      step = if (accepted) step else NA_real_, accepted = accepted,
      anchor_weight = iter.anchor.weight, smooth_weight = 0,
      bend_weight = iter.bend.weight, stringsAsFactors = FALSE
    )
    if (!accepted) break
    current <- candidate
    state <- candidate.state
    accepted.frames[[length(accepted.frames) + 1L]] <- current
  }

  trace.df <- do.call(rbind, trace.rows[seq_len(used)])
  if (!isTRUE(return_trace)) {
    trace.df <- trace.df[, c("iteration", "energy", "gmds_energy", "anchor_energy", "smooth_energy", "bend_energy", "gradient_norm", "step", "accepted", "anchor_weight", "smooth_weight", "bend_weight"), drop = FALSE]
    accepted.frames <- list(current)
  }
  list(coords = current, trace = trace.df, frames = accepted.frames)
}
