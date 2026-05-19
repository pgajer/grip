# Final edge-KK polish reference for grip.layout.misf.edge.kk().
#
# This development-only layer validates the final auxiliary-free polish stage.
# It uses only original edge stress and the same Armijo reference machinery.

mek_profile_edge_scale <- function(coords,
                                   edges,
                                   edge_lengths,
                                   stiffness = rep(1, nrow(edges)),
                                   eps = 1e-8) {
  coords <- mek_validate_coords(coords)
  edges <- mek_validate_pair_matrix(edges, nrow(coords), "edges")
  edge_lengths <- mek_validate_positive_vector(
    edge_lengths, nrow(edges), "edge_lengths"
  )
  stiffness <- mek_validate_positive_vector(stiffness, nrow(edges), "stiffness")
  lengths <- vapply(seq_len(nrow(edges)), function(e) {
    u <- edges[e, 1L]
    v <- edges[e, 2L]
    sqrt(sum((coords[u, ] - coords[v, ])^2) + eps^2)
  }, numeric(1L))
  denom <- sum(stiffness * edge_lengths^2)
  if (!is.finite(denom) || denom <= 0) {
    return(1)
  }
  scale <- sum(stiffness * edge_lengths * lengths) / denom
  if (is.finite(scale) && scale > 0) {
    scale
  } else {
    1
  }
}

mek_polish_edge_scale <- function(coords,
                                  edges,
                                  edge_lengths,
                                  edge_stiffness,
                                  scale_mode,
                                  fixed_scale) {
  if (identical(scale_mode, "profiled")) {
    mek_profile_edge_scale(
      coords = coords,
      edges = edges,
      edge_lengths = edge_lengths,
      stiffness = edge_stiffness
    )
  } else if (identical(scale_mode, "identity")) {
    1
  } else {
    fixed_scale
  }
}

mek_final_edge_kk_polish_reference <- function(coords,
                                               edges,
                                               edge_lengths,
                                               max_iter = 25L,
                                               initial_step = 0.25,
                                               scale_mode = c(
                                                 "profiled",
                                                 "identity",
                                                 "global_fixed",
                                                 "level_fixed",
                                                 "user"
                                               ),
                                               fixed_scale = 1) {
  coords <- mek_validate_coords(coords)
  edges <- mek_validate_pair_matrix(edges, nrow(coords), "edges")
  edge_lengths <- mek_validate_positive_vector(
    edge_lengths, nrow(edges), "edge_lengths"
  )
  scale_mode <- match.arg(scale_mode)
  if (!is.finite(fixed_scale) || fixed_scale <= 0) {
    stop("fixed_scale must be finite and > 0", call. = FALSE)
  }
  mek_validate_misf_scale_policy(scale_mode, rho = 0, lambda = 0)

  active <- seq_len(nrow(coords))
  edge_stiffness <- rep(1, nrow(edges))
  empty_pairs <- matrix(integer(0L), ncol = 2L)
  anchor_targets <- coords
  anchor_weights <- numeric(nrow(coords))

  state_fn <- function(x) {
    edge_scale <- mek_polish_edge_scale(
      coords = x,
      edges = edges,
      edge_lengths = edge_lengths,
      edge_stiffness = edge_stiffness,
      scale_mode = scale_mode,
      fixed_scale = fixed_scale
    )
    mek_total_state(
      coords = x,
      edges = edges,
      edge_lengths = edge_lengths,
      edge_stiffness = edge_stiffness,
      metric_pairs = empty_pairs,
      metric_distances = numeric(0L),
      metric_stiffness = numeric(0L),
      repulsion_pairs = empty_pairs,
      repulsion_weights = numeric(0L),
      anchor_targets = anchor_targets,
      anchor_weights = anchor_weights,
      rho = 0,
      lambda = 0,
      edge_scale = edge_scale,
      metric_scale = 1,
      edge_eps = 1e-8,
      repulsion_eps = 1e-4
    )
  }

  initial_state <- state_fn(coords)
  start_time <- proc.time()[["elapsed"]]
  fit <- mek_armijo_refine(
    coords = coords,
    active = active,
    state_fn = state_fn,
    max_iter = max_iter,
    initial_step = initial_step,
    recenter = TRUE
  )
  elapsed <- proc.time()[["elapsed"]] - start_time
  final_state <- fit$final_state
  last_accepted <- which(fit$trace$accepted)
  accepted_step <- if (length(last_accepted)) {
    fit$trace$step[[tail(last_accepted, 1L)]]
  } else {
    NA_real_
  }

  polish_trace <- data.frame(
    stage = "final_polish",
    active_count = length(active),
    edge_count = nrow(edges),
    metric_count = 0L,
    repulsion_count = 0L,
    rho = 0,
    lambda = 0,
    anchor_weight_mean = 0,
    anchor_weight_max = 0,
    scale_mode = scale_mode,
    scale_value = final_state$edge$scale,
    initial_total_energy = initial_state$energy,
    total_energy = final_state$energy,
    edge_energy = final_state$edge$energy,
    metric_energy = final_state$metric$energy,
    repulsion_energy = final_state$repulsion$energy,
    anchor_energy = final_state$anchor$energy,
    gradient_norm = mek_active_gradient_norm(final_state$gradient, active),
    accepted_step = accepted_step,
    edge_rmse = final_state$edge$rel_rmse,
    metric_rmse = final_state$metric$rel_rmse,
    armijo_iterations = nrow(fit$trace),
    armijo_accepted = sum(fit$trace$accepted),
    stop_reason = fit$stop_reason,
    elapsed_time = elapsed,
    stringsAsFactors = FALSE
  )

  list(
    coords = fit$coords,
    polish_trace = polish_trace,
    armijo_trace = fit$trace,
    active = active,
    edges = edges,
    metric_pairs = empty_pairs,
    repulsion_pairs = empty_pairs,
    anchor_weights = anchor_weights,
    initial_state = initial_state,
    final_state = final_state
  )
}
