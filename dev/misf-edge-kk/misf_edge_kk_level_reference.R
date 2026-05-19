# Single-level reference refinement for grip.layout.misf.edge.kk().
#
# This development-only layer integrates the already validated objective,
# constraint, scale-policy, and Armijo pieces for one active MISF level.

mek_metric_pairs_matrix <- function(metric_constraints) {
  if (nrow(metric_constraints) == 0L) {
    return(matrix(integer(0L), ncol = 2L))
  }
  pairs <- cbind(metric_constraints$i, metric_constraints$j)
  storage.mode(pairs) <- "integer"
  pairs
}

mek_full_anchor_targets <- function(coords, anchors) {
  targets <- coords
  targets[anchors$vertices, ] <- anchors$coords
  targets
}

mek_full_anchor_weights <- function(n, anchors) {
  weights <- numeric(n)
  weights[anchors$vertices] <- anchors$weights
  weights
}

mek_repulsion_weights_by_mode <- function(pairs,
                                          active_count,
                                          mode = c("active_count",
                                                   "unit",
                                                   "sqrt_active")) {
  mode <- match.arg(mode)
  if (nrow(pairs) == 0L) {
    return(numeric(0L))
  }
  total_mass <- switch(
    mode,
    active_count = active_count,
    unit = 1,
    sqrt_active = sqrt(active_count)
  )
  rep(total_mass / nrow(pairs), nrow(pairs))
}

mek_refine_level_reference <- function(n,
                                       edges,
                                       edge_lengths,
                                       coords_post_insertion,
                                       misf_order,
                                       level_sizes,
                                       level_index,
                                       candidate_metric_pairs,
                                       rho = 0.5,
                                       lambda = 0.1,
                                       anchor_weight = 0.25,
                                       exact_repulsion_below = 64L,
                                       repulsion_sample_count = 128L,
                                       repulsion_seed = 1L,
                                       max_iter = 10L,
                                       initial_step = 0.5,
                                       edge_scale = 1,
                                       metric_scale = 1,
                                       scale_mode = "identity",
                                       repulsion_weight_mode = c("active_count",
                                                                 "unit",
                                                                 "sqrt_active")) {
  coords_post_insertion <- mek_validate_coords(coords_post_insertion)
  if (nrow(coords_post_insertion) != n) {
    stop("coords_post_insertion row count must equal n", call. = FALSE)
  }
  repulsion_weight_mode <- match.arg(repulsion_weight_mode)
  mek_validate_misf_scale_policy(scale_mode, rho = rho, lambda = lambda)

  active <- mek_misf_active_set(misf_order, level_sizes, level_index)
  inserted <- mek_misf_new_vertices(misf_order, level_sizes, level_index)
  active_edges <- mek_active_edges(edges, active)
  active_edge_keys <- mek_pair_key(edges)
  active_edge_index <- match(mek_pair_key(active_edges), active_edge_keys)
  active_edge_lengths <- edge_lengths[active_edge_index]
  active_edge_stiffness <- rep(1, length(active_edge_lengths))

  metric <- mek_metric_constraints(
    n = n,
    edges = edges,
    edge_lengths = edge_lengths,
    active = active,
    candidate_pairs = candidate_metric_pairs
  )
  metric_pairs <- mek_metric_pairs_matrix(metric)
  metric_stiffness <- rep(1, nrow(metric_pairs))

  repulsion_pairs <- mek_repulsion_pairs(
    active,
    exact_below = exact_repulsion_below,
    sample_count = repulsion_sample_count,
    seed = repulsion_seed,
    level_index = level_index
  )
  repulsion_weights <- mek_repulsion_weights_by_mode(
    repulsion_pairs,
    active_count = length(active),
    mode = repulsion_weight_mode
  )
  repulsion_key <- mek_pair_key(repulsion_pairs)

  anchors <- mek_capture_anchors(
    coords = coords_post_insertion,
    active = active,
    weights = rep(anchor_weight, length(active))
  )
  anchor_targets <- mek_full_anchor_targets(coords_post_insertion, anchors)
  anchor_weights <- mek_full_anchor_weights(n, anchors)

  state_calls <- 0L
  state_fn <- function(coords) {
    state_calls <<- state_calls + 1L
    if (!identical(repulsion_key, mek_pair_key(repulsion_pairs))) {
      stop("repulsion pairs changed during level refinement", call. = FALSE)
    }
    mek_total_state(
      coords = coords,
      edges = active_edges,
      edge_lengths = active_edge_lengths,
      edge_stiffness = active_edge_stiffness,
      metric_pairs = metric_pairs,
      metric_distances = metric$distance,
      metric_stiffness = metric_stiffness,
      repulsion_pairs = repulsion_pairs,
      repulsion_weights = repulsion_weights,
      anchor_targets = anchor_targets,
      anchor_weights = anchor_weights,
      rho = rho,
      lambda = lambda,
      edge_scale = edge_scale,
      metric_scale = metric_scale,
      edge_eps = 1e-8,
      repulsion_eps = 1e-4
    )
  }

  start_time <- proc.time()[["elapsed"]]
  fit <- mek_armijo_refine(
    coords = coords_post_insertion,
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

  level_trace <- data.frame(
    level = level_index,
    active_count = length(active),
    inserted_count = length(inserted),
    edge_count = nrow(active_edges),
    metric_count = nrow(metric_pairs),
    repulsion_count = nrow(repulsion_pairs),
    repulsion_weight_mode = repulsion_weight_mode,
    repulsion_weight_mass = sum(repulsion_weights),
    rho = rho,
    lambda = lambda,
    anchor_weight_mean = mean(anchor_weights[active]),
    anchor_weight_max = max(anchor_weights[active]),
    scale_mode = scale_mode,
    edge_scale = edge_scale,
    metric_scale = metric_scale,
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
    state_calls = state_calls,
    stop_reason = fit$stop_reason,
    elapsed_time = elapsed,
    stringsAsFactors = FALSE
  )

  list(
    coords = fit$coords,
    level_trace = level_trace,
    armijo_trace = fit$trace,
    active = active,
    inserted = inserted,
    active_edges = active_edges,
    metric_constraints = metric,
    repulsion_pairs = repulsion_pairs,
    repulsion_weights = repulsion_weights,
    anchors = anchors,
    final_state = final_state
  )
}

mek_level_fixture <- function(dim = 2L) {
  fixture <- mek_contract_fixture()
  coords <- fixture$coords_post_insertion
  if (dim > 2L) {
    extra <- matrix(seq(0.05, 0.25, length.out = fixture$n * (dim - 2L)),
                    nrow = fixture$n)
    coords <- cbind(coords, extra)
  }
  fixture$coords_post_insertion <- coords
  fixture
}
