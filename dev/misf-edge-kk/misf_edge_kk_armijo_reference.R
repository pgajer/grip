# Reference Armijo refinement for grip.layout.misf.edge.kk().
#
# Development-only implementation used to validate active-level optimizer
# contracts before adding any package API or compiled backend.

mek_recenter_active <- function(coords, active) {
  coords <- mek_validate_coords(coords)
  active <- as.integer(active)
  if (any(active < 1L | active > nrow(coords))) {
    stop("active contains vertices outside coords", call. = FALSE)
  }
  center <- colMeans(coords[active, , drop = FALSE])
  coords[active, ] <- sweep(coords[active, , drop = FALSE], 2L, center, "-")
  coords
}

mek_active_gradient_norm <- function(gradient, active) {
  sqrt(sum(gradient[active, , drop = FALSE]^2))
}

mek_armijo_refine <- function(coords,
                              active,
                              state_fn,
                              max_iter = 25L,
                              initial_step = 1,
                              step_shrink = 0.5,
                              armijo_factor = 1e-4,
                              grad_tol = 1e-8,
                              min_step = 1e-8,
                              recenter = TRUE) {
  coords <- mek_validate_coords(coords)
  active <- as.integer(active)
  if (any(active < 1L | active > nrow(coords))) {
    stop("active contains vertices outside coords", call. = FALSE)
  }
  if (!is.function(state_fn)) {
    stop("state_fn must be a function", call. = FALSE)
  }
  if (max_iter < 0L) {
    stop("max_iter must be non-negative", call. = FALSE)
  }
  if (!is.finite(initial_step) || initial_step <= 0) {
    stop("initial_step must be finite and > 0", call. = FALSE)
  }
  if (!is.finite(step_shrink) || step_shrink <= 0 || step_shrink >= 1) {
    stop("step_shrink must be finite and in (0, 1)", call. = FALSE)
  }
  if (!is.finite(armijo_factor) || armijo_factor < 0) {
    stop("armijo_factor must be finite and >= 0", call. = FALSE)
  }
  if (!is.finite(grad_tol) || grad_tol < 0) {
    stop("grad_tol must be finite and >= 0", call. = FALSE)
  }
  if (!is.finite(min_step) || min_step <= 0) {
    stop("min_step must be finite and > 0", call. = FALSE)
  }

  trace <- list()
  stop_reason <- "iteration_budget"
  current <- coords

  for (iter in seq_len(max_iter)) {
    state <- state_fn(current)
    grad_norm <- mek_active_gradient_norm(state$gradient, active)
    if (!is.finite(state$energy) || !is.finite(grad_norm)) {
      stop_reason <- "nonfinite_state"
      break
    }
    if (grad_norm <= grad_tol) {
      stop_reason <- "gradient_tolerance"
      trace[[length(trace) + 1L]] <- data.frame(
        iteration = iter,
        energy_before = state$energy,
        energy_trial = state$energy,
        energy_after_recenter = state$energy,
        gradient_norm = grad_norm,
        step = NA_real_,
        shrink_attempts = 0L,
        accepted = FALSE,
        armijo_target = state$energy,
        armijo_satisfied = TRUE,
        stopping_reason = stop_reason,
        stringsAsFactors = FALSE
      )
      break
    }

    step <- initial_step
    shrink_attempts <- 0L
    accepted <- FALSE
    trial_energy <- NA_real_
    trial_state <- NULL
    target <- NA_real_
    trial <- current

    while (is.finite(step) && step >= min_step) {
      trial <- current
      trial[active, ] <- trial[active, , drop = FALSE] -
        step * state$gradient[active, , drop = FALSE]
      trial_state <- state_fn(trial)
      trial_energy <- trial_state$energy
      target <- state$energy - armijo_factor * step * grad_norm^2
      if (is.finite(trial_energy) && trial_energy <= target) {
        accepted <- TRUE
        break
      }
      step <- step * step_shrink
      shrink_attempts <- shrink_attempts + 1L
    }

    if (!accepted) {
      stop_reason <- "no_accepted_step"
      trace[[length(trace) + 1L]] <- data.frame(
        iteration = iter,
        energy_before = state$energy,
        energy_trial = trial_energy,
        energy_after_recenter = NA_real_,
        gradient_norm = grad_norm,
        step = step,
        shrink_attempts = shrink_attempts,
        accepted = FALSE,
        armijo_target = target,
        armijo_satisfied = FALSE,
        stopping_reason = stop_reason,
        stringsAsFactors = FALSE
      )
      break
    }

    current <- trial
    if (isTRUE(recenter)) {
      current <- mek_recenter_active(current, active)
    }
    after_state <- state_fn(current)
    trace[[length(trace) + 1L]] <- data.frame(
      iteration = iter,
      energy_before = state$energy,
      energy_trial = trial_energy,
      energy_after_recenter = after_state$energy,
      gradient_norm = grad_norm,
      step = step,
      shrink_attempts = shrink_attempts,
      accepted = TRUE,
      armijo_target = target,
      armijo_satisfied = trial_energy <= target,
      stopping_reason = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  if (max_iter == 0L) {
    stop_reason <- "iteration_budget"
  }
  trace_df <- if (length(trace)) {
    do.call(rbind, trace)
  } else {
    data.frame(
      iteration = integer(),
      energy_before = numeric(),
      energy_trial = numeric(),
      energy_after_recenter = numeric(),
      gradient_norm = numeric(),
      step = numeric(),
      shrink_attempts = integer(),
      accepted = logical(),
      armijo_target = numeric(),
      armijo_satisfied = logical(),
      stopping_reason = character(),
      stringsAsFactors = FALSE
    )
  }
  if (nrow(trace_df) && all(is.na(trace_df$stopping_reason))) {
    trace_df$stopping_reason[[nrow(trace_df)]] <- stop_reason
  }
  list(
    coords = current,
    trace = trace_df,
    stop_reason = stop_reason,
    final_state = state_fn(current)
  )
}

mek_armijo_fixture <- function(dim = 3L,
                               repeated = FALSE,
                               include_regularization = TRUE) {
  fixture <- mek_objective_fixture(dim = dim, repeated = repeated)
  active <- seq_len(nrow(fixture$coords))
  rho <- if (include_regularization) 0.25 else 0
  lambda <- if (include_regularization) 0.1 else 0
  anchor_weights <- if (include_regularization) {
    fixture$anchor_weights
  } else {
    rep(0, nrow(fixture$coords))
  }
  state_fn <- function(x) {
    mek_total_state(
      coords = x,
      edges = fixture$edges,
      edge_lengths = fixture$edge_lengths,
      edge_stiffness = fixture$edge_stiffness,
      metric_pairs = fixture$metric_pairs,
      metric_distances = fixture$metric_distances,
      metric_stiffness = fixture$metric_stiffness,
      repulsion_pairs = fixture$repulsion_pairs,
      repulsion_weights = fixture$repulsion_weights,
      anchor_targets = fixture$anchor_targets,
      anchor_weights = anchor_weights,
      rho = rho,
      lambda = lambda,
      edge_scale = 1,
      metric_scale = 1,
      edge_eps = 1e-8,
      repulsion_eps = 1e-4
    )
  }
  list(
    coords = fixture$coords,
    active = active,
    state_fn = state_fn,
    repulsion_pairs = fixture$repulsion_pairs,
    repulsion_weights = fixture$repulsion_weights,
    rho = rho,
    lambda = lambda
  )
}
