.grip.edge.repulsive.center.coords <- function(Z) {
  sweep(as.matrix(Z), 2L, colMeans(Z), "-", check.margin = FALSE)
}

.grip.edge.repulsive.potential.state <- function(s,
                                                 family = c("quadratic", "upper_barrier"),
                                                 eps.plus = 0.35,
                                                 beta = 0,
                                                 wall.margin = 1e-10) {
  family <- match.arg(family)
  s <- as.numeric(s)
  if (identical(family, "quadratic") || beta <= 0) {
    return(list(value = 0.5 * (s - 1)^2, deriv = s - 1, feasible = rep(TRUE, length(s))))
  }
  u <- (s - 1) / eps.plus
  feasible <- u < 1 - wall.margin
  value <- rep(Inf, length(s))
  deriv <- rep(NA_real_, length(s))
  ok <- feasible & is.finite(u)
  value[ok] <- 0.5 * (s[ok] - 1)^2 + beta * (-log1p(-u[ok]) - u[ok])
  deriv[ok] <- (s[ok] - 1) + beta * (1 / (1 - u[ok]) - 1) / eps.plus
  list(value = value, deriv = deriv, feasible = feasible)
}

.grip.edge.repulsive.repulsion.state <- function(r,
                                                 family = c("log", "inverse_power"),
                                                 delta = 1e-3,
                                                 power = 1) {
  family <- match.arg(family)
  r <- as.numeric(r)
  rr <- r^2 + delta^2
  if (identical(family, "log")) {
    value <- -0.5 * log(rr)
    deriv.r <- -r / rr
  } else {
    value <- rr^(-power / 2)
    deriv.r <- -power * r * rr^(-power / 2 - 1)
  }
  list(value = value, deriv.r = deriv.r)
}

.grip.repulsive.state.R <- function(coords,
                                    lambda = 1,
                                    pair.index = NULL,
                                    pair.weights = NULL,
                                    repulsion.family = "log",
                                    repulsion.delta = 1e-3,
                                    repulsion.power = 1,
                                    distance.eps = 1e-10) {
  Z <- as.matrix(coords)
  n <- nrow(Z)
  dim <- ncol(Z)
  pair.index <- .grip.edge.repulsive.normalize.pairs(pair.index, n, lambda)
  if (is.null(pair.weights)) pair.weights <- rep(1, nrow(pair.index))
  if (length(pair.weights) != nrow(pair.index)) {
    stop("pair.weights must be parallel to pair.index", call. = FALSE)
  }
  pdiff <- if (nrow(pair.index)) {
    Z[pair.index[, 1L], , drop = FALSE] - Z[pair.index[, 2L], , drop = FALSE]
  } else {
    matrix(numeric(), nrow = 0L, ncol = dim)
  }
  pr <- if (nrow(pair.index)) sqrt(rowSums(pdiff^2) + distance.eps^2) else numeric()
  repel.energy <- 0
  grad <- matrix(0, n, dim)
  if (lambda > 0 && nrow(pair.index)) {
    rp <- .grip.edge.repulsive.repulsion.state(
      pr,
      family = repulsion.family,
      delta = repulsion.delta,
      power = repulsion.power
    )
    repel.energy <- sum(pair.weights * rp$value)
    coeff <- lambda * pair.weights * rp$deriv.r / pr
    for (p in seq_len(nrow(pair.index))) {
      g <- coeff[[p]] * pdiff[p, ]
      grad[pair.index[p, 1L], ] <- grad[pair.index[p, 1L], ] + g
      grad[pair.index[p, 2L], ] <- grad[pair.index[p, 2L], ] - g
    }
  }
  list(
    energy = lambda * repel.energy,
    repel.energy = repel.energy,
    gradient = grad,
    gradient.norm = sqrt(sum(grad^2)),
    pair.embedded.lengths = pr
  )
}

.grip.edge.repulsive.normalize.edges <- function(edges) {
  edges <- as.matrix(edges)
  if (!is.numeric(edges) || ncol(edges) != 2L || nrow(edges) < 1L) {
    stop("edges must be a numeric or integer matrix with two columns and at least one row", call. = FALSE)
  }
  if (any(!is.finite(edges)) || any(edges != floor(edges)) || any(edges < 1L)) {
    stop("edges must contain finite positive 1-based vertex indices", call. = FALSE)
  }
  storage.mode(edges) <- "integer"
  edges
}

.grip.edge.repulsive.normalize.pairs <- function(pair.index, n, lambda) {
  if (is.null(pair.index)) {
    if (lambda > 0) {
      pair.index <- t(utils::combn(n, 2L))
    } else {
      pair.index <- matrix(integer(), ncol = 2L)
    }
  }
  pair.index <- as.matrix(pair.index)
  if (ncol(pair.index) != 2L) {
    stop("pair.index must be a matrix with two columns", call. = FALSE)
  }
  if (length(pair.index) &&
      (any(!is.finite(pair.index)) || any(pair.index != floor(pair.index)) ||
        any(pair.index < 1L) || any(pair.index > n))) {
    stop("pair.index must contain valid 1-based vertex indices", call. = FALSE)
  }
  storage.mode(pair.index) <- "integer"
  pair.index
}

grip.edge.repulsive.state.R <- function(coords,
                                        edges,
                                        edge.lengths,
                                        edge.weights = NULL,
                                        edge.family = "quadratic",
                                        eps.plus = 0.35,
                                        beta = 0,
                                        lambda = 0,
                                        pair.index = NULL,
                                        pair.weights = NULL,
                                        repulsion.family = "log",
                                        repulsion.delta = 1e-3,
                                        repulsion.power = 1,
                                        distance.eps = 1e-10) {
  Z <- as.matrix(coords)
  n <- nrow(Z)
  dim <- ncol(Z)
  if (is.null(edge.weights)) edge.weights <- rep(1, nrow(edges))
  if (!nrow(edges)) {
    stop("edges must contain at least one edge", call. = FALSE)
  }
  if (length(edge.lengths) != nrow(edges)) {
    stop("edge.lengths must be parallel to edges", call. = FALSE)
  }
  diffs <- Z[edges[, 1L], , drop = FALSE] - Z[edges[, 2L], , drop = FALSE]
  d <- sqrt(rowSums(diffs^2) + distance.eps^2)
  s <- d / edge.lengths
  ep <- .grip.edge.repulsive.potential.state(s, family = edge.family, eps.plus = eps.plus, beta = beta)
  feasible <- all(ep$feasible)
  edge.energy <- sum(edge.weights * ep$value)
  grad <- matrix(0, n, dim)
  if (is.finite(edge.energy)) {
    coeff <- edge.weights * ep$deriv / edge.lengths / d
    for (e in seq_len(nrow(edges))) {
      g <- coeff[[e]] * diffs[e, ]
      grad[edges[e, 1L], ] <- grad[edges[e, 1L], ] + g
      grad[edges[e, 2L], ] <- grad[edges[e, 2L], ] - g
    }
  }

  repel.energy <- 0
  if (lambda > 0) {
    if (is.null(pair.index)) {
      pair.index <- t(utils::combn(n, 2L))
    }
    if (is.null(pair.weights)) pair.weights <- rep(1, nrow(pair.index))
    pdiff <- Z[pair.index[, 1L], , drop = FALSE] - Z[pair.index[, 2L], , drop = FALSE]
    pr <- sqrt(rowSums(pdiff^2) + distance.eps^2)
    rp <- .grip.edge.repulsive.repulsion.state(
      pr,
      family = repulsion.family,
      delta = repulsion.delta,
      power = repulsion.power
    )
    repel.energy <- sum(pair.weights * rp$value)
    coeff <- lambda * pair.weights * rp$deriv.r / pr
    for (p in seq_len(nrow(pair.index))) {
      g <- coeff[[p]] * pdiff[p, ]
      grad[pair.index[p, 1L], ] <- grad[pair.index[p, 1L], ] + g
      grad[pair.index[p, 2L], ] <- grad[pair.index[p, 2L], ] - g
    }
  }

  list(
    energy = edge.energy + lambda * repel.energy,
    edge.energy = edge.energy,
    repel.energy = repel.energy,
    gradient = grad,
    gradient.norm = sqrt(sum(grad^2)),
    feasible = feasible,
    edge.embedded.lengths = d,
    edge.relative.lengths = s,
    edge.residuals = d - edge.lengths,
    n.wall.violations = sum(!ep$feasible)
  )
}

grip.optimize.edge.repulsive.stage.R <- function(coords,
                                                 edges,
                                                 edge.lengths,
                                                 edge.weights = NULL,
                                                 lambda = 0,
                                                 edge.family = "quadratic",
                                                 eps.plus = 0.35,
                                                 beta = 0,
                                                 pair.index = NULL,
                                                 pair.weights = NULL,
                                                 repulsion.family = "log",
                                                 repulsion.delta = 1e-3,
                                                 repulsion.power = 1,
                                                 max.iter = 80L,
                                                 initial.step = 0.02,
                                                 step.shrink = 0.5,
                                                 armijo = 1e-4,
                                                 min.step = 1e-8,
                                                 grad.tol = 1e-7,
                                                 recenter = TRUE,
                                                 distance.eps = 1e-10,
                                                 return.frames = FALSE) {
  Z <- .grip.edge.repulsive.center.coords(coords)
  trace <- list()
  frames <- if (isTRUE(return.frames)) list(Z) else NULL
  state <- grip.edge.repulsive.state.R(
    Z, edges, edge.lengths, edge.weights = edge.weights,
    edge.family = edge.family, eps.plus = eps.plus, beta = beta,
    lambda = lambda, pair.index = pair.index, pair.weights = pair.weights,
    repulsion.family = repulsion.family, repulsion.delta = repulsion.delta,
    repulsion.power = repulsion.power, distance.eps = distance.eps
  )
  frame.index <- if (isTRUE(return.frames)) 0L else NA_integer_
  for (iter in seq_len(max.iter)) {
    trace[[length(trace) + 1L]] <- data.frame(
      iteration = iter - 1L,
      energy = state$energy,
      edge.energy = state$edge.energy,
      repel.energy = state$repel.energy,
      gradient.norm = state$gradient.norm,
      step = NA_real_,
      accepted = TRUE,
      n.wall.violations = state$n.wall.violations,
      frame.index = frame.index,
      stringsAsFactors = FALSE
    )
    if (!is.finite(state$energy) || !is.finite(state$gradient.norm) ||
        state$gradient.norm <= grad.tol) {
      break
    }
    step <- initial.step
    accepted <- FALSE
    candidate <- Z
    candidate.state <- state
    while (is.finite(step) && step >= min.step) {
      proposal <- Z - step * state$gradient
      if (isTRUE(recenter)) proposal <- .grip.edge.repulsive.center.coords(proposal)
      proposal.state <- grip.edge.repulsive.state.R(
        proposal, edges, edge.lengths, edge.weights = edge.weights,
        edge.family = edge.family, eps.plus = eps.plus, beta = beta,
        lambda = lambda, pair.index = pair.index, pair.weights = pair.weights,
        repulsion.family = repulsion.family, repulsion.delta = repulsion.delta,
        repulsion.power = repulsion.power, distance.eps = distance.eps
      )
      target <- state$energy - armijo * step * state$gradient.norm^2
      if (proposal.state$feasible && is.finite(proposal.state$energy) &&
          proposal.state$energy <= target) {
        accepted <- TRUE
        candidate <- proposal
        candidate.state <- proposal.state
        break
      }
      step <- step * step.shrink
    }
    if (!accepted) break
    Z <- candidate
    state <- candidate.state
    if (isTRUE(return.frames)) {
      frames[[length(frames) + 1L]] <- Z
      frame.index <- length(frames) - 1L
    }
    trace[[length(trace) + 1L]] <- data.frame(
      iteration = iter,
      energy = state$energy,
      edge.energy = state$edge.energy,
      repel.energy = state$repel.energy,
      gradient.norm = state$gradient.norm,
      step = step,
      accepted = TRUE,
      n.wall.violations = state$n.wall.violations,
      frame.index = frame.index,
      stringsAsFactors = FALSE
    )
  }
  list(coords = Z, state = state, trace = do.call(rbind, trace), frames = frames)
}

.grip.edge.repulsive.cpp.available <- function() {
  exists("grip_edge_repulsive_state_cpp", mode = "function", inherits = TRUE) &&
    exists("grip_optimize_edge_repulsive_stage_cpp", mode = "function", inherits = TRUE) &&
    exists("grip_repulsive_state_cpp", mode = "function", inherits = TRUE) &&
    exists("grip_optimize_repulsive_stage_cpp", mode = "function", inherits = TRUE)
}

#' Evaluate a pure repulsive layout objective
#'
#' `grip.repulsive.state()` evaluates only the repulsion term used by the
#' experimental GMDS repulsive-unfolding operators. Unlike
#' [grip.edge.repulsive.state()], this function has no graph-edge term and does
#' not require placeholder edges or edge lengths.
#'
#' @param coords Numeric `n` by `dim` coordinate matrix.
#' @param lambda Repulsion strength. The returned total energy is
#'   `lambda * repel.energy`; the gradient is scaled by `lambda`.
#' @param pair.index Optional two-column matrix of 1-based vertex pairs. If
#'   `NULL` and `lambda > 0`, all unordered pairs are used.
#' @param pair.weights Optional non-negative pair weights parallel to
#'   `pair.index`.
#' @param repulsion.family Repulsion potential family, either `"log"` or
#'   `"inverse_power"`.
#' @param repulsion.delta Small positive softening parameter for pair distances.
#' @param repulsion.power Power used by the `"inverse_power"` repulsion.
#' @param distance.eps Small positive distance floor used in derivatives.
#' @param engine Backend engine. `"cpp"` is the default; `"R"` uses the
#'   reference implementation.
#'
#' @return A list containing energy, unscaled repulsion energy, gradient,
#'   gradient norm, and embedded pair lengths.
#' @export
grip.repulsive.state <- function(coords,
                                 lambda = 1,
                                 pair.index = NULL,
                                 pair.weights = NULL,
                                 repulsion.family = c("log", "inverse_power"),
                                 repulsion.delta = 1e-3,
                                 repulsion.power = 1,
                                 distance.eps = 1e-10,
                                 engine = c("cpp", "R")) {
  engine <- match.arg(engine)
  repulsion.family <- match.arg(repulsion.family)
  Z <- as.matrix(coords)
  if (!is.numeric(Z) || nrow(Z) < 1L || ncol(Z) < 1L) {
    stop("coords must be a non-empty numeric matrix", call. = FALSE)
  }
  lambda <- as.numeric(lambda)
  pair.index <- .grip.edge.repulsive.normalize.pairs(pair.index, nrow(Z), lambda)
  if (is.null(pair.weights)) pair.weights <- rep(1, nrow(pair.index))
  pair.weights <- as.numeric(pair.weights)
  if (length(pair.weights) != nrow(pair.index) || any(!is.finite(pair.weights)) ||
      any(pair.weights < 0)) {
    stop("pair.weights must be finite non-negative values parallel to pair.index", call. = FALSE)
  }
  if (identical(engine, "cpp") && .grip.edge.repulsive.cpp.available()) {
    return(grip_repulsive_state_cpp(
      coords = Z,
      lambda = lambda,
      pair_index = pair.index,
      pair_weights = pair.weights,
      repulsion_family = repulsion.family,
      repulsion_delta = repulsion.delta,
      repulsion_power = repulsion.power,
      distance_eps = distance.eps
    ))
  }
  .grip.repulsive.state.R(
    coords = Z, lambda = lambda, pair.index = pair.index,
    pair.weights = pair.weights, repulsion.family = repulsion.family,
    repulsion.delta = repulsion.delta, repulsion.power = repulsion.power,
    distance.eps = distance.eps
  )
}

#' Optimize one pure repulsive layout stage
#'
#' `grip.optimize.repulsive.stage()` performs a gradient-descent stage for the
#' pure repulsion objective evaluated by [grip.repulsive.state()]. It is intended
#' for GMDS experiments that need to inspect repulsion separately from
#' edge-length repair.
#'
#' @inheritParams grip.repulsive.state
#' @param max.iter Maximum number of gradient-descent iterations.
#' @param initial.step Initial gradient step size.
#' @param step.shrink Multiplicative shrink factor used by backtracking.
#' @param armijo Armijo sufficient-decrease coefficient.
#' @param min.step Minimum allowable step before the stage stops.
#' @param grad.tol Gradient-norm stopping tolerance.
#' @param recenter Logical; whether to recenter coordinates after each proposal.
#' @param return.frames Logical; whether to return accepted coordinate frames.
#'   Frame 0 is the starting coordinate matrix and later frames are accepted
#'   updates.
#'
#' @return A list with `coords`, final `state`, a data-frame `trace`, and,
#'   when requested, a list of coordinate `frames`.
#' @export
grip.optimize.repulsive.stage <- function(coords,
                                          lambda = 1,
                                          pair.index = NULL,
                                          pair.weights = NULL,
                                          repulsion.family = c("log", "inverse_power"),
                                          repulsion.delta = 1e-3,
                                          repulsion.power = 1,
                                          max.iter = 80L,
                                          initial.step = 0.02,
                                          step.shrink = 0.5,
                                          armijo = 1e-4,
                                          min.step = 1e-8,
                                          grad.tol = 1e-7,
                                          recenter = TRUE,
                                          distance.eps = 1e-10,
                                          return.frames = FALSE,
                                          engine = c("cpp", "R")) {
  engine <- match.arg(engine)
  repulsion.family <- match.arg(repulsion.family)
  Z <- as.matrix(coords)
  if (!is.numeric(Z) || nrow(Z) < 1L || ncol(Z) < 1L) {
    stop("coords must be a non-empty numeric matrix", call. = FALSE)
  }
  lambda <- as.numeric(lambda)
  pair.index <- .grip.edge.repulsive.normalize.pairs(pair.index, nrow(Z), lambda)
  if (is.null(pair.weights)) pair.weights <- rep(1, nrow(pair.index))
  pair.weights <- as.numeric(pair.weights)
  if (length(pair.weights) != nrow(pair.index) || any(!is.finite(pair.weights)) ||
      any(pair.weights < 0)) {
    stop("pair.weights must be finite non-negative values parallel to pair.index", call. = FALSE)
  }
  if (identical(engine, "cpp") && .grip.edge.repulsive.cpp.available()) {
    return(grip_optimize_repulsive_stage_cpp(
      coords = Z,
      lambda = lambda,
      pair_index = pair.index,
      pair_weights = pair.weights,
      repulsion_family = repulsion.family,
      repulsion_delta = repulsion.delta,
      repulsion_power = repulsion.power,
      max_iter = as.integer(max.iter),
      initial_step = initial.step,
      step_shrink = step.shrink,
      armijo = armijo,
      min_step = min.step,
      grad_tol = grad.tol,
      recenter = recenter,
      distance_eps = distance.eps,
      return_frames = isTRUE(return.frames)
    ))
  }
  Z <- .grip.edge.repulsive.center.coords(Z)
  trace <- list()
  frames <- if (isTRUE(return.frames)) list(Z) else NULL
  frame.index <- if (isTRUE(return.frames)) 0L else NA_integer_
  state <- .grip.repulsive.state.R(
    Z, lambda = lambda, pair.index = pair.index, pair.weights = pair.weights,
    repulsion.family = repulsion.family, repulsion.delta = repulsion.delta,
    repulsion.power = repulsion.power, distance.eps = distance.eps
  )
  for (iter in seq_len(max.iter)) {
    trace[[length(trace) + 1L]] <- data.frame(
      iteration = iter - 1L,
      energy = state$energy,
      repel.energy = state$repel.energy,
      gradient.norm = state$gradient.norm,
      step = NA_real_,
      accepted = TRUE,
      frame.index = frame.index,
      stringsAsFactors = FALSE
    )
    if (!is.finite(state$energy) || !is.finite(state$gradient.norm) ||
        state$gradient.norm <= grad.tol) {
      break
    }
    step <- initial.step
    accepted <- FALSE
    candidate <- Z
    candidate.state <- state
    while (is.finite(step) && step >= min.step) {
      proposal <- Z - step * state$gradient
      if (isTRUE(recenter)) proposal <- .grip.edge.repulsive.center.coords(proposal)
      proposal.state <- .grip.repulsive.state.R(
        proposal, lambda = lambda, pair.index = pair.index,
        pair.weights = pair.weights, repulsion.family = repulsion.family,
        repulsion.delta = repulsion.delta, repulsion.power = repulsion.power,
        distance.eps = distance.eps
      )
      target <- state$energy - armijo * step * state$gradient.norm^2
      if (is.finite(proposal.state$energy) && proposal.state$energy <= target) {
        accepted <- TRUE
        candidate <- proposal
        candidate.state <- proposal.state
        break
      }
      step <- step * step.shrink
    }
    if (!accepted) break
    Z <- candidate
    state <- candidate.state
    if (isTRUE(return.frames)) {
      frames[[length(frames) + 1L]] <- Z
      frame.index <- length(frames) - 1L
    }
    trace[[length(trace) + 1L]] <- data.frame(
      iteration = iter,
      energy = state$energy,
      repel.energy = state$repel.energy,
      gradient.norm = state$gradient.norm,
      step = step,
      accepted = TRUE,
      frame.index = frame.index,
      stringsAsFactors = FALSE
    )
  }
  list(coords = Z, state = state, trace = do.call(rbind, trace), frames = frames)
}

#' Evaluate an edge-isometric repulsive unfolding objective
#'
#' `grip.edge.repulsive.state()` evaluates an experimental GMDS objective that
#' combines an edge-isometric term with a configurable repulsion term. The edge
#' term penalizes deviation of embedded edge lengths from graph edge lengths.
#' The repulsion term can be evaluated on all vertex pairs or on a selected pair
#' set.
#'
#' The default `engine = "cpp"` uses the compiled backend. `engine = "R"` keeps
#' a reference implementation available for diagnostics and regression tests.
#'
#' @param coords Numeric `n` by `dim` coordinate matrix.
#' @param edges Integer or numeric matrix with two columns containing 1-based
#'   graph edge endpoints.
#' @param edge.lengths Numeric vector of target edge lengths, parallel to
#'   `edges`.
#' @param edge.weights Optional non-negative edge weights. Defaults to one.
#' @param edge.family Edge potential family, either `"quadratic"` or
#'   `"upper_barrier"`.
#' @param eps.plus Upper-barrier slack parameter.
#' @param beta Upper-barrier strength. When `beta <= 0`, the edge potential is
#'   quadratic.
#' @param lambda Repulsion strength.
#' @param pair.index Optional two-column matrix of 1-based vertex pairs for the
#'   repulsion term. If `NULL` and `lambda > 0`, all unordered pairs are used.
#' @param pair.weights Optional pair weights, parallel to `pair.index`.
#' @param repulsion.family Repulsion potential family, either `"log"` or
#'   `"inverse_power"`.
#' @param repulsion.delta Small positive softening parameter for pair distances.
#' @param repulsion.power Power used by the `"inverse_power"` repulsion.
#' @param distance.eps Small positive distance floor used in derivatives.
#' @param engine Backend engine. `"cpp"` is the default; `"R"` uses the
#'   reference implementation.
#'
#' @return A list containing total energy, edge and repulsion energies,
#'   gradient, gradient norm, feasibility flag, embedded edge lengths, relative
#'   edge lengths, edge residuals, and number of upper-barrier wall violations.
#' @export
grip.edge.repulsive.state <- function(coords,
                                      edges,
                                      edge.lengths,
                                      edge.weights = NULL,
                                      edge.family = c("quadratic", "upper_barrier"),
                                      eps.plus = 0.35,
                                      beta = 0,
                                      lambda = 0,
                                      pair.index = NULL,
                                      pair.weights = NULL,
                                      repulsion.family = c("log", "inverse_power"),
                                      repulsion.delta = 1e-3,
                                      repulsion.power = 1,
                                      distance.eps = 1e-10,
                                      engine = c("cpp", "R")) {
  engine <- match.arg(engine)
  edge.family <- match.arg(edge.family)
  repulsion.family <- match.arg(repulsion.family)
  Z <- as.matrix(coords)
  if (!is.numeric(Z) || nrow(Z) < 1L || ncol(Z) < 1L) {
    stop("coords must be a non-empty numeric matrix", call. = FALSE)
  }
  edges <- .grip.edge.repulsive.normalize.edges(edges)
  if (max(edges) > nrow(Z)) {
    stop("edges contain vertex indices outside coords", call. = FALSE)
  }
  edge.lengths <- as.numeric(edge.lengths)
  if (length(edge.lengths) != nrow(edges) || any(!is.finite(edge.lengths)) ||
      any(edge.lengths <= 0)) {
    stop("edge.lengths must be positive finite values parallel to edges", call. = FALSE)
  }
  if (is.null(edge.weights)) edge.weights <- rep(1, nrow(edges))
  edge.weights <- as.numeric(edge.weights)
  if (length(edge.weights) != nrow(edges) || any(!is.finite(edge.weights)) ||
      any(edge.weights < 0)) {
    stop("edge.weights must be finite non-negative values parallel to edges", call. = FALSE)
  }
  lambda <- as.numeric(lambda)
  pair.index <- .grip.edge.repulsive.normalize.pairs(pair.index, nrow(Z), lambda)
  if (is.null(pair.weights)) pair.weights <- rep(1, nrow(pair.index))
  pair.weights <- as.numeric(pair.weights)
  if (length(pair.weights) != nrow(pair.index) || any(!is.finite(pair.weights)) ||
      any(pair.weights < 0)) {
    stop("pair.weights must be finite non-negative values parallel to pair.index", call. = FALSE)
  }
  if (identical(engine, "cpp") && .grip.edge.repulsive.cpp.available()) {
    return(grip_edge_repulsive_state_cpp(
      coords = Z,
      edges = edges,
      edge_lengths = edge.lengths,
      edge_weights = edge.weights,
      edge_family = edge.family,
      eps_plus = eps.plus,
      beta = beta,
      lambda = lambda,
      pair_index = pair.index,
      pair_weights = pair.weights,
      repulsion_family = repulsion.family,
      repulsion_delta = repulsion.delta,
      repulsion_power = repulsion.power,
      distance_eps = distance.eps
    ))
  }
  grip.edge.repulsive.state.R(
    coords = Z, edges = edges, edge.lengths = edge.lengths,
    edge.weights = edge.weights, edge.family = edge.family, eps.plus = eps.plus,
    beta = beta, lambda = lambda, pair.index = pair.index,
    pair.weights = pair.weights, repulsion.family = repulsion.family,
    repulsion.delta = repulsion.delta, repulsion.power = repulsion.power,
    distance.eps = distance.eps
  )
}

#' Optimize one edge-isometric repulsive unfolding stage
#'
#' `grip.optimize.edge.repulsive.stage()` performs one gradient-descent stage
#' for the objective evaluated by `grip.edge.repulsive.state()`. It uses Armijo
#' backtracking, optionally recenters coordinates after each proposal, and
#' returns the final coordinates plus a per-iteration trace.
#'
#' @inheritParams grip.edge.repulsive.state
#' @param max.iter Maximum number of gradient-descent iterations.
#' @param initial.step Initial gradient step size.
#' @param step.shrink Multiplicative shrink factor used by backtracking.
#' @param armijo Armijo sufficient-decrease coefficient.
#' @param min.step Minimum allowable step before the stage stops.
#' @param grad.tol Gradient-norm stopping tolerance.
#' @param recenter Logical; whether to recenter coordinates after each proposal.
#' @param return.frames Logical; whether to return accepted coordinate frames.
#'   Frame 0 is the starting coordinate matrix and later frames are accepted
#'   updates.
#'
#' @return A list with `coords`, final `state`, a data-frame `trace`, and,
#'   when requested, a list of coordinate `frames`.
#' @export
grip.optimize.edge.repulsive.stage <- function(coords,
                                               edges,
                                               edge.lengths,
                                               edge.weights = NULL,
                                               lambda = 0,
                                               edge.family = c("quadratic", "upper_barrier"),
                                               eps.plus = 0.35,
                                               beta = 0,
                                               pair.index = NULL,
                                               pair.weights = NULL,
                                               repulsion.family = c("log", "inverse_power"),
                                               repulsion.delta = 1e-3,
                                               repulsion.power = 1,
                                               max.iter = 80L,
                                               initial.step = 0.02,
                                               step.shrink = 0.5,
                                               armijo = 1e-4,
                                               min.step = 1e-8,
                                               grad.tol = 1e-7,
                                               recenter = TRUE,
                                               distance.eps = 1e-10,
                                               return.frames = FALSE,
                                               engine = c("cpp", "R")) {
  engine <- match.arg(engine)
  edge.family <- match.arg(edge.family)
  repulsion.family <- match.arg(repulsion.family)
  Z <- as.matrix(coords)
  if (!is.numeric(Z) || nrow(Z) < 1L || ncol(Z) < 1L) {
    stop("coords must be a non-empty numeric matrix", call. = FALSE)
  }
  edges <- .grip.edge.repulsive.normalize.edges(edges)
  if (max(edges) > nrow(Z)) {
    stop("edges contain vertex indices outside coords", call. = FALSE)
  }
  edge.lengths <- as.numeric(edge.lengths)
  if (length(edge.lengths) != nrow(edges) || any(!is.finite(edge.lengths)) ||
      any(edge.lengths <= 0)) {
    stop("edge.lengths must be positive finite values parallel to edges", call. = FALSE)
  }
  if (is.null(edge.weights)) edge.weights <- rep(1, nrow(edges))
  edge.weights <- as.numeric(edge.weights)
  if (length(edge.weights) != nrow(edges) || any(!is.finite(edge.weights)) ||
      any(edge.weights < 0)) {
    stop("edge.weights must be finite non-negative values parallel to edges", call. = FALSE)
  }
  lambda <- as.numeric(lambda)
  pair.index <- .grip.edge.repulsive.normalize.pairs(pair.index, nrow(Z), lambda)
  if (is.null(pair.weights)) pair.weights <- rep(1, nrow(pair.index))
  pair.weights <- as.numeric(pair.weights)
  if (length(pair.weights) != nrow(pair.index) || any(!is.finite(pair.weights)) ||
      any(pair.weights < 0)) {
    stop("pair.weights must be finite non-negative values parallel to pair.index", call. = FALSE)
  }
  if (identical(engine, "cpp") && .grip.edge.repulsive.cpp.available()) {
    return(grip_optimize_edge_repulsive_stage_cpp(
      coords = Z,
      edges = edges,
      edge_lengths = edge.lengths,
      edge_weights = edge.weights,
      lambda = lambda,
      edge_family = edge.family,
      eps_plus = eps.plus,
      beta = beta,
      pair_index = pair.index,
      pair_weights = pair.weights,
      repulsion_family = repulsion.family,
      repulsion_delta = repulsion.delta,
      repulsion_power = repulsion.power,
      max_iter = as.integer(max.iter),
      initial_step = initial.step,
      step_shrink = step.shrink,
      armijo = armijo,
      min_step = min.step,
      grad_tol = grad.tol,
      recenter = recenter,
      distance_eps = distance.eps,
      return_frames = isTRUE(return.frames)
    ))
  }
  grip.optimize.edge.repulsive.stage.R(
    coords = Z, edges = edges, edge.lengths = edge.lengths,
    edge.weights = edge.weights, lambda = lambda, edge.family = edge.family,
    eps.plus = eps.plus, beta = beta, pair.index = pair.index,
    pair.weights = pair.weights, repulsion.family = repulsion.family,
    repulsion.delta = repulsion.delta, repulsion.power = repulsion.power,
    max.iter = max.iter, initial.step = initial.step,
    step.shrink = step.shrink, armijo = armijo, min.step = min.step,
    grad.tol = grad.tol, recenter = recenter, distance.eps = distance.eps,
    return.frames = isTRUE(return.frames)
  )
}

.grip.edge.repulsive.finite.difference.gradient.check <- function() {
  set.seed(1)
  Z <- matrix(stats::rnorm(15), ncol = 3)
  edges <- rbind(c(1L, 2L), c(2L, 3L), c(3L, 4L), c(4L, 5L))
  ell <- rep(1, 4)
  st <- grip.edge.repulsive.state.R(Z, edges, ell, lambda = 0.05)
  i <- 2L
  j <- 1L
  eps <- 1e-6
  Zp <- Z
  Zm <- Z
  Zp[i, j] <- Zp[i, j] + eps
  Zm[i, j] <- Zm[i, j] - eps
  fp <- grip.edge.repulsive.state.R(Zp, edges, ell, lambda = 0.05)$energy
  fm <- grip.edge.repulsive.state.R(Zm, edges, ell, lambda = 0.05)$energy
  fd <- (fp - fm) / (2 * eps)
  c(analytic = st$gradient[i, j], finite_difference = fd, abs_error = abs(fd - st$gradient[i, j]))
}
