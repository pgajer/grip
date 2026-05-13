test_that("edge-length density stiffness handles uniform and continuation cases", {
  equal <- grip.edge.length.density.stiffness(rep(2, 5), method = "density")
  expect_s3_class(equal, "grip_edge_length_stiffness")
  expect_equal(equal$stiffness, rep(1, 5), tolerance = 1e-12)
  expect_equal(mean(equal$stiffness), 1, tolerance = 1e-12)

  mixed <- grip.edge.length.density.stiffness(
    c(rep(1, 20), rep(4, 4)),
    method = "density",
    mix = 1
  )
  expect_equal(mixed$stiffness, rep(1, 24), tolerance = 1e-12)
})

test_that("edge-length density stiffness emphasizes modal edge scale", {
  weights <- c(rep(1, 40), rep(4, 8), rep(8, 4))
  stiff <- grip.edge.length.density.stiffness(
    weights,
    method = "density",
    mix = 0,
    transform = "identity"
  )

  expect_equal(mean(stiff$stiffness), 1, tolerance = 1e-12)
  expect_gt(mean(stiff$stiffness[weights == 1]), mean(stiff$stiffness[weights == 8]))
  expect_lt(abs(stiff$mode - 1), 0.5)
})

test_that("edge-length stiffness clipping is respected after normalization", {
  weights <- c(rep(1, 20), rep(5, 5))
  stiff <- grip.edge.length.density.stiffness(
    weights,
    method = "density",
    stiffness_floor = 0.5,
    stiffness_ceiling = 1.5
  )

  expect_equal(mean(stiff$stiffness), 1, tolerance = 1e-12)
  expect_gte(min(stiff$mixed_signal), 0.5)
  expect_lte(max(stiff$mixed_signal), 1.5)
})

test_that("edge-isometric energy gradient matches finite differences", {
  edges <- rbind(c(1L, 2L), c(2L, 3L), c(3L, 4L))
  coords <- matrix(c(
    0.0, 0.0,
    1.2, 0.1,
    1.9, 0.8,
    3.1, 0.2
  ), ncol = 2, byrow = TRUE)
  edge.weights <- c(1, 1.2, 1.1)
  stiffness <- c(0.7, 1.4, 0.9)
  state <- grip:::grip.edge.isometric.energy.gradient(
    coords = coords,
    edges = edges,
    edge_weights = edge.weights,
    stiffness = stiffness,
    scale = 1.1,
    edge_length_epsilon = 1e-8
  )

  eps <- 1e-6
  fd <- matrix(0, nrow(coords), ncol(coords))
  for (i in seq_len(nrow(coords))) {
    for (j in seq_len(ncol(coords))) {
      plus <- coords
      minus <- coords
      plus[i, j] <- plus[i, j] + eps
      minus[i, j] <- minus[i, j] - eps
      e.plus <- grip:::grip.edge.isometric.energy.gradient(
        plus, edges, edge.weights, stiffness, scale = 1.1
      )$energy
      e.minus <- grip:::grip.edge.isometric.energy.gradient(
        minus, edges, edge.weights, stiffness, scale = 1.1
      )$energy
      fd[i, j] <- (e.plus - e.minus) / (2 * eps)
    }
  }

  expect_equal(state$gradient, fd, tolerance = 1e-5)
})

test_that("edge-isometric optimizer preserves exact weighted path layouts", {
  prepared <- grip.prepare.graph.geodesic.mds(
    edges = edges.path(4L),
    n = 4L,
    edge_weights = c(1, 2, 1)
  )
  coords <- cbind(c(0, 1, 3, 4), 0)
  fit <- grip.optimize.edge.isometric.layout(
    coords = coords,
    prepared = prepared,
    dim = 2L,
    stiffness_method = "uniform",
    density_mix_schedule = 1,
    scale_mode = "identity",
    max_iter = 5L
  )

  expect_s3_class(fit, "grip_gmds_layout")
  expect_equal(fit$method, "edge_isometric_gkk")
  expect_lt(fit$diagnostics$edge.rel.rmse[[1L]], 1e-8)
  expect_lt(fit$diagnostics$gmds.stress[[1L]], 1e-8)
})

test_that("edge-isometric optimizer decreases edge error from perturbed layout", {
  prepared <- grip.prepare.graph.geodesic.mds(
    edges = edges.path(5L),
    n = 5L,
    edge_weights = rep(1, 4L)
  )
  start <- cbind(c(0, 0.7, 1.9, 2.4, 4.2), c(0, 0.4, -0.2, 0.5, -0.1))
  before <- grip.score.gmds.layout(
    start,
    prepared = prepared,
    scale_mode = "identity"
  )
  fit <- grip.optimize.edge.isometric.layout(
    coords = start,
    prepared = prepared,
    dim = 2L,
    stiffness_method = "uniform",
    density_mix_schedule = 1,
    scale_mode = "identity",
    max_iter = 80L,
    initial_step = 0.25,
    return_trace = TRUE
  )

  expect_lt(fit$diagnostics$edge.rel.rmse[[1L]], before$edge.rel.rmse[[1L]])
  expect_true(is.data.frame(fit$trace))
  expect_true(all(c("stage", "mix", "energy", "edge.rel.rmse") %in% names(fit$trace)))
  expect_true(is.data.frame(fit$metadata$stage_summaries))
})

test_that("C++ edge-isometric optimizer matches R reference engine", {
  prepared <- grip.prepare.graph.geodesic.mds(
    edges = rbind(c(1L, 2L), c(2L, 3L), c(3L, 4L), c(4L, 1L), c(1L, 3L)),
    n = 4L,
    edge_weights = c(1, 1.4, 1, 1.3, 1.8)
  )
  start <- matrix(c(
    0.0, 0.0,
    0.8, 0.4,
    1.7, 1.0,
    -0.2, 1.1
  ), ncol = 2, byrow = TRUE)
  args <- list(
    coords = start,
    prepared = prepared,
    dim = 2L,
    stiffness_method = "density",
    stiffness_transform = "sqrt",
    density_mix_schedule = c(0, 0.5, 1),
    scale_mode = "profiled",
    max_iter = 12L,
    initial_step = 0.2,
    return_trace = TRUE
  )
  fit.cpp <- do.call(grip.optimize.edge.isometric.layout, c(args, list(engine = "cpp")))
  fit.r <- do.call(grip.optimize.edge.isometric.layout, c(args, list(engine = "R")))

  expect_equal(fit.cpp$coords, fit.r$coords, tolerance = 1e-10)
  expect_equal(fit.cpp$trace$energy, fit.r$trace$energy, tolerance = 1e-10)
  expect_equal(fit.cpp$trace$edge.rel.rmse, fit.r$trace$edge.rel.rmse, tolerance = 1e-10)
  expect_equal(fit.cpp$metadata$engine, "cpp_gradient_descent_armijo")
  expect_equal(fit.r$metadata$engine, "r_gradient_descent_armijo")
})

test_that("C++ edge-isometric optimizer supports fixed and user scale modes", {
  prepared <- grip.prepare.graph.geodesic.mds(
    edges = edges.path(4L),
    n = 4L,
    edge_weights = c(1, 2, 1)
  )
  start <- 2 * cbind(c(0, 1.1, 3.2, 4.0), c(0, 0.2, -0.1, 0.1))
  fixed <- grip.optimize.edge.isometric.layout(
    coords = start,
    prepared = prepared,
    dim = 2L,
    stiffness_method = "uniform",
    density_mix_schedule = c(0, 1),
    scale_mode = "fixed_initial",
    max_iter = 3L,
    engine = "cpp"
  )
  user <- grip.optimize.edge.isometric.layout(
    coords = start,
    prepared = prepared,
    dim = 2L,
    stiffness_method = "uniform",
    density_mix_schedule = 1,
    scale_mode = "user",
    scale = 2,
    max_iter = 3L,
    engine = "cpp"
  )

  expect_s3_class(fixed, "grip_gmds_layout")
  expect_s3_class(user, "grip_gmds_layout")
  expect_true(all(is.finite(fixed$coords)))
  expect_true(all(is.finite(user$trace$edge.scale)))
  expect_equal(unique(user$trace$edge.scale), 2)
})
