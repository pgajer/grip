test_that("edge-length density stiffness handles uniform and continuation cases", {
  equal <- edge.length.density.stiffness(rep(2, 5), method = "density")
  expect_s3_class(equal, "grip_edge_length_stiffness")
  expect_equal(equal$stiffness, rep(1, 5), tolerance = 1e-12)
  expect_equal(mean(equal$stiffness), 1, tolerance = 1e-12)

  mixed <- edge.length.density.stiffness(
    c(rep(1, 20), rep(4, 4)),
    method = "density",
    mix = 1
  )
  expect_equal(mixed$stiffness, rep(1, 24), tolerance = 1e-12)
})

test_that("edge-length density stiffness emphasizes modal edge scale", {
  weights <- c(rep(1, 40), rep(4, 8), rep(8, 4))
  stiff <- edge.length.density.stiffness(
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
  stiff <- edge.length.density.stiffness(
    weights,
    method = "density",
    stiffness_floor = 0.5,
    stiffness_ceiling = 1.5
  )

  expect_equal(mean(stiff$stiffness), 1, tolerance = 1e-12)
  expect_gte(min(stiff$mixed_signal), 0.5)
  expect_lte(max(stiff$mixed_signal), 1.5)
})

test_that("edge-only edge-KK preparation avoids all-pairs caches", {
  edges <- rbind(c(1L, 2L), c(2L, 3L), c(3L, 4L), c(1L, 4L))
  edge.weights <- c(1, 1.3, 0.9, 1.7)
  prepared <- prepare.edge.kk(
    edges = edges,
    n = 4L,
    edge_weights = edge.weights
  )

  expect_s3_class(prepared, "grip_edge_kk_prepared")
  expect_s3_class(prepared, "grip_gmds_prepared")
  expect_equal(prepared$n, 4L)
  expect_equal(
    prepared$edges,
    rbind(c(1L, 2L), c(1L, 4L), c(2L, 3L), c(3L, 4L))
  )
  expect_equal(prepared$edge_targets, c(1, 1.7, 1.3, 0.9))
  expect_equal(prepared$pair_mode, "edge_only")
  expect_equal(nrow(prepared$pair_matrix), 0L)
  expect_length(prepared$pair_graph_distance, 0L)
  expect_null(prepared$distance_matrix)
  expect_equal(prepared$n_components, 1L)
})

test_that("edge-only prepared objects support weighted-GRIP to edge-KK repair", {
  edges <- rbind(c(1L, 2L), c(2L, 3L), c(3L, 4L), c(1L, 4L), c(1L, 3L))
  edge.weights <- c(1, 1.4, 1.1, 1.2, 1.8)
  prepared <- prepare.edge.kk(
    edges = edges,
    n = 4L,
    edge_weights = edge.weights
  )
  init <- grip(metric = "edge_length",
    edges = edges,
    n = 4L,
    edge_weights = edge.weights,
    dim = 3L,
    rounds = 4L,
    final_rounds = 4L,
    seed = 11L
  )
  before <- score.gmds(init, prepared = prepared)
  fit <- edge.kk(
    coords = init,
    prepared = prepared,
    dim = 3L,
    stiffness_method = "uniform",
    density_mix_schedule = 1,
    scale_mode = "profiled",
    max_iter = 4L,
    diagnostics = TRUE,
    engine = "cpp"
  )

  expect_s3_class(fit, "grip_gmds_layout")
  expect_equal(fit$method, "edge_kk")
  expect_equal(dim(fit$coords), c(4L, 3L))
  expect_s3_class(fit$prepared, "grip_edge_kk_prepared")
  expect_true(is.na(before$gmds.stress[[1L]]))
  expect_true(is.na(fit$diagnostics$gmds.stress[[1L]]))
  expect_lte(fit$diagnostics$edge.rel.rmse[[1L]], before$edge.rel.rmse[[1L]])
})

test_that("edge-KK weighted-GRIP initialization uses edge-only preparation", {
  edges <- rbind(c(1L, 2L), c(2L, 3L), c(3L, 4L), c(1L, 4L), c(1L, 3L))
  edge.weights <- c(1, 1.4, 1.1, 1.2, 1.8)
  prepared <- prepare.edge.kk(
    edges = edges,
    n = 4L,
    edge_weights = edge.weights
  )
  weighted.args <- list(rounds = 4L, final_rounds = 4L, num_init = 3L)
  init <- do.call(grip, c(
    list(
      edges = prepared$edges,
      n = prepared$n,
      edge_weights = prepared$edge_targets,
      dim = 3L,
      seed = 11L,
      metric = "edge_length"
    ),
    weighted.args
  ))
  manual <- edge.kk(
    coords = init,
    prepared = prepared,
    dim = 3L,
    stiffness_method = "uniform",
    density_mix_schedule = 1,
    scale_mode = "profiled",
    max_iter = 3L,
    diagnostics = FALSE,
    return_trace = FALSE,
    engine = "cpp"
  )
  direct <- edge.kk(
    prepared = prepared,
    dim = 3L,
    init = "weighted_grip",
    weighted.grip.args = weighted.args,
    stiffness_method = "uniform",
    density_mix_schedule = 1,
    scale_mode = "profiled",
    max_iter = 3L,
    diagnostics = FALSE,
    return_trace = FALSE,
    seed = 11L,
    engine = "cpp"
  )

  expect_s3_class(direct$prepared, "grip_edge_kk_prepared")
  expect_equal(direct$prepared$pair_mode, "edge_only")
  expect_equal(direct$coords, manual$coords, tolerance = 1e-12)
  expect_null(direct$trace)
  expect_null(direct$metadata$frames)
})

test_that("edge-KK raw graph input uses edge-only preparation when coordinates are supplied", {
  edges <- edges.mesh(4L, 4L, connectivity = "orthogonal")
  coords <- cbind(
    rep(seq_len(4L), times = 4L),
    rep(seq_len(4L), each = 4L)
  )
  fit <- edge.kk(
    coords = coords,
    edges = edges,
    n = nrow(coords),
    edge_weights = rep(1, nrow(edges)),
    stiffness_method = "uniform",
    density_mix_schedule = 1,
    max_iter = 1L,
    diagnostics = TRUE,
    return_trace = FALSE,
    engine = "cpp"
  )

  expect_s3_class(fit$prepared, "grip_edge_kk_prepared")
  expect_equal(fit$prepared$pair_mode, "edge_only")
  expect_null(fit$prepared$distance_matrix)
  expect_true(is.na(fit$diagnostics$gmds.stress[[1L]]))
  expect_null(fit$trace)
  expect_null(fit$metadata$frames)
  expect_s3_class(fit$metadata$stage_summaries, "data.frame")
  expect_equal(nrow(fit$coords), nrow(coords))
})

test_that("edge-KK random initialization uses edge-only preparation from raw graph input", {
  edges <- edges.path(6L)
  fit <- edge.kk(
    edges = edges,
    n = 6L,
    init = "random",
    stiffness_method = "uniform",
    density_mix_schedule = 1,
    max_iter = 1L,
    diagnostics = FALSE,
    return_trace = FALSE,
    seed = 5L,
    engine = "cpp"
  )

  expect_s3_class(fit$prepared, "grip_edge_kk_prepared")
  expect_equal(fit$prepared$pair_mode, "edge_only")
  expect_null(fit$prepared$distance_matrix)
  expect_null(fit$trace)
  expect_null(fit$metadata$frames)
  expect_equal(dim(fit$coords), c(6L, 2L))
})

test_that("edge-KK weighted-GRIP initialization uses edge-only raw graph preparation", {
  edges <- edges.path(6L)
  fit <- edge.kk(
    edges = edges,
    n = 6L,
    edge_weights = rep(1, nrow(edges)),
    init = "weighted_grip",
    weighted.grip.args = list(rounds = 3L, final_rounds = 3L, num_init = 3L),
    stiffness_method = "uniform",
    density_mix_schedule = 1,
    max_iter = 1L,
    diagnostics = FALSE,
    return_trace = FALSE,
    seed = 5L,
    engine = "cpp"
  )

  expect_s3_class(fit$prepared, "grip_edge_kk_prepared")
  expect_equal(fit$prepared$pair_mode, "edge_only")
  expect_null(fit$prepared$distance_matrix)
  expect_equal(dim(fit$coords), c(6L, 2L))
})

test_that("edge-KK weighted-GRIP initialization validates companion arguments", {
  prepared <- prepare.edge.kk(
    edges = edges.path(4L),
    n = 4L,
    edge_weights = rep(1, 3L)
  )

  expect_error(
    edge.kk(
      prepared = prepared,
      init = "weighted_grip",
      weighted.grip.args = list(edges = edges.path(4L)),
      max_iter = 1L
    ),
    "weighted.grip.args must not include"
  )
  expect_error(
    edge.kk(
      prepared = prepared,
      init = "weighted_grip",
      dim = 4L,
      max_iter = 1L
    ),
    "init = \"weighted_grip\" requires dim to be 2 or 3"
  )
})

test_that("edge-KK omits trace rows and frames when tracing is disabled in R engine", {
  edges <- edges.path(5L)
  start <- cbind(seq_len(5L), c(0.1, -0.2, 0.3, -0.1, 0.2))
  fit <- edge.kk(
    coords = start,
    edges = edges,
    n = 5L,
    stiffness_method = "uniform",
    density_mix_schedule = c(0, 1),
    max_iter = 2L,
    diagnostics = FALSE,
    return_trace = FALSE,
    engine = "R"
  )

  expect_null(fit$trace)
  expect_null(fit$metadata$frames)
  expect_s3_class(fit$metadata$stage_summaries, "data.frame")
  expect_equal(nrow(fit$metadata$stage_summaries), 2L)
})

test_that("edge-only preparation rejects duplicate undirected edges", {
  expect_error(
    prepare.edge.kk(
      edges = rbind(c(1L, 2L), c(2L, 1L)),
      n = 2L,
      edge_weights = c(1, 1)
    ),
    "duplicate undirected edges"
  )
})

test_that("edge-only preparation does not support classical-MDS initialization", {
  prepared <- prepare.edge.kk(
    edges = edges.path(4L),
    n = 4L,
    edge_weights = c(1, 1, 1)
  )

  expect_error(
    classical.mds(prepared = prepared),
    "classical MDS requires an all-pairs prepared object"
  )
  expect_error(
    edge.kk(
      prepared = prepared,
      dim = 2L,
      init = "classical_mds",
      max_iter = 1L
    ),
    "MDS initialization requires an all-pairs prepared object"
  )
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

test_that("edge-isometric energy gradient supports higher-dimensional coordinates", {
  edges <- rbind(c(1L, 2L), c(2L, 3L), c(3L, 4L), c(4L, 5L))
  coords <- matrix(c(
    0.0, 0.0, 0.2, -0.1, 0.3,
    1.1, 0.1, 0.0,  0.4, -0.2,
    1.8, 0.6, 0.3, -0.3, 0.5,
    2.7, 0.4, 0.8,  0.2, 0.1,
    3.5, 0.9, 0.4, -0.4, 0.6
  ), ncol = 5, byrow = TRUE)
  edge.weights <- c(1.2, 1.1, 1.4, 1.3)
  stiffness <- c(0.8, 1.2, 0.7, 1.5)
  state <- grip:::grip.edge.isometric.energy.gradient(
    coords = coords,
    edges = edges,
    edge_weights = edge.weights,
    stiffness = stiffness,
    scale = 1.05,
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
        plus, edges, edge.weights, stiffness, scale = 1.05
      )$energy
      e.minus <- grip:::grip.edge.isometric.energy.gradient(
        minus, edges, edge.weights, stiffness, scale = 1.05
      )$energy
      fd[i, j] <- (e.plus - e.minus) / (2 * eps)
    }
  }

  expect_equal(dim(state$gradient), c(5L, 5L))
  expect_equal(state$gradient, fd, tolerance = 1e-5)
})

test_that("edge-KK optimizer preserves exact weighted path layouts", {
  prepared <- prepare.graph.geodesic.mds(
    edges = edges.path(4L),
    n = 4L,
    edge_weights = c(1, 2, 1)
  )
  coords <- cbind(c(0, 1, 3, 4), 0)
  fit <- edge.kk(
    coords = coords,
    prepared = prepared,
    dim = 2L,
    stiffness_method = "uniform",
    density_mix_schedule = 1,
    scale_mode = "identity",
    max_iter = 5L
  )

  expect_s3_class(fit, "grip_gmds_layout")
  expect_equal(fit$method, "edge_kk")
  expect_lt(fit$diagnostics$edge.rel.rmse[[1L]], 1e-8)
  expect_lt(fit$diagnostics$gmds.stress[[1L]], 1e-8)
})

test_that("edge-KK optimizer preserves exact higher-dimensional weighted path layouts", {
  prepared <- prepare.graph.geodesic.mds(
    edges = edges.path(5L),
    n = 5L,
    edge_weights = c(1, 2, 1.5, 0.75)
  )
  coords <- cbind(c(0, 1, 3, 4.5, 5.25), matrix(0, nrow = 5L, ncol = 3L))
  fit <- edge.kk(
    coords = coords,
    prepared = prepared,
    dim = 4L,
    stiffness_method = "uniform",
    density_mix_schedule = 1,
    scale_mode = "identity",
    max_iter = 5L,
    engine = "cpp"
  )

  expect_s3_class(fit, "grip_gmds_layout")
  expect_equal(dim(fit$coords), c(5L, 4L))
  expect_equal(fit$method, "edge_kk")
  expect_lt(fit$diagnostics$edge.rel.rmse[[1L]], 1e-8)
  expect_lt(fit$diagnostics$gmds.stress[[1L]], 1e-8)
})

test_that("edge-KK optimizer decreases edge error from perturbed layout", {
  prepared <- prepare.graph.geodesic.mds(
    edges = edges.path(5L),
    n = 5L,
    edge_weights = rep(1, 4L)
  )
  start <- cbind(c(0, 0.7, 1.9, 2.4, 4.2), c(0, 0.4, -0.2, 0.5, -0.1))
  before <- score.gmds(
    start,
    prepared = prepared,
    scale_mode = "identity"
  )
  fit <- edge.kk(
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

test_that("C++ edge-KK optimizer matches R reference engine", {
  prepared <- prepare.graph.geodesic.mds(
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
  fit.cpp <- do.call(edge.kk, c(args, list(engine = "cpp")))
  fit.r <- do.call(edge.kk, c(args, list(engine = "R")))

  expect_equal(fit.cpp$coords, fit.r$coords, tolerance = 1e-10)
  expect_equal(fit.cpp$trace$energy, fit.r$trace$energy, tolerance = 1e-10)
  expect_equal(fit.cpp$trace$edge.rel.rmse, fit.r$trace$edge.rel.rmse, tolerance = 1e-10)
  expect_equal(fit.cpp$metadata$engine, "cpp_gradient_descent_armijo")
  expect_equal(fit.r$metadata$engine, "r_gradient_descent_armijo")
})

test_that("C++ edge-KK optimizer matches R reference engine in higher dimensions", {
  prepared <- prepare.graph.geodesic.mds(
    edges = rbind(
      c(1L, 2L), c(2L, 3L), c(3L, 4L),
      c(4L, 5L), c(1L, 5L), c(2L, 5L)
    ),
    n = 5L,
    edge_weights = c(1, 1.4, 1, 1.3, 2.1, 1.7)
  )
  start <- matrix(c(
    0.0, 0.0, 0.2, -0.1,
    0.8, 0.4, 0.0,  0.3,
    1.7, 1.0, 0.5, -0.2,
    2.4, 0.9, 1.0,  0.1,
    0.2, 1.3, 0.4,  0.7
  ), ncol = 4, byrow = TRUE)
  args <- list(
    coords = start,
    prepared = prepared,
    dim = 4L,
    stiffness_method = "density",
    stiffness_transform = "sqrt",
    density_mix_schedule = c(0, 0.5, 1),
    scale_mode = "profiled",
    max_iter = 8L,
    initial_step = 0.2,
    return_trace = TRUE,
    diagnostics = FALSE
  )
  fit.cpp <- do.call(edge.kk, c(args, list(engine = "cpp")))
  fit.r <- do.call(edge.kk, c(args, list(engine = "R")))

  expect_equal(dim(fit.cpp$coords), c(5L, 4L))
  expect_equal(fit.cpp$coords, fit.r$coords, tolerance = 1e-10)
  expect_equal(fit.cpp$trace$energy, fit.r$trace$energy, tolerance = 1e-10)
  expect_equal(fit.cpp$trace$edge.rel.rmse, fit.r$trace$edge.rel.rmse, tolerance = 1e-10)
})

test_that("C++ edge-KK optimizer supports fixed and user scale modes", {
  prepared <- prepare.graph.geodesic.mds(
    edges = edges.path(4L),
    n = 4L,
    edge_weights = c(1, 2, 1)
  )
  start <- 2 * cbind(c(0, 1.1, 3.2, 4.0), c(0, 0.2, -0.1, 0.1))
  fixed <- edge.kk(
    coords = start,
    prepared = prepared,
    dim = 2L,
    stiffness_method = "uniform",
    density_mix_schedule = c(0, 1),
    scale_mode = "fixed_initial",
    max_iter = 3L,
    engine = "cpp"
  )
  user <- edge.kk(
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

test_that("classical-MDS initialization supports higher-dimensional edge-KK layouts", {
  prepared <- prepare.graph.geodesic.mds(
    edges = edges.path(6L),
    n = 6L,
    edge_weights = c(1, 1.5, 0.75, 1.25, 1)
  )
  init <- suppressWarnings(
    classical.mds(
      prepared = prepared,
      dim = 4L,
      diagnostics = TRUE
    )
  )
  fit <- suppressWarnings(
    edge.kk(
      prepared = prepared,
      dim = 4L,
      init = "classical_mds",
      stiffness_method = "uniform",
      density_mix_schedule = 1,
      scale_mode = "identity",
      max_iter = 3L,
      diagnostics = FALSE,
      engine = "cpp"
    )
  )

  expect_s3_class(init, "grip_gmds_layout")
  expect_equal(dim(init$coords), c(6L, 4L))
  expect_equal(dim(fit$coords), c(6L, 4L))
  expect_true(all(is.finite(fit$coords)))
})
