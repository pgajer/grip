test_that("geodesic-MDS scoring exposes smoothness energy", {
  bundle <- mesh.surface.graph(
    4,
    4,
    surface = "paraboloid",
    amplitude = 0.2,
    connectivity = "orthogonal",
    normalize = "median"
  )
  prepared <- prepare.geodesic.kk(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average"
  )
  coords <- grip:::grip.classical.mds.embedding(prepared, dim = 3L)$coords

  plain <- grip.score.geodesic.mds(coords, prepared = prepared)
  smooth <- grip.score.geodesic.mds(coords, prepared = prepared, smoothness_weight = 0.1)

  expect_equal(plain$smooth.energy[[1L]], 0)
  expect_equal(plain$smooth.raw.penalty[[1L]], 0)
  expect_gt(smooth$smooth.energy[[1L]], 0)
  expect_gt(smooth$smooth.raw.penalty[[1L]], 0)
})

test_that("C++ and R smoothness-regularized GMDS agree for one iteration", {
  bundle <- mesh.surface.graph(
    4,
    4,
    surface = "paraboloid",
    amplitude = 0.2,
    connectivity = "orthogonal",
    normalize = "median"
  )
  prepared <- prepare.geodesic.kk(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average"
  )
  coords <- grip:::grip.classical.mds.embedding(prepared, dim = 3L)$coords

  fit.r <- grip.optimize.geodesic.mds(
    coords = coords,
    prepared = prepared,
    engine = "r",
    max_iter = 1L,
    initial_step = 0.05,
    smoothness_weight = 0.1,
    return_trace = TRUE,
    recenter = TRUE
  )
  fit.cpp <- grip.optimize.geodesic.mds(
    coords = coords,
    prepared = prepared,
    engine = "cpp",
    max_iter = 1L,
    initial_step = 0.05,
    smoothness_weight = 0.1,
    return_trace = TRUE,
    recenter = TRUE,
    n_threads = 1L
  )

  expect_true("smooth_energy" %in% names(fit.cpp$trace))
  expect_true("smooth_weight" %in% names(fit.cpp$trace))
  expect_equal(fit.r$trace$smooth_weight[[1L]], 0.1)
  expect_equal(fit.cpp$trace$smooth_weight[[1L]], 0.1)
  expect_equal(unname(fit.cpp$coords), unname(fit.r$coords), tolerance = 1e-7)
  expect_equal(
    fit.cpp$score$smooth.energy[[1L]],
    fit.r$score$smooth.energy[[1L]],
    tolerance = 1e-8
  )
})
