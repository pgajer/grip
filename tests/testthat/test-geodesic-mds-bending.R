test_that("geodesic-MDS scoring exposes bending energy", {
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
  bending <- grip.rectangular.grid.bending.stencils(4L, 4L)

  plain <- grip.score.geodesic.mds(coords, prepared = prepared)
  bent <- grip.score.geodesic.mds(
    coords,
    prepared = prepared,
    bending_stencils = bending,
    bending_weight = 0.1
  )

  expect_equal(plain$bend.energy[[1L]], 0)
  expect_equal(plain$bend.raw.penalty[[1L]], 0)
  expect_gt(bent$bend.energy[[1L]], 0)
  expect_gt(bent$bend.raw.penalty[[1L]], 0)
})

test_that("C++ and R bending-regularized GMDS agree for one iteration", {
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
  bending <- grip.rectangular.grid.bending.stencils(4L, 4L)

  fit.r <- grip.optimize.geodesic.mds(
    coords = coords,
    prepared = prepared,
    engine = "r",
    max_iter = 1L,
    initial_step = 0.05,
    bending_stencils = bending,
    bending_weight = 0.1,
    return_trace = TRUE,
    recenter = TRUE
  )
  fit.cpp <- grip.optimize.geodesic.mds(
    coords = coords,
    prepared = prepared,
    engine = "cpp",
    max_iter = 1L,
    initial_step = 0.05,
    bending_stencils = bending,
    bending_weight = 0.1,
    return_trace = TRUE,
    recenter = TRUE,
    n_threads = 1L
  )

  expect_true("bend_energy" %in% names(fit.cpp$trace))
  expect_true("bend_weight" %in% names(fit.cpp$trace))
  expect_equal(fit.r$trace$bend_weight[[1L]], 0.1)
  expect_equal(fit.cpp$trace$bend_weight[[1L]], 0.1)
  expect_equal(unname(fit.cpp$coords), unname(fit.r$coords), tolerance = 1e-7)
  expect_equal(
    fit.cpp$score$bend.energy[[1L]],
    fit.r$score$bend.energy[[1L]],
    tolerance = 1e-8
  )
})
