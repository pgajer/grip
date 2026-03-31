phase1.run.gmds.mesh.case <- function(bundle,
                                      truth_coords,
                                      dim = 2L,
                                      max_iter = 8L) {
  prepared <- grip.prepare.geodesic.kk(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights
  )

  cmd <- grip:::grip.classical.mds.embedding(prepared, dim = dim, eig = TRUE)
  before <- grip.score.geodesic.mds(cmd$coords, prepared = prepared)
  opt <- grip.optimize.geodesic.mds(
    coords = cmd$coords,
    prepared = prepared,
    max_iter = max_iter,
    engine = "cpp",
    return_trace = TRUE
  )
  fit.before <- grip:::grip.align.to.target.nd(
    cmd$coords,
    truth_coords,
    allow.reflection = TRUE
  )
  fit.after <- grip:::grip.align.to.target.nd(
    opt$coords,
    truth_coords,
    allow.reflection = TRUE
  )

  list(
    prepared = prepared,
    cmd = cmd,
    before = before,
    opt = opt,
    fit_before = fit.before,
    fit_after = fit.after
  )
}

test_that("Phase 1 classical MDS diagnostics are finite on a flat mesh", {
  bundle <- mesh.surface.graph(
    6, 6,
    surface = "saddle",
    amplitude = 0,
    normalize = "median"
  )
  prepared <- grip.prepare.geodesic.kk(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights
  )
  fit <- grip:::grip.classical.mds.embedding(prepared, dim = 2L, eig = TRUE)
  score <- grip:::grip.classical.mds.score.stats(fit$coords, prepared)

  expect_equal(dim(fit$coords), c(bundle$n, 2L))
  expect_true(all(is.finite(fit$coords)))
  expect_true(is.finite(score$raw_stress))
  expect_true(is.finite(score$stress))
  expect_true(is.finite(fit$positive_eigen_fraction))
  expect_true(is.finite(fit$negative_eigen_fraction))
  expect_equal(
    fit$positive_eigen_fraction + fit$negative_eigen_fraction,
    1,
    tolerance = 1e-10
  )
})

test_that("Phase 1 GMDS refinement sharply lowers path stress on flat meshes", {
  bundle <- mesh.surface.graph(
    6, 6,
    surface = "saddle",
    amplitude = 0,
    normalize = "median"
  )
  run <- phase1.run.gmds.mesh.case(
    bundle = bundle,
    truth_coords = bundle$coords_param,
    dim = 2L
  )

  expect_lt(
    run$opt$score$gmds.stress[[1L]],
    run$before$gmds.stress[[1L]] * 0.2
  )
  expect_true(nrow(run$opt$trace) >= 2L)
  expect_true(length(run$opt$frames) >= 2L)
  expect_true(run$fit_after$rmse > run$fit_before$rmse)
})

test_that("Phase 1 GMDS refinement lowers path stress on slit-channel meshes", {
  bundle <- occupied.mesh.surface.graph(
    keep.slit.channels(
      10, 10,
      orientation = "vertical",
      slit_period = 4,
      slit_width = 1,
      bridge_spacing = 3,
      bridge_size = 1,
      offset = 2
    ),
    surface = "saddle",
    amplitude = 0,
    normalize = "median"
  )
  run <- phase1.run.gmds.mesh.case(
    bundle = bundle,
    truth_coords = bundle$coords_param,
    dim = 2L
  )

  expect_lt(
    run$opt$score$gmds.stress[[1L]],
    run$before$gmds.stress[[1L]] * 0.3
  )
  expect_lte(run$fit_after$rmse, run$fit_before$rmse * 1.1)
  expect_true(all(is.finite(run$opt$coords)))
})

test_that("Phase 1 GMDS refinement lowers path stress for 3D ripple meshes", {
  bundle <- mesh.surface.graph(
    6, 6,
    surface = "ripple",
    amplitude = 0.5,
    freq_u = 2,
    freq_v = 2,
    normalize = "median"
  )
  run <- phase1.run.gmds.mesh.case(
    bundle = bundle,
    truth_coords = bundle$coords_surface,
    dim = 3L
  )

  expect_equal(dim(run$opt$coords), c(bundle$n, 3L))
  expect_lt(
    run$opt$score$gmds.stress[[1L]],
    run$before$gmds.stress[[1L]] * 0.2
  )
  expect_true(is.finite(run$fit_after$rmse))
})

test_that("diagonal mesh connectivity improves GMDS geometry on simple flat meshes", {
  orth <- mesh.surface.graph(
    6, 6,
    surface = "saddle",
    amplitude = 0,
    connectivity = "orthogonal",
    normalize = "median"
  )
  diag <- mesh.surface.graph(
    6, 6,
    surface = "saddle",
    amplitude = 0,
    connectivity = "diagonal",
    normalize = "median"
  )

  orth.run <- phase1.run.gmds.mesh.case(
    bundle = orth,
    truth_coords = orth$coords_param,
    dim = 2L
  )
  diag.run <- phase1.run.gmds.mesh.case(
    bundle = diag,
    truth_coords = diag$coords_param,
    dim = 2L
  )

  expect_lt(diag.run$fit_after$rmse, orth.run$fit_after$rmse)
})

test_that("diagonal mesh connectivity improves GMDS geometry on simple ripple meshes", {
  orth <- mesh.surface.graph(
    6, 6,
    surface = "ripple",
    amplitude = 0.5,
    freq_u = 2,
    freq_v = 2,
    connectivity = "orthogonal",
    normalize = "median"
  )
  diag <- mesh.surface.graph(
    6, 6,
    surface = "ripple",
    amplitude = 0.5,
    freq_u = 2,
    freq_v = 2,
    connectivity = "diagonal",
    normalize = "median"
  )

  orth.run <- phase1.run.gmds.mesh.case(
    bundle = orth,
    truth_coords = orth$coords_param,
    dim = 2L
  )
  diag.run <- phase1.run.gmds.mesh.case(
    bundle = diag,
    truth_coords = diag$coords_param,
    dim = 2L
  )

  expect_lt(diag.run$fit_after$rmse, orth.run$fit_after$rmse)
})
