test_that("GMDS layout result constructor preserves common shape", {
  prepared <- grip.prepare.graph.geodesic.mds(
    edges = edges.path(4L),
    n = 4L,
    edge_weights = c(1, 2, 1)
  )
  coords <- cbind(c(0, 1, 3, 4), 0)
  diagnostics <- grip.score.gmds.layout(coords, prepared = prepared)
  layout <- grip.gmds.layout.result(
    coords = coords,
    method = "fixture",
    prepared = prepared,
    diagnostics = diagnostics,
    metadata = list(source = "unit-test")
  )

  expect_s3_class(layout, "grip_gmds_layout")
  expect_equal(layout$method, "fixture")
  expect_equal(layout$coords, coords)
  expect_equal(layout$metadata$source, "unit-test")
  expect_true(all(c(
    "edge.rel.rmse",
    "gmds.stress",
    "gmds.short.stress",
    "gmds.mid.stress",
    "gmds.long.stress",
    "metric.chord.stress",
    "spread.score"
  ) %in% names(layout$diagnostics)))
})

test_that("common diagnostics are exact on a weighted path realization", {
  prepared <- grip.prepare.graph.geodesic.mds(
    edges = edges.path(4L),
    n = 4L,
    edge_weights = c(1, 2, 1)
  )
  coords <- cbind(c(0, 1, 3, 4), 0)
  diagnostics <- grip.score.gmds.layout(
    coords = coords,
    prepared = prepared,
    scale_mode = "identity"
  )

  expect_lt(diagnostics$edge.rel.rmse[[1L]], 1e-8)
  expect_lt(diagnostics$gmds.stress[[1L]], 1e-8)
  expect_lt(diagnostics$gmds.short.stress[[1L]], 1e-8)
  expect_lt(diagnostics$gmds.mid.stress[[1L]], 1e-8)
  expect_lt(diagnostics$gmds.long.stress[[1L]], 1e-8)
  expect_lt(diagnostics$metric.chord.stress[[1L]], 1e-8)
  expect_equal(diagnostics$shortcut.fraction[[1L]], 0)
})

test_that("metric MDS baseline returns common GMDS layout result", {
  prepared <- grip.prepare.graph.geodesic.mds(
    edges = edges.path(5L),
    n = 5L,
    edge_weights = rep(1, 4L)
  )
  layout <- grip.metric.mds.layout(prepared = prepared, dim = 2L)

  expect_s3_class(layout, "grip_gmds_layout")
  expect_equal(layout$method, "metric_mds")
  expect_equal(dim(layout$coords), c(5L, 2L))
  expect_true(all(is.finite(layout$coords)))
  expect_true(is.data.frame(layout$diagnostics))
  expect_lt(layout$diagnostics$edge.rel.rmse[[1L]], 1e-8)
  expect_lt(layout$diagnostics$gmds.stress[[1L]], 1e-8)
  expect_equal(layout$metadata$engine, "cmdscale")
  expect_true("positive_eigen_fraction" %in% names(layout$metadata))
})

test_that("profiled diagnostics factor out global scale", {
  prepared <- grip.prepare.graph.geodesic.mds(
    edges = edges.path(4L),
    n = 4L,
    edge_weights = c(1, 2, 1)
  )
  coords <- 3 * cbind(c(0, 1, 3, 4), 0)
  profiled <- grip.score.gmds.layout(
    coords = coords,
    prepared = prepared,
    scale_mode = "profiled"
  )
  identity <- grip.score.gmds.layout(
    coords = coords,
    prepared = prepared,
    scale_mode = "identity"
  )

  expect_lt(profiled$edge.rel.rmse[[1L]], 1e-8)
  expect_lt(profiled$gmds.stress[[1L]], 1e-8)
  expect_gt(identity$edge.rel.rmse[[1L]], 1)
  expect_gt(identity$gmds.stress[[1L]], 1)
  expect_equal(profiled$edge.scale[[1L]], 3, tolerance = 1e-10)
  expect_equal(profiled$gmds.scale[[1L]], 3, tolerance = 1e-10)
})
