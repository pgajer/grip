test_that("GMDS layout result constructor preserves common shape", {
  prepared <- prepare.graph.geodesic.mds(
    edges = edges.path(4L),
    n = 4L,
    edge_weights = c(1, 2, 1)
  )
  coords <- cbind(c(0, 1, 3, 4), 0)
  diagnostics <- score.gmds(coords, prepared = prepared)
  layout <- grip:::gmds.result(
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
  prepared <- prepare.graph.geodesic.mds(
    edges = edges.path(4L),
    n = 4L,
    edge_weights = c(1, 2, 1)
  )
  coords <- cbind(c(0, 1, 3, 4), 0)
  diagnostics <- score.gmds(
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

test_that("flat path-length diagnostics match the R fallback", {
  prepared <- prepare.graph.geodesic.mds(
    edges = edges.mesh(4L, 4L),
    n = 16L
  )
  coords <- cbind(
    seq_len(16L),
    rev(seq_len(16L)),
    sqrt(seq_len(16L)),
    sin(seq_len(16L))
  )
  fast <- grip:::grip.geodesic.kk.path.lengths(coords, prepared)
  fallback <- grip:::grip.geodesic.kk.path.lengths.r(coords, prepared)

  expect_equal(fast, fallback, tolerance = 1e-12)
})

test_that("flat path-length diagnostics handle tie-averaged path coefficients", {
  edges <- rbind(c(1L, 2L), c(2L, 4L), c(1L, 3L), c(3L, 4L))
  prepared <- prepare.graph.geodesic.mds(
    edges = edges,
    n = 4L,
    tie_mode = "average"
  )
  coords <- rbind(
    c(0, 0),
    c(1, 0),
    c(0, 2),
    c(2, 2)
  )
  fast <- grip:::grip.geodesic.kk.path.lengths(coords, prepared)
  fallback <- grip:::grip.geodesic.kk.path.lengths.r(coords, prepared)

  expect_equal(fast, fallback, tolerance = 1e-12)
})

test_that("classical MDS baseline returns common GMDS layout result", {
  prepared <- prepare.graph.geodesic.mds(
    edges = edges.path(5L),
    n = 5L,
    edge_weights = rep(1, 4L)
  )
  layout <- classical.mds(prepared = prepared, dim = 2L)

  expect_s3_class(layout, "grip_gmds_layout")
  expect_equal(layout$method, "classical_mds")
  expect_equal(dim(layout$coords), c(5L, 2L))
  expect_true(all(is.finite(layout$coords)))
  expect_true(is.data.frame(layout$diagnostics))
  expect_lt(layout$diagnostics$edge.rel.rmse[[1L]], 1e-8)
  expect_lt(layout$diagnostics$gmds.stress[[1L]], 1e-8)
  expect_equal(layout$metadata$engine, "cmdscale")
  expect_true("positive_eigen_fraction" %in% names(layout$metadata))
})

test_that("classical MDS skips full path cache for cold no-diagnostic calls", {
  edges <- edges.mesh(6L, 6L)
  prepared <- prepare.graph.geodesic.mds(edges = edges, n = 36L)
  full <- classical.mds(prepared = prepared, dim = 2L, diagnostics = FALSE)
  fast <- classical.mds(edges = edges, n = 36L, dim = 2L, diagnostics = FALSE)

  expect_s3_class(fast, "grip_gmds_layout")
  expect_s3_class(fast$prepared, "grip_metric_mds_prepared")
  expect_equal(fast$prepared$pair_mode, "distance_matrix_only")
  expect_equal(length(fast$prepared$path_edges), 0L)
  expect_equal(fast$coords, full$coords, tolerance = 1e-10)
  expect_null(fast$diagnostics)
})

test_that("classical MDS keeps full path cache when diagnostics are requested", {
  edges <- edges.mesh(4L, 4L)
  layout <- classical.mds(edges = edges, n = 16L, dim = 2L, diagnostics = TRUE)

  expect_s3_class(layout$prepared, "grip_gmds_prepared")
  expect_false(inherits(layout$prepared, "grip_metric_mds_prepared"))
  expect_equal(layout$prepared$pair_mode, "all_pairs")
  expect_gt(length(layout$prepared$path_edges), 0L)
  expect_s3_class(layout$diagnostics, "data.frame")
})

test_that("profiled diagnostics factor out global scale", {
  prepared <- prepare.graph.geodesic.mds(
    edges = edges.path(4L),
    n = 4L,
    edge_weights = c(1, 2, 1)
  )
  coords <- 3 * cbind(c(0, 1, 3, 4), 0)
  profiled <- score.gmds(
    coords = coords,
    prepared = prepared,
    scale_mode = "profiled"
  )
  identity <- score.gmds(
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
