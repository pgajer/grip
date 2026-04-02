test_that("weighted trace requires edge weights", {
  edges <- edges.mesh(4, 4)
  expect_error(
    grip.layout.trace.weighted(edges = edges, n = 16, dim = 2, seed = 1),
    "requires edge weights"
  )
})

test_that("weighted round trace returns frames, metadata, and inactive NA rows", {
  graph <- mesh.surface.graph(4, 5, surface = "saddle", amplitude = 0.6)
  tr <- grip.layout.trace.weighted(
    edges = graph$edges,
    edge_weights = graph$edge_weights,
    n = graph$n,
    dim = 2,
    rounds = 6,
    final_rounds = 4,
    num_init = 5,
    num_nbrs = 6,
    coarse_repulsion_factor = 0.3,
    coarse_repulsion_sample = 8,
    coarse_repulsion_exact_below = 32,
    trace = "round",
    trace.every = 1,
    seed = 123
  )

  expect_s3_class(tr, "grip_layout_trace")
  expect_equal(dim(tr$final), c(graph$n, 2))
  expect_true(all(is.finite(tr$final)))
  expect_true(length(tr$frames) >= 3L)
  expect_equal(names(tr$meta),
               c("frame", "phase", "level_index", "misf_level", "round_in_level", "active_vertices"))
  expect_identical(tr$meta$phase[[1]], "init")
  expect_identical(tail(tr$meta$phase, 1), "final")
  expect_true(any(tr$meta$phase == "round"))
  expect_true(all(diff(tr$meta$active_vertices) >= 0))

  first_frame <- tr$frames[[1L]]
  expect_equal(dim(first_frame), c(graph$n, 2))
  expect_true(any(rowSums(is.na(first_frame)) == 2L))
  expect_equal(sum(rowSums(is.na(first_frame)) == 2L),
               graph$n - tr$meta$active_vertices[[1L]])
  expect_equal(tr$frames[[length(tr$frames)]], tr$final)
})

test_that("weighted trace final matches weighted layout", {
  graph <- mesh.surface.graph(5, 5, surface = "ripple", amplitude = 0.5)
  tr <- grip.layout.trace.weighted(
    edges = graph$edges,
    edge_weights = graph$edge_weights,
    n = graph$n,
    dim = 2,
    rounds = 8,
    final_rounds = 8,
    num_init = 6,
    num_nbrs = 8,
    coarse_repulsion_factor = 0.3,
    coarse_repulsion_sample = 8,
    coarse_repulsion_exact_below = 32,
    trace = "level",
    trace.every = 1,
    seed = 41
  )
  coords <- grip.layout.weighted(
    edges = graph$edges,
    edge_weights = graph$edge_weights,
    n = graph$n,
    dim = 2,
    rounds = 8,
    final_rounds = 8,
    num_init = 6,
    num_nbrs = 8,
    coarse_repulsion_factor = 0.3,
    coarse_repulsion_sample = 8,
    coarse_repulsion_exact_below = 32,
    seed = 41
  )
  expect_identical(tr$final, coords)
})

test_that("weighted trace final matches weighted layout for LGKK polish", {
  graph <- mesh.surface.graph(5, 5, surface = "paraboloid", amplitude = 0.4)
  tr <- grip.layout.trace.weighted(
    edges = graph$edges,
    edge_weights = graph$edge_weights,
    n = graph$n,
    dim = 2,
    rounds = 8,
    final_rounds = 8,
    num_init = 6,
    num_nbrs = 8,
    coarse_repulsion_factor = 0.3,
    coarse_repulsion_sample = 8,
    coarse_repulsion_exact_below = 32,
    lgkk_polish_rounds = 4,
    lgkk_local_nbrs = 6,
    lgkk_landmark_count = 4,
    trace = "level",
    trace.every = 1,
    seed = 47
  )
  coords <- grip.layout.weighted(
    edges = graph$edges,
    edge_weights = graph$edge_weights,
    n = graph$n,
    dim = 2,
    rounds = 8,
    final_rounds = 8,
    num_init = 6,
    num_nbrs = 8,
    coarse_repulsion_factor = 0.3,
    coarse_repulsion_sample = 8,
    coarse_repulsion_exact_below = 32,
    lgkk_polish_rounds = 4,
    lgkk_local_nbrs = 6,
    lgkk_landmark_count = 4,
    seed = 47
  )
  expect_identical(tr$final, coords)
  expect_true(nrow(tr$lgkk.polish) >= 1L)
})

test_that("weighted trace diagnostics are produced", {
  graph <- mesh.surface.graph(4, 5, surface = "saddle", amplitude = 0.6)
  tr <- grip.layout.trace.weighted(
    edges = graph$edges,
    edge_weights = graph$edge_weights,
    n = graph$n,
    dim = 2,
    rounds = 6,
    final_rounds = 4,
    num_init = 5,
    num_nbrs = 6,
    coarse_repulsion_factor = 0.3,
    coarse_repulsion_sample = 8,
    coarse_repulsion_exact_below = 32,
    trace = "level",
    trace.every = 1,
    diagnostics = "light",
    seed = 131
  )

  expect_equal(nrow(tr$diagnostics), nrow(tr$meta))
  expect_true(all(c("edge.length.cv", "sampled.nonedge.sep.ratio", "sampled.stress") %in%
                    names(tr$diagnostics)))
})
