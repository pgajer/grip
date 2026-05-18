test_that("weighted ND layout returns finite deterministic coordinates above 3D", {
  edges <- edges.mesh(3, 4)
  weights <- rep(1, nrow(edges))

  coords1 <- grip.layout.weighted.nd(
    edges = edges,
    edge_weights = weights,
    n = 12,
    dim = 5,
    rounds = 8,
    final_rounds = 8,
    num_init = 6,
    num_nbrs = 8,
    repulsion_factor = 0.4,
    seed = 71
  )
  coords2 <- grip.layout.weighted.nd(
    edges = edges,
    edge_weights = weights,
    n = 12,
    dim = 5,
    rounds = 8,
    final_rounds = 8,
    num_init = 6,
    num_nbrs = 8,
    repulsion_factor = 0.4,
    seed = 71
  )

  expect_equal(dim(coords1), c(12, 5))
  expect_identical(colnames(coords1), paste0("Dim", 1:5))
  expect_true(all(is.finite(coords1)))
  expect_identical(coords1, coords2)
})

test_that("weighted ND layout remains available in 2D and 3D without using legacy validators", {
  edges <- edges.cycle(8)
  weights <- seq_len(nrow(edges)) / nrow(edges) + 1

  coords2 <- grip.layout.weighted.nd(
    edges = edges,
    edge_weights = weights,
    n = 8,
    dim = 2,
    rounds = 4,
    final_rounds = 4,
    num_init = 4,
    seed = 73
  )
  coords3 <- grip.layout.weighted.nd(
    edges = edges,
    edge_weights = weights,
    n = 8,
    dim = 3,
    rounds = 4,
    final_rounds = 4,
    num_init = 5,
    seed = 73
  )

  expect_equal(dim(coords2), c(8, 2))
  expect_equal(dim(coords3), c(8, 3))
  expect_true(all(is.finite(coords2)))
  expect_true(all(is.finite(coords3)))
})

test_that("weighted ND layout responds to nonuniform edge lengths", {
  edges <- edges.path(7)
  weights_flat <- rep(1, nrow(edges))
  weights_shaped <- c(1, 1, 4, 4, 1, 1)

  coords_flat <- grip.layout.weighted.nd(
    edges = edges,
    edge_weights = weights_flat,
    n = 7,
    dim = 4,
    rounds = 12,
    final_rounds = 12,
    num_init = 5,
    repulsion_factor = 0,
    seed = 79
  )
  coords_shaped <- grip.layout.weighted.nd(
    edges = edges,
    edge_weights = weights_shaped,
    n = 7,
    dim = 4,
    rounds = 12,
    final_rounds = 12,
    num_init = 5,
    repulsion_factor = 0,
    seed = 79
  )

  expect_gt(max(abs(coords_flat - coords_shaped)), 1e-6)
})

test_that("weighted ND layout packs disconnected components", {
  edges <- matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE)
  coords <- grip.layout.weighted.nd(
    edges = edges,
    edge_weights = c(1, 2),
    n = 4,
    dim = 4,
    rounds = 4,
    final_rounds = 4,
    num_init = 5,
    disconnected = "components",
    seed = 83
  )

  expect_equal(dim(coords), c(4, 4))
  expect_true(all(is.finite(coords)))
  expect_error(
    grip.layout.weighted.nd(
      edges = edges,
      edge_weights = c(1, 2),
      n = 4,
      dim = 4,
      rounds = 1,
      final_rounds = 1,
      num_init = 5,
      disconnected = "error",
      seed = 83
    ),
    "graph must be connected"
  )
})

test_that("weighted ND layout keeps weighted presets usable with ND constraints", {
  edges <- edges.kary.tree(k = 2, depth = 3)
  coords <- grip.layout.weighted.nd(
    edges = edges,
    edge_weights = rep(1, nrow(edges)),
    n = max(edges),
    dim = 4,
    preset = "tree",
    rounds = 4,
    final_rounds = 4,
    seed = 87
  )

  expect_equal(dim(coords), c(max(edges), 4))
  expect_true(all(is.finite(coords)))
  expect_error(
    grip.layout.weighted.nd(
      edges = edges,
      edge_weights = rep(1, nrow(edges)),
      n = max(edges),
      dim = 4,
      placement = "circle",
      rounds = 1,
      final_rounds = 1,
      seed = 87
    ),
    "placement = 'barycenter'"
  )
})

test_that("weighted ND layout validates dimensions and required weights", {
  edges <- edges.path(5)

  expect_error(
    grip.layout.weighted.nd(edges = edges, n = 5, dim = 4),
    "requires edge weights"
  )
  expect_error(
    grip.layout.weighted.nd(
      edges = edges,
      edge_weights = rep(1, nrow(edges)),
      n = 5,
      dim = 1
    ),
    "dim must be >= 2"
  )
  expect_error(
    grip.layout.weighted.nd(
      edges = edges,
      edge_weights = rep(1, nrow(edges)),
      n = 5,
      dim = 3.5
    ),
    "dim must be an integer"
  )
  expect_error(
    grip.layout.weighted.nd(
      edges = edges,
      edge_weights = rep(1, nrow(edges)),
      n = 5,
      dim = 6,
      num_init = 6
    ),
    "num_init must be >= 7"
  )
})

test_that("legacy weighted GRIP dimensionality remains capped at 3D", {
  edges <- edges.path(6)
  expect_error(
    grip.layout.globalrep.weighted(
      edges = edges,
      edge_weights = rep(1, nrow(edges)),
      n = 6,
      dim = 4,
      rounds = 2,
      final_rounds = 2,
      num_init = 5,
      seed = 89
    ),
    "dim must be 2 or 3"
  )
})
