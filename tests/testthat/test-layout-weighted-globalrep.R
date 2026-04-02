test_that("weighted GRIP requires edge weights", {
  edges <- edges.mesh(4, 4)
  expect_error(
    grip.layout.globalrep.weighted(edges = edges, n = 16, dim = 2, seed = 1),
    "requires edge weights"
  )
  expect_error(
    grip.build.misf.weighted(edges = edges, n = 16, seed = 1),
    "requires edge weights"
  )
})

test_that("weighted globalrep layout returns a finite deterministic matrix", {
  graph <- mesh.surface.graph(4, 5, surface = "saddle", amplitude = 0.6)
  coords1 <- grip.layout.globalrep.weighted(
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
  coords2 <- grip.layout.globalrep.weighted(
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

  expect_equal(dim(coords1), c(graph$n, 2))
  expect_true(all(is.finite(coords1)))
  expect_identical(coords1, coords2)
})

test_that("grip.layout.weighted is an alias of weighted globalrep", {
  graph <- mesh.surface.graph(4, 5, surface = "ripple", amplitude = 0.5)
  coords_primary <- grip.layout.weighted(
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
    seed = 17
  )
  coords_alias <- grip.layout.globalrep.weighted(
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
    seed = 17
  )
  expect_identical(coords_primary, coords_alias)
})

test_that("weighted layout is invariant to global weight rescaling under normalization", {
  graph <- mesh.surface.graph(4, 5, surface = "paraboloid", amplitude = 0.4)
  coords_base <- grip.layout.globalrep.weighted(
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
    length_normalization = "median",
    seed = 19
  )
  coords_scaled <- grip.layout.globalrep.weighted(
    edges = graph$edges,
    edge_weights = graph$edge_weights * 7,
    n = graph$n,
    dim = 2,
    rounds = 8,
    final_rounds = 8,
    num_init = 6,
    num_nbrs = 8,
    coarse_repulsion_factor = 0.3,
    coarse_repulsion_sample = 8,
    coarse_repulsion_exact_below = 32,
    length_normalization = "median",
    seed = 19
  )

  expect_equal(coords_base, coords_scaled, tolerance = 1e-10)
})

test_that("weighted GRIP responds to nontrivial geometry differently than combinatorial GRIP", {
  graph <- mesh.surface.graph(5, 5, surface = "saddle", amplitude = 0.9)
  coords_comb <- grip.layout.globalrep(
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
    seed = 23
  )
  coords_weighted <- grip.layout.globalrep.weighted(
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
    seed = 23
  )

  expect_gt(max(abs(coords_comb - coords_weighted)), 1e-6)
})

test_that("weighted MISF helper is deterministic and normalization-aware", {
  graph <- mesh.surface.graph(5, 5, surface = "ripple", amplitude = 0.6)
  misf1 <- grip.build.misf.weighted(
    edges = graph$edges,
    edge_weights = graph$edge_weights,
    n = graph$n,
    num_init = 6,
    num_nbrs = 8,
    seed = 29
  )
  misf2 <- grip.build.misf.weighted(
    edges = graph$edges,
    edge_weights = graph$edge_weights,
    n = graph$n,
    num_init = 6,
    num_nbrs = 8,
    seed = 29
  )
  misf_scaled <- grip.build.misf.weighted(
    edges = graph$edges,
    edge_weights = graph$edge_weights * 11,
    n = graph$n,
    num_init = 6,
    num_nbrs = 8,
    seed = 29
  )

  expect_identical(misf1$levels, misf2$levels)
  expect_identical(misf1$mish_order, misf2$mish_order)
  expect_identical(misf1$vertex_depth, misf2$vertex_depth)
  expect_identical(misf1$levels, misf_scaled$levels)
  expect_true(is.finite(misf1$weight_scale))
  expect_identical(misf1$length_normalization, "median")
})
