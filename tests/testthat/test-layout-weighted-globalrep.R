test_that("weighted GRIP requires edge weights", {
  edges <- edges.mesh(4, 4)
  expect_error(
    globalrep.weighted.grip(edges = edges, n = 16, dim = 2, seed = 1),
    "requires edge weights"
  )
  expect_error(
    build.weighted.misf(edges = edges, n = 16, seed = 1),
    "requires edge weights"
  )
})

test_that("weighted globalrep layout returns a finite deterministic matrix", {
  graph <- mesh.surface.graph(4, 5, surface = "saddle", amplitude = 0.6)
  coords1 <- globalrep.weighted.grip(
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 2,
    rounds = 8,
    final.rounds = 8,
    num.init = 6,
    num.nbrs = 8,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    seed = 41
  )
  coords2 <- globalrep.weighted.grip(
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 2,
    rounds = 8,
    final.rounds = 8,
    num.init = 6,
    num.nbrs = 8,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    seed = 41
  )

  expect_equal(dim(coords1), c(graph$n, 2))
  expect_true(all(is.finite(coords1)))
  expect_identical(coords1, coords2)
})

test_that("grip edge-length mode matches the weighted backend", {
  graph <- mesh.surface.graph(4, 5, surface = "ripple", amplitude = 0.5)
  coords_primary <- grip(metric = "edge_length",
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 2,
    rounds = 8,
    final.rounds = 8,
    num.init = 6,
    num.nbrs = 8,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    seed = 17
  )
  coords_alias <- globalrep.weighted.grip(
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 2,
    rounds = 8,
    final.rounds = 8,
    num.init = 6,
    num.nbrs = 8,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    seed = 17
  )
  expect_identical(coords_primary, coords_alias)
})

test_that("weighted layout is invariant to global weight rescaling under normalization", {
  graph <- mesh.surface.graph(4, 5, surface = "paraboloid", amplitude = 0.4)
  coords_base <- globalrep.weighted.grip(
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 2,
    rounds = 8,
    final.rounds = 8,
    num.init = 6,
    num.nbrs = 8,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    length.normalization = "median",
    seed = 19
  )
  coords_scaled <- globalrep.weighted.grip(
    edges = graph$edges,
    edge.weights = graph$edge_weights * 7,
    n = graph$n,
    dim = 2,
    rounds = 8,
    final.rounds = 8,
    num.init = 6,
    num.nbrs = 8,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    length.normalization = "median",
    seed = 19
  )

  expect_equal(coords_base, coords_scaled, tolerance = 1e-10)
})

test_that("weighted GRIP responds to nontrivial geometry differently than combinatorial GRIP", {
  graph <- mesh.surface.graph(5, 5, surface = "saddle", amplitude = 0.9)
  coords_comb <- globalrep.grip(
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 2,
    rounds = 8,
    final.rounds = 8,
    num.init = 6,
    num.nbrs = 8,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    seed = 23
  )
  coords_weighted <- globalrep.weighted.grip(
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 2,
    rounds = 8,
    final.rounds = 8,
    num.init = 6,
    num.nbrs = 8,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    seed = 23
  )

  expect_gt(max(abs(coords_comb - coords_weighted)), 1e-6)
})

test_that("weighted MISF helper is deterministic and normalization-aware", {
  graph <- mesh.surface.graph(5, 5, surface = "ripple", amplitude = 0.6)
  misf1 <- build.weighted.misf(
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    num.init = 6,
    num.nbrs = 8,
    seed = 29
  )
  misf2 <- build.weighted.misf(
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    num.init = 6,
    num.nbrs = 8,
    seed = 29
  )
  misf_scaled <- build.weighted.misf(
    edges = graph$edges,
    edge.weights = graph$edge_weights * 11,
    n = graph$n,
    num.init = 6,
    num.nbrs = 8,
    seed = 29
  )

  expect_identical(misf1$levels, misf2$levels)
  expect_identical(misf1$mish_order, misf2$mish_order)
  expect_identical(misf1$vertex_depth, misf2$vertex_depth)
  expect_identical(misf1$levels, misf_scaled$levels)
  expect_true(is.finite(misf1$weight_scale))
  expect_identical(misf1$length_normalization, "median")
})

test_that("weighted mesh preset matches the explicit tuning profile", {
  graph <- mesh.surface.graph(4, 4, surface = "saddle", amplitude = 0.5)
  coords_preset <- grip(metric = "edge_length",
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 2,
    preset = "mesh",
    seed = 101
  )
  coords_explicit <- grip(metric = "edge_length",
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 2,
    placement = "barycenter",
    rounds = 128,
    final.rounds = 128,
    num.init = 12,
    num.nbrs = 20,
    r = 0.10,
    s = 4.5,
    repulsion.factor = 1.5,
    seed = 101
  )
  expect_identical(coords_preset, coords_explicit)
})

test_that("weighted cylinder preset matches the explicit tuning profile", {
  graph <- cylinder.surface.graph(4, 6, surface = "hourglass", amplitude = 0.25)
  coords_preset <- grip(metric = "edge_length",
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 3,
    preset = "cylinder",
    seed = 103
  )
  coords_explicit <- grip(metric = "edge_length",
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 3,
    placement = "barycenter",
    rounds = 160,
    final.rounds = 224,
    num.init = 14,
    num.nbrs = 22,
    r = 0.08,
    s = 5.8,
    repulsion.factor = 1.10,
    seed = 103
  )
  expect_identical(coords_preset, coords_explicit)
})

test_that("weighted sphere and irregular presets match explicit tuning profiles", {
  sphere <- sphere.surface.graph(4, 5, surface = "ellipsoid", amplitude = 0.15)
  irregular <- irregular.annulus.surface.graph(
    rings = 4,
    outer.count = 16,
    surface = "folded",
    amplitude = 0.25
  )

  sphere_preset <- grip(metric = "edge_length",
    edges = sphere$edges,
    edge.weights = sphere$edge_weights,
    n = sphere$n,
    dim = 3,
    preset = "sphere",
    seed = 107
  )
  sphere_explicit <- grip(metric = "edge_length",
    edges = sphere$edges,
    edge.weights = sphere$edge_weights,
    n = sphere$n,
    dim = 3,
    placement = "barycenter",
    rounds = 176,
    final.rounds = 240,
    num.init = 14,
    num.nbrs = 24,
    r = 0.06,
    s = 6.5,
    repulsion.factor = 0.90,
    seed = 107
  )

  irregular_preset <- grip(metric = "edge_length",
    edges = irregular$edges,
    edge.weights = irregular$edge_weights,
    n = irregular$n,
    dim = 3,
    preset = "irregular",
    seed = 109
  )
  irregular_explicit <- grip(metric = "edge_length",
    edges = irregular$edges,
    edge.weights = irregular$edge_weights,
    n = irregular$n,
    dim = 3,
    placement = "barycenter",
    rounds = 192,
    final.rounds = 256,
    num.init = 18,
    num.nbrs = 24,
    r = 0.05,
    s = 6.5,
    repulsion.factor = 1.10,
    seed = 109
  )

  expect_identical(sphere_preset, sphere_explicit)
  expect_identical(irregular_preset, irregular_explicit)
})

test_that("weighted tree preset matches explicit tuning profile and overrides cleanly", {
  tree <- kary.tree.weighted.graph(k = 2, depth = 3)
  coords_preset <- grip(metric = "edge_length",
    edges = tree$edges,
    edge.weights = tree$edge_weights,
    n = tree$n,
    dim = 2,
    preset = "tree",
    seed = 113
  )
  coords_explicit <- grip(metric = "edge_length",
    edges = tree$edges,
    edge.weights = tree$edge_weights,
    n = tree$n,
    dim = 2,
    placement = "circle",
    rounds = 64,
    final.rounds = 160,
    num.init = 28,
    num.nbrs = 8,
    r = 0.05,
    s = 7.5,
    repulsion.factor = 0.0,
    seed = 113
  )
  coords_override <- grip(metric = "edge_length",
    edges = tree$edges,
    edge.weights = tree$edge_weights,
    n = tree$n,
    dim = 2,
    preset = "tree",
    repulsion.factor = 0.5,
    seed = 127
  )
  coords_override_explicit <- grip(metric = "edge_length",
    edges = tree$edges,
    edge.weights = tree$edge_weights,
    n = tree$n,
    dim = 2,
    placement = "circle",
    rounds = 64,
    final.rounds = 160,
    num.init = 28,
    num.nbrs = 8,
    r = 0.05,
    s = 7.5,
    repulsion.factor = 0.5,
    seed = 127
  )

  expect_identical(coords_preset, coords_explicit)
  expect_identical(coords_override, coords_override_explicit)
})

test_that("invalid weighted preset is rejected", {
  graph <- mesh.surface.graph(4, 4, surface = "saddle", amplitude = 0.4)
  expect_error(
    grip(metric = "edge_length",
      edges = graph$edges,
      edge.weights = graph$edge_weights,
      n = graph$n,
      dim = 2,
      preset = "bogus",
      seed = 131
    ),
    "preset for globalrep.weighted.grip must be NULL"
  )
})

test_that("weighted metric neighbor cap preserves exact results when sufficiently large", {
  graph <- torus.surface.graph(5, 5, surface = "pinched", amplitude = 0.18)
  coords_exact <- grip(metric = "edge_length",
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 3,
    preset = "torus",
    seed = 151
  )
  coords_capped <- grip(metric = "edge_length",
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 3,
    preset = "torus",
    metric.neighbor.cap = 128,
    seed = 151
  )

  expect_identical(coords_exact, coords_capped)
})

test_that("weighted metric neighbor cap yields deterministic approximate layouts", {
  graph <- irregular.annulus.surface.graph(
    rings = 5,
    outer.count = 18,
    surface = "folded",
    amplitude = 0.3
  )
  coords1 <- grip(metric = "edge_length",
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 3,
    preset = "irregular",
    metric.neighbor.cap = 8,
    seed = 157
  )
  coords2 <- grip(metric = "edge_length",
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 3,
    preset = "irregular",
    metric.neighbor.cap = 8,
    seed = 157
  )
  coords_exact <- grip(metric = "edge_length",
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 3,
    preset = "irregular",
    seed = 157
  )

  expect_true(all(is.finite(coords1)))
  expect_identical(coords1, coords2)
  expect_gt(max(abs(coords1 - coords_exact)), 1e-8)
})

test_that("invalid weighted metric neighbor cap is rejected", {
  graph <- mesh.surface.graph(4, 4, surface = "saddle", amplitude = 0.4)
  expect_error(
    grip(metric = "edge_length",
      edges = graph$edges,
      edge.weights = graph$edge_weights,
      n = graph$n,
      dim = 2,
      metric.neighbor.cap = 0,
      seed = 163
    ),
    "metric.neighbor.cap must be a positive integer"
  )
})

test_that("weighted globalrep multiscale LGKK knobs can change the layout", {
  graph <- mesh.surface.graph(5, 5, surface = "saddle", amplitude = 0.8)
  coords_base <- globalrep.weighted.grip(
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 2,
    rounds = 8,
    final.rounds = 8,
    num.init = 6,
    num.nbrs = 8,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    lgkk.multiscale.rounds = 0,
    seed = 167
  )
  coords_lgkk <- globalrep.weighted.grip(
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 2,
    rounds = 8,
    final.rounds = 8,
    num.init = 6,
    num.nbrs = 8,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    lgkk.multiscale.rounds = 2,
    lgkk.local.nbrs = 6,
    lgkk.landmark.count = 4,
    lgkk.multiscale.scope = "all",
    lgkk.active.limit = 512,
    seed = 167
  )

  expect_gt(max(abs(coords_base - coords_lgkk)), 1e-6)
})

test_that("weighted globalrep staged LGKK budgets can change layouts", {
  graph <- cylinder.surface.graph(5, 6, surface = "hourglass", amplitude = 0.25)
  coords_shared <- globalrep.weighted.grip(
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 3,
    rounds = 8,
    final.rounds = 8,
    num.init = 6,
    num.nbrs = 8,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    lgkk.multiscale.rounds = 3,
    lgkk.local.nbrs = 6,
    lgkk.landmark.count = 6,
    lgkk.multiscale.scope = "all",
    lgkk.active.limit = 4096,
    seed = 173
  )
  coords_staged <- globalrep.weighted.grip(
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 3,
    rounds = 8,
    final.rounds = 8,
    num.init = 6,
    num.nbrs = 8,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    lgkk.multiscale.rounds = 0,
    lgkk.rounds.coarse = 1,
    lgkk.rounds.pre.final = 2,
    lgkk.rounds.final = 4,
    lgkk.local.nbrs = 6,
    lgkk.landmark.count = 6,
    lgkk.multiscale.scope = "all",
    lgkk.active.limit = 4096,
    seed = 173
  )

  expect_gt(max(abs(coords_shared - coords_staged)), 1e-6)
})

test_that("weighted globalrep validates multiscale LGKK round budgets", {
  graph <- mesh.surface.graph(4, 4, surface = "saddle", amplitude = 0.4)
  expect_error(
    globalrep.weighted.grip(
      edges = graph$edges,
      edge.weights = graph$edge_weights,
      n = graph$n,
      dim = 2,
      lgkk.rounds.coarse = -1,
      seed = 179
    ),
    "lgkk.rounds.coarse must be a non-negative integer"
  )
  expect_error(
    globalrep.weighted.grip(
      edges = graph$edges,
      edge.weights = graph$edge_weights,
      n = graph$n,
      dim = 2,
      lgkk.rounds.pre.final = -1,
      seed = 181
    ),
    "lgkk.rounds.pre.final must be a non-negative integer"
  )
  expect_error(
    globalrep.weighted.grip(
      edges = graph$edges,
      edge.weights = graph$edge_weights,
      n = graph$n,
      dim = 2,
      lgkk.rounds.final = -1,
      seed = 183
    ),
    "lgkk.rounds.final must be a non-negative integer"
  )
})
