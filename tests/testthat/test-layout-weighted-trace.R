test_that("weighted trace requires edge weights", {
  edges <- edges.mesh(4, 4)
  expect_error(
    trace.grip(metric = "edge_length", edges = edges, n = 16, dim = 2, seed = 1),
    "requires edge weights"
  )
})

test_that("weighted round trace returns frames, metadata, and inactive NA rows", {
  graph <- mesh.surface.graph(4, 5, surface = "saddle", amplitude = 0.6)
  tr <- trace.grip(metric = "edge_length",
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 2,
    rounds = 6,
    final.rounds = 4,
    num.init = 5,
    num.nbrs = 6,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
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
  tr <- trace.grip(metric = "edge_length",
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
    trace = "level",
    trace.every = 1,
    seed = 41
  )
  coords <- grip(metric = "edge_length",
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
  expect_identical(tr$final, coords)
})

test_that("weighted trace final matches weighted layout for LGKK polish", {
  graph <- mesh.surface.graph(5, 5, surface = "paraboloid", amplitude = 0.4)
  tr <- trace.grip(metric = "edge_length",
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
    lgkk.polish.rounds = 4,
    lgkk.local.nbrs = 6,
    lgkk.landmark.count = 4,
    trace = "level",
    trace.every = 1,
    seed = 47
  )
  coords <- grip(metric = "edge_length",
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
    lgkk.polish.rounds = 4,
    lgkk.local.nbrs = 6,
    lgkk.landmark.count = 4,
    seed = 47
  )
  expect_identical(tr$final, coords)
  expect_true(nrow(tr$lgkk.polish) >= 1L)
})

test_that("weighted trace diagnostics are produced", {
  graph <- mesh.surface.graph(4, 5, surface = "saddle", amplitude = 0.6)
  tr <- trace.grip(metric = "edge_length",
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 2,
    rounds = 6,
    final.rounds = 4,
    num.init = 5,
    num.nbrs = 6,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    trace = "level",
    trace.every = 1,
    diagnostics = "light",
    seed = 131
  )

  expect_equal(nrow(tr$diagnostics), nrow(tr$meta))
  expect_true(all(c("edge.length.cv", "sampled.nonedge.sep.ratio", "sampled.stress") %in%
                    names(tr$diagnostics)))
})

test_that("weighted trace presets match explicit tuning profiles", {
  mesh <- mesh.surface.graph(4, 4, surface = "saddle", amplitude = 0.45)
  irregular <- irregular.annulus.surface.graph(
    rings = 4,
    outer.count = 16,
    surface = "folded",
    amplitude = 0.25
  )

  tr_mesh_preset <- trace.grip(metric = "edge_length",
    edges = mesh$edges,
    edge.weights = mesh$edge_weights,
    n = mesh$n,
    dim = 2,
    preset = "mesh",
    trace = "level",
    trace.every = 1,
    seed = 211
  )
  tr_mesh_explicit <- trace.grip(metric = "edge_length",
    edges = mesh$edges,
    edge.weights = mesh$edge_weights,
    n = mesh$n,
    dim = 2,
    placement = "barycenter",
    rounds = 128,
    final.rounds = 128,
    num.init = 12,
    num.nbrs = 20,
    r = 0.10,
    s = 4.5,
    repulsion.factor = 1.5,
    trace = "level",
    trace.every = 1,
    seed = 211
  )

  tr_irregular_preset <- trace.grip(metric = "edge_length",
    edges = irregular$edges,
    edge.weights = irregular$edge_weights,
    n = irregular$n,
    dim = 3,
    preset = "irregular",
    trace = "level",
    trace.every = 1,
    seed = 223
  )
  tr_irregular_explicit <- trace.grip(metric = "edge_length",
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
    trace = "level",
    trace.every = 1,
    seed = 223
  )

  expect_identical(tr_mesh_preset$final, tr_mesh_explicit$final)
  expect_identical(tr_mesh_preset$meta, tr_mesh_explicit$meta)
  expect_identical(tr_irregular_preset$final, tr_irregular_explicit$final)
  expect_identical(tr_irregular_preset$meta, tr_irregular_explicit$meta)
})

test_that("weighted trace rejects invalid presets", {
  graph <- mesh.surface.graph(4, 4, surface = "saddle", amplitude = 0.4)
  expect_error(
    trace.grip(metric = "edge_length",
      edges = graph$edges,
      edge.weights = graph$edge_weights,
      n = graph$n,
      dim = 2,
      preset = "bogus",
      seed = 227
    ),
    "preset for trace.grip must be NULL"
  )
})

test_that("weighted trace respects metric neighbor cap and matches weighted layout", {
  graph <- torus.surface.graph(5, 5, surface = "pinched", amplitude = 0.18)
  tr <- trace.grip(metric = "edge_length",
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 3,
    preset = "torus",
    metric.neighbor.cap = 16,
    trace = "level",
    trace.every = 1,
    seed = 229
  )
  coords <- grip(metric = "edge_length",
    edges = graph$edges,
    edge.weights = graph$edge_weights,
    n = graph$n,
    dim = 3,
    preset = "torus",
    metric.neighbor.cap = 16,
    seed = 229
  )

  expect_identical(tr$final, coords)
  expect_true(nrow(tr$meta) >= 2L)
})

test_that("weighted trace rejects invalid metric neighbor caps", {
  graph <- mesh.surface.graph(4, 4, surface = "saddle", amplitude = 0.4)
  expect_error(
    trace.grip(metric = "edge_length",
      edges = graph$edges,
      edge.weights = graph$edge_weights,
      n = graph$n,
      dim = 2,
      metric.neighbor.cap = -1,
      seed = 233
    ),
    "metric.neighbor.cap must be a positive integer"
  )
})

test_that("weighted trace final matches weighted layout for multiscale LGKK", {
  graph <- mesh.surface.graph(5, 5, surface = "saddle", amplitude = 0.8)
  tr <- trace.grip(metric = "edge_length",
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
    trace = "round",
    trace.every = 1,
    seed = 239
  )
  coords <- grip(metric = "edge_length",
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
    seed = 239
  )

  expect_identical(tr$final, coords)
  expect_true(any(tr$meta$phase == "lgkk"))
})

test_that("weighted trace final matches weighted layout for staged multiscale LGKK", {
  graph <- cylinder.surface.graph(5, 6, surface = "hourglass", amplitude = 0.25)
  tr <- trace.grip(metric = "edge_length",
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
    trace = "round",
    trace.every = 1,
    seed = 241
  )
  coords <- grip(metric = "edge_length",
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
    seed = 241
  )

  expect_identical(tr$final, coords)
  expect_true(any(tr$meta$phase == "lgkk"))
})
