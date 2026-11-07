test_that("MISF-GMDS and MISF-GKK expose the same canonical stage schema", {
  edges <- edges.mesh(4, 4)

  gmds.fit <- grip.optimize.misf.geodesic.mds(
    edges = edges,
    n = 16L,
    num.init = 4L,
    dim = 2L,
    top.level.restarts = 2L,
    top.level.max.iter = 2L,
    insertion.max.iter = 6L,
    refinement.local.nbrs = 3L,
    refinement.landmark.count = 2L,
    refinement.max.iter = 2L,
    refinement.engine = "cpp",
    final.polish.max.iter = 2L,
    final.polish.engine = "cpp",
    n.threads = 1L,
    return.trace = TRUE,
    return.frames = TRUE,
    seed = 101L
  )

  gkk.fit <- misf.geodesic.kk(
    edges = edges,
    n = 16L,
    num.init = 4L,
    dim = 2L,
    top.level.restarts = 2L,
    top.level.max.iter = 2L,
    insertion.mode = "geodesic",
    insertion.max.iter = 6L,
    refinement.pair.mode = "auto",
    refinement.full.limit = 4L,
    refinement.max.iter = 2L,
    final.pair.mode = "landmark",
    final.full.limit = 4L,
    final.max.iter = 2L,
    return.trace = TRUE,
    return.frames = TRUE,
    seed = 202L
  )

  expect_identical(names(gmds.fit$stage_trace), names(gkk.fit$stage_trace))
  expect_equal(unique(gmds.fit$stage_trace$trace_schema_version), 1L)
  expect_equal(unique(gkk.fit$stage_trace$trace_schema_version), 1L)
  expect_true(all(c("seed", "initial_placement", "top_level", "final_polish") %in% gmds.fit$stage_trace$stage))
  expect_true(all(c("seed", "initial_placement", "top_level", "final_polish") %in% gkk.fit$stage_trace$stage))
})

test_that("canonical stage lookup returns per-stage data records", {
  edges <- edges.mesh(4, 4)

  fit <- grip.optimize.misf.geodesic.mds(
    edges = edges,
    n = 16L,
    num.init = 4L,
    dim = 2L,
    top.level.restarts = 2L,
    top.level.max.iter = 2L,
    insertion.max.iter = 6L,
    refinement.local.nbrs = 3L,
    refinement.landmark.count = 2L,
    refinement.max.iter = 2L,
    refinement.engine = "cpp",
    final.polish.max.iter = 2L,
    final.polish.engine = "cpp",
    n.threads = 1L,
    return.trace = TRUE,
    return.frames = TRUE,
    seed = 303L
  )

  seed.stage <- grip.geodesic.misf.trace.stage.lookup(
    fit,
    stage = "seed",
    level = fit$prepared$top_level_level
  )
  top.stage <- grip.geodesic.misf.trace.stage.lookup(
    fit,
    stage = "top_level",
    level = fit$prepared$top_level_level
  )
  final.stage <- grip.geodesic.misf.trace.stage.lookup(
    fit,
    stage = "final_polish",
    level = 0L
  )

  expect_true(is.list(seed.stage))
  expect_true(is.matrix(seed.stage$coords_full))
  expect_true(is.list(top.stage$frames))
  expect_true(is.matrix(top.stage$coords_full))
  expect_true(is.matrix(final.stage$coords_full))
  expect_equal(final.stage$label, "Final full-graph polish")
})

test_that("stage payload helper returns app-facing boundary states", {
  edges <- edges.mesh(4, 4)

  fit <- grip.optimize.misf.geodesic.mds(
    edges = edges,
    n = 16L,
    num.init = 4L,
    dim = 2L,
    top.level.restarts = 2L,
    top.level.max.iter = 2L,
    insertion.max.iter = 6L,
    refinement.local.nbrs = 3L,
    refinement.landmark.count = 2L,
    refinement.max.iter = 2L,
    refinement.engine = "cpp",
    final.polish.max.iter = 2L,
    final.polish.engine = "cpp",
    n.threads = 1L,
    return.trace = TRUE,
    return.frames = TRUE,
    seed = 404L
  )

  payloads <- grip.geodesic.misf.trace.stage.payloads(
    fit,
    target = fit$coords,
    states = c("top_level", "after_insertion", "after_refinement", "final_polish")
  )

  expect_identical(names(payloads), c("top_level", "after_insertion", "after_refinement", "final_polish"))
  expect_true(all(vapply(payloads, is.list, logical(1L))))
  expect_true(all(vapply(payloads, function(x) is.matrix(x$coords), logical(1L))))
  expect_true(all(vapply(payloads, function(x) is.matrix(x$display_coords), logical(1L))))
  expect_true(is.matrix(payloads$final_polish$active_edges))
  expect_equal(payloads$final_polish$stage, "final_polish")
})
