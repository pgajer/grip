test_that("MISF-GMDS and MISF-GKK expose the same canonical stage schema", {
  edges <- edges.mesh(4, 4)

  gmds.fit <- grip.optimize.misf.geodesic.mds(
    edges = edges,
    n = 16L,
    num_init = 4L,
    dim = 2L,
    top_level_restarts = 2L,
    top_level_max_iter = 2L,
    insertion_max_iter = 6L,
    refinement_local_nbrs = 3L,
    refinement_landmark_count = 2L,
    refinement_max_iter = 2L,
    refinement_engine = "cpp",
    final_polish_max_iter = 2L,
    final_polish_engine = "cpp",
    n_threads = 1L,
    return_trace = TRUE,
    return_frames = TRUE,
    seed = 101L
  )

  gkk.fit <- grip.optimize.misf.geodesic.kk(
    edges = edges,
    n = 16L,
    num_init = 4L,
    dim = 2L,
    top_level_restarts = 2L,
    top_level_max_iter = 2L,
    insertion_mode = "geodesic",
    insertion_max_iter = 6L,
    refinement_pair_mode = "auto",
    refinement_full_limit = 4L,
    refinement_max_iter = 2L,
    final_pair_mode = "landmark",
    final_full_limit = 4L,
    final_max_iter = 2L,
    return_trace = TRUE,
    return_frames = TRUE,
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
    num_init = 4L,
    dim = 2L,
    top_level_restarts = 2L,
    top_level_max_iter = 2L,
    insertion_max_iter = 6L,
    refinement_local_nbrs = 3L,
    refinement_landmark_count = 2L,
    refinement_max_iter = 2L,
    refinement_engine = "cpp",
    final_polish_max_iter = 2L,
    final_polish_engine = "cpp",
    n_threads = 1L,
    return_trace = TRUE,
    return_frames = TRUE,
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
    num_init = 4L,
    dim = 2L,
    top_level_restarts = 2L,
    top_level_max_iter = 2L,
    insertion_max_iter = 6L,
    refinement_local_nbrs = 3L,
    refinement_landmark_count = 2L,
    refinement_max_iter = 2L,
    refinement_engine = "cpp",
    final_polish_max_iter = 2L,
    final_polish_engine = "cpp",
    n_threads = 1L,
    return_trace = TRUE,
    return_frames = TRUE,
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
