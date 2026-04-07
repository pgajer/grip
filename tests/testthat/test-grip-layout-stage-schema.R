test_that("grip.layout.trace exposes the canonical MISF stage schema", {
  edges <- edges.mesh(6, 6)

  tr <- grip.layout.trace(
    edges = edges,
    n = 36L,
    dim = 3L,
    placement = "barycenter",
    rounds = 2L,
    final_rounds = 2L,
    num_init = 6L,
    num_nbrs = 8L,
    trace = "round",
    trace.every = 1L,
    seed = 77L
  )

  expect_s3_class(tr, "grip_layout_trace")
  expect_true(is.data.frame(tr$stage_trace))
  expect_true(length(tr$stage_data) >= 4L)
  expect_equal(unique(tr$stage_trace$trace_schema_version), 1L)
  expect_true(all(c("seed", "initial_placement", "top_level", "insertion", "refinement", "final_polish") %in% tr$stage_trace$stage))

  top.level <- max(tr$stage_trace$level)
  seed.stage <- grip.geodesic.misf.trace.stage.lookup(
    tr,
    stage = "seed",
    level = top.level
  )
  top.stage <- grip.geodesic.misf.trace.stage.lookup(
    tr,
    stage = "top_level",
    level = top.level
  )
  final.stage <- grip.geodesic.misf.trace.state.record(
    tr,
    state = "final_polish"
  )

  expect_equal(seed.stage$label, sprintf("Recorded start of V_%d", top.level))
  expect_equal(top.stage$label, "Top-level GRIP solve")
  expect_true(is.matrix(seed.stage$coords_full))
  expect_true(is.list(top.stage$frames))
  expect_true(is.matrix(final.stage$coords_full))
})
