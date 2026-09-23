test_that("paraboloid slice runs on regular and irregular meshes", {
  regular <- mesh.surface.graph(
    6, 6,
    surface = "paraboloid",
    amplitude = 0.25,
    connectivity = "orthogonal",
    normalize = "median"
  )
  irregular <- occupied.mesh.surface.graph(
    keep = keep.asymmetric.notches(6, 6, notch.depth = 2, notch.width = 1),
    surface = "paraboloid",
    amplitude = 0.25,
    connectivity = "orthogonal",
    normalize = "median"
  )

  cases <- list(regular = regular, irregular = irregular)
  for (case_name in names(cases)) {
    bundle <- cases[[case_name]]
    prepared <- prepare.graph.geodesic.mds(
      edges = bundle$edges,
      n = bundle$n,
      edge.weights = bundle$edge_weights,
      tie.mode = "average"
    )

    cmd <- grip:::grip.classical.mds.embedding(prepared, dim = 3L, eig = TRUE)
    pure_fit <- grip.optimize.geodesic.mds(
      coords = cmd$coords,
      prepared = prepared,
      engine = "cpp",
      max.iter = 4L,
      n.threads = 1L,
      return.trace = TRUE
    )
    misf_fit <- grip.optimize.misf.geodesic.mds(
      prepared = prepared,
      dim = 3L,
      top.level.restarts = 2L,
      top.level.max.iter = 3L,
      insertion.max.iter = 24L,
      refinement.local.nbrs = 4L,
      refinement.landmark.count = 2L,
      refinement.pair.mode = "sparse",
      refinement.anchor.weight = 0.05,
      refinement.anchor.weight.end = 0.01,
      refinement.continuation = "linear",
      refinement.max.iter = 3L,
      refinement.engine = "cpp",
      final.polish.max.iter = 3L,
      final.polish.engine = "cpp",
      n.threads = 1L,
      return.trace = TRUE,
      seed = 31L
    )

    expect_true(all(is.finite(cmd$coords)), info = case_name)
    expect_true(all(is.finite(pure_fit$coords)), info = case_name)
    expect_true(all(is.finite(misf_fit$coords)), info = case_name)
    expect_true(is.data.frame(misf_fit$stage_trace), info = case_name)
    expect_true(all(c("top_level", "refinement", "final_polish") %in% misf_fit$stage_trace$stage), info = case_name)
    expect_true(is.finite(misf_fit$score$final.gmds.energy[[1L]]), info = case_name)
    expect_true(is.finite(misf_fit$score$final.gmds.stress[[1L]]), info = case_name)
    expect_true(misf_fit$score$top.level.n[[1L]] >= 1L, info = case_name)
    expect_equal(nrow(misf_fit$coords), bundle$n, info = case_name)
    expect_equal(ncol(misf_fit$coords), 3L, info = case_name)
  }
})
