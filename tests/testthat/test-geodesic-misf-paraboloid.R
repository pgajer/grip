test_that("paraboloid slice runs on regular and irregular meshes", {
  regular <- mesh.surface.graph(
    6, 6,
    surface = "paraboloid",
    amplitude = 0.25,
    connectivity = "orthogonal",
    normalize = "median"
  )
  irregular <- occupied.mesh.surface.graph(
    keep = keep.asymmetric.notches(6, 6, notch_depth = 2, notch_width = 1),
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
      edge_weights = bundle$edge_weights,
      tie_mode = "average"
    )

    cmd <- grip:::grip.classical.mds.embedding(prepared, dim = 3L, eig = TRUE)
    pure_fit <- grip.optimize.geodesic.mds(
      coords = cmd$coords,
      prepared = prepared,
      engine = "cpp",
      max_iter = 4L,
      n_threads = 1L,
      return_trace = TRUE
    )
    misf_fit <- grip.optimize.misf.geodesic.mds(
      prepared = prepared,
      dim = 3L,
      top_level_restarts = 2L,
      top_level_max_iter = 3L,
      insertion_max_iter = 24L,
      refinement_local_nbrs = 4L,
      refinement_landmark_count = 2L,
      refinement_pair_mode = "sparse",
      refinement_anchor_weight = 0.05,
      refinement_anchor_weight_end = 0.01,
      refinement_continuation = "linear",
      refinement_max_iter = 3L,
      refinement_engine = "cpp",
      final_polish_max_iter = 3L,
      final_polish_engine = "cpp",
      n_threads = 1L,
      return_trace = TRUE,
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
