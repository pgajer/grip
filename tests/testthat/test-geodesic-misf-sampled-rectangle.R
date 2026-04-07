test_that("MISF-GMDS runs on sampled rectangle paraboloid k-sequences", {
  sequences <- list(
    sampled.rectangle.surface.graphs(
      n = 50,
      k = c(6, 7, 8),
      xmin = -1.6,
      xmax = 1.6,
      ymin = -1,
      ymax = 1,
      seed = 1050,
      surface = "paraboloid",
      amplitude = 0.35,
      graph_space = "surface",
      normalize = "median"
    ),
    sampled.rectangle.surface.graphs(
      n = 75,
      k = c(6, 7, 8),
      xmin = -1.6,
      xmax = 1.6,
      ymin = -1,
      ymax = 1,
      seed = 1075,
      surface = "paraboloid",
      amplitude = 0.35,
      graph_space = "surface",
      normalize = "median"
    )
  )

  for (seq_spec in sequences) {
    shared_param <- seq_spec$coords_param
    shared_surface <- seq_spec$coords_surface
    for (graph_name in names(seq_spec$graphs)) {
      graph <- seq_spec$graphs[[graph_name]]
      info <- sprintf("n=%d, k=%d", graph$n, graph$k)

      expect_equal(graph$coords_param, shared_param, info = info)
      expect_equal(graph$coords_surface, shared_surface, info = info)

      prepared <- grip.prepare.misf.geodesic.mds(
        edges = graph$edges,
        n = graph$n,
        edge_weights = graph$edge_weights,
        tie_mode = "average",
        num_init = 6L,
        dim = 3L,
        top_level_mode = "skip",
        seed = 3000L + graph$n + graph$k
      )
      fit <- grip.optimize.misf.geodesic.mds(
        prepared = prepared,
        dim = 3L,
        top_level_restarts = 1L,
        top_level_max_iter = 2L,
        top_level_engine = "cpp",
        insertion_anchor_policy = "prev_level_spread",
        insertion_max_iter = 12L,
        refinement_local_nbrs = 4L,
        refinement_landmark_count = 2L,
        refinement_pair_mode = "sparse",
        refinement_anchor_weight = 0.05,
        refinement_anchor_weight_end = 0.01,
        refinement_continuation = "linear",
        refinement_max_iter = 2L,
        refinement_engine = "cpp",
        final_polish_max_iter = 2L,
        final_polish_engine = "cpp",
        n_threads = 1L,
        return_trace = TRUE,
        return_frames = FALSE,
        seed = 3000L + graph$n + graph$k
      )

      expect_true(length(prepared$misf$levels) >= 2L, info = info)
      expect_true(prepared$top_level_level >= 1L, info = info)
      expect_equal(nrow(fit$coords), graph$n, info = info)
      expect_equal(ncol(fit$coords), 3L, info = info)
      expect_true(all(is.finite(fit$coords)), info = info)
      expect_true(is.data.frame(fit$stage_trace), info = info)
      expect_true(all(c("top_level", "refinement", "final_polish") %in% fit$stage_trace$stage), info = info)
      expect_true(is.finite(fit$score$final.gmds.energy[[1L]]), info = info)
      expect_true(is.finite(fit$score$final.gmds.stress[[1L]]), info = info)
    }
  }
})
