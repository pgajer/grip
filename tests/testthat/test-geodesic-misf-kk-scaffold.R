test_that("MISF-GKK prepare returns a layered prepared object", {
  edges <- edges.mesh(4, 4)

  prepared <- grip.prepare.misf.geodesic.kk(
    edges = edges,
    n = 16L,
    tie_mode = "average",
    num_init = 4L,
    top_level_mode = "skip",
    seed = 1L
  )

  expect_s3_class(prepared, "grip_misf_gkk_prepared")
  expect_s3_class(prepared, "grip_gkk_prepared")
  expect_equal(prepared$top_level_pair_mode, "auto")
  expect_true(prepared$top_level_effective_pair_mode %in% c("full", "landmark"))
  expect_true(is.list(prepared$top_level_prepared_full))
  expect_true(is.list(prepared$top_level_prepared_sparse))
})

test_that("MISF-GKK prepare can solve the top MISF level immediately", {
  edges <- edges.mesh(4, 4)

  prepared <- grip.prepare.misf.geodesic.kk(
    edges = edges,
    n = 16L,
    num_init = 4L,
    top_level_mode = "solve",
    top_level_restarts = 2L,
    top_level_max_iter = 3L,
    seed = 1L
  )

  expect_s3_class(prepared, "grip_misf_gkk_prepared")
  expect_true(is.list(prepared$top_level_fit))
  expect_true(all(is.finite(prepared$top_level_fit$coords)))
  expect_true(all(is.finite(prepared$top_level_fit$coords_full[prepared$top_level_vertices, , drop = FALSE])))
})

test_that("MISF-GKK scorer summarizes external coordinates against the prepared object", {
  edges <- edges.mesh(4, 4)
  prepared <- grip.prepare.misf.geodesic.kk(
    edges = edges,
    n = 16L,
    num_init = 4L,
    top_level_mode = "skip",
    seed = 1L
  )
  coords <- matrix(runif(32L), nrow = 16L, ncol = 2L)

  score <- grip.score.misf.geodesic.kk(
    coords = coords,
    prepared = prepared,
    return_trace = TRUE
  )

  expect_true(is.data.frame(score))
  expect_equal(nrow(score), 1L)
  expect_true(all(c(
    "multiscale.mode",
    "top.level",
    "final.gkk.energy",
    "final.gkk.weighted.rel.rmse",
    "stage.trace"
  ) %in% names(score)))
  expect_equal(score$multiscale.mode[[1L]], "misf")
  expect_true(is.list(score$stage.trace))
})

test_that("MISF-GKK scorer supports landmark scoring mode", {
  edges <- edges.mesh(4, 4)
  prepared <- grip.prepare.misf.geodesic.kk(
    edges = edges,
    n = 16L,
    num_init = 4L,
    top_level_mode = "skip",
    seed = 1L
  )
  coords <- matrix(runif(32L), nrow = 16L, ncol = 2L)

  score <- grip.score.misf.geodesic.kk(
    coords = coords,
    prepared = prepared,
    score_pair_mode = "landmark",
    score_full_limit = 4L,
    score_local_nbrs = 4L,
    score_landmark_count = 2L
  )

  expect_equal(score$final.score.mode[[1L]], "landmark")
  expect_equal(score$final.score.requested.mode[[1L]], "landmark")
  expect_true("final.lgkk.energy" %in% names(score))
  expect_true(is.finite(score$final.energy[[1L]]))
})

test_that("MISF-GKK optimizer returns a multistage fit with traces", {
  edges <- edges.mesh(4, 4)
  fit <- grip.optimize.misf.geodesic.kk(
    edges = edges,
    n = 16L,
    num_init = 4L,
    dim = 2L,
    top_level_restarts = 2L,
    top_level_max_iter = 3L,
    insertion_mode = "geodesic",
    insertion_max_iter = 6L,
    refinement_pair_mode = "auto",
    refinement_full_limit = 4L,
    refinement_max_iter = 3L,
    final_pair_mode = "landmark",
    final_full_limit = 4L,
    final_max_iter = 3L,
    return_trace = TRUE,
    seed = 1L
  )

  expect_s3_class(fit, "grip_misf_gkk_fit")
  expect_true(all(is.finite(fit$coords)))
  expect_true(is.data.frame(fit$stage_trace))
  expect_true(all(c("top_level", "refinement", "final_polish") %in% fit$stage_trace$stage))
  expect_true(all(c("trace_schema_version", "stage_key", "pair_mode", "weighted_rel_rmse") %in% names(fit$stage_trace)))
  expect_equal(fit$stage_trace$pair_mode[fit$stage_trace$stage == "final_polish"][1L], "landmark")
  expect_true(is.list(fit$trace))
  expect_true(is.list(fit$trace$stage_data))
  expect_false(is.null(grip.geodesic.misf.trace.stage.lookup(fit, stage = "top_level", level = fit$prepared$top_level_level)))
  expect_true(is.data.frame(fit$trace$top_restart_summary))
  expect_true(is.data.frame(fit$trace$final_polish_trace))
  expect_true(is.data.frame(fit$score))
  expect_equal(nrow(fit$score), 1L)
  expect_true(is.finite(fit$score$final.gkk.energy[[1L]]))
  expect_true(any(vapply(
    fit$refinement$level_results,
    function(result) {
      !is.null(result$fit$final_anchor_weight) && is.finite(result$fit$final_anchor_weight) &&
        result$fit$final_anchor_weight > 0
    },
    logical(1L)
  )))
})
