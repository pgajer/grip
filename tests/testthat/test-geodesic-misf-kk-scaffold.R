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

test_that("MISF-GKK optimizer stub fails with an informative message", {
  edges <- edges.mesh(4, 4)

  expect_error(
    grip.optimize.misf.geodesic.kk(
      edges = edges,
      n = 16L,
      num_init = 4L,
      seed = 1L
    ),
    "scaffolded but not implemented yet"
  )
})
