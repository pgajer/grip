test_that("gripui paraboloid GMDS payload computes expected outputs", {
  skip_if_not_installed("rgl")

  compute_payload <- getFromNamespace("gripui.paraboloid.gmds.compute.payload", "grip")
  payload <- compute_payload(
    side = 4L,
    amplitude = 0.25,
    lambda_anchor = 0.10,
    lambda_edge = 0.20,
    lambda_repulsion = 0.30,
    max_iter = 1L,
    n_threads = 0L
  )

  expect_true(is.list(payload))
  expect_equal(payload$case$side, 4L)
  expect_equal(nrow(payload$metric_table), 3L)
  expect_equal(
    payload$metric_table$Method,
    c("Reference paraboloid", "Anchor + Repulsion", "Anchor + Edge Spring")
  )
  expect_equal(dim(payload$reference$display_coords), c(16L, 3L))
  expect_equal(dim(payload$anchor_repulsion$display_coords), c(16L, 3L))
  expect_equal(dim(payload$anchor_edge_spring$display_coords), c(16L, 3L))
})

test_that("gripui paraboloid GMDS app builds", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("rgl")

  old <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old), add = TRUE)

  app <- gripui_paraboloid_gmds_app()
  expect_s3_class(app, "shiny.appobj")
})
