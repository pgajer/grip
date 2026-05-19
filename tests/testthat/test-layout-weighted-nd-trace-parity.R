test_that("weighted ND trace parity diagnostics run when explicitly enabled", {
  skip_if_not(
    identical(Sys.getenv("GRIP_RUN_TRACE_PARITY_TESTS"), "true"),
    "Set GRIP_RUN_TRACE_PARITY_TESTS=true to run weighted ND trace parity diagnostics."
  )

  metrics <- grip_weighted_nd_trace_parity_run(dims = c(2L, 3L))
  first <- grip_weighted_nd_trace_parity_first_divergence(metrics)
  neighbors <- grip_weighted_nd_trace_parity_neighbor_run(dims = c(2L, 3L))

  expect_s3_class(metrics, "data.frame")
  expect_s3_class(first, "data.frame")
  expect_s3_class(neighbors, "data.frame")
  expect_true(nrow(metrics) > 0L)
  expect_true(nrow(first) > 0L)
  expect_true(nrow(neighbors) > 0L)
  expect_true(all(metrics$metadata_match))
  expect_true(all(metrics$frame_count_match))
  expect_true(all(is.finite(metrics$direct_rmse)))
  expect_true(all(is.finite(metrics$direct_max_abs)))
  expect_true(all(is.finite(metrics$procrustes_rmse)))
  expect_true(all(neighbors$neighbor_match))
  expect_true(all(neighbors$distance_match))

  if (identical(Sys.getenv("GRIP_ENFORCE_TRACE_PARITY"), "true")) {
    expect_true(all(first$status == "within_tolerance"))
    expect_lte(max(metrics$direct_max_abs), 1e-10)
    expect_lte(max(metrics$procrustes_rmse), 1e-8)
  }
})

test_that("weighted ND final-anchor trace parity runs when explicitly enabled", {
  skip_if_not(
    identical(Sys.getenv("GRIP_RUN_FINAL_ANCHOR_TRACE_PARITY_TESTS"), "true"),
    "Set GRIP_RUN_FINAL_ANCHOR_TRACE_PARITY_TESTS=true to run final-anchor trace parity diagnostics."
  )

  tuning_fun <- function(dim) {
    tuning <- grip_weighted_nd_trace_parity_tuning(dim)
    tuning$final_anchor_factor <- 0.7
    tuning$final_move_scale_after_first <- 1
    tuning
  }
  metrics <- grip_weighted_nd_trace_parity_run(
    dims = c(2L, 3L),
    tuning_fun = tuning_fun
  )
  first <- grip_weighted_nd_trace_parity_first_divergence(metrics)

  expect_s3_class(metrics, "data.frame")
  expect_s3_class(first, "data.frame")
  expect_true(all(metrics$metadata_match))
  expect_true(all(metrics$frame_count_match))
  expect_true(all(is.finite(metrics$direct_rmse)))
  expect_true(all(is.finite(metrics$procrustes_rmse)))

  if (identical(Sys.getenv("GRIP_ENFORCE_FINAL_ANCHOR_TRACE_PARITY"), "true")) {
    expect_true(all(first$status == "within_tolerance"))
    expect_lte(max(metrics$direct_max_abs), 1e-10)
    expect_lte(max(metrics$procrustes_rmse), 1e-8)
  }
})
