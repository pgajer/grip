test_that("gflow-backed weighted ND parity diagnostics run when explicitly enabled", {
  skip_if_not(
    identical(Sys.getenv("GRIP_RUN_GFLOW_PARITY_TESTS"), "true"),
    "Set GRIP_RUN_GFLOW_PARITY_TESTS=true to run gflow-backed parity diagnostics."
  )
  skip_if_not_installed("gflow")

  result <- tryCatch(
    grip_weighted_nd_gflow_parity_run(mode = "smoke", dims = c(2L, 3L)),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    skip(paste("gflow-backed parity fixtures unavailable:", conditionMessage(result)))
  }

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 12L)
  expect_true(all(result$finite))
  expect_true(all(is.finite(result$direct_rmse)))
  expect_true(all(is.finite(result$procrustes_rmse)))

  if (identical(Sys.getenv("GRIP_ENFORCE_GFLOW_PARITY"), "true")) {
    thresholds <- grip_weighted_nd_gflow_parity_thresholds()
    expect_lte(max(result$direct_max_abs), thresholds$direct_max_abs)
    expect_lte(max(result$procrustes_rmse), thresholds$procrustes_rmse)
    expect_lte(max(result$procrustes_max_abs), thresholds$procrustes_max_abs)
  }
})

test_that("expanded gflow-backed weighted ND parity diagnostics run when explicitly enabled", {
  skip_if_not(
    identical(Sys.getenv("GRIP_RUN_GFLOW_FULL_PARITY_TESTS"), "true"),
    "Set GRIP_RUN_GFLOW_FULL_PARITY_TESTS=true to run expanded gflow-backed parity diagnostics."
  )
  skip_if_not_installed("gflow")

  result <- tryCatch(
    grip_weighted_nd_gflow_parity_run(mode = "full", dims = c(2L, 3L)),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    skip(paste("expanded gflow-backed parity fixtures unavailable:", conditionMessage(result)))
  }

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 32L)
  expect_equal(length(unique(result$case_id)), 16L)
  expect_true(all(result$finite))
  expect_true(all(is.finite(result$direct_rmse)))
  expect_true(all(is.finite(result$procrustes_rmse)))
  expect_true(any(grepl("quadform_2d_adaptive_radius_index_0", result$case_id)))
  expect_true(any(grepl("quadform_2d_cknn_index_0", result$case_id)))
  expect_true(any(grepl("quadform_3d_adaptive_radius_index_1", result$case_id)))
  expect_true(any(grepl("quadform_3d_cknn_index_1", result$case_id)))

  if (identical(Sys.getenv("GRIP_ENFORCE_GFLOW_FULL_PARITY"), "true")) {
    thresholds <- grip_weighted_nd_gflow_parity_thresholds()
    expect_lte(max(result$direct_max_abs), thresholds$direct_max_abs)
    expect_lte(max(result$procrustes_rmse), thresholds$procrustes_rmse)
    expect_lte(max(result$procrustes_max_abs), thresholds$procrustes_max_abs)
  }
})

test_that("stress gflow-backed weighted ND parity diagnostics run in parallel when explicitly enabled", {
  skip_if_not(
    identical(Sys.getenv("GRIP_RUN_GFLOW_STRESS_PARITY_TESTS"), "true"),
    "Set GRIP_RUN_GFLOW_STRESS_PARITY_TESTS=true to run stress gflow-backed parity diagnostics."
  )
  skip_if_not_installed("gflow")

  output_dir <- file.path(
    tempdir(),
    "weighted-grip-nd-parity",
    "stress-parity-parts"
  )
  result <- tryCatch(
    grip_weighted_nd_gflow_parity_run_parallel(
      mode = "stress",
      dims = c(2L, 3L),
      output_dir = output_dir
    ),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    skip(paste("stress gflow-backed parity fixtures unavailable:", conditionMessage(result)))
  }

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 72L)
  expect_equal(length(unique(result$case_id)), 36L)
  expect_true(all(result$finite))
  expect_true(all(is.finite(result$direct_rmse)))
  expect_true(all(is.finite(result$procrustes_rmse)))
  expect_true(any(grepl("quadform_2d_.*index_0_curv_low", result$case_id)))
  expect_true(any(grepl("quadform_2d_.*index_1_curv_high", result$case_id)))
  expect_true(any(grepl("quadform_3d_.*index_0_curv_low", result$case_id)))
  expect_true(any(grepl("quadform_3d_.*index_1_curv_high", result$case_id)))

  if (identical(Sys.getenv("GRIP_ENFORCE_GFLOW_STRESS_PARITY"), "true")) {
    thresholds <- grip_weighted_nd_gflow_parity_thresholds()
    expect_lte(max(result$direct_max_abs), thresholds$direct_max_abs)
    expect_lte(max(result$procrustes_rmse), thresholds$procrustes_rmse)
    expect_lte(max(result$procrustes_max_abs), thresholds$procrustes_max_abs)
  }
})

test_that("gflow-backed weighted ND trace ratchet runs when explicitly enabled", {
  skip_if_not(
    identical(Sys.getenv("GRIP_RUN_GFLOW_TRACE_PARITY_TESTS"), "true"),
    "Set GRIP_RUN_GFLOW_TRACE_PARITY_TESTS=true to run gflow-backed trace diagnostics."
  )
  skip_if_not_installed("gflow")

  result <- tryCatch(
    grip_weighted_nd_gflow_trace_run(mode = "smoke", dims = "intrinsic"),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    skip(paste("gflow-backed trace fixtures unavailable:", conditionMessage(result)))
  }
  first <- grip_weighted_nd_gflow_trace_first_divergence(result)

  expect_s3_class(result, "data.frame")
  expect_s3_class(first, "data.frame")
  expect_equal(length(unique(result$case_id)), 4L)
  expect_true(all(result$metadata_match))
  expect_true(all(result$frame_count_match))
  expect_true(all(is.finite(result$direct_rmse)))
  expect_true(all(is.finite(result$procrustes_rmse)))

  if (identical(Sys.getenv("GRIP_ENFORCE_GFLOW_TRACE_PARITY"), "true")) {
    expect_true(all(first$status == "within_tolerance"))
    expect_lte(max(result$direct_max_abs), 1e-10)
    expect_lte(max(result$procrustes_rmse), 1e-8)
  }
})

test_that("gflow-backed weighted ND cross-dimensional trace ratchet runs when explicitly enabled", {
  skip_if_not(
    identical(Sys.getenv("GRIP_RUN_GFLOW_TRACE_ALL_TESTS"), "true"),
    "Set GRIP_RUN_GFLOW_TRACE_ALL_TESTS=true to run cross-dimensional gflow trace diagnostics."
  )
  skip_if_not_installed("gflow")

  result <- tryCatch(
    grip_weighted_nd_gflow_trace_run(mode = "smoke", dims = "all"),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    skip(paste("gflow-backed trace fixtures unavailable:", conditionMessage(result)))
  }
  first <- grip_weighted_nd_gflow_trace_first_divergence(result)

  expect_s3_class(result, "data.frame")
  expect_s3_class(first, "data.frame")
  expect_equal(length(unique(interaction(result$case_id, result$layout_dim))), 8L)
  expect_true(all(result$metadata_match))
  expect_true(all(result$frame_count_match))
  expect_true(all(is.finite(result$direct_rmse)))
  expect_true(all(is.finite(result$procrustes_rmse)))

  if (identical(Sys.getenv("GRIP_ENFORCE_GFLOW_TRACE_ALL_PARITY"), "true")) {
    expect_true(all(first$status == "within_tolerance"))
    expect_lte(max(result$direct_max_abs), 1e-10)
    expect_lte(max(result$procrustes_rmse), 1e-8)
  }
})

test_that("gflow-backed weighted ND refinement-step trace runs when explicitly enabled", {
  skip_if_not(
    identical(Sys.getenv("GRIP_RUN_GFLOW_REFINEMENT_STEP_TRACE_TESTS"), "true"),
    "Set GRIP_RUN_GFLOW_REFINEMENT_STEP_TRACE_TESTS=true to run refinement-step trace diagnostics."
  )
  skip_if_not_installed("gflow")

  result <- tryCatch(
    grip_weighted_nd_gflow_refinement_step_trace_run(mode = "smoke"),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    skip(paste("gflow-backed refinement-step trace unavailable:", conditionMessage(result)))
  }
  summary <- grip_weighted_nd_gflow_refinement_step_summary(result)
  term_summary <- grip_weighted_nd_gflow_refinement_step_term_summary(result)
  source <- grip_weighted_nd_gflow_refinement_source_run(mode = "smoke")
  terms <- attr(result, "attraction_terms")

  expect_s3_class(result, "data.frame")
  expect_s3_class(summary, "data.frame")
  expect_s3_class(term_summary, "data.frame")
  expect_s3_class(terms, "data.frame")
  expect_type(source, "list")
  expect_s3_class(source$target, "data.frame")
  expect_s3_class(source$summary, "data.frame")
  expect_s3_class(source$vertex_trace, "data.frame")
  expect_s3_class(source$term_trace, "data.frame")
  expect_gt(nrow(result), 0L)
  expect_gt(nrow(terms), 0L)
  expect_gt(nrow(source$term_trace), 0L)
  expect_true(all(result$layout_dim == 2L))
  expect_true(all(result$misf_level == 0L))
  expect_true(all(result$round_in_level >= 1L))
  expect_true(all(is.finite(result$heat_before)))
  expect_true(all(is.finite(result$heat_after)))
  expect_true(all(is.finite(result$pre_temp_disp_norm)))
  expect_true(all(is.finite(result$attraction_disp_max_abs)))
  expect_true(all(is.finite(result$repulsion_disp_max_abs)))
  expect_true(all(result$attraction_edges_match))
  expect_true(all(result$repulsion_neighbors_match))
  expect_true(all(is.finite(terms$norm2_legacy)))
  expect_true(all(is.finite(terms$scale_legacy)))
  expect_true(any(term_summary$max_abs > 0))
  expect_true(any(summary$max_abs > 0))
  expect_true(any(source$term_trace$step_max_abs > 0))
})
