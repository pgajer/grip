test_that("GMDS app bundle builds a canonical MISF trace bundle", {
  compute_bundle <- getFromNamespace("gripui.gmds.compute.bundle", "grip")
  merge_values <- getFromNamespace(".gripui.family.merge.values", "grip")

  catalog <- gripui_graph_family_catalog()
  desc <- catalog$mesh
  values <- merge_values(desc, preset_id = "default")
  values$h <- 4L
  values$w <- 4L
  values$surface <- "paraboloid"
  values$amplitude <- 0.25

  bundle <- compute_bundle(
    desc = desc,
    values = values,
    dim = 3L,
    num_init = 6L,
    prepare_seed = 1101L,
    optimizer_seed = 2101L,
    top_level_max_iter = 1L,
    insertion_max_iter = 4L,
    refinement_max_iter = 1L,
    final_polish_max_iter = 1L,
    n_threads = 0L
  )

  expect_true(is.list(bundle))
  expect_equal(bundle$payload$family_id, "mesh")
  expect_s3_class(bundle$fit, "grip_misf_gmds_fit")
  expect_true(nrow(bundle$stage_trace) >= 3L)
  expect_true(length(bundle$stage_data) >= 3L)
  expect_true(bundle$prepared$top_level_level >= 0L)
  expect_true(all(c("seed", "top_level", "final_polish") %in% bundle$stage_trace$stage))
  expect_true(is.list(bundle$stage_payloads))
  expect_true(all(c("seed", "initial_placement", "top_level") %in% names(bundle$stage_payloads)))
  expect_true(length(bundle$stage_payloads$seed$active_vertices) >= 3L)
  expect_equal(bundle$stage_payloads$initial_placement$level, bundle$prepared$top_level_level)
  expect_equal(bundle$stage_payloads$top_level$level, bundle$prepared$top_level_level)

  expansion_levels <- getFromNamespace("gripui.gmds.expansion.levels", "grip")(bundle)
  expect_true(length(expansion_levels) >= 1L)

  level_stage_payload <- getFromNamespace("gripui.gmds.level.stage.payload", "grip")
  insertion_payload <- level_stage_payload(bundle, "insertion", expansion_levels[[1L]])
  refinement_payload <- level_stage_payload(bundle, "refinement", expansion_levels[[1L]])
  expect_false(is.null(insertion_payload))
  expect_false(is.null(refinement_payload))
  expect_equal(insertion_payload$level, expansion_levels[[1L]])
  expect_equal(refinement_payload$level, expansion_levels[[1L]])
})

test_that("GMDS stage explorer app builds", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("rgl")

  old <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old), add = TRUE)

  app <- gripui_gmds_app(catalog = gripui_graph_family_catalog()[c("mesh", "sampled_rectangle")])
  expect_s3_class(app, "shiny.appobj")
})
