test_that("graph family catalog covers the implemented explorer families", {
  catalog <- gripui_graph_family_catalog()

  expected <- c(
    "mesh",
    "irregular_rectangle",
    "sampled_rectangle",
    "cylinder",
    "torus",
    "sphere",
    "recursive_mask_grid",
    "sierpinski_carpet",
    "vicsek",
    "occupied_mesh",
    "recursive_triangle_mask",
    "sierpinski_triangle",
    "recursive_tetrahedron_mask",
    "sierpinski_tetrahedron",
    "recursive_cube_mask",
    "menger_sponge",
    "cube_periodic_tunnels",
    "cube_asymmetric_cavities",
    "cube_channel_network",
    "triangulated_polyhedron",
    "triangulated_annulus",
    "triangulated_pair_of_pants",
    "irregular_annulus",
    "irregular_sphere",
    "irregular_pair_of_pants",
    "irregular_torus",
    "irregular_double_torus",
    "irregular_ball",
    "irregular_shell",
    "kary_tree"
  )

  expect_true(is.list(catalog))
  expect_setequal(names(catalog), expected)
})

test_that("default catalog families build normalized payloads", {
  catalog <- gripui_graph_family_catalog()

  for (id in names(catalog)) {
    desc <- catalog[[id]]
    values <- grip:::.gripui.family.param.defaults(desc)
    payload <- grip:::gripui.family.build.payload(desc, values)

    expect_identical(payload$family_id, desc$id, info = id)
    expect_identical(nrow(payload$coords_display), payload$n, info = id)
    expect_identical(nrow(payload$coords_plot), payload$n, info = id)
    expect_equal(length(payload$edge_weights), nrow(payload$edges), info = id)
    expect_true(is.list(payload$graph$adj_list), info = id)
    expect_true(is.character(payload$code) && nzchar(payload$code), info = id)
  }
})

test_that("generic family code paths show helper calls without quoting them", {
  catalog <- gripui_graph_family_catalog()
  ids <- c(
    "recursive_mask_grid",
    "occupied_mesh",
    "recursive_triangle_mask",
    "recursive_tetrahedron_mask",
    "recursive_cube_mask"
  )

  for (id in ids) {
    desc <- catalog[[id]]
    values <- grip:::.gripui.family.param.defaults(desc)
    code <- desc$code(values)
    expect_false(grepl('mask = "', code, fixed = TRUE), info = id)
    expect_false(grepl('keep = "', code, fixed = TRUE), info = id)
  }
})

test_that("compare helpers build shared color choices and summaries", {
  catalog <- gripui_graph_family_catalog()

  desc1 <- catalog$mesh
  desc2 <- catalog$kary_tree
  payload1 <- grip:::gripui.family.build.payload(desc1, grip:::.gripui.family.param.defaults(desc1))
  payload2 <- grip:::gripui.family.build.payload(desc2, grip:::.gripui.family.param.defaults(desc2))

  payload1$compare_slot <- 1L
  payload1$compare_label <- "Current"
  payload1$compare_preset <- "current_controls"

  payload2$compare_slot <- 2L
  payload2$compare_label <- "Variant 2"
  payload2$compare_preset <- "default"

  payloads <- list(payload1, payload2)
  choices <- grip:::gripui.family.compare.color.choices(payloads)
  summary <- grip:::gripui.family.compare.summary(payloads)

  expect_true("degree" %in% unname(choices))
  expect_true("surface_z" %in% unname(choices))
  expect_identical(nrow(summary), 2L)
  expect_identical(summary$label[[1L]], "Current")
  expect_identical(summary$preset[[2L]], "default")
})

test_that("compare slot selection uses current slot and valid preset fallbacks", {
  catalog <- gripui_graph_family_catalog()
  mock_input <- list(
    compare_family_2 = "mesh",
    compare_preset_2 = "not_a_real_preset"
  )

  current <- grip:::gripui.family.compare.slot.selection(
    idx = 1L,
    input = mock_input,
    catalog = catalog,
    current_family_id = "mesh",
    include_current = TRUE,
    lock_family = FALSE
  )
  fallback <- grip:::gripui.family.compare.slot.selection(
    idx = 2L,
    input = mock_input,
    catalog = catalog,
    current_family_id = "mesh",
    include_current = TRUE,
    lock_family = FALSE
  )

  expect_identical(current$source, "current")
  expect_identical(fallback$family_id, "mesh")
  expect_true(fallback$preset_id %in% c("default", names(catalog$mesh$presets)))
})

test_that("initial explore state keeps family and category aligned", {
  catalog <- gripui_graph_family_catalog()
  initial <- grip:::gripui.family.initial.explore.state(catalog)

  expect_identical(initial$family_id, "mesh")
  expect_identical(initial$category, catalog[[initial$family_id]]$category)
  expect_true(initial$family_id %in% unname(initial$choices))
})

test_that("sampled rectangle family is marked stochastic", {
  catalog <- gripui_graph_family_catalog()

  expect_true(grip:::gripui.family.is.stochastic(catalog$sampled_rectangle))
  expect_false(grip:::gripui.family.is.stochastic(catalog$mesh))
  expect_identical(grip:::gripui.family.seed.spec(catalog$sampled_rectangle)$id, "seed")
})

test_that("resample seed helper always returns a valid nonnegative integer", {
  seed <- grip:::gripui.family.resample.seed(current_seed = 1L, max_seed = 1000000L)
  repeated <- grip:::gripui.family.resample.seed(current_seed = seed, max_seed = 1000000L)
  zero_seed <- grip:::gripui.family.resample.seed(current_seed = 0L, max_seed = 0L)

  expect_true(is.integer(seed))
  expect_length(seed, 1L)
  expect_false(is.na(seed))
  expect_gte(seed, 0L)
  expect_lte(seed, 1000000L)
  expect_true(is.integer(repeated))
  expect_false(is.na(repeated))
  expect_gte(repeated, 0L)
  expect_lte(repeated, 1000000L)
  expect_identical(zero_seed, 0L)
})

test_that("family save helpers create readable paths and complete bundles", {
  catalog <- gripui_graph_family_catalog()
  desc <- catalog$sampled_rectangle
  values <- grip:::.gripui.family.param.defaults(desc)
  payload <- grip:::gripui.family.build.payload(desc, values)
  root <- tempfile("gripui-family-save-")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  stamp <- as.POSIXct("2026-04-06 12:34:56", tz = "America/New_York")

  path <- grip:::gripui.family.save.path(payload, root = root, timestamp = stamp)
  saved <- grip:::gripui.family.save.bundle(payload, path, saved_at = stamp)
  bundle <- readRDS(path)
  unique_path <- grip:::gripui.family.unique.save.path(path)

  expect_identical(saved, path)
  expect_true(file.exists(path))
  expect_match(path, "tmp/gripui-family-graphs/sampled_rectangle", fixed = TRUE)
  expect_match(basename(path), "sampled_rectangle__n-80__k-6__seed-1")
  expect_lte(nchar(basename(path), type = "bytes"), 120L)
  expect_match(basename(path), "20260406-123456\\.rds$")
  expect_match(basename(unique_path), "__01\\.rds$")
  expect_identical(bundle$family_id, payload$family_id)
  expect_identical(bundle$values, payload$values)
  expect_identical(bundle$code, payload$code)
  expect_identical(bundle$payload$raw$edges, payload$raw$edges)
})

test_that("saved bundle helpers list files newest-first and expose load modes", {
  catalog <- gripui_graph_family_catalog()
  desc <- catalog$sampled_rectangle
  payload <- grip:::gripui.family.build.payload(desc, grip:::.gripui.family.param.defaults(desc))
  root <- tempfile("gripui-family-bundles-")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)

  path1 <- grip:::gripui.family.save.path(
    payload,
    root = root,
    timestamp = as.POSIXct("2026-04-06 08:00:00", tz = "America/New_York")
  )
  path2 <- grip:::gripui.family.save.path(
    payload,
    root = root,
    timestamp = as.POSIXct("2026-04-06 09:00:00", tz = "America/New_York")
  )
  grip:::gripui.family.save.bundle(payload, path1)
  grip:::gripui.family.save.bundle(payload, path2)
  Sys.setFileTime(
    c(path1, path2),
    as.POSIXct("2026-04-06 10:00:00", tz = "America/New_York")
  )

  files <- grip:::gripui.family.saved.bundle.files(desc$id, root = root)
  choices <- grip:::gripui.family.saved.bundle.choices(desc, root = root)
  modes <- grip:::gripui.family.saved.load.mode.choices(desc)

  expect_identical(files, c(path2, path1))
  expect_identical(unname(choices), c(path2, path1))
  expect_true(all(nzchar(names(choices))))
  expect_identical(unname(modes), c("exact", "sample_topology", "sample_rebuild"))
  expect_match(grip:::gripui.family.saved.load.mode.help(desc), "Rebuild iKNN", fixed = TRUE)
})

test_that("sample-topology loads keep saved topology controls while reusing the sample", {
  catalog <- gripui_graph_family_catalog()
  desc <- catalog$sampled_rectangle
  saved_values <- grip:::.gripui.family.param.defaults(desc)
  bundle <- list(
    family_id = desc$id,
    values = saved_values,
    payload = list(
      raw = desc$builder(saved_values)
    )
  )
  current_values <- saved_values
  current_values$k <- 11L
  current_values$graph_space <- "param"
  current_values$surface <- "saddle"

  loaded <- grip:::gripui.family.sampled.rectangle.raw.from.bundle(
    desc = desc,
    current_values = current_values,
    bundle = bundle,
    mode = "sample_topology"
  )

  expect_identical(loaded$values$k, saved_values$k)
  expect_identical(loaded$values$graph_space, saved_values$graph_space)
  expect_identical(loaded$values$surface, "saddle")
  expect_equal(loaded$raw$edges, bundle$payload$raw$edges)
})

test_that("family explorer app constructs when optional packages are available", {
  old <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old), add = TRUE)

  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("rgl")

  app <- gripui_family_app()
  expect_s3_class(app, "shiny.appobj")
})
