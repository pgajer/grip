test_that("graph family catalog covers the implemented explorer families", {
  catalog <- gripui_graph_family_catalog()

  expected <- c(
    "mesh",
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
