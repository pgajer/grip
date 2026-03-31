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
