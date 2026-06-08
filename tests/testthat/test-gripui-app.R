test_that("gripui_app constructs a Shiny app object", {
  old.rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old.rgl), add = TRUE)

  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("htmltools")
  skip_if_not_installed("rgl")

  graph <- list(
    adj_list = list(c(2L, 3L), c(1L, 4L), c(1L, 4L), c(2L, 3L)),
    vertex_data = data.frame(
      cluster = c("A", "A", "B", "B"),
      stringsAsFactors = FALSE
    )
  )
  layouts <- data.frame(
    candidate = "toy.3d",
    stage = "layout",
    seed = 1L,
    status = "ok",
    stringsAsFactors = FALSE
  )
  layouts$coords <- list(
    matrix(c(
      0, 0, 0,
      1, 0, 0,
      0, 1, 0,
      1, 1, 0
    ), ncol = 3, byrow = TRUE)
  )

  project <- gripui_project(graph = graph, layouts = layouts, title = "Toy app")
  app <- gripui_app(project)

  expect_s3_class(app, "shiny.appobj")
  expect_true(is.function(app$httpHandler))
  expect_true(is.function(app$serverFuncSource))
})

test_that("gripui_app constructs from compare.layouts output", {
  old.rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old.rgl), add = TRUE)

  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("htmltools")
  skip_if_not_installed("rgl")

  edges <- edges.path(5)
  cmp <- compare.layouts(
    edges = edges,
    n = 5,
    candidates = "default",
    seeds = 1L,
    return.layouts = TRUE
  )
  graph <- list(adj_list = list(2L, c(1L, 3L), c(2L, 4L), c(3L, 5L), 4L))

  project <- gripui_project_from_compare(cmp, graph = graph, title = "Compare app")
  app <- gripui_app(project)

  expect_s3_class(app, "shiny.appobj")
  expect_true(is.function(app$httpHandler))
  expect_true(is.function(app$serverFuncSource))
})

test_that("run_gripui can auto-stop for automated checks", {
  old.rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old.rgl), add = TRUE)

  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("htmltools")
  skip_if_not_installed("later")
  skip_if_not_installed("rgl")

  graph <- list(adj_list = list(2L, c(1L, 3L), 2L))
  layouts <- data.frame(
    candidate = "toy.layout",
    stage = "layout",
    seed = 1L,
    status = "ok",
    stringsAsFactors = FALSE
  )
  project <- gripui_project(graph = graph, layouts = layouts, title = "Toy project")

  expect_invisible(
    run_gripui(
      project,
      launch.browser = FALSE,
      quiet = TRUE,
      auto.stop.after = 0.1
    )
  )
})
