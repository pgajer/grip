test_that("public functions and parameters use the package naming convention", {
  exports <- getNamespaceExports("grip")
  public <- exports[vapply(exports, function(nm) is.function(getExportedValue("grip", nm)), logical(1))]
  expect_gt(length(public), 100L)
  expect_true(all(grepl("^[a-z][a-z0-9]*(\\.[a-z0-9]+)*$", public)))
  for (nm in public) {
    arguments <- names(formals(getExportedValue("grip", nm)))
    exceptions <- c("...", "X", "scale.L0")
    expect_true(all(arguments %in% exceptions |
      grepl("^[a-z][a-z0-9]*(\\.[a-z0-9]+)*$", arguments)), info = nm)
  }
  # These methods must retain the names of the established S3 classes.
  expect_true(is.function(getS3method("print", "grip_gmds_layout")))
  expect_true(is.function(getS3method("print", "grip_misf")))
})

test_that("obsolete public names are rejected rather than aliased", {
  old.functions <- c("gripui_app", "gripui_family_app", "gripui_graph_family_catalog",
    "gripui_project", "gripui_project_from_compare", "gripui_project_from_dir",
    "gripui_validate_project", "run_gripui", "run_gripui_family")
  expect_false(any(old.functions %in% getNamespaceExports("grip")))
  expect_true(all(gsub("_", ".", old.functions, fixed = TRUE) %in% getNamespaceExports("grip")))
  expect_error(metric.mds(max_iter = 1), "unused argument")
  expect_error(metric.mds(pair_weights = "uniform"), "unused argument")
  expect_error(grip(edge_weights = 1), "unused argument")
  expect_error(gripui.app(layout_source = "x"), "unused argument")
  expect_error(metric.mds(edges = edges.path(4), sgd.control = list(learning_rate = .5)),
    "documented controls")
  expect_error(metric.mds(edges = edges.path(4), approximation = "sparse",
    sparse.control = list(n_pivots = 2)), "unique names")
  expect_error(compare.layouts(edges = edges.path(4), candidates = list(old = list(num_init = 3))),
    "unsupported layout arguments")
  expect_error(compare.layouts(edges = edges.path(4), search = list(final_rounds = 2)),
    "unsupported layout arguments")
  expect_error(edge.kk(edges = edges.path(4), init = "weighted_grip",
    weighted.grip.args = list(final_rounds = 2)), "unknown arguments")
})

test_that("new controls work and saved graph fields keep their schema", {
  graph <- mesh.surface.graph(3, 3, surface = "saddle")
  expect_true(all(c("edges", "edge_weights", "coords_surface") %in% names(graph)))
  fit <- suppressWarnings(metric.mds(edges = graph$edges, edge.weights = graph$edge_weights,
    dim = 3, max.iter = 2, sgd.control = list(learning.rate = .4, checkpoint.every = 1)))
  expect_equal(dim(fit$coords), c(9L, 3L))
  expect_equal(fit$metadata$settings$sgd_control$learning.rate, .4)
  controls <- grip.mds.sgd.control(list(), 2L)
  expect_false(any(grepl("_", names(controls), fixed = TRUE)))
  expect_false(any(grepl("_", names(params.from.summary(data.frame(
    preset = NA_character_, placement = "barycenter", rounds = 2L, final.rounds = 3L, num.init = 4L,
    num.nbrs = 2L, r = .1, s = 3, repulsion.factor = 1, tinit.factor = 6))), fixed = TRUE)))
})
