test_that("estimates distinguish dense distances from sparse retained routes", {
  full <- estimate.preparation(10000, 20000)
  expect_equal(full$dense.distance.bytes.lower.bound, 800000000)
  expect_equal(full$pair.count.upper.bound, 49995000)
  expect_equal(estimate.preparation(10000, pair_mode = "landmark_sparse")$dense.distance.bytes.lower.bound,
               full$dense.distance.bytes.lower.bound)
  expect_equal(estimate.preparation(10000, pair_mode = "edge_only")$dense.distance.bytes.lower.bound, 0)
  expect_gt(as.numeric(object.size(matrix(0, 10, 10))),
            estimate.preparation(10)$dense.distance.bytes.lower.bound)
  expect_error(estimate.preparation(4.5), "integer")
  expect_error(estimate.preparation(4, n.edges = 3.5), "integer")
})

test_that("large-work warning precedes graph searches and has an explicit override", {
  old <- options(grip.preparation.warn.bytes = 1)
  on.exit(options(old), add = TRUE)
  local_mocked_bindings(grip.shortest.path.tree = function(...) stop("SEARCH"))
  events <- character()
  tryCatch(withCallingHandlers(prepare.geodesic.kk(edges.path(4), n = 4),
    warning = function(w) { events <<- c(events, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) events <<- c(events, conditionMessage(e)))
  expect_match(events[[1]], "before routes, copies, and workspace")
  expect_identical(events[[2]], "SEARCH")
  options(grip.preparation.warn.bytes = Inf)
  expect_error(prepare.geodesic.kk(edges.path(4), n = 4), "SEARCH")
  expect_silent(prepare.edge.kk(edges.path(4), n = 4))
})
