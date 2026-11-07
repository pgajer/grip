test_that("coordinates survive extraction without transformation", {
  e <- edges.path(6)
  z <- grip(e, n = 6, dim = 2, seed = 11)
  dimnames(z) <- list(letters[1:6], c("horizontal", "vertical"))
  attr(z, "units") <- "drawing units"
  expect_identical(layout.coords(z), z)
  trace <- trace.grip(e, n = 6, dim = 2, seed = 11)
  fit <- classical.mds(edges = e, n = 6)
  refinement <- geodesic.kk(z, edges = e, n = 6, max.iter = 0)
  landmark <- landmark.geodesic.kk(z, edges = e, n = 6, max.iter = 0)
  for (obj in list(trace, fit, refinement, landmark)) {
    expected <- if (inherits(obj, "grip_layout_trace")) obj$final else obj$coords
    expect_identical(layout.coords(obj), expected)
    expect_identical(score.layout(layout.coords(obj), edges = e, n = 6),
                     score.layout(expected, edges = e, n = 6))
  }
  expect_identical(layout.coords(matrix(1:20, 5, 4)), matrix(1:20, 5, 4))
})

test_that("ambiguous and malformed objects cannot silently choose coordinates", {
  z <- matrix(1:12, 6, 2)
  expect_error(layout.coords(list(coords = z)), "Unsupported")
  expect_error(layout.coords(structure(list(coords = z, final = z), class = "grip_gmds_layout")), "Ambiguous")
  expect_error(layout.coords(structure(list(coords = z[, 1, drop = FALSE]), class = "grip_gmds_layout")), "at least two")
  expect_error(layout.coords(structure(list(final = z * NA), class = "grip_layout_trace")), "finite")
  expect_error(layout.coords(structure(list(), class = "grip_layout_trace")), "Malformed")
  expect_error(layout.coords(structure(list(coords = z, prepared = list(n = 4)), class = "grip_gmds_layout")), "rows")
  expect_error(layout.coords(1:4), "Unsupported")
})
