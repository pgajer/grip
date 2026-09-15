test_that("static region overrides and positional edges work", {
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({ grDevices::dev.off(); unlink(file) }, add = TRUE)
  xy <- cbind(0:3, c(0, 1, 0, 1))
  edges <- edges.path(4)
  expect_invisible(plot.layout(xy, edges, asp = 2, xlim = c(-5, 5),
                              xlab = "Horizontal", col = "blue"))
  expect_invisible(plot.layout(cbind(xy, 0:3), edges, projection = "ortho",
                              asp = 2, xlab = "Projected coordinate",
                              xlim = c(-5, 5), ylim = c(-4, 4), pch = 2, cex = 2))
})

test_that("orthographic options reach the vertex renderer", {
  calls <- list()
  local_mocked_bindings(points = function(...) calls[[length(calls) + 1L]] <<- list(...),
                       .package = "graphics")
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({ grDevices::dev.off(); unlink(file) }, add = TRUE)
  xyz <- cbind(0:3, c(0, 1, 0, 1), 0:3)
  plot.layout(xyz, projection = "ortho", pch = 2, cex = 2, col = "red", bg = "white")
  expect_equal(calls[[1]]$pch, 2)
  expect_equal(calls[[1]]$cex, 2)
  expect_identical(calls[[1]]$col, "red")
  expect_identical(calls[[1]]$bg, "white")
  plot.layout(xyz, projection = "ortho")
  expect_equal(calls[[2]]$pch, 16)
  expect_equal(calls[[2]]$cex, 0.55)
})
