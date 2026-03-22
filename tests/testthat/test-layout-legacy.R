test_that("mesh example runs", {
  edges <- edges.mesh(5, 5)
  coords <- grip.layout(edges, n = 25, dim = 3,
                        placement = "barycenter",
                        rounds = 5, final_rounds = 3,
                        num_init = 6, num_nbrs = 8,
                        seed = 101)
  expect_equal(dim(coords), c(25, 3))
  expect_true(all(is.finite(coords)))
})

test_that("cylinder example runs", {
  edges <- edges.cylinder(4, 6)
  coords <- grip.layout(edges, n = 24, dim = 3,
                        placement = "barycenter",
                        rounds = 5, final_rounds = 3,
                        num_init = 6, num_nbrs = 8,
                        seed = 202)
  expect_equal(dim(coords), c(24, 3))
  expect_true(all(is.finite(coords)))
})

test_that("torus example runs", {
  edges <- edges.torus(4, 4)
  coords <- grip.layout(edges, n = 16, dim = 3,
                        placement = "barycenter",
                        rounds = 5, final_rounds = 3,
                        num_init = 5, num_nbrs = 7,
                        seed = 303)
  expect_equal(dim(coords), c(16, 3))
  expect_true(all(is.finite(coords)))
})

test_that("sierpinski example runs", {
  edges <- edges.sierpinski.triangle(2)
  n <- max(edges)
  coords <- grip.layout(edges, n = n, dim = 2,
                        placement = "circle",
                        rounds = 5, final_rounds = 3,
                        num_init = 5, num_nbrs = 7,
                        seed = 404)
  expect_equal(dim(coords), c(n, 2))
  expect_true(all(is.finite(coords)))
})
