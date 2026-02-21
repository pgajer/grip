test_that("basic layout returns finite matrix", {
  edges <- edges_path(10)
  coords <- grip.layout(edges, n = 10, dim = 2,
                        engine = "mish_v5",
                        placement = "barycenter",
                        rounds = 5, final_rounds = 3,
                        num_init = 5, num_nbrs = 6,
                        seed = 123)
  expect_equal(dim(coords), c(10, 2))
  expect_true(all(is.finite(coords)))
})

test_that("seeded runs are deterministic", {
  edges <- edges_cycle(12)
  coords1 <- grip.layout(edges, n = 12, dim = 2,
                         engine = "mish_v5",
                         placement = "barycenter",
                         rounds = 4, final_rounds = 2,
                         num_init = 4, num_nbrs = 5,
                         seed = 42)
  coords2 <- grip.layout(edges, n = 12, dim = 2,
                         engine = "mish_v5",
                         placement = "barycenter",
                         rounds = 4, final_rounds = 2,
                         num_init = 4, num_nbrs = 5,
                         seed = 42)
  expect_identical(coords1, coords2)
})

test_that("mish_v6 runs and returns expected shape", {
  edges <- edges_cycle(15)
  coords <- grip.layout(edges, n = 15, dim = 2,
                        engine = "mish_v6",
                        placement = "barycenter",
                        rounds = 5, final_rounds = 3,
                        num_init = 5, num_nbrs = 6,
                        seed = 7)
  expect_equal(dim(coords), c(15, 2))
  expect_true(all(is.finite(coords)))
})

test_that("circle placement works in 2D", {
  edges <- edges_path(8)
  coords <- grip.layout(edges, n = 8, dim = 2,
                        engine = "mish_v5",
                        placement = "circle",
                        rounds = 4, final_rounds = 2,
                        num_init = 4, num_nbrs = 5,
                        seed = 9)
  expect_equal(dim(coords), c(8, 2))
  expect_true(all(is.finite(coords)))
})

test_that("circle placement falls back in 3D with warning", {
  edges <- edges_mesh(4, 4)
  expect_warning({
    coords <- grip.layout(edges, n = 16, dim = 3,
                          engine = "mish_v5",
                          placement = "circle",
                          rounds = 4, final_rounds = 2,
                          num_init = 5, num_nbrs = 6,
                          seed = 11)
    expect_equal(dim(coords), c(16, 3))
    expect_true(all(is.finite(coords)))
  }, "circle placement is only used for 2D")
})

test_that("tree example runs", {
  edges <- edges_kary_tree(k = 2, depth = 2)
  n <- max(edges)
  coords <- grip.layout(edges, n = n, dim = 2,
                        engine = "mish_v6",
                        placement = "barycenter",
                        rounds = 4, final_rounds = 2,
                        num_init = 4, num_nbrs = 5,
                        seed = 5)
  expect_equal(dim(coords), c(n, 2))
  expect_true(all(is.finite(coords)))
})

test_that("adj_list + weight_list input works", {
  adj_list <- list(c(2), c(1, 3), c(2, 4), c(3))
  weight_list <- list(c(1.0), c(1.0, 2.0), c(2.0, 1.5), c(1.5))
  coords <- grip.layout(adj_list = adj_list,
                        weight_list = weight_list,
                        n = 4,
                        dim = 2,
                        engine = "mish_v5",
                        placement = "barycenter",
                        rounds = 4, final_rounds = 2,
                        num_init = 3, num_nbrs = 3,
                        seed = 12)
  expect_equal(dim(coords), c(4, 2))
  expect_true(all(is.finite(coords)))
})

test_that("adj_list input works without weights", {
  adj_list <- list(c(2), c(1, 3), c(2, 4), c(3))
  coords <- grip.layout(adj_list = adj_list,
                        n = 4,
                        dim = 2,
                        engine = "mish_v5",
                        placement = "barycenter",
                        rounds = 4, final_rounds = 2,
                        num_init = 3, num_nbrs = 3,
                        seed = 13)
  expect_equal(dim(coords), c(4, 2))
  expect_true(all(is.finite(coords)))
})

test_that("disconnected graph defaults to safe component layouts", {
  edges <- rbind(
    cbind(1:4, 2:5),
    cbind(7:9, 8:10)
  )
  expect_warning({
    coords <- grip.layout(edges = edges,
                          n = 10,
                          dim = 3,
                          engine = "mish_v5",
                          placement = "barycenter",
                          rounds = 4,
                          final_rounds = 2,
                          num_init = 4,
                          num_nbrs = 5,
                          seed = 21)
    expect_equal(dim(coords), c(10, 3))
    expect_true(all(is.finite(coords)))
  }, "connected components")
})

test_that("disconnected mode can be set to error", {
  edges <- rbind(
    cbind(1:2, 2:3),
    cbind(5:6, 6:7)
  )
  expect_error(
    grip.layout(edges = edges,
                n = 7,
                dim = 2,
                disconnected = "error",
                seed = 7),
    "connected components"
  )
})

test_that("disconnected adj_list + weights is handled safely", {
  adj_list <- list(
    c(2), c(1),           # component 1
    c(4), c(3),           # component 2
    integer(0)            # isolated
  )
  weight_list <- list(
    c(1.0), c(1.0),
    c(2.0), c(2.0),
    numeric(0)
  )
  expect_warning({
    coords <- grip.layout(adj_list = adj_list,
                          weight_list = weight_list,
                          n = 5,
                          dim = 2,
                          engine = "mish_v6",
                          placement = "barycenter",
                          rounds = 4,
                          final_rounds = 2,
                          num_init = 3,
                          num_nbrs = 3,
                          seed = 31)
    expect_equal(dim(coords), c(5, 2))
    expect_true(all(is.finite(coords)))
  }, "connected components")
})
