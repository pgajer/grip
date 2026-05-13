test_that("edge-repulsive R gradient matches finite differences", {
  grad <- grip:::.grip.edge.repulsive.finite.difference.gradient.check()
  expect_true(is.finite(grad[["abs_error"]]))
  expect_lt(grad[["abs_error"]], 1e-5)
})

test_that("edge-repulsive objective recognizes exact edge lengths", {
  Z <- rbind(c(0, 0, 0), c(1, 0, 0), c(2, 0, 0))
  edges <- rbind(c(1L, 2L), c(2L, 3L))
  st <- grip.edge.repulsive.state(Z, edges, c(1, 1), lambda = 0)
  expect_lt(abs(st$edge.energy), 1e-8)
  expect_lt(st$gradient.norm, 1e-6)
})

test_that("upper barrier reports infeasible stretched edges", {
  Z <- rbind(c(0, 0, 0), c(1, 0, 0), c(2, 0, 0))
  edges <- rbind(c(1L, 2L), c(2L, 3L))
  bad <- grip.edge.repulsive.state(
    Z * 2, edges, c(1, 1),
    edge.family = "upper_barrier",
    eps.plus = 0.25,
    beta = 1
  )
  expect_false(bad$feasible)
  expect_equal(bad$n.wall.violations, 2L)
})

test_that("C++ backend matches R reference for edge-repulsive state and optimizer", {
  skip_if_not(grip:::.grip.edge.repulsive.cpp.available(), "C++ backend not loaded")
  set.seed(12)
  Z <- matrix(rnorm(18), ncol = 3)
  edges <- rbind(c(1L, 2L), c(1L, 4L), c(2L, 3L), c(3L, 5L), c(5L, 6L))
  edge.lengths <- c(0.8, 1.1, 0.9, 1.3, 0.7)
  edge.weights <- c(1.0, 1.3, 0.8, 1.1, 1.4)
  pair.index <- rbind(c(1L, 3L), c(1L, 6L), c(2L, 5L), c(4L, 6L))
  pair.weights <- c(0.7, 1.2, 0.9, 1.5)

  st.r <- grip.edge.repulsive.state(
    Z, edges, edge.lengths,
    edge.weights = edge.weights,
    lambda = 0.025,
    pair.index = pair.index,
    pair.weights = pair.weights,
    repulsion.family = "log",
    engine = "R"
  )
  st.cpp <- grip.edge.repulsive.state(
    Z, edges, edge.lengths,
    edge.weights = edge.weights,
    lambda = 0.025,
    pair.index = pair.index,
    pair.weights = pair.weights,
    repulsion.family = "log",
    engine = "cpp"
  )
  expect_equal(st.cpp$energy, st.r$energy, tolerance = 1e-10)
  expect_equal(st.cpp$edge.energy, st.r$edge.energy, tolerance = 1e-10)
  expect_equal(st.cpp$repel.energy, st.r$repel.energy, tolerance = 1e-10)
  expect_lt(max(abs(st.cpp$gradient - st.r$gradient)), 1e-10)
  expect_lt(max(abs(st.cpp$edge.embedded.lengths - st.r$edge.embedded.lengths)), 1e-10)

  st.r.inv <- grip.edge.repulsive.state(
    Z, edges, edge.lengths,
    edge.weights = edge.weights,
    lambda = 0.01,
    pair.index = pair.index,
    pair.weights = pair.weights,
    repulsion.family = "inverse_power",
    repulsion.power = 1.5,
    engine = "R"
  )
  st.cpp.inv <- grip.edge.repulsive.state(
    Z, edges, edge.lengths,
    edge.weights = edge.weights,
    lambda = 0.01,
    pair.index = pair.index,
    pair.weights = pair.weights,
    repulsion.family = "inverse_power",
    repulsion.power = 1.5,
    engine = "cpp"
  )
  expect_equal(st.cpp.inv$energy, st.r.inv$energy, tolerance = 1e-10)
  expect_lt(max(abs(st.cpp.inv$gradient - st.r.inv$gradient)), 1e-10)

  bad.r <- grip.edge.repulsive.state(
    Z * 4, edges, edge.lengths,
    edge.family = "upper_barrier",
    eps.plus = 0.2,
    beta = 0.5,
    lambda = 0,
    pair.index = pair.index,
    engine = "R"
  )
  bad.cpp <- grip.edge.repulsive.state(
    Z * 4, edges, edge.lengths,
    edge.family = "upper_barrier",
    eps.plus = 0.2,
    beta = 0.5,
    lambda = 0,
    pair.index = pair.index,
    engine = "cpp"
  )
  expect_identical(bad.cpp$feasible, bad.r$feasible)
  expect_identical(as.integer(bad.cpp$n.wall.violations), as.integer(bad.r$n.wall.violations))
  expect_lt(max(abs(bad.cpp$gradient - bad.r$gradient)), 1e-10)

  fit.r <- grip.optimize.edge.repulsive.stage(
    Z, edges, edge.lengths,
    edge.weights = edge.weights,
    lambda = 0.01,
    pair.index = pair.index,
    pair.weights = pair.weights,
    max.iter = 6L,
    initial.step = 0.01,
    engine = "R"
  )
  fit.cpp <- grip.optimize.edge.repulsive.stage(
    Z, edges, edge.lengths,
    edge.weights = edge.weights,
    lambda = 0.01,
    pair.index = pair.index,
    pair.weights = pair.weights,
    max.iter = 6L,
    initial.step = 0.01,
    engine = "cpp"
  )
  expect_equal(fit.cpp$coords, fit.r$coords, tolerance = 1e-9)
  expect_equal(fit.cpp$state$energy, fit.r$state$energy, tolerance = 1e-10)
  expect_gte(nrow(fit.cpp$trace), 1L)
})
