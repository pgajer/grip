edge_gkk_barrier_fixture <- function(kind) {
  switch(kind,
    triangle = list(
      coords = matrix(c(
        0.00, 0.00,
        1.25, 0.25,
        0.35, 1.40
      ), ncol = 2, byrow = TRUE),
      edges = rbind(c(1L, 2L), c(2L, 3L), c(1L, 3L)),
      ell = c(1.00, 1.35, 1.10),
      k = c(0.8, 1.7, 1.2)
    ),
    path = list(
      coords = matrix(c(
        0.0, 0.0, 0.1,
        1.2, 0.3, -0.1,
        2.1, -0.4, 0.2,
        3.4, 0.1, 0.5,
        4.2, -0.2, -0.3
      ), ncol = 3, byrow = TRUE),
      edges = rbind(c(1L, 2L), c(2L, 3L), c(3L, 4L), c(4L, 5L)),
      ell = c(0.9, 1.4, 1.1, 1.6),
      k = c(1.0, 0.7, 1.8, 1.2)
    ),
    irregular = list(
      coords = matrix(c(
        0.0, 0.0, 0.0,
        1.0, 0.2, 0.1,
        0.4, 1.2, -0.2,
        1.5, 1.0, 0.4,
        2.1, 0.1, -0.3,
        2.3, 1.4, 0.2
      ), ncol = 3, byrow = TRUE),
      edges = rbind(
        c(1L, 2L), c(1L, 3L), c(2L, 3L), c(2L, 4L),
        c(3L, 4L), c(4L, 5L), c(4L, 6L), c(5L, 6L)
      ),
      ell = c(1.0, 1.15, 1.3, 0.95, 1.05, 1.4, 1.1, 1.25),
      k = c(0.9, 1.1, 0.6, 1.8, 1.4, 0.75, 1.25, 1.6)
    ),
    paraboloid = {
      x <- matrix(c(
        -0.85, -0.35,
        -0.35, 0.60,
        0.25, -0.75,
        0.70, 0.15,
        0.95, 0.80,
        -0.10, 0.05,
        0.45, 0.55
      ), ncol = 2, byrow = TRUE)
      coords <- cbind(x, x[, 1]^2 + x[, 2]^2)
      edges <- rbind(
        c(1L, 2L), c(1L, 6L), c(2L, 6L), c(2L, 7L),
        c(3L, 4L), c(3L, 6L), c(4L, 6L), c(4L, 7L),
        c(5L, 7L), c(6L, 7L)
      )
      edge.d <- sqrt(rowSums((coords[edges[, 1L], , drop = FALSE] -
        coords[edges[, 2L], , drop = FALSE])^2))
      list(
        coords = coords,
        edges = edges,
        ell = edge.d * c(1.05, 0.92, 1.08, 0.97, 1.12, 0.90, 1.04, 1.09, 0.94, 1.02),
        k = seq(0.7, 1.6, length.out = nrow(edges))
      )
    },
    stop("unknown fixture")
  )
}

check_edge_gkk_barrier_parity <- function(fixture, distance.eps = 1e-10) {
  gkk <- grip:::grip.edge.isometric.energy.gradient(
    coords = fixture$coords,
    edges = fixture$edges,
    edge_weights = fixture$ell,
    stiffness = fixture$k,
    scale = 1,
    edge_length_epsilon = distance.eps
  )
  barrier <- edge.repulsive.state(
    coords = fixture$coords,
    edges = fixture$edges,
    edge.lengths = fixture$ell,
    edge.weights = fixture$k * fixture$ell^2,
    edge.family = "quadratic",
    lambda = 0,
    distance.eps = distance.eps,
    engine = "cpp"
  )
  barrier.r <- edge.repulsive.state(
    coords = fixture$coords,
    edges = fixture$edges,
    edge.lengths = fixture$ell,
    edge.weights = fixture$k * fixture$ell^2,
    edge.family = "quadratic",
    lambda = 0,
    distance.eps = distance.eps,
    engine = "R"
  )
  list(gkk = gkk, barrier = barrier, barrier.r = barrier.r)
}

test_that("edge-KK objective equals quadratic edge barrier with ell-squared weights", {
  for (kind in c("triangle", "path", "irregular", "paraboloid")) {
    fixture <- edge_gkk_barrier_fixture(kind)
    states <- check_edge_gkk_barrier_parity(fixture)

    expect_equal(states$barrier$energy, states$gkk$energy, tolerance = 1e-10,
                 info = kind)
    expect_equal(states$barrier$edge.energy, states$gkk$energy, tolerance = 1e-10,
                 info = kind)
    expect_equal(states$barrier$repel.energy, 0, tolerance = 1e-14,
                 info = kind)
    expect_equal(states$barrier.r$energy, states$gkk$energy, tolerance = 1e-10,
                 info = kind)
  }
})

test_that("edge-KK gradient equals quadratic edge-barrier gradient", {
  for (kind in c("triangle", "path", "irregular", "paraboloid")) {
    fixture <- edge_gkk_barrier_fixture(kind)
    states <- check_edge_gkk_barrier_parity(fixture)

    expect_equal(states$barrier$gradient, states$gkk$gradient, tolerance = 1e-10,
                 info = kind)
    expect_equal(states$barrier.r$gradient, states$gkk$gradient, tolerance = 1e-10,
                 info = kind)
    expect_equal(states$barrier$gradient.norm, states$gkk$gradient_norm,
                 tolerance = 1e-10, info = kind)
  }
})

test_that("ell-squared conversion is necessary for non-unit edge lengths", {
  fixture <- edge_gkk_barrier_fixture("irregular")
  gkk <- grip:::grip.edge.isometric.energy.gradient(
    coords = fixture$coords,
    edges = fixture$edges,
    edge_weights = fixture$ell,
    stiffness = fixture$k,
    scale = 1,
    edge_length_epsilon = 1e-10
  )
  relative.weight.barrier <- edge.repulsive.state(
    coords = fixture$coords,
    edges = fixture$edges,
    edge.lengths = fixture$ell,
    edge.weights = fixture$k,
    edge.family = "quadratic",
    lambda = 0,
    distance.eps = 1e-10,
    engine = "cpp"
  )

  expect_gt(abs(relative.weight.barrier$energy - gkk$energy), 1e-3)
  expect_gt(max(abs(relative.weight.barrier$gradient - gkk$gradient)), 1e-3)
})
