test_that("Riemannian star structure builds antipodal-weighted center pairs", {
  prepared <- grip.prepare.graph.geodesic.mds(
    edges = rbind(c(1L, 2L), c(2L, 3L)),
    n = 3L,
    edge_weights = c(1, 1)
  )
  X <- cbind(c(0, 1, 2), c(0, 0, 0))
  star <- graph.riemannian.star.structure(
    prepared = prepared,
    X = X,
    angle.power = 4,
    reliability = "length.balance"
  )

  expect_s3_class(star, "grip_riemannian_star")
  expect_equal(nrow(star$pairs), 1L)
  expect_equal(star$pairs$center, 2L)
  expect_equal(star$pairs$cos_angle, -1, tolerance = 1e-12)
  expect_equal(star$pairs$angle_weight, 1, tolerance = 1e-12)
})

test_that("Riemannian star angle.power downweights non-antipodal pairs", {
  prepared <- grip.prepare.graph.geodesic.mds(
    edges = rbind(c(1L, 2L), c(1L, 3L), c(1L, 4L)),
    n = 4L,
    edge_weights = c(1, 1, 1)
  )
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(-1, 0),
    c(0, 1)
  )
  star <- graph.riemannian.star.structure(prepared = prepared, X = X, angle.power = 8)

  antipodal <- star$pairs$angle_weight[star$pairs$cos_angle < -0.9]
  orthogonal <- star$pairs$angle_weight[abs(star$pairs$cos_angle) < 0.1]
  expect_equal(antipodal, 1, tolerance = 1e-12)
  expect_true(all(orthogonal < antipodal))
})

test_that("Riemannian star quantile filter keeps strongest star pairs", {
  prepared <- grip.prepare.graph.geodesic.mds(
    edges = rbind(c(1L, 2L), c(1L, 3L), c(1L, 4L), c(1L, 5L)),
    n = 5L,
    edge_weights = rep(1, 4)
  )
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(-1, 0),
    c(0, 1),
    c(0.7, 0.7)
  )

  full <- graph.riemannian.star.structure(prepared = prepared, X = X, angle.power = 2)
  filtered <- graph.riemannian.star.structure(
    prepared = prepared,
    X = X,
    angle.power = 2,
    star.quantile = 0.75
  )
  cutoff <- as.double(stats::quantile(full$pairs$angle_weight, 0.75, names = FALSE))

  expect_equal(filtered$star.quantile, 0.75)
  expect_lt(nrow(filtered$pairs), nrow(full$pairs))
  expect_true(all(filtered$pairs$angle_weight >= cutoff))
})

test_that("kernel Gram energy gradient matches finite differences", {
  edges <- rbind(c(1L, 2L), c(2L, 3L), c(2L, 4L))
  prepared <- grip.prepare.graph.geodesic.mds(edges = edges, n = 4L, edge_weights = c(1, 1, 1))
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(2, 0),
    c(1, 1)
  )
  star <- graph.riemannian.star.structure(prepared = prepared, X = X, angle.power = 2)
  coords <- matrix(c(
    0.1, 0.0,
    1.0, 0.2,
    1.8, 0.8,
    0.8, 1.2
  ), ncol = 2, byrow = TRUE)
  state <- grip:::grip.kernel.gram.energy.gradient(
    coords = coords,
    edges = prepared$edges,
    edge_weights = prepared$edge_targets,
    edge_stiffness = c(0.8, 1.2, 1.1),
    star = star,
    edge_scale = 1.05,
    lambda.edge = 0.7,
    lambda.gram = 0.4
  )

  eps <- 1e-6
  fd <- matrix(0, nrow(coords), ncol(coords))
  for (i in seq_len(nrow(coords))) {
    for (j in seq_len(ncol(coords))) {
      plus <- coords
      minus <- coords
      plus[i, j] <- plus[i, j] + eps
      minus[i, j] <- minus[i, j] - eps
      e.plus <- grip:::grip.kernel.gram.energy.gradient(
        plus, prepared$edges, prepared$edge_targets, c(0.8, 1.2, 1.1),
        star, edge_scale = 1.05, lambda.edge = 0.7, lambda.gram = 0.4
      )$energy
      e.minus <- grip:::grip.kernel.gram.energy.gradient(
        minus, prepared$edges, prepared$edge_targets, c(0.8, 1.2, 1.1),
        star, edge_scale = 1.05, lambda.edge = 0.7, lambda.gram = 0.4
      )$energy
      fd[i, j] <- (e.plus - e.minus) / (2 * eps)
    }
  }

  expect_equal(state$gradient, fd, tolerance = 1e-5)
})

test_that("C++ kernel Gram-gKK matches R reference engine", {
  prepared <- grip.prepare.graph.geodesic.mds(
    edges = rbind(c(1L, 2L), c(2L, 3L), c(2L, 4L), c(4L, 5L)),
    n = 5L,
    edge_weights = c(1, 1.1, 1, 1.2)
  )
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(2, 0),
    c(1, 1),
    c(1, 2)
  )
  start <- matrix(c(
    0.0, 0.2,
    0.9, 0.1,
    1.7, 0.7,
    0.7, 1.0,
    0.5, 1.9
  ), ncol = 2, byrow = TRUE)
  args <- list(
    coords = start,
    prepared = prepared,
    X = X,
    dim = 2L,
    angle.power = 4,
    lambda.edge = 1,
    lambda.gram = 0.5,
    stiffness_method = "uniform",
    density_mix = 1,
    scale_mode = "profiled",
    max_iter = 10L,
    initial_step = 0.05,
    return_trace = TRUE
  )
  fit.cpp <- do.call(grip.optimize.kernel.gram.gkk.layout, c(args, list(engine = "cpp")))
  fit.r <- do.call(grip.optimize.kernel.gram.gkk.layout, c(args, list(engine = "R")))

  expect_equal(fit.cpp$coords, fit.r$coords, tolerance = 1e-10)
  expect_equal(fit.cpp$trace$energy, fit.r$trace$energy, tolerance = 1e-10)
  expect_equal(fit.cpp$trace$gram.rel.rmse, fit.r$trace$gram.rel.rmse, tolerance = 1e-10)
  expect_equal(fit.cpp$metadata$engine, "cpp_gradient_descent_armijo")
})

test_that("kernel Gram-gKK improves folded antipodal star Gram error", {
  prepared <- grip.prepare.graph.geodesic.mds(
    edges = rbind(c(1L, 2L), c(2L, 3L)),
    n = 3L,
    edge_weights = c(1, 1)
  )
  X <- cbind(c(0, 1, 2), c(0, 0, 0))
  folded <- rbind(
    c(0, 0),
    c(1, 0),
    c(1, 1)
  )
  star <- graph.riemannian.star.structure(prepared = prepared, X = X, angle.power = 4)
  before <- grip:::grip.kernel.gram.score(folded, star, edge_scale = 1)
  fit <- grip.optimize.kernel.gram.gkk.layout(
    coords = folded,
    prepared = prepared,
    star = star,
    dim = 2L,
    stiffness_method = "uniform",
    density_mix = 1,
    scale_mode = "identity",
    lambda.edge = 1,
    lambda.gram = 1,
    max_iter = 80L,
    initial_step = 0.05,
    engine = "cpp"
  )

  expect_s3_class(fit, "grip_gmds_layout")
  expect_lt(fit$diagnostics$gram.rel.rmse[[1L]], before$gram.rel.rmse[[1L]])
  expect_lt(fit$diagnostics$edge.rel.rmse[[1L]], 0.2)
})
