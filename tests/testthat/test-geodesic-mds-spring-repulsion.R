finite_difference_gradient <- function(fn, x, eps = 1e-6) {
  grad <- numeric(length(x))
  for (i in seq_along(x)) {
    plus <- x
    minus <- x
    plus[[i]] <- plus[[i]] + eps
    minus[[i]] <- minus[[i]] - eps
    grad[[i]] <- (fn(plus) - fn(minus)) / (2 * eps)
  }
  grad
}

test_that("prepared geodesic objects retain graph edge targets and repulsion caches", {
  prepared <- grip.prepare.geodesic.kk(edges = edges.path(4L), n = 4L)
  expect_equal(prepared$edge_targets, rep(1, 3))

  cached <- grip:::grip.geodesic.mds.ensure.graph.term.cache(
    prepared = prepared,
    repulsion_weight = 0.2,
    repulsion_quantile = 0,
    repulsion_scale = 0.5,
    repulsion_cap_quantile = 1,
    repulsion_hop_min = 2L
  )

  expect_equal(cached$graph_edge_target, rep(1, 3))
  expect_true(is.matrix(cached$repulsion_pair_matrix))
  expect_true(nrow(cached$repulsion_pair_matrix) >= 1L)
  expect_true(all(cached$repulsion_target > 0))
})

test_that("edge spring energy gradient matches finite differences", {
  prepared <- grip.prepare.geodesic.kk(edges = edges.path(4L), n = 4L)
  coords <- rbind(
    c(0.0, 0.0),
    c(0.8, 0.2),
    c(1.4, -0.1),
    c(2.9, 0.3)
  )
  weight <- 0.7

  analytic <- grip:::grip.geodesic.mds.edge.spring.stats(
    coords = coords,
    prepared = prepared,
    edge_spring_weight = weight
  )

  numeric <- finite_difference_gradient(function(vec) {
    pts <- matrix(vec, ncol = 2L, byrow = FALSE)
    grip:::grip.geodesic.mds.edge.spring.stats(
      coords = pts,
      prepared = prepared,
      edge_spring_weight = weight
    )$energy
  }, as.double(coords))

  expect_equal(as.double(analytic$gradient), numeric, tolerance = 1e-5)
})

test_that("graph-aware repulsion energy gradient matches finite differences", {
  prepared <- grip.prepare.geodesic.kk(edges = edges.path(4L), n = 4L)
  prepared <- grip:::grip.geodesic.mds.ensure.graph.term.cache(
    prepared = prepared,
    repulsion_weight = 0.9,
    repulsion_quantile = 0,
    repulsion_scale = 0.8,
    repulsion_cap_quantile = 1,
    repulsion_hop_min = 2L
  )
  coords <- rbind(
    c(0.0, 0.0),
    c(0.3, 0.0),
    c(0.6, 0.0),
    c(0.8, 0.1)
  )
  weight <- 0.9

  analytic <- grip:::grip.geodesic.mds.repulsion.stats(
    coords = coords,
    prepared = prepared,
    repulsion_weight = weight
  )

  numeric <- finite_difference_gradient(function(vec) {
    pts <- matrix(vec, ncol = 2L, byrow = FALSE)
    grip:::grip.geodesic.mds.repulsion.stats(
      coords = pts,
      prepared = prepared,
      repulsion_weight = weight
    )$energy
  }, as.double(coords))

  expect_true(analytic$active_pair_count >= 1L)
  expect_equal(as.double(analytic$gradient), numeric, tolerance = 1e-5)
})

test_that("geodesic MDS scoring reports spring and repulsion contributions", {
  prepared <- grip.prepare.geodesic.kk(edges = edges.path(4L), n = 4L)
  coords <- rbind(
    c(0.0, 0.0),
    c(0.5, 0.0),
    c(0.9, 0.0),
    c(1.2, 0.0)
  )

  score <- grip.score.geodesic.mds(
    coords = coords,
    prepared = prepared,
    edge_spring_weight = 0.25,
    repulsion_weight = 0.5,
    repulsion_quantile = 0,
    repulsion_scale = 0.7,
    repulsion_cap_quantile = 1,
    repulsion_hop_min = 2L
  )

  expect_true(score$edge.spring.energy[[1L]] > 0)
  expect_true(score$repulsion.energy[[1L]] > 0)
  expect_true(score$repulsion.active.pair.count[[1L]] >= 1L)
  expect_equal(
    score$gmds.energy[[1L]],
    score$gmds.base.energy[[1L]] +
      score$anchor.energy[[1L]] +
      score$edge.spring.energy[[1L]] +
      score$repulsion.energy[[1L]] +
      score$smooth.energy[[1L]],
    tolerance = 1e-8
  )
})

test_that("optimizer falls back to the R engine and reduces repulsion energy", {
  prepared <- grip.prepare.geodesic.kk(edges = edges.path(4L), n = 4L)
  coords <- rbind(
    c(0.0, 0.0),
    c(0.2, 0.0),
    c(0.4, 0.0),
    c(0.6, 0.0)
  )
  args <- list(
    prepared = prepared,
    edge_spring_weight = 0.1,
    repulsion_weight = 0.4,
    repulsion_quantile = 0,
    repulsion_scale = 0.6,
    repulsion_cap_quantile = 1,
    repulsion_hop_min = 2L
  )

  before <- do.call(grip.score.geodesic.mds, c(list(coords = coords), args))
  opt <- NULL
  expect_warning(
    opt <- do.call(
      grip.optimize.geodesic.mds,
      c(
        list(coords = coords, engine = "cpp", max_iter = 6L, return_trace = TRUE),
        args
      )
    ),
    "falling back to the R engine"
  )
  after <- do.call(grip.score.geodesic.mds, c(list(coords = opt$coords), args))

  expect_true(all(c("edge_spring_energy", "repulsion_energy") %in% names(opt$trace)))
  expect_lt(after$gmds.energy[[1L]], before$gmds.energy[[1L]])
  expect_lt(after$repulsion.energy[[1L]], before$repulsion.energy[[1L]])
  expect_gt(
    sqrt(sum((opt$coords[1L, ] - opt$coords[4L, ])^2)),
    sqrt(sum((coords[1L, ] - coords[4L, ])^2))
  )
})
