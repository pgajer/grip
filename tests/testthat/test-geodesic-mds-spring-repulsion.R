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
  prepared <- prepare.geodesic.kk(edges = edges.path(4L), n = 4L)
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
  expect_equal(cached$flat_repulsion_u, cached$repulsion_pair_matrix[, 1L] - 1L)
  expect_equal(cached$flat_repulsion_v, cached$repulsion_pair_matrix[, 2L] - 1L)
  expect_equal(cached$flat_repulsion_target, cached$repulsion_target)
})

test_that("repulsion pair selection is deterministic and matches the flattened cache", {
  prepared <- prepare.geodesic.kk(edges = edges.path(4L), n = 4L)
  cached <- grip:::grip.geodesic.mds.ensure.graph.term.cache(
    prepared = prepared,
    repulsion_weight = 0.2,
    repulsion_quantile = 0,
    repulsion_scale = 0.5,
    repulsion_cap_quantile = 1,
    repulsion_hop_min = 2L
  )

  expect_equal(
    cached$repulsion_pair_matrix,
    matrix(c(1L, 3L,
             1L, 4L,
             2L, 4L), ncol = 2L, byrow = TRUE)
  )
  expect_equal(cached$repulsion_target, c(1, 1.5, 1))
  expect_equal(cached$flat_repulsion_u, c(0L, 0L, 1L))
  expect_equal(cached$flat_repulsion_v, c(2L, 3L, 3L))
  expect_equal(cached$flat_repulsion_target, c(1, 1.5, 1))
})

test_that("edge spring energy gradient matches finite differences", {
  prepared <- prepare.geodesic.kk(edges = edges.path(4L), n = 4L)
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
  prepared <- prepare.geodesic.kk(edges = edges.path(4L), n = 4L)
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

test_that("compiled repulsion stats match the R implementation exactly on a deterministic cache", {
  prepared <- prepare.geodesic.kk(edges = edges.path(4L), n = 4L)
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

  r.stats <- grip:::grip.geodesic.mds.repulsion.stats(
    coords = coords,
    prepared = prepared,
    repulsion_weight = weight
  )
  cpp.stats <- grip:::grip_geodesic_mds_flat_repulsion_stats_cpp(
    repulsion_u = prepared$flat_repulsion_u,
    repulsion_v = prepared$flat_repulsion_v,
    repulsion_target = prepared$flat_repulsion_target,
    coords = coords,
    edge_length_epsilon = 1e-8,
    repulsion_weight = weight
  )

  expect_equal(cpp.stats$pair_count, r.stats$pair_count)
  expect_equal(cpp.stats$active_pair_count, r.stats$active_pair_count)
  expect_equal(cpp.stats$energy, r.stats$energy, tolerance = 1e-10)
  expect_equal(as.double(cpp.stats$gradient), as.double(r.stats$gradient), tolerance = 1e-10)
})

test_that("geodesic MDS scoring reports spring and repulsion contributions", {
  prepared <- prepare.geodesic.kk(edges = edges.path(4L), n = 4L)
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

test_that("compiled optimizer matches the R engine for spring-repulsion GMDS", {
  prepared <- prepare.geodesic.kk(edges = edges.path(4L), n = 4L)
  coords <- rbind(
    c(0.0, 0.0),
    c(0.2, 0.0),
    c(0.4, 0.0),
    c(0.6, 0.0)
  )
  score.args <- list(
    prepared = prepared,
    edge_spring_weight = 0.1,
    repulsion_weight = 0.4,
    repulsion_quantile = 0,
    repulsion_scale = 0.6,
    repulsion_cap_quantile = 1,
    repulsion_hop_min = 2L
  )
  opt.args <- c(
    score.args,
    list(
    initial_step = 0.5,
    max_iter = 6L,
    return_trace = TRUE,
    n_threads = 1L
    )
  )

  before <- do.call(grip.score.geodesic.mds, c(list(coords = coords), score.args))
  opt.r <- do.call(
    grip.optimize.geodesic.mds,
    c(list(coords = coords, engine = "r"), opt.args)
  )
  opt.cpp <- NULL
  expect_no_warning(
    opt.cpp <- do.call(
      grip.optimize.geodesic.mds,
      c(list(coords = coords, engine = "cpp"), opt.args)
    )
  )
  after <- do.call(grip.score.geodesic.mds, c(list(coords = opt.cpp$coords), score.args))

  expect_true(all(c("edge_spring_energy", "repulsion_energy") %in% names(opt.cpp$trace)))
  expect_true(all(c("edge_spring_weight", "repulsion_weight", "repulsion_pair_count", "repulsion_active_pair_count") %in% names(opt.cpp$trace)))
  expect_equal(opt.cpp$coords, opt.r$coords, tolerance = 1e-6)
  expect_equal(opt.cpp$trace$energy, opt.r$trace$energy, tolerance = 1e-6)
  expect_equal(opt.cpp$trace$edge_spring_energy, opt.r$trace$edge_spring_energy, tolerance = 1e-6)
  expect_equal(opt.cpp$trace$repulsion_energy, opt.r$trace$repulsion_energy, tolerance = 1e-6)
  expect_equal(opt.cpp$trace$repulsion_pair_count, opt.r$trace$repulsion_pair_count)
  expect_equal(opt.cpp$trace$repulsion_active_pair_count, opt.r$trace$repulsion_active_pair_count)
  expect_equal(opt.cpp$score$gmds.energy[[1L]], opt.r$score$gmds.energy[[1L]], tolerance = 1e-6)
  expect_lt(after$gmds.energy[[1L]], before$gmds.energy[[1L]])
  expect_lt(after$repulsion.energy[[1L]], before$repulsion.energy[[1L]])
  expect_gt(
    sqrt(sum((opt.cpp$coords[1L, ] - opt.cpp$coords[4L, ])^2)),
    sqrt(sum((coords[1L, ] - coords[4L, ])^2))
  )
  expect_equal(opt.cpp$n_threads_used, 1L)
})
