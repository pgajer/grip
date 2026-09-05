test_that("classical MDS retains cmdscale distances and metadata", {
  p <- prepare.graph.geodesic.mds(edges.cycle(7), n = 7)
  old <- stats::cmdscale(stats::as.dist(p$distance_matrix), k = 3, eig = TRUE)
  fit <- classical.mds(prepared = p, dim = 3, diagnostics = FALSE)
  expect_equal(as.matrix(dist(fit$coords)), as.matrix(dist(old$points)), tolerance = 1e-12)
  expect_identical(fit$method, "classical_mds")
  expect_identical(fit$metadata$objective, "classical_strain")
  expect_equal(fit$metadata$eig, old$eig)
})

test_that("ratio SMACOF wrapper restores input units and matches the backend", {
  skip_if_not_installed("smacof", "2.1-7")
  p <- prepare.graph.geodesic.mds(edges.cycle(7), n = 7, edge_weights = rep(3, 7))
  initial <- classical.mds(prepared = p, diagnostics = FALSE)$coords
  target <- as.double(as.dist(p$distance_matrix))
  rms <- sqrt(mean(target^2))
  x <- initial / rms
  d <- as.double(dist(x))
  x <- x * sum(d * target / rms) / sum(d^2)
  backend <- smacof::mds(p$distance_matrix/rms, type = "ratio", ndim = 2,
                         init = x, eps = 1e-10, itmax = 2000)
  bd <- as.double(dist(backend$conf))
  multiplier <- sum(bd * target) / sum(bd^2)
  fit <- metric.mds(prepared = p, init = initial, eps = 1e-10, max_iter = 2000)
  observed <- as.double(dist(fit$coords))
  expect_equal(observed, multiplier * bd, tolerance = 1e-9)
  expect_equal(fit$metadata$raw_stress, sum((observed-target)^2), tolerance = 1e-10)
  expect_equal(fit$metadata$target_normalized_rmse,
               sqrt(sum((observed-target)^2)/sum(target^2)), tolerance = 1e-12)
  expect_equal(fit$metadata$stress1_identity,
               sqrt(sum((observed-target)^2)/sum(observed^2)), tolerance = 1e-12)
  a <- sum(observed*target)/sum(target^2)
  expect_equal(fit$metadata$stress1_profiled,
               sqrt(sum((observed-a*target)^2)/sum(observed^2)), tolerance = 1e-12)
  expect_equal(fit$metadata$stress1_profiled, fit$metadata$target_normalized_rmse,
               tolerance = 1e-10)
  expect_lt(fit$metadata$raw_stress, sum((as.double(dist(initial))-target)^2))
  expect_identical(fit$metadata$engine, "smacof")
  expect_identical(fit$metadata$objective, "raw_distance_stress")
  expect_true(is.data.frame(fit$diagnostics))
})

test_that("exact Euclidean targets are recovered and changes of units are equivariant", {
  skip_if_not_installed("smacof", "2.1-7")
  x <- rbind(c(0,0,0), c(1,0,0), c(0,2,0), c(0,0,3), c(1,1,1))
  edges <- t(combn(nrow(x), 2))
  lengths <- sqrt(rowSums((x[edges[,1],]-x[edges[,2],])^2))
  fit <- metric.mds(edges = edges, n = nrow(x), edge_weights = lengths,
                    dim = 3, diagnostics = FALSE, eps = 1e-12)
  expect_equal(as.matrix(dist(fit$coords)), as.matrix(dist(x)), tolerance = 1e-8)
  expect_lt(fit$metadata$target_normalized_rmse, 1e-10)
  expect_identical(fit$prepared$pair_mode, "distance_matrix_only")
  large <- metric.mds(edges = edges, n = nrow(x), edge_weights = 1e5*lengths,
                      dim = 3, diagnostics = FALSE, eps = 1e-12)
  expect_equal(as.matrix(dist(large$coords))/1e5, as.matrix(dist(fit$coords)),
               tolerance = 1e-8)
})

test_that("multiple starts are reproducible, preserve the RNG, and select achieved stress", {
  skip_if_not_installed("smacof", "2.1-7")
  p <- prepare.graph.geodesic.mds(edges.cycle(8), n = 8)
  set.seed(817)
  before <- .Random.seed
  a <- suppressWarnings(metric.mds(prepared = p, init = "random", n_init = 3,
                                   seed = 19, diagnostics = FALSE))
  expect_identical(.Random.seed, before)
  b <- suppressWarnings(metric.mds(prepared = p, init = "random", n_init = 3,
                                   seed = 19, diagnostics = FALSE))
  expect_equal(a$coords, b$coords)
  expect_equal(a$metadata$starts, b$metadata$starts)
  expect_equal(nrow(a$metadata$starts), 3)
  expect_equal(a$metadata$raw_stress, min(a$metadata$starts$raw_stress))
  expect_true(all(a$metadata$starts$raw_stress <=
                   a$metadata$starts$initial_raw_stress + 1e-8))
})

test_that("iteration limits remain visible and do not claim convergence", {
  skip_if_not_installed("smacof", "2.1-7")
  expect_warning(fit <- metric.mds(edges = edges.cycle(8), n = 8, init = "random",
                                   max_iter = 1, diagnostics = FALSE), "iteration_limit")
  expect_false(fit$metadata$converged)
  expect_identical(fit$metadata$termination, "iteration_limit")
  expect_equal(fit$metadata$starts$iterations, 1)
})

test_that("MDS initializers dispatch explicitly and defaults remain classical", {
  skip_if_not_installed("smacof", "2.1-7")
  p <- prepare.graph.geodesic.mds(edges.cycle(6), n = 6)
  c <- classical.mds(prepared = p, diagnostics = FALSE)
  m <- metric.mds(prepared = p, diagnostics = FALSE)
  common <- list(prepared = p, max_iter = 0, diagnostics = FALSE,
                  density_mix_schedule = 1, scale_mode = "identity")
  default <- do.call(edge.kk, common)
  explicit <- do.call(edge.kk, c(common, list(init = "classical_mds")))
  stress <- do.call(edge.kk, c(common, list(init = "metric_mds")))
  expect_equal(as.matrix(dist(default$coords)), as.matrix(dist(c$coords)))
  expect_equal(default$coords, explicit$coords)
  expect_equal(as.matrix(dist(stress$coords)), as.matrix(dist(m$coords)))
  cold <- edge.kk(edges = edges.cycle(6), n = 6, init = "metric_mds",
                  max_iter = 0, density_mix_schedule = 1, diagnostics = FALSE)
  expect_equal(as.matrix(dist(cold$coords)), as.matrix(dist(m$coords)))
})

test_that("stress MDS rejects unsuitable inputs instead of silently changing the problem", {
  skip_if_not_installed("smacof", "2.1-7")
  p <- prepare.graph.geodesic.mds(edges.cycle(5), n = 5)
  expect_error(metric.mds(prepared = p, add = TRUE), "unused argument")
  expect_error(metric.mds(prepared = p, n_init = 1.5), "positive integer")
  expect_error(metric.mds(prepared = p, max_iter = 0), "positive integer")
  expect_error(metric.mds(prepared = p, scale_mode = "user"), "arg")
  expect_error(metric.mds(prepared = p, dim = 5), "less than")
  expect_error(metric.mds(prepared = p, init = matrix(0,5,2)), "collapsed")
  expect_error(metric.mds(prepared = p, init = matrix(0,4,2)), "n by dim")
  expect_error(metric.mds(prepared = prepare.edge.kk(edges.cycle(5), n = 5)), "all-pairs")
  expect_error(metric.mds(edges = rbind(c(1,2),c(3,4)), n = 4,
                          diagnostics = FALSE), "connected graph")
  for (bad in c(NA_real_, Inf, -1)) {
    q <- p
    q$distance_matrix[1,2] <- q$distance_matrix[2,1] <- bad
    expect_error(metric.mds(prepared = q, diagnostics = FALSE), "finite, symmetric")
  }
  q <- p
  q$distance_matrix[,] <- 0
  expect_error(metric.mds(prepared = q, diagnostics = FALSE), "positive")
})

test_that("zero target distances and coincident starting points are handled explicitly", {
  skip_if_not_installed("smacof", "2.1-7")
  p <- prepare.graph.geodesic.mds(edges.cycle(5), n = 5)
  x <- rbind(c(0,0),c(0,0),c(1,0),c(0,1),c(1,1))
  p$distance_matrix <- as.matrix(dist(x))
  fit <- metric.mds(prepared = p, init = x, diagnostics = FALSE)
  expect_equal(as.matrix(dist(fit$coords)), as.matrix(dist(x)), tolerance = 1e-8)
})

test_that("backend failures are recorded and never become silent classical fallbacks", {
  skip_if_not_installed("smacof", "2.1-7")
  original <- smacof::mds
  calls <- 0L
  testthat::local_mocked_bindings(mds = function(...) {
    calls <<- calls + 1L
    if (calls == 1L) stop("controlled backend failure")
    original(...)
  }, .package = "smacof")
  fit <- suppressWarnings(metric.mds(edges = edges.cycle(6), n = 6,
                                     n_init = 2, diagnostics = FALSE))
  expect_identical(fit$metadata$starts$termination[1], "backend_error")
  expect_match(fit$metadata$starts$error[1], "controlled backend failure")
  expect_equal(fit$metadata$selected_start, 2)
  calls <- 0L
  expect_error(metric.mds(edges = edges.cycle(6), n = 6, diagnostics = FALSE),
                "All SMACOF starts failed")
})

test_that("a missing optional backend gives an actionable error", {
  testthat::local_mocked_bindings(grip.mds.has.smacof = function() FALSE)
  expect_error(metric.mds(edges = edges.cycle(6), n = 6), "optional 'smacof'")
})

test_that("near-tie routes and strict MDS targets remain separate", {
  skip_if_not_installed("smacof", "2.1-7")
  e <- rbind(c(1L, 2L), c(1L, 3L), c(2L, 3L))
  w <- c(1, 2 + 1e-8, 1)
  p <- prepare.graph.geodesic.mds(e, n = 3L, edge_weights = w)
  expect_identical(p$distance_matrix, t(p$distance_matrix))
  expect_equal(unname(p$distance_matrix), matrix(c(0,1,2,1,0,1,2,1,0), 3), tolerance = 1e-14)
  # The established near-tie convention chooses the direct 1--3 edge.
  i <- which(p$pair_matrix[,1] == 1L & p$pair_matrix[,2] == 3L)
  expect_equal(p$pair_graph_distance[i], 2 + 1e-8, tolerance = 1e-14)
  expect_equal(unname(p$path_edges[[i]]), matrix(c(1L,3L), ncol=2L))
  m <- metric.mds(prepared = p, dim = 2L, diagnostics = FALSE)
  strict <- metric.mds(edges = e, n = 3L, edge_weights = w,
                       dim = 2L, diagnostics = FALSE)
  expect_equal(as.matrix(dist(m$coords)), as.matrix(dist(strict$coords)), tolerance=1e-10)
  expect_lt(m$metadata$raw_stress, 1e-12)
  expect_identical(p$pair_graph_distance[i], 2 + 1e-8)
  q <- p
  q$distance_matrix[1,2] <- 1.1
  expect_error(metric.mds(prepared=q, diagnostics=FALSE), "finite, symmetric")
})
