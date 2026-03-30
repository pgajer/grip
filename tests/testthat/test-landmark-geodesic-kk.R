build_test_sierpinski_carpet_lgkk <- function(level) {
  side <- 3L^level
  grid <- expand.grid(x = 0:(side - 1L), y = 0:(side - 1L))

  keep_cell <- function(x, y) {
    while (x > 0L || y > 0L) {
      if ((x %% 3L) == 1L && (y %% 3L) == 1L) {
        return(FALSE)
      }
      x <- x %/% 3L
      y <- y %/% 3L
    }
    TRUE
  }

  keep <- mapply(keep_cell, grid$x, grid$y)
  cells <- grid[keep, , drop = FALSE]
  coords <- cbind(
    x = cells$x + 0.5,
    y = (side - 1L - cells$y) + 0.5
  )

  id_map <- matrix(0L, nrow = side, ncol = side)
  for (i in seq_len(nrow(cells))) {
    id_map[cells$x[i] + 1L, cells$y[i] + 1L] <- i
  }

  edges <- list()
  for (i in seq_len(nrow(cells))) {
    x <- cells$x[i]
    y <- cells$y[i]
    if (x + 1L < side) {
      nbr <- id_map[x + 2L, y + 1L]
      if (nbr > 0L) edges[[length(edges) + 1L]] <- c(i, nbr)
    }
    if (y + 1L < side) {
      nbr <- id_map[x + 1L, y + 2L]
      if (nbr > 0L) edges[[length(edges) + 1L]] <- c(i, nbr)
    }
  }

  list(edges = do.call(rbind, edges), coords = coords)
}

test_that("landmark geodesic KK scoring is exact on a weighted path realization", {
  edges <- matrix(c(
    1L, 2L,
    2L, 3L
  ), ncol = 2L, byrow = TRUE)
  coords <- cbind(
    x = c(0, 1, 3),
    y = c(0, 0, 0)
  )
  prepared <- grip.prepare.landmark.geodesic.kk(
    edges = edges,
    n = 3L,
    edge_weights = c(1, 2),
    local_nbrs = 1L,
    landmark_count = 1L
  )
  score <- grip.score.landmark.geodesic.kk(coords, prepared = prepared)

  expect_equal(score$n.pairs[[1L]], 3L)
  expect_equal(score$scale.L0[[1L]], 1, tolerance = 1e-10)
  expect_lt(score$lgkk.energy[[1L]], 1e-10)
  expect_lt(score$lgkk.weighted.rmse[[1L]], 1e-10)
  expect_lt(score$lgkk.weighted.rel.rmse[[1L]], 1e-10)
})

test_that("prepared landmark geodesic KK caches are deterministic and reusable", {
  edges <- edges.mesh(3, 3)
  prepared1 <- grip.prepare.landmark.geodesic.kk(
    edges = edges,
    n = 9L,
    local_nbrs = 2L,
    landmark_count = 2L
  )
  prepared2 <- grip.prepare.landmark.geodesic.kk(
    edges = edges,
    n = 9L,
    local_nbrs = 2L,
    landmark_count = 2L
  )

  expect_equal(prepared1$pair_matrix, prepared2$pair_matrix)
  expect_equal(prepared1$pair_graph_distance, prepared2$pair_graph_distance)
  expect_equal(prepared1$path_vertices, prepared2$path_vertices)

  set.seed(22)
  coords <- matrix(rnorm(18), ncol = 2L)
  direct <- grip.score.landmark.geodesic.kk(
    coords = coords,
    edges = edges,
    n = 9L,
    local_nbrs = 2L,
    landmark_count = 2L
  )
  cached <- grip.score.landmark.geodesic.kk(coords = coords, prepared = prepared1)

  expect_equal(direct$lgkk.energy[[1L]], cached$lgkk.energy[[1L]], tolerance = 1e-10)
  expect_equal(direct$scale.L0[[1L]], cached$scale.L0[[1L]], tolerance = 1e-10)
})

test_that("landmark geodesic KK prefers the canonical carpet over a perturbed carpet", {
  built <- build_test_sierpinski_carpet_lgkk(2L)
  prepared <- grip.prepare.landmark.geodesic.kk(
    edges = built$edges,
    n = nrow(built$coords),
    local_nbrs = 6L,
    landmark_count = 4L
  )

  perturbed <- built$coords
  set.seed(11)
  perturbed <- perturbed + matrix(rnorm(length(perturbed), sd = 0.15), ncol = 2L)

  canonical_score <- grip.score.landmark.geodesic.kk(built$coords, prepared = prepared)
  perturbed_score <- grip.score.landmark.geodesic.kk(perturbed, prepared = prepared)

  expect_lt(canonical_score$lgkk.energy[[1L]], perturbed_score$lgkk.energy[[1L]])
  expect_lt(canonical_score$lgkk.weighted.rmse[[1L]], perturbed_score$lgkk.weighted.rmse[[1L]])
})

test_that("landmark geodesic KK optimizer decreases the prototype energy", {
  built <- build_test_sierpinski_carpet_lgkk(2L)
  prepared <- grip.prepare.landmark.geodesic.kk(
    edges = built$edges,
    n = nrow(built$coords),
    local_nbrs = 6L,
    landmark_count = 4L
  )

  perturbed <- built$coords
  set.seed(19)
  perturbed <- perturbed + matrix(rnorm(length(perturbed), sd = 0.2), ncol = 2L)

  before <- grip.score.landmark.geodesic.kk(perturbed, prepared = prepared)
  opt <- grip.optimize.landmark.geodesic.kk(
    coords = perturbed,
    prepared = prepared,
    max_iter = 8L,
    return_trace = TRUE
  )
  after <- grip.score.landmark.geodesic.kk(opt$coords, prepared = prepared)

  expect_lt(after$lgkk.energy[[1L]], before$lgkk.energy[[1L]])
  expect_lt(after$lgkk.weighted.rel.rmse[[1L]], before$lgkk.weighted.rel.rmse[[1L]])
  expect_true(nrow(opt$trace) >= 2L)
  expect_true(length(opt$frames) >= 2L)
})
