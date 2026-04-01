path_raw_stress <- function(score_df) {
  as.double(score_df$gmds.raw_stress[[1L]])
}

test_that("path graphs admit exact zero-stress realizations", {
  prepared <- grip.prepare.graph.geodesic.mds(
    edges = edges.path(4L),
    n = 4L,
    edge_weights = c(1, 1.5, 0.7),
    tie_mode = "single"
  )

  coords_exact <- cbind(c(0, 1, 2.5, 3.2), 0)
  score_exact <- grip.score.geodesic.mds(coords = coords_exact, prepared = prepared)
  expect_lt(path_raw_stress(score_exact), 1e-12)

  coords_perturbed <- coords_exact + rbind(
    c(0.0, 0.0),
    c(0.1, 0.2),
    c(-0.1, 0.1),
    c(0.2, -0.1)
  )
  fit <- grip.optimize.geodesic.mds(
    coords = coords_perturbed,
    prepared = prepared,
    engine = "cpp",
    max_iter = 25L,
    return_trace = TRUE,
    n_threads = 1L
  )
  expect_lt(path_raw_stress(fit$score), 1e-6)
})

test_that("shared-edge triangles show nonunique exact GMDS minima", {
  s3 <- sqrt(3) / 2
  prepared <- grip.prepare.graph.geodesic.mds(
    edges = matrix(
      c(
        1L, 2L,
        1L, 3L,
        2L, 3L,
        2L, 4L,
        3L, 4L
      ),
      ncol = 2L,
      byrow = TRUE
    ),
    n = 4L,
    edge_weights = rep(1, 5L),
    tie_mode = "single"
  )

  coords_separated <- rbind(
    c(0.5, s3),
    c(0.0, 0.0),
    c(1.0, 0.0),
    c(0.5, -s3)
  )
  coords_overlapped <- rbind(
    c(0.5, s3),
    c(0.0, 0.0),
    c(1.0, 0.0),
    c(0.5, s3)
  )

  score_separated <- grip.score.geodesic.mds(coords = coords_separated, prepared = prepared)
  score_overlapped <- grip.score.geodesic.mds(coords = coords_overlapped, prepared = prepared)

  expect_lt(path_raw_stress(score_separated), 1e-12)
  expect_lt(path_raw_stress(score_overlapped), 1e-12)
  expect_gt(sum((coords_separated - coords_overlapped)^2), 1)
})

test_that("tie_mode average records tied shortest paths on a square", {
  square_edges <- matrix(
    c(
      1L, 2L,
      2L, 3L,
      3L, 4L,
      4L, 1L
    ),
    ncol = 2L,
    byrow = TRUE
  )
  prepared_single <- grip.prepare.graph.geodesic.mds(
    edges = square_edges,
    n = 4L,
    edge_weights = rep(1, 4L),
    tie_mode = "single"
  )
  prepared_average <- grip.prepare.graph.geodesic.mds(
    edges = square_edges,
    n = 4L,
    edge_weights = rep(1, 4L),
    tie_mode = "average"
  )

  pair_idx <- which(
    prepared_single$pair_matrix[, 1L] == 1L &
      prepared_single$pair_matrix[, 2L] == 3L
  )
  expect_length(pair_idx, 1L)

  single_edges <- prepared_single$path_edges[[pair_idx]]
  avg_edges <- prepared_average$path_edges[[pair_idx]]
  avg_coeffs <- grip:::grip.path.edge.coefficients(prepared_average, pair_idx, nrow(avg_edges))

  expect_equal(nrow(single_edges), 2L)
  expect_equal(nrow(avg_edges), 4L)
  expect_equal(prepared_average$pair_path_count_log[[pair_idx]], log(2), tolerance = 1e-10)
  expect_equal(sort(avg_coeffs), rep(0.5, 4L))
})

test_that("tie_mode average is relabeling-invariant on the square objective", {
  square_edges <- matrix(
    c(
      1L, 2L,
      2L, 3L,
      3L, 4L,
      4L, 1L
    ),
    ncol = 2L,
    byrow = TRUE
  )
  perm <- c(2L, 3L, 4L, 1L)
  square_edges_perm <- matrix(perm[square_edges], ncol = 2L)

  prepared_a <- grip.prepare.graph.geodesic.mds(
    edges = square_edges,
    n = 4L,
    edge_weights = rep(1, 4L),
    tie_mode = "average"
  )
  prepared_b <- grip.prepare.graph.geodesic.mds(
    edges = square_edges_perm,
    n = 4L,
    edge_weights = rep(1, 4L),
    tie_mode = "average"
  )

  square_coords <- rbind(
    c(0, 0),
    c(1, 0),
    c(1, 1),
    c(0, 1)
  )
  score_a <- grip.score.geodesic.mds(coords = square_coords, prepared = prepared_a)
  score_b <- grip.score.geodesic.mds(coords = square_coords[perm, , drop = FALSE], prepared = prepared_b)

  expect_equal(path_raw_stress(score_a), path_raw_stress(score_b), tolerance = 1e-12)
})
