path_raw_stress <- function(score_df) {
  as.double(score_df$gmds.raw_stress[[1L]])
}

edge_lengths_from_coords <- function(edges, coords) {
  diffs <- coords[edges[, 1L], , drop = FALSE] - coords[edges[, 2L], , drop = FALSE]
  sqrt(rowSums(diffs^2))
}

test_that("path graphs admit exact zero-stress realizations", {
  prepared <- prepare.graph.geodesic.mds(
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
  prepared <- prepare.graph.geodesic.mds(
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

test_that("selected tree examples admit exact zero-stress realizations", {
  star.edges <- matrix(
    c(
      1L, 2L,
      1L, 3L,
      1L, 4L,
      1L, 5L
    ),
    ncol = 2L,
    byrow = TRUE
  )
  star.coords <- rbind(
    c(0.0, 0.0),
    c(1.0, 0.0),
    c(-2.0, 0.0),
    c(0.0, 1.5),
    c(0.0, -0.75)
  )
  star.prepared <- prepare.graph.geodesic.mds(
    edges = star.edges,
    n = 5L,
    edge_weights = edge_lengths_from_coords(star.edges, star.coords),
    tie_mode = "single"
  )
  star.score <- grip.score.geodesic.mds(coords = star.coords, prepared = star.prepared)
  expect_lt(path_raw_stress(star.score), 1e-12)

  tree.edges <- edges.kary.tree(k = 2L, depth = 2L)
  tree.coords <- rbind(
    c(0.0, 0.0),
    c(-1.0, 1.0),
    c(1.25, 0.8),
    c(-2.0, 2.0),
    c(-0.1, 2.3),
    c(0.4, 2.1),
    c(2.2, 1.9)
  )
  tree.prepared <- prepare.graph.geodesic.mds(
    edges = tree.edges,
    n = nrow(tree.coords),
    edge_weights = edge_lengths_from_coords(tree.edges, tree.coords),
    tie_mode = "single"
  )
  tree.score <- grip.score.geodesic.mds(coords = tree.coords, prepared = tree.prepared)
  expect_lt(path_raw_stress(tree.score), 1e-12)
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
  prepared_single <- prepare.graph.geodesic.mds(
    edges = square_edges,
    n = 4L,
    edge_weights = rep(1, 4L),
    tie_mode = "single"
  )
  prepared_average <- prepare.graph.geodesic.mds(
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

  prepared_a <- prepare.graph.geodesic.mds(
    edges = square_edges,
    n = 4L,
    edge_weights = rep(1, 4L),
    tie_mode = "average"
  )
  prepared_b <- prepare.graph.geodesic.mds(
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

test_that("diamond ties document single-mode label dependence and average-mode invariance", {
  diamond_edges <- matrix(
    c(
      1L, 2L,
      1L, 3L,
      2L, 4L,
      3L, 4L
    ),
    ncol = 2L,
    byrow = TRUE
  )
  perm <- c(1L, 3L, 2L, 4L)
  diamond_edges_perm <- matrix(perm[diamond_edges], ncol = 2L)
  coords <- rbind(
    c(0.0, 0.0),
    c(1.0, 0.0),
    c(0.0, 3.0),
    c(1.0, 1.0)
  )

  prepared_single <- prepare.graph.geodesic.mds(
    edges = diamond_edges,
    n = 4L,
    edge_weights = rep(1, 4L),
    tie_mode = "single"
  )
  prepared_single_perm <- prepare.graph.geodesic.mds(
    edges = diamond_edges_perm,
    n = 4L,
    edge_weights = rep(1, 4L),
    tie_mode = "single"
  )
  prepared_average <- prepare.graph.geodesic.mds(
    edges = diamond_edges,
    n = 4L,
    edge_weights = rep(1, 4L),
    tie_mode = "average"
  )
  prepared_average_perm <- prepare.graph.geodesic.mds(
    edges = diamond_edges_perm,
    n = 4L,
    edge_weights = rep(1, 4L),
    tie_mode = "average"
  )

  pair_idx <- which(prepared_average$pair_matrix[, 1L] == 1L & prepared_average$pair_matrix[, 2L] == 4L)
  expect_length(pair_idx, 1L)
  expect_equal(exp(prepared_average$pair_path_count_log[[pair_idx]]), 2, tolerance = 1e-12)
  expect_equal(sort(grip:::grip.path.edge.coefficients(prepared_average, pair_idx, nrow(prepared_average$path_edges[[pair_idx]]))), rep(0.5, 4L))

  score_single_a <- grip.score.geodesic.mds(coords = coords, prepared = prepared_single)
  score_single_b <- grip.score.geodesic.mds(coords = coords[perm, , drop = FALSE], prepared = prepared_single_perm)
  score_average_a <- grip.score.geodesic.mds(coords = coords, prepared = prepared_average)
  score_average_b <- grip.score.geodesic.mds(coords = coords[perm, , drop = FALSE], prepared = prepared_average_perm)

  expect_gt(abs(path_raw_stress(score_single_a) - path_raw_stress(score_single_b)), 1)
  expect_equal(path_raw_stress(score_average_a), path_raw_stress(score_average_b), tolerance = 1e-12)
})

test_that("3x3 lattice patches restore symmetry under tie_mode average", {
  lattice_edges <- edges.mesh(3L, 3L)
  perm <- c(1L, 4L, 7L, 2L, 5L, 8L, 3L, 6L, 9L)
  lattice_edges_perm <- matrix(perm[lattice_edges], ncol = 2L)
  coords <- cbind(
    c(0, 1, 2,
      0, 1, 2,
      0, 1, 2),
    c(0.0, 0.4, 0.0,
      1.5, 1.8, 1.3,
      2.7, 2.9, 2.5)
  )

  prepared_single <- prepare.graph.geodesic.mds(
    edges = lattice_edges,
    n = 9L,
    edge_weights = rep(1, nrow(lattice_edges)),
    tie_mode = "single"
  )
  prepared_single_perm <- prepare.graph.geodesic.mds(
    edges = lattice_edges_perm,
    n = 9L,
    edge_weights = rep(1, nrow(lattice_edges_perm)),
    tie_mode = "single"
  )
  prepared_average <- prepare.graph.geodesic.mds(
    edges = lattice_edges,
    n = 9L,
    edge_weights = rep(1, nrow(lattice_edges)),
    tie_mode = "average"
  )
  prepared_average_perm <- prepare.graph.geodesic.mds(
    edges = lattice_edges_perm,
    n = 9L,
    edge_weights = rep(1, nrow(lattice_edges_perm)),
    tie_mode = "average"
  )

  pair_idx <- which(prepared_average$pair_matrix[, 1L] == 1L & prepared_average$pair_matrix[, 2L] == 9L)
  expect_length(pair_idx, 1L)
  expect_equal(exp(prepared_average$pair_path_count_log[[pair_idx]]), 6, tolerance = 1e-12)
  expect_equal(sum(prepared_average$path_edge_weights[[pair_idx]]), 4, tolerance = 1e-12)
  expect_equal(nrow(prepared_average$path_edges[[pair_idx]]), 12L)

  score_single_a <- grip.score.geodesic.mds(coords = coords, prepared = prepared_single)
  score_single_b <- grip.score.geodesic.mds(coords = coords[perm, , drop = FALSE], prepared = prepared_single_perm)
  score_average_a <- grip.score.geodesic.mds(coords = coords, prepared = prepared_average)
  score_average_b <- grip.score.geodesic.mds(coords = coords[perm, , drop = FALSE], prepared = prepared_average_perm)

  expect_gt(abs(path_raw_stress(score_single_a) - path_raw_stress(score_single_b)), 0.1)
  expect_equal(path_raw_stress(score_average_a), path_raw_stress(score_average_b), tolerance = 1e-12)
})
