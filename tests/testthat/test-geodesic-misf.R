test_that("MISF-GMDS preparation carries the coarse level induced graph", {
  edges <- edges.mesh(4, 4)

  prepared <- grip.prepare.misf.geodesic.mds(
    edges = edges,
    n = 16L,
    tie_mode = "average",
    num_init = 4L,
    num_nbrs = 6L,
    top_level_mode = "skip",
    seed = 9L
  )

  expect_s3_class(prepared, "grip_misf_gmds_prepared")
  expect_s3_class(prepared, "grip_gmds_prepared")
  expect_identical(prepared$multiscale_mode, "misf")
  expect_identical(prepared$top_level_vertices, prepared$misf$levels[[length(prepared$misf$levels)]])
  expect_identical(prepared$active_levels, prepared$misf$levels)
  expect_identical(prepared$insertion_order, prepared$misf$mish_order)
  expect_null(prepared$top_level_fit)

  coarse <- prepared$top_level_graph
  expect_equal(coarse$n, length(prepared$top_level_vertices))
  expect_equal(sort(unique(as.vector(coarse$global_edge_matrix))), sort(prepared$top_level_vertices))

  sub.dist <- prepared$distance_matrix[prepared$top_level_vertices, prepared$top_level_vertices, drop = FALSE]
  expect_equal(coarse$distance_matrix, sub.dist)
  expect_equal(
    coarse$edge_weights,
    as.double(sub.dist[cbind(coarse$edges[, 1L], coarse$edges[, 2L])])
  )
  expect_equal(prepared$top_level_prepared$distance_matrix, sub.dist)
})

test_that("top-level MISF-GMDS solve improves over its random restart family", {
  edges <- edges.mesh(4, 4)

  prepared <- grip.prepare.misf.geodesic.mds(
    edges = edges,
    n = 16L,
    tie_mode = "average",
    num_init = 4L,
    num_nbrs = 6L,
    dim = 2L,
    top_level_mode = "solve",
    top_level_restarts = 3L,
    top_level_max_iter = 4L,
    top_level_engine = "cpp",
    seed = 11L
  )

  fit <- prepared$top_level_fit
  expect_true(is.list(fit))
  expect_equal(nrow(fit$coords), length(prepared$top_level_vertices))
  expect_equal(ncol(fit$coords), 2L)
  expect_true(all(is.finite(fit$coords)))
  expect_true(is.data.frame(fit$restart_summary))
  expect_equal(nrow(fit$restart_summary), 3L)
  expect_true(any(fit$restart_summary$improved))
  expect_lte(
    fit$score$gmds.energy[[1L]],
    min(fit$restart_summary$initial.energy) + 1e-10
  )
  expect_equal(fit$vertex_ids, prepared$top_level_vertices)
  expect_equal(nrow(fit$coords_full), prepared$n)
  expect_true(all(is.na(fit$coords_full[-prepared$top_level_vertices, , drop = FALSE])))
  expect_true(all(is.finite(fit$coords_full[prepared$top_level_vertices, , drop = FALSE])))
})

test_that("grip.geodesic.misf.induced_level_graph accepts explicit level selection", {
  edges <- edges.mesh(5, 5)
  prepared <- grip.prepare.misf.geodesic.mds(
    edges = edges,
    n = 25L,
    num_init = 5L,
    num_nbrs = 8L,
    top_level_mode = "skip",
    seed = 5L
  )

  coarse.level1 <- grip:::grip.geodesic.misf.induced_level_graph(prepared, level = 1L)
  level.vertices <- prepared$misf$levels[[2L]]

  expect_equal(coarse.level1$level, 1L)
  expect_equal(coarse.level1$vertex_ids, level.vertices)
  expect_equal(coarse.level1$n, length(level.vertices))
  expect_equal(
    coarse.level1$distance_matrix,
    prepared$distance_matrix[level.vertices, level.vertices, drop = FALSE]
  )
})
