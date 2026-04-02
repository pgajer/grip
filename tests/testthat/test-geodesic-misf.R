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

test_that("compiled MISF insertion recovers a realizable point from anchor distances", {
  anchor.coords <- rbind(
    c(0, 0),
    c(2, 0),
    c(0, 2),
    c(2, 2)
  )
  target <- c(1.25, 0.75)
  anchor.dist <- sqrt(rowSums((anchor.coords - matrix(target, nrow = 4L, ncol = 2L, byrow = TRUE))^2))

  fit <- grip:::grip_geodesic_misf_insert_vertex_cpp(
    anchor_coords = anchor.coords,
    anchor_distance = anchor.dist,
    anchor_weights = rep(1, 4L),
    max_iter = 80L
  )

  expect_equal(as.double(fit$coord), target, tolerance = 1e-6)
  expect_lte(fit$objective, 1e-10)
})

test_that("MISF anchor selection is deterministic and previous-level restricted", {
  bundle <- mesh.surface.graph(
    5, 5,
    surface = "paraboloid",
    amplitude = 0.2,
    connectivity = "orthogonal",
    normalize = "median"
  )
  prepared <- grip.prepare.misf.geodesic.mds(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average",
    num_init = 5L,
    num_nbrs = 8L,
    dim = 3L,
    top_level_mode = "solve",
    top_level_restarts = 2L,
    top_level_max_iter = 2L,
    seed = 4L
  )
  coords <- prepared$top_level_fit$coords_full
  level <- prepared$top_level_level - 1L
  vertex <- grip:::grip.geodesic.misf.level.insert.vertices(prepared, level)[[1L]]

  first1 <- grip:::grip.geodesic.misf.select.anchors(
    prepared = prepared,
    coords = coords,
    vertex = vertex,
    level = level,
    anchor_policy = "prev_level_first",
    anchor_count = 4L
  )
  first2 <- grip:::grip.geodesic.misf.select.anchors(
    prepared = prepared,
    coords = coords,
    vertex = vertex,
    level = level,
    anchor_policy = "prev_level_first",
    anchor_count = 4L
  )
  spread1 <- grip:::grip.geodesic.misf.select.anchors(
    prepared = prepared,
    coords = coords,
    vertex = vertex,
    level = level,
    anchor_policy = "prev_level_spread",
    anchor_count = 4L
  )
  spread2 <- grip:::grip.geodesic.misf.select.anchors(
    prepared = prepared,
    coords = coords,
    vertex = vertex,
    level = level,
    anchor_policy = "prev_level_spread",
    anchor_count = 4L
  )
  prev.level <- prepared$misf$levels[[level + 2L]]

  expect_identical(first1$anchor_ids, first2$anchor_ids)
  expect_identical(spread1$anchor_ids, spread2$anchor_ids)
  expect_true(all(first1$anchor_ids %in% prev.level))
  expect_true(all(spread1$anchor_ids %in% prev.level))
})

test_that("MISF level insertion yields finite layouts on regular and irregular paraboloids", {
  regular <- mesh.surface.graph(
    6, 6,
    surface = "paraboloid",
    amplitude = 0.25,
    connectivity = "orthogonal",
    normalize = "median"
  )
  irregular <- occupied.mesh.surface.graph(
    keep = keep.asymmetric.notches(6, 6, notch_depth = 2, notch_width = 1),
    surface = "paraboloid",
    amplitude = 0.25,
    connectivity = "orthogonal",
    normalize = "median"
  )

  cases <- list(regular = regular, irregular = irregular)
  for (case.name in names(cases)) {
    bundle <- cases[[case.name]]
    prepared <- grip.prepare.misf.geodesic.mds(
      edges = bundle$edges,
      n = bundle$n,
      edge_weights = bundle$edge_weights,
      tie_mode = "average",
      num_init = 6L,
      num_nbrs = 8L,
      dim = 3L,
      top_level_mode = "solve",
      top_level_restarts = 2L,
      top_level_max_iter = 3L,
      seed = 8L
    )
    inserted <- grip:::grip.geodesic.misf.insert.all.levels(
      prepared = prepared,
      anchor_policy = "prev_level_spread",
      max_iter = 40L
    )

    expect_equal(dim(inserted$coords), c(bundle$n, 3L), info = case.name)
    expect_true(all(is.finite(inserted$coords)), info = case.name)
    expect_true(nrow(inserted$level_trace) >= 1L, info = case.name)
  }
})

test_that("MISF sparse level-pair builder respects the active level", {
  bundle <- mesh.surface.graph(
    6, 6,
    surface = "paraboloid",
    amplitude = 0.2,
    connectivity = "orthogonal",
    normalize = "median"
  )
  prepared <- grip.prepare.misf.geodesic.mds(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average",
    num_init = 6L,
    num_nbrs = 8L,
    dim = 3L,
    top_level_mode = "solve",
    top_level_restarts = 2L,
    top_level_max_iter = 3L,
    seed = 3L
  )
  built <- grip:::grip.geodesic.misf.build.level.pairs(
    prepared = prepared,
    level = 1L,
    local_nbrs = 3L,
    landmark_count = 2L,
    pair_mode = "sparse"
  )

  expect_equal(built$active_vertices, prepared$misf$levels[[2L]])
  expect_s3_class(built$active_prepared, "grip_gmds_prepared")
  expect_true(nrow(built$pair_matrix) > 0L)
  expect_true(all(as.vector(built$global_pair_matrix) %in% built$active_vertices))
  expect_equal(built$active_prepared$active_vertex_ids, built$active_vertices)
})

test_that("MISF sparse refinement lowers active-level energy and final polish lowers global energy", {
  bundle <- mesh.surface.graph(
    6, 6,
    surface = "paraboloid",
    amplitude = 0.25,
    connectivity = "orthogonal",
    normalize = "median"
  )
  prepared <- grip.prepare.misf.geodesic.mds(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average",
    num_init = 6L,
    num_nbrs = 8L,
    dim = 3L,
    top_level_mode = "solve",
    top_level_restarts = 2L,
    top_level_max_iter = 3L,
    seed = 12L
  )
  inserted <- grip:::grip.geodesic.misf.insert.all.levels(
    prepared = prepared,
    anchor_policy = "prev_level_spread",
    max_iter = 40L
  )
  before.global <- grip.score.geodesic.mds(inserted$coords, prepared = prepared)
  refined <- grip:::grip.geodesic.misf.refine.level(
    prepared = prepared,
    coords = inserted$coords,
    level = 0L,
    local_nbrs = 4L,
    landmark_count = 2L,
    pair_mode = "sparse",
    anchor_weight = 0.05,
    max_iter = 4L,
    engine = "cpp",
    n_threads = 1L,
    return_trace = TRUE
  )
  polished <- grip:::grip.geodesic.misf.final.polish(
    prepared = prepared,
    coords = refined$coords,
    max_iter = 4L,
    engine = "cpp",
    n_threads = 1L,
    return_trace = TRUE
  )

  expect_lte(refined$after$gmds.energy[[1L]], refined$before$gmds.energy[[1L]] + 1e-8)
  expect_true(length(refined$pinned_vertices) > 0L)
  expect_true(any(refined$anchor_vertex_weight > 0))
  expect_true(all(is.finite(refined$coords)))
  expect_true(all(is.finite(polished$coords)))
  expect_lte(polished$score$gmds.energy[[1L]], before.global$gmds.energy[[1L]] + 1e-8)
})
