test_that("geometric coarse seed builds finite top-level placements", {
  edges <- edges.mesh(6, 6)
  prepared2 <- grip.prepare.misf.geodesic.mds(
    edges = edges,
    n = 36,
    num_init = 6L,
    dim = 2L,
    top_level_mode = "skip",
    seed = 17L
  )
  init2 <- grip:::grip.geodesic.misf.build.geometric.seed.coords(
    distance_matrix = prepared2$top_level_graph$distance_matrix,
    dim = 2L,
    vertex_ids = prepared2$top_level_vertices,
    insertion_order = prepared2$insertion_order[prepared2$insertion_order %in% prepared2$top_level_vertices],
    anchor_count = prepared2$insertion_anchor_count,
    anchor_weight_mode = prepared2$insertion_anchor_weight_mode
  )

  expect_true(all(is.finite(init2$coords)))
  expect_equal(nrow(init2$coords), length(prepared2$top_level_vertices))
  expect_equal(ncol(init2$coords), 2L)
  expect_equal(length(init2$seed_vertices), min(3L, length(prepared2$top_level_vertices)))

  prepared3 <- grip.prepare.misf.geodesic.mds(
    edges = edges,
    n = 36,
    num_init = 6L,
    dim = 3L,
    top_level_mode = "skip",
    seed = 17L
  )
  init3 <- grip:::grip.geodesic.misf.build.geometric.seed.coords(
    distance_matrix = prepared3$top_level_graph$distance_matrix,
    dim = 3L,
    vertex_ids = prepared3$top_level_vertices,
    insertion_order = prepared3$insertion_order[prepared3$insertion_order %in% prepared3$top_level_vertices],
    anchor_count = prepared3$insertion_anchor_count,
    anchor_weight_mode = prepared3$insertion_anchor_weight_mode
  )

  expect_true(all(is.finite(init3$coords)))
  expect_equal(nrow(init3$coords), length(prepared3$top_level_vertices))
  expect_equal(ncol(init3$coords), 3L)
  expect_equal(length(init3$seed_vertices), min(4L, length(prepared3$top_level_vertices)))
})

test_that("geometric coarse seed handles top levels smaller than d plus one", {
  edges <- cbind(1:2, 2:3)
  prepared <- grip.prepare.misf.geodesic.mds(
    edges = edges,
    n = 3,
    num_init = 2L,
    dim = 3L,
    top_level_mode = "skip",
    seed = 19L
  )
  init <- grip:::grip.geodesic.misf.build.geometric.seed.coords(
    distance_matrix = prepared$top_level_graph$distance_matrix,
    dim = 3L,
    vertex_ids = prepared$top_level_vertices,
    insertion_order = prepared$insertion_order[prepared$insertion_order %in% prepared$top_level_vertices],
    anchor_count = prepared$insertion_anchor_count,
    anchor_weight_mode = prepared$insertion_anchor_weight_mode
  )

  expect_true(all(is.finite(init$coords)))
  expect_equal(length(init$seed_vertices), length(prepared$top_level_vertices))
})

test_that("GMDS and GKK top-level solvers retain geometric initialization metadata", {
  edges <- edges.mesh(5, 5)

  prepared.gmds <- grip.prepare.misf.geodesic.mds(
    edges = edges,
    n = 25,
    num_init = 6L,
    dim = 2L,
    top_level_mode = "solve",
    top_level_init = "geometric",
    top_level_restarts = 2L,
    top_level_max_iter = 1L,
    seed = 23L
  )
  expect_identical(prepared.gmds$top_level_fit$top_level_init, "geometric")
  expect_true(all(is.finite(prepared.gmds$top_level_fit$coords)))
  expect_false(is.null(prepared.gmds$top_level_fit$initial_placement))
  expect_equal(
    length(prepared.gmds$top_level_fit$initial_placement$seed_vertices),
    min(3L, length(prepared.gmds$top_level_vertices))
  )

  prepared.gkk <- grip.prepare.misf.geodesic.kk(
    edges = edges,
    n = 25,
    num_init = 6L,
    dim = 2L,
    top_level_mode = "solve",
    top_level_init = "geometric",
    top_level_restarts = 2L,
    top_level_max_iter = 1L,
    seed = 23L
  )
  expect_identical(prepared.gkk$top_level_fit$top_level_init, "geometric")
  expect_true(all(is.finite(prepared.gkk$top_level_fit$coords)))
  expect_false(is.null(prepared.gkk$top_level_fit$initial_placement))
  expect_equal(
    length(prepared.gkk$top_level_fit$initial_placement$seed_vertices),
    min(3L, length(prepared.gkk$top_level_vertices))
  )
})

test_that("3D geometric coarse seed prefers a full-rank top-level seed when available", {
  bundle <- mesh.surface.graph(
    12, 12,
    surface = "paraboloid",
    amplitude = 0.35,
    connectivity = "orthogonal",
    normalize = "median"
  )
  prepared <- grip.prepare.misf.geodesic.mds(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average",
    dim = 3L,
    top_level_mode = "skip",
    seed = 12L
  )
  init <- grip:::grip.geodesic.misf.build.geometric.seed.coords(
    distance_matrix = prepared$top_level_graph$distance_matrix,
    dim = 3L,
    vertex_ids = prepared$top_level_vertices,
    insertion_order = prepared$insertion_order[prepared$insertion_order %in% prepared$top_level_vertices],
    anchor_count = prepared$insertion_anchor_count,
    anchor_weight_mode = prepared$insertion_anchor_weight_mode
  )

  seed.local <- match(init$seed_vertices, init$vertex_ids)
  seed.coords <- init$coords[seed.local, , drop = FALSE]
  expect_equal(length(init$seed_vertices), 4L)
  expect_true(qr(scale(seed.coords, scale = FALSE))$rank == 3L)
  expect_true(qr(scale(init$coords, scale = FALSE))$rank == 3L)
})
