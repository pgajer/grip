test_that("geodesic MDS scoring is exact on a weighted path realization", {
  edges <- matrix(c(
    1L, 2L,
    2L, 3L
  ), ncol = 2L, byrow = TRUE)
  coords <- cbind(
    x = c(0, 1, 3),
    y = c(0, 0, 0)
  )
  prepared <- grip.prepare.geodesic.kk(
    edges = edges,
    n = 3L,
    edge_weights = c(1, 2)
  )
  score <- grip.score.geodesic.mds(coords, prepared = prepared)

  expect_equal(score$n.pairs[[1L]], 3L)
  expect_lt(score$gmds.energy[[1L]], 1e-10)
  expect_lt(score$gmds.raw_stress[[1L]], 1e-10)
  expect_lt(score$gmds.stress[[1L]], 1e-10)
})

test_that("geodesic MDS preparation is deterministic and augments disconnected k-NN graphs", {
  data <- rbind(
    c(0, 0),
    c(0, 1),
    c(10, 0),
    c(10, 1)
  )

  prepared1 <- grip.prepare.geodesic.mds(data = data, k = 1L, connect = "mst")
  prepared2 <- grip.prepare.geodesic.mds(data = data, k = 1L, connect = "mst")

  expect_s3_class(prepared1, "grip_gmds_prepared")
  expect_equal(prepared1$edges, prepared2$edges)
  expect_equal(prepared1$mst_added_edges, prepared2$mst_added_edges)
  expect_equal(nrow(prepared1$mst_added_edges), 1L)
  expect_true(all(is.finite(prepared1$pair_graph_distance)))

  adj <- grip:::grip.build.adj.from.edges(prepared1$edges, n = nrow(data))$adj_list
  comp <- grip:::grip.connected.components(adj, nrow(data))
  expect_equal(length(unique(comp)), 1L)
})

test_that("graph-first geodesic MDS preparation matches the data-native wrapper", {
  data <- rbind(
    c(0, 0),
    c(1, 0),
    c(2, 0.2),
    c(0.2, 1.1),
    c(1.3, 1.2)
  )

  built <- grip:::grip.prepare.geodesic.mds.graph(data = data, k = 2L, connect = "mst")
  prepared.graph <- grip.prepare.graph.geodesic.mds(
    edges = built$edges,
    n = nrow(data),
    edge_weights = built$edge_weights,
    tie_mode = "average"
  )
  prepared.data <- grip.prepare.geodesic.mds(
    data = data,
    k = 2L,
    connect = "mst",
    tie_mode = "average"
  )

  expect_s3_class(prepared.graph, "grip_gmds_prepared")
  expect_identical(prepared.graph$graph_build_mode, "graph_input")
  expect_identical(prepared.data$graph_build_mode, "symmetric_knn")

  common.fields <- c(
    "n", "edges", "edge_targets", "pair_matrix", "pair_graph_distance",
    "path_vertices", "path_edges", "path_edge_weights", "pair_path_count_log",
    "flat_pair_edge_offsets", "flat_edge_u", "flat_edge_v", "flat_edge_coeff",
    "graph_diameter", "distance_matrix", "pair_mode", "tie_mode"
  )
  for (field in common.fields) {
    expect_equal(prepared.graph[[field]], prepared.data[[field]], info = field)
  }

  coords <- cbind(
    x = seq_len(nrow(data)),
    y = c(0, 1, -0.5, 0.75, -1.25)
  )
  score.graph <- grip.score.geodesic.mds(coords, prepared = prepared.graph)
  score.data <- grip.score.geodesic.mds(coords, prepared = prepared.data)
  expect_equal(score.graph$gmds.raw_stress[[1L]], score.data$gmds.raw_stress[[1L]])
  expect_equal(score.graph$gmds.energy[[1L]], score.data$gmds.energy[[1L]])
})

test_that("compiled geodesic MDS optimizer decreases path stress on a perturbed path", {
  edges <- edges.path(5L)
  canonical <- cbind(
    x = 0:4,
    y = rep(0, 5)
  )
  prepared <- grip.prepare.geodesic.kk(edges = edges, n = 5L)

  perturbed <- canonical
  set.seed(41)
  perturbed <- perturbed + matrix(rnorm(length(perturbed), sd = 0.3), ncol = 2L)

  before <- grip.score.geodesic.mds(perturbed, prepared = prepared)
  opt <- grip.optimize.geodesic.mds(
    coords = perturbed,
    prepared = prepared,
    max_iter = 8L,
    return_trace = TRUE,
    engine = "cpp"
  )
  after <- grip.score.geodesic.mds(opt$coords, prepared = prepared)

  expect_lt(after$gmds.raw_stress[[1L]], before$gmds.raw_stress[[1L]])
  expect_true(nrow(opt$trace) >= 2L)
  expect_true(length(opt$frames) >= 2L)
})

test_that("geodesic MDS can optimize directly from data with cmdscale initialization", {
  data <- cbind(
    x = c(0, 1, 2, 3, 4),
    y = c(0, 0.1, -0.1, 0.1, 0)
  )

  opt <- grip.optimize.geodesic.mds(
    data = data,
    k = 1L,
    dim = 2L,
    max_iter = 3L,
    return_trace = TRUE
  )

  expect_s3_class(opt$prepared, "grip_gmds_prepared")
  expect_true(is.data.frame(opt$trace))
  expect_true(length(opt$frames) >= 1L)
  expect_true(is.finite(opt$score$gmds.stress[[1L]]))
})
