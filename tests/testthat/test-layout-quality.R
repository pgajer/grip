test_that("grip.score.layout returns the expected metrics", {
  edges <- edges.mesh(4, 4)
  coords <- grip.layout(edges, n = 16, dim = 2, preset = "mesh", seed = 1)
  score <- grip.score.layout(coords, edges = edges, n = 16)

  expect_s3_class(score, "data.frame")
  expect_equal(nrow(score), 1)
  expect_true(all(c(
    "n.vertices",
    "n.edges",
    "dim",
    "sampled.stress",
    "edge.length.cv",
    "median.edge.length",
    "sampled.nonedge.sep.ratio",
    "edge.crossings",
    "cluster.separation"
  ) %in% names(score)))
  expect_true(all(is.finite(unlist(score[1, c(
    "sampled.stress",
    "edge.length.cv",
    "median.edge.length",
    "sampled.nonedge.sep.ratio"
  )]))))
})

test_that("grip.score.layout counts 2D edge crossings exactly", {
  coords <- matrix(c(
    0, 0,
    1, 1,
    0, 1,
    1, 0
  ), ncol = 2, byrow = TRUE)
  edges <- rbind(
    c(1, 2),
    c(3, 4)
  )

  score <- grip.score.layout(
    coords = coords,
    edges = edges,
    n = 4,
    edge.crossings = "always"
  )

  expect_equal(score$edge.crossings[[1L]], 1)
})

test_that("grip.score.layout computes cluster separation when clusters are supplied", {
  coords <- matrix(c(
    0.0, 0.0,
    0.2, 0.1,
    5.0, 0.0,
    5.1, 0.2
  ), ncol = 2, byrow = TRUE)
  edges <- rbind(
    c(1, 2),
    c(3, 4)
  )

  score <- grip.score.layout(
    coords = coords,
    edges = edges,
    n = 4,
    clusters = c("A", "A", "B", "B"),
    edge.crossings = "always"
  )

  expect_true(is.finite(score$cluster.separation[[1L]]))
  expect_gt(score$cluster.separation[[1L]], 1)
})

test_that("grip.compare.layouts summarizes runs and seed stability", {
  edges <- edges.mesh(4, 4)
  cmp <- grip.compare.layouts(
    edges = edges,
    n = 16,
    dim = 2,
    candidates = c("default", "mesh"),
    seeds = 1:2,
    return.layouts = TRUE
  )

  expect_s3_class(cmp$runs, "data.frame")
  expect_s3_class(cmp$summary, "data.frame")
  expect_equal(sort(cmp$summary$candidate), c("default", "mesh"))
  expect_true(all(c(
    "sampled.stress.mean",
    "edge.length.cv.mean",
    "sampled.nonedge.sep.ratio.mean",
    "stability.procrustes.mean",
    "score.composite"
  ) %in% names(cmp$summary)))
  expect_true(all(vapply(cmp$layouts, length, integer(1L)) == 2L))
  mesh.row <- cmp$summary[cmp$summary$candidate == "mesh", , drop = FALSE]
  default.row <- cmp$summary[cmp$summary$candidate == "default", , drop = FALSE]
  expect_lt(mesh.row$score.composite[[1L]], default.row$score.composite[[1L]])
  expect_true(is.finite(mesh.row$stability.procrustes.mean[[1L]]))
})

test_that("grip.compare.layouts expands search grids into candidates", {
  edges <- edges.mesh(4, 4)
  cmp <- grip.compare.layouts(
    edges = edges,
    n = 16,
    dim = 2,
    search = list(
      candidate.prefix = "mesh.search",
      preset = "mesh",
      final_rounds = c(96L, 128L),
      num_nbrs = c(16L, 20L),
      repulsion_factor = c(1.0, 1.5)
    ),
    seeds = 1
  )

  expect_equal(nrow(cmp$summary), 8L)
  expect_true(all(grepl("^mesh.search", cmp$summary$candidate)))
  expect_true(all(c(
    "placement",
    "rounds",
    "final.rounds",
    "num.init",
    "num.nbrs",
    "r",
    "s",
    "repulsion.factor",
    "tinit.factor"
  ) %in% names(cmp$summary)))
  expect_true(all(sort(unique(cmp$summary$repulsion.factor)) == c(1.0, 1.5)))
  expect_true(all(sort(unique(cmp$summary$num.nbrs)) == c(16L, 20L)))
  expect_true(all(sort(unique(cmp$summary$final.rounds)) == c(96L, 128L)))
})
