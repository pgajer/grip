test_that("weighted ND layout returns finite deterministic coordinates above 3D", {
  edges <- edges.mesh(3, 4)
  weights <- rep(1, nrow(edges))

  coords1 <- weighted.grip.nd(
    edges = edges,
    edge_weights = weights,
    n = 12,
    dim = 5,
    rounds = 8,
    final_rounds = 8,
    num_init = 6,
    num_nbrs = 8,
    repulsion_factor = 0.4,
    seed = 71
  )
  coords2 <- weighted.grip.nd(
    edges = edges,
    edge_weights = weights,
    n = 12,
    dim = 5,
    rounds = 8,
    final_rounds = 8,
    num_init = 6,
    num_nbrs = 8,
    repulsion_factor = 0.4,
    seed = 71
  )

  expect_equal(dim(coords1), c(12, 5))
  expect_identical(colnames(coords1), paste0("Dim", 1:5))
  expect_true(all(is.finite(coords1)))
  expect_identical(coords1, coords2)
})

test_that("weighted ND layout remains available in 2D and 3D without using legacy validators", {
  edges <- edges.cycle(8)
  weights <- seq_len(nrow(edges)) / nrow(edges) + 1

  coords2 <- weighted.grip.nd(
    edges = edges,
    edge_weights = weights,
    n = 8,
    dim = 2,
    rounds = 4,
    final_rounds = 4,
    num_init = 4,
    seed = 73
  )
  coords3 <- weighted.grip.nd(
    edges = edges,
    edge_weights = weights,
    n = 8,
    dim = 3,
    rounds = 4,
    final_rounds = 4,
    num_init = 5,
    seed = 73
  )

  expect_equal(dim(coords2), c(8, 2))
  expect_equal(dim(coords3), c(8, 3))
  expect_true(all(is.finite(coords2)))
  expect_true(all(is.finite(coords3)))
})

test_that("weighted ND layout responds to nonuniform edge lengths", {
  edges <- edges.path(7)
  weights_flat <- rep(1, nrow(edges))
  weights_shaped <- c(1, 1, 4, 4, 1, 1)

  coords_flat <- weighted.grip.nd(
    edges = edges,
    edge_weights = weights_flat,
    n = 7,
    dim = 4,
    rounds = 12,
    final_rounds = 12,
    num_init = 5,
    repulsion_factor = 0,
    seed = 79
  )
  coords_shaped <- weighted.grip.nd(
    edges = edges,
    edge_weights = weights_shaped,
    n = 7,
    dim = 4,
    rounds = 12,
    final_rounds = 12,
    num_init = 5,
    repulsion_factor = 0,
    seed = 79
  )

  expect_gt(max(abs(coords_flat - coords_shaped)), 1e-6)
})

test_that("weighted ND MISF skeleton matches legacy weighted MISF ordering", {
  edges <- edges.mesh(4, 4)
  weights <- seq_len(nrow(edges)) / nrow(edges) + 1

  legacy <- build.weighted.misf(
    edges = edges,
    edge_weights = weights,
    n = 16,
    num_init = 6,
    num_nbrs = 5,
    seed = 83
  )
  nd <- grip:::grip.build.misf.weighted.nd(
    edges = edges,
    edge_weights = weights,
    n = 16,
    num_init = 6,
    num_nbrs = 5,
    seed = 83
  )

  expect_identical(nd$mish_order, legacy$mish_order)
  expect_identical(nd$misf_size, legacy$misf_size)
  expect_identical(nd$vertex_depth, as.integer(legacy$vertex_depth))
  expect_identical(nd$num_nbrs_schedule, as.integer(legacy$num_nbrs_schedule))
  expect_equal(nd$levels, lapply(legacy$levels, as.integer))
})

test_that("weighted ND trace records initialization, refinement, and final frames", {
  edges <- edges.mesh(3, 3)
  trace <- grip:::grip.layout.weighted.nd.trace(
    edges = edges,
    edge_weights = rep(1, nrow(edges)),
    n = 9,
    dim = 4,
    rounds = 3,
    final_rounds = 2,
    num_init = 5,
    num_nbrs = 6,
    trace.every = 2,
    seed = 81
  )

  expect_s3_class(trace, "grip_layout_weighted_nd_trace")
  expect_equal(dim(trace$final), c(9, 4))
  expect_identical(colnames(trace$final), paste0("Dim", 1:4))
  expect_identical(trace$meta$phase[[1]], "init")
  expect_identical(tail(trace$meta$phase, 1), "final")
  expect_true(any(trace$meta$phase == "level_start"))
  expect_true(any(trace$meta$phase == "round"))
  expect_equal(length(trace$frames), nrow(trace$meta))
  expect_equal(trace$frames[[length(trace$frames)]], trace$final)
  expect_true(all(is.finite(trace$final)))
  expect_true(all(vapply(trace$frames, function(x) {
    all(is.finite(x) | is.na(x))
  }, logical(1))))
  expect_equal(trace$meta$active_vertices, c(5L, 5L, 9L, 9L, 9L))
  expect_equal(trace$meta$round_in_level, c(0L, 2L, 0L, 2L, 2L))
})

test_that("weighted ND trace mirrors legacy weighted-GRIP active-level phases", {
  edges <- edges.mesh(3, 3)
  weights <- rep(1, nrow(edges))
  args <- list(
    edges = edges,
    edge_weights = weights,
    n = 9,
    dim = 2,
    rounds = 3,
    final_rounds = 2,
    num_init = 5,
    num_nbrs = 6,
    trace.every = 1,
    seed = 81
  )
  legacy <- do.call(trace.legacy.grip, c(args, list(trace = "round")))
  nd <- do.call(grip:::grip.layout.weighted.nd.trace, args)

  expect_identical(nd$meta$phase, legacy$meta$phase)
  expect_identical(nd$meta$level_index, legacy$meta$level_index)
  expect_identical(nd$meta$misf_level, legacy$meta$misf_level)
  expect_identical(nd$meta$round_in_level, legacy$meta$round_in_level)
  expect_identical(nd$meta$active_vertices, legacy$meta$active_vertices)
  expect_equal(
    vapply(nd$frames, function(x) sum(stats::complete.cases(x)), integer(1)),
    legacy$meta$active_vertices
  )
})

test_that("weighted ND trace initializes only top vertices and inserts at level boundaries", {
  edges <- edges.mesh(3, 3)
  trace <- grip:::grip.layout.weighted.nd.trace(
    edges = edges,
    edge_weights = rep(1, nrow(edges)),
    n = 9,
    dim = 3,
    rounds = 3,
    final_rounds = 2,
    num_init = 5,
    num_nbrs = 6,
    trace.every = 1,
    seed = 81
  )

  init <- trace$frames[[1L]]
  before.level <- trace$frames[[which(trace$meta$phase == "round" &
                                        trace$meta$level_index == 1L)[3L]]]
  level.start <- trace$frames[[which(trace$meta$phase == "level_start")[1L]]]
  init.active <- stats::complete.cases(init)
  opened.active <- stats::complete.cases(level.start)

  expect_equal(sum(init.active), 5L)
  expect_equal(sum(opened.active), 9L)
  expect_true(all(is.na(init[!init.active, , drop = FALSE])))
  expect_true(all(is.finite(level.start[opened.active, , drop = FALSE])))
  expect_equal(
    level.start[init.active, , drop = FALSE],
    before.level[init.active, , drop = FALSE],
    tolerance = 1e-12
  )
  expect_true(any(abs(level.start[!init.active, , drop = FALSE]) > 0))
})

test_that("weighted ND insertion polish moves new vertices beyond raw anchor barycenters", {
  edges <- edges.mesh(3, 3)
  weights <- rep(1, nrow(edges))
  trace <- grip:::grip.layout.weighted.nd.trace(
    edges = edges,
    edge_weights = weights,
    n = 9,
    dim = 2,
    rounds = 3,
    final_rounds = 2,
    num_init = 5,
    num_nbrs = 6,
    trace.every = 1,
    seed = 81
  )
  misf <- grip:::grip.build.misf.weighted.nd(
    edges = edges,
    edge_weights = weights,
    n = 9,
    num_init = 5,
    num_nbrs = 6,
    seed = 81
  )

  before.level <- trace$frames[[which(trace$meta$phase == "round" &
                                        trace$meta$level_index == 1L)[3L]]]
  level.start <- trace$frames[[which(trace$meta$phase == "level_start")[1L]]]
  previous.active <- stats::complete.cases(before.level)
  inserted <- which(!previous.active & stats::complete.cases(level.start))

  adj <- vector("list", 9L)
  for (i in seq_len(nrow(edges))) {
    u <- edges[i, 1L]
    v <- edges[i, 2L]
    adj[[u]] <- c(adj[[u]], v)
    adj[[v]] <- c(adj[[v]], u)
  }
  graph.dist <- function(root) {
    dist <- rep(Inf, 9L)
    dist[root] <- 0
    queue <- root
    while (length(queue)) {
      u <- queue[[1L]]
      queue <- queue[-1L]
      for (v in adj[[u]]) {
        if (!is.finite(dist[v])) {
          dist[v] <- dist[u] + 1
          queue <- c(queue, v)
        }
      }
    }
    dist
  }

  moved <- vapply(inserted, function(root) {
    eligible <- which(previous.active & misf$vertex_depth > misf$vertex_depth[root])
    d <- graph.dist(root)
    anchors <- eligible[order(d[eligible], eligible)][seq_len(min(3L, length(eligible)))]
    raw.barycenter <- colMeans(before.level[anchors, , drop = FALSE])
    max(abs(level.start[root, ] - raw.barycenter)) > 1e-8
  }, logical(1))

  expect_true(any(moved))
})

test_that("weighted ND insertion placement and anchor controls are usable", {
  edges <- edges.mesh(3, 3)
  weights <- rep(1, nrow(edges))

  circle <- weighted.grip.nd(
    edges = edges,
    edge_weights = weights,
    n = 9,
    dim = 2,
    placement = "circle",
    rounds = 3,
    final_rounds = 2,
    num_init = 5,
    num_nbrs = 6,
    seed = 91
  )
  tuned <- weighted.grip.nd(
    edges = edges,
    edge_weights = weights,
    n = 9,
    dim = 3,
    rounds = 3,
    final_rounds = 2,
    num_init = 5,
    num_nbrs = 6,
    insertion_anchor_count = 4,
    insertion_anchor_scope = "prev_misf",
    insertion_anchor_strategy = "balanced_band",
    level0_insertion_mode = "least_squares",
    level0_anchor_count = 4,
    level0_local_kk_steps = 4,
    seed = 91
  )
  spread <- grip:::grip.layout.weighted.nd.trace(
    edges = edges,
    edge_weights = weights,
    n = 9,
    dim = 2,
    rounds = 3,
    final_rounds = 2,
    num_init = 5,
    num_nbrs = 6,
    insertion_anchor_scope = "prev_misf",
    insertion_anchor_strategy = "spread_prev",
    level0_insertion_mode = "barycenter",
    level0_local_kk_steps = 0,
    trace.every = 1,
    seed = 91
  )

  expect_equal(dim(circle), c(9, 2))
  expect_equal(dim(tuned), c(9, 3))
  expect_s3_class(spread, "grip_layout_weighted_nd_trace")
  expect_true(all(is.finite(circle)))
  expect_true(all(is.finite(tuned)))
  expect_true(all(is.finite(spread$final)))
  expect_error(
    weighted.grip.nd(
      edges = edges,
      edge_weights = weights,
      n = 9,
      dim = 3,
      placement = "circle",
      rounds = 1,
      final_rounds = 1,
      num_init = 5,
      seed = 91
    ),
    "placement = 'circle' requires dim = 2"
  )
})

test_that("weighted ND legacy final-stage and metric-search controls are usable", {
  edges <- edges.mesh(3, 4)
  weights <- seq_len(nrow(edges)) / nrow(edges) + 1

  baseline <- weighted.grip.nd(
    edges = edges,
    edge_weights = weights,
    n = 12,
    dim = 3,
    rounds = 4,
    final_rounds = 4,
    num_init = 5,
    num_nbrs = 6,
    seed = 93
  )
  scaled <- weighted.grip.nd(
    edges = edges,
    edge_weights = weights,
    n = 12,
    dim = 3,
    rounds = 4,
    final_rounds = 4,
    num_init = 5,
    num_nbrs = 6,
    final_move_scale_after_first = 0.5,
    metric_neighbor_cap = 8,
    seed = 93
  )
  kk <- weighted.grip.nd(
    edges = edges,
    edge_weights = weights,
    n = 12,
    dim = 3,
    rounds = 4,
    final_rounds = 4,
    num_init = 5,
    num_nbrs = 6,
    final_mode = "kk_repulse",
    seed = 93
  )

  expect_equal(dim(scaled), c(12, 3))
  expect_equal(dim(kk), c(12, 3))
  expect_true(all(is.finite(scaled)))
  expect_true(all(is.finite(kk)))
  expect_gt(max(abs(baseline - scaled)), 1e-8)
  expect_gt(max(abs(baseline - kk)), 1e-8)
})

test_that("weighted ND layout packs disconnected components", {
  edges <- matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE)
  coords <- weighted.grip.nd(
    edges = edges,
    edge_weights = c(1, 2),
    n = 4,
    dim = 4,
    rounds = 4,
    final_rounds = 4,
    num_init = 5,
    disconnected = "components",
    seed = 83
  )

  expect_equal(dim(coords), c(4, 4))
  expect_true(all(is.finite(coords)))
  expect_error(
    weighted.grip.nd(
      edges = edges,
      edge_weights = c(1, 2),
      n = 4,
      dim = 4,
      rounds = 1,
      final_rounds = 1,
      num_init = 5,
      disconnected = "error",
      seed = 83
    ),
    "graph must be connected"
  )
})

test_that("weighted ND layout keeps weighted presets usable with ND constraints", {
  edges <- edges.kary.tree(k = 2, depth = 3)
  coords <- weighted.grip.nd(
    edges = edges,
    edge_weights = rep(1, nrow(edges)),
    n = max(edges),
    dim = 4,
    preset = "tree",
    rounds = 4,
    final_rounds = 4,
    seed = 87
  )

  expect_equal(dim(coords), c(max(edges), 4))
  expect_true(all(is.finite(coords)))
  expect_error(
    weighted.grip.nd(
      edges = edges,
      edge_weights = rep(1, nrow(edges)),
      n = max(edges),
      dim = 4,
      placement = "circle",
      rounds = 1,
      final_rounds = 1,
      seed = 87
    ),
    "placement = 'circle' requires dim = 2"
  )
})

test_that("weighted ND layout validates dimensions and required weights", {
  edges <- edges.path(5)

  expect_error(
    weighted.grip.nd(edges = edges, n = 5, dim = 4),
    "requires edge weights"
  )
  expect_error(
    weighted.grip.nd(
      edges = edges,
      edge_weights = rep(1, nrow(edges)),
      n = 5,
      dim = 1
    ),
    "dim must be >= 2"
  )
  expect_error(
    weighted.grip.nd(
      edges = edges,
      edge_weights = rep(1, nrow(edges)),
      n = 5,
      dim = 3.5
    ),
    "dim must be an integer"
  )
  expect_error(
    weighted.grip.nd(
      edges = edges,
      edge_weights = rep(1, nrow(edges)),
      n = 5,
      dim = 6,
      num_init = 6
    ),
    "num_init must be >= 7"
  )
  expect_error(
    weighted.grip.nd(
      edges = edges,
      edge_weights = rep(1, nrow(edges)),
      n = 5,
      dim = 3,
      final_move_scale_after_first = 1.1
    ),
    "final_move_scale_after_first must be in \\[0, 1\\]"
  )
  expect_error(
    weighted.grip.nd(
      edges = edges,
      edge_weights = rep(1, nrow(edges)),
      n = 5,
      dim = 3,
      final_mode = "banana"
    ),
    "'arg' should be one of"
  )
  expect_error(
    weighted.grip.nd(
      edges = edges,
      edge_weights = rep(1, nrow(edges)),
      n = 5,
      dim = 3,
      metric_neighbor_cap = 0
    ),
    "metric_neighbor_cap must be a positive integer"
  )
})

test_that("legacy weighted GRIP dimensionality remains capped at 3D", {
  edges <- edges.path(6)
  expect_error(
    globalrep.weighted.grip(
      edges = edges,
      edge_weights = rep(1, nrow(edges)),
      n = 6,
      dim = 4,
      rounds = 2,
      final_rounds = 2,
      num_init = 5,
      seed = 89
    ),
    "dim must be 2 or 3"
  )
})
