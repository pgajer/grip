test_that("legacy grip arguments preserve both metrics and preset overrides", {
  graph <- list(2L, c(1L, 3L), c(2L, 4L), 3L)
  lengths <- list(1, c(1, 2), c(2, 1), 1)
  for (metric in c("hop", "edge_length")) {
    current <- grip(adj.list = graph, weight.list = lengths, metric = metric,
      dim = 2, preset = "mesh", rounds = 2, final.rounds = 3, num.init = 4,
      level0.anchor.count = 2, seed = 71)
    expect_no_warning(previous <- grip(adj_list = graph, weight_list = lengths,
      metric = metric, dim = 2, preset = "mesh", rounds = 2, final_rounds = 3,
      num_init = 4, level0_anchor_count = 2, seed = 71))
    expect_identical(previous, current)
  }
  expect_error(grip(edges.path(4), metric_neighbor_cap = 2), "only available")
  expect_identical(
    grip(edges.path(4), n = 4, metric = "edge_length", edge_weights = c(1, 2, 1),
         metric_neighbor_cap = 2, length_normalization = "none",
         rounds = 1, final_rounds = 1),
    grip(edges.path(4), n = 4, metric = "edge_length", edge.weights = c(1, 2, 1),
         metric.neighbor.cap = 2, length.normalization = "none",
         rounds = 1, final.rounds = 1))
})

test_that("legacy edge-KK arguments survive downstream formal filtering", {
  graph <- list(2L, c(1L, 3L), c(2L, 4L), 3L)
  lengths <- list(1, c(1, 2), c(2, 1), 1)
  coords <- cbind(c(0, 1, 3, 4), c(0, 0.1, -0.1, 0))
  # Published geosmooth filters its settings by formals unless ... is present.
  old <- list(coords = coords, adj_list = graph, weight_list = lengths, dim = 2L,
    stiffness_method = "uniform", density_mix_schedule = 1, scale_mode = "identity",
    max_iter = 2L, return_trace = FALSE, diagnostics = FALSE, engine = "cpp")
  fmls <- names(formals(edge.kk))
  forwarded <- if ("..." %in% fmls) old else old[names(old) %in% fmls]
  expect_no_warning(previous <- do.call(edge.kk, forwarded))
  current <- edge.kk(coords = coords, adj.list = graph, weight.list = lengths,
    dim = 2L, stiffness.method = "uniform", density.mix.schedule = 1,
    scale.mode = "identity", max.iter = 2L, return.trace = FALSE,
    diagnostics = FALSE, engine = "cpp")
  expect_identical(previous, current)
  expect_identical(previous$method, "edge_kk")
})

test_that("legacy weighted initializer settings are translated safely", {
  old <- edge.kk(edges = edges.path(5), init = "weighted_grip",
    weighted_grip_args = list(rounds = 1, final_rounds = 2, num_init = 3),
    max_iter = 0, diagnostics = FALSE)
  current <- edge.kk(edges = edges.path(5), init = "weighted_grip",
    weighted.grip.args = list(rounds = 1, final.rounds = 2, num.init = 3),
    max.iter = 0, diagnostics = FALSE)
  expect_identical(old, current)
  expect_error(edge.kk(edges = edges.path(5), init = "weighted_grip",
    weighted_grip_args = list(edge_weights = rep(1, 4))), "must not include graph")
  expect_error(edge.kk(edges = edges.path(5), init = "weighted_grip",
    weighted_grip_args = list(final_rounds = 2, final.rounds = 2)), "Conflicting")
})

test_that("transition calls reject ambiguity and unknown names", {
  expect_error(grip(final.rounds = 1, final_rounds = 1), "Conflicting")
  expect_error(grip(final.r = 1, final_rounds = 1), "Conflicting")
  expect_error(grip(final_rounds = 1, final_rounds = 2), "Conflicting")
  expect_error(edge.kk(max.iter = 1, max_iter = 1), "Conflicting")
  expect_error(grip(final_round = 1), "unused argument")
  expect_error(edge.kk(coords = stop("must not force"), unknown_control = 1),
               "unused argument")
})

test_that("transition calls keep caller scope and do not force unused values", {
  run <- function(...) edge.kk(...)
  calls <- 0L
  coords <- cbind(0:3, c(0, 1, 0, 1))
  expected <- edge.kk(coords = coords, edges = edges.path(4), max.iter = 0,
                     diagnostics = FALSE)
  result <- run(coords = { calls <- calls + 1L; coords }, edges = edges.path(4),
                max_iter = 0, diagnostics = FALSE,
                weighted_grip_args = stop("unused when coords are supplied"))
  expect_identical(calls, 1L)
  expect_identical(result, expected)
})
