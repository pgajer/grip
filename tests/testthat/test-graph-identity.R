graph.identity.entry.points <- function(prepared = FALSE) {
  exports <- getNamespaceExports("grip")
  required <- c("edges", "n", "adj_list", "weight_list", "edge_weights")
  if (prepared) required <- c(required, "prepared")
  exports[vapply(exports, function(name) {
    all(required %in% names(formals(getExportedValue("grip", name))))
  }, logical(1))]
}

test_that("public graph entry points reject fractional identity before dispatch", {
  for (name in graph.identity.entry.points()) {
    fun <- getExportedValue("grip", name)
    for (id in c(1.9, 1 + 1e-10, NA_real_, NaN, Inf, 0, -1, 2^31)) {
      expect_error(do.call(fun, list(edges = rbind(c(id, 2), c(2, 3)), n = 3)),
                   "edges\\[1, 1\\].*integer vertex id", info = name)
    }
    for (n in list(3.9, 3 + 1e-10, NA_real_, Inf, 0, -1, 2^31, c(3, 4))) {
      expect_error(do.call(fun, list(edges = cbind(1:2, 2:3), n = n)),
                   "n must be a single finite positive integer", info = name)
    }
    expect_error(do.call(fun, list(adj_list = list(2.5, c(1, 3), 2), n = 3)),
                 "adj_list\\[\\[1\\]\\]\\[1\\].*integer vertex id", info = name)
  }
})

test_that("public graph entry points reject competing representations and misplaced lengths", {
  edges <- cbind(1:2, 2:3)
  adj <- list(2L, c(1L, 3L), 2L)
  for (name in graph.identity.entry.points()) {
    fun <- getExportedValue("grip", name)
    expect_error(do.call(fun, list(edges = edges, adj_list = adj, n = 3)),
                 "either edges or adj_list, not both", info = name)
    expect_error(do.call(fun, list(edges = edges, weight_list = list(1, c(1, 1), 1), n = 3)),
                 "weight_list requires adj_list", info = name)
    expect_error(do.call(fun, list(adj_list = adj, edge_weights = c(1, 1), n = 3)),
                 "edge_weights requires edges", info = name)
  }
})

test_that("inferred counts cannot hide fractional endpoint ids", {
  for (fun in list(prepare.edge.kk, prepare.geodesic.kk, prepare.graph.geodesic.mds,
                   prepare.landmark.geodesic.kk, classical.mds, metric.mds,
                   weighted.grip.nd)) {
    expect_error(fun(edges = rbind(c(1.5, 2), c(2, 3))), "edges\\[1, 1\\]")
  }
  expect_error(classical.mds(edges = rbind(c(1.5, 2), c(2, 3)), diagnostics = FALSE),
               "edges\\[1, 1\\]")
  expect_error(metric.mds(edges = rbind(c(1, 2), c(2, 3.5)), diagnostics = FALSE),
               "edges\\[2, 2\\]")
})

test_that("prepared graphs reject raw graph inputs and contradictory counts", {
  edges <- cbind(1:2, 2:3)
  prepared <- prepare.graph.geodesic.mds(edges, n = 3)
  for (name in graph.identity.entry.points(prepared = TRUE)) {
    fun <- getExportedValue("grip", name)
    for (raw in list(list(edges = edges), list(adj_list = list(2, c(1, 3), 2)),
                     list(edge_weights = c(1, 1)), list(weight_list = list(1, c(1, 1), 1)))) {
      expect_error(do.call(fun, c(list(prepared = prepared), raw)),
                   "prepared cannot be combined", info = name)
    }
    expect_error(do.call(fun, list(prepared = prepared, n = 4)),
                 "n must match the graph size stored in prepared", info = name)
    expect_error(do.call(fun, list(prepared = prepared, n = 3.1)),
                 "n must be a single finite positive integer", info = name)
  }
  coords <- cbind(0:2, 0)
  expect_identical(score.gmds(coords, prepared = prepared, n = 3),
                   score.gmds(coords, prepared = prepared))
  expect_identical(classical.mds(prepared = prepared, n = 3),
                   classical.mds(prepared = prepared))
})

test_that("valid numeric ids retain graph order, isolates, and layout behavior", {
  edges <- cbind(1:2, 2:3)
  numeric.edges <- matrix(as.double(edges), ncol = 2)
  adj <- list(2L, c(1L, 3L), 2L, integer())
  args <- list(n = 4, dim = 2, rounds = 3, final_rounds = 3, seed = 19)
  expect_warning(from.edges <- do.call(grip, c(list(edges = edges), args)), "2 connected components")
  expect_warning(from.numeric <- do.call(grip, c(list(edges = numeric.edges), args)), "2 connected components")
  expect_warning(from.adj <- do.call(grip, c(list(adj_list = adj), args)), "2 connected components")
  expect_identical(from.numeric, from.edges)
  expect_identical(from.adj, from.edges)
  expect_equal(dim(from.edges), c(4, 2))
  expect_true(all(is.finite(from.edges)))
  expect_warning(empty <- grip(matrix(integer(), ncol = 2), n = 2, dim = 2), "2 connected components")
  expect_equal(dim(empty), c(2, 2))

  weighted <- list(n = 4, edge_weights = c(1, 2), dim = 2, seed = 19,
                   rounds = 3, final_rounds = 3, metric = "edge_length")
  expect_warning(weighted.numeric <- do.call(grip, c(list(edges = numeric.edges), weighted)), "2 connected components")
  expect_warning(weighted.integer <- do.call(grip, c(list(edges = edges), weighted)), "2 connected components")
  expect_identical(weighted.numeric, weighted.integer)
  expect_identical(prepare.edge.kk(numeric.edges), prepare.edge.kk(edges))
  expect_identical(prepare.graph.geodesic.mds(numeric.edges),
                   prepare.graph.geodesic.mds(edges))
})

test_that("star graph bundles cannot override another graph or carry fractional ids", {
  edges <- cbind(1:2, 2:3)
  coords <- cbind(0:2, 0)
  graph <- list(n = 3, edges = edges, edge_weights = c(1, 1))
  expect_error(graph.riemannian.star.structure(graph, coords, edges = edges),
               "graph cannot be combined")
  expect_error(graph.riemannian.star.structure(graph, coords, n = 3.5), "n must")
  expect_error(graph.riemannian.star.structure(graph, coords, n = 4), "n must match")
  graph$edges[1, 1] <- 1.5
  expect_error(graph.riemannian.star.structure(graph, coords), "edges\\[1, 1\\]")
  graph$edges <- edges
  graph$adj_list <- list(2L, c(1L, 3L), 2L)
  graph$weight_list <- list(1, c(1, 1), 1)
  expect_equal(graph.riemannian.star.structure(graph, coords)$n, 3)
  graph$edges <- edges.cycle(3)
  graph$edge_weights <- rep(1, 3)
  expect_error(graph.riemannian.star.structure(graph, coords), "contradictory edge and adjacency")
})
