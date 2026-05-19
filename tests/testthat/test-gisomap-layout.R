test_that("GISOMAP wrapper supports metric MDS and weighted-GRIP initializers in dimensions 2, 3, and 4", {
  edges <- edges.mesh(3L, 3L)
  edge_weights <- rep(1, nrow(edges))
  prepared <- grip.prepare.graph.geodesic.mds(
    edges = edges,
    n = 9L,
    edge_weights = edge_weights
  )

  for (initializer in c("metric_mds", "weighted_grip")) {
    for (dim in 2:4) {
      init_args <- if (identical(initializer, "weighted_grip")) {
        list(
          placement = "barycenter",
          rounds = 3L,
          final_rounds = 3L,
          num_init = max(5L, dim + 1L),
          num_nbrs = 5L,
          repulsion_factor = 0.4
        )
      } else {
        list()
      }
      fit <- grip.layout.gisomap(
        prepared = prepared,
        dim = dim,
        init = initializer,
        init_args = init_args,
        edge_kk_args = list(
          stiffness_method = "uniform",
          density_mix_schedule = 1,
          scale_mode = "identity",
          max_iter = 3L,
          return_trace = TRUE,
          engine = "cpp"
        ),
        diagnostics = TRUE,
        seed = 91L
      )

      expect_s3_class(fit, "grip_gisomap_layout")
      expect_s3_class(fit, "grip_gmds_layout")
      expect_equal(fit$method, "gisomap")
      expect_equal(dim(fit$coords), c(9L, dim))
      expect_true(all(is.finite(fit$coords)))
      expect_equal(fit$metadata$initializer, initializer)
      expect_equal(fit$metadata$polish, "edge_kk")
      expect_equal(fit$metadata$edge_kk_method, "edge_kk")
      expect_equal(fit$metadata$dim, dim)
      expect_true(is.data.frame(fit$diagnostics))
      expect_true(is.data.frame(fit$metadata$initial_diagnostics))
      expect_true(is.data.frame(fit$trace))
      expect_true(all(c("stage", "mix", "energy", "edge.rel.rmse") %in% names(fit$trace)))
    }
  }
})

test_that("GISOMAP wrapper is deterministic for weighted-GRIP initialization with a fixed seed", {
  edges <- edges.mesh(3L, 3L)
  edge_weights <- rep(1, nrow(edges))
  args <- list(
    edges = edges,
    n = 9L,
    edge_weights = edge_weights,
    dim = 4L,
    init = "weighted_grip",
    init_args = list(
      placement = "barycenter",
      rounds = 3L,
      final_rounds = 3L,
      num_init = 5L,
      num_nbrs = 5L,
      repulsion_factor = 0.4
    ),
    edge_kk_args = list(
      stiffness_method = "uniform",
      density_mix_schedule = 1,
      scale_mode = "identity",
      max_iter = 3L,
      return_trace = TRUE,
      engine = "cpp"
    ),
    seed = 92L
  )

  fit1 <- do.call(grip.layout.gisomap, args)
  fit2 <- do.call(grip.layout.gisomap, args)

  expect_equal(fit1$coords, fit2$coords, tolerance = 1e-12)
  expect_equal(fit1$trace$energy, fit2$trace$energy, tolerance = 1e-12)
})

test_that("GISOMAP wrapper validates dimensions and delegated argument lists", {
  edges <- edges.path(4L)
  edge_weights <- rep(1, nrow(edges))

  expect_error(
    grip.layout.gisomap(edges = edges, n = 4L, edge_weights = edge_weights, dim = 1L),
    "dim must be >= 2"
  )
  expect_error(
    grip.layout.gisomap(
      edges = edges,
      n = 4L,
      edge_weights = edge_weights,
      init_args = list(dim = 3L)
    ),
    "init_args must not override"
  )
  expect_error(
    grip.layout.gisomap(
      edges = edges,
      n = 4L,
      edge_weights = edge_weights,
      edge_kk_args = list(coords = matrix(0, 4L, 2L))
    ),
    "edge_kk_args must not override"
  )
  expect_error(
    grip.layout.gisomap(
      edges = edges,
      n = 4L,
      edge_weights = edge_weights,
      init_args = list(1L)
    ),
    "init_args must be a named list"
  )
})
