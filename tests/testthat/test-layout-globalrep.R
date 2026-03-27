test_that("globalrep layout returns a finite matrix", {
  edges <- edges.path(10)
  coords <- grip.layout.globalrep(edges, n = 10, dim = 2,
                                  placement = "barycenter",
                                  rounds = 5, final_rounds = 3,
                                  num_init = 5, num_nbrs = 6,
                                  coarse_repulsion_factor = 0.2,
                                  coarse_repulsion_sample = 8,
                                  coarse_repulsion_exact_below = 32,
                                  seed = 123)
  expect_equal(dim(coords), c(10, 2))
  expect_true(all(is.finite(coords)))
})

test_that("globalrep seeded runs are deterministic", {
  edges <- edges.cycle(12)
  coords1 <- grip.layout.globalrep(edges, n = 12, dim = 2,
                                   placement = "barycenter",
                                   rounds = 4, final_rounds = 2,
                                   num_init = 4, num_nbrs = 5,
                                   coarse_repulsion_factor = 0.2,
                                   coarse_repulsion_sample = 8,
                                   coarse_repulsion_exact_below = 32,
                                   seed = 42)
  coords2 <- grip.layout.globalrep(edges, n = 12, dim = 2,
                                   placement = "barycenter",
                                   rounds = 4, final_rounds = 2,
                                   num_init = 4, num_nbrs = 5,
                                   coarse_repulsion_factor = 0.2,
                                   coarse_repulsion_sample = 8,
                                   coarse_repulsion_exact_below = 32,
                                   seed = 42)
  expect_identical(coords1, coords2)
})

test_that("globalrep with zero coarse repulsion matches grip.layout", {
  edges <- edges.mesh(5, 5)
  coords_base <- grip.layout(edges, n = 25, dim = 2,
                             placement = "barycenter",
                             rounds = 8, final_rounds = 8,
                             num_init = 6, num_nbrs = 8,
                             repulsion_factor = 1.5,
                             seed = 29)
  coords_globalrep <- grip.layout.globalrep(edges, n = 25, dim = 2,
                                            placement = "barycenter",
                                            rounds = 8, final_rounds = 8,
                                            num_init = 6, num_nbrs = 8,
                                            repulsion_factor = 1.5,
                                            coarse_repulsion_factor = 0,
                                            coarse_repulsion_sample = 8,
                                            coarse_repulsion_exact_below = 32,
                                            seed = 29)
  expect_identical(coords_base, coords_globalrep)
})

test_that("globalrep coarse repulsion changes the layout with a fixed seed", {
  edges <- edges.mesh(5, 5)
  coords_none <- grip.layout.globalrep(edges, n = 25, dim = 2,
                                       placement = "barycenter",
                                       rounds = 8, final_rounds = 8,
                                       num_init = 6, num_nbrs = 8,
                                       coarse_repulsion_factor = 0,
                                       coarse_repulsion_sample = 8,
                                       coarse_repulsion_exact_below = 32,
                                       seed = 31)
  coords_more <- grip.layout.globalrep(edges, n = 25, dim = 2,
                                       placement = "barycenter",
                                       rounds = 8, final_rounds = 8,
                                       num_init = 6, num_nbrs = 8,
                                       coarse_repulsion_factor = 0.3,
                                       coarse_repulsion_sample = 8,
                                       coarse_repulsion_exact_below = 32,
                                       seed = 31)
  expect_gt(max(abs(coords_none - coords_more)), 1e-6)
})

test_that("globalrep validates the new tuning parameters", {
  edges <- edges.cycle(10)
  expect_error(
    grip.layout.globalrep(edges, n = 10, dim = 2,
                          coarse_repulsion_factor = -0.1,
                          seed = 1),
    "coarse_repulsion_factor must be >= 0"
  )
  expect_error(
    grip.layout.globalrep(edges, n = 10, dim = 2,
                          coarse_repulsion_sample = 0,
                          seed = 1),
    "coarse_repulsion_sample must be a positive integer"
  )
  expect_error(
    grip.layout.globalrep(edges, n = 10, dim = 2,
                          coarse_repulsion_exact_below = 0,
                          seed = 1),
    "coarse_repulsion_exact_below must be a positive integer"
  )
})

test_that("globalrep disconnected handling matches grip.layout when disabled", {
  edges <- rbind(
    cbind(1:2, 2:3),
    cbind(5:6, 6:7)
  )
  expect_warning(
    coords_base <- grip.layout(edges, n = 7, dim = 2, seed = 11),
    "laying out components separately"
  )
  expect_warning(
    coords_globalrep <- grip.layout.globalrep(edges, n = 7, dim = 2,
                                              coarse_repulsion_factor = 0,
                                              coarse_repulsion_sample = 8,
                                              coarse_repulsion_exact_below = 32,
                                              seed = 11),
    "laying out components separately"
  )
  expect_identical(coords_base, coords_globalrep)
})
