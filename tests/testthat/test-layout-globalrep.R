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

test_that("globalrep with zero coarse repulsion matches grip.layout.legacy", {
  edges <- edges.mesh(5, 5)
  coords_base <- grip.layout.legacy(edges, n = 25, dim = 2,
                                    placement = "barycenter",
                                    rounds = 8, final_rounds = 8,
                                    num_init = 6, num_nbrs = 8,
                                    r = 0.15, s = 3.0,
                                    repulsion_factor = 1.5,
                                    seed = 29)
  coords_globalrep <- grip.layout.globalrep(edges, n = 25, dim = 2,
                                            placement = "barycenter",
                                            rounds = 8, final_rounds = 8,
                                            num_init = 6, num_nbrs = 8,
                                            r = 0.15, s = 3.0,
                                            repulsion_factor = 1.5,
                                            coarse_repulsion_factor = 0,
                                            coarse_repulsion_sample = 8,
                                            coarse_repulsion_exact_below = 32,
                                            seed = 29)
  expect_identical(coords_base, coords_globalrep)
})

test_that("grip.layout is an alias of grip.layout.globalrep", {
  edges <- edges.mesh(5, 5)
  coords_primary <- grip.layout(edges, n = 25, dim = 2, seed = 17)
  coords_alias <- grip.layout.globalrep(edges, n = 25, dim = 2, seed = 17)
  expect_identical(coords_primary, coords_alias)
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

test_that("globalrep final_mode can switch the final stage", {
  edges <- edges.mesh(5, 5)
  coords_fr <- grip.layout.globalrep(edges, n = 25, dim = 2,
                                     placement = "barycenter",
                                     rounds = 8, final_rounds = 8,
                                     num_init = 6, num_nbrs = 8,
                                     coarse_repulsion_factor = 0.3,
                                     coarse_repulsion_sample = 8,
                                     coarse_repulsion_exact_below = 32,
                                     final_mode = "fr",
                                     seed = 31)
  coords_kkr <- grip.layout.globalrep(edges, n = 25, dim = 2,
                                      placement = "barycenter",
                                      rounds = 8, final_rounds = 8,
                                      num_init = 6, num_nbrs = 8,
                                      coarse_repulsion_factor = 0.3,
                                      coarse_repulsion_sample = 8,
                                      coarse_repulsion_exact_below = 32,
                                      final_mode = "kk_repulse",
                                      seed = 31)
  expect_gt(max(abs(coords_fr - coords_kkr)), 1e-6)
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
  expect_error(
    grip.layout.globalrep(edges, n = 10, dim = 2,
                          final_mode = "banana",
                          seed = 1),
    "'arg' should be one of"
  )
})

test_that("globalrep adaptive default final_rounds schedule is stable", {
  expect_equal(grip:::grip.globalrep.default.final_rounds(1000L), 384L)
  expect_equal(grip:::grip.globalrep.default.final_rounds(1001L), 320L)
  expect_equal(grip:::grip.globalrep.default.final_rounds(5001L), 256L)
  expect_equal(grip:::grip.globalrep.default.final_rounds(20001L), 200L)
  expect_equal(grip:::grip.globalrep.default.final_rounds(50001L), 128L)
})

test_that("globalrep small-graph defaults match the new fixed candidate profile", {
  edges <- edges.path(12)
  coords_default <- grip.layout.globalrep(edges, n = 12, dim = 2, seed = 7)
  coords_explicit <- grip.layout.globalrep(
    edges, n = 12, dim = 2,
    placement = "barycenter",
    rounds = 160, final_rounds = 384,
    num_init = 24, num_nbrs = 20,
    r = 0.03, s = 7.5,
    repulsion_factor = 2.5,
    coarse_repulsion_factor = 1.5,
    coarse_repulsion_sample = 16,
    coarse_repulsion_exact_below = 64,
    seed = 7
  )
  expect_identical(coords_default, coords_explicit)
})

test_that("globalrep larger-graph defaults taper final_rounds only", {
  n <- 1001L
  edges <- edges.path(n)
  coords_default <- grip.layout.globalrep(edges, n = n, dim = 2, seed = 9)
  coords_explicit <- grip.layout.globalrep(
    edges, n = n, dim = 2,
    placement = "barycenter",
    rounds = 160, final_rounds = 320,
    num_init = 24, num_nbrs = 20,
    r = 0.03, s = 7.5,
    repulsion_factor = 2.5,
    coarse_repulsion_factor = 1.5,
    coarse_repulsion_sample = 16,
    coarse_repulsion_exact_below = 64,
    seed = 9
  )
  expect_identical(coords_default, coords_explicit)
})

test_that("globalrep disconnected handling matches grip.layout.legacy when disabled", {
  edges <- rbind(
    cbind(1:2, 2:3),
    cbind(5:6, 6:7)
  )
  expect_warning(
    coords_base <- grip.layout.legacy(edges, n = 7, dim = 2, seed = 11),
    "laying out components separately"
  )
  expect_warning(
    coords_globalrep <- grip.layout.globalrep(edges, n = 7, dim = 2,
                                              rounds = 20,
                                              final_rounds = 25,
                                              num_init = 36,
                                              num_nbrs = 10,
                                              r = 0.15,
                                              s = 3.0,
                                              repulsion_factor = 1.0,
                                              coarse_repulsion_factor = 0,
                                              coarse_repulsion_sample = 8,
                                              coarse_repulsion_exact_below = 32,
                                              seed = 11),
    "laying out components separately"
  )
  expect_identical(coords_base, coords_globalrep)
})
