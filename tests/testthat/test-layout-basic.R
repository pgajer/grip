test_that("basic layout returns finite matrix", {
  edges <- edges.path(10)
  coords <- grip(edges, n = 10, dim = 2,
                        placement = "barycenter",
                        rounds = 5, final_rounds = 3,
                        num_init = 5, num_nbrs = 6,
                        seed = 123)
  expect_equal(dim(coords), c(10, 2))
  expect_true(all(is.finite(coords)))
})

test_that("vignette path example returns finite coordinates", {
  edges <- edges.path(12)
  coords <- grip(edges, n = 12, dim = 2,
                        placement = "barycenter",
                        rounds = 25, final_rounds = 25,
                        num_init = 5, num_nbrs = 6,
                        seed = 1)
  expect_equal(dim(coords), c(12, 2))
  expect_true(all(is.finite(coords)))
})

test_that("seeded runs are deterministic", {
  edges <- edges.cycle(12)
  coords1 <- grip(edges, n = 12, dim = 2,
                         placement = "barycenter",
                         rounds = 4, final_rounds = 2,
                         num_init = 4, num_nbrs = 5,
                         seed = 42)
  coords2 <- grip(edges, n = 12, dim = 2,
                         placement = "barycenter",
                         rounds = 4, final_rounds = 2,
                         num_init = 4, num_nbrs = 5,
                         seed = 42)
  expect_identical(coords1, coords2)
})

test_that("num_nbrs changes the layout with a fixed seed", {
  edges <- edges.mesh(5, 5)
  coords_small <- grip(edges, n = 25, dim = 2,
                              placement = "barycenter",
                              rounds = 8, final_rounds = 8,
                              num_init = 6, num_nbrs = 2,
                              seed = 19)
  coords_large <- grip(edges, n = 25, dim = 2,
                              placement = "barycenter",
                              rounds = 8, final_rounds = 8,
                              num_init = 6, num_nbrs = 10,
                              seed = 19)
  expect_gt(max(abs(coords_small - coords_large)), 1e-6)
})

test_that("r and s change the layout with a fixed seed", {
  edges <- edges.mesh(5, 5)
  coords_cool <- grip(edges, n = 25, dim = 2,
                             placement = "barycenter",
                             rounds = 8, final_rounds = 8,
                             num_init = 6, num_nbrs = 8,
                             r = 0.00, s = 0.00,
                             seed = 23)
  coords_adaptive <- grip(edges, n = 25, dim = 2,
                                 placement = "barycenter",
                                 rounds = 8, final_rounds = 8,
                                 num_init = 6, num_nbrs = 8,
                                 r = 0.30, s = 6.00,
                                 seed = 23)
  expect_gt(max(abs(coords_cool - coords_adaptive)), 1e-6)
})

test_that("repulsion_factor changes the layout with a fixed seed", {
  edges <- edges.mesh(5, 5)
  coords_none <- grip(edges, n = 25, dim = 2,
                             placement = "barycenter",
                             rounds = 8, final_rounds = 8,
                             num_init = 6, num_nbrs = 8,
                             repulsion_factor = 0,
                             seed = 29)
  coords_more <- grip(edges, n = 25, dim = 2,
                             placement = "barycenter",
                             rounds = 8, final_rounds = 8,
                             num_init = 6, num_nbrs = 8,
                             repulsion_factor = 2,
                             seed = 29)
  expect_gt(max(abs(coords_none - coords_more)), 1e-6)
})

test_that("carpet preset matches the explicit carpet tuning profile", {
  edges <- edges.sierpinski.carpet(2)
  n <- max(edges)
  coords_preset <- grip(edges, n = n, dim = 2,
                               preset = "carpet",
                               seed = 41)
  coords_explicit <- grip(edges, n = n, dim = 2,
                                 placement = "barycenter",
                                 rounds = 160, final_rounds = 288,
                                 num_init = 28, num_nbrs = 24,
                                 r = 0.03, s = 6.0,
                                 repulsion_factor = 2.5,
                                 coarse_repulsion_factor = 1.5,
                                 coarse_repulsion_sample = 16,
                                 coarse_repulsion_exact_below = 64,
                                 seed = 41)
  expect_identical(coords_preset, coords_explicit)
})

test_that("mesh preset matches the explicit mesh tuning profile", {
  edges <- edges.mesh(4, 4)
  n <- max(edges)
  coords_preset <- grip(edges, n = n, dim = 2,
                               preset = "mesh",
                               seed = 45)
  coords_explicit <- grip(edges, n = n, dim = 2,
                                 placement = "barycenter",
                                 rounds = 128, final_rounds = 128,
                                 num_init = 12, num_nbrs = 20,
                                 r = 0.10, s = 4.5,
                                 repulsion_factor = 1.5,
                                 coarse_repulsion_factor = 1.5,
                                 coarse_repulsion_sample = 16,
                                 coarse_repulsion_exact_below = 64,
                                 seed = 45)
  expect_identical(coords_preset, coords_explicit)
})

test_that("torus preset matches the explicit torus tuning profile", {
  edges <- edges.torus(4, 4)
  n <- max(edges)
  coords_preset <- grip(edges, n = n, dim = 3,
                               preset = "torus",
                               seed = 47)
  coords_explicit <- grip(edges, n = n, dim = 3,
                                 placement = "barycenter",
                                 rounds = 192, final_rounds = 288,
                                 num_init = 12, num_nbrs = 28,
                                 r = 0.05, s = 7.5,
                                 repulsion_factor = 0.75,
                                 coarse_repulsion_factor = 1.5,
                                 coarse_repulsion_sample = 16,
                                 coarse_repulsion_exact_below = 64,
                                 seed = 47)
  expect_identical(coords_preset, coords_explicit)
})

test_that("tree preset matches the explicit tree tuning profile", {
  edges <- edges.kary.tree(k = 2, depth = 4)
  n <- max(edges)
  coords_preset <- grip(edges, n = n, dim = 2,
                               preset = "tree",
                               seed = 61)
  coords_explicit <- grip(edges, n = n, dim = 2,
                                 placement = "circle",
                                 rounds = 64, final_rounds = 160,
                                 num_init = 28, num_nbrs = 8,
                                 r = 0.05, s = 7.5,
                                 repulsion_factor = 0,
                                 coarse_repulsion_factor = 1.5,
                                 coarse_repulsion_sample = 16,
                                 coarse_repulsion_exact_below = 64,
                                 seed = 61)
  expect_identical(coords_preset, coords_explicit)
})

test_that("tree preset uses barycenter placement in 3D without warning", {
  edges <- edges.kary.tree(k = 2, depth = 4)
  n <- max(edges)
  expect_no_warning({
    coords_preset <- grip(edges, n = n, dim = 3,
                                 preset = "tree",
                                 seed = 63)
    coords_explicit <- grip(edges, n = n, dim = 3,
                                   placement = "barycenter",
                                   rounds = 64, final_rounds = 160,
                                   num_init = 28, num_nbrs = 8,
                                   r = 0.05, s = 7.5,
                                   repulsion_factor = 0,
                                   coarse_repulsion_factor = 1.5,
                                   coarse_repulsion_sample = 16,
                                   coarse_repulsion_exact_below = 64,
                                   seed = 63)
    expect_identical(coords_preset, coords_explicit)
  })
})

test_that("explicit tuning args override the carpet preset", {
  edges <- edges.sierpinski.carpet(2)
  n <- max(edges)
  coords_preset <- grip(edges, n = n, dim = 2,
                               preset = "carpet",
                               repulsion_factor = 1.75,
                               seed = 43)
  coords_explicit <- grip(edges, n = n, dim = 2,
                                 placement = "barycenter",
                                 rounds = 160, final_rounds = 288,
                                 num_init = 28, num_nbrs = 24,
                                 r = 0.03, s = 6.0,
                                 repulsion_factor = 1.75,
                                 coarse_repulsion_factor = 1.5,
                                 coarse_repulsion_sample = 16,
                                 coarse_repulsion_exact_below = 64,
                                 seed = 43)
  expect_identical(coords_preset, coords_explicit)
})

test_that("explicit tuning args override the mesh preset", {
  edges <- edges.mesh(4, 4)
  n <- max(edges)
  coords_preset <- grip(edges, n = n, dim = 2,
                               preset = "mesh",
                               repulsion_factor = 0.75,
                               seed = 49)
  coords_explicit <- grip(edges, n = n, dim = 2,
                                 placement = "barycenter",
                                 rounds = 128, final_rounds = 128,
                                 num_init = 12, num_nbrs = 20,
                                 r = 0.10, s = 4.5,
                                 repulsion_factor = 0.75,
                                 coarse_repulsion_factor = 1.5,
                                 coarse_repulsion_sample = 16,
                                 coarse_repulsion_exact_below = 64,
                                 seed = 49)
  expect_identical(coords_preset, coords_explicit)
})

test_that("explicit tuning args override the torus preset", {
  edges <- edges.torus(4, 4)
  n <- max(edges)
  coords_preset <- grip(edges, n = n, dim = 3,
                               preset = "torus",
                               final_rounds = 320,
                               seed = 53)
  coords_explicit <- grip(edges, n = n, dim = 3,
                                 placement = "barycenter",
                                 rounds = 192, final_rounds = 320,
                                 num_init = 12, num_nbrs = 28,
                                 r = 0.05, s = 7.5,
                                 repulsion_factor = 0.75,
                                 coarse_repulsion_factor = 1.5,
                                 coarse_repulsion_sample = 16,
                                 coarse_repulsion_exact_below = 64,
                                 seed = 53)
  expect_identical(coords_preset, coords_explicit)
})

test_that("explicit tuning args override the tree preset", {
  edges <- edges.kary.tree(k = 2, depth = 4)
  n <- max(edges)
  coords_preset <- grip(edges, n = n, dim = 2,
                               preset = "tree",
                               repulsion_factor = 0.5,
                               seed = 67)
  coords_explicit <- grip(edges, n = n, dim = 2,
                                 placement = "circle",
                                 rounds = 64, final_rounds = 160,
                                 num_init = 28, num_nbrs = 8,
                                 r = 0.05, s = 7.5,
                                 repulsion_factor = 0.5,
                                 coarse_repulsion_factor = 1.5,
                                 coarse_repulsion_sample = 16,
                                 coarse_repulsion_exact_below = 64,
                                 seed = 67)
  expect_identical(coords_preset, coords_explicit)
})

test_that("invalid tuning parameters are rejected", {
  edges <- edges.cycle(10)
  expect_error(
    grip(edges, n = 10, dim = 2, num_nbrs = 0, seed = 1),
    "num_nbrs must be a positive integer"
  )
  expect_error(
    grip(edges, n = 10, dim = 2, r = -0.1, seed = 1),
    "r must be in \\[0, 1\\]"
  )
  expect_error(
    grip(edges, n = 10, dim = 2, s = -1, seed = 1),
    "s must be >= 0"
  )
  expect_error(
    grip(edges, n = 10, dim = 2, repulsion_factor = -0.1, seed = 1),
    "repulsion_factor must be >= 0"
  )
  expect_error(
    grip(edges, n = 10, dim = 2, preset = "bogus", seed = 1),
    "preset must be NULL, 'carpet', 'mesh', 'torus', or 'tree'"
  )
})

test_that("circle placement works in 2D", {
  edges <- edges.path(8)
  coords <- grip(edges, n = 8, dim = 2,
                        placement = "circle",
                        rounds = 4, final_rounds = 2,
                        num_init = 4, num_nbrs = 5,
                        seed = 9)
  expect_equal(dim(coords), c(8, 2))
  expect_true(all(is.finite(coords)))
})

test_that("circle placement falls back in 3D with warning", {
  edges <- edges.mesh(4, 4)
  expect_warning({
    coords <- grip(edges, n = 16, dim = 3,
                          placement = "circle",
                          rounds = 4, final_rounds = 2,
                          num_init = 5, num_nbrs = 6,
                          seed = 11)
    expect_equal(dim(coords), c(16, 3))
    expect_true(all(is.finite(coords)))
  }, "circle placement is only used for 2D")
})

test_that("tree example runs", {
  edges <- edges.kary.tree(k = 2, depth = 2)
  n <- max(edges)
  coords <- grip(edges, n = n, dim = 2,
                        placement = "barycenter",
                        rounds = 4, final_rounds = 2,
                        num_init = 4, num_nbrs = 5,
                        seed = 5)
  expect_equal(dim(coords), c(n, 2))
  expect_true(all(is.finite(coords)))
})

test_that("adj_list + weight_list input works", {
  adj_list <- list(c(2), c(1, 3), c(2, 4), c(3))
  weight_list <- list(c(1.0), c(1.0, 2.0), c(2.0, 1.5), c(1.5))
  coords <- grip(adj_list = adj_list,
                        weight_list = weight_list,
                        n = 4,
                        dim = 2,
                        placement = "barycenter",
                        rounds = 4, final_rounds = 2,
                        num_init = 3, num_nbrs = 3,
                        seed = 12)
  expect_equal(dim(coords), c(4, 2))
  expect_true(all(is.finite(coords)))
})

test_that("adj_list input works without weights", {
  adj_list <- list(c(2), c(1, 3), c(2, 4), c(3))
  coords <- grip(adj_list = adj_list,
                        n = 4,
                        dim = 2,
                        placement = "barycenter",
                        rounds = 4, final_rounds = 2,
                        num_init = 3, num_nbrs = 3,
                        seed = 13)
  expect_equal(dim(coords), c(4, 2))
  expect_true(all(is.finite(coords)))
})

test_that("disconnected graph defaults to safe component layouts", {
  edges <- rbind(
    cbind(1:4, 2:5),
    cbind(7:9, 8:10)
  )
  expect_warning({
    coords <- grip(edges = edges,
                          n = 10,
                          dim = 3,
                          placement = "barycenter",
                          rounds = 4,
                          final_rounds = 2,
                          num_init = 4,
                          num_nbrs = 5,
                          seed = 21)
    expect_equal(dim(coords), c(10, 3))
    expect_true(all(is.finite(coords)))
  }, "connected components")
})

test_that("disconnected mode can be set to error", {
  edges <- rbind(
    cbind(1:2, 2:3),
    cbind(5:6, 6:7)
  )
  expect_error(
    grip(edges = edges,
                n = 7,
                dim = 2,
                disconnected = "error",
                seed = 7),
    "connected components"
  )
})

test_that("disconnected adj_list + weights is handled safely", {
  adj_list <- list(
    c(2), c(1),           # component 1
    c(4), c(3),           # component 2
    integer(0)            # isolated
  )
  weight_list <- list(
    c(1.0), c(1.0),
    c(2.0), c(2.0),
    numeric(0)
  )
  expect_warning({
    coords <- grip(adj_list = adj_list,
                          weight_list = weight_list,
                          n = 5,
                          dim = 2,
                          placement = "barycenter",
                          rounds = 4,
                          final_rounds = 2,
                          num_init = 3,
                          num_nbrs = 3,
                          seed = 31)
    expect_equal(dim(coords), c(5, 2))
    expect_true(all(is.finite(coords)))
  }, "connected components")
})
