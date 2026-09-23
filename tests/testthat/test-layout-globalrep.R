test_that("globalrep layout returns a finite matrix", {
  edges <- edges.path(10)
  coords <- globalrep.grip(edges, n = 10, dim = 2,
                                  placement = "barycenter",
                                  rounds = 5, final.rounds = 3,
                                  num.init = 5, num.nbrs = 6,
                                  coarse.repulsion.factor = 0.2,
                                  coarse.repulsion.sample = 8,
                                  coarse.repulsion.exact.below = 32,
                                  seed = 123)
  expect_equal(dim(coords), c(10, 2))
  expect_true(all(is.finite(coords)))
})

test_that("globalrep seeded runs are deterministic", {
  edges <- edges.cycle(12)
  coords1 <- globalrep.grip(edges, n = 12, dim = 2,
                                   placement = "barycenter",
                                   rounds = 4, final.rounds = 2,
                                   num.init = 4, num.nbrs = 5,
                                   coarse.repulsion.factor = 0.2,
                                   coarse.repulsion.sample = 8,
                                   coarse.repulsion.exact.below = 32,
                                   seed = 42)
  coords2 <- globalrep.grip(edges, n = 12, dim = 2,
                                   placement = "barycenter",
                                   rounds = 4, final.rounds = 2,
                                   num.init = 4, num.nbrs = 5,
                                   coarse.repulsion.factor = 0.2,
                                   coarse.repulsion.sample = 8,
                                   coarse.repulsion.exact.below = 32,
                                   seed = 42)
  expect_identical(coords1, coords2)
})

test_that("globalrep with zero coarse repulsion matches legacy.grip", {
  edges <- edges.mesh(5, 5)
  coords_base <- legacy.grip(edges, n = 25, dim = 2,
                                    placement = "barycenter",
                                    rounds = 8, final.rounds = 8,
                                    num.init = 6, num.nbrs = 8,
                                    r = 0.15, s = 3.0,
                                    repulsion.factor = 1.5,
                                    seed = 29)
  coords_globalrep <- globalrep.grip(edges, n = 25, dim = 2,
                                            placement = "barycenter",
                                            rounds = 8, final.rounds = 8,
                                            num.init = 6, num.nbrs = 8,
                                            r = 0.15, s = 3.0,
                                            repulsion.factor = 1.5,
                                            coarse.repulsion.factor = 0,
                                            coarse.repulsion.sample = 8,
                                            coarse.repulsion.exact.below = 32,
                                            seed = 29)
  expect_identical(coords_base, coords_globalrep)
})

test_that("grip is an alias of globalrep.grip", {
  edges <- edges.mesh(5, 5)
  coords_primary <- grip(edges, n = 25, dim = 2, seed = 17)
  coords_alias <- globalrep.grip(edges, n = 25, dim = 2, seed = 17)
  expect_identical(coords_primary, coords_alias)
})

test_that("globalrep coarse repulsion changes the layout with a fixed seed", {
  edges <- edges.mesh(5, 5)
  coords_none <- globalrep.grip(edges, n = 25, dim = 2,
                                       placement = "barycenter",
                                       rounds = 8, final.rounds = 8,
                                       num.init = 6, num.nbrs = 8,
                                       coarse.repulsion.factor = 0,
                                       coarse.repulsion.sample = 8,
                                       coarse.repulsion.exact.below = 32,
                                       seed = 31)
  coords_more <- globalrep.grip(edges, n = 25, dim = 2,
                                       placement = "barycenter",
                                       rounds = 8, final.rounds = 8,
                                       num.init = 6, num.nbrs = 8,
                                       coarse.repulsion.factor = 0.3,
                                       coarse.repulsion.sample = 8,
                                       coarse.repulsion.exact.below = 32,
                                       seed = 31)
  expect_gt(max(abs(coords_none - coords_more)), 1e-6)
})

test_that("globalrep final_mode can switch the final stage", {
  edges <- edges.mesh(5, 5)
  coords_fr <- globalrep.grip(edges, n = 25, dim = 2,
                                     placement = "barycenter",
                                     rounds = 8, final.rounds = 8,
                                     num.init = 6, num.nbrs = 8,
                                     coarse.repulsion.factor = 0.3,
                                     coarse.repulsion.sample = 8,
                                     coarse.repulsion.exact.below = 32,
                                     final.mode = "fr",
                                     seed = 31)
  coords_kkr <- globalrep.grip(edges, n = 25, dim = 2,
                                      placement = "barycenter",
                                      rounds = 8, final.rounds = 8,
                                      num.init = 6, num.nbrs = 8,
                                      coarse.repulsion.factor = 0.3,
                                      coarse.repulsion.sample = 8,
                                      coarse.repulsion.exact.below = 32,
                                      final.mode = "kk_repulse",
                                      seed = 31)
  expect_gt(max(abs(coords_fr - coords_kkr)), 1e-6)
})

test_that("globalrep structural final-stage knobs can change the layout", {
  edges <- edges.mesh(5, 5)
  coords_base <- globalrep.grip(edges, n = 25, dim = 2,
                                       placement = "barycenter",
                                       rounds = 8, final.rounds = 8,
                                       num.init = 6, num.nbrs = 8,
                                       coarse.repulsion.factor = 0.3,
                                       coarse.repulsion.sample = 8,
                                       coarse.repulsion.exact.below = 32,
                                       final.mode = "fr",
                                       seed = 31)
  coords_struct <- globalrep.grip(edges, n = 25, dim = 2,
                                         placement = "barycenter",
                                         rounds = 8, final.rounds = 8,
                                         num.init = 6, num.nbrs = 8,
                                         coarse.repulsion.factor = 0.3,
                                         coarse.repulsion.sample = 8,
                                         coarse.repulsion.exact.below = 32,
                                         final.anchor.factor = 1,
                                         final.move.scale.after.first = 0.5,
                                         final.mode = "fr",
                                         seed = 31)
  expect_gt(max(abs(coords_base - coords_struct)), 1e-6)
})

test_that("globalrep level-0 insertion knobs are wired and can change the layout", {
  edges <- edges.sierpinski.carpet(2)
  n <- max(edges)

  coords_base <- globalrep.grip(edges, n = n, dim = 2,
                                       placement = "barycenter",
                                       rounds = 8, final.rounds = 8,
                                       num.init = 6, num.nbrs = 8,
                                       coarse.repulsion.factor = 0.3,
                                       coarse.repulsion.sample = 8,
                                       coarse.repulsion.exact.below = 32,
                                       level0.insertion.mode = "inherit",
                                       seed = 41)
  coords_bary <- globalrep.grip(edges, n = n, dim = 2,
                                       placement = "barycenter",
                                       rounds = 8, final.rounds = 8,
                                       num.init = 6, num.nbrs = 8,
                                       coarse.repulsion.factor = 0.3,
                                       coarse.repulsion.sample = 8,
                                       coarse.repulsion.exact.below = 32,
                                       level0.insertion.mode = "barycenter",
                                       seed = 41)
  coords_no_local_kk <- globalrep.grip(edges, n = n, dim = 2,
                                              placement = "barycenter",
                                              rounds = 8, final.rounds = 8,
                                              num.init = 6, num.nbrs = 8,
                                              coarse.repulsion.factor = 0.3,
                                              coarse.repulsion.sample = 8,
                                              coarse.repulsion.exact.below = 32,
                                              level0.local.kk.steps = 0,
                                              seed = 41)
  coords_ls <- globalrep.grip(edges, n = n, dim = 2,
                                     placement = "barycenter",
                                     rounds = 8, final.rounds = 8,
                                     num.init = 6, num.nbrs = 8,
                                     coarse.repulsion.factor = 0.3,
                                     coarse.repulsion.sample = 8,
                                     coarse.repulsion.exact.below = 32,
                                     level0.insertion.mode = "least_squares",
                                     level0.anchor.count = 6,
                                     seed = 41)

  expect_identical(coords_base, coords_bary)
  expect_gt(max(abs(coords_bary - coords_no_local_kk)), 1e-6)
  expect_gt(max(abs(coords_no_local_kk - coords_ls)), 1e-6)
})

test_that("global insertion anchor knobs can change the layout", {
  edges <- edges.sierpinski.carpet(2)
  n <- max(edges)

  coords_base <- globalrep.grip(edges, n = n, dim = 2,
                                       placement = "barycenter",
                                       rounds = 8, final.rounds = 8,
                                       num.init = 6, num.nbrs = 8,
                                       coarse.repulsion.factor = 0.3,
                                       coarse.repulsion.sample = 8,
                                       coarse.repulsion.exact.below = 32,
                                       insertion.anchor.count = 3,
                                       insertion.anchor.scope = "any_higher",
                                       seed = 53)
  coords_more <- globalrep.grip(edges, n = n, dim = 2,
                                       placement = "barycenter",
                                       rounds = 8, final.rounds = 8,
                                       num.init = 6, num.nbrs = 8,
                                       coarse.repulsion.factor = 0.3,
                                       coarse.repulsion.sample = 8,
                                       coarse.repulsion.exact.below = 32,
                                       insertion.anchor.count = 6,
                                       insertion.anchor.scope = "any_higher",
                                       seed = 53)
  coords_prev <- globalrep.grip(edges, n = n, dim = 2,
                                       placement = "barycenter",
                                       rounds = 8, final.rounds = 8,
                                       num.init = 6, num.nbrs = 8,
                                       coarse.repulsion.factor = 0.3,
                                       coarse.repulsion.sample = 8,
                                       coarse.repulsion.exact.below = 32,
                                       insertion.anchor.count = 6,
                                       insertion.anchor.scope = "prev_misf",
                                       seed = 53)

  expect_gt(max(abs(coords_base - coords_more)), 1e-6)
  expect_gt(max(abs(coords_more - coords_prev)), 1e-6)
})

test_that("insertion anchor strategy can change the layout", {
  edges <- edges.sierpinski.carpet(2)
  n <- max(edges)

  coords_first <- globalrep.grip(
    edges, n = n, dim = 2,
    rounds = 8, final.rounds = 8,
    num.init = 6, num.nbrs = 8,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    insertion.anchor.count = 6,
    insertion.anchor.scope = "prev_misf",
    insertion.anchor.strategy = "first",
    seed = 71
  )
  coords_band <- globalrep.grip(
    edges, n = n, dim = 2,
    rounds = 8, final.rounds = 8,
    num.init = 6, num.nbrs = 8,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    insertion.anchor.count = 6,
    insertion.anchor.scope = "prev_misf",
    insertion.anchor.strategy = "distance_band",
    seed = 71
  )
  coords_balanced <- globalrep.grip(
    edges, n = n, dim = 2,
    rounds = 8, final.rounds = 8,
    num.init = 6, num.nbrs = 8,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    insertion.anchor.count = 6,
    insertion.anchor.scope = "prev_misf",
    insertion.anchor.strategy = "balanced_band",
    seed = 71
  )
  coords_spread_prev <- globalrep.grip(
    edges, n = n, dim = 2,
    rounds = 8, final.rounds = 8,
    num.init = 6, num.nbrs = 8,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    insertion.anchor.count = 6,
    insertion.anchor.scope = "prev_misf",
    insertion.anchor.strategy = "spread_prev",
    level0.insertion.mode = "least_squares",
    level0.local.kk.steps = 1,
    seed = 71
  )

  expect_gt(max(abs(coords_first - coords_band)), 1e-6)
  expect_gt(max(abs(coords_band - coords_balanced)), 1e-6)
  expect_gt(max(abs(coords_balanced - coords_spread_prev)), 1e-6)
})

test_that("globalrep LGKK polish knobs can change the layout", {
  edges <- edges.mesh(5, 5)
  coords_base <- globalrep.grip(
    edges, n = 25, dim = 2,
    rounds = 8, final.rounds = 8,
    num.init = 6, num.nbrs = 8,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    lgkk.polish.rounds = 0,
    seed = 61
  )
  coords_polish <- globalrep.grip(
    edges, n = 25, dim = 2,
    rounds = 8, final.rounds = 8,
    num.init = 6, num.nbrs = 8,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    lgkk.polish.rounds = 4,
    lgkk.local.nbrs = 6,
    lgkk.landmark.count = 4,
    seed = 61
  )
  expect_gt(max(abs(coords_base - coords_polish)), 1e-6)
})

test_that("globalrep multiscale LGKK knobs can change the layout", {
  edges <- edges.mesh(5, 5)
  coords_base <- globalrep.grip(
    edges, n = 25, dim = 2,
    rounds = 8, final.rounds = 8,
    num.init = 6, num.nbrs = 8,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    lgkk.multiscale.rounds = 0,
    seed = 67
  )
  coords_lgkk <- globalrep.grip(
    edges, n = 25, dim = 2,
    rounds = 8, final.rounds = 8,
    num.init = 6, num.nbrs = 8,
    coarse.repulsion.factor = 0.3,
    coarse.repulsion.sample = 8,
    coarse.repulsion.exact.below = 32,
    lgkk.multiscale.rounds = 2,
    lgkk.local.nbrs = 6,
    lgkk.landmark.count = 4,
    lgkk.multiscale.scope = "all",
    lgkk.active.limit = 512,
    seed = 67
  )
  expect_gt(max(abs(coords_base - coords_lgkk)), 1e-6)
})

test_that("globalrep staged LGKK budgets can change weighted layouts", {
  h <- 6L
  w <- 6L
  edges <- edges.mesh(h, w)
  row_of <- function(v) ((v - 1L) %/% w) + 1L
  col_of <- function(v) ((v - 1L) %% w) + 1L
  edge_weights <- apply(edges, 1L, function(e) {
    r1 <- row_of(e[[1L]])
    r2 <- row_of(e[[2L]])
    c1 <- col_of(e[[1L]])
    c2 <- col_of(e[[2L]])
    if (r1 == r2 && abs(c1 - c2) == 1L) 1 else 2
  })

  coords_shared <- globalrep.grip(
    edges, n = h * w, dim = 2,
    edge.weights = edge_weights,
    lgkk.multiscale.rounds = 4,
    lgkk.local.nbrs = 6,
    lgkk.landmark.count = 8,
    lgkk.multiscale.scope = "all",
    lgkk.active.limit = 4096,
    seed = 1
  )
  coords_a3 <- globalrep.grip(
    edges, n = h * w, dim = 2,
    edge.weights = edge_weights,
    lgkk.multiscale.rounds = 0,
    lgkk.rounds.coarse = 1,
    lgkk.rounds.pre.final = 2,
    lgkk.rounds.final = 4,
    lgkk.local.nbrs = 6,
    lgkk.landmark.count = 8,
    lgkk.multiscale.scope = "all",
    lgkk.active.limit = 4096,
    seed = 1
  )

  expect_gt(max(abs(coords_shared - coords_a3)), 1e-6)
})

test_that("globalrep validates the new tuning parameters", {
  edges <- edges.cycle(10)
  expect_error(
    globalrep.grip(edges, n = 10, dim = 2,
                          coarse.repulsion.factor = -0.1,
                          seed = 1),
    "coarse.repulsion.factor must be >= 0"
  )
  expect_error(
    globalrep.grip(edges, n = 10, dim = 2,
                          coarse.repulsion.sample = 0,
                          seed = 1),
    "coarse.repulsion.sample must be a positive integer"
  )
  expect_error(
    globalrep.grip(edges, n = 10, dim = 2,
                          coarse.repulsion.exact.below = 0,
                          seed = 1),
    "coarse.repulsion.exact.below must be a positive integer"
  )
  expect_error(
    globalrep.grip(edges, n = 10, dim = 2,
                          final.anchor.factor = -0.1,
                          seed = 1),
    "final.anchor.factor must be >= 0"
  )
  expect_error(
    globalrep.grip(edges, n = 10, dim = 2,
                          final.move.scale.after.first = 1.1,
                          seed = 1),
    "final.move.scale.after.first must be in \\[0, 1\\]"
  )
  expect_error(
    globalrep.grip(edges, n = 10, dim = 2,
                          insertion.anchor.count = 0,
                          seed = 1),
    "insertion.anchor.count must be a positive integer"
  )
  expect_error(
    globalrep.grip(edges, n = 10, dim = 2,
                          insertion.anchor.scope = "banana",
                          seed = 1),
    "'arg' should be one of"
  )
  expect_error(
    globalrep.grip(edges, n = 10, dim = 2,
                          final.mode = "banana",
                          seed = 1),
    "'arg' should be one of"
  )
  expect_error(
    globalrep.grip(edges, n = 10, dim = 2,
                          level0.insertion.mode = "banana",
                          seed = 1),
    "'arg' should be one of"
  )
  expect_error(
    globalrep.grip(edges, n = 10, dim = 2,
                          level0.anchor.count = 0,
                          seed = 1),
    "level0.anchor.count must be a positive integer"
  )
  expect_error(
    globalrep.grip(edges, n = 10, dim = 2,
                          level0.local.kk.steps = -1,
                          seed = 1),
    "level0.local.kk.steps must be a non-negative integer"
  )
  expect_error(
    globalrep.grip(edges, n = 10, dim = 2,
                          lgkk.polish.rounds = -1,
                          seed = 1),
    "lgkk.polish.rounds must be a non-negative integer"
  )
  expect_error(
    globalrep.grip(edges, n = 10, dim = 2,
                          lgkk.local.nbrs = -1,
                          seed = 1),
    "lgkk.local.nbrs must be a non-negative integer"
  )
  expect_error(
    globalrep.grip(edges, n = 10, dim = 2,
                          lgkk.landmark.count = -1,
                          seed = 1),
    "lgkk.landmark.count must be a non-negative integer"
  )
  expect_error(
    globalrep.grip(edges, n = 10, dim = 2,
                          lgkk.rounds.coarse = -1,
                          seed = 1),
    "lgkk.rounds.coarse must be a non-negative integer"
  )
  expect_error(
    globalrep.grip(edges, n = 10, dim = 2,
                          lgkk.rounds.pre.final = -1,
                          seed = 1),
    "lgkk.rounds.pre.final must be a non-negative integer"
  )
  expect_error(
    globalrep.grip(edges, n = 10, dim = 2,
                          lgkk.rounds.final = -1,
                          seed = 1),
    "lgkk.rounds.final must be a non-negative integer"
  )
})

test_that("globalrep adaptive default final_rounds schedule is stable", {
  expect_equal(grip:::grip.globalrep.default.final.rounds(1000L), 384L)
  expect_equal(grip:::grip.globalrep.default.final.rounds(1001L), 320L)
  expect_equal(grip:::grip.globalrep.default.final.rounds(5001L), 256L)
  expect_equal(grip:::grip.globalrep.default.final.rounds(20001L), 200L)
  expect_equal(grip:::grip.globalrep.default.final.rounds(50001L), 128L)
})

test_that("globalrep small-graph defaults match the new fixed candidate profile", {
  edges <- edges.path(12)
  coords_default <- globalrep.grip(edges, n = 12, dim = 2, seed = 7)
  coords_explicit <- globalrep.grip(
    edges, n = 12, dim = 2,
    placement = "barycenter",
    rounds = 160, final.rounds = 384,
    num.init = 24, num.nbrs = 20,
    r = 0.03, s = 7.5,
    repulsion.factor = 2.5,
    coarse.repulsion.factor = 1.5,
    coarse.repulsion.sample = 16,
    coarse.repulsion.exact.below = 64,
    insertion.anchor.count = 3,
    insertion.anchor.scope = "any_higher",
    level0.insertion.mode = "inherit",
    level0.anchor.count = 3,
    level0.local.kk.steps = 3,
    lgkk.polish.rounds = 0,
    lgkk.local.nbrs = 20,
    lgkk.landmark.count = 8,
    seed = 7
  )
  expect_identical(coords_default, coords_explicit)
})

test_that("globalrep larger-graph defaults taper final_rounds only", {
  n <- 1001L
  edges <- edges.path(n)
  coords_default <- globalrep.grip(edges, n = n, dim = 2, seed = 9)
  coords_explicit <- globalrep.grip(
    edges, n = n, dim = 2,
    placement = "barycenter",
    rounds = 160, final.rounds = 320,
    num.init = 24, num.nbrs = 20,
    r = 0.03, s = 7.5,
    repulsion.factor = 2.5,
    coarse.repulsion.factor = 1.5,
    coarse.repulsion.sample = 16,
    coarse.repulsion.exact.below = 64,
    insertion.anchor.count = 3,
    insertion.anchor.scope = "any_higher",
    level0.insertion.mode = "inherit",
    level0.anchor.count = 3,
    level0.local.kk.steps = 3,
    lgkk.polish.rounds = 0,
    lgkk.local.nbrs = 20,
    lgkk.landmark.count = 8,
    seed = 9
  )
  expect_identical(coords_default, coords_explicit)
})

test_that("globalrep disconnected handling matches legacy.grip when disabled", {
  edges <- rbind(
    cbind(1:2, 2:3),
    cbind(5:6, 6:7)
  )
  expect_warning(
    coords_base <- legacy.grip(edges, n = 7, dim = 2, seed = 11),
    "laying out components separately"
  )
  expect_warning(
    coords_globalrep <- globalrep.grip(edges, n = 7, dim = 2,
                                              rounds = 20,
                                              final.rounds = 25,
                                              num.init = 36,
                                              num.nbrs = 10,
                                              r = 0.15,
                                              s = 3.0,
                                              repulsion.factor = 1.0,
                                              coarse.repulsion.factor = 0,
                                              coarse.repulsion.sample = 8,
                                              coarse.repulsion.exact.below = 32,
                                              seed = 11),
    "laying out components separately"
  )
  expect_identical(coords_base, coords_globalrep)
})
