test_that("round trace returns frames, metadata, and inactive NA rows", {
  edges <- edges.path(8)
  tr <- grip.layout.trace(edges = edges,
                          n = 8,
                          dim = 2,
                          placement = "barycenter",
                          rounds = 4,
                          final_rounds = 3,
                          num_init = 3,
                          num_nbrs = 4,
                          trace = "round",
                          trace.every = 1,
                          seed = 123)

  expect_s3_class(tr, "grip_layout_trace")
  expect_equal(dim(tr$final), c(8, 2))
  expect_true(all(is.finite(tr$final)))
  expect_true(length(tr$frames) >= 3L)
  expect_equal(names(tr$meta),
               c("frame", "phase", "level_index", "misf_level", "round_in_level", "active_vertices"))
  expect_identical(tr$meta$phase[[1]], "init")
  expect_identical(tail(tr$meta$phase, 1), "final")
  expect_true(any(tr$meta$phase == "round"))
  expect_true(all(diff(tr$meta$active_vertices) >= 0))

  first_frame <- tr$frames[[1L]]
  expect_equal(dim(first_frame), c(8, 2))
  expect_true(any(rowSums(is.na(first_frame)) == 2L))
  expect_equal(sum(rowSums(is.na(first_frame)) == 2L),
               8L - tr$meta$active_vertices[[1L]])
  expect_equal(tr$frames[[length(tr$frames)]], tr$final)
})

test_that("grip.layout.trace final matches grip.layout", {
  edges <- edges.path(8)
  tr_primary <- grip.layout.trace(edges = edges,
                                  n = 8,
                                  dim = 2,
                                  trace = "level",
                                  trace.every = 1,
                                  seed = 41)
  coords <- grip.layout(edges = edges,
                        n = 8,
                        dim = 2,
                        seed = 41)
  expect_identical(tr_primary$final, coords)
})

test_that("grip.layout.trace final matches grip.layout for kk_repulse mode", {
  edges <- edges.mesh(5, 5)
  tr_primary <- grip.layout.trace(edges = edges,
                                  n = 25,
                                  dim = 2,
                                  rounds = 8,
                                  final_rounds = 8,
                                  num_init = 6,
                                  num_nbrs = 8,
                                  coarse_repulsion_factor = 0.3,
                                  coarse_repulsion_sample = 8,
                                  coarse_repulsion_exact_below = 32,
                                  final_mode = "kk_repulse",
                                  trace = "level",
                                  trace.every = 1,
                                  seed = 41)
  coords <- grip.layout(edges = edges,
                        n = 25,
                        dim = 2,
                        rounds = 8,
                        final_rounds = 8,
                        num_init = 6,
                        num_nbrs = 8,
                        coarse_repulsion_factor = 0.3,
                        coarse_repulsion_sample = 8,
                        coarse_repulsion_exact_below = 32,
                        final_mode = "kk_repulse",
                        seed = 41)
  expect_identical(tr_primary$final, coords)
})

test_that("grip.layout.trace final matches grip.layout for structural FR knobs", {
  edges <- edges.mesh(5, 5)
  tr_primary <- grip.layout.trace(edges = edges,
                                  n = 25,
                                  dim = 2,
                                  rounds = 8,
                                  final_rounds = 8,
                                  num_init = 6,
                                  num_nbrs = 8,
                                  coarse_repulsion_factor = 0.3,
                                  coarse_repulsion_sample = 8,
                                  coarse_repulsion_exact_below = 32,
                                  final_anchor_factor = 1,
                                  final_move_scale_after_first = 0.5,
                                  final_mode = "fr",
                                  trace = "level",
                                  trace.every = 1,
                                  seed = 41)
  coords <- grip.layout(edges = edges,
                        n = 25,
                        dim = 2,
                        rounds = 8,
                        final_rounds = 8,
                        num_init = 6,
                        num_nbrs = 8,
                        coarse_repulsion_factor = 0.3,
                        coarse_repulsion_sample = 8,
                        coarse_repulsion_exact_below = 32,
                        final_anchor_factor = 1,
                        final_move_scale_after_first = 0.5,
                        final_mode = "fr",
                        seed = 41)
  expect_identical(tr_primary$final, coords)
})

test_that("grip.layout.trace final matches grip.layout for level-0 insertion knobs", {
  edges <- edges.sierpinski.carpet(2)
  n <- max(edges)
  tr_primary <- grip.layout.trace(edges = edges,
                                  n = n,
                                  dim = 2,
                                  rounds = 8,
                                  final_rounds = 8,
                                  num_init = 6,
                                  num_nbrs = 8,
                                  coarse_repulsion_factor = 0.3,
                                  coarse_repulsion_sample = 8,
                                  coarse_repulsion_exact_below = 32,
                                  level0_insertion_mode = "least_squares",
                                  level0_anchor_count = 6,
                                  level0_local_kk_steps = 0,
                                  trace = "level",
                                  trace.every = 1,
                                  seed = 41)
  coords <- grip.layout(edges = edges,
                        n = n,
                        dim = 2,
                        rounds = 8,
                        final_rounds = 8,
                        num_init = 6,
                        num_nbrs = 8,
                        coarse_repulsion_factor = 0.3,
                        coarse_repulsion_sample = 8,
                        coarse_repulsion_exact_below = 32,
                        level0_insertion_mode = "least_squares",
                        level0_anchor_count = 6,
                        level0_local_kk_steps = 0,
                        seed = 41)
  expect_identical(tr_primary$final, coords)
})

test_that("trace can append per-frame diagnostics against a target", {
  edges <- edges.mesh(4, 4)
  n <- max(edges)
  target <- grip.layout(edges = edges, n = n, dim = 2, seed = 17)
  tr <- grip.layout.trace(edges = edges,
                          n = n,
                          dim = 2,
                          trace = "level",
                          trace.every = 1,
                          diagnostics = "light",
                          target_coords = target,
                          diagnostic_sample_size_nonedge = 50,
                          seed = 17)

  expect_true(is.data.frame(tr$diagnostics))
  expect_equal(nrow(tr$diagnostics), length(tr$frames))
  expect_true(all(c(
    "frame", "phase", "level_index", "misf_level", "round_in_level",
    "active_vertices", "active.edges", "edge.length.cv",
    "median.edge.length", "sampled.nonedge.sep.ratio",
    "sampled.stress", "procrustes.rmse"
  ) %in% names(tr$diagnostics)))
  expect_lt(tail(tr$diagnostics$procrustes.rmse, 1), 1e-8)
  expect_true(all(is.na(tr$diagnostics$sampled.stress)))
})

test_that("level trace thins level-start snapshots and keeps endpoints", {
  edges <- edges.cycle(12)
  tr_dense <- grip.layout.trace(edges = edges,
                                n = 12,
                                dim = 2,
                                placement = "barycenter",
                                rounds = 5,
                                final_rounds = 3,
                                num_init = 4,
                                num_nbrs = 5,
                                trace = "level",
                                trace.every = 1,
                                seed = 77)
  tr_sparse <- grip.layout.trace(edges = edges,
                                 n = 12,
                                 dim = 2,
                                 placement = "barycenter",
                                 rounds = 5,
                                 final_rounds = 3,
                                 num_init = 4,
                                 num_nbrs = 5,
                                 trace = "level",
                                 trace.every = 2,
                                 seed = 77)

  expect_true(all(tr_dense$meta$phase %in% c("init", "level_start", "final")))
  expect_true(all(tr_sparse$meta$phase %in% c("init", "level_start", "final")))
  expect_identical(tr_dense$meta$phase[[1]], "init")
  expect_identical(tail(tr_dense$meta$phase, 1), "final")
  expect_identical(tr_sparse$meta$phase[[1]], "init")
  expect_identical(tail(tr_sparse$meta$phase, 1), "final")
  expect_true(length(tr_sparse$frames) <= length(tr_dense$frames))
  expect_equal(tr_dense$final, tr_sparse$final)
})

test_that("trace rejects disconnected graphs for now", {
  edges <- rbind(
    cbind(1:2, 2:3),
    cbind(5:6, 6:7)
  )
  expect_error(
    grip.layout.trace(edges = edges,
                      n = 7,
                      dim = 2,
                      trace = "round",
                      seed = 7),
    "currently supports only connected graphs"
  )
})

test_that("trace validates tuning parameters", {
  edges <- edges.path(8)
  expect_error(
    grip.layout.trace(edges = edges,
                      n = 8,
                      dim = 2,
                      num_nbrs = 0,
                      trace = "round",
                      seed = 123),
    "num_nbrs must be a positive integer"
  )
  expect_error(
    grip.layout.trace(edges = edges,
                      n = 8,
                      dim = 2,
                      repulsion_factor = -0.1,
                      trace = "round",
                      seed = 123),
    "repulsion_factor must be >= 0"
  )
  expect_error(
    grip.layout.trace(edges = edges,
                      n = 8,
                      dim = 2,
                      coarse_repulsion_factor = -0.1,
                      trace = "round",
                      seed = 123),
    "coarse_repulsion_factor must be >= 0"
  )
  expect_error(
    grip.layout.trace(edges = edges,
                      n = 8,
                      dim = 2,
                      final_anchor_factor = -0.1,
                      trace = "round",
                      seed = 123),
    "final_anchor_factor must be >= 0"
  )
  expect_error(
    grip.layout.trace(edges = edges,
                      n = 8,
                      dim = 2,
                      final_move_scale_after_first = 1.1,
                      trace = "round",
                      seed = 123),
    "final_move_scale_after_first must be in \\[0, 1\\]"
  )
  expect_error(
    grip.layout.trace(edges = edges,
                      n = 8,
                      dim = 2,
                      diagnostics = "light",
                      target_coords = matrix(0, nrow = 7, ncol = 2),
                      trace = "round",
                      seed = 123),
    "nrow\\(target_coords\\) must match the graph size"
  )
  expect_error(
    grip.layout.trace(edges = edges,
                      n = 8,
                      dim = 2,
                      final_mode = "banana",
                      trace = "round",
                      seed = 123),
    "'arg' should be one of"
  )
  expect_error(
    grip.layout.trace(edges = edges,
                      n = 8,
                      dim = 2,
                      level0_insertion_mode = "banana",
                      trace = "round",
                      seed = 123),
    "'arg' should be one of"
  )
  expect_error(
    grip.layout.trace(edges = edges,
                      n = 8,
                      dim = 2,
                      level0_anchor_count = 0,
                      trace = "round",
                      seed = 123),
    "level0_anchor_count must be a positive integer"
  )
  expect_error(
    grip.layout.trace(edges = edges,
                      n = 8,
                      dim = 2,
                      level0_local_kk_steps = -1,
                      trace = "round",
                      seed = 123),
    "level0_local_kk_steps must be a non-negative integer"
  )
})

test_that("trace repulsion_factor changes the final layout with a fixed seed", {
  edges <- edges.mesh(5, 5)
  tr_none <- grip.layout.trace(edges = edges,
                               n = 25,
                               dim = 2,
                               placement = "barycenter",
                               rounds = 8,
                               final_rounds = 8,
                               num_init = 6,
                               num_nbrs = 8,
                               repulsion_factor = 0,
                               trace = "level",
                               trace.every = 1,
                               seed = 29)
  tr_more <- grip.layout.trace(edges = edges,
                               n = 25,
                               dim = 2,
                               placement = "barycenter",
                               rounds = 8,
                               final_rounds = 8,
                               num_init = 6,
                               num_nbrs = 8,
                               repulsion_factor = 2,
                               trace = "level",
                               trace.every = 1,
                               seed = 29)
  expect_gt(max(abs(tr_none$final - tr_more$final)), 1e-6)
})

test_that("trace carpet preset matches the explicit carpet tuning profile", {
  edges <- edges.sierpinski.carpet(2)
  n <- max(edges)
  tr_preset <- grip.layout.trace(edges = edges,
                                 n = n,
                                 dim = 2,
                                 preset = "carpet",
                                 trace = "level",
                                 trace.every = 1,
                                 seed = 41)
  tr_explicit <- grip.layout.trace(edges = edges,
                                   n = n,
                                   dim = 2,
                                   placement = "barycenter",
                                   rounds = 160,
                                   final_rounds = 288,
                                   num_init = 28,
                                   num_nbrs = 24,
                                   r = 0.03,
                                   s = 6.0,
                                   repulsion_factor = 2.5,
                                   trace = "level",
                                   trace.every = 1,
                                   seed = 41)
  expect_identical(tr_preset$final, tr_explicit$final)
  expect_identical(tr_preset$meta, tr_explicit$meta)
})

test_that("trace mesh preset matches the explicit mesh tuning profile", {
  edges <- edges.mesh(4, 4)
  n <- max(edges)
  tr_preset <- grip.layout.trace(edges = edges,
                                 n = n,
                                 dim = 2,
                                 preset = "mesh",
                                 trace = "level",
                                 trace.every = 1,
                                 seed = 51)
  tr_explicit <- grip.layout.trace(edges = edges,
                                   n = n,
                                   dim = 2,
                                   placement = "barycenter",
                                   rounds = 128,
                                   final_rounds = 128,
                                   num_init = 12,
                                   num_nbrs = 20,
                                   r = 0.10,
                                   s = 4.5,
                                   repulsion_factor = 1.5,
                                   trace = "level",
                                   trace.every = 1,
                                   seed = 51)
  expect_identical(tr_preset$final, tr_explicit$final)
  expect_identical(tr_preset$meta, tr_explicit$meta)
})

test_that("trace torus preset matches the explicit torus tuning profile", {
  edges <- edges.torus(4, 4)
  n <- max(edges)
  tr_preset <- grip.layout.trace(edges = edges,
                                 n = n,
                                 dim = 3,
                                 preset = "torus",
                                 trace = "level",
                                 trace.every = 1,
                                 seed = 59)
  tr_explicit <- grip.layout.trace(edges = edges,
                                   n = n,
                                   dim = 3,
                                   placement = "barycenter",
                                   rounds = 192,
                                   final_rounds = 288,
                                   num_init = 12,
                                   num_nbrs = 28,
                                   r = 0.05,
                                   s = 7.5,
                                   repulsion_factor = 0.75,
                                   trace = "level",
                                   trace.every = 1,
                                   seed = 59)
  expect_identical(tr_preset$final, tr_explicit$final)
  expect_identical(tr_preset$meta, tr_explicit$meta)
})

test_that("trace tree preset matches the explicit tree tuning profile", {
  edges <- edges.kary.tree(k = 2, depth = 4)
  n <- max(edges)
  tr_preset <- grip.layout.trace(edges = edges,
                                 n = n,
                                 dim = 2,
                                 preset = "tree",
                                 trace = "level",
                                 trace.every = 1,
                                 seed = 71)
  tr_explicit <- grip.layout.trace(edges = edges,
                                   n = n,
                                   dim = 2,
                                   placement = "circle",
                                   rounds = 64,
                                   final_rounds = 160,
                                   num_init = 28,
                                   num_nbrs = 8,
                                   r = 0.05,
                                   s = 7.5,
                                   repulsion_factor = 0,
                                   trace = "level",
                                   trace.every = 1,
                                   seed = 71)
  expect_identical(tr_preset$final, tr_explicit$final)
  expect_identical(tr_preset$meta, tr_explicit$meta)
})

test_that("trace tree preset uses barycenter placement in 3D without warning", {
  edges <- edges.kary.tree(k = 2, depth = 4)
  n <- max(edges)
  expect_no_warning({
    tr_preset <- grip.layout.trace(edges = edges,
                                   n = n,
                                   dim = 3,
                                   preset = "tree",
                                   trace = "level",
                                   trace.every = 1,
                                   seed = 73)
    tr_explicit <- grip.layout.trace(edges = edges,
                                     n = n,
                                     dim = 3,
                                     placement = "barycenter",
                                     rounds = 64,
                                     final_rounds = 160,
                                     num_init = 28,
                                     num_nbrs = 8,
                                     r = 0.05,
                                     s = 7.5,
                                     repulsion_factor = 0,
                                     trace = "level",
                                     trace.every = 1,
                                     seed = 73)
    expect_identical(tr_preset$final, tr_explicit$final)
    expect_identical(tr_preset$meta, tr_explicit$meta)
  })
})
