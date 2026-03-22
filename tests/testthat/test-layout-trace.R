test_that("round trace returns frames, metadata, and inactive NA rows", {
  edges <- edges.path(8)
  tr <- grip.layout.trace(edges = edges,
                          n = 8,
                          dim = 2,
                          engine = "mish_v5",
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

test_that("level trace thins level-start snapshots and keeps endpoints", {
  edges <- edges.cycle(12)
  tr_dense <- grip.layout.trace(edges = edges,
                                n = 12,
                                dim = 2,
                                engine = "mish_v6",
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
                                 engine = "mish_v6",
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
