test_that("streamed continuation fingerprints preserve the legacy byte hash", {
  withr::local_preserve_seed()
  set.seed(73)
  before <- .Random.seed
  shared <- new.env(parent = emptyenv())
  shared$self <- shared
  shared$value <- c(NA_real_, NaN, Inf, -Inf, -0, 0, pi)
  objects <- list(NULL, logical(), integer(), raw(), shared$value,
    c("plain", enc2utf8("\u03bc"), NA_character_),
    structure(matrix(seq_len(12), 3), note = "attributes survive"),
    list(first = shared, second = shared), data.frame(x = 1:7, y = letters[1:7]))
  for (object in objects) {
    expected <- grip_sgd_state_hash_cpp(serialize(object, NULL, version = 2))
    expect_identical(grip_sgd_state_hash_object_cpp(object), expected)
  }
  expect_identical(.Random.seed, before)

  fit <- metric.mds(edges = edges.cycle(12), n = 12,
    approximation = "sparse", backend = "sgd", sparse.control = list(n.pivots = 4),
    dim = 3, max.iter = 3, seed = 19,
    sgd.control = list(scheduler = "zheng", schedule.epochs = 10, keep.state = TRUE))
  state <- fit$continuation$state
  expect_identical(fit$continuation$checksum,
    grip_sgd_state_hash_cpp(serialize(state, NULL, version = 2)))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  restored <- readRDS(path)
  expect_identical(.metric.mds.state.hash(restored$continuation$state),
    fit$continuation$checksum)
})

test_that("bounded weight extrema preserve complete Zheng rate sequences", {
  first <- c(0, 1e-200, .125, 7, NA_real_, Inf)
  second <- c(NaN, 0, 1e200, .25, 3)
  range <- grip_sgd_weight_range_cpp(first, second)
  expect_identical(range, c(1e-200, 1e200))
  expect_identical(grip_sgd_weight_range_cpp(c(0, NA_real_), c(Inf, NaN)), numeric())
  for (epochs in c(20L, 500L, 1000L)) {
    control <- grip.mds.sgd.control(list(scheduler = "zheng", schedule.epochs = 500L), epochs)
    expect_identical(grip.mds.sgd.rates(control, epochs, range),
      grip.mds.sgd.rates(control, epochs, c(first, second)))
  }
})
