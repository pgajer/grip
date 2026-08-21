metric.dispatch.fixture <- function() {
  list(
    edges = cbind(1:7, 2:8),
    edge_lengths = c(1, 1.5, 0.75, 2, 1.25, 0.8, 1.1)
  )
}

metric.dispatch.layout.args <- function() {
  list(
    n = 8,
    dim = 2,
    rounds = 3,
    final_rounds = 3,
    num_init = 4,
    num_nbrs = 4,
    seed = 19
  )
}

test_that("grip dispatches hop metric without changing existing behavior", {
  graph <- metric.dispatch.fixture()
  args <- c(list(edges = graph$edges), metric.dispatch.layout.args())

  coords.default <- do.call(grip, args)
  coords.hop <- do.call(grip, c(args, list(metric = "hop")))
  coords.backend <- do.call(globalrep.grip, args)

  expect_identical(coords.default, coords.hop)
  expect_identical(coords.hop, coords.backend)
})

test_that("new metric arguments do not shift the existing positional API", {
  layout.formals <- names(formals(globalrep.grip))
  trace.formals <- names(formals(getFromNamespace("grip.trace.hop", "grip")))

  expect_identical(
    head(names(formals(grip)), length(layout.formals)),
    layout.formals
  )
  expect_identical(
    head(names(formals(trace.grip)), length(trace.formals)),
    trace.formals
  )
})

test_that("grip dispatches edge-length metric to the weighted backend", {
  graph <- metric.dispatch.fixture()
  args <- c(
    list(edges = graph$edges, edge_weights = graph$edge_lengths),
    metric.dispatch.layout.args()
  )

  coords.unified <- do.call(grip, c(args, list(metric = "edge_length")))
  coords.backend <- do.call(globalrep.weighted.grip, args)

  expect_identical(coords.unified, coords.backend)
  expect_false(exists("weighted.grip", envir = asNamespace("grip"), inherits = FALSE))
})

test_that("metric-specific arguments are validated at the unified interface", {
  graph <- metric.dispatch.fixture()

  expect_error(
    grip(
      graph$edges,
      n = 8,
      metric = "hop",
      metric_neighbor_cap = 20
    ),
    "metric_neighbor_cap is only available"
  )
  expect_error(
    grip(
      graph$edges,
      n = 8,
      metric = "hop",
      length_normalization = "none"
    ),
    "length_normalization is only available"
  )
  expect_error(
    grip(graph$edges, n = 8, metric = "edge_length"),
    "requires edge weights|edge_weights is required|weight_list is required"
  )
  expect_error(
    grip(graph$edges, n = 8, metric = "unknown"),
    "should be one of"
  )
})

test_that("trace.grip mirrors layout metric dispatch", {
  graph <- metric.dispatch.fixture()
  args <- c(
    list(edges = graph$edges, edge_weights = graph$edge_lengths),
    metric.dispatch.layout.args(),
    list(trace = "level", diagnostics = "none")
  )

  trace.unified <- do.call(
    trace.grip,
    c(args, list(metric = "edge_length"))
  )
  expect_s3_class(trace.unified$stage_trace, "data.frame")
  expect_type(trace.unified$stage_data, "list")
  expect_false(exists(
    "trace.weighted.grip",
    envir = asNamespace("grip"),
    inherits = FALSE
  ))

  hop.args <- args
  hop.args$edge_weights <- NULL
  trace.default <- do.call(trace.grip, hop.args)
  trace.hop <- do.call(trace.grip, c(hop.args, list(metric = "hop")))
  expect_identical(trace.default$final, trace.hop$final)
})
