recipe <- new.env(parent = globalenv())
sys.source(system.file('scripts', 'metric-mds-comparison.R', package = 'grip'), envir = recipe)

test_that('comparison ambient fixtures recover their intended public-API targets', {
  skip_if_not_installed('withr')
  cases <- recipe$comparison_cases(3L)
  expect_length(cases, 7L)
  for (case in cases) {
    p <- prepare.graph.geodesic.mds(edges = case$edges, n = case$n, edge_weights = case$weights)
    expect_true(all(is.finite(p$distance_matrix)))
    expect_equal(p$distance_matrix, t(p$distance_matrix))
    if (case$target == 'Ambient Euclidean') {
      expect_equal(unname(p$distance_matrix), unname(as.matrix(dist(case$X))), tolerance = 1e-10)
    } else {
      expect_true(all(p$distance_matrix >= as.matrix(dist(case$X)) - 1e-10))
    }
    expect_true(all(case$weights > 0))
    expect_equal(ncol(recipe$comparison_start(case, 3L)), case$dimension)
  }
  expect_equal(cases[['saddle-9']]$X[,3], .8 * (cases[['saddle-9']]$X[,1]^2 - cases[['saddle-9']]$X[,2]^2))
})

test_that('overlay alignment preserves distances and does not scale fits', {
  X <- matrix(c(0,0,0, 1,0,0, 0,2,0, 0,0,3), byrow = TRUE, ncol = 3)
  reflection <- diag(c(-1,1,1))
  moved <- sweep(2 * X %*% reflection, 2, c(4,5,6), '+')
  aligned <- recipe$comparison_align(moved, X)
  expect_equal(as.vector(dist(aligned)), as.vector(dist(moved)), tolerance = 1e-12)
  expect_equal(colMeans(aligned), colMeans(X), tolerance = 1e-12)
  expect_gt(sum((as.vector(dist(aligned)) - as.vector(dist(X)))^2), 0)
})

test_that('comparison failures and warnings remain explicit and do not stop later fits', {
  skip_if_not_installed('withr')
  case <- recipe$comparison_cases(3L)[['saddle-9']]
  p <- prepare.graph.geodesic.mds(edges = case$edges, n = case$n, edge_weights = case$weights)
  failed <- recipe$comparison_fit(case, p, 1L, 10L, 'sgd', fit_fun = function(...) {
    warning('diagnostic warning'); stop('deliberate failed fit')
  })
  expect_identical(failed$row$status, 'error')
  expect_match(failed$row$warnings, 'diagnostic warning')
  expect_match(failed$row$failure, 'deliberate failed fit')
  expect_null(failed$coords)
  ok <- recipe$comparison_fit(case, p, 1L, 10L, 'sgd')
  expect_identical(ok$row$status, 'ok')
  expect_identical(ok$row$termination, 'iteration_limit')
  expect_match(ok$row$warnings, 'iteration_limit')
  expect_equal(recipe$comparison_score(ok$coords, p$distance_matrix)[['raw_stress']], ok$row$raw_stress)
})
