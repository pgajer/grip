test_that('surface samples integrate triangle area and centroid', {
  e <- new.env()
  sys.source(system.file('scripts', 'metric-mds-surface-alignment.R', package='grip'), e)
  X <- rbind(c(0,0,0), c(2,0,0), c(0,1,0))
  for (m in c(1L,3L,6L)) {
    s <- e$comparison.surface.samples(X, matrix(1:3,1), m)
    expect_equal(nrow(s$points), m^2)
    expect_equal(sum(s$weights), 1)
    expect_equal(colSums(s$points*s$weights), colMeans(X))
  }
})

test_that('surface registration improves relabeled surfaces without deforming them', {
  skip_if_not_installed('FNN')
  e <- new.env()
  sys.source(system.file('scripts', 'metric-mds-comparison.R', package='grip'), e)
  sys.source(system.file('scripts', 'metric-mds-surface-alignment.R', package='grip'), e)
  X <- rbind(c(0,0,0), c(2,0,0), c(0,1,0), c(.2,.3,1.5))
  triangles <- t(combn(4L,3L))
  # All four faces are present, so relabeling leaves the same tetrahedral surface.
  Y <- sweep(X[c(2,1,3,4), ] %*% diag(c(-1,1,1)), 2, c(3,-2,1), '+')
  fit <- e$comparison.surface.align(Y, X, triangles, max.iter=120L)
  expect_gt(fit$initial.rmse, .01)
  # This finite multistart local search need not find the exact global match.
  expect_lt(fit$rmse, fit$initial.rmse / 2)
  expect_equal(as.vector(dist(fit$coords)), as.vector(dist(Y)), tolerance=1e-10)
  expect_true(all(is.finite(fit$coords)))
  exact <- e$comparison.surface.align(sweep(X %*% diag(c(-1,1,1)), 2,
    c(3,-2,1), '+'), X, triangles, angles=0)
  expect_lt(exact$rmse, 1e-10)
  expect_equal(exact$coords, X, tolerance=1e-10)
  audit <- e$comparison.surface.audit(Y, X, triangles, fit)
  expect_true(all(audit$surface.alignment < audit$vertex.alignment))
})
