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

test_that('the large sparse saddle mesh follows the saved grid and graph', {
  e <- new.env()
  sys.source(system.file('scripts', 'sparse-mds-comparison.R', package='grip'), e)
  bundle <- readRDS(system.file('extdata','sparse-mds-comparison','benchmark.rds',package='grip'))
  case <- bundle$cases[['saddle-4096']]
  triangles <- e$sparse.comparison.surface.triangles(case)
  expect_equal(dim(triangles), c(7938L,3L))
  expect_equal(sort(unique(as.vector(triangles))), seq_len(case$n))
  key <- function(E) paste(pmin(E[,1],E[,2]),pmax(E[,1],E[,2]),sep=':')
  sides <- rbind(triangles[,1:2],triangles[,2:3],triangles[,c(3,1)])
  expect_true(all(key(sides) %in% key(case$edges)))
  # Grid recovery must not depend on the stored vertex order.
  reordered <- case; reordered$X <- case$X[case$n:1L,]
  reindexed <- e$sparse.comparison.surface.triangles(reordered)
  expect_equal(reordered$X[as.vector(reindexed),],case$X[as.vector(triangles),])
})

test_that('sparse saddle displays align surfaces after applying the saved scale', {
  skip_if_not_installed('FNN'); skip_if_not_installed('withr')
  e <- new.env()
  for (file in c('metric-mds-comparison.R','metric-mds-surface-alignment.R',
                 'sparse-mds-comparison.R'))
    sys.source(system.file('scripts',file,package='grip'), e)
  # Capture the view inputs without requiring WebGL or optional ivue.
  e$comparison_view <- function(case, fits, ...) c(list(fits=fits),list(...))
  case <- e$comparison_cases(3L)[['saddle_graph-9']]
  coords <- sweep(2*case$X %*% diag(c(-1,1,1)),2,c(3,4,5),'+')
  fits <- list(full=list(coords=coords,row=list(relative_scale=.5)))
  before <- serialize(fits,NULL)
  view <- e$sparse_comparison_view(case,fits)
  expect_identical(serialize(fits,NULL),before)
  expect_match(view$description,'Surface-shape alignment',fixed=TRUE)
  expect_equal(view$fits$full$coords,coords*.5)
  expect_equal(as.vector(dist(view$surface.alignments$full$coords)),
    as.vector(dist(coords*.5)),tolerance=1e-10)
  unscaled <- e$sparse_comparison_view(case,fits,adjust_scale=FALSE)
  expect_equal(as.vector(dist(unscaled$surface.alignments$full$coords)),
    as.vector(dist(coords)),tolerance=1e-10)
})
