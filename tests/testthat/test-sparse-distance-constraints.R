test_that('explicit constraints preserve distances rather than taking graph paths', {
  E <- rbind(c(1,2),c(2,3),c(1,3))
  cst <- list(pairs=E,targets=c(1,1,9))
  fit <- metric.mds(n=3,approximation='sparse',constraints=cst,max.iter=2,seed=19)
  expect_identical(fit$metadata$constraint_source,'explicit_distances')
  expect_equal(fit$metadata$sparse$targets,c(1,9,1))
  expect_null(fit$prepared$distance_matrix)
  expect_equal(fit$prepared$graph_build_mode,'distance_constraints')
  expect_equal(fit$coords,metric.mds(n=3,approximation='sparse',constraints=cst,max.iter=2,seed=19)$coords)
  d <- fit$metadata$sparse$targets
  P <- fit$metadata$sparse$pairs
  stress <- sum((sqrt(rowSums((fit$coords[P[,1],]-fit$coords[P[,2],])^2))/d-1)^2)
  expect_equal(fit$metadata$sparse_proxy_stress,stress,tolerance=1e-10)
})

test_that('explicit constraints agree with existing full-pivot sparse updates', {
  E <- edges.path(5)
  ref <- metric.mds(edges=E,n=5,approximation='sparse',sparse.control=list(pivots=1:5),seed=77,max.iter=4)
  S <- ref$metadata$sparse
  fit <- metric.mds(n=5,approximation='sparse',constraints=S[c('pairs','targets','count_i','count_j')],seed=77,max.iter=4)
  expect_equal(fit$coords,ref$coords,tolerance=1e-12)
  # Canonicalizing reversed pairs must also reverse endpoint multiplicities.
  cst <- list(pairs=rbind(c(3,1),c(2,3)),targets=c(4,2),count_i=c(3,0),count_j=c(0,2))
  x <- .metric.mds.validate.constraints(cst,3,1e7)
  expect_equal(x$count_i,c(0,0)); expect_equal(x$count_j,c(3,2))
})

test_that('invalid and ambiguous explicit constraints fail before fitting', {
  cst <- list(pairs=rbind(c(1,2),c(2,3)),targets=c(1,2))
  run <- function(x=cst,...) metric.mds(n=3,approximation='sparse',constraints=x,...)
  expect_error(metric.mds(n=3,constraints=cst),'sparse approximation')
  expect_error(run(edges=edges.path(3)),'cannot accompany')
  expect_error(run(sparse.control=list(pivots=1)),'without pivot')
  expect_error(run(list(pairs=matrix(c(1,2),1),targets=1)),'connected')
  expect_error(run(modifyList(cst,list(targets=c(0,1)))),'positive finite')
  expect_error(run(modifyList(cst,list(targets=c(NA,1)))),'positive finite')
  expect_error(run(list(pairs=rbind(c(1,2),c(2,1)),targets=c(1,1))),'unique')
  expect_error(run(modifyList(cst,list(count_i=c(0,0),count_j=c(0,1)))),'positive finite total')
  expect_error(run(modifyList(cst,list(count_i=rep(8e307,2),count_j=rep(1e307,2)))),'positive finite total')
  expect_error(run(sgd.control=list(max.workspace.bytes=1)),'workspace')
  expect_error(run(modifyList(cst,list(pairs=rbind(c(1.5,2),c(2,3))))),'integer vertex')
})
