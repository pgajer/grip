continued_fit <- function(...) suppressWarnings(metric.mds(...))
continuation_fixture <- function() {
  x <- cbind(seq_len(9)/4,sin(seq_len(9)),cos(seq_len(9)*.7))
  list(x=x,d=as.matrix(dist(x)),start=x[,1:2])
}

test_that("full continuation preserves RNG, terminal state and schedule across RDS", {
  z<-continuation_fixture()
  for (weights in c("uniform","inverse_squared")) for (scheduler in c("zheng","hybrid","exponential")) {
    control<-list(keep.state=TRUE,schedule.epochs=17,checkpoint.every=3,scheduler=scheduler)
    a<-continued_fit(distance.matrix=z$d,init=z$start,max.iter=7,pair.weights=weights,sgd.control=control)
    full<-continued_fit(distance.matrix=z$d,init=z$start,max.iter=28,pair.weights=weights,sgd.control=control)
    path<-tempfile(fileext=".rds");saveRDS(a,path);a<-readRDS(path);unlink(path)
    set.seed(19);before<-.Random.seed
    b<-metric.mds(resume=a,max.iter=5)
    result<-metric.mds(resume=b,max.iter=16)
    expect_identical(.Random.seed,before)
    expect_equal(result$coords,full$coords,tolerance=0)
    expect_identical(result$continuation$state$native,full$continuation$state$native)
    expect_equal(result$metadata$raw_stress,full$metadata$raw_stress,tolerance=0)
    expect_identical(result$metadata$history,result$trace)
    # A pause at epoch 7 is recorded, but never changes selection/stopping.
    canonical<-function(x) {t<-x$trace[x$trace$scheduled==1L,,drop=FALSE];rownames(t)<-NULL;t}
    expect_equal(canonical(result),canonical(full),tolerance=0)
    expect_equal(result$metadata$stopping$epochs_completed,28)
    expect_true(result$metadata$stopping$schedule_completed)
    expect_false(result$metadata$stopping$assessed)
    expect_false(result$metadata$converged)
    expect_equal(result$metadata$stopping$returned_iterate,"best_profiled_checkpoint")
  }
})

test_that("sparse graph continuation does not reprepare pivots or reset the seed", {
  control<-list(keep.state=TRUE,schedule.epochs=13,checkpoint.every=4,scheduler="zheng")
  args<-list(edges=edges.cycle(12),n=12,approximation="sparse",dim=3,
    sparse.control=list(n.pivots=4),seed=91,sgd.control=control)
  a<-do.call(metric.mds,c(args,list(max.iter=5)))
  full<-do.call(metric.mds,c(args,list(max.iter=22)))
  local_mocked_bindings(prepare.edge.kk=function(...)stop("reprepared"),.package="grip")
  b<-metric.mds(resume=a,max.iter=17)
  expect_equal(b$coords,full$coords,tolerance=0)
  expect_identical(b$continuation$state$native,full$continuation$state$native)
  expect_identical(b$metadata$sparse,full$metadata$sparse)
  expect_equal(b$metadata$sparse_proxy_stress,full$metadata$sparse_proxy_stress,tolerance=0)
  expect_identical(b$metadata$history,b$trace)
  expect_equal(b$metadata$stopping$returned_iterate,"terminal")
})

test_that("observation preparation metadata survives continuation", {
  skip_if_not_installed("dgraphs")
  if(!"search.data" %in% names(formals(dgraphs::create.sknn.graph)))skip("new dgraphs required")
  z<-continuation_fixture()
  p<-prepare.sparse.mds(z$x,n.pivots=3,q=2,neighbor.variance.target=NULL,seed=7)
  control<-list(keep.state=TRUE,schedule.epochs=11)
  a<-metric.mds(prepared=p,approximation="sparse",max.iter=4,init=z$start,sgd.control=control)
  full<-metric.mds(prepared=p,approximation="sparse",max.iter=16,init=z$start,sgd.control=control)
  b<-metric.mds(resume=a,max.iter=12)
  expect_equal(b$coords,full$coords,tolerance=0)
  expect_identical(b$metadata$data_preparation,a$metadata$data_preparation)
  expect_identical(b$metadata$sparse$pivots,p$pivots)
  expect_identical(b$metadata$sparse$region,p$region)
})

test_that("practical stopping requires patience and survives chunk boundaries", {
  x<-rbind(c(0,0),c(1,0),c(0,1),c(2,3))
  for(cp in c(1,4)) {
    control<-list(keep.state=TRUE,schedule.epochs=30,checkpoint.every=cp,
      stopping=list(min.epochs=5,patience=3))
    full<-continued_fit(distance.matrix=as.matrix(dist(x)),init=x,max.iter=30,sgd.control=control)
    a<-continued_fit(distance.matrix=as.matrix(dist(x)),init=x,max.iter=3,sgd.control=control)
    b<-metric.mds(resume=a,max.iter=27)
    expect_identical(b$continuation$state$native,full$continuation$state$native)
    expect_equal(b$coords,full$coords,tolerance=0)
    expect_equal(full$metadata$stopping$epochs_completed,if(cp==1)7 else 16)
    expect_equal(full$metadata$termination,"practical_stabilization")
    expect_true(full$metadata$stopping$assessed)
    expect_true(full$metadata$stopping$stabilized)
    expect_false(full$metadata$stopping$certified)
    expect_false(full$metadata$converged)
    expect_error(metric.mds(resume=b,max.iter=1),"already met")
  }
})

test_that("tiny learning rates alone cannot manufacture stabilization", {
  x<-rbind(c(0,0),c(3,0),c(0,4),c(5,6))
  rule<-list(relative.tol=1e-4,motion.tol=1e-4,residual.tol=1e-4,patience=2L,min.epochs=1L)
  for (retain in c(TRUE,FALSE)) {
    fit<-grip_sgd_mds_cpp(x,c(1,2,3,1,2,1),rep(1e-30,8),1L,1L,1e6,
      retain_best=retain,stopping=rule)
    expect_false(fit$stabilized)
    expect_true(all(fit$trace$centered_motion[-1]<1e-4))
    expect_true(all(fit$trace$relative_change[-1]<1e-4))
    expect_true(all(fit$trace$update_residual[-1]>1e-4))
    expect_equal(fit$niter,8)
  }
})

test_that("state and control errors fail before numerical continuation", {
  z<-continuation_fixture()
  a<-continued_fit(distance.matrix=z$d,max.iter=2,sgd.control=list(keep.state=TRUE))
  old<-a;old$continuation<-NULL
  expect_error(metric.mds(resume=old),"historical fits")
  expect_error(metric.mds(resume=a,seed=4),"frozen")
  expect_error(metric.mds(resume=a,distance.matrix=z$d),"frozen")
  expect_error(metric.mds(resume=a,sgd.control=list()),"frozen")
  expect_error(metric.mds(resume=a,max.iter=0),"positive")
  bad<-a;bad$continuation$state$arguments$targets[1]<-9
  expect_error(metric.mds(resume=bad),"integrity")
  bad<-a;bad$continuation$state$version<-"future"
  expect_error(metric.mds(resume=bad),"version")
  expect_error(metric.mds(distance.matrix=z$d,n.init=2,sgd.control=list(keep.state=TRUE)),"n.init=1")
  expect_error(metric.mds(distance.matrix=z$d,sgd.control=list(stopping=list(patience=0))),"positive")
  expect_error(metric.mds(distance.matrix=z$d,sgd.control=list(stopping=list(residual.tol=-1))),"residual.tol")
  expect_error(metric.mds(distance.matrix=z$d,sgd.control=list(stopping=list(typo=1))),"stopping")
})

test_that("default runs retain histories and explicitly decline convergence assessment", {
  z<-continuation_fixture()
  a<-continued_fit(distance.matrix=z$d,max.iter=4,n.init=2)
  expect_null(a$continuation)
  expect_identical(a$trace,a$metadata$history)
  expect_equal(length(a$metadata$sgd),2)
  expect_equal(a$metadata$stopping$reason,"epoch_budget")
  expect_false(a$metadata$stopping$assessed)
  expect_true(a$metadata$stopping$schedule_completed)
  expect_equal(a$metadata$termination,"iteration_limit")
})

test_that("endpoint residual matches an independent asymmetric calculation", {
  x <- rbind(c(0,0),c(2,1),c(-1,3))
  pairs <- rbind(c(1L,2L),c(1L,3L),c(2L,3L))
  target <- c(1,2,1); wi <- c(1,3,0); wj <- c(2,0,4)
  rule <- list(relative.tol=0,motion.tol=0,residual.tol=0,patience=2L,min.epochs=1L)
  fit <- grip_sgd_mds_cpp(x,target,1e-30,1L,1L,1e6,weights=wi,
    endpoints=pairs,reverse_weights=wj,retain_best=FALSE,stopping=rule)
  force <- matrix(0,3,2); total <- numeric(3)
  for(k in seq_along(target)) {
    i <- pairs[k,1];j <- pairs[k,2];v <- x[i,]-x[j,];d <- sqrt(sum(v*v))
    force[i,] <- force[i,]-wi[k]*(d-target[k])*v/d
    force[j,] <- force[j,]+wj[k]*(d-target[k])*v/d
    total[i] <- total[i]+wi[k];total[j] <- total[j]+wj[k]
  }
  expect_equal(tail(fit$trace$update_residual,1),max(sqrt(rowSums((force/total)^2))))
  expect_false(fit$stabilized)
})

test_that("a failed checkpoint resets the practical stopping streak", {
  rule <- list(relative.tol=.01,motion.tol=.01,residual.tol=10,patience=2L,min.epochs=1L)
  fit <- grip_sgd_mds_cpp(rbind(c(0,0),c(4,0)),1,
    c(.5,.5,1e-4,.5,1e-4,1e-4,1e-4),1L,1L,1e6,retain_best=FALSE,stopping=rule)
  expect_equal(fit$trace$stable_checkpoints,c(0,0,0,1,0,1,2))
  expect_equal(fit$niter,6)
  expect_true(fit$stabilized)
})

test_that("graph diagnostics and sparse collision RNG survive continuation", {
  args <- list(edges=edges.cycle(8),n=8,init="random",diagnostics=TRUE,
    sgd.control=list(keep.state=TRUE,schedule.epochs=12))
  a <- do.call(continued_fit,c(args,list(max.iter=3)))
  full <- do.call(continued_fit,c(args,list(max.iter=12)))
  set.seed(71); before <- .Random.seed
  b <- metric.mds(resume=a,max.iter=9)
  expect_identical(.Random.seed,before)
  expect_equal(b$diagnostics,full$diagnostics,tolerance=0)
  expect_equal(b$coords,full$coords,tolerance=0)
  args <- list(edges=edges.cycle(8),n=8,approximation="sparse",init=matrix(0,8,2),
    sparse.control=list(n.pivots=3),sgd.control=list(keep.state=TRUE,schedule.epochs=9))
  a <- do.call(metric.mds,c(args,list(max.iter=2)))
  full <- do.call(metric.mds,c(args,list(max.iter=13)))
  b <- metric.mds(resume=a,max.iter=11)
  expect_identical(b$continuation$state$native,full$continuation$state$native)
  expect_equal(b$coords,full$coords,tolerance=0)
})

test_that("public sparse stopping is scale invariant and explicitly uncertified", {
  x <- rbind(c(0,0),c(1,0),c(0,1),c(2,3))
  pairs <- t(combn(4,2)); storage.mode(pairs) <- "integer"
  control <- list(keep.state=TRUE,stopping=list(min.epochs=2,patience=2))
  for(scale in c(1,1e5)) {
    fit <- metric.mds(n=4,approximation="sparse",init=x*scale,max.iter=20,
      constraints=list(pairs=pairs,targets=as.double(dist(x))*scale),sgd.control=control)
    expect_equal(fit$metadata$stopping$epochs_completed,3)
    expect_true(fit$metadata$stopping$stabilized)
    expect_false(fit$metadata$converged)
    expect_output(print(fit),"practical stabilization rule met")
    expect_output(print(fit),"convergence: not certified")
  }
  expect_error(grip.mds.sgd.control(list(schedule.epochs=1000000,max.workspace.bytes=10000),1),"rate vector")
  control <- grip.mds.sgd.control(list(schedule.epochs=4),4)
  rates <- grip.mds.sgd.rates(control,7)
  expect_equal(rates[5:7],rates[4]/(1+(1:3)/4))
})
