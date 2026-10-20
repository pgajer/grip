trace_native <- function(x, targets, ...) {
  grip_sgd_mds_cpp(x, targets, c(.3,.15,.05), 17L, 1L, 1e7, ...)
}
without_recording <- function(x) {
  x$update_trace <- NULL
  x$workspace_bytes <- NULL
  x
}

test_that("tracing preserves every native result and RNG for dense and sparse fits", {
  x <- rbind(c(0,0),c(0,0),c(1,2),c(-1,1))
  target <- c(0,2,1,1,2,3)
  pairs <- t(utils::combn(4,2))
  for (shuffle in c(FALSE,TRUE)) for (sparse in c(FALSE,TRUE)) {
    args <- list(x=x,targets=target,shuffle=shuffle,keep_state=TRUE,
      weights=c(1,2,1,3,2,1),reverse_weights=if(sparse)c(0,1,2,1,0,2) else NULL,
      endpoints=if(sparse)pairs else NULL,retain_best=!sparse)
    off <- do.call(trace_native,args)
    on <- do.call(trace_native,c(args,list(update_trace=list(every=1,max.bytes=1e6))))
    expect_identical(without_recording(on),without_recording(off))
    expect_null(off$update_trace)
    expect_length(on$update_trace$coordinates,19)
    expect_identical(on$update_trace$coordinates[[1]],x)
    expect_identical(tail(on$update_trace$coordinates,1)[[1]],on$terminal_conf)
    expect_equal(on$update_trace$steps$pair_updates,0:18)
    for(i in 2:19) expect_lte(sum(rowSums(abs(on$update_trace$coordinates[[i]]-
      on$update_trace$coordinates[[i-1]]))>0),2L)
  }
})

test_that("fixed-order individual states match an independent pair-update calculation", {
  x <- rbind(c(0,0),c(3,0),c(0,4));targets<-c(2,3,4);pairs<-t(utils::combn(3,2))
  on <- trace_native(x,targets,shuffle=FALSE,update_trace=list(every=1,max.bytes=1e6))
  index<-1L
  for(rate in c(.3,.15,.05)) for(k in 1:3) {
    i<-pairs[k,1];j<-pairs[k,2];delta<-x[i,]-x[j,];distance<-sqrt(sum(delta^2))
    step<-(distance-targets[k])*.5*min(rate,1)*delta/distance
    x[i,]<-x[i,]-step;x[j,]<-x[j,]+step;index<-index+1L
    expect_equal(on$update_trace$coordinates[[index]],x,tolerance=1e-14)
  }
})

test_that("cadence includes epoch endpoints and no-op pairs", {
  x<-rbind(c(0,0),c(0,0),c(1,1));target<-c(0,2,2)
  on<-trace_native(x,target,shuffle=FALSE,update_trace=list(every=2,max.bytes=1e6))
  expect_equal(on$update_trace$steps$pair_updates,c(0,2,3,5,6,8,9))
  expect_equal(on$update_trace$steps$epoch,c(0,1,1,2,2,3,3))
  one<-grip_sgd_mds_cpp(x,target,.1,1L,1L,1e6,FALSE,
    update_trace=list(every=1,max.bytes=1e5))
  expect_identical(one$update_trace$coordinates[[1]],one$update_trace$coordinates[[2]])
  # A positive-distance collision consumes RNG; snapshots must not alter it.
  a<-trace_native(x,c(1,2,2),keep_state=TRUE)
  b<-trace_native(x,c(1,2,2),keep_state=TRUE,update_trace=list(every=9,max.bytes=1e6))
  expect_identical(without_recording(a),without_recording(b))
})

test_that("trace budgets and malformed controls fail before recording", {
  x<-matrix(1:8,4,2);targets<-as.double(dist(x))
  expect_error(trace_native(x,targets,update_trace=list(every=1,max.bytes=1)),"trace.*max.bytes")
  expect_error(grip_sgd_mds_cpp(x,targets,.1,1L,1L,1000,
    update_trace=list(every=1,max.bytes=1e7)),"workspace")
  for(v in list(0,-1,NA_real_,Inf,1.5))
    expect_error(trace_native(x,targets,update_trace=list(every=v,max.bytes=1e6)),"tracing")
  for(control in list(TRUE,list(unknown=1),list(every=0),list(every=NA),list(max.bytes=-1)))
    expect_error(grip.mds.sgd.control(list(trace=control),3))
  expect_null(grip.mds.sgd.control(list(),3)$trace)
})

test_that("public tracing preserves fits and continues with absolute counts", {
  set.seed(19);x<-matrix(rnorm(24),8,3);D<-as.matrix(dist(x))
  base<-list(distance.matrix=D,dim=2,init="random",n.init=1,max.iter=3,seed=42,
    backend="sgd",diagnostics=FALSE)
  ctl<-list(keep.state=TRUE,schedule.epochs=6)
  off<-suppressWarnings(do.call(metric.mds,c(base,list(sgd.control=ctl))))
  ctl$trace<-list(every=3)
  on<-suppressWarnings(do.call(metric.mds,c(base,list(sgd.control=ctl))))
  expect_identical(on$coords,off$coords)
  expect_identical(on$trace,off$trace)
  expect_identical(on$continuation$state$native,off$continuation$state$native)
  expect_identical(on$update_trace$scope,"current_call")
  expect_equal(tail(on$update_trace$coordinates,1)[[1]],
    on$continuation$state$native$coords*on$metadata$input_rms_distance)
  resumed<-metric.mds(resume=on,max.iter=3)
  base$max.iter<-6
  full<-suppressWarnings(do.call(metric.mds,c(base,list(sgd.control=ctl))))
  expect_identical(resumed$coords,full$coords)
  expect_identical(resumed$continuation$state$native,full$continuation$state$native)
  expect_equal(resumed$update_trace$steps$pair_updates[1],3*28)
  expect_equal(resumed$update_trace$steps$epoch[1],3)
  expect_equal(resumed$update_trace$coordinates,c(tail(on$update_trace$coordinates,1),
    full$update_trace$coordinates[(length(on$update_trace$coordinates)+1):length(full$update_trace$coordinates)]))
  base$n.init<-2
  expect_error(do.call(metric.mds,c(base,list(sgd.control=ctl))),"tracing requires n.init=1")
})

test_that("sparse tracing uses input units and practical stopping trims frames", {
  set.seed(28);x<-matrix(rnorm(36),12,3)
  prep<-prepare.sparse.mds(x,method="euclidean",n.pivots=3,q=3,seed=1)
  args<-list(prepared=prep,approximation="sparse",backend="sgd",dim=2,
    max.iter=3,seed=4,diagnostics=FALSE)
  off<-do.call(metric.mds,c(args,list(sgd.control=list(keep.state=TRUE))))
  on<-do.call(metric.mds,c(args,list(sgd.control=list(keep.state=TRUE,trace=list(every=7)))))
  expect_identical(on$coords,off$coords)
  expect_identical(on$trace,off$trace)
  expect_identical(on$continuation$state$native,off$continuation$state$native)
  final<-tail(on$update_trace$coordinates,1)[[1]]
  expect_equal(sweep(final,2,colMeans(final)),on$coords)
  later<-metric.mds(resume=on,max.iter=2)
  expect_equal(later$update_trace$steps$epoch[1],3)
  expect_equal(tail(later$update_trace$steps$epoch,1),5)
  targets<-as.double(dist(x))
  rule<-list(relative.tol=1e10,motion.tol=1e10,residual.tol=1e10,patience=1L,min.epochs=1L)
  early<-grip_sgd_mds_cpp(x,targets,rep(.1,20),1L,1L,1e7,stopping=rule,
    update_trace=list(every=20,max.bytes=1e7))
  expect_true(early$stabilized)
  expect_equal(early$niter,1)
  expect_equal(tail(early$update_trace$steps$epoch,1),1)
  expect_length(early$update_trace$coordinates,5)
})
