test_that('weighted native updates and checkpoint scores match an independent loop', {
  X <- rbind(c(0,0),c(2,1),c(-1,2),c(1,-2))
  D <- c(1,3,2,2,1,4); W <- c(8,.25,2,9,.5,.125); rates <- c(.7,.3,.1)
  pairs <- rbind(c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4))
  score <- function(Z) {
    d <- as.vector(dist(Z)); a <- sum(W*d*D)/sum(W*d^2)
    c(raw=sum(W*(d-D)^2), profiled=sum(W*(a*d-D)^2))
  }
  path <- list(X); losses <- list(score(X)); Z <- X
  for (rate in rates) {
    for (k in seq_len(nrow(pairs))) {
      i <- pairs[k,1]; j <- pairs[k,2]; v <- Z[i,]-Z[j,]; d <- sqrt(sum(v^2))
      step <- min(1,rate*W[k])*(d-D[k])*v/(2*d)
      Z[i,] <- Z[i,]-step; Z[j,] <- Z[j,]+step
    }
    path[[length(path)+1L]] <- Z; losses[[length(losses)+1L]] <- score(Z)
  }
  losses <- do.call(rbind,losses)
  native <- grip_sgd_mds_cpp(X,D,rates,7L,1L,1e6,FALSE,W)
  expect_equal(native$terminal_conf,Z,tolerance=1e-12)
  expect_equal(native$trace$raw_stress,unname(losses[,'raw']))
  expect_equal(native$trace$profiled_stress,unname(losses[,'profiled']))
  expect_equal(native$conf,path[[which.min(losses[,'profiled'])]])
  expect_equal(colMeans(native$terminal_conf),colMeans(X),tolerance=1e-12)
  expect_error(grip_sgd_mds_cpp(X,D,rates,7L,1L,1e6,FALSE,1), 'weight count')
  for (bad in c(0,-1,NA_real_,Inf)) {
    w <- W; w[2] <- bad
    expect_error(grip_sgd_mds_cpp(X,D,rates,7L,1L,1e6,FALSE,w), 'weights')
  }
})

test_that('both backends fit and report the inverse-squared objective consistently', {
  skip_if_not_installed('smacof')
  for (backend in c('sgd','smacof')) {
    args <- list(edges=edges.cycle(8),n=8,edge.weights=seq_len(8),dim=3,
                 init='random',n.init=3,seed=92,max.iter=80,backend=backend,
                 pair.weights='inverse_squared',diagnostics=FALSE)
    set.seed(79); saved <- .Random.seed
    fit <- suppressWarnings(do.call(metric.mds,args))
    expect_identical(.Random.seed,saved)
    d <- as.vector(dist(fit$coords)); target <- as.vector(as.dist(fit$prepared$distance_matrix))
    w <- 1/target^2; raw <- sum(((d-target)/target)^2)
    expect_identical(fit$metadata$pair_weights,'inverse_squared')
    expect_equal(fit$metadata$raw_stress,raw,tolerance=1e-10)
    expect_equal(fit$metadata$target_normalized_rmse,sqrt(mean(((d-target)/target)^2)),tolerance=1e-10)
    expect_equal(fit$metadata$stress1_identity,sqrt(raw/sum(w*d^2)),tolerance=1e-10)
    a <- sum(w*d*target)/sum(w*target^2)
    expect_equal(fit$metadata$stress1_profiled,sqrt(sum(w*(d-a*target)^2)/sum(w*d^2)),tolerance=1e-10)
    expect_equal(sum(w*d*target)/sum(w*d^2),1,tolerance=1e-10)
    expect_equal(raw,min(fit$metadata$starts$raw_stress),tolerance=1e-10)
    expect_true(all(fit$metadata$starts$raw_stress <= fit$metadata$starts$initial_raw_stress+1e-9))
    for (unit in c(1e-5,1e5)) {
      changed <- args; changed$edge.weights <- unit*args$edge.weights
      scaled <- suppressWarnings(do.call(metric.mds,changed))
      expect_equal(as.vector(dist(scaled$coords))/unit,d,tolerance=1e-7)
      expect_equal(scaled$metadata$raw_stress,raw,tolerance=1e-7)
    }
    if (backend=='sgd') {
      trace <- fit$metadata$sgd[[fit$metadata$selected_start]]$trace
      expect_equal(min(trace$profiled_stress),raw,tolerance=1e-9)
    }
    args$pair.weights <- 'uniform'
    uniform <- suppressWarnings(do.call(metric.mds,args))
    args$pair.weights <- NULL
    default <- suppressWarnings(do.call(metric.mds,args[!vapply(args,is.null,logical(1))]))
    expect_equal(uniform$coords,default$coords,tolerance=0)
    expect_equal(uniform$metadata$raw_stress,default$metadata$raw_stress,tolerance=0)
    expect_gt(max(abs(as.vector(dist(uniform$coords))-d)),1e-3)
  }
})

test_that('weighted SMACOF matches an explicit weight matrix and weighted rescaling', {
  skip_if_not_installed('smacof')
  p <- prepare.graph.geodesic.mds(edges.cycle(6),n=6,edge.weights=seq_len(6))
  target <- as.vector(as.dist(p$distance_matrix)); rms <- sqrt(mean(target^2))
  t <- target/rms; w <- 1/t^2
  start <- cbind(1:6,c(2,-1,3,0,1,-2)); X <- scale(start,scale=FALSE)/rms
  d <- as.vector(dist(X)); X <- X*sum(w*d*t)/sum(w*d^2)
  W <- 1/(p$distance_matrix/rms)^2; diag(W) <- 0
  direct <- suppressWarnings(smacof::mds(p$distance_matrix/rms,ndim=2,type='ratio',weightmat=W,
                        init=X,itmax=100,eps=1e-9,principal=FALSE))
  d <- as.vector(dist(direct$conf)); expected <- d*sum(w*d*target)/sum(w*d^2)
  fit <- suppressWarnings(metric.mds(prepared=p,dim=2,init=start,max.iter=100,
    eps=1e-9,backend='smacof',pair.weights='inverse_squared',diagnostics=FALSE))
  expect_equal(as.vector(dist(fit$coords)),expected,tolerance=1e-9)
})

test_that('weighted zero-distance and invalid-choice policies are explicit', {
  X <- rbind(c(0,0),c(0,0),c(1,0),c(0,1))
  p <- prepare.graph.geodesic.mds(edges.cycle(4),n=4)
  p$distance_matrix <- as.matrix(dist(X))
  for (backend in c('sgd','smacof')) {
    if (backend=='smacof' && !requireNamespace('smacof',quietly=TRUE)) next
    expect_error(metric.mds(prepared=p,backend=backend,pair.weights='inverse_squared'), 'strictly positive')
    fit <- suppressWarnings(metric.mds(prepared=p,backend=backend,init=X,max.iter=5,diagnostics=FALSE))
    expect_equal(as.vector(dist(fit$coords)),as.vector(dist(X)),tolerance=1e-8)
    expect_error(metric.mds(prepared=p,backend=backend,pair.weights='bad'), 'arg')
    p2 <- p; p2$distance_matrix[1,2] <- p2$distance_matrix[2,1] <- 1e-200
    expect_error(metric.mds(prepared=p2,backend=backend,pair.weights='inverse_squared'), 'numeric range')
  }
})
