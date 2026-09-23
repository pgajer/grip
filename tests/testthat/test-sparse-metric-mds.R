# Independent dense reference used only on tiny test graphs.
sparse_reference <- function(n, E, lengths, pivots) {
  D <- matrix(Inf,n,n); diag(D) <- 0
  for (k in seq_len(nrow(E))) D[E[k,1],E[k,2]] <- D[E[k,2],E[k,1]] <- lengths[k]
  for (k in seq_len(n)) for (i in seq_len(n)) for (j in seq_len(n))
    D[i,j] <- min(D[i,j],D[i,k]+D[k,j])
  region <- pivots[apply(D[,pivots,drop=FALSE],1,which.min)]
  pairs <- t(combn(n,2)); rows <- list()
  for (k in seq_len(nrow(pairs))) {
    i <- pairs[k,1]; j <- pairs[k,2]
    edge <- which((E[,1]==i & E[,2]==j) | (E[,1]==j & E[,2]==i))
    ci <- if (j %in% pivots) sum(D[j,region==j] <= D[i,j]/2) else 0
    cj <- if (i %in% pivots) sum(D[i,region==i] <= D[i,j]/2) else 0
    if (length(edge)) { ci <- cj <- 1; d <- lengths[edge] } else d <- D[i,j]
    if (ci+cj) rows[[length(rows)+1L]] <- c(i,j,d,ci,cj)
  }
  list(constraints=do.call(rbind,rows),region=region,D=D)
}

test_that('sparse preparation agrees with independent regions, counts, and edge targets', {
  examples <- list(list(n=7,E=edges.path(7),w=c(1,2,1,3,2,1),p=c(6L,2L)),
                   list(n=4,E=edges.cycle(4),w=rep(1,4),p=c(3L,1L)),
                   list(n=3,E=rbind(c(1,2),c(1,3),c(2,3)),w=c(1,9,1),p=1:3))
  for (a in examples) {
    ref <- sparse_reference(a$n,a$E,a$w,a$p)
    x <- grip_sparse_prepare_cpp(a$n,a$E,a$w,length(a$p),a$p,1L,1e7)
    expect_equal(cbind(x$pairs,x$targets,x$count_i,x$count_j),ref$constraints)
    expect_equal(x$region,ref$region)
    expect_equal(x$pivots,a$p)
  }
  # Earliest selected pivot owns exact equidistance ties.
  expect_equal(grip_sparse_prepare_cpp(4L,edges.cycle(4),rep(1,4),2L,c(3L,1L),1L,1e7)$region,
               c(1L,3L,3L,3L))
})

test_that('asymmetric updates preserve an inactive endpoint and match hand calculation', {
  X <- rbind(c(0,0,1),c(4,2,3),c(1,-2,4))
  E <- rbind(c(1L,2L),c(2L,3L)); d <- c(2,3); wi <- c(0,3); wj <- c(.2,1)
  rates <- c(.5,.1); Z <- X; losses <- numeric(3)
  proxy <- function(z) sum((wi+wj)/2*(sqrt(rowSums((z[E[,1],]-z[E[,2],])^2))-d)^2)
  losses[1] <- proxy(Z)
  for (e in seq_along(rates)) {
    for (k in 1:2) {
      i <- E[k,1]; j <- E[k,2]; v <- Z[i,]-Z[j,]; r <- sqrt(sum(v^2))
      change <- (r-d[k])/2*v/r
      Z[i,] <- Z[i,]-min(1,rates[e]*wi[k])*change
      Z[j,] <- Z[j,]+min(1,rates[e]*wj[k])*change
    }
    losses[e+1] <- proxy(Z)
  }
  fit <- grip_sgd_mds_cpp(X,d,rates,9L,1L,1e7,FALSE,wi,E,wj,FALSE)
  expect_equal(fit$terminal_conf,Z,tolerance=1e-12)
  expect_equal(fit$conf,Z,tolerance=1e-12)
  expect_equal(fit$terminal_conf[1,],X[1,],tolerance=0)
  expect_equal(fit$trace$raw_stress,losses)
  # A collision must not move the zero-weight endpoint.
  Y <- rbind(c(0,0,0),c(0,0,0),c(1,1,1))
  coll <- grip_sgd_mds_cpp(Y,c(1,1),1,9L,1L,1e7,FALSE,c(0,1),
                          rbind(c(1L,2L),c(2L,3L)),c(1,1),FALSE)
  expect_equal(coll$terminal_conf[1,],Y[1,],tolerance=0)
  expect_true(all(is.finite(coll$terminal_conf)))
})

test_that('all-pivot updates reduce to full updates with matched order, scale and schedule', {
  n <- 6L; E <- edges.cycle(n); w <- rep(1,n)
  for (dim in 2:3) {
    X <- matrix(sin(seq_len(n*dim)),n,dim)
    x <- grip_sparse_prepare_cpp(n,E,w,n,seq_len(n),1L,1e7)
    target <- x$targets; target <- target/sqrt(mean(target^2)); weights <- 1/target^2
    full <- grip_sgd_mds_cpp(X,target,c(.5,.2,.1),10L,1L,1e7,TRUE,weights)
    sparse <- grip_sgd_mds_cpp(X,target,c(.5,.2,.1),10L,1L,1e7,TRUE,
                              x$count_i*weights,x$pairs,x$count_j*weights,FALSE)
    expect_equal(sparse$terminal_conf,full$terminal_conf,tolerance=0)
    expect_equal(sparse$trace$raw_stress,full$trace$raw_stress,tolerance=0)
  }
})

test_that('public sparse fits are repeatable, scale equivariant, and preserve RNG', {
  E <- edges.cycle(12); X <- cbind(sin(1:12),cos(1:12),sin((1:12)*2))
  set.seed(17); saved <- .Random.seed
  a <- sparse.metric.mds(E,edge_weights=sqrt(1:12)+.137,pivots=c(1,5,9),init=X,max_iter=8)
  expect_identical(.Random.seed,saved)
  expect_identical(a$metadata$sparse$pivots,c(1L,5L,9L))
  b <- sparse.metric.mds(E,edge_weights=sqrt(1:12)+.137,pivots=c(1,5,9),init=X,max_iter=8)
  expect_equal(a$coords,b$coords,tolerance=0)
  for (scale in c(1e-150,1e150)) {
    scaled <- sparse.metric.mds(E,edge_weights=scale*(sqrt(1:12)+.137),pivots=c(1,5,9),init=X*scale,max_iter=8)
    expect_equal(scaled$coords/scale,a$coords,tolerance=1e-10)
    expect_equal(scaled$metadata$sparse_proxy_stress,a$metadata$sparse_proxy_stress,tolerance=1e-10)
  }
  expect_null(a$prepared$distance_matrix)
  expect_equal(layout.coords(a),a$coords)
  expect_equal(colMeans(a$coords),rep(0,3),tolerance=1e-12)
  s <- a$metadata$sparse
  d <- sqrt(rowSums((a$coords[s$pairs[,1],]-a$coords[s$pairs[,2],])^2))
  expected <- sum((s$count_i+s$count_j)/2*((d-s$targets)/s$targets)^2)
  expect_equal(a$metadata$sparse_proxy_stress,expected)
  expect_null(a$metadata$raw_stress)
  expect_false(a$metadata$converged)
  expect_identical(a$metadata$termination,'iteration_limit')
})

test_that('sparse API never invokes dense preparation or classical initialization', {
  local_mocked_bindings(prepare.graph.geodesic.mds=function(...) stop('DENSE'),
    grip.metric.mds.distance.prepared=function(...) stop('DENSE'),
    classical.mds=function(...) stop('DENSE'))
  for (h in c(1L,3L,20L)) {
    fit <- sparse.metric.mds(edges.path(20),n_pivots=h,max_iter=2)
    expect_equal(length(fit$metadata$sparse$pivots),h)
    expect_lte(nrow(fit$metadata$sparse$pairs),20*h+19)
    expect_null(fit$prepared$distance_matrix)
  }
  for (dim in 2:3) expect_equal(ncol(sparse.metric.mds(edges.path(5),dim=dim,max_iter=2)$coords),dim)
})

test_that('sparse invalid inputs and resource limits fail explicitly', {
  E <- edges.path(5)
  expect_error(sparse.metric.mds(E,n=6),'connected')
  expect_error(sparse.metric.mds(E,edge_weights=c(0,1,1,1)),'positive|> 0')
  expect_error(sparse.metric.mds(E,pivots=c(1,1)),'distinct')
  expect_error(sparse.metric.mds(E,pivots=c(1,6)),'vertex')
  expect_error(sparse.metric.mds(E,pivots=c(1,3),n_pivots=3),'match')
  expect_error(sparse.metric.mds(E,dim=4),'dim')
  expect_error(sparse.metric.mds(E,max_iter=0),'positive')
  expect_error(sparse.metric.mds(E,init='classical'),'init')
  expect_true(all(is.finite(sparse.metric.mds(E,init=matrix(0,5,3),max_iter=2)$coords)))
  expect_error(sparse.metric.mds(E,seed=NA),'seed')
  expect_error(sparse.metric.mds(E,sgd_control=list(max_workspace_bytes=1000)),'workspace')
  expect_error(grip_sparse_prepare_cpp(5L,E,rep(1,4),2L,c(1L,1L),1L,1e7),'distinct')
  expect_error(grip_sparse_prepare_cpp(5L,E,rep(1,4),2L,integer(),1L,1),'workspace')
})


test_that('unit scaling does not alter exact half-distance boundaries on a path', {
  E <- edges.path(12)
  a <- sparse.metric.mds(E,n_pivots=3,seed=125,max_iter=3)
  for (scale in c(1e-6,1e-150,1e150)) {
    b <- sparse.metric.mds(E,edge_weights=rep(scale,11),n_pivots=3,seed=125,max_iter=3)
    expect_identical(a$metadata$sparse$count_i,b$metadata$sparse$count_i)
    expect_identical(a$metadata$sparse$count_j,b$metadata$sparse$count_j)
    expect_equal(a$coords,b$coords/scale,tolerance=1e-12)
  }
})

test_that('adjacency input agrees and printed sparse scores are labeled', {
  E <- edges.path(8); prep <- prepare.edge.kk(E,n=8,edge_weights=sqrt(1:7))
  a <- sparse.metric.mds(E,n=8,edge_weights=sqrt(1:7),n_pivots=3,max_iter=3)
  b <- sparse.metric.mds(adj_list=prep$adj_list,weight_list=prep$weight_list,n_pivots=3,max_iter=3)
  expect_equal(a$coords,b$coords,tolerance=0)
  expect_output(print(a),'sparse proxy stress [(]not full stress[)]')
})

test_that('independent scale diagnostic retains coordinates and original errors', {
  e <- new.env(parent=globalenv())
  sys.source(system.file('scripts','sparse-mds-comparison.R',package='grip'),envir=e)
  X <- rbind(c(0,0,0),c(1,0,0),c(0,1,0))
  b <- list(evaluation=list(example=list(sources=1:3,targets=as.matrix(dist(X)))),
    fits=list(fit=list(coords=X*4,row=data.frame(case='example',status='ok',relative_error=3))))
  adjusted <- e$sparse_comparison_profile(b)
  expect_equal(adjusted$fits$fit$coords,b$fits$fit$coords)
  expect_equal(adjusted$fits$fit$row$relative_error,3)
  expect_equal(adjusted$fits$fit$row$relative_scale,.25)
  expect_equal(adjusted$fits$fit$row$profiled_relative_error,0)
})
