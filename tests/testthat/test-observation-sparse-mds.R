observation_fixture <- function() {
  u <- matrix(sin(seq_len(96)*1.37),32,3)
  cbind(u,.8*(u[,1]^2-u[,2]^2+u[,3]^2))
}

test_that("on-demand distances and region counts match independent dense calculations", {
  skip_if_not_installed("dgraphs")
  X <- observation_fixture(); D <- as.matrix(dist(X))
  for (method in c("uniform","euclidean","geodesic")) {
    a <- prepare.sparse.mds(X,method=method,n.pivots=7,q=5,seed=19)
    if (method=="geodesic") {
      skip_if_not_installed("igraph")
      ed <- dgraphs::graph.edges(a$graph)
      g <- igraph::graph_from_data_frame(ed[,1:2],directed=FALSE,vertices=seq_len(nrow(X)))
      expected <- igraph::distances(g,weights=ed$length)
    } else expected <- D
    expect_equal(unname(a$distance_rows),unname(expected[a$pivots,,drop=FALSE]),tolerance=1e-12)
    owner <- a$pivots[max.col(-t(expected[a$pivots,,drop=FALSE]),ties.method="first")]
    expect_equal(a$region,owner)
    C <- a$constraints; ed <- if(method=="uniform") NULL else dgraphs::graph.edges(a$graph)
    local <- if(is.null(ed)) character() else paste(ed$from,ed$to)
    for (row in seq_len(nrow(C$pairs))) {
      i <- C$pairs[row,1]; j <- C$pairs[row,2]
      if (method=="uniform" || paste(i,j) %in% local) { ai<-1; bj<-1 } else {
        ai <- if(j %in% a$pivots) sum(expected[j,owner==j] <= expected[j,i]/2) else 0
        bj <- if(i %in% a$pivots) sum(expected[i,owner==i] <= expected[i,j]/2) else 0
      }
      expect_equal(c(C$count_i[row],C$count_j[row]),c(ai,bj))
      expect_equal(C$targets[row],expected[i,j],tolerance=1e-12)
    }
    expect_false(a$metadata$dense_distance_matrix_allocated)
  }
})

test_that("selection is reproducible, spread-out, unit invariant and RNG preserving", {
  skip_if_not_installed("dgraphs")
  X <- observation_fixture()
  set.seed(122); state <- .Random.seed
  a <- prepare.sparse.mds(X,n.pivots=6,q=5,seed=13,pivot.selection="farthest")
  expect_identical(.Random.seed,state)
  b <- prepare.sparse.mds(X,n.pivots=6,q=5,seed=13,pivot.selection="farthest")
  expect_identical(a$pivots,b$pivots); expect_equal(a$constraints,b$constraints)
  D <- as.matrix(dist(X))
  for (j in 2:6) {
    nearest <- apply(D[a$pivots[seq_len(j-1)],,drop=FALSE],2,min)
    expect_equal(a$pivots[j],unname(which.max(nearest)))
  }
  for (method in c("uniform","euclidean","geodesic")) {
    a <- prepare.sparse.mds(X,method=method,n.pivots=6,q=5,seed=13)
    b <- prepare.sparse.mds(X*100,method=method,n.pivots=6,q=5,seed=13)
    expect_identical(a$pivots,b$pivots)
    expect_equal(a$constraints$count_i,b$constraints$count_i)
    expect_equal(a$constraints$count_j,b$constraints$count_j)
    expect_equal(a$constraints$targets,b$constraints$targets/100,tolerance=1e-11)
  }
})

test_that("PCA selection is minimal and local targets retain original geometry", {
  skip_if_not_installed("dgraphs")
  X <- observation_fixture()
  a <- prepare.sparse.mds(X,n.pivots=6,q=5,neighbor.variance.target=.9)
  cv <- cumsum(a$pca$variance)/sum(a$pca$variance); r <- a$pca$dimension
  expect_gte(cv[r],.9); if(r>1) expect_lt(cv[r-1],.9)
  U <- sweep(X,2,a$pca$center) %*% a$pca$loadings
  expect_equal(sum(as.numeric(dist(U))^2)/sum(as.numeric(dist(X))^2),cv[r],tolerance=1e-12)
  E <- dgraphs::graph.edges(a$graph)
  expect_equal(E$length,sqrt(rowSums((X[E$from,]-X[E$to,])^2)))
})

test_that("all-pivot limit and exact-coordinate stress are correct", {
  X <- observation_fixture(); n <- nrow(X)
  C <- grip_sparse_prepare_cpp(n,matrix(integer(),ncol=2),numeric(),n,seq_len(n),1L,
    1e8,observations=X,save_distances=FALSE)
  expect_equal(nrow(C$pairs),choose(n,2))
  expect_true(all(C$count_i==1 & C$count_j==1))
  D <- as.matrix(dist(X)); expect_equal(C$targets,D[C$pairs],tolerance=1e-12)
  fit <- metric.mds(n=n,constraints=C[c("pairs","targets","count_i","count_j")],
    approximation="sparse",dim=4,init=X,max.iter=3,sgd.control=list(scheduler="zheng"))
  expect_lt(fit$metadata$sparse_proxy_stress,1e-20)
})

test_that("published schedule includes endpoints and handles zero influence", {
  control <- grip.mds.sgd.control(list(scheduler="zheng",schedule.epsilon=.1),5)
  rates <- grip.mds.sgd.rates(control,5,c(0,2,4,8))
  expect_equal(rates[c(1,5)],c(.5,.1/8))
  expect_equal(diff(log(rates)),rep(diff(range(log(rates)))/-4,4))
  expect_equal(grip.mds.sgd.rates(control,1,c(0,2,8)),.5)
  expect_error(grip.mds.sgd.control(list(scheduler="zheng",learning.rate=1),5),"omit")
})

test_that("observation preparations fit both dimensions and refuse invalid inputs early", {
  X <- observation_fixture()
  a <- prepare.sparse.mds(X,method="uniform",n.pivots=6)
  for (k in 2:3) {
    f <- metric.mds(prepared=a,approximation="sparse",dim=k,max.iter=5,seed=14)
    expect_equal(dim(f$coords),c(nrow(X),k))
    expect_true(all(is.finite(f$coords))); expect_false(f$metadata$converged)
    expect_equal(f$metadata$settings$sgd_control$scheduler,"zheng")
    expect_equal(f$metadata$sparse$pivots,a$pivots)
  }
  expect_error(prepare.sparse.mds(X,max.workspace.bytes=1),"before allocation")
  expect_error(prepare.sparse.mds(rbind(X,X[1,])),"duplicate")
  expect_error(prepare.sparse.mds(X,q=nrow(X)),"q must")
  expect_error(prepare.sparse.mds(X,neighbor.variance.target=0),"neighbor.variance.target")
  expect_error(prepare.sparse.mds(X,method="geodesic",pivot.selection="farthest"),"only")
  expect_error(metric.mds(prepared=a),"requires approximation")
  expect_error(metric.mds(prepared=a,approximation="sparse",distance.matrix=as.matrix(dist(X))),"other data")
})

test_that("PCA search handles ambient lifting and collapsed projected pairs", {
  skip_if_not_installed("dgraphs")
  set.seed(88)
  U <- matrix(runif(80,-1,1),40,2)
  X <- cbind(U,.2*rowSums(U^2))
  Q <- qr.Q(qr(matrix(rnorm(150),50,3)))
  lifted <- X %*% t(Q)
  a <- prepare.sparse.mds(lifted,n.pivots=8,q=5)
  expect_lt(a$pca$dimension,ncol(lifted))
  expect_equal(as.numeric(dist(lifted)),as.numeric(dist(X)),tolerance=1e-12)
  expect_true(all(a$constraints$targets>0))
  # Even coincident search rows retain separate original vertices and lengths.
  g <- dgraphs::create.sknn.graph(X,3,search.data=matrix(0,40,1),neighbor.method="ann",graph.detail="minimal")
  expect_true(all(dgraphs::graph.edges(g)$length>0))
  expect_false(g$metadata$dense_distance_matrix_allocated)
})
