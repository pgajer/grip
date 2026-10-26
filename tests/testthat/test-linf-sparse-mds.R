linf_fixture <- function() as.matrix(expand.grid(x=0:3, y=0:2))

test_that("L-infinity neighbors, ties and regional constraints match dense oracles", {
  X <- linf_fixture(); n <- nrow(X); D <- as.matrix(dist(X, method="maximum"))
  for (q in c(1L, 3L, n-1L)) {
    prep <- prepare.sparse.mds(X, method="linf", n.pivots=4, q=q, seed=19)
    neighbors <- t(vapply(seq_len(n), function(i) {
      ids <- setdiff(seq_len(n),i); ids[order(D[i,ids],ids)][seq_len(q)]
    }, integer(q)))
    a <- rep(seq_len(n),each=q); b <- as.vector(t(neighbors))
    expected <- unique(data.frame(from=pmin(a,b),to=pmax(a,b)))
    expected <- expected[order(expected$from,expected$to),]; rownames(expected) <- NULL
    edges <- prep$graph$edges
    expect_equal(edges[,1:2],expected)
    expect_equal(edges$length,D[as.matrix(expected)])
    expect_equal(prep$distance_rows,unname(D[prep$pivots,,drop=FALSE]))
    owner <- prep$pivots[max.col(-t(D[prep$pivots,,drop=FALSE]), ties.method="first")]
    expect_identical(prep$region,owner)
    oracle <- landmark.mds.constraints(D[prep$pivots,,drop=FALSE],prep$pivots,
      as.matrix(edges[,1:2]),edges$length,weighting="region")
    expect_equal(prep$constraints,oracle[c("pairs","targets","count_i","count_j")])
    expect_null(prep$pca); expect_null(prep$metadata$neighbor_variance_target)
    expect_identical(prep$metadata$target_metric,"linf")
    expect_false(prep$metadata$dense_distance_matrix_allocated)
  }
})

test_that("L-infinity selection is reproducible, scale invariant and RNG preserving", {
  X <- linf_fixture(); D <- as.matrix(dist(X,method="maximum"))
  set.seed(2026); state <- .Random.seed
  for (selection in c("randomized","farthest")) {
    args <- list(X=X,method="linf",n.pivots=5,q=3,seed=42,pivot.selection=selection)
    a <- do.call(prepare.sparse.mds,args); b <- do.call(prepare.sparse.mds,args)
    expect_identical(a$pivots,b$pivots); expect_identical(a$constraints,b$constraints)
    args$X <- X*100
    scaled <- do.call(prepare.sparse.mds,args)
    expect_identical(a$pivots,scaled$pivots)
    expect_identical(a$constraints$count_i,scaled$constraints$count_i)
    expect_identical(a$constraints$count_j,scaled$constraints$count_j)
    expect_equal(a$constraints$targets,scaled$constraints$targets/100)
    if(selection=="farthest") for(k in 2:5) {
      nearest <- apply(D[a$pivots[seq_len(k-1)],,drop=FALSE],2,min)
      expect_equal(a$pivots[k],unname(which.max(nearest)))
    }
  }
  expect_identical(.Random.seed,state)
})

test_that("L-infinity preparation integrates with SGD without changing constraints", {
  X <- linf_fixture()
  prep <- prepare.sparse.mds(X,method="linf",n.pivots=5,q=3)
  a <- metric.mds(prepared=prep,approximation="sparse",backend="sgd",
    dim=3,max.iter=10,seed=17)
  b <- metric.mds(n=nrow(X),constraints=prep$constraints,approximation="sparse",
    backend="sgd",dim=3,max.iter=10,seed=17,sgd.control=list(scheduler="zheng"))
  expect_identical(a$coords,b$coords)
  expect_identical(a$metadata$data_preparation$target_metric,"linf")
  expect_true(all(is.finite(a$coords)))
  expect_false(a$metadata$converged)
})

test_that("L-infinity rejects inappropriate PCA, invalid input and insufficient workspace", {
  X <- linf_fixture()
  expect_error(prepare.sparse.mds(X,method="linf",neighbor.variance.target=.9),"PCA changes")
  expect_error(prepare.sparse.mds(X,method="linf",q=nrow(X)),"q must")
  expect_error(prepare.sparse.mds(rbind(X,X[1,]),method="linf"),"duplicate")
  expect_error(prepare.sparse.mds(X,method="linf",max.workspace.bytes=1),"before allocation")
  expect_error(grip_linf_knn_cpp(X,3L,1),"workspace")
  bad <- X; bad[1,1] <- Inf
  expect_error(prepare.sparse.mds(bad,method="linf"),"finite")
  expect_error(grip_linf_knn_cpp(bad,3L,1e8),"Nonfinite")
  extreme <- matrix(c(-1e308,0,1e308),ncol=1)
  expect_error(prepare.sparse.mds(extreme,method="linf",q=2),"numeric range")
})
