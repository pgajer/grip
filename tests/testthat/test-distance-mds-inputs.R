test_that("direct dense distances are unchanged and bypass shortest paths", {
  D <- as.matrix(dist(matrix(sin(seq_len(24)), 8, 3)))
  D[1, 2] <- D[2, 1] <- 20 # deliberate triangle inequality violation
  rownames(D) <- colnames(D) <- letters[1:8]
  fit <- suppressWarnings(metric.mds(distance.matrix = D, dim = 3,
    init = "random", max.iter = 3, seed = 11))
  expect_identical(fit$prepared$distance_matrix, D)
  expect_equal(dim(fit$coords), c(8L, 3L))
  expect_equal(rownames(fit$coords), letters[1:8])
  expect_null(fit$diagnostics)
  expect_equal(fit$prepared$graph_build_mode, "supplied_distances")
  expect_error(metric.mds(distance.matrix = D, n = 7), "n must match")
  expect_error(metric.mds(distance.matrix = D, diagnostics = TRUE), "diagnostics")
  expect_error(metric.mds(distance.matrix = D, approximation = "sparse"), "full approximation")
  expect_error(metric.mds(distance.matrix = D, edges = matrix(c(1,2), 1)), "cannot accompany")
  bad <- D; bad[1,2] <- NA
  expect_error(metric.mds(distance.matrix = bad), "finite, symmetric")
  bad <- D; bad[1,2] <- 1
  expect_error(metric.mds(distance.matrix = bad), "finite, symmetric")
  classical <- suppressWarnings(metric.mds(distance.matrix = as.matrix(dist(matrix(1:24,8))),
    dim = 3, max.iter = 2))
  expect_true(all(is.finite(classical$coords)))
})

test_that("landmark adapter deduplicates exact targets with explicit weights", {
  D <- abs(outer(0:4, 0:4, "-"))
  C <- landmark.mds.constraints(D[c(1,5),], c(1,5))
  expect_equal(nrow(C$pairs), 7)
  expect_equal(C$targets, D[C$pairs])
  expect_true(all(C$count_i == 1 & C$count_j == 1))
  # For landmark1, region {1,2,3}; sample5 target4 includes {1,2,3}.
  R <- landmark.mds.constraints(D[c(1,5),], c(1,5), weighting = "region")
  pair <- which(R$pairs[,1] == 1 & R$pairs[,2] == 5)
  expect_equal(R$count_i[pair], 2) # landmark5 region {4,5}
  expect_equal(R$count_j[pair], 3)
  L <- landmark.mds.constraints(D[c(1,5),], c(1,5),
    matrix(c(5,1,2,3), byrow = TRUE, ncol = 2), c(4,1), weighting = "region")
  pair <- which(L$pairs[,1] == 1 & L$pairs[,2] == 5)
  expect_equal(L$count_i[pair], 1)
  expect_equal(L$count_j[pair], 1)
  expect_equal(nrow(L$pairs), 8)
  fit <- metric.mds(n = 5, constraints = L, approximation = "sparse", dim = 3, max.iter = 2)
  expect_true(all(is.finite(fit$coords)))
  expect_error(landmark.mds.constraints(D[c(1,5),], c(1,1)), "distinct")
  expect_error(landmark.mds.constraints(D[c(1,5),], c(1,5), matrix(c(1,5),1), 3), "agree")
  bad <- D[c(1,5),]; bad[1,2] <- 0
  expect_error(landmark.mds.constraints(bad,c(1,5)), "consolidate")
})

test_that("uniform landmark adapter matches the general path and supports 10D refinement", {
 set.seed(27);D<-as.matrix(dist(matrix(rnorm(90),30,3)));landmarks<-c(23L,1L,9L,17L)
 a<-landmark.mds.constraints(D[landmarks,],landmarks)
 b<-landmark.mds.constraints(D[landmarks,],landmarks,matrix(c(1,23),1),D[1,23])
 expect_identical(a,b)
 fit<-metric.mds(n=30L,constraints=a,approximation="sparse",dim=10L,max.iter=3,seed=29)
 expect_equal(dim(fit$coords),c(30L,10L));expect_true(all(is.finite(fit$coords)))
 init<-prcomp(fit$coords,rank.=3)$x[,1:3]
 refined<-metric.mds(n=30L,constraints=a,approximation="sparse",dim=3L,init=init,max.iter=3,seed=29)
 expect_equal(refined$metadata$sparse$targets,fit$metadata$sparse$targets)
 expect_equal(dim(refined$coords),c(30L,3L))
 expect_error(metric.mds(n=30L,constraints=a,approximation="sparse",dim=2.5),"dim")
})
