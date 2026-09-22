# A deliberately simple reference: no native helpers or shared scoring code.
sgd_reference <- function(x, targets, rates) {
  pairs <- do.call(rbind, lapply(seq_len(nrow(x) - 1L), function(i)
    cbind(i, seq.int(i + 1L, nrow(x)))))
  history <- list()
  configurations <- list(x)
  evaluate <- function(y) {
    observed <- as.double(dist(y))
    scale <- sum(observed * targets) / sum(observed^2)
    c(raw = sum((observed - targets)^2),
      profiled = sum((scale * observed - targets)^2))
  }
  history[[1]] <- evaluate(x)
  for (epoch in seq_along(rates)) {
    for (k in seq_len(nrow(pairs))) {
      i <- pairs[k, 1]; j <- pairs[k, 2]
      v <- x[i, ] - x[j, ]
      r <- sqrt(sum(v^2))
      stopifnot(r > 0)
      step <- min(rates[epoch], 1) * (r - targets[k]) / 2 * v / r
      x[i, ] <- x[i, ] - step
      x[j, ] <- x[j, ] + step
    }
    configurations[[epoch + 1L]] <- x
    history[[epoch + 1L]] <- evaluate(x)
  }
  history <- do.call(rbind, history)
  best <- which.min(history[, "profiled"])
  list(terminal = x, best = configurations[[best]], history = history, best_epoch = best - 1L)
}

sgd_call <- function(...) suppressWarnings(metric.mds(..., backend = "sgd"))

test_that("fixed-order updates and independently scored checkpoints match R", {
  x <- rbind(c(0, 0, 0), c(2, 1, -1), c(-1, 2, 1), c(1, -2, 3))
  targets <- c(1, 3, 2, 2, 1, 4)
  rates <- c(0.5, 0.2, 0.1)
  before <- x
  native <- grip_sgd_mds_cpp(x, targets, rates, 9L, 1L, 1e6, FALSE)
  reference <- sgd_reference(x, targets, rates)
  expect_equal(x, before)
  expect_equal(native$terminal_conf, reference$terminal, tolerance = 1e-12)
  expect_equal(native$conf, reference$best, tolerance = 1e-12)
  expect_equal(native$trace$raw_stress, reference$history[, "raw"], ignore_attr = TRUE)
  expect_equal(native$trace$profiled_stress, reference$history[, "profiled"], ignore_attr = TRUE)
  expect_equal(native$best_epoch, reference$best_epoch)
  expect_equal(native$pair_updates, 18)
  expect_equal(colMeans(native$terminal_conf), colMeans(x), tolerance = 1e-12)
})

test_that("pair movement clips overshoot and zero targets are not floored", {
  x <- rbind(c(0, 0), c(4, 0), c(0, 3))
  for (targets in list(c(2, 3, 4), c(0, 3, 4))) {
    native <- grip_sgd_mds_cpp(x, targets, 4, 1L, 1L, 1e6, FALSE)
    ref <- sgd_reference(x, targets, 4)
    expect_equal(native$terminal_conf, ref$terminal)
  }
})

test_that("schedule values have explicit boundary and small-budget conventions", {
  control <- grip.mds.sgd.control(list(), 10)
  rate <- grip.mds.sgd.rates(control, 10)
  expect_equal(rate[1], .5)
  expect_equal(rate[5], .05)
  expect_equal(rate[10], .05 / (1 + 4 * 5/6))
  expect_true(all(diff(rate) <= 0))
  control$scheduler <- "exponential"
  expect_equal(grip.mds.sgd.rates(control, 10), .5 * (.01/.5)^((0:9)/10))
  control$scheduler <- "hybrid"
  expect_equal(grip.mds.sgd.rates(control, 1), .5)
  control$switch_ratio <- 0
  expect_equal(grip.mds.sgd.rates(control, 2), c(.5, .5/(1+49/2)))
  control$switch_ratio <- 1
  expect_equal(grip.mds.sgd.rates(control, 2), .5 * 10^(-c(0,1)/2))
  control$final_rate <- .4
  expect_equal(grip.mds.sgd.rates(control, 2), .5 * 10^(-c(0,1)/2))
  control$switch_ratio <- .4
  expect_equal(grip.mds.sgd.rates(control, 10)[5:10], rep(.05,6))
})

test_that("initial and terminal checkpoints are mandatory and best is retained", {
  x <- rbind(c(0, 0), c(1, 0), c(0, 1), c(2, 3))
  target <- as.double(dist(x))
  exact <- grip_sgd_mds_cpp(x, target, rep(.5, 5), 19L, 3L, 1e6)
  expect_equal(exact$trace$epoch, c(0,3,5))
  expect_equal(exact$best_epoch, 0L)
  expect_equal(exact$conf, x)
  expect_equal(exact$stress, 0)
  expect_equal(exact$trace$pair_updates, c(0,18,30))
  # Selection is by recomputed scale-profiled stress, not raw stress.
  expect_equal(exact$stress, min(exact$trace$profiled_stress))
})

test_that("native shuffling is deterministic and independent of R RNG", {
  x <- matrix(seq_len(20), 5, 4)
  targets <- seq_len(10) / 3
  set.seed(991); before <- .Random.seed
  a <- grip_sgd_mds_cpp(x, targets, c(.5,.2), 4L, 1L, 1e6)
  expect_identical(.Random.seed, before)
  b <- grip_sgd_mds_cpp(x, targets, c(.5,.2), 4L, 1L, 1e6)
  c <- grip_sgd_mds_cpp(x, targets, c(.5,.2), 5L, 1L, 1e6)
  expect_equal(a, b)
  expect_false(isTRUE(all.equal(a$terminal_conf, c$terminal_conf)))
})

test_that("positive-target collisions separate reproducibly and stay finite", {
  x <- rbind(c(0,0,0), c(0,0,0), c(1,1,1), c(-1,1,0))
  d <- c(1,2,3,2,1,2)
  a <- grip_sgd_mds_cpp(x, d, .5, 13L, 1L, 1e6, FALSE)
  b <- grip_sgd_mds_cpp(x, d, .5, 13L, 1L, 1e6, FALSE)
  expect_equal(a, b)
  expect_true(all(is.finite(a$terminal_conf)))
  expect_gt(sum((a$terminal_conf[1,]-a$terminal_conf[2,])^2), 0)
  expect_equal(colMeans(a$terminal_conf), colMeans(x), tolerance = 1e-12)
  d[1] <- 0
  zero <- grip_sgd_mds_cpp(x, d, .5, 13L, 1L, 1e6, FALSE)
  expect_true(all(is.finite(zero$terminal_conf)))
})

test_that("near-coincident pairs do not overflow the update ratio", {
  x <- rbind(c(0,0), c(1e-300,0), c(1,1))
  a <- grip_sgd_mds_cpp(x, c(1,2,3), .5, 13L, 1L, 1e6, FALSE)
  expect_true(all(is.finite(a$terminal_conf)))
})

test_that("workspace and native input checks fail before unsafe allocation", {
  x <- matrix(1:12, 4, 3)
  expect_error(grip_sgd_mds_cpp(x, rep(1,6), .5, 1L, 1L, 1), "workspace")
  expect_error(grip_sgd_mds_cpp(x, rep(1,5), .5, 1L, 1L, 1e6), "target count")
  expect_error(grip_sgd_mds_cpp(x, rep(NA_real_,6), .5, 1L, 1L, 1e6), "target")
  expect_error(grip_sgd_mds_cpp(x, rep(-1,6), .5, 1L, 1L, 1e6), "target")
  expect_error(grip_sgd_mds_cpp(x, rep(1,6), 0, 1L, 1L, 1e6), "rates")
  expect_error(grip_sgd_mds_cpp(x * Inf, rep(1,6), .5, 1L, 1L, 1e6), "start")
})

test_that("SGD is the default and controls reject backend ambiguity", {
  args <- list(edges = edges.cycle(5), n = 5, diagnostics = FALSE, max_iter = 2)
  expect_error(do.call(metric.mds, c(args, list(backend="sgd", eps=1e-8))), "SMACOF tolerance")
  expect_error(do.call(metric.mds, c(args, list(backend="smacof", sgd_control=list(learning_rate=.5)))), "requires")
  for (bad in list(list(unknown=1), list(1), list(scheduler="bad"),
                  list(final_rate=1), list(checkpoint_every=1.5),
                  list(learning_rate=NA_real_), list(switch_ratio=2),
                  list(learning_rate=.5, learning_rate=.4), list(final_rate=NULL))) {
    expect_error(do.call(metric.mds, c(args, list(backend="sgd", sgd_control=bad))))
  }
  expect_error(grip.mds.sgd.control(list(max_workspace_bytes=1), 10), "rate vector")
  a <- suppressWarnings(do.call(metric.mds, args))
  b <- suppressWarnings(do.call(metric.mds, c(args, list(backend="sgd"))))
  expect_identical(a$metadata$engine, "sgd")
  expect_equal(a$coords, b$coords)
  a$metadata$starts$elapsed_seconds <- NULL
  b$metadata$starts$elapsed_seconds <- NULL
  expect_equal(a, b)
})

test_that("SGD recovers exact Euclidean distances and improves imperfect starts", {
  truth <- rbind(c(0,0,0), c(1,0,0), c(0,2,0), c(0,0,3), c(1,1,1))
  p <- prepare.graph.geodesic.mds(edges.cycle(5), n=5)
  p$distance_matrix <- as.matrix(dist(truth))
  exact <- sgd_call(prepared=p, dim=3, init=truth, max_iter=30, diagnostics=FALSE)
  expect_equal(as.matrix(dist(exact$coords)), p$distance_matrix, tolerance=1e-10)
  random <- sgd_call(prepared=p, dim=3, init="random", seed=41,
                     max_iter=200, diagnostics=FALSE)
  expect_lt(random$metadata$raw_stress, random$metadata$starts$initial_raw_stress / 10)
  expect_lt(random$metadata$target_normalized_rmse, .02)
  expect_equal(random$metadata$raw_stress,
               sum((as.double(dist(random$coords))-as.double(as.dist(p$distance_matrix)))^2))
  expect_equal(random$metadata$engine, "sgd")
  expect_equal(random$metadata$backend_version, "grip-sgd-mds-v1")
})

test_that("distance units and all supported dimensions are respected", {
  for (dimension in 2:4) {
    a <- sgd_call(edges=edges.cycle(7), n=7, dim=dimension,
                  edge_weights=rep(1,7), init="random", max_iter=10, diagnostics=FALSE)
    b <- sgd_call(edges=edges.cycle(7), n=7, dim=dimension,
                  edge_weights=rep(1e5,7), init="random", max_iter=10, diagnostics=FALSE)
    expect_equal(b$coords / 1e5, a$coords, tolerance=1e-10)
    expect_equal(b$metadata$raw_stress / 1e10, a$metadata$raw_stress, tolerance=1e-10)
    expect_equal(ncol(a$coords), dimension)
  }
})

test_that("rank-deficient supplied starts are not silently perturbed", {
  p <- prepare.graph.geodesic.mds(edges.cycle(6), n=6)
  start <- cbind(cos((1:6)*pi/3), sin((1:6)*pi/3), 0)
  fit <- sgd_call(prepared=p, dim=3, init=start, max_iter=10, diagnostics=FALSE)
  expect_equal(fit$coords[,3], rep(0,6))
})

test_that("SGD multiple starts preserve RNG and keep reproducible histories", {
  set.seed(129); before <- .Random.seed
  a <- sgd_call(edges=edges.cycle(7), n=7, dim=3, init="random", n_init=3,
                seed=5, max_iter=10, diagnostics=FALSE)
  expect_identical(.Random.seed, before)
  b <- sgd_call(edges=edges.cycle(7), n=7, dim=3, init="random", n_init=3,
                seed=5, max_iter=10, diagnostics=FALSE)
  expect_equal(a$coords, b$coords)
  expect_equal(a$metadata$sgd, b$metadata$sgd)
  expect_equal(a$metadata$starts$native_seed, b$metadata$starts$native_seed)
  expect_equal(a$metadata$raw_stress, min(a$metadata$starts$raw_stress))
  expect_true(all(a$metadata$starts$raw_stress <= a$metadata$starts$initial_raw_stress + 1e-10))
  expect_equal(a$metadata$starts$pair_updates, rep(210,3))
})

test_that("seed scope includes supplied starts, absent RNG, and errors", {
  withr::local_preserve_seed()
  if (exists(".Random.seed", envir=.GlobalEnv)) rm(".Random.seed", envir=.GlobalEnv)
  p <- prepare.graph.geodesic.mds(edges.cycle(5), n=5)
  start <- cbind(1:5, c(0,1,0,-1,2))
  fit <- sgd_call(prepared=p, init=start, max_iter=3, diagnostics=FALSE)
  expect_false(exists(".Random.seed", envir=.GlobalEnv))
  expect_error(sgd_call(prepared=p, init=start, max_iter=3, diagnostics=FALSE,
                       sgd_control=list(max_workspace_bytes=30)), "All SGD starts failed")
  expect_false(exists(".Random.seed", envir=.GlobalEnv))
  set.seed(88); before <- .Random.seed
  sgd_call(prepared=p, init=start, max_iter=3, seed=NULL, diagnostics=FALSE)
  expect_false(identical(.Random.seed,before))
})

test_that("SGD works without SMACOF and reports schedule completion honestly", {
  testthat::local_mocked_bindings(grip.mds.has.smacof=function() FALSE)
  expect_warning(fit <- metric.mds(edges=edges.cycle(5), n=5,
                                   max_iter=1, diagnostics=FALSE), "iteration_limit")
  expect_false(fit$metadata$converged)
  expect_identical(fit$metadata$termination, "iteration_limit")
  expect_equal(fit$metadata$sgd[[1]]$trace$epoch, c(0,1))
  expect_error(metric.mds(edges=edges.cycle(5), n=5, backend="smacof"), "optional 'smacof'")
})

test_that("one failed SGD start does not erase other starts or hide failure", {
  original <- grip.mds.sgd.fit
  count <- 0L
  testthat::local_mocked_bindings(grip.mds.sgd.fit=function(...) {
    count <<- count+1L
    if (count==1L) stop("controlled SGD error")
    original(...)
  })
  a <- sgd_call(edges=edges.cycle(6), n=6, n_init=2, max_iter=4, diagnostics=FALSE)
  expect_equal(a$metadata$selected_start, 2)
  expect_equal(a$metadata$starts$termination[1], "backend_error")
  expect_match(a$metadata$starts$error[1], "controlled SGD error")
  expect_match(a$metadata$sgd[[1]]$error, "controlled SGD error")
  count <- 0L
  expect_error(sgd_call(edges=edges.cycle(6), n=6, max_iter=4, diagnostics=FALSE),
               "All SGD starts failed")
})

test_that("interrupts propagate rather than masquerading as numerical failures", {
  testthat::local_mocked_bindings(grip.mds.sgd.fit=function(...) {
    signalCondition(structure(list(message="user interrupt",call=NULL),
                              class=c("interrupt","condition")))
    stop("interrupt was swallowed")
  })
  caught <- tryCatch(sgd_call(edges=edges.cycle(6),n=6,max_iter=2,diagnostics=FALSE),
                    interrupt=function(e) "interrupted")
  expect_identical(caught, "interrupted")
})

test_that("SGD preserves diagnostics and integrates by explicit edge-KK coordinates", {
  p <- prepare.graph.geodesic.mds(edges.cycle(6),n=6)
  a <- sgd_call(prepared=p,dim=3,max_iter=5)
  expect_s3_class(a,"grip_gmds_layout")
  expect_true(is.data.frame(a$diagnostics))
  refined <- edge.kk(prepared=p,dim=3,coords=a$coords,max_iter=0,diagnostics=FALSE)
  expect_equal(as.matrix(dist(refined$coords)),as.matrix(dist(a$coords)),tolerance=1e-8)
})

test_that("SGD uses common invalid-input and duplicate-observation policies", {
  p <- prepare.graph.geodesic.mds(edges.cycle(5),n=5)
  for (bad in c(NA_real_,Inf,-1)) {
    q <- p; q$distance_matrix[1,2] <- q$distance_matrix[2,1] <- bad
    expect_error(sgd_call(prepared=q,max_iter=2), "finite, symmetric")
  }
  expect_error(sgd_call(prepared=p,init=matrix(0,5,2),max_iter=2), "collapsed")
  expect_error(sgd_call(edges=rbind(c(1,2),c(3,4)),n=4,max_iter=2,diagnostics=FALSE), "connected")
  x <- rbind(c(0,0),c(0,0),c(1,0),c(0,1),c(1,1))
  p$distance_matrix <- as.matrix(dist(x))
  fit <- sgd_call(prepared=p,init=x,max_iter=5,diagnostics=FALSE)
  expect_equal(as.matrix(dist(fit$coords)),p$distance_matrix,tolerance=1e-10)
})
