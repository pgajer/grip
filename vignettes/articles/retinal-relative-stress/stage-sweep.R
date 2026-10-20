# One independently schedulable preparation or fit; all targets remain PC1814 Euclidean.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3L)
mode <- args[1L]; root <- normalizePath(args[2L]); out <- normalizePath(args[3L])
pkgload::load_all(file.path(dirname(root), "dgraphs"), recompile = FALSE, quiet = TRUE)
pkgload::load_all(root, recompile = FALSE, quiet = TRUE)
hash.file <- function(f) digest::digest(file = f, algo = "sha256")
atomic.rds <- function(x, f) {
  saveRDS(x, paste0(f, ".tmp")); stopifnot(file.rename(paste0(f, ".tmp"), f))
}
files <- unlist(lapply(c(root, file.path(dirname(root), "dgraphs")), function(r)
  c(list.files(file.path(r, "R"), "\\.R$", full.names = TRUE),
    list.files(file.path(r, "src"), "\\.(cpp|h|hpp|so)$", full.names = TRUE))))
files <- c(files, file.path(root, "vignettes/articles/retinal-relative-stress",
                            c("stage-sweep.R", "stage-sweep.py")))
fingerprints <- vapply(files, hash.file, "")
expected <- read.csv(file.path(dirname(out), "source-manifest.csv"))
stopifnot(identical(files, expected$path), identical(unname(fingerprints), expected$sha256))
config <- jsonlite::read_json(file.path(out, "input.json"), simplifyVector = TRUE)
stopifnot(hash.file(file.path(out, "input.bin")) == config$input_sha256,
          hash.file(file.path(out, "cells.csv")) == config$cells_sha256)
pair.dist <- function(X, pairs) {
  ans <- numeric(nrow(pairs))
  for (lo in seq.int(1L, nrow(pairs), 2000L)) {
    ix <- lo:min(lo + 1999L, nrow(pairs))
    ans[ix] <- sqrt(rowSums((X[pairs[ix, 1L], , drop = FALSE] -
                             X[pairs[ix, 2L], , drop = FALSE])^2))
  }
  ans
}
score <- function(Y, target) 100 * sqrt(mean((as.numeric(dist(Y)) / target - 1)^2))
if (mode == "prepare") {
  stopifnot(!file.exists(file.path(out, "prepared.rds")))
  con <- file(file.path(out, "input.bin"), "rb")
  v <- readBin(con, "double", 1500L * 1814L, endian = "little"); close(con)
  stopifnot(length(v) == 1500L * 1814L, all(is.finite(v)))
  X <- matrix(v, 1500L, 1814L, byrow = TRUE); rm(v)
  target <- as.numeric(dist(X)); stopifnot(all(target > 0))
  began <- proc.time()[["elapsed"]]
  prep <- grip::prepare.sparse.mds(X, method = "euclidean", n.pivots = 200L,
    pivot.selection = "randomized", q = 10L, neighbor.variance.target = .9,
    seed = config$preparation_seed, max.workspace.bytes = 4 * 1024^3)
  preparation.seconds <- proc.time()[["elapsed"]] - began
  stopifnot(max(abs(pair.dist(X, prep$constraints$pairs) - prep$constraints$targets)) < 1e-8)
  Z <- sweep(X, 2L, prep$pca$center, "-") %*% prep$pca$loadings
  neighbors <- function(A) {
    D <- as.matrix(dist(A)); diag(D) <- Inf
    t(apply(D, 1L, function(row) order(row)[1:10]))
  }
  original <- neighbors(X); reduced <- neighbors(Z)
  retention <- vapply(seq_len(nrow(X)), function(i)
    length(intersect(original[i, ], reduced[i, ])) / 10, 0.0)
  protocol <- list(input = config, source.sha256 = fingerprints, n = 1500L,
    dimensions = c(50L, 25L, 20L, 10L, 5L, 4L, 3L), input.dimension = 1814L,
    pivots = 200L, pivot.selection = "randomized", q = 10L, neighbor.variance.target = .9,
    backend = "sgd", fit.seed = 20190619L, starts.per.dimension = 1L,
    schedule = "zheng", schedule.epochs = 100L, schedule.epsilon = .1, epochs = 100L,
    initialization = "nested normal columns, centered and scale-fitted on all pairs before fitting",
    target = "ambient Euclidean in the same full-population 1814-PC coordinates",
    evaluation = "all 1124250 within-stage unordered pairs; terminal coordinates; no rescaling",
    native.workspace.bytes = 4 * 1024^3, threads.per.worker = 1L)
  key <- digest::digest(protocol, algo = "sha256", serializeVersion = 2)
  atomic.rds(list(protocol = protocol, cache.key = key, X = X, target = target, prep = prep),
             file.path(out, "prepared.rds"))
  diagnostic <- data.frame(stage = config$stage, n = 1500L, pairs = nrow(prep$constraints$pairs),
    search.dimension = prep$pca$dimension, search.variance = prep$pca$achieved,
    neighbor.retention = mean(retention), target.q25 = unname(quantile(target, .25)),
    target.median = median(target), target.q75 = unname(quantile(target, .75)),
    preparation.seconds = preparation.seconds)
  write.csv(diagnostic, file.path(out, "preparation.csv"), row.names = FALSE)
  jsonlite::write_json(list(cache.key = key, protocol = protocol, diagnostic = diagnostic),
    file.path(out, "protocol.json"), pretty = TRUE, auto_unbox = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))
  cat("Prepared", config$stage, "\n")
} else if (mode == "fit") {
  d <- as.integer(args[4L]); shared <- readRDS(file.path(out, "prepared.rds"))
  p <- shared$protocol; z <- shared$prep$constraints
  stopifnot(d %in% p$dimensions, identical(fingerprints, p$source.sha256), identical(config, p$input))
  set.seed(p$fit.seed)
  G <- matrix(rnorm(p$n * 50L), p$n, 50L)
  Y <- scale(G[, seq_len(d), drop = FALSE], center = TRUE, scale = FALSE)
  relative <- as.numeric(dist(Y)) / shared$target
  Y <- Y * sum(relative) / sum(relative^2)
  start <- Y; fit <- NULL; seconds <- 0
  rows <- list(data.frame(stage = config$stage, dimension = d, epoch = 0L,
    rms.relative.percent = score(Y, shared$target), proxy.rms.percent = NA_real_, fit.seconds = 0))
  for (epoch in seq.int(10L, 100L, 10L)) {
    began <- proc.time()[["elapsed"]]
    if (is.null(fit)) fit <- grip::metric.mds(prepared = shared$prep, approximation = "sparse",
      backend = "sgd", dim = d, init = Y, n.init = 1L, seed = p$fit.seed, max.iter = 10L,
      diagnostics = FALSE, pair.weights = "inverse_squared", sgd.control = list(
        scheduler = "zheng", schedule.epsilon = .1, schedule.epochs = 100L,
        checkpoint.every = 10L, keep.state = TRUE, max.workspace.bytes = p$native.workspace.bytes))
    else fit <- grip::metric.mds(resume = fit, max.iter = 10L)
    seconds <- seconds + proc.time()[["elapsed"]] - began
    stopifnot(fit$metadata$engine == "sgd", ncol(fit$coords) == d,
      fit$metadata$stopping$epochs_completed == epoch,
      fit$metadata$settings$seed == p$fit.seed,
      fit$metadata$settings$sgd_control$schedule.epochs == 100L,
      all(is.finite(fit$coords)))
    w <- (z$count_i + z$count_j) / 2
    proxy <- 100 * sqrt(sum(w * (pair.dist(fit$coords, z$pairs) / z$targets - 1)^2) / sum(w))
    stopifnot(abs(proxy - 100 * fit$metadata$sparse_proxy_normalized_rmse) < 1e-8)
    rows[[length(rows) + 1L]] <- data.frame(stage = config$stage, dimension = d, epoch = epoch,
      rms.relative.percent = score(fit$coords, shared$target), proxy.rms.percent = proxy,
      fit.seconds = seconds)
    atomic.rds(list(cache.key = shared$cache.key, dimension = d, start = start,
      fit = fit, scores = do.call(rbind, rows)), file.path(out, sprintf("dimension-%03d.rds", d)))
    write.csv(do.call(rbind, rows), file.path(out, sprintf("scores-%03d.csv", d)), row.names = FALSE)
    cat(config$stage, d, epoch, tail(rows, 1L)[[1L]]$rms.relative.percent, "\n"); flush.console()
  }
} else stop("Unknown task")
stopifnot(identical(fingerprints, vapply(files, hash.file, "")))
