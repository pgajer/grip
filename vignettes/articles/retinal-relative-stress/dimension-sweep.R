# Small-subset dimensionality experiment; never resumes the full retinal run.
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1L]; root <- normalizePath(args[2L]); out <- normalizePath(args[3L])
hash.file <- function(f) digest::digest(file = f, algo = "sha256")
atomic.rds <- function(x, f) {
  saveRDS(x, paste0(f, ".tmp")); stopifnot(file.rename(paste0(f, ".tmp"), f))
}
if (mode == "export") {
  prior <- readRDS(args[4L])
  stopifnot(nrow(prior$X) == 1500L, length(prior$cells) == 1500L)
  write.csv(data.frame(index = prior$indices, cell = prior$cells),
            file.path(out, "subset.csv"), row.names = FALSE)
  quit(save = "no")
}
pkgload::load_all(file.path(dirname(root), "dgraphs"), recompile = FALSE, quiet = TRUE)
pkgload::load_all(root, recompile = FALSE, quiet = TRUE)
source.files <- unlist(lapply(c(root, file.path(dirname(root), "dgraphs")), function(r)
  c(list.files(file.path(r, "R"), "\\.R$", full.names = TRUE),
    list.files(file.path(r, "src"), "\\.(cpp|h|hpp|so)$", full.names = TRUE))))
source.files <- c(source.files, file.path(root, "vignettes/articles/retinal-relative-stress",
                                         c("dimension-sweep.R", "dimension-sweep.py")))
fingerprints <- vapply(source.files, hash.file, "")
score <- function(Y, target) 100 * sqrt(mean((as.numeric(dist(Y)) / target - 1)^2))
pair.dist <- function(Y, pairs) {
  ans <- numeric(nrow(pairs))
  for (lo in seq.int(1L, nrow(pairs), 2000L)) {
    ix <- lo:min(nrow(pairs), lo + 1999L)
    ans[ix] <- sqrt(rowSums((Y[pairs[ix, 1L], , drop = FALSE] -
                             Y[pairs[ix, 2L], , drop = FALSE])^2))
  }
  ans
}
if (mode == "prepare") {
  stopifnot(!file.exists(file.path(out, "prepared.rds")))
  con <- file(file.path(out, "subset-pc1814.bin"), "rb")
  X <- matrix(readBin(con, "double", 1500L * 1814L, endian = "little"),
              nrow = 1500L, ncol = 1814L, byrow = TRUE)
  close(con)
  stopifnot(all(is.finite(X)))
  target <- as.numeric(dist(X))
  stopifnot(length(target) == choose(1500L, 2L), all(target > 0))
  began <- proc.time()[["elapsed"]]
  prep <- grip::prepare.sparse.mds(X, method = "euclidean", n.pivots = 200L,
    pivot.selection = "randomized", q = 10L, neighbor.variance.target = .9,
    seed = 20261021L, max.workspace.bytes = 32 * 1024^3)
  preparation.seconds <- proc.time()[["elapsed"]] - began
  stopifnot(max(abs(pair.dist(X, prep$constraints$pairs) - prep$constraints$targets)) < 1e-8)
  dimensions <- c(300L, 200L, 100L, 50L, 25L, 20L, 10L, 5L, 4L, 3L)
  # This subset PCA is a reference, distinct from the full-population input PCA.
  pca <- prcomp(X, center = TRUE, scale. = FALSE, rank. = max(dimensions))
  pca.rows <- lapply(dimensions, function(d) data.frame(dimension = d,
    rms.relative.percent = score(pca$x[, seq_len(d), drop = FALSE], target),
    subset.variance.percent = 100 * sum(pca$sdev[seq_len(d)]^2) / sum(pca$sdev^2)))
  write.csv(do.call(rbind, pca.rows), file.path(out, "pca-reference.csv"), row.names = FALSE)
  search <- sweep(X, 2L, prep$pca$center, "-") %*% prep$pca$loadings
  neighbor.ids <- function(A) {
    D <- as.matrix(dist(A)); diag(D) <- Inf
    t(apply(D, 1L, function(x) order(x)[seq_len(10L)]))
  }
  original <- neighbor.ids(X); reduced <- neighbor.ids(search)
  retention <- vapply(seq_len(nrow(X)), function(i)
    length(intersect(original[i, ], reduced[i, ])) / 10, 0.0)
  protocol <- list(n = nrow(X), input.dimension = ncol(X), dimensions = dimensions,
    backend = "sgd", epochs = 100L, schedule.epochs = 100L,
    schedule = "zheng", schedule.epsilon = .1, fit.seed = 20190619L,
    starts.per.dimension = 1L, preparation.seed = 20261021L, pivots = 200L,
    pivot.selection = "randomized", q = 10L, neighbor.variance.target = .9,
    target = "ambient Euclidean distances in the original 1814-PC input",
    evaluation = "all 1124250 unordered pairs; terminal coordinates; no post-fit rescaling",
    initialization = "shared nested normal coordinates, centered and scale-fitted once per dimension",
    source.sha256 = fingerprints, input.sha256 = hash.file(file.path(out, "subset-pc1814.bin")),
    subset.sha256 = hash.file(file.path(out, "subset.csv")),
    workers = 10L, threads.per.worker = 1L, total.rss.limit.bytes = 50 * 1024^3,
    native.workspace.bytes.per.fit = 4 * 1024^3)
  key <- digest::digest(protocol, algo = "sha256", serializeVersion = 2)
  saveRDS(list(protocol = protocol, cache.key = key, X = X, target = target, prep = prep),
          file.path(out, "prepared.rds"))
  jsonlite::write_json(list(cache.key = key, protocol = protocol,
    preparation.seconds = preparation.seconds, constraints = nrow(prep$constraints$pairs),
    search.dimension = prep$pca$dimension, search.variance = prep$pca$achieved,
    neighbor.retention.mean = mean(retention), neighbor.retention.quantiles = quantile(retention)),
    file.path(out, "protocol.json"), pretty = TRUE, auto_unbox = TRUE)
  write.csv(data.frame(path = source.files, sha256 = fingerprints),
            file.path(out, "source-manifest.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))
  cat("Shared preparation complete\n")
} else if (mode == "fit") {
  d <- as.integer(args[4L]); shared <- readRDS(file.path(out, "prepared.rds"))
  p <- shared$protocol; z <- shared$prep$constraints
  stopifnot(d %in% p$dimensions, identical(fingerprints, p$source.sha256),
    hash.file(file.path(out, "subset-pc1814.bin")) == p$input.sha256)
  set.seed(p$fit.seed)
  G <- matrix(rnorm(p$n * max(p$dimensions)), p$n, max(p$dimensions))
  Y <- scale(G[, seq_len(d), drop = FALSE], center = TRUE, scale = FALSE)
  r <- as.numeric(dist(Y)) / shared$target
  Y <- Y * sum(r) / sum(r^2)
  start <- Y; fit <- NULL; completed <- 0L; elapsed <- 0
  rows <- list(data.frame(dimension = d, epoch = 0L,
    rms.relative.percent = score(Y, shared$target), proxy.rms.percent = NA_real_, fit.seconds = 0))
  for (epoch in seq.int(10L, 100L, 10L)) {
    began <- proc.time()[["elapsed"]]
    if (is.null(fit)) fit <- grip::metric.mds(prepared = shared$prep,
      approximation = "sparse", backend = "sgd", dim = d, init = Y, n.init = 1L,
      seed = p$fit.seed, max.iter = 10L, diagnostics = FALSE,
      pair.weights = "inverse_squared", sgd.control = list(scheduler = "zheng",
        schedule.epsilon = .1, schedule.epochs = 100L, checkpoint.every = 10L,
        keep.state = TRUE, max.workspace.bytes = p$native.workspace.bytes.per.fit))
    else fit <- grip::metric.mds(resume = fit, max.iter = 10L)
    elapsed <- elapsed + proc.time()[["elapsed"]] - began
    stopifnot(fit$metadata$engine == "sgd", ncol(fit$coords) == d,
      fit$metadata$stopping$epochs_completed == epoch,
      fit$metadata$settings$seed == p$fit.seed,
      fit$metadata$settings$sgd_control$schedule.epochs == 100L,
      all(is.finite(fit$coords)))
    a <- (z$count_i + z$count_j) / 2
    proxy <- 100 * sqrt(sum(a * (pair.dist(fit$coords, z$pairs) / z$targets - 1)^2) / sum(a))
    stopifnot(abs(proxy - 100 * fit$metadata$sparse_proxy_normalized_rmse) < 1e-8)
    rows[[length(rows) + 1L]] <- data.frame(dimension = d, epoch = epoch,
      rms.relative.percent = score(fit$coords, shared$target),
      proxy.rms.percent = proxy, fit.seconds = elapsed)
    atomic.rds(list(cache.key = shared$cache.key, dimension = d, start = start,
      fit = fit, scores = do.call(rbind, rows)),
      file.path(out, sprintf("dimension-%03d.rds", d)))
    write.csv(do.call(rbind, rows), file.path(out, sprintf("scores-%03d.csv", d)), row.names = FALSE)
    cat("dimension", d, "epoch", epoch, "RMS", tail(rows, 1L)[[1L]]$rms.relative.percent, "\n")
    flush.console(); completed <- epoch
  }
  stopifnot(identical(fingerprints, vapply(source.files, hash.file, "")))
} else stop("Unknown mode")
