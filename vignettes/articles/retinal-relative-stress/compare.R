# Deliberately small dense reference; never loads or resumes the full retinal fit.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
root <- normalizePath(args[1])
input <- normalizePath(args[2])
out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(out, "protocol.rds"))) stop("Use a fresh output directory")
Sys.setenv(TZ = "America/New_York")
pkgload::load_all(file.path(dirname(root), "dgraphs"), recompile = FALSE, quiet = TRUE)
pkgload::load_all(root, recompile = FALSE, quiet = TRUE)
hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2)
files <- c(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE),
           list.files(file.path(root, "src"), "\\.(cpp|h|hpp|so)$", full.names = TRUE),
           list.files(file.path(dirname(root), "dgraphs/R"), "\\.R$", full.names = TRUE),
           list.files(file.path(dirname(root), "dgraphs/src"), "\\.(cpp|h|hpp|so)$",
                      full.names = TRUE),
           file.path(root, "vignettes/articles/retinal-relative-stress/compare.R"))
fingerprints <- vapply(files, function(f) digest::digest(file = f, algo = "sha256"), "")
write.csv(data.frame(path = files, sha256 = fingerprints),
          file.path(out, "source-manifest.csv"), row.names = FALSE)
n <- 1500L
epochs <- c(1L, 10L, 25L, 50L, 80L, 100L, 110L, 150L, 200L, 300L, 400L, 500L)
seeds <- c(20190619L, 20190620L, 20190621L)
pc <- readRDS(input)
stopifnot(nrow(pc$scores) == 120804L, ncol(pc$scores) == 20L)
set.seed(20261020)
indices <- sort(sample.int(nrow(pc$scores), n))
X <- pc$scores[indices, , drop = FALSE]
pairs <- t(combn(n, 2L))
storage.mode(pairs) <- "integer"
D <- as.matrix(dist(X))
targets <- D[pairs]
stopifnot(all(is.finite(X)), all(targets > 0))
prepared <- lapply(c(25L, 200L), function(h) {
  grip::prepare.sparse.mds(X, method = "euclidean", n.pivots = h,
    pivot.selection = "randomized", q = 10L, neighbor.variance.target = .9,
    seed = 20261021L, max.workspace.bytes = 32 * 1024^3)
})
names(prepared) <- c("regional_25", "regional_200")
constraints <- list(pairs = pairs, targets = targets,
                    count_i = rep(1, length(targets)), count_j = rep(1, length(targets)))
score <- function(Y) {
  d <- as.numeric(dist(Y))
  target <- as.numeric(as.dist(D))
  100 * sqrt(mean((d / target - 1)^2))
}
proxy <- function(Y, z) {
  d <- sqrt(rowSums((Y[z$pairs[, 1], , drop = FALSE] -
                     Y[z$pairs[, 2], , drop = FALSE])^2))
  a <- (z$count_i + z$count_j) / 2
  100 * sqrt(sum(a * (d / z$targets - 1)^2) / sum(a))
}
starts <- lapply(seeds, function(s) {
  set.seed(s)
  Y <- scale(matrix(rnorm(n * 3L), n, 3L), center = TRUE, scale = FALSE)
  # Shared starting scale is fitted once on all pairs, then used by every method.
  d <- as.numeric(dist(Y)); t <- as.numeric(as.dist(D))
  Y * sum(d / t) / sum((d / t)^2)
})
protocol <- list(n = n, input.sha256 = digest::digest(file = input, algo = "sha256"),
  subset.seed = 20261020L, preparation.seed = 20261021L, fit.seeds = seeds,
  epochs = epochs, schedule.epochs = 500L, schedule = "zheng", epsilon = .1,
  q = 10L, neighbor.variance.target = .9, pivots = c(25L, 200L),
  backend = "sgd", starts.per.call = 1L, target = "ambient_euclidean_pc20",
  evaluation = "all unordered pairs; no fitted evaluation scale",
  methods = c(names(prepared), "all_pairs"), source.hash = hash(fingerprints))
key <- hash(protocol)
saveRDS(list(protocol = protocol, cache.key = key, X = X, cells = pc$cells[indices],
             indices = indices, starts = starts), file.path(out, "protocol.rds"))
saveRDS(prepared, file.path(out, "preparations.rds"))
prep.rows <- lapply(names(prepared), function(m) {
  p <- prepared[[m]]; z <- p$constraints
  stopifnot(max(abs(z$targets - D[z$pairs])) < 1e-8)
  data.frame(method = m, pairs = nrow(z$pairs),
    asymmetric.fraction = mean(z$count_i != z$count_j),
    zero.endpoint.fraction = mean(z$count_i == 0 | z$count_j == 0),
    pca.dimension = p$pca$dimension, pca.fraction = p$pca$achieved)
})
write.csv(do.call(rbind, prep.rows), file.path(out, "preparations.csv"), row.names = FALSE)

# Exact objective identity and continuation checks on a tiny synthetic reference.
set.seed(17)
T <- matrix(rnorm(20 * 4), 20, 4); P <- t(combn(20L, 2L))
C <- list(pairs = P, targets = as.matrix(dist(T))[P])
Y <- matrix(rnorm(20 * 3), 20, 3)
control <- list(scheduler = "zheng", schedule.epsilon = .1,
  schedule.epochs = 10L, checkpoint.every = 1L, keep.state = TRUE,
  max.workspace.bytes = 32 * 1024^3)
testfit <- function(k) grip::metric.mds(n = 20L, constraints = C,
  approximation = "sparse", backend = "sgd", dim = 3L, init = Y, seed = 7L,
  n.init = 1L, max.iter = k, diagnostics = FALSE, sgd.control = control)
a <- testfit(10L); b <- grip::metric.mds(resume = testfit(4L), max.iter = 6L)
stopifnot(identical(a$coords, b$coords),
  abs(a$metadata$sparse_proxy_normalized_rmse -
      sqrt(mean((as.numeric(dist(a$coords)) / as.numeric(dist(T)) - 1)^2))) < 1e-12)
writeLines(c("PASS: all-pair proxy equals independently computed relative stress",
             "PASS: 4+6 continuation gives identical coordinates to 10 epochs",
             "PASS: regional targets equal original PC20 Euclidean distances"),
           file.path(out, "validation.txt"))
rm(a, b)
rows <- list(); status.rows <- list()
for (s in seq_along(seeds)) for (method in protocol$methods) {
  cat(format(Sys.time()), "starting", method, "seed", seeds[s], "\n"); flush.console()
  z <- if (method == "all_pairs") constraints else prepared[[method]]$constraints
  fit <- NULL; completed <- 0L; seconds <- 0
  coords <- list(`0` = starts[[s]])
  one <- function(epoch, Y, proxy.value, rate) data.frame(
    method = method, seed = seeds[s], epoch = epoch,
    rms.relative.percent = score(Y), proxy.rms.percent = proxy.value,
    learning.rate = rate, fit.seconds = seconds, pairs.per.epoch = nrow(z$pairs))
  rows[[length(rows) + 1L]] <- one(0L, starts[[s]], proxy(starts[[s]], z), NA_real_)
  result <- tryCatch({
    for (epoch in epochs) {
      began <- proc.time()[["elapsed"]]
      if (is.null(fit)) {
        common <- list(approximation = "sparse", backend = "sgd", dim = 3L,
          init = starts[[s]], n.init = 1L, seed = seeds[s], max.iter = epoch,
          diagnostics = FALSE, pair.weights = "inverse_squared",
          sgd.control = list(scheduler = "zheng", schedule.epsilon = .1,
            schedule.epochs = 500L, checkpoint.every = 10L, keep.state = TRUE,
            max.workspace.bytes = 32 * 1024^3))
        if (method == "all_pairs") common <- c(common, list(n = n, constraints = constraints))
        else common <- c(common, list(prepared = prepared[[method]]))
        fit <- do.call(grip::metric.mds, common)
        initial <- fit$trace$sparse_proxy_stress[fit$trace$epoch == 0L]
        energy <- sum((z$count_i + z$count_j) / 2)
        stopifnot(abs(100 * sqrt(initial / energy) - proxy(starts[[s]], z)) < 1e-8)
      } else fit <- grip::metric.mds(resume = fit, max.iter = epoch - completed)
      seconds <- seconds + proc.time()[["elapsed"]] - began
      stopifnot(fit$metadata$engine == "sgd",
        fit$metadata$stopping$epochs_completed == epoch,
        fit$metadata$settings$seed == seeds[s],
        fit$metadata$settings$sgd_control$schedule.epochs == 500L,
        all(is.finite(fit$coords)))
      native <- fit$continuation$state$native$coords *
        fit$metadata$input_rms_distance
      native <- sweep(native, 2L, colMeans(native), "-")
      stopifnot(max(abs(native - fit$coords)) < 1e-8)
      q <- 100 * fit$metadata$sparse_proxy_normalized_rmse
      stopifnot(abs(q - proxy(fit$coords, z)) < 1e-8)
      if (method == "all_pairs") stopifnot(abs(q - score(fit$coords)) < 1e-8)
      rows[[length(rows) + 1L]] <- one(epoch, fit$coords, q,
        tail(fit$trace$learning_rate, 1L))
      coords[[as.character(epoch)]] <- fit$coords
      completed <- epoch
      write.csv(do.call(rbind, rows), file.path(out, "scores.csv"), row.names = FALSE)
      cat(method, seeds[s], epoch, tail(rows, 1L)[[1]]$rms.relative.percent, "\n")
      flush.console()
    }
    saveRDS(list(cache.key = key, method = method, seed = seeds[s], coords = coords,
      trace = fit$trace, metadata = fit$metadata[setdiff(names(fit$metadata), c("sparse", "history"))]),
      file.path(out, sprintf("%s-%s.rds", method, seeds[s])))
    NULL
  }, error = function(e) conditionMessage(e))
  status.rows[[length(status.rows) + 1L]] <- data.frame(method = method, seed = seeds[s],
    status = if (is.null(result)) "ok" else "failed", completed.epoch = completed,
    message = if (is.null(result)) "" else result)
  write.csv(do.call(rbind, status.rows), file.path(out, "fit-status.csv"), row.names = FALSE)
  rm(fit); gc()
}
stopifnot(identical(fingerprints, vapply(files,
  function(f) digest::digest(file = f, algo = "sha256"), "")))
writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))
writeLines(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), file.path(out, "generated-at.txt"))
if (any(do.call(rbind, status.rows)$status != "ok")) stop("Fit failures preserved")
cat("All nine fits completed; source fingerprints unchanged.\n")
