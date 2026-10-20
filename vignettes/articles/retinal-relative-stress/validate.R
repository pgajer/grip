args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- normalizePath(args[1])
p <- readRDS(file.path(out, "protocol.rds"))
preps <- readRDS(file.path(out, "preparations.rds"))
d <- read.csv(file.path(out, "scores.csv"))
D <- as.matrix(dist(p$X)); diag(D) <- Inf
neighbors <- function(D) t(vapply(seq_len(nrow(D)),
  function(i) order(D[i, ])[seq_len(10L)], integer(10L)))
original <- neighbors(D)
prep <- preps[[1]]
Z <- sweep(p$X, 2, prep$pca$center, "-") %*% prep$pca$loadings
projected <- as.matrix(dist(Z)); diag(projected) <- Inf
projected <- neighbors(projected)
overlap <- vapply(seq_len(nrow(p$X)),
  function(i) length(intersect(original[i, ], projected[i, ])) / 10, 0)
write.csv(data.frame(mean.retention = mean(overlap), minimum = min(overlap),
  q05 = unname(quantile(overlap, .05)), median = median(overlap)),
  file.path(out, "neighbor-retention.csv"), row.names = FALSE)
for (s in p$protocol$fit.seeds) {
  initial <- d[d$seed == s & d$epoch == 0, "rms.relative.percent"]
  stopifnot(diff(range(initial)) < 1e-12)
}
for (m in p$protocol$methods) for (s in p$protocol$fit.seeds) {
  f <- readRDS(file.path(out, sprintf("%s-%s.rds", m, s)))
  stopifnot(identical(f$cache.key, p$cache.key),
            identical(f$seed, s), identical(f$method, m))
  for (epoch in p$protocol$epochs) {
    Y <- f$coords[[as.character(epoch)]]
    # Independent row-wise evaluation, not stats::dist() used by the runner.
    loss <- 0; count <- 0L
    for (i in seq_len(nrow(Y)-1L)) {
      j <- seq.int(i+1L, nrow(Y))
      embedded <- sqrt(rowSums(sweep(Y[j, , drop = FALSE], 2, Y[i, ], "-")^2))
      loss <- loss + sum((embedded / D[i, j] - 1)^2)
      count <- count + length(j)
    }
    expected <- d[d$method == m & d$seed == s & d$epoch == epoch,
                  "rms.relative.percent"]
    stopifnot(abs(100 * sqrt(loss / count) - expected) < 1e-9)
  }
}
writeLines(c("PASS: identical initial evaluation for each matched start",
  "PASS: independent row-wise recomputation of all 108 post-fit evaluations",
  "PASS: saved result cache keys, method labels and seeds match protocol",
  "Neighbor-retention diagnostic is descriptive; it does not change Euclidean targets."),
  file.path(out, "independent-validation.txt"))
cat("Independent evaluation checks passed; neighbor retention:\n")
print(read.csv(file.path(out, "neighbor-retention.csv")))
