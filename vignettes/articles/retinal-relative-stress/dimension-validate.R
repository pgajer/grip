args <- commandArgs(trailingOnly = TRUE)
out <- normalizePath(args[1L])
shared <- readRDS(file.path(out, "prepared.rds"))
manifest <- read.csv(file.path(out, "source-manifest.csv"))
stopifnot(identical(manifest$sha256, vapply(manifest$path,
  function(f) digest::digest(file = f, algo = "sha256"), "", USE.NAMES = FALSE)))
# Recompute distances through a Gram matrix, independently of stats::dist().
distances <- function(X) {
  g <- tcrossprod(X); s <- diag(g)
  D <- sqrt(pmax(outer(s, s, "+") - 2 * g, 0))
  D[lower.tri(D)]
}
targets <- distances(shared$X)
stopifnot(max(abs(targets - shared$target)) < 1e-8)
rows <- lapply(shared$protocol$dimensions, function(d) {
  result <- readRDS(file.path(out, sprintf("dimension-%03d.rds", d)))
  fit <- result$fit
  stopifnot(identical(result$cache.key, shared$cache.key),
    fit$metadata$engine == "sgd", ncol(fit$coords) == d,
    fit$metadata$stopping$epochs_completed == 100L,
    fit$metadata$settings$sgd_control$schedule.epochs == 100L)
  actual <- 100 * sqrt(mean((distances(fit$coords) / targets - 1)^2))
  row <- tail(result$scores, 1L)
  stopifnot(abs(row$rms.relative.percent - actual) < 1e-8)
  native <- fit$continuation$state$native$coords * fit$metadata$input_rms_distance
  native <- sweep(native, 2L, colMeans(native), "-")
  stopifnot(max(abs(native - fit$coords)) < 1e-8)
  row
})
summary <- do.call(rbind, rows)
pca <- read.csv(file.path(out, "pca-reference.csv"))
names(pca)[2:3] <- c("pca.rms.relative.percent", "pca.subset.variance.percent")
summary <- merge(summary, pca, by = "dimension", sort = TRUE)
write.csv(summary, file.path(out, "dimension-summary.csv"), row.names = FALSE)
trace <- do.call(rbind, lapply(shared$protocol$dimensions, function(d)
  read.csv(file.path(out, sprintf("scores-%03d.csv", d)))))
write.csv(trace, file.path(out, "dimension-epoch-scores.csv"), row.names = FALSE)
png(file.path(out, "dimension-error.png"), width = 1600, height = 1000, res = 160)
par(mar = c(5.1, 5.2, 4.5, 1.2), mgp = c(3, .8, 0), las = 1)
yr <- c(0, max(summary$rms.relative.percent, summary$proxy.rms.percent) * 1.08)
plot(summary$dimension, summary$rms.relative.percent, type = "n", log = "x", axes = FALSE,
  ylim = yr, xlab = "Embedding dimension (log scale)", ylab = "RMS relative distance error (%)")
abline(h = pretty(yr), col = "#eeeeee", lwd = 1)
axis(1, at = summary$dimension, labels = summary$dimension, cex.axis = .85)
axis(2); box(bty = "l")
lines(summary$dimension, summary$rms.relative.percent, col = "#087f8c", lwd = 2.5,
      type = "b", pch = 16)
lines(summary$dimension, summary$proxy.rms.percent, col = "#555555", lwd = 2,
      type = "b", pch = 1, lty = 2)
legend("topright", c("All-pairs relative error", "Sparse fitting proxy"),
       col = c("#087f8c", "#555555"), lty = c(1, 2), pch = c(16, 1), lwd = 2, bty = "n")
title(main = "Retinal distance error after 100 SGD epochs", line = 2.7, adj = 0)
mtext("1,500 cells; targets in 1,814 PCs; all 1,124,250 pairs; no post-fit rescaling",
      side = 3, line = 1.2, adj = 0, cex = .85)
dev.off()
writeLines(c(
  "PASS: all ten dimensions completed 100 epochs with a 100-epoch schedule",
  "PASS: backend SGD, matching preparation cache keys and unchanged source files",
  "PASS: all final errors independently recomputed through Gram matrices (tolerance 1e-8 percentage points)",
  "PASS: returned coordinates match terminal optimizer state in original distance units",
  "CAVEAT: one start per dimension, one shared 200-landmark preparation; not a convergence certificate",
  "CAVEAT: PCA minimizes squared reconstruction error, not this relative-distance objective",
  "CAVEAT: results concern this 1500-cell subset, not the paused 120804-cell fit"),
  file.path(out, "validation.txt"))
print(summary, row.names = FALSE)
