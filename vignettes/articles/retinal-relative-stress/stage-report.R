# Independently validate saved fits, then generate compact publication artifacts.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
out <- normalizePath(args[1L]); baseline <- normalizePath(args[2L]); article <- normalizePath(args[3L])
stages <- c("E11", "E12", "E14", "E16", "E18", "P0", "P2", "P5", "P8", "P14")
dims <- c(3L, 4L, 5L, 10L, 20L, 25L, 50L)
status <- jsonlite::read_json(file.path(out, "status.json"), simplifyVector = TRUE)
stopifnot(status$status == "completed", status$completed_fits == 70L, status$full_run_unchanged)
manifest <- read.csv(file.path(out, "source-manifest.csv"))
stopifnot(identical(manifest$sha256, vapply(manifest$path,
  function(f) digest::digest(file = f, algo = "sha256"), "", USE.NAMES = FALSE)))
pooled <- read.csv(file.path(baseline, "dimension-summary.csv"))
pooled <- pooled[match(dims, pooled$dimension), ]
pooled.protocol <- readRDS(file.path(baseline, "prepared.rds"))$protocol
stopifnot(pooled.protocol$backend == "sgd", pooled.protocol$n == 1500L,
  pooled.protocol$input.dimension == 1814L, pooled.protocol$pivots == 200L,
  pooled.protocol$schedule.epochs == 100L, pooled.protocol$fit.seed == 20190619L)
gram.distances <- function(X) {
  G <- tcrossprod(X); v <- diag(G)
  sqrt(pmax(outer(v, v, "+") - 2 * G, 0))
}
all.rows <- list(); all.traces <- list(); preparations <- list(); cells.used <- integer()
for (stage in stages) {
  home <- file.path(out, stage)
  shared <- readRDS(file.path(home, "prepared.rds")); p <- shared$protocol
  cells <- read.csv(file.path(home, "cells.csv"))
  stopifnot(nrow(cells) == 1500L, !anyDuplicated(cells$index),
    all(cells$stage == stage), identical(sort(p$dimensions), dims),
    p$n == 1500L, p$input.dimension == 1814L, p$pivots == 200L,
    p$q == 10L, p$backend == "sgd", p$fit.seed == 20190619L,
    p$schedule.epochs == 100L, p$schedule.epsilon == .1)
  cells.used <- c(cells.used, cells$index)
  D <- gram.distances(shared$X); lower <- lower.tri(D)
  stopifnot(max(abs(D[lower] - shared$target)) < 1e-8)
  z <- shared$prep$constraints; w <- (z$count_i + z$count_j) / 2
  selected <- matrix(FALSE, p$n, p$n)
  selected[z$pairs[w > 0, , drop = FALSE]] <- TRUE
  selected <- selected | t(selected)
  F <- lower & selected; U <- lower & !selected; fraction <- sum(F) / sum(lower)
  for (d in dims) {
    result <- readRDS(file.path(home, sprintf("dimension-%03d.rds", d)))
    fit <- result$fit
    stopifnot(result$cache.key == shared$cache.key, fit$metadata$engine == "sgd",
      fit$metadata$stopping$epochs_completed == 100L,
      fit$metadata$settings$seed == p$fit.seed,
      fit$metadata$settings$sgd_control$schedule.epochs == 100L,
      identical(dim(fit$coords), c(1500L, d)))
    errors <- (gram.distances(fit$coords) / D - 1)^2
    actual <- 100 * sqrt(mean(errors[lower]))
    proxy <- 100 * sqrt(sum(w * errors[z$pairs]) / sum(w))
    fitted <- 100 * sqrt(mean(errors[F])); omitted <- 100 * sqrt(mean(errors[U]))
    row <- tail(result$scores, 1L)
    stopifnot(abs(row$rms.relative.percent - actual) < 1e-8,
      abs(row$proxy.rms.percent - proxy) < 1e-8,
      abs(actual^2 - (fraction * fitted^2 + (1 - fraction) * omitted^2)) < 1e-8)
    native <- fit$continuation$state$native$coords * fit$metadata$input_rms_distance
    native <- sweep(native, 2L, colMeans(native), "-")
    stopifnot(max(abs(native - fit$coords)) < 1e-8)
    row$fitted.rms.percent <- fitted; row$remaining.rms.percent <- omitted
    row$fitted.fraction <- fraction
    row$pooled.rms.percent <- pooled$rms.relative.percent[pooled$dimension == d]
    row$excess.over.pooled.pp <- actual - row$pooled.rms.percent
    all.rows[[length(all.rows) + 1L]] <- row
    all.traces[[length(all.traces) + 1L]] <- result$scores
  }
  preparations[[length(preparations) + 1L]] <- read.csv(file.path(home, "preparation.csv"))
  cat("Validated", stage, "seven final fits\n"); flush.console()
}
stopifnot(length(cells.used) == 15000L, !anyDuplicated(cells.used))
scores <- do.call(rbind, all.rows); prep <- do.call(rbind, preparations)
stopifnot(nrow(scores) == 70L)
write.csv(do.call(rbind, all.traces), file.path(out, "stage-epoch-scores.csv"), row.names = FALSE)
for (dest in c(out, article)) {
  write.csv(scores, file.path(dest, "stage-summary.csv"), row.names = FALSE)
  write.csv(prep, file.path(dest, "stage-preparations.csv"), row.names = FALSE)
}
table.rows <- vapply(stages, function(s) {
  r <- scores[scores$stage == s, ]; values <- r$rms.relative.percent[match(dims, r$dimension)]
  paste0(s, " & ", paste(sprintf("%.2f", values), collapse = " & "), " \\\\")
}, "")
writeLines(c("\\begin{tabular}{lrrrrrrr}\\toprule",
  "Stage & 3D & 4D & 5D & 10D & 20D & 25D & 50D \\\\ \\midrule", table.rows,
  "\\midrule", paste0("Pooled & ", paste(sprintf("%.2f", pooled$rms.relative.percent), collapse = " & "),
                        " \\\\"), "\\bottomrule\\end{tabular}"), file.path(article, "stage-results-table.tex"))
writeLines(c("\\begin{tabular}{lrrr}\\toprule",
  "Stage & Search PCs & Neighbor retention (\\%) & Median target distance \\\\ \\midrule",
  sprintf("%s & %d & %.1f & %.2f \\\\", prep$stage, prep$search.dimension,
          100 * prep$neighbor.retention, prep$target.median),
  "\\bottomrule\\end{tabular}"), file.path(article, "stage-preparations.tex"))
stamp <- format(as.POSIXct(status$updated, origin = "1970-01-01", tz = "UTC"),
                "%Y-%m-%d %H:%M:%S %Z", tz = "America/New_York")
writeLines(c(sprintf("\\newcommand{\\stagegenerated}{%s}", stamp),
  sprintf("\\newcommand{\\stagepeakmemory}{%.2f}", status$peak_rss_bytes / 1024^3),
  sprintf("\\newcommand{\\stageelapsedminutes}{%.2f}", status$elapsed_seconds / 60)),
  file.path(article, "stage-info.tex"))
draw <- function() {
  par(mfrow = c(5, 2), mar = c(2.0, 2.7, 1.7, .6), oma = c(3.4, 2.2, 2.6, .3),
      mgp = c(1.6, .5, 0), las = 1, cex = .85)
  ymax <- ceiling(max(scores$rms.relative.percent, pooled$rms.relative.percent) / 5) * 5
  for (stage in stages) {
    r <- scores[scores$stage == stage, ]
    plot(r$dimension, r$rms.relative.percent, type = "n", log = "x", axes = FALSE,
      ylim = c(0, ymax), xlim = range(dims), xlab = "", ylab = "")
    abline(h = seq(0, ymax, 10), col = "#eeeeee")
    axis(1, at = c(3, 5, 10, 25, 50), labels = c(3, 5, 10, 25, 50))
    axis(2, at = seq(0, ymax, 10)); box(bty = "l")
    lines(pooled$dimension, pooled$rms.relative.percent, col = "#222222", lty = 3, lwd = 1.6)
    lines(r$dimension, r$proxy.rms.percent, col = "#999999", lty = 2, lwd = 1.2, type = "b", pch = 1, cex = .6)
    lines(r$dimension, r$rms.relative.percent, col = "#087f8c", lwd = 1.7, type = "b", pch = 16, cex = .65)
    title(stage, adj = 0, line = .4, cex.main = 1)
  }
  mtext("Embedding dimension (log scale)", side = 1, outer = TRUE, line = 1.7, cex = .8)
  mtext("RMS relative distance error (%)", side = 2, outer = TRUE, line = .7, las = 0, cex = .8)
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = rep(0, 4), oma = rep(0, 4))
  plot.new()
  legend("top", inset = c(0, .005), horiz = TRUE, bty = "n", cex = .8,
    legend = c("Stage: all pairs", "Stage: sparse proxy", "Pooled: all pairs"),
    col = c("#087f8c", "#999999", "#222222"), lty = c(1, 2, 3), lwd = c(1.7, 1.2, 1.6))
}
pdf(file.path(article, "stage-error-curves.pdf"), width = 7, height = 7.5, pointsize = 10)
draw(); dev.off()
draw.contrasts <- function() {
  par(mfrow = c(4, 2), mar = c(3.0, 2.7, 1.6, .6), oma = c(2.5, 2.4, .4, .3),
      mgp = c(1.6, .5, 0), cex = .85)
  yr <- range(c(0, scores$excess.over.pooled.pp)) + c(-.2, .2)
  for (d in dims) {
    r <- scores[scores$dimension == d, ]
    r <- r[match(stages, r$stage), ]
    plot(seq_along(stages), r$excess.over.pooled.pp, type = "n", axes = FALSE,
      xlim = c(.6, 10.4), ylim = yr, xlab = "", ylab = "")
    abline(h = 0, col = "#555555", lty = 2)
    axis(1, at = seq_along(stages), labels = stages, las = 2, cex.axis = .85)
    axis(2, las = 1); box(bty = "l")
    points(seq_along(stages), r$excess.over.pooled.pp, pch = 16, col = "#087f8c")
    title(paste(d, "dimensions"), adj = 0, line = .3, cex.main = 1)
  }
  plot.new()
  text(.05, .9, "Above zero: higher stage error\nBelow zero: lower stage error\nOne fit per condition; no intervals.",
       adj = c(0, 1), cex = .75)
  mtext("Developmental stage (categorical spacing)", side = 1, outer = TRUE, line = 1, cex = .8)
  mtext("Stage minus pooled RMS error (percentage points)", side = 2, outer = TRUE,
        line = .8, cex = .8)
}
pdf(file.path(article, "stage-pooled-differences.pdf"), width = 7, height = 7, pointsize = 10)
draw.contrasts(); dev.off()
png(file.path(out, "stage-pooled-differences.png"), width = 1400, height = 1400, res = 200)
draw.contrasts(); dev.off()
png(file.path(out, "stage-error-curves.png"), width = 1400, height = 1500, res = 200)
draw(); dev.off()
writeLines(c(
  "PASS: 70 completed 100-epoch SGD fits, each on 1500 cells and 1814 input PCs",
  "PASS: all final all-pairs errors independently recomputed through Gram matrices",
  "PASS: sparse proxies, pair-set decompositions, terminal-state units and cache identities verified",
  "PASS: ten nonoverlapping stage samples, with 1500 distinct cells per stage",
  "PASS: source files unchanged and full-run pause/status records preserved by coordinator",
  "CAVEAT: one sample, one preparation and one initialization per stage; no uncertainty intervals",
  "CAVEAT: pooled and within-stage scores evaluate different cell-pair populations"),
  file.path(out, "validation.txt"))
print(reshape(scores[c("stage", "dimension", "rms.relative.percent")],
              idvar = "stage", timevar = "dimension", direction = "wide"), row.names = FALSE)
