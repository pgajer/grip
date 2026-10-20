# Render figures and compact tables from completed fits; never rerun optimization.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
out <- normalizePath(args[1]); article <- normalizePath(args[2])
d <- read.csv(file.path(out, "scores.csv"))
status <- read.csv(file.path(out, "fit-status.csv"))
stopifnot(nrow(status) == 9L, all(status$status == "ok"),
          all(status$completed.epoch == 500L),
          !anyDuplicated(d[c("method", "seed", "epoch")]), nrow(d) == 117L)
p <- readRDS(file.path(out, "protocol.rds"))
methods <- p$protocol$methods
labels <- c("Regional: 25 landmarks", "Regional: 200 landmarks", "All pairs")
colors <- c("#B44A31", "#247C83", "#333333")
makeplot <- function(name, fun, width = 7, height = 4.1) {
  pdf(file.path(article, name), width = width, height = height, useDingbats = FALSE)
  fun()
  dev.off()
}
makeplot("error-trajectories.pdf", function() {
  par(mar = c(4.2, 4.4, .8, .6), las = 1)
  x <- d[d$epoch > 0, ]
  plot(range(x$epoch), range(x$rms.relative.percent), type = "n",
       xlab = "Completed SGD epochs", ylab = "RMS relative distance error (%)")
  for (i in seq_along(methods)) {
    z <- x[x$method == methods[i], ]
    for (s in unique(z$seed)) {
      zz <- z[z$seed == s, ]
      lines(zz$epoch, zz$rms.relative.percent,
            col = adjustcolor(colors[i], .35), lty = i, lwd = 1)
    }
    a <- aggregate(rms.relative.percent ~ epoch, z, mean)
    lines(a$epoch, a$rms.relative.percent, col = colors[i], lty = i, lwd = 2.5)
  }
  legend("topright", labels, col = colors, lty = 1:3, lwd = 2.5, bty = "n", cex = .85)
})
makeplot("proxy-trajectories.pdf", function() {
  par(mfrow = c(1, 2), mar = c(4, 4.1, 2, .5), las = 1)
  all.means <- aggregate(cbind(rms.relative.percent, proxy.rms.percent) ~ method + epoch,
    d[d$method != "all_pairs" & d$epoch > 0, ], mean)
  common.ylim <- range(all.means[c("rms.relative.percent", "proxy.rms.percent")])
  for (i in 1:2) {
    z <- d[d$method == methods[i] & d$epoch > 0, ]
    a <- aggregate(cbind(rms.relative.percent, proxy.rms.percent) ~ epoch, z, mean)
    plot(a$epoch, a$rms.relative.percent, type = "l", lwd = 2, col = colors[i],
      ylim = common.ylim,
      main = labels[i], cex.main = .9, xlab = "SGD epochs", ylab = "RMS error (%)")
    lines(a$epoch, a$proxy.rms.percent, lwd = 2, lty = 2, col = "#333333")
    legend("topright", c("All-pairs evaluation", "Sparse proxy"),
      col = c(colors[i], "#333333"), lty = 1:2, lwd = 2, bty = "n", cex = .75)
  }
}, height = 3.5)
final <- d[d$epoch == 500L, ]
ref <- final[final$method == "all_pairs", c("seed", "rms.relative.percent")]
names(ref)[2] <- "reference"
delta <- merge(final[final$method != "all_pairs", ], ref, by = "seed")
delta$difference <- delta$rms.relative.percent - delta$reference
summaries <- lapply(methods[1:2], function(m) {
  z <- delta[delta$method == m, ]
  v <- z$difference; half <- qt(.975, length(v)-1) * sd(v)/sqrt(length(v))
  data.frame(method = m, mean.difference = mean(v), lower = mean(v)-half,
    upper = mean(v)+half, better = sum(v < -1e-8), tied = sum(abs(v) <= 1e-8),
    worse = sum(v > 1e-8))
})
delta.summary <- do.call(rbind, summaries)
makeplot("paired-differences.pdf", function() {
  par(mar = c(4.3, 4.5, .5, 1), las = 1)
  lim <- range(c(0, delta$difference, delta.summary$lower, delta.summary$upper))
  plot(c(1, 2), lim, type = "n", xaxt = "n", xlim = c(.6, 2.4),
    xlab = "", ylab = "Excess RMS error (percentage points)")
  abline(h = 0, col = "gray50", lty = 3)
  axis(1, at = 1:2, labels = c("25 landmarks", "200 landmarks"))
  for (i in 1:2) {
    z <- delta[delta$method == methods[i], ]; a <- delta.summary[i, ]
    points(i + c(-.08, 0, .08), z$difference, pch = 1:3, col = colors[i])
    arrows(i+.15, a$lower, i+.15, a$upper, angle = 90, code = 3, length = .06)
    points(i+.15, a$mean.difference, pch = 18, cex = 1.4)
  }
}, height = 3.7)
write.csv(delta.summary, file.path(out, "paired-summary.csv"), row.names = FALSE)
write.csv(d, file.path(article, "scores.csv"), row.names = FALSE)
write.csv(status, file.path(article, "fit-status.csv"), row.names = FALSE)
write.csv(read.csv(file.path(out, "preparations.csv")),
          file.path(article, "preparations.csv"), row.names = FALSE)
write.csv(delta.summary, file.path(article, "paired-summary.csv"), row.names = FALSE)
small <- d[d$epoch %in% c(10L, 50L, 80L, 100L, 110L, 500L), ]
table <- aggregate(rms.relative.percent ~ epoch + method, small, mean)
lines <- c("\\begin{tabular}{rrrr}", "\\toprule",
  "Epoch & Regional, 25 & Regional, 200 & All pairs \\\\", "\\midrule")
for (ep in sort(unique(table$epoch))) {
  z <- table[table$epoch == ep, ]
  vals <- z$rms.relative.percent[match(methods, z$method)]
  lines <- c(lines, sprintf("%d & %.3f & %.3f & %.3f \\\\", ep, vals[1], vals[2], vals[3]))
}
writeLines(c(lines, "\\bottomrule", "\\end{tabular}"), file.path(article, "error-table.tex"))
lines <- c("\\begin{tabular}{lrrr}", "\\toprule",
  "Method & Start 1 & Start 2 & Start 3 \\\\", "\\midrule")
for (i in seq_along(methods)) {
  z <- final[final$method == methods[i], ]
  z <- z[order(z$seed), ]
  lines <- c(lines, sprintf("%s & %.3f & %.3f & %.3f \\\\",
    labels[i], z$rms.relative.percent[1], z$rms.relative.percent[2], z$rms.relative.percent[3]))
}
writeLines(c(lines, "\\bottomrule", "\\end{tabular}"), file.path(article, "final-table.tex"))
metadata <- c(paste0("\\newcommand{\\resultgenerated}{",
  readLines(file.path(out, "generated-at.txt")), "}"),
  paste0("\\newcommand{\\resultkey}{\\texttt{", substr(p$cache.key, 1, 16), "}}"))
writeLines(metadata, file.path(article, "results-info.tex"))
means <- vapply(methods, function(m) mean(final$rms.relative.percent[final$method == m]), 0)
proxy.means <- vapply(methods, function(m) mean(final$proxy.rms.percent[final$method == m]), 0)
stopifnot(all(delta.summary$worse == 3L), means[1] > means[2], means[2] > means[3],
          proxy.means[1] < proxy.means[2])
interpretation <- sprintf(paste0(
  "At epoch 500, mean error is %.3f\\%% for 25 landmarks, %.3f\\%% for ",
  "200 landmarks, and %.3f\\%% for the all-pairs reference. ",
  "Both regional methods have higher error than the reference in all three ",
  "matched starts (zero better, zero tied, three worse for each comparison). ",
  "For 25 landmarks the mean excess is %.3f percentage points, with a ",
  "paired 95-percent interval [%.3f, %.3f]. For 200 landmarks it is %.3f ",
  "percentage points, with interval [%.3f, %.3f]. ",
  "The very narrow latter interval reflects similar final results across ",
  "these three starts; it does not include cell-sampling or preparation uncertainty.\n\n",
  "Increasing the landmark count reduces the observed approximation gap ",
  "in this subset, but does not eliminate it. The experiment does not show ",
  "that 200 landmarks are sufficient for another population, nor that the ",
  "2,000-landmark full fit will achieve these errors."),
  means[1], means[2], means[3],
  delta.summary$mean.difference[1], delta.summary$lower[1], delta.summary$upper[1],
  delta.summary$mean.difference[2], delta.summary$lower[2], delta.summary$upper[2])
interpretation <- c(interpretation, sprintf(paste0(
  "\n\nThe proxy illustrates why the distinction matters: the 25-landmark ",
  "method has a mean proxy RMS of %.3f\\%%, below %.3f\\%% for 200 landmarks, ",
  "even though its all-pairs error is higher. Ranking these two fits by ",
  "their own proxies would reverse the ranking given by the common ",
  "evaluation metric. The all-pairs reference has equal proxy and evaluation ",
  "values, %.3f\\%%, as required by the unit-multiplier identity."),
  proxy.means[1], proxy.means[2], proxy.means[3]))
writeLines(interpretation, file.path(article, "interpretation.tex"))
retention <- read.csv(file.path(out, "neighbor-retention.csv"))
writeLines(sprintf(paste0("The mean retention of original 20-PC ten-neighbor sets ",
  "under the 13-dimensional search projection is %.2f\\%%; the median is ",
  "%.0f\\%% and the fifth percentile is %.0f\\%%. Thus high retained variance ",
  "does not imply identical local interactions."),
  100 * retention$mean.retention, 100 * retention$median, 100 * retention$q05),
  file.path(article, "neighbor-check.tex"))
write.csv(retention, file.path(article, "neighbor-retention.csv"), row.names = FALSE)
cat("Final scores:\n"); print(final[c("method", "seed", "rms.relative.percent", "fit.seconds")])
cat("Paired differences and 95% t intervals across three starts:\n"); print(delta.summary)
