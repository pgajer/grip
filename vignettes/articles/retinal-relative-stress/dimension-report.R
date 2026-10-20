# Report-only diagnostics: consume completed fits without any optimization.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
out <- normalizePath(args[1L]); article <- normalizePath(args[2L])
shared <- readRDS(file.path(out, "prepared.rds"))
summary <- read.csv(file.path(out, "dimension-summary.csv"))
stopifnot(identical(summary, read.csv(file.path(article, "dimension-summary.csv"))))
z <- shared$prep$constraints
w <- (z$count_i + z$count_j) / 2
active <- w > 0
pairs <- z$pairs[active, , drop = FALSE]
D <- as.matrix(dist(shared$X))
selected <- matrix(FALSE, nrow(D), ncol(D))
selected[pairs] <- TRUE; selected <- selected | t(selected)
lower <- lower.tri(D)
F <- lower & selected; U <- lower & !selected
stopifnot(sum(F) == nrow(pairs), sum(F) + sum(U) == choose(nrow(D), 2))
decomposition <- do.call(rbind, lapply(c(100L, 200L), function(d) {
  fit <- readRDS(file.path(out, sprintf("dimension-%03d.rds", d)))
  stopifnot(fit$cache.key == shared$cache.key,
            fit$fit$metadata$stopping$epochs_completed == 100L)
  E <- (as.matrix(dist(fit$fit$coords)) / D - 1)^2
  proxy <- sum(w[active] * E[pairs]) / sum(w[active])
  all <- mean(E[lower]); fitted <- mean(E[F]); omitted <- mean(E[U])
  fraction <- sum(F) / sum(lower)
  stopifnot(abs(all - (fraction * fitted + (1 - fraction) * omitted)) < 1e-12,
    abs(100 * sqrt(all) - summary$rms.relative.percent[summary$dimension == d]) < 1e-8,
    abs(100 * sqrt(proxy) - summary$proxy.rms.percent[summary$dimension == d]) < 1e-8)
  data.frame(dimension = d, all.rms.percent = 100 * sqrt(all),
    fitted.rms.percent = 100 * sqrt(fitted), omitted.rms.percent = 100 * sqrt(omitted),
    proxy.rms.percent = 100 * sqrt(proxy), fitted.pairs = sum(F),
    omitted.pairs = sum(U), fitted.fraction = fraction,
    fitted.contribution = fraction * fitted, omitted.contribution = (1 - fraction) * omitted)
}))
write.csv(decomposition, file.path(article, "dimension-decomposition.csv"), row.names = FALSE)
write.csv(decomposition, file.path(out, "dimension-decomposition.csv"), row.names = FALSE)
writeLines(c("\\begin{tabular}{rrr}\\toprule",
  "$d$ & All-pairs RMS (\\%) & Sparse proxy RMS (\\%) \\\\ \\midrule",
  sprintf("%d & %.3f & %.3f \\\\", summary$dimension, summary$rms.relative.percent,
          summary$proxy.rms.percent), "\\bottomrule\\end{tabular}"),
  file.path(article, "dimension-table.tex"))
writeLines(c("\\begin{tabular}{rrrrr}\\toprule",
  "$d$ & $R$ & $R_F$ & $R_U$ & Proxy \\\\ \\midrule",
  sprintf("%d & %.3f & %.3f & %.3f & %.3f \\\\", decomposition$dimension,
    decomposition$all.rms.percent, decomposition$fitted.rms.percent,
    decomposition$omitted.rms.percent, decomposition$proxy.rms.percent),
  "\\bottomrule\\end{tabular}"), file.path(article, "dimension-decomposition.tex"))
status <- jsonlite::read_json(file.path(out, "status.json"), simplifyVector = TRUE)
stamp <- format(as.POSIXct(status$updated, origin = "1970-01-01", tz = "UTC"),
                "%Y-%m-%d %H:%M:%S %Z", tz = "America/New_York")
writeLines(c(sprintf("\\newcommand{\\dimensiongenerated}{%s}", stamp),
  sprintf("\\newcommand{\\dimensionfittedpairs}{%s}", format(sum(F), big.mark = ",")),
  sprintf("\\newcommand{\\dimensionomittedpairs}{%s}", format(sum(U), big.mark = ",")),
  sprintf("\\newcommand{\\dimensionfittedfraction}{%.2f}", 100 * sum(F) / sum(lower))),
  file.path(article, "dimension-info.tex"))
pdf(file.path(article, "dimension-error.pdf"), width = 7, height = 4.5, pointsize = 10)
par(mar = c(4.2, 4.5, 1, .5), mgp = c(2.8, .7, 0), las = 1)
plot(summary$dimension, summary$rms.relative.percent, type = "n", log = "x", axes = FALSE,
  ylim = c(0, 33), xlab = "Embedding dimension (log scale)",
  ylab = "RMS relative distance error (%)")
abline(h = seq(0, 30, 5), col = "#eeeeee")
axis(1, at = summary$dimension, labels = summary$dimension, cex.axis = .85)
axis(2); box(bty = "l")
lines(summary$dimension, summary$rms.relative.percent, col = "#087f8c", lwd = 2,
      type = "b", pch = 16)
lines(summary$dimension, summary$proxy.rms.percent, col = "#555555", lwd = 1.5,
      type = "b", pch = 1, lty = 2)
legend("topright", c("All-pairs relative error", "Sparse fitting proxy"),
       col = c("#087f8c", "#555555"), lty = c(1, 2), pch = c(16, 1), lwd = 2, bty = "n")
dev.off()
annotations <- read.csv(args[3L])
stages <- c("E11", "E12", "E14", "E16", "E18", "P0", "P2", "P5", "P8", "P14")
stopifnot(nrow(annotations) == 120804L, setequal(unique(annotations$stage), stages))
counts <- data.frame(stage = stages, available = as.integer(table(annotations$stage)[stages]),
                     sample = 1500L)
stopifnot(all(counts$available >= counts$sample))
write.csv(counts, file.path(article, "stage-counts.csv"), row.names = FALSE)
writeLines(c("\\begin{tabular}{lrr}\\toprule",
  "Stage & Available cells & Stage sample \\\\ \\midrule",
  sprintf("%s & %s & 1,500 \\\\", counts$stage, format(counts$available, big.mark = ",", trim = TRUE)),
  "\\bottomrule\\end{tabular}"), file.path(article, "stage-counts.tex"))
print(decomposition, row.names = FALSE)
