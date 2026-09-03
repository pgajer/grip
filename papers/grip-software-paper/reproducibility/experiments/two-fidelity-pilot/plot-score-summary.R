#!/usr/bin/env Rscript
# Reproduce the paired-dot summary from the validated, Git-tracked pilot tables:
# Rscript plot-score-summary.R pilot-summary path/to/build/two-fidelity-pilot
# No package dependencies, fitted objects, or new optimization runs are needed.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
input <- normalizePath(args[1], mustWork = TRUE)
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(args[2], mustWork = TRUE)
scores <- read.csv(file.path(input, "layout-scores.csv"))
reference <- read.csv(file.path(input, "reference-validation.csv"))
status <- read.csv(file.path(input, "fit-status.csv"))

methods <- c("Original saddle", "Metric MDS", "MDS + edge-KK")
metrics <- c("path_rel", "edge_rel", "stress1")
clouds <- 1:5
stopifnot(nrow(scores) == 15L, setequal(scores$replicate, clouds),
          setequal(scores$method, methods),
          all(table(scores$replicate, scores$method) == 1L),
          all(is.finite(as.matrix(scores[metrics]))),
          all(as.matrix(scores[metrics]) >= 0),
          nrow(reference) == 5L, setequal(reference$replicate, clouds))
for (r in clouds) {
  stopifnot(all(scores$k[scores$replicate == r] ==
                reference$best_k[reference$replicate == r]))
}
primary <- status[status$run == "primary", ]
stopifnot(nrow(primary) == 25L, all(primary$budget == 200L),
          all(primary$hit_budget))

# Fixed method colors; the same small offset identifies each cloud across groups.
method_colors <- c("#737373", "#286EAB", "#B66027")
offsets <- seq(-0.16, 0.16, length.out = length(clouds))
values <- lapply(metrics, function(metric) {
  vapply(methods, function(method) {
    d <- scores[scores$method == method, ]
    100 * d[[metric]][match(clouds, d$replicate)]
  }, numeric(length(clouds)))
})
names(values) <- metrics
medians <- vapply(values, function(y) apply(y, 2, median), numeric(3))
write.csv(data.frame(method = methods, medians, check.names = FALSE),
          file.path(out, "pilot-score-summary-medians-percent.csv"), row.names = FALSE)

figure_caption <- paste(
  "Figure 1. Each small dot represents one independently sampled cloud; thin lines join results for the same cloud.",
  "Open circles and numeric labels show medians. Each selected graph is shared by all three configurations;",
  sprintf("selected k values are %s.", paste(reference$best_k[match(clouds, reference$replicate)], collapse = ", ")),
  "Fixed-path and chord scores use all 499,500 unordered pairs; edge scores use all graph edges.",
  "All scores use unprojected 3D coordinates. The three statistics have different normalizations, and vertical scales differ.",
  "All five clouds are included. All 25 primary edge-KK stages reached their 200-iteration limit;",
  "these are achieved errors, not certified minima. No confidence intervals are shown for this five-cloud pilot."
)
writeLines(figure_caption, file.path(out, "pilot-score-summary-caption.txt"))

draw <- function() {
  par(mfrow = c(1, 3), mar = c(4.3, 4.25, 2.15, 0.7),
      oma = c(5.9, 0.15, 3.7, 0.25), family = "sans", cex = 0.96,
      mgp = c(2.6, 0.7, 0), tcl = -0.25, las = 1)
  titles <- c("A   Fixed-path relative RMSE", "B   Edge relative RMSE",
              "C   MDS Stress-1")
  for (m in seq_along(values)) {
    y <- values[[m]]
    ymax <- max(y)
    plot(NA, xlim = c(0.65, 3.35), ylim = c(-0.045 * ymax, 1.14 * ymax),
         xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
         xlab = "", ylab = if (m < 3L) "Relative RMSE (%)" else "Stress-1 (%)")
    ticks <- pretty(c(0, ymax), n = 4)
    ticks <- ticks[ticks >= 0 & ticks <= 1.1 * ymax]
    abline(h = ticks, col = "#E9E9E9", lwd = 0.7)
    axis(2, at = ticks, col = "#555555", col.axis = "#333333", lwd = 0.75)
    axis(1, at = 1:3,
         labels = c("Original saddle\ncoordinates", "Metric\nMDS", "Metric MDS\n+ edge-KK"),
         tick = FALSE, line = 0.35, cex.axis = 0.93, col.axis = "#333333")
    box(bty = "l", col = "#555555", lwd = 0.75)
    for (r in seq_along(clouds)) {
      lines(1:3 + offsets[r], y[r, ], col = "#C2C2C2", lwd = 0.8)
    }
    for (j in seq_along(methods)) {
      points(j + offsets, y[, j], pch = 16, cex = 0.88,
             col = adjustcolor(method_colors[j], alpha.f = 0.7))
      # Transparent open rings preserve the individual dot at a shared median.
      points(j, medians[j, m], pch = 1, cex = 1.6, lwd = 1.65,
             col = method_colors[j])
      text(j, max(y[, j]) + 0.064 * ymax, sprintf("%.3f", medians[j, m]),
           col = method_colors[j], cex = 0.85)
    }
    title(titles[m], font.main = 1, cex.main = 1.03, line = 0.8, col.main = "#202020")
  }
  mtext("Three-dimensional embedding fidelity across saddle samples",
        side = 3, outer = TRUE, line = 2.0, cex = 1.25, col = "#202020")
  mtext("Five independent clouds; n = 1,000; graph calibrated separately for each cloud",
        side = 3, outer = TRUE, line = 0.7, cex = 0.94, col = "#555555")
  lines <- strwrap(figure_caption, width = 153)
  for (i in seq_along(lines)) {
    mtext(lines[i], side = 1, outer = TRUE, line = 0.4 + 1.02 * (i - 1),
          cex = 0.77, col = "#333333")
  }
}

pdf(file.path(out, "pilot-score-summary.pdf"), width = 12, height = 6.1,
    useDingbats = FALSE, title = "Three-dimensional embedding fidelity across saddle samples")
draw()
dev.off()
png(file.path(out, "pilot-score-summary.png"), width = 12, height = 6.1,
    units = "in", res = 200, type = "cairo")
draw()
dev.off()
print(data.frame(method = methods, medians, check.names = FALSE), row.names = FALSE)
