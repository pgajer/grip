#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  file.path("output", "gkk_lgkk_paper", "tmp", "sampled-final-rounds-suite-2026-03-28")
}

overall_path <- file.path(root, "sampled-final-rounds-overall-summary.csv")
best_family_path <- file.path(root, "sampled-final-rounds-best-by-family.csv")
best_graph_path <- file.path(root, "sampled-final-rounds-best-by-graph.csv")
summary_path <- file.path(root, "sampled-final-rounds-summary.csv")

if (!file.exists(overall_path) || !file.exists(best_graph_path) || !file.exists(summary_path)) {
  stop("Missing required sampled-final-rounds CSV outputs")
}

overall <- utils::read.csv(overall_path, stringsAsFactors = FALSE)
best_graph <- utils::read.csv(best_graph_path, stringsAsFactors = FALSE)
summary_df <- utils::read.csv(summary_path, stringsAsFactors = FALSE)

family_levels_sorted <- sort(unique(best_graph$family))
best_family <- do.call(
  rbind,
  lapply(family_levels_sorted, function(fam) {
    fam_df <- best_graph[best_graph$family == fam, , drop = FALSE]
    counts <- table(fam_df$final_rounds)
    mode_vals <- as.integer(names(counts)[counts == max(counts)])
    data.frame(
      family = fam,
      median_best_final_rounds = stats::median(fam_df$final_rounds),
      mean_best_final_rounds = mean(fam_df$final_rounds),
      modal_best_final_rounds = paste(mode_vals, collapse = "/"),
      modal_count = max(counts),
      graphs = nrow(fam_df),
      best_values = paste(fam_df$final_rounds[order(fam_df$graph_label)], collapse = ","),
      stringsAsFactors = FALSE
    )
  })
)
detailed_family_path <- file.path(root, "sampled-final-rounds-best-by-family-detailed.csv")
utils::write.csv(best_family, detailed_family_path, row.names = FALSE)

overall_png <- file.path(root, "sampled-final-rounds-overall-curves.png")
family_png <- file.path(root, "sampled-final-rounds-family-modal-best.png")
heatmap_png <- file.path(root, "sampled-final-rounds-family-heatmap.png")

grDevices::png(overall_png, width = 1800, height = 1200, res = 170, bg = "#f7f3ea")
graphics::par(mfrow = c(2, 1), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
graphics::plot(
  overall$final_rounds,
  overall$geometry_score_mean,
  type = "b",
  pch = 19,
  lwd = 2,
  col = "#16324f",
  xlab = "sampled final_rounds",
  ylab = "mean rank score",
  main = "Global mean scores across sampled_f0..32"
)
graphics::lines(
  overall$final_rounds,
  overall$full_score_mean,
  type = "b",
  pch = 17,
  lwd = 2,
  col = "#b23a48"
)
graphics::legend(
  "topright",
  legend = c("geometry score", "full score"),
  col = c("#16324f", "#b23a48"),
  pch = c(19, 17),
  bty = "n"
)

graphics::plot(
  overall$final_rounds,
  overall$global_symmetry_score_mean,
  type = "b",
  pch = 19,
  lwd = 2,
  col = "#355070",
  xlab = "sampled final_rounds",
  ylab = "global mean metric",
  ylim = range(
    c(
      overall$global_symmetry_score_mean,
      overall$local_angle_deviation_mean,
      overall$edge_axis_deviation_mean
    ),
    na.rm = TRUE
  ),
  main = "Global mean component diagnostics"
)
graphics::lines(
  overall$final_rounds,
  overall$local_angle_deviation_mean,
  type = "b",
  pch = 17,
  lwd = 2,
  col = "#b23a48"
)
graphics::lines(
  overall$final_rounds,
  overall$edge_axis_deviation_mean,
  type = "b",
  pch = 15,
  lwd = 2,
  col = "#3d7a57"
)
graphics::legend(
  "topright",
  legend = c("symmetry", "angle deviation", "axis deviation"),
  col = c("#355070", "#b23a48", "#3d7a57"),
  pch = c(19, 17, 15),
  bty = "n"
)
graphics::mtext("Sampled final-rounds suite", side = 3, outer = TRUE, cex = 1.2, font = 2)
grDevices::dev.off()

best_family <- best_family[order(best_family$median_best_final_rounds, best_family$mean_best_final_rounds), , drop = FALSE]
grDevices::png(family_png, width = 1800, height = 900, res = 170, bg = "#f7f3ea")
graphics::par(mar = c(8, 4, 3, 1))
bar_cols <- grDevices::colorRampPalette(c("#c9d6ea", "#355070"))(nrow(best_family))
graphics::barplot(
  best_family$median_best_final_rounds,
  names.arg = best_family$family,
  las = 2,
  col = bar_cols,
  border = NA,
  ylim = c(0, max(best_family$median_best_final_rounds) + 4),
  ylab = "median best sampled final_rounds",
  main = "Family-level median winners in the sampled_f0..32 sweep"
)
graphics::abline(h = c(8, 16, 24, 32), col = "#d8cfc2", lty = 3)
grDevices::dev.off()

family_levels <- unique(summary_df$family)
final_levels <- sort(unique(summary_df$final_rounds))
mat <- matrix(NA_real_, nrow = length(family_levels), ncol = length(final_levels))
rownames(mat) <- family_levels
colnames(mat) <- final_levels
for (i in seq_along(family_levels)) {
  fam <- family_levels[[i]]
  fam_df <- summary_df[summary_df$family == fam, , drop = FALSE]
  fam_mean <- stats::aggregate(full_score ~ final_rounds, data = fam_df, FUN = mean)
  mat[i, match(fam_mean$final_rounds, final_levels)] <- fam_mean$full_score
}

grDevices::png(heatmap_png, width = 2000, height = 1000, res = 170, bg = "#f7f3ea")
graphics::par(mar = c(6, 9, 3, 2))
palette_vals <- grDevices::colorRampPalette(c("#1f1f1f", "#355070", "#9ec1cf", "#f7f3ea"))(100)
filled <- mat
filled[!is.finite(filled)] <- max(filled, na.rm = TRUE)
heatmap_z <- t(filled[rev(seq_len(nrow(filled))), , drop = FALSE])
graphics::image(
  x = seq_along(final_levels),
  y = seq_along(family_levels),
  z = heatmap_z,
  col = rev(palette_vals),
  axes = FALSE,
  xlab = "sampled final_rounds",
  ylab = "",
  main = "Family mean full-score heatmap across sampled_f0..32"
)
graphics::axis(1, at = seq_along(final_levels), labels = final_levels, las = 2, cex.axis = 0.7)
graphics::axis(2, at = seq_along(family_levels), labels = rev(family_levels), las = 2, cex.axis = 0.9)
grDevices::dev.off()

cat(sprintf("Wrote %s\n", overall_png))
cat(sprintf("Wrote %s\n", family_png))
cat(sprintf("Wrote %s\n", heatmap_png))
cat(sprintf("Wrote %s\n", detailed_family_path))
