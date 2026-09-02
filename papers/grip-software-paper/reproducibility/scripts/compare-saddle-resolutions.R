#!/usr/bin/env Rscript
# Generate both requested resolutions, all scores, and matched figure exports.
# Usage: Rscript scripts/compare-saddle-resolutions.R OUTPUT_DIRECTORY [INPUT.rds OUTPUT.rds]
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
paper_library <- Sys.getenv("GRIP_RJOURNAL_PACKAGE_LIBRARY")
if (nzchar(paper_library)) .libPaths(c(paper_library, .libPaths()))
source(file.path(script_dir, "weighted-saddle-comparison.R"))
source(file.path(script_dir, "plot-weighted-saddle.R"))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(1L, 3L)) {
  stop("Provide OUTPUT_DIRECTORY, optionally followed by benchmark INPUT.rds OUTPUT.rds")
}
dir.create(args[[1L]], recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(args[[1L]], mustWork = TRUE)
cases <- setNames(lapply(c(10L, 15L), function(size) {
  cat("Generating", size, "x", size, "comparison...\n")
  weighted_saddle_comparison(grid_size = size)
}), c("10x10", "15x15"))
limits <- weighted_saddle_limits(cases)
for (name in names(cases)) {
  saddle <- cases[[name]]
  stem <- file.path(out_dir, paste0("weighted-saddle-", name))
  saveRDS(saddle, paste0(stem, ".rds"), compress = "xz")
  pdf(paste0(stem, ".pdf"), width = 11.2, height = 6.8, useDingbats = FALSE)
  plot_weighted_saddle(saddle, limits)
  dev.off()
  png(paste0(stem, ".png"), width = 11.2, height = 6.8, units = "in", res = 180)
  plot_weighted_saddle(saddle, limits)
  dev.off()
  print(saddle$scores[, c("method", "edge.rel.rmse", "gmds.stress",
                          "gkk.weighted.rel.rmse", "gkk.mean.rel.path.error")])
  print(saddle$validation)
}
scores <- do.call(rbind, lapply(names(cases), function(name) {
  cbind(mesh = name, cases[[name]]$scores)
}))
write.csv(scores, file.path(out_dir, "weighted-saddle-resolution-scores.csv"), row.names = FALSE)
panel_scores <- do.call(rbind, lapply(names(cases), function(name) {
  cbind(mesh = name, weighted_saddle_panel_scores(cases[[name]]))
}))
write.csv(panel_scores, file.path(out_dir, "weighted-saddle-panel-scores.csv"), row.names = FALSE)
comparison <- list(cases = cases, display_limits = limits, selected_grid = "10x10")
saveRDS(comparison,
        file.path(out_dir, "weighted-saddle-resolutions.rds"), compress = "xz")
if (length(args) == 3L) {
  original <- readRDS(args[[2L]])
  updated <- original
  updated$weighted_saddle <- cases[[comparison$selected_grid]]
  updated$weighted_saddle_resolutions <- comparison
  untouched <- setdiff(names(original), c("weighted_saddle", "weighted_saddle_resolutions"))
  stopifnot(identical(original[untouched], updated[untouched]))
  saveRDS(updated, args[[3L]], compress = "xz")
  stopifnot(identical(updated, readRDS(args[[3L]])))
  cat("Updated saddle components only; other benchmark objects are unchanged.\n")
}
cat("Figures, scores, validation, and coordinates saved to", out_dir, "\n")
