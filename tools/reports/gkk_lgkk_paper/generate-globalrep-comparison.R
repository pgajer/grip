#!/usr/bin/env Rscript

sys.source(file.path("tools", "benchmarks", "gkk_lgkk_paper", "benchmark-sierpinski-baseline.R"), envir = environment())

parse_named_args <- function(args) {
  out <- list()
  for (arg in args) {
    parts <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) {
      stop("Arguments must use key=value format")
    }
    out[[parts[[1L]]]] <- parts[[2L]]
  }
  out
}

parse_int_scalar <- function(x, name) {
  val <- suppressWarnings(as.integer(x))
  if (length(val) != 1L || is.na(val)) {
    stop(sprintf("%s must be a single integer", name))
  }
  val
}

plot_triptych <- function(canonical_coords,
                          baseline_coords,
                          globalrep_coords,
                          edges,
                          title_text,
                          subtitle_text) {
  can_norm <- normalize_coords(canonical_coords)
  base_norm <- normalize_coords(baseline_coords)
  global_norm <- normalize_coords(globalrep_coords)
  xs <- c(can_norm[, 1L], base_norm[, 1L], global_norm[, 1L])
  ys <- c(can_norm[, 2L], base_norm[, 2L], global_norm[, 2L])
  xr <- range(xs)
  yr <- range(ys)
  xpad <- max(diff(xr) * 0.12, 0.12)
  ypad <- max(diff(yr) * 0.12, 0.12)
  xlim <- xr + c(-xpad, xpad)
  ylim <- yr + c(-ypad, ypad)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, 3), mar = c(0, 0, 2.8, 0), xaxs = "i", yaxs = "i")
  plot_layout_panel(canonical_coords, edges, paste(title_text, "- canonical"), subtitle_text,
                    xlim = xlim, ylim = ylim)
  plot_layout_panel(baseline_coords, edges, paste(title_text, "- baseline"), subtitle_text,
                    xlim = xlim, ylim = ylim)
  plot_layout_panel(globalrep_coords, edges, paste(title_text, "- globalrep"), subtitle_text,
                    xlim = xlim, ylim = ylim)
}

write_triptych_outputs <- function(pdf_path,
                                   png_path,
                                   canonical_coords,
                                   baseline_coords,
                                   globalrep_coords,
                                   edges,
                                   title_text,
                                   subtitle_text) {
  grDevices::pdf(pdf_path, width = 21, height = 8.5,
                 paper = "special", bg = "#f7f3ea", useDingbats = FALSE)
  plot_triptych(canonical_coords, baseline_coords, globalrep_coords, edges,
                title_text, subtitle_text)
  grDevices::dev.off()

  grDevices::png(png_path, width = 3024, height = 1224, res = 144,
                 bg = "#f7f3ea")
  plot_triptych(canonical_coords, baseline_coords, globalrep_coords, edges,
                title_text, subtitle_text)
  grDevices::dev.off()
}

compute_metrics <- function(coords, canonical_coords, edges, seed) {
  adj <- make_adj_list(edges, nrow(coords))
  aligned <- align_to_target(coords, canonical_coords)
  edge_stats <- edge_length_stats(coords, edges)
  list(
    rmse = aligned$rmse,
    edge_cv = edge_stats$cv,
    stress = sampled_stress(coords, adj, sample_size = 2000L, rng_seed = 1000L + seed),
    sep = sampled_nonedge_separation_ratio(coords, edges, sample_size = 5000L, rng_seed = 2000L + seed)
  )
}

write_summary <- function(path,
                          family,
                          level,
                          seed,
                          baseline_metrics,
                          globalrep_metrics,
                          pdf_path,
                          png_path) {
  lines <- c(
    sprintf("# %s level %d default-vs-globalrep comparison", tools::toTitleCase(family), level),
    "",
    sprintf("- seed: `%d`", seed),
    "- baseline call: `grip.layout.legacy(edges, n = n, dim = 2)`",
    "- current call: `grip.layout(edges, n = n, dim = 2)`",
    "",
    "## Metrics",
    "",
    "| Layout | Procrustes RMSE | Edge CV | Sampled stress | Non-edge sep ratio |",
    "| --- | ---: | ---: | ---: | ---: |",
    sprintf("| baseline | %.4f | %.4f | %.4f | %.4f |",
            baseline_metrics$rmse,
            baseline_metrics$edge_cv,
            baseline_metrics$stress,
            baseline_metrics$sep),
    sprintf("| current | %.4f | %.4f | %.4f | %.4f |",
            globalrep_metrics$rmse,
            globalrep_metrics$edge_cv,
            globalrep_metrics$stress,
            globalrep_metrics$sep),
    "",
    "## Outputs",
    "",
    sprintf("- PDF: `%s`", pdf_path),
    sprintf("- PNG: `%s`", png_path)
  )
  writeLines(lines, con = path)
}

args <- parse_named_args(commandArgs(trailingOnly = TRUE))

family <- if (!is.null(args$family)) args$family else "carpet"
level <- if (!is.null(args$level)) parse_int_scalar(args$level, "level") else 3L
seed <- if (!is.null(args$seed)) parse_int_scalar(args$seed, "seed") else 6L
tag <- if (!is.null(args$tag)) args$tag else sprintf("globalrep-%s-level%d-default", family, level)

if (!family %in% c("triangle", "carpet")) {
  stop("family must be one of: triangle, carpet")
}

if (identical(family, "triangle")) {
  built <- build_sierpinski_triangle(level)
  edges <- edges.sierpinski.triangle(level)
} else {
  built <- build_sierpinski_carpet(level)
  edges <- edges.sierpinski.carpet(level)
}

if (!identical(unname(built$edges), unname(edges))) {
  stop("Canonical builder does not match package edge generator")
}

canonical_coords <- built$coords
n <- max(edges)

baseline_coords <- grip.layout.legacy(
  edges = edges,
  n = n,
  dim = 2,
  seed = seed
)

globalrep_coords <- grip.layout(
  edges = edges,
  n = n,
  dim = 2,
  seed = seed
)

baseline_metrics <- compute_metrics(baseline_coords, canonical_coords, edges, seed)
globalrep_metrics <- compute_metrics(globalrep_coords, canonical_coords, edges, seed)

pdf_dir_local <- file.path("dev", "manual", "pdf", tag, family)
tmp_dir_local <- file.path("dev", "manual", "tmp", tag)
preview_dir_local <- file.path(tmp_dir_local, "pdf-previews")
dir.create(pdf_dir_local, recursive = TRUE, showWarnings = FALSE)
dir.create(preview_dir_local, recursive = TRUE, showWarnings = FALSE)

base_name <- sprintf("sierpinski-%s-level-%d-canonical-baseline-globalrep", family, level)
pdf_path <- file.path(pdf_dir_local, paste0(base_name, ".pdf"))
png_path <- file.path(preview_dir_local, paste0(base_name, ".png"))
summary_path <- file.path(tmp_dir_local, sprintf("%s-summary.md", base_name))

subtitle <- sprintf(
  "seed=%d; baseline=grip.layout.legacy defaults; current=grip.layout defaults",
  seed
)
title_text <- sprintf("Sierpinski %s level %d", family, level)

write_triptych_outputs(
  pdf_path = pdf_path,
  png_path = png_path,
  canonical_coords = canonical_coords,
  baseline_coords = baseline_coords,
  globalrep_coords = globalrep_coords,
  edges = edges,
  title_text = title_text,
  subtitle_text = subtitle
)

write_summary(
  path = summary_path,
  family = family,
  level = level,
  seed = seed,
  baseline_metrics = baseline_metrics,
  globalrep_metrics = globalrep_metrics,
  pdf_path = pdf_path,
  png_path = png_path
)

message(sprintf("Summary written to %s", summary_path))
message(sprintf("PDF written to %s", pdf_path))
message(sprintf("PNG written to %s", png_path))
