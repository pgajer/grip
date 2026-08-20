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

parse_int_vector <- function(x, name) {
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  val <- suppressWarnings(as.integer(parts))
  if (length(val) == 0L || any(is.na(val))) {
    stop(sprintf("%s must be a comma-separated integer list", name))
  }
  val
}

format_num <- function(x, digits = 4L) {
  if (!is.finite(x)) {
    return("NA")
  }
  format(round(x, digits), nsmall = digits, trim = TRUE)
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

best_candidate_cfg <- list(
  placement = "barycenter",
  rounds = 160L,
  num_init = 24L,
  num_nbrs = 20L,
  r = 0.03,
  s = 7.5,
  repulsion_factor = 2.5,
  coarse_repulsion_factor = 1.5,
  coarse_repulsion_sample = 16L,
  coarse_repulsion_exact_below = 64L
)

run_layout <- function(edges, n, seed, final_rounds) {
  grip.layout(
    edges = edges,
    n = n,
    dim = 2,
    placement = best_candidate_cfg$placement,
    rounds = best_candidate_cfg$rounds,
    final_rounds = as.integer(final_rounds),
    num_init = best_candidate_cfg$num_init,
    num_nbrs = best_candidate_cfg$num_nbrs,
    r = best_candidate_cfg$r,
    s = best_candidate_cfg$s,
    repulsion_factor = best_candidate_cfg$repulsion_factor,
    coarse_repulsion_factor = best_candidate_cfg$coarse_repulsion_factor,
    coarse_repulsion_sample = best_candidate_cfg$coarse_repulsion_sample,
    coarse_repulsion_exact_below = best_candidate_cfg$coarse_repulsion_exact_below,
    seed = as.integer(seed)
  )
}

plot_seed_grid <- function(results,
                           edges,
                           title_text,
                           subtitle_text) {
  all_coords <- lapply(results$coords, normalize_coords)
  xs <- unlist(lapply(all_coords, function(x) x[, 1L]), use.names = FALSE)
  ys <- unlist(lapply(all_coords, function(x) x[, 2L]), use.names = FALSE)
  xr <- range(xs)
  yr <- range(ys)
  xpad <- max(diff(xr) * 0.12, 0.12)
  ypad <- max(diff(yr) * 0.12, 0.12)
  xlim <- xr + c(-xpad, xpad)
  ylim <- yr + c(-ypad, ypad)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(length(unique(results$seed)), length(unique(results$variant_label))),
    mar = c(0, 0, 2.3, 0),
    oma = c(0.5, 0, 3.2, 0),
    xaxs = "i",
    yaxs = "i"
  )

  for (i in seq_len(nrow(results))) {
    row <- results[i, , drop = FALSE]
    metric_text <- sprintf(
      "RMSE=%s | sep=%s",
      format_num(row$rmse[[1L]], digits = 3L),
      format_num(row$sep[[1L]], digits = 3L)
    )
    plot_layout_panel(
      row$coords[[1L]],
      edges,
      sprintf("seed=%d | final_rounds=%d", row$seed[[1L]], row$final_rounds[[1L]]),
      metric_text,
      xlim = xlim,
      ylim = ylim
    )
  }

  graphics::mtext(title_text, outer = TRUE, side = 3, line = 1.8,
                  cex = 1.3, font = 2, col = "#16324f")
  graphics::mtext(subtitle_text, outer = TRUE, side = 3, line = 0.4,
                  cex = 0.95, col = "#466074")
}

write_outputs <- function(png_path,
                          pdf_path,
                          results,
                          edges,
                          title_text,
                          subtitle_text) {
  grDevices::pdf(pdf_path, width = 15, height = 22,
                 paper = "special", bg = "#f7f3ea", useDingbats = FALSE)
  plot_seed_grid(results, edges, title_text, subtitle_text)
  grDevices::dev.off()

  grDevices::png(png_path, width = 2160, height = 3168, res = 144,
                 bg = "#f7f3ea")
  plot_seed_grid(results, edges, title_text, subtitle_text)
  grDevices::dev.off()
}

write_summary <- function(path,
                          level,
                          seeds,
                          baseline_final_rounds,
                          variant_final_rounds,
                          results,
                          png_path,
                          pdf_path,
                          csv_path) {
  summary_rows <- do.call(rbind, lapply(split(results, results$variant_label), function(df) {
    data.frame(
      variant_label = df$variant_label[[1L]],
      final_rounds = df$final_rounds[[1L]],
      rmse_mean = mean(df$rmse),
      edge_cv_mean = mean(df$edge_cv),
      stress_mean = mean(df$stress),
      sep_mean = mean(df$sep),
      stringsAsFactors = FALSE
    )
  }))
  summary_rows <- summary_rows[match(c("best_384", "variant_200"), summary_rows$variant_label), , drop = FALSE]

  lines <- c(
    sprintf("# Sierpinski carpet level %d final-rounds seed comparison", level),
    "",
    sprintf("- seeds: `%s`", paste(seeds, collapse = ", ")),
    sprintf("- left column: best tuned candidate with `final_rounds = %d`", baseline_final_rounds),
    sprintf("- right column: same candidate with `final_rounds = %d`", variant_final_rounds),
    sprintf(
      "- shared parameters: `placement=%s`, `rounds=%d`, `num_init=%d`, `num_nbrs=%d`, `r=%.2f`, `s=%.1f`, `repulsion_factor=%.1f`, `coarse_repulsion_factor=%.1f`, `coarse_repulsion_sample=%d`, `coarse_repulsion_exact_below=%d`",
      best_candidate_cfg$placement,
      best_candidate_cfg$rounds,
      best_candidate_cfg$num_init,
      best_candidate_cfg$num_nbrs,
      best_candidate_cfg$r,
      best_candidate_cfg$s,
      best_candidate_cfg$repulsion_factor,
      best_candidate_cfg$coarse_repulsion_factor,
      best_candidate_cfg$coarse_repulsion_sample,
      best_candidate_cfg$coarse_repulsion_exact_below
    ),
    "",
    "## Mean metrics across seeds",
    "",
    "| Variant | final_rounds | RMSE | Edge CV | Sampled stress | Non-edge sep ratio |",
    "| --- | ---: | ---: | ---: | ---: | ---: |",
    sprintf(
      "| best_384 | %d | %.4f | %.4f | %.4f | %.4f |",
      summary_rows$final_rounds[[1L]],
      summary_rows$rmse_mean[[1L]],
      summary_rows$edge_cv_mean[[1L]],
      summary_rows$stress_mean[[1L]],
      summary_rows$sep_mean[[1L]]
    ),
    sprintf(
      "| variant_200 | %d | %.4f | %.4f | %.4f | %.4f |",
      summary_rows$final_rounds[[2L]],
      summary_rows$rmse_mean[[2L]],
      summary_rows$edge_cv_mean[[2L]],
      summary_rows$stress_mean[[2L]],
      summary_rows$sep_mean[[2L]]
    ),
    "",
    "## Outputs",
    "",
    sprintf("- PNG: `%s`", png_path),
    sprintf("- PDF: `%s`", pdf_path),
    sprintf("- per-seed metrics CSV: `%s`", csv_path)
  )
  writeLines(lines, con = path)
}

args <- parse_named_args(commandArgs(trailingOnly = TRUE))

level <- if (!is.null(args$level)) parse_int_scalar(args$level, "level") else 3L
seeds <- if (!is.null(args$seeds)) parse_int_vector(args$seeds, "seeds") else 1:6
baseline_final_rounds <- if (!is.null(args$baseline_final_rounds)) parse_int_scalar(args$baseline_final_rounds, "baseline_final_rounds") else 384L
variant_final_rounds <- if (!is.null(args$variant_final_rounds)) parse_int_scalar(args$variant_final_rounds, "variant_final_rounds") else 200L
tag <- if (!is.null(args$tag)) args$tag else sprintf(
  "globalrep-carpet-level%d-finalrounds-%d-vs-%d-seeds-%s",
  level,
  baseline_final_rounds,
  variant_final_rounds,
  paste(seeds, collapse = "-")
)

built <- build_sierpinski_carpet(level)
edges <- edges.sierpinski.carpet(level)
if (!identical(unname(built$edges), unname(edges))) {
  stop("Canonical builder does not match package edge generator")
}

canonical_coords <- built$coords
n <- max(edges)

variants <- data.frame(
  variant_label = c("best_384", "variant_200"),
  final_rounds = c(baseline_final_rounds, variant_final_rounds),
  stringsAsFactors = FALSE
)

results <- do.call(rbind, lapply(seeds, function(seed) {
  do.call(rbind, lapply(seq_len(nrow(variants)), function(i) {
    variant <- variants[i, , drop = FALSE]
    coords <- run_layout(edges, n, seed = seed, final_rounds = variant$final_rounds[[1L]])
    metrics <- compute_metrics(coords, canonical_coords, edges, seed)
    data.frame(
      seed = as.integer(seed),
      variant_label = variant$variant_label[[1L]],
      final_rounds = as.integer(variant$final_rounds[[1L]]),
      rmse = metrics$rmse,
      edge_cv = metrics$edge_cv,
      stress = metrics$stress,
      sep = metrics$sep,
      stringsAsFactors = FALSE
    )
  }))
}))

coords_list <- vector("list", nrow(results))
for (i in seq_len(nrow(results))) {
  coords_list[[i]] <- run_layout(
    edges,
    n,
    seed = results$seed[[i]],
    final_rounds = results$final_rounds[[i]]
  )
}
results$coords <- coords_list

results <- results[order(results$seed, match(results$variant_label, c("best_384", "variant_200"))), , drop = FALSE]

tmp_dir_local <- file.path("output", "gkk_lgkk_paper", "tmp", tag)
preview_dir_local <- file.path(tmp_dir_local, "pdf-previews")
pdf_dir_local <- file.path("output", "gkk_lgkk_paper", "reports", tag, "carpet")
dir.create(preview_dir_local, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir_local, recursive = TRUE, showWarnings = FALSE)

base_name <- sprintf(
  "sierpinski-carpet-level-%d-globalrep-finalrounds-%d-vs-%d-seed-grid",
  level,
  baseline_final_rounds,
  variant_final_rounds
)
png_path <- file.path(preview_dir_local, paste0(base_name, ".png"))
pdf_path <- file.path(pdf_dir_local, paste0(base_name, ".pdf"))
csv_path <- file.path(tmp_dir_local, paste0(base_name, "-metrics.csv"))
summary_path <- file.path(tmp_dir_local, paste0(base_name, "-summary.md"))

title_text <- sprintf("Sierpinski carpet level %d globalrep comparison", level)
subtitle_text <- sprintf(
  "Left: best tuned final_rounds=%d | Right: same parameters with final_rounds=%d",
  baseline_final_rounds,
  variant_final_rounds
)

write_outputs(
  png_path = png_path,
  pdf_path = pdf_path,
  results = results,
  edges = edges,
  title_text = title_text,
  subtitle_text = subtitle_text
)

results_out <- results
results_out$coords <- NULL
utils::write.csv(results_out, csv_path, row.names = FALSE)
write_summary(
  path = summary_path,
  level = level,
  seeds = seeds,
  baseline_final_rounds = baseline_final_rounds,
  variant_final_rounds = variant_final_rounds,
  results = results_out,
  png_path = png_path,
  pdf_path = pdf_path,
  csv_path = csv_path
)

message(sprintf("Summary written to %s", summary_path))
message(sprintf("CSV written to %s", csv_path))
message(sprintf("PDF written to %s", pdf_path))
message(sprintf("PNG written to %s", png_path))
