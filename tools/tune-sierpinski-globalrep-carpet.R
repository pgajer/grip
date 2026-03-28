#!/usr/bin/env Rscript

sys.source(file.path("tools", "benchmark-sierpinski-baseline.R"), envir = environment())

score_weights <- c(
  procrustes_rmse = 0.45,
  edge_length_cv = 0.20,
  sampled_stress = 0.20,
  sampled_nonedge_sep_ratio = 0.15
)

candidate_fields <- c(
  "placement", "rounds", "final_rounds", "num_init",
  "num_nbrs", "r", "s", "repulsion_factor",
  "coarse_repulsion_factor", "coarse_repulsion_sample",
  "coarse_repulsion_exact_below"
)

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

parse_num_vector <- function(x, name) {
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  val <- suppressWarnings(as.double(parts))
  if (length(val) == 0L || any(is.na(val) | !is.finite(val))) {
    stop(sprintf("%s must be a comma-separated numeric list", name))
  }
  val
}

validate_run_tag <- function(x) {
  if (!grepl("^[A-Za-z0-9._-]+$", x)) {
    stop("tag must match ^[A-Za-z0-9._-]+$")
  }
  x
}

search_space <- list(
  placement = c("barycenter"),
  rounds = c(128L, 160L, 192L),
  final_rounds = c(256L, 320L, 384L),
  num_init = c(20L, 24L, 28L, 32L),
  num_nbrs = c(16L, 20L, 24L, 28L),
  r = c(0.02, 0.03, 0.05, 0.07),
  s = c(3.0, 4.5, 6.0, 7.5),
  repulsion_factor = c(2.0, 2.25, 2.5, 2.75, 3.0),
  coarse_repulsion_factor = c(0.0, 0.1, 0.2, 0.4, 0.8, 1.5, 3.0, 6.0),
  coarse_repulsion_sample = c(8L, 16L, 32L, 64L),
  coarse_repulsion_exact_below = c(64L, 96L, 128L, 192L, 256L)
)

standard_default_cfg <- list(
  placement = "barycenter",
  rounds = 20L,
  final_rounds = 25L,
  num_init = 36L,
  num_nbrs = 10L,
  r = 0.15,
  s = 3.0,
  repulsion_factor = 1.0
)

globalrep_default_cfg <- c(
  standard_default_cfg,
  list(
    coarse_repulsion_factor = 0.2,
    coarse_repulsion_sample = 16L,
    coarse_repulsion_exact_below = 128L
  )
)

anchor_cfg <- list(
  placement = "barycenter",
  rounds = 160L,
  final_rounds = 320L,
  num_init = 24L,
  num_nbrs = 20L,
  r = 0.03,
  s = 4.5,
  repulsion_factor = 2.5,
  coarse_repulsion_factor = 0.2,
  coarse_repulsion_sample = 16L,
  coarse_repulsion_exact_below = 128L
)

candidate_key <- function(candidate) {
  paste(
    candidate$placement,
    candidate$rounds,
    candidate$final_rounds,
    candidate$num_init,
    candidate$num_nbrs,
    sprintf("%.4f", candidate$r),
    sprintf("%.4f", candidate$s),
    sprintf("%.4f", candidate$repulsion_factor),
    sprintf("%.4f", candidate$coarse_repulsion_factor),
    candidate$coarse_repulsion_sample,
    candidate$coarse_repulsion_exact_below,
    sep = "|"
  )
}

candidate_row <- function(candidate) {
  out <- as.data.frame(candidate[c("candidate_id", "candidate_source", candidate_fields)],
                       stringsAsFactors = FALSE)
  out$rounds <- as.integer(out$rounds)
  out$final_rounds <- as.integer(out$final_rounds)
  out$num_init <- as.integer(out$num_init)
  out$num_nbrs <- as.integer(out$num_nbrs)
  out$r <- as.double(out$r)
  out$s <- as.double(out$s)
  out$repulsion_factor <- as.double(out$repulsion_factor)
  out$coarse_repulsion_factor <- as.double(out$coarse_repulsion_factor)
  out$coarse_repulsion_sample <- as.integer(out$coarse_repulsion_sample)
  out$coarse_repulsion_exact_below <- as.integer(out$coarse_repulsion_exact_below)
  out
}

format_candidate <- function(row) {
  sprintf(
    paste(
      "placement=%s, rounds=%d, final_rounds=%d, num_init=%d, num_nbrs=%d,",
      "r=%.2f, s=%.1f, repulsion_factor=%.2f, coarse_repulsion_factor=%.2f,",
      "coarse_repulsion_sample=%d, coarse_repulsion_exact_below=%d"
    ),
    row$placement,
    row$rounds,
    row$final_rounds,
    row$num_init,
    row$num_nbrs,
    row$r,
    row$s,
    row$repulsion_factor,
    row$coarse_repulsion_factor,
    row$coarse_repulsion_sample,
    row$coarse_repulsion_exact_below
  )
}

sample_candidate <- function(space) {
  rounds <- sample(space$rounds, 1L)
  final_choices <- space$final_rounds[space$final_rounds >= rounds]
  if (length(final_choices) == 0L) {
    final_choices <- space$final_rounds
  }
  list(
    placement = sample(space$placement, 1L),
    rounds = as.integer(rounds),
    final_rounds = as.integer(sample(final_choices, 1L)),
    num_init = as.integer(sample(space$num_init, 1L)),
    num_nbrs = as.integer(sample(space$num_nbrs, 1L)),
    r = as.double(sample(space$r, 1L)),
    s = as.double(sample(space$s, 1L)),
    repulsion_factor = as.double(sample(space$repulsion_factor, 1L)),
    coarse_repulsion_factor = as.double(sample(space$coarse_repulsion_factor, 1L)),
    coarse_repulsion_sample = as.integer(sample(space$coarse_repulsion_sample, 1L)),
    coarse_repulsion_exact_below = as.integer(sample(space$coarse_repulsion_exact_below, 1L))
  )
}

generate_candidates <- function(n_random, search_seed) {
  candidates <- list(
    c(list(candidate_id = "globalrep_default", candidate_source = "baseline"), globalrep_default_cfg),
    c(list(candidate_id = "globalrep_anchor", candidate_source = "anchor"), anchor_cfg)
  )
  seen <- vapply(candidates, candidate_key, character(1L))

  set.seed(search_seed)
  attempts <- 0L
  max_attempts <- max(1000L, n_random * 300L)
  while (length(candidates) < n_random + 2L && attempts < max_attempts) {
    attempts <- attempts + 1L
    sampled <- sample_candidate(search_space)
    key <- candidate_key(sampled)
    if (key %in% seen) {
      next
    }
    idx <- length(candidates) - 1L
    candidates[[length(candidates) + 1L]] <- c(
      list(candidate_id = sprintf("globalrep_rand_%03d", idx), candidate_source = "random"),
      sampled
    )
    seen <- c(seen, key)
  }

  if (length(candidates) < n_random + 2L) {
    warning(sprintf(
      "Requested %d random candidates but found only %d unique combinations before hitting the attempt limit.",
      n_random,
      length(candidates) - 2L
    ), call. = FALSE)
  }

  candidates
}

rank01 <- function(x, higher_better = FALSE) {
  n <- length(x)
  out <- rep(1, n)
  ok <- is.finite(x)
  if (!any(ok)) {
    return(out)
  }
  vals <- if (higher_better) -x[ok] else x[ok]
  if (length(vals) == 1L || all(abs(vals - vals[[1L]]) <= sqrt(.Machine$double.eps))) {
    out[ok] <- 0
    return(out)
  }
  ranks <- rank(vals, ties.method = "average")
  out[ok] <- (ranks - 1) / (length(ranks) - 1)
  out
}

format_num <- function(x, digits = 4) {
  ifelse(is.finite(x), format(round(x, digits), nsmall = digits, trim = TRUE), "NA")
}

run_candidate <- function(spec, candidate, seed) {
  n <- max(spec$edges)
  adj <- make_adj_list(spec$edges, n)
  started <- proc.time()[["elapsed"]]
  coords <- grip.layout(
    edges = spec$edges,
    n = n,
    dim = 2,
    placement = candidate$placement,
    rounds = candidate$rounds,
    final_rounds = candidate$final_rounds,
    num_init = candidate$num_init,
    num_nbrs = candidate$num_nbrs,
    r = candidate$r,
    s = candidate$s,
    repulsion_factor = candidate$repulsion_factor,
    coarse_repulsion_factor = candidate$coarse_repulsion_factor,
    coarse_repulsion_sample = candidate$coarse_repulsion_sample,
    coarse_repulsion_exact_below = candidate$coarse_repulsion_exact_below,
    seed = seed
  )
  elapsed <- proc.time()[["elapsed"]] - started
  aligned <- align_to_target(coords, spec$canonical)
  edge_stats <- edge_length_stats(coords, spec$edges)
  data.frame(
    candidate_id = candidate$candidate_id,
    candidate_source = candidate$candidate_source,
    seed = seed,
    placement = candidate$placement,
    rounds = candidate$rounds,
    final_rounds = candidate$final_rounds,
    num_init = candidate$num_init,
    num_nbrs = candidate$num_nbrs,
    r = candidate$r,
    s = candidate$s,
    repulsion_factor = candidate$repulsion_factor,
    coarse_repulsion_factor = candidate$coarse_repulsion_factor,
    coarse_repulsion_sample = candidate$coarse_repulsion_sample,
    coarse_repulsion_exact_below = candidate$coarse_repulsion_exact_below,
    procrustes_rmse = aligned$rmse,
    edge_length_cv = edge_stats$cv,
    sampled_stress = sampled_stress(coords, adj, sample_size = 2000L, rng_seed = 1000L + seed),
    sampled_nonedge_sep_ratio = sampled_nonedge_separation_ratio(
      coords, spec$edges, sample_size = 5000L, rng_seed = 2000L + seed
    ),
    elapsed_sec = elapsed,
    stringsAsFactors = FALSE
  )
}

summarize_candidates <- function(raw_metrics) {
  groups <- split(raw_metrics, raw_metrics$candidate_id)
  do.call(rbind, lapply(groups, function(df) {
    best_seed <- df$seed[[which.min(df$procrustes_rmse)]]
    data.frame(
      candidate_id = df$candidate_id[[1L]],
      candidate_source = df$candidate_source[[1L]],
      placement = df$placement[[1L]],
      rounds = df$rounds[[1L]],
      final_rounds = df$final_rounds[[1L]],
      num_init = df$num_init[[1L]],
      num_nbrs = df$num_nbrs[[1L]],
      r = df$r[[1L]],
      s = df$s[[1L]],
      repulsion_factor = df$repulsion_factor[[1L]],
      coarse_repulsion_factor = df$coarse_repulsion_factor[[1L]],
      coarse_repulsion_sample = df$coarse_repulsion_sample[[1L]],
      coarse_repulsion_exact_below = df$coarse_repulsion_exact_below[[1L]],
      n_runs = nrow(df),
      procrustes_rmse_mean = mean(df$procrustes_rmse),
      procrustes_rmse_sd = if (nrow(df) > 1L) stats::sd(df$procrustes_rmse) else 0,
      edge_length_cv_mean = mean(df$edge_length_cv),
      edge_length_cv_sd = if (nrow(df) > 1L) stats::sd(df$edge_length_cv) else 0,
      sampled_stress_mean = mean(df$sampled_stress),
      sampled_stress_sd = if (nrow(df) > 1L) stats::sd(df$sampled_stress) else 0,
      sampled_nonedge_sep_ratio_mean = mean(df$sampled_nonedge_sep_ratio),
      sampled_nonedge_sep_ratio_sd = if (nrow(df) > 1L) stats::sd(df$sampled_nonedge_sep_ratio) else 0,
      elapsed_sec_mean = mean(df$elapsed_sec),
      best_seed = best_seed,
      stringsAsFactors = FALSE
    )
  }))
}

score_candidates <- function(summary_metrics) {
  out <- summary_metrics
  out$rank_rmse <- rank01(out$procrustes_rmse_mean, higher_better = FALSE)
  out$rank_edge_cv <- rank01(out$edge_length_cv_mean, higher_better = FALSE)
  out$rank_stress <- rank01(out$sampled_stress_mean, higher_better = FALSE)
  out$rank_sep <- rank01(out$sampled_nonedge_sep_ratio_mean, higher_better = TRUE)
  out$score <- score_weights[["procrustes_rmse"]] * out$rank_rmse +
    score_weights[["edge_length_cv"]] * out$rank_edge_cv +
    score_weights[["sampled_stress"]] * out$rank_stress +
    score_weights[["sampled_nonedge_sep_ratio"]] * out$rank_sep
  out[order(out$score, out$procrustes_rmse_mean), , drop = FALSE]
}

plot_triptych <- function(canonical_coords,
                          baseline_coords,
                          tuned_coords,
                          edges,
                          title_text,
                          subtitle_text) {
  can_norm <- normalize_coords(canonical_coords)
  base_norm <- normalize_coords(baseline_coords)
  tuned_norm <- normalize_coords(tuned_coords)
  xs <- c(can_norm[, 1L], base_norm[, 1L], tuned_norm[, 1L])
  ys <- c(can_norm[, 2L], base_norm[, 2L], tuned_norm[, 2L])
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
  plot_layout_panel(baseline_coords, edges, paste(title_text, "- legacy default"), subtitle_text,
                    xlim = xlim, ylim = ylim)
  plot_layout_panel(tuned_coords, edges, paste(title_text, "- current tuned"), subtitle_text,
                    xlim = xlim, ylim = ylim)
}

write_triptych_outputs <- function(pdf_path,
                                   png_path,
                                   canonical_coords,
                                   baseline_coords,
                                   tuned_coords,
                                   edges,
                                   title_text,
                                   subtitle_text) {
  grDevices::pdf(pdf_path, width = 21, height = 8.5,
                 paper = "special", bg = "#f7f3ea", useDingbats = FALSE)
  plot_triptych(canonical_coords, baseline_coords, tuned_coords, edges,
                title_text, subtitle_text)
  grDevices::dev.off()

  grDevices::png(png_path, width = 3024, height = 1224, res = 144,
                 bg = "#f7f3ea")
  plot_triptych(canonical_coords, baseline_coords, tuned_coords, edges,
                title_text, subtitle_text)
  grDevices::dev.off()
}

write_summary_markdown <- function(path,
                                   config,
                                   family_rankings,
                                   top_n,
                                   standard_default_metrics,
                                   tuned_metrics,
                                   tuned_candidate,
                                   comparison_pdf,
                                   comparison_png) {
  top_rows <- head(family_rankings, top_n)
  lines <- c(
    "# Sierpinski Global-Repulsion Carpet Tuning",
    "",
    "Search setup:",
    sprintf("- level: `%d`", config$level),
    sprintf("- random candidates: `%d`", config$n_random),
    sprintf("- seeds: `%s`", paste(config$seeds, collapse = ", ")),
    sprintf("- comparison seed: `%d`", config$comparison_seed),
    sprintf("- search seed: `%d`", config$search_seed),
    "",
    "Score definition (lower is better):",
    sprintf(
      "- `score = %.2f * rank(RMSE) + %.2f * rank(edge_length_cv) + %.2f * rank(sampled_stress) + %.2f * rank(-nonedge_sep)`",
      score_weights[["procrustes_rmse"]],
      score_weights[["edge_length_cv"]],
      score_weights[["sampled_stress"]],
      score_weights[["sampled_nonedge_sep_ratio"]]
    ),
    "",
    "Search space:",
    sprintf("- rounds: `%s`", paste(search_space$rounds, collapse = ", ")),
    sprintf("- final_rounds: `%s`", paste(search_space$final_rounds, collapse = ", ")),
    sprintf("- num_init: `%s`", paste(search_space$num_init, collapse = ", ")),
    sprintf("- num_nbrs: `%s`", paste(search_space$num_nbrs, collapse = ", ")),
    sprintf("- r: `%s`", paste(search_space$r, collapse = ", ")),
    sprintf("- s: `%s`", paste(search_space$s, collapse = ", ")),
    sprintf("- repulsion_factor: `%s`", paste(search_space$repulsion_factor, collapse = ", ")),
    sprintf("- coarse_repulsion_factor: `%s`", paste(search_space$coarse_repulsion_factor, collapse = ", ")),
    sprintf("- coarse_repulsion_sample: `%s`", paste(search_space$coarse_repulsion_sample, collapse = ", ")),
    sprintf("- coarse_repulsion_exact_below: `%s`", paste(search_space$coarse_repulsion_exact_below, collapse = ", ")),
    "",
    "Best candidate:",
    sprintf("- `%s`", tuned_candidate$candidate_id[[1L]]),
    sprintf("- %s", format_candidate(tuned_candidate[1L, ])),
    sprintf("- mean score across seeds: %s", format_num(tuned_candidate$score[[1L]], 4)),
    "",
    "Top candidates:",
    "",
    "| Candidate | Source | Mean score | RMSE | Edge CV | Stress | Non-edge sep |",
    "| --- | --- | ---: | ---: | ---: | ---: | ---: |"
  )

  for (i in seq_len(nrow(top_rows))) {
    row <- top_rows[i, ]
    lines <- c(lines, sprintf(
      "| %s | %s | %s | %s | %s | %s | %s |",
      row$candidate_id,
      row$candidate_source,
      format_num(row$score, 4),
      format_num(row$procrustes_rmse_mean, 4),
      format_num(row$edge_length_cv_mean, 4),
      format_num(row$sampled_stress_mean, 4),
      format_num(row$sampled_nonedge_sep_ratio_mean, 4)
    ))
  }

  lines <- c(
    lines,
    "",
    "Comparison on the shared comparison seed:",
    "",
    "| Layout | Procrustes RMSE | Edge CV | Sampled stress | Non-edge sep ratio |",
    "| --- | ---: | ---: | ---: | ---: |",
    sprintf("| legacy default | %s | %s | %s | %s |",
            format_num(standard_default_metrics$rmse, 4),
            format_num(standard_default_metrics$edge_cv, 4),
            format_num(standard_default_metrics$stress, 4),
            format_num(standard_default_metrics$sep, 4)),
    sprintf("| current tuned | %s | %s | %s | %s |",
            format_num(tuned_metrics$rmse, 4),
            format_num(tuned_metrics$edge_cv, 4),
            format_num(tuned_metrics$stress, 4),
            format_num(tuned_metrics$sep, 4)),
    "",
    "Outputs:",
    "",
    sprintf("- comparison PDF: `%s`", comparison_pdf),
    sprintf("- comparison PNG: `%s`", comparison_png)
  )

  writeLines(lines, con = path)
}

args <- parse_named_args(commandArgs(trailingOnly = TRUE))

level <- if (!is.null(args$level)) parse_int_scalar(args$level, "level") else 3L
n_random <- if (!is.null(args$n_random)) parse_int_scalar(args$n_random, "n_random") else 120L
search_seed <- if (!is.null(args$search_seed)) parse_int_scalar(args$search_seed, "search_seed") else 20260327L
comparison_seed <- if (!is.null(args$comparison_seed)) parse_int_scalar(args$comparison_seed, "comparison_seed") else 6L
seeds <- if (!is.null(args$seeds)) parse_int_vector(args$seeds, "seeds") else 1:10
top_n <- if (!is.null(args$top_n)) parse_int_scalar(args$top_n, "top_n") else 12L
run_tag <- validate_run_tag(if (!is.null(args$tag)) args$tag else sprintf("globalrep-carpet-level%d-tuning", level))

tuning_root <- file.path("dev", "manual")
tuning_pdf_dir <- file.path(tuning_root, "pdf", run_tag, "carpet")
tuning_tmp_dir <- file.path(tuning_root, "tmp", run_tag)
tuning_preview_dir <- file.path(tuning_tmp_dir, "pdf-previews")
dir.create(tuning_pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tuning_tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tuning_preview_dir, recursive = TRUE, showWarnings = FALSE)

built <- build_sierpinski_carpet(level)
edges <- edges.sierpinski.carpet(level)
if (!identical(unname(built$edges), unname(edges))) {
  stop("Canonical carpet builder does not match package edge generator")
}
spec <- list(family = "carpet", level = level, edges = edges, canonical = built$coords)

candidates <- generate_candidates(n_random = n_random, search_seed = search_seed)

pb_total <- length(candidates) * length(seeds)
pb <- utils::txtProgressBar(min = 0, max = pb_total, style = 3)
pb_idx <- 0L
raw_metrics <- tryCatch(
  do.call(
    rbind,
    lapply(candidates, function(candidate) {
      do.call(
        rbind,
        lapply(seeds, function(seed) {
          pb_idx <<- pb_idx + 1L
          utils::setTxtProgressBar(pb, pb_idx)
          run_candidate(spec, candidate, seed)
        })
      )
    })
  ),
  finally = close(pb)
)

summary_metrics <- summarize_candidates(raw_metrics)
family_rankings <- score_candidates(summary_metrics)
best_candidate <- family_rankings[1L, , drop = FALSE]

n <- max(edges)
standard_default_coords <- grip.layout.legacy(
  edges = edges,
  n = n,
  dim = 2,
  seed = comparison_seed
)
tuned_coords <- grip.layout(
  edges = edges,
  n = n,
  dim = 2,
  placement = best_candidate$placement[[1L]],
  rounds = best_candidate$rounds[[1L]],
  final_rounds = best_candidate$final_rounds[[1L]],
  num_init = best_candidate$num_init[[1L]],
  num_nbrs = best_candidate$num_nbrs[[1L]],
  r = best_candidate$r[[1L]],
  s = best_candidate$s[[1L]],
  repulsion_factor = best_candidate$repulsion_factor[[1L]],
  coarse_repulsion_factor = best_candidate$coarse_repulsion_factor[[1L]],
  coarse_repulsion_sample = best_candidate$coarse_repulsion_sample[[1L]],
  coarse_repulsion_exact_below = best_candidate$coarse_repulsion_exact_below[[1L]],
  seed = comparison_seed
)

standard_default_metrics <- {
  adj <- make_adj_list(edges, n)
  aligned <- align_to_target(standard_default_coords, spec$canonical)
  edge_stats <- edge_length_stats(standard_default_coords, edges)
  list(
    rmse = aligned$rmse,
    edge_cv = edge_stats$cv,
    stress = sampled_stress(standard_default_coords, adj, sample_size = 2000L,
                            rng_seed = 1000L + comparison_seed),
    sep = sampled_nonedge_separation_ratio(standard_default_coords, edges,
                                           sample_size = 5000L,
                                           rng_seed = 2000L + comparison_seed)
  )
}

tuned_metrics <- {
  adj <- make_adj_list(edges, n)
  aligned <- align_to_target(tuned_coords, spec$canonical)
  edge_stats <- edge_length_stats(tuned_coords, edges)
  list(
    rmse = aligned$rmse,
    edge_cv = edge_stats$cv,
    stress = sampled_stress(tuned_coords, adj, sample_size = 2000L,
                            rng_seed = 1000L + comparison_seed),
    sep = sampled_nonedge_separation_ratio(tuned_coords, edges,
                                           sample_size = 5000L,
                                           rng_seed = 2000L + comparison_seed)
  )
}

comparison_pdf <- file.path(
  tuning_pdf_dir,
  sprintf("sierpinski-carpet-level-%d-standard-default-vs-globalrep-best.pdf", level)
)
comparison_png <- file.path(
  tuning_preview_dir,
  sprintf("sierpinski-carpet-level-%d-standard-default-vs-globalrep-best.png", level)
)
summary_md <- file.path(tuning_tmp_dir, "sierpinski-globalrep-carpet-summary.md")
raw_csv <- file.path(tuning_tmp_dir, "sierpinski-globalrep-carpet-raw-metrics.csv")
ranking_csv <- file.path(tuning_tmp_dir, "sierpinski-globalrep-carpet-family-ranking.csv")

subtitle <- sprintf(
  "comparison seed=%d | tuned candidate=%s",
  comparison_seed,
  best_candidate$candidate_id[[1L]]
)
write_triptych_outputs(
  pdf_path = comparison_pdf,
  png_path = comparison_png,
  canonical_coords = spec$canonical,
  baseline_coords = standard_default_coords,
  tuned_coords = tuned_coords,
  edges = edges,
  title_text = sprintf("Sierpinski carpet level %d", level),
  subtitle_text = subtitle
)

utils::write.csv(raw_metrics, raw_csv, row.names = FALSE)
utils::write.csv(family_rankings, ranking_csv, row.names = FALSE)
write_summary_markdown(
  path = summary_md,
  config = list(
    level = level,
    n_random = n_random,
    seeds = seeds,
    comparison_seed = comparison_seed,
    search_seed = search_seed
  ),
  family_rankings = family_rankings,
  top_n = top_n,
  standard_default_metrics = standard_default_metrics,
  tuned_metrics = tuned_metrics,
  tuned_candidate = best_candidate,
  comparison_pdf = comparison_pdf,
  comparison_png = comparison_png
)

message(sprintf("Raw tuning metrics written to %s", raw_csv))
message(sprintf("Family ranking written to %s", ranking_csv))
message(sprintf("Markdown summary written to %s", summary_md))
message(sprintf("Comparison PDF written to %s", comparison_pdf))
message(sprintf("Comparison PNG written to %s", comparison_png))
