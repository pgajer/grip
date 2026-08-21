#!/usr/bin/env Rscript

sys.source(file.path("tools", "benchmarks", "gkk_lgkk_paper", "benchmark-sierpinski-baseline.R"), envir = environment())
sys.source(file.path("tools", "reports", "gkk_lgkk_paper", "generate-sierpinski-diagnostics.R"), envir = environment())

search_space <- list(
  placement = c("barycenter"),
  rounds = c(64L, 96L, 128L, 160L, 192L, 224L),
  final_rounds = c(128L, 160L, 192L, 224L, 256L, 288L, 320L, 384L),
  num_init = c(12L, 20L, 28L, 36L, 49L, 64L),
  num_nbrs = c(8L, 12L, 16L, 20L, 24L, 28L, 32L),
  r = c(0.01, 0.03, 0.05, 0.07, 0.10, 0.15),
  s = c(1.5, 3.0, 4.5, 6.0, 7.5, 9.0),
  repulsion_factor = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0)
)

score_weights <- c(
  procrustes_rmse = 0.45,
  edge_length_cv = 0.20,
  sampled_stress = 0.20,
  sampled_nonedge_sep_ratio = 0.15
)

candidate_fields <- c(
  "placement", "rounds", "final_rounds", "num_init",
  "num_nbrs", "r", "s", "repulsion_factor"
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

parse_char_vector <- function(x, name) {
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0L) {
    stop(sprintf("%s must be a comma-separated list", name))
  }
  parts
}

validate_run_tag <- function(x) {
  if (!grepl("^[A-Za-z0-9._-]+$", x)) {
    stop("tag must match ^[A-Za-z0-9._-]+$")
  }
  x
}

apply_search_overrides <- function(space, args) {
  for (field in names(space)) {
    if (is.null(args[[field]])) {
      next
    }
    value <- switch(
      field,
      placement = parse_char_vector(args[[field]], field),
      rounds = parse_int_vector(args[[field]], field),
      final_rounds = parse_int_vector(args[[field]], field),
      num_init = parse_int_vector(args[[field]], field),
      num_nbrs = parse_int_vector(args[[field]], field),
      r = parse_num_vector(args[[field]], field),
      s = parse_num_vector(args[[field]], field),
      repulsion_factor = parse_num_vector(args[[field]], field)
    )
    if (field == "placement") {
      bad <- setdiff(value, c("barycenter"))
      if (length(bad) > 0L) {
        stop("placement values for tetrahedron tuning must be 'barycenter'")
      }
    }
    if (field %in% c("rounds", "final_rounds", "num_init", "num_nbrs") && any(value <= 0L)) {
      stop(sprintf("%s values must be positive integers", field))
    }
    if (field == "r" && any(value < 0 | value > 1)) {
      stop("r values must be in [0, 1]")
    }
    if (field == "s" && any(value < 0)) {
      stop("s values must be >= 0")
    }
    if (field == "repulsion_factor" && any(value < 0)) {
      stop("repulsion_factor values must be >= 0")
    }
    space[[field]] <- value
  }
  space
}

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
  out
}

format_candidate <- function(row) {
  sprintf(
    "placement=%s, rounds=%d, final_rounds=%d, num_init=%d, num_nbrs=%d, r=%.2f, s=%.1f, repulsion_factor=%.2f",
    row$placement,
    row$rounds,
    row$final_rounds,
    row$num_init,
    row$num_nbrs,
    row$r,
    row$s,
    row$repulsion_factor
  )
}

sample_candidate <- function(space) {
  rounds <- sample(space$rounds, 1L)
  final_choices <- space$final_rounds[space$final_rounds >= rounds]
  if (length(final_choices) == 0L) {
    final_choices <- space$final_rounds
  }
  list(
    placement = "barycenter",
    rounds = as.integer(rounds),
    final_rounds = as.integer(sample(final_choices, 1L)),
    num_init = as.integer(sample(space$num_init, 1L)),
    num_nbrs = as.integer(sample(space$num_nbrs, 1L)),
    r = as.double(sample(space$r, 1L)),
    s = as.double(sample(space$s, 1L)),
    repulsion_factor = as.double(sample(space$repulsion_factor, 1L))
  )
}

generate_candidates <- function(space, n_random, search_seed, include_default_reference = TRUE) {
  carpet_candidate <- c(
    list(candidate_id = "tetrahedron_carpet_preset", candidate_source = "baseline"),
    carpet_preset_profile()[candidate_fields]
  )

  candidates <- list(carpet_candidate)
  seen <- candidate_key(carpet_candidate)

  if (isTRUE(include_default_reference)) {
    default_candidate <- c(
      list(candidate_id = "tetrahedron_default_reference", candidate_source = "reference"),
      default_profile()[candidate_fields]
    )
    candidates[[length(candidates) + 1L]] <- default_candidate
    seen <- c(seen, candidate_key(default_candidate))
  }

  set.seed(search_seed)
  attempts <- 0L
  max_attempts <- max(1000L, n_random * 250L)
  while (sum(vapply(candidates, function(x) identical(x$candidate_source, "random"), logical(1L))) < n_random &&
         attempts < max_attempts) {
    attempts <- attempts + 1L
    sampled <- sample_candidate(space)
    key <- candidate_key(sampled)
    if (key %in% seen) {
      next
    }
    idx <- sum(vapply(candidates, function(x) identical(x$candidate_source, "random"), logical(1L))) + 1L
    candidates[[length(candidates) + 1L]] <- c(
      list(candidate_id = sprintf("tetrahedron_rand_%03d", idx), candidate_source = "random"),
      sampled
    )
    seen <- c(seen, key)
  }

  candidates
}

build_graph_specs_for_levels <- function(levels) {
  lapply(levels, function(level) {
    built <- build_sierpinski_tetrahedron(level)
    package_edges <- edges.sierpinski.tetrahedron(level)
    if (!identical(unname(built$edges), unname(package_edges))) {
      stop(sprintf("Canonical tetrahedron builder does not match edges.sierpinski.tetrahedron(%d)", level))
    }
    list(
      family = "tetrahedron",
      level = level,
      edges = package_edges,
      canonical = built$coords
    )
  })
}

run_one_layout_tetra <- function(spec, cfg, seed, stress_sample = 3000L, sep_sample = 6000L) {
  n <- max(spec$edges)
  adj <- make_adj_list(spec$edges, n)
  started <- proc.time()[["elapsed"]]
  coords <- legacy.grip(
    edges = spec$edges,
    n = n,
    dim = 3,
    placement = cfg$placement,
    rounds = cfg$rounds,
    final_rounds = cfg$final_rounds,
    num_init = cfg$num_init,
    num_nbrs = cfg$num_nbrs,
    r = cfg$r,
    s = cfg$s,
    repulsion_factor = cfg$repulsion_factor,
    seed = seed
  )
  elapsed <- proc.time()[["elapsed"]] - started
  aligned <- align_to_target_nd(coords, spec$canonical)
  edge_stats <- edge_length_stats(coords, spec$edges)

  data.frame(
    family = spec$family,
    level = spec$level,
    seed = seed,
    vertices = n,
    edges = nrow(spec$edges),
    placement = cfg$placement,
    rounds = cfg$rounds,
    final_rounds = cfg$final_rounds,
    num_init = cfg$num_init,
    num_nbrs = cfg$num_nbrs,
    r = cfg$r,
    s = cfg$s,
    repulsion_factor = cfg$repulsion_factor,
    procrustes_rmse = sqrt(mean(rowSums((aligned - normalize_coords(spec$canonical))^2))),
    edge_length_cv = edge_stats$cv,
    median_edge_length = edge_stats$median,
    sampled_stress = sampled_stress(coords, adj, sample_size = stress_sample, rng_seed = 1000L + seed),
    sampled_nonedge_sep_ratio = sampled_nonedge_separation_ratio(
      coords, spec$edges, sample_size = sep_sample, rng_seed = 2000L + seed
    ),
    elapsed_sec = elapsed,
    stringsAsFactors = FALSE
  )
}

run_one_layout_safe <- function(spec, candidate, seed) {
  started <- proc.time()[["elapsed"]]
  metrics <- tryCatch(
    run_one_layout_tetra(spec, candidate, seed),
    error = function(e) {
      n <- max(spec$edges)
      data.frame(
        family = spec$family,
        level = spec$level,
        seed = seed,
        vertices = n,
        edges = nrow(spec$edges),
        placement = candidate$placement,
        rounds = candidate$rounds,
        final_rounds = candidate$final_rounds,
        num_init = candidate$num_init,
        num_nbrs = candidate$num_nbrs,
        r = candidate$r,
        s = candidate$s,
        repulsion_factor = candidate$repulsion_factor,
        procrustes_rmse = NA_real_,
        edge_length_cv = NA_real_,
        median_edge_length = NA_real_,
        sampled_stress = NA_real_,
        sampled_nonedge_sep_ratio = NA_real_,
        elapsed_sec = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  )
  if (!"elapsed_sec" %in% names(metrics) || !is.finite(metrics$elapsed_sec[[1L]])) {
    metrics$elapsed_sec <- proc.time()[["elapsed"]] - started
  }
  metrics$candidate_id <- candidate$candidate_id
  metrics$candidate_source <- candidate$candidate_source
  metrics$status <- if (is.finite(metrics$procrustes_rmse[[1L]])) "ok" else "error"
  metrics$error_message <- if (identical(metrics$status[[1L]], "ok")) "" else "layout failed"
  metrics[, c(
    "family", "level", "candidate_id", "candidate_source", "seed",
    "vertices", "edges", "placement", "rounds", "final_rounds", "num_init",
    "num_nbrs", "r", "s", "repulsion_factor", "status", "error_message",
    "procrustes_rmse", "edge_length_cv", "median_edge_length", "sampled_stress",
    "sampled_nonedge_sep_ratio", "elapsed_sec"
  )]
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

summarize_candidate_levels <- function(raw_metrics) {
  groups <- split(raw_metrics, paste(raw_metrics$family, raw_metrics$level, raw_metrics$candidate_id, sep = "|"))
  do.call(rbind, lapply(groups, function(df) {
    ok <- df$status == "ok"
    good <- df[ok, , drop = FALSE]
    best_seed <- if (nrow(good) > 0L) good$seed[[which.min(good$procrustes_rmse)]] else NA_integer_
    data.frame(
      family = df$family[[1L]],
      level = df$level[[1L]],
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
      vertices = df$vertices[[1L]],
      edges = df$edges[[1L]],
      n_runs = nrow(df),
      n_ok = sum(ok),
      n_fail = sum(!ok),
      procrustes_rmse_mean = if (nrow(good) > 0L) mean(good$procrustes_rmse) else NA_real_,
      procrustes_rmse_sd = if (nrow(good) > 1L) stats::sd(good$procrustes_rmse) else 0,
      edge_length_cv_mean = if (nrow(good) > 0L) mean(good$edge_length_cv) else NA_real_,
      edge_length_cv_sd = if (nrow(good) > 1L) stats::sd(good$edge_length_cv) else 0,
      sampled_stress_mean = if (nrow(good) > 0L) mean(good$sampled_stress) else NA_real_,
      sampled_stress_sd = if (nrow(good) > 1L) stats::sd(good$sampled_stress) else 0,
      sampled_nonedge_sep_ratio_mean = if (nrow(good) > 0L) mean(good$sampled_nonedge_sep_ratio) else NA_real_,
      sampled_nonedge_sep_ratio_sd = if (nrow(good) > 1L) stats::sd(good$sampled_nonedge_sep_ratio) else 0,
      elapsed_sec_mean = if (nrow(good) > 0L) mean(good$elapsed_sec) else NA_real_,
      best_seed = best_seed,
      stringsAsFactors = FALSE
    )
  }))
}

score_candidate_levels <- function(level_summary) {
  split_groups <- split(level_summary, level_summary$level)
  scored <- lapply(split_groups, function(df) {
    df$rank_rmse <- rank01(df$procrustes_rmse_mean, higher_better = FALSE)
    df$rank_edge_cv <- rank01(df$edge_length_cv_mean, higher_better = FALSE)
    df$rank_stress <- rank01(df$sampled_stress_mean, higher_better = FALSE)
    df$rank_sep <- rank01(df$sampled_nonedge_sep_ratio_mean, higher_better = TRUE)
    df$score_level <-
      score_weights[["procrustes_rmse"]] * df$rank_rmse +
      score_weights[["edge_length_cv"]] * df$rank_edge_cv +
      score_weights[["sampled_stress"]] * df$rank_stress +
      score_weights[["sampled_nonedge_sep_ratio"]] * df$rank_sep
    df$score_level[!is.finite(df$procrustes_rmse_mean)] <- Inf
    df[order(df$score_level, df$procrustes_rmse_mean), , drop = FALSE]
  })
  do.call(rbind, scored)
}

summarize_family_rankings <- function(scored_levels) {
  groups <- split(scored_levels, scored_levels$candidate_id)
  out <- do.call(rbind, lapply(groups, function(df) {
    data.frame(
      family = df$family[[1L]],
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
      n_levels = nrow(df),
      n_fail_levels = sum(!is.finite(df$score_level)),
      score_mean = if (all(!is.finite(df$score_level))) Inf else mean(df$score_level[is.finite(df$score_level)]),
      score_sd = if (sum(is.finite(df$score_level)) > 1L) stats::sd(df$score_level[is.finite(df$score_level)]) else 0,
      procrustes_rmse_mean = mean(df$procrustes_rmse_mean, na.rm = TRUE),
      edge_length_cv_mean = mean(df$edge_length_cv_mean, na.rm = TRUE),
      sampled_stress_mean = mean(df$sampled_stress_mean, na.rm = TRUE),
      sampled_nonedge_sep_ratio_mean = mean(df$sampled_nonedge_sep_ratio_mean, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  out[order(out$score_mean, out$procrustes_rmse_mean), , drop = FALSE]
}

format_num <- function(x, digits = 4) {
  ifelse(is.finite(x), format(round(x, digits), nsmall = digits, trim = TRUE), "NA")
}

build_level_comparison <- function(scored_levels, best_candidate_id, baseline_candidate_id) {
  levels <- sort(unique(scored_levels$level))
  do.call(rbind, lapply(levels, function(level) {
    level_df <- scored_levels[scored_levels$level == level, , drop = FALSE]
    baseline_row <- level_df[level_df$candidate_id == baseline_candidate_id, , drop = FALSE]
    tuned_row <- level_df[level_df$candidate_id == best_candidate_id, , drop = FALSE]
    data.frame(
      family = "tetrahedron",
      level = level,
      baseline_score = baseline_row$score_level[[1L]],
      tuned_score = tuned_row$score_level[[1L]],
      baseline_rmse = baseline_row$procrustes_rmse_mean[[1L]],
      tuned_rmse = tuned_row$procrustes_rmse_mean[[1L]],
      baseline_edge_cv = baseline_row$edge_length_cv_mean[[1L]],
      tuned_edge_cv = tuned_row$edge_length_cv_mean[[1L]],
      baseline_stress = baseline_row$sampled_stress_mean[[1L]],
      tuned_stress = tuned_row$sampled_stress_mean[[1L]],
      baseline_sep = baseline_row$sampled_nonedge_sep_ratio_mean[[1L]],
      tuned_sep = tuned_row$sampled_nonedge_sep_ratio_mean[[1L]],
      stringsAsFactors = FALSE
    )
  }))
}

write_tuning_summary <- function(path,
                                 config,
                                 family_rankings,
                                 scored_levels,
                                 level_comparison,
                                 pdf_paths) {
  top_rows <- head(family_rankings, config$top_n)
  best_overall <- top_rows[1L, , drop = FALSE]
  best_by_level <- do.call(rbind, lapply(split(scored_levels, scored_levels$level), function(df) {
    df[order(df$score_level, df$procrustes_rmse_mean), , drop = FALSE][1L, , drop = FALSE]
  }))
  best_by_level <- best_by_level[order(best_by_level$level), , drop = FALSE]

  lines <- c(
    "# Sierpinski Tetrahedron Parameter Tuning",
    "",
    "Search setup:",
    sprintf("- levels: `%s`", paste(config$levels, collapse = ", ")),
    sprintf("- random candidates: `%d`", config$n_random),
    sprintf("- seeds: `%s`", paste(config$seeds, collapse = ", ")),
    sprintf("- search seed: `%d`", config$search_seed),
    sprintf("- baseline candidate: `%s`", config$baseline_candidate_id),
    sprintf("- reference candidate: `%s`", config$reference_candidate_id),
    "",
    "Score definition (lower is better):",
    sprintf(
      "- `score_level = %.2f * rank(RMSE) + %.2f * rank(edge_length_cv) + %.2f * rank(sampled_stress) + %.2f * rank(-nonedge_sep)`",
      score_weights[["procrustes_rmse"]],
      score_weights[["edge_length_cv"]],
      score_weights[["sampled_stress"]],
      score_weights[["sampled_nonedge_sep_ratio"]]
    ),
    "",
    "Search space:",
    sprintf("- placement: `%s`", paste(config$search_space$placement, collapse = ", ")),
    sprintf("- rounds: `%s`", paste(config$search_space$rounds, collapse = ", ")),
    sprintf("- final_rounds: `%s`", paste(config$search_space$final_rounds, collapse = ", ")),
    sprintf("- num_init: `%s`", paste(config$search_space$num_init, collapse = ", ")),
    sprintf("- num_nbrs: `%s`", paste(config$search_space$num_nbrs, collapse = ", ")),
    sprintf("- r: `%s`", paste(config$search_space$r, collapse = ", ")),
    sprintf("- s: `%s`", paste(config$search_space$s, collapse = ", ")),
    sprintf("- repulsion_factor: `%s`", paste(config$search_space$repulsion_factor, collapse = ", ")),
    "",
    "Best overall candidate:",
    sprintf("- `%s`", best_overall$candidate_id[[1L]]),
    sprintf("- %s", format_candidate(best_overall[1L, ])),
    sprintf("- mean score across levels: %s", format_num(best_overall$score_mean[[1L]], digits = 4)),
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
      format_num(row$score_mean, 4),
      format_num(row$procrustes_rmse_mean, 4),
      format_num(row$edge_length_cv_mean, 4),
      format_num(row$sampled_stress_mean, 4),
      format_num(row$sampled_nonedge_sep_ratio_mean, 4)
    ))
  }

  lines <- c(lines, "", "Best candidate by level:", "", "| Level | Candidate | Score | Params |", "| ---: | --- | ---: | --- |")
  for (i in seq_len(nrow(best_by_level))) {
    row <- best_by_level[i, ]
    lines <- c(lines, sprintf(
      "| %d | %s | %s | %s |",
      row$level,
      row$candidate_id,
      format_num(row$score_level, 4),
      format_candidate(row)
    ))
  }

  lines <- c(lines, "", "Baseline vs best overall:", "", "| Level | Baseline score | Tuned score | Baseline RMSE | Tuned RMSE | Baseline edge CV | Tuned edge CV | Baseline stress | Tuned stress | Baseline sep | Tuned sep |", "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
  for (i in seq_len(nrow(level_comparison))) {
    row <- level_comparison[i, ]
    lines <- c(lines, sprintf(
      "| %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
      row$level,
      format_num(row$baseline_score, 4),
      format_num(row$tuned_score, 4),
      format_num(row$baseline_rmse, 4),
      format_num(row$tuned_rmse, 4),
      format_num(row$baseline_edge_cv, 4),
      format_num(row$tuned_edge_cv, 4),
      format_num(row$baseline_stress, 4),
      format_num(row$tuned_stress, 4),
      format_num(row$baseline_sep, 4),
      format_num(row$tuned_sep, 4)
    ))
  }

  lines <- c(lines, "", "Comparison PDFs:", "")
  lines <- c(lines, sprintf("- `%s`", pdf_paths))
  writeLines(lines, con = path)
}

args <- parse_named_args(commandArgs(trailingOnly = TRUE))
run_tag <- validate_run_tag(if (!is.null(args$tag)) args$tag else "sierpinski-tetrahedron-tuning")
levels <- if (!is.null(args$levels)) parse_int_vector(args$levels, "levels") else 3L
n_random <- if (!is.null(args$n_random)) parse_int_scalar(args$n_random, "n_random") else 120L
top_n <- if (!is.null(args$top_n)) parse_int_scalar(args$top_n, "top_n") else 8L
search_seed <- if (!is.null(args$search_seed)) parse_int_scalar(args$search_seed, "search_seed") else 20260322L
seeds <- if (!is.null(args$seeds)) parse_int_vector(args$seeds, "seeds") else 1:10
include_default_reference <- if (!is.null(args$include_default_reference)) {
  tolower(args$include_default_reference) %in% c("1", "true", "yes")
} else {
  TRUE
}

if (n_random < 0L) stop("n_random must be >= 0")
if (top_n <= 0L) stop("top_n must be >= 1")

search_space <- apply_search_overrides(search_space, args)

tuning_root <- file.path("output", "gkk_lgkk_paper")
tuning_pdf_dir <- file.path(tuning_root, "reports", run_tag, "tetrahedron")
tuning_tmp_dir <- file.path(tuning_root, "tmp", run_tag)
tuning_preview_dir <- file.path(tuning_tmp_dir, "pdf-previews")
dir.create(tuning_pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tuning_tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tuning_preview_dir, recursive = TRUE, showWarnings = FALSE)

candidates <- generate_candidates(search_space, n_random, search_seed, include_default_reference)
graphs <- build_graph_specs_for_levels(levels)

raw_metrics <- do.call(
  rbind,
  lapply(graphs, function(spec) {
    do.call(
      rbind,
      lapply(candidates, function(candidate) {
        do.call(rbind, lapply(seeds, function(seed) run_one_layout_safe(spec, candidate, seed)))
      })
    )
  })
)

level_summary <- summarize_candidate_levels(raw_metrics)
scored_levels <- score_candidate_levels(level_summary)
family_rankings <- summarize_family_rankings(scored_levels)

baseline_candidate_id <- "tetrahedron_carpet_preset"
reference_candidate_id <- if (include_default_reference) "tetrahedron_default_reference" else ""
best_candidate_id <- family_rankings$candidate_id[[1L]]
level_comparison <- build_level_comparison(scored_levels, best_candidate_id, baseline_candidate_id)

comparison_pdfs <- character()
for (spec in graphs) {
  level_df <- scored_levels[scored_levels$level == spec$level, , drop = FALSE]
  baseline_row <- level_df[level_df$candidate_id == baseline_candidate_id, , drop = FALSE]
  tuned_row <- level_df[level_df$candidate_id == best_candidate_id, , drop = FALSE]

  baseline_coords <- legacy.grip(
    edges = spec$edges,
    n = max(spec$edges),
    dim = 3,
    placement = baseline_row$placement[[1L]],
    rounds = baseline_row$rounds[[1L]],
    final_rounds = baseline_row$final_rounds[[1L]],
    num_init = baseline_row$num_init[[1L]],
    num_nbrs = baseline_row$num_nbrs[[1L]],
    r = baseline_row$r[[1L]],
    s = baseline_row$s[[1L]],
    repulsion_factor = baseline_row$repulsion_factor[[1L]],
    seed = baseline_row$best_seed[[1L]]
  )
  tuned_coords <- legacy.grip(
    edges = spec$edges,
    n = max(spec$edges),
    dim = 3,
    placement = tuned_row$placement[[1L]],
    rounds = tuned_row$rounds[[1L]],
    final_rounds = tuned_row$final_rounds[[1L]],
    num_init = tuned_row$num_init[[1L]],
    num_nbrs = tuned_row$num_nbrs[[1L]],
    r = tuned_row$r[[1L]],
    s = tuned_row$s[[1L]],
    repulsion_factor = tuned_row$repulsion_factor[[1L]],
    seed = tuned_row$best_seed[[1L]]
  )

  baseline_aligned <- align_to_target_nd(baseline_coords, spec$canonical)
  tuned_aligned <- align_to_target_nd(tuned_coords, spec$canonical)
  pdf_path <- file.path(
    tuning_pdf_dir,
    sprintf("sierpinski-tetrahedron-level-%d-baseline-vs-best-overall-3d.pdf", spec$level)
  )
  subtitle <- sprintf(
    "baseline=%s (seed=%d) | tuned=%s (seed=%d)",
    baseline_row$candidate_id[[1L]], baseline_row$best_seed[[1L]],
    tuned_row$candidate_id[[1L]], tuned_row$best_seed[[1L]]
  )
  write_3d_diagnostic_pdf(
    path = pdf_path,
    canonical_coords = spec$canonical,
    baseline_coords = baseline_aligned,
    tuned_coords = tuned_aligned,
    edges = spec$edges,
    title_text = sprintf("Sierpinski tetrahedron level %d", spec$level),
    subtitle_text = subtitle,
    tuned_label = "tuned"
  )
  comparison_pdfs <- c(comparison_pdfs, pdf_path)
}

raw_csv_path <- file.path(tuning_tmp_dir, "sierpinski-tetrahedron-tuning-raw-metrics.csv")
level_csv_path <- file.path(tuning_tmp_dir, "sierpinski-tetrahedron-tuning-level-summary.csv")
family_csv_path <- file.path(tuning_tmp_dir, "sierpinski-tetrahedron-tuning-family-ranking.csv")
comparison_csv_path <- file.path(tuning_tmp_dir, "sierpinski-tetrahedron-tuning-baseline-vs-best.csv")
summary_md_path <- file.path(tuning_tmp_dir, "sierpinski-tetrahedron-tuning-summary.md")

utils::write.csv(raw_metrics, raw_csv_path, row.names = FALSE)
utils::write.csv(scored_levels, level_csv_path, row.names = FALSE)
utils::write.csv(family_rankings, family_csv_path, row.names = FALSE)
utils::write.csv(level_comparison, comparison_csv_path, row.names = FALSE)
write_tuning_summary(
  path = summary_md_path,
  config = list(
    levels = levels,
    n_random = n_random,
    top_n = top_n,
    seeds = seeds,
    search_seed = search_seed,
    baseline_candidate_id = baseline_candidate_id,
    reference_candidate_id = reference_candidate_id,
    search_space = search_space
  ),
  family_rankings = family_rankings,
  scored_levels = scored_levels,
  level_comparison = level_comparison,
  pdf_paths = comparison_pdfs
)
render_pdf_previews(comparison_pdfs, tuning_preview_dir)

message(sprintf("Raw tuning metrics written to %s", raw_csv_path))
message(sprintf("Level summary written to %s", level_csv_path))
message(sprintf("Family ranking written to %s", family_csv_path))
message(sprintf("Baseline-vs-best comparison written to %s", comparison_csv_path))
message(sprintf("Markdown summary written to %s", summary_md_path))
message(sprintf("Comparison PDFs written to %s", tuning_pdf_dir))
