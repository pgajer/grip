#!/usr/bin/env Rscript

sys.source(file.path("tools", "benchmarks", "gkk_lgkk_paper", "tune-torus-parameters.R"), envir = environment())

parse_named_args_local <- function(args) {
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

collapse_arg_values <- function(values) {
  paste(values, collapse = ",")
}

escape_tag_part <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "-", x)
}

build_space_args <- function(space) {
  c(
    sprintf("placement=%s", collapse_arg_values(space$placement)),
    sprintf("rounds=%s", collapse_arg_values(space$rounds)),
    sprintf("final_rounds=%s", collapse_arg_values(space$final_rounds)),
    sprintf("num_init=%s", collapse_arg_values(space$num_init)),
    sprintf("num_nbrs=%s", collapse_arg_values(space$num_nbrs)),
    sprintf("r=%s", collapse_arg_values(space$r)),
    sprintf("s=%s", collapse_arg_values(space$s)),
    sprintf("repulsion_factor=%s", collapse_arg_values(space$repulsion_factor))
  )
}

build_baseline_args <- function(profile, baseline_id) {
  c(
    sprintf("baseline_id=%s", baseline_id),
    sprintf("baseline_placement=%s", profile$placement),
    sprintf("baseline_rounds=%d", profile$rounds),
    sprintf("baseline_final_rounds=%d", profile$final_rounds),
    sprintf("baseline_num_init=%d", profile$num_init),
    sprintf("baseline_num_nbrs=%d", profile$num_nbrs),
    sprintf("baseline_r=%s", profile$r),
    sprintf("baseline_s=%s", profile$s),
    sprintf("baseline_repulsion_factor=%s", profile$repulsion_factor)
  )
}

extract_profile_from_row <- function(row) {
  list(
    placement = row$placement[[1L]],
    rounds = as.integer(row$rounds[[1L]]),
    final_rounds = as.integer(row$final_rounds[[1L]]),
    num_init = as.integer(row$num_init[[1L]]),
    num_nbrs = as.integer(row$num_nbrs[[1L]]),
    r = as.double(row$r[[1L]]),
    s = as.double(row$s[[1L]]),
    repulsion_factor = as.double(row$repulsion_factor[[1L]])
  )
}

profile_to_text <- function(profile) {
  sprintf(
    "placement=%s, rounds=%d, final_rounds=%d, num_init=%d, num_nbrs=%d, r=%.2f, s=%.1f, repulsion_factor=%.2f",
    profile$placement,
    profile$rounds,
    profile$final_rounds,
    profile$num_init,
    profile$num_nbrs,
    profile$r,
    profile$s,
    profile$repulsion_factor
  )
}

grid_match_indices <- function(grid, values) {
  idx <- integer()
  for (value in values) {
    if (is.numeric(grid)) {
      distances <- abs(as.double(grid) - as.double(value))
      matches <- which(distances == min(distances))
    } else {
      matches <- which(grid == value)
    }
    idx <- c(idx, matches)
  }
  sort(unique(idx))
}

derive_local_space <- function(top_rows, global_space, radius = 1L) {
  out <- list()
  out$placement <- global_space$placement
  for (field in setdiff(candidate_fields, "placement")) {
    grid <- global_space[[field]]
    selected <- unique(top_rows[[field]])
    idx <- grid_match_indices(grid, selected)
    if (length(idx) == 0L) {
      stop(sprintf("Could not map top-candidate values for field '%s' back onto the master grid", field))
    }
    expanded <- sort(unique(unlist(lapply(idx, function(i) {
      seq.int(max(1L, i - radius), min(length(grid), i + radius))
    }))))
    out[[field]] <- grid[expanded]
  }
  out
}

run_torus_tuner <- function(tag,
                            sizes,
                            n_random,
                            top_n,
                            seeds,
                            search_seed,
                            space,
                            baseline_profile = NULL,
                            baseline_id = NULL,
                            compare_candidate_id = NULL,
                            include_default_reference = TRUE,
                            include_carpet_reference = TRUE) {
  args <- c(
    sprintf("tag=%s", tag),
    sprintf("sizes=%s", collapse_arg_values(vapply(sizes, `[[`, "", "label"))),
    sprintf("n_random=%d", n_random),
    sprintf("top_n=%d", top_n),
    sprintf("seeds=%s", collapse_arg_values(seeds)),
    sprintf("search_seed=%d", search_seed),
    sprintf("include_default_reference=%s", if (include_default_reference) "true" else "false"),
    sprintf("include_carpet_reference=%s", if (include_carpet_reference) "true" else "false")
  )
  args <- c(args, build_space_args(space))
  if (!is.null(baseline_profile)) {
    args <- c(args, build_baseline_args(baseline_profile, baseline_id))
  }
  if (!is.null(compare_candidate_id)) {
    args <- c(args, sprintf("compare_candidate_id=%s", compare_candidate_id))
  }

  status <- system2("Rscript", c(file.path("tools", "benchmarks", "gkk_lgkk_paper", "tune-torus-parameters.R"), args))
  if (!identical(status, 0L)) {
    stop(sprintf("torus tuning run '%s' failed with exit code %s", tag, status))
  }

  tmp_dir <- file.path("dev", "manual", "tmp", tag)
  list(
    tag = tag,
    tmp_dir = tmp_dir,
    pdf_dir = file.path("dev", "manual", "pdf", tag, "torus"),
    summary_path = file.path(tmp_dir, "torus-tuning-summary.md"),
    family_csv = file.path(tmp_dir, "torus-tuning-family-ranking.csv"),
    graph_csv = file.path(tmp_dir, "torus-tuning-graph-summary.csv"),
    comparison_csv = file.path(tmp_dir, "torus-tuning-baseline-vs-best.csv")
  )
}

read_top_rows <- function(family_csv, top_k) {
  df <- utils::read.csv(family_csv, stringsAsFactors = FALSE)
  df[seq_len(min(top_k, nrow(df))), , drop = FALSE]
}

read_graph_summary <- function(graph_csv) {
  utils::read.csv(graph_csv, stringsAsFactors = FALSE)
}

write_master_summary <- function(path,
                                 config,
                                 stage1_result,
                                 stage1_best_row,
                                 local_space,
                                 stage2_8_result,
                                 dense_results) {
  stage1_best_profile <- extract_profile_from_row(stage1_best_row)
  lines <- c(
    "# Torus Second-Pass Tuning",
    "",
    "Workflow:",
    sprintf("- Stage 1 broad search on `%s`", config$stage1_size$label),
    sprintf("- Stage 2 local refinement on `%s`", config$stage1_size$label),
    sprintf("- Stage 2 local transfer searches on `%s`", paste(vapply(config$dense_sizes, `[[`, "", "label"), collapse = ", ")),
    "",
    "Stage 1 best candidate:",
    sprintf("- `%s`", stage1_best_row$candidate_id[[1L]]),
    sprintf("- %s", profile_to_text(stage1_best_profile)),
    sprintf("- stage 1 summary: `%s`", stage1_result$summary_path),
    "",
    "Derived second-pass local search space:",
    sprintf("- rounds: `%s`", collapse_arg_values(local_space$rounds)),
    sprintf("- final_rounds: `%s`", collapse_arg_values(local_space$final_rounds)),
    sprintf("- num_init: `%s`", collapse_arg_values(local_space$num_init)),
    sprintf("- num_nbrs: `%s`", collapse_arg_values(local_space$num_nbrs)),
    sprintf("- r: `%s`", collapse_arg_values(local_space$r)),
    sprintf("- s: `%s`", collapse_arg_values(local_space$s)),
    sprintf("- repulsion_factor: `%s`", collapse_arg_values(local_space$repulsion_factor)),
    "",
    "Stage 2 local refinement on stage-1 size:",
    sprintf("- summary: `%s`", stage2_8_result$summary_path),
    sprintf("- comparison table: `%s`", stage2_8_result$comparison_csv),
    ""
  )

  lines <- c(lines, "Stage 2 dense-size runs:", "")
  for (res in dense_results) {
    graph_df <- read_graph_summary(res$graph_csv)
    best_rows <- graph_df[graph_df$candidate_id != config$baseline_id, , drop = FALSE]
    best_by_graph <- do.call(rbind, lapply(split(graph_df, graph_df$graph_label), function(df) {
      df[order(df$score_graph, df$procrustes_rmse_mean), , drop = FALSE][1L, , drop = FALSE]
    }))
    best_by_graph <- best_by_graph[order(best_by_graph$h, best_by_graph$w), , drop = FALSE]
    lines <- c(lines, sprintf("- `%s`", res$tag))
    lines <- c(lines, sprintf("  summary: `%s`", res$summary_path))
    lines <- c(lines, sprintf("  comparison: `%s`", res$comparison_csv))
    if (nrow(best_by_graph) > 0L) {
      for (i in seq_len(nrow(best_by_graph))) {
        row <- best_by_graph[i, , drop = FALSE]
        lines <- c(lines, sprintf(
          "  best for %s: `%s` with %s",
          row$graph_label[[1L]],
          row$candidate_id[[1L]],
          format_candidate(row)
        ))
      }
    }
  }

  writeLines(lines, con = path)
}

args <- parse_named_args_local(commandArgs(trailingOnly = TRUE))
run_tag <- validate_run_tag(if (!is.null(args$tag)) args$tag else "torus-second-pass")
stage1_size <- if (!is.null(args$stage1_size)) parse_size_specs(args$stage1_size)[[1L]] else parse_size_specs("8x8")[[1L]]
dense_sizes <- if (!is.null(args$dense_sizes)) parse_size_specs(args$dense_sizes) else parse_size_specs("12x12,16x16,20x20")
stage1_n_random <- if (!is.null(args$stage1_n_random)) parse_int_scalar(args$stage1_n_random, "stage1_n_random") else 80L
stage2_n_random <- if (!is.null(args$stage2_n_random)) parse_int_scalar(args$stage2_n_random, "stage2_n_random") else 40L
stage1_top_n <- if (!is.null(args$stage1_top_n)) parse_int_scalar(args$stage1_top_n, "stage1_top_n") else 10L
stage2_top_n <- if (!is.null(args$stage2_top_n)) parse_int_scalar(args$stage2_top_n, "stage2_top_n") else 8L
top_k <- if (!is.null(args$top_k)) parse_int_scalar(args$top_k, "top_k") else 5L
radius <- if (!is.null(args$radius)) parse_int_scalar(args$radius, "radius") else 1L
stage1_search_seed <- if (!is.null(args$stage1_search_seed)) parse_int_scalar(args$stage1_search_seed, "stage1_search_seed") else 20260322L
stage2_search_seed <- if (!is.null(args$stage2_search_seed)) parse_int_scalar(args$stage2_search_seed, "stage2_search_seed") else 20260323L
stage1_seeds <- if (!is.null(args$stage1_seeds)) parse_int_vector(args$stage1_seeds, "stage1_seeds") else 1:6
stage2_seeds <- if (!is.null(args$stage2_seeds)) parse_int_vector(args$stage2_seeds, "stage2_seeds") else 1:4

if (top_k <= 0L) stop("top_k must be >= 1")
if (radius < 0L) stop("radius must be >= 0")

master_search_space <- search_space
stage1_tag <- validate_run_tag(sprintf("%s-stage1-%s", run_tag, escape_tag_part(stage1_size$label)))
stage1_result <- run_torus_tuner(
  tag = stage1_tag,
  sizes = list(stage1_size),
  n_random = stage1_n_random,
  top_n = stage1_top_n,
  seeds = stage1_seeds,
  search_seed = stage1_search_seed,
  space = master_search_space
)

stage1_top_rows <- read_top_rows(stage1_result$family_csv, top_k)
stage1_best_row <- stage1_top_rows[1L, , drop = FALSE]
stage1_best_profile <- extract_profile_from_row(stage1_best_row)
baseline_id <- sprintf("%s-baseline-%s", run_tag, escape_tag_part(stage1_size$label))
local_space <- derive_local_space(stage1_top_rows, master_search_space, radius = radius)

stage2_8_tag <- validate_run_tag(sprintf("%s-stage2-%s-local", run_tag, escape_tag_part(stage1_size$label)))
stage2_8_result <- run_torus_tuner(
  tag = stage2_8_tag,
  sizes = list(stage1_size),
  n_random = stage2_n_random,
  top_n = stage2_top_n,
  seeds = stage2_seeds,
  search_seed = stage2_search_seed,
  space = local_space,
  baseline_profile = stage1_best_profile,
  baseline_id = baseline_id
)

dense_results <- list()
for (i in seq_along(dense_sizes)) {
  size <- dense_sizes[[i]]
  dense_tag <- validate_run_tag(sprintf("%s-stage2-%s-local", run_tag, escape_tag_part(size$label)))
  dense_results[[i]] <- run_torus_tuner(
    tag = dense_tag,
    sizes = list(size),
    n_random = stage2_n_random,
    top_n = stage2_top_n,
    seeds = stage2_seeds,
    search_seed = stage2_search_seed + i,
    space = local_space,
    baseline_profile = stage1_best_profile,
    baseline_id = baseline_id
  )
}

combined_tag <- validate_run_tag(sprintf("%s-stage2-dense-combined", run_tag))
dense_results[[length(dense_results) + 1L]] <- run_torus_tuner(
  tag = combined_tag,
  sizes = dense_sizes,
  n_random = stage2_n_random,
  top_n = stage2_top_n,
  seeds = stage2_seeds,
  search_seed = stage2_search_seed + 100L,
  space = local_space,
  baseline_profile = stage1_best_profile,
  baseline_id = baseline_id
)

master_tmp_dir <- file.path("dev", "manual", "tmp", run_tag)
dir.create(master_tmp_dir, recursive = TRUE, showWarnings = FALSE)
master_summary_path <- file.path(master_tmp_dir, "torus-second-pass-summary.md")
write_master_summary(
  path = master_summary_path,
  config = list(
    stage1_size = stage1_size,
    dense_sizes = dense_sizes,
    baseline_id = baseline_id
  ),
  stage1_result = stage1_result,
  stage1_best_row = stage1_best_row,
  local_space = local_space,
  stage2_8_result = stage2_8_result,
  dense_results = dense_results
)

cat(master_summary_path, "\n")
