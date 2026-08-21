#!/usr/bin/env Rscript

sys.source(file.path("tools", "benchmarks", "gkk_lgkk_paper", "benchmark-sierpinski-baseline.R"), envir = environment())
sys.source(file.path("tools", "reports", "gkk_lgkk_paper", "generate-sierpinski-diagnostics.R"), envir = environment())
sys.source(file.path("tools", "benchmarks", "gkk_lgkk_paper", "tune-torus-parameters.R"), envir = environment())

score_weights_generic <- c(
  procrustes_rmse = 0.45,
  edge_length_cv = 0.20,
  sampled_stress = 0.20,
  sampled_nonedge_sep_ratio = 0.15
)

candidate_fields_generic <- c(
  "placement", "rounds", "final_rounds", "num_init",
  "num_nbrs", "r", "s", "repulsion_factor"
)

mesh_preset_profile <- function() {
  get("grip.mesh.preset.defaults", envir = asNamespace("grip"))()
}

torus_preset_profile <- function() {
  get("grip.torus.preset.defaults", envir = asNamespace("grip"))()
}

tree_preset_profile <- function() {
  get("grip.tree.preset.defaults", envir = asNamespace("grip"))()
}

align_to_canonical_nd <- function(source, target) {
  aligned <- align_to_target_nd(source, target)
  target_norm <- normalize_coords(target)
  list(
    aligned = aligned,
    rmse = sqrt(mean(rowSums((aligned - target_norm)^2)))
  )
}

build_mesh_spec <- function(h, w = h) {
  edges <- edges.mesh(h, w)
  coords <- matrix(0, nrow = h * w, ncol = 2L)
  idx <- 1L
  for (i in seq_len(h)) {
    for (j in seq_len(w)) {
      coords[idx, ] <- c(j - 1L, h - i)
      idx <- idx + 1L
    }
  }
  list(
    family = "mesh",
    graph_label = sprintf("%dx%d", h, w),
    title = sprintf("Mesh %dx%d", h, w),
    dim = 2L,
    edges = edges,
    canonical = coords,
    aligner = "nd"
  )
}

build_cylinder_spec <- function(h, w) {
  edges <- edges.cylinder(h, w)
  coords <- matrix(0, nrow = h * w, ncol = 3L)
  idx <- 1L
  for (i in seq_len(h)) {
    z <- seq(1, -1, length.out = h)[[i]]
    for (j in seq_len(w)) {
      phi <- 2 * pi * (j - 1L) / w
      coords[idx, ] <- c(cos(phi), sin(phi), z)
      idx <- idx + 1L
    }
  }
  arr <- coords_to_grid_array(coords, h, w)
  list(
    family = "cylinder",
    graph_label = sprintf("%dx%d", h, w),
    title = sprintf("Cylinder %dx%d", h, w),
    dim = 3L,
    edges = edges,
    canonical = coords,
    canonical_arr = arr,
    aligner = "cylinder"
  )
}

best_cylinder_alignment <- function(source_coords, spec) {
  variants <- list(
    spec$canonical_arr,
    spec$canonical_arr[seq(dim(spec$canonical_arr)[1L], 1L), , , drop = FALSE],
    spec$canonical_arr[, seq(dim(spec$canonical_arr)[2L], 1L), , drop = FALSE],
    spec$canonical_arr[seq(dim(spec$canonical_arr)[1L], 1L),
                       seq(dim(spec$canonical_arr)[2L], 1L), , drop = FALSE]
  )

  best <- NULL
  for (arr in variants) {
    for (shift in 0:(dim(arr)[2L] - 1L)) {
      target <- grid_array_to_coords(shift_grid_array(arr, shift, 2L))
      fit <- align_to_target_orthogonal(source_coords, target, allow_reflection = TRUE)
      if (is.null(best) || fit$rmse < best$rmse) {
        best <- fit
      }
    }
  }
  best
}

build_binary_tree_spec <- function(depth, k = 2L) {
  edges <- edges.kary.tree(k = k, depth = depth)
  n <- max(edges)
  adj <- make_adj_list(edges, n)
  parent <- integer(n)
  level <- integer(n)
  order <- integer(n)
  head <- 1L
  tail <- 1L
  order[[tail]] <- 1L
  parent[[1L]] <- 0L
  level[[1L]] <- 0L
  while (head <= tail) {
    v <- order[[head]]
    head <- head + 1L
    kids <- sort(adj[[v]][adj[[v]] != parent[[v]]])
    if (length(kids) == 0L) {
      next
    }
    for (u in kids) {
      if (u == parent[[v]]) next
      tail <- tail + 1L
      order[[tail]] <- u
      parent[[u]] <- v
      level[[u]] <- level[[v]] + 1L
    }
  }

  children <- vector("list", n)
  for (v in 2:n) {
    children[[parent[[v]]]] <- c(children[[parent[[v]]]], v)
  }
  children <- lapply(children, sort)

  x <- numeric(n)
  next_leaf <- 1L
  place <- function(v) {
    kids <- children[[v]]
    if (length(kids) == 0L) {
      x[[v]] <<- next_leaf
      next_leaf <<- next_leaf + 1L
      return(invisible(NULL))
    }
    for (u in kids) place(u)
    x[[v]] <<- mean(x[kids])
  }
  place(1L)

  coords <- cbind(
    x = x,
    y = -level
  )
  list(
    family = "kary.tree",
    graph_label = sprintf("k%d-depth%d", k, depth),
    title = sprintf("Binary tree depth %d", depth),
    dim = 2L,
    edges = edges,
    canonical = coords,
    aligner = "nd"
  )
}

sample_candidate_generic <- function(space) {
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
    repulsion_factor = as.double(sample(space$repulsion_factor, 1L))
  )
}

generic_candidate_key <- function(candidate) {
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

generate_candidates_generic <- function(space,
                                        n_random,
                                        search_seed,
                                        baseline_profile,
                                        baseline_id,
                                        reference_profiles = list()) {
  baseline_candidate <- c(
    list(candidate_id = baseline_id, candidate_source = "baseline"),
    baseline_profile[candidate_fields_generic]
  )
  candidates <- list(baseline_candidate)
  seen <- generic_candidate_key(baseline_candidate)

  if (length(reference_profiles) > 0L) {
    for (ref_id in names(reference_profiles)) {
      ref_candidate <- c(
        list(candidate_id = ref_id, candidate_source = "reference"),
        reference_profiles[[ref_id]][candidate_fields_generic]
      )
      candidates[[length(candidates) + 1L]] <- ref_candidate
      seen <- c(seen, generic_candidate_key(ref_candidate))
    }
  }

  set.seed(search_seed)
  attempts <- 0L
  max_attempts <- max(1000L, n_random * 250L)
  while (sum(vapply(candidates, function(x) identical(x$candidate_source, "random"), logical(1L))) < n_random &&
         attempts < max_attempts) {
    attempts <- attempts + 1L
    sampled <- sample_candidate_generic(space)
    key <- generic_candidate_key(sampled)
    if (key %in% seen) next
    idx <- sum(vapply(candidates, function(x) identical(x$candidate_source, "random"), logical(1L))) + 1L
    candidates[[length(candidates) + 1L]] <- c(
      list(candidate_id = sprintf("%s_rand_%03d", space$family_label, idx), candidate_source = "random"),
      sampled
    )
    seen <- c(seen, key)
  }
  candidates
}

run_one_layout_family <- function(spec, cfg, seed, stress_sample = 3000L, sep_sample = 6000L) {
  n <- max(spec$edges)
  adj <- make_adj_list(spec$edges, n)
  started <- proc.time()[["elapsed"]]
  coords <- legacy.grip(
    edges = spec$edges,
    n = n,
    dim = spec$dim,
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
  fit <- if (identical(spec$aligner, "cylinder")) {
    best_cylinder_alignment(coords, spec)
  } else {
    align_to_canonical_nd(coords, spec$canonical)
  }
  edge_stats <- edge_length_stats(coords, spec$edges)

  data.frame(
    family = spec$family,
    graph_label = spec$graph_label,
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
    procrustes_rmse = fit$rmse,
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

run_one_layout_family_safe <- function(spec, candidate, seed) {
  started <- proc.time()[["elapsed"]]
  metrics <- tryCatch(
    run_one_layout_family(spec, candidate, seed),
    error = function(e) {
      n <- max(spec$edges)
      data.frame(
        family = spec$family,
        graph_label = spec$graph_label,
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
    "family", "graph_label", "candidate_id", "candidate_source", "seed",
    "vertices", "edges", "placement", "rounds", "final_rounds", "num_init",
    "num_nbrs", "r", "s", "repulsion_factor", "status", "error_message",
    "procrustes_rmse", "edge_length_cv", "median_edge_length", "sampled_stress",
    "sampled_nonedge_sep_ratio", "elapsed_sec"
  )]
}

summarize_candidate_graphs_generic <- function(raw_metrics) {
  groups <- split(raw_metrics, paste(raw_metrics$graph_label, raw_metrics$candidate_id, sep = "|"))
  do.call(rbind, lapply(groups, function(df) {
    ok <- df$status == "ok"
    good <- df[ok, , drop = FALSE]
    best_seed <- if (nrow(good) > 0L) good$seed[[which.min(good$procrustes_rmse)]] else NA_integer_
    data.frame(
      family = df$family[[1L]],
      graph_label = df$graph_label[[1L]],
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
      edge_length_cv_mean = if (nrow(good) > 0L) mean(good$edge_length_cv) else NA_real_,
      sampled_stress_mean = if (nrow(good) > 0L) mean(good$sampled_stress) else NA_real_,
      sampled_nonedge_sep_ratio_mean = if (nrow(good) > 0L) mean(good$sampled_nonedge_sep_ratio) else NA_real_,
      elapsed_sec_mean = if (nrow(good) > 0L) mean(good$elapsed_sec) else NA_real_,
      best_seed = best_seed,
      stringsAsFactors = FALSE
    )
  }))
}

score_candidate_graphs_generic <- function(graph_summary) {
  split_groups <- split(graph_summary, graph_summary$graph_label)
  scored <- lapply(split_groups, function(df) {
    df$rank_rmse <- rank01(df$procrustes_rmse_mean, higher_better = FALSE)
    df$rank_edge_cv <- rank01(df$edge_length_cv_mean, higher_better = FALSE)
    df$rank_stress <- rank01(df$sampled_stress_mean, higher_better = FALSE)
    df$rank_sep <- rank01(df$sampled_nonedge_sep_ratio_mean, higher_better = TRUE)
    df$score_graph <-
      score_weights_generic[["procrustes_rmse"]] * df$rank_rmse +
      score_weights_generic[["edge_length_cv"]] * df$rank_edge_cv +
      score_weights_generic[["sampled_stress"]] * df$rank_stress +
      score_weights_generic[["sampled_nonedge_sep_ratio"]] * df$rank_sep
    df$score_graph[!is.finite(df$procrustes_rmse_mean)] <- Inf
    df[order(df$score_graph, df$procrustes_rmse_mean), , drop = FALSE]
  })
  do.call(rbind, scored)
}

summarize_rankings_generic <- function(scored_graphs) {
  groups <- split(scored_graphs, scored_graphs$candidate_id)
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
      n_graphs = nrow(df),
      score_mean = if (all(!is.finite(df$score_graph))) Inf else mean(df$score_graph[is.finite(df$score_graph)]),
      procrustes_rmse_mean = mean(df$procrustes_rmse_mean, na.rm = TRUE),
      edge_length_cv_mean = mean(df$edge_length_cv_mean, na.rm = TRUE),
      sampled_stress_mean = mean(df$sampled_stress_mean, na.rm = TRUE),
      sampled_nonedge_sep_ratio_mean = mean(df$sampled_nonedge_sep_ratio_mean, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  out[order(out$score_mean, out$procrustes_rmse_mean), , drop = FALSE]
}

grid_match_indices_local <- function(grid, values) {
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
  for (field in names(global_space)) {
    grid <- global_space[[field]]
    selected <- unique(top_rows[[field]])
    idx <- grid_match_indices_local(grid, selected)
    expanded <- sort(unique(unlist(lapply(idx, function(i) {
      seq.int(max(1L, i - radius), min(length(grid), i + radius))
    }))))
    out[[field]] <- grid[expanded]
  }
  out
}

write_stage_summary <- function(path, result, pdf_paths) {
  top_rows <- head(result$rankings, result$top_n)
  best_row <- result$best_row
  comparison <- result$comparison
  lines <- c(
    sprintf("# %s tuning", result$spec$title),
    "",
    sprintf("- stage tag: `%s`", result$tag),
    sprintf("- baseline candidate: `%s`", result$baseline_id),
    sprintf("- comparison candidate: `%s`", result$comparison_candidate_id),
    sprintf("- search seed: `%d`", result$search_seed),
    sprintf("- random candidates: `%d`", result$n_random),
    sprintf("- seeds: `%s`", paste(result$seeds, collapse = ", ")),
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
  lines <- c(
    lines,
    "",
    "Best candidate:",
    sprintf("- `%s`", best_row$candidate_id[[1L]]),
    sprintf("- %s", format_candidate(best_row[1L, ])),
    "",
    "Baseline vs comparison candidate:",
    sprintf("- baseline RMSE `%s` vs tuned `%s`", format_num(comparison$baseline_rmse[[1L]], 4), format_num(comparison$tuned_rmse[[1L]], 4)),
    sprintf("- baseline edge CV `%s` vs tuned `%s`", format_num(comparison$baseline_edge_cv[[1L]], 4), format_num(comparison$tuned_edge_cv[[1L]], 4)),
    sprintf("- baseline stress `%s` vs tuned `%s`", format_num(comparison$baseline_stress[[1L]], 4), format_num(comparison$tuned_stress[[1L]], 4)),
    sprintf("- baseline non-edge sep `%s` vs tuned `%s`", format_num(comparison$baseline_sep[[1L]], 4), format_num(comparison$tuned_sep[[1L]], 4)),
    "",
    "Diagnostic PDFs:",
    sprintf("- `%s`", pdf_paths)
  )
  writeLines(lines, con = path)
}

run_family_stage <- function(tag,
                             spec,
                             search_space,
                             baseline_profile,
                             baseline_id,
                             reference_profiles,
                             n_random,
                             top_n,
                             seeds,
                             search_seed,
                             out_root) {
  stage_tmp_dir <- file.path(out_root, "tmp", tag)
  stage_pdf_dir <- file.path(out_root, "reports", tag, spec$family)
  stage_preview_dir <- file.path(stage_tmp_dir, "pdf-previews")
  dir.create(stage_tmp_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(stage_pdf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(stage_preview_dir, recursive = TRUE, showWarnings = FALSE)

  space <- c(list(family_label = spec$family), search_space)
  candidates <- generate_candidates_generic(
    space = space,
    n_random = n_random,
    search_seed = search_seed,
    baseline_profile = baseline_profile,
    baseline_id = baseline_id,
    reference_profiles = reference_profiles
  )

  raw_metrics <- do.call(
    rbind,
    lapply(candidates, function(candidate) {
      do.call(rbind, lapply(seeds, function(seed) run_one_layout_family_safe(spec, candidate, seed)))
    })
  )
  graph_summary <- summarize_candidate_graphs_generic(raw_metrics)
  scored_graphs <- score_candidate_graphs_generic(graph_summary)
  rankings <- summarize_rankings_generic(scored_graphs)
  best_row <- rankings[1L, , drop = FALSE]
  comparison_candidate_id <- if (!identical(best_row$candidate_id[[1L]], baseline_id)) {
    best_row$candidate_id[[1L]]
  } else if (nrow(rankings) > 1L) {
    rankings$candidate_id[[2L]]
  } else {
    baseline_id
  }

  comparison <- data.frame(
    baseline_rmse = scored_graphs$procrustes_rmse_mean[scored_graphs$candidate_id == baseline_id],
    tuned_rmse = scored_graphs$procrustes_rmse_mean[scored_graphs$candidate_id == comparison_candidate_id],
    baseline_edge_cv = scored_graphs$edge_length_cv_mean[scored_graphs$candidate_id == baseline_id],
    tuned_edge_cv = scored_graphs$edge_length_cv_mean[scored_graphs$candidate_id == comparison_candidate_id],
    baseline_stress = scored_graphs$sampled_stress_mean[scored_graphs$candidate_id == baseline_id],
    tuned_stress = scored_graphs$sampled_stress_mean[scored_graphs$candidate_id == comparison_candidate_id],
    baseline_sep = scored_graphs$sampled_nonedge_sep_ratio_mean[scored_graphs$candidate_id == baseline_id],
    tuned_sep = scored_graphs$sampled_nonedge_sep_ratio_mean[scored_graphs$candidate_id == comparison_candidate_id]
  )

  best_seed_baseline <- graph_summary$best_seed[graph_summary$candidate_id == baseline_id][[1L]]
  best_seed_tuned <- graph_summary$best_seed[graph_summary$candidate_id == comparison_candidate_id][[1L]]
  baseline_profile_resolved <- as.list(scored_graphs[scored_graphs$candidate_id == baseline_id, candidate_fields_generic][1L, ])
  tuned_profile_resolved <- as.list(scored_graphs[scored_graphs$candidate_id == comparison_candidate_id, candidate_fields_generic][1L, ])

  baseline_coords <- legacy.grip(
    edges = spec$edges,
    n = max(spec$edges),
    dim = spec$dim,
    placement = baseline_profile_resolved$placement,
    rounds = baseline_profile_resolved$rounds,
    final_rounds = baseline_profile_resolved$final_rounds,
    num_init = baseline_profile_resolved$num_init,
    num_nbrs = baseline_profile_resolved$num_nbrs,
    r = baseline_profile_resolved$r,
    s = baseline_profile_resolved$s,
    repulsion_factor = baseline_profile_resolved$repulsion_factor,
    seed = best_seed_baseline
  )
  tuned_coords <- legacy.grip(
    edges = spec$edges,
    n = max(spec$edges),
    dim = spec$dim,
    placement = tuned_profile_resolved$placement,
    rounds = tuned_profile_resolved$rounds,
    final_rounds = tuned_profile_resolved$final_rounds,
    num_init = tuned_profile_resolved$num_init,
    num_nbrs = tuned_profile_resolved$num_nbrs,
    r = tuned_profile_resolved$r,
    s = tuned_profile_resolved$s,
    repulsion_factor = tuned_profile_resolved$repulsion_factor,
    seed = best_seed_tuned
  )

  pdf_path <- file.path(stage_pdf_dir, sprintf("%s-baseline-vs-tuned.pdf", spec$graph_label))
  subtitle <- sprintf(
    "baseline=%s (seed=%d) | tuned=%s (seed=%d)",
    baseline_id, best_seed_baseline, comparison_candidate_id, best_seed_tuned
  )
  if (spec$dim == 2L) {
    baseline_aligned <- align_to_canonical_nd(baseline_coords, spec$canonical)$aligned
    tuned_aligned <- align_to_canonical_nd(tuned_coords, spec$canonical)$aligned
    plot_2d_triptych(
      path = pdf_path,
      canonical_coords = spec$canonical,
      baseline_coords = baseline_aligned,
      tuned_coords = tuned_aligned,
      edges = spec$edges,
      title_text = spec$title,
      subtitle_text = subtitle,
      tuned_label = "tuned"
    )
  } else {
    baseline_aligned <- best_cylinder_alignment(baseline_coords, spec)$aligned
    tuned_aligned <- best_cylinder_alignment(tuned_coords, spec)$aligned
    write_3d_diagnostic_pdf(
      path = pdf_path,
      canonical_coords = spec$canonical,
      baseline_coords = baseline_aligned,
      tuned_coords = tuned_aligned,
      edges = spec$edges,
      title_text = spec$title,
      subtitle_text = subtitle,
      tuned_label = "tuned"
    )
  }
  render_pdf_previews(pdf_path, stage_preview_dir)

  raw_csv <- file.path(stage_tmp_dir, "tuning-raw-metrics.csv")
  graph_csv <- file.path(stage_tmp_dir, "tuning-graph-summary.csv")
  rank_csv <- file.path(stage_tmp_dir, "tuning-family-ranking.csv")
  comparison_csv <- file.path(stage_tmp_dir, "tuning-baseline-vs-best.csv")
  summary_md <- file.path(stage_tmp_dir, "tuning-summary.md")
  utils::write.csv(raw_metrics, raw_csv, row.names = FALSE)
  utils::write.csv(scored_graphs, graph_csv, row.names = FALSE)
  utils::write.csv(rankings, rank_csv, row.names = FALSE)
  utils::write.csv(comparison, comparison_csv, row.names = FALSE)

  result <- list(
    tag = tag,
    spec = spec,
    rankings = rankings,
    best_row = best_row,
    top_n = top_n,
    baseline_id = baseline_id,
    comparison_candidate_id = comparison_candidate_id,
    comparison = comparison,
    seeds = seeds,
    search_seed = search_seed,
    n_random = n_random,
    summary_path = summary_md,
    pdf_path = pdf_path,
    ranking_csv = rank_csv
  )
  write_stage_summary(summary_md, result, pdf_path)
  result
}

extract_profile_from_row_generic <- function(row) {
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

family_configs <- list(
  list(
    key = "mesh",
    broad_spec = build_mesh_spec(8, 8),
    validation_spec = build_mesh_spec(12, 12),
    baseline_profile = mesh_preset_profile(),
    baseline_id = "mesh_preset_baseline",
    reference_profiles = list(
      mesh_default_reference = default_profile(),
      mesh_carpet_reference = carpet_preset_profile()
    ),
    search_space = list(
      placement = c("barycenter", "circle"),
      rounds = c(32L, 64L, 96L, 128L, 160L),
      final_rounds = c(64L, 96L, 128L, 160L, 192L, 224L, 288L),
      num_init = c(7L, 12L, 20L, 28L, 36L),
      num_nbrs = c(6L, 8L, 10L, 12L, 16L, 20L, 24L),
      r = c(0.03, 0.05, 0.07, 0.10, 0.15),
      s = c(1.5, 3.0, 4.5, 6.0, 7.5),
      repulsion_factor = c(0.5, 0.75, 1.0, 1.5, 2.0)
    )
  ),
  list(
    key = "cylinder",
    broad_spec = build_cylinder_spec(8, 12),
    validation_spec = build_cylinder_spec(12, 16),
    baseline_profile = torus_preset_profile(),
    baseline_id = "cylinder_torus_baseline",
    reference_profiles = list(cylinder_default_reference = default_profile()),
    search_space = list(
      placement = c("barycenter"),
      rounds = c(96L, 128L, 160L, 192L, 224L),
      final_rounds = c(160L, 192L, 224L, 256L, 288L, 320L),
      num_init = c(12L, 20L, 28L, 36L),
      num_nbrs = c(16L, 20L, 24L, 28L, 32L),
      r = c(0.03, 0.05, 0.07, 0.10),
      s = c(4.5, 6.0, 7.5, 9.0),
      repulsion_factor = c(0.5, 0.75, 1.0, 1.5)
    )
  ),
  list(
    key = "kary-tree",
    broad_spec = build_binary_tree_spec(5),
    validation_spec = build_binary_tree_spec(6),
    baseline_profile = tree_preset_profile(),
    baseline_id = "tree_preset_baseline",
    reference_profiles = list(tree_default_reference = default_profile()),
    search_space = list(
      placement = c("barycenter", "circle"),
      rounds = c(32L, 64L, 96L, 128L),
      final_rounds = c(64L, 96L, 128L, 160L, 192L),
      num_init = c(7L, 12L, 20L, 28L),
      num_nbrs = c(6L, 8L, 10L, 12L, 16L),
      r = c(0.05, 0.07, 0.10, 0.15),
      s = c(3.0, 4.5, 6.0, 7.5, 9.0),
      repulsion_factor = c(0.0, 0.5, 1.0, 1.5)
    )
  )
)

write_master_summary <- function(path, family_results) {
  lines <- c(
    "# Next Standard Families Tuning",
    "",
    "Action plan executed:",
    "1. Broad search on one representative graph per family.",
    "2. Derive a local second-pass search space from the best broad candidates.",
    "3. Validate on a larger graph for the same family.",
    "4. Review the generated diagnostic PDFs visually.",
    ""
  )
  for (entry in family_results) {
    lines <- c(
      lines,
      sprintf("## %s", entry$config$key),
      sprintf("- broad summary: `%s`", entry$broad$summary_path),
      sprintf("- validation summary: `%s`", entry$validation$summary_path),
      sprintf("- validation PDF: `%s`", entry$validation$pdf_path),
      sprintf("- broad best: %s", format_candidate(entry$broad$best_row[1L, ])),
      sprintf("- validation best: %s", format_candidate(entry$validation$best_row[1L, ])),
      ""
    )
  }
  writeLines(lines, con = path)
}

args <- parse_named_args(commandArgs(trailingOnly = TRUE))
run_tag <- validate_run_tag(if (!is.null(args$tag)) args$tag else "next-standard-families")
broad_n_random <- if (!is.null(args$broad_n_random)) parse_int_scalar(args$broad_n_random, "broad_n_random") else 40L
validation_n_random <- if (!is.null(args$validation_n_random)) parse_int_scalar(args$validation_n_random, "validation_n_random") else 24L
broad_top_n <- if (!is.null(args$broad_top_n)) parse_int_scalar(args$broad_top_n, "broad_top_n") else 8L
validation_top_n <- if (!is.null(args$validation_top_n)) parse_int_scalar(args$validation_top_n, "validation_top_n") else 8L
top_k <- if (!is.null(args$top_k)) parse_int_scalar(args$top_k, "top_k") else 5L
radius <- if (!is.null(args$radius)) parse_int_scalar(args$radius, "radius") else 1L
broad_seeds <- if (!is.null(args$broad_seeds)) parse_int_vector(args$broad_seeds, "broad_seeds") else 1:4
validation_seeds <- if (!is.null(args$validation_seeds)) parse_int_vector(args$validation_seeds, "validation_seeds") else 1:4
base_search_seed <- if (!is.null(args$search_seed)) parse_int_scalar(args$search_seed, "search_seed") else 20260331L

out_root <- file.path("output", "gkk_lgkk_paper")
family_results <- list()
for (i in seq_along(family_configs)) {
  cfg <- family_configs[[i]]
  broad_tag <- validate_run_tag(sprintf("%s-%s-broad", run_tag, cfg$key))
  broad_result <- run_family_stage(
    tag = broad_tag,
    spec = cfg$broad_spec,
    search_space = cfg$search_space,
    baseline_profile = cfg$baseline_profile,
    baseline_id = cfg$baseline_id,
    reference_profiles = cfg$reference_profiles,
    n_random = broad_n_random,
    top_n = broad_top_n,
    seeds = broad_seeds,
    search_seed = base_search_seed + i,
    out_root = out_root
  )

  top_rows <- head(broad_result$rankings, top_k)
  local_space <- derive_local_space(top_rows, cfg$search_space, radius = radius)
  validation_tag <- validate_run_tag(sprintf("%s-%s-validation", run_tag, cfg$key))
  validation_result <- run_family_stage(
    tag = validation_tag,
    spec = cfg$validation_spec,
    search_space = local_space,
    baseline_profile = extract_profile_from_row_generic(broad_result$best_row),
    baseline_id = sprintf("%s-broad-best", cfg$key),
    reference_profiles = cfg$reference_profiles,
    n_random = validation_n_random,
    top_n = validation_top_n,
    seeds = validation_seeds,
    search_seed = base_search_seed + 100L + i,
    out_root = out_root
  )

  family_results[[i]] <- list(
    config = cfg,
    broad = broad_result,
    validation = validation_result
  )
}

master_tmp_dir <- file.path(out_root, "tmp", run_tag)
dir.create(master_tmp_dir, recursive = TRUE, showWarnings = FALSE)
master_summary_path <- file.path(master_tmp_dir, "next-standard-families-summary.md")
write_master_summary(master_summary_path, family_results)
cat(master_summary_path, "\n")
