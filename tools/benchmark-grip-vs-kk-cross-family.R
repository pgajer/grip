#!/usr/bin/env Rscript

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Package 'pkgload' is required to run the cross-family GRIP vs KK benchmark.")
}

pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)

sys.source(file.path("tools", "benchmark-globalrep-fixed-candidate.R"), envir = environment())

if (!requireNamespace("callr", quietly = TRUE)) {
  stop("Package 'callr' is required to run the cross-family GRIP vs KK benchmark.")
}

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
grip_namespace <- asNamespace("grip")
large_graph_vertex_threshold <- 10000L

package_internal <- function(name) {
  get(name, envir = grip_namespace)
}

runtime_weight <- 0.20

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

validate_run_tag <- function(x) {
  if (!grepl("^[A-Za-z0-9._-]+$", x)) {
    stop("tag must match ^[A-Za-z0-9._-]+$")
  }
  x
}

current_generic_best_cfg <- list(
  placement = "barycenter",
  rounds = 160L,
  final_rounds = 32L,
  num_init = 24L,
  num_nbrs = 20L,
  r = 0.03,
  s = 7.5,
  repulsion_factor = 2.5,
  coarse_repulsion_factor = 1.5,
  coarse_repulsion_sample = 16L,
  coarse_repulsion_exact_below = 64L
)

fixed_final_cfg <- function(final_rounds, exact = FALSE, n = NULL) {
  cfg <- package_internal("grip.globalrep.base.defaults")(if (is.null(n)) 1000L else n)
  cfg$final_rounds <- as.integer(final_rounds)
  if (isTRUE(exact)) {
    cfg$coarse_repulsion_exact_below <- if (is.null(n)) 4096L else as.integer(n)
  }
  cfg
}

mesh_preset_reference_cfg <- function(n = NULL) {
  cfg <- package_internal("grip.globalrep.base.defaults")(if (is.null(n)) 1000L else n)
  preset <- package_internal("grip.mesh.preset.defaults")()
  cfg[names(preset)] <- preset
  cfg
}

torus_preset_reference_cfg <- function(n = NULL) {
  cfg <- package_internal("grip.globalrep.base.defaults")(if (is.null(n)) 1000L else n)
  preset <- package_internal("grip.torus.preset.defaults")()
  cfg[names(preset)] <- preset
  cfg
}

tree_preset_reference_cfg <- function(n = NULL, dim = 2L) {
  cfg <- package_internal("grip.globalrep.base.defaults")(if (is.null(n)) 1000L else n)
  preset <- package_internal("grip.tree.preset.defaults")(dim = dim)
  cfg[names(preset)] <- preset
  cfg
}

all_universal_candidate_ids <- c(
  "grip_default_adaptive",
  "grip_f32",
  "grip_f64",
  "grip_f96",
  "grip_f128",
  "grip_f192",
  "grip_f384",
  "grip_exact_f64",
  "grip_exact_f96"
)

candidate_label <- function(candidate_id) {
  switch(
    candidate_id,
    grip_default_adaptive = "GRIP adaptive default",
    grip_f32 = "GRIP f32",
    grip_f64 = "GRIP f64",
    grip_f96 = "GRIP f96",
    grip_f128 = "GRIP f128",
    grip_f192 = "GRIP f192",
    grip_f384 = "GRIP f384",
    grip_exact_f64 = "GRIP exact f64",
    grip_exact_f96 = "GRIP exact f96",
    grip_mesh_preset_reference = "GRIP mesh preset",
    grip_torus_preset_reference = "GRIP torus preset",
    grip_tree_preset_reference = "GRIP tree preset",
    igraph_kk_default = "igraph KK",
    candidate_id
  )
}

is_grip_candidate <- function(candidate_id) {
  startsWith(candidate_id, "grip_")
}

candidate_applicable <- function(candidate_id, spec) {
  if (candidate_id == "grip_mesh_preset_reference") {
    return(identical(spec$family, "mesh"))
  }
  if (candidate_id == "grip_torus_preset_reference") {
    return(identical(spec$family, "torus"))
  }
  if (candidate_id == "grip_tree_preset_reference") {
    return(identical(spec$family, "kary.tree"))
  }
  TRUE
}

resolve_candidate_cfg <- function(candidate_id, spec) {
  n <- max(spec$edges)
  switch(
    candidate_id,
    grip_default_adaptive = package_internal("grip.globalrep.base.defaults")(n),
    grip_f32 = fixed_final_cfg(32L, exact = FALSE, n = n),
    grip_f64 = fixed_final_cfg(64L, exact = FALSE, n = n),
    grip_f96 = fixed_final_cfg(96L, exact = FALSE, n = n),
    grip_f128 = fixed_final_cfg(128L, exact = FALSE, n = n),
    grip_f192 = fixed_final_cfg(192L, exact = FALSE, n = n),
    grip_f384 = fixed_final_cfg(384L, exact = FALSE, n = n),
    grip_exact_f64 = fixed_final_cfg(64L, exact = TRUE, n = n),
    grip_exact_f96 = fixed_final_cfg(96L, exact = TRUE, n = n),
    grip_mesh_preset_reference = mesh_preset_reference_cfg(n = n),
    grip_torus_preset_reference = torus_preset_reference_cfg(n = n),
    grip_tree_preset_reference = tree_preset_reference_cfg(n = n, dim = spec$dim),
    stop(sprintf("Unknown GRIP candidate_id: %s", candidate_id))
  )
}

seed_budget_grip <- function(n_vertices) {
  if (n_vertices <= 500L) {
    return(1:3)
  }
  if (n_vertices <= 5000L) {
    return(1:2)
  }
  1L
}

seed_budget_kk <- function(n_vertices) {
  1L
}

grip_time_limit_budget <- function(n_vertices) {
  if (n_vertices <= 1000L) {
    return(60L)
  }
  if (n_vertices <= 10000L) {
    return(180L)
  }
  180L
}

kk_time_limit_budget <- function(n_vertices) {
  if (n_vertices <= 1000L) {
    return(60L)
  }
  if (n_vertices <= 5000L) {
    return(180L)
  }
  240L
}

is_ok_status <- function(x) {
  identical(x, "ok")
}

grip_failure_row <- function(spec, candidate_id, seed, cfg, status, error_message) {
  n <- max(spec$edges)
  data.frame(
    family = spec$family,
    graph_label = spec$graph_label,
    title = spec$title,
    dim = spec$dim,
    candidate_id = candidate_id,
    candidate_label = candidate_label(candidate_id),
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
    coarse_repulsion_factor = cfg$coarse_repulsion_factor,
    coarse_repulsion_sample = cfg$coarse_repulsion_sample,
    coarse_repulsion_exact_below = cfg$coarse_repulsion_exact_below,
    status = status,
    error_message = error_message,
    procrustes_rmse = NA_real_,
    edge_length_cv = NA_real_,
    median_edge_length = NA_real_,
    sampled_stress = NA_real_,
    sampled_nonedge_sep_ratio = NA_real_,
    elapsed_sec = NA_real_,
    stringsAsFactors = FALSE
  )
}

run_one_layout_grip_external_safe <- function(spec, candidate_id, seed) {
  cfg <- resolve_candidate_cfg(candidate_id, spec)
  n <- max(spec$edges)
  adj <- make_adj_list(spec$edges, n)
  samples <- sample_budget(n)
  time_limit_sec <- grip_time_limit_budget(n)

  result <- tryCatch({
    callr::r(
      function(edges, dim, cfg, seed) {
        pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)
        started <- proc.time()[["elapsed"]]
        coords <- grip.layout(
          edges = edges,
          n = max(edges),
          dim = dim,
          placement = cfg$placement,
          rounds = cfg$rounds,
          final_rounds = cfg$final_rounds,
          num_init = cfg$num_init,
          num_nbrs = cfg$num_nbrs,
          r = cfg$r,
          s = cfg$s,
          repulsion_factor = cfg$repulsion_factor,
          coarse_repulsion_factor = cfg$coarse_repulsion_factor,
          coarse_repulsion_sample = cfg$coarse_repulsion_sample,
          coarse_repulsion_exact_below = cfg$coarse_repulsion_exact_below,
          seed = seed
        )
        elapsed <- proc.time()[["elapsed"]] - started
        list(coords = coords, elapsed_sec = elapsed)
      },
      args = list(edges = spec$edges, dim = spec$dim, cfg = cfg, seed = seed),
      timeout = time_limit_sec,
      stdout = "|",
      stderr = "|"
    )
  }, error = function(e) e)

  if (inherits(result, "error")) {
    status <- if (inherits(result, "callr_timeout_error") ||
      grepl("timeout|timed out|elapsed time limit", conditionMessage(result), ignore.case = TRUE)) {
      "timeout"
    } else {
      "error"
    }
    return(grip_failure_row(spec, candidate_id, seed, cfg, status, conditionMessage(result)))
  }

  coords <- result$coords
  elapsed <- result$elapsed_sec
  fit <- align_layout_to_spec(coords, spec)
  edge_stats <- edge_length_stats(coords, spec$edges)

  data.frame(
    family = spec$family,
    graph_label = spec$graph_label,
    title = spec$title,
    dim = spec$dim,
    candidate_id = candidate_id,
    candidate_label = candidate_label(candidate_id),
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
    coarse_repulsion_factor = cfg$coarse_repulsion_factor,
    coarse_repulsion_sample = cfg$coarse_repulsion_sample,
    coarse_repulsion_exact_below = cfg$coarse_repulsion_exact_below,
    status = "ok",
    error_message = "",
    procrustes_rmse = fit$rmse,
    edge_length_cv = edge_stats$cv,
    median_edge_length = edge_stats$median,
    sampled_stress = sampled_stress(coords, adj, sample_size = samples$stress, rng_seed = 1000L + seed),
    sampled_nonedge_sep_ratio = sampled_nonedge_separation_ratio(
      coords,
      spec$edges,
      sample_size = samples$sep,
      rng_seed = 2000L + seed
    ),
    elapsed_sec = elapsed,
    stringsAsFactors = FALSE
  )
}

run_one_layout_grip_safe <- function(spec, candidate_id, seed) {
  n <- max(spec$edges)
  if (n > large_graph_vertex_threshold) {
    return(run_one_layout_grip_external_safe(spec, candidate_id, seed))
  }
  cfg <- resolve_candidate_cfg(candidate_id, spec)
  run_one_layout_globalrep_safe(
    spec = spec,
    cfg = cfg,
    candidate_id = candidate_id,
    seed = seed,
    base_time_limit_sec = grip_time_limit_budget(n)
  )
}

run_one_layout_kk_safe <- function(spec, seed, kk_max_vertices) {
  n <- max(spec$edges)
  if (n > kk_max_vertices) {
    return(data.frame(
      family = spec$family,
      graph_label = spec$graph_label,
      title = spec$title,
      dim = spec$dim,
      candidate_id = "igraph_kk_default",
      candidate_label = candidate_label("igraph_kk_default"),
      seed = seed,
      vertices = n,
      edges = nrow(spec$edges),
      placement = NA_character_,
      rounds = NA_integer_,
      final_rounds = NA_integer_,
      num_init = NA_integer_,
      num_nbrs = NA_integer_,
      r = NA_real_,
      s = NA_real_,
      repulsion_factor = NA_real_,
      coarse_repulsion_factor = NA_real_,
      coarse_repulsion_sample = NA_integer_,
      coarse_repulsion_exact_below = NA_integer_,
      status = "skipped",
      error_message = sprintf("Skipped: %d vertices exceeds kk_max_vertices=%d", n, kk_max_vertices),
      procrustes_rmse = NA_real_,
      edge_length_cv = NA_real_,
      median_edge_length = NA_real_,
      sampled_stress = NA_real_,
      sampled_nonedge_sep_ratio = NA_real_,
      elapsed_sec = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  adj <- make_adj_list(spec$edges, n)
  samples <- sample_budget(n)
  time_limit_sec <- kk_time_limit_budget(n)

  result <- tryCatch(
    {
      callr::r(
        function(edges, dim, seed) {
          library(igraph)
          graph <- igraph::graph_from_edgelist(as.matrix(edges), directed = FALSE)
          set.seed(seed)
          started <- proc.time()[["elapsed"]]
          coords <- igraph::layout_with_kk(graph, dim = dim)
          elapsed <- proc.time()[["elapsed"]] - started
          list(coords = coords, elapsed_sec = elapsed)
        },
        args = list(edges = spec$edges, dim = spec$dim, seed = seed),
        timeout = time_limit_sec,
        spinner = FALSE
      )
    },
    error = function(e) e
  )

  if (inherits(result, "error")) {
    msg <- conditionMessage(result)
    status <- if (inherits(result, "callr_timeout_error") ||
                  grepl("timed out|timeout", msg, ignore.case = TRUE)) "timeout" else "error"
    return(data.frame(
      family = spec$family,
      graph_label = spec$graph_label,
      title = spec$title,
      dim = spec$dim,
      candidate_id = "igraph_kk_default",
      candidate_label = candidate_label("igraph_kk_default"),
      seed = seed,
      vertices = n,
      edges = nrow(spec$edges),
      placement = NA_character_,
      rounds = NA_integer_,
      final_rounds = NA_integer_,
      num_init = NA_integer_,
      num_nbrs = NA_integer_,
      r = NA_real_,
      s = NA_real_,
      repulsion_factor = NA_real_,
      coarse_repulsion_factor = NA_real_,
      coarse_repulsion_sample = NA_integer_,
      coarse_repulsion_exact_below = NA_integer_,
      status = status,
      error_message = msg,
      procrustes_rmse = NA_real_,
      edge_length_cv = NA_real_,
      median_edge_length = NA_real_,
      sampled_stress = NA_real_,
      sampled_nonedge_sep_ratio = NA_real_,
      elapsed_sec = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  coords <- result$coords
  fit <- align_layout_to_spec(coords, spec)
  edge_stats <- edge_length_stats(coords, spec$edges)

  data.frame(
    family = spec$family,
    graph_label = spec$graph_label,
    title = spec$title,
    dim = spec$dim,
    candidate_id = "igraph_kk_default",
    candidate_label = candidate_label("igraph_kk_default"),
    seed = seed,
    vertices = n,
    edges = nrow(spec$edges),
    placement = NA_character_,
    rounds = NA_integer_,
    final_rounds = NA_integer_,
    num_init = NA_integer_,
    num_nbrs = NA_integer_,
    r = NA_real_,
    s = NA_real_,
    repulsion_factor = NA_real_,
    coarse_repulsion_factor = NA_real_,
    coarse_repulsion_sample = NA_integer_,
    coarse_repulsion_exact_below = NA_integer_,
    status = "ok",
    error_message = "",
    procrustes_rmse = fit$rmse,
    edge_length_cv = edge_stats$cv,
    median_edge_length = edge_stats$median,
    sampled_stress = sampled_stress(coords, adj, sample_size = samples$stress, rng_seed = 5000L + seed),
    sampled_nonedge_sep_ratio = sampled_nonedge_separation_ratio(
      coords, spec$edges, sample_size = samples$sep, rng_seed = 6000L + seed
    ),
    elapsed_sec = result$elapsed_sec,
    stringsAsFactors = FALSE
  )
}

stage1_representative_specs <- function(specs, representative_max_vertices = 10000L) {
  non_sphere <- Filter(function(spec) !identical(spec$family, "sphere"), specs)
  families <- unique(vapply(non_sphere, `[[`, character(1L), "family"))
  reps <- lapply(families, function(family) {
    candidates <- Filter(function(spec) identical(spec$family, family) &&
                           max(spec$edges) <= representative_max_vertices, non_sphere)
    if (length(candidates) == 0L) {
      candidates <- Filter(function(spec) identical(spec$family, family), non_sphere)
    }
    candidates[[which.max(vapply(candidates, function(spec) max(spec$edges), numeric(1L)))]]
  })
  reps
}

stage2_specs <- function(specs, max_vertices = 100000L) {
  Filter(function(spec) !identical(spec$family, "sphere") && max(spec$edges) <= max_vertices, specs)
}

summarize_graph_results_general <- function(raw_metrics) {
  key <- paste(raw_metrics$family, raw_metrics$graph_label, raw_metrics$candidate_id, sep = "||")
  pieces <- split(raw_metrics, key)
  out <- do.call(rbind, lapply(pieces, function(df) {
    ok <- df$status == "ok"
    good <- df[ok, , drop = FALSE]
    best_seed <- if (nrow(good) > 0L) good$seed[[which.min(good$procrustes_rmse)]] else NA_integer_
    data.frame(
      family = df$family[[1L]],
      graph_label = df$graph_label[[1L]],
      title = df$title[[1L]],
      dim = df$dim[[1L]],
      candidate_id = df$candidate_id[[1L]],
      candidate_label = df$candidate_label[[1L]],
      vertices = df$vertices[[1L]],
      edges = df$edges[[1L]],
      n_runs = nrow(df),
      n_ok = sum(ok),
      n_fail = sum(!ok),
      seeds = paste(df$seed, collapse = ","),
      best_seed = best_seed,
      procrustes_rmse_mean = if (nrow(good) > 0L) mean(good$procrustes_rmse) else NA_real_,
      edge_length_cv_mean = if (nrow(good) > 0L) mean(good$edge_length_cv) else NA_real_,
      sampled_stress_mean = if (nrow(good) > 0L) mean(good$sampled_stress) else NA_real_,
      sampled_nonedge_sep_ratio_mean = if (nrow(good) > 0L) mean(good$sampled_nonedge_sep_ratio) else NA_real_,
      elapsed_sec_mean = if (nrow(good) > 0L) mean(good$elapsed_sec) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

score_graph_results_general <- function(graph_summary) {
  pieces <- split(graph_summary, paste(graph_summary$family, graph_summary$graph_label, sep = "||"))
  scored <- lapply(pieces, function(df) {
    df$rank_rmse <- rank01(df$procrustes_rmse_mean, higher_better = FALSE)
    df$rank_edge_cv <- rank01(df$edge_length_cv_mean, higher_better = FALSE)
    df$rank_stress <- rank01(df$sampled_stress_mean, higher_better = FALSE)
    df$rank_sep <- rank01(df$sampled_nonedge_sep_ratio_mean, higher_better = TRUE)
    df$rank_runtime <- rank01(df$elapsed_sec_mean, higher_better = FALSE)
    df$quality_score <- score_weights[["procrustes_rmse"]] * df$rank_rmse +
      score_weights[["edge_length_cv"]] * df$rank_edge_cv +
      score_weights[["sampled_stress"]] * df$rank_stress +
      score_weights[["sampled_nonedge_sep_ratio"]] * df$rank_sep
    df$value_score <- (1 - runtime_weight) * df$quality_score + runtime_weight * df$rank_runtime
    df[order(df$quality_score, df$procrustes_rmse_mean), , drop = FALSE]
  })
  out <- do.call(rbind, scored)
  rownames(out) <- NULL
  out
}

summarize_candidate_results_general <- function(scored_graphs) {
  pieces <- split(scored_graphs, scored_graphs$candidate_id)
  out <- do.call(rbind, lapply(pieces, function(df) {
    data.frame(
      candidate_id = df$candidate_id[[1L]],
      candidate_label = df$candidate_label[[1L]],
      n_graphs = nrow(df),
      n_fail_graphs = sum(df$n_ok == 0L),
      quality_score_mean = mean(df$quality_score, na.rm = TRUE),
      value_score_mean = mean(df$value_score, na.rm = TRUE),
      procrustes_rmse_mean = mean(df$procrustes_rmse_mean, na.rm = TRUE),
      edge_length_cv_mean = mean(df$edge_length_cv_mean, na.rm = TRUE),
      sampled_stress_mean = mean(df$sampled_stress_mean, na.rm = TRUE),
      sampled_nonedge_sep_ratio_mean = mean(df$sampled_nonedge_sep_ratio_mean, na.rm = TRUE),
      elapsed_sec_mean = mean(df$elapsed_sec_mean, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out[order(out$quality_score_mean), , drop = FALSE]
}

summarize_family_results_general <- function(scored_graphs) {
  pieces <- split(scored_graphs, paste(scored_graphs$family, scored_graphs$candidate_id, sep = "||"))
  out <- do.call(rbind, lapply(pieces, function(df) {
    data.frame(
      family = df$family[[1L]],
      candidate_id = df$candidate_id[[1L]],
      candidate_label = df$candidate_label[[1L]],
      n_graphs = nrow(df),
      n_fail_graphs = sum(df$n_ok == 0L),
      quality_score_mean = mean(df$quality_score, na.rm = TRUE),
      value_score_mean = mean(df$value_score, na.rm = TRUE),
      procrustes_rmse_mean = mean(df$procrustes_rmse_mean, na.rm = TRUE),
      edge_length_cv_mean = mean(df$edge_length_cv_mean, na.rm = TRUE),
      sampled_stress_mean = mean(df$sampled_stress_mean, na.rm = TRUE),
      sampled_nonedge_sep_ratio_mean = mean(df$sampled_nonedge_sep_ratio_mean, na.rm = TRUE),
      elapsed_sec_mean = mean(df$elapsed_sec_mean, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out[order(out$family, out$quality_score_mean), , drop = FALSE]
}

winner_by_family_general <- function(family_summary, score_col = "quality_score_mean") {
  pieces <- split(family_summary, family_summary$family)
  do.call(rbind, lapply(pieces, function(df) df[which.min(df[[score_col]]), , drop = FALSE]))
}

run_stage1 <- function(specs) {
  raw_rows <- list()
  idx <- 0L
  for (spec in specs) {
    seeds <- seed_budget_grip(max(spec$edges))
    message(sprintf("Stage 1: %s %s (%d vertices) with seeds %s",
                    spec$family, spec$graph_label, max(spec$edges), paste(seeds, collapse = ",")))
    for (candidate_id in all_universal_candidate_ids) {
      if (!candidate_applicable(candidate_id, spec)) next
      for (seed in seeds) {
        idx <- idx + 1L
        raw_rows[[idx]] <- run_one_layout_grip_safe(spec, candidate_id, seed)
      }
    }
  }
  raw_metrics <- do.call(rbind, raw_rows)
  graph_summary <- summarize_graph_results_general(raw_metrics)
  scored_graphs <- score_graph_results_general(graph_summary)
  candidate_summary <- summarize_candidate_results_general(scored_graphs)
  list(raw_metrics = raw_metrics, graph_summary = graph_summary, scored_graphs = scored_graphs, candidate_summary = candidate_summary)
}

select_stage2_universal_candidates <- function(stage1_candidate_summary, top_k = 4L) {
  ranked <- stage1_candidate_summary$candidate_id
  ranked <- ranked[ranked %in% all_universal_candidate_ids]
  selected <- unique(c(
    head(ranked, top_k),
    "grip_f32",
    "grip_default_adaptive"
  ))
  selected
}

stage2_candidate_ids_for_spec <- function(spec, stage2_universal_candidates) {
  n_vertices <- max(spec$edges)
  universal_candidates <- stage2_universal_candidates
  if (n_vertices > large_graph_vertex_threshold) {
    universal_candidates <- unique(c(
      head(stage2_universal_candidates, 3L),
      "grip_f32",
      "grip_default_adaptive"
    ))
  }
  unique(c(
    universal_candidates,
    "grip_mesh_preset_reference",
    "grip_torus_preset_reference",
    "grip_tree_preset_reference",
    "igraph_kk_default"
  ))
}

run_stage2 <- function(specs, stage2_universal_candidates, kk_max_vertices = 10000L) {
  raw_rows <- list()
  idx <- 0L
  for (spec in specs) {
    n_vertices <- max(spec$edges)
    message(sprintf("Stage 2: %s %s (%d vertices)", spec$family, spec$graph_label, n_vertices))
    candidate_ids <- stage2_candidate_ids_for_spec(spec, stage2_universal_candidates)
    for (candidate_id in candidate_ids) {
      if (!candidate_applicable(candidate_id, spec)) next
      seeds <- if (candidate_id == "igraph_kk_default") seed_budget_kk(n_vertices) else seed_budget_grip(n_vertices)
      for (seed in seeds) {
        idx <- idx + 1L
        raw_rows[[idx]] <- if (candidate_id == "igraph_kk_default") {
          run_one_layout_kk_safe(spec, seed, kk_max_vertices = kk_max_vertices)
        } else {
          run_one_layout_grip_safe(spec, candidate_id, seed)
        }
      }
    }
  }
  raw_metrics <- do.call(rbind, raw_rows)
  graph_summary <- summarize_graph_results_general(raw_metrics)
  scored_graphs <- score_graph_results_general(graph_summary)
  candidate_summary <- summarize_candidate_results_general(scored_graphs)
  family_summary <- summarize_family_results_general(scored_graphs)
  family_quality_winners <- winner_by_family_general(family_summary, score_col = "quality_score_mean")
  family_value_winners <- winner_by_family_general(family_summary, score_col = "value_score_mean")
  list(
    raw_metrics = raw_metrics,
    graph_summary = graph_summary,
    scored_graphs = scored_graphs,
    candidate_summary = candidate_summary,
    family_summary = family_summary,
    family_quality_winners = family_quality_winners,
    family_value_winners = family_value_winners
  )
}

draw_stage1_barplot <- function(path, candidate_summary) {
  df <- candidate_summary
  df <- df[order(df$quality_score_mean), , drop = FALSE]
  grDevices::png(path, width = 1800, height = 1000, res = 170, bg = "#f7f3ea")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(8, 5, 4, 1))
  graphics::barplot(
    height = df$quality_score_mean,
    names.arg = df$candidate_id,
    las = 2,
    col = "#355070",
    border = "#16324f",
    ylab = "Mean quality score (lower is better)",
    main = "Stage 1 broad screen: universal GRIP candidates"
  )
}

draw_stage2_scatter <- function(path, candidate_summary) {
  df <- candidate_summary
  cols <- ifelse(df$candidate_id == "igraph_kk_default", "#b23a48", "#355070")
  pchs <- ifelse(grepl("preset", df$candidate_id), 17, 19)
  grDevices::png(path, width = 1600, height = 1100, res = 170, bg = "#f7f3ea")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(5, 5, 4, 1))
  graphics::plot(
    df$elapsed_sec_mean, df$quality_score_mean,
    log = "x",
    pch = pchs,
    col = cols,
    cex = 1.4,
    xlab = "Mean elapsed sec (log scale)",
    ylab = "Mean quality score (lower is better)",
    main = "Stage 2 overall candidate tradeoff"
  )
  graphics::text(df$elapsed_sec_mean, df$quality_score_mean,
                 labels = df$candidate_id, pos = 3, cex = 0.85, col = "#16324f")
  graphics::legend("topright",
                   legend = c("GRIP candidate", "igraph KK", "family preset reference"),
                   col = c("#355070", "#b23a48", "#355070"),
                   pch = c(19, 19, 17),
                   bty = "n")
}

draw_family_heatmap <- function(path, family_summary) {
  fams <- unique(family_summary$family)
  cands <- unique(family_summary$candidate_id)
  mat <- matrix(NA_real_, nrow = length(fams), ncol = length(cands),
                dimnames = list(fams, cands))
  for (i in seq_len(nrow(family_summary))) {
    row <- family_summary[i, , drop = FALSE]
    mat[row$family, row$candidate_id] <- row$quality_score_mean
  }
  z <- mat
  finite_vals <- z[is.finite(z)]
  z_norm <- matrix(NA_real_, nrow = nrow(z), ncol = ncol(z), dimnames = dimnames(z))
  if (length(finite_vals) > 0L) {
    rng <- range(finite_vals)
    if (diff(rng) <= sqrt(.Machine$double.eps)) {
      z_norm[is.finite(z)] <- 0.5
    } else {
      z_norm[is.finite(z)] <- (z[is.finite(z)] - rng[[1L]]) / diff(rng)
    }
  }
  palette_fun <- grDevices::colorRampPalette(c("#f7f3ea", "#c5d6e8", "#355070"))
  pal <- palette_fun(100)

  grDevices::png(path, width = 1800, height = 1200, res = 170, bg = "#f7f3ea")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(9, 8, 4, 2))
  graphics::plot(c(0, ncol(z_norm)), c(0, nrow(z_norm)),
                 type = "n", xlab = "", ylab = "", axes = FALSE,
                 main = "Stage 2 family-by-candidate quality heatmap")
  for (i in seq_len(nrow(z_norm))) {
    for (j in seq_len(ncol(z_norm))) {
      val <- z_norm[i, j]
      col <- if (is.na(val)) "#d9d9d9" else pal[[max(1L, min(100L, floor(val * 99L) + 1L))]]
      graphics::rect(j - 1, nrow(z_norm) - i, j, nrow(z_norm) - i + 1, col = col, border = "#f7f3ea")
      if (is.finite(mat[i, j])) {
        graphics::text(j - 0.5, nrow(z_norm) - i + 0.5,
                       labels = sprintf("%.3f", mat[i, j]), cex = 0.7, col = "#16324f")
      }
    }
  }
  graphics::axis(1, at = seq_len(ncol(z_norm)) - 0.5, labels = cands, las = 2)
  graphics::axis(2, at = rev(seq_len(nrow(z_norm)) - 0.5), labels = fams, las = 2)
  graphics::box()
}

write_stage1_summary <- function(path, representative_specs, candidate_summary, selected_candidates) {
  lines <- c(
    "# Stage 1 broad screen",
    "",
    "Representative graph set:",
    vapply(representative_specs, function(spec) {
      sprintf("- `%s %s` (%d vertices, dim=%d)", spec$family, spec$graph_label, max(spec$edges), spec$dim)
    }, character(1L)),
    "",
    "| Candidate | Mean quality score | RMSE | Edge CV | Stress | Non-edge sep | Mean sec |",
    "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"
  )
  for (i in seq_len(nrow(candidate_summary))) {
    row <- candidate_summary[i, , drop = FALSE]
    lines <- c(lines, sprintf(
      "| %s | %.4f | %.4f | %.4f | %.4f | %.4f | %.3f |",
      row$candidate_id,
      row$quality_score_mean,
      row$procrustes_rmse_mean,
      row$edge_length_cv_mean,
      row$sampled_stress_mean,
      row$sampled_nonedge_sep_ratio_mean,
      row$elapsed_sec_mean
    ))
  }
  lines <- c(
    lines,
    "",
    sprintf("Selected universal candidates for Stage 2: `%s`", paste(selected_candidates, collapse = "`, `"))
  )
  writeLines(lines, con = path)
}

write_stage2_summary <- function(path,
                                 config,
                                 stage2_result,
                                 global_quality_winner,
                                 global_value_winner,
                                 comparable_subset_summary) {
  cs <- comparable_subset_summary
  lines <- c(
    "# Stage 2 full benchmark",
    "",
    sprintf("- stage-2 graphs benchmarked: `%d`", nrow(unique(stage2_result$graph_summary[c("family", "graph_label")]))),
    sprintf("- graph vertex cap: `%d`", config$max_vertices),
    sprintf("- KK vertex cap: `%d`", config$kk_max_vertices),
    "",
    sprintf("- global quality winner: `%s` (mean quality score %.4f, mean elapsed %.3f sec)",
            global_quality_winner$candidate_id,
            global_quality_winner$quality_score_mean,
            global_quality_winner$elapsed_sec_mean),
    sprintf("- global value winner: `%s` (mean value score %.4f, mean elapsed %.3f sec)",
            global_value_winner$candidate_id,
            global_value_winner$value_score_mean,
            global_value_winner$elapsed_sec_mean),
    "",
    "## Candidate summary",
    "",
    "| Candidate | Mean quality score | Mean value score | RMSE | Edge CV | Stress | Non-edge sep | Mean sec | Failed graphs |",
    "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"
  )
  for (i in seq_len(nrow(stage2_result$candidate_summary))) {
    row <- stage2_result$candidate_summary[i, , drop = FALSE]
    lines <- c(lines, sprintf(
      "| %s | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.3f | %d |",
      row$candidate_id,
      row$quality_score_mean,
      row$value_score_mean,
      row$procrustes_rmse_mean,
      row$edge_length_cv_mean,
      row$sampled_stress_mean,
      row$sampled_nonedge_sep_ratio_mean,
      row$elapsed_sec_mean,
      row$n_fail_graphs
    ))
  }

  lines <- c(lines, "", "## Family quality winners", "")
  for (i in seq_len(nrow(stage2_result$family_quality_winners))) {
    row <- stage2_result$family_quality_winners[i, , drop = FALSE]
    lines <- c(lines, sprintf(
      "- `%s`: `%s` (mean quality score %.4f, mean elapsed %.3f sec)",
      row$family, row$candidate_id, row$quality_score_mean, row$elapsed_sec_mean
    ))
  }

  lines <- c(lines, "", "## Comparable subset: grip_f32 vs igraph KK", "")
  if (nrow(cs) == 0L) {
    lines <- c(lines, "- no graphs where both methods completed successfully")
  } else {
    lines <- c(
      lines,
      sprintf("- graphs compared: `%d`", nrow(cs)),
      sprintf("- grip_f32 mean RMSE: `%.4f`", mean(cs$grip_f32_rmse)),
      sprintf("- igraph KK mean RMSE: `%.4f`", mean(cs$igraph_kk_rmse)),
      sprintf("- grip_f32 mean elapsed sec: `%.3f`", mean(cs$grip_f32_elapsed)),
      sprintf("- igraph KK mean elapsed sec: `%.3f`", mean(cs$igraph_kk_elapsed))
    )
  }

  writeLines(lines, con = path)
}

build_comparable_subset_summary <- function(graph_summary) {
  pieces <- split(graph_summary, paste(graph_summary$family, graph_summary$graph_label, sep = "||"))
  rows <- lapply(pieces, function(df) {
    grip_row <- df[df$candidate_id == "grip_f32" & df$n_ok > 0L, , drop = FALSE]
    kk_row <- df[df$candidate_id == "igraph_kk_default" & df$n_ok > 0L, , drop = FALSE]
    if (nrow(grip_row) == 0L || nrow(kk_row) == 0L) {
      return(NULL)
    }
    data.frame(
      family = grip_row$family,
      graph_label = grip_row$graph_label,
      grip_f32_rmse = grip_row$procrustes_rmse_mean,
      igraph_kk_rmse = kk_row$procrustes_rmse_mean,
      grip_f32_elapsed = grip_row$elapsed_sec_mean,
      igraph_kk_elapsed = kk_row$elapsed_sec_mean,
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

write_tracked_summary <- function(path,
                                  stage1_result,
                                  selected_candidates,
                                  stage2_result,
                                  global_quality_winner,
                                  global_value_winner,
                                  comparable_subset_summary,
                                  run_tag) {
  cs <- comparable_subset_summary
  lines <- c(
    "# GRIP vs KK Cross-Family Benchmark Summary",
    "",
    sprintf("Run tag: `%s`", run_tag),
    "",
    "## Stage 1",
    "",
    sprintf("- screened universal candidates: `%s`", paste(all_universal_candidate_ids, collapse = "`, `")),
    sprintf("- selected universal candidates for Stage 2: `%s`", paste(selected_candidates, collapse = "`, `")),
    sprintf("- Stage-1 winner: `%s` (mean quality score %.4f)", stage1_result$candidate_summary$candidate_id[[1L]], stage1_result$candidate_summary$quality_score_mean[[1L]]),
    "",
    "## Stage 2",
    "",
    sprintf("- global quality winner: `%s` (mean quality score %.4f)", global_quality_winner$candidate_id, global_quality_winner$quality_score_mean),
    sprintf("- global value winner: `%s` (mean value score %.4f)", global_value_winner$candidate_id, global_value_winner$value_score_mean),
    "",
    "Family quality winners:",
    vapply(seq_len(nrow(stage2_result$family_quality_winners)), function(i) {
      row <- stage2_result$family_quality_winners[i, , drop = FALSE]
      sprintf("- `%s`: `%s`", row$family, row$candidate_id)
    }, character(1L))
  )

  if (nrow(cs) > 0L) {
    lines <- c(
      lines,
      "",
      "Comparable subset (`grip_f32` vs `igraph_kk_default`):",
      sprintf("- graphs where both completed: `%d`", nrow(cs)),
      sprintf("- mean RMSE: grip_f32 `%.4f`, igraph KK `%.4f`", mean(cs$grip_f32_rmse), mean(cs$igraph_kk_rmse)),
      sprintf("- mean elapsed sec: grip_f32 `%.3f`, igraph KK `%.3f`", mean(cs$grip_f32_elapsed), mean(cs$igraph_kk_elapsed))
    )
  }

  lines <- c(
    lines,
    "",
    "Primary outputs:",
    sprintf("- Stage-1 summary: `dev/design/tmp/%s/stage1-summary.md`", run_tag),
    sprintf("- Stage-2 summary: `dev/design/tmp/%s/stage2-summary.md`", run_tag),
    sprintf("- Stage-2 candidate CSV: `dev/design/tmp/%s/stage2-candidate-summary.csv`", run_tag),
    sprintf("- Stage-2 family CSV: `dev/design/tmp/%s/stage2-family-summary.csv`", run_tag)
  )
  writeLines(lines, con = path)
}

if (sys.nframe() == 0L) {
  args <- parse_named_args(commandArgs(trailingOnly = TRUE))
  run_tag <- validate_run_tag(if (!is.null(args$tag)) args$tag else "grip-kk-cross-family-2026-03-28")
  max_vertices <- if (!is.null(args$max_vertices)) parse_int_scalar(args$max_vertices, "max_vertices") else 100000L
  kk_max_vertices <- if (!is.null(args$kk_max_vertices)) parse_int_scalar(args$kk_max_vertices, "kk_max_vertices") else 10000L
  representative_max_vertices <- if (!is.null(args$representative_max_vertices)) parse_int_scalar(args$representative_max_vertices, "representative_max_vertices") else 10000L

  out_root <- file.path("dev", "manual")
  out_tmp_dir <- file.path(out_root, "tmp", run_tag)
  dir.create(out_tmp_dir, recursive = TRUE, showWarnings = FALSE)

  all_specs <- build_graph_specs()
  rep_specs <- stage1_representative_specs(all_specs, representative_max_vertices = representative_max_vertices)
  stage1_result <- run_stage1(rep_specs)
  selected_candidates <- select_stage2_universal_candidates(stage1_result$candidate_summary, top_k = 4L)

  stage2_all_specs <- stage2_specs(all_specs, max_vertices = max_vertices)
  stage2_result <- run_stage2(stage2_all_specs, selected_candidates, kk_max_vertices = kk_max_vertices)

  global_quality_winner <- stage2_result$candidate_summary[which.min(stage2_result$candidate_summary$quality_score_mean), , drop = FALSE]
  global_value_winner <- stage2_result$candidate_summary[which.min(stage2_result$candidate_summary$value_score_mean), , drop = FALSE]
  comparable_subset_summary <- build_comparable_subset_summary(stage2_result$graph_summary)

  stage1_raw_csv <- file.path(out_tmp_dir, "stage1-raw-metrics.csv")
  stage1_graph_csv <- file.path(out_tmp_dir, "stage1-graph-summary.csv")
  stage1_candidate_csv <- file.path(out_tmp_dir, "stage1-candidate-summary.csv")
  stage1_summary_md <- file.path(out_tmp_dir, "stage1-summary.md")
  stage1_plot_png <- file.path(out_tmp_dir, "stage1-candidate-ranking.png")
  stage2_raw_csv <- file.path(out_tmp_dir, "stage2-raw-metrics.csv")
  stage2_graph_csv <- file.path(out_tmp_dir, "stage2-graph-summary.csv")
  stage2_candidate_csv <- file.path(out_tmp_dir, "stage2-candidate-summary.csv")
  stage2_family_csv <- file.path(out_tmp_dir, "stage2-family-summary.csv")
  stage2_quality_winners_csv <- file.path(out_tmp_dir, "stage2-family-quality-winners.csv")
  stage2_value_winners_csv <- file.path(out_tmp_dir, "stage2-family-value-winners.csv")
  stage2_summary_md <- file.path(out_tmp_dir, "stage2-summary.md")
  stage2_scatter_png <- file.path(out_tmp_dir, "stage2-candidate-tradeoff.png")
  stage2_heatmap_png <- file.path(out_tmp_dir, "stage2-family-heatmap.png")
  comparable_csv <- file.path(out_tmp_dir, "stage2-grip-f32-vs-kk-comparable.csv")

  utils::write.csv(stage1_result$raw_metrics, stage1_raw_csv, row.names = FALSE)
  utils::write.csv(stage1_result$graph_summary, stage1_graph_csv, row.names = FALSE)
  utils::write.csv(stage1_result$candidate_summary, stage1_candidate_csv, row.names = FALSE)
  write_stage1_summary(stage1_summary_md, rep_specs, stage1_result$candidate_summary, selected_candidates)
  draw_stage1_barplot(stage1_plot_png, stage1_result$candidate_summary)

  utils::write.csv(stage2_result$raw_metrics, stage2_raw_csv, row.names = FALSE)
  utils::write.csv(stage2_result$graph_summary, stage2_graph_csv, row.names = FALSE)
  utils::write.csv(stage2_result$candidate_summary, stage2_candidate_csv, row.names = FALSE)
  utils::write.csv(stage2_result$family_summary, stage2_family_csv, row.names = FALSE)
  utils::write.csv(stage2_result$family_quality_winners, stage2_quality_winners_csv, row.names = FALSE)
  utils::write.csv(stage2_result$family_value_winners, stage2_value_winners_csv, row.names = FALSE)
  utils::write.csv(comparable_subset_summary, comparable_csv, row.names = FALSE)
  write_stage2_summary(
    path = stage2_summary_md,
    config = list(max_vertices = max_vertices, kk_max_vertices = kk_max_vertices),
    stage2_result = stage2_result,
    global_quality_winner = global_quality_winner,
    global_value_winner = global_value_winner,
    comparable_subset_summary = comparable_subset_summary
  )
  draw_stage2_scatter(stage2_scatter_png, stage2_result$candidate_summary)
  draw_family_heatmap(stage2_heatmap_png, stage2_result$family_summary)

  tracked_summary_path <- file.path("dev", "manual", "grip_kk_cross_family_summary_2026-03-28.md")
  write_tracked_summary(
    path = tracked_summary_path,
    stage1_result = stage1_result,
    selected_candidates = selected_candidates,
    stage2_result = stage2_result,
    global_quality_winner = global_quality_winner,
    global_value_winner = global_value_winner,
    comparable_subset_summary = comparable_subset_summary,
    run_tag = run_tag
  )

  message(sprintf("Stage 1 summary written to %s", stage1_summary_md))
  message(sprintf("Stage 2 summary written to %s", stage2_summary_md))
  message(sprintf("Tracked summary written to %s", tracked_summary_path))
}
