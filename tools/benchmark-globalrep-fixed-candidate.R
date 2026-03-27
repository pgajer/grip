#!/usr/bin/env Rscript

sys.source(file.path("tools", "benchmark-sierpinski-baseline.R"), envir = environment())
sys.source(file.path("tools", "generate-sierpinski-diagnostics.R"), envir = environment())
sys.source(file.path("tools", "tune-torus-parameters.R"), envir = environment())

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

score_weights <- c(
  procrustes_rmse = 0.45,
  edge_length_cv = 0.20,
  sampled_stress = 0.20,
  sampled_nonedge_sep_ratio = 0.15
)

package_globalrep_default_cfg <- function(n) {
  get("grip.globalrep.base.defaults", envir = asNamespace("grip"))(n = n)
}

fixed_candidate_cfg <- list(
  placement = "barycenter",
  rounds = 160L,
  final_rounds = 384L,
  num_init = 24L,
  num_nbrs = 20L,
  r = 0.03,
  s = 7.5,
  repulsion_factor = 2.5,
  coarse_repulsion_factor = 1.5,
  coarse_repulsion_sample = 16L,
  coarse_repulsion_exact_below = 64L
)

candidate_ids <- c("globalrep_default", "globalrep_fixed_candidate")

candidate_label <- function(candidate_id) {
  switch(
    candidate_id,
    globalrep_default = "globalrep defaults",
    globalrep_fixed_candidate = "fixed candidate",
    candidate_id
  )
}

resolve_candidate_cfg <- function(candidate_id, n) {
  switch(
    candidate_id,
    globalrep_default = package_globalrep_default_cfg(n),
    globalrep_fixed_candidate = fixed_candidate_cfg,
    stop(sprintf("Unknown candidate_id: %s", candidate_id))
  )
}

format_cfg <- function(cfg) {
  paste(names(cfg), unlist(cfg), sep = "=", collapse = ", ")
}

align_to_canonical_nd <- function(source, target) {
  fit <- align_to_target_orthogonal(source, target, allow_reflection = TRUE)
  list(
    aligned = fit$aligned,
    rmse = fit$rmse
  )
}

graph_seed_budget <- function(n_vertices, base_seeds) {
  if (n_vertices <= 2000L) {
    return(base_seeds)
  }
  if (n_vertices <= 20000L) {
    return(base_seeds[seq_len(min(2L, length(base_seeds)))])
  }
  base_seeds[[1L]]
}

sample_budget <- function(n_vertices) {
  if (n_vertices <= 5000L) {
    return(list(stress = 3000L, sep = 6000L))
  }
  if (n_vertices <= 30000L) {
    return(list(stress = 1500L, sep = 3000L))
  }
  list(stress = 750L, sep = 1500L)
}

time_limit_budget <- function(n_vertices, base_limit_sec) {
  if (n_vertices <= 5000L) {
    return(base_limit_sec)
  }
  if (n_vertices <= 30000L) {
    return(max(base_limit_sec, 180L))
  }
  max(base_limit_sec, 300L)
}

build_path_spec <- function(n) {
  edges <- edges.path(n)
  coords <- cbind(
    x = seq_len(n) - 1L,
    y = rep(0, n)
  )
  list(
    family = "path",
    graph_label = sprintf("n%d", n),
    title = sprintf("Path n=%d", n),
    dim = 2L,
    edges = edges,
    canonical = coords,
    aligner = "nd"
  )
}

build_cycle_spec <- function(n) {
  edges <- edges.cycle(n)
  theta <- seq(0, 2 * pi, length.out = n + 1L)[seq_len(n)]
  coords <- cbind(
    x = cos(theta),
    y = sin(theta)
  )
  list(
    family = "cycle",
    graph_label = sprintf("n%d", n),
    title = sprintf("Cycle n=%d", n),
    dim = 2L,
    edges = edges,
    canonical = coords,
    aligner = "nd"
  )
}

build_triangle_spec <- function(level) {
  built <- build_sierpinski_triangle(level)
  list(
    family = "sierpinski.triangle",
    graph_label = sprintf("level%d", level),
    title = sprintf("Sierpinski triangle level %d", level),
    dim = 2L,
    edges = built$edges,
    canonical = built$coords,
    aligner = "nd"
  )
}

build_carpet_spec <- function(level) {
  built <- build_sierpinski_carpet(level)
  list(
    family = "sierpinski.carpet",
    graph_label = sprintf("level%d", level),
    title = sprintf("Sierpinski carpet level %d", level),
    dim = 2L,
    edges = built$edges,
    canonical = built$coords,
    aligner = "nd"
  )
}

build_tetrahedron_spec <- function(level) {
  built <- build_sierpinski_tetrahedron(level)
  list(
    family = "sierpinski.tetrahedron",
    graph_label = sprintf("level%d", level),
    title = sprintf("Sierpinski tetrahedron level %d", level),
    dim = 3L,
    edges = built$edges,
    canonical = built$coords,
    aligner = "nd"
  )
}

build_mesh_spec_local <- function(h, w = h) {
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

build_cylinder_spec_local <- function(h, w) {
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
  list(
    family = "cylinder",
    graph_label = sprintf("%dx%d", h, w),
    title = sprintf("Cylinder %dx%d", h, w),
    dim = 3L,
    edges = edges,
    canonical = coords,
    aligner = "nd"
  )
}

build_sphere_graph <- function(h, w) {
  h <- as.integer(h)
  w <- as.integer(w)
  edges <- edges.sphere(h, w)
  ring.count <- h - 2L
  coords <- matrix(0, nrow = 2L + ring.count * w, ncol = 3L)
  coords[1L, ] <- c(0, 0, 1)
  idx <- 2L
  for (i in seq_len(ring.count)) {
    theta <- pi * i / (h - 1L)
    for (j in seq_len(w)) {
      phi <- 2 * pi * (j - 1L) / w
      coords[idx, ] <- c(
        sin(theta) * cos(phi),
        sin(theta) * sin(phi),
        cos(theta)
      )
      idx <- idx + 1L
    }
  }
  coords[idx, ] <- c(0, 0, -1)
  list(edges = edges, coords = coords)
}

build_sphere_spec <- function(h, w) {
  built <- build_sphere_graph(h, w)
  list(
    family = "sphere",
    graph_label = sprintf("%dx%d", h, w),
    title = sprintf("Sphere %dx%d", h, w),
    dim = 3L,
    edges = built$edges,
    canonical = built$coords,
    aligner = "nd"
  )
}

build_cube_graph <- function(side) {
  side <- as.integer(side)
  edges <- edges.cube(side)
  grid <- expand.grid(
    x = seq_len(side),
    y = seq_len(side),
    z = seq_len(side)
  )
  keep <- grid$x %in% c(1L, side) |
    grid$y %in% c(1L, side) |
    grid$z %in% c(1L, side)
  coords <- as.matrix(grid[keep, , drop = FALSE])
  coords <- scale(coords, center = TRUE, scale = FALSE)
  list(edges = edges, coords = coords)
}

build_cube_spec <- function(side) {
  built <- build_cube_graph(side)
  list(
    family = "cube",
    graph_label = sprintf("side%d", side),
    title = sprintf("Cube side %d", side),
    dim = 3L,
    edges = built$edges,
    canonical = built$coords,
    aligner = "nd"
  )
}

build_kary_tree_spec_local <- function(k, depth) {
  edges <- edges.kary.tree(k = k, depth = depth)
  n <- max(edges)
  adj <- make_adj_list(edges, n)
  parent <- integer(n)
  level <- integer(n)
  order <- integer(n)
  head <- 1L
  tail <- 1L
  order[[tail]] <- 1L
  while (head <= tail) {
    v <- order[[head]]
    head <- head + 1L
    kids <- sort(adj[[v]][adj[[v]] != parent[[v]]])
    if (length(kids) == 0L) {
      next
    }
    for (u in kids) {
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
    for (u in kids) {
      place(u)
    }
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
    title = sprintf("%d-ary tree depth %d", k, depth),
    dim = 2L,
    edges = edges,
    canonical = coords,
    aligner = "nd"
  )
}

build_graph_specs <- function() {
  specs <- list()

  for (level in 2:6) {
    specs[[length(specs) + 1L]] <- build_carpet_spec(level)
    specs[[length(specs) + 1L]] <- build_triangle_spec(level)
    specs[[length(specs) + 1L]] <- build_tetrahedron_spec(level)
  }

  for (size in c(8L, 10L, 12L)) {
    specs[[length(specs) + 1L]] <- build_mesh_spec_local(size, size)
    specs[[length(specs) + 1L]] <- build_cylinder_spec_local(size, size)
  }

  for (size in default_torus_sizes()) {
    built <- build_graph_specs_for_sizes(list(size))[[1L]]
    built$title <- sprintf("Torus %s", built$graph_label)
    built$dim <- 3L
    built$aligner <- "torus"
    specs[[length(specs) + 1L]] <- built
  }

  for (size in c(32L, 64L, 128L)) {
    specs[[length(specs) + 1L]] <- build_path_spec(size)
    specs[[length(specs) + 1L]] <- build_cycle_spec(size)
  }

  sphere_sizes <- list(
    list(h = 6L, w = 12L),
    list(h = 8L, w = 16L),
    list(h = 10L, w = 20L)
  )
  for (size in sphere_sizes) {
    specs[[length(specs) + 1L]] <- build_sphere_spec(size$h, size$w)
  }

  for (side in c(4L, 6L, 8L)) {
    specs[[length(specs) + 1L]] <- build_cube_spec(side)
  }

  specs[[length(specs) + 1L]] <- build_kary_tree_spec_local(2L, 4L)
  specs[[length(specs) + 1L]] <- build_kary_tree_spec_local(2L, 5L)
  specs[[length(specs) + 1L]] <- build_kary_tree_spec_local(2L, 6L)
  specs[[length(specs) + 1L]] <- build_kary_tree_spec_local(3L, 4L)

  specs
}

align_layout_to_spec <- function(coords, spec) {
  if (identical(spec$aligner, "torus")) {
    return(best_torus_alignment(coords, spec))
  }
  align_to_canonical_nd(coords, spec$canonical)
}

run_one_layout_globalrep <- function(spec,
                                     cfg,
                                     candidate_id,
                                     seed,
                                     base_time_limit_sec) {
  n <- max(spec$edges)
  adj <- make_adj_list(spec$edges, n)
  samples <- sample_budget(n)
  time_limit_sec <- time_limit_budget(n, base_time_limit_sec)
  started <- proc.time()[["elapsed"]]
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE), add = TRUE)
  setTimeLimit(elapsed = time_limit_sec, transient = TRUE)

  coords <- grip.layout.globalrep(
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
    coarse_repulsion_factor = cfg$coarse_repulsion_factor,
    coarse_repulsion_sample = cfg$coarse_repulsion_sample,
    coarse_repulsion_exact_below = cfg$coarse_repulsion_exact_below,
    seed = seed
  )

  elapsed <- proc.time()[["elapsed"]] - started
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

run_one_layout_globalrep_safe <- function(spec,
                                          cfg,
                                          candidate_id,
                                          seed,
                                          base_time_limit_sec) {
  tryCatch(
    run_one_layout_globalrep(spec, cfg, candidate_id, seed, base_time_limit_sec),
    error = function(e) {
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
        status = "error",
        error_message = conditionMessage(e),
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
}

make_skipped_row <- function(spec, cfg, candidate_id, seed, reason) {
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
    status = "skipped",
    error_message = reason,
    procrustes_rmse = NA_real_,
    edge_length_cv = NA_real_,
    median_edge_length = NA_real_,
    sampled_stress = NA_real_,
    sampled_nonedge_sep_ratio = NA_real_,
    elapsed_sec = NA_real_,
    stringsAsFactors = FALSE
  )
}

summarize_graph_results <- function(raw_metrics) {
  key <- paste(raw_metrics$family, raw_metrics$graph_label, raw_metrics$candidate_id, sep = "||")
  rows <- split(raw_metrics, key)
  out <- do.call(rbind, lapply(rows, function(df) {
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

score_graph_results <- function(graph_summary) {
  pieces <- split(graph_summary, paste(graph_summary$family, graph_summary$graph_label, sep = "||"))
  scored <- lapply(pieces, function(df) {
    df$rank_rmse <- rank01(df$procrustes_rmse_mean, higher_better = FALSE)
    df$rank_edge_cv <- rank01(df$edge_length_cv_mean, higher_better = FALSE)
    df$rank_stress <- rank01(df$sampled_stress_mean, higher_better = FALSE)
    df$rank_sep <- rank01(df$sampled_nonedge_sep_ratio_mean, higher_better = TRUE)
    df$score_graph <- score_weights[["procrustes_rmse"]] * df$rank_rmse +
      score_weights[["edge_length_cv"]] * df$rank_edge_cv +
      score_weights[["sampled_stress"]] * df$rank_stress +
      score_weights[["sampled_nonedge_sep_ratio"]] * df$rank_sep
    df[order(df$score_graph, df$procrustes_rmse_mean), , drop = FALSE]
  })
  out <- do.call(rbind, scored)
  rownames(out) <- NULL
  out
}

summarize_family_results <- function(scored_graphs) {
  pieces <- split(scored_graphs, paste(scored_graphs$family, scored_graphs$candidate_id, sep = "||"))
  out <- do.call(rbind, lapply(pieces, function(df) {
    data.frame(
      family = df$family[[1L]],
      candidate_id = df$candidate_id[[1L]],
      candidate_label = df$candidate_label[[1L]],
      n_graphs = nrow(df),
      n_fail_graphs = sum(df$n_ok == 0L),
      score_graph_mean = mean(df$score_graph, na.rm = TRUE),
      procrustes_rmse_mean = mean(df$procrustes_rmse_mean, na.rm = TRUE),
      edge_length_cv_mean = mean(df$edge_length_cv_mean, na.rm = TRUE),
      sampled_stress_mean = mean(df$sampled_stress_mean, na.rm = TRUE),
      sampled_nonedge_sep_ratio_mean = mean(df$sampled_nonedge_sep_ratio_mean, na.rm = TRUE),
      elapsed_sec_mean = mean(df$elapsed_sec_mean, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out[order(out$family, out$score_graph_mean), , drop = FALSE]
}

summarize_candidate_results <- function(scored_graphs) {
  graph_winner_ids <- vapply(
    split(scored_graphs, paste(scored_graphs$family, scored_graphs$graph_label, sep = "||")),
    function(df) {
      ok <- is.finite(df$score_graph)
      if (!any(ok)) {
        return(NA_character_)
      }
      df$candidate_id[ok][[which.min(df$score_graph[ok])]]
    },
    character(1L)
  )
  graph_win_table <- table(graph_winner_ids[!is.na(graph_winner_ids)])
  pieces <- split(scored_graphs, scored_graphs$candidate_id)
  out <- do.call(rbind, lapply(pieces, function(df) {
    candidate_id <- df$candidate_id[[1L]]
    data.frame(
      candidate_id = candidate_id,
      candidate_label = df$candidate_label[[1L]],
      n_graphs = nrow(df),
      n_fail_graphs = sum(df$n_ok == 0L),
      graph_wins = if (candidate_id %in% names(graph_win_table)) as.integer(graph_win_table[[candidate_id]]) else 0L,
      score_graph_mean = mean(df$score_graph, na.rm = TRUE),
      procrustes_rmse_mean = mean(df$procrustes_rmse_mean, na.rm = TRUE),
      edge_length_cv_mean = mean(df$edge_length_cv_mean, na.rm = TRUE),
      sampled_stress_mean = mean(df$sampled_stress_mean, na.rm = TRUE),
      sampled_nonedge_sep_ratio_mean = mean(df$sampled_nonedge_sep_ratio_mean, na.rm = TRUE),
      elapsed_sec_mean = mean(df$elapsed_sec_mean, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out[order(out$score_graph_mean), , drop = FALSE]
}

winner_by_family <- function(family_summary) {
  pieces <- split(family_summary, family_summary$family)
  do.call(rbind, lapply(pieces, function(df) df[which.min(df$score_graph_mean), , drop = FALSE]))
}

render_representative_figures <- function(specs,
                                          raw_metrics,
                                          pdf_dir,
                                          preview_dir) {
  pdf_paths <- character()
  families <- unique(vapply(specs, `[[`, character(1L), "family"))
  for (family in families) {
    family_rows <- raw_metrics[raw_metrics$family == family & raw_metrics$status == "ok", , drop = FALSE]
    if (nrow(family_rows) == 0L) {
      next
    }
    graph_sizes <- aggregate(vertices ~ graph_label + title, data = family_rows, FUN = max)
    graph_sizes <- graph_sizes[order(-graph_sizes$vertices), , drop = FALSE]
    graph_label <- graph_sizes$graph_label[[1L]]
    spec_idx <- which(vapply(specs, function(x) identical(x$family, family) && identical(x$graph_label, graph_label), logical(1L)))[[1L]]
    spec <- specs[[spec_idx]]
    candidate_rows <- family_rows[family_rows$graph_label == graph_label, , drop = FALSE]
    seed <- sort(unique(candidate_rows$seed))[[1L]]
    n <- max(spec$edges)
    default_cfg <- resolve_candidate_cfg("globalrep_default", n)
    fixed_cfg <- resolve_candidate_cfg("globalrep_fixed_candidate", n)

    default_coords <- grip.layout.globalrep(
      edges = spec$edges,
      n = n,
      dim = spec$dim,
      placement = default_cfg$placement,
      rounds = default_cfg$rounds,
      final_rounds = default_cfg$final_rounds,
      num_init = default_cfg$num_init,
      num_nbrs = default_cfg$num_nbrs,
      r = default_cfg$r,
      s = default_cfg$s,
      repulsion_factor = default_cfg$repulsion_factor,
      coarse_repulsion_factor = default_cfg$coarse_repulsion_factor,
      coarse_repulsion_sample = default_cfg$coarse_repulsion_sample,
      coarse_repulsion_exact_below = default_cfg$coarse_repulsion_exact_below,
      seed = seed
    )
    fixed_coords <- grip.layout.globalrep(
      edges = spec$edges,
      n = n,
      dim = spec$dim,
      placement = fixed_cfg$placement,
      rounds = fixed_cfg$rounds,
      final_rounds = fixed_cfg$final_rounds,
      num_init = fixed_cfg$num_init,
      num_nbrs = fixed_cfg$num_nbrs,
      r = fixed_cfg$r,
      s = fixed_cfg$s,
      repulsion_factor = fixed_cfg$repulsion_factor,
      coarse_repulsion_factor = fixed_cfg$coarse_repulsion_factor,
      coarse_repulsion_sample = fixed_cfg$coarse_repulsion_sample,
      coarse_repulsion_exact_below = fixed_cfg$coarse_repulsion_exact_below,
      seed = seed
    )

    default_fit <- align_layout_to_spec(default_coords, spec)
    fixed_fit <- align_layout_to_spec(fixed_coords, spec)
    pdf_path <- file.path(pdf_dir, family, sprintf("%s-default-vs-fixed-candidate.pdf", graph_label))
    dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)

    subtitle <- sprintf(
      "seed=%d | default vs fixed candidate | %d vertices",
      seed,
      n
    )
    if (spec$dim == 2L) {
      plot_2d_triptych(
        path = pdf_path,
        canonical_coords = spec$canonical,
        baseline_coords = default_fit$aligned,
        tuned_coords = fixed_fit$aligned,
        edges = spec$edges,
        title_text = spec$title,
        subtitle_text = subtitle,
        tuned_label = "fixed candidate"
      )
    } else {
      write_3d_diagnostic_pdf(
        path = pdf_path,
        canonical_coords = spec$canonical,
        baseline_coords = default_fit$aligned,
        tuned_coords = fixed_fit$aligned,
        edges = spec$edges,
        title_text = spec$title,
        subtitle_text = subtitle,
        tuned_label = "fixed candidate"
      )
    }
    pdf_paths <- c(pdf_paths, pdf_path)
  }

  if (length(pdf_paths) > 0L) {
    render_pdf_previews(pdf_paths, preview_dir)
  }
  pdf_paths
}

write_summary <- function(path,
                          config,
                          candidate_summary,
                          family_summary,
                          family_winners,
                          representative_pdfs) {
  lines <- c(
    "# Globalrep fixed-candidate benchmark",
    "",
    "Benchmark setup:",
    sprintf("- base seeds: `%s`", paste(config$base_seeds, collapse = ", ")),
    sprintf("- base per-run time limit: `%d` seconds", config$base_time_limit_sec),
    sprintf("- max benchmarked vertices per graph: `%d`", config$max_vertices),
    "- adaptive policy: graphs with <= 2,000 vertices use all base seeds; <= 20,000 use the first two seeds; larger graphs use the first seed only",
    "- adaptive policy: stress/non-edge sample sizes and time limits grow looser for larger graphs",
    "- graphs above the max-vertex budget are recorded as skipped so one ultra-large case cannot stall the full report",
    "",
    "Candidates:",
    sprintf("- globalrep defaults: quality-first package defaults with adaptive `final_rounds` schedule `<=1000 -> 384`, `<=5000 -> 320`, `<=20000 -> 256`, `<=50000 -> 200`, `>50000 -> 128`; shared base profile `%s`", format_cfg(package_globalrep_default_cfg(1000L))),
    sprintf("- fixed candidate: `%s`", format_cfg(fixed_candidate_cfg)),
    "",
    "## Overall summary",
    "",
    "| Candidate | Mean graph score | RMSE | Edge CV | Stress | Non-edge sep | Mean elapsed sec | Failed graphs |",
    "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"
  )

  for (i in seq_len(nrow(candidate_summary))) {
    row <- candidate_summary[i, , drop = FALSE]
    lines <- c(lines, sprintf(
      "| %s | %.4f | %.4f | %.4f | %.4f | %.4f | %.3f | %d |",
      row$candidate_label[[1L]],
      row$score_graph_mean[[1L]],
      row$procrustes_rmse_mean[[1L]],
      row$edge_length_cv_mean[[1L]],
      row$sampled_stress_mean[[1L]],
      row$sampled_nonedge_sep_ratio_mean[[1L]],
      row$elapsed_sec_mean[[1L]],
      row$n_fail_graphs[[1L]]
    ))
  }

  lines <- c(lines, "", "## Family winners", "")
  for (i in seq_len(nrow(family_winners))) {
    row <- family_winners[i, , drop = FALSE]
    lines <- c(lines, sprintf(
      "- `%s`: %s (mean score %.4f, mean elapsed %.3f sec)",
      row$family[[1L]],
      row$candidate_label[[1L]],
      row$score_graph_mean[[1L]],
      row$elapsed_sec_mean[[1L]]
    ))
  }

  lines <- c(
    lines,
    "",
    "## Family summary",
    "",
    "| Family | Candidate | Mean graph score | RMSE | Edge CV | Stress | Non-edge sep | Mean elapsed sec | Failed graphs |",
    "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"
  )
  for (i in seq_len(nrow(family_summary))) {
    row <- family_summary[i, , drop = FALSE]
    lines <- c(lines, sprintf(
      "| %s | %s | %.4f | %.4f | %.4f | %.4f | %.4f | %.3f | %d |",
      row$family[[1L]],
      row$candidate_label[[1L]],
      row$score_graph_mean[[1L]],
      row$procrustes_rmse_mean[[1L]],
      row$edge_length_cv_mean[[1L]],
      row$sampled_stress_mean[[1L]],
      row$sampled_nonedge_sep_ratio_mean[[1L]],
      row$elapsed_sec_mean[[1L]],
      row$n_fail_graphs[[1L]]
    ))
  }

  lines <- c(lines, "", "## Representative PDFs", "")
  if (length(representative_pdfs) == 0L) {
    lines <- c(lines, "- none generated")
  } else {
    for (pdf in representative_pdfs) {
      lines <- c(lines, sprintf("- `%s`", pdf))
    }
  }

  writeLines(lines, con = path)
}

if (sys.nframe() == 0L) {
  args <- parse_named_args(commandArgs(trailingOnly = TRUE))
  run_tag <- validate_run_tag(if (!is.null(args$tag)) args$tag else "globalrep-fixed-candidate-benchmark")
  base_seeds <- if (!is.null(args$seeds)) parse_int_vector(args$seeds, "seeds") else 1:3
  base_time_limit_sec <- if (!is.null(args$time_limit_sec)) parse_int_scalar(args$time_limit_sec, "time_limit_sec") else 90L
  max_vertices <- if (!is.null(args$max_vertices)) parse_int_scalar(args$max_vertices, "max_vertices") else 100000L
  families_filter <- if (!is.null(args$families)) parse_char_vector(args$families, "families") else NULL

  specs <- build_graph_specs()
  if (!is.null(families_filter)) {
    specs <- Filter(function(x) x$family %in% families_filter, specs)
  }
  if (length(specs) == 0L) {
    stop("No graph specs selected")
  }

  out_root <- file.path("dev", "manual")
  out_tmp_dir <- file.path(out_root, "tmp", run_tag)
  out_pdf_dir <- file.path(out_root, "pdf", run_tag)
  out_preview_dir <- file.path(out_tmp_dir, "pdf-previews")
  dir.create(out_tmp_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_pdf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_preview_dir, recursive = TRUE, showWarnings = FALSE)

  raw_rows <- list()
  idx <- 0L
  for (spec in specs) {
    n_vertices <- max(spec$edges)
    seeds <- graph_seed_budget(n_vertices, base_seeds)
    if (n_vertices > max_vertices) {
      message(sprintf("Skipping %s %s (%d vertices) because it exceeds max_vertices=%d",
                      spec$family, spec$graph_label, n_vertices, max_vertices))
      reason <- sprintf("Skipped: %d vertices exceeds max_vertices=%d", n_vertices, max_vertices)
      for (candidate_id in candidate_ids) {
        cfg <- resolve_candidate_cfg(candidate_id, n_vertices)
        for (seed in seeds) {
          idx <- idx + 1L
          raw_rows[[idx]] <- make_skipped_row(spec, cfg, candidate_id, seed, reason)
        }
      }
      next
    }
    message(sprintf("Benchmarking %s %s (%d vertices) on seeds %s",
                    spec$family,
                    spec$graph_label,
                    n_vertices,
                    paste(seeds, collapse = ",")))
    for (candidate_id in candidate_ids) {
      cfg <- resolve_candidate_cfg(candidate_id, n_vertices)
      for (seed in seeds) {
        idx <- idx + 1L
        raw_rows[[idx]] <- run_one_layout_globalrep_safe(
          spec = spec,
          cfg = cfg,
          candidate_id = candidate_id,
          seed = seed,
          base_time_limit_sec = base_time_limit_sec
        )
      }
    }
  }

  raw_metrics <- do.call(rbind, raw_rows)
  graph_summary <- summarize_graph_results(raw_metrics)
  scored_graphs <- score_graph_results(graph_summary)
  family_summary <- summarize_family_results(scored_graphs)
  family_winners <- winner_by_family(family_summary)
  candidate_summary <- summarize_candidate_results(scored_graphs)
  representative_pdfs <- render_representative_figures(specs, raw_metrics, out_pdf_dir, out_preview_dir)

  raw_csv <- file.path(out_tmp_dir, "globalrep-fixed-candidate-raw-metrics.csv")
  graph_csv <- file.path(out_tmp_dir, "globalrep-fixed-candidate-graph-summary.csv")
  family_csv <- file.path(out_tmp_dir, "globalrep-fixed-candidate-family-summary.csv")
  candidate_csv <- file.path(out_tmp_dir, "globalrep-fixed-candidate-candidate-summary.csv")
  summary_md <- file.path(out_tmp_dir, "globalrep-fixed-candidate-summary.md")

  utils::write.csv(raw_metrics, raw_csv, row.names = FALSE)
  utils::write.csv(scored_graphs, graph_csv, row.names = FALSE)
  utils::write.csv(family_summary, family_csv, row.names = FALSE)
  utils::write.csv(candidate_summary, candidate_csv, row.names = FALSE)
  write_summary(
    path = summary_md,
    config = list(base_seeds = base_seeds, base_time_limit_sec = base_time_limit_sec, max_vertices = max_vertices),
    candidate_summary = candidate_summary,
    family_summary = family_summary,
    family_winners = family_winners,
    representative_pdfs = representative_pdfs
  )

  message(sprintf("Raw metrics written to %s", raw_csv))
  message(sprintf("Graph summary written to %s", graph_csv))
  message(sprintf("Family summary written to %s", family_csv))
  message(sprintf("Candidate summary written to %s", candidate_csv))
  message(sprintf("Markdown summary written to %s", summary_md))
}
