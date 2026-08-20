#!/usr/bin/env Rscript

sys.source(file.path("tools", "benchmarks", "gkk_lgkk_paper", "benchmark-globalrep-fixed-candidate.R"), envir = environment())

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

build_requested_specs <- function() {
  specs <- list()

  for (level in 3:6) {
    specs[[length(specs) + 1L]] <- build_carpet_spec(level)
    specs[[length(specs) + 1L]] <- build_triangle_spec(level)
    specs[[length(specs) + 1L]] <- build_tetrahedron_spec(level)
  }

  for (size in c(8L, 10L, 12L)) {
    specs[[length(specs) + 1L]] <- build_mesh_spec_local(size, size)
    specs[[length(specs) + 1L]] <- build_cylinder_spec_local(size, size)
  }

  for (size in list(
    list(h = 8L, w = 8L, label = "8x8"),
    list(h = 12L, w = 12L, label = "12x12"),
    list(h = 16L, w = 16L, label = "16x16")
  )) {
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

  specs[[length(specs) + 1L]] <- build_kary_tree_spec_local(2L, 4L)
  specs[[length(specs) + 1L]] <- build_kary_tree_spec_local(2L, 5L)
  specs[[length(specs) + 1L]] <- build_kary_tree_spec_local(2L, 6L)

  specs
}

compute_layout_metrics <- function(coords, spec, seed) {
  n <- max(spec$edges)
  adj <- make_adj_list(spec$edges, n)
  fit <- align_layout_to_spec(coords, spec)
  edge_stats <- edge_length_stats(coords, spec$edges)
  list(
    aligned = fit$aligned,
    rmse = fit$rmse,
    edge_cv = edge_stats$cv,
    stress = sampled_stress(coords, adj, sample_size = 3000L, rng_seed = 1000L + seed),
    sep = sampled_nonedge_separation_ratio(coords, spec$edges, sample_size = 6000L, rng_seed = 2000L + seed)
  )
}

run_layout_profile <- function(spec, cfg, seed) {
  n <- max(spec$edges)
  started <- proc.time()[["elapsed"]]
  coords <- grip.layout(
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
  metrics <- compute_layout_metrics(coords, spec, seed)
  c(list(coords = coords, elapsed_sec = elapsed), metrics)
}

write_summary <- function(path, seed, rows) {
  lines <- c(
    "# Globalrep diagnostic suite",
    "",
    sprintf("- seed: `%d`", seed),
    "- figures: canonical vs current globalrep defaults vs fixed candidate",
    "- note: the package defaults now use the tuned quality-first profile and only taper `final_rounds` on larger graphs",
    "",
    "| Family | Graph | Dim | Default final_rounds | Fixed final_rounds | Default RMSE | Fixed RMSE | Default stress | Fixed stress | Default elapsed sec | Fixed elapsed sec | PNG |",
    "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |"
  )
  for (i in seq_len(nrow(rows))) {
    row <- rows[i, , drop = FALSE]
    lines <- c(lines, sprintf(
      "| %s | %s | %d | %d | %d | %.4f | %.4f | %.4f | %.4f | %.3f | %.3f | `%s` |",
      row$family[[1L]],
      row$graph_label[[1L]],
      row$dim[[1L]],
      row$default_final_rounds[[1L]],
      row$fixed_final_rounds[[1L]],
      row$default_rmse[[1L]],
      row$fixed_rmse[[1L]],
      row$default_stress[[1L]],
      row$fixed_stress[[1L]],
      row$default_elapsed_sec[[1L]],
      row$fixed_elapsed_sec[[1L]],
      row$png_path[[1L]]
    ))
  }
  writeLines(lines, con = path)
}

render_one_diagnostic <- function(spec, seed, pdf_root, preview_root) {
  n <- max(spec$edges)
  default_cfg <- resolve_candidate_cfg("globalrep_default", n)
  fixed_cfg <- resolve_candidate_cfg("globalrep_fixed_candidate", n)
  default_res <- run_layout_profile(spec, default_cfg, seed)
  fixed_res <- run_layout_profile(spec, fixed_cfg, seed)

  pdf_dir <- file.path(pdf_root, spec$family)
  preview_dir <- file.path(preview_root, spec$family)
  dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)

  base_name <- sprintf("%s-default-vs-fixed-candidate", spec$graph_label)
  pdf_path <- file.path(pdf_dir, paste0(base_name, ".pdf"))
  png_path <- file.path(preview_dir, paste0(base_name, ".png"))
  subtitle <- sprintf(
    "seed=%d | default final_rounds=%d RMSE=%.4f | fixed final_rounds=%d RMSE=%.4f",
    seed,
    default_cfg$final_rounds,
    default_res$rmse,
    fixed_cfg$final_rounds,
    fixed_res$rmse
  )

  if (spec$dim == 2L) {
    plot_2d_triptych(
      path = pdf_path,
      canonical_coords = spec$canonical,
      baseline_coords = default_res$aligned,
      tuned_coords = fixed_res$aligned,
      edges = spec$edges,
      title_text = spec$title,
      subtitle_text = subtitle,
      tuned_label = "fixed candidate"
    )
  } else {
    write_3d_diagnostic_pdf(
      path = pdf_path,
      canonical_coords = spec$canonical,
      baseline_coords = default_res$aligned,
      tuned_coords = fixed_res$aligned,
      edges = spec$edges,
      title_text = spec$title,
      subtitle_text = subtitle,
      tuned_label = "fixed candidate"
    )
  }
  render_pdf_previews(pdf_path, preview_dir)

  data.frame(
    family = spec$family,
    graph_label = spec$graph_label,
    title = spec$title,
    dim = spec$dim,
    seed = seed,
    vertices = n,
    edges = nrow(spec$edges),
    default_final_rounds = default_cfg$final_rounds,
    fixed_final_rounds = fixed_cfg$final_rounds,
    default_rmse = default_res$rmse,
    fixed_rmse = fixed_res$rmse,
    default_edge_cv = default_res$edge_cv,
    fixed_edge_cv = fixed_res$edge_cv,
    default_stress = default_res$stress,
    fixed_stress = fixed_res$stress,
    default_sep = default_res$sep,
    fixed_sep = fixed_res$sep,
    default_elapsed_sec = default_res$elapsed_sec,
    fixed_elapsed_sec = fixed_res$elapsed_sec,
    pdf_path = pdf_path,
    png_path = png_path,
    stringsAsFactors = FALSE
  )
}

if (sys.nframe() == 0L) {
  args <- parse_named_args(commandArgs(trailingOnly = TRUE))
  run_tag <- validate_run_tag(if (!is.null(args$tag)) args$tag else "globalrep-diagnostic-suite")
  seed <- if (!is.null(args$seed)) parse_int_scalar(args$seed, "seed") else 1L
  families_filter <- if (!is.null(args$families)) parse_char_vector(args$families, "families") else NULL
  graph_filter <- if (!is.null(args$graph_labels)) parse_char_vector(args$graph_labels, "graph_labels") else NULL
  exclude_graph_filter <- if (!is.null(args$exclude_graph_labels)) parse_char_vector(args$exclude_graph_labels, "exclude_graph_labels") else NULL

  specs <- build_requested_specs()
  if (!is.null(families_filter)) {
    specs <- Filter(function(x) x$family %in% families_filter, specs)
  }
  if (!is.null(graph_filter)) {
    specs <- Filter(function(x) x$graph_label %in% graph_filter, specs)
  }
  if (!is.null(exclude_graph_filter)) {
    specs <- Filter(function(x) !(x$graph_label %in% exclude_graph_filter), specs)
  }
  if (length(specs) == 0L) {
    stop("No graph specs selected")
  }

  out_root <- file.path("output", "gkk_lgkk_paper")
  out_tmp_dir <- file.path(out_root, "tmp", run_tag)
  out_pdf_dir <- file.path(out_root, "reports", run_tag)
  out_preview_dir <- file.path(out_tmp_dir, "pdf-previews")
  dir.create(out_tmp_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_pdf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_preview_dir, recursive = TRUE, showWarnings = FALSE)

  rows <- list()
  for (i in seq_along(specs)) {
    spec <- specs[[i]]
    message(sprintf("Rendering %s %s (%d vertices)", spec$family, spec$graph_label, max(spec$edges)))
    rows[[i]] <- render_one_diagnostic(spec, seed, out_pdf_dir, out_preview_dir)
  }
  summary_df <- do.call(rbind, rows)

  csv_path <- file.path(out_tmp_dir, "globalrep-diagnostic-suite-summary.csv")
  md_path <- file.path(out_tmp_dir, "globalrep-diagnostic-suite-summary.md")
  utils::write.csv(summary_df, csv_path, row.names = FALSE)
  write_summary(md_path, seed, summary_df)

  message(sprintf("Summary CSV written to %s", csv_path))
  message(sprintf("Summary markdown written to %s", md_path))
}
