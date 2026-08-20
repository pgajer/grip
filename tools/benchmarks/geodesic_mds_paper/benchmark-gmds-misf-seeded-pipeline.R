#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) == 0L) {
  stop("Unable to determine script path from commandArgs().")
}
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

args <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% args

run_tag <- if (smoke) {
  sprintf("gmds-misf-seeded-pipeline-smoke-%s", format(Sys.Date(), "%Y-%m-%d"))
} else {
  "gmds-misf-seeded-pipeline-2026-04-02"
}

design_root <- file.path(repo_root, "output", "geodesic_mds_paper")
tmp_dir <- file.path(design_root, "tmp", run_tag)
pdf_dir <- file.path(design_root, "reports", run_tag)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

tex_path <- file.path(design_root, "reports", "gmds_misf_seeded_pipeline_report_2026-04-02.tex")
pdf_path <- file.path(design_root, "reports", "gmds_misf_seeded_pipeline_report_2026-04-02.pdf")
rds_path <- file.path(tmp_dir, "gmds_misf_seeded_pipeline_results.rds")
metrics_csv <- file.path(tmp_dir, "gmds_misf_seeded_pipeline_metrics.csv")
stage_trace_csv <- file.path(tmp_dir, "gmds_misf_seeded_pipeline_stage_traces.csv")
phase_b_rds <- file.path(
  design_root,
  "tmp",
  "gmds-misf-top-level-corrections-2026-04-02",
  "gmds_misf_top_level_corrections_results.rds"
)
phase_a_rds <- file.path(
  design_root,
  "tmp",
  "gmds-misf-top-level-initializers-2026-04-02",
  "gmds_misf_top_level_initializer_results.rds"
)
phase_a_script <- file.path(repo_root, "tools", "benchmark-gmds-misf-top-level-initializers.R")
phase_b_script <- file.path(repo_root, "tools", "benchmark-gmds-misf-top-level-corrections.R")

ensure_phase_b_bundle <- function() {
  if (file.exists(phase_b_rds)) {
    return(invisible(phase_b_rds))
  }
  message("Phase B bundle missing; bootstrapping prerequisite artifacts.")
  if (!file.exists(phase_a_rds)) {
    if (!file.exists(phase_a_script)) {
      stop("Phase A bundle is missing and bootstrap script not found: ", phase_a_script)
    }
    status <- system2("Rscript", phase_a_script)
    if (!identical(status, 0L) || !file.exists(phase_a_rds)) {
      stop("Failed to regenerate the Phase A bundle needed for Phase C.")
    }
  }
  if (!file.exists(phase_b_script)) {
    stop("Phase B bundle is missing and bootstrap script not found: ", phase_b_script)
  }
  status <- system2("Rscript", phase_b_script)
  if (!identical(status, 0L) || !file.exists(phase_b_rds)) {
    stop("Failed to regenerate the Phase B bundle needed for Phase C.")
  }
  invisible(phase_b_rds)
}

ensure_phase_b_bundle()

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo_root, quiet = TRUE, helpers = FALSE)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(repo_root, quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'pkgload', 'devtools', or the 'grip' package to run this benchmark.")
}

ns <- asNamespace("grip")
align_to_target_nd <- get("grip.align.to.target.nd", envir = ns)
classical_mds_embedding <- get("grip.classical.mds.embedding", envir = ns)
trace_stage_payloads <- get("grip.geodesic.misf.trace.stage.payloads", envir = ns)

cfg <- list(
  phase_seed = 20260402L,
  amplitude = 0.35,
  dim = 3L,
  top_level_restarts = if (smoke) 1L else 2L,
  top_level_max_iter = if (smoke) 3L else 6L,
  insertion_max_iter = if (smoke) 10L else 20L,
  refinement_local_nbrs = if (smoke) 3L else 4L,
  refinement_landmark_count = if (smoke) 2L else 2L,
  refinement_pair_mode = "sparse",
  refinement_anchor_weight = 0.05,
  refinement_anchor_weight_end = 0.01,
  refinement_continuation = "linear",
  refinement_max_iter = if (smoke) 2L else 3L,
  final_polish_max_iter = if (smoke) 2L else 4L,
  n_threads = 0L
)

fmt_num <- function(x, digits = 4L) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}

fmt_time <- function(x) {
  ifelse(
    is.finite(x),
    ifelse(x < 1, formatC(x, format = "f", digits = 3L), formatC(x, format = "f", digits = 2L)),
    "--"
  )
}

tex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x, perl = TRUE)
  x
}

grid_mesh_triangles <- function(h, w) {
  index <- matrix(seq_len(h * w), nrow = h, ncol = w, byrow = TRUE)
  triangles <- vector("list", 2L * (h - 1L) * (w - 1L))
  k <- 1L
  for (r in seq_len(h - 1L)) {
    for (c in seq_len(w - 1L)) {
      triangles[[k]] <- c(index[r, c], index[r + 1L, c], index[r, c + 1L])
      k <- k + 1L
      triangles[[k]] <- c(index[r + 1L, c], index[r + 1L, c + 1L], index[r, c + 1L])
      k <- k + 1L
    }
  }
  do.call(rbind, triangles)
}

triangle_areas <- function(coords, triangles) {
  coords <- as.matrix(coords)
  if (!nrow(triangles)) {
    return(numeric(0L))
  }
  if (ncol(coords) == 2L) {
    coords <- cbind(coords, 0)
  }
  v1 <- coords[triangles[, 2L], , drop = FALSE] - coords[triangles[, 1L], , drop = FALSE]
  v2 <- coords[triangles[, 3L], , drop = FALSE] - coords[triangles[, 1L], , drop = FALSE]
  cross <- cbind(
    v1[, 2L] * v2[, 3L] - v1[, 3L] * v2[, 2L],
    v1[, 3L] * v2[, 1L] - v1[, 1L] * v2[, 3L],
    v1[, 1L] * v2[, 2L] - v1[, 2L] * v2[, 1L]
  )
  0.5 * sqrt(rowSums(cross^2))
}

area_floor_ratio <- function(coords, triangles) {
  areas <- triangle_areas(coords, triangles)
  if (!length(areas)) {
    return(NA_real_)
  }
  med <- stats::median(areas)
  if (!is.finite(med) || med <= 0) {
    return(NA_real_)
  }
  as.double(stats::quantile(areas, probs = 0.05, names = FALSE)) / med
}

mesh_roughness <- function(coords, adj_list, edges) {
  centered <- sweep(coords, 2L, colMeans(coords), FUN = "-", check.margin = FALSE)
  median.edge <- stats::median(sqrt(rowSums(
    (centered[edges[, 1L], , drop = FALSE] - centered[edges[, 2L], , drop = FALSE])^2
  )))
  if (!is.finite(median.edge) || median.edge <= 0) {
    return(NA_real_)
  }
  residuals <- vapply(seq_len(nrow(centered)), function(i) {
    nbrs <- adj_list[[i]]
    if (length(nbrs) == 0L) {
      return(0)
    }
    delta <- centered[i, ] - colMeans(centered[nbrs, , drop = FALSE])
    sum(delta^2)
  }, numeric(1L))
  sqrt(mean(residuals)) / median.edge
}

make_regular_case <- function(side) {
  bundle <- mesh.surface.graph(
    side,
    side,
    surface = "paraboloid",
    amplitude = cfg$amplitude,
    connectivity = "orthogonal",
    normalize = "median"
  )
  misf.seed <- cfg$phase_seed + side
  prepared <- grip.prepare.misf.geodesic.mds(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average",
    dim = cfg$dim,
    top_level_mode = "skip",
    seed = misf.seed
  )
  started <- proc.time()[["elapsed"]]
  cmd <- classical_mds_embedding(prepared, dim = cfg$dim, eig = TRUE)
  cmd.elapsed <- proc.time()[["elapsed"]] - started
  list(
    id = sprintf("paraboloid_regular_%dx%d", side, side),
    top_case_id = sprintf("paraboloid_regular_%dx%d_top", side, side),
    label = sprintf("Regular paraboloid mesh %dx%d", side, side),
    family = "regular",
    side = as.integer(side),
    n = bundle$n,
    edges = bundle$edges,
    triangles = grid_mesh_triangles(side, side),
    truth = bundle$coords_surface,
    prepared = prepared,
    adj_list = prepared$adj_list,
    cmd = cmd,
    cmd_elapsed = as.double(cmd.elapsed),
    misf_seed = misf.seed
  )
}

make_irregular_rectangle_case <- function(side) {
  bundle <- irregular.rectangle.surface.graph(
    side,
    side,
    surface = "paraboloid",
    amplitude = cfg$amplitude,
    connectivity = "orthogonal",
    normalize = "median"
  )
  misf.seed <- cfg$phase_seed + side + 100L
  prepared <- grip.prepare.misf.geodesic.mds(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average",
    dim = cfg$dim,
    top_level_mode = "skip",
    seed = misf.seed
  )
  started <- proc.time()[["elapsed"]]
  cmd <- classical_mds_embedding(prepared, dim = cfg$dim, eig = TRUE)
  cmd.elapsed <- proc.time()[["elapsed"]] - started
  list(
    id = sprintf("paraboloid_irregular_rectangle_%dx%d", side, side),
    top_case_id = sprintf("paraboloid_irregular_rectangle_%dx%d_top", side, side),
    label = sprintf("Irregular rectangle paraboloid %dx%d", side, side),
    family = "irregular_rectangle",
    side = as.integer(side),
    n = bundle$n,
    edges = bundle$edges,
    triangles = grid_mesh_triangles(side, side),
    truth = bundle$coords_surface,
    prepared = prepared,
    adj_list = prepared$adj_list,
    cmd = cmd,
    cmd_elapsed = as.double(cmd.elapsed),
    misf_seed = misf.seed
  )
}

compute_metrics <- function(case,
                            spec,
                            coords,
                            pipeline_elapsed_sec = NA_real_,
                            seed_elapsed_sec = 0,
                            stage_trace = NULL,
                            note = NULL) {
  score_df <- grip.score.geodesic.mds(coords = coords, prepared = case$prepared)
  aligned <- if (identical(spec$id, "reference")) {
    list(aligned = as.matrix(coords), rmse = 0)
  } else {
    align_to_target_nd(coords, case$truth, allow.reflection = TRUE)
  }
  total.elapsed <- if (is.finite(pipeline_elapsed_sec) || is.finite(seed_elapsed_sec)) {
    sum(c(pipeline_elapsed_sec, seed_elapsed_sec), na.rm = TRUE)
  } else {
    NA_real_
  }
  data.frame(
    case_id = case$id,
    case_label = case$label,
    family = case$family,
    side = case$side,
    n = case$n,
    method_id = spec$id,
    method_label = spec$label,
    short_label = spec$short,
    seed_source = if (!is.null(spec$seed_source)) spec$seed_source else "",
    seed_kind = if (!is.null(spec$seed_kind)) spec$seed_kind else "",
    seed_elapsed_sec = as.double(seed_elapsed_sec),
    pipeline_elapsed_sec = as.double(pipeline_elapsed_sec),
    elapsed_sec = as.double(total.elapsed),
    gmds_energy = score_df$gmds.energy[[1L]],
    gmds_stress = score_df$gmds.stress[[1L]],
    gmds_raw_stress = score_df$gmds.raw_stress[[1L]],
    procrustes_rmse = aligned$rmse,
    roughness = mesh_roughness(coords, case$adj_list, case$edges),
    area_q05_ratio = area_floor_ratio(coords, case$triangles),
    stage_rows = if (is.null(stage_trace) || !nrow(stage_trace)) 0L else nrow(stage_trace),
    note = if (is.null(note)) "" else as.character(note),
    stringsAsFactors = FALSE
  )
}

align_stage_coords <- function(coords, target) {
  coords <- as.matrix(coords)
  keep <- stats::complete.cases(coords)
  if (sum(keep) < 3L) {
    return(coords)
  }
  aligned <- align_to_target_nd(coords[keep, , drop = FALSE], target[keep, , drop = FALSE], allow.reflection = TRUE)
  out <- matrix(NA_real_, nrow = nrow(coords), ncol = ncol(coords))
  out[keep, ] <- aligned$aligned
  out
}

phase_b_bundle <- readRDS(phase_b_rds)
phase_b_map <- setNames(phase_b_bundle$case_results, vapply(phase_b_bundle$case_results, function(x) x$case$id, character(1L)))

pick_phase_b_winners <- function(phase_b_case_result) {
  methods <- phase_b_case_result$methods
  pure.metrics <- lapply(methods, function(method) method$pure$metrics)
  sigma.vals <- vapply(pure.metrics, function(df) df$sigma_geo[[1L]], numeric(1L))
  rho.vals <- vapply(pure.metrics, function(df) df$rho[[1L]], numeric(1L))
  list(
    sigma = methods[[which.min(sigma.vals)]],
    rho = methods[[which.min(rho.vals)]]
  )
}

partial_coords <- function(coords, vertex_ids, n) {
  out <- matrix(NA_real_, nrow = n, ncol = ncol(coords))
  out[vertex_ids, ] <- coords
  out
}

build_injected_top_fit <- function(case, phase_b_case_result, phase_b_method_result, seed_kind) {
  source.ids <- as.integer(phase_b_case_result$case$top_vertex_ids)
  target.ids <- as.integer(case$prepared$top_level_vertices)
  idx <- match(target.ids, source.ids)
  if (anyNA(idx)) {
    stop("Phase B top-level vertices do not align with the full-case MISF top level for ", case$id)
  }
  top.coords <- as.matrix(phase_b_method_result$pure$coords[idx, , drop = FALSE])
  top.prepared <- case$prepared$top_level_prepared
  top.score <- grip.score.geodesic.mds(coords = top.coords, prepared = top.prepared)
  trace.df <- data.frame(
    iteration = 0L,
    energy = top.score$gmds.energy[[1L]],
    gmds_energy = top.score$gmds.energy[[1L]],
    anchor_energy = 0,
    edge_spring_energy = 0,
    repulsion_energy = 0,
    repulsion_pair_count = 0L,
    repulsion_active_pair_count = 0L,
    smooth_energy = 0,
    gradient_norm = NA_real_,
    step = NA_real_,
    accepted = TRUE,
    anchor_weight = 0,
    edge_spring_weight = 0,
    repulsion_weight = 0,
    smooth_weight = 0,
    stringsAsFactors = FALSE
  )
  seed.elapsed <- sum(c(
    phase_b_method_result$initializer$metrics$elapsed_sec[[1L]],
    phase_b_method_result$anchor$metrics$elapsed_sec[[1L]],
    phase_b_method_result$pure$metrics$elapsed_sec[[1L]]
  ), na.rm = TRUE)
  list(
    fit = list(
      coords = top.coords,
      trace = trace.df,
      frames = list(top.coords),
      prepared = top.prepared,
      score = top.score,
      restart_summary = data.frame(
        restart = 1L,
        seed = NA_integer_,
        initial.energy = top.score$gmds.energy[[1L]],
        initial.stress = top.score$gmds.stress[[1L]],
        final.energy = top.score$gmds.energy[[1L]],
        final.stress = top.score$gmds.stress[[1L]],
        improved = TRUE,
        trace.rows = 1L,
        stringsAsFactors = FALSE
      ),
      best_restart = 1L,
      best_restart_row = data.frame(
        restart = 1L,
        seed = NA_integer_,
        initial.energy = top.score$gmds.energy[[1L]],
        initial.stress = top.score$gmds.stress[[1L]],
        final.energy = top.score$gmds.energy[[1L]],
        final.stress = top.score$gmds.stress[[1L]],
        improved = TRUE,
        trace.rows = 1L,
        stringsAsFactors = FALSE
      ),
      vertex_ids = target.ids,
      coords_full = partial_coords(top.coords, target.ids, case$prepared$n),
      injected = TRUE,
      injected_seed_kind = seed_kind,
      injected_seed_label = phase_b_method_result$method_label,
      injected_seed_elapsed_sec = seed.elapsed
    ),
    seed_label = phase_b_method_result$method_label,
    seed_elapsed_sec = as.double(seed.elapsed)
  )
}

base_misf_args <- function() {
  list(
    dim = cfg$dim,
    top_level_restarts = cfg$top_level_restarts,
    top_level_max_iter = cfg$top_level_max_iter,
    top_level_engine = "cpp",
    insertion_anchor_policy = "prev_level_spread",
    insertion_max_iter = cfg$insertion_max_iter,
    refinement_local_nbrs = cfg$refinement_local_nbrs,
    refinement_landmark_count = cfg$refinement_landmark_count,
    refinement_pair_mode = cfg$refinement_pair_mode,
    refinement_anchor_weight = cfg$refinement_anchor_weight,
    refinement_anchor_weight_end = cfg$refinement_anchor_weight_end,
    refinement_continuation = cfg$refinement_continuation,
    refinement_max_iter = cfg$refinement_max_iter,
    refinement_engine = "cpp",
    final_polish_max_iter = cfg$final_polish_max_iter,
    final_polish_engine = "cpp",
    n_threads = cfg$n_threads,
    return_trace = TRUE,
    return_frames = FALSE
  )
}

make_method_specs <- function(case) {
  phase_b_case <- phase_b_map[[case$top_case_id]]
  if (is.null(phase_b_case)) {
    stop("No matching Phase B coarse case for ", case$top_case_id)
  }
  winners <- pick_phase_b_winners(phase_b_case)
  list(
    list(id = "reference", label = "Reference surface", short = "Ref", kind = "reference"),
    list(id = "cmd_pure_gmds", label = "cmdscale -> pure GMDS", short = "CMD->GMDS", kind = "pure_gmds"),
    list(
      id = "misf_default",
      label = "MISF-GMDS baseline",
      short = "MISF",
      kind = "misf_default"
    ),
    list(
      id = "misf_seed_sigma",
      label = sprintf("MISF seeded by %s (best coarse sigma)", winners$sigma$method_label),
      short = "MISF+sig",
      kind = "misf_seeded",
      seed_kind = "best_sigma",
      seed_method = winners$sigma,
      seed_source = winners$sigma$method_label
    ),
    list(
      id = "misf_seed_rho",
      label = sprintf("MISF seeded by %s (best coarse rho)", winners$rho$method_label),
      short = "MISF+rho",
      kind = "misf_seeded",
      seed_kind = "best_rho",
      seed_method = winners$rho,
      seed_source = winners$rho$method_label
    )
  )
}

run_method <- function(case, spec, phase_b_case_result) {
  if (identical(spec$kind, "reference")) {
    row <- compute_metrics(case = case, spec = spec, coords = case$truth)
    return(list(
      coords = case$truth,
      display_coords = case$truth,
      metrics = row,
      stage_trace = NULL,
      stage_coords = NULL,
      fit = NULL
    ))
  }

  if (identical(spec$kind, "pure_gmds")) {
    started <- proc.time()[["elapsed"]]
    fit <- grip.optimize.geodesic.mds(
      coords = case$cmd$coords,
      prepared = case$prepared,
      init = "user",
      anchor_mode = "none",
      engine = "cpp",
      max_iter = cfg$final_polish_max_iter,
      n_threads = cfg$n_threads,
      return_trace = TRUE,
      recenter = TRUE
    )
    elapsed <- proc.time()[["elapsed"]] - started
    row <- compute_metrics(
      case = case,
      spec = spec,
      coords = fit$coords,
      pipeline_elapsed_sec = elapsed,
      seed_elapsed_sec = case$cmd_elapsed,
      stage_trace = NULL,
      note = "Full-graph cMDS start followed by pure GMDS"
    )
    return(list(
      coords = fit$coords,
      display_coords = align_to_target_nd(fit$coords, case$truth, allow.reflection = TRUE)$aligned,
      metrics = row,
      stage_trace = NULL,
      stage_coords = NULL,
      fit = fit
    ))
  }

  method.args <- base_misf_args()
  prepared <- case$prepared
  seed.elapsed <- 0

  if (identical(spec$kind, "misf_seeded")) {
    injected <- build_injected_top_fit(case, phase_b_case_result, spec$seed_method, seed_kind = spec$seed_kind)
    prepared <- injected$fit$prepared
    prepared <- case$prepared
    prepared$top_level_fit <- injected$fit
    seed.elapsed <- injected$seed_elapsed_sec
  }

  fit <- do.call(
    grip.optimize.misf.geodesic.mds,
    c(list(prepared = prepared, seed = case$misf_seed), method.args)
  )

  note <- if (identical(spec$kind, "misf_seeded")) {
    sprintf("Top MISF level injected from Phase B %s winner: %s", spec$seed_kind, spec$seed_source)
  } else {
    "Current MISF-GMDS baseline with internal top-level pure-GMDS solve"
  }
  row <- compute_metrics(
    case = case,
    spec = spec,
    coords = fit$coords,
    pipeline_elapsed_sec = fit$timing$total,
    seed_elapsed_sec = seed.elapsed,
    stage_trace = fit$stage_trace,
    note = note
  )
  stage.payloads <- trace_stage_payloads(
    fit,
    target = case$truth,
    states = c("top_level", "after_insertion", "after_refinement", "final_polish")
  )
  stage.coords <- lapply(stage.payloads, function(payload) payload$coords)
  stage.display <- lapply(stage.payloads, function(payload) payload$display_coords)
  list(
    coords = fit$coords,
    display_coords = align_to_target_nd(fit$coords, case$truth, allow.reflection = TRUE)$aligned,
    metrics = row,
    stage_trace = fit$stage_trace,
    stage_coords = stage.coords,
    stage_display = stage.display,
    fit = fit
  )
}

stage_trace_rows <- function(case_result) {
  rows <- list()
  for (method in case_result$methods) {
    if (is.null(method$stage_trace) || !nrow(method$stage_trace)) {
      next
    }
    trace_df <- method$stage_trace
    trace_df <- trace_df[is.finite(trace_df$energy), , drop = FALSE]
    if (!nrow(trace_df)) {
      next
    }
    trace_df$case_id <- case_result$case$id
    trace_df$case_label <- case_result$case$label
    trace_df$family <- case_result$case$family
    trace_df$method_id <- method$metrics$method_id[[1L]]
    trace_df$method_label <- method$metrics$method_label[[1L]]
    trace_df$stage_index <- seq_len(nrow(trace_df))
    trace_df$stage_label <- ifelse(
      trace_df$stage == "refinement",
      sprintf("ref L%d", trace_df$level),
      ifelse(trace_df$stage == "top_level", "top", "final")
    )
    rows[[length(rows) + 1L]] <- trace_df
  }
  if (!length(rows)) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

save_case_panel_grid <- function(case_result, output_path) {
  methods <- case_result$methods
  grDevices::png(output_path, width = 2400L, height = 1600L, res = 180, bg = "#ffffff")
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(2L, 3L), mar = c(1.2, 1.2, 2.9, 0.4), oma = c(0, 0, 1.2, 0))

  for (i in seq_len(6L)) {
    if (i > length(methods)) {
      graphics::plot.new()
      next
    }
    method <- methods[[i]]
    grip.plot(
      coords = method$display_coords,
      edges = case_result$case$edges,
      projection = "ortho",
      azimuth = 35,
      elevation = 24,
      vertex.col = if (identical(method$metrics$method_id[[1L]], "reference")) "#bc6c25" else "#3a5a40",
      edge.col = "#adb5bd",
      main = ""
    )
    row <- method$metrics[1L, , drop = FALSE]
    ttl <- if (identical(row$method_id[[1L]], "reference")) {
      "Reference surface"
    } else {
      sprintf(
        "%s\nsigma %s, rho %s\ntotal %ss",
        row$method_label[[1L]],
        fmt_num(row$gmds_stress[[1L]], 4L),
        fmt_num(row$procrustes_rmse[[1L]], 4L),
        fmt_time(row$elapsed_sec[[1L]])
      )
    }
    graphics::mtext(ttl, side = 3L, line = 0.3, cex = 0.74)
  }
  graphics::mtext(case_result$case$label, side = 3L, outer = TRUE, line = -0.3, cex = 1.1, font = 2L)
}

save_stage_grid <- function(case_result, output_path) {
  stage.methods <- Filter(function(method) !is.null(method$stage_display), case_result$methods)
  if (!length(stage.methods)) {
    return(invisible(NULL))
  }
  stage.names <- c("top_level", "after_insertion", "after_refinement", "final_polish")
  grDevices::png(output_path, width = 3200L, height = max(1000L, 800L * length(stage.methods)), res = 180, bg = "#ffffff")
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(length(stage.methods), length(stage.names)), mar = c(1.0, 1.0, 2.7, 0.4), oma = c(0, 0, 1.2, 0))

  for (method in stage.methods) {
    for (stage.name in stage.names) {
      coords <- method$stage_display[[stage.name]]
      keep <- stats::complete.cases(coords)
      stage.edges <- case_result$case$edges
      if (!all(keep)) {
        good.edges <- keep[stage.edges[, 1L]] & keep[stage.edges[, 2L]]
        stage.edges <- stage.edges[good.edges, , drop = FALSE]
      }
      if (sum(keep) >= 2L && nrow(stage.edges) > 0L) {
        grip.plot(
          coords = coords,
          edges = stage.edges,
          projection = "ortho",
          azimuth = 35,
          elevation = 24,
          vertex.col = "#355070",
          edge.col = "#c6d1db",
          main = ""
        )
      } else if (sum(keep) >= 1L) {
        graphics::plot.new()
        graphics::points(
          coords[keep, 1L],
          coords[keep, 2L],
          pch = 16,
          col = "#355070"
        )
      } else {
        graphics::plot.new()
      }
      ttl <- sprintf(
        "%s\n%s",
        method$metrics$short_label[[1L]],
        switch(stage.name,
          top_level = "top level",
          after_insertion = "after insertion",
          after_refinement = "after refinement",
          final_polish = "final polish"
        )
      )
      graphics::mtext(ttl, side = 3L, line = 0.3, cex = 0.75)
    }
  }
  graphics::mtext(sprintf("%s: MISF stage layouts", case_result$case$label), side = 3L, outer = TRUE, line = -0.3, cex = 1.05, font = 2L)
}

save_stage_trace_grid <- function(stage_df, cases, output_path) {
  case_ids <- vapply(cases, `[[`, character(1L), "id")
  grDevices::png(output_path, width = 2400L, height = 1500L, res = 180, bg = "#ffffff")
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(2L, 2L), mar = c(3.5, 3.7, 2.2, 0.8), oma = c(0, 0, 1.0, 0))

  cols <- c("misf_default" = "#3a5a40", "misf_seed_sigma" = "#9c6644", "misf_seed_rho" = "#355070")
  for (case_id in case_ids) {
    df <- stage_df[stage_df$case_id == case_id, , drop = FALSE]
    if (!nrow(df)) {
      graphics::plot.new()
      next
    }
    yvals <- log10(pmax(df$energy, 1e-12))
    xlim <- c(1, max(df$stage_index))
    ylim <- range(yvals)
    graphics::plot(NA,
      xlim = xlim,
      ylim = ylim,
      xlab = "Stage step",
      ylab = "log10 energy",
      main = unique(df$case_label),
      xaxt = "n"
    )
    for (method_id in names(cols)) {
      sub <- df[df$method_id == method_id, , drop = FALSE]
      if (!nrow(sub)) {
        next
      }
      graphics::lines(
        x = sub$stage_index,
        y = log10(pmax(sub$energy, 1e-12)),
        type = "b",
        pch = 16,
        lwd = 2,
        col = cols[[method_id]]
      )
      graphics::text(
        x = sub$stage_index,
        y = log10(pmax(sub$energy, 1e-12)),
        labels = sub$stage_label,
        pos = 3,
        cex = 0.63,
        col = cols[[method_id]]
      )
    }
    graphics::axis(1, at = sort(unique(df$stage_index)))
    graphics::legend(
      "topright",
      legend = c("MISF baseline", "MISF seeded (sigma)", "MISF seeded (rho)"),
      col = cols,
      lty = 1,
      lwd = 2,
      pch = 16,
      bty = "n",
      cex = 0.8
    )
  }
  graphics::mtext("Phase C per-stage MISF energy traces", side = 3L, outer = TRUE, line = -0.3, cex = 1.1, font = 2L)
}

write_case_table <- function(df) {
  rows <- vapply(seq_len(nrow(df)), function(i) {
    sprintf(
      "%s & %s & %s & %s & %s & %s & %s & %s \\\\",
      tex_escape(df$method_label[[i]]),
      fmt_time(df$seed_elapsed_sec[[i]]),
      fmt_time(df$pipeline_elapsed_sec[[i]]),
      fmt_time(df$elapsed_sec[[i]]),
      fmt_num(df$gmds_stress[[i]], 4L),
      fmt_num(df$procrustes_rmse[[i]], 4L),
      fmt_num(df$roughness[[i]], 4L),
      fmt_num(df$area_q05_ratio[[i]], 4L)
    )
  }, character(1L))
  paste(
    "\\begin{table}[H]",
    "\\centering",
    "\\scriptsize",
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{lrrrrrrr}",
    "\\toprule",
    "Method & $t_{seed}$ (s) & $t_{pipe}$ (s) & $t$ (s) & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ \\\\",
    "\\midrule",
    paste(rows, collapse = "\n"),
    "\\bottomrule",
    "\\end{tabular}}",
    "\\end{table}",
    sep = "\n"
  )
}

case_summary_paragraph <- function(df) {
  baseline <- df[df$method_id == "misf_default", , drop = FALSE]
  seeded.sigma <- df[df$method_id == "misf_seed_sigma", , drop = FALSE]
  seeded.rho <- df[df$method_id == "misf_seed_rho", , drop = FALSE]
  cmd <- df[df$method_id == "cmd_pure_gmds", , drop = FALSE]
  sprintf(
    "On %s, the direct cmdscale-to-GMDS control finishes at $\\sigma=%s$ and $\\rho=%s$. The current MISF-GMDS baseline lands at $\\sigma=%s$ and $\\rho=%s$. Seeding the full MISF pipeline with the best Phase~B coarse $\\sigma$ winner (%s) yields $\\sigma=%s$ and $\\rho=%s$, while seeding with the best Phase~B coarse $\\rho$ winner (%s) yields $\\sigma=%s$ and $\\rho=%s$.",
    tex_escape(df$case_label[[1L]]),
    fmt_num(cmd$gmds_stress[[1L]], 4L),
    fmt_num(cmd$procrustes_rmse[[1L]], 4L),
    fmt_num(baseline$gmds_stress[[1L]], 4L),
    fmt_num(baseline$procrustes_rmse[[1L]], 4L),
    tex_escape(seeded.sigma$seed_source[[1L]]),
    fmt_num(seeded.sigma$gmds_stress[[1L]], 4L),
    fmt_num(seeded.sigma$procrustes_rmse[[1L]], 4L),
    tex_escape(seeded.rho$seed_source[[1L]]),
    fmt_num(seeded.rho$gmds_stress[[1L]], 4L),
    fmt_num(seeded.rho$procrustes_rmse[[1L]], 4L)
  )
}

case_figure_rel <- function(case_id) {
  file.path(pdf_dir, sprintf("%s_final_grid.png", case_id))
}

stage_figure_rel <- function(case_id) {
  file.path(pdf_dir, sprintf("%s_stage_grid.png", case_id))
}

stage_trace_figure_rel <- file.path(pdf_dir, "gmds_misf_seeded_pipeline_stage_trace_grid.png")

cases <- if (smoke) {
  list(
    make_regular_case(12L),
    make_irregular_rectangle_case(15L)
  )
} else {
  list(
    make_regular_case(12L),
    make_regular_case(15L),
    make_irregular_rectangle_case(15L)
  )
}

case_results <- lapply(cases, function(case) {
  message("Running Phase C full seeded MISF panel for: ", case$label)
  phase_b_case_result <- phase_b_map[[case$top_case_id]]
  if (is.null(phase_b_case_result)) {
    stop("Missing Phase B result for top case ", case$top_case_id)
  }
  if (!setequal(case$prepared$top_level_vertices, phase_b_case_result$case$top_vertex_ids)) {
    stop("Top-level vertex sets differ between full case and Phase B result for ", case$id)
  }
  method_specs <- make_method_specs(case)
  methods <- lapply(method_specs, function(spec) run_method(case, spec, phase_b_case_result = phase_b_case_result))
  metrics <- do.call(rbind, lapply(methods, `[[`, "metrics"))
  out <- list(case = case, methods = methods, metrics = metrics)
  save_case_panel_grid(out, case_figure_rel(case$id))
  save_stage_grid(out, stage_figure_rel(case$id))
  out
})

compact_case_result <- function(case_result) {
  list(
    case = list(
      id = case_result$case$id,
      label = case_result$case$label,
      family = case_result$case$family,
      side = case_result$case$side,
      n = case_result$case$n,
      edges = case_result$case$edges
    ),
    methods = lapply(case_result$methods, function(method) {
      list(
        display_coords = method$display_coords,
        metrics = method$metrics,
        stage_display = method$stage_display
      )
    }),
    metrics = case_result$metrics
  )
}

bundle_case_results <- lapply(case_results, compact_case_result)

metrics_df <- do.call(rbind, lapply(case_results, `[[`, "metrics"))
stage_trace_df <- do.call(
  rbind,
  Filter(function(x) nrow(x) > 0L, lapply(case_results, stage_trace_rows))
)

if (nrow(stage_trace_df) > 0L) {
  save_stage_trace_grid(stage_trace_df, cases = lapply(case_results, `[[`, "case"), output_path = stage_trace_figure_rel)
}

utils::write.csv(metrics_df, metrics_csv, row.names = FALSE)
utils::write.csv(stage_trace_df, stage_trace_csv, row.names = FALSE)
saveRDS(
  list(
    run_tag = run_tag,
    cfg = cfg,
    case_results = bundle_case_results,
    metrics = metrics_df,
    stage_traces = stage_trace_df
  ),
  rds_path
)

seeded_summary <- subset(metrics_df, method_id %in% c("misf_default", "misf_seed_sigma", "misf_seed_rho"))
best_seeded_by_case <- do.call(rbind, lapply(split(seeded_summary, seeded_summary$case_id), function(df) {
  seeded <- df[df$method_id %in% c("misf_seed_sigma", "misf_seed_rho"), , drop = FALSE]
  baseline <- df[df$method_id == "misf_default", , drop = FALSE]
  best.rho <- seeded$method_label[[which.min(seeded$procrustes_rmse)]]
  best.sigma <- seeded$method_label[[which.min(seeded$gmds_stress)]]
  data.frame(
    case_id = baseline$case_id[[1L]],
    case_label = baseline$case_label[[1L]],
    baseline_sigma = baseline$gmds_stress[[1L]],
    baseline_rho = baseline$procrustes_rmse[[1L]],
    best_seeded_sigma_label = best.sigma,
    best_seeded_sigma = min(seeded$gmds_stress),
    best_seeded_rho_label = best.rho,
    best_seeded_rho = min(seeded$procrustes_rmse),
    stringsAsFactors = FALSE
  )
}))

seeded_better_rho <- sum(best_seeded_by_case$best_seeded_rho < best_seeded_by_case$baseline_rho)
seeded_better_sigma <- sum(best_seeded_by_case$best_seeded_sigma < best_seeded_by_case$baseline_sigma)

case_sections <- vapply(case_results, function(case_result) {
  df <- subset(case_result$metrics, method_id != "reference")
  paste(
    sprintf("\\section*{%s}", tex_escape(case_result$case$label)),
    case_summary_paragraph(df),
    "",
    sprintf(
      "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{%s. Final full-graph layouts for the direct cmdscale control, the current MISF baseline, and the two Phase~C seeded MISF pipelines.}\\end{figure}",
      case_figure_rel(case_result$case$id),
      tex_escape(case_result$case$label)
    ),
    sprintf(
      "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{%s. Stagewise geometry for the MISF pipelines: top-level scaffold, after insertion, after sparse level refinement, and final layout.}\\end{figure}",
      stage_figure_rel(case_result$case$id),
      tex_escape(case_result$case$label)
    ),
    write_case_table(df),
    sep = "\n\n"
  )
}, character(1L))

tex_lines <- c(
  "\\documentclass[11pt]{article}",
  "\\usepackage[margin=1in]{geometry}",
  "\\usepackage{graphicx}",
  "\\usepackage{booktabs}",
  "\\usepackage{float}",
  "\\usepackage{amsmath}",
  "\\usepackage[T1]{fontenc}",
  "\\usepackage[utf8]{inputenc}",
  "\\title{Phase C: Seeded MISF-GMDS Full-Pipeline Benchmark}",
  "\\author{MISF-GMDS top-level initializer plan}",
  "\\date{2026-04-02}",
  "\\begin{document}",
  "\\maketitle",
  "\\section{Purpose}",
  "This report implements Group~C of the MISF-GMDS top-level initializer plan. The question is no longer which coarse top-level initializer looks best in isolation, but whether reusing the best Phase~B coarse corrected top-level embeddings improves the \\emph{entire} multiscale MISF-GMDS pipeline once insertion, sparse refinement, and final polish are all run on the full graph.",
  "\\section{Experimental design}",
  "For each family, we compare four optimization paths: a direct full-graph cmdscale-to-GMDS control, the current MISF-GMDS baseline with its internal top-level pure-GMDS solve, a seeded MISF-GMDS run initialized from the best Phase~B coarse $\\sigma$ winner, and a seeded MISF-GMDS run initialized from the best Phase~B coarse $\\rho$ winner. The seeded pipelines keep the same insertion, refinement, and final-polish settings as the baseline; only the top MISF level is replaced.",
  sprintf(
    "Across the %d current cases, the best seeded variant improves over the current MISF baseline in $\\rho$ on %d cases and in $\\sigma$ on %d cases.",
    nrow(best_seeded_by_case),
    seeded_better_rho,
    seeded_better_sigma
  ),
  paste(case_sections, collapse = "\n\n"),
  "\\section*{Per-stage energy ladders}",
  "Figure below condenses the MISF stage traces for the baseline and seeded multiscale runs. It shows the coarse top-level stage, the per-level refinement endpoints, and the final polish, making it easier to see whether the seeded pipelines descend into a better full-graph basin than the baseline multiscale solve.",
  sprintf(
    "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{Phase~C per-stage MISF energy traces.}\\end{figure}",
    stage_trace_figure_rel
  ),
  "\\section*{Interactive companion}",
  paste(
    "The companion HTML gallery is generated by",
    "\\texttt{tools/reports/geodesic_mds_paper/render-gmds-misf-seeded-pipeline-html.R}.",
    "It shows all final layouts and all saved stage layouts for the MISF pipelines as interactive \\texttt{rglwidget} panels."
  ),
  "\\end{document}"
)

writeLines(tex_lines, tex_path)

message("Wrote Phase C metrics: ", metrics_csv)
message("Wrote Phase C stage traces: ", stage_trace_csv)
message("Wrote Phase C bundle: ", rds_path)
message("Wrote Phase C report TeX: ", tex_path)
