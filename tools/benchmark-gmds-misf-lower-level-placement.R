#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) == 0L) {
  stop("Unable to determine script path from commandArgs().")
}
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

args <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% args

run_tag <- if (smoke) {
  sprintf("gmds-misf-lower-level-placement-smoke-%s", format(Sys.Date(), "%Y-%m-%d"))
} else {
  "gmds-misf-lower-level-placement-2026-04-02"
}

design_root <- file.path(repo_root, "dev", "design")
tmp_dir <- file.path(design_root, "tmp", run_tag)
pdf_dir <- file.path(design_root, "pdf", run_tag)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

tex_path <- file.path(design_root, "pdf", "gmds_misf_lower_level_placement_report_2026-04-02.tex")
pdf_path <- file.path(design_root, "pdf", "gmds_misf_lower_level_placement_report_2026-04-02.pdf")
rds_path <- file.path(tmp_dir, "gmds_misf_lower_level_placement_results.rds")
metrics_csv <- file.path(tmp_dir, "gmds_misf_lower_level_placement_metrics.csv")
stage_csv <- file.path(tmp_dir, "gmds_misf_lower_level_placement_stage_traces.csv")

phase_a_rds <- file.path(
  design_root,
  "tmp",
  "gmds-misf-top-level-initializers-2026-04-02",
  "gmds_misf_top_level_initializer_results.rds"
)
phase_b_rds <- file.path(
  design_root,
  "tmp",
  "gmds-misf-top-level-corrections-2026-04-02",
  "gmds_misf_top_level_corrections_results.rds"
)
phase_c_rds <- file.path(
  design_root,
  "tmp",
  "gmds-misf-seeded-pipeline-2026-04-02",
  "gmds_misf_seeded_pipeline_results.rds"
)
phase_a_script <- file.path(repo_root, "tools", "benchmark-gmds-misf-top-level-initializers.R")
phase_b_script <- file.path(repo_root, "tools", "benchmark-gmds-misf-top-level-corrections.R")
phase_c_script <- file.path(repo_root, "tools", "benchmark-gmds-misf-seeded-pipeline.R")

ensure_phase_c_bundle <- function() {
  if (file.exists(phase_c_rds)) {
    return(invisible(phase_c_rds))
  }
  message("Phase C bundle missing; bootstrapping prerequisite artifacts.")
  if (!file.exists(phase_a_rds)) {
    status <- system2("Rscript", phase_a_script)
    if (!identical(status, 0L) || !file.exists(phase_a_rds)) {
      stop("Failed to regenerate the Phase A bundle needed for Phase D.")
    }
  }
  if (!file.exists(phase_b_rds)) {
    status <- system2("Rscript", phase_b_script)
    if (!identical(status, 0L) || !file.exists(phase_b_rds)) {
      stop("Failed to regenerate the Phase B bundle needed for Phase D.")
    }
  }
  status <- system2("Rscript", phase_c_script)
  if (!identical(status, 0L) || !file.exists(phase_c_rds)) {
    stop("Failed to regenerate the Phase C bundle needed for Phase D.")
  }
  invisible(phase_c_rds)
}

ensure_phase_c_bundle()

for (pkg in c("igraph")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required to run this benchmark.", pkg))
  }
}

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
place_level_with_layout <- get("grip.geodesic.misf.place.level.with.layout", envir = ns)
active_level_vertices <- get("grip.geodesic.misf.active.level.vertices", envir = ns)
complete_edge_matrix <- get("grip.geodesic.misf.complete.edge.matrix", envir = ns)

cfg <- list(
  phase_seed = 20260402L,
  amplitude = 0.35,
  dim = 3L,
  insertion_max_iter = if (smoke) 10L else 20L,
  refinement_local_nbrs = if (smoke) 3L else 4L,
  refinement_landmark_count = if (smoke) 2L else 2L,
  refinement_pair_mode = "sparse",
  refinement_anchor_weight = 0.05,
  refinement_anchor_weight_end = 0.01,
  refinement_continuation = "linear",
  refinement_max_iter = if (smoke) 2L else 3L,
  final_polish_max_iter = if (smoke) 2L else 4L,
  n_threads = 0L,
  layout_k = 6L,
  fr_niter = if (smoke) 160L else 480L
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

resolve_unweighted_grip_args <- function(family) {
  if (identical(family, "irregular_rectangle")) {
    return(list(
      placement = "barycenter",
      rounds = 96L,
      final_rounds = 128L,
      num_init = 18L,
      num_nbrs = 24L,
      r = 0.05,
      s = 6.5,
      repulsion_factor = 1.10
    ))
  }
  list(
    placement = "barycenter",
    rounds = 72L,
    final_rounds = 96L,
    num_init = 12L,
    num_nbrs = 18L,
    r = 0.08,
    s = 5.0,
    repulsion_factor = 1.25
  )
}

resolve_weighted_grip_args <- function(top_n) {
  list(
    rounds = 96L,
    final_rounds = 128L,
    num_init = min(18L, max(4L, top_n - 1L)),
    num_nbrs = min(20L, max(4L, top_n - 1L))
  )
}

align_partial_to_truth <- function(coords, target) {
  coords <- as.matrix(coords)
  keep <- stats::complete.cases(coords)
  out <- matrix(NA_real_, nrow = nrow(coords), ncol = ncol(coords))
  if (sum(keep) >= 3L) {
    aligned <- align_to_target_nd(coords[keep, , drop = FALSE], target[keep, , drop = FALSE], allow.reflection = TRUE)
    out[keep, ] <- aligned$aligned
    return(list(aligned = out, rmse = aligned$rmse))
  }
  out[keep, ] <- coords[keep, , drop = FALSE]
  list(aligned = out, rmse = NA_real_)
}

build_active_stage_prepared <- function(prepared, level) {
  active.vertices <- active_level_vertices(prepared, level)
  active.distance <- prepared$distance_matrix[active.vertices, active.vertices, drop = FALSE]
  active.edges <- complete_edge_matrix(length(active.vertices))
  active.weights <- if (nrow(active.edges) == 0L) {
    numeric(0L)
  } else {
    as.double(active.distance[cbind(active.edges[, 1L], active.edges[, 2L])])
  }
  list(
    active_vertices = active.vertices,
    prepared = grip.prepare.graph.geodesic.mds(
      edges = active.edges,
      n = length(active.vertices),
      edge_weights = active.weights,
      tie_mode = prepared$tie_mode
    )
  )
}

make_case_regular <- function(side) {
  bundle <- mesh.surface.graph(
    side,
    side,
    surface = "paraboloid",
    amplitude = cfg$amplitude,
    connectivity = "orthogonal",
    normalize = "median"
  )
  prepared <- grip.prepare.misf.geodesic.mds(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average",
    dim = cfg$dim,
    top_level_mode = "skip",
    seed = cfg$phase_seed + side
  )
  level.ids <- seq.int(from = prepared$top_level_level, to = 0L, by = -1L)
  active.prepared <- setNames(lapply(level.ids, function(level) {
    build_active_stage_prepared(prepared, level)
  }), as.character(level.ids))
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
    weighted_preset = "mesh",
    unweighted_grip_args = resolve_unweighted_grip_args("regular"),
    weighted_grip_args = resolve_weighted_grip_args(bundle$n),
    active_prepared = active.prepared,
    misf_seed = cfg$phase_seed + side
  )
}

make_case_irregular_rectangle <- function(side) {
  bundle <- irregular.rectangle.surface.graph(
    side,
    side,
    surface = "paraboloid",
    amplitude = cfg$amplitude,
    connectivity = "orthogonal",
    normalize = "median"
  )
  prepared <- grip.prepare.misf.geodesic.mds(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average",
    dim = cfg$dim,
    top_level_mode = "skip",
    seed = cfg$phase_seed + side + 100L
  )
  level.ids <- seq.int(from = prepared$top_level_level, to = 0L, by = -1L)
  active.prepared <- setNames(lapply(level.ids, function(level) {
    build_active_stage_prepared(prepared, level)
  }), as.character(level.ids))
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
    weighted_preset = "irregular",
    unweighted_grip_args = resolve_unweighted_grip_args("irregular_rectangle"),
    weighted_grip_args = resolve_weighted_grip_args(bundle$n),
    active_prepared = active.prepared,
    misf_seed = cfg$phase_seed + side + 100L
  )
}

compute_metrics <- function(case,
                            spec,
                            coords,
                            top_seed_source,
                            top_seed_kind,
                            seed_elapsed_sec,
                            placement_elapsed_sec,
                            refinement_elapsed_sec,
                            final_elapsed_sec,
                            note = NULL) {
  score_df <- grip.score.geodesic.mds(coords = coords, prepared = case$prepared)
  aligned <- if (identical(spec$id, "reference")) {
    list(aligned = as.matrix(coords), rmse = 0)
  } else {
    align_to_target_nd(coords, case$truth, allow.reflection = TRUE)
  }
  total.elapsed <- sum(c(seed_elapsed_sec, placement_elapsed_sec, refinement_elapsed_sec, final_elapsed_sec), na.rm = TRUE)
  data.frame(
    case_id = case$id,
    case_label = case$label,
    family = case$family,
    side = case$side,
    n = case$n,
    method_id = spec$id,
    method_label = spec$label,
    short_label = spec$short,
    top_seed_source = top_seed_source,
    top_seed_kind = top_seed_kind,
    seed_elapsed_sec = as.double(seed_elapsed_sec),
    placement_elapsed_sec = as.double(placement_elapsed_sec),
    refinement_elapsed_sec = as.double(refinement_elapsed_sec),
    final_elapsed_sec = as.double(final_elapsed_sec),
    elapsed_sec = as.double(total.elapsed),
    gmds_energy = score_df$gmds.energy[[1L]],
    gmds_stress = score_df$gmds.stress[[1L]],
    gmds_raw_stress = score_df$gmds.raw_stress[[1L]],
    procrustes_rmse = aligned$rmse,
    roughness = mesh_roughness(coords, case$adj_list, case$edges),
    area_q05_ratio = area_floor_ratio(coords, case$triangles),
    note = if (is.null(note)) "" else as.character(note),
    stringsAsFactors = FALSE
  )
}

record_stage <- function(case,
                         spec,
                         coords,
                         level,
                         stage_kind,
                         stage_label,
                         elapsed_sec,
                         cumulative_elapsed,
                         rows,
                         snapshots) {
  active.info <- case$active_prepared[[as.character(level)]]
  active.vertices <- active.info$active_vertices
  active.coords <- coords[active.vertices, , drop = FALSE]
  active.score <- grip.score.geodesic.mds(coords = active.coords, prepared = active.info$prepared)
  aligned <- align_partial_to_truth(coords, case$truth)
  aligned.active <- align_to_target_nd(active.coords, case$truth[active.vertices, , drop = FALSE], allow.reflection = TRUE)

  rows[[length(rows) + 1L]] <- data.frame(
    case_id = case$id,
    case_label = case$label,
    family = case$family,
    method_id = spec$id,
    method_label = spec$label,
    stage_index = length(rows) + 1L,
    stage_kind = stage_kind,
    stage_label = stage_label,
    level = as.integer(level),
    active_n = length(active.vertices),
    energy = active.score$gmds.energy[[1L]],
    stress = active.score$gmds.stress[[1L]],
    rho = aligned.active$rmse,
    elapsed_sec = as.double(elapsed_sec),
    cumulative_elapsed = as.double(cumulative_elapsed),
    stringsAsFactors = FALSE
  )
  snapshots[[length(snapshots) + 1L]] <- list(
    stage_index = length(snapshots) + 1L,
    stage_kind = stage_kind,
    stage_label = stage_label,
    level = as.integer(level),
    active_n = length(active.vertices),
    display_coords = aligned$aligned
  )
  list(rows = rows, snapshots = snapshots)
}

phase_b_bundle <- readRDS(phase_b_rds)
phase_b_map <- setNames(phase_b_bundle$case_results, vapply(phase_b_bundle$case_results, function(x) x$case$id, character(1L)))
phase_c_bundle <- readRDS(phase_c_rds)

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

pick_phase_c_seed_kind <- function(case_id) {
  df <- phase_c_bundle$metrics
  df <- df[df$case_id == case_id & df$method_id %in% c("misf_seed_sigma", "misf_seed_rho"), , drop = FALSE]
  df <- df[order(df$procrustes_rmse, df$gmds_stress), , drop = FALSE]
  df$method_id[[1L]]
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
  seed.elapsed <- sum(c(
    phase_b_method_result$initializer$metrics$elapsed_sec[[1L]],
    phase_b_method_result$anchor$metrics$elapsed_sec[[1L]],
    phase_b_method_result$pure$metrics$elapsed_sec[[1L]]
  ), na.rm = TRUE)
  list(
    fit = list(
      coords = top.coords,
      trace = data.frame(),
      frames = list(top.coords),
      prepared = top.prepared,
      score = top.score,
      restart_summary = data.frame(),
      best_restart = 1L,
      best_restart_row = data.frame(),
      vertex_ids = target.ids,
      coords_full = partial_coords(top.coords, target.ids, case$prepared$n),
      injected = TRUE,
      injected_seed_kind = seed_kind,
      injected_seed_label = phase_b_method_result$method_label,
      injected_seed_elapsed_sec = seed.elapsed
    ),
    seed_label = phase_b_method_result$method_label,
    seed_elapsed_sec = as.double(seed.elapsed),
    seed_kind = seed_kind
  )
}

method_specs <- list(
  list(id = "geodesic_insert", label = "Geodesic insertion", short = "Geo"),
  list(id = "weighted_kk_insert", label = "Weighted KK seed", short = "W-KK"),
  list(id = "kk_insert", label = "KK seed", short = "KK"),
  list(id = "fr_insert", label = "FR seed", short = "FR"),
  list(id = "grip_insert", label = "GRIP seed", short = "GRIP"),
  list(id = "weighted_grip_insert", label = "Weighted GRIP seed", short = "W-GRIP")
)

run_method <- function(case, spec, top_seed) {
  prepared <- case$prepared
  prepared$top_level_fit <- top_seed$fit
  coords <- prepared$top_level_fit$coords_full
  placement.elapsed.total <- 0
  refinement.elapsed.total <- 0
  final.elapsed <- 0
  stage.rows <- list()
  stage.snapshots <- list()
  cumulative <- top_seed$seed_elapsed_sec

  rec <- record_stage(
    case = case,
    spec = spec,
    coords = coords,
    level = prepared$top_level_level,
    stage_kind = "top_level",
    stage_label = sprintf("top L%d", prepared$top_level_level),
    elapsed_sec = 0,
    cumulative_elapsed = cumulative,
    rows = stage.rows,
    snapshots = stage.snapshots
  )
  stage.rows <- rec$rows
  stage.snapshots <- rec$snapshots

  level.ids <- seq.int(from = prepared$top_level_level - 1L, to = 0L, by = -1L)
  for (level in level.ids) {
    place.seed <- cfg$phase_seed + case$side * 100L + level * 10L + match(spec$id, vapply(method_specs, `[[`, character(1L), "id"))
    placement.start <- proc.time()[["elapsed"]]
    placement <- switch(
      spec$id,
      geodesic_insert = grip.geodesic.misf.insert.level(
        prepared = prepared,
        coords = coords,
        level = level,
        anchor_policy = "prev_level_spread",
        max_iter = cfg$insertion_max_iter
      ),
      weighted_kk_insert = place_level_with_layout(
        prepared = prepared,
        coords = coords,
        level = level,
        method = "weighted_kk",
        layout_k = cfg$layout_k,
        seed = place.seed
      ),
      kk_insert = place_level_with_layout(
        prepared = prepared,
        coords = coords,
        level = level,
        method = "kk",
        layout_k = cfg$layout_k,
        seed = place.seed
      ),
      fr_insert = place_level_with_layout(
        prepared = prepared,
        coords = coords,
        level = level,
        method = "fr",
        layout_k = cfg$layout_k,
        fr_niter = cfg$fr_niter,
        seed = place.seed
      ),
      grip_insert = place_level_with_layout(
        prepared = prepared,
        coords = coords,
        level = level,
        method = "grip",
        layout_k = cfg$layout_k,
        grip_args = case$unweighted_grip_args,
        seed = place.seed
      ),
      weighted_grip_insert = place_level_with_layout(
        prepared = prepared,
        coords = coords,
        level = level,
        method = "weighted_grip",
        layout_k = cfg$layout_k,
        weighted_preset = case$weighted_preset,
        weighted_args = case$weighted_grip_args,
        seed = place.seed
      )
    )
    placement.elapsed <- proc.time()[["elapsed"]] - placement.start
    placement.elapsed.total <- placement.elapsed.total + placement.elapsed
    cumulative <- cumulative + placement.elapsed
    coords <- placement$coords
    rec <- record_stage(
      case = case,
      spec = spec,
      coords = coords,
      level = level,
      stage_kind = "placement",
      stage_label = sprintf("place L%d", level),
      elapsed_sec = placement.elapsed,
      cumulative_elapsed = cumulative,
      rows = stage.rows,
      snapshots = stage.snapshots
    )
    stage.rows <- rec$rows
    stage.snapshots <- rec$snapshots
  }

  for (level in seq.int(from = prepared$top_level_level, to = 0L, by = -1L)) {
    refinement.start <- proc.time()[["elapsed"]]
    refinement <- grip.geodesic.misf.refine.level(
      prepared = prepared,
      coords = coords,
      level = level,
      local_nbrs = cfg$refinement_local_nbrs,
      landmark_count = cfg$refinement_landmark_count,
      pair_mode = cfg$refinement_pair_mode,
      anchor_weight = cfg$refinement_anchor_weight,
      anchor_weight_end = cfg$refinement_anchor_weight_end,
      continuation = cfg$refinement_continuation,
      max_iter = cfg$refinement_max_iter,
      engine = "cpp",
      n_threads = cfg$n_threads,
      return_trace = FALSE
    )
    refinement.elapsed <- proc.time()[["elapsed"]] - refinement.start
    refinement.elapsed.total <- refinement.elapsed.total + refinement.elapsed
    cumulative <- cumulative + refinement.elapsed
    coords <- refinement$coords
    rec <- record_stage(
      case = case,
      spec = spec,
      coords = coords,
      level = level,
      stage_kind = "refinement",
      stage_label = sprintf("refine L%d", level),
      elapsed_sec = refinement.elapsed,
      cumulative_elapsed = cumulative,
      rows = stage.rows,
      snapshots = stage.snapshots
    )
    stage.rows <- rec$rows
    stage.snapshots <- rec$snapshots
  }

  final.start <- proc.time()[["elapsed"]]
  final.fit <- grip.geodesic.misf.final.polish(
    prepared = prepared,
    coords = coords,
    max_iter = cfg$final_polish_max_iter,
    engine = "cpp",
    n_threads = cfg$n_threads,
    return_trace = FALSE
  )
  final.elapsed <- proc.time()[["elapsed"]] - final.start
  cumulative <- cumulative + final.elapsed
  coords <- final.fit$coords
  rec <- record_stage(
    case = case,
    spec = spec,
    coords = coords,
    level = 0L,
    stage_kind = "final",
    stage_label = "final",
    elapsed_sec = final.elapsed,
    cumulative_elapsed = cumulative,
    rows = stage.rows,
    snapshots = stage.snapshots
  )
  stage.rows <- rec$rows
  stage.snapshots <- rec$snapshots

  metrics <- compute_metrics(
    case = case,
    spec = spec,
    coords = coords,
    top_seed_source = top_seed$seed_label,
    top_seed_kind = top_seed$seed_kind,
    seed_elapsed_sec = top_seed$seed_elapsed_sec,
    placement_elapsed_sec = placement.elapsed.total,
    refinement_elapsed_sec = refinement.elapsed.total,
    final_elapsed_sec = final.elapsed,
    note = sprintf("Top-level scaffold reused from Phase C (%s)", top_seed$seed_label)
  )
  aligned.final <- align_to_target_nd(coords, case$truth, allow.reflection = TRUE)

  list(
    coords = coords,
    display_coords = aligned.final$aligned,
    metrics = metrics,
    stage_trace = do.call(rbind, stage.rows),
    stage_display = stage.snapshots
  )
}

safe_run_method <- function(case, spec, top_seed) {
  message("  lower-level placement: ", spec$label)
  tryCatch(
    run_method(case, spec, top_seed = top_seed),
    error = function(e) {
      coords <- matrix(NA_real_, nrow = case$n, ncol = cfg$dim)
      list(
        coords = coords,
        display_coords = coords,
        metrics = compute_metrics(
          case = case,
          spec = spec,
          coords = case$truth,
          top_seed_source = top_seed$seed_label,
          top_seed_kind = top_seed$seed_kind,
          seed_elapsed_sec = top_seed$seed_elapsed_sec,
          placement_elapsed_sec = NA_real_,
          refinement_elapsed_sec = NA_real_,
          final_elapsed_sec = NA_real_,
          note = paste("failed:", conditionMessage(e))
        ),
        stage_trace = data.frame(),
        stage_display = list()
      )
    }
  )
}

save_final_panel <- function(case_result, file_path) {
  entries <- c(
    list(list(
      display_coords = case_result$case$truth,
      metrics = data.frame(
        method_id = "reference",
        method_label = "Reference surface",
        gmds_stress = NA_real_,
        procrustes_rmse = 0,
        stringsAsFactors = FALSE
      )
    )),
    case_result$methods
  )

  grDevices::png(file_path, width = 2400, height = 2400, res = 220, bg = "#ffffff", type = "cairo")
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(3L, 3L), mar = c(1.1, 1.1, 3.0, 0.4), oma = c(0, 0, 1.0, 0))

  for (i in seq_len(9L)) {
    if (i > length(entries)) {
      graphics::plot.new()
      next
    }
    entry <- entries[[i]]
    row <- entry$metrics[1L, , drop = FALSE]
    if (all(is.finite(entry$display_coords))) {
      grip.plot(
        coords = entry$display_coords,
        edges = case_result$case$edges,
        projection = "ortho",
        azimuth = 35,
        elevation = 24,
        vertex.col = if (identical(row$method_id[[1L]], "reference")) "#bc6c25" else "#3a5a40",
        edge.col = "#adb5bd",
        main = ""
      )
    } else {
      graphics::plot.new()
      graphics::text(0.5, 0.55, labels = row$method_label[[1L]], cex = 1.0, font = 2L)
      graphics::text(0.5, 0.40, labels = "layout failed", cex = 0.95, col = "#8d0801")
    }
    ttl <- if (identical(row$method_id[[1L]], "reference")) {
      "Reference surface"
    } else {
      sprintf(
        "%s\nsigma %s, rho %s",
        row$method_label[[1L]],
        fmt_num(row$gmds_stress[[1L]], 4L),
        fmt_num(row$procrustes_rmse[[1L]], 4L)
      )
    }
    graphics::mtext(ttl, side = 3L, line = 0.3, cex = 0.82)
  }
  graphics::mtext(case_result$case$label, side = 3L, outer = TRUE, line = -0.3, cex = 1.15, font = 2L)
}

save_stage_panel <- function(case_result, file_path) {
  df <- case_result$metrics
  best.sigma.id <- df$method_id[[which.min(df$gmds_stress)]]
  best.rho.id <- df$method_id[[which.min(df$procrustes_rmse)]]
  chosen.ids <- unique(c("geodesic_insert", best.sigma.id, best.rho.id))
  methods <- Filter(function(x) x$metrics$method_id[[1L]] %in% chosen.ids, case_result$methods)
  if (!length(methods)) {
    return(invisible(NULL))
  }
  ncol.panel <- max(vapply(methods, function(method) length(method$stage_display), integer(1L)))
  if (!is.finite(ncol.panel) || ncol.panel < 1L) {
    return(invisible(NULL))
  }

  grDevices::png(file_path, width = 420L * ncol.panel, height = 320L * length(methods), res = 180, bg = "#ffffff", type = "cairo")
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(length(methods), ncol.panel), mar = c(0.8, 0.8, 2.7, 0.3), oma = c(0, 0, 1.0, 0))

  for (method in methods) {
    for (j in seq_len(ncol.panel)) {
      if (j <= length(method$stage_display)) {
        snap <- method$stage_display[[j]]
        keep <- stats::complete.cases(snap$display_coords)
        stage.edges <- case_result$case$edges
        if (!all(keep)) {
          good.edges <- keep[stage.edges[, 1L]] & keep[stage.edges[, 2L]]
          stage.edges <- stage.edges[good.edges, , drop = FALSE]
        }
        if (sum(keep) >= 2L && nrow(stage.edges) > 0L) {
          grip.plot(
            coords = snap$display_coords,
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
            snap$display_coords[keep, 1L],
            snap$display_coords[keep, 2L],
            pch = 16,
            col = "#355070"
          )
        } else {
          graphics::plot.new()
        }
        graphics::mtext(
          sprintf("%s\n%s", method$metrics$short_label[[1L]], snap$stage_label),
          side = 3L,
          line = 0.3,
          cex = 0.74
        )
      } else {
        graphics::plot.new()
      }
    }
  }
  graphics::mtext(sprintf("%s: lower-level placement trajectories", case_result$case$label), side = 3L, outer = TRUE, line = -0.3, cex = 1.05, font = 2L)
}

save_stage_trace_grid <- function(stage_df, cases, output_path) {
  case.ids <- vapply(cases, `[[`, character(1L), "id")
  cols <- c(
    geodesic_insert = "#1b4332",
    weighted_kk_insert = "#355070",
    kk_insert = "#6d597a",
    fr_insert = "#9c6644",
    grip_insert = "#bc6c25",
    weighted_grip_insert = "#3a5a40"
  )
  grDevices::png(output_path, width = 2400L, height = 1500L, res = 180, bg = "#ffffff", type = "cairo")
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(2L, 2L), mar = c(3.5, 3.7, 2.2, 0.8), oma = c(0, 0, 1.0, 0))

  for (case.id in case.ids) {
    df <- stage_df[stage_df$case_id == case.id, , drop = FALSE]
    if (!nrow(df)) {
      graphics::plot.new()
      next
    }
    yvals <- log10(pmax(df$energy, 1e-12))
    graphics::plot(
      NA,
      xlim = c(1, max(df$stage_index)),
      ylim = range(yvals),
      xlab = "Stage step",
      ylab = "log10 active energy",
      main = unique(df$case_label),
      xaxt = "n"
    )
    for (method.id in names(cols)) {
      sub <- df[df$method_id == method.id, , drop = FALSE]
      if (!nrow(sub)) {
        next
      }
      graphics::lines(
        x = sub$stage_index,
        y = log10(pmax(sub$energy, 1e-12)),
        type = "b",
        pch = 16,
        lwd = 2,
        col = cols[[method.id]]
      )
      graphics::text(
        x = sub$stage_index,
        y = log10(pmax(sub$energy, 1e-12)),
        labels = sub$stage_label,
        pos = 3,
        cex = 0.58,
        col = cols[[method.id]]
      )
    }
    graphics::axis(1, at = sort(unique(df$stage_index)))
    graphics::legend(
      "topright",
      legend = c("Geo", "W-KK", "KK", "FR", "GRIP", "W-GRIP"),
      col = cols,
      lty = 1,
      lwd = 2,
      pch = 16,
      bty = "n",
      cex = 0.75
    )
  }
  graphics::mtext("Phase D lower-level placement active-stage energy traces", side = 3L, outer = TRUE, line = -0.3, cex = 1.1, font = 2L)
}

write_case_table <- function(df) {
  rows <- vapply(seq_len(nrow(df)), function(i) {
    sprintf(
      "%s & %s & %s & %s & %s & %s & %s & %s \\\\",
      tex_escape(df$method_label[[i]]),
      fmt_time(df$placement_elapsed_sec[[i]]),
      fmt_time(df$refinement_elapsed_sec[[i]]),
      fmt_time(df$final_elapsed_sec[[i]]),
      fmt_time(df$elapsed_sec[[i]]),
      fmt_num(df$gmds_stress[[i]], 4L),
      fmt_num(df$procrustes_rmse[[i]], 4L),
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
    "Method & $t_{place}$ (s) & $t_{ref}$ (s) & $t_{final}$ (s) & $t$ (s) & $\\sigma$ & $\\rho$ & $\\alpha_{0.05}$ \\\\",
    "\\midrule",
    paste(rows, collapse = "\n"),
    "\\bottomrule",
    "\\end{tabular}}",
    "\\end{table}",
    sep = "\n"
  )
}

case_summary_paragraph <- function(case_result) {
  df <- case_result$metrics
  base <- df[df$method_id == "geodesic_insert", , drop = FALSE]
  best.sigma <- df[which.min(df$gmds_stress), , drop = FALSE]
  best.rho <- df[which.min(df$procrustes_rmse), , drop = FALSE]
  sprintf(
    "For %s, the top MISF scaffold is fixed to the best Phase~C seeded pipeline by final $\\rho$ (%s from %s). The current geodesic insertion baseline finishes at $\\sigma=%s$ and $\\rho=%s$. The best lower-level placement by final geodesic stress is %s with $\\sigma=%s$ and $\\rho=%s$, while the best placement by final shape fidelity is %s with $\\sigma=%s$ and $\\rho=%s$.",
    tex_escape(case_result$case$label),
    tex_escape(unique(df$top_seed_kind)[1L]),
    tex_escape(unique(df$top_seed_source)[1L]),
    fmt_num(base$gmds_stress[[1L]], 4L),
    fmt_num(base$procrustes_rmse[[1L]], 4L),
    tex_escape(best.sigma$method_label[[1L]]),
    fmt_num(best.sigma$gmds_stress[[1L]], 4L),
    fmt_num(best.sigma$procrustes_rmse[[1L]], 4L),
    tex_escape(best.rho$method_label[[1L]]),
    fmt_num(best.rho$gmds_stress[[1L]], 4L),
    fmt_num(best.rho$procrustes_rmse[[1L]], 4L)
  )
}

cases <- if (smoke) {
  list(
    make_case_regular(12L),
    make_case_irregular_rectangle(15L)
  )
} else {
  list(
    make_case_regular(12L),
    make_case_regular(15L),
    make_case_irregular_rectangle(15L)
  )
}

case_results <- lapply(cases, function(case) {
  message("Running Phase D lower-level placement panel for: ", case$label)
  phase_b_case <- phase_b_map[[case$top_case_id]]
  winners <- pick_phase_b_winners(phase_b_case)
  phase_c_best <- pick_phase_c_seed_kind(case$id)
  chosen <- if (identical(phase_c_best, "misf_seed_sigma")) winners$sigma else winners$rho
  top_seed <- build_injected_top_fit(case, phase_b_case, chosen, seed_kind = phase_c_best)

  methods <- lapply(method_specs, function(spec) safe_run_method(case, spec, top_seed = top_seed))
  metrics <- do.call(rbind, lapply(methods, `[[`, "metrics"))
  stage.list <- Filter(function(x) nrow(x) > 0L, lapply(methods, `[[`, "stage_trace"))
  stage_df <- if (length(stage.list)) do.call(rbind, stage.list) else data.frame()

  out <- list(
    case = case,
    methods = methods,
    metrics = metrics,
    stage_trace = stage_df,
    top_seed = list(
      seed_label = top_seed$seed_label,
      seed_kind = top_seed$seed_kind,
      seed_elapsed_sec = top_seed$seed_elapsed_sec
    )
  )
  save_final_panel(out, file.path(pdf_dir, sprintf("%s_final_grid.png", case$id)))
  save_stage_panel(out, file.path(pdf_dir, sprintf("%s_stage_grid.png", case$id)))
  out
})

metrics_df <- do.call(rbind, lapply(case_results, `[[`, "metrics"))
stage.trace.list <- Filter(function(x) nrow(x) > 0L, lapply(case_results, `[[`, "stage_trace"))
stage_trace_df <- if (length(stage.trace.list)) do.call(rbind, stage.trace.list) else data.frame()
trace_grid_path <- file.path(pdf_dir, "gmds_misf_lower_level_placement_stage_trace_grid.png")
if (nrow(stage_trace_df) > 0L) {
  save_stage_trace_grid(stage_trace_df, lapply(case_results, `[[`, "case"), trace_grid_path)
}

compact_case_result <- function(case_result) {
  list(
    case = list(
      id = case_result$case$id,
      label = case_result$case$label,
      family = case_result$case$family,
      side = case_result$case$side,
      n = case_result$case$n,
      edges = case_result$case$edges,
      truth = case_result$case$truth
    ),
    metrics = case_result$metrics,
    methods = lapply(case_result$methods, function(method) {
      list(
        display_coords = method$display_coords,
        metrics = method$metrics,
        stage_display = method$stage_display
      )
    }),
    top_seed = case_result$top_seed
  )
}

saveRDS(
  list(
    run_tag = run_tag,
    cfg = cfg,
    case_results = lapply(case_results, compact_case_result),
    metrics = metrics_df,
    stage_traces = stage_trace_df
  ),
  rds_path
)
utils::write.csv(metrics_df, metrics_csv, row.names = FALSE)
utils::write.csv(stage_trace_df, stage_csv, row.names = FALSE)

baseline <- metrics_df[metrics_df$method_id == "geodesic_insert", , drop = FALSE]
nonbaseline <- metrics_df[metrics_df$method_id != "geodesic_insert", , drop = FALSE]
better.rho <- sum(vapply(split(nonbaseline, nonbaseline$case_id), function(df) {
  base <- baseline[baseline$case_id == df$case_id[[1L]], , drop = FALSE]
  min(df$procrustes_rmse) < base$procrustes_rmse[[1L]]
}, logical(1L)))
better.sigma <- sum(vapply(split(nonbaseline, nonbaseline$case_id), function(df) {
  base <- baseline[baseline$case_id == df$case_id[[1L]], , drop = FALSE]
  min(df$gmds_stress) < base$gmds_stress[[1L]]
}, logical(1L)))

case_sections <- vapply(case_results, function(case_result) {
  df <- case_result$metrics
  paste(
    sprintf("\\section*{%s}", tex_escape(case_result$case$label)),
    case_summary_paragraph(case_result),
    "",
    sprintf(
      "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{%s. Final layouts for the fixed-top-seed lower-level placement comparison.}\\end{figure}",
      file.path(pdf_dir, sprintf("%s_final_grid.png", case_result$case$id)),
      tex_escape(case_result$case$label)
    ),
    sprintf(
      "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{%s. Levelwise trajectories for the baseline and the best lower-level placement alternatives.}\\end{figure}",
      file.path(pdf_dir, sprintf("%s_stage_grid.png", case_result$case$id)),
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
  "\\title{Phase D: MISF-GMDS Lower-Level Placement Benchmark}",
  "\\author{MISF-GMDS top-level initializer plan}",
  "\\date{2026-04-02}",
  "\\begin{document}",
  "\\maketitle",
  "\\section{Purpose}",
  "This report implements Phase~D of the MISF-GMDS initializer plan. The top MISF level is no longer the variable under study: for each case we reuse the best Phase~C seeded top-level scaffold, then compare lower-level placement policies before the same sparse level refinement and final full-graph GMDS polish.",
  "\\section{Experimental design}",
  "Each method shares the same fixed top-level scaffold, the same sparse active-level refinement settings, and the same final full-graph pure-GMDS polish. The only changed component is how newly activated lower-level vertices are placed before refinement: current geodesic anchor insertion, unweighted KK, weighted KK, FR, GRIP, or weighted GRIP.",
  sprintf(
    "Across the %d current cases, at least one alternative lower-level placement beats the current geodesic insertion baseline in final $\\rho$ on %d cases and in final $\\sigma$ on %d cases.",
    length(case_results),
    better.rho,
    better.sigma
  ),
  paste(case_sections, collapse = "\n\n"),
  "\\section*{Per-stage active-energy traces}",
  "The figure below condenses the active-subgraph GMDS energy after the top-level stage, after each lower-level placement, after each lower-level refinement, and after the final polish. It makes it easier to see whether a placement policy helps the multiscale descent enter a better basin before full-graph polish begins.",
  sprintf(
    "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{Phase~D active-stage energy traces.}\\end{figure}",
    trace_grid_path
  ),
  "\\section*{Interactive companion}",
  paste(
    "The companion HTML gallery is generated by",
    "\\texttt{tools/render-gmds-misf-lower-level-placement-html.R}.",
    "It includes all saved stage layouts for every Phase~D method as interactive \\texttt{rglwidget} panels."
  ),
  "\\end{document}"
)

writeLines(tex_lines, tex_path)

message("Wrote Phase D metrics: ", metrics_csv)
message("Wrote Phase D stage traces: ", stage_csv)
message("Wrote Phase D bundle: ", rds_path)
message("Wrote Phase D report TeX: ", tex_path)
