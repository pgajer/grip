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
  sprintf("gmds-misf-phase-e-integrated-smoke-%s", format(Sys.Date(), "%Y-%m-%d"))
} else {
  "gmds-misf-phase-e-integrated-2026-04-02"
}

design_root <- file.path(repo_root, "dev", "design")
tmp_dir <- file.path(design_root, "tmp", run_tag)
pdf_dir <- file.path(design_root, "pdf", run_tag)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

tex_path <- file.path(design_root, "pdf", "gmds_misf_phase_e_integrated_report_2026-04-02.tex")
pdf_path <- file.path(design_root, "pdf", "gmds_misf_phase_e_integrated_report_2026-04-02.pdf")
rds_path <- file.path(tmp_dir, "gmds_misf_phase_e_integrated_results.rds")
selection_csv <- file.path(tmp_dir, "gmds_misf_phase_e_integrated_selection_metrics.csv")
final_csv <- file.path(tmp_dir, "gmds_misf_phase_e_integrated_final_metrics.csv")

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
phase_b_script <- file.path(repo_root, "tools", "benchmark-gmds-misf-top-level-corrections.R")
phase_c_script <- file.path(repo_root, "tools", "benchmark-gmds-misf-seeded-pipeline.R")

ensure_bundle <- function(path, script) {
  if (file.exists(path)) {
    return(invisible(path))
  }
  status <- system2("Rscript", script)
  if (!identical(status, 0L) || !file.exists(path)) {
    stop("Failed to generate required bundle: ", path)
  }
  invisible(path)
}

ensure_bundle(phase_b_rds, phase_b_script)
ensure_bundle(phase_c_rds, phase_c_script)

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
active_level_vertices <- get("grip.geodesic.misf.active.level.vertices", envir = ns)
induced_active_graph <- get("grip.geodesic.misf.induced_active_graph", envir = ns)
place_level_with_layout <- get("grip.geodesic.misf.place.level.with.layout", envir = ns)
prepare_active_level <- get("grip.geodesic.misf.prepare.active.level", envir = ns)

cfg <- list(
  phase_seed = 20260402L,
  amplitude = 0.35,
  dim = 3L,
  candidate_ids = c("cmdscale", "kk", "grip", "weighted_grip", "weighted_grip_polish_lgkk"),
  insertion_layout_k = 6L,
  lookahead_refine_local_nbrs = if (smoke) 3L else 4L,
  lookahead_refine_landmarks = 2L,
  lookahead_refine_max_iter = if (smoke) 2L else 3L,
  full_refine_local_nbrs = if (smoke) 3L else 4L,
  full_refine_landmarks = 2L,
  full_refine_max_iter = if (smoke) 2L else 3L,
  final_polish_max_iter = if (smoke) 2L else 4L,
  refinement_anchor_weight = 0.05,
  refinement_anchor_weight_end = 0.01,
  refinement_continuation = "linear",
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

edge_floor_ratio <- function(coords, edges) {
  coords <- as.matrix(coords)
  edges <- as.matrix(edges)
  if (!nrow(edges)) {
    return(NA_real_)
  }
  lengths <- sqrt(rowSums(
    (coords[edges[, 1L], , drop = FALSE] - coords[edges[, 2L], , drop = FALSE])^2
  ))
  med <- stats::median(lengths)
  if (!is.finite(med) || med <= 0) {
    return(NA_real_)
  }
  as.double(stats::quantile(lengths, probs = 0.05, names = FALSE)) / med
}

mesh_roughness <- function(coords, adj_list, edges) {
  centered <- sweep(coords, 2L, colMeans(coords), FUN = "-", check.margin = FALSE)
  if (!nrow(edges)) {
    return(NA_real_)
  }
  median.edge <- stats::median(sqrt(rowSums(
    (centered[edges[, 1L], , drop = FALSE] - centered[edges[, 2L], , drop = FALSE])^2
  )))
  if (!is.finite(median.edge) || median.edge <= 0) {
    return(NA_real_)
  }
  residuals <- vapply(seq_len(nrow(centered)), function(i) {
    nbrs <- adj_list[[i]]
    if (!length(nbrs)) {
      return(0)
    }
    delta <- centered[i, ] - colMeans(centered[nbrs, , drop = FALSE])
    sum(delta^2)
  }, numeric(1L))
  sqrt(mean(residuals)) / median.edge
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

complete_partial_coords <- function(coords) {
  coords <- as.matrix(coords)
  keep <- stats::complete.cases(coords)
  if (all(keep)) {
    return(list(coords = coords, missing = which(!keep)))
  }
  if (!any(keep)) {
    stop("Cannot complete a partial coordinate matrix with no finite rows.")
  }
  fill.center <- colMeans(coords[keep, , drop = FALSE])
  out <- coords
  out[!keep, ] <- matrix(fill.center, nrow = sum(!keep), ncol = ncol(coords), byrow = TRUE)
  list(coords = out, missing = which(!keep))
}

plot_partial_layout <- function(coords,
                                edges,
                                projection = "ortho",
                                azimuth = 35,
                                elevation = 24,
                                vertex.col = "#355070",
                                edge.col = "#adb5bd",
                                main = "") {
  coords <- as.matrix(coords)
  keep <- stats::complete.cases(coords)
  if (!any(keep)) {
    graphics::plot.new()
    return(invisible(NULL))
  }
  index.map <- integer(nrow(coords))
  index.map[keep] <- seq_len(sum(keep))
  display.edges <- as.matrix(edges)
  if (nrow(display.edges)) {
    good.edges <- keep[display.edges[, 1L]] & keep[display.edges[, 2L]]
    display.edges <- display.edges[good.edges, , drop = FALSE]
    if (nrow(display.edges)) {
      display.edges <- cbind(
        as.integer(index.map[display.edges[, 1L]]),
        as.integer(index.map[display.edges[, 2L]])
      )
    }
  }
  grip.plot(
    coords = coords[keep, , drop = FALSE],
    edges = display.edges,
    projection = projection,
    azimuth = azimuth,
    elevation = elevation,
    vertex.col = vertex.col,
    edge.col = edge.col,
    main = main
  )
  invisible(NULL)
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
    prepared = grip.prepare.misf.geodesic.mds(
      edges = bundle$edges,
      n = bundle$n,
      edge_weights = bundle$edge_weights,
      tie_mode = "average",
      dim = cfg$dim,
      top_level_mode = "skip",
      seed = cfg$phase_seed + side
    )
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
    prepared = grip.prepare.misf.geodesic.mds(
      edges = bundle$edges,
      n = bundle$n,
      edge_weights = bundle$edge_weights,
      tie_mode = "average",
      dim = cfg$dim,
      top_level_mode = "skip",
      seed = cfg$phase_seed + side + 100L
    )
  )
}

phase_b_bundle <- readRDS(phase_b_rds)
phase_b_map <- setNames(phase_b_bundle$case_results, vapply(phase_b_bundle$case_results, function(x) x$case$id, character(1L)))
phase_c_bundle <- readRDS(phase_c_rds)
phase_c_map <- setNames(phase_c_bundle$case_results, vapply(phase_c_bundle$case_results, function(x) x$case$id, character(1L)))

phase_c_best_seed_method <- function(case_result) {
  seeded <- case_result$metrics[case_result$metrics$method_id %in% c("misf_seed_sigma", "misf_seed_rho"), , drop = FALSE]
  seeded <- seeded[order(seeded$procrustes_rmse, seeded$gmds_stress), , drop = FALSE]
  seeded$method_id[[1L]]
}

phase_c_seed_source_to_candidate_id <- function(seed_source) {
  switch(
    seed_source,
    "cMDS" = "cmdscale",
    "KK" = "kk",
    "GRIP" = "grip",
    "Weighted GRIP" = "weighted_grip",
    "Weighted GRIP + polish LGKK" = "weighted_grip_polish_lgkk",
    stop("Unsupported Phase C seed source: ", seed_source)
  )
}

partial_coords <- function(coords, vertex_ids, n) {
  out <- matrix(NA_real_, nrow = n, ncol = ncol(coords))
  out[vertex_ids, ] <- coords
  out
}

build_injected_top_fit <- function(case,
                                   phase_b_case_result,
                                   phase_b_method_result,
                                   with_frame_stub = FALSE) {
  source.ids <- as.integer(phase_b_case_result$case$top_vertex_ids)
  target.ids <- as.integer(case$prepared$top_level_vertices)
  idx <- match(target.ids, source.ids)
  if (anyNA(idx)) {
    stop("Phase B top-level vertices do not align with full-case top vertices for ", case$id)
  }
  top.coords <- as.matrix(phase_b_method_result$pure$coords[idx, , drop = FALSE])
  top.prepared <- case$prepared$top_level_prepared
  top.score <- grip.score.geodesic.mds(coords = top.coords, prepared = top.prepared)
  seed.elapsed <- sum(c(
    phase_b_method_result$initializer$metrics$elapsed_sec[[1L]],
    phase_b_method_result$anchor$metrics$elapsed_sec[[1L]],
    phase_b_method_result$pure$metrics$elapsed_sec[[1L]]
  ), na.rm = TRUE)
  frames <- if (isTRUE(with_frame_stub)) list(top.coords, top.coords) else list(top.coords)
  list(
    fit = list(
      coords = top.coords,
      trace = data.frame(iter = 0L, stringsAsFactors = FALSE),
      frames = frames,
      prepared = top.prepared,
      score = top.score,
      restart_summary = data.frame(),
      best_restart = 1L,
      best_restart_row = data.frame(),
      vertex_ids = target.ids,
      coords_full = partial_coords(top.coords, target.ids, case$prepared$n),
      injected = TRUE,
      injected_seed_label = phase_b_method_result$method_label,
      injected_seed_elapsed_sec = seed.elapsed
    ),
    seed_label = phase_b_method_result$method_label,
    seed_elapsed_sec = as.double(seed.elapsed),
    method_id = phase_b_method_result$method_id
  )
}

compute_final_metrics <- function(case, method_id, method_label, coords, elapsed_sec, note = "") {
  score_df <- grip.score.geodesic.mds(coords = coords, prepared = case$prepared)
  aligned <- align_to_target_nd(coords, case$truth, allow.reflection = TRUE)
  data.frame(
    case_id = case$id,
    case_label = case$label,
    family = case$family,
    side = case$side,
    n = case$n,
    method_id = method_id,
    method_label = method_label,
    elapsed_sec = as.double(elapsed_sec),
    gmds_stress = score_df$gmds.stress[[1L]],
    gmds_energy = score_df$gmds.energy[[1L]],
    procrustes_rmse = aligned$rmse,
    roughness = mesh_roughness(coords, case$prepared$adj_list, case$edges),
    area_q05_ratio = area_floor_ratio(coords, case$triangles),
    note = note,
    stringsAsFactors = FALSE
  )
}

compute_lookahead_metrics <- function(case, prepared, coords, level) {
  active.vertices <- active_level_vertices(prepared, level)
  active.full.coords <- as.matrix(coords)
  active.coords <- active.full.coords[active.vertices, , drop = FALSE]
  active.prepared <- prepare_active_level(
    prepared = prepared,
    active_vertices = active.vertices,
    local_nbrs = cfg$lookahead_refine_local_nbrs,
    landmark_count = cfg$lookahead_refine_landmarks,
    pair_mode = "full"
  )
  active.graph <- induced_active_graph(prepared, active.vertices)
  active.score <- grip.score.geodesic.mds(active.full.coords, prepared = active.prepared)
  active.align <- align_to_target_nd(active.coords, case$truth[active.vertices, , drop = FALSE], allow.reflection = TRUE)
  data.frame(
    active_level = level,
    active_n = length(active.vertices),
    lookahead_sigma = active.score$gmds.stress[[1L]],
    lookahead_energy = active.score$gmds.energy[[1L]],
    lookahead_rho = active.align$rmse,
    lookahead_roughness = mesh_roughness(active.coords, active.graph$adj_list, active.graph$edges),
    lookahead_edge_floor = edge_floor_ratio(active.coords, active.graph$edges),
    stringsAsFactors = FALSE
  )
}

run_candidate <- function(case, phase_b_case_result, phase_b_method_result) {
  top.seed <- build_injected_top_fit(case, phase_b_case_result, phase_b_method_result, with_frame_stub = FALSE)
  prepared <- case$prepared
  prepared$top_level_fit <- top.seed$fit
  top.display <- align_partial_to_truth(top.seed$fit$coords_full, case$truth)$aligned

  next.level <- prepared$top_level_level - 1L
  lookahead.coords <- top.seed$fit$coords_full
  if (next.level >= 0L) {
    placement <- place_level_with_layout(
      prepared = prepared,
      coords = lookahead.coords,
      level = next.level,
      method = "weighted_kk",
      layout_k = cfg$insertion_layout_k,
      seed = cfg$phase_seed + case$side + match(phase_b_method_result$method_id, cfg$candidate_ids)
    )
    lookahead.partial <- placement$coords
    lookahead.completed <- complete_partial_coords(lookahead.partial)
    refined <- grip:::grip.geodesic.misf.refine.level(
      prepared = prepared,
      coords = lookahead.completed$coords,
      level = next.level,
      local_nbrs = cfg$lookahead_refine_local_nbrs,
      landmark_count = cfg$lookahead_refine_landmarks,
      pair_mode = "sparse",
      anchor_weight = cfg$refinement_anchor_weight,
      anchor_weight_end = cfg$refinement_anchor_weight_end,
      continuation = cfg$refinement_continuation,
      max_iter = cfg$lookahead_refine_max_iter,
      engine = "cpp",
      n_threads = cfg$n_threads,
      recenter = FALSE,
      return_trace = FALSE
    )
    lookahead.metrics <- compute_lookahead_metrics(case, prepared, refined$coords, next.level)
    lookahead.coords <- refined$coords
    if (length(lookahead.completed$missing)) {
      lookahead.coords[lookahead.completed$missing, ] <- NA_real_
    }
  } else {
    lookahead.metrics <- data.frame(
      active_level = prepared$top_level_level,
      active_n = length(prepared$top_level_vertices),
      lookahead_sigma = top.seed$fit$score$gmds.stress[[1L]],
      lookahead_energy = top.seed$fit$score$gmds.energy[[1L]],
      lookahead_rho = align_to_target_nd(top.seed$fit$coords, case$truth[prepared$top_level_vertices, , drop = FALSE], allow.reflection = TRUE)$rmse,
      lookahead_roughness = NA_real_,
      lookahead_edge_floor = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  lookahead.display <- align_partial_to_truth(lookahead.coords, case$truth)$aligned

  full.start <- proc.time()[["elapsed"]]
  full.fit <- grip.optimize.misf.geodesic.mds(
    prepared = prepared,
    insertion_mode = "weighted_kk",
    insertion_layout_k = cfg$insertion_layout_k,
    refinement_local_nbrs = cfg$full_refine_local_nbrs,
    refinement_landmark_count = cfg$full_refine_landmarks,
    refinement_pair_mode = "sparse",
    refinement_anchor_weight = cfg$refinement_anchor_weight,
    refinement_anchor_weight_end = cfg$refinement_anchor_weight_end,
    refinement_continuation = cfg$refinement_continuation,
    refinement_max_iter = cfg$full_refine_max_iter,
    refinement_engine = "cpp",
    final_polish_max_iter = cfg$final_polish_max_iter,
    final_polish_engine = "cpp",
    n_threads = cfg$n_threads,
    return_trace = FALSE,
    return_frames = FALSE,
    seed = cfg$phase_seed + case$side
  )
  full.elapsed <- proc.time()[["elapsed"]] - full.start
  final.metrics <- compute_final_metrics(
    case = case,
    method_id = phase_b_method_result$method_id,
    method_label = phase_b_method_result$method_label,
    coords = full.fit$coords,
    elapsed_sec = full.elapsed,
    note = "Phase E full weighted-KK lower-level placement pipeline"
  )

  list(
    method_id = phase_b_method_result$method_id,
    method_label = phase_b_method_result$method_label,
    top_seed = top.seed,
    top_display = top.display,
    lookahead_display = lookahead.display,
    final_display = align_to_target_nd(full.fit$coords, case$truth, allow.reflection = TRUE)$aligned,
    lookahead_metrics = lookahead.metrics,
    final_metrics = final.metrics
  )
}

run_candidate_with_frames <- function(case, phase_b_case_result, phase_b_method_result) {
  top.seed <- build_injected_top_fit(case, phase_b_case_result, phase_b_method_result, with_frame_stub = TRUE)
  prepared <- case$prepared
  prepared$top_level_fit <- top.seed$fit
  fit <- grip.optimize.misf.geodesic.mds(
    prepared = prepared,
    insertion_mode = "weighted_kk",
    insertion_layout_k = cfg$insertion_layout_k,
    refinement_local_nbrs = cfg$full_refine_local_nbrs,
    refinement_landmark_count = cfg$full_refine_landmarks,
    refinement_pair_mode = "sparse",
    refinement_anchor_weight = cfg$refinement_anchor_weight,
    refinement_anchor_weight_end = cfg$refinement_anchor_weight_end,
    refinement_continuation = cfg$refinement_continuation,
    refinement_max_iter = cfg$full_refine_max_iter,
    refinement_engine = "cpp",
    final_polish_max_iter = cfg$final_polish_max_iter,
    final_polish_engine = "cpp",
    n_threads = cfg$n_threads,
    return_trace = FALSE,
    return_frames = TRUE,
    seed = cfg$phase_seed + case$side
  )
  list(
    top_level = align_partial_to_truth(fit$frames$after_top_level, case$truth)$aligned,
    after_insertion = align_partial_to_truth(fit$frames$after_insertion, case$truth)$aligned,
    after_refinement = align_partial_to_truth(fit$frames$after_refinement, case$truth)$aligned,
    final = align_partial_to_truth(fit$frames$final, case$truth)$aligned
  )
}

save_candidate_grid <- function(case_result, stage = c("top", "lookahead"), output_path) {
  stage <- match.arg(stage)
  entries <- c(
    list(list(
      display = case_result$case$truth,
      title = "Reference surface",
      subtitle = case_result$case$label
    )),
    lapply(case_result$candidates, function(candidate) {
      look <- candidate$lookahead_metrics[1L, , drop = FALSE]
      list(
        display = if (identical(stage, "top")) candidate$top_display else candidate$lookahead_display,
        title = candidate$method_label,
        subtitle = if (identical(stage, "top")) {
          sprintf("top seed")
        } else {
          sprintf("look sigma %s, rho %s", fmt_num(look$lookahead_sigma[[1L]], 4L), fmt_num(look$lookahead_rho[[1L]], 4L))
        }
      )
    })
  )
  n <- length(entries)
  ncol.panel <- 3L
  nrow.panel <- ceiling(n / ncol.panel)
  grDevices::png(output_path, width = 2400L, height = 820L * nrow.panel, res = 180, bg = "#ffffff", type = "cairo")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(nrow.panel, ncol.panel), mar = c(1.1, 1.1, 3.0, 0.4), oma = c(0, 0, 1.0, 0))
  for (i in seq_len(nrow.panel * ncol.panel)) {
    if (i > length(entries)) {
      graphics::plot.new()
      next
    }
    entry <- entries[[i]]
    plot_partial_layout(
      coords = entry$display,
      edges = case_result$case$edges,
      projection = "ortho",
      azimuth = 35,
      elevation = 24,
      vertex.col = if (i == 1L) "#bc6c25" else "#3a5a40",
      edge.col = "#adb5bd",
      main = ""
    )
    graphics::mtext(sprintf("%s\n%s", entry$title, entry$subtitle), side = 3L, line = 0.3, cex = 0.82)
  }
  ttl <- if (identical(stage, "top")) "Phase E corrected top-level seeds" else "Phase E one-level weighted-KK lookahead"
  graphics::mtext(sprintf("%s: %s", case_result$case$label, ttl), side = 3L, outer = TRUE, line = -0.3, cex = 1.1, font = 2L)
}

save_final_comparison_grid <- function(case_result, output_path) {
  entries <- case_result$comparison_entries
  n <- length(entries)
  ncol.panel <- 3L
  nrow.panel <- ceiling(n / ncol.panel)
  grDevices::png(output_path, width = 2400L, height = 820L * nrow.panel, res = 180, bg = "#ffffff", type = "cairo")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(nrow.panel, ncol.panel), mar = c(1.1, 1.1, 3.0, 0.4), oma = c(0, 0, 1.0, 0))
  for (i in seq_len(nrow.panel * ncol.panel)) {
    if (i > length(entries)) {
      graphics::plot.new()
      next
    }
    entry <- entries[[i]]
    plot_partial_layout(
      coords = entry$display_coords,
      edges = case_result$case$edges,
      projection = "ortho",
      azimuth = 35,
      elevation = 24,
      vertex.col = if (entry$method_id == "reference") "#bc6c25" else "#355070",
      edge.col = "#adb5bd",
      main = ""
    )
    graphics::mtext(
      sprintf("%s\nsigma %s, rho %s", entry$method_label, fmt_num(entry$gmds_stress, 4L), fmt_num(entry$procrustes_rmse, 4L)),
      side = 3L,
      line = 0.3,
      cex = 0.78
    )
  }
  graphics::mtext(sprintf("%s: final pipeline comparison", case_result$case$label), side = 3L, outer = TRUE, line = -0.3, cex = 1.1, font = 2L)
}

save_selected_stage_grid <- function(case_result, output_path) {
  methods <- case_result$selected_stage_methods
  stage.names <- c("top_level", "after_insertion", "after_refinement", "final")
  grDevices::png(output_path, width = 420L * length(stage.names), height = 320L * length(methods), res = 180, bg = "#ffffff", type = "cairo")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(length(methods), length(stage.names)), mar = c(0.8, 0.8, 2.7, 0.3), oma = c(0, 0, 1.0, 0))
  for (method in methods) {
    for (stage.name in stage.names) {
      coords <- method$stage_display[[stage.name]]
      plot_partial_layout(
        coords = coords,
        edges = case_result$case$edges,
        projection = "ortho",
        azimuth = 35,
        elevation = 24,
        vertex.col = "#355070",
        edge.col = "#c6d1db",
        main = ""
      )
      label <- switch(
        stage.name,
        top_level = "top",
        after_insertion = "after insertion",
        after_refinement = "after refinement",
        final = "final",
        stage.name
      )
      graphics::mtext(sprintf("%s\n%s", method$method_label, label), side = 3L, line = 0.3, cex = 0.72)
    }
  }
  graphics::mtext(sprintf("%s: selected pipeline stage layouts", case_result$case$label), side = 3L, outer = TRUE, line = -0.3, cex = 1.05, font = 2L)
}

write_candidate_table <- function(df) {
  rows <- vapply(seq_len(nrow(df)), function(i) {
    sprintf(
      "%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
      tex_escape(df$method_label[[i]]),
      fmt_num(df$lookahead_sigma[[i]], 4L),
      fmt_num(df$lookahead_rho[[i]], 4L),
      fmt_num(df$lookahead_roughness[[i]], 4L),
      fmt_num(df$lookahead_edge_floor[[i]], 4L),
      as.integer(df$proxy_rank_sum[[i]]),
      fmt_num(df$final_sigma[[i]], 4L),
      fmt_num(df$final_rho[[i]], 4L),
      fmt_time(df$elapsed_sec[[i]])
    )
  }, character(1L))
  paste(
    "\\begin{table}[H]",
    "\\centering",
    "\\scriptsize",
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{lrrrrrrrr}",
    "\\toprule",
    "Candidate & $\\sigma_{look}$ & $\\rho_{look}$ & $\\eta_{look}$ & $q_{0.05}/q_{0.50}$ & proxy rank & $\\sigma$ & $\\rho$ & $t$ (s) \\\\",
    "\\midrule",
    paste(rows, collapse = "\n"),
    "\\bottomrule",
    "\\end{tabular}}",
    "\\end{table}",
    sep = "\n"
  )
}

cases <- if (smoke) {
  list(make_regular_case(12L), make_irregular_rectangle_case(15L))
} else {
  list(make_regular_case(12L), make_regular_case(15L), make_irregular_rectangle_case(15L))
}

case_results <- lapply(cases, function(case) {
  message("Running Phase E integrated panel for: ", case$label)
  phase_b_case <- phase_b_map[[case$top_case_id]]
  if (is.null(phase_b_case)) {
    stop("Missing Phase B result for ", case$top_case_id)
  }
  phase_c_case <- phase_c_map[[case$id]]
  if (is.null(phase_c_case)) {
    stop("Missing Phase C result for ", case$id)
  }

  methods_all <- setNames(phase_b_case$methods, vapply(phase_b_case$methods, `[[`, character(1L), "method_id"))
  candidates <- lapply(cfg$candidate_ids, function(method_id) {
    message("  candidate: ", method_id)
    run_candidate(case, phase_b_case, methods_all[[method_id]])
  })
  names(candidates) <- cfg$candidate_ids

  selection_df <- do.call(rbind, lapply(candidates, function(candidate) {
    look <- candidate$lookahead_metrics[1L, , drop = FALSE]
    fin <- candidate$final_metrics[1L, , drop = FALSE]
    data.frame(
      case_id = case$id,
      case_label = case$label,
      method_id = candidate$method_id,
      method_label = candidate$method_label,
      lookahead_sigma = look$lookahead_sigma[[1L]],
      lookahead_rho = look$lookahead_rho[[1L]],
      lookahead_roughness = look$lookahead_roughness[[1L]],
      lookahead_edge_floor = look$lookahead_edge_floor[[1L]],
      final_sigma = fin$gmds_stress[[1L]],
      final_rho = fin$procrustes_rmse[[1L]],
      final_area_q05_ratio = fin$area_q05_ratio[[1L]],
      elapsed_sec = fin$elapsed_sec[[1L]],
      stringsAsFactors = FALSE
    )
  }))
  selection_df$rank_sigma <- rank(selection_df$lookahead_sigma, ties.method = "min")
  selection_df$rank_roughness <- rank(selection_df$lookahead_roughness, ties.method = "min")
  selection_df$rank_edge_floor <- rank(-selection_df$lookahead_edge_floor, ties.method = "min")
  selection_df$proxy_rank_sum <- selection_df$rank_sigma + selection_df$rank_roughness + selection_df$rank_edge_floor

  proxy.method.id <- selection_df$method_id[[which.min(selection_df$proxy_rank_sum)]]
  oracle.sigma.id <- selection_df$method_id[[which.min(selection_df$final_sigma)]]
  oracle.rho.id <- selection_df$method_id[[which.min(selection_df$final_rho)]]

  phase_c_best_id <- phase_c_best_seed_method(phase_c_case)
  phase_c_best_row <- phase_c_case$metrics[phase_c_case$metrics$method_id == phase_c_best_id, , drop = FALSE]
  fixed_top_method_id <- phase_c_seed_source_to_candidate_id(phase_c_best_row$seed_source[[1L]])

  selected.ids <- unique(c(proxy.method.id, fixed_top_method_id, oracle.rho.id))
  selected.stage.methods <- lapply(selected.ids, function(method_id) {
    method <- methods_all[[method_id]]
    stage.display <- run_candidate_with_frames(case, phase_b_case, method)
    list(method_id = method_id, method_label = method$method_label, stage_display = stage.display)
  })

  candidate_grid_top <- file.path(pdf_dir, sprintf("%s_top_seed_grid.png", case$id))
  candidate_grid_lookahead <- file.path(pdf_dir, sprintf("%s_lookahead_grid.png", case$id))
  final_grid <- file.path(pdf_dir, sprintf("%s_final_comparison_grid.png", case$id))
  stage_grid <- file.path(pdf_dir, sprintf("%s_selected_stage_grid.png", case$id))

  comparison_entries <- list(
    list(
      method_id = "cmd_pure_gmds",
      method_label = phase_c_case$methods[[2L]]$metrics$method_label[[1L]],
      display_coords = phase_c_case$methods[[2L]]$display_coords,
      gmds_stress = phase_c_case$methods[[2L]]$metrics$gmds_stress[[1L]],
      procrustes_rmse = phase_c_case$methods[[2L]]$metrics$procrustes_rmse[[1L]]
    ),
    list(
      method_id = phase_c_best_id,
      method_label = paste0("Phase C best seeded (", phase_c_best_row$seed_source[[1L]], ")"),
      display_coords = phase_c_case$methods[[match(phase_c_best_id, vapply(phase_c_case$methods, function(x) x$metrics$method_id[[1L]], character(1L)))]]$display_coords,
      gmds_stress = phase_c_best_row$gmds_stress[[1L]],
      procrustes_rmse = phase_c_best_row$procrustes_rmse[[1L]]
    ),
    list(
      method_id = fixed_top_method_id,
      method_label = paste0("Fixed-top W-KK (", methods_all[[fixed_top_method_id]]$method_label, ")"),
      display_coords = candidates[[fixed_top_method_id]]$final_display,
      gmds_stress = candidates[[fixed_top_method_id]]$final_metrics$gmds_stress[[1L]],
      procrustes_rmse = candidates[[fixed_top_method_id]]$final_metrics$procrustes_rmse[[1L]]
    ),
    list(
      method_id = proxy.method.id,
      method_label = paste0("Phase E proxy winner (", methods_all[[proxy.method.id]]$method_label, ")"),
      display_coords = candidates[[proxy.method.id]]$final_display,
      gmds_stress = candidates[[proxy.method.id]]$final_metrics$gmds_stress[[1L]],
      procrustes_rmse = candidates[[proxy.method.id]]$final_metrics$procrustes_rmse[[1L]]
    )
  )
  if (!oracle.rho.id %in% vapply(comparison_entries, `[[`, character(1L), "method_id")) {
    comparison_entries[[length(comparison_entries) + 1L]] <- list(
      method_id = oracle.rho.id,
      method_label = paste0("Oracle best final rho (", methods_all[[oracle.rho.id]]$method_label, ")"),
      display_coords = candidates[[oracle.rho.id]]$final_display,
      gmds_stress = candidates[[oracle.rho.id]]$final_metrics$gmds_stress[[1L]],
      procrustes_rmse = candidates[[oracle.rho.id]]$final_metrics$procrustes_rmse[[1L]]
    )
  }
  comparison_entries <- c(list(list(
    method_id = "reference",
    method_label = "Reference surface",
    display_coords = case$truth,
    gmds_stress = 0,
    procrustes_rmse = 0
  )), comparison_entries)

  case.out <- list(
    case = case,
    candidates = candidates,
    selection_metrics = selection_df,
    proxy_method_id = proxy.method.id,
    oracle_sigma_method_id = oracle.sigma.id,
    oracle_rho_method_id = oracle.rho.id,
    fixed_top_method_id = fixed_top_method_id,
    phase_c_best_id = phase_c_best_id,
    comparison_entries = comparison_entries,
    selected_stage_methods = selected.stage.methods
  )

  save_candidate_grid(case.out, stage = "top", output_path = candidate_grid_top)
  save_candidate_grid(case.out, stage = "lookahead", output_path = candidate_grid_lookahead)
  save_final_comparison_grid(case.out, output_path = final_grid)
  save_selected_stage_grid(case.out, output_path = stage_grid)

  case.out
})

selection_df <- do.call(rbind, lapply(case_results, `[[`, "selection_metrics"))
final_rows <- do.call(rbind, lapply(case_results, function(case_result) {
  rbind(
    data.frame(
      case_id = case_result$case$id,
      case_label = case_result$case$label,
      baseline = "phase_c_best_seeded",
      method_id = case_result$phase_c_best_id,
      method_label = case_result$comparison_entries[[3L]]$method_label,
      gmds_stress = case_result$comparison_entries[[3L]]$gmds_stress,
      procrustes_rmse = case_result$comparison_entries[[3L]]$procrustes_rmse,
      stringsAsFactors = FALSE
    ),
    data.frame(
      case_id = case_result$case$id,
      case_label = case_result$case$label,
      baseline = "phase_d_fixed_top_weighted_kk",
      method_id = case_result$fixed_top_method_id,
      method_label = paste0("Fixed-top W-KK (", case_result$fixed_top_method_id, ")"),
      gmds_stress = case_result$candidates[[case_result$fixed_top_method_id]]$final_metrics$gmds_stress[[1L]],
      procrustes_rmse = case_result$candidates[[case_result$fixed_top_method_id]]$final_metrics$procrustes_rmse[[1L]],
      stringsAsFactors = FALSE
    ),
    data.frame(
      case_id = case_result$case$id,
      case_label = case_result$case$label,
      baseline = "phase_e_proxy",
      method_id = case_result$proxy_method_id,
      method_label = paste0("Phase E proxy winner (", case_result$proxy_method_id, ")"),
      gmds_stress = case_result$candidates[[case_result$proxy_method_id]]$final_metrics$gmds_stress[[1L]],
      procrustes_rmse = case_result$candidates[[case_result$proxy_method_id]]$final_metrics$procrustes_rmse[[1L]],
      stringsAsFactors = FALSE
    )
  )
}))

utils::write.csv(selection_df, selection_csv, row.names = FALSE)
utils::write.csv(final_rows, final_csv, row.names = FALSE)

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
    candidates = lapply(case_result$candidates, function(candidate) {
      list(
        method_id = candidate$method_id,
        method_label = candidate$method_label,
        top_display = candidate$top_display,
        lookahead_display = candidate$lookahead_display,
        final_display = candidate$final_display,
        lookahead_metrics = candidate$lookahead_metrics,
        final_metrics = candidate$final_metrics
      )
    }),
    selection_metrics = case_result$selection_metrics,
    proxy_method_id = case_result$proxy_method_id,
    oracle_sigma_method_id = case_result$oracle_sigma_method_id,
    oracle_rho_method_id = case_result$oracle_rho_method_id,
    fixed_top_method_id = case_result$fixed_top_method_id,
    phase_c_best_id = case_result$phase_c_best_id,
    comparison_entries = case_result$comparison_entries,
    selected_stage_methods = case_result$selected_stage_methods
  )
}

saveRDS(
  list(
    run_tag = run_tag,
    cfg = cfg,
    case_results = lapply(case_results, compact_case_result),
    selection_metrics = selection_df,
    final_metrics = final_rows
  ),
  rds_path
)

proxy_matches_oracle_rho <- sum(vapply(case_results, function(x) identical(x$proxy_method_id, x$oracle_rho_method_id), logical(1L)))
proxy_beats_fixed_rho <- sum(vapply(case_results, function(x) {
  x$candidates[[x$proxy_method_id]]$final_metrics$procrustes_rmse[[1L]] <
    x$candidates[[x$fixed_top_method_id]]$final_metrics$procrustes_rmse[[1L]]
}, logical(1L)))
proxy_beats_fixed_sigma <- sum(vapply(case_results, function(x) {
  x$candidates[[x$proxy_method_id]]$final_metrics$gmds_stress[[1L]] <
    x$candidates[[x$fixed_top_method_id]]$final_metrics$gmds_stress[[1L]]
}, logical(1L)))

case_sections <- vapply(case_results, function(case_result) {
  candidate.table <- case_result$selection_metrics[order(case_result$selection_metrics$proxy_rank_sum, case_result$selection_metrics$final_rho), , drop = FALSE]
  proxy.row <- candidate.table[candidate.table$method_id == case_result$proxy_method_id, , drop = FALSE]
  fixed.row <- candidate.table[candidate.table$method_id == case_result$fixed_top_method_id, , drop = FALSE]
  oracle.row <- candidate.table[candidate.table$method_id == case_result$oracle_rho_method_id, , drop = FALSE]
  paste(
    sprintf("\\section*{%s}", tex_escape(case_result$case$label)),
    sprintf(
      paste(
        "Phase~E compares a small corrected top-level seed pool, scores each candidate after a one-level weighted-KK lookahead,",
        "and then runs the full weighted-KK lower-level pipeline from each candidate.",
        "The proxy winner is \\textbf{%s}; the oracle best final $\\\\rho$ candidate is \\textbf{%s};",
        "the fixed-top weighted-KK baseline inherited from Phases~C/D is \\textbf{%s}.",
        "On this case the proxy winner finishes at $\\\\sigma=%s$, $\\\\rho=%s$, compared with the fixed-top baseline at $\\\\sigma=%s$, $\\\\rho=%s$."
      ),
      tex_escape(proxy.row$method_label[[1L]]),
      tex_escape(oracle.row$method_label[[1L]]),
      tex_escape(fixed.row$method_label[[1L]]),
      fmt_num(proxy.row$final_sigma[[1L]], 4L),
      fmt_num(proxy.row$final_rho[[1L]], 4L),
      fmt_num(fixed.row$final_sigma[[1L]], 4L),
      fmt_num(fixed.row$final_rho[[1L]], 4L)
    ),
    sprintf(
      "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{%s. Corrected top-level seed candidates reused in Phase~E.}\\end{figure}",
      file.path(pdf_dir, sprintf("%s_top_seed_grid.png", case_result$case$id)),
      tex_escape(case_result$case$label)
    ),
    sprintf(
      "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{%s. One-level weighted-KK lookahead used to rank the candidate top seeds.}\\end{figure}",
      file.path(pdf_dir, sprintf("%s_lookahead_grid.png", case_result$case$id)),
      tex_escape(case_result$case$label)
    ),
    write_candidate_table(candidate.table),
    sprintf(
      "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{%s. Final Phase~E comparison: direct control, previous seeded baselines, and the new integrated Phase~E pipelines.}\\end{figure}",
      file.path(pdf_dir, sprintf("%s_final_comparison_grid.png", case_result$case$id)),
      tex_escape(case_result$case$label)
    ),
    sprintf(
      "\\begin{figure}[p]\\centering\\includegraphics[width=0.98\\linewidth]{%s}\\caption{%s. Stage layouts for the selected integrated pipelines.}\\end{figure}",
      file.path(pdf_dir, sprintf("%s_selected_stage_grid.png", case_result$case$id)),
      tex_escape(case_result$case$label)
    ),
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
  "\\title{Phase E: Integrated Seed Selection with Weighted-KK Lower-Level Placement}",
  "\\author{MISF-GMDS follow-up experiment}",
  "\\date{2026-04-02}",
  "\\begin{document}",
  "\\maketitle",
  "\\section{Purpose}",
  "Phase~E combines the main lessons of Phases~A--D. Instead of fixing the top seed in advance, it compares a small corrected top-level seed pool, scores those candidates after a one-level weighted-KK lookahead, selects a top seed by a truth-free proxy rank, and then runs the full MISF-GMDS pipeline with weighted-KK lower-level placement.",
  "\\section{Design}",
  "The candidate pool is \\{cMDS, KK, GRIP, Weighted GRIP, Weighted GRIP + polish LGKK\\}. Each candidate inherits the common Phase~B corrected top-level geometry. The one-level lookahead inserts the next MISF level with weighted KK, refines that active subgraph, and records active geodesic stress, active roughness, and an edge-floor collapse proxy. The proxy score is the rank-sum of those three quantities.",
  sprintf(
    "Across the %d current cases, the Phase~E proxy winner matches the oracle best final $\\\\rho$ candidate on %d cases, beats the Phase~D fixed-top weighted-KK baseline in final $\\\\rho$ on %d cases, and beats it in final $\\\\sigma$ on %d cases.",
    length(case_results),
    proxy_matches_oracle_rho,
    proxy_beats_fixed_rho,
    proxy_beats_fixed_sigma
  ),
  paste(case_sections, collapse = "\n\n"),
  "\\section*{Interactive companion}",
  "The companion HTML gallery is generated by \\texttt{tools/reports/geodesic_mds_paper/render-gmds-misf-phase-e-integrated-html.R}. It includes all corrected top-level seeds, all weighted-KK lookahead layouts, the final integrated pipelines, and the selected stagewise trajectories as interactive \\texttt{rglwidget} panels.",
  "\\end{document}"
)

writeLines(tex_lines, tex_path)

message("Wrote Phase E selection metrics: ", selection_csv)
message("Wrote Phase E final metrics: ", final_csv)
message("Wrote Phase E bundle: ", rds_path)
message("Wrote Phase E report TeX: ", tex_path)
