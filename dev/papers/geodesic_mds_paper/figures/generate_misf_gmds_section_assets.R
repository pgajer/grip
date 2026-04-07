#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (!length(file_arg)) {
  stop("Unable to determine script path from commandArgs().")
}
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
figure_dir <- normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE)
manuscript_dir <- normalizePath(file.path(figure_dir, ".."), winslash = "/", mustWork = TRUE)
interactive_dir <- file.path(manuscript_dir, "interactive")
tmp_dir <- file.path(manuscript_dir, "tmp", "misf-gmds-manuscript-2026-04-02")
dir.create(interactive_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

grip_root <- "/Users/pgajer/current_projects/grip"

options(rgl.useNULL = TRUE)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(grip_root, quiet = TRUE, helpers = FALSE)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(grip_root, quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'pkgload', 'devtools', or the 'grip' package to generate the MISF manuscript assets.")
}

for (pkg in c("htmltools", "htmlwidgets", "rgl")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required to generate the HTML supplement.", pkg))
  }
}

ns <- asNamespace("grip")
align_to_target_nd <- get("grip.align.to.target.nd", envir = ns)
normalize_coords_with_meta <- get("grip.normalize.coords.with.meta", envir = ns)
trace_stage_data <- get("grip.geodesic.misf.trace.stage.data", envir = ns)
trace_stage_lookup <- get("grip.geodesic.misf.trace.stage.lookup", envir = ns)

fmt_num <- function(x, digits = 4L) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}

filter_edges_to_vertices <- function(edges, vertex_ids) {
  edges <- as.matrix(edges)
  vertex_ids <- as.integer(vertex_ids)
  if (!length(vertex_ids) || !nrow(edges)) {
    return(matrix(integer(), ncol = 2L))
  }
  map <- integer(max(max(edges), max(vertex_ids)))
  map[vertex_ids] <- seq_along(vertex_ids)
  keep <- map[edges[, 1L]] > 0L & map[edges[, 2L]] > 0L
  if (!any(keep)) {
    return(matrix(integer(), ncol = 2L))
  }
  cbind(
    as.integer(map[edges[keep, 1L]]),
    as.integer(map[edges[keep, 2L]])
  )
}

project3d <- function(coords, azimuth = 35, elevation = 24) {
  grip::grip.project.3d(coords, azimuth = azimuth, elevation = elevation)
}

align_stage_to_active_reference <- function(stage_coords,
                                            truth_coords,
                                            active_vertices,
                                            allow_reflection = TRUE) {
  active_vertices <- as.integer(active_vertices)
  target_active <- as.matrix(truth_coords[active_vertices, , drop = FALSE])
  target_meta <- normalize_coords_with_meta(target_active)
  target_active_norm <- target_meta$normalized
  target_full_norm <- sweep(as.matrix(truth_coords), 2L, target_meta$center, FUN = "-", check.margin = FALSE) / target_meta$radius

  if (is.null(stage_coords)) {
    out <- matrix(NA_real_, nrow = nrow(truth_coords), ncol = ncol(truth_coords))
    out[active_vertices, ] <- target_active_norm
    return(list(
      aligned_full = out,
      aligned_active = target_active_norm,
      target_active = target_active_norm,
      target_full = target_full_norm,
      rmse = 0
    ))
  }

  source_active <- as.matrix(stage_coords[active_vertices, , drop = FALSE])
  ok <- rowSums(is.finite(source_active)) == ncol(source_active)
  if (!all(ok)) {
    source_active <- source_active[ok, , drop = FALSE]
    target_active_sub <- target_active[ok, , drop = FALSE]
    active_sub <- active_vertices[ok]
  } else {
    target_active_sub <- target_active
    active_sub <- active_vertices
  }
  aligned <- align_to_target_nd(source_active, target_active_sub, allow.reflection = allow_reflection)
  out <- matrix(NA_real_, nrow = nrow(truth_coords), ncol = ncol(truth_coords))
  out[active_sub, ] <- aligned$aligned
  list(
    aligned_full = out,
    aligned_active = aligned$aligned,
    target_active = target_active_norm,
    target_full = target_full_norm,
    rmse = aligned$rmse
  )
}

make_regular_case <- function(side, amplitude = 0.35, seed = side) {
  bundle <- grip::mesh.surface.graph(
    side,
    side,
    surface = "paraboloid",
    amplitude = amplitude,
    connectivity = "orthogonal",
    normalize = "median"
  )
  prepared <- grip::grip.prepare.misf.geodesic.mds(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average",
    dim = 3L,
    top_level_mode = "skip",
    seed = seed
  )
  list(
    id = sprintf("regular_%dx%d", side, side),
    label = sprintf("Regular paraboloid mesh %dx%d", side, side),
    short_label = sprintf("Regular %dx%d", side, side),
    family = "regular",
    side = side,
    bundle = bundle,
    prepared = prepared,
    seed = seed
  )
}

make_irregular_case <- function(side, amplitude = 0.35, seed = side + 100L) {
  bundle <- grip::irregular.rectangle.surface.graph(
    side,
    side,
    surface = "paraboloid",
    amplitude = amplitude,
    connectivity = "orthogonal",
    normalize = "median"
  )
  prepared <- grip::grip.prepare.misf.geodesic.mds(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average",
    dim = 3L,
    top_level_mode = "skip",
    seed = seed
  )
  list(
    id = sprintf("irregular_rectangle_%dx%d", side, side),
    label = sprintf("Irregular rectangle paraboloid %dx%d", side, side),
    short_label = sprintf("Irregular %dx%d", side, side),
    family = "irregular_rectangle",
    side = side,
    bundle = bundle,
    prepared = prepared,
    seed = seed
  )
}

make_sampled_rectangle_cases <- function(n,
                                         k_values = c(6L, 7L, 8L),
                                         xmin = -1.6,
                                         xmax = 1.6,
                                         ymin = -1,
                                         ymax = 1,
                                         amplitude = 0.35,
                                         sample_seed = 1000L + n,
                                         num_init = 6L) {
  seq_spec <- grip::sampled.rectangle.surface.graphs(
    n = n,
    k = as.integer(k_values),
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    seed = sample_seed,
    surface = "paraboloid",
    amplitude = amplitude,
    graph_space = "surface",
    normalize = "median"
  )

  lapply(names(seq_spec$graphs), function(k_name) {
    bundle <- seq_spec$graphs[[k_name]]
    fit_seed <- as.integer(2000L + bundle$n + bundle$k)
    prepared <- grip::grip.prepare.misf.geodesic.mds(
      edges = bundle$edges,
      n = bundle$n,
      edge_weights = bundle$edge_weights,
      tie_mode = "average",
      num_init = num_init,
      dim = 3L,
      top_level_mode = "skip",
      seed = fit_seed
    )
    list(
      id = sprintf("sampled_rectangle_n%d_k%d", n, bundle$k),
      label = sprintf("Sampled rectangle paraboloid n=%d, k=%d", n, bundle$k),
      short_label = sprintf("Sampled n=%d, k=%d", n, bundle$k),
      family = "sampled_rectangle",
      sample_n = n,
      iknn_k = as.integer(bundle$k),
      sample_seed = as.integer(sample_seed),
      num_init = as.integer(num_init),
      bundle = bundle,
      prepared = prepared,
      seed = fit_seed
    )
  })
}

run_case_fit <- function(case) {
  fit <- grip::grip.optimize.misf.geodesic.mds(
    prepared = case$prepared,
    dim = 3L,
    top_level_restarts = 2L,
    top_level_max_iter = 6L,
    top_level_engine = "cpp",
    insertion_anchor_policy = "prev_level_spread",
    insertion_max_iter = 20L,
    refinement_local_nbrs = 4L,
    refinement_landmark_count = 2L,
    refinement_pair_mode = "sparse",
    refinement_anchor_weight = 0.05,
    refinement_anchor_weight_end = 0.01,
    refinement_continuation = "linear",
    refinement_max_iter = 3L,
    refinement_engine = "cpp",
    final_polish_max_iter = 4L,
    final_polish_engine = "cpp",
    n_threads = 0L,
    return_trace = TRUE,
    return_frames = TRUE,
    seed = case$seed
  )
  fit
}

build_stage_rows <- function(case, fit) {
  prepared <- fit$prepared
  truth <- as.matrix(case$bundle$coords_surface)
  stage_data <- trace_stage_data(fit)
  levels <- seq.int(from = prepared$top_level_level, to = 0L, by = -1L)
  rows <- vector("list", length(levels))

  for (i in seq_along(levels)) {
    level <- levels[[i]]
    active_vertices <- as.integer(prepared$misf$levels[[level + 1L]])
    active_edges <- filter_edges_to_vertices(case$bundle$edges, active_vertices)
    ref_align <- align_stage_to_active_reference(NULL, truth, active_vertices)

    geometric_seed <- NULL
    coarse_initial <- NULL

    if (identical(level, prepared$top_level_level)) {
      init_record <- trace_stage_lookup(stage_data, stage = "top_level", level = level)
      init_frame <- if (!is.null(init_record)) init_record$coords_full else fit$frames$after_top_level
      init_label <- if (!is.null(init_record)) init_record$label else "Top-level GMDS solve"
      top_init <- trace_stage_lookup(stage_data, stage = "initial_placement", level = level)
      if (!is.null(top_init) && !is.null(top_init$coords_full)) {
        top_init_align <- align_stage_to_active_reference(
          stage_coords = top_init$coords_full,
          truth_coords = truth,
          active_vertices = active_vertices
        )
        coarse_initial <- list(
          label = top_init$label,
          caption = sprintf("|V_%d| = %d, rho = %s", level, length(active_vertices), fmt_num(top_init_align$rmse, 4L)),
          coords_active = top_init_align$aligned_active,
          coords_full = top_init_align$target_full,
          rmse = top_init_align$rmse
        )

        seed_stage <- trace_stage_lookup(stage_data, stage = "seed", level = level)
        if (!is.null(seed_stage) && !is.null(seed_stage$coords_full) && length(seed_stage$active_vertices) > 0L) {
          seed_vertices <- as.integer(seed_stage$active_vertices)
          seed_align <- align_stage_to_active_reference(
            stage_coords = seed_stage$coords_full,
            truth_coords = truth,
            active_vertices = seed_vertices
          )
          geometric_seed <- list(
            label = seed_stage$label,
            caption = sprintf("|S_%d| = %d, rho = %s", level, length(seed_vertices), fmt_num(seed_align$rmse, 4L)),
            coords_active = seed_align$aligned_active,
            coords_full = seed_align$target_full,
            rmse = seed_align$rmse,
            active_vertices = seed_vertices,
            active_edges = matrix(integer(), ncol = 2L)
          )
        }
      }
    } else {
      init_record <- trace_stage_lookup(stage_data, stage = "insertion", level = level)
      init_frame <- if (!is.null(init_record)) init_record$coords_full else fit$frames$insertion_levels[[paste0("level_", level)]]
      init_label <- if (!is.null(init_record)) init_record$label else sprintf("Insertion of V_%d", level)
    }
    refine_record <- trace_stage_lookup(stage_data, stage = "refinement", level = level)
    refine_frame <- if (!is.null(refine_record)) refine_record$coords_full else fit$frames$refinement_levels[[paste0("level_", level)]]

    init_align <- align_stage_to_active_reference(init_frame, truth, active_vertices)
    refine_align <- align_stage_to_active_reference(refine_frame, truth, active_vertices)

    rows[[i]] <- list(
      level = level,
      active_vertices = active_vertices,
      active_edges = active_edges,
      reference = list(
        label = sprintf("Reference V_%d", level),
        caption = sprintf("|V_%d| = %d", level, length(active_vertices)),
        coords_active = ref_align$target_active,
        coords_full = ref_align$target_full,
        rmse = 0
      ),
      geometric_seed = geometric_seed,
      coarse_initial = coarse_initial,
      initial = list(
        label = init_label,
        caption = sprintf("|V_%d| = %d, rho = %s", level, length(active_vertices), fmt_num(init_align$rmse, 4L)),
        coords_active = init_align$aligned_active,
        coords_full = init_align$target_full,
        rmse = init_align$rmse
      ),
      refinement = list(
        label = if (!is.null(refine_record)) refine_record$label else sprintf("Refinement of V_%d", level),
        caption = sprintf("|V_%d| = %d, rho = %s", level, length(active_vertices), fmt_num(refine_align$rmse, 4L)),
        coords_active = refine_align$aligned_active,
        coords_full = refine_align$target_full,
        rmse = refine_align$rmse
      )
    )
  }

  final_record <- trace_stage_lookup(stage_data, stage = "final_polish", level = 0L)
  final_coords <- if (!is.null(final_record) && !is.null(final_record$coords_full)) {
    final_record$coords_full
  } else {
    fit$coords
  }
  final_align <- align_to_target_nd(final_coords, truth, allow.reflection = TRUE)
  list(
    levels = rows,
    final_polish = list(
      label = if (!is.null(final_record)) final_record$label else "Final full-graph polish",
      caption = sprintf("n = %d, rho = %s", nrow(truth), fmt_num(final_align$rmse, 4L)),
      coords_active = final_align$aligned,
      coords_full = final_align$target,
      edges = as.matrix(case$bundle$edges),
      rmse = final_align$rmse
    )
  )
}

plot_wireframe_panel <- function(stage,
                                 active_edges,
                                 full_edges,
                                 main_title,
                                 sub_title = NULL,
                                 vertex_col = "#24577a",
                                 background_edge_col = "#dddddd",
                                 active_edge_col = "#64748b",
                                 background_point_col = "#e5e7eb",
                                 azimuth = 35,
                                 elevation = 24) {
  xy_bg <- project3d(stage$coords_full, azimuth = azimuth, elevation = elevation)
  xy_active <- project3d(stage$coords_active, azimuth = azimuth, elevation = elevation)
  xlim <- range(c(xy_bg[, 1L], xy_active[, 1L]), finite = TRUE)
  ylim <- range(c(xy_bg[, 2L], xy_active[, 2L]), finite = TRUE)
  xpad <- 0.08 * diff(xlim)
  ypad <- 0.08 * diff(ylim)
  if (!is.finite(xpad) || xpad == 0) xpad <- 0.2
  if (!is.finite(ypad) || ypad == 0) ypad <- 0.2

  graphics::plot(
    xy_bg[, 1L], xy_bg[, 2L],
    type = "n",
    asp = 1,
    axes = FALSE,
    xlab = "",
    ylab = "",
    xlim = xlim + c(-xpad, xpad),
    ylim = ylim + c(-ypad, ypad),
    main = ""
  )
  if (nrow(full_edges) > 0L) {
    apply(full_edges, 1L, function(e) {
      graphics::segments(
        xy_bg[e[1L], 1L], xy_bg[e[1L], 2L],
        xy_bg[e[2L], 1L], xy_bg[e[2L], 2L],
        col = background_edge_col
      )
    })
  }
  graphics::points(xy_bg[, 1L], xy_bg[, 2L], pch = 16, cex = 0.25, col = background_point_col)
  if (nrow(active_edges) > 0L) {
    apply(active_edges, 1L, function(e) {
      graphics::segments(
        xy_active[e[1L], 1L], xy_active[e[1L], 2L],
        xy_active[e[2L], 1L], xy_active[e[2L], 2L],
        col = active_edge_col,
        lwd = 1.2
      )
    })
  }
  graphics::points(xy_active[, 1L], xy_active[, 2L], pch = 16, cex = 0.75, col = vertex_col)
  graphics::mtext(main_title, side = 3L, line = 0.3, cex = 0.82, font = 2L)
  if (!is.null(sub_title) && nzchar(sub_title)) {
    graphics::mtext(sub_title, side = 3L, line = -0.7, cex = 0.72, col = "#4b5563")
  }
}

save_filtration_grid <- function(case_results,
                                 output_pdf,
                                 output_png,
                                 azimuth = 35,
                                 elevation = 24) {
  max_cols <- max(vapply(case_results, function(x) x$fit$prepared$top_level_level + 1L, integer(1L)))

  draw_once <- function() {
    old.par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old.par), add = TRUE)
    graphics::par(
      mfrow = c(length(case_results), max_cols),
      mar = c(0.7, 0.7, 2.7, 0.2),
      oma = c(0, 0, 0.4, 0)
    )

    for (case_idx in seq_along(case_results)) {
      case_result <- case_results[[case_idx]]
      prepared <- case_result$fit$prepared
      coords2d <- as.matrix(case_result$case$bundle$coords_param)
      edges <- as.matrix(case_result$case$bundle$edges)
      levels <- seq.int(from = prepared$top_level_level, to = 0L, by = -1L)
      panel_title <- if (case_result$case$family == "regular") {
        sprintf("(%s) Regular %dx%d mesh graph", LETTERS[case_idx], case_result$case$side, case_result$case$side)
      } else {
        sprintf("(%s) Irregular %dx%d mesh graph", LETTERS[case_idx], case_result$case$side, case_result$case$side)
      }

      blanks <- max_cols - length(levels)
      if (blanks > 0L) {
        for (i in seq_len(blanks)) graphics::plot.new()
      }

      xlim <- range(coords2d[, 1L], finite = TRUE)
      ylim <- range(coords2d[, 2L], finite = TRUE)
      xpad <- 0.08 * diff(xlim)
      ypad <- 0.08 * diff(ylim)
      if (!is.finite(xpad) || xpad == 0) xpad <- 0.2
      if (!is.finite(ypad) || ypad == 0) ypad <- 0.2

      for (level_idx in seq_along(levels)) {
        level <- levels[[level_idx]]
        active <- as.integer(prepared$misf$levels[[level + 1L]])
        graphics::plot(
          coords2d[, 1L], coords2d[, 2L],
          type = "n",
          asp = 1,
          axes = FALSE,
          xlab = "",
          ylab = "",
          xlim = xlim + c(-xpad, xpad),
          ylim = ylim + c(-ypad, ypad),
          main = ""
        )
        if (level_idx == 1L) {
          graphics::mtext(
            panel_title,
            side = 3L,
            line = 1.35,
            cex = 0.88,
            font = 2L,
            adj = 0
          )
        }
        if (nrow(edges) > 0L) {
          apply(edges, 1L, function(e) {
            graphics::segments(
              coords2d[e[1L], 1L], coords2d[e[1L], 2L],
              coords2d[e[2L], 1L], coords2d[e[2L], 2L],
              col = "#d7d7d7"
            )
          })
        }
        graphics::points(coords2d[, 1L], coords2d[, 2L], pch = 16, cex = 0.22, col = "#cbd5e1")
        graphics::points(coords2d[active, 1L], coords2d[active, 2L], pch = 16, cex = 0.75, col = "#9c6644")
        graphics::mtext(
          sprintf("V_%d", level),
          side = 3L,
          line = 0.2,
          cex = 0.86,
          font = 2L
        )
        graphics::mtext(
          sprintf("|V_%d| = %d", level, length(active)),
          side = 3L,
          line = -0.7,
          cex = 0.72,
          col = "#4b5563"
        )
      }
    }
  }

  grDevices::pdf(output_pdf, width = 12, height = 8.8, useDingbats = FALSE)
  draw_once()
  grDevices::dev.off()

  grDevices::png(output_png, width = 2400L, height = 1760L, res = 180, bg = "#ffffff")
  draw_once()
  grDevices::dev.off()
}

save_stage_grid <- function(case_result, output_pdf, output_png) {
  rows <- case_result$stages$levels
  n_rows <- length(rows)

  draw_once <- function() {
    old.par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old.par), add = TRUE)
    graphics::par(
      mfrow = c(n_rows, 3L),
      mar = c(0.7, 0.7, 2.1, 0.2),
      oma = c(0, 0, 1.2, 0)
    )

    full_edges <- as.matrix(case_result$case$bundle$edges)
    for (row in rows) {
      plot_wireframe_panel(
        stage = row$reference,
        active_edges = row$active_edges,
        full_edges = full_edges,
        main_title = row$reference$label,
        sub_title = row$reference$caption,
        vertex_col = "#bc6c25",
        active_edge_col = "#c48a55"
      )
      plot_wireframe_panel(
        stage = row$initial,
        active_edges = row$active_edges,
        full_edges = full_edges,
        main_title = row$initial$label,
        sub_title = row$initial$caption,
        vertex_col = "#24577a",
        active_edge_col = "#4d7ea8"
      )
      plot_wireframe_panel(
        stage = row$refinement,
        active_edges = row$active_edges,
        full_edges = full_edges,
        main_title = row$refinement$label,
        sub_title = row$refinement$caption,
        vertex_col = "#3a5a40",
        active_edge_col = "#588157"
      )
    }
    graphics::mtext(case_result$case$label, side = 3L, outer = TRUE, line = -0.3, cex = 1.08, font = 2L)
  }

  height_in <- max(3.1 * n_rows, 6.5)
  grDevices::pdf(output_pdf, width = 11.5, height = height_in, useDingbats = FALSE)
  draw_once()
  grDevices::dev.off()

  grDevices::png(output_png, width = 2200L, height = as.integer(round(520 * n_rows)), res = 180, bg = "#ffffff")
  draw_once()
  grDevices::dev.off()
}

build_metrics_table <- function(case_results) {
  rows <- list()
  for (case_result in case_results) {
    for (row in case_result$stages$levels) {
      rows[[length(rows) + 1L]] <- data.frame(
        case_id = case_result$case$id,
        case_label = case_result$case$label,
        level = row$level,
        active_n = length(row$active_vertices),
        initial_rho = row$initial$rmse,
        refinement_rho = row$refinement$rmse,
        stringsAsFactors = FALSE
      )
    }
    rows[[length(rows) + 1L]] <- data.frame(
      case_id = case_result$case$id,
      case_label = case_result$case$label,
      level = -1L,
      active_n = case_result$case$bundle$n,
      initial_rho = NA_real_,
      refinement_rho = case_result$stages$final_polish$rmse,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

make_widget <- function(stage,
                        active_edges,
                        full_edges,
                        stage_type = c("reference", "initial", "refinement", "final"),
                        width = 420,
                        height = 320) {
  stage_type <- match.arg(stage_type)
  active_coords <- as.matrix(stage$coords_active)
  full_coords <- as.matrix(stage$coords_full)
  full_edges <- as.matrix(full_edges)
  active_edges <- as.matrix(active_edges)

  rgl::open3d(useNULL = TRUE)
  on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)
  rgl::bg3d(color = "#ffffff")
  rgl::material3d(specular = "#555555")
  rgl::light3d(theta = 35, phi = 20, viewpoint.rel = TRUE)
  rgl::light3d(theta = -55, phi = 25, viewpoint.rel = TRUE)

  if (nrow(full_edges) > 0L) {
    seg_bg <- matrix(NA_real_, nrow = 3L * nrow(full_edges), ncol = 3L)
    seg_bg[seq(1L, nrow(seg_bg), by = 3L), ] <- full_coords[full_edges[, 1L], , drop = FALSE]
    seg_bg[seq(2L, nrow(seg_bg), by = 3L), ] <- full_coords[full_edges[, 2L], , drop = FALSE]
    rgl::lines3d(seg_bg, color = grDevices::adjustcolor("#d7d7d7", alpha.f = 0.9), lwd = 1.0)
  }

  if (nrow(active_edges) > 0L) {
    seg_active <- matrix(NA_real_, nrow = 3L * nrow(active_edges), ncol = 3L)
    seg_active[seq(1L, nrow(seg_active), by = 3L), ] <- active_coords[active_edges[, 1L], , drop = FALSE]
    seg_active[seq(2L, nrow(seg_active), by = 3L), ] <- active_coords[active_edges[, 2L], , drop = FALSE]
    edge_col <- switch(
      stage_type,
      reference = "#c48a55",
      initial = "#4d7ea8",
      refinement = "#588157",
      final = "#7c6a9f"
    )
    rgl::lines3d(seg_active, color = grDevices::adjustcolor(edge_col, alpha.f = 0.92), lwd = 2.0)
  }

  point_col <- switch(
    stage_type,
    reference = "#bc6c25",
    initial = "#24577a",
    refinement = "#3a5a40",
    final = "#6d597a"
  )
  rgl::points3d(active_coords[, 1L], active_coords[, 2L], active_coords[, 3L], color = point_col, size = 7, alpha = 0.98)
  rgl::aspect3d(1, 1, 1)
  rgl::view3d(theta = 35, phi = 24, zoom = 0.82, fov = 30)
  rgl::rglwidget(width = width, height = height)
}

card_div <- function(title_text, widget, note = NULL) {
  lines <- strsplit(title_text, "\n", fixed = TRUE)[[1L]]
  htmltools::tags$div(
    class = "panel-card",
    htmltools::tags$div(
      class = "panel-card-title",
      htmltools::tags$div(class = "panel-card-line panel-card-line-main", lines[[1L]]),
      if (length(lines) > 1L) {
        lapply(lines[-1L], function(x) htmltools::tags$div(class = "panel-card-line panel-card-line-metric", x))
      }
    ),
    htmltools::tags$div(class = "panel-widget", widget),
    if (!is.null(note)) htmltools::tags$p(class = "panel-card-note", note)
  )
}

save_interactive_html <- function(case_results, output_html) {
  full_reference_note <- function(case) {
    if (identical(case$family, "sampled_rectangle")) {
      "Full reference sampled paraboloid graph, normalized for display."
    } else {
      "Full reference paraboloid mesh, normalized for display."
    }
  }

  active_reference_note <- function(case) {
    if (identical(case$family, "sampled_rectangle")) {
      "Reference active-set sampled subgraph used to judge whether the multiscale stage still resembles the target paraboloid."
    } else {
      "Reference active-set sample used to judge whether the multiscale stage still resembles the target paraboloid."
    }
  }

  case_sections <- lapply(case_results, function(case_result) {
    full_edges <- as.matrix(case_result$case$bundle$edges)
    cards <- list(
      card_div(
        title_text = sprintf("Reference full surface\nn = %d", case_result$case$bundle$n),
        widget = make_widget(
          stage = list(coords_active = case_result$stages$final_polish$coords_full, coords_full = case_result$stages$final_polish$coords_full),
          active_edges = full_edges,
          full_edges = full_edges,
          stage_type = "reference"
        ),
        note = full_reference_note(case_result$case)
      ),
      card_div(
        title_text = sprintf("%s\n%s", case_result$stages$final_polish$label, case_result$stages$final_polish$caption),
        widget = make_widget(
          stage = case_result$stages$final_polish,
          active_edges = case_result$stages$final_polish$edges,
          full_edges = full_edges,
          stage_type = "final"
        ),
        note = "Optional full-graph polish after all insertion and levelwise refinement stages."
      )
    )

    for (row in case_result$stages$levels) {
      cards[[length(cards) + 1L]] <- card_div(
        title_text = sprintf("%s\n%s", row$reference$label, row$reference$caption),
        widget = make_widget(
          stage = row$reference,
          active_edges = row$active_edges,
          full_edges = full_edges,
          stage_type = "reference"
        ),
        note = active_reference_note(case_result$case)
      )
      if (!is.null(row$geometric_seed)) {
        cards[[length(cards) + 1L]] <- card_div(
          title_text = sprintf("%s\n%s", row$geometric_seed$label, row$geometric_seed$caption),
          widget = make_widget(
            stage = row$geometric_seed,
            active_edges = row$geometric_seed$active_edges,
            full_edges = full_edges,
            stage_type = "initial"
          ),
          note = "The spread d+1 seed placed first on the coarsest MISF level before inserting the remaining coarse vertices."
        )
      }
      if (!is.null(row$coarse_initial)) {
        cards[[length(cards) + 1L]] <- card_div(
          title_text = sprintf("%s\n%s", row$coarse_initial$label, row$coarse_initial$caption),
          widget = make_widget(
            stage = row$coarse_initial,
            active_edges = row$active_edges,
            full_edges = full_edges,
            stage_type = "initial"
          ),
          note = "All coarsest-level vertices after geometric seed placement and anchor-based insertion, before the pure-GMDS top-level solve."
        )
      }
      cards[[length(cards) + 1L]] <- card_div(
        title_text = sprintf("%s\n%s", row$initial$label, row$initial$caption),
        widget = make_widget(
          stage = row$initial,
          active_edges = row$active_edges,
          full_edges = full_edges,
          stage_type = "initial"
        ),
        note = if (identical(row$level, case_result$fit$prepared$top_level_level)) {
          "Coarsest MISF level after the restartable pure-GMDS solve that starts from the geometric initialization shown above."
        } else {
          "Active-set embedding immediately after geodesic anchor insertion of this level."
        }
      )
      cards[[length(cards) + 1L]] <- card_div(
        title_text = sprintf("%s\n%s", row$refinement$label, row$refinement$caption),
        widget = make_widget(
          stage = row$refinement,
          active_edges = row$active_edges,
          full_edges = full_edges,
          stage_type = "refinement"
        ),
        note = "Sparse GMDS refinement on the active MISF level, with the current scaffold pinned by anchor weights."
      )
    }

    htmltools::tags$section(
      class = "case-section",
      htmltools::tags$h2(case_result$case$label),
        htmltools::tags$p(
          class = "case-lead",
          sprintf(
            "Levels: %s. All non-reference panels are rigidly aligned to the corresponding active-set reference before rendering.",
            paste(
              rev(vapply(case_result$fit$prepared$misf$levels, length, integer(1L))),
              collapse = " \u2192 "
            )
          )
      ),
      htmltools::tags$div(class = "panel-grid", cards)
    )
  })

  page <- htmltools::browsable(
    htmltools::tagList(
      htmltools::tags$head(
        htmltools::tags$title("MISF-GMDS multiscale embeddings"),
        htmltools::tags$style(htmltools::HTML(
          "
          body {
            margin: 0;
            background: #f5f1ea;
            color: #1f2933;
            font-family: 'Iowan Old Style', 'Palatino Linotype', 'Book Antiqua', Palatino, serif;
          }
          .page {
            max-width: 1600px;
            margin: 0 auto;
            padding: 28px 28px 36px;
          }
          .hero {
            background: linear-gradient(135deg, #fffdf9 0%, #f3ece1 100%);
            border: 1px solid #d8c9b7;
            border-radius: 18px;
            padding: 22px 24px;
            box-shadow: 0 10px 28px rgba(58, 36, 14, 0.08);
            margin-bottom: 24px;
          }
          .hero h1 {
            margin: 0 0 10px;
            font-size: 30px;
            line-height: 1.2;
            color: #6b3e26;
          }
          .hero p {
            margin: 10px 0 0;
            font-size: 17px;
            line-height: 1.55;
            max-width: 1180px;
          }
          .case-section {
            margin-top: 26px;
          }
          .case-section h2 {
            margin: 0 0 8px;
            font-size: 28px;
            color: #4a3728;
          }
          .case-lead {
            margin: 0 0 14px;
            color: #5f5445;
            font-size: 15px;
            line-height: 1.5;
          }
          .panel-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(420px, 1fr));
            gap: 18px;
          }
          .panel-card {
            background: #fffdfb;
            border: 1px solid #d9d3c7;
            border-radius: 16px;
            box-shadow: 0 8px 22px rgba(51, 41, 28, 0.08);
            overflow: hidden;
          }
          .panel-card-title {
            padding: 16px 18px 8px;
            border-bottom: 1px solid #ede6da;
            background: linear-gradient(180deg, #fffdfa 0%, #faf5ed 100%);
            min-height: 86px;
          }
          .panel-card-line-main {
            font-size: 20px;
            color: #4a3728;
            font-weight: 700;
            line-height: 1.22;
          }
          .panel-card-line-metric {
            margin-top: 4px;
            font-size: 14px;
            color: #5f5445;
          }
          .panel-widget {
            padding: 8px 10px 4px;
          }
          .panel-card-note {
            margin: 0;
            padding: 0 16px 16px;
            font-size: 13px;
            line-height: 1.45;
            color: #5f5445;
          }
          "
        ))
      ),
      htmltools::tags$div(
        class = "page",
        htmltools::tags$section(
          class = "hero",
          htmltools::tags$h1("MISF-GMDS multiscale embeddings"),
          htmltools::tags$p(
            "This supplement shows the multiscale states used in the MISF-GMDS sections of the GMDS manuscript: the coarsest-level pure-GMDS solve, each finer-level insertion state, each active-level refinement state, and the final full-graph polish. Light-gray wireframes show the normalized reference graph geometry; colored points and edges show the current active embedding."
          )
        ),
        case_sections
      )
    )
  )

  htmltools::save_html(
    page,
    file = output_html,
    libdir = paste0(tools::file_path_sans_ext(output_html), "_files")
  )
}

mesh_cases <- list(
  make_regular_case(12L),
  make_regular_case(15L),
  make_irregular_case(15L)
)

sampled_cases <- c(
  make_sampled_rectangle_cases(50L),
  make_sampled_rectangle_cases(75L)
)

mesh_case_results <- lapply(mesh_cases, function(case) {
  fit <- run_case_fit(case)
  list(
    case = case,
    fit = fit,
    stages = build_stage_rows(case, fit)
  )
})

sampled_case_results <- lapply(sampled_cases, function(case) {
  fit <- run_case_fit(case)
  list(
    case = case,
    fit = fit,
    stages = build_stage_rows(case, fit)
  )
})

case_results <- c(mesh_case_results, sampled_case_results)

filtration_pdf <- file.path(figure_dir, "fig_misf_filtration_cases.pdf")
filtration_png <- file.path(figure_dir, "fig_misf_filtration_cases.png")
save_filtration_grid(mesh_case_results, filtration_pdf, filtration_png)

stage_files <- lapply(case_results, function(case_result) {
  stem <- switch(
    case_result$case$id,
    regular_12x12 = "fig_misf_stages_regular_12x12",
    regular_15x15 = "fig_misf_stages_regular_15x15",
    irregular_rectangle_15x15 = "fig_misf_stages_irregular_rectangle_15x15",
    sampled_rectangle_n50_k6 = "fig_misf_stages_sampled_rectangle_n50_k6",
    sampled_rectangle_n50_k7 = "fig_misf_stages_sampled_rectangle_n50_k7",
    sampled_rectangle_n50_k8 = "fig_misf_stages_sampled_rectangle_n50_k8",
    sampled_rectangle_n75_k6 = "fig_misf_stages_sampled_rectangle_n75_k6",
    sampled_rectangle_n75_k7 = "fig_misf_stages_sampled_rectangle_n75_k7",
    sampled_rectangle_n75_k8 = "fig_misf_stages_sampled_rectangle_n75_k8",
    paste0("fig_", case_result$case$id)
  )
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  save_stage_grid(case_result, pdf_path, png_path)
  list(pdf = pdf_path, png = png_path)
})

metrics_df <- build_metrics_table(case_results)
metrics_csv <- file.path(tmp_dir, "misf_gmds_stage_metrics.csv")
utils::write.csv(metrics_df, metrics_csv, row.names = FALSE)

bundle_rds <- file.path(tmp_dir, "misf_gmds_stage_bundle.rds")
saveRDS(
  list(
    generated_at = Sys.time(),
    case_results = case_results,
    metrics = metrics_df,
    filtration_pdf = filtration_pdf,
    filtration_png = filtration_png,
    stage_files = stage_files
  ),
  bundle_rds
)

interactive_html <- file.path(interactive_dir, "misf_gmds_multiscale_embeddings_2026-04-02.html")
save_interactive_html(case_results, interactive_html)

message("Wrote MIS filtration figure: ", filtration_pdf)
message("Wrote stage figures: ", paste(vapply(stage_files, `[[`, character(1L), "pdf"), collapse = ", "))
message("Wrote stage metrics CSV: ", metrics_csv)
message("Wrote stage bundle: ", bundle_rds)
message("Wrote interactive HTML: ", interactive_html)
