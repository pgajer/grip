#!/usr/bin/env Rscript

source("tools/benchmark-sierpinski-baseline.R")

args <- commandArgs(trailingOnly = TRUE)
arg_map <- list()
for (arg in args) {
  if (!grepl("=", arg, fixed = TRUE)) {
    next
  }
  parts <- strsplit(arg, "=", fixed = TRUE)[[1L]]
  key <- parts[[1L]]
  value <- paste(parts[-1L], collapse = "=")
  arg_map[[key]] <- value
}

arg_value <- function(name, default = NULL) {
  if (name %in% names(arg_map)) {
    arg_map[[name]]
  } else {
    default
  }
}

parse_levels <- function(text, default) {
  if (is.null(text) || !nzchar(text)) {
    return(as.integer(default))
  }
  as.integer(strsplit(text, ",", fixed = TRUE)[[1L]])
}

family <- arg_value("family", "all")
seed <- as.integer(arg_value("seed", "1"))
triangle_levels <- parse_levels(arg_value("triangle_levels", NULL), 3:7)
tetrahedron_levels <- parse_levels(arg_value("tetrahedron_levels", NULL), 3:6)
manual_root <- file.path("dev", "manual")
pdf_out_dir <- file.path(manual_root, "pdf", "diagnostics")
preview_dir <- file.path(manual_root, "tmp", "diagnostics-previews")
dir.create(pdf_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)

normalize_undirected_edges_local <- function(edges) {
  if (is.null(edges) || length(edges) == 0L || nrow(edges) == 0L) {
    return(matrix(integer(), ncol = 2L))
  }
  edges <- cbind(pmin(edges[, 1L], edges[, 2L]),
                 pmax(edges[, 1L], edges[, 2L]))
  edges <- edges[edges[, 1L] != edges[, 2L], , drop = FALSE]
  storage.mode(edges) <- "integer"
  unique(edges)
}

align_to_target_nd <- function(source, target) {
  src <- normalize_coords(source)
  dst <- normalize_coords(target)
  cross <- t(src) %*% dst
  sv <- svd(cross)
  rot <- sv$u %*% t(sv$v)
  if (det(rot) < 0) {
    fix <- diag(ncol(rot))
    fix[ncol(fix), ncol(fix)] <- -1
    rot <- sv$u %*% fix %*% t(sv$v)
  }
  src %*% rot
}

rotate_xyz <- function(coords, yaw = 0, pitch = 0, roll = 0) {
  cy <- cos(yaw)
  sy <- sin(yaw)
  cp <- cos(pitch)
  sp <- sin(pitch)
  cr <- cos(roll)
  sr <- sin(roll)

  ry <- matrix(c(cy, 0, sy,
                 0, 1, 0,
                 -sy, 0, cy), nrow = 3L, byrow = TRUE)
  rx <- matrix(c(1, 0, 0,
                 0, cp, -sp,
                 0, sp, cp), nrow = 3L, byrow = TRUE)
  rz <- matrix(c(cr, -sr, 0,
                 sr, cr, 0,
                 0, 0, 1), nrow = 3L, byrow = TRUE)

  coords %*% ry %*% rx %*% rz
}

project_perspective <- function(coords, camera = 4.6) {
  z_camera <- coords[, 3L] + camera
  scale <- camera / pmax(z_camera, 1e-6)
  cbind(x = coords[, 1L] * scale, y = coords[, 2L] * scale, z = coords[, 3L])
}

format_profile <- function(profile) {
  sprintf(
    "%s %d/%d n0=%d k=%d r=%.2f s=%.1f rep=%.1f",
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

default_profile <- function() {
  list(
    placement = "barycenter",
    rounds = 20L,
    final_rounds = 25L,
    num_init = 36L,
    num_nbrs = 10L,
    r = 0.15,
    s = 3.0,
    repulsion_factor = 1.0
  )
}

carpet_preset_profile <- function() {
  get("grip.carpet.preset.defaults", envir = asNamespace("grip"))()
}

build_sierpinski_tetrahedron <- function(level) {
  stopifnot(length(level) == 1L, is.finite(level), level == as.integer(level), level >= 0L)
  level <- as.integer(level)

  coords_list <- vector("list", 4L)
  coords_list[[1L]] <- c(0, 0, 0)
  coords_list[[2L]] <- c(1, 0, 0)
  coords_list[[3L]] <- c(0.5, sqrt(3) / 2, 0)
  coords_list[[4L]] <- c(0.5, sqrt(3) / 6, sqrt(2 / 3))

  state <- new.env(parent = emptyenv())
  state$edges <- list()
  state$next_id <- 5L

  add_tetrahedron <- function(a, b, c, d) {
    state$edges[[length(state$edges) + 1L]] <- rbind(
      c(a, b), c(a, c), c(a, d),
      c(b, c), c(b, d), c(c, d)
    )
  }

  midpoint_id <- function(u, v) {
    id <- state$next_id
    state$next_id <- state$next_id + 1L
    coords_list[[id]] <<- (coords_list[[u]] + coords_list[[v]]) / 2
    id
  }

  recurse <- function(current_level, a, b, c, d) {
    if (current_level >= level) {
      add_tetrahedron(a, b, c, d)
      return(invisible(NULL))
    }

    e <- midpoint_id(a, b)
    f <- midpoint_id(a, c)
    g <- midpoint_id(a, d)
    h <- midpoint_id(b, c)
    i <- midpoint_id(b, d)
    j <- midpoint_id(c, d)

    recurse(current_level + 1L, a, g, e, f)
    recurse(current_level + 1L, e, b, i, h)
    recurse(current_level + 1L, f, c, j, h)
    recurse(current_level + 1L, g, d, i, j)
  }

  recurse(0L, 1L, 2L, 3L, 4L)

  edges <- normalize_undirected_edges_local(do.call(rbind, state$edges))
  coords <- do.call(rbind, coords_list[seq_len(state$next_id - 1L)])
  storage.mode(coords) <- "double"

  package_edges <- edges.sierpinski.tetrahedron(level)
  if (!identical(package_edges, edges)) {
    stop(sprintf("Canonical tetrahedron builder does not match package edges at level %d", level))
  }

  list(edges = edges, coords = coords)
}

plot_2d_triptych <- function(path,
                             canonical_coords,
                             baseline_coords,
                             tuned_coords,
                             edges,
                             title_text,
                             subtitle_text,
                             tuned_label) {
  can_norm <- normalize_coords(canonical_coords)
  base_norm <- normalize_coords(baseline_coords)
  tuned_norm <- normalize_coords(tuned_coords)
  xs <- c(can_norm[, 1L], base_norm[, 1L], tuned_norm[, 1L])
  ys <- c(can_norm[, 2L], base_norm[, 2L], tuned_norm[, 2L])
  xr <- range(xs)
  yr <- range(ys)
  xpad <- max(diff(xr) * 0.12, 0.12)
  ypad <- max(diff(yr) * 0.12, 0.12)
  xlim <- xr + c(-xpad, xpad)
  ylim <- yr + c(-ypad, ypad)

  grDevices::pdf(path, width = 21, height = 8.5,
                 paper = "special", bg = "#f7f3ea", useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, 3), mar = c(0, 0, 2.8, 0), xaxs = "i", yaxs = "i")
  plot_layout_panel(canonical_coords, edges, paste(title_text, "- canonical"), subtitle_text,
                    xlim = xlim, ylim = ylim)
  plot_layout_panel(baseline_coords, edges, paste(title_text, "- baseline"), subtitle_text,
                    xlim = xlim, ylim = ylim)
  plot_layout_panel(tuned_coords, edges, paste(title_text, "-", tuned_label), subtitle_text,
                    xlim = xlim, ylim = ylim)
}

plot_3d_panel <- function(projected_coords, edges, title_text, subtitle_text,
                          xlim, ylim) {
  n <- nrow(projected_coords)
  node_cex <- if (n <= 200L) {
    0.95
  } else if (n <= 2500L) {
    0.42
  } else {
    0.12
  }
  edge_lwd <- if (n <= 200L) {
    1.1
  } else if (n <= 2500L) {
    0.45
  } else {
    0.12
  }

  plot(NA,
       xlim = xlim,
       ylim = ylim,
       asp = 1,
       axes = FALSE,
       xlab = "",
       ylab = "")

  edge_depth <- rowMeans(cbind(projected_coords[edges[, 1L], 3L],
                               projected_coords[edges[, 2L], 3L]))
  edge_order <- order(edge_depth)
  for (idx in edge_order) {
    e <- edges[idx, ]
    depth_alpha <- 0.12 + 0.20 * (projected_coords[e[1L], 3L] + projected_coords[e[2L], 3L] -
                                  min(projected_coords[, 3L])) /
      max(diff(range(projected_coords[, 3L])), 1e-6)
    graphics::segments(projected_coords[e[1L], 1L], projected_coords[e[1L], 2L],
                       projected_coords[e[2L], 1L], projected_coords[e[2L], 2L],
                       col = grDevices::adjustcolor("#0f3b5f", alpha.f = min(depth_alpha, 0.45)),
                       lwd = edge_lwd)
  }

  point_order <- order(projected_coords[, 3L])
  for (idx in point_order) {
    depth_scale <- 0.70 + 0.55 * (projected_coords[idx, 3L] - min(projected_coords[, 3L])) /
      max(diff(range(projected_coords[, 3L])), 1e-6)
    graphics::points(projected_coords[idx, 1L], projected_coords[idx, 2L],
                     pch = 21,
                     bg = "#f05a28",
                     col = "#16324f",
                     cex = node_cex * depth_scale,
                     lwd = 0.45)
  }

  graphics::title(main = title_text, sub = subtitle_text,
                  col.main = "#16324f", cex.main = 1.05,
                  col.sub = "#466074", cex.sub = 0.82)
}

write_3d_diagnostic_pdf <- function(path,
                                    canonical_coords,
                                    baseline_coords,
                                    tuned_coords,
                                    edges,
                                    title_text,
                                    subtitle_text,
                                    tuned_label) {
  views <- list(
    list(name = "view A", yaw = 0.82, pitch = 0.55, roll = 0.00, camera = 4.8),
    list(name = "view B", yaw = -0.95, pitch = 0.28, roll = 0.18, camera = 4.8)
  )

  grDevices::pdf(path, width = 21, height = 13,
                 paper = "special", bg = "#f7f3ea", useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(length(views), 3), mar = c(0, 0, 2.4, 0), xaxs = "i", yaxs = "i")

  for (view in views) {
    can_view <- project_perspective(
      rotate_xyz(normalize_coords(canonical_coords), yaw = view$yaw, pitch = view$pitch, roll = view$roll),
      camera = view$camera
    )
    base_view <- project_perspective(
      rotate_xyz(normalize_coords(baseline_coords), yaw = view$yaw, pitch = view$pitch, roll = view$roll),
      camera = view$camera
    )
    tuned_view <- project_perspective(
      rotate_xyz(normalize_coords(tuned_coords), yaw = view$yaw, pitch = view$pitch, roll = view$roll),
      camera = view$camera
    )

    xs <- c(can_view[, 1L], base_view[, 1L], tuned_view[, 1L])
    ys <- c(can_view[, 2L], base_view[, 2L], tuned_view[, 2L])
    xr <- range(xs)
    yr <- range(ys)
    xpad <- max(diff(xr) * 0.14, 0.12)
    ypad <- max(diff(yr) * 0.14, 0.12)
    xlim <- xr + c(-xpad, xpad)
    ylim <- yr + c(-ypad, ypad)

    plot_3d_panel(can_view, edges,
                  paste(title_text, "- canonical", sprintf("(%s)", view$name)),
                  subtitle_text, xlim, ylim)
    plot_3d_panel(base_view, edges,
                  paste(title_text, "- baseline", sprintf("(%s)", view$name)),
                  subtitle_text, xlim, ylim)
    plot_3d_panel(tuned_view, edges,
                  paste(title_text, "-", tuned_label, sprintf("(%s)", view$name)),
                  subtitle_text, xlim, ylim)
  }
}

run_triangle_diagnostics <- function(levels, seed) {
  cfg_baseline <- baseline_profile$triangle
  pdfs <- character()
  for (level in levels) {
    message(sprintf("Building level %d triangle canonical graph...", level))
    built <- build_sierpinski_triangle(level)
    edges <- built$edges
    n <- max(edges)
    message(sprintf("Level %d: %d vertices, %d edges", level, n, nrow(edges)))

    baseline_coords <- grip.layout(
      edges = edges,
      n = n,
      dim = 2,
      placement = cfg_baseline$placement,
      rounds = cfg_baseline$rounds,
      final_rounds = cfg_baseline$final_rounds,
      num_init = cfg_baseline$num_init,
      num_nbrs = cfg_baseline$num_nbrs,
      r = cfg_baseline$r,
      s = cfg_baseline$s,
      repulsion_factor = cfg_baseline$repulsion_factor,
      seed = seed
    )
    tuned_coords <- grip.layout(
      edges = edges,
      n = n,
      dim = 2,
      preset = "carpet",
      seed = seed
    )

    baseline_aligned <- align_to_target_nd(baseline_coords, built$coords)
    tuned_aligned <- align_to_target_nd(tuned_coords, built$coords)
    pdf_path <- file.path(
      pdf_out_dir,
      sprintf("sierpinski-triangle-level-%d-canonical-baseline-carpet-preset.pdf", level)
    )
    subtitle <- paste(
      sprintf("seed=%d", seed),
      "| baseline", format_profile(cfg_baseline),
      "| preset=carpet", format_profile(carpet_preset_profile())
    )

    message(sprintf("Writing %s", pdf_path))
    plot_2d_triptych(
      path = pdf_path,
      canonical_coords = built$coords,
      baseline_coords = baseline_aligned,
      tuned_coords = tuned_aligned,
      edges = edges,
      title_text = sprintf("Sierpinski triangle level %d", level),
      subtitle_text = subtitle,
      tuned_label = "preset=carpet"
    )
    render_pdf_previews(pdf_path, preview_dir)
    pdfs <- c(pdfs, pdf_path)
    message(sprintf("Finished level %d", level))
  }
  invisible(pdfs)
}

run_tetrahedron_diagnostics <- function(levels, seed) {
  cfg_baseline <- default_profile()
  cfg_tuned <- carpet_preset_profile()
  pdfs <- character()
  for (level in levels) {
    message(sprintf("Building level %d tetrahedron canonical graph...", level))
    built <- build_sierpinski_tetrahedron(level)
    edges <- built$edges
    n <- max(edges)
    message(sprintf("Level %d: %d vertices, %d edges", level, n, nrow(edges)))

    baseline_coords <- grip.layout(
      edges = edges,
      n = n,
      dim = 3,
      seed = seed
    )
    tuned_coords <- grip.layout(
      edges = edges,
      n = n,
      dim = 3,
      preset = "carpet",
      seed = seed
    )

    baseline_aligned <- align_to_target_nd(baseline_coords, built$coords)
    tuned_aligned <- align_to_target_nd(tuned_coords, built$coords)

    pdf_path <- file.path(
      pdf_out_dir,
      sprintf("sierpinski-tetrahedron-level-%d-canonical-default-carpet-preset-3d.pdf", level)
    )
    subtitle <- paste(
      sprintf("seed=%d", seed),
      "| defaults", format_profile(cfg_baseline),
      "| preset=carpet", format_profile(cfg_tuned)
    )

    message(sprintf("Writing %s", pdf_path))
    write_3d_diagnostic_pdf(
      path = pdf_path,
      canonical_coords = built$coords,
      baseline_coords = baseline_aligned,
      tuned_coords = tuned_aligned,
      edges = edges,
      title_text = sprintf("Sierpinski tetrahedron level %d", level),
      subtitle_text = subtitle,
      tuned_label = "preset=carpet"
    )
    render_pdf_previews(pdf_path, preview_dir)
    pdfs <- c(pdfs, pdf_path)
    message(sprintf("Finished level %d", level))
  }
  invisible(pdfs)
}

pdfs <- character()
if (family %in% c("triangle", "all")) {
  pdfs <- c(pdfs, run_triangle_diagnostics(triangle_levels, seed))
}
if (family %in% c("tetrahedron", "all")) {
  pdfs <- c(pdfs, run_tetrahedron_diagnostics(tetrahedron_levels, seed))
}

if (length(pdfs) == 0L) {
  stop("No diagnostics requested. Use family=triangle, family=tetrahedron, or family=all.")
}

cat(paste(pdfs, collapse = "\n"), "\n")
