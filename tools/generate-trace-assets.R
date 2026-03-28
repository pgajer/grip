#!/usr/bin/env Rscript

fig_dir <- file.path("dev", "manual", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'devtools' or the 'grip' package to generate trace assets.")
}

edges_sierpinski_carpet <- function(level) {
  stopifnot(level >= 0)
  side <- 3L^level
  grid <- expand.grid(x = 0:(side - 1L), y = 0:(side - 1L))

  keep_cell <- function(x, y) {
    while (x > 0L || y > 0L) {
      if ((x %% 3L) == 1L && (y %% 3L) == 1L) {
        return(FALSE)
      }
      x <- x %/% 3L
      y <- y %/% 3L
    }
    TRUE
  }

  keep <- mapply(keep_cell, grid$x, grid$y)
  cells <- grid[keep, , drop = FALSE]

  id_map <- matrix(0L, nrow = side, ncol = side)
  for (i in seq_len(nrow(cells))) {
    id_map[cells$x[i] + 1L, cells$y[i] + 1L] <- i
  }

  edges <- list()
  for (i in seq_len(nrow(cells))) {
    x <- cells$x[i]
    y <- cells$y[i]
    if (x + 1L < side) {
      nbr <- id_map[x + 2L, y + 1L]
      if (nbr > 0L) edges[[length(edges) + 1L]] <- c(i, nbr)
    }
    if (y + 1L < side) {
      nbr <- id_map[x + 1L, y + 2L]
      if (nbr > 0L) edges[[length(edges) + 1L]] <- c(i, nbr)
    }
  }

  do.call(rbind, edges)
}

edges_sierpinski_triangle <- function(level) {
  stopifnot(level >= 0)

  merge_nodes <- function(edges, from, to) {
    edges[edges == from] <- to
    edges
  }

  build <- function(k) {
    if (k == 0) {
      edges <- rbind(c(1L, 2L), c(2L, 3L), c(3L, 1L))
      return(list(edges = edges, corners = c(1L, 2L, 3L), n = 3L))
    }

    left <- build(k - 1)
    right <- build(k - 1)
    top <- build(k - 1)

    off1 <- left$n
    off2 <- left$n + right$n

    edges <- rbind(left$edges,
                   right$edges + off1,
                   top$edges + off2)

    L <- left$corners
    R <- right$corners + off1
    T <- top$corners + off2

    edges <- merge_nodes(edges, R[1], L[2])
    edges <- merge_nodes(edges, T[1], L[3])
    edges <- merge_nodes(edges, T[2], R[3])

    ids <- sort(unique(c(edges)))
    map <- seq_along(ids)
    names(map) <- ids
    edges <- cbind(map[as.character(edges[, 1])],
                   map[as.character(edges[, 2])])
    edges <- t(apply(edges, 1, sort))
    edges <- unique(edges)

    corners <- c(map[as.character(L[1])],
                 map[as.character(R[2])],
                 map[as.character(T[3])])

    list(edges = edges, corners = corners, n = length(ids))
  }

  build(level)$edges
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
                 -sy, 0, cy), nrow = 3, byrow = TRUE)
  rx <- matrix(c(1, 0, 0,
                 0, cp, -sp,
                 0, sp, cp), nrow = 3, byrow = TRUE)
  rz <- matrix(c(cr, -sr, 0,
                 sr, cr, 0,
                 0, 0, 1), nrow = 3, byrow = TRUE)

  coords %*% ry %*% rx %*% rz
}

center_trace_frame <- function(coords) {
  active <- stats::complete.cases(coords)
  if (!any(active)) {
    return(coords)
  }
  centered <- coords
  centered[active, ] <- sweep(coords[active, , drop = FALSE],
                              2L,
                              colMeans(coords[active, , drop = FALSE]),
                              "-")
  centered
}

normalize_trace_frames <- function(frames) {
  centered <- lapply(frames, center_trace_frame)
  radius <- max(vapply(centered, function(z) {
    active <- stats::complete.cases(z)
    if (!any(active)) return(0)
    max(sqrt(rowSums(z[active, , drop = FALSE]^2)))
  }, numeric(1L)))
  if (!is.finite(radius) || radius <= 0) {
    return(centered)
  }
  lapply(centered, function(z) {
    active <- stats::complete.cases(z)
    if (any(active)) {
      z[active, ] <- z[active, , drop = FALSE] / radius
    }
    z
  })
}

trace_subtitle <- function(meta_row) {
  phase <- switch(meta_row$phase,
                  init = "coarsest initialization",
                  level_start = "new level introduced",
                  round = "force refinement",
                  final = "final layout",
                  meta_row$phase)
  sprintf("%s | level %d | round %d | active %d",
          phase,
          meta_row$level_index,
          meta_row$round_in_level,
          meta_row$active_vertices)
}

compute_2d_limits <- function(frames) {
  xs <- c()
  ys <- c()
  for (coords in frames) {
    active <- stats::complete.cases(coords)
    if (!any(active)) next
    xs <- c(xs, coords[active, 1])
    ys <- c(ys, coords[active, 2])
  }
  xr <- range(xs)
  yr <- range(ys)
  xpad <- max(diff(xr) * 0.12, 0.12)
  ypad <- max(diff(yr) * 0.12, 0.12)
  list(xlim = xr + c(-xpad, xpad), ylim = yr + c(-ypad, ypad))
}

edge_mask <- function(coords, edges) {
  active <- stats::complete.cases(coords)
  active[edges[, 1]] & active[edges[, 2]]
}

draw_2d_trace_frame <- function(path, coords, edges, title_text, subtitle_text,
                                xlim, ylim, width = 900, height = 900) {
  n <- nrow(coords)
  node_cex <- if (n <= 80) {
    1.5
  } else if (n <= 180) {
    1.0
  } else if (n <= 600) {
    0.56
  } else if (n <= 1500) {
    0.30
  } else {
    0.12
  }
  edge_lwd <- if (n <= 80) {
    2.0
  } else if (n <= 180) {
    1.35
  } else if (n <= 600) {
    0.75
  } else if (n <= 1500) {
    0.38
  } else {
    0.16
  }

  keep_edges <- edge_mask(coords, edges)
  active <- stats::complete.cases(coords)

  png(path, width = width, height = height, res = 160, bg = "#f7f3ea")
  on.exit(dev.off(), add = TRUE)
  par(mar = c(0, 0, 3.4, 0), xaxs = "i", yaxs = "i")
  plot(NA,
       xlim = xlim,
       ylim = ylim,
       asp = 1,
       axes = FALSE,
       xlab = "",
       ylab = "")
  if (any(keep_edges)) {
    apply(edges[keep_edges, , drop = FALSE], 1, function(e) {
      segments(coords[e[1], 1], coords[e[1], 2],
               coords[e[2], 1], coords[e[2], 2],
               col = grDevices::adjustcolor("#0f3b5f", alpha.f = 0.18),
               lwd = edge_lwd)
    })
  }
  if (any(active)) {
    points(coords[active, 1], coords[active, 2],
           pch = 21,
           bg = "#f05a28",
           col = "#16324f",
           cex = node_cex,
           lwd = 0.7)
  }
  title(main = title_text, sub = subtitle_text,
        col.main = "#16324f", cex.main = 1.28,
        col.sub = "#466074", cex.sub = 0.85)
}

plot_2d_trace_gif <- function(path, trace_obj, edges, title_text,
                              width = 900, height = 900, fps = 10,
                              final_hold_sec = 0) {
  if (!requireNamespace("magick", quietly = TRUE)) {
    message("Skipping GIF generation: package 'magick' is not installed.")
    return(invisible(FALSE))
  }

  frames <- normalize_trace_frames(trace_obj$frames)
  limits <- compute_2d_limits(frames)
  frame_dir <- tempfile("grip-trace-2d-")
  dir.create(frame_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(frame_dir, recursive = TRUE, force = TRUE), add = TRUE)

  frame_paths <- file.path(frame_dir, sprintf("frame-%03d.png", seq_along(frames)))
  for (i in seq_along(frames)) {
    draw_2d_trace_frame(frame_paths[i],
                        frames[[i]],
                        edges,
                        title_text = title_text,
                        subtitle_text = trace_subtitle(trace_obj$meta[i, , drop = FALSE]),
                        xlim = limits$xlim,
                        ylim = limits$ylim,
                        width = width,
                        height = height)
  }

  if (final_hold_sec > 0) {
    hold_frames <- max(1L, as.integer(round(fps * final_hold_sec)))
    frame_paths <- c(frame_paths, rep(tail(frame_paths, 1L), hold_frames))
  }

  anim <- magick::image_animate(magick::image_read(frame_paths),
                                fps = fps,
                                optimize = TRUE)
  magick::image_write(anim, path)
  invisible(TRUE)
}

project_frame_3d <- function(coords, yaw, pitch, camera = 4.2) {
  projected <- matrix(NA_real_, nrow = nrow(coords), ncol = 3)
  colnames(projected) <- c("x", "y", "z")
  active <- stats::complete.cases(coords)
  if (!any(active)) {
    return(projected)
  }
  rotated <- rotate_xyz(coords[active, , drop = FALSE], yaw = yaw, pitch = pitch)
  z_camera <- rotated[, 3] + camera
  scale <- camera / pmax(z_camera, 1e-6)
  projected[active, 1] <- rotated[, 1] * scale
  projected[active, 2] <- rotated[, 2] * scale
  projected[active, 3] <- rotated[, 3]
  projected
}

compute_3d_limits <- function(frames, yaw, pitch, camera = 4.2) {
  xs <- c()
  ys <- c()
  for (coords in frames) {
    projected <- project_frame_3d(coords, yaw = yaw, pitch = pitch, camera = camera)
    active <- stats::complete.cases(projected)
    if (!any(active)) next
    xs <- c(xs, projected[active, 1])
    ys <- c(ys, projected[active, 2])
  }
  xr <- range(xs)
  yr <- range(ys)
  xpad <- max(diff(xr) * 0.14, 0.12)
  ypad <- max(diff(yr) * 0.14, 0.12)
  list(xlim = xr + c(-xpad, xpad), ylim = yr + c(-ypad, ypad))
}

draw_3d_trace_frame <- function(path, coords, edges, yaw, pitch,
                                title_text, subtitle_text,
                                xlim, ylim,
                                width = 900, height = 900,
                                camera = 4.2) {
  projected <- project_frame_3d(coords, yaw = yaw, pitch = pitch, camera = camera)
  active <- stats::complete.cases(projected)
  keep_edges <- active[edges[, 1]] & active[edges[, 2]]
  n <- nrow(coords)
  edge_lwd <- if (n <= 100) 2.0 else if (n <= 600) 1.2 else if (n <= 1500) 0.55 else 0.24
  base_cex <- if (n <= 100) 0.9 else if (n <= 600) 0.42 else if (n <= 1500) 0.22 else 0.12

  png(path, width = width, height = height, res = 160, bg = "#f7f3ea")
  on.exit(dev.off(), add = TRUE)
  par(mar = c(0, 0, 3.4, 0), xaxs = "i", yaxs = "i")
  plot(NA,
       xlim = xlim,
       ylim = ylim,
       asp = 1,
       axes = FALSE,
       xlab = "",
       ylab = "")

  if (any(keep_edges)) {
    edge_subset <- edges[keep_edges, , drop = FALSE]
    edge_depth <- rowMeans(matrix(projected[as.vector(t(edge_subset)), 3], ncol = 2))
    edge_order <- order(edge_depth)
    for (idx in edge_order) {
      e <- edge_subset[idx, ]
      alpha <- 0.10 + 0.36 * ((edge_depth[idx] - min(edge_depth)) /
        max(diff(range(edge_depth)), 1e-6))
      segments(projected[e[1], 1], projected[e[1], 2],
               projected[e[2], 1], projected[e[2], 2],
               col = grDevices::adjustcolor("#1f4d6b", alpha.f = alpha),
               lwd = edge_lwd)
    }
  }

  if (any(active)) {
    point_order <- order(projected[active, 3])
    z_vals <- projected[active, 3]
    z_range <- range(z_vals)
    z_span <- max(diff(z_range), 1e-6)
    point_level <- (z_vals - z_range[1]) / z_span
    palette <- grDevices::colorRampPalette(c("#f7b267", "#f4845f", "#d1495b", "#7b2d26"))(100)
    cols <- palette[pmax(1L, pmin(100L, floor(point_level * 99) + 1L))]
    cex_vals <- base_cex + base_cex * 1.2 * point_level
    xy <- projected[active, 1:2, drop = FALSE]
    points(xy[point_order, 1], xy[point_order, 2],
           pch = 21,
           bg = cols[point_order],
           col = "#16324f",
           cex = cex_vals[point_order],
           lwd = 0.55)
  }

  title(main = title_text, sub = subtitle_text,
        col.main = "#16324f", cex.main = 1.28,
        col.sub = "#466074", cex.sub = 0.85)
}

plot_3d_trace_gif <- function(path, trace_obj, edges, title_text,
                              width = 900, height = 900, fps = 10,
                              final_hold_sec = 0,
                              yaw = 0.9, pitch = 0.68, camera = 4.2) {
  if (!requireNamespace("magick", quietly = TRUE)) {
    message("Skipping GIF generation: package 'magick' is not installed.")
    return(invisible(FALSE))
  }

  frames <- normalize_trace_frames(trace_obj$frames)
  limits <- compute_3d_limits(frames, yaw = yaw, pitch = pitch, camera = camera)
  frame_dir <- tempfile("grip-trace-3d-")
  dir.create(frame_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(frame_dir, recursive = TRUE, force = TRUE), add = TRUE)

  frame_paths <- file.path(frame_dir, sprintf("frame-%03d.png", seq_along(frames)))
  for (i in seq_along(frames)) {
    draw_3d_trace_frame(frame_paths[i],
                        frames[[i]],
                        edges,
                        yaw = yaw,
                        pitch = pitch,
                        title_text = title_text,
                        subtitle_text = trace_subtitle(trace_obj$meta[i, , drop = FALSE]),
                        xlim = limits$xlim,
                        ylim = limits$ylim,
                        width = width,
                        height = height,
                        camera = camera)
  }

  if (final_hold_sec > 0) {
    hold_frames <- max(1L, as.integer(round(fps * final_hold_sec)))
    frame_paths <- c(frame_paths, rep(tail(frame_paths, 1L), hold_frames))
  }

  anim <- magick::image_animate(magick::image_read(frame_paths),
                                fps = fps,
                                optimize = TRUE)
  magick::image_write(anim, path)
  invisible(TRUE)
}

message("Computing trace animations...")

carpet4 <- edges_sierpinski_carpet(4)
triangle4 <- edges_sierpinski_triangle(4)

# Keep the animation scripts on the primary layout defaults and only
# downsample the saved trace frames for a manageable GIF size.
carpet4_trace_2d <- grip.layout.trace(carpet4,
                                      n = max(carpet4),
                                      dim = 2,
                                      trace = "round",
                                      trace.every = 8,
                                      seed = 24)

carpet4_trace_3d <- grip.layout.trace(carpet4,
                                      n = max(carpet4),
                                      dim = 3,
                                      trace = "round",
                                      trace.every = 12,
                                      seed = 24)

triangle4_trace_2d <- grip.layout.trace(triangle4,
                                        n = max(triangle4),
                                        dim = 2,
                                        placement = "circle",
                                        trace = "round",
                                        trace.every = 6,
                                        seed = 25)

triangle4_trace_3d <- grip.layout.trace(triangle4,
                                        n = max(triangle4),
                                        dim = 3,
                                        trace = "round",
                                        trace.every = 8,
                                        seed = 25)

plot_2d_trace_gif(file.path(fig_dir, "trace-sierpinski-carpet-level-4-2d.gif"),
                  carpet4_trace_2d,
                  carpet4,
                  title_text = "Sierpinski Carpet (Level 4, 2D Trace)",
                  fps = 5,
                  final_hold_sec = 5)
plot_3d_trace_gif(file.path(fig_dir, "trace-sierpinski-carpet-level-4-3d.gif"),
                  carpet4_trace_3d,
                  carpet4,
                  title_text = "Sierpinski Carpet (Level 4, 3D Trace)")
plot_2d_trace_gif(file.path(fig_dir, "trace-sierpinski-triangle-level-4-2d.gif"),
                  triangle4_trace_2d,
                  triangle4,
                  title_text = "Sierpinski Triangle (Level 4, 2D Trace)",
                  fps = 5,
                  final_hold_sec = 5)
plot_3d_trace_gif(file.path(fig_dir, "trace-sierpinski-triangle-level-4-3d.gif"),
                  triangle4_trace_3d,
                  triangle4,
                  title_text = "Sierpinski Triangle (Level 4, 3D Trace)")

message(sprintf("Trace assets written to %s", fig_dir))
