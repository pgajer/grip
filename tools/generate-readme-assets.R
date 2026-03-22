#!/usr/bin/env Rscript

dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'devtools' or the 'grip' package to generate README assets.")
}

edges_cycle <- function(n) {
  if (n < 2) {
    return(matrix(integer(), ncol = 2))
  }
  rbind(cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L), c(n, 1L))
}

edges_mesh <- function(h, w = h) {
  stopifnot(h >= 1, w >= 1)
  idx <- function(i, j) (i - 1L) * w + j
  edges <- list()
  for (i in seq_len(h)) {
    for (j in seq_len(w)) {
      v <- idx(i, j)
      if (i < h) edges[[length(edges) + 1L]] <- c(v, idx(i + 1L, j))
      if (j < w) edges[[length(edges) + 1L]] <- c(v, idx(i, j + 1L))
    }
  }
  do.call(rbind, edges)
}

edges_grid3d <- function(nx, ny = nx, nz = nx) {
  stopifnot(nx >= 1, ny >= 1, nz >= 1)
  idx <- function(i, j, k) (k - 1L) * nx * ny + (i - 1L) * ny + j
  edges <- list()
  for (k in seq_len(nz)) {
    for (i in seq_len(nx)) {
      for (j in seq_len(ny)) {
        v <- idx(i, j, k)
        if (i < nx) edges[[length(edges) + 1L]] <- c(v, idx(i + 1L, j, k))
        if (j < ny) edges[[length(edges) + 1L]] <- c(v, idx(i, j + 1L, k))
        if (k < nz) edges[[length(edges) + 1L]] <- c(v, idx(i, j, k + 1L))
      }
    }
  }
  do.call(rbind, edges)
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
  coords <- cbind(
    x = cells$x + 0.5,
    y = (side - 1L - cells$y) + 0.5
  )

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

  list(
    edges = do.call(rbind, edges),
    coords = coords
  )
}

edges_sierpinski_triangle <- function(level) {
  stopifnot(level >= 0)

  merge_nodes <- function(edges, from, to) {
    edges[edges == from] <- to
    edges
  }

  build <- function(k) {
    if (k == 0) {
      coords <- rbind(
        c(0, 0),
        c(1, 0),
        c(0.5, sqrt(3) / 2)
      )
      edges <- rbind(c(1L, 2L), c(2L, 3L), c(3L, 1L))
      return(list(edges = edges, coords = coords, corners = c(1L, 2L, 3L), n = 3L))
    }

    left <- build(k - 1)
    right <- build(k - 1)
    top <- build(k - 1)

    left_coords <- left$coords / 2
    right_coords <- right$coords / 2 + matrix(c(0.5, 0), nrow(right$coords), 2, byrow = TRUE)
    top_coords <- top$coords / 2 + matrix(c(0.25, sqrt(3) / 4), nrow(top$coords), 2, byrow = TRUE)

    off1 <- left$n
    off2 <- left$n + right$n

    edges <- rbind(left$edges,
                   right$edges + off1,
                   top$edges + off2)
    coords <- rbind(left_coords, right_coords, top_coords)

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

    new_coords <- matrix(NA_real_, nrow = length(ids), ncol = 2)
    for (i in seq_along(ids)) {
      new_coords[i, ] <- coords[ids[i], ]
    }

    corners <- c(map[as.character(L[1])],
                 map[as.character(R[2])],
                 map[as.character(T[3])])

    list(edges = edges, coords = new_coords, corners = corners, n = length(ids))
  }

  build(level)
}

normalize_coords <- function(coords) {
  centered <- scale(coords, center = TRUE, scale = FALSE)
  radius <- max(sqrt(rowSums(centered^2)))
  if (!is.finite(radius) || radius <= 0) {
    return(centered)
  }
  centered / radius
}

align_2d <- function(target, source) {
  target_norm <- normalize_coords(target)
  source_norm <- normalize_coords(source)
  fit <- svd(t(source_norm) %*% target_norm)
  rotation <- fit$u %*% t(fit$v)
  source_norm %*% rotation
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

project_perspective <- function(coords, camera = 4.2) {
  z_camera <- coords[, 3] + camera
  scale <- camera / pmax(z_camera, 1e-6)
  cbind(x = coords[, 1] * scale, y = coords[, 2] * scale, z = coords[, 3])
}

plot_2d_layout <- function(path) {
  edges <- edges_cycle(18)
  coords <- grip.layout(edges, n = 18, dim = 2,
                        engine = "mish_v6",
                        placement = "circle",
                        rounds = 10, final_rounds = 10,
                        num_init = 6, num_nbrs = 6,
                        seed = 2)

  x <- coords[, 1]
  y <- coords[, 2]
  xr <- range(x)
  yr <- range(y)
  xpad <- diff(xr) * 0.15
  ypad <- diff(yr) * 0.15

  png(path, width = 1200, height = 1200, res = 160, bg = "#f7f3ea")
  on.exit(dev.off(), add = TRUE)
  par(mar = c(0, 0, 2.6, 0), xaxs = "i", yaxs = "i")
  plot(NA,
       xlim = xr + c(-xpad, xpad),
       ylim = yr + c(-ypad, ypad),
       asp = 1,
       axes = FALSE,
       xlab = "",
       ylab = "")
  apply(edges, 1, function(e) {
    segments(coords[e[1], 1], coords[e[1], 2],
             coords[e[2], 1], coords[e[2], 2],
             col = grDevices::adjustcolor("#0f3b5f", alpha.f = 0.34),
             lwd = 2.7)
  })
  points(coords[, 1], coords[, 2],
         pch = 21,
         bg = "#f05a28",
         col = "#16324f",
         cex = 4.2,
         lwd = 1.5)
  title("2D Cycle Layout", col.main = "#16324f", cex.main = 1.4)
}

draw_2d_frame <- function(path, coords, edges, title_text,
                          width = 900, height = 900) {
  xr <- range(coords[, 1])
  yr <- range(coords[, 2])
  xpad <- diff(xr) * 0.16
  ypad <- diff(yr) * 0.16
  if (!is.finite(xpad) || xpad == 0) xpad <- 0.2
  if (!is.finite(ypad) || ypad == 0) ypad <- 0.2

  n <- nrow(coords)
  node_cex <- if (n <= 80) {
    1.5
  } else if (n <= 180) {
    1.0
  } else if (n <= 600) {
    0.56
  } else if (n <= 1500) {
    0.34
  } else {
    0.22
  }
  edge_lwd <- if (n <= 80) {
    2.0
  } else if (n <= 180) {
    1.35
  } else if (n <= 600) {
    0.75
  } else if (n <= 1500) {
    0.46
  } else {
    0.30
  }

  png(path, width = width, height = height, res = 160, bg = "#f7f3ea")
  on.exit(dev.off(), add = TRUE)
  par(mar = c(0, 0, 2.4, 0), xaxs = "i", yaxs = "i")
  plot(NA,
       xlim = xr + c(-xpad, xpad),
       ylim = yr + c(-ypad, ypad),
       asp = 1,
       axes = FALSE,
       xlab = "",
       ylab = "")
  apply(edges, 1, function(e) {
    segments(coords[e[1], 1], coords[e[1], 2],
             coords[e[2], 1], coords[e[2], 2],
             col = grDevices::adjustcolor("#0f3b5f", alpha.f = 0.22),
             lwd = edge_lwd)
  })
  points(coords[, 1], coords[, 2],
         pch = 21,
         bg = "#f05a28",
         col = "#16324f",
         cex = node_cex,
         lwd = 0.9)
  title(title_text, col.main = "#16324f", cex.main = 1.3)
}

plot_2d_morph_gif <- function(path, from_coords, to_coords, edges, title_text,
                              width = 800, height = 800) {
  if (!requireNamespace("magick", quietly = TRUE)) {
    message("Skipping GIF generation: package 'magick' is not installed.")
    return(invisible(FALSE))
  }

  from_norm <- normalize_coords(from_coords)
  to_norm <- align_2d(from_norm, to_coords)

  frame_dir <- tempfile("grip-sierpinski-frames-")
  dir.create(frame_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(frame_dir, recursive = TRUE, force = TRUE), add = TRUE)

  phases <- c(rep(0, 2),
              seq(0, 1, length.out = 7)[-1],
              rep(1, 2),
              seq(1, 0, length.out = 7)[-1],
              rep(0, 2))
  frame_paths <- file.path(frame_dir, sprintf("frame-%02d.png", seq_along(phases)))

  for (i in seq_along(phases)) {
    t <- 0.5 - 0.5 * cos(pi * phases[i])
    coords <- (1 - t) * from_norm + t * to_norm
    draw_2d_frame(frame_paths[i], coords, edges, title_text, width = width, height = height)
  }

  frames <- magick::image_read(frame_paths)
  anim <- magick::image_animate(frames, fps = 10, optimize = TRUE)
  magick::image_write(anim, path)
  invisible(TRUE)
}

draw_3d_projection <- function(path, coords, edges, yaw, pitch,
                               title_text,
                               width = 1200, height = 1200) {
  rotated <- rotate_xyz(normalize_coords(coords), yaw = yaw, pitch = pitch)
  projected <- project_perspective(rotated)

  xr <- range(projected[, 1])
  yr <- range(projected[, 2])
  xpad <- diff(xr) * 0.18
  ypad <- diff(yr) * 0.18

  edge_depth <- rowMeans(matrix(rotated[as.vector(t(edges)), 3], ncol = 2))
  edge_order <- order(edge_depth)
  point_order <- order(rotated[, 3])
  z_range <- range(rotated[, 3])
  z_span <- diff(z_range)
  if (!is.finite(z_span) || z_span == 0) {
    z_span <- 1
  }
  point_level <- (rotated[, 3] - z_range[1]) / z_span

  png(path, width = width, height = height, res = 160, bg = "#f7f3ea")
  on.exit(dev.off(), add = TRUE)
  par(mar = c(0, 0, 2.6, 0), xaxs = "i", yaxs = "i")
  plot(NA,
       xlim = xr + c(-xpad, xpad),
       ylim = yr + c(-ypad, ypad),
       asp = 1,
       axes = FALSE,
       xlab = "",
       ylab = "")

  for (idx in edge_order) {
    e <- edges[idx, ]
    alpha <- 0.14 + 0.42 * ((edge_depth[idx] - min(edge_depth)) /
      max(diff(range(edge_depth)), 1e-6))
    segments(projected[e[1], 1], projected[e[1], 2],
             projected[e[2], 1], projected[e[2], 2],
             col = grDevices::adjustcolor("#1f4d6b", alpha.f = alpha),
             lwd = 2.5)
  }

  palette <- grDevices::colorRampPalette(c("#f7b267", "#f4845f", "#d1495b", "#7b2d26"))(100)
  cols <- palette[pmax(1L, pmin(100L, floor(point_level * 99) + 1L))]
  cex_vals <- 1.0 + 1.4 * point_level
  points(projected[point_order, 1], projected[point_order, 2],
         pch = 21,
         bg = cols[point_order],
         col = "#16324f",
         cex = cex_vals[point_order],
         lwd = 1.2)
  title(title_text, col.main = "#16324f", cex.main = 1.4)
}

plot_3d_layout <- function(path) {
  edges <- edges_grid3d(3, 3, 3)
  coords <- grip.layout(edges, n = 27, dim = 3,
                        engine = "mish_v6",
                        placement = "barycenter",
                        rounds = 12, final_rounds = 12,
                        num_init = 6, num_nbrs = 7,
                        seed = 3)
  draw_3d_projection(path, coords, edges,
                     yaw = 0.8, pitch = 0.7,
                     title_text = "3D Lattice Layout")
  invisible(list(coords = coords, edges = edges))
}

plot_3d_gif <- function(path, coords, edges) {
  if (!requireNamespace("magick", quietly = TRUE)) {
    message("Skipping GIF generation: package 'magick' is not installed.")
    return(invisible(FALSE))
  }

  frame_dir <- tempfile("grip-readme-frames-")
  dir.create(frame_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(frame_dir, recursive = TRUE, force = TRUE), add = TRUE)

  angles <- seq(0, 2 * pi, length.out = 25L)[-25L]
  frame_paths <- file.path(frame_dir, sprintf("frame-%02d.png", seq_along(angles)))
  for (i in seq_along(angles)) {
    draw_3d_projection(frame_paths[i], coords, edges,
                       yaw = angles[i], pitch = 0.62,
                       title_text = "Animated 3D Rotation",
                       width = 900, height = 900)
  }

  frames <- magick::image_read(frame_paths)
  anim <- magick::image_animate(frames, fps = 10, optimize = TRUE)
  magick::image_write(anim, path)
  invisible(TRUE)
}

sierpinski2 <- edges_sierpinski_carpet(2)
sierpinski2_layout <- grip.layout(sierpinski2$edges,
                                  n = nrow(sierpinski2$coords),
                                  dim = 2,
                                  engine = "mish_v6",
                                  placement = "barycenter",
                                  rounds = 12, final_rounds = 12,
                                  num_init = 10, num_nbrs = 8,
                                  seed = 22)

sierpinski3 <- edges_sierpinski_carpet(3)
sierpinski3_layout <- grip.layout(sierpinski3$edges,
                                  n = nrow(sierpinski3$coords),
                                  dim = 2,
                                  engine = "mish_v6",
                                  placement = "barycenter",
                                  rounds = 14, final_rounds = 14,
                                  num_init = 18, num_nbrs = 10,
                                  seed = 23)

sierpinski4 <- edges_sierpinski_carpet(4)
sierpinski4_layout <- grip.layout(sierpinski4$edges,
                                  n = nrow(sierpinski4$coords),
                                  dim = 2,
                                  engine = "mish_v6",
                                  placement = "barycenter",
                                  rounds = 16, final_rounds = 16,
                                  num_init = 24, num_nbrs = 12,
                                  seed = 24)

sierpinski_triangle2 <- edges_sierpinski_triangle(2)
sierpinski_triangle2_layout <- grip.layout(sierpinski_triangle2$edges,
                                           n = nrow(sierpinski_triangle2$coords),
                                           dim = 2,
                                           engine = "mish_v6",
                                           placement = "circle",
                                           rounds = 25, final_rounds = 25,
                                           num_init = 5, num_nbrs = 7,
                                           seed = 4)

sierpinski_triangle3 <- edges_sierpinski_triangle(3)
sierpinski_triangle3_layout <- grip.layout(sierpinski_triangle3$edges,
                                           n = nrow(sierpinski_triangle3$coords),
                                           dim = 2,
                                           engine = "mish_v6",
                                           placement = "circle",
                                           rounds = 28, final_rounds = 28,
                                           num_init = 6, num_nbrs = 8,
                                           seed = 24)

sierpinski_triangle4 <- edges_sierpinski_triangle(4)
sierpinski_triangle4_layout <- grip.layout(sierpinski_triangle4$edges,
                                           n = nrow(sierpinski_triangle4$coords),
                                           dim = 2,
                                           engine = "mish_v6",
                                           placement = "circle",
                                           rounds = 32, final_rounds = 32,
                                           num_init = 7, num_nbrs = 9,
                                           seed = 25)

plot_2d_layout("man/figures/readme-layout-2d-cycle.png")
mesh <- plot_3d_layout("man/figures/readme-layout-3d-mesh.png")
plot_3d_gif("man/figures/readme-layout-3d-rotation.gif", mesh$coords, mesh$edges)
plot_2d_morph_gif("man/figures/readme-sierpinski-carpet-level-2.gif",
                  sierpinski2$coords, sierpinski2_layout, sierpinski2$edges,
                  "Sierpinski Carpet (Level 2)")
plot_2d_morph_gif("man/figures/readme-sierpinski-carpet-level-3.gif",
                  sierpinski3$coords, sierpinski3_layout, sierpinski3$edges,
                  "Sierpinski Carpet (Level 3)")
plot_2d_morph_gif("man/figures/readme-sierpinski-carpet-level-4.gif",
                  sierpinski4$coords, sierpinski4_layout, sierpinski4$edges,
                  "Sierpinski Carpet (Level 4)",
                  width = 900, height = 900)
plot_2d_morph_gif("man/figures/readme-sierpinski-triangle-level-2.gif",
                  sierpinski_triangle2$coords, sierpinski_triangle2_layout, sierpinski_triangle2$edges,
                  "Sierpinski Triangle (Level 2)")
plot_2d_morph_gif("man/figures/readme-sierpinski-triangle-level-3.gif",
                  sierpinski_triangle3$coords, sierpinski_triangle3_layout, sierpinski_triangle3$edges,
                  "Sierpinski Triangle (Level 3)")
plot_2d_morph_gif("man/figures/readme-sierpinski-triangle-level-4.gif",
                  sierpinski_triangle4$coords, sierpinski_triangle4_layout, sierpinski_triangle4$edges,
                  "Sierpinski Triangle (Level 4)",
                  width = 900, height = 900)

message("README assets written to man/figures/")
