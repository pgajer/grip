#!/usr/bin/env Rscript

dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'devtools' or the 'grip' package to generate README assets.")
}

if (!requireNamespace("magick", quietly = TRUE)) {
  stop("Package 'magick' is required to generate README GIF assets.")
}

build_sierpinski_carpet <- function(level) {
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

  edges <- do.call(rbind, edges)
  package_edges <- edges.sierpinski.carpet(level)
  if (!identical(package_edges, edges)) {
    stop(sprintf(
      "Canonical carpet builder does not match package edges at level %d",
      level
    ))
  }
  list(edges = edges, coords = coords)
}

build_sierpinski_triangle <- function(level) {
  stopifnot(level >= 0)

  merge_nodes <- function(edges, from, to) {
    edges[edges == from] <- to
    edges
  }

  build <- function(k) {
    if (k == 0L) {
      coords <- rbind(
        c(0, 0),
        c(1, 0),
        c(0.5, sqrt(3) / 2)
      )
      edges <- rbind(c(1L, 2L), c(2L, 3L), c(3L, 1L))
      return(list(edges = edges, coords = coords, corners = c(1L, 2L, 3L), n = 3L))
    }

    left <- build(k - 1L)
    right <- build(k - 1L)
    top <- build(k - 1L)

    left_coords <- left$coords / 2
    right_coords <- right$coords / 2 +
      matrix(c(0.5, 0), nrow(right$coords), 2, byrow = TRUE)
    top_coords <- top$coords / 2 +
      matrix(c(0.25, sqrt(3) / 4), nrow(top$coords), 2, byrow = TRUE)

    off1 <- left$n
    off2 <- left$n + right$n

    edges <- rbind(
      left$edges,
      right$edges + off1,
      top$edges + off2
    )
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
    edges <- cbind(map[as.character(edges[, 1])], map[as.character(edges[, 2])])
    coords <- coords[ids, , drop = FALSE]
    edges <- t(apply(edges, 1, sort))
    edges <- unique(edges)

    corners <- c(
      map[as.character(L[1])],
      map[as.character(R[2])],
      map[as.character(T[3])]
    )

    list(edges = edges, coords = coords, corners = corners, n = length(ids))
  }

  out <- build(level)
  package_edges <- edges.sierpinski.triangle(level)
  if (!identical(package_edges, out$edges)) {
    stop(sprintf(
      "Canonical triangle builder does not match package edges at level %d",
      level
    ))
  }
  list(edges = out$edges, coords = out$coords, corners = out$corners)
}

center_trace_frame <- function(coords) {
  active <- stats::complete.cases(coords)
  if (!any(active)) {
    return(coords)
  }
  centered <- coords
  centered[active, ] <- sweep(
    coords[active, , drop = FALSE],
    2L,
    colMeans(coords[active, , drop = FALSE]),
    "-"
  )
  centered
}

normalize_trace_frames <- function(frames) {
  centered <- lapply(frames, center_trace_frame)
  radius <- max(vapply(centered, function(z) {
    active <- stats::complete.cases(z)
    if (!any(active)) {
      return(0)
    }
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

align_trace_frame_to_target <- function(coords, target_coords) {
  active <- stats::complete.cases(coords) & stats::complete.cases(target_coords)
  if (!any(active)) {
    return(coords)
  }

  aligned <- coords
  src <- coords[active, , drop = FALSE]
  dst <- target_coords[active, , drop = FALSE]

  src <- sweep(src, 2L, colMeans(src), "-")
  dst <- sweep(dst, 2L, colMeans(dst), "-")

  src_radius <- max(sqrt(rowSums(src^2)))
  dst_radius <- max(sqrt(rowSums(dst^2)))
  if (is.finite(src_radius) && src_radius > 0) {
    src <- src / src_radius
  }
  if (is.finite(dst_radius) && dst_radius > 0) {
    dst <- dst / dst_radius
  }

  if (sum(active) >= 2L) {
    cross <- t(src) %*% dst
    sv <- svd(cross)
    rot <- sv$u %*% t(sv$v)
    if (det(rot) < 0) {
      fix <- diag(ncol(rot))
      fix[ncol(fix), ncol(fix)] <- -1
      rot <- sv$u %*% fix %*% t(sv$v)
    }
    src <- src %*% rot
  }

  aligned[active, ] <- src
  aligned
}

align_trace_frames_to_target <- function(frames, target_coords) {
  lapply(frames, align_trace_frame_to_target, target_coords = target_coords)
}

canonicalize_triangle_symmetry <- function(frames, target_coords, target_corners) {
  if (is.null(target_coords) || is.null(target_corners) ||
      length(target_corners) < 3L || length(frames) == 0L) {
    return(frames)
  }

  final_frame <- frames[[length(frames)]]
  active <- stats::complete.cases(final_frame) & stats::complete.cases(target_coords)
  if (!any(active)) {
    return(frames)
  }
  if (!all(stats::complete.cases(final_frame[target_corners, , drop = FALSE]))) {
    return(frames)
  }

  normalize_coords_local <- function(coords) {
    centered <- sweep(coords, 2L, colMeans(coords), "-")
    radius <- max(sqrt(rowSums(centered^2)))
    if (is.finite(radius) && radius > 0) {
      centered <- centered / radius
    }
    centered
  }

  fit_orthogonal_transform <- function(src, dst) {
    cross <- t(src) %*% dst
    sv <- svd(cross)
    sv$u %*% t(sv$v)
  }

  current_coords <- normalize_coords_local(final_frame[active, , drop = FALSE])
  target_coords_norm <- normalize_coords_local(target_coords[active, , drop = FALSE])
  current_corner_coords <- normalize_coords_local(final_frame[target_corners, , drop = FALSE])
  target_corner_coords <- normalize_coords_local(target_coords[target_corners, , drop = FALSE])

  perms <- rbind(
    c(1L, 2L, 3L),
    c(1L, 3L, 2L),
    c(2L, 1L, 3L),
    c(2L, 3L, 1L),
    c(3L, 1L, 2L),
    c(3L, 2L, 1L)
  )

  transforms <- lapply(seq_len(nrow(perms)), function(i) {
    fit_orthogonal_transform(
      current_corner_coords[perms[i, ], , drop = FALSE],
      target_corner_coords
    )
  })

  best_idx <- 1L
  best_err <- Inf
  for (i in seq_along(transforms)) {
    transformed <- current_coords %*% transforms[[i]]
    err <- sum((transformed - target_coords_norm)^2)
    if (err < best_err) {
      best_err <- err
      best_idx <- i
    }
  }

  transform <- transforms[[best_idx]]
  lapply(frames, function(z) {
    active <- stats::complete.cases(z)
    if (any(active)) {
      z[active, ] <- z[active, , drop = FALSE] %*% transform
    }
    z
  })
}

subsample_trace <- function(trace_obj, motion_frames) {
  total <- length(trace_obj$frames)
  if (is.null(motion_frames) || motion_frames >= total) {
    return(trace_obj)
  }

  idx <- floor(seq(0, total - 1L, length.out = motion_frames)) + 1L
  idx <- pmax.int(1L, pmin.int(total, idx))
  idx <- unique(idx)
  if (length(idx) < motion_frames) {
    extras <- setdiff(seq_len(total), idx)
    idx <- sort(c(idx, extras[seq_len(motion_frames - length(idx))]))
  }

  list(
    final = trace_obj$final,
    frames = trace_obj$frames[idx],
    meta = trace_obj$meta[idx, , drop = FALSE]
  )
}

compute_2d_limits <- function(frames, pad_frac = 0.06, min_pad = 0.06) {
  xs <- c()
  ys <- c()
  for (coords in frames) {
    active <- stats::complete.cases(coords)
    if (!any(active)) {
      next
    }
    xs <- c(xs, coords[active, 1L])
    ys <- c(ys, coords[active, 2L])
  }
  xr <- range(xs)
  yr <- range(ys)
  xpad <- max(diff(xr) * pad_frac, min_pad)
  ypad <- max(diff(yr) * pad_frac, min_pad)
  list(xlim = xr + c(-xpad, xpad), ylim = yr + c(-ypad, ypad))
}

edge_mask <- function(coords, edges) {
  active <- stats::complete.cases(coords)
  active[edges[, 1L]] & active[edges[, 2L]]
}

draw_readme_trace_frame <- function(path, coords, edges, xlim, ylim,
                                    width = 1200, height = 1200) {
  n <- nrow(coords)
  node_cex <- if (n <= 80) {
    1.65
  } else if (n <= 180) {
    1.15
  } else if (n <= 600) {
    0.62
  } else if (n <= 1500) {
    0.32
  } else {
    0.14
  }
  edge_lwd <- if (n <= 80) {
    2.0
  } else if (n <= 180) {
    1.35
  } else if (n <= 600) {
    0.78
  } else if (n <= 1500) {
    0.40
  } else {
    0.17
  }

  keep_edges <- edge_mask(coords, edges)
  active <- stats::complete.cases(coords)

  png(path, width = width, height = height, res = 180, bg = "#f7f3ea")
  on.exit(dev.off(), add = TRUE)
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot(
    NA,
    xlim = xlim,
    ylim = ylim,
    asp = 1,
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  if (any(keep_edges)) {
    apply(edges[keep_edges, , drop = FALSE], 1, function(e) {
      graphics::segments(
        coords[e[1L], 1L], coords[e[1L], 2L],
        coords[e[2L], 1L], coords[e[2L], 2L],
        col = grDevices::adjustcolor("#0f3b5f", alpha.f = 0.18),
        lwd = edge_lwd
      )
    })
  }
  if (any(active)) {
    graphics::points(
      coords[active, 1L], coords[active, 2L],
      pch = 21,
      bg = "#f05a28",
      col = "#16324f",
      cex = node_cex,
      lwd = 0.75
    )
  }
}

plot_readme_trace_gif <- function(path, trace_obj, edges,
                                  width = 1200, height = 1200,
                                  fps = 5, motion_frames = 112,
                                  final_hold_sec = 8,
                                  target_coords = NULL,
                                  target_corners = NULL,
                                  pad_frac = 0.06,
                                  min_pad = 0.06) {
  trace_obj <- subsample_trace(trace_obj, motion_frames = motion_frames)
  frames <- trace_obj$frames
  if (!is.null(target_coords)) {
    frames <- align_trace_frames_to_target(frames, target_coords)
    frames <- canonicalize_triangle_symmetry(frames, target_coords, target_corners)
  }
  frames <- normalize_trace_frames(frames)
  limits <- compute_2d_limits(frames, pad_frac = pad_frac, min_pad = min_pad)
  frame_dir <- tempfile("grip-readme-trace-")
  dir.create(frame_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(frame_dir, recursive = TRUE, force = TRUE), add = TRUE)

  frame_paths <- file.path(frame_dir, sprintf("frame-%03d.png", seq_along(frames)))
  for (i in seq_along(frames)) {
    draw_readme_trace_frame(
      frame_paths[i],
      frames[[i]],
      edges,
      xlim = limits$xlim,
      ylim = limits$ylim,
      width = width,
      height = height
    )
  }

  hold_frames <- max(1L, as.integer(round(fps * final_hold_sec)))
  frame_paths <- c(frame_paths, rep(tail(frame_paths, 1L), hold_frames))
  anim <- magick::image_animate(
    magick::image_read(frame_paths),
    fps = fps,
    optimize = TRUE
  )
  magick::image_write(anim, path)
  invisible(path)
}

carpet4 <- build_sierpinski_carpet(4)
triangle6 <- build_sierpinski_triangle(6)

carpet4_trace_2d <- grip.layout.trace(
  carpet4$edges,
  n = max(carpet4$edges),
  dim = 2,
  preset = "carpet",
  trace = "round",
  trace.every = 1,
  seed = 24
)

# Match the diagnostic triangle panels by using the carpet preset and aligning
# each trace frame to the canonical recursive triangle before rendering.
triangle6_trace_2d <- grip.layout.trace(
  triangle6$edges,
  n = max(triangle6$edges),
  dim = 2,
  preset = "carpet",
  trace = "round",
  trace.every = 1,
  seed = 1
)

plot_readme_trace_gif(
  "man/figures/readme-sierpinski-carpet-level-4-trace.gif",
  carpet4_trace_2d,
  carpet4$edges,
  target_coords = carpet4$coords,
  pad_frac = 0.025,
  min_pad = 0.025
)

plot_readme_trace_gif(
  "man/figures/readme-sierpinski-triangle-level-6-trace.gif",
  triangle6_trace_2d,
  triangle6$edges,
  target_coords = triangle6$coords,
  target_corners = triangle6$corners
)

message("README assets written to man/figures/")
