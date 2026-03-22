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

if (!requireNamespace("magick", quietly = TRUE)) {
  stop("Package 'magick' is required to generate GIF trace assets.")
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

subsample_trace <- function(trace_obj, motion_frames) {
  total <- length(trace_obj$frames)
  if (is.null(motion_frames) || motion_frames >= total) {
    return(trace_obj)
  }

  idx <- floor(seq(0, total - 1, length.out = motion_frames)) + 1L
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
                              width = 900, height = 900, fps = 5,
                              motion_frames = 112,
                              final_hold_sec = 25) {
  trace_obj <- subsample_trace(trace_obj, motion_frames = motion_frames)
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

  hold_frames <- max(1L, as.integer(round(fps * final_hold_sec)))
  frame_paths <- c(frame_paths, rep(tail(frame_paths, 1L), hold_frames))

  anim <- magick::image_animate(magick::image_read(frame_paths),
                                fps = fps,
                                optimize = TRUE)
  magick::image_write(anim, path)

  list(
    motion_frames = length(frames),
    hold_frames = hold_frames,
    total_frames = length(frame_paths)
  )
}

message("Computing requested mish_v6 trace animations...")

carpet4 <- edges_sierpinski_carpet(4)
triangle5 <- edges_sierpinski_triangle(5)

carpet4_trace_2d <- grip.layout.trace(carpet4,
                                      n = max(carpet4),
                                      dim = 2,
                                      engine = "mish_v6",
                                      placement = "barycenter",
                                      rounds = 128,
                                      final_rounds = 128,
                                      num_init = 24,
                                      num_nbrs = 12,
                                      trace = "round",
                                      trace.every = 1,
                                      seed = 24)

triangle5_trace_2d <- grip.layout.trace(triangle5,
                                        n = max(triangle5),
                                        dim = 2,
                                        engine = "mish_v6",
                                        placement = "circle",
                                        rounds = 128,
                                        final_rounds = 128,
                                        num_init = 7,
                                        num_nbrs = 9,
                                        trace = "round",
                                        trace.every = 1,
                                        seed = 25)

carpet_info <- plot_2d_trace_gif(
  file.path(fig_dir, "trace-sierpinski-carpet-level-4-2d-v6-112f.gif"),
  carpet4_trace_2d,
  carpet4,
  title_text = "Sierpinski Carpet (Level 4, 2D Trace, mish_v6)"
)

triangle_info <- plot_2d_trace_gif(
  file.path(fig_dir, "trace-sierpinski-triangle-level-5-2d-v6-112f.gif"),
  triangle5_trace_2d,
  triangle5,
  title_text = "Sierpinski Triangle (Level 5, 2D Trace, mish_v6)"
)

message(sprintf(
  "carpet: motion=%d, hold=%d, total=%d",
  carpet_info$motion_frames,
  carpet_info$hold_frames,
  carpet_info$total_frames
))
message(sprintf(
  "triangle: motion=%d, hold=%d, total=%d",
  triangle_info$motion_frames,
  triangle_info$hold_frames,
  triangle_info$total_frames
))
message(sprintf("Requested trace assets written to %s", fig_dir))
