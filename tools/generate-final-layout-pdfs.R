#!/usr/bin/env Rscript

pdf_dir <- file.path("dev", "manual", "pdf")
preview_dir <- file.path("dev", "manual", "tmp", "pdfs")
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'devtools' or the 'grip' package to generate final layout PDFs.")
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

normalize_coords <- function(coords) {
  centered <- scale(coords, center = TRUE, scale = FALSE)
  radius <- max(sqrt(rowSums(centered^2)))
  if (!is.finite(radius) || radius <= 0) {
    return(centered)
  }
  centered / radius
}

plot_layout_panel <- function(coords, edges, title_text, subtitle_text,
                              xlim = NULL, ylim = NULL) {
  coords <- normalize_coords(coords)
  if (is.null(xlim) || is.null(ylim)) {
    xr <- range(coords[, 1])
    yr <- range(coords[, 2])
    xpad <- max(diff(xr) * 0.12, 0.12)
    ypad <- max(diff(yr) * 0.12, 0.12)
    xlim <- xr + c(-xpad, xpad)
    ylim <- yr + c(-ypad, ypad)
  }

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

  plot(NA,
       xlim = xlim,
       ylim = ylim,
       asp = 1,
       axes = FALSE,
       xlab = "",
       ylab = "")
  apply(edges, 1, function(e) {
    segments(coords[e[1], 1], coords[e[1], 2],
             coords[e[2], 1], coords[e[2], 2],
             col = grDevices::adjustcolor("#0f3b5f", alpha.f = 0.18),
             lwd = edge_lwd)
  })
  points(coords[, 1], coords[, 2],
         pch = 21,
         bg = "#f05a28",
         col = "#16324f",
         cex = node_cex,
         lwd = 0.7)
  title(main = title_text, sub = subtitle_text,
        col.main = "#16324f", cex.main = 1.2,
        col.sub = "#466074", cex.sub = 0.9)
}

write_single_pdf <- function(path, coords, edges, title_text, subtitle_text) {
  grDevices::pdf(path, width = 8.5, height = 8.5,
                 paper = "special", bg = "#f7f3ea", useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  par(mar = c(0, 0, 2.8, 0), xaxs = "i", yaxs = "i")
  plot_layout_panel(coords, edges, title_text, subtitle_text)
}

write_compare_pdf <- function(path, left_coords, right_coords, edges,
                              left_title, right_title, subtitle_text) {
  left_norm <- normalize_coords(left_coords)
  right_norm <- normalize_coords(right_coords)
  xs <- c(left_norm[, 1], right_norm[, 1])
  ys <- c(left_norm[, 2], right_norm[, 2])
  xr <- range(xs)
  yr <- range(ys)
  xpad <- max(diff(xr) * 0.12, 0.12)
  ypad <- max(diff(yr) * 0.12, 0.12)
  xlim <- xr + c(-xpad, xpad)
  ylim <- yr + c(-ypad, ypad)

  grDevices::pdf(path, width = 15.5, height = 8.5,
                 paper = "special", bg = "#f7f3ea", useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  par(mfrow = c(1, 2), mar = c(0, 0, 2.8, 0), xaxs = "i", yaxs = "i")
  plot_layout_panel(left_coords, edges, left_title, subtitle_text, xlim = xlim, ylim = ylim)
  plot_layout_panel(right_coords, edges, right_title, subtitle_text, xlim = xlim, ylim = ylim)
}

compute_layout <- function(edges, placement, rounds, final_rounds,
                           num_init, num_nbrs, seed) {
  grip.layout(edges = edges,
              n = max(edges),
              dim = 2,
              placement = placement,
              rounds = rounds,
              final_rounds = final_rounds,
              num_init = num_init,
              num_nbrs = num_nbrs,
              seed = seed)
}

message("Computing final 2D layouts for PDF export...")

carpet4 <- edges_sierpinski_carpet(4)
triangle4 <- edges_sierpinski_triangle(4)

carpet128 <- compute_layout(carpet4, "barycenter", 128, 128, 24, 12, 24)
triangle128 <- compute_layout(triangle4, "circle", 128, 128, 7, 9, 25)

carpet256 <- compute_layout(carpet4, "barycenter", 256, 256, 24, 12, 24)
triangle256 <- compute_layout(triangle4, "circle", 256, 256, 7, 9, 25)

write_single_pdf(file.path(pdf_dir, "sierpinski-carpet-level-4-final-128.pdf"),
                 carpet128, carpet4,
                 "Sierpinski Carpet Level 4",
                 "rounds=128, final_rounds=128")
write_compare_pdf(file.path(pdf_dir, "sierpinski-carpet-level-4-final-128-compare.pdf"),
                  carpet128, carpet256, carpet4,
                  "128 rounds", "256 rounds",
                  "current layout implementation")

write_single_pdf(file.path(pdf_dir, "sierpinski-triangle-level-4-final-128.pdf"),
                 triangle128, triangle4,
                 "Sierpinski Triangle Level 4",
                 "rounds=128, final_rounds=128")
write_compare_pdf(file.path(pdf_dir, "sierpinski-triangle-level-4-final-128-compare.pdf"),
                  triangle128, triangle256, triangle4,
                  "128 rounds", "256 rounds",
                  "current layout implementation")

write_single_pdf(file.path(pdf_dir, "sierpinski-carpet-level-4-final-256.pdf"),
                 carpet256, carpet4,
                 "Sierpinski Carpet Level 4",
                 "rounds=256, final_rounds=256")

write_single_pdf(file.path(pdf_dir, "sierpinski-triangle-level-4-final-256.pdf"),
                 triangle256, triangle4,
                 "Sierpinski Triangle Level 4",
                 "rounds=256, final_rounds=256")

message(sprintf("PDF assets written to %s", pdf_dir))
