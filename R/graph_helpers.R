# Helper utilities shared by the exported graph generators.
.empty_edge_matrix <- function() {
  matrix(integer(), ncol = 2L)
}

.as_whole_number <- function(x, name, min = 0L) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x != as.integer(x) || x < min) {
    stop(sprintf("%s must be a single integer >= %d", name, min), call. = FALSE)
  }
  as.integer(x)
}

.bind_edges <- function(edges) {
  if (length(edges) == 0L) {
    return(.empty_edge_matrix())
  }
  out <- do.call(rbind, edges)
  storage.mode(out) <- "integer"
  out
}

.normalize_undirected_edges <- function(edges) {
  if (is.null(edges) || length(edges) == 0L || nrow(edges) == 0L) {
    return(.empty_edge_matrix())
  }
  edges <- cbind(pmin(edges[, 1L], edges[, 2L]),
                 pmax(edges[, 1L], edges[, 2L]))
  edges <- edges[edges[, 1L] != edges[, 2L], , drop = FALSE]
  if (nrow(edges) == 0L) {
    return(.empty_edge_matrix())
  }
  storage.mode(edges) <- "integer"
  unique(edges)
}

.as_finite_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(sprintf("%s must be a single finite numeric value", name), call. = FALSE)
  }
  as.double(x)
}

.as_positive_scalar <- function(x, name) {
  x <- .as_finite_scalar(x, name)
  if (x <= 0) {
    stop(sprintf("%s must be > 0", name), call. = FALSE)
  }
  x
}

.grid_axis <- function(n) {
  n <- .as_whole_number(n, "n", min = 1L)
  if (n == 1L) {
    return(0)
  }
  seq(-1, 1, length.out = n)
}

.mesh.param.coords <- function(h, w, x_scale = 1, y_scale = 1) {
  h <- .as_whole_number(h, "h", min = 1L)
  w <- .as_whole_number(w, "w", min = 1L)
  x_scale <- .as_positive_scalar(x_scale, "x_scale")
  y_scale <- .as_positive_scalar(y_scale, "y_scale")

  x_vals <- .grid_axis(w) * x_scale
  y_vals <- rev(.grid_axis(h) * y_scale)
  coords <- matrix(0, nrow = h * w, ncol = 2L)
  for (i in seq_len(h)) {
    for (j in seq_len(w)) {
      id <- (i - 1L) * w + j
      coords[id, ] <- c(x_vals[[j]], y_vals[[i]])
    }
  }
  colnames(coords) <- c("u", "v")
  coords
}

.cylinder.param.coords <- function(h, w, radius = 1, height = 2) {
  h <- .as_whole_number(h, "h", min = 1L)
  w <- .as_whole_number(w, "w", min = 3L)
  radius <- .as_positive_scalar(radius, "radius")
  height <- .as_positive_scalar(height, "height")

  theta_vals <- seq(0, 2 * pi * (1 - 1 / w), length.out = w)
  z_vals <- rev(.grid_axis(h) * (height / 2))
  coords <- matrix(0, nrow = h * w, ncol = 2L)
  for (i in seq_len(h)) {
    for (j in seq_len(w)) {
      id <- (i - 1L) * w + j
      coords[id, ] <- c(radius * theta_vals[[j]], z_vals[[i]])
    }
  }
  colnames(coords) <- c("theta_arc", "z")
  coords
}

.torus.param.coords <- function(h, w, major_radius = 2, minor_radius = 0.75) {
  h <- .as_whole_number(h, "h", min = 3L)
  w <- .as_whole_number(w, "w", min = 3L)
  major_radius <- .as_positive_scalar(major_radius, "major_radius")
  minor_radius <- .as_positive_scalar(minor_radius, "minor_radius")

  theta_vals <- seq(0, 2 * pi * (1 - 1 / w), length.out = w)
  phi_vals <- seq(0, 2 * pi * (1 - 1 / h), length.out = h)
  coords <- matrix(0, nrow = h * w, ncol = 2L)
  for (i in seq_len(h)) {
    for (j in seq_len(w)) {
      id <- (i - 1L) * w + j
      coords[id, ] <- c(major_radius * theta_vals[[j]], minor_radius * phi_vals[[i]])
    }
  }
  colnames(coords) <- c("theta_arc", "phi_arc")
  coords
}

.sphere.param.coords <- function(h, w, radius = 1) {
  h <- .as_whole_number(h, "h", min = 3L)
  w <- .as_whole_number(w, "w", min = 3L)
  radius <- .as_positive_scalar(radius, "radius")

  ring.count <- h - 2L
  n <- 2L + ring.count * w
  theta_vals <- seq(0, 2 * pi * (1 - 1 / w), length.out = w)
  lat_vals <- seq(pi / 2, -pi / 2, length.out = h)
  coords <- matrix(0, nrow = n, ncol = 2L)
  coords[1L, ] <- c(0, radius * lat_vals[[1L]])
  if (ring.count > 0L) {
    for (i in seq_len(ring.count)) {
      lat <- lat_vals[[i + 1L]]
      for (j in seq_len(w)) {
        id <- 1L + (i - 1L) * w + j
        coords[id, ] <- c(radius * theta_vals[[j]], radius * lat)
      }
    }
  }
  coords[n, ] <- c(0, radius * lat_vals[[h]])
  colnames(coords) <- c("theta_arc", "latitude_arc")
  coords
}

.as_square_keep_mask <- function(mask, name = "mask", min_size = 2L) {
  if (!is.matrix(mask) || nrow(mask) != ncol(mask) || nrow(mask) < min_size) {
    stop(sprintf("%s must be a square matrix with at least %d rows and columns",
                 name, min_size),
         call. = FALSE)
  }
  if (!(is.logical(mask) || is.numeric(mask) || is.integer(mask))) {
    stop(sprintf("%s must be a logical or numeric square matrix", name),
         call. = FALSE)
  }
  if (any(is.na(mask)) || (is.numeric(mask) && any(!is.finite(mask)))) {
    stop(sprintf("%s must not contain NA or non-finite values", name),
         call. = FALSE)
  }
  keep <- mask != 0
  storage.mode(keep) <- "logical"
  if (!any(keep)) {
    stop(sprintf("%s must keep at least one cell", name), call. = FALSE)
  }
  keep
}

.as_keep_grid <- function(keep, name = "keep", min_h = 1L, min_w = 1L) {
  if (!is.matrix(keep) || nrow(keep) < min_h || ncol(keep) < min_w) {
    stop(sprintf("%s must be a matrix with at least %d rows and %d columns",
                 name, min_h, min_w),
         call. = FALSE)
  }
  if (!(is.logical(keep) || is.numeric(keep) || is.integer(keep))) {
    stop(sprintf("%s must be a logical or numeric matrix", name), call. = FALSE)
  }
  if (any(is.na(keep)) || (is.numeric(keep) && any(!is.finite(keep)))) {
    stop(sprintf("%s must not contain NA or non-finite values", name), call. = FALSE)
  }
  out <- keep != 0
  storage.mode(out) <- "logical"
  if (!any(out)) {
    stop(sprintf("%s must keep at least one cell", name), call. = FALSE)
  }
  out
}

.as_named_choice <- function(x, choices, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("%s must be a single character value", name), call. = FALSE)
  }
  if (!x %in% choices) {
    stop(sprintf("%s must be one of %s",
                 name,
                 paste(sprintf('"%s"', choices), collapse = ", ")),
         call. = FALSE)
  }
  x
}

.as_odd_whole_number <- function(x, name, min = 1L) {
  x <- .as_whole_number(x, name, min = min)
  if ((x %% 2L) != 1L) {
    stop(sprintf("%s must be odd", name), call. = FALSE)
  }
  x
}

.surface.z.from_uv <- function(u,
                               v,
                               surface,
                               amplitude,
                               freq_u,
                               freq_v) {
  switch(
    surface,
    saddle = amplitude * (u^2 - v^2),
    paraboloid = amplitude * (u^2 + v^2),
    ripple = amplitude * sin(pi * freq_u * u) * cos(pi * freq_v * v)
  )
}

.occupied.mesh.cells <- function(keep) {
  keep <- .as_keep_grid(keep, "keep", min_h = 1L, min_w = 1L)
  h <- nrow(keep)
  w <- ncol(keep)
  grid <- expand.grid(row = seq_len(h), col = seq_len(w))
  keep_idx <- keep[cbind(grid$row, grid$col)]
  cells <- data.frame(
    row = grid$row[keep_idx],
    col = grid$col[keep_idx]
  )
  rownames(cells) <- NULL
  list(
    keep = keep,
    h = h,
    w = w,
    cells = cells
  )
}

.occupied.mesh.param.coords <- function(keep, x_scale = 1, y_scale = 1) {
  spec <- .occupied.mesh.cells(keep)
  x_scale <- .as_positive_scalar(x_scale, "x_scale")
  y_scale <- .as_positive_scalar(y_scale, "y_scale")

  x_vals <- .grid_axis(spec$w) * x_scale
  y_vals <- rev(.grid_axis(spec$h) * y_scale)
  coords <- cbind(
    u = x_vals[spec$cells$col],
    v = y_vals[spec$cells$row]
  )
  storage.mode(coords) <- "double"
  coords
}

.occupied.mesh.edges <- function(keep) {
  spec <- .occupied.mesh.cells(keep)
  cells <- spec$cells
  id_map <- matrix(0L, nrow = spec$h, ncol = spec$w)
  for (i in seq_len(nrow(cells))) {
    id_map[cells$row[i], cells$col[i]] <- i
  }

  edges <- list()
  for (i in seq_len(nrow(cells))) {
    row <- cells$row[i]
    col <- cells$col[i]
    if (col < spec$w) {
      nbr <- id_map[row, col + 1L]
      if (nbr > 0L) edges[[length(edges) + 1L]] <- c(i, nbr)
    }
    if (row < spec$h) {
      nbr <- id_map[row + 1L, col]
      if (nbr > 0L) edges[[length(edges) + 1L]] <- c(i, nbr)
    }
  }

  .bind_edges(edges)
}

.recursive.mask.grid.cells <- function(mask, level) {
  mask <- .as_square_keep_mask(mask, "mask", min_size = 2L)
  level <- .as_whole_number(level, "level", min = 1L)
  base <- as.integer(nrow(mask))
  side <- as.integer(base^level)
  grid <- expand.grid(x = 0:(side - 1L), y = 0:(side - 1L))

  keep_cell <- function(x, y) {
    for (step in seq_len(level)) {
      x_digit <- (x %% base) + 1L
      y_digit <- (y %% base) + 1L
      if (!mask[y_digit, x_digit]) {
        return(FALSE)
      }
      x <- x %/% base
      y <- y %/% base
    }
    TRUE
  }

  keep <- mapply(keep_cell, grid$x, grid$y)
  cells <- grid[keep, , drop = FALSE]
  rownames(cells) <- NULL
  list(
    level = level,
    base = base,
    side = side,
    mask = mask,
    cells = cells
  )
}

.recursive.mask.grid.param.coords <- function(mask,
                                              level,
                                              x_scale = 1,
                                              y_scale = 1) {
  spec <- .recursive.mask.grid.cells(mask, level)
  x_scale <- .as_positive_scalar(x_scale, "x_scale")
  y_scale <- .as_positive_scalar(y_scale, "y_scale")

  x_vals <- .grid_axis(spec$side) * x_scale
  y_vals <- rev(.grid_axis(spec$side) * y_scale)
  coords <- cbind(
    u = x_vals[spec$cells$x + 1L],
    v = y_vals[spec$cells$y + 1L]
  )
  storage.mode(coords) <- "double"
  coords
}

.recursive.mask.grid.edges <- function(mask, level) {
  spec <- .recursive.mask.grid.cells(mask, level)
  side <- spec$side
  cells <- spec$cells

  id_map <- matrix(0L, nrow = side, ncol = side)
  for (i in seq_len(nrow(cells))) {
    id_map[cells$y[i] + 1L, cells$x[i] + 1L] <- i
  }

  edges <- list()
  for (i in seq_len(nrow(cells))) {
    x <- cells$x[i]
    y <- cells$y[i]
    if (x + 1L < side) {
      nbr <- id_map[y + 1L, x + 2L]
      if (nbr > 0L) edges[[length(edges) + 1L]] <- c(i, nbr)
    }
    if (y + 1L < side) {
      nbr <- id_map[y + 2L, x + 1L]
      if (nbr > 0L) edges[[length(edges) + 1L]] <- c(i, nbr)
    }
  }

  .bind_edges(edges)
}

.sierpinski.carpet.mask <- function() {
  matrix(
    c(
      1, 1, 1,
      1, 0, 1,
      1, 1, 1
    ),
    nrow = 3L,
    byrow = TRUE
  )
}

.vicsek.mask <- function() {
  matrix(
    c(
      0, 1, 0,
      1, 1, 1,
      0, 1, 0
    ),
    nrow = 3L,
    byrow = TRUE
  )
}

.mask.corner.keep <- function(k,
                              width = 2,
                              corner = c("top_left", "top_right",
                                         "bottom_left", "bottom_right")) {
  k <- .as_whole_number(k, "k", min = 2L)
  width <- .as_whole_number(width, "width", min = 1L)
  if (width > k) {
    stop("width must be <= k", call. = FALSE)
  }
  corner <- .as_named_choice(corner[[1L]],
                             c("top_left", "top_right", "bottom_left", "bottom_right"),
                             "corner")
  keep <- matrix(FALSE, nrow = k, ncol = k)
  rows <- switch(corner,
                 top_left = seq_len(width),
                 top_right = seq_len(width),
                 bottom_left = (k - width + 1L):k,
                 bottom_right = (k - width + 1L):k)
  cols <- switch(corner,
                 top_left = seq_len(width),
                 top_right = (k - width + 1L):k,
                 bottom_left = seq_len(width),
                 bottom_right = (k - width + 1L):k)
  keep[rows, cols] <- TRUE
  keep
}

.mask.asymmetric.holes.keep <- function(k, hole_size = 1) {
  k <- .as_odd_whole_number(k, "k", min = 5L)
  hole_size <- .as_whole_number(hole_size, "hole_size", min = 1L)
  max_hole <- max(1L, (k - 3L) %/% 2L)
  if (hole_size > max_hole) {
    stop(sprintf("hole_size must be <= %d for the chosen k", max_hole), call. = FALSE)
  }

  keep <- matrix(TRUE, nrow = k, ncol = k)
  max_start <- k - hole_size
  mid_start <- max(2L, min(max_start, ((k + 1L) %/% 2L) - (hole_size %/% 2L)))
  hole_blocks <- unique(
    list(
      c(2L, 2L),
      c(2L, max_start),
      c(max_start, mid_start)
    )
  )
  for (block in hole_blocks) {
    rows <- block[[1L]]:(block[[1L]] + hole_size - 1L)
    cols <- block[[2L]]:(block[[2L]] + hole_size - 1L)
    if (min(rows) >= 2L && max(rows) <= k - 1L &&
        min(cols) >= 2L && max(cols) <= k - 1L) {
      keep[rows, cols] <- FALSE
    }
  }
  keep
}

.keep.periodic.holes <- function(h,
                                 w = h,
                                 hole_period = 4,
                                 hole_height = 1,
                                 hole_width = hole_height,
                                 row_offset = 2,
                                 col_offset = 2) {
  h <- .as_whole_number(h, "h", min = 1L)
  w <- .as_whole_number(w, "w", min = 1L)
  hole_period <- .as_whole_number(hole_period, "hole_period", min = 2L)
  hole_height <- .as_whole_number(hole_height, "hole_height", min = 1L)
  hole_width <- .as_whole_number(hole_width, "hole_width", min = 1L)
  row_offset <- .as_whole_number(row_offset, "row_offset", min = 1L)
  col_offset <- .as_whole_number(col_offset, "col_offset", min = 1L)

  keep <- matrix(TRUE, nrow = h, ncol = w)
  row_starts <- seq(row_offset, h - hole_height + 1L, by = hole_period)
  col_starts <- seq(col_offset, w - hole_width + 1L, by = hole_period)
  for (r in row_starts) {
    for (c in col_starts) {
      keep[r:(r + hole_height - 1L), c:(c + hole_width - 1L)] <- FALSE
    }
  }
  keep
}

.keep.staggered.windows <- function(h,
                                    w = h,
                                    window_height = 1,
                                    window_width = 2,
                                    row_period = 4,
                                    col_period = 5,
                                    row_offset = 2,
                                    col_offset = 2) {
  h <- .as_whole_number(h, "h", min = 1L)
  w <- .as_whole_number(w, "w", min = 1L)
  window_height <- .as_whole_number(window_height, "window_height", min = 1L)
  window_width <- .as_whole_number(window_width, "window_width", min = 1L)
  row_period <- .as_whole_number(row_period, "row_period", min = 2L)
  col_period <- .as_whole_number(col_period, "col_period", min = 2L)
  row_offset <- .as_whole_number(row_offset, "row_offset", min = 1L)
  col_offset <- .as_whole_number(col_offset, "col_offset", min = 1L)

  keep <- matrix(TRUE, nrow = h, ncol = w)
  row_starts <- seq(row_offset, h - window_height + 1L, by = row_period)
  for (band_idx in seq_along(row_starts)) {
    r <- row_starts[[band_idx]]
    band_shift <- if ((band_idx %% 2L) == 0L) col_period %/% 2L else 0L
    col_starts <- seq(col_offset + band_shift, w - window_width + 1L, by = col_period)
    for (c in col_starts) {
      if (c >= 1L && c + window_width - 1L <= w) {
        keep[r:(r + window_height - 1L), c:(c + window_width - 1L)] <- FALSE
      }
    }
  }
  keep
}

.keep.slit.channels <- function(h,
                                w = h,
                                orientation = c("vertical", "horizontal"),
                                slit_period = 5,
                                slit_width = 1,
                                bridge_spacing = 4,
                                bridge_size = 1,
                                offset = 2) {
  h <- .as_whole_number(h, "h", min = 1L)
  w <- .as_whole_number(w, "w", min = 1L)
  orientation <- match.arg(orientation)
  slit_period <- .as_whole_number(slit_period, "slit_period", min = 2L)
  slit_width <- .as_whole_number(slit_width, "slit_width", min = 1L)
  bridge_spacing <- .as_whole_number(bridge_spacing, "bridge_spacing", min = 2L)
  bridge_size <- .as_whole_number(bridge_size, "bridge_size", min = 1L)
  offset <- .as_whole_number(offset, "offset", min = 1L)

  keep <- matrix(TRUE, nrow = h, ncol = w)
  if (orientation == "vertical") {
    col_starts <- seq(offset, w - slit_width + 1L, by = slit_period)
    for (c in col_starts) {
      for (r in seq_len(h)) {
        in_bridge <- ((r - 1L) %% bridge_spacing) < bridge_size
        if (!in_bridge) {
          keep[r, c:(c + slit_width - 1L)] <- FALSE
        }
      }
    }
  } else {
    row_starts <- seq(offset, h - slit_width + 1L, by = slit_period)
    for (r in row_starts) {
      for (c in seq_len(w)) {
        in_bridge <- ((c - 1L) %% bridge_spacing) < bridge_size
        if (!in_bridge) {
          keep[r:(r + slit_width - 1L), c] <- FALSE
        }
      }
    }
  }
  keep
}

.keep.asymmetric.notches <- function(h,
                                     w = h,
                                     notch_depth = 3,
                                     notch_width = 2) {
  h <- .as_whole_number(h, "h", min = 3L)
  w <- .as_whole_number(w, "w", min = 3L)
  notch_depth <- .as_whole_number(notch_depth, "notch_depth", min = 1L)
  notch_width <- .as_whole_number(notch_width, "notch_width", min = 1L)
  keep <- matrix(TRUE, nrow = h, ncol = w)

  top_depth <- min(notch_depth, max(1L, h - 2L))
  left_depth <- min(notch_depth, max(1L, w - 2L))
  bottom_depth <- min(max(1L, notch_depth - 1L), max(1L, h - 2L))
  width1 <- min(notch_width, max(1L, w - 2L))
  width2 <- min(max(1L, notch_width + 1L), max(1L, h - 2L))
  width3 <- min(max(1L, notch_width), max(1L, w - 3L))

  c1_start <- max(2L, (w %/% 3L))
  c1_end <- min(w - 1L, c1_start + width1 - 1L)
  keep[1L:top_depth, c1_start:c1_end] <- FALSE

  r2_start <- max(2L, (h %/% 2L) - width2 %/% 2L)
  r2_end <- min(h - 1L, r2_start + width2 - 1L)
  keep[r2_start:r2_end, (w - left_depth + 1L):w] <- FALSE

  c3_end <- min(w - 2L, max(2L, (2L * w) %/% 3L))
  c3_start <- max(2L, c3_end - width3 + 1L)
  keep[(h - bottom_depth + 1L):h, c3_start:c3_end] <- FALSE

  keep
}

#' Mask pattern helpers for recursive grid families
#'
#' Convenience constructors for connected \eqn{k \times k} keep-masks that can
#' be passed to \code{\link{edges.recursive.mask.grid}()} or
#' \code{\link{recursive.mask.grid.surface.graph}()} to build generalized mask
#' carpets. These helpers are meant to provide explicit benchmark motifs such
#' as cross, border, corner, and asymmetric-hole patterns.
#'
#' The returned masks use standard matrix display orientation: rows run from top
#' to bottom and columns run from left to right.
#'
#' @param k Mask side length.
#' @param arm_width Width of the retained cross arms.
#' @param thickness Border thickness for \code{mask.border()}.
#' @param width Side length of the retained corner block for
#'   \code{mask.corner()}.
#' @param corner Which corner to retain in \code{mask.corner()}.
#' @param hole_size Side length of each removed interior hole block for
#'   \code{mask.asymmetric.holes()}.
#'
#' @return A logical \eqn{k \times k} keep-mask.
#' @name mask_pattern_helpers
NULL

#' @rdname mask_pattern_helpers
#' @export
mask.cross <- function(k = 5, arm_width = 1) {
  k <- .as_odd_whole_number(k, "k", min = 3L)
  arm_width <- .as_whole_number(arm_width, "arm_width", min = 1L)
  if (arm_width > k) {
    stop("arm_width must be <= k", call. = FALSE)
  }
  start <- ((k - arm_width) %/% 2L) + 1L
  stop_idx <- start + arm_width - 1L
  keep <- matrix(FALSE, nrow = k, ncol = k)
  keep[start:stop_idx, ] <- TRUE
  keep[, start:stop_idx] <- TRUE
  keep
}

#' @rdname mask_pattern_helpers
#' @export
mask.border <- function(k = 5, thickness = 1) {
  k <- .as_whole_number(k, "k", min = 3L)
  thickness <- .as_whole_number(thickness, "thickness", min = 1L)
  if (thickness > k) {
    stop("thickness must be <= k", call. = FALSE)
  }
  keep <- matrix(TRUE, nrow = k, ncol = k)
  if ((2L * thickness) < k) {
    keep[(thickness + 1L):(k - thickness), (thickness + 1L):(k - thickness)] <- FALSE
  }
  keep
}

#' @rdname mask_pattern_helpers
#' @export
mask.corner <- function(k = 5,
                        width = 2,
                        corner = c("top_left", "top_right", "bottom_left", "bottom_right")) {
  .mask.corner.keep(
    k = k,
    width = width,
    corner = match.arg(corner)
  )
}

#' @rdname mask_pattern_helpers
#' @export
mask.asymmetric.holes <- function(k = 5, hole_size = 1) {
  .mask.asymmetric.holes.keep(
    k = k,
    hole_size = hole_size
  )
}

#' Occupied-mesh and perforated-grid helpers
#'
#' \code{edges.occupied.mesh()} builds the orthogonal adjacency graph on the
#' occupied cells of a finite rectangular grid. The corresponding surface
#' helpers lift those occupied cells into \eqn{\mathbb{R}^3} using the same
#' \code{"saddle"}, \code{"paraboloid"}, and \code{"ripple"} families used
#' for regular meshes.
#'
#' @param keep Logical or numeric occupancy matrix. Non-zero entries are kept.
#' @param surface Surface family used for the lift. One of \code{"saddle"},
#'   \code{"paraboloid"}, or \code{"ripple"}.
#' @param amplitude Finite numeric amplitude controlling the non-flat
#'   displacement.
#' @param freq_u Positive ripple frequency in the horizontal parameter
#'   direction. Used only when \code{surface = "ripple"}.
#' @param freq_v Positive ripple frequency in the vertical parameter direction.
#'   Used only when \code{surface = "ripple"}.
#' @param x_scale Positive horizontal scaling of the parameter domain.
#' @param y_scale Positive vertical scaling of the parameter domain.
#' @param normalize Normalization applied to the induced edge lengths. One of
#'   \code{"median"}, \code{"mean"}, or \code{"none"}.
#'
#' @return
#' \code{occupied.mesh.surface.embedding()} returns an \code{n x 3} numeric
#' matrix with columns \code{x}, \code{y}, and \code{z}.
#'
#' \code{occupied.mesh.surface.graph()} returns a list with components:
#' \itemize{
#'   \item \code{edges}: the occupied-cell mesh edges,
#'   \item \code{n}: number of vertices,
#'   \item \code{edge_weights}: induced positive edge lengths,
#'   \item \code{coords_surface}: the 3D surface embedding,
#'   \item \code{coords_param}: the 2D parameter coordinates,
#'   \item \code{weight_scale}: the normalization constant applied to the raw
#'     edge lengths,
#'   \item \code{family}: always \code{"occupied.mesh"},
#'   \item \code{surface}: the chosen surface name,
#'   \item \code{keep}: the logical occupancy matrix,
#'   \item \code{label}: a human-readable family label.
#' }
#'
#' @name occupied_mesh_surface_helpers
NULL

#' @rdname occupied_mesh_surface_helpers
#' @export
occupied.mesh.surface.embedding <- function(keep,
                                            surface = c("saddle", "paraboloid", "ripple"),
                                            amplitude = 0.75,
                                            freq_u = 1,
                                            freq_v = 1,
                                            x_scale = 1,
                                            y_scale = 1) {
  keep <- .as_keep_grid(keep, "keep", min_h = 1L, min_w = 1L)
  surface <- match.arg(surface)
  amplitude <- .as_finite_scalar(amplitude, "amplitude")
  freq_u <- .as_positive_scalar(freq_u, "freq_u")
  freq_v <- .as_positive_scalar(freq_v, "freq_v")

  coords_param <- .occupied.mesh.param.coords(
    keep = keep,
    x_scale = x_scale,
    y_scale = y_scale
  )
  u <- coords_param[, 1L]
  v <- coords_param[, 2L]
  z <- .surface.z.from_uv(
    u = u,
    v = v,
    surface = surface,
    amplitude = amplitude,
    freq_u = freq_u,
    freq_v = freq_v
  )

  coords <- cbind(x = u, y = v, z = z)
  storage.mode(coords) <- "double"
  coords
}

#' @rdname occupied_mesh_surface_helpers
#' @export
occupied.mesh.surface.graph <- function(keep,
                                        surface = c("saddle", "paraboloid", "ripple"),
                                        amplitude = 0.75,
                                        freq_u = 1,
                                        freq_v = 1,
                                        x_scale = 1,
                                        y_scale = 1,
                                        normalize = c("median", "mean", "none")) {
  keep <- .as_keep_grid(keep, "keep", min_h = 1L, min_w = 1L)
  surface <- match.arg(surface)
  normalize <- match.arg(normalize)

  edges <- .occupied.mesh.edges(keep)
  coords_param <- .occupied.mesh.param.coords(
    keep = keep,
    x_scale = x_scale,
    y_scale = y_scale
  )
  coords_surface <- occupied.mesh.surface.embedding(
    keep = keep,
    surface = surface,
    amplitude = amplitude,
    freq_u = freq_u,
    freq_v = freq_v,
    x_scale = x_scale,
    y_scale = y_scale
  )
  weights <- .edge.weights.from.embedding(
    edges = edges,
    coords = coords_surface,
    normalize = normalize
  )

  out <- list(
    edges = edges,
    n = as.integer(nrow(coords_surface)),
    edge_weights = weights$edge_weights,
    coords_surface = coords_surface,
    coords_param = coords_param,
    weight_scale = weights$weight_scale,
    family = "occupied.mesh",
    surface = surface,
    keep = keep,
    normalize = normalize,
    label = sprintf("%s occupied mesh %dx%d",
                    tools::toTitleCase(surface),
                    nrow(keep),
                    ncol(keep))
  )
  class(out) <- c("grip_occupied_mesh_surface_graph", "list")
  out
}

#' Deterministic perforated-grid occupancy helpers
#'
#' Convenience constructors for finite occupied grids with non-random holes or
#' channels. These matrices can be passed to \code{\link{edges.occupied.mesh}()}
#' or \code{\link{occupied.mesh.surface.graph}()}.
#'
#' @param h Number of rows.
#' @param w Number of columns. Defaults to \code{h}.
#' @param hole_period Spacing between periodic holes.
#' @param hole_height Height of each removed rectangular hole.
#' @param hole_width Width of each removed rectangular hole.
#' @param row_offset Starting row index for the first removed block.
#' @param col_offset Starting column index for the first removed block.
#' @param window_height Height of each staggered removed window.
#' @param window_width Width of each staggered removed window.
#' @param row_period Vertical spacing between staggered window bands.
#' @param col_period Horizontal spacing between staggered windows within a band.
#' @param orientation Whether slit channels run \code{"vertical"} or
#'   \code{"horizontal"}.
#' @param slit_period Spacing between repeated slit channels.
#' @param slit_width Width of each slit channel.
#' @param bridge_spacing Spacing between preserved bridge segments.
#' @param bridge_size Size of each preserved bridge segment along a slit.
#' @param offset Starting row or column index for the first slit channel.
#' @param notch_depth Depth of each asymmetric notch cut from a boundary.
#' @param notch_width Width of the notched opening along that boundary.
#'
#' @return A logical occupancy matrix whose \code{TRUE} entries represent
#'   retained cells.
#' @name perforated_grid_helpers
NULL

#' @rdname perforated_grid_helpers
#' @export
keep.periodic.holes <- function(h,
                                w = h,
                                hole_period = 4,
                                hole_height = 1,
                                hole_width = hole_height,
                                row_offset = 2,
                                col_offset = 2) {
  .keep.periodic.holes(
    h = h,
    w = w,
    hole_period = hole_period,
    hole_height = hole_height,
    hole_width = hole_width,
    row_offset = row_offset,
    col_offset = col_offset
  )
}

#' @rdname perforated_grid_helpers
#' @export
keep.staggered.windows <- function(h,
                                   w = h,
                                   window_height = 1,
                                   window_width = 2,
                                   row_period = 4,
                                   col_period = 5,
                                   row_offset = 2,
                                   col_offset = 2) {
  .keep.staggered.windows(
    h = h,
    w = w,
    window_height = window_height,
    window_width = window_width,
    row_period = row_period,
    col_period = col_period,
    row_offset = row_offset,
    col_offset = col_offset
  )
}

#' @rdname perforated_grid_helpers
#' @export
keep.slit.channels <- function(h,
                               w = h,
                               orientation = c("vertical", "horizontal"),
                               slit_period = 5,
                               slit_width = 1,
                               bridge_spacing = 4,
                               bridge_size = 1,
                               offset = 2) {
  .keep.slit.channels(
    h = h,
    w = w,
    orientation = match.arg(orientation),
    slit_period = slit_period,
    slit_width = slit_width,
    bridge_spacing = bridge_spacing,
    bridge_size = bridge_size,
    offset = offset
  )
}

#' @rdname perforated_grid_helpers
#' @export
keep.asymmetric.notches <- function(h,
                                    w = h,
                                    notch_depth = 3,
                                    notch_width = 2) {
  .keep.asymmetric.notches(
    h = h,
    w = w,
    notch_depth = notch_depth,
    notch_width = notch_width
  )
}

.edge.weights.from.embedding <- function(edges,
                                         coords,
                                         normalize = c("median", "mean", "none")) {
  normalize <- match.arg(normalize)
  if (!is.matrix(coords) || !is.numeric(coords) || ncol(coords) < 2L) {
    stop("coords must be a numeric matrix with at least 2 columns", call. = FALSE)
  }
  if (!is.matrix(edges) || ncol(edges) != 2L) {
    stop("edges must be a two-column matrix", call. = FALSE)
  }
  if (nrow(edges) == 0L) {
    return(list(edge_weights = numeric(0L), weight_scale = 1))
  }

  diffs <- coords[edges[, 1L], , drop = FALSE] - coords[edges[, 2L], , drop = FALSE]
  edge_weights <- sqrt(rowSums(diffs^2))
  if (any(!is.finite(edge_weights) | edge_weights <= 0)) {
    stop("embedded edge weights must all be finite and > 0", call. = FALSE)
  }

  weight_scale <- switch(
    normalize,
    median = stats::median(edge_weights),
    mean = mean(edge_weights),
    none = 1
  )
  if (!is.finite(weight_scale) || weight_scale <= 0) {
    weight_scale <- 1
  }

  list(
    edge_weights = as.double(edge_weights / weight_scale),
    weight_scale = as.double(weight_scale)
  )
}

#' Sample graph generators
#'
#' Convenience helpers that build small undirected graph families as two-column
#' integer edge matrices suitable for \code{\link{grip.layout}()}. These
#' helpers are meant for examples, experiments, and reproducible tests.
#'
#' The occupied-grid, recursive masked-grid, and Sierpinski families are
#' exposed explicitly rather than overloading a single generator with
#' layout-dimension-dependent behavior: \code{\link{edges.occupied.mesh}()}
#' builds a finite perforated-mesh family from an occupancy matrix,
#' \code{\link{edges.recursive.mask.grid}()} builds a generic square-mask
#' family, \code{\link{edges.vicsek}()} builds the connected axial-cross
#' variant, \code{\link{edges.sierpinski.triangle}()} builds the 2-simplex
#' family, \code{\link{edges.sierpinski.tetrahedron}()} builds the 3-simplex
#' family, and \code{\link{edges.sierpinski.carpet}()} builds a 2D
#' cell-adjacency carpet graph.
#'
#' @return A two-column integer matrix of undirected edges. Vertex labels are
#'   consecutive integers starting at 1.
#' @name graph_generators
NULL

#' @describeIn graph_generators Path graph on \code{n} vertices.
#' @param n Number of vertices.
#' @examples
#' edges <- edges.path(6)
#' coords <- grip.layout(edges, n = 6, dim = 2, seed = 1)
#' grip.plot(coords, edges, main = "Path graph", pch = 16, cex = 0.8)
#' @export
edges.path <- function(n) {
  n <- .as_whole_number(n, "n")
  if (n < 2L) {
    return(.empty_edge_matrix())
  }
  cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L)
}

#' @describeIn graph_generators Cycle graph on \code{n} vertices.
#' @export
edges.cycle <- function(n) {
  n <- .as_whole_number(n, "n")
  if (n < 2L) {
    return(.empty_edge_matrix())
  }
  .normalize_undirected_edges(rbind(edges.path(n), c(n, 1L)))
}

#' @describeIn graph_generators Rectangular grid graph with \code{h} rows and \code{w}
#'   columns.
#' @param h Number of rows.
#' @param w Number of columns. Defaults to \code{h}.
#' @param keep Logical or numeric occupancy matrix. Non-zero entries are kept.
#' @export
edges.mesh <- function(h, w = h) {
  h <- .as_whole_number(h, "h", min = 1L)
  w <- .as_whole_number(w, "w", min = 1L)
  idx <- function(i, j) (i - 1L) * w + j
  edges <- list()
  for (i in seq_len(h)) {
    for (j in seq_len(w)) {
      v <- idx(i, j)
      if (i < h) edges[[length(edges) + 1L]] <- c(v, idx(i + 1L, j))
      if (j < w) edges[[length(edges) + 1L]] <- c(v, idx(i, j + 1L))
    }
  }
  .normalize_undirected_edges(.bind_edges(edges))
}

#' @describeIn graph_generators Rectangular occupied-grid graph whose vertices
#'   are kept cells and whose edges connect orthogonally adjacent kept cells.
#' @export
edges.occupied.mesh <- function(keep) {
  .occupied.mesh.edges(keep)
}

#' Weighted mesh surface helpers
#'
#' Convenience helpers that lift a rectangular mesh parameter grid into
#' \eqn{\mathbb{R}^3} and use the induced Euclidean edge lengths as positive
#' graph weights. These helpers are intended for benchmark families where the
#' topology is a plain mesh but the intended metric comes from a curved ambient
#' geometry.
#'
#' `mesh.surface.embedding()` returns the 3D coordinates of the lifted grid.
#' `mesh.surface.graph()` returns a reusable weighted-graph bundle containing the
#' mesh edges, induced edge weights, the 3D surface coordinates, and the 2D
#' parameter coordinates.
#'
#' @param h Number of rows.
#' @param w Number of columns. Defaults to \code{h}.
#' @param surface Surface family used for the lift. One of \code{"saddle"},
#'   \code{"paraboloid"}, or \code{"ripple"}.
#' @param amplitude Finite numeric amplitude controlling the non-flat
#'   displacement.
#' @param freq_u Positive ripple frequency in the horizontal parameter
#'   direction. Used only when \code{surface = "ripple"}.
#' @param freq_v Positive ripple frequency in the vertical parameter direction.
#'   Used only when \code{surface = "ripple"}.
#' @param x_scale Positive horizontal scaling of the parameter domain.
#' @param y_scale Positive vertical scaling of the parameter domain.
#' @param normalize Normalization applied to the induced edge lengths. One of
#'   \code{"median"}, \code{"mean"}, or \code{"none"}.
#'
#' @return
#' \code{mesh.surface.embedding()} returns an \code{n x 3} numeric matrix with
#' columns \code{x}, \code{y}, and \code{z}.
#'
#' \code{mesh.surface.graph()} returns a list with components:
#' \itemize{
#'   \item \code{edges}: the undirected mesh edges,
#'   \item \code{n}: number of vertices,
#'   \item \code{edge_weights}: induced positive edge lengths,
#'   \item \code{coords_surface}: the 3D surface embedding,
#'   \item \code{coords_param}: the 2D parameter-grid coordinates,
#'   \item \code{weight_scale}: the normalization constant applied to the raw
#'     edge lengths,
#'   \item \code{family}: always \code{"mesh"},
#'   \item \code{surface}: the chosen surface name,
#'   \item \code{label}: a human-readable family label.
#' }
#'
#' @name mesh_surface_helpers
NULL

#' @rdname mesh_surface_helpers
#' @export
mesh.surface.embedding <- function(h,
                                   w = h,
                                   surface = c("saddle", "paraboloid", "ripple"),
                                   amplitude = 0.75,
                                   freq_u = 1,
                                   freq_v = 1,
                                   x_scale = 1,
                                   y_scale = 1) {
  h <- .as_whole_number(h, "h", min = 1L)
  w <- .as_whole_number(w, "w", min = 1L)
  surface <- match.arg(surface)
  amplitude <- .as_finite_scalar(amplitude, "amplitude")
  freq_u <- .as_positive_scalar(freq_u, "freq_u")
  freq_v <- .as_positive_scalar(freq_v, "freq_v")

  coords_param <- .mesh.param.coords(
    h = h,
    w = w,
    x_scale = x_scale,
    y_scale = y_scale
  )
  u <- coords_param[, 1L]
  v <- coords_param[, 2L]
  z <- switch(
    surface,
    saddle = amplitude * (u^2 - v^2),
    paraboloid = amplitude * (u^2 + v^2),
    ripple = amplitude * sin(pi * freq_u * u) * cos(pi * freq_v * v)
  )

  coords <- cbind(x = u, y = v, z = z)
  storage.mode(coords) <- "double"
  coords
}

#' @rdname mesh_surface_helpers
#' @export
mesh.surface.graph <- function(h,
                               w = h,
                               surface = c("saddle", "paraboloid", "ripple"),
                               amplitude = 0.75,
                               freq_u = 1,
                               freq_v = 1,
                               x_scale = 1,
                               y_scale = 1,
                               normalize = c("median", "mean", "none")) {
  h <- .as_whole_number(h, "h", min = 1L)
  w <- .as_whole_number(w, "w", min = 1L)
  surface <- match.arg(surface)
  normalize <- match.arg(normalize)

  edges <- edges.mesh(h, w)
  coords_param <- .mesh.param.coords(
    h = h,
    w = w,
    x_scale = x_scale,
    y_scale = y_scale
  )
  coords_surface <- mesh.surface.embedding(
    h = h,
    w = w,
    surface = surface,
    amplitude = amplitude,
    freq_u = freq_u,
    freq_v = freq_v,
    x_scale = x_scale,
    y_scale = y_scale
  )
  weights <- .edge.weights.from.embedding(
    edges = edges,
    coords = coords_surface,
    normalize = normalize
  )

  out <- list(
    edges = edges,
    n = as.integer(h * w),
    edge_weights = weights$edge_weights,
    coords_surface = coords_surface,
    coords_param = coords_param,
    weight_scale = weights$weight_scale,
    family = "mesh",
    surface = surface,
    normalize = normalize,
    label = sprintf("%s mesh %dx%d", tools::toTitleCase(surface), h, w)
  )
  class(out) <- c("grip_mesh_surface_graph", "list")
  out
}

#' Weighted cylinder surface helpers
#'
#' Convenience helpers that embed a cylindrical grid graph into
#' \eqn{\mathbb{R}^3} and use the induced Euclidean edge lengths as positive
#' graph weights. These helpers are intended for benchmark families where the
#' graph topology is cylindrical but the intended metric comes from a curved or
#' spatially varying 3D realization.
#'
#' `cylinder.surface.embedding()` returns the 3D coordinates of the embedded
#' cylindrical grid. `cylinder.surface.graph()` returns a reusable weighted-graph
#' bundle containing the cylinder edges, induced edge weights, the 3D surface
#' coordinates, and a 2D unwrapped parameterization.
#'
#' @param h Number of rows.
#' @param w Number of columns. Defaults to \code{h}.
#' @param surface Cylinder surface family. One of \code{"barrel"},
#'   \code{"hourglass"}, or \code{"wavy"}.
#' @param radius Positive baseline cylinder radius.
#' @param height Positive cylinder height.
#' @param amplitude Finite numeric deformation amplitude. The resulting radius
#'   profile must stay positive everywhere.
#' @param freq_theta Positive angular frequency used only when
#'   \code{surface = "wavy"}.
#' @param freq_z Positive vertical frequency used only when
#'   \code{surface = "wavy"}.
#' @param twist Finite angular twist applied linearly with height.
#' @param normalize Normalization applied to the induced edge lengths. One of
#'   \code{"median"}, \code{"mean"}, or \code{"none"}.
#'
#' @return
#' \code{cylinder.surface.embedding()} returns an \code{n x 3} numeric matrix
#' with columns \code{x}, \code{y}, and \code{z}.
#'
#' \code{cylinder.surface.graph()} returns a list with components:
#' \itemize{
#'   \item \code{edges}: the undirected cylindrical-grid edges,
#'   \item \code{n}: number of vertices,
#'   \item \code{edge_weights}: induced positive edge lengths,
#'   \item \code{coords_surface}: the 3D surface embedding,
#'   \item \code{coords_param}: the 2D unwrapped parameter coordinates,
#'   \item \code{weight_scale}: the normalization constant applied to the raw
#'     edge lengths,
#'   \item \code{family}: always \code{"cylinder"},
#'   \item \code{surface}: the chosen surface name,
#'   \item \code{label}: a human-readable family label.
#' }
#'
#' @name cylinder_surface_helpers
NULL

#' @rdname cylinder_surface_helpers
#' @export
cylinder.surface.embedding <- function(h,
                                       w = h,
                                       surface = c("barrel", "hourglass", "wavy"),
                                       radius = 1,
                                       height = 2,
                                       amplitude = 0.3,
                                       freq_theta = 2,
                                       freq_z = 1,
                                       twist = 0.25) {
  h <- .as_whole_number(h, "h", min = 1L)
  w <- .as_whole_number(w, "w", min = 3L)
  surface <- match.arg(surface)
  radius <- .as_positive_scalar(radius, "radius")
  height <- .as_positive_scalar(height, "height")
  amplitude <- .as_finite_scalar(amplitude, "amplitude")
  freq_theta <- .as_positive_scalar(freq_theta, "freq_theta")
  freq_z <- .as_positive_scalar(freq_z, "freq_z")
  twist <- .as_finite_scalar(twist, "twist")

  theta_vals <- seq(0, 2 * pi * (1 - 1 / w), length.out = w)
  s_vals <- rev(.grid_axis(h))
  coords <- matrix(0, nrow = h * w, ncol = 3L)
  for (i in seq_len(h)) {
    s <- s_vals[[i]]
    z <- (height / 2) * s
    for (j in seq_len(w)) {
      theta <- theta_vals[[j]]
      local_radius <- switch(
        surface,
        barrel = radius * (1 + amplitude * (1 - s^2)),
        hourglass = radius * (1 - amplitude * (1 - s^2)),
        wavy = radius * (1 + amplitude * sin(freq_theta * theta) * cos(pi * freq_z * s))
      )
      if (!is.finite(local_radius) || local_radius <= 0) {
        stop("cylinder surface parameters produce a non-positive radius; adjust amplitude or surface settings",
             call. = FALSE)
      }
      angle <- theta + twist * s
      id <- (i - 1L) * w + j
      coords[id, ] <- c(
        local_radius * cos(angle),
        local_radius * sin(angle),
        z
      )
    }
  }
  colnames(coords) <- c("x", "y", "z")
  storage.mode(coords) <- "double"
  coords
}

#' @rdname cylinder_surface_helpers
#' @export
cylinder.surface.graph <- function(h,
                                   w = h,
                                   surface = c("barrel", "hourglass", "wavy"),
                                   radius = 1,
                                   height = 2,
                                   amplitude = 0.3,
                                   freq_theta = 2,
                                   freq_z = 1,
                                   twist = 0.25,
                                   normalize = c("median", "mean", "none")) {
  h <- .as_whole_number(h, "h", min = 1L)
  w <- .as_whole_number(w, "w", min = 3L)
  surface <- match.arg(surface)
  normalize <- match.arg(normalize)

  edges <- edges.cylinder(h, w)
  coords_param <- .cylinder.param.coords(
    h = h,
    w = w,
    radius = radius,
    height = height
  )
  coords_surface <- cylinder.surface.embedding(
    h = h,
    w = w,
    surface = surface,
    radius = radius,
    height = height,
    amplitude = amplitude,
    freq_theta = freq_theta,
    freq_z = freq_z,
    twist = twist
  )
  weights <- .edge.weights.from.embedding(
    edges = edges,
    coords = coords_surface,
    normalize = normalize
  )

  out <- list(
    edges = edges,
    n = as.integer(h * w),
    edge_weights = weights$edge_weights,
    coords_surface = coords_surface,
    coords_param = coords_param,
    weight_scale = weights$weight_scale,
    family = "cylinder",
    surface = surface,
    normalize = normalize,
    label = sprintf("%s cylinder %dx%d", tools::toTitleCase(surface), h, w)
  )
  class(out) <- c("grip_cylinder_surface_graph", "list")
  out
}

#' Weighted torus surface helpers
#'
#' Convenience helpers that embed a toroidal grid graph into
#' \eqn{\mathbb{R}^3} and use the induced Euclidean edge lengths as positive
#' graph weights. These helpers are intended for benchmark families where the
#' graph topology is toroidal but the intended metric comes from a curved or
#' spatially varying 3D realization.
#'
#' `torus.surface.embedding()` returns the 3D coordinates of the embedded
#' toroidal grid. `torus.surface.graph()` returns a reusable weighted-graph
#' bundle containing the torus edges, induced edge weights, the 3D surface
#' coordinates, and a 2D unwrapped parameterization.
#'
#' @param h Number of wrapped rows. Must be at least \code{3}.
#' @param w Number of wrapped columns. Defaults to \code{h}; must be at least
#'   \code{3}.
#' @param surface Torus surface family. One of \code{"standard"},
#'   \code{"pinched"}, or \code{"wavy"}.
#' @param major_radius Positive distance from the torus center to the center of
#'   the tube.
#' @param minor_radius Positive baseline radius of the torus tube. The local
#'   tube radius must remain strictly smaller than \code{major_radius}
#'   everywhere.
#' @param amplitude Finite numeric deformation amplitude.
#' @param freq_major Positive angular frequency around the major cycle, used by
#'   \code{"wavy"} and the periodic twist.
#' @param freq_minor Positive angular frequency around the minor cycle, used by
#'   \code{"wavy"}.
#' @param twist Finite phase twist applied periodically to the minor angle as a
#'   function of the major angle.
#' @param normalize Normalization applied to the induced edge lengths. One of
#'   \code{"median"}, \code{"mean"}, or \code{"none"}.
#'
#' @return
#' \code{torus.surface.embedding()} returns an \code{n x 3} numeric matrix with
#' columns \code{x}, \code{y}, and \code{z}.
#'
#' \code{torus.surface.graph()} returns a list with components:
#' \itemize{
#'   \item \code{edges}: the undirected toroidal-grid edges,
#'   \item \code{n}: number of vertices,
#'   \item \code{edge_weights}: induced positive edge lengths,
#'   \item \code{coords_surface}: the 3D surface embedding,
#'   \item \code{coords_param}: the 2D unwrapped parameter coordinates,
#'   \item \code{weight_scale}: the normalization constant applied to the raw
#'     edge lengths,
#'   \item \code{family}: always \code{"torus"},
#'   \item \code{surface}: the chosen surface name,
#'   \item \code{label}: a human-readable family label.
#' }
#'
#' @name torus_surface_helpers
NULL

#' @rdname torus_surface_helpers
#' @export
torus.surface.embedding <- function(h,
                                    w = h,
                                    surface = c("standard", "pinched", "wavy"),
                                    major_radius = 2,
                                    minor_radius = 0.75,
                                    amplitude = 0.2,
                                    freq_major = 2,
                                    freq_minor = 1,
                                    twist = 0.25) {
  h <- .as_whole_number(h, "h", min = 3L)
  w <- .as_whole_number(w, "w", min = 3L)
  surface <- match.arg(surface)
  major_radius <- .as_positive_scalar(major_radius, "major_radius")
  minor_radius <- .as_positive_scalar(minor_radius, "minor_radius")
  amplitude <- .as_finite_scalar(amplitude, "amplitude")
  freq_major <- .as_positive_scalar(freq_major, "freq_major")
  freq_minor <- .as_positive_scalar(freq_minor, "freq_minor")
  twist <- .as_finite_scalar(twist, "twist")

  theta_vals <- seq(0, 2 * pi * (1 - 1 / w), length.out = w)
  phi_vals <- seq(0, 2 * pi * (1 - 1 / h), length.out = h)
  coords <- matrix(0, nrow = h * w, ncol = 3L)
  for (i in seq_len(h)) {
    phi <- phi_vals[[i]]
    for (j in seq_len(w)) {
      theta <- theta_vals[[j]]
      local_minor_radius <- switch(
        surface,
        standard = minor_radius,
        pinched = minor_radius * (1 - amplitude * cos(theta)),
        wavy = minor_radius * (1 + amplitude * sin(freq_major * theta) * cos(freq_minor * phi))
      )
      if (!is.finite(local_minor_radius) || local_minor_radius <= 0) {
        stop("torus surface parameters produce a non-positive tube radius; adjust amplitude or surface settings",
             call. = FALSE)
      }
      if (local_minor_radius >= major_radius) {
        stop("torus surface parameters require major_radius to exceed the local tube radius everywhere",
             call. = FALSE)
      }
      phi_eff <- phi + twist * sin(freq_major * theta)
      ring_radius <- major_radius + local_minor_radius * cos(phi_eff)
      id <- (i - 1L) * w + j
      coords[id, ] <- c(
        ring_radius * cos(theta),
        ring_radius * sin(theta),
        local_minor_radius * sin(phi_eff)
      )
    }
  }
  colnames(coords) <- c("x", "y", "z")
  storage.mode(coords) <- "double"
  coords
}

#' @rdname torus_surface_helpers
#' @export
torus.surface.graph <- function(h,
                                w = h,
                                surface = c("standard", "pinched", "wavy"),
                                major_radius = 2,
                                minor_radius = 0.75,
                                amplitude = 0.2,
                                freq_major = 2,
                                freq_minor = 1,
                                twist = 0.25,
                                normalize = c("median", "mean", "none")) {
  h <- .as_whole_number(h, "h", min = 3L)
  w <- .as_whole_number(w, "w", min = 3L)
  surface <- match.arg(surface)
  normalize <- match.arg(normalize)

  edges <- edges.torus(h, w)
  coords_param <- .torus.param.coords(
    h = h,
    w = w,
    major_radius = major_radius,
    minor_radius = minor_radius
  )
  coords_surface <- torus.surface.embedding(
    h = h,
    w = w,
    surface = surface,
    major_radius = major_radius,
    minor_radius = minor_radius,
    amplitude = amplitude,
    freq_major = freq_major,
    freq_minor = freq_minor,
    twist = twist
  )
  weights <- .edge.weights.from.embedding(
    edges = edges,
    coords = coords_surface,
    normalize = normalize
  )

  out <- list(
    edges = edges,
    n = as.integer(h * w),
    edge_weights = weights$edge_weights,
    coords_surface = coords_surface,
    coords_param = coords_param,
    weight_scale = weights$weight_scale,
    family = "torus",
    surface = surface,
    normalize = normalize,
    label = sprintf("%s torus %dx%d", tools::toTitleCase(surface), h, w)
  )
  class(out) <- c("grip_torus_surface_graph", "list")
  out
}

#' Weighted sphere surface helpers
#'
#' Convenience helpers that embed the sampled sphere graph into
#' \eqn{\mathbb{R}^3} and use the induced Euclidean edge lengths as positive
#' graph weights. These helpers are intended for benchmark families where the
#' graph topology follows the current pole-plus-latitude-rings sphere graph but
#' the intended metric comes from a curved or spatially varying 3D realization.
#'
#' `sphere.surface.embedding()` returns the 3D coordinates of the sampled
#' surface in the same vertex order as \code{edges.sphere()}: north pole first,
#' then latitude rings from north to south, then the south pole.
#' `sphere.surface.graph()` returns a reusable weighted-graph bundle containing
#' the sphere edges, induced edge weights, the 3D surface coordinates, and a 2D
#' longitude-latitude parameterization.
#'
#' @param h Number of latitude levels including the poles. Must be at least
#'   \code{3}.
#' @param w Wrapped longitude count. Defaults to \code{h}; must be at least
#'   \code{3}.
#' @param surface Sphere surface family. One of \code{"standard"},
#'   \code{"ellipsoid"}, or \code{"wavy"}.
#' @param radius Positive baseline radius.
#' @param amplitude Finite deformation amplitude. For \code{"ellipsoid"},
#'   positive values make the shape oblate and negative values make it prolate,
#'   while keeping all axis lengths positive. For \code{"wavy"}, the local
#'   radius must remain positive everywhere.
#' @param freq_theta Positive longitudinal frequency used only when
#'   \code{surface = "wavy"}.
#' @param freq_lat Positive latitudinal frequency used only when
#'   \code{surface = "wavy"}.
#' @param twist Finite longitude twist applied smoothly by latitude. The twist
#'   vanishes at the poles.
#' @param normalize Normalization applied to the induced edge lengths. One of
#'   \code{"median"}, \code{"mean"}, or \code{"none"}.
#'
#' @return
#' \code{sphere.surface.embedding()} returns an \code{n x 3} numeric matrix
#' with columns \code{x}, \code{y}, and \code{z}, where
#' \code{n = 2 + (h - 2) * w}.
#'
#' \code{sphere.surface.graph()} returns a list with components:
#' \itemize{
#'   \item \code{edges}: the undirected sphere-graph edges,
#'   \item \code{n}: number of vertices,
#'   \item \code{edge_weights}: induced positive edge lengths,
#'   \item \code{coords_surface}: the 3D surface embedding,
#'   \item \code{coords_param}: the 2D longitude-latitude parameter coordinates,
#'   \item \code{weight_scale}: the normalization constant applied to the raw
#'     edge lengths,
#'   \item \code{family}: always \code{"sphere"},
#'   \item \code{surface}: the chosen surface name,
#'   \item \code{label}: a human-readable family label.
#' }
#'
#' @name sphere_surface_helpers
NULL

#' @rdname sphere_surface_helpers
#' @export
sphere.surface.embedding <- function(h,
                                     w = h,
                                     surface = c("standard", "ellipsoid", "wavy"),
                                     radius = 1,
                                     amplitude = 0.2,
                                     freq_theta = 3,
                                     freq_lat = 2,
                                     twist = 0.25) {
  h <- .as_whole_number(h, "h", min = 3L)
  w <- .as_whole_number(w, "w", min = 3L)
  surface <- match.arg(surface)
  radius <- .as_positive_scalar(radius, "radius")
  amplitude <- .as_finite_scalar(amplitude, "amplitude")
  freq_theta <- .as_positive_scalar(freq_theta, "freq_theta")
  freq_lat <- .as_positive_scalar(freq_lat, "freq_lat")
  twist <- .as_finite_scalar(twist, "twist")

  ring.count <- h - 2L
  n <- 2L + ring.count * w
  theta_vals <- seq(0, 2 * pi * (1 - 1 / w), length.out = w)
  lat_vals <- seq(pi / 2, -pi / 2, length.out = h)
  coords <- matrix(0, nrow = n, ncol = 3L)

  if (surface == "ellipsoid") {
    equatorial_radius <- radius * (1 + amplitude)
    polar_radius <- radius * (1 - amplitude)
    if (!is.finite(equatorial_radius) || !is.finite(polar_radius) ||
        equatorial_radius <= 0 || polar_radius <= 0) {
      stop("sphere ellipsoid parameters require positive equatorial and polar radii; adjust amplitude",
           call. = FALSE)
    }
    coords[1L, ] <- c(0, 0, polar_radius)
    if (ring.count > 0L) {
      for (i in seq_len(ring.count)) {
        lat <- lat_vals[[i + 1L]]
        theta_shift <- twist * sin(2 * lat)
        for (j in seq_len(w)) {
          theta_eff <- theta_vals[[j]] + theta_shift
          id <- 1L + (i - 1L) * w + j
          coords[id, ] <- c(
            equatorial_radius * cos(lat) * cos(theta_eff),
            equatorial_radius * cos(lat) * sin(theta_eff),
            polar_radius * sin(lat)
          )
        }
      }
    }
    coords[n, ] <- c(0, 0, -polar_radius)
  } else {
    coords[1L, ] <- c(0, 0, radius)
    if (ring.count > 0L) {
      for (i in seq_len(ring.count)) {
        lat <- lat_vals[[i + 1L]]
        theta_shift <- twist * sin(2 * lat)
        pole_taper <- cos(lat)^2
        for (j in seq_len(w)) {
          theta <- theta_vals[[j]]
          theta_eff <- theta + theta_shift
          local_radius <- switch(
            surface,
            standard = radius,
            wavy = radius * (1 + amplitude * pole_taper *
                               sin(freq_theta * theta) * cos(freq_lat * lat))
          )
          if (!is.finite(local_radius) || local_radius <= 0) {
            stop("sphere surface parameters produce a non-positive local radius; adjust amplitude or surface settings",
                 call. = FALSE)
          }
          id <- 1L + (i - 1L) * w + j
          coords[id, ] <- c(
            local_radius * cos(lat) * cos(theta_eff),
            local_radius * cos(lat) * sin(theta_eff),
            local_radius * sin(lat)
          )
        }
      }
    }
    coords[n, ] <- c(0, 0, -radius)
  }

  colnames(coords) <- c("x", "y", "z")
  storage.mode(coords) <- "double"
  coords
}

#' @rdname sphere_surface_helpers
#' @export
sphere.surface.graph <- function(h,
                                 w = h,
                                 surface = c("standard", "ellipsoid", "wavy"),
                                 radius = 1,
                                 amplitude = 0.2,
                                 freq_theta = 3,
                                 freq_lat = 2,
                                 twist = 0.25,
                                 normalize = c("median", "mean", "none")) {
  h <- .as_whole_number(h, "h", min = 3L)
  w <- .as_whole_number(w, "w", min = 3L)
  surface <- match.arg(surface)
  normalize <- match.arg(normalize)

  edges <- edges.sphere(h, w)
  coords_param <- .sphere.param.coords(
    h = h,
    w = w,
    radius = radius
  )
  coords_surface <- sphere.surface.embedding(
    h = h,
    w = w,
    surface = surface,
    radius = radius,
    amplitude = amplitude,
    freq_theta = freq_theta,
    freq_lat = freq_lat,
    twist = twist
  )
  weights <- .edge.weights.from.embedding(
    edges = edges,
    coords = coords_surface,
    normalize = normalize
  )

  out <- list(
    edges = edges,
    n = as.integer(2L + (h - 2L) * w),
    edge_weights = weights$edge_weights,
    coords_surface = coords_surface,
    coords_param = coords_param,
    weight_scale = weights$weight_scale,
    family = "sphere",
    surface = surface,
    normalize = normalize,
    label = sprintf("%s sphere %dx%d", tools::toTitleCase(surface), h, w)
  )
  class(out) <- c("grip_sphere_surface_graph", "list")
  out
}

#' Weighted recursive mask-grid surface helpers
#'
#' Convenience helpers that recursively subdivide a square grid according to a
#' square keep-mask, retain the occupied cells, and then lift those cells into
#' \eqn{\mathbb{R}^3}. The induced Euclidean edge lengths become positive graph
#' weights. This provides a generic mesh-derived family that includes the
#' Sierpinski carpet and related masked-grid fractals such as Vicsek-style
#' cross families.
#'
#' The \code{mask} is interpreted in standard matrix display orientation: rows
#' run from top to bottom and columns run from left to right. Non-zero entries
#' are retained at each recursive subdivision step.
#'
#' `recursive.mask.grid.surface.embedding()` returns the 3D coordinates of the
#' occupied cells in the same vertex order as \code{edges.recursive.mask.grid()}.
#' `recursive.mask.grid.surface.graph()` returns a reusable weighted-graph
#' bundle containing the masked-grid edges, induced edge weights, the 3D
#' surface coordinates, and the 2D parameter coordinates of the occupied cells.
#'
#' @param mask Square logical or numeric keep-mask with at least one retained
#'   cell. Non-zero entries are kept at each recursive step.
#' @param level Recursion depth. Must be at least \code{1}.
#' @param surface Surface family used for the lift. One of \code{"saddle"},
#'   \code{"paraboloid"}, or \code{"ripple"}.
#' @param amplitude Finite numeric amplitude controlling the non-flat
#'   displacement.
#' @param freq_u Positive ripple frequency in the horizontal parameter
#'   direction. Used only when \code{surface = "ripple"}.
#' @param freq_v Positive ripple frequency in the vertical parameter direction.
#'   Used only when \code{surface = "ripple"}.
#' @param x_scale Positive horizontal scaling of the parameter domain.
#' @param y_scale Positive vertical scaling of the parameter domain.
#' @param normalize Normalization applied to the induced edge lengths. One of
#'   \code{"median"}, \code{"mean"}, or \code{"none"}.
#'
#' @return
#' \code{recursive.mask.grid.surface.embedding()} returns an \code{n x 3}
#' numeric matrix with columns \code{x}, \code{y}, and \code{z}.
#'
#' \code{recursive.mask.grid.surface.graph()} returns a list with components:
#' \itemize{
#'   \item \code{edges}: the occupied-cell masked-grid edges,
#'   \item \code{n}: number of vertices,
#'   \item \code{edge_weights}: induced positive edge lengths,
#'   \item \code{coords_surface}: the 3D surface embedding,
#'   \item \code{coords_param}: the 2D parameter coordinates,
#'   \item \code{weight_scale}: the normalization constant applied to the raw
#'     edge lengths,
#'   \item \code{family}: always \code{"recursive.mask.grid"},
#'   \item \code{surface}: the chosen surface name,
#'   \item \code{level}: the recursion depth,
#'   \item \code{mask}: the logical keep-mask,
#'   \item \code{mask_size}: the side length of the mask,
#'   \item \code{side}: the side length of the fully refined grid,
#'   \item \code{label}: a human-readable family label.
#' }
#'
#' @name recursive_mask_grid_surface_helpers
NULL

#' @rdname recursive_mask_grid_surface_helpers
#' @export
recursive.mask.grid.surface.embedding <- function(mask,
                                                  level = 2,
                                                  surface = c("saddle", "paraboloid", "ripple"),
                                                  amplitude = 0.75,
                                                  freq_u = 1,
                                                  freq_v = 1,
                                                  x_scale = 1,
                                                  y_scale = 1) {
  mask <- .as_square_keep_mask(mask, "mask", min_size = 2L)
  level <- .as_whole_number(level, "level", min = 1L)
  surface <- match.arg(surface)
  amplitude <- .as_finite_scalar(amplitude, "amplitude")
  freq_u <- .as_positive_scalar(freq_u, "freq_u")
  freq_v <- .as_positive_scalar(freq_v, "freq_v")

  coords_param <- .recursive.mask.grid.param.coords(
    mask = mask,
    level = level,
    x_scale = x_scale,
    y_scale = y_scale
  )
  u <- coords_param[, 1L]
  v <- coords_param[, 2L]
  z <- .surface.z.from_uv(
    u = u,
    v = v,
    surface = surface,
    amplitude = amplitude,
    freq_u = freq_u,
    freq_v = freq_v
  )

  coords <- cbind(x = u, y = v, z = z)
  storage.mode(coords) <- "double"
  coords
}

#' @rdname recursive_mask_grid_surface_helpers
#' @export
recursive.mask.grid.surface.graph <- function(mask,
                                              level = 2,
                                              surface = c("saddle", "paraboloid", "ripple"),
                                              amplitude = 0.75,
                                              freq_u = 1,
                                              freq_v = 1,
                                              x_scale = 1,
                                              y_scale = 1,
                                              normalize = c("median", "mean", "none")) {
  mask <- .as_square_keep_mask(mask, "mask", min_size = 2L)
  level <- .as_whole_number(level, "level", min = 1L)
  surface <- match.arg(surface)
  normalize <- match.arg(normalize)

  spec <- .recursive.mask.grid.cells(mask, level)
  edges <- .recursive.mask.grid.edges(mask, level)
  coords_param <- .recursive.mask.grid.param.coords(
    mask = mask,
    level = level,
    x_scale = x_scale,
    y_scale = y_scale
  )
  coords_surface <- recursive.mask.grid.surface.embedding(
    mask = mask,
    level = level,
    surface = surface,
    amplitude = amplitude,
    freq_u = freq_u,
    freq_v = freq_v,
    x_scale = x_scale,
    y_scale = y_scale
  )
  weights <- .edge.weights.from.embedding(
    edges = edges,
    coords = coords_surface,
    normalize = normalize
  )

  out <- list(
    edges = edges,
    n = as.integer(nrow(coords_surface)),
    edge_weights = weights$edge_weights,
    coords_surface = coords_surface,
    coords_param = coords_param,
    weight_scale = weights$weight_scale,
    family = "recursive.mask.grid",
    surface = surface,
    level = level,
    mask = mask,
    mask_size = as.integer(nrow(mask)),
    side = spec$side,
    normalize = normalize,
    label = sprintf("Recursive %dx%d mask level %d",
                    nrow(mask), ncol(mask), level)
  )
  class(out) <- c("grip_recursive_mask_grid_surface_graph", "list")
  out
}

#' Weighted Sierpinski carpet surface helpers
#'
#' Convenience wrappers around \code{recursive.mask.grid.surface.*()} using the
#' classic \eqn{3 \times 3} carpet mask with the center cell removed.
#'
#' `sierpinski.carpet.surface.embedding()` returns the 3D coordinates of the
#' occupied cells in the same vertex order as \code{edges.sierpinski.carpet()}.
#' `sierpinski.carpet.surface.graph()` returns a reusable weighted-graph bundle
#' for the named Sierpinski carpet family.
#'
#' @inheritParams recursive.mask.grid.surface.embedding
#'
#' @return
#' \code{sierpinski.carpet.surface.embedding()} returns an \code{n x 3} numeric
#' matrix with columns \code{x}, \code{y}, and \code{z}, where
#' \code{n = 8^level}.
#'
#' \code{sierpinski.carpet.surface.graph()} returns the same components as
#' \code{recursive.mask.grid.surface.graph()}, with \code{family} set to
#' \code{"sierpinski.carpet"} and a family-specific class and label.
#'
#' @name sierpinski_carpet_surface_helpers
NULL

#' @rdname sierpinski_carpet_surface_helpers
#' @export
sierpinski.carpet.surface.embedding <- function(level = 2,
                                                surface = c("saddle", "paraboloid", "ripple"),
                                                amplitude = 0.75,
                                                freq_u = 1,
                                                freq_v = 1,
                                                x_scale = 1,
                                                y_scale = 1) {
  recursive.mask.grid.surface.embedding(
    mask = .sierpinski.carpet.mask(),
    level = level,
    surface = surface,
    amplitude = amplitude,
    freq_u = freq_u,
    freq_v = freq_v,
    x_scale = x_scale,
    y_scale = y_scale
  )
}

#' @rdname sierpinski_carpet_surface_helpers
#' @export
sierpinski.carpet.surface.graph <- function(level = 2,
                                            surface = c("saddle", "paraboloid", "ripple"),
                                            amplitude = 0.75,
                                            freq_u = 1,
                                            freq_v = 1,
                                            x_scale = 1,
                                            y_scale = 1,
                                            normalize = c("median", "mean", "none")) {
  out <- recursive.mask.grid.surface.graph(
    mask = .sierpinski.carpet.mask(),
    level = level,
    surface = surface,
    amplitude = amplitude,
    freq_u = freq_u,
    freq_v = freq_v,
    x_scale = x_scale,
    y_scale = y_scale,
    normalize = normalize
  )
  out$family <- "sierpinski.carpet"
  out$label <- sprintf("%s Sierpinski carpet level %d",
                       tools::toTitleCase(out$surface),
                       out$level)
  class(out) <- c("grip_sierpinski_carpet_surface_graph", class(out))
  out
}

#' Weighted Vicsek surface helpers
#'
#' Convenience wrappers around \code{recursive.mask.grid.surface.*()} using the
#' connected \eqn{3 \times 3} axial-cross mask: center plus the four cardinal
#' neighbors. This produces a mesh-derived fractal family with strong
#' bottlenecks while remaining orthogonally connected at every level.
#'
#' `vicsek.surface.embedding()` returns the 3D coordinates of the occupied
#' cells in the same vertex order as \code{edges.vicsek()}.
#' `vicsek.surface.graph()` returns a reusable weighted-graph bundle for the
#' named Vicsek family.
#'
#' @inheritParams recursive.mask.grid.surface.embedding
#'
#' @return
#' \code{vicsek.surface.embedding()} returns an \code{n x 3} numeric matrix
#' with columns \code{x}, \code{y}, and \code{z}, where \code{n = 5^level}.
#'
#' \code{vicsek.surface.graph()} returns the same components as
#' \code{recursive.mask.grid.surface.graph()}, with \code{family} set to
#' \code{"vicsek"} and a family-specific class and label.
#'
#' @name vicsek_surface_helpers
NULL

#' @rdname vicsek_surface_helpers
#' @export
vicsek.surface.embedding <- function(level = 2,
                                     surface = c("saddle", "paraboloid", "ripple"),
                                     amplitude = 0.75,
                                     freq_u = 1,
                                     freq_v = 1,
                                     x_scale = 1,
                                     y_scale = 1) {
  recursive.mask.grid.surface.embedding(
    mask = .vicsek.mask(),
    level = level,
    surface = surface,
    amplitude = amplitude,
    freq_u = freq_u,
    freq_v = freq_v,
    x_scale = x_scale,
    y_scale = y_scale
  )
}

#' @rdname vicsek_surface_helpers
#' @export
vicsek.surface.graph <- function(level = 2,
                                 surface = c("saddle", "paraboloid", "ripple"),
                                 amplitude = 0.75,
                                 freq_u = 1,
                                 freq_v = 1,
                                 x_scale = 1,
                                 y_scale = 1,
                                 normalize = c("median", "mean", "none")) {
  out <- recursive.mask.grid.surface.graph(
    mask = .vicsek.mask(),
    level = level,
    surface = surface,
    amplitude = amplitude,
    freq_u = freq_u,
    freq_v = freq_v,
    x_scale = x_scale,
    y_scale = y_scale,
    normalize = normalize
  )
  out$family <- "vicsek"
  out$label <- sprintf("%s Vicsek level %d",
                       tools::toTitleCase(out$surface),
                       out$level)
  class(out) <- c("grip_vicsek_surface_graph", class(out))
  out
}

#' @describeIn graph_generators Cylindrical grid graph with \code{h} rows and wrapped
#'   width \code{w}.
#' @export
edges.cylinder <- function(h, w = h) {
  h <- .as_whole_number(h, "h", min = 1L)
  w <- .as_whole_number(w, "w", min = 1L)
  idx <- function(i, j) (i - 1L) * w + j
  edges <- list()
  for (i in seq_len(h)) {
    for (j in seq_len(w)) {
      v <- idx(i, j)
      if (i < h) edges[[length(edges) + 1L]] <- c(v, idx(i + 1L, j))
      edges[[length(edges) + 1L]] <- c(v, idx(i, (j %% w) + 1L))
    }
  }
  .normalize_undirected_edges(.bind_edges(edges))
}

#' @describeIn graph_generators Toroidal grid graph with wrapped height and
#'   width.
#' @export
edges.torus <- function(h, w = h) {
  h <- .as_whole_number(h, "h", min = 1L)
  w <- .as_whole_number(w, "w", min = 1L)
  idx <- function(i, j) (i - 1L) * w + j
  edges <- list()
  for (i in seq_len(h)) {
    for (j in seq_len(w)) {
      v <- idx(i, j)
      edges[[length(edges) + 1L]] <- c(v, idx((i %% h) + 1L, j))
      edges[[length(edges) + 1L]] <- c(v, idx(i, (j %% w) + 1L))
    }
  }
  .normalize_undirected_edges(.bind_edges(edges))
}

#' @describeIn graph_generators Sphere surface graph with \code{h} latitude
#'   levels (including the poles) and wrapped longitude \code{w}.
#' @export
edges.sphere <- function(h, w = h) {
  h <- .as_whole_number(h, "h", min = 3L)
  w <- .as_whole_number(w, "w", min = 3L)
  ring.count <- h - 2L
  north <- 1L
  south <- 2L + ring.count * w
  idx <- function(i, j) 1L + (i - 1L) * w + j

  edges <- list()
  for (j in seq_len(w)) {
    edges[[length(edges) + 1L]] <- c(north, idx(1L, j))
  }
  for (i in seq_len(ring.count)) {
    for (j in seq_len(w)) {
      v <- idx(i, j)
      edges[[length(edges) + 1L]] <- c(v, idx(i, (j %% w) + 1L))
      if (i < ring.count) {
        edges[[length(edges) + 1L]] <- c(v, idx(i + 1L, j))
      } else {
        edges[[length(edges) + 1L]] <- c(v, south)
      }
    }
  }
  .normalize_undirected_edges(.bind_edges(edges))
}

#' @describeIn graph_generators Cube surface graph on the boundary of a
#'   \code{side x side x side} lattice.
#' @param side Number of lattice points along each cube edge.
#' @export
edges.cube <- function(side = 2) {
  side <- .as_whole_number(side, "side", min = 2L)
  grid <- expand.grid(
    x = seq_len(side),
    y = seq_len(side),
    z = seq_len(side)
  )
  keep <- grid$x %in% c(1L, side) |
    grid$y %in% c(1L, side) |
    grid$z %in% c(1L, side)
  points <- grid[keep, , drop = FALSE]
  ids <- seq_len(nrow(points))
  names(ids) <- paste(points$x, points$y, points$z, sep = ":")
  dirs <- rbind(c(1L, 0L, 0L), c(0L, 1L, 0L), c(0L, 0L, 1L))

  edges <- list()
  for (i in seq_len(nrow(points))) {
    p <- unlist(points[i, , drop = TRUE], use.names = FALSE)
    for (k in seq_len(nrow(dirs))) {
      q <- p + dirs[k, ]
      if (any(q < 1L | q > side)) {
        next
      }
      key <- paste(q[[1L]], q[[2L]], q[[3L]], sep = ":")
      if (!key %in% names(ids)) {
        next
      }
      edges[[length(edges) + 1L]] <- c(i, ids[[key]])
    }
  }
  .normalize_undirected_edges(.bind_edges(edges))
}

#' @describeIn graph_generators Full \code{k}-ary tree of depth \code{depth}.
#' @param k Branching factor.
#' @param depth Number of levels below the root.
#' @export
edges.kary.tree <- function(k = 2, depth = 2) {
  k <- .as_whole_number(k, "k", min = 1L)
  depth <- .as_whole_number(depth, "depth")
  if (depth == 0L) {
    return(.empty_edge_matrix())
  }

  edges <- list()
  current <- 1L
  next_id <- 2L
  for (d in seq_len(depth)) {
    new_parents <- integer()
    for (p in current) {
      for (i in seq_len(k)) {
        edges[[length(edges) + 1L]] <- c(p, next_id)
        new_parents <- c(new_parents, next_id)
        next_id <- next_id + 1L
      }
    }
    current <- new_parents
  }
  .normalize_undirected_edges(.bind_edges(edges))
}

#' @describeIn graph_generators Recursively refined square-mask grid graph whose
#'   vertices are occupied cells and whose edges connect orthogonally adjacent
#'   cells.
#' @param mask Square logical or numeric keep-mask. Non-zero entries are kept at
#'   each recursive subdivision step. Rows run from top to bottom and columns
#'   run from left to right in the motif.
#' @param level Recursion depth. For \code{edges.recursive.mask.grid()},
#'   \code{edges.vicsek()}, and \code{edges.sierpinski.carpet()},
#'   \code{level} must be at least 1.
#' @export
edges.recursive.mask.grid <- function(mask, level = 2) {
  .recursive.mask.grid.edges(mask, level)
}

#' @describeIn graph_generators Connected Vicsek-style cross family derived from
#'   a \code{3 x 3} axial-cross keep-mask.
#' @export
edges.vicsek <- function(level = 2) {
  edges.recursive.mask.grid(.vicsek.mask(), level = level)
}

#' @describeIn graph_generators Two-dimensional Sierpinski triangle graph at
#'   recursion depth \code{level}.
#' @examples
#' edges <- edges.sierpinski.triangle(2)
#' n <- max(edges)
#' coords <- grip.layout(edges, n = n, dim = 2,
#'                       placement = "circle",
#'                       seed = 1)
#' grip.plot(coords, edges, main = "Sierpinski triangle", pch = 16, cex = 0.7)
#' @export
edges.sierpinski.triangle <- function(level = 2) {
  level <- .as_whole_number(level, "level")

  merge_nodes <- function(edges, from, to) {
    edges[edges == from] <- to
    edges
  }

  build <- function(k) {
    if (k == 0L) {
      edges <- rbind(c(1L, 2L), c(2L, 3L), c(3L, 1L))
      return(list(edges = edges, corners = c(1L, 2L, 3L), n = 3L))
    }

    left <- build(k - 1L)
    right <- build(k - 1L)
    top <- build(k - 1L)

    off1 <- left$n
    off2 <- left$n + right$n

    edges <- rbind(left$edges,
                   right$edges + off1,
                   top$edges + off2)

    L <- left$corners
    R <- right$corners + off1
    T <- top$corners + off2

    edges <- merge_nodes(edges, R[1L], L[2L])
    edges <- merge_nodes(edges, T[1L], L[3L])
    edges <- merge_nodes(edges, T[2L], R[3L])

    ids <- sort(unique(c(edges)))
    map <- seq_along(ids)
    names(map) <- ids
    edges <- cbind(map[as.character(edges[, 1L])],
                   map[as.character(edges[, 2L])])
    edges <- .normalize_undirected_edges(edges)

    corners <- c(map[as.character(L[1L])],
                 map[as.character(R[2L])],
                 map[as.character(T[3L])])

    list(edges = edges, corners = corners, n = length(ids))
  }

  build(level)$edges
}

#' @describeIn graph_generators Three-dimensional tetrahedral Sierpinski graph
#'   at recursion depth \code{level}.
#' @export
edges.sierpinski.tetrahedron <- function(level = 2) {
  level <- .as_whole_number(level, "level")
  state <- new.env(parent = emptyenv())
  state$edges <- list()
  state$next_id <- 5L

  add_tetrahedron <- function(a, b, c, d) {
    state$edges[[length(state$edges) + 1L]] <<- rbind(
      c(a, b), c(a, c), c(a, d),
      c(b, c), c(b, d), c(c, d)
    )
  }

  recurse <- function(current_level, a, b, c, d) {
    if (current_level >= level) {
      add_tetrahedron(a, b, c, d)
    } else {
      e <- state$next_id
      f <- e + 1L
      g <- e + 2L
      h <- e + 3L
      i <- e + 4L
      j <- e + 5L
      state$next_id <- j + 1L

      recurse(current_level + 1L, a, g, e, f)
      recurse(current_level + 1L, e, b, i, h)
      recurse(current_level + 1L, f, c, j, h)
      recurse(current_level + 1L, g, d, i, j)
    }
  }

  recurse(0L, 1L, 2L, 3L, 4L)
  .normalize_undirected_edges(.bind_edges(state$edges))
}

#' @describeIn graph_generators Two-dimensional Sierpinski carpet graph whose
#'   vertices are occupied cells and whose edges connect orthogonally adjacent
#'   cells.
#' @export
edges.sierpinski.carpet <- function(level = 2) {
  edges.recursive.mask.grid(.sierpinski.carpet.mask(), level = level)
}
