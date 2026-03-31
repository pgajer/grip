# Helper utilities shared by the exported graph generators.
.empty_edge_matrix <- function() {
  matrix(integer(), ncol = 2L)
}

.empty_face_matrix <- function() {
  matrix(integer(), ncol = 3L)
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

.bind_faces <- function(faces) {
  if (length(faces) == 0L) {
    return(.empty_face_matrix())
  }
  out <- do.call(rbind, faces)
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

.as_cube_keep_mask <- function(mask, name = "mask", min_size = 2L) {
  if (!is.array(mask) || length(dim(mask)) != 3L) {
    stop(sprintf("%s must be a 3D array", name), call. = FALSE)
  }
  dims <- dim(mask)
  if (any(dims < min_size) || length(unique(dims)) != 1L) {
    stop(sprintf("%s must be a cubic 3D array with side length at least %d",
                 name, min_size),
         call. = FALSE)
  }
  if (!(is.logical(mask) || is.numeric(mask) || is.integer(mask))) {
    stop(sprintf("%s must be a logical or numeric 3D array", name),
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

.as_triangle_keep_mask <- function(mask, name = "mask") {
  slots <- c("left", "right", "top", "center")
  if (!(is.logical(mask) || is.numeric(mask) || is.integer(mask))) {
    stop(sprintf("%s must be a logical or numeric vector", name), call. = FALSE)
  }
  if (is.matrix(mask)) {
    if (length(mask) != 4L) {
      stop(sprintf("%s must have exactly 4 entries", name), call. = FALSE)
    }
    mask <- as.vector(mask)
  }
  if (length(mask) != 4L) {
    stop(sprintf("%s must have exactly 4 entries", name), call. = FALSE)
  }
  if (any(is.na(mask)) || (is.numeric(mask) && any(!is.finite(mask)))) {
    stop(sprintf("%s must not contain NA or non-finite values", name), call. = FALSE)
  }
  keep <- mask != 0
  storage.mode(keep) <- "logical"
  if (!is.null(names(mask))) {
    mask_names <- names(mask)
    if (length(unique(mask_names)) != 4L || !all(slots %in% mask_names)) {
      stop(sprintf("%s names must be exactly %s",
                   name,
                   paste(sprintf('"%s"', slots), collapse = ", ")),
           call. = FALSE)
    }
    keep <- keep[match(slots, mask_names)]
  }
  names(keep) <- slots
  if (!any(keep)) {
    stop(sprintf("%s must keep at least one subtriangle", name), call. = FALSE)
  }
  keep
}

.as_tetrahedron_keep_mask <- function(mask, name = "mask") {
  slots <- c("base_left", "base_right", "base_back", "apex")
  if (!(is.logical(mask) || is.numeric(mask) || is.integer(mask))) {
    stop(sprintf("%s must be a logical or numeric vector", name), call. = FALSE)
  }
  if (is.matrix(mask)) {
    if (length(mask) != 4L) {
      stop(sprintf("%s must have exactly 4 entries", name), call. = FALSE)
    }
    mask <- as.vector(mask)
  }
  if (length(mask) != 4L) {
    stop(sprintf("%s must have exactly 4 entries", name), call. = FALSE)
  }
  if (any(is.na(mask)) || (is.numeric(mask) && any(!is.finite(mask)))) {
    stop(sprintf("%s must not contain NA or non-finite values", name), call. = FALSE)
  }
  keep <- mask != 0
  storage.mode(keep) <- "logical"
  if (!is.null(names(mask))) {
    mask_names <- names(mask)
    if (length(unique(mask_names)) != 4L || !all(slots %in% mask_names)) {
      stop(sprintf("%s names must be exactly %s",
                   name,
                   paste(sprintf('"%s"', slots), collapse = ", ")),
           call. = FALSE)
    }
    keep <- keep[match(slots, mask_names)]
  }
  names(keep) <- slots
  if (!any(keep)) {
    stop(sprintf("%s must keep at least one subtetrahedron", name), call. = FALSE)
  }
  keep
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

.normalize.center_radius.coords <- function(coords) {
  coords <- as.matrix(coords)
  centered <- scale(coords, center = TRUE, scale = FALSE)
  radius <- max(sqrt(rowSums(centered^2)))
  if (!is.finite(radius) || radius <= 0) {
    return(centered)
  }
  centered / radius
}

.normalize.row.coords <- function(coords) {
  coords <- as.matrix(coords)
  norms <- sqrt(rowSums(coords^2))
  if (any(!is.finite(norms) | norms <= 0)) {
    stop("coords must have positive finite row norms", call. = FALSE)
  }
  sweep(coords, 1L, norms, "/")
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

.recursive.cube.mask.cells <- function(mask, level) {
  mask <- .as_cube_keep_mask(mask, "mask", min_size = 2L)
  level <- .as_whole_number(level, "level", min = 1L)
  base <- as.integer(dim(mask)[1L])
  side <- as.integer(base^level)
  grid <- expand.grid(col = 0:(side - 1L),
                      row = 0:(side - 1L),
                      layer = 0:(side - 1L))

  keep_cell <- function(col, row, layer) {
    for (step in seq_len(level)) {
      col_digit <- (col %% base) + 1L
      row_digit <- (row %% base) + 1L
      layer_digit <- (layer %% base) + 1L
      if (!mask[row_digit, col_digit, layer_digit]) {
        return(FALSE)
      }
      col <- col %/% base
      row <- row %/% base
      layer <- layer %/% base
    }
    TRUE
  }

  keep <- mapply(keep_cell, grid$col, grid$row, grid$layer)
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

.recursive.cube.mask.param.coords <- function(mask,
                                              level,
                                              x_scale = 1,
                                              y_scale = 1,
                                              z_scale = 1) {
  spec <- .recursive.cube.mask.cells(mask, level)
  x_scale <- .as_positive_scalar(x_scale, "x_scale")
  y_scale <- .as_positive_scalar(y_scale, "y_scale")
  z_scale <- .as_positive_scalar(z_scale, "z_scale")

  x_vals <- .grid_axis(spec$side) * x_scale
  y_vals <- rev(.grid_axis(spec$side) * y_scale)
  z_vals <- .grid_axis(spec$side) * z_scale
  coords <- cbind(
    x = x_vals[spec$cells$col + 1L],
    y = y_vals[spec$cells$row + 1L],
    z = z_vals[spec$cells$layer + 1L]
  )
  storage.mode(coords) <- "double"
  coords
}

.recursive.cube.mask.edges <- function(mask, level) {
  spec <- .recursive.cube.mask.cells(mask, level)
  side <- spec$side
  cells <- spec$cells

  id_map <- array(0L, dim = c(side, side, side))
  for (i in seq_len(nrow(cells))) {
    id_map[cells$row[i] + 1L, cells$col[i] + 1L, cells$layer[i] + 1L] <- i
  }

  edges <- list()
  for (i in seq_len(nrow(cells))) {
    col <- cells$col[i]
    row <- cells$row[i]
    layer <- cells$layer[i]
    if (col + 1L < side) {
      nbr <- id_map[row + 1L, col + 2L, layer + 1L]
      if (nbr > 0L) edges[[length(edges) + 1L]] <- c(i, nbr)
    }
    if (row + 1L < side) {
      nbr <- id_map[row + 2L, col + 1L, layer + 1L]
      if (nbr > 0L) edges[[length(edges) + 1L]] <- c(i, nbr)
    }
    if (layer + 1L < side) {
      nbr <- id_map[row + 1L, col + 1L, layer + 2L]
      if (nbr > 0L) edges[[length(edges) + 1L]] <- c(i, nbr)
    }
  }

  .bind_edges(edges)
}

.triangulated.polyhedron.base <- function(base = c("tetrahedron", "octahedron", "icosahedron")) {
  base <- match.arg(base)
  out <- switch(
    base,
    tetrahedron = list(
      vertices = rbind(
        c(1, 1, 1),
        c(-1, -1, 1),
        c(-1, 1, -1),
        c(1, -1, -1)
      ),
      faces = matrix(
        c(
          1, 2, 3,
          1, 4, 2,
          1, 3, 4,
          2, 4, 3
        ),
        ncol = 3L,
        byrow = TRUE
      )
    ),
    octahedron = list(
      vertices = rbind(
        c(1, 0, 0),
        c(-1, 0, 0),
        c(0, 1, 0),
        c(0, -1, 0),
        c(0, 0, 1),
        c(0, 0, -1)
      ),
      faces = matrix(
        c(
          1, 3, 5,
          3, 2, 5,
          2, 4, 5,
          4, 1, 5,
          3, 1, 6,
          2, 3, 6,
          4, 2, 6,
          1, 4, 6
        ),
        ncol = 3L,
        byrow = TRUE
      )
    ),
    icosahedron = {
      phi <- (1 + sqrt(5)) / 2
      list(
        vertices = rbind(
          c(-1, phi, 0),
          c(1, phi, 0),
          c(-1, -phi, 0),
          c(1, -phi, 0),
          c(0, -1, phi),
          c(0, 1, phi),
          c(0, -1, -phi),
          c(0, 1, -phi),
          c(phi, 0, -1),
          c(phi, 0, 1),
          c(-phi, 0, -1),
          c(-phi, 0, 1)
        ),
        faces = matrix(
          c(
            1, 12, 6,
            1, 6, 2,
            1, 2, 8,
            1, 8, 11,
            1, 11, 12,
            2, 6, 10,
            6, 12, 5,
            12, 11, 3,
            11, 8, 7,
            8, 2, 9,
            4, 10, 5,
            4, 5, 3,
            4, 3, 7,
            4, 7, 9,
            4, 9, 10,
            5, 10, 6,
            3, 5, 12,
            7, 3, 11,
            9, 7, 8,
            10, 9, 2
          ),
          ncol = 3L,
          byrow = TRUE
        )
      )
    }
  )
  out$vertices <- .normalize.row.coords(out$vertices)
  storage.mode(out$vertices) <- "double"
  storage.mode(out$faces) <- "integer"
  out
}

.triangulated.polyhedron.canonical <- function(base = c("tetrahedron", "octahedron", "icosahedron"),
                                               level = 1) {
  base <- match.arg(base)
  level <- .as_whole_number(level, "level")
  spec <- .triangulated.polyhedron.base(base)
  vertices <- spec$vertices
  faces <- spec$faces
  subdiv <- as.integer(2^level)

  face_build <- function(a, b, c) {
    id_map <- matrix(NA_integer_, nrow = subdiv + 1L, ncol = subdiv + 1L)
    coords <- vector("list", ((subdiv + 1L) * (subdiv + 2L)) %/% 2L)
    next_id <- 1L
    for (i in 0:subdiv) {
      for (j in 0:(subdiv - i)) {
        id_map[i + 1L, j + 1L] <- next_id
        point <- a + (i / subdiv) * (b - a) + (j / subdiv) * (c - a)
        coords[[next_id]] <- point
        next_id <- next_id + 1L
      }
    }

    edges <- list()
    faces <- list()
    for (i in 0:subdiv) {
      for (j in 0:(subdiv - i)) {
        cur <- id_map[i + 1L, j + 1L]
        if ((i + 1L + j) <= subdiv) {
          edges[[length(edges) + 1L]] <- c(cur, id_map[i + 2L, j + 1L])
        }
        if ((i + j + 1L) <= subdiv) {
          edges[[length(edges) + 1L]] <- c(cur, id_map[i + 1L, j + 2L])
        }
        if (j >= 1L && (i + 1L + j - 1L) <= subdiv) {
          edges[[length(edges) + 1L]] <- c(cur, id_map[i + 2L, j])
        }
      }
    }
    if (subdiv >= 1L) {
      for (i in 0:(subdiv - 1L)) {
        for (j in 0:(subdiv - 1L - i)) {
          v00 <- id_map[i + 1L, j + 1L]
          v10 <- id_map[i + 2L, j + 1L]
          v01 <- id_map[i + 1L, j + 2L]
          faces[[length(faces) + 1L]] <- c(v00, v10, v01)
          if (j <= (subdiv - i - 2L)) {
            v11 <- id_map[i + 2L, j + 2L]
            faces[[length(faces) + 1L]] <- c(v10, v11, v01)
          }
        }
      }
    }

    list(
      coords = do.call(rbind, coords),
      edges = .bind_edges(edges),
      faces = .bind_faces(faces)
    )
  }

  coords_list <- list()
  edges_list <- list()
  faces_list <- list()
  offset <- 0L
  for (face_idx in seq_len(nrow(faces))) {
    face <- faces[face_idx, ]
    built <- face_build(vertices[face[1L], ],
                        vertices[face[2L], ],
                        vertices[face[3L], ])
    coords_list[[face_idx]] <- built$coords
    edges_list[[face_idx]] <- built$edges + offset
    faces_list[[face_idx]] <- built$faces + offset
    offset <- offset + nrow(built$coords)
  }

  merged <- .deduplicate.coordinate.graph(
    edges = .bind_edges(edges_list),
    coords = do.call(rbind, coords_list),
    return_map = TRUE
  )
  coords <- merged$coords
  colnames(coords) <- c("x", "y", "z")
  storage.mode(coords) <- "double"
  faces_all <- .bind_faces(faces_list)
  faces_new <- if (nrow(faces_all) == 0L) {
    .empty_face_matrix()
  } else {
    mapped <- cbind(
      merged$old_to_new[faces_all[, 1L]],
      merged$old_to_new[faces_all[, 2L]],
      merged$old_to_new[faces_all[, 3L]]
    )
    keep_face <- apply(mapped, 1L, function(row) length(unique(row)) == 3L)
    mapped <- mapped[keep_face, , drop = FALSE]
    if (nrow(mapped) == 0L) {
      .empty_face_matrix()
    } else {
      unique(t(apply(mapped, 1L, sort)))
    }
  }
  storage.mode(faces_new) <- "integer"
  list(
    edges = merged$edges,
    faces = faces_new,
    coords = coords,
    n = as.integer(nrow(coords)),
    base = base,
    subdivision = subdiv
  )
}

.largest.connected.vertex.set <- function(edges, n) {
  n <- .as_whole_number(n, "n", min = 1L)
  if (n == 1L || nrow(edges) == 0L) {
    return(seq_len(n))
  }

  neighbors <- vector("list", n)
  for (i in seq_len(nrow(edges))) {
    a <- edges[i, 1L]
    b <- edges[i, 2L]
    neighbors[[a]] <- c(neighbors[[a]], b)
    neighbors[[b]] <- c(neighbors[[b]], a)
  }

  seen <- rep(FALSE, n)
  best <- integer(0L)
  for (start in seq_len(n)) {
    if (seen[[start]]) {
      next
    }
    queue <- start
    seen[[start]] <- TRUE
    component <- integer(0L)
    head <- 1L
    while (head <= length(queue)) {
      v <- queue[[head]]
      head <- head + 1L
      component <- c(component, v)
      for (nbr in neighbors[[v]]) {
        if (!seen[[nbr]]) {
          seen[[nbr]] <- TRUE
          queue <- c(queue, nbr)
        }
      }
    }
    if (length(component) > length(best)) {
      best <- component
    }
  }
  sort(best)
}

.induce.subgraph.coords <- function(edges, coords, keep_idx) {
  keep_idx <- sort(unique(as.integer(keep_idx)))
  coords <- as.matrix(coords)
  if (length(keep_idx) == 0L) {
    stop("keep_idx must contain at least one vertex", call. = FALSE)
  }

  map <- integer(nrow(coords))
  map[keep_idx] <- seq_along(keep_idx)
  keep_edge <- edges[, 1L] %in% keep_idx & edges[, 2L] %in% keep_idx
  edges_new <- if (!any(keep_edge)) {
    .empty_edge_matrix()
  } else {
    cbind(
      map[edges[keep_edge, 1L]],
      map[edges[keep_edge, 2L]]
    )
  }
  storage.mode(edges_new) <- "integer"
  coords_new <- coords[keep_idx, , drop = FALSE]
  rownames(coords_new) <- NULL
  list(
    edges = .normalize_undirected_edges(edges_new),
    coords = coords_new
  )
}

.triangular.lattice.region.canonical <- function(bbox, spacing, keep_fn) {
  if (!is.numeric(bbox) || length(bbox) != 4L || any(!is.finite(bbox))) {
    stop("bbox must be a finite numeric vector c(xmin, xmax, ymin, ymax)",
         call. = FALSE)
  }
  spacing <- .as_positive_scalar(spacing, "spacing")
  if (!is.function(keep_fn)) {
    stop("keep_fn must be a function", call. = FALSE)
  }
  xmin <- bbox[[1L]]
  xmax <- bbox[[2L]]
  ymin <- bbox[[3L]]
  ymax <- bbox[[4L]]
  if (!(xmin < xmax && ymin < ymax)) {
    stop("bbox must satisfy xmin < xmax and ymin < ymax", call. = FALSE)
  }

  dy <- spacing * sqrt(3) / 2
  row_vals <- seq(ymin - dy, ymax + dy, by = dy)
  col_count <- .as_whole_number(ceiling((xmax - xmin) / spacing) + 4L,
                                "col_count",
                                min = 3L)
  n_rows <- length(row_vals)
  x_mat <- matrix(0, nrow = n_rows, ncol = col_count)
  y_mat <- matrix(0, nrow = n_rows, ncol = col_count)
  keep <- matrix(FALSE, nrow = n_rows, ncol = col_count)
  for (row in seq_len(n_rows)) {
    x_offset <- 0.5 * ((row - 1L) %% 2L)
    x_vals <- xmin - spacing + spacing * ((0:(col_count - 1L)) + x_offset)
    y_val <- row_vals[[row]]
    x_mat[row, ] <- x_vals
    y_mat[row, ] <- y_val
    keep[row, ] <- keep_fn(x_vals, rep(y_val, col_count))
  }

  positions <- which(keep, arr.ind = TRUE)
  if (nrow(positions) == 0L) {
    stop("shape and resolution keep no triangulation vertices", call. = FALSE)
  }
  positions <- positions[order(positions[, 1L], positions[, 2L]), , drop = FALSE]
  id_map <- matrix(0L, nrow = n_rows, ncol = col_count)
  id_map[positions] <- seq_len(nrow(positions))
  coords <- cbind(
    u = x_mat[positions],
    v = y_mat[positions]
  )
  storage.mode(coords) <- "double"

  edges <- list()
  for (i in seq_len(nrow(positions))) {
    row <- positions[i, 1L]
    col <- positions[i, 2L]
    if (col < col_count) {
      nbr <- id_map[row, col + 1L]
      if (nbr > 0L) {
        edges[[length(edges) + 1L]] <- c(i, nbr)
      }
    }
    if (row < n_rows) {
      down_same <- id_map[row + 1L, col]
      if (down_same > 0L) {
        edges[[length(edges) + 1L]] <- c(i, down_same)
      }
      if (((row - 1L) %% 2L) == 0L) {
        if (col > 1L) {
          down_diag <- id_map[row + 1L, col - 1L]
          if (down_diag > 0L) {
            edges[[length(edges) + 1L]] <- c(i, down_diag)
          }
        }
      } else if (col < col_count) {
        down_diag <- id_map[row + 1L, col + 1L]
        if (down_diag > 0L) {
          edges[[length(edges) + 1L]] <- c(i, down_diag)
        }
      }
    }
  }

  edges <- .normalize_undirected_edges(.bind_edges(edges))
  if (nrow(edges) == 0L) {
    stop("shape and resolution produce no triangulation edges", call. = FALSE)
  }
  keep_idx <- .largest.connected.vertex.set(edges, nrow(coords))
  induced <- .induce.subgraph.coords(edges, coords, keep_idx)
  colnames(induced$coords) <- c("u", "v")
  list(
    edges = induced$edges,
    coords = induced$coords,
    n = as.integer(nrow(induced$coords))
  )
}

.triangulated.annulus.canonical <- function(resolution = 12,
                                            outer_radius = 1,
                                            inner_radius = 0.45) {
  resolution <- .as_whole_number(resolution, "resolution", min = 4L)
  outer_radius <- .as_positive_scalar(outer_radius, "outer_radius")
  inner_radius <- .as_positive_scalar(inner_radius, "inner_radius")
  if (inner_radius >= outer_radius) {
    stop("inner_radius must be < outer_radius", call. = FALSE)
  }
  spacing <- outer_radius / resolution
  if ((outer_radius - inner_radius) < (2 * spacing)) {
    stop("annulus thickness is too small for the chosen resolution", call. = FALSE)
  }

  keep_fn <- function(x, y) {
    r2 <- x^2 + y^2
    r2 <= outer_radius^2 & r2 >= inner_radius^2
  }
  built <- .triangular.lattice.region.canonical(
    bbox = c(-outer_radius, outer_radius, -outer_radius, outer_radius),
    spacing = spacing,
    keep_fn = keep_fn
  )
  built$resolution <- resolution
  built$outer_radius <- outer_radius
  built$inner_radius <- inner_radius
  built
}

.triangulated.pair.of.pants.canonical <- function(resolution = 12,
                                                  outer_radius = 1.1,
                                                  hole_radius = 0.24,
                                                  hole_offset = 0.38,
                                                  hole_height = 0.18) {
  resolution <- .as_whole_number(resolution, "resolution", min = 4L)
  outer_radius <- .as_positive_scalar(outer_radius, "outer_radius")
  hole_radius <- .as_positive_scalar(hole_radius, "hole_radius")
  hole_offset <- .as_positive_scalar(hole_offset, "hole_offset")
  hole_height <- .as_finite_scalar(hole_height, "hole_height")
  if (hole_offset <= hole_radius) {
    stop("hole_offset must exceed hole_radius so the two holes do not overlap",
         call. = FALSE)
  }
  hole_center_radius <- sqrt(hole_offset^2 + hole_height^2)
  if ((hole_center_radius + hole_radius) >= outer_radius) {
    stop("holes must lie strictly inside the outer boundary", call. = FALSE)
  }
  spacing <- outer_radius / resolution
  if ((outer_radius - hole_center_radius - hole_radius) < spacing) {
    stop("pair-of-pants neck is too thin for the chosen resolution", call. = FALSE)
  }

  keep_fn <- function(x, y) {
    outer <- (x^2 + y^2) <= outer_radius^2
    left_hole <- ((x + hole_offset)^2 + (y - hole_height)^2) < hole_radius^2
    right_hole <- ((x - hole_offset)^2 + (y - hole_height)^2) < hole_radius^2
    outer & !left_hole & !right_hole
  }
  built <- .triangular.lattice.region.canonical(
    bbox = c(-outer_radius, outer_radius, -outer_radius, outer_radius),
    spacing = spacing,
    keep_fn = keep_fn
  )
  built$resolution <- resolution
  built$outer_radius <- outer_radius
  built$hole_radius <- hole_radius
  built$hole_offset <- hole_offset
  built$hole_height <- hole_height
  built
}

.normalize.angles <- function(theta) {
  theta <- theta %% (2 * pi)
  theta[theta < 0] <- theta[theta < 0] + 2 * pi
  theta
}

.cyclic.ring.edges <- function(ids) {
  ids <- as.integer(ids)
  if (length(ids) < 3L) {
    stop("cyclic rings require at least 3 vertices", call. = FALSE)
  }
  edges <- cbind(ids, c(ids[-1L], ids[1L]))
  storage.mode(edges) <- "integer"
  .normalize_undirected_edges(edges)
}

.cyclic.bracketing.indices <- function(theta, ring_angles) {
  ring_angles <- .normalize.angles(ring_angles)
  order_idx <- order(ring_angles)
  sorted_angles <- ring_angles[order_idx]
  theta <- .normalize.angles(theta)
  if (length(sorted_angles) == 1L) {
    return(rep(order_idx[[1L]], 2L))
  }
  aug <- c(sorted_angles, sorted_angles[1L] + 2 * pi)
  theta_aug <- theta
  if (theta_aug < sorted_angles[1L]) {
    theta_aug <- theta_aug + 2 * pi
  }
  lower_sorted <- findInterval(theta_aug, aug, rightmost.closed = TRUE)
  lower_sorted <- min(max(lower_sorted, 1L), length(sorted_angles))
  upper_sorted <- (lower_sorted %% length(sorted_angles)) + 1L
  c(order_idx[lower_sorted], order_idx[upper_sorted])
}

.connect.cyclic.rings <- function(ids_a, angles_a, ids_b, angles_b) {
  ids_a <- as.integer(ids_a)
  ids_b <- as.integer(ids_b)
  angles_a <- .normalize.angles(angles_a)
  angles_b <- .normalize.angles(angles_b)
  edges <- list()

  for (j in seq_along(ids_b)) {
    nbr <- .cyclic.bracketing.indices(angles_b[[j]], angles_a)
    edges[[length(edges) + 1L]] <- c(ids_b[[j]], ids_a[nbr[[1L]]])
    edges[[length(edges) + 1L]] <- c(ids_b[[j]], ids_a[nbr[[2L]]])
  }
  for (i in seq_along(ids_a)) {
    nbr <- .cyclic.bracketing.indices(angles_a[[i]], angles_b)
    edges[[length(edges) + 1L]] <- c(ids_a[[i]], ids_b[nbr[[1L]]])
    edges[[length(edges) + 1L]] <- c(ids_a[[i]], ids_b[nbr[[2L]]])
  }

  .normalize_undirected_edges(.bind_edges(edges))
}

.connect.cyclic.component_sets <- function(comps_a,
                                           comps_b,
                                           overlap_scale = 1.05) {
  overlap_scale <- .as_positive_scalar(overlap_scale, "overlap_scale")
  if (length(comps_a) == 0L || length(comps_b) == 0L) {
    return(matrix(integer(), ncol = 2L))
  }

  centers_a <- do.call(rbind, lapply(comps_a, `[[`, "center"))
  centers_b <- do.call(rbind, lapply(comps_b, `[[`, "center"))
  radii_a <- vapply(comps_a, `[[`, numeric(1L), "radius")
  radii_b <- vapply(comps_b, `[[`, numeric(1L), "radius")
  dmat <- outer(
    seq_len(nrow(centers_a)),
    seq_len(nrow(centers_b)),
    Vectorize(function(i, j) sqrt(sum((centers_a[i, ] - centers_b[j, ])^2)))
  )
  thresh <- outer(radii_a, radii_b, "+") * overlap_scale
  overlap_pairs <- which(dmat <= thresh, arr.ind = TRUE)

  pair_keys <- character()
  pair_idx <- matrix(integer(), ncol = 2L)
  add_pair <- function(i, j) {
    key <- sprintf("%d-%d", i, j)
    if (!(key %in% pair_keys)) {
      pair_keys <<- c(pair_keys, key)
      pair_idx <<- rbind(pair_idx, c(i, j))
    }
  }

  if (nrow(overlap_pairs) > 0L) {
    for (k in seq_len(nrow(overlap_pairs))) {
      add_pair(overlap_pairs[k, 1L], overlap_pairs[k, 2L])
    }
  }
  for (i in seq_along(comps_a)) {
    if (!any(pair_idx[, 1L] == i)) {
      add_pair(i, which.min(dmat[i, ]))
    }
  }
  for (j in seq_along(comps_b)) {
    if (!any(pair_idx[, 2L] == j)) {
      add_pair(which.min(dmat[, j]), j)
    }
  }

  edges <- vector("list", nrow(pair_idx))
  for (k in seq_len(nrow(pair_idx))) {
    i <- pair_idx[k, 1L]
    j <- pair_idx[k, 2L]
    edges[[k]] <- .connect.cyclic.rings(
      ids_a = comps_a[[i]]$ids,
      angles_a = comps_a[[i]]$angles,
      ids_b = comps_b[[j]]$ids,
      angles_b = comps_b[[j]]$angles
    )
  }

  .normalize_undirected_edges(.bind_edges(edges))
}

.irregular.annulus.canonical <- function(rings = 6,
                                         outer_count = 28,
                                         outer_radius = 1,
                                         inner_radius = 0.45,
                                         count_irregularity = 0.2,
                                         radial_irregularity = 0.35,
                                         phase_twist = 0.35) {
  rings <- .as_whole_number(rings, "rings", min = 2L)
  outer_count <- .as_whole_number(outer_count, "outer_count", min = 6L)
  outer_radius <- .as_positive_scalar(outer_radius, "outer_radius")
  inner_radius <- .as_positive_scalar(inner_radius, "inner_radius")
  if (inner_radius >= outer_radius) {
    stop("inner_radius must be < outer_radius", call. = FALSE)
  }
  count_irregularity <- .as_finite_scalar(count_irregularity, "count_irregularity")
  radial_irregularity <- .as_finite_scalar(radial_irregularity, "radial_irregularity")
  phase_twist <- .as_finite_scalar(phase_twist, "phase_twist")
  if (count_irregularity < 0 || count_irregularity >= 1) {
    stop("count_irregularity must be in [0, 1)", call. = FALSE)
  }
  if (radial_irregularity < 0 || radial_irregularity > 1) {
    stop("radial_irregularity must be in [0, 1]", call. = FALSE)
  }

  radii <- seq(inner_radius, outer_radius, length.out = rings)
  ring_gap <- if (rings > 1L) min(diff(radii)) else (outer_radius - inner_radius)
  radial_amp <- radial_irregularity * 0.45 * ring_gap
  coords_list <- list()
  ring_ids <- vector("list", rings)
  ring_angles <- vector("list", rings)
  ring_sizes <- integer(rings)
  next_id <- 1L

  for (ring in seq_len(rings)) {
    t <- if (rings == 1L) 0 else (ring - 1L) / (rings - 1L)
    target_count <- outer_count *
      (radii[[ring]] / outer_radius) *
      (1 + count_irregularity * cos(2 * pi * t + 0.4))
    count <- max(5L, as.integer(round(target_count)))
    theta <- seq(0, 2 * pi * (1 - 1 / count), length.out = count)
    theta <- .normalize.angles(theta + phase_twist * sin(pi * (t - 0.5)))
    theta <- sort(theta)
    radius_local <- radii[[ring]] + radial_amp * sin(3 * theta + 2 * pi * t)
    if (any(!is.finite(radius_local) | radius_local <= 0)) {
      stop("irregular annulus parameters produce non-positive local radii",
           call. = FALSE)
    }

    ids <- next_id:(next_id + count - 1L)
    next_id <- next_id + count
    ring_ids[[ring]] <- as.integer(ids)
    ring_angles[[ring]] <- theta
    ring_sizes[[ring]] <- count
    coords_list[[ring]] <- cbind(
      u = radius_local * cos(theta),
      v = radius_local * sin(theta)
    )
  }

  coords <- do.call(rbind, coords_list)
  storage.mode(coords) <- "double"
  edges <- list()
  for (ring in seq_len(rings)) {
    edges[[length(edges) + 1L]] <- .cyclic.ring.edges(ring_ids[[ring]])
    if (ring < rings) {
      edges[[length(edges) + 1L]] <- .connect.cyclic.rings(
        ids_a = ring_ids[[ring]],
        angles_a = ring_angles[[ring]],
        ids_b = ring_ids[[ring + 1L]],
        angles_b = ring_angles[[ring + 1L]]
      )
    }
  }

  list(
    edges = .normalize_undirected_edges(.bind_edges(edges)),
    coords = coords,
    n = as.integer(nrow(coords)),
    rings = rings,
    ring_sizes = ring_sizes,
    outer_radius = outer_radius,
    inner_radius = inner_radius
  )
}

.irregular.sphere.canonical <- function(bands = 6,
                                        equator_count = 28,
                                        count_irregularity = 0.2,
                                        lat_irregularity = 0.35,
                                        phase_twist = 0.35) {
  bands <- .as_whole_number(bands, "bands", min = 2L)
  equator_count <- .as_whole_number(equator_count, "equator_count", min = 6L)
  count_irregularity <- .as_finite_scalar(count_irregularity, "count_irregularity")
  lat_irregularity <- .as_finite_scalar(lat_irregularity, "lat_irregularity")
  phase_twist <- .as_finite_scalar(phase_twist, "phase_twist")
  if (count_irregularity < 0 || count_irregularity >= 1) {
    stop("count_irregularity must be in [0, 1)", call. = FALSE)
  }
  if (lat_irregularity < 0 || lat_irregularity > 1) {
    stop("lat_irregularity must be in [0, 1]", call. = FALSE)
  }

  base_gap <- pi / (bands + 1L)
  lat_base <- seq(pi / 2 - base_gap, -pi / 2 + base_gap, length.out = bands)
  lat_shift <- lat_irregularity * 0.22 * base_gap *
    sin(seq_len(bands) * (2 * pi / (bands + 1L)))
  lat_vals <- sort(lat_base + lat_shift, decreasing = TRUE)

  north_id <- 1L
  coords_param_list <- list(matrix(c(0, pi / 2), nrow = 1L,
                                   dimnames = list(NULL, c("theta", "latitude"))))
  ring_ids <- vector("list", bands)
  ring_angles <- vector("list", bands)
  band_sizes <- integer(bands)
  next_id <- 2L
  edges <- list()

  for (band in seq_len(bands)) {
    lat <- lat_vals[[band]]
    t <- (band - 1L) / max(1L, bands - 1L)
    target_count <- equator_count *
      max(cos(lat), 0.2) *
      (1 + count_irregularity * sin(3 * pi * t + 0.3))
    count <- max(3L, as.integer(round(target_count)))
    theta <- seq(0, 2 * pi * (1 - 1 / count), length.out = count)
    theta <- .normalize.angles(theta + phase_twist * sin(2 * lat))
    theta <- sort(theta)
    ids <- next_id:(next_id + count - 1L)
    next_id <- next_id + count
    ring_ids[[band]] <- as.integer(ids)
    ring_angles[[band]] <- theta
    band_sizes[[band]] <- count
    coords_param_list[[band + 1L]] <- cbind(theta = theta, latitude = rep(lat, count))
    edges[[length(edges) + 1L]] <- .cyclic.ring.edges(ids)
  }

  south_id <- next_id
  coords_param_list[[bands + 2L]] <- matrix(c(0, -pi / 2), nrow = 1L,
                                            dimnames = list(NULL, c("theta", "latitude")))
  if (bands >= 1L) {
    first_ring <- ring_ids[[1L]]
    last_ring <- ring_ids[[bands]]
    for (id in first_ring) {
      edges[[length(edges) + 1L]] <- c(north_id, id)
    }
    for (id in last_ring) {
      edges[[length(edges) + 1L]] <- c(south_id, id)
    }
  }
  if (bands >= 2L) {
    for (band in seq_len(bands - 1L)) {
      edges[[length(edges) + 1L]] <- .connect.cyclic.rings(
        ids_a = ring_ids[[band]],
        angles_a = ring_angles[[band]],
        ids_b = ring_ids[[band + 1L]],
        angles_b = ring_angles[[band + 1L]]
      )
    }
  }

  coords_param <- do.call(rbind, coords_param_list)
  colnames(coords_param) <- c("theta", "latitude")
  storage.mode(coords_param) <- "double"
  list(
    edges = .normalize_undirected_edges(.bind_edges(edges)),
    coords_param = coords_param,
    n = as.integer(nrow(coords_param)),
    bands = bands,
    band_sizes = band_sizes
  )
}

.irregular.torus.canonical <- function(major_rings = 8,
                                       tube_count = 16,
                                       count_irregularity = 0.2,
                                       major_irregularity = 0.25,
                                       phase_twist = 0.35) {
  major_rings <- .as_whole_number(major_rings, "major_rings", min = 4L)
  tube_count <- .as_whole_number(tube_count, "tube_count", min = 6L)
  count_irregularity <- .as_finite_scalar(count_irregularity, "count_irregularity")
  major_irregularity <- .as_finite_scalar(major_irregularity, "major_irregularity")
  phase_twist <- .as_finite_scalar(phase_twist, "phase_twist")
  if (count_irregularity < 0 || count_irregularity >= 1) {
    stop("count_irregularity must be in [0, 1)", call. = FALSE)
  }
  if (major_irregularity < 0 || major_irregularity > 1) {
    stop("major_irregularity must be in [0, 1]", call. = FALSE)
  }

  theta_gap <- 2 * pi / major_rings
  theta_base <- seq(0, 2 * pi * (1 - 1 / major_rings), length.out = major_rings)
  theta_shift <- major_irregularity * 0.22 * theta_gap *
    sin(seq_len(major_rings) * (2 * pi / major_rings))
  theta_vals <- sort(.normalize.angles(theta_base + theta_shift))

  coords_param_list <- vector("list", major_rings)
  ring_ids <- vector("list", major_rings)
  ring_angles <- vector("list", major_rings)
  ring_sizes <- integer(major_rings)
  edges <- list()
  next_id <- 1L

  for (ring in seq_len(major_rings)) {
    theta <- theta_vals[[ring]]
    t <- (ring - 1L) / major_rings
    target_count <- tube_count * (1 + count_irregularity * sin(2 * pi * t + 0.4))
    count <- max(5L, as.integer(round(target_count)))
    phi <- seq(0, 2 * pi * (1 - 1 / count), length.out = count)
    phi <- sort(.normalize.angles(phi + phase_twist * sin(2 * theta)))

    ids <- next_id:(next_id + count - 1L)
    next_id <- next_id + count
    ring_ids[[ring]] <- as.integer(ids)
    ring_angles[[ring]] <- phi
    ring_sizes[[ring]] <- count
    coords_param_list[[ring]] <- cbind(theta = rep(theta, count), phi = phi)
    edges[[length(edges) + 1L]] <- .cyclic.ring.edges(ids)
  }

  for (ring in seq_len(major_rings)) {
    next_ring <- (ring %% major_rings) + 1L
    edges[[length(edges) + 1L]] <- .connect.cyclic.rings(
      ids_a = ring_ids[[ring]],
      angles_a = ring_angles[[ring]],
      ids_b = ring_ids[[next_ring]],
      angles_b = ring_angles[[next_ring]]
    )
  }

  coords_param <- do.call(rbind, coords_param_list)
  colnames(coords_param) <- c("theta", "phi")
  storage.mode(coords_param) <- "double"
  list(
    edges = .normalize_undirected_edges(.bind_edges(edges)),
    coords_param = coords_param,
    n = as.integer(nrow(coords_param)),
    major_rings = major_rings,
    ring_sizes = ring_sizes
  )
}

.double.torus.slice.components <- function(x,
                                           branch_length,
                                           transition_width,
                                           branch_offset,
                                           tube_radius) {
  branch_length <- .as_positive_scalar(branch_length, "branch_length")
  transition_width <- .as_positive_scalar(transition_width, "transition_width")
  branch_offset <- .as_positive_scalar(branch_offset, "branch_offset")
  tube_radius <- .as_positive_scalar(tube_radius, "tube_radius")
  x <- .as_finite_scalar(x, "x")

  u <- abs(x)
  branch_weight <- (branch_length + transition_width - u) / transition_width
  branch_weight <- min(max(branch_weight, 0), 1)
  if (branch_weight <= 0.12) {
    radius <- tube_radius * (0.92 + 0.55 * branch_weight)
    return(list(
      centers = matrix(c(0, 0), nrow = 1L, dimnames = list(NULL, c("y", "z"))),
      radii = radius,
      branch_weight = branch_weight
    ))
  }

  spread <- branch_offset * (0.35 + 0.65 * branch_weight)
  centers <- cbind(
    y = c(-spread, 0, spread),
    z = c(0, 0, 0)
  )
  radii <- tube_radius * c(0.96, 1.08, 0.96) * (0.95 + 0.08 * cos(pi * x))
  list(
    centers = centers,
    radii = radii,
    branch_weight = branch_weight
  )
}

.irregular.double.torus.canonical <- function(slices = 11,
                                              tube_count = 14,
                                              branch_length = 0.85,
                                              branch_offset = 0.72,
                                              tube_radius = 0.28,
                                              transition_width = 0.42,
                                              count_irregularity = 0.2,
                                              axial_irregularity = 0.3,
                                              phase_twist = 0.35) {
  slices <- .as_whole_number(slices, "slices", min = 7L)
  tube_count <- .as_whole_number(tube_count, "tube_count", min = 6L)
  branch_length <- .as_positive_scalar(branch_length, "branch_length")
  branch_offset <- .as_positive_scalar(branch_offset, "branch_offset")
  tube_radius <- .as_positive_scalar(tube_radius, "tube_radius")
  transition_width <- .as_positive_scalar(transition_width, "transition_width")
  count_irregularity <- .as_finite_scalar(count_irregularity, "count_irregularity")
  axial_irregularity <- .as_finite_scalar(axial_irregularity, "axial_irregularity")
  phase_twist <- .as_finite_scalar(phase_twist, "phase_twist")
  if (count_irregularity < 0 || count_irregularity >= 1) {
    stop("count_irregularity must be in [0, 1)", call. = FALSE)
  }
  if (axial_irregularity < 0 || axial_irregularity > 1) {
    stop("axial_irregularity must be in [0, 1]", call. = FALSE)
  }

  x_extent <- branch_length + transition_width + tube_radius
  gap <- (2 * x_extent) / (slices + 1L)
  x_base <- seq(-x_extent + gap / 2, x_extent - gap / 2, length.out = slices)
  x_shift <- axial_irregularity * 0.22 * gap *
    sin(seq_len(slices) * (2 * pi / (slices + 1L)))
  x_vals <- sort(x_base + x_shift)

  west_id <- 1L
  coords_list <- list(matrix(c(-x_extent, 0, 0), nrow = 1L,
                             dimnames = list(NULL, c("x", "y", "z"))))
  slice_specs <- vector("list", slices)
  slice_sizes <- integer(slices)
  slice_components <- integer(slices)
  edges <- list()
  next_id <- 2L

  for (slice in seq_len(slices)) {
    x <- x_vals[[slice]]
    spec <- .double.torus.slice.components(
      x = x,
      branch_length = branch_length,
      transition_width = transition_width,
      branch_offset = branch_offset,
      tube_radius = tube_radius
    )
    interval_specs <- vector("list", nrow(spec$centers))
    t <- (slice - 1L) / max(1L, slices - 1L)

    for (comp in seq_len(nrow(spec$centers))) {
      radius <- spec$radii[[comp]]
      phase <- phase_twist * (slice + 0.7 * comp)
      target_count <- tube_count * (radius / tube_radius) *
        (1 + count_irregularity * sin(2 * pi * t + 0.8 * comp))
      count <- max(6L, as.integer(round(target_count)))
      angles <- seq(0, 2 * pi * (1 - 1 / count), length.out = count)
      angles <- sort(.normalize.angles(
        angles + phase * sin(pi * (t - 0.5))
      ))

      ids <- next_id:(next_id + count - 1L)
      next_id <- next_id + count
      center <- spec$centers[comp, ]
      coords_list[[length(coords_list) + 1L]] <- cbind(
        x = rep(x, count),
        y = center[[1L]] + radius * cos(angles),
        z = center[[2L]] + radius * sin(angles)
      )
      interval_specs[[comp]] <- list(
        ids = as.integer(ids),
        angles = angles,
        center = as.double(center),
        radius = radius
      )
      edges[[length(edges) + 1L]] <- .cyclic.ring.edges(ids)
    }

    slice_specs[[slice]] <- list(
      x = x,
      components = interval_specs
    )
    slice_sizes[[slice]] <- sum(vapply(interval_specs,
                                       function(comp) length(comp$ids),
                                       integer(1L)))
    slice_components[[slice]] <- length(interval_specs)
  }

  if (slice_components[[1L]] != 1L || slice_components[[slices]] != 1L) {
    stop("double torus end slices must collapse to a single component; adjust parameters",
         call. = FALSE)
  }
  for (id in slice_specs[[1L]]$components[[1L]]$ids) {
    edges[[length(edges) + 1L]] <- c(west_id, id)
  }

  if (slices >= 2L) {
    for (slice in seq_len(slices - 1L)) {
      edges[[length(edges) + 1L]] <- .connect.cyclic.component_sets(
        comps_a = slice_specs[[slice]]$components,
        comps_b = slice_specs[[slice + 1L]]$components,
        overlap_scale = 1.1
      )
    }
  }

  east_id <- next_id
  coords_list[[length(coords_list) + 1L]] <- matrix(c(x_extent, 0, 0), nrow = 1L,
                                                    dimnames = list(NULL, c("x", "y", "z")))
  for (id in slice_specs[[slices]]$components[[1L]]$ids) {
    edges[[length(edges) + 1L]] <- c(id, east_id)
  }

  coords <- do.call(rbind, coords_list)
  storage.mode(coords) <- "double"
  list(
    edges = .normalize_undirected_edges(.bind_edges(edges)),
    coords = coords,
    n = as.integer(nrow(coords)),
    slices = slices,
    slice_sizes = slice_sizes,
    slice_components = slice_components,
    x_extent = x_extent
  )
}

.layered.radial.values <- function(layers,
                                   inner_radius,
                                   outer_radius,
                                   irregularity,
                                   include_center = FALSE) {
  layers <- .as_whole_number(layers, "layers", min = 1L)
  inner_radius <- .as_finite_scalar(inner_radius, "inner_radius")
  outer_radius <- .as_positive_scalar(outer_radius, "outer_radius")
  irregularity <- .as_finite_scalar(irregularity, "irregularity")
  if (inner_radius < 0 || inner_radius >= outer_radius) {
    stop("inner_radius must satisfy 0 <= inner_radius < outer_radius", call. = FALSE)
  }
  if (irregularity < 0 || irregularity > 1) {
    stop("irregularity must be in [0, 1]", call. = FALSE)
  }

  base <- if (include_center) {
    seq(outer_radius / layers, outer_radius, length.out = layers)
  } else {
    seq(inner_radius, outer_radius, length.out = layers)
  }
  if (layers >= 3L) {
    gap <- min(diff(base))
    shift <- irregularity * 0.22 * gap *
      sin(seq_len(layers) * (2 * pi / (layers + 1L)))
    if (include_center) {
      shift[[layers]] <- 0
      base <- sort(base + shift)
      if (base[[1L]] <= 0) {
        base[[1L]] <- max(outer_radius / (2 * layers), 1e-8)
      }
    } else {
      shift[[1L]] <- 0
      shift[[layers]] <- 0
      base <- sort(base + shift)
    }
  }
  as.double(base)
}

.rotate.z.coords <- function(coords, angle) {
  coords <- as.matrix(coords)
  angle <- .as_finite_scalar(angle, "angle")
  rot <- matrix(
    c(cos(angle), -sin(angle), 0,
      sin(angle),  cos(angle), 0,
      0,           0,          1),
    nrow = 3L,
    byrow = TRUE
  )
  coords %*% t(rot)
}

.layered.triangulated.solid.canonical <- function(base = c("tetrahedron", "octahedron", "icosahedron"),
                                                  level = 1,
                                                  layers = 3,
                                                  inner_radius = 0,
                                                  outer_radius = 1,
                                                  radial_irregularity = 0.25,
                                                  layer_twist = 0.35,
                                                  include_center = FALSE) {
  base <- match.arg(base)
  level <- .as_whole_number(level, "level")
  layers <- .as_whole_number(layers, "layers", min = 1L)
  layer_twist <- .as_finite_scalar(layer_twist, "layer_twist")

  surf <- .triangulated.polyhedron.canonical(base = base, level = level)
  dirs <- .normalize.row.coords(surf$coords)
  n_surface <- nrow(dirs)
  radii <- .layered.radial.values(
    layers = layers,
    inner_radius = inner_radius,
    outer_radius = outer_radius,
    irregularity = radial_irregularity,
    include_center = include_center
  )

  coords_list <- list()
  edges <- list()
  layer_offsets <- integer(layers)
  next_id <- 1L

  if (include_center) {
    coords_list[[1L]] <- matrix(c(0, 0, 0), nrow = 1L,
                                dimnames = list(NULL, c("x", "y", "z")))
    next_id <- 2L
  }

  for (layer in seq_len(layers)) {
    t <- if (layers == 1L) 1 else layer / layers
    angle <- layer_twist * sin(pi * (t - 0.5))
    dirs_layer <- .rotate.z.coords(dirs, angle)
    coords_layer <- dirs_layer * radii[[layer]]
    layer_offsets[[layer]] <- next_id - 1L
    coords_list[[length(coords_list) + 1L]] <- coords_layer
    edges[[length(edges) + 1L]] <- surf$edges + layer_offsets[[layer]]
    next_id <- next_id + n_surface
  }

  if (include_center) {
    first_ids <- layer_offsets[[1L]] + seq_len(n_surface)
    for (id in first_ids) {
      edges[[length(edges) + 1L]] <- c(1L, id)
    }
  }

  if (layers >= 2L) {
    face_template <- surf$faces
    for (layer in seq_len(layers - 1L)) {
      lower_ids <- layer_offsets[[layer]] + seq_len(n_surface)
      upper_ids <- layer_offsets[[layer + 1L]] + seq_len(n_surface)
      edges[[length(edges) + 1L]] <- cbind(lower_ids, upper_ids)
      for (face_idx in seq_len(nrow(face_template))) {
        face <- face_template[face_idx, ]
        lower <- layer_offsets[[layer]] + face
        upper <- layer_offsets[[layer + 1L]] + face
        edges[[length(edges) + 1L]] <- rbind(
          c(lower[[1L]], upper[[2L]]),
          c(lower[[1L]], upper[[3L]]),
          c(lower[[2L]], upper[[3L]])
        )
      }
    }
  }

  coords <- do.call(rbind, coords_list)
  storage.mode(coords) <- "double"
  list(
    edges = .normalize_undirected_edges(.bind_edges(edges)),
    coords = coords,
    n = as.integer(nrow(coords)),
    base = base,
    level = level,
    layers = layers,
    radii = radii,
    subdivision = surf$subdivision,
    n_surface = n_surface,
    include_center = include_center
  )
}

.irregular.ball.canonical <- function(base = c("tetrahedron", "octahedron", "icosahedron"),
                                      level = 1,
                                      layers = 3,
                                      outer_radius = 1,
                                      radial_irregularity = 0.25,
                                      layer_twist = 0.35) {
  .layered.triangulated.solid.canonical(
    base = match.arg(base),
    level = level,
    layers = layers,
    inner_radius = 0,
    outer_radius = outer_radius,
    radial_irregularity = radial_irregularity,
    layer_twist = layer_twist,
    include_center = TRUE
  )
}

.irregular.shell.canonical <- function(base = c("tetrahedron", "octahedron", "icosahedron"),
                                       level = 1,
                                       layers = 3,
                                       inner_radius = 0.45,
                                       outer_radius = 1,
                                       radial_irregularity = 0.25,
                                       layer_twist = 0.35) {
  .layered.triangulated.solid.canonical(
    base = match.arg(base),
    level = level,
    layers = layers,
    inner_radius = inner_radius,
    outer_radius = outer_radius,
    radial_irregularity = radial_irregularity,
    layer_twist = layer_twist,
    include_center = FALSE
  )
}

.path.bracketing.indices <- function(x, x_vals) {
  x_vals <- as.double(x_vals)
  n <- length(x_vals)
  if (n == 0L) {
    stop("x_vals must contain at least one point", call. = FALSE)
  }
  if (n == 1L) {
    return(c(1L, 1L))
  }
  if (x <= x_vals[[1L]]) {
    return(c(1L, 2L))
  }
  if (x >= x_vals[[n]]) {
    return(c(n - 1L, n))
  }
  lower <- max(which(x_vals <= x))
  if (isTRUE(all.equal(x, x_vals[[lower]], tolerance = 1e-12))) {
    return(c(lower, lower))
  }
  c(lower, lower + 1L)
}

.connect_open_paths <- function(ids_a, x_a, ids_b, x_b) {
  ids_a <- as.integer(ids_a)
  ids_b <- as.integer(ids_b)
  ord_a <- order(x_a)
  ord_b <- order(x_b)
  ids_a <- ids_a[ord_a]
  ids_b <- ids_b[ord_b]
  x_a <- as.double(x_a[ord_a])
  x_b <- as.double(x_b[ord_b])
  edges <- list()

  for (j in seq_along(ids_b)) {
    nbr <- .path.bracketing.indices(x_b[[j]], x_a)
    edges[[length(edges) + 1L]] <- c(ids_b[[j]], ids_a[nbr[[1L]]])
    if (nbr[[2L]] != nbr[[1L]]) {
      edges[[length(edges) + 1L]] <- c(ids_b[[j]], ids_a[nbr[[2L]]])
    }
  }
  for (i in seq_along(ids_a)) {
    nbr <- .path.bracketing.indices(x_a[[i]], x_b)
    edges[[length(edges) + 1L]] <- c(ids_a[[i]], ids_b[nbr[[1L]]])
    if (nbr[[2L]] != nbr[[1L]]) {
      edges[[length(edges) + 1L]] <- c(ids_a[[i]], ids_b[nbr[[2L]]])
    }
  }

  .normalize_undirected_edges(.bind_edges(edges))
}

.pair.of.pants.slice.intervals <- function(y,
                                           outer_radius,
                                           hole_radius,
                                           hole_offset,
                                           hole_height,
                                           tol = 1e-10) {
  outer_dx_sq <- outer_radius^2 - y^2
  if (outer_dx_sq <= 0) {
    return(matrix(numeric(), ncol = 2L))
  }
  intervals <- matrix(c(-sqrt(outer_dx_sq), sqrt(outer_dx_sq)), ncol = 2L)
  holes <- list()

  hole_y <- y - hole_height
  hole_dx_sq <- hole_radius^2 - hole_y^2
  if (hole_dx_sq > 0) {
    hole_dx <- sqrt(hole_dx_sq)
    holes <- list(
      c(-hole_offset - hole_dx, -hole_offset + hole_dx),
      c(hole_offset - hole_dx, hole_offset + hole_dx)
    )
  }

  subtract_interval <- function(current, hole) {
    out <- matrix(numeric(), ncol = 2L)
    for (i in seq_len(nrow(current))) {
      start <- current[i, 1L]
      end <- current[i, 2L]
      h_start <- hole[[1L]]
      h_end <- hole[[2L]]
      if (h_end <= start + tol || h_start >= end - tol) {
        out <- rbind(out, c(start, end))
      } else {
        if (h_start > start + tol) {
          out <- rbind(out, c(start, h_start))
        }
        if (h_end < end - tol) {
          out <- rbind(out, c(h_end, end))
        }
      }
    }
    out
  }

  for (hole in holes) {
    intervals <- subtract_interval(intervals, hole)
    if (nrow(intervals) == 0L) {
      break
    }
  }
  intervals[order(intervals[, 1L]), , drop = FALSE]
}

.sample.open.interval <- function(start,
                                  end,
                                  count,
                                  phase,
                                  irregularity) {
  count <- .as_whole_number(count, "count", min = 1L)
  phase <- .as_finite_scalar(phase, "phase")
  irregularity <- .as_finite_scalar(irregularity, "irregularity")
  if (count == 1L) {
    return((start + end) / 2)
  }
  s <- seq(0, 1, length.out = count)
  if (count > 2L && irregularity > 0) {
    delta <- 0.18 * irregularity * sin(pi * s) * cos(phase + 2 * pi * s)
    s <- pmin(pmax(s + delta / max(count - 1L, 1L), 0), 1)
    s[[1L]] <- 0
    s[[count]] <- 1
    s <- sort(s)
  }
  start + (end - start) * s
}

.irregular.pair.of.pants.canonical <- function(slices = 11,
                                               outer_count = 28,
                                               outer_radius = 1.1,
                                               hole_radius = 0.24,
                                               hole_offset = 0.38,
                                               hole_height = 0.18,
                                               count_irregularity = 0.2,
                                               vertical_irregularity = 0.35,
                                               phase_twist = 0.35) {
  slices <- .as_whole_number(slices, "slices", min = 5L)
  outer_count <- .as_whole_number(outer_count, "outer_count", min = 8L)
  outer_radius <- .as_positive_scalar(outer_radius, "outer_radius")
  hole_radius <- .as_positive_scalar(hole_radius, "hole_radius")
  hole_offset <- .as_positive_scalar(hole_offset, "hole_offset")
  hole_height <- .as_finite_scalar(hole_height, "hole_height")
  count_irregularity <- .as_finite_scalar(count_irregularity, "count_irregularity")
  vertical_irregularity <- .as_finite_scalar(vertical_irregularity, "vertical_irregularity")
  phase_twist <- .as_finite_scalar(phase_twist, "phase_twist")
  if (count_irregularity < 0 || count_irregularity >= 1) {
    stop("count_irregularity must be in [0, 1)", call. = FALSE)
  }
  if (vertical_irregularity < 0 || vertical_irregularity > 1) {
    stop("vertical_irregularity must be in [0, 1]", call. = FALSE)
  }
  if (hole_offset <= hole_radius) {
    stop("hole_offset must exceed hole_radius so the two holes remain disjoint",
         call. = FALSE)
  }
  hole_center_radius <- sqrt(hole_offset^2 + hole_height^2)
  if ((hole_center_radius + hole_radius) >= outer_radius) {
    stop("holes must lie strictly inside the outer boundary", call. = FALSE)
  }

  gap <- (2 * outer_radius) / (slices + 1L)
  y_base <- seq(outer_radius - gap / 2, -outer_radius + gap / 2, length.out = slices)
  y_shift <- vertical_irregularity * 0.22 * gap *
    sin(seq_len(slices) * (2 * pi / (slices + 1L)))
  y_vals <- sort(y_base + y_shift, decreasing = TRUE)
  base_spacing <- (2 * outer_radius) / outer_count

  coords_list <- list()
  slice_specs <- vector("list", slices)
  slice_sizes <- integer(slices)
  slice_components <- integer(slices)
  next_id <- 1L
  edges <- list()

  for (slice in seq_len(slices)) {
    y <- y_vals[[slice]]
    intervals <- .pair.of.pants.slice.intervals(
      y = y,
      outer_radius = outer_radius,
      hole_radius = hole_radius,
      hole_offset = hole_offset,
      hole_height = hole_height
    )
    if (nrow(intervals) == 0L) {
      stop("slice construction produced an empty cross-section; adjust parameters",
           call. = FALSE)
    }

    interval_specs <- vector("list", nrow(intervals))
    t <- (slice - 1L) / max(1L, slices - 1L)
    for (comp in seq_len(nrow(intervals))) {
      start <- intervals[comp, 1L]
      end <- intervals[comp, 2L]
      len <- end - start
      target <- outer_count * (len / (2 * outer_radius)) *
        (1 + count_irregularity * sin(2 * pi * t + 0.9 * comp))
      min_count <- if (len > 0.75 * base_spacing) 2L else 1L
      count <- max(min_count, as.integer(round(target)))
      x_vals <- .sample.open.interval(
        start = start,
        end = end,
        count = count,
        phase = phase_twist * (slice + comp),
        irregularity = count_irregularity
      )
      ids <- next_id:(next_id + count - 1L)
      next_id <- next_id + count
      coords_list[[length(coords_list) + 1L]] <- cbind(
        u = x_vals,
        v = rep(y, count)
      )
      interval_specs[[comp]] <- list(
        ids = as.integer(ids),
        x = as.double(x_vals),
        start = start,
        end = end
      )
      if (count >= 2L) {
        edges[[length(edges) + 1L]] <- cbind(ids[-count], ids[-1L])
      }
    }
    slice_specs[[slice]] <- list(
      y = y,
      intervals = interval_specs
    )
    slice_sizes[[slice]] <- sum(vapply(interval_specs,
                                       function(spec) length(spec$ids),
                                       integer(1L)))
    slice_components[[slice]] <- length(interval_specs)
  }

  for (slice in seq_len(slices - 1L)) {
    current <- slice_specs[[slice]]$intervals
    next_slice <- slice_specs[[slice + 1L]]$intervals
    for (i in seq_along(current)) {
      for (j in seq_along(next_slice)) {
        overlap <- min(current[[i]]$end, next_slice[[j]]$end) -
          max(current[[i]]$start, next_slice[[j]]$start)
        if (overlap >= (-0.15 * base_spacing)) {
          edges[[length(edges) + 1L]] <- .connect_open_paths(
            ids_a = current[[i]]$ids,
            x_a = current[[i]]$x,
            ids_b = next_slice[[j]]$ids,
            x_b = next_slice[[j]]$x
          )
        }
      }
    }
  }

  coords <- do.call(rbind, coords_list)
  colnames(coords) <- c("u", "v")
  storage.mode(coords) <- "double"
  list(
    edges = .normalize_undirected_edges(.bind_edges(edges)),
    coords = coords,
    n = as.integer(nrow(coords)),
    slices = slices,
    slice_sizes = slice_sizes,
    slice_components = slice_components
  )
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

.menger.sponge.mask <- function() {
  mask <- array(TRUE, dim = c(3L, 3L, 3L))
  for (row in seq_len(3L)) {
    for (col in seq_len(3L)) {
      for (layer in seq_len(3L)) {
        middle_count <- sum(c(row, col, layer) == 2L)
        if (middle_count >= 2L) {
          mask[row, col, layer] <- FALSE
        }
      }
    }
  }
  mask
}

.cube.band.starts <- function(side,
                              width,
                              period,
                              offset,
                              width_name = "width",
                              period_name = "period",
                              offset_name = "offset") {
  width <- .as_whole_number(width, width_name, min = 1L)
  if (width >= side) {
    stop(sprintf("%s must be < side", width_name), call. = FALSE)
  }
  period <- .as_whole_number(period, period_name, min = 1L)
  offset <- .as_whole_number(offset, offset_name, min = 1L)
  max_start <- side - width
  if (offset > max_start) {
    stop(sprintf("%s must be <= %d for the chosen side and width",
                 offset_name,
                 max_start),
         call. = FALSE)
  }
  as.integer(seq(offset, max_start, by = period))
}

.cube.band.mask <- function(side, width, starts) {
  side <- .as_whole_number(side, "side", min = 2L)
  width <- .as_whole_number(width, "width", min = 1L)
  starts <- as.integer(starts)
  out <- rep(FALSE, side)
  for (start in starts) {
    out[start:(start + width - 1L)] <- TRUE
  }
  out
}

.set_cube_block <- function(arr,
                            row_start,
                            col_start,
                            layer_start,
                            row_size,
                            col_size = row_size,
                            layer_size = row_size,
                            value = FALSE) {
  dims <- dim(arr)
  row_start <- .as_whole_number(row_start, "row_start", min = 1L)
  col_start <- .as_whole_number(col_start, "col_start", min = 1L)
  layer_start <- .as_whole_number(layer_start, "layer_start", min = 1L)
  row_size <- .as_whole_number(row_size, "row_size", min = 1L)
  col_size <- .as_whole_number(col_size, "col_size", min = 1L)
  layer_size <- .as_whole_number(layer_size, "layer_size", min = 1L)

  row_end <- row_start + row_size - 1L
  col_end <- col_start + col_size - 1L
  layer_end <- layer_start + layer_size - 1L
  if (row_end > dims[1L] || col_end > dims[2L] || layer_end > dims[3L]) {
    stop("requested cube block lies outside the keep-array bounds", call. = FALSE)
  }

  arr[row_start:row_end, col_start:col_end, layer_start:layer_end] <- value
  arr
}

.cube.periodic.tunnels.mask <- function(side = 5,
                                        tunnel_width = 1,
                                        tunnel_period = 2,
                                        tunnel_offset = 2) {
  side <- .as_whole_number(side, "side", min = 3L)
  starts <- .cube.band.starts(
    side = side,
    width = tunnel_width,
    period = tunnel_period,
    offset = tunnel_offset,
    width_name = "tunnel_width",
    period_name = "tunnel_period",
    offset_name = "tunnel_offset"
  )
  bands <- .cube.band.mask(side, tunnel_width, starts)
  grid <- expand.grid(
    row = seq_len(side),
    col = seq_len(side),
    layer = seq_len(side)
  )
  middle_count <- bands[grid$row] + bands[grid$col] + bands[grid$layer]
  keep <- array(TRUE, dim = c(side, side, side))
  keep[cbind(grid$row[middle_count >= 2L],
             grid$col[middle_count >= 2L],
             grid$layer[middle_count >= 2L])] <- FALSE
  keep
}

.cube.asymmetric.cavities.mask <- function(side = 5,
                                           cavity_size = 2,
                                           pocket_size = max(1L, cavity_size - 1L)) {
  side <- .as_whole_number(side, "side", min = 5L)
  cavity_size <- .as_whole_number(cavity_size, "cavity_size", min = 1L)
  pocket_size <- .as_whole_number(pocket_size, "pocket_size", min = 1L)
  max_interior <- side - 2L
  if (cavity_size > max_interior) {
    stop(sprintf("cavity_size must be <= %d for the chosen side", max_interior),
         call. = FALSE)
  }
  if (pocket_size > max_interior) {
    stop(sprintf("pocket_size must be <= %d for the chosen side", max_interior),
         call. = FALSE)
  }

  keep <- array(TRUE, dim = c(side, side, side))
  keep <- .set_cube_block(
    keep,
    row_start = 2L,
    col_start = max(2L, side - cavity_size),
    layer_start = 2L,
    row_size = cavity_size,
    col_size = cavity_size,
    layer_size = cavity_size,
    value = FALSE
  )
  keep <- .set_cube_block(
    keep,
    row_start = max(2L, side - pocket_size),
    col_start = 2L,
    layer_start = max(2L, side - pocket_size),
    row_size = pocket_size,
    col_size = pocket_size,
    layer_size = pocket_size,
    value = FALSE
  )
  keep
}

.cube.channel.network.mask <- function(side = 5,
                                       channel_width = 1,
                                       branch_offset = 2) {
  side <- .as_whole_number(side, "side", min = 5L)
  channel_width <- .as_whole_number(channel_width, "channel_width", min = 1L)
  if (channel_width >= side) {
    stop("channel_width must be < side", call. = FALSE)
  }
  max_start <- side - channel_width
  branch_offset <- .as_whole_number(branch_offset, "branch_offset", min = 2L)
  if (branch_offset > max_start) {
    stop(sprintf("branch_offset must be <= %d for the chosen side and channel_width",
                 max_start),
         call. = FALSE)
  }

  keep <- array(TRUE, dim = c(side, side, side))
  mid_start <- ((side - channel_width) %/% 2L) + 1L
  high_start <- side - channel_width

  keep <- .set_cube_block(
    keep,
    row_start = mid_start,
    col_start = 1L,
    layer_start = mid_start,
    row_size = channel_width,
    col_size = side,
    layer_size = channel_width,
    value = FALSE
  )
  keep <- .set_cube_block(
    keep,
    row_start = 1L,
    col_start = mid_start,
    layer_start = mid_start,
    row_size = side,
    col_size = channel_width,
    layer_size = channel_width,
    value = FALSE
  )
  keep <- .set_cube_block(
    keep,
    row_start = mid_start,
    col_start = mid_start,
    layer_start = 1L,
    row_size = channel_width,
    col_size = channel_width,
    layer_size = side,
    value = FALSE
  )
  keep <- .set_cube_block(
    keep,
    row_start = branch_offset,
    col_start = 1L,
    layer_start = mid_start,
    row_size = channel_width,
    col_size = side,
    layer_size = channel_width,
    value = FALSE
  )
  keep <- .set_cube_block(
    keep,
    row_start = mid_start,
    col_start = high_start,
    layer_start = 1L,
    row_size = channel_width,
    col_size = channel_width,
    layer_size = side,
    value = FALSE
  )
  keep
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

.sierpinski.triangle.canonical <- function(level) {
  level <- .as_whole_number(level, "level")

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
      matrix(c(0.5, 0), nrow(right$coords), 2L, byrow = TRUE)
    top_coords <- top$coords / 2 +
      matrix(c(0.25, sqrt(3) / 4), nrow(top$coords), 2L, byrow = TRUE)

    off1 <- left$n
    off2 <- left$n + right$n

    edges <- rbind(left$edges,
                   right$edges + off1,
                   top$edges + off2)
    coords <- rbind(left_coords, right_coords, top_coords)

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
    coords <- coords[ids, , drop = FALSE]
    edges <- .normalize_undirected_edges(edges)

    corners <- c(map[as.character(L[1L])],
                 map[as.character(R[2L])],
                 map[as.character(T[3L])])

    list(edges = edges, coords = coords, corners = corners, n = length(ids))
  }

  out <- build(level)
  coords <- out$coords
  colnames(coords) <- c("u", "v")
  storage.mode(coords) <- "double"
  list(
    edges = out$edges,
    coords = coords,
    corners = out$corners,
    n = as.integer(nrow(coords))
  )
}

.sierpinski.tetrahedron.canonical <- function(level) {
  level <- .as_whole_number(level, "level")

  coords_list <- vector("list", 4L)
  coords_list[[1L]] <- c(0, 0, 0)
  coords_list[[2L]] <- c(1, 0, 0)
  coords_list[[3L]] <- c(0.5, sqrt(3) / 2, 0)
  coords_list[[4L]] <- c(0.5, sqrt(3) / 6, sqrt(2 / 3))

  state <- new.env(parent = emptyenv())
  state$edges <- list()
  state$next_id <- 5L

  add_tetrahedron <- function(a, b, c, d) {
    state$edges[[length(state$edges) + 1L]] <<- rbind(
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

  edges <- .normalize_undirected_edges(.bind_edges(state$edges))
  coords <- do.call(rbind, coords_list[seq_len(state$next_id - 1L)])
  colnames(coords) <- c("x", "y", "z")
  storage.mode(coords) <- "double"
  list(
    edges = edges,
    coords = coords,
    n = as.integer(nrow(coords))
  )
}

.tetrahedron.classic.mask <- function() {
  stats::setNames(rep(TRUE, 4L), c("base_left", "base_right", "base_back", "apex"))
}

.tetrahedron.corner.missing.mask <- function(omit = c("apex", "base_left",
                                                      "base_right", "base_back")) {
  omit <- .as_named_choice(
    omit[[1L]],
    c("apex", "base_left", "base_right", "base_back"),
    "omit"
  )
  keep <- .tetrahedron.classic.mask()
  keep[[omit]] <- FALSE
  keep
}

.tetrahedron.mask.is_classic <- function(mask) {
  identical(
    unname(.as_tetrahedron_keep_mask(mask, "mask")),
    unname(.tetrahedron.classic.mask())
  )
}

.tetrahedron.mask.label <- function(mask) {
  keep <- .as_tetrahedron_keep_mask(mask, "mask")
  if (identical(unname(keep), unname(.tetrahedron.classic.mask()))) {
    return("Sierpinski tetrahedron mask")
  }
  if (sum(keep) == 3L) {
    omitted <- names(keep)[!keep][1L]
    return(sprintf("Corner-mask tetrahedron (%s omitted)", omitted))
  }
  sprintf("Tetrahedron mask (%s)", paste(names(keep)[keep], collapse = ", "))
}

.recursive.tetrahedron.mask.canonical <- function(mask, level) {
  mask <- .as_tetrahedron_keep_mask(mask, "mask")
  level <- .as_whole_number(level, "level")
  if (.tetrahedron.mask.is_classic(mask)) {
    out <- .sierpinski.tetrahedron.canonical(level)
    out$mask <- mask
    return(out)
  }

  base_vertices <- rbind(
    c(0, 0, 0),
    c(1, 0, 0),
    c(0.5, sqrt(3) / 2, 0),
    c(0.5, sqrt(3) / 6, sqrt(2 / 3))
  )
  base_edges <- rbind(
    c(1L, 2L), c(1L, 3L), c(1L, 4L),
    c(2L, 3L), c(2L, 4L), c(3L, 4L)
  )

  build <- function(current_level, vertices) {
    if (current_level == 0L) {
      coords <- as.matrix(vertices)
      storage.mode(coords) <- "double"
      return(list(edges = base_edges, coords = coords, n = 4L))
    }

    a <- vertices[1L, ]
    b <- vertices[2L, ]
    c <- vertices[3L, ]
    d <- vertices[4L, ]
    ab <- (a + b) / 2
    ac <- (a + c) / 2
    ad <- (a + d) / 2
    bc <- (b + c) / 2
    bd <- (b + d) / 2
    cd <- (c + d) / 2
    child_vertices <- list(
      base_left = rbind(a, ab, ac, ad),
      base_right = rbind(ab, b, bc, bd),
      base_back = rbind(ac, bc, c, cd),
      apex = rbind(ad, bd, cd, d)
    )

    edges <- list()
    coords_list <- list()
    offset <- 0L
    active_slots <- names(mask)[mask]
    for (slot in active_slots) {
      child <- build(current_level - 1L, child_vertices[[slot]])
      edges[[length(edges) + 1L]] <- child$edges + offset
      coords_list[[length(coords_list) + 1L]] <- child$coords
      offset <- offset + child$n
    }

    merged <- .deduplicate.coordinate.graph(
      edges = .bind_edges(edges),
      coords = do.call(rbind, coords_list)
    )
    list(
      edges = merged$edges,
      coords = merged$coords,
      n = as.integer(nrow(merged$coords))
    )
  }

  out <- build(level, base_vertices)
  colnames(out$coords) <- c("x", "y", "z")
  storage.mode(out$coords) <- "double"
  out$mask <- mask
  out
}

.triangle.classic.mask <- function() {
  stats::setNames(c(TRUE, TRUE, TRUE, FALSE), c("left", "right", "top", "center"))
}

.triangle.bridge.mask <- function(missing = c("top", "left", "right")) {
  missing <- .as_named_choice(missing[[1L]], c("top", "left", "right"), "missing")
  keep <- stats::setNames(rep(TRUE, 4L), c("left", "right", "top", "center"))
  keep[[missing]] <- FALSE
  keep
}

.triangle.mask.is_classic <- function(mask) {
  identical(
    unname(.as_triangle_keep_mask(mask, "mask")),
    unname(.triangle.classic.mask())
  )
}

.triangle.mask.label <- function(mask) {
  keep <- .as_triangle_keep_mask(mask, "mask")
  if (identical(unname(keep), unname(.triangle.classic.mask()))) {
    return("Sierpinski triangle mask")
  }
  if (isTRUE(keep[["center"]]) && sum(keep) == 3L) {
    omitted <- names(keep)[!keep][1L]
    return(sprintf("Bridge triangle mask (%s omitted)", omitted))
  }
  sprintf("Triangle mask (%s)", paste(names(keep)[keep], collapse = ", "))
}

.deduplicate.coordinate.graph <- function(edges,
                                         coords,
                                         digits = 12L,
                                         return_map = FALSE) {
  digits <- .as_whole_number(digits, "digits", min = 1L)
  coords <- as.matrix(coords)
  if (nrow(coords) == 0L) {
    out <- list(edges = .empty_edge_matrix(), coords = coords)
    if (return_map) {
      out$old_to_new <- integer(0L)
      out$keep_rows <- logical(0L)
    }
    return(out)
  }
  rounded <- round(coords, digits)
  rounded[abs(rounded) < (10^(-digits))] <- 0
  key_parts <- apply(
    rounded,
    2L,
    function(col) formatC(col, digits = digits, format = "f")
  )
  keys <- apply(key_parts, 1L, paste, collapse = ":")
  keep_rows <- !duplicated(keys)
  old_to_new <- match(keys, keys[keep_rows])
  coords_dedup <- coords[keep_rows, , drop = FALSE]
  edges_dedup <- if (nrow(edges) == 0L) {
    .empty_edge_matrix()
  } else {
    cbind(old_to_new[edges[, 1L]], old_to_new[edges[, 2L]])
  }
  storage.mode(coords_dedup) <- "double"
  out <- list(
    edges = .normalize_undirected_edges(edges_dedup),
    coords = coords_dedup
  )
  if (return_map) {
    out$old_to_new <- as.integer(old_to_new)
    out$keep_rows <- keep_rows
  }
  out
}

.recursive.triangle.mask.canonical <- function(mask, level) {
  mask <- .as_triangle_keep_mask(mask, "mask")
  level <- .as_whole_number(level, "level")
  if (.triangle.mask.is_classic(mask)) {
    out <- .sierpinski.triangle.canonical(level)
    out$mask <- mask
    return(out)
  }

  height <- sqrt(3) / 2
  base_vertices <- rbind(
    c(0, 0),
    c(1, 0),
    c(0.5, height)
  )
  base_edges <- rbind(c(1L, 2L), c(2L, 3L), c(3L, 1L))

  build <- function(current_level, vertices) {
    if (current_level == 0L) {
      coords <- as.matrix(vertices)
      storage.mode(coords) <- "double"
      return(list(edges = base_edges, coords = coords, n = 3L))
    }

    a <- vertices[1L, ]
    b <- vertices[2L, ]
    c <- vertices[3L, ]
    ab <- (a + b) / 2
    bc <- (b + c) / 2
    ac <- (a + c) / 2
    child_vertices <- list(
      left = rbind(a, ab, ac),
      right = rbind(ab, b, bc),
      top = rbind(ac, bc, c),
      center = rbind(ab, bc, ac)
    )

    edges <- list()
    coords_list <- list()
    offset <- 0L
    active_slots <- names(mask)[mask]
    for (slot in active_slots) {
      child <- build(current_level - 1L, child_vertices[[slot]])
      edges[[length(edges) + 1L]] <- child$edges + offset
      coords_list[[length(coords_list) + 1L]] <- child$coords
      offset <- offset + child$n
    }

    merged <- .deduplicate.coordinate.graph(
      edges = .bind_edges(edges),
      coords = do.call(rbind, coords_list)
    )
    list(
      edges = merged$edges,
      coords = merged$coords,
      n = as.integer(nrow(merged$coords))
    )
  }

  out <- build(level, base_vertices)
  colnames(out$coords) <- c("u", "v")
  storage.mode(out$coords) <- "double"
  out$mask <- mask
  out
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

#' Cube mask pattern helpers for recursive porous families
#'
#' Convenience constructors for cubic keep-arrays that can be passed to
#' \code{\link{edges.recursive.cube.mask}()} or
#' \code{\link{recursive.cube.mask.surface.graph}()} to build cubical porous
#' benchmark families. These helpers encode three qualitatively different void
#' patterns: repeated tunnel lattices, interior asymmetric cavities, and
#' deterministic channel networks.
#'
#' The returned arrays use the same orientation as the recursive cube-mask
#' helpers: the first dimension runs from top to bottom, the second from left
#' to right, and the third from front to back.
#'
#' @param side Side length of the cubic keep-array.
#' @param tunnel_width Width of each removed tunnel band.
#' @param tunnel_period Spacing between successive tunnel bands.
#' @param tunnel_offset Starting index of the first tunnel band.
#' @param cavity_size Side length of the larger interior cavity block.
#' @param pocket_size Side length of the smaller secondary cavity block.
#' @param channel_width Width of each removed channel in the channel-network
#'   family.
#' @param branch_offset Interior offset of the extra branch channel in
#'   \code{mask.cube.channel.network()}.
#'
#' @return A logical cubic keep-array.
#' @name cube_mask_pattern_helpers
NULL

#' @rdname cube_mask_pattern_helpers
#' @export
mask.cube.periodic.tunnels <- function(side = 5,
                                       tunnel_width = 1,
                                       tunnel_period = 2,
                                       tunnel_offset = 2) {
  .cube.periodic.tunnels.mask(
    side = side,
    tunnel_width = tunnel_width,
    tunnel_period = tunnel_period,
    tunnel_offset = tunnel_offset
  )
}

#' @rdname cube_mask_pattern_helpers
#' @export
mask.cube.asymmetric.cavities <- function(side = 5,
                                          cavity_size = 2,
                                          pocket_size = max(1L, cavity_size - 1L)) {
  .cube.asymmetric.cavities.mask(
    side = side,
    cavity_size = cavity_size,
    pocket_size = pocket_size
  )
}

#' @rdname cube_mask_pattern_helpers
#' @export
mask.cube.channel.network <- function(side = 5,
                                      channel_width = 1,
                                      branch_offset = 2) {
  .cube.channel.network.mask(
    side = side,
    channel_width = channel_width,
    branch_offset = branch_offset
  )
}

#' Triangle mask helpers for recursive gasket families
#'
#' Convenience constructors for the four-slot triangle masks used by
#' \code{\link{edges.recursive.triangle.mask}()} and
#' \code{\link{recursive.triangle.mask.surface.graph}()}. The slots correspond
#' to the three corner subtriangles and the central inverted subtriangle created
#' by one barycentric refinement of an equilateral triangle.
#'
#' The returned mask entries are ordered and named as \code{left},
#' \code{right}, \code{top}, and \code{center}.
#'
#' @param missing Which outer corner is omitted in
#'   \code{mask.triangle.bridge()}.
#'
#' @return A named logical vector of length 4.
#' @name triangle_mask_helpers
NULL

#' @rdname triangle_mask_helpers
#' @export
mask.triangle.classic <- function() {
  .triangle.classic.mask()
}

#' @rdname triangle_mask_helpers
#' @export
mask.triangle.bridge <- function(missing = c("top", "left", "right")) {
  .triangle.bridge.mask(match.arg(missing))
}

#' Tetrahedron mask helpers for recursive gasket families
#'
#' Convenience constructors for the four-corner tetrahedron masks used by
#' \code{\link{edges.recursive.tetrahedron.mask}()} and
#' \code{\link{recursive.tetrahedron.mask.surface.graph}()}. The slots
#' correspond to the four corner subtetrahedra created by one barycentric
#' refinement of a tetrahedron.
#'
#' The returned mask entries are ordered and named as \code{base_left},
#' \code{base_right}, \code{base_back}, and \code{apex}.
#'
#' @param omit Which tetrahedron corner is omitted in
#'   \code{mask.tetrahedron.corner.missing()}.
#'
#' @return A named logical vector of length 4.
#' @name tetrahedron_mask_helpers
NULL

#' @rdname tetrahedron_mask_helpers
#' @export
mask.tetrahedron.classic <- function() {
  .tetrahedron.classic.mask()
}

#' @rdname tetrahedron_mask_helpers
#' @export
mask.tetrahedron.corner.missing <- function(
    omit = c("apex", "base_left", "base_right", "base_back")) {
  .tetrahedron.corner.missing.mask(match.arg(omit))
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

.normalize.edge.weights <- function(edge_weights,
                                    normalize = c("median", "mean", "none")) {
  normalize <- match.arg(normalize)
  if (!is.numeric(edge_weights)) {
    stop("edge_weights must be numeric", call. = FALSE)
  }
  edge_weights <- as.double(edge_weights)
  if (length(edge_weights) == 0L) {
    return(list(edge_weights = numeric(0L), weight_scale = 1))
  }
  if (any(!is.finite(edge_weights) | edge_weights <= 0)) {
    stop("edge_weights must all be finite and > 0", call. = FALSE)
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

.kary.tree.vertex.count <- function(k, depth) {
  k <- .as_whole_number(k, "k", min = 1L)
  depth <- .as_whole_number(depth, "depth")
  if (depth == 0L) {
    return(1L)
  }
  if (k == 1L) {
    return(as.integer(depth + 1L))
  }
  count <- (k^(depth + 1L) - 1) / (k - 1L)
  if (!is.finite(count) || count > .Machine$integer.max) {
    stop("k and depth produce too many tree vertices", call. = FALSE)
  }
  as.integer(round(count))
}

.kary.tree.structure <- function(k = 2, depth = 2) {
  k <- .as_whole_number(k, "k", min = 1L)
  depth <- .as_whole_number(depth, "depth")
  n <- .kary.tree.vertex.count(k, depth)
  if (depth == 0L) {
    return(list(
      edges = .empty_edge_matrix(),
      n = 1L,
      parent = 0L,
      vertex_depth = 0L,
      edge_parent = integer(0L),
      edge_child = integer(0L),
      branch_index = integer(0L),
      parent_depth = integer(0L),
      child_depth = integer(0L)
    ))
  }

  parent <- integer(n)
  vertex_depth <- integer(n)
  edge_count <- n - 1L
  edge_parent <- integer(edge_count)
  edge_child <- integer(edge_count)
  branch_index <- integer(edge_count)
  parent_depth <- integer(edge_count)
  child_depth <- integer(edge_count)

  current <- 1L
  next_id <- 2L
  edge_id <- 1L
  for (d in seq_len(depth)) {
    new_parents <- integer(length(current) * k)
    new_idx <- 1L
    for (p in current) {
      for (slot in seq_len(k)) {
        child <- next_id
        next_id <- next_id + 1L
        parent[[child]] <- p
        vertex_depth[[child]] <- d
        edge_parent[[edge_id]] <- p
        edge_child[[edge_id]] <- child
        branch_index[[edge_id]] <- slot
        parent_depth[[edge_id]] <- d - 1L
        child_depth[[edge_id]] <- d
        new_parents[[new_idx]] <- child
        new_idx <- new_idx + 1L
        edge_id <- edge_id + 1L
      }
    }
    current <- new_parents
  }

  edges <- cbind(edge_parent, edge_child)
  storage.mode(edges) <- "integer"
  list(
    edges = edges,
    n = n,
    parent = parent,
    vertex_depth = vertex_depth,
    edge_parent = edge_parent,
    edge_child = edge_child,
    branch_index = branch_index,
    parent_depth = parent_depth,
    child_depth = child_depth
  )
}

.kary.tree.depth.factors <- function(depth,
                                     depth_rule = c("geometric", "constant", "custom"),
                                     depth_decay = 0.85,
                                     depth_factors = NULL) {
  depth <- .as_whole_number(depth, "depth")
  depth_rule <- .as_named_choice(depth_rule[[1L]],
                                 c("geometric", "constant", "custom"),
                                 "depth_rule")
  if (depth == 0L) {
    return(list(rule = depth_rule, factors = numeric(0L)))
  }

  factors <- switch(
    depth_rule,
    geometric = {
      depth_decay <- .as_positive_scalar(depth_decay, "depth_decay")
      depth_decay^(seq_len(depth) - 1L)
    },
    constant = rep(1, depth),
    custom = {
      if (is.null(depth_factors)) {
        stop("depth_factors must be supplied when depth_rule = \"custom\"",
             call. = FALSE)
      }
      if (!is.numeric(depth_factors) || any(is.na(depth_factors)) ||
          any(!is.finite(depth_factors)) || any(depth_factors <= 0)) {
        stop("depth_factors must be a positive finite numeric vector",
             call. = FALSE)
      }
      if (length(depth_factors) == 1L) {
        rep(as.double(depth_factors), depth)
      } else if (length(depth_factors) == depth) {
        as.double(depth_factors)
      } else {
        stop(sprintf("depth_factors must have length 1 or %d", depth),
             call. = FALSE)
      }
    }
  )
  list(rule = depth_rule, factors = as.double(factors))
}

.kary.tree.branch.factors <- function(k,
                                      branch_rule = c("linear", "uniform", "custom"),
                                      branch_spread = 0.3,
                                      branch_factors = NULL) {
  k <- .as_whole_number(k, "k", min = 1L)
  branch_rule <- .as_named_choice(branch_rule[[1L]],
                                  c("linear", "uniform", "custom"),
                                  "branch_rule")

  factors <- switch(
    branch_rule,
    linear = {
      branch_spread <- .as_finite_scalar(branch_spread, "branch_spread")
      if (branch_spread < 0) {
        stop("branch_spread must be >= 0", call. = FALSE)
      }
      vals <- seq(1 - branch_spread / 2, 1 + branch_spread / 2, length.out = k)
      if (any(vals <= 0)) {
        stop("branch_spread is too large and produces non-positive branch factors",
             call. = FALSE)
      }
      vals
    },
    uniform = rep(1, k),
    custom = {
      if (is.null(branch_factors)) {
        stop("branch_factors must be supplied when branch_rule = \"custom\"",
             call. = FALSE)
      }
      if (!is.numeric(branch_factors) || any(is.na(branch_factors)) ||
          any(!is.finite(branch_factors)) || any(branch_factors <= 0)) {
        stop("branch_factors must be a positive finite numeric vector",
             call. = FALSE)
      }
      if (length(branch_factors) == 1L) {
        rep(as.double(branch_factors), k)
      } else if (length(branch_factors) == k) {
        as.double(branch_factors)
      } else {
        stop(sprintf("branch_factors must have length 1 or %d", k),
             call. = FALSE)
      }
    }
  )
  list(rule = branch_rule, factors = as.double(factors))
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
#' family, \code{\link{edges.recursive.triangle.mask}()} builds a generic
#' triangle-mask family, \code{\link{edges.recursive.tetrahedron.mask}()}
#' builds a generic tetrahedron-mask family,
#' \code{\link{edges.recursive.cube.mask}()} builds a generic cube-mask
#' family, \code{\link{edges.vicsek}()} builds the connected axial-cross
#' variant, \code{\link{edges.menger.sponge}()} builds the classic cubical
#' sponge variant, \code{\link{edges.triangulated.polyhedron}()} builds a
#' generic irregular triangulated-surface family,
#' \code{\link{edges.sierpinski.triangle}()} builds the 2-simplex family,
#' \code{\link{edges.sierpinski.tetrahedron}()} builds the 3-simplex family,
#' and \code{\link{edges.sierpinski.carpet}()} builds a 2D cell-adjacency
#' carpet graph.
#'
#' @return A two-column integer matrix of undirected edges. Vertex labels are
#'   consecutive integers starting at 1.
#' @param tunnel_width Width of each removed tunnel band in
#'   \code{edges.cube.periodic.tunnels()}.
#' @param tunnel_period Spacing between successive tunnel bands in
#'   \code{edges.cube.periodic.tunnels()}.
#' @param tunnel_offset Starting index of the first tunnel band in
#'   \code{edges.cube.periodic.tunnels()}.
#' @param cavity_size Side length of the larger interior cavity block in
#'   \code{edges.cube.asymmetric.cavities()}.
#' @param pocket_size Side length of the smaller secondary cavity block in
#'   \code{edges.cube.asymmetric.cavities()}.
#' @param channel_width Width of each removed channel in
#'   \code{edges.cube.channel.network()}.
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

#' Weighted irregular torus surface helpers
#'
#' Convenience helpers that build a torus from deterministically irregular
#' major-cycle rings with varying sample counts, stitch adjacent rings into a
#' locally triangulated toroidal graph, and then realize that graph in
#' \eqn{\mathbb{R}^3} using either a standard, pinched, or wavy tube geometry.
#'
#' @param major_rings Number of cyclic major rings around the torus.
#' @param tube_count Approximate number of vertices around each minor cycle.
#' @param count_irregularity Irregularity level for the per-ring sample counts.
#'   Must lie in \code{[0, 1)}.
#' @param major_irregularity Irregularity level for the major-angle ring
#'   spacing. Must lie in \code{[0, 1]}.
#' @param phase_twist Finite phase offset used to desynchronize neighboring
#'   minor cycles.
#' @param surface Torus geometry family. One of \code{"standard"},
#'   \code{"pinched"}, or \code{"wavy"}.
#' @param major_radius Positive distance from the torus center to the center of
#'   the tube.
#' @param minor_radius Positive baseline radius of the torus tube.
#' @param amplitude Finite deformation amplitude.
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
#' \code{irregular.torus.surface.embedding()} returns an \code{n x 3} numeric
#' matrix with columns \code{x}, \code{y}, and \code{z}.
#'
#' \code{irregular.torus.surface.graph()} returns a list with components:
#' \itemize{
#'   \item \code{edges}: the irregular torus edges,
#'   \item \code{n}: number of vertices,
#'   \item \code{edge_weights}: induced positive edge lengths,
#'   \item \code{coords_surface}: the 3D embedding,
#'   \item \code{coords_param}: the irregular angular parameter coordinates,
#'   \item \code{weight_scale}: the normalization constant applied to the raw
#'     edge lengths,
#'   \item \code{family}: always \code{"irregular.torus"},
#'   \item \code{surface}: the chosen surface name,
#'   \item \code{major_rings}: number of major-cycle rings,
#'   \item \code{ring_sizes}: sample counts around each ring,
#'   \item \code{label}: a human-readable family label.
#' }
#'
#' @name irregular_torus_surface_helpers
NULL

#' @rdname irregular_torus_surface_helpers
#' @export
irregular.torus.surface.embedding <- function(
    major_rings = 8,
    tube_count = 16,
    count_irregularity = 0.2,
    major_irregularity = 0.25,
    phase_twist = 0.35,
    surface = c("standard", "pinched", "wavy"),
    major_radius = 2,
    minor_radius = 0.75,
    amplitude = 0.2,
    freq_major = 2,
    freq_minor = 1,
    twist = 0.25) {
  surface <- match.arg(surface)
  major_radius <- .as_positive_scalar(major_radius, "major_radius")
  minor_radius <- .as_positive_scalar(minor_radius, "minor_radius")
  amplitude <- .as_finite_scalar(amplitude, "amplitude")
  freq_major <- .as_positive_scalar(freq_major, "freq_major")
  freq_minor <- .as_positive_scalar(freq_minor, "freq_minor")
  twist <- .as_finite_scalar(twist, "twist")

  built <- .irregular.torus.canonical(
    major_rings = major_rings,
    tube_count = tube_count,
    count_irregularity = count_irregularity,
    major_irregularity = major_irregularity,
    phase_twist = phase_twist
  )
  theta <- built$coords_param[, 1L]
  phi <- built$coords_param[, 2L]
  phi_eff <- phi + twist * sin(freq_major * theta)
  local_minor_radius <- switch(
    surface,
    standard = rep(minor_radius, length(phi)),
    pinched = minor_radius * (1 - amplitude * cos(theta)),
    wavy = minor_radius * (1 + amplitude *
                             sin(freq_major * theta) * cos(freq_minor * phi))
  )
  if (any(!is.finite(local_minor_radius) | local_minor_radius <= 0)) {
    stop("irregular torus parameters produce a non-positive tube radius; adjust amplitude or frequencies",
         call. = FALSE)
  }
  if (any(local_minor_radius >= major_radius)) {
    stop("irregular torus parameters require major_radius to exceed the local tube radius everywhere",
         call. = FALSE)
  }

  ring_radius <- major_radius + local_minor_radius * cos(phi_eff)
  coords <- cbind(
    x = ring_radius * cos(theta),
    y = ring_radius * sin(theta),
    z = local_minor_radius * sin(phi_eff)
  )
  storage.mode(coords) <- "double"
  coords
}

#' @rdname irregular_torus_surface_helpers
#' @export
irregular.torus.surface.graph <- function(
    major_rings = 8,
    tube_count = 16,
    count_irregularity = 0.2,
    major_irregularity = 0.25,
    phase_twist = 0.35,
    surface = c("standard", "pinched", "wavy"),
    major_radius = 2,
    minor_radius = 0.75,
    amplitude = 0.2,
    freq_major = 2,
    freq_minor = 1,
    twist = 0.25,
    normalize = c("median", "mean", "none")) {
  surface <- match.arg(surface)
  normalize <- match.arg(normalize)
  built <- .irregular.torus.canonical(
    major_rings = major_rings,
    tube_count = tube_count,
    count_irregularity = count_irregularity,
    major_irregularity = major_irregularity,
    phase_twist = phase_twist
  )
  coords_surface <- irregular.torus.surface.embedding(
    major_rings = major_rings,
    tube_count = tube_count,
    count_irregularity = count_irregularity,
    major_irregularity = major_irregularity,
    phase_twist = phase_twist,
    surface = surface,
    major_radius = major_radius,
    minor_radius = minor_radius,
    amplitude = amplitude,
    freq_major = freq_major,
    freq_minor = freq_minor,
    twist = twist
  )
  weights <- .edge.weights.from.embedding(
    edges = built$edges,
    coords = coords_surface,
    normalize = normalize
  )

  out <- list(
    edges = built$edges,
    n = built$n,
    edge_weights = weights$edge_weights,
    coords_surface = coords_surface,
    coords_param = built$coords_param,
    weight_scale = weights$weight_scale,
    family = "irregular.torus",
    surface = surface,
    major_rings = built$major_rings,
    ring_sizes = built$ring_sizes,
    normalize = normalize,
    label = sprintf("%s irregular torus rings %d",
                    tools::toTitleCase(surface),
                    built$major_rings)
  )
  class(out) <- c("grip_irregular_torus_surface_graph", "list")
  out
}

.irregular.double.torus.surface.coords <- function(coords_param,
                                                   surface,
                                                   amplitude,
                                                   freq_x,
                                                   freq_theta,
                                                   twist) {
  surface <- .as_named_choice(surface,
                              c("standard", "bulged", "twisted", "wavy"),
                              "surface")
  amplitude <- .as_finite_scalar(amplitude, "amplitude")
  freq_x <- .as_positive_scalar(freq_x, "freq_x")
  freq_theta <- .as_positive_scalar(freq_theta, "freq_theta")
  twist <- .as_finite_scalar(twist, "twist")

  coords_param <- as.matrix(coords_param)
  x_extent <- max(abs(coords_param[, 1L]))
  x_norm <- if (x_extent > 0) coords_param[, 1L] / x_extent else rep(0, nrow(coords_param))
  theta <- atan2(coords_param[, 3L], coords_param[, 2L])

  coords <- switch(
    surface,
    standard = coords_param,
    bulged = {
      scale <- 1 + amplitude * (1 - x_norm^2)
      if (any(!is.finite(scale) | scale <= 0)) {
        stop("double torus bulged parameters require positive radial scale; adjust amplitude",
             call. = FALSE)
      }
      cbind(
        x = coords_param[, 1L],
        y = coords_param[, 2L] * scale,
        z = coords_param[, 3L] * scale
      )
    },
    twisted = {
      angle <- twist * sin(pi * x_norm)
      cbind(
        x = coords_param[, 1L],
        y = coords_param[, 2L] * cos(angle) - coords_param[, 3L] * sin(angle),
        z = coords_param[, 2L] * sin(angle) + coords_param[, 3L] * cos(angle)
      )
    },
    wavy = {
      scale <- 1 + amplitude * sin(freq_x * pi * x_norm) * cos(freq_theta * theta)
      if (any(!is.finite(scale) | scale <= 0)) {
        stop("double torus wavy parameters require positive radial scale; adjust amplitude or frequencies",
             call. = FALSE)
      }
      cbind(
        x = coords_param[, 1L],
        y = coords_param[, 2L] * scale,
        z = coords_param[, 3L] * scale
      )
    }
  )

  storage.mode(coords) <- "double"
  coords
}

#' Weighted irregular double-torus surface helpers
#'
#' Convenience helpers that build a closed genus-2 surface from a
#' deterministically irregular slice family whose cyclic cross-sections follow a
#' `1 -> 3 -> 1` loop transition between two poles. This gives a non-lattice,
#' point-sampled double-torus graph together with either the canonical 3D
#' embedding or simple bulged, twisted, and wavy deformations of that geometry.
#'
#' @param slices Number of non-pole slices through the double torus.
#' @param tube_count Approximate number of vertices around each tube-like loop.
#' @param branch_length Half-length of the three-loop central region.
#' @param branch_offset Offset of the outer loop centers from the central loop.
#' @param tube_radius Baseline radius of each tube-like loop.
#' @param transition_width Width of the left and right transition regions
#'   between the single-loop and three-loop slices.
#' @param count_irregularity Irregularity level for the per-component sample
#'   counts. Must lie in \code{[0, 1)}.
#' @param axial_irregularity Irregularity level for the slice spacing. Must lie
#'   in \code{[0, 1]}.
#' @param phase_twist Finite phase offset used to desynchronize neighboring
#'   cyclic slices.
#' @param surface Double-torus geometry family. One of \code{"standard"},
#'   \code{"bulged"}, \code{"twisted"}, or \code{"wavy"}.
#' @param amplitude Finite deformation amplitude.
#' @param freq_x Positive modulation frequency along the axial direction. Used
#'   only when \code{surface = "wavy"}.
#' @param freq_theta Positive angular modulation frequency around the local
#'   tube direction. Used only when \code{surface = "wavy"}.
#' @param twist Finite twist strength used only when
#'   \code{surface = "twisted"}.
#' @param normalize Normalization applied to the induced edge lengths. One of
#'   \code{"median"}, \code{"mean"}, or \code{"none"}.
#'
#' @return
#' \code{irregular.double.torus.surface.embedding()} returns an \code{n x 3}
#' numeric matrix with columns \code{x}, \code{y}, and \code{z}.
#'
#' \code{irregular.double.torus.surface.graph()} returns a list with
#' components:
#' \itemize{
#'   \item \code{edges}: the irregular double-torus edges,
#'   \item \code{n}: number of vertices,
#'   \item \code{edge_weights}: induced positive edge lengths,
#'   \item \code{coords_surface}: the 3D embedding,
#'   \item \code{coords_param}: the canonical 3D double-torus coordinates,
#'   \item \code{weight_scale}: the normalization constant applied to the raw
#'     edge lengths,
#'   \item \code{family}: always \code{"irregular.double.torus"},
#'   \item \code{surface}: the chosen surface name,
#'   \item \code{slices}: number of non-pole slices,
#'   \item \code{slice_sizes}: vertex counts in each slice,
#'   \item \code{slice_components}: number of cyclic components in each slice,
#'   \item \code{label}: a human-readable family label.
#' }
#'
#' @name irregular_double_torus_surface_helpers
NULL

#' @rdname irregular_double_torus_surface_helpers
#' @export
irregular.double.torus.surface.embedding <- function(
    slices = 11,
    tube_count = 14,
    branch_length = 0.85,
    branch_offset = 0.72,
    tube_radius = 0.28,
    transition_width = 0.42,
    count_irregularity = 0.2,
    axial_irregularity = 0.3,
    phase_twist = 0.35,
    surface = c("standard", "bulged", "twisted", "wavy"),
    amplitude = 0.25,
    freq_x = 2,
    freq_theta = 2,
    twist = 0.6) {
  surface <- match.arg(surface)
  built <- .irregular.double.torus.canonical(
    slices = slices,
    tube_count = tube_count,
    branch_length = branch_length,
    branch_offset = branch_offset,
    tube_radius = tube_radius,
    transition_width = transition_width,
    count_irregularity = count_irregularity,
    axial_irregularity = axial_irregularity,
    phase_twist = phase_twist
  )
  .irregular.double.torus.surface.coords(
    coords_param = built$coords,
    surface = surface,
    amplitude = amplitude,
    freq_x = freq_x,
    freq_theta = freq_theta,
    twist = twist
  )
}

#' @rdname irregular_double_torus_surface_helpers
#' @export
irregular.double.torus.surface.graph <- function(
    slices = 11,
    tube_count = 14,
    branch_length = 0.85,
    branch_offset = 0.72,
    tube_radius = 0.28,
    transition_width = 0.42,
    count_irregularity = 0.2,
    axial_irregularity = 0.3,
    phase_twist = 0.35,
    surface = c("standard", "bulged", "twisted", "wavy"),
    amplitude = 0.25,
    freq_x = 2,
    freq_theta = 2,
    twist = 0.6,
    normalize = c("median", "mean", "none")) {
  surface <- match.arg(surface)
  normalize <- match.arg(normalize)
  built <- .irregular.double.torus.canonical(
    slices = slices,
    tube_count = tube_count,
    branch_length = branch_length,
    branch_offset = branch_offset,
    tube_radius = tube_radius,
    transition_width = transition_width,
    count_irregularity = count_irregularity,
    axial_irregularity = axial_irregularity,
    phase_twist = phase_twist
  )
  coords_surface <- irregular.double.torus.surface.embedding(
    slices = slices,
    tube_count = tube_count,
    branch_length = branch_length,
    branch_offset = branch_offset,
    tube_radius = tube_radius,
    transition_width = transition_width,
    count_irregularity = count_irregularity,
    axial_irregularity = axial_irregularity,
    phase_twist = phase_twist,
    surface = surface,
    amplitude = amplitude,
    freq_x = freq_x,
    freq_theta = freq_theta,
    twist = twist
  )
  weights <- .edge.weights.from.embedding(
    edges = built$edges,
    coords = coords_surface,
    normalize = normalize
  )

  out <- list(
    edges = built$edges,
    n = built$n,
    edge_weights = weights$edge_weights,
    coords_surface = coords_surface,
    coords_param = built$coords,
    weight_scale = weights$weight_scale,
    family = "irregular.double.torus",
    surface = surface,
    slices = built$slices,
    slice_sizes = built$slice_sizes,
    slice_components = built$slice_components,
    normalize = normalize,
    label = sprintf("%s irregular double torus slices %d",
                    tools::toTitleCase(surface),
                    built$slices)
  )
  class(out) <- c("grip_irregular_double_torus_surface_graph", "list")
  out
}

.irregular.radial.solid.surface.coords <- function(coords_param,
                                                   surface,
                                                   amplitude,
                                                   freq_theta,
                                                   freq_phi,
                                                   twist) {
  surface <- .as_named_choice(surface,
                              c("standard", "bulged", "twisted", "wavy"),
                              "surface")
  amplitude <- .as_finite_scalar(amplitude, "amplitude")
  freq_theta <- .as_positive_scalar(freq_theta, "freq_theta")
  freq_phi <- .as_positive_scalar(freq_phi, "freq_phi")
  twist <- .as_finite_scalar(twist, "twist")

  coords_param <- as.matrix(coords_param)
  norms <- sqrt(rowSums(coords_param^2))
  max_norm <- max(norms)
  r_norm <- if (is.finite(max_norm) && max_norm > 0) norms / max_norm else rep(0, length(norms))
  theta <- atan2(coords_param[, 2L], coords_param[, 1L])
  phi <- ifelse(norms > 0, acos(pmin(pmax(coords_param[, 3L] / norms, -1), 1)), 0)
  unit <- coords_param
  keep <- norms > 0
  if (any(keep)) {
    unit[keep, ] <- sweep(coords_param[keep, , drop = FALSE], 1L, norms[keep], "/")
  }

  coords <- switch(
    surface,
    standard = coords_param,
    bulged = {
      scale <- 1 + amplitude * (1 - (2 * r_norm - 1)^2)
      if (any(!is.finite(scale) | scale <= 0)) {
        stop("radial solid bulged parameters require positive radial scale; adjust amplitude",
             call. = FALSE)
      }
      unit * (norms * scale)
    },
    twisted = {
      angle <- twist * (2 * r_norm - 1)
      x <- coords_param[, 1L] * cos(angle) - coords_param[, 2L] * sin(angle)
      y <- coords_param[, 1L] * sin(angle) + coords_param[, 2L] * cos(angle)
      cbind(x = x, y = y, z = coords_param[, 3L])
    },
    wavy = {
      scale <- 1 + amplitude * (0.35 + 0.65 * r_norm) *
        sin(freq_theta * theta) * cos(freq_phi * (phi - pi / 2))
      if (any(!is.finite(scale) | scale <= 0)) {
        stop("radial solid wavy parameters require positive radial scale; adjust amplitude or frequencies",
             call. = FALSE)
      }
      unit * (norms * scale)
    }
  )

  storage.mode(coords) <- "double"
  coords
}

#' Weighted irregular ball solid helpers
#'
#' Convenience helpers that build a layered simplicial ball by taking a
#' subdivided triangulated polyhedron, normalizing its vertices to directions on
#' the sphere, placing those directions on nested radial layers, and connecting
#' adjacent layers by a deterministic prism-to-tetrahedra edge pattern. The
#' center vertex is included explicitly, so the resulting graph is genuinely
#' volumetric rather than just a closed surface.
#'
#' @param base Base polyhedron. One of \code{"tetrahedron"},
#'   \code{"octahedron"}, or \code{"icosahedron"}.
#' @param level Surface subdivision depth used for each radial layer.
#' @param layers Number of non-center radial layers.
#' @param outer_radius Positive outer radius of the ball.
#' @param radial_irregularity Irregularity level for the radial layer spacing.
#'   Must lie in \code{[0, 1]}.
#' @param layer_twist Finite z-axis twist applied smoothly across radial
#'   layers.
#' @param surface Solid geometry family. One of \code{"standard"},
#'   \code{"bulged"}, \code{"twisted"}, or \code{"wavy"}.
#' @param amplitude Finite deformation amplitude.
#' @param freq_theta Positive azimuthal modulation frequency used only when
#'   \code{surface = "wavy"}.
#' @param freq_phi Positive polar modulation frequency used only when
#'   \code{surface = "wavy"}.
#' @param twist Finite twist strength used only when
#'   \code{surface = "twisted"}.
#' @param normalize Normalization applied to the induced edge lengths. One of
#'   \code{"median"}, \code{"mean"}, or \code{"none"}.
#'
#' @return
#' \code{irregular.ball.solid.embedding()} returns an \code{n x 3} numeric
#' matrix with columns \code{x}, \code{y}, and \code{z}.
#'
#' \code{irregular.ball.solid.graph()} returns a list with components:
#' \itemize{
#'   \item \code{edges}: the irregular ball edges,
#'   \item \code{n}: number of vertices,
#'   \item \code{edge_weights}: induced positive edge lengths,
#'   \item \code{coords_surface}: the 3D embedding,
#'   \item \code{coords_param}: the canonical 3D ball coordinates,
#'   \item \code{weight_scale}: the normalization constant applied to the raw
#'     edge lengths,
#'   \item \code{family}: always \code{"irregular.ball"},
#'   \item \code{surface}: the chosen solid geometry family,
#'   \item \code{base}: the chosen base polyhedron,
#'   \item \code{level}: the surface subdivision depth,
#'   \item \code{layers}: number of non-center radial layers,
#'   \item \code{radii}: the radial layer values,
#'   \item \code{label}: a human-readable family label.
#' }
#'
#' @name irregular_ball_solid_helpers
NULL

#' @rdname irregular_ball_solid_helpers
#' @export
irregular.ball.solid.embedding <- function(
    base = c("tetrahedron", "octahedron", "icosahedron"),
    level = 1,
    layers = 3,
    outer_radius = 1,
    radial_irregularity = 0.25,
    layer_twist = 0.35,
    surface = c("standard", "bulged", "twisted", "wavy"),
    amplitude = 0.2,
    freq_theta = 2,
    freq_phi = 2,
    twist = 0.6) {
  base <- match.arg(base)
  surface <- match.arg(surface)
  built <- .irregular.ball.canonical(
    base = base,
    level = level,
    layers = layers,
    outer_radius = outer_radius,
    radial_irregularity = radial_irregularity,
    layer_twist = layer_twist
  )
  .irregular.radial.solid.surface.coords(
    coords_param = built$coords,
    surface = surface,
    amplitude = amplitude,
    freq_theta = freq_theta,
    freq_phi = freq_phi,
    twist = twist
  )
}

#' @rdname irregular_ball_solid_helpers
#' @export
irregular.ball.solid.graph <- function(
    base = c("tetrahedron", "octahedron", "icosahedron"),
    level = 1,
    layers = 3,
    outer_radius = 1,
    radial_irregularity = 0.25,
    layer_twist = 0.35,
    surface = c("standard", "bulged", "twisted", "wavy"),
    amplitude = 0.2,
    freq_theta = 2,
    freq_phi = 2,
    twist = 0.6,
    normalize = c("median", "mean", "none")) {
  base <- match.arg(base)
  surface <- match.arg(surface)
  normalize <- match.arg(normalize)
  built <- .irregular.ball.canonical(
    base = base,
    level = level,
    layers = layers,
    outer_radius = outer_radius,
    radial_irregularity = radial_irregularity,
    layer_twist = layer_twist
  )
  coords_surface <- irregular.ball.solid.embedding(
    base = base,
    level = level,
    layers = layers,
    outer_radius = outer_radius,
    radial_irregularity = radial_irregularity,
    layer_twist = layer_twist,
    surface = surface,
    amplitude = amplitude,
    freq_theta = freq_theta,
    freq_phi = freq_phi,
    twist = twist
  )
  weights <- .edge.weights.from.embedding(
    edges = built$edges,
    coords = coords_surface,
    normalize = normalize
  )

  out <- list(
    edges = built$edges,
    n = built$n,
    edge_weights = weights$edge_weights,
    coords_surface = coords_surface,
    coords_param = built$coords,
    weight_scale = weights$weight_scale,
    family = "irregular.ball",
    surface = surface,
    base = built$base,
    level = built$level,
    layers = built$layers,
    radii = built$radii,
    subdivision = built$subdivision,
    normalize = normalize,
    label = sprintf("%s irregular ball %s level %d layers %d",
                    tools::toTitleCase(surface),
                    built$base,
                    built$level,
                    built$layers)
  )
  class(out) <- c("grip_irregular_ball_solid_graph", "list")
  out
}

#' Weighted irregular shell solid helpers
#'
#' Convenience helpers that build a layered simplicial shell by placing the
#' vertices of a subdivided triangulated polyhedron on nested inner-to-outer
#' radial layers and connecting adjacent layers by the same deterministic
#' prism-to-tetrahedra edge pattern used for \code{irregular.ball.solid.*()}.
#' The result is a genuinely volumetric shell graph with a hollow interior.
#'
#' @param inner_radius Positive inner radius of the shell.
#' @inheritParams irregular_ball_solid_helpers
#'
#' @return
#' \code{irregular.shell.solid.embedding()} returns an \code{n x 3} numeric
#' matrix with columns \code{x}, \code{y}, and \code{z}.
#'
#' \code{irregular.shell.solid.graph()} returns the same components as
#' \code{irregular.ball.solid.graph()}, with \code{family} set to
#' \code{"irregular.shell"}.
#'
#' @name irregular_shell_solid_helpers
NULL

#' @rdname irregular_shell_solid_helpers
#' @export
irregular.shell.solid.embedding <- function(
    base = c("tetrahedron", "octahedron", "icosahedron"),
    level = 1,
    layers = 3,
    inner_radius = 0.45,
    outer_radius = 1,
    radial_irregularity = 0.25,
    layer_twist = 0.35,
    surface = c("standard", "bulged", "twisted", "wavy"),
    amplitude = 0.2,
    freq_theta = 2,
    freq_phi = 2,
    twist = 0.6) {
  base <- match.arg(base)
  surface <- match.arg(surface)
  built <- .irregular.shell.canonical(
    base = base,
    level = level,
    layers = layers,
    inner_radius = inner_radius,
    outer_radius = outer_radius,
    radial_irregularity = radial_irregularity,
    layer_twist = layer_twist
  )
  .irregular.radial.solid.surface.coords(
    coords_param = built$coords,
    surface = surface,
    amplitude = amplitude,
    freq_theta = freq_theta,
    freq_phi = freq_phi,
    twist = twist
  )
}

#' @rdname irregular_shell_solid_helpers
#' @export
irregular.shell.solid.graph <- function(
    base = c("tetrahedron", "octahedron", "icosahedron"),
    level = 1,
    layers = 3,
    inner_radius = 0.45,
    outer_radius = 1,
    radial_irregularity = 0.25,
    layer_twist = 0.35,
    surface = c("standard", "bulged", "twisted", "wavy"),
    amplitude = 0.2,
    freq_theta = 2,
    freq_phi = 2,
    twist = 0.6,
    normalize = c("median", "mean", "none")) {
  base <- match.arg(base)
  surface <- match.arg(surface)
  normalize <- match.arg(normalize)
  built <- .irregular.shell.canonical(
    base = base,
    level = level,
    layers = layers,
    inner_radius = inner_radius,
    outer_radius = outer_radius,
    radial_irregularity = radial_irregularity,
    layer_twist = layer_twist
  )
  coords_surface <- irregular.shell.solid.embedding(
    base = base,
    level = level,
    layers = layers,
    inner_radius = inner_radius,
    outer_radius = outer_radius,
    radial_irregularity = radial_irregularity,
    layer_twist = layer_twist,
    surface = surface,
    amplitude = amplitude,
    freq_theta = freq_theta,
    freq_phi = freq_phi,
    twist = twist
  )
  weights <- .edge.weights.from.embedding(
    edges = built$edges,
    coords = coords_surface,
    normalize = normalize
  )

  out <- list(
    edges = built$edges,
    n = built$n,
    edge_weights = weights$edge_weights,
    coords_surface = coords_surface,
    coords_param = built$coords,
    weight_scale = weights$weight_scale,
    family = "irregular.shell",
    surface = surface,
    base = built$base,
    level = built$level,
    layers = built$layers,
    radii = built$radii,
    subdivision = built$subdivision,
    normalize = normalize,
    label = sprintf("%s irregular shell %s level %d layers %d",
                    tools::toTitleCase(surface),
                    built$base,
                    built$level,
                    built$layers)
  )
  class(out) <- c("grip_irregular_shell_solid_graph", "list")
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

#' Weighted irregular sphere surface helpers
#'
#' Convenience helpers that build a sphere from deterministically irregular
#' latitude bands with varying sample counts, stitch adjacent bands into a
#' locally triangulated graph, and then realize that graph in \eqn{\mathbb{R}^3}
#' using either a standard sphere, an ellipsoid, or a wavy radial modulation.
#'
#' @param bands Number of non-pole latitude bands.
#' @param equator_count Approximate number of vertices near the equator.
#' @param count_irregularity Irregularity level for per-band sample counts.
#'   Must lie in \code{[0, 1)}.
#' @param lat_irregularity Irregularity level for latitude-band spacing. Must
#'   lie in \code{[0, 1]}.
#' @param phase_twist Finite angular phase offset used to desynchronize
#'   neighboring latitude bands.
#' @param surface Sphere geometry family. One of \code{"standard"},
#'   \code{"ellipsoid"}, or \code{"wavy"}.
#' @param radius Positive baseline radius.
#' @param amplitude Finite deformation amplitude. For \code{"ellipsoid"},
#'   positive values make the shape oblate and negative values make it prolate.
#' @param freq_theta Positive longitudinal frequency used only when
#'   \code{surface = "wavy"}.
#' @param freq_lat Positive latitudinal frequency used only when
#'   \code{surface = "wavy"}.
#' @param twist Finite longitude twist applied smoothly by latitude.
#' @param normalize Normalization applied to the induced edge lengths. One of
#'   \code{"median"}, \code{"mean"}, or \code{"none"}.
#'
#' @return
#' \code{irregular.sphere.surface.embedding()} returns an \code{n x 3} numeric
#' matrix with columns \code{x}, \code{y}, and \code{z}.
#'
#' \code{irregular.sphere.surface.graph()} returns a list with components:
#' \itemize{
#'   \item \code{edges}: the irregular sphere edges,
#'   \item \code{n}: number of vertices,
#'   \item \code{edge_weights}: induced positive edge lengths,
#'   \item \code{coords_surface}: the 3D embedding,
#'   \item \code{coords_param}: the irregular angular parameter coordinates,
#'   \item \code{weight_scale}: the normalization constant applied to the raw
#'     edge lengths,
#'   \item \code{family}: always \code{"irregular.sphere"},
#'   \item \code{surface}: the chosen surface name,
#'   \item \code{bands}: number of non-pole latitude bands,
#'   \item \code{band_sizes}: sample counts in each latitude band,
#'   \item \code{label}: a human-readable family label.
#' }
#'
#' @name irregular_sphere_surface_helpers
NULL

#' @rdname irregular_sphere_surface_helpers
#' @export
irregular.sphere.surface.embedding <- function(
    bands = 6,
    equator_count = 28,
    count_irregularity = 0.2,
    lat_irregularity = 0.35,
    phase_twist = 0.35,
    surface = c("standard", "ellipsoid", "wavy"),
    radius = 1,
    amplitude = 0.2,
    freq_theta = 3,
    freq_lat = 2,
    twist = 0.25) {
  surface <- match.arg(surface)
  radius <- .as_positive_scalar(radius, "radius")
  amplitude <- .as_finite_scalar(amplitude, "amplitude")
  freq_theta <- .as_positive_scalar(freq_theta, "freq_theta")
  freq_lat <- .as_positive_scalar(freq_lat, "freq_lat")
  twist <- .as_finite_scalar(twist, "twist")

  built <- .irregular.sphere.canonical(
    bands = bands,
    equator_count = equator_count,
    count_irregularity = count_irregularity,
    lat_irregularity = lat_irregularity,
    phase_twist = phase_twist
  )
  theta <- built$coords_param[, 1L]
  lat <- built$coords_param[, 2L]
  theta_eff <- theta + twist * sin(2 * lat)

  if (surface == "ellipsoid") {
    equatorial_radius <- radius * (1 + amplitude)
    polar_radius <- radius * (1 - amplitude)
    if (!is.finite(equatorial_radius) || !is.finite(polar_radius) ||
        equatorial_radius <= 0 || polar_radius <= 0) {
      stop("irregular sphere ellipsoid parameters require positive equatorial and polar radii; adjust amplitude",
           call. = FALSE)
    }
    coords <- cbind(
      x = equatorial_radius * cos(lat) * cos(theta_eff),
      y = equatorial_radius * cos(lat) * sin(theta_eff),
      z = polar_radius * sin(lat)
    )
  } else {
    pole_taper <- cos(lat)^2
    local_radius <- switch(
      surface,
      standard = rep(radius, length(lat)),
      wavy = radius * (1 + amplitude * pole_taper *
                         sin(freq_theta * theta) * cos(freq_lat * lat))
    )
    if (any(!is.finite(local_radius) | local_radius <= 0)) {
      stop("irregular sphere parameters produce a non-positive local radius; adjust amplitude or frequencies",
           call. = FALSE)
    }
    coords <- cbind(
      x = local_radius * cos(lat) * cos(theta_eff),
      y = local_radius * cos(lat) * sin(theta_eff),
      z = local_radius * sin(lat)
    )
  }

  colnames(coords) <- c("x", "y", "z")
  storage.mode(coords) <- "double"
  coords
}

#' @rdname irregular_sphere_surface_helpers
#' @export
irregular.sphere.surface.graph <- function(
    bands = 6,
    equator_count = 28,
    count_irregularity = 0.2,
    lat_irregularity = 0.35,
    phase_twist = 0.35,
    surface = c("standard", "ellipsoid", "wavy"),
    radius = 1,
    amplitude = 0.2,
    freq_theta = 3,
    freq_lat = 2,
    twist = 0.25,
    normalize = c("median", "mean", "none")) {
  surface <- match.arg(surface)
  normalize <- match.arg(normalize)
  built <- .irregular.sphere.canonical(
    bands = bands,
    equator_count = equator_count,
    count_irregularity = count_irregularity,
    lat_irregularity = lat_irregularity,
    phase_twist = phase_twist
  )
  coords_surface <- irregular.sphere.surface.embedding(
    bands = bands,
    equator_count = equator_count,
    count_irregularity = count_irregularity,
    lat_irregularity = lat_irregularity,
    phase_twist = phase_twist,
    surface = surface,
    radius = radius,
    amplitude = amplitude,
    freq_theta = freq_theta,
    freq_lat = freq_lat,
    twist = twist
  )
  weights <- .edge.weights.from.embedding(
    edges = built$edges,
    coords = coords_surface,
    normalize = normalize
  )

  out <- list(
    edges = built$edges,
    n = built$n,
    edge_weights = weights$edge_weights,
    coords_surface = coords_surface,
    coords_param = built$coords_param,
    weight_scale = weights$weight_scale,
    family = "irregular.sphere",
    surface = surface,
    bands = built$bands,
    band_sizes = built$band_sizes,
    normalize = normalize,
    label = sprintf("%s irregular sphere bands %d",
                    tools::toTitleCase(surface),
                    built$bands)
  )
  class(out) <- c("grip_irregular_sphere_surface_graph", "list")
  out
}

.triangulated.polyhedron.surface.coords <- function(coords_param,
                                                    surface,
                                                    amplitude,
                                                    freq,
                                                    twist) {
  surface <- .as_named_choice(surface,
                              c("standard", "inflated", "twisted", "wavy"),
                              "surface")
  amplitude <- .as_finite_scalar(amplitude, "amplitude")
  freq <- .as_positive_scalar(freq, "freq")
  twist <- .as_finite_scalar(twist, "twist")

  coords_param <- as.matrix(coords_param)
  norms <- sqrt(rowSums(coords_param^2))
  unit <- .normalize.row.coords(coords_param)
  coords_norm <- .normalize.center_radius.coords(coords_param)
  coords <- switch(
    surface,
    standard = coords_param,
    inflated = {
      target_radius <- norms + amplitude * (1 - norms)
      if (any(!is.finite(target_radius) | target_radius <= 0)) {
        stop("inflated polyhedron parameters produce a non-positive radius; adjust amplitude",
             call. = FALSE)
      }
      unit * target_radius
    },
    twisted = {
      theta <- twist * coords_norm[, 3L]
      x <- coords_param[, 1L] * cos(theta) - coords_param[, 2L] * sin(theta)
      y <- coords_param[, 1L] * sin(theta) + coords_param[, 2L] * cos(theta)
      cbind(x = x, y = y, z = coords_param[, 3L])
    },
    wavy = {
      radial_scale <- 1 + amplitude *
        sin(pi * freq * coords_norm[, 1L]) *
        cos(pi * freq * coords_norm[, 2L]) *
        sin(pi * freq * coords_norm[, 3L])
      if (any(!is.finite(radial_scale) | radial_scale <= 0)) {
        stop("wavy polyhedron parameters produce a non-positive local radius; adjust amplitude or frequency",
             call. = FALSE)
      }
      unit * (norms * radial_scale)
    }
  )

  storage.mode(coords) <- "double"
  coords
}

#' Weighted triangulated polyhedron surface helpers
#'
#' Convenience helpers that take a triangular-face polyhedron, repeatedly split
#' each face into four subtriangles, merge shared-edge vertices, and use the
#' resulting triangulated-surface graph as a reusable benchmark family. The
#' canonical geometry is a piecewise-flat triangulated surface, while the
#' weighted variants can inflate, twist, or radially modulate that surface in
#' \eqn{\mathbb{R}^3}.
#'
#' The current supported bases are the regular tetrahedron, octahedron, and
#' icosahedron. These give closed triangulated surfaces with non-grid topology
#' and a small number of extraordinary vertices, making them useful complements
#' to the mesh-derived families already in the package.
#'
#' @param base Base polyhedron. One of \code{"tetrahedron"},
#'   \code{"octahedron"}, or \code{"icosahedron"}.
#' @param level Subdivision depth. \code{level = 0} returns the base
#'   triangulation, \code{level = 1} splits each face into four, and so on.
#' @param surface Geometry family used for the 3D embedding. One of
#'   \code{"standard"}, \code{"inflated"}, \code{"twisted"}, or
#'   \code{"wavy"}.
#' @param amplitude Finite deformation amplitude.
#' @param freq Positive modulation frequency used only when
#'   \code{surface = "wavy"}.
#' @param twist Finite twist strength used only when
#'   \code{surface = "twisted"}.
#' @param normalize Normalization applied to the induced edge lengths. One of
#'   \code{"median"}, \code{"mean"}, or \code{"none"}.
#'
#' @return
#' \code{triangulated.polyhedron.surface.embedding()} returns an \code{n x 3}
#' numeric matrix with columns \code{x}, \code{y}, and \code{z}.
#'
#' \code{triangulated.polyhedron.surface.graph()} returns a list with
#' components:
#' \itemize{
#'   \item \code{edges}: the triangulated-surface edges,
#'   \item \code{n}: number of vertices,
#'   \item \code{edge_weights}: induced positive edge lengths,
#'   \item \code{coords_surface}: the 3D embedding,
#'   \item \code{coords_param}: the canonical piecewise-flat coordinates,
#'   \item \code{weight_scale}: the normalization constant applied to the raw
#'     edge lengths,
#'   \item \code{family}: always \code{"triangulated.polyhedron"},
#'   \item \code{base}: the chosen base polyhedron,
#'   \item \code{surface}: the chosen surface name,
#'   \item \code{level}: the subdivision depth,
#'   \item \code{subdivision}: the linear face subdivision factor
#'     \code{2^level},
#'   \item \code{label}: a human-readable family label.
#' }
#'
#' @name triangulated_polyhedron_surface_helpers
NULL

#' @rdname triangulated_polyhedron_surface_helpers
#' @export
triangulated.polyhedron.surface.embedding <- function(
    base = c("tetrahedron", "octahedron", "icosahedron"),
    level = 1,
    surface = c("standard", "inflated", "twisted", "wavy"),
    amplitude = 0.25,
    freq = 2,
    twist = 0.6) {
  base <- match.arg(base)
  level <- .as_whole_number(level, "level")
  surface <- match.arg(surface)

  built <- .triangulated.polyhedron.canonical(base = base, level = level)
  .triangulated.polyhedron.surface.coords(
    coords_param = built$coords,
    surface = surface,
    amplitude = amplitude,
    freq = freq,
    twist = twist
  )
}

#' @rdname triangulated_polyhedron_surface_helpers
#' @export
triangulated.polyhedron.surface.graph <- function(
    base = c("tetrahedron", "octahedron", "icosahedron"),
    level = 1,
    surface = c("standard", "inflated", "twisted", "wavy"),
    amplitude = 0.25,
    freq = 2,
    twist = 0.6,
    normalize = c("median", "mean", "none")) {
  base <- match.arg(base)
  level <- .as_whole_number(level, "level")
  surface <- match.arg(surface)
  normalize <- match.arg(normalize)

  built <- .triangulated.polyhedron.canonical(base = base, level = level)
  coords_surface <- triangulated.polyhedron.surface.embedding(
    base = base,
    level = level,
    surface = surface,
    amplitude = amplitude,
    freq = freq,
    twist = twist
  )
  weights <- .edge.weights.from.embedding(
    edges = built$edges,
    coords = coords_surface,
    normalize = normalize
  )

  out <- list(
    edges = built$edges,
    n = built$n,
    edge_weights = weights$edge_weights,
    coords_surface = coords_surface,
    coords_param = built$coords,
    weight_scale = weights$weight_scale,
    family = "triangulated.polyhedron",
    base = base,
    surface = surface,
    level = level,
    subdivision = built$subdivision,
    normalize = normalize,
    label = sprintf("%s triangulated %s level %d",
                    tools::toTitleCase(surface),
                    base,
                    level)
  )
  class(out) <- c("grip_triangulated_polyhedron_surface_graph", "list")
  out
}

.triangulated.planar.surface.coords <- function(coords_param,
                                                surface,
                                                amplitude,
                                                freq_u,
                                                freq_v) {
  surface <- .as_named_choice(surface,
                              c("flat", "saddle", "paraboloid", "ripple", "folded"),
                              "surface")
  amplitude <- .as_finite_scalar(amplitude, "amplitude")
  freq_u <- .as_positive_scalar(freq_u, "freq_u")
  freq_v <- .as_positive_scalar(freq_v, "freq_v")

  coords_param <- as.matrix(coords_param)
  coords_norm <- .normalize.center_radius.coords(coords_param)
  z <- switch(
    surface,
    flat = rep(0, nrow(coords_param)),
    saddle = .surface.z.from_uv(coords_norm[, 1L], coords_norm[, 2L],
                                "saddle", amplitude, freq_u, freq_v),
    paraboloid = .surface.z.from_uv(coords_norm[, 1L], coords_norm[, 2L],
                                    "paraboloid", amplitude, freq_u, freq_v),
    ripple = .surface.z.from_uv(coords_norm[, 1L], coords_norm[, 2L],
                                "ripple", amplitude, freq_u, freq_v),
    folded = amplitude * abs(coords_norm[, 1L])
  )
  coords <- cbind(x = coords_param[, 1L], y = coords_param[, 2L], z = z)
  storage.mode(coords) <- "double"
  coords
}

#' Weighted triangulated annulus surface helpers
#'
#' Convenience helpers that triangulate a planar annulus by clipping a regular
#' triangular lattice to the region between two concentric circles. The
#' resulting graph is a triangulated manifold with boundary and a single hole,
#' making it a useful complement to the closed triangulated-polyhedron
#' families.
#'
#' \code{triangulated.annulus.surface.embedding()} returns the 3D coordinates
#' of the clipped lattice vertices in the same vertex order as
#' \code{edges.triangulated.annulus()}. \code{triangulated.annulus.surface.graph()}
#' returns a reusable weighted-graph bundle with the annulus edges, induced edge
#' weights, the 3D embedding, and the 2D annulus parameter coordinates.
#'
#' @param resolution Positive lattice-resolution control. Larger values produce
#'   finer triangulations.
#' @param outer_radius Positive outer annulus radius.
#' @param inner_radius Positive inner annulus radius. Must be strictly smaller
#'   than \code{outer_radius}.
#' @param surface Geometry family used for the 3D lift. One of \code{"flat"},
#'   \code{"saddle"}, \code{"paraboloid"}, \code{"ripple"}, or
#'   \code{"folded"}.
#' @param amplitude Finite deformation amplitude.
#' @param freq_u Positive ripple frequency in the first planar coordinate. Used
#'   only when \code{surface = "ripple"}.
#' @param freq_v Positive ripple frequency in the second planar coordinate. Used
#'   only when \code{surface = "ripple"}.
#' @param normalize Normalization applied to the induced edge lengths. One of
#'   \code{"median"}, \code{"mean"}, or \code{"none"}.
#'
#' @return
#' \code{triangulated.annulus.surface.embedding()} returns an \code{n x 3}
#' numeric matrix with columns \code{x}, \code{y}, and \code{z}.
#'
#' \code{triangulated.annulus.surface.graph()} returns a list with components:
#' \itemize{
#'   \item \code{edges}: the triangulated-annulus edges,
#'   \item \code{n}: number of vertices,
#'   \item \code{edge_weights}: induced positive edge lengths,
#'   \item \code{coords_surface}: the 3D embedding,
#'   \item \code{coords_param}: the canonical 2D annulus coordinates,
#'   \item \code{weight_scale}: the normalization constant applied to the raw
#'     edge lengths,
#'   \item \code{family}: always \code{"triangulated.annulus"},
#'   \item \code{surface}: the chosen surface name,
#'   \item \code{resolution}: the lattice-resolution parameter,
#'   \item \code{label}: a human-readable family label.
#' }
#'
#' @name triangulated_annulus_surface_helpers
NULL

#' @rdname triangulated_annulus_surface_helpers
#' @export
triangulated.annulus.surface.embedding <- function(
    resolution = 12,
    outer_radius = 1,
    inner_radius = 0.45,
    surface = c("flat", "saddle", "paraboloid", "ripple", "folded"),
    amplitude = 0.6,
    freq_u = 1,
    freq_v = 1) {
  surface <- match.arg(surface)
  built <- .triangulated.annulus.canonical(
    resolution = resolution,
    outer_radius = outer_radius,
    inner_radius = inner_radius
  )
  .triangulated.planar.surface.coords(
    coords_param = built$coords,
    surface = surface,
    amplitude = amplitude,
    freq_u = freq_u,
    freq_v = freq_v
  )
}

#' @rdname triangulated_annulus_surface_helpers
#' @export
triangulated.annulus.surface.graph <- function(
    resolution = 12,
    outer_radius = 1,
    inner_radius = 0.45,
    surface = c("flat", "saddle", "paraboloid", "ripple", "folded"),
    amplitude = 0.6,
    freq_u = 1,
    freq_v = 1,
    normalize = c("median", "mean", "none")) {
  surface <- match.arg(surface)
  normalize <- match.arg(normalize)
  built <- .triangulated.annulus.canonical(
    resolution = resolution,
    outer_radius = outer_radius,
    inner_radius = inner_radius
  )
  coords_surface <- triangulated.annulus.surface.embedding(
    resolution = resolution,
    outer_radius = outer_radius,
    inner_radius = inner_radius,
    surface = surface,
    amplitude = amplitude,
    freq_u = freq_u,
    freq_v = freq_v
  )
  weights <- .edge.weights.from.embedding(
    edges = built$edges,
    coords = coords_surface,
    normalize = normalize
  )

  out <- list(
    edges = built$edges,
    n = built$n,
    edge_weights = weights$edge_weights,
    coords_surface = coords_surface,
    coords_param = built$coords,
    weight_scale = weights$weight_scale,
    family = "triangulated.annulus",
    surface = surface,
    resolution = built$resolution,
    normalize = normalize,
    label = sprintf("%s triangulated annulus res %d",
                    tools::toTitleCase(surface),
                    built$resolution)
  )
  class(out) <- c("grip_triangulated_annulus_surface_graph", "list")
  out
}

#' Weighted triangulated pair-of-pants surface helpers
#'
#' Convenience helpers that triangulate a planar pair-of-pants domain by
#' clipping a regular triangular lattice to a disk with two interior circular
#' holes. The resulting graph is a triangulated surface with three boundary
#' components and a deterministic non-grid topology.
#'
#' \code{triangulated.pair.of.pants.surface.embedding()} returns the 3D
#' coordinates of the clipped lattice vertices in the same vertex order as
#' \code{edges.triangulated.pair.of.pants()}.
#' \code{triangulated.pair.of.pants.surface.graph()} returns a reusable
#' weighted-graph bundle with the pair-of-pants edges, induced edge weights, 3D
#' embedding, and the 2D parameter coordinates of the clipped domain.
#'
#' @param resolution Positive lattice-resolution control. Larger values produce
#'   finer triangulations.
#' @param outer_radius Positive radius of the outer boundary.
#' @param hole_radius Positive radius of each interior hole.
#' @param hole_offset Positive horizontal offset of the two hole centers from
#'   the vertical axis.
#' @param hole_height Vertical coordinate shared by the two hole centers.
#' @param surface Geometry family used for the 3D lift. One of \code{"flat"},
#'   \code{"saddle"}, \code{"paraboloid"}, \code{"ripple"}, or
#'   \code{"folded"}.
#' @param amplitude Finite deformation amplitude.
#' @param freq_u Positive ripple frequency in the first planar coordinate. Used
#'   only when \code{surface = "ripple"}.
#' @param freq_v Positive ripple frequency in the second planar coordinate. Used
#'   only when \code{surface = "ripple"}.
#' @param normalize Normalization applied to the induced edge lengths. One of
#'   \code{"median"}, \code{"mean"}, or \code{"none"}.
#'
#' @return
#' \code{triangulated.pair.of.pants.surface.embedding()} returns an \code{n x
#' 3} numeric matrix with columns \code{x}, \code{y}, and \code{z}.
#'
#' \code{triangulated.pair.of.pants.surface.graph()} returns a list with
#' components:
#' \itemize{
#'   \item \code{edges}: the triangulated pair-of-pants edges,
#'   \item \code{n}: number of vertices,
#'   \item \code{edge_weights}: induced positive edge lengths,
#'   \item \code{coords_surface}: the 3D embedding,
#'   \item \code{coords_param}: the canonical 2D pair-of-pants coordinates,
#'   \item \code{weight_scale}: the normalization constant applied to the raw
#'     edge lengths,
#'   \item \code{family}: always \code{"triangulated.pair.of.pants"},
#'   \item \code{surface}: the chosen surface name,
#'   \item \code{resolution}: the lattice-resolution parameter,
#'   \item \code{label}: a human-readable family label.
#' }
#'
#' @name triangulated_pair_of_pants_surface_helpers
NULL

#' @rdname triangulated_pair_of_pants_surface_helpers
#' @export
triangulated.pair.of.pants.surface.embedding <- function(
    resolution = 12,
    outer_radius = 1.1,
    hole_radius = 0.24,
    hole_offset = 0.38,
    hole_height = 0.18,
    surface = c("flat", "saddle", "paraboloid", "ripple", "folded"),
    amplitude = 0.6,
    freq_u = 1,
    freq_v = 1) {
  surface <- match.arg(surface)
  built <- .triangulated.pair.of.pants.canonical(
    resolution = resolution,
    outer_radius = outer_radius,
    hole_radius = hole_radius,
    hole_offset = hole_offset,
    hole_height = hole_height
  )
  .triangulated.planar.surface.coords(
    coords_param = built$coords,
    surface = surface,
    amplitude = amplitude,
    freq_u = freq_u,
    freq_v = freq_v
  )
}

#' @rdname triangulated_pair_of_pants_surface_helpers
#' @export
triangulated.pair.of.pants.surface.graph <- function(
    resolution = 12,
    outer_radius = 1.1,
    hole_radius = 0.24,
    hole_offset = 0.38,
    hole_height = 0.18,
    surface = c("flat", "saddle", "paraboloid", "ripple", "folded"),
    amplitude = 0.6,
    freq_u = 1,
    freq_v = 1,
    normalize = c("median", "mean", "none")) {
  surface <- match.arg(surface)
  normalize <- match.arg(normalize)
  built <- .triangulated.pair.of.pants.canonical(
    resolution = resolution,
    outer_radius = outer_radius,
    hole_radius = hole_radius,
    hole_offset = hole_offset,
    hole_height = hole_height
  )
  coords_surface <- triangulated.pair.of.pants.surface.embedding(
    resolution = resolution,
    outer_radius = outer_radius,
    hole_radius = hole_radius,
    hole_offset = hole_offset,
    hole_height = hole_height,
    surface = surface,
    amplitude = amplitude,
    freq_u = freq_u,
    freq_v = freq_v
  )
  weights <- .edge.weights.from.embedding(
    edges = built$edges,
    coords = coords_surface,
    normalize = normalize
  )

  out <- list(
    edges = built$edges,
    n = built$n,
    edge_weights = weights$edge_weights,
    coords_surface = coords_surface,
    coords_param = built$coords,
    weight_scale = weights$weight_scale,
    family = "triangulated.pair.of.pants",
    surface = surface,
    resolution = built$resolution,
    normalize = normalize,
    label = sprintf("%s triangulated pair-of-pants res %d",
                    tools::toTitleCase(surface),
                    built$resolution)
  )
  class(out) <- c("grip_triangulated_pair_of_pants_surface_graph", "list")
  out
}

#' Weighted irregular annulus surface helpers
#'
#' Convenience helpers that build an annulus from deterministically irregular
#' concentric rings with varying sample counts, stitch adjacent rings into a
#' locally triangulated graph, and then optionally lift the resulting irregular
#' planar parameterization into \eqn{\mathbb{R}^3}. This provides a non-lattice
#' annulus family with boundary and nonuniform local valence.
#'
#' @param rings Number of concentric sample rings, including the inner and
#'   outer boundary cycles.
#' @param outer_count Approximate number of vertices on the outer boundary.
#' @param outer_radius Positive outer annulus radius.
#' @param inner_radius Positive inner annulus radius.
#' @param count_irregularity Irregularity level for the per-ring sample counts.
#'   Must lie in \code{[0, 1)}.
#' @param radial_irregularity Irregularity level for within-ring radial
#'   perturbations. Must lie in \code{[0, 1]}.
#' @param phase_twist Finite angular phase offset used to desynchronize
#'   neighboring rings.
#' @param surface Geometry family used for the 3D lift. One of \code{"flat"},
#'   \code{"saddle"}, \code{"paraboloid"}, \code{"ripple"}, or
#'   \code{"folded"}.
#' @param amplitude Finite deformation amplitude.
#' @param freq_u Positive ripple frequency in the first planar coordinate. Used
#'   only when \code{surface = "ripple"}.
#' @param freq_v Positive ripple frequency in the second planar coordinate. Used
#'   only when \code{surface = "ripple"}.
#' @param normalize Normalization applied to the induced edge lengths. One of
#'   \code{"median"}, \code{"mean"}, or \code{"none"}.
#'
#' @return
#' \code{irregular.annulus.surface.embedding()} returns an \code{n x 3} numeric
#' matrix with columns \code{x}, \code{y}, and \code{z}.
#'
#' \code{irregular.annulus.surface.graph()} returns a list with components:
#' \itemize{
#'   \item \code{edges}: the irregular annulus edges,
#'   \item \code{n}: number of vertices,
#'   \item \code{edge_weights}: induced positive edge lengths,
#'   \item \code{coords_surface}: the 3D embedding,
#'   \item \code{coords_param}: the irregular planar annulus coordinates,
#'   \item \code{weight_scale}: the normalization constant applied to the raw
#'     edge lengths,
#'   \item \code{family}: always \code{"irregular.annulus"},
#'   \item \code{surface}: the chosen surface name,
#'   \item \code{rings}: number of concentric rings,
#'   \item \code{ring_sizes}: sample counts on each ring,
#'   \item \code{label}: a human-readable family label.
#' }
#'
#' @name irregular_annulus_surface_helpers
NULL

#' @rdname irregular_annulus_surface_helpers
#' @export
irregular.annulus.surface.embedding <- function(
    rings = 6,
    outer_count = 28,
    outer_radius = 1,
    inner_radius = 0.45,
    count_irregularity = 0.2,
    radial_irregularity = 0.35,
    phase_twist = 0.35,
    surface = c("flat", "saddle", "paraboloid", "ripple", "folded"),
    amplitude = 0.6,
    freq_u = 1,
    freq_v = 1) {
  surface <- match.arg(surface)
  built <- .irregular.annulus.canonical(
    rings = rings,
    outer_count = outer_count,
    outer_radius = outer_radius,
    inner_radius = inner_radius,
    count_irregularity = count_irregularity,
    radial_irregularity = radial_irregularity,
    phase_twist = phase_twist
  )
  .triangulated.planar.surface.coords(
    coords_param = built$coords,
    surface = surface,
    amplitude = amplitude,
    freq_u = freq_u,
    freq_v = freq_v
  )
}

#' @rdname irregular_annulus_surface_helpers
#' @export
irregular.annulus.surface.graph <- function(
    rings = 6,
    outer_count = 28,
    outer_radius = 1,
    inner_radius = 0.45,
    count_irregularity = 0.2,
    radial_irregularity = 0.35,
    phase_twist = 0.35,
    surface = c("flat", "saddle", "paraboloid", "ripple", "folded"),
    amplitude = 0.6,
    freq_u = 1,
    freq_v = 1,
    normalize = c("median", "mean", "none")) {
  surface <- match.arg(surface)
  normalize <- match.arg(normalize)
  built <- .irregular.annulus.canonical(
    rings = rings,
    outer_count = outer_count,
    outer_radius = outer_radius,
    inner_radius = inner_radius,
    count_irregularity = count_irregularity,
    radial_irregularity = radial_irregularity,
    phase_twist = phase_twist
  )
  coords_surface <- irregular.annulus.surface.embedding(
    rings = rings,
    outer_count = outer_count,
    outer_radius = outer_radius,
    inner_radius = inner_radius,
    count_irregularity = count_irregularity,
    radial_irregularity = radial_irregularity,
    phase_twist = phase_twist,
    surface = surface,
    amplitude = amplitude,
    freq_u = freq_u,
    freq_v = freq_v
  )
  weights <- .edge.weights.from.embedding(
    edges = built$edges,
    coords = coords_surface,
    normalize = normalize
  )

  out <- list(
    edges = built$edges,
    n = built$n,
    edge_weights = weights$edge_weights,
    coords_surface = coords_surface,
    coords_param = built$coords,
    weight_scale = weights$weight_scale,
    family = "irregular.annulus",
    surface = surface,
    rings = built$rings,
    ring_sizes = built$ring_sizes,
    normalize = normalize,
    label = sprintf("%s irregular annulus rings %d",
                    tools::toTitleCase(surface),
                    built$rings)
  )
  class(out) <- c("grip_irregular_annulus_surface_graph", "list")
  out
}

#' Weighted irregular pair-of-pants surface helpers
#'
#' Convenience helpers that build a pair-of-pants domain from deterministically
#' irregular horizontal slices, sample each surviving slice interval with
#' varying point counts, and stitch adjacent slices into a locally triangulated
#' surface-with-boundary graph. The resulting irregular planar parameterization
#' can then be lifted into \eqn{\mathbb{R}^3} using the same planar surface
#' families as the triangulated annulus and irregular annulus helpers.
#'
#' @param slices Number of horizontal sample slices through the pair-of-pants
#'   domain.
#' @param outer_count Approximate number of vertices across the widest slices.
#' @param outer_radius Positive outer boundary radius.
#' @param hole_radius Positive radius of each interior hole.
#' @param hole_offset Positive horizontal offset of the two hole centers.
#' @param hole_height Shared vertical coordinate of the two hole centers.
#' @param count_irregularity Irregularity level for the per-slice sample counts.
#'   Must lie in \code{[0, 1)}.
#' @param vertical_irregularity Irregularity level for slice spacing. Must lie
#'   in \code{[0, 1]}.
#' @param phase_twist Finite phase offset used to desynchronize neighboring
#'   slice samples.
#' @param surface Geometry family used for the 3D lift. One of \code{"flat"},
#'   \code{"saddle"}, \code{"paraboloid"}, \code{"ripple"}, or
#'   \code{"folded"}.
#' @param amplitude Finite deformation amplitude.
#' @param freq_u Positive ripple frequency in the first planar coordinate. Used
#'   only when \code{surface = "ripple"}.
#' @param freq_v Positive ripple frequency in the second planar coordinate. Used
#'   only when \code{surface = "ripple"}.
#' @param normalize Normalization applied to the induced edge lengths. One of
#'   \code{"median"}, \code{"mean"}, or \code{"none"}.
#'
#' @return
#' \code{irregular.pair.of.pants.surface.embedding()} returns an \code{n x 3}
#' numeric matrix with columns \code{x}, \code{y}, and \code{z}.
#'
#' \code{irregular.pair.of.pants.surface.graph()} returns a list with
#' components:
#' \itemize{
#'   \item \code{edges}: the irregular pair-of-pants edges,
#'   \item \code{n}: number of vertices,
#'   \item \code{edge_weights}: induced positive edge lengths,
#'   \item \code{coords_surface}: the 3D embedding,
#'   \item \code{coords_param}: the irregular planar pair-of-pants coordinates,
#'   \item \code{weight_scale}: the normalization constant applied to the raw
#'     edge lengths,
#'   \item \code{family}: always \code{"irregular.pair.of.pants"},
#'   \item \code{surface}: the chosen surface name,
#'   \item \code{slices}: number of horizontal sample slices,
#'   \item \code{slice_sizes}: vertex counts in each slice,
#'   \item \code{slice_components}: number of connected slice intervals in each
#'     slice,
#'   \item \code{label}: a human-readable family label.
#' }
#'
#' @name irregular_pair_of_pants_surface_helpers
NULL

#' @rdname irregular_pair_of_pants_surface_helpers
#' @export
irregular.pair.of.pants.surface.embedding <- function(
    slices = 11,
    outer_count = 28,
    outer_radius = 1.1,
    hole_radius = 0.24,
    hole_offset = 0.38,
    hole_height = 0.18,
    count_irregularity = 0.2,
    vertical_irregularity = 0.35,
    phase_twist = 0.35,
    surface = c("flat", "saddle", "paraboloid", "ripple", "folded"),
    amplitude = 0.6,
    freq_u = 1,
    freq_v = 1) {
  surface <- match.arg(surface)
  built <- .irregular.pair.of.pants.canonical(
    slices = slices,
    outer_count = outer_count,
    outer_radius = outer_radius,
    hole_radius = hole_radius,
    hole_offset = hole_offset,
    hole_height = hole_height,
    count_irregularity = count_irregularity,
    vertical_irregularity = vertical_irregularity,
    phase_twist = phase_twist
  )
  .triangulated.planar.surface.coords(
    coords_param = built$coords,
    surface = surface,
    amplitude = amplitude,
    freq_u = freq_u,
    freq_v = freq_v
  )
}

#' @rdname irregular_pair_of_pants_surface_helpers
#' @export
irregular.pair.of.pants.surface.graph <- function(
    slices = 11,
    outer_count = 28,
    outer_radius = 1.1,
    hole_radius = 0.24,
    hole_offset = 0.38,
    hole_height = 0.18,
    count_irregularity = 0.2,
    vertical_irregularity = 0.35,
    phase_twist = 0.35,
    surface = c("flat", "saddle", "paraboloid", "ripple", "folded"),
    amplitude = 0.6,
    freq_u = 1,
    freq_v = 1,
    normalize = c("median", "mean", "none")) {
  surface <- match.arg(surface)
  normalize <- match.arg(normalize)
  built <- .irregular.pair.of.pants.canonical(
    slices = slices,
    outer_count = outer_count,
    outer_radius = outer_radius,
    hole_radius = hole_radius,
    hole_offset = hole_offset,
    hole_height = hole_height,
    count_irregularity = count_irregularity,
    vertical_irregularity = vertical_irregularity,
    phase_twist = phase_twist
  )
  coords_surface <- irregular.pair.of.pants.surface.embedding(
    slices = slices,
    outer_count = outer_count,
    outer_radius = outer_radius,
    hole_radius = hole_radius,
    hole_offset = hole_offset,
    hole_height = hole_height,
    count_irregularity = count_irregularity,
    vertical_irregularity = vertical_irregularity,
    phase_twist = phase_twist,
    surface = surface,
    amplitude = amplitude,
    freq_u = freq_u,
    freq_v = freq_v
  )
  weights <- .edge.weights.from.embedding(
    edges = built$edges,
    coords = coords_surface,
    normalize = normalize
  )

  out <- list(
    edges = built$edges,
    n = built$n,
    edge_weights = weights$edge_weights,
    coords_surface = coords_surface,
    coords_param = built$coords,
    weight_scale = weights$weight_scale,
    family = "irregular.pair.of.pants",
    surface = surface,
    slices = built$slices,
    slice_sizes = built$slice_sizes,
    slice_components = built$slice_components,
    normalize = normalize,
    label = sprintf("%s irregular pair-of-pants slices %d",
                    tools::toTitleCase(surface),
                    built$slices)
  )
  class(out) <- c("grip_irregular_pair_of_pants_surface_graph", "list")
  out
}

#' Weighted Sierpinski triangle surface helpers
#'
#' Convenience helpers that take the canonical recursive embedding of the
#' Sierpinski triangle graph and lift it into \eqn{\mathbb{R}^3}. These helpers
#' are intended for benchmark families where the topology is the usual
#' triangular gasket but the intended metric comes from a curved or folded
#' geometry rather than the flat equilateral reference drawing.
#'
#' @param level Recursion depth. May be \code{0} or larger.
#' @param surface Triangle surface family. One of \code{"flat"},
#'   \code{"saddle"}, \code{"paraboloid"}, \code{"ripple"}, or
#'   \code{"folded"}.
#' @param amplitude Finite deformation amplitude.
#' @param freq_u Positive ripple frequency in the first canonical triangle
#'   coordinate. Used only when \code{surface = "ripple"}.
#' @param freq_v Positive ripple frequency in the second canonical triangle
#'   coordinate. Used only when \code{surface = "ripple"}.
#' @param x_scale Positive horizontal scaling applied to the canonical
#'   triangle coordinates.
#' @param y_scale Positive vertical scaling applied to the canonical triangle
#'   coordinates.
#' @param normalize Normalization applied to the induced edge lengths. One of
#'   \code{"median"}, \code{"mean"}, or \code{"none"}.
#'
#' @return
#' \code{sierpinski.triangle.surface.embedding()} returns an \code{n x 3}
#' numeric matrix with columns \code{x}, \code{y}, and \code{z}.
#'
#' \code{sierpinski.triangle.surface.graph()} returns a list with components:
#' \itemize{
#'   \item \code{edges}: the Sierpinski triangle edges,
#'   \item \code{n}: number of vertices,
#'   \item \code{edge_weights}: induced positive edge lengths,
#'   \item \code{coords_surface}: the 3D surface embedding,
#'   \item \code{coords_param}: the canonical 2D triangle coordinates,
#'   \item \code{weight_scale}: the normalization constant applied to the raw
#'     edge lengths,
#'   \item \code{family}: always \code{"sierpinski.triangle"},
#'   \item \code{surface}: the chosen surface name,
#'   \item \code{level}: the recursion depth,
#'   \item \code{label}: a human-readable family label.
#' }
#'
#' @name sierpinski_triangle_surface_helpers
NULL

#' @rdname sierpinski_triangle_surface_helpers
#' @export
sierpinski.triangle.surface.embedding <- function(level = 2,
                                                  surface = c("flat", "saddle",
                                                              "paraboloid", "ripple",
                                                              "folded"),
                                                  amplitude = 0.75,
                                                  freq_u = 1,
                                                  freq_v = 1,
                                                  x_scale = 1,
                                                  y_scale = 1) {
  level <- .as_whole_number(level, "level")
  surface <- match.arg(surface)
  amplitude <- .as_finite_scalar(amplitude, "amplitude")
  freq_u <- .as_positive_scalar(freq_u, "freq_u")
  freq_v <- .as_positive_scalar(freq_v, "freq_v")
  x_scale <- .as_positive_scalar(x_scale, "x_scale")
  y_scale <- .as_positive_scalar(y_scale, "y_scale")

  built <- .sierpinski.triangle.canonical(level)
  coords_param <- cbind(
    u = built$coords[, 1L] * x_scale,
    v = built$coords[, 2L] * y_scale
  )
  coords_norm <- .normalize.center_radius.coords(coords_param)
  z <- switch(
    surface,
    flat = rep(0, nrow(coords_param)),
    saddle = .surface.z.from_uv(coords_norm[, 1L], coords_norm[, 2L],
                                "saddle", amplitude, freq_u, freq_v),
    paraboloid = .surface.z.from_uv(coords_norm[, 1L], coords_norm[, 2L],
                                    "paraboloid", amplitude, freq_u, freq_v),
    ripple = .surface.z.from_uv(coords_norm[, 1L], coords_norm[, 2L],
                                "ripple", amplitude, freq_u, freq_v),
    folded = amplitude * abs(coords_norm[, 1L])
  )

  coords <- cbind(x = coords_param[, 1L], y = coords_param[, 2L], z = z)
  storage.mode(coords) <- "double"
  coords
}

#' @rdname sierpinski_triangle_surface_helpers
#' @export
sierpinski.triangle.surface.graph <- function(level = 2,
                                              surface = c("flat", "saddle",
                                                          "paraboloid", "ripple",
                                                          "folded"),
                                              amplitude = 0.75,
                                              freq_u = 1,
                                              freq_v = 1,
                                              x_scale = 1,
                                              y_scale = 1,
                                              normalize = c("median", "mean", "none")) {
  level <- .as_whole_number(level, "level")
  surface <- match.arg(surface)
  normalize <- match.arg(normalize)

  built <- .sierpinski.triangle.canonical(level)
  coords_param <- cbind(
    u = built$coords[, 1L] * x_scale,
    v = built$coords[, 2L] * y_scale
  )
  coords_surface <- sierpinski.triangle.surface.embedding(
    level = level,
    surface = surface,
    amplitude = amplitude,
    freq_u = freq_u,
    freq_v = freq_v,
    x_scale = x_scale,
    y_scale = y_scale
  )
  weights <- .edge.weights.from.embedding(
    edges = built$edges,
    coords = coords_surface,
    normalize = normalize
  )

  out <- list(
    edges = built$edges,
    n = built$n,
    edge_weights = weights$edge_weights,
    coords_surface = coords_surface,
    coords_param = coords_param,
    weight_scale = weights$weight_scale,
    family = "sierpinski.triangle",
    surface = surface,
    level = level,
    normalize = normalize,
    label = sprintf("%s Sierpinski triangle level %d",
                    tools::toTitleCase(surface),
                    level)
  )
  class(out) <- c("grip_sierpinski_triangle_surface_graph", "list")
  out
}

.tetrahedron.surface.coords <- function(coords_param,
                                        surface,
                                        amplitude,
                                        freq,
                                        twist) {
  surface <- .as_named_choice(surface, c("standard", "squashed", "twisted", "wavy"),
                              "surface")
  amplitude <- .as_finite_scalar(amplitude, "amplitude")
  freq <- .as_positive_scalar(freq, "freq")
  twist <- .as_finite_scalar(twist, "twist")

  coords_param <- as.matrix(coords_param)
  coords_norm <- .normalize.center_radius.coords(coords_param)
  coords <- switch(
    surface,
    standard = coords_param,
    squashed = {
      xy_scale <- 1 - amplitude / 2
      z_scale <- 1 + amplitude
      if (xy_scale <= 0 || z_scale <= 0) {
        stop("squashed tetrahedron parameters require positive coordinate scales; adjust amplitude",
             call. = FALSE)
      }
      cbind(
        x = coords_param[, 1L] * xy_scale,
        y = coords_param[, 2L] * xy_scale,
        z = coords_param[, 3L] * z_scale
      )
    },
    twisted = {
      theta <- twist * coords_norm[, 3L]
      x <- coords_param[, 1L] * cos(theta) - coords_param[, 2L] * sin(theta)
      y <- coords_param[, 1L] * sin(theta) + coords_param[, 2L] * cos(theta)
      cbind(x = x, y = y, z = coords_param[, 3L])
    },
    wavy = {
      local_scale <- 1 + amplitude * sin(pi * freq *
                                           (coords_norm[, 1L] +
                                              0.5 * coords_norm[, 2L] -
                                              0.25 * coords_norm[, 3L]))
      if (any(!is.finite(local_scale) | local_scale <= 0)) {
        stop("wavy tetrahedron parameters produce a non-positive local scale; adjust amplitude or frequency",
             call. = FALSE)
      }
      coords_param * local_scale
    }
  )

  storage.mode(coords) <- "double"
  coords
}

#' Weighted recursive tetrahedron-mask surface helpers
#'
#' Convenience helpers that recursively subdivide a tetrahedron into four
#' corner subtetrahedra according to a four-corner keep-mask, retain the
#' selected subtetrahedra at each recursion step, and deform the resulting
#' simplicial graph inside \eqn{\mathbb{R}^3}. This provides a generic
#' tetrahedral gasket family that includes the classic Sierpinski tetrahedron
#' together with corner-omission variants.
#'
#' The \code{mask} entries are interpreted in the order
#' \code{base_left}, \code{base_right}, \code{base_back}, and \code{apex}.
#' Named masks are reordered automatically to match that convention.
#'
#' @param mask Four-entry logical or numeric keep-mask. Non-zero entries are
#'   retained at each recursive step.
#' @param level Recursion depth. May be \code{0} or larger.
#' @param surface Tetrahedron surface family. One of \code{"standard"},
#'   \code{"squashed"}, \code{"twisted"}, or \code{"wavy"}.
#' @param amplitude Finite deformation amplitude.
#' @param freq Positive modulation frequency used only when
#'   \code{surface = "wavy"}.
#' @param twist Finite twist strength used only when
#'   \code{surface = "twisted"}.
#' @param normalize Normalization applied to the induced edge lengths. One of
#'   \code{"median"}, \code{"mean"}, or \code{"none"}.
#'
#' @return
#' \code{recursive.tetrahedron.mask.surface.embedding()} returns an \code{n x 3}
#' numeric matrix with columns \code{x}, \code{y}, and \code{z}.
#'
#' \code{recursive.tetrahedron.mask.surface.graph()} returns a list with
#' components:
#' \itemize{
#'   \item \code{edges}: the retained recursive tetrahedron-mask edges,
#'   \item \code{n}: number of vertices,
#'   \item \code{edge_weights}: induced positive edge lengths,
#'   \item \code{coords_surface}: the 3D surface embedding,
#'   \item \code{coords_param}: the canonical 3D tetrahedron coordinates,
#'   \item \code{weight_scale}: the normalization constant applied to the raw
#'     edge lengths,
#'   \item \code{family}: always \code{"recursive.tetrahedron.mask"},
#'   \item \code{surface}: the chosen surface name,
#'   \item \code{level}: the recursion depth,
#'   \item \code{mask}: the logical keep-mask,
#'   \item \code{label}: a human-readable family label.
#' }
#'
#' @name recursive_tetrahedron_mask_surface_helpers
NULL

#' @rdname recursive_tetrahedron_mask_surface_helpers
#' @export
recursive.tetrahedron.mask.surface.embedding <- function(
    mask = mask.tetrahedron.classic(),
    level = 2,
    surface = c("standard", "squashed", "twisted", "wavy"),
    amplitude = 0.3,
    freq = 2,
    twist = 0.6) {
  mask <- .as_tetrahedron_keep_mask(mask, "mask")
  level <- .as_whole_number(level, "level")
  surface <- match.arg(surface)

  built <- .recursive.tetrahedron.mask.canonical(mask, level)
  .tetrahedron.surface.coords(
    coords_param = built$coords,
    surface = surface,
    amplitude = amplitude,
    freq = freq,
    twist = twist
  )
}

#' @rdname recursive_tetrahedron_mask_surface_helpers
#' @export
recursive.tetrahedron.mask.surface.graph <- function(
    mask = mask.tetrahedron.classic(),
    level = 2,
    surface = c("standard", "squashed", "twisted", "wavy"),
    amplitude = 0.3,
    freq = 2,
    twist = 0.6,
    normalize = c("median", "mean", "none")) {
  mask <- .as_tetrahedron_keep_mask(mask, "mask")
  level <- .as_whole_number(level, "level")
  surface <- match.arg(surface)
  normalize <- match.arg(normalize)

  built <- .recursive.tetrahedron.mask.canonical(mask, level)
  coords_surface <- recursive.tetrahedron.mask.surface.embedding(
    mask = mask,
    level = level,
    surface = surface,
    amplitude = amplitude,
    freq = freq,
    twist = twist
  )
  weights <- .edge.weights.from.embedding(
    edges = built$edges,
    coords = coords_surface,
    normalize = normalize
  )

  out <- list(
    edges = built$edges,
    n = built$n,
    edge_weights = weights$edge_weights,
    coords_surface = coords_surface,
    coords_param = built$coords,
    weight_scale = weights$weight_scale,
    family = "recursive.tetrahedron.mask",
    surface = surface,
    level = level,
    mask = mask,
    normalize = normalize,
    label = sprintf("%s %s level %d",
                    tools::toTitleCase(surface),
                    .tetrahedron.mask.label(mask),
                    level)
  )
  class(out) <- c("grip_recursive_tetrahedron_mask_surface_graph", "list")
  out
}

#' Weighted Sierpinski tetrahedron surface helpers
#'
#' Convenience wrappers around \code{recursive.tetrahedron.mask.surface.*()}
#' using the classic all-corner tetrahedron mask. These helpers are intended
#' for benchmark families where the topology is the standard simplicial
#' tetrahedral gasket but the intended metric comes from a squashed, twisted,
#' or spatially varying realization.
#'
#' @inheritParams recursive.tetrahedron.mask.surface.embedding
#'
#' @return
#' \code{sierpinski.tetrahedron.surface.embedding()} returns an \code{n x 3}
#' numeric matrix with columns \code{x}, \code{y}, and \code{z}.
#'
#' \code{sierpinski.tetrahedron.surface.graph()} returns the same components as
#' \code{recursive.tetrahedron.mask.surface.graph()}, with \code{family} set to
#' \code{"sierpinski.tetrahedron"} and a family-specific class and label.
#'
#' @name sierpinski_tetrahedron_surface_helpers
NULL

#' @rdname sierpinski_tetrahedron_surface_helpers
#' @export
sierpinski.tetrahedron.surface.embedding <- function(level = 2,
                                                     surface = c("standard", "squashed",
                                                                 "twisted", "wavy"),
                                                     amplitude = 0.3,
                                                     freq = 2,
                                                     twist = 0.6) {
  recursive.tetrahedron.mask.surface.embedding(
    mask = .tetrahedron.classic.mask(),
    level = level,
    surface = surface,
    amplitude = amplitude,
    freq = freq,
    twist = twist
  )
}

#' @rdname sierpinski_tetrahedron_surface_helpers
#' @export
sierpinski.tetrahedron.surface.graph <- function(level = 2,
                                                 surface = c("standard", "squashed",
                                                             "twisted", "wavy"),
                                                 amplitude = 0.3,
                                                 freq = 2,
                                                 twist = 0.6,
                                                 normalize = c("median", "mean", "none")) {
  out <- recursive.tetrahedron.mask.surface.graph(
    mask = .tetrahedron.classic.mask(),
    level = level,
    surface = surface,
    amplitude = amplitude,
    freq = freq,
    twist = twist,
    normalize = normalize
  )
  out$family <- "sierpinski.tetrahedron"
  out$label <- sprintf("%s Sierpinski tetrahedron level %d",
                       tools::toTitleCase(out$surface),
                       out$level)
  class(out) <- c("grip_sierpinski_tetrahedron_surface_graph", class(out))
  out
}

.cube.mask.surface.coords <- function(coords_param,
                                      surface,
                                      amplitude,
                                      freq,
                                      twist) {
  surface <- .as_named_choice(surface, c("standard", "bulged", "twisted", "wavy"),
                              "surface")
  amplitude <- .as_finite_scalar(amplitude, "amplitude")
  freq <- .as_positive_scalar(freq, "freq")
  twist <- .as_finite_scalar(twist, "twist")

  coords_param <- as.matrix(coords_param)
  coords_norm <- .normalize.center_radius.coords(coords_param)
  coords <- switch(
    surface,
    standard = coords_param,
    bulged = {
      local_scale <- 1 + amplitude * rowSums(coords_norm^2)
      if (any(!is.finite(local_scale) | local_scale <= 0)) {
        stop("bulged cube-mask parameters produce a non-positive local scale; adjust amplitude",
             call. = FALSE)
      }
      coords_param * local_scale
    },
    twisted = {
      theta <- twist * coords_norm[, 3L]
      x <- coords_param[, 1L] * cos(theta) - coords_param[, 2L] * sin(theta)
      y <- coords_param[, 1L] * sin(theta) + coords_param[, 2L] * cos(theta)
      cbind(x = x, y = y, z = coords_param[, 3L])
    },
    wavy = {
      dx <- amplitude * sin(pi * freq * coords_norm[, 2L]) *
        cos(pi * freq * coords_norm[, 3L])
      dy <- amplitude * sin(pi * freq * coords_norm[, 3L]) *
        cos(pi * freq * coords_norm[, 1L])
      dz <- amplitude * sin(pi * freq * coords_norm[, 1L]) *
        cos(pi * freq * coords_norm[, 2L])
      coords_param + cbind(dx, dy, dz)
    }
  )

  storage.mode(coords) <- "double"
  coords
}

#' Weighted recursive cube-mask surface helpers
#'
#' Convenience helpers that recursively subdivide a cube into a regular
#' \eqn{k \times k \times k} array of subcubes according to a cubic keep-mask,
#' retain the selected cells at each recursion step, and deform the resulting
#' face-adjacency graph inside \eqn{\mathbb{R}^3}. This provides a generic
#' cubical 3D family that includes Menger-sponge-style cell adjacency as a
#' named special case.
#'
#' The \code{mask} is interpreted as a cubic array whose first dimension runs
#' from top to bottom, second from left to right, and third from front to back.
#' Non-zero entries are retained at each recursive step.
#'
#' @param mask Cubic logical or numeric keep-array. Non-zero entries are
#'   retained at each recursive step.
#' @param level Recursion depth. Must be at least \code{1}.
#' @param surface Cube-mask geometry family. One of \code{"standard"},
#'   \code{"bulged"}, \code{"twisted"}, or \code{"wavy"}.
#' @param amplitude Finite deformation amplitude.
#' @param freq Positive modulation frequency used only when
#'   \code{surface = "wavy"}.
#' @param twist Finite twist strength used only when
#'   \code{surface = "twisted"}.
#' @param x_scale Positive horizontal scaling applied to the canonical cube
#'   coordinates.
#' @param y_scale Positive vertical scaling applied to the canonical cube
#'   coordinates.
#' @param z_scale Positive depth scaling applied to the canonical cube
#'   coordinates.
#' @param normalize Normalization applied to the induced edge lengths. One of
#'   \code{"median"}, \code{"mean"}, or \code{"none"}.
#'
#' @return
#' \code{recursive.cube.mask.surface.embedding()} returns an \code{n x 3}
#' numeric matrix with columns \code{x}, \code{y}, and \code{z}.
#'
#' \code{recursive.cube.mask.surface.graph()} returns a list with components:
#' \itemize{
#'   \item \code{edges}: the retained recursive cube-mask face-adjacency edges,
#'   \item \code{n}: number of vertices,
#'   \item \code{edge_weights}: induced positive edge lengths,
#'   \item \code{coords_surface}: the 3D embedding,
#'   \item \code{coords_param}: the canonical 3D voxel-center coordinates,
#'   \item \code{weight_scale}: the normalization constant applied to the raw
#'     edge lengths,
#'   \item \code{family}: always \code{"recursive.cube.mask"},
#'   \item \code{surface}: the chosen surface name,
#'   \item \code{level}: the recursion depth,
#'   \item \code{mask}: the logical keep-array,
#'   \item \code{mask_side}: side length of the base keep-mask,
#'   \item \code{side}: side length of the fully refined cube grid,
#'   \item \code{label}: a human-readable family label.
#' }
#'
#' @name recursive_cube_mask_surface_helpers
NULL

#' @rdname recursive_cube_mask_surface_helpers
#' @export
recursive.cube.mask.surface.embedding <- function(mask,
                                                  level = 2,
                                                  surface = c("standard", "bulged",
                                                              "twisted", "wavy"),
                                                  amplitude = 0.2,
                                                  freq = 2,
                                                  twist = 0.6,
                                                  x_scale = 1,
                                                  y_scale = 1,
                                                  z_scale = 1) {
  mask <- .as_cube_keep_mask(mask, "mask", min_size = 2L)
  level <- .as_whole_number(level, "level", min = 1L)
  surface <- match.arg(surface)
  x_scale <- .as_positive_scalar(x_scale, "x_scale")
  y_scale <- .as_positive_scalar(y_scale, "y_scale")
  z_scale <- .as_positive_scalar(z_scale, "z_scale")

  coords_param <- .recursive.cube.mask.param.coords(
    mask = mask,
    level = level,
    x_scale = x_scale,
    y_scale = y_scale,
    z_scale = z_scale
  )
  .cube.mask.surface.coords(
    coords_param = coords_param,
    surface = surface,
    amplitude = amplitude,
    freq = freq,
    twist = twist
  )
}

#' @rdname recursive_cube_mask_surface_helpers
#' @export
recursive.cube.mask.surface.graph <- function(mask,
                                              level = 2,
                                              surface = c("standard", "bulged",
                                                          "twisted", "wavy"),
                                              amplitude = 0.2,
                                              freq = 2,
                                              twist = 0.6,
                                              x_scale = 1,
                                              y_scale = 1,
                                              z_scale = 1,
                                              normalize = c("median", "mean", "none")) {
  mask <- .as_cube_keep_mask(mask, "mask", min_size = 2L)
  level <- .as_whole_number(level, "level", min = 1L)
  surface <- match.arg(surface)
  normalize <- match.arg(normalize)

  spec <- .recursive.cube.mask.cells(mask, level)
  edges <- .recursive.cube.mask.edges(mask, level)
  coords_param <- .recursive.cube.mask.param.coords(
    mask = mask,
    level = level,
    x_scale = x_scale,
    y_scale = y_scale,
    z_scale = z_scale
  )
  coords_surface <- recursive.cube.mask.surface.embedding(
    mask = mask,
    level = level,
    surface = surface,
    amplitude = amplitude,
    freq = freq,
    twist = twist,
    x_scale = x_scale,
    y_scale = y_scale,
    z_scale = z_scale
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
    family = "recursive.cube.mask",
    surface = surface,
    level = level,
    mask = mask,
    mask_side = as.integer(dim(mask)[1L]),
    side = spec$side,
    normalize = normalize,
    label = sprintf("%s recursive cube mask level %d",
                    tools::toTitleCase(surface),
                    level)
  )
  class(out) <- c("grip_recursive_cube_mask_surface_graph", "list")
  out
}

#' Weighted Menger sponge surface helpers
#'
#' Convenience wrappers around \code{recursive.cube.mask.surface.*()} using the
#' classic \eqn{3 \times 3 \times 3} Menger-sponge keep-mask. These helpers
#' provide a named cubical 3D benchmark family whose vertices are occupied
#' subcubes and whose edges connect face-adjacent occupied cells.
#'
#' @inheritParams recursive.cube.mask.surface.embedding
#'
#' @return
#' \code{menger.sponge.surface.embedding()} returns an \code{n x 3} numeric
#' matrix with columns \code{x}, \code{y}, and \code{z}, where
#' \code{n = 20^level}.
#'
#' \code{menger.sponge.surface.graph()} returns the same components as
#' \code{recursive.cube.mask.surface.graph()}, with \code{family} set to
#' \code{"menger.sponge"} and a family-specific class and label.
#'
#' @name menger_sponge_surface_helpers
NULL

#' @rdname menger_sponge_surface_helpers
#' @export
menger.sponge.surface.embedding <- function(level = 2,
                                            surface = c("standard", "bulged",
                                                        "twisted", "wavy"),
                                            amplitude = 0.2,
                                            freq = 2,
                                            twist = 0.6,
                                            x_scale = 1,
                                            y_scale = 1,
                                            z_scale = 1) {
  recursive.cube.mask.surface.embedding(
    mask = .menger.sponge.mask(),
    level = level,
    surface = surface,
    amplitude = amplitude,
    freq = freq,
    twist = twist,
    x_scale = x_scale,
    y_scale = y_scale,
    z_scale = z_scale
  )
}

#' @rdname menger_sponge_surface_helpers
#' @export
menger.sponge.surface.graph <- function(level = 2,
                                        surface = c("standard", "bulged",
                                                    "twisted", "wavy"),
                                        amplitude = 0.2,
                                        freq = 2,
                                        twist = 0.6,
                                        x_scale = 1,
                                        y_scale = 1,
                                        z_scale = 1,
                                        normalize = c("median", "mean", "none")) {
  out <- recursive.cube.mask.surface.graph(
    mask = .menger.sponge.mask(),
    level = level,
    surface = surface,
    amplitude = amplitude,
    freq = freq,
    twist = twist,
    x_scale = x_scale,
    y_scale = y_scale,
    z_scale = z_scale,
    normalize = normalize
  )
  out$family <- "menger.sponge"
  out$label <- sprintf("%s Menger sponge level %d",
                       tools::toTitleCase(out$surface),
                       out$level)
  class(out) <- c("grip_menger_sponge_surface_graph", class(out))
  out
}

#' Weighted porous cube-mask surface helpers
#'
#' Convenience wrappers around \code{recursive.cube.mask.surface.*()} for three
#' named porous cubical families: a periodic tunnel lattice, an asymmetric
#' cavity family, and a deterministic channel network. Each family uses a
#' reusable base keep-array and then recurses it in exactly the same way as the
#' generic cube-mask helpers.
#'
#' @param level Recursion depth. Must be at least \code{1}.
#' @param side Side length of the base cubic keep-array.
#' @param surface Cube-mask geometry family. One of \code{"standard"},
#'   \code{"bulged"}, \code{"twisted"}, or \code{"wavy"}.
#' @param amplitude Finite deformation amplitude.
#' @param freq Positive modulation frequency used only when
#'   \code{surface = "wavy"}.
#' @param twist Finite twist strength used only when
#'   \code{surface = "twisted"}.
#' @param x_scale Positive horizontal scaling applied to the canonical cube
#'   coordinates.
#' @param y_scale Positive vertical scaling applied to the canonical cube
#'   coordinates.
#' @param z_scale Positive depth scaling applied to the canonical cube
#'   coordinates.
#' @param normalize Normalization applied to the induced edge lengths. One of
#'   \code{"median"}, \code{"mean"}, or \code{"none"}.
#' @param tunnel_width Width of each removed tunnel band.
#' @param tunnel_period Spacing between successive tunnel bands.
#' @param tunnel_offset Starting index of the first tunnel band.
#' @param cavity_size Side length of the larger interior cavity block.
#' @param pocket_size Side length of the smaller secondary cavity block.
#' @param channel_width Width of each removed channel in the channel-network
#'   family.
#' @param branch_offset Interior offset of the extra branch channel in the
#'   channel-network family.
#'
#' @return
#' Each \code{*.surface.embedding()} helper returns an \code{n x 3} numeric
#' matrix with columns \code{x}, \code{y}, and \code{z}.
#'
#' Each \code{*.surface.graph()} helper returns the same components as
#' \code{recursive.cube.mask.surface.graph()}, with a family-specific
#' \code{family} field, class, and label.
#'
#' @name porous_cube_surface_helpers
NULL

#' @rdname porous_cube_surface_helpers
#' @export
cube.periodic.tunnels.surface.embedding <- function(level = 2,
                                                    side = 5,
                                                    tunnel_width = 1,
                                                    tunnel_period = 2,
                                                    tunnel_offset = 2,
                                                    surface = c("standard", "bulged",
                                                                "twisted", "wavy"),
                                                    amplitude = 0.2,
                                                    freq = 2,
                                                    twist = 0.6,
                                                    x_scale = 1,
                                                    y_scale = 1,
                                                    z_scale = 1) {
  recursive.cube.mask.surface.embedding(
    mask = mask.cube.periodic.tunnels(
      side = side,
      tunnel_width = tunnel_width,
      tunnel_period = tunnel_period,
      tunnel_offset = tunnel_offset
    ),
    level = level,
    surface = surface,
    amplitude = amplitude,
    freq = freq,
    twist = twist,
    x_scale = x_scale,
    y_scale = y_scale,
    z_scale = z_scale
  )
}

#' @rdname porous_cube_surface_helpers
#' @export
cube.periodic.tunnels.surface.graph <- function(level = 2,
                                                side = 5,
                                                tunnel_width = 1,
                                                tunnel_period = 2,
                                                tunnel_offset = 2,
                                                surface = c("standard", "bulged",
                                                            "twisted", "wavy"),
                                                amplitude = 0.2,
                                                freq = 2,
                                                twist = 0.6,
                                                x_scale = 1,
                                                y_scale = 1,
                                                z_scale = 1,
                                                normalize = c("median", "mean", "none")) {
  out <- recursive.cube.mask.surface.graph(
    mask = mask.cube.periodic.tunnels(
      side = side,
      tunnel_width = tunnel_width,
      tunnel_period = tunnel_period,
      tunnel_offset = tunnel_offset
    ),
    level = level,
    surface = surface,
    amplitude = amplitude,
    freq = freq,
    twist = twist,
    x_scale = x_scale,
    y_scale = y_scale,
    z_scale = z_scale,
    normalize = normalize
  )
  out$family <- "cube.periodic.tunnels"
  out$label <- sprintf("%s periodic tunnels cube level %d",
                       tools::toTitleCase(out$surface),
                       out$level)
  class(out) <- c("grip_cube_periodic_tunnels_surface_graph", class(out))
  out
}

#' @rdname porous_cube_surface_helpers
#' @export
cube.asymmetric.cavities.surface.embedding <- function(level = 2,
                                                       side = 5,
                                                       cavity_size = 2,
                                                       pocket_size = max(1L, cavity_size - 1L),
                                                       surface = c("standard", "bulged",
                                                                   "twisted", "wavy"),
                                                       amplitude = 0.2,
                                                       freq = 2,
                                                       twist = 0.6,
                                                       x_scale = 1,
                                                       y_scale = 1,
                                                       z_scale = 1) {
  recursive.cube.mask.surface.embedding(
    mask = mask.cube.asymmetric.cavities(
      side = side,
      cavity_size = cavity_size,
      pocket_size = pocket_size
    ),
    level = level,
    surface = surface,
    amplitude = amplitude,
    freq = freq,
    twist = twist,
    x_scale = x_scale,
    y_scale = y_scale,
    z_scale = z_scale
  )
}

#' @rdname porous_cube_surface_helpers
#' @export
cube.asymmetric.cavities.surface.graph <- function(level = 2,
                                                   side = 5,
                                                   cavity_size = 2,
                                                   pocket_size = max(1L, cavity_size - 1L),
                                                   surface = c("standard", "bulged",
                                                               "twisted", "wavy"),
                                                   amplitude = 0.2,
                                                   freq = 2,
                                                   twist = 0.6,
                                                   x_scale = 1,
                                                   y_scale = 1,
                                                   z_scale = 1,
                                                   normalize = c("median", "mean", "none")) {
  out <- recursive.cube.mask.surface.graph(
    mask = mask.cube.asymmetric.cavities(
      side = side,
      cavity_size = cavity_size,
      pocket_size = pocket_size
    ),
    level = level,
    surface = surface,
    amplitude = amplitude,
    freq = freq,
    twist = twist,
    x_scale = x_scale,
    y_scale = y_scale,
    z_scale = z_scale,
    normalize = normalize
  )
  out$family <- "cube.asymmetric.cavities"
  out$label <- sprintf("%s asymmetric cavities cube level %d",
                       tools::toTitleCase(out$surface),
                       out$level)
  class(out) <- c("grip_cube_asymmetric_cavities_surface_graph", class(out))
  out
}

#' @rdname porous_cube_surface_helpers
#' @export
cube.channel.network.surface.embedding <- function(level = 2,
                                                   side = 5,
                                                   channel_width = 1,
                                                   branch_offset = 2,
                                                   surface = c("standard", "bulged",
                                                               "twisted", "wavy"),
                                                   amplitude = 0.2,
                                                   freq = 2,
                                                   twist = 0.6,
                                                   x_scale = 1,
                                                   y_scale = 1,
                                                   z_scale = 1) {
  recursive.cube.mask.surface.embedding(
    mask = mask.cube.channel.network(
      side = side,
      channel_width = channel_width,
      branch_offset = branch_offset
    ),
    level = level,
    surface = surface,
    amplitude = amplitude,
    freq = freq,
    twist = twist,
    x_scale = x_scale,
    y_scale = y_scale,
    z_scale = z_scale
  )
}

#' @rdname porous_cube_surface_helpers
#' @export
cube.channel.network.surface.graph <- function(level = 2,
                                               side = 5,
                                               channel_width = 1,
                                               branch_offset = 2,
                                               surface = c("standard", "bulged",
                                                           "twisted", "wavy"),
                                               amplitude = 0.2,
                                               freq = 2,
                                               twist = 0.6,
                                               x_scale = 1,
                                               y_scale = 1,
                                               z_scale = 1,
                                               normalize = c("median", "mean", "none")) {
  out <- recursive.cube.mask.surface.graph(
    mask = mask.cube.channel.network(
      side = side,
      channel_width = channel_width,
      branch_offset = branch_offset
    ),
    level = level,
    surface = surface,
    amplitude = amplitude,
    freq = freq,
    twist = twist,
    x_scale = x_scale,
    y_scale = y_scale,
    z_scale = z_scale,
    normalize = normalize
  )
  out$family <- "cube.channel.network"
  out$label <- sprintf("%s cube channel network level %d",
                       tools::toTitleCase(out$surface),
                       out$level)
  class(out) <- c("grip_cube_channel_network_surface_graph", class(out))
  out
}

#' Weighted recursive triangle-mask surface helpers
#'
#' Convenience helpers that recursively subdivide an equilateral triangle into
#' three corner subtriangles plus one central inverted subtriangle according to
#' a four-slot keep-mask, retain the selected subtriangles at each recursion
#' step, and lift the resulting vertex graph into \eqn{\mathbb{R}^3}. This
#' provides a generic gasket-style triangular family that includes the classic
#' Sierpinski triangle together with bridge and asymmetric variants.
#'
#' The \code{mask} entries are interpreted in the order
#' \code{left}, \code{right}, \code{top}, and \code{center}. Named masks are
#' reordered automatically to match that convention.
#'
#' @param mask Four-entry logical or numeric keep-mask. Non-zero entries are
#'   retained at each recursive step.
#' @param level Recursion depth. May be \code{0} or larger.
#' @param surface Triangle surface family. One of \code{"flat"},
#'   \code{"saddle"}, \code{"paraboloid"}, \code{"ripple"}, or
#'   \code{"folded"}.
#' @param amplitude Finite deformation amplitude.
#' @param freq_u Positive ripple frequency in the first canonical triangle
#'   coordinate. Used only when \code{surface = "ripple"}.
#' @param freq_v Positive ripple frequency in the second canonical triangle
#'   coordinate. Used only when \code{surface = "ripple"}.
#' @param x_scale Positive horizontal scaling applied to the canonical
#'   triangle coordinates.
#' @param y_scale Positive vertical scaling applied to the canonical triangle
#'   coordinates.
#' @param normalize Normalization applied to the induced edge lengths. One of
#'   \code{"median"}, \code{"mean"}, or \code{"none"}.
#'
#' @return
#' \code{recursive.triangle.mask.surface.embedding()} returns an \code{n x 3}
#' numeric matrix with columns \code{x}, \code{y}, and \code{z}.
#'
#' \code{recursive.triangle.mask.surface.graph()} returns a list with
#' components:
#' \itemize{
#'   \item \code{edges}: the retained recursive triangle-mask edges,
#'   \item \code{n}: number of vertices,
#'   \item \code{edge_weights}: induced positive edge lengths,
#'   \item \code{coords_surface}: the 3D surface embedding,
#'   \item \code{coords_param}: the canonical 2D triangle coordinates,
#'   \item \code{weight_scale}: the normalization constant applied to the raw
#'     edge lengths,
#'   \item \code{family}: always \code{"recursive.triangle.mask"},
#'   \item \code{surface}: the chosen surface name,
#'   \item \code{level}: the recursion depth,
#'   \item \code{mask}: the logical keep-mask,
#'   \item \code{label}: a human-readable family label.
#' }
#'
#' @name recursive_triangle_mask_surface_helpers
NULL

#' @rdname recursive_triangle_mask_surface_helpers
#' @export
recursive.triangle.mask.surface.embedding <- function(mask = mask.triangle.classic(),
                                                      level = 2,
                                                      surface = c("flat", "saddle",
                                                                  "paraboloid", "ripple",
                                                                  "folded"),
                                                      amplitude = 0.75,
                                                      freq_u = 1,
                                                      freq_v = 1,
                                                      x_scale = 1,
                                                      y_scale = 1) {
  mask <- .as_triangle_keep_mask(mask, "mask")
  level <- .as_whole_number(level, "level")
  surface <- match.arg(surface)
  amplitude <- .as_finite_scalar(amplitude, "amplitude")
  freq_u <- .as_positive_scalar(freq_u, "freq_u")
  freq_v <- .as_positive_scalar(freq_v, "freq_v")
  x_scale <- .as_positive_scalar(x_scale, "x_scale")
  y_scale <- .as_positive_scalar(y_scale, "y_scale")

  built <- .recursive.triangle.mask.canonical(mask, level)
  coords_param <- cbind(
    u = built$coords[, 1L] * x_scale,
    v = built$coords[, 2L] * y_scale
  )
  coords_norm <- .normalize.center_radius.coords(coords_param)
  z <- switch(
    surface,
    flat = rep(0, nrow(coords_param)),
    saddle = .surface.z.from_uv(coords_norm[, 1L], coords_norm[, 2L],
                                "saddle", amplitude, freq_u, freq_v),
    paraboloid = .surface.z.from_uv(coords_norm[, 1L], coords_norm[, 2L],
                                    "paraboloid", amplitude, freq_u, freq_v),
    ripple = .surface.z.from_uv(coords_norm[, 1L], coords_norm[, 2L],
                                "ripple", amplitude, freq_u, freq_v),
    folded = amplitude * abs(coords_norm[, 1L])
  )

  coords <- cbind(x = coords_param[, 1L], y = coords_param[, 2L], z = z)
  storage.mode(coords) <- "double"
  coords
}

#' @rdname recursive_triangle_mask_surface_helpers
#' @export
recursive.triangle.mask.surface.graph <- function(mask = mask.triangle.classic(),
                                                  level = 2,
                                                  surface = c("flat", "saddle",
                                                              "paraboloid", "ripple",
                                                              "folded"),
                                                  amplitude = 0.75,
                                                  freq_u = 1,
                                                  freq_v = 1,
                                                  x_scale = 1,
                                                  y_scale = 1,
                                                  normalize = c("median", "mean", "none")) {
  mask <- .as_triangle_keep_mask(mask, "mask")
  level <- .as_whole_number(level, "level")
  surface <- match.arg(surface)
  normalize <- match.arg(normalize)

  built <- .recursive.triangle.mask.canonical(mask, level)
  coords_param <- cbind(
    u = built$coords[, 1L] * x_scale,
    v = built$coords[, 2L] * y_scale
  )
  coords_surface <- recursive.triangle.mask.surface.embedding(
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
    edges = built$edges,
    coords = coords_surface,
    normalize = normalize
  )

  out <- list(
    edges = built$edges,
    n = built$n,
    edge_weights = weights$edge_weights,
    coords_surface = coords_surface,
    coords_param = coords_param,
    weight_scale = weights$weight_scale,
    family = "recursive.triangle.mask",
    surface = surface,
    level = level,
    mask = mask,
    normalize = normalize,
    label = sprintf("%s %s level %d",
                    tools::toTitleCase(surface),
                    .triangle.mask.label(mask),
                    level)
  )
  class(out) <- c("grip_recursive_triangle_mask_surface_graph", "list")
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

#' @describeIn graph_generators Deterministically irregular tetrahedralized ball
#'   graph built from nested subdivided polyhedral shells connected by a
#'   layered prism-to-tetrahedra edge pattern.
#' @param base Base polyhedron for \code{edges.triangulated.polyhedron()},
#'   \code{edges.irregular.ball()}, and \code{edges.irregular.shell()}. One of
#'   \code{"tetrahedron"}, \code{"octahedron"}, or \code{"icosahedron"}.
#' @param level Surface subdivision depth used by \code{edges.irregular.ball()}
#'   and \code{edges.irregular.shell()}.
#' @param layers Number of non-center radial layers for
#'   \code{edges.irregular.ball()} and number of inner-to-outer layers for
#'   \code{edges.irregular.shell()}.
#' @param outer_radius Positive outer radius for \code{edges.irregular.ball()}
#'   and \code{edges.irregular.shell()}.
#' @param radial_irregularity Radial layer-spacing irregularity level for
#'   \code{edges.irregular.ball()} and \code{edges.irregular.shell()}.
#' @param layer_twist Finite z-axis twist applied across radial layers in
#'   \code{edges.irregular.ball()} and \code{edges.irregular.shell()}.
#' @export
edges.irregular.ball <- function(base = c("tetrahedron", "octahedron", "icosahedron"),
                                 level = 1,
                                 layers = 3,
                                 outer_radius = 1,
                                 radial_irregularity = 0.25,
                                 layer_twist = 0.35) {
  .irregular.ball.canonical(
    base = match.arg(base),
    level = level,
    layers = layers,
    outer_radius = outer_radius,
    radial_irregularity = radial_irregularity,
    layer_twist = layer_twist
  )$edges
}

#' @describeIn graph_generators Deterministically irregular tetrahedralized shell
#'   graph built from nested subdivided polyhedral shells connected by a
#'   layered prism-to-tetrahedra edge pattern.
#' @param inner_radius Positive inner radius for \code{edges.irregular.shell()}.
#' @export
edges.irregular.shell <- function(base = c("tetrahedron", "octahedron", "icosahedron"),
                                  level = 1,
                                  layers = 3,
                                  inner_radius = 0.45,
                                  outer_radius = 1,
                                  radial_irregularity = 0.25,
                                  layer_twist = 0.35) {
  .irregular.shell.canonical(
    base = match.arg(base),
    level = level,
    layers = layers,
    inner_radius = inner_radius,
    outer_radius = outer_radius,
    radial_irregularity = radial_irregularity,
    layer_twist = layer_twist
  )$edges
}

#' @describeIn graph_generators Deterministically irregular torus graph built
#'   from major-cycle rings with varying sample counts and stitched into a
#'   locally triangulated closed surface.
#' @param major_rings Number of cyclic major rings for
#'   \code{edges.irregular.torus()}.
#' @param tube_count Approximate number of vertices around each minor cycle for
#'   \code{edges.irregular.torus()} and around each tube-like loop for
#'   \code{edges.irregular.double.torus()}.
#' @param major_irregularity Major-angle ring-spacing irregularity level for
#'   \code{edges.irregular.torus()}.
#' @export
edges.irregular.torus <- function(major_rings = 8,
                                  tube_count = 16,
                                  count_irregularity = 0.2,
                                  major_irregularity = 0.25,
                                  phase_twist = 0.35) {
  .irregular.torus.canonical(
    major_rings = major_rings,
    tube_count = tube_count,
    count_irregularity = count_irregularity,
    major_irregularity = major_irregularity,
    phase_twist = phase_twist
  )$edges
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

#' @describeIn graph_generators Deterministically irregular annulus graph built
#'   from concentric sample rings with varying sample counts and stitched into a
#'   locally triangulated surface-with-boundary graph.
#' @param rings Number of concentric sample rings for
#'   \code{edges.irregular.annulus()}.
#' @param outer_count Approximate number of vertices on the outer boundary for
#'   \code{edges.irregular.annulus()} and across the widest slices for
#'   \code{edges.irregular.pair.of.pants()}.
#' @param outer_radius Positive outer annulus radius for
#'   \code{edges.irregular.annulus()} and positive outer boundary radius for
#'   \code{edges.irregular.pair.of.pants()}.
#' @param inner_radius Positive inner annulus radius for
#'   \code{edges.irregular.annulus()}.
#' @param equator_count Approximate number of vertices near the equator for
#'   \code{edges.irregular.sphere()}.
#' @param count_irregularity Irregularity level for sample counts in
#'   \code{edges.irregular.torus()}, \code{edges.irregular.annulus()},
#'   \code{edges.irregular.pair.of.pants()},
#'   \code{edges.irregular.double.torus()}, and \code{edges.irregular.sphere()}.
#' @param radial_irregularity Within-ring radial irregularity level for
#'   \code{edges.irregular.annulus()}.
#' @param phase_twist Angular phase offset used to desynchronize neighboring
#'   rings, slice samples, or latitude bands in the irregular torus,
#'   irregular annulus, irregular pair-of-pants, irregular double torus,
#'   and irregular sphere families.
#' @param bands Number of non-pole latitude bands for
#'   \code{edges.irregular.sphere()}.
#' @param lat_irregularity Latitude-band spacing irregularity level for
#'   \code{edges.irregular.sphere()}.
#' @param slices Number of horizontal sample slices for
#'   \code{edges.irregular.pair.of.pants()}.
#' @param hole_radius Positive radius of each interior hole for
#'   \code{edges.irregular.pair.of.pants()}.
#' @param hole_offset Positive horizontal offset of the two hole centers for
#'   \code{edges.irregular.pair.of.pants()}.
#' @param hole_height Shared vertical coordinate of the two hole centers for
#'   \code{edges.irregular.pair.of.pants()}.
#' @param vertical_irregularity Slice-spacing irregularity level for
#'   \code{edges.irregular.pair.of.pants()}.
#' @param branch_length Half-length of the three-loop central region for
#'   \code{edges.irregular.double.torus()}.
#' @param branch_offset Offset of the outer loop centers from the middle loop
#'   for \code{edges.irregular.double.torus()}.
#' @param tube_radius Baseline radius of each tube-like loop for
#'   \code{edges.irregular.double.torus()}.
#' @param transition_width Width of the single-loop to three-loop transition
#'   regions for \code{edges.irregular.double.torus()}.
#' @param axial_irregularity Slice-spacing irregularity level for
#'   \code{edges.irregular.double.torus()}.
#' @export
edges.irregular.annulus <- function(rings = 6,
                                    outer_count = 28,
                                    outer_radius = 1,
                                    inner_radius = 0.45,
                                    count_irregularity = 0.2,
                                    radial_irregularity = 0.35,
                                    phase_twist = 0.35) {
  .irregular.annulus.canonical(
    rings = rings,
    outer_count = outer_count,
    outer_radius = outer_radius,
    inner_radius = inner_radius,
    count_irregularity = count_irregularity,
    radial_irregularity = radial_irregularity,
    phase_twist = phase_twist
  )$edges
}

#' @describeIn graph_generators Deterministically irregular pair-of-pants graph
#'   built from horizontal slice intervals with varying sample counts and
#'   stitched into a locally triangulated surface-with-boundary graph.
#' @export
edges.irregular.pair.of.pants <- function(slices = 11,
                                          outer_count = 28,
                                          outer_radius = 1.1,
                                          hole_radius = 0.24,
                                          hole_offset = 0.38,
                                          hole_height = 0.18,
                                          count_irregularity = 0.2,
                                          vertical_irregularity = 0.35,
                                          phase_twist = 0.35) {
  .irregular.pair.of.pants.canonical(
    slices = slices,
    outer_count = outer_count,
    outer_radius = outer_radius,
    hole_radius = hole_radius,
    hole_offset = hole_offset,
    hole_height = hole_height,
    count_irregularity = count_irregularity,
    vertical_irregularity = vertical_irregularity,
    phase_twist = phase_twist
  )$edges
}

#' @describeIn graph_generators Deterministically irregular double-torus graph
#'   built from cyclic slices whose cross-sections follow a `1 -> 3 -> 1` loop
#'   transition between two poles.
#' @export
edges.irregular.double.torus <- function(slices = 11,
                                         tube_count = 14,
                                         branch_length = 0.85,
                                         branch_offset = 0.72,
                                         tube_radius = 0.28,
                                         transition_width = 0.42,
                                         count_irregularity = 0.2,
                                         axial_irregularity = 0.3,
                                         phase_twist = 0.35) {
  .irregular.double.torus.canonical(
    slices = slices,
    tube_count = tube_count,
    branch_length = branch_length,
    branch_offset = branch_offset,
    tube_radius = tube_radius,
    transition_width = transition_width,
    count_irregularity = count_irregularity,
    axial_irregularity = axial_irregularity,
    phase_twist = phase_twist
  )$edges
}

#' @describeIn graph_generators Deterministically irregular sphere graph built
#'   from latitude bands with varying sample counts and stitched into a locally
#'   triangulated closed surface.
#' @export
edges.irregular.sphere <- function(bands = 6,
                                   equator_count = 28,
                                   count_irregularity = 0.2,
                                   lat_irregularity = 0.35,
                                   phase_twist = 0.35) {
  .irregular.sphere.canonical(
    bands = bands,
    equator_count = equator_count,
    count_irregularity = count_irregularity,
    lat_irregularity = lat_irregularity,
    phase_twist = phase_twist
  )$edges
}

#' @describeIn graph_generators Cube surface graph on the boundary of a
#'   \code{side x side x side} lattice.
#' @param side Number of lattice points along each cube edge.
#' @param base Base polyhedron for \code{edges.triangulated.polyhedron()}. One
#'   of \code{"tetrahedron"}, \code{"octahedron"}, or \code{"icosahedron"}.
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
  .kary.tree.structure(k = k, depth = depth)$edges
}

#' Intrinsically weighted k-ary tree helpers
#'
#' Convenience helpers that keep the exact \code{\link{edges.kary.tree}()}
#' topology but assign edge lengths directly from combinatorial tree metadata
#' rather than from an ambient Euclidean embedding. Each edge weight is
#' constructed as
#' \deqn{\mathrm{base\_length} \times \mathrm{depth\_factor(child\ depth)} \times
#' \mathrm{branch\_factor(child\ slot)}}
#' so the user can control tapering with depth separately from child-slot
#' asymmetry.
#'
#' The returned object is intended as an intrinsic weighted-tree benchmark
#' family. It exposes the same edge matrix as \code{edges.kary.tree()} together
#' with edge weights, vertex depths, parent indices, and an edge table that
#' records the depth and branch slot associated with each tree edge.
#'
#' @param k Branching factor. Must be at least \code{1}.
#' @param depth Number of levels below the root. May be \code{0}.
#' @param base_length Positive global edge-length multiplier before
#'   normalization.
#' @param depth_rule Rule used to build the per-depth multipliers. One of
#'   \code{"geometric"}, \code{"constant"}, or \code{"custom"}.
#' @param depth_decay Positive decay factor used when
#'   \code{depth_rule = "geometric"}.
#' @param depth_factors Positive custom depth multipliers. Used only when
#'   \code{depth_rule = "custom"}. Must have length \code{1} or \code{depth}.
#' @param branch_rule Rule used to build the per-child-slot multipliers. One of
#'   \code{"linear"}, \code{"uniform"}, or \code{"custom"}.
#' @param branch_spread Non-negative spread used when
#'   \code{branch_rule = "linear"}.
#' @param branch_factors Positive custom branch multipliers. Used only when
#'   \code{branch_rule = "custom"}. Must have length \code{1} or \code{k}.
#' @param normalize Normalization applied to the raw intrinsic edge lengths.
#'   One of \code{"median"}, \code{"mean"}, or \code{"none"}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{edges}: the \code{k}-ary tree edges,
#'   \item \code{n}: number of vertices,
#'   \item \code{edge_weights}: normalized positive intrinsic edge lengths,
#'   \item \code{weight_scale}: normalization constant applied to the raw edge
#'     lengths,
#'   \item \code{vertex_depth}: depth of each vertex, with the root at depth
#'     \code{0},
#'   \item \code{parent}: parent vertex index for each vertex, with the root
#'     parent recorded as \code{0},
#'   \item \code{edge_table}: data frame containing \code{parent},
#'     \code{child}, \code{parent_depth}, \code{child_depth},
#'     \code{branch_index}, \code{raw_weight}, and normalized
#'     \code{edge_weight},
#'   \item \code{depth_factors}: the applied per-depth multipliers,
#'   \item \code{branch_factors}: the applied per-child-slot multipliers,
#'   \item \code{family}: always \code{"kary.tree.weighted"},
#'   \item \code{k}: the branching factor,
#'   \item \code{depth}: the requested tree depth,
#'   \item \code{label}: a human-readable family label.
#' }
#'
#' @name kary_tree_weighted_graph_helpers
NULL

#' @rdname kary_tree_weighted_graph_helpers
#' @export
kary.tree.weighted.graph <- function(
    k = 2,
    depth = 2,
    base_length = 1,
    depth_rule = c("geometric", "constant", "custom"),
    depth_decay = 0.85,
    depth_factors = NULL,
    branch_rule = c("linear", "uniform", "custom"),
    branch_spread = 0.3,
    branch_factors = NULL,
    normalize = c("median", "mean", "none")) {
  k <- .as_whole_number(k, "k", min = 1L)
  depth <- .as_whole_number(depth, "depth")
  structure <- .kary.tree.structure(k = k, depth = depth)
  base_length <- .as_positive_scalar(base_length, "base_length")
  normalize <- match.arg(normalize)
  depth_spec <- .kary.tree.depth.factors(
    depth = depth,
    depth_rule = depth_rule,
    depth_decay = depth_decay,
    depth_factors = depth_factors
  )
  branch_spec <- .kary.tree.branch.factors(
    k = k,
    branch_rule = branch_rule,
    branch_spread = branch_spread,
    branch_factors = branch_factors
  )

  raw_weights <- if (nrow(structure$edges) == 0L) {
    numeric(0L)
  } else {
    base_length *
      depth_spec$factors[structure$child_depth] *
      branch_spec$factors[structure$branch_index]
  }
  weights <- .normalize.edge.weights(raw_weights, normalize = normalize)

  edge_table <- data.frame(
    parent = structure$edge_parent,
    child = structure$edge_child,
    parent_depth = structure$parent_depth,
    child_depth = structure$child_depth,
    branch_index = structure$branch_index,
    raw_weight = as.double(raw_weights),
    edge_weight = weights$edge_weights,
    stringsAsFactors = FALSE
  )

  out <- list(
    edges = structure$edges,
    n = structure$n,
    edge_weights = weights$edge_weights,
    weight_scale = weights$weight_scale,
    vertex_depth = structure$vertex_depth,
    parent = structure$parent,
    edge_table = edge_table,
    depth_factors = depth_spec$factors,
    branch_factors = branch_spec$factors,
    family = "kary.tree.weighted",
    k = k,
    depth = depth,
    base_length = as.double(base_length),
    depth_rule = depth_spec$rule,
    branch_rule = branch_spec$rule,
    normalize = normalize,
    label = sprintf("%s/%s intrinsic %d-ary tree depth %d",
                    tools::toTitleCase(depth_spec$rule),
                    branch_spec$rule,
                    k,
                    depth)
  )
  class(out) <- c("grip_kary_tree_weighted_graph", "list")
  out
}

#' @describeIn graph_generators Recursively refined square-mask grid graph whose
#'   vertices are occupied cells and whose edges connect orthogonally adjacent
#'   cells.
#' @param mask Keep-mask describing which recursive cells are retained. For
#'   \code{edges.recursive.mask.grid()}, \code{mask} must be a square logical
#'   or numeric keep-matrix whose non-zero entries are kept at each recursive
#'   subdivision step. For
#'   \code{edges.recursive.triangle.mask()}, \code{mask} must instead be a
#'   four-entry vector in \code{left}, \code{right}, \code{top},
#'   \code{center} order, and for \code{edges.recursive.tetrahedron.mask()},
#'   \code{mask} must be a four-entry vector in \code{base_left},
#'   \code{base_right}, \code{base_back}, \code{apex} order. For
#'   \code{edges.recursive.cube.mask()}, \code{mask} must be a cubic logical
#'   or numeric keep-array whose non-zero entries are kept at each recursive
#'   subdivision step.
#' @param level Recursion depth. For \code{edges.recursive.mask.grid()},
#'   \code{edges.recursive.cube.mask()}, \code{edges.vicsek()},
#'   \code{edges.menger.sponge()}, and \code{edges.sierpinski.carpet()},
#'   \code{level} must be at least 1; \code{edges.recursive.triangle.mask()},
#'   \code{edges.recursive.tetrahedron.mask()}, and
#'   \code{edges.triangulated.polyhedron()} also allow \code{level = 0}.
#' @export
edges.recursive.mask.grid <- function(mask, level = 2) {
  .recursive.mask.grid.edges(mask, level)
}

#' @describeIn graph_generators Recursively refined triangle-mask graph whose
#'   vertices are the retained subdivision vertices of an equilateral triangle.
#' @export
edges.recursive.triangle.mask <- function(mask = mask.triangle.classic(), level = 2) {
  .recursive.triangle.mask.canonical(mask, level)$edges
}

#' @describeIn graph_generators Recursively refined tetrahedron-mask graph whose
#'   vertices are the retained subdivision vertices of a tetrahedron.
#' @export
edges.recursive.tetrahedron.mask <- function(mask = mask.tetrahedron.classic(),
                                             level = 2) {
  .recursive.tetrahedron.mask.canonical(mask, level)$edges
}

#' @describeIn graph_generators Recursively refined cube-mask graph whose
#'   vertices are occupied subcubes and whose edges connect face-adjacent
#'   occupied cells.
#' @export
edges.recursive.cube.mask <- function(mask, level = 2) {
  .recursive.cube.mask.edges(mask, level)
}

#' @describeIn graph_generators Triangulated closed-surface graph obtained by
#'   repeatedly splitting the triangular faces of a tetrahedron, octahedron, or
#'   icosahedron.
#' @export
edges.triangulated.polyhedron <- function(
    base = c("tetrahedron", "octahedron", "icosahedron"),
    level = 1) {
  base <- match.arg(base)
  .triangulated.polyhedron.canonical(base = base, level = level)$edges
}

#' @describeIn graph_generators Triangulated annulus graph obtained by clipping
#'   a regular triangular lattice to the region between two concentric circles.
#' @param resolution Positive lattice-resolution control used by
#'   \code{edges.triangulated.annulus()} and
#'   \code{edges.triangulated.pair.of.pants()}.
#' @param outer_radius Positive outer boundary radius for
#'   \code{edges.triangulated.annulus()} and
#'   \code{edges.triangulated.pair.of.pants()}.
#' @param inner_radius Positive inner annulus radius for
#'   \code{edges.triangulated.annulus()}.
#' @param hole_radius Positive radius of each interior hole for
#'   \code{edges.triangulated.pair.of.pants()}.
#' @param hole_offset Positive horizontal offset of the two hole centers for
#'   \code{edges.triangulated.pair.of.pants()}.
#' @param hole_height Shared vertical coordinate of the two hole centers for
#'   \code{edges.triangulated.pair.of.pants()}.
#' @export
edges.triangulated.annulus <- function(resolution = 12,
                                       outer_radius = 1,
                                       inner_radius = 0.45) {
  .triangulated.annulus.canonical(
    resolution = resolution,
    outer_radius = outer_radius,
    inner_radius = inner_radius
  )$edges
}

#' @describeIn graph_generators Triangulated pair-of-pants graph obtained by
#'   clipping a regular triangular lattice to a disk with two interior circular
#'   holes.
#' @export
edges.triangulated.pair.of.pants <- function(resolution = 12,
                                             outer_radius = 1.1,
                                             hole_radius = 0.24,
                                             hole_offset = 0.38,
                                             hole_height = 0.18) {
  .triangulated.pair.of.pants.canonical(
    resolution = resolution,
    outer_radius = outer_radius,
    hole_radius = hole_radius,
    hole_offset = hole_offset,
    hole_height = hole_height
  )$edges
}

#' @describeIn graph_generators Connected Vicsek-style cross family derived from
#'   a \code{3 x 3} axial-cross keep-mask.
#' @export
edges.vicsek <- function(level = 2) {
  edges.recursive.mask.grid(.vicsek.mask(), level = level)
}

#' @describeIn graph_generators Classic Menger-sponge cubical cell-adjacency
#'   graph derived from the \code{3 x 3 x 3} keep-mask that removes the center
#'   cube and the six face-center cubes at each recursion step.
#' @export
edges.menger.sponge <- function(level = 2) {
  edges.recursive.cube.mask(.menger.sponge.mask(), level = level)
}

#' @describeIn graph_generators Periodic cubical tunnel family derived from a
#'   repeated tunnel-band keep-mask. The classic Menger sponge appears as the
#'   \code{side = 3}, \code{tunnel_width = 1} special case.
#' @export
edges.cube.periodic.tunnels <- function(level = 2,
                                        side = 5,
                                        tunnel_width = 1,
                                        tunnel_period = 2,
                                        tunnel_offset = 2) {
  edges.recursive.cube.mask(
    mask.cube.periodic.tunnels(
      side = side,
      tunnel_width = tunnel_width,
      tunnel_period = tunnel_period,
      tunnel_offset = tunnel_offset
    ),
    level = level
  )
}

#' @describeIn graph_generators Cubical porous family with two offset interior
#'   cavity blocks repeated recursively.
#' @export
edges.cube.asymmetric.cavities <- function(level = 2,
                                           side = 5,
                                           cavity_size = 2,
                                           pocket_size = max(1L, cavity_size - 1L)) {
  edges.recursive.cube.mask(
    mask.cube.asymmetric.cavities(
      side = side,
      cavity_size = cavity_size,
      pocket_size = pocket_size
    ),
    level = level
  )
}

#' @describeIn graph_generators Cubical porous family with a deterministic
#'   branched channel network carved through each recursive block.
#' @export
edges.cube.channel.network <- function(level = 2,
                                       side = 5,
                                       channel_width = 1,
                                       branch_offset = 2) {
  edges.recursive.cube.mask(
    mask.cube.channel.network(
      side = side,
      channel_width = channel_width,
      branch_offset = branch_offset
    ),
    level = level
  )
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
  .sierpinski.triangle.canonical(level)$edges
}

#' @describeIn graph_generators Three-dimensional tetrahedral Sierpinski graph
#'   at recursion depth \code{level}.
#' @export
edges.sierpinski.tetrahedron <- function(level = 2) {
  .sierpinski.tetrahedron.canonical(level)$edges
}

#' @describeIn graph_generators Two-dimensional Sierpinski carpet graph whose
#'   vertices are occupied cells and whose edges connect orthogonally adjacent
#'   cells.
#' @export
edges.sierpinski.carpet <- function(level = 2) {
  edges.recursive.mask.grid(.sierpinski.carpet.mask(), level = level)
}
