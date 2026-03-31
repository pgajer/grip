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
#' The Sierpinski families are exposed explicitly rather than overloading a
#' single generator with layout-dimension-dependent behavior:
#' \code{\link{edges.sierpinski.triangle}()} builds the 2-simplex family,
#' \code{\link{edges.sierpinski.tetrahedron}()} builds the 3-simplex family,
#' and \code{\link{edges.sierpinski.carpet}()} builds a 2D cell-adjacency
#' carpet graph.
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

#' @describeIn graph_generators Two-dimensional Sierpinski triangle graph at
#'   recursion depth \code{level}.
#' @param level Recursion depth. For \code{edges.sierpinski.carpet()},
#'   \code{level} must be at least 1.
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
  level <- .as_whole_number(level, "level", min = 1L)
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

  .bind_edges(edges)
}
