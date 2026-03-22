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
