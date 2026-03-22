## GRIP layout examples (R script version of the vignette)

library(grip)

edges_path <- function(n) {
  if (n < 2) return(matrix(integer(), ncol = 2))
  cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L)
}

edges_cycle <- function(n) {
  if (n < 2) return(matrix(integer(), ncol = 2))
  rbind(edges_path(n), c(n, 1L))
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

edges_cylinder <- function(h, w = h) {
  stopifnot(h >= 1, w >= 1)
  idx <- function(i, j) (i - 1L) * w + j
  edges <- list()
  for (i in seq_len(h)) {
    for (j in seq_len(w)) {
      v <- idx(i, j)
      if (i < h) edges[[length(edges) + 1L]] <- c(v, idx(i + 1L, j))
      edges[[length(edges) + 1L]] <- c(v, idx(i, (j %% w) + 1L))
    }
  }
  do.call(rbind, edges)
}

edges_torus <- function(h, w = h) {
  stopifnot(h >= 1, w >= 1)
  idx <- function(i, j) (i - 1L) * w + j
  edges <- list()
  for (i in seq_len(h)) {
    for (j in seq_len(w)) {
      v <- idx(i, j)
      edges[[length(edges) + 1L]] <- c(v, idx((i %% h) + 1L, j))
      edges[[length(edges) + 1L]] <- c(v, idx(i, (j %% w) + 1L))
    }
  }
  do.call(rbind, edges)
}

edges_sierpinski_triangle <- function(level = 2) {
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

## 2D examples
edges <- edges_path(12)
coords <- grip.layout(edges, n = 12, dim = 2,
                      engine = "mish_v5",
                      placement = "barycenter",
                      rounds = 5, final_rounds = 3,
                      num_init = 5, num_nbrs = 6,
                      seed = 1)
grip.plot(coords, edges)

edges <- edges_cycle(16)
coords <- grip.layout(edges, n = 16, dim = 2,
                      engine = "mish_v5",
                      placement = "barycenter",
                      rounds = 5, final_rounds = 3,
                      num_init = 5, num_nbrs = 6,
                      seed = 2)
grip.plot(coords, edges)

edges <- edges_mesh(5, 5)
coords <- grip.layout(edges, n = 25, dim = 2,
                      engine = "mish_v6",
                      placement = "barycenter",
                      rounds = 5, final_rounds = 3,
                      num_init = 6, num_nbrs = 8,
                      seed = 3)
grip.plot(coords, edges)

edges <- edges_sierpinski_triangle(2)
n <- max(edges)
coords <- grip.layout(edges, n = n, dim = 2,
                      engine = "mish_v5",
                      placement = "circle",
                      rounds = 5, final_rounds = 3,
                      num_init = 5, num_nbrs = 7,
                      seed = 4)
grip.plot(coords, edges)

## 3D examples
edges <- edges_mesh(4, 4)
coords <- grip.layout(edges, n = 16, dim = 3,
                      engine = "mish_v5",
                      placement = "barycenter",
                      rounds = 5, final_rounds = 3,
                      num_init = 5, num_nbrs = 7,
                      seed = 5)
head(coords)

edges <- edges_cylinder(4, 6)
coords <- grip.layout(edges, n = 24, dim = 3,
                      engine = "mish_v5",
                      placement = "barycenter",
                      rounds = 5, final_rounds = 3,
                      num_init = 6, num_nbrs = 8,
                      seed = 6)
head(coords)

edges <- edges_torus(4, 4)
coords <- grip.layout(edges, n = 16, dim = 3,
                      engine = "mish_v6",
                      placement = "barycenter",
                      rounds = 5, final_rounds = 3,
                      num_init = 5, num_nbrs = 7,
                      seed = 7)
head(coords)

## Adjacency list + weights (weight_list is optional; omit for all 1s)
adj_list <- list(c(2), c(1, 3), c(2, 4), c(3))
weight_list <- list(c(1.0), c(1.0, 2.0), c(2.0, 1.5), c(1.5))
coords <- grip.layout(adj_list = adj_list,
                      weight_list = weight_list,
                      n = 4,
                      dim = 2,
                      engine = "mish_v5",
                      placement = "barycenter",
                      rounds = 4, final_rounds = 2,
                      num_init = 3, num_nbrs = 3,
                      seed = 12)
grip.plot(coords)
