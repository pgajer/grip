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
      # wrap around in the width direction
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

edges_kary_tree <- function(k = 2, depth = 2) {
  stopifnot(k >= 1, depth >= 0)
  if (depth == 0) return(matrix(integer(), ncol = 2))
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
  do.call(rbind, edges)
}
