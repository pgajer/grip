misf_level_min_hop <- function(levels, edges, n) {
  built <- grip:::grip.build.adj.from.edges(edges, n = n)
  adj <- built$adj_list

  bfs_dist <- function(src) {
    dist <- rep.int(Inf, n)
    queue <- integer(n)
    head <- 1L
    tail <- 1L
    queue[[tail]] <- src
    dist[[src]] <- 0
    while (head <= tail) {
      vert <- queue[[head]]
      head <- head + 1L
      for (nbr in adj[[vert]]) {
        if (!is.finite(dist[[nbr]])) {
          tail <- tail + 1L
          queue[[tail]] <- nbr
          dist[[nbr]] <- dist[[vert]] + 1L
        }
      }
    }
    dist
  }

  all.dist <- lapply(seq_len(n), bfs_dist)
  vapply(seq_along(levels), function(idx) {
    verts <- levels[[idx]]
    if (length(verts) < 2L) {
      return(Inf)
    }
    min.hop <- Inf
    for (i in seq_len(length(verts) - 1L)) {
      min.hop <- min(min.hop, min(all.dist[[verts[[i]]]][verts[(i + 1L):length(verts)]]))
    }
    min.hop
  }, numeric(1L))
}

test_that("grip.build.misf returns nested deterministic levels on a mesh", {
  edges <- edges.mesh(4, 4)

  misf1 <- grip.build.misf(edges = edges, n = 16L, num_init = 6L, num_nbrs = 8L, seed = 11L)
  misf2 <- grip.build.misf(edges = edges, n = 16L, num_init = 6L, num_nbrs = 8L, seed = 11L)

  expect_s3_class(misf1, "grip_misf")
  expect_identical(misf1$levels, misf2$levels)
  expect_identical(misf1$mish_order, misf2$mish_order)
  expect_identical(misf1$vertex_depth, misf2$vertex_depth)
  expect_identical(misf1$misf_size, misf2$misf_size)
  expect_identical(misf1$num_nbrs_schedule, misf2$num_nbrs_schedule)

  expect_equal(sort(misf1$levels[[1L]]), seq_len(16L))
  expect_equal(length(misf1$levels), misf1$misf_height + 1L)
  expect_equal(unname(vapply(misf1$levels, length, integer(1))), misf1$misf_size)

  for (i in seq_len(length(misf1$levels) - 1L)) {
    expect_true(all(misf1$levels[[i + 1L]] %in% misf1$levels[[i]]), info = i)
  }

  expected.depth <- integer(misf1$n)
  for (i in seq_along(misf1$levels)) {
    expected.depth[misf1$levels[[i]]] <- i - 1L
  }
  expect_identical(misf1$vertex_depth, expected.depth)
  expect_equal(sort(misf1$mish_order), seq_len(16L))
  expect_lte(tail(misf1$misf_size, 1L), 6L)
})

test_that("grip.build.misf preserves the MIS hop bound on every returned level", {
  edges <- edges.mesh(12L, 12L)
  misf <- grip.build.misf(
    edges = edges,
    n = 144L,
    num_init = 24L,
    num_nbrs = 20L,
    seed = 12L
  )

  min.hops <- misf_level_min_hop(misf$levels, edges, 144L)
  required <- c(1, vapply(seq.int(1L, length(misf$levels) - 1L), function(level) {
    as.integer(2^(level - 1L) + 1L)
  }, integer(1L)))

  expect_true(all(min.hops >= required))
})

test_that("grip.build.misf handles small graphs with a single MISF level", {
  edges <- edges.path(5L)

  misf <- grip.build.misf(edges = edges, n = 5L, num_init = 8L, num_nbrs = 4L, seed = 3L)

  expect_equal(length(misf$levels), 1L)
  expect_equal(misf$misf_height, 0L)
  expect_equal(misf$top_level_size, 5L)
  expect_equal(sort(misf$levels[[1L]]), seq_len(5L))
  expect_equal(misf$vertex_depth, integer(5L))
  expect_equal(misf$misf_size, 5L)
})

test_that("grip.build.misf accepts weighted graph input and matches adj-list mode", {
  edges <- matrix(c(
    1L, 2L,
    2L, 3L,
    3L, 4L,
    4L, 5L,
    2L, 5L
  ), byrow = TRUE, ncol = 2L)
  edge.weights <- c(1, 2, 3, 4, 5)
  built <- grip:::grip.build.adj.from.edges(edges, n = 5L, edge_weights = edge.weights)

  misf.edges <- grip.build.misf(
    edges = edges,
    n = 5L,
    edge_weights = edge.weights,
    num_init = 3L,
    num_nbrs = 4L,
    seed = 7L
  )
  misf.adj <- grip.build.misf(
    adj_list = built$adj_list,
    weight_list = built$weight_list,
    n = 5L,
    num_init = 3L,
    num_nbrs = 4L,
    seed = 7L
  )
  misf.unweighted <- grip.build.misf(
    edges = edges,
    n = 5L,
    num_init = 3L,
    num_nbrs = 4L,
    seed = 7L
  )

  expect_identical(misf.edges$levels, misf.adj$levels)
  expect_identical(misf.edges$vertex_depth, misf.adj$vertex_depth)
  expect_identical(misf.edges$mish_order, misf.adj$mish_order)
  expect_identical(misf.edges$levels, misf.unweighted$levels)
})
