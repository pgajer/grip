test_that('catalogue and all bundled graphs agree, including disconnected graphs', {
  catalogue <- suitesparse.graphs()
  expect_equal(nrow(catalogue), 85L)
  expect_equal(sum(catalogue$gallery), 80L)
  expect_equal(sum(catalogue$original_six), 6L)
  expect_false(anyDuplicated(catalogue$id) > 0L)
  for (i in seq_len(nrow(catalogue))) {
    row <- catalogue[i, ]; g <- suitesparse.graph(row$id)
    expect_identical(g, suitesparse.graph(row$name))
    expect_identical(g$n, row$vertices)
    expect_identical(nrow(g$edges), row$edges)
    expect_type(g$edges, 'integer')
    expect_true(all(g$edges >= 1L & g$edges <= g$n))
    expect_true(all(g$edges[, 1] < g$edges[, 2]))
    expect_identical(order(g$edges[, 1], g$edges[, 2]), seq_len(nrow(g$edges)))
    expect_equal(anyDuplicated(data.frame(g$edges)), 0L)
    expect_identical(g$edge_weights, rep(1, nrow(g$edges)))
    expect_identical(g$vertex_data$vertex, seq_len(g$n))
    expect_length(g$vertex_data$label, g$n)
    expect_identical(g$provenance$graph_id, row$id)
    expect_equal(g$n / row$original_vertices, row$retained_fraction)
    # Independent union-find check, counting isolated vertices too.
    parent <- seq_len(g$n)
    root <- function(v) { while (parent[v] != v) v <- parent[v]; v }
    for (j in seq_len(nrow(g$edges))) parent[root(g$edges[j, 2])] <- root(g$edges[j, 1])
    expect_equal(length(unique(vapply(seq_len(g$n), root, 1L))), row$components)
    expect_true(length(g$provenance$source_header) > 0)
    expect_false(any(grepl('/Users/|/home/', unlist(g$provenance), fixed = FALSE)))
  }
  expect_gt(suitesparse.graphs()$components[suitesparse.graphs()$name == 'lock1074'], 1)
  g <- suitesparse.graph('lock1074')
  expect_gt(sum(tabulate(g$edges, nbins=g$n) == 0), 0)
})

test_that('bipartite and component vertex identities remain distinct', {
  g <- suitesparse.graph('illc1033')
  expect_equal(g$n, 1353L)
  expect_equal(table(g$vertex_data$source_partition)[c('row','column')], c(row=1033L,column=320L), ignore_attr=TRUE)
  expect_identical(g$vertex_data$source_index, c(seq_len(1033L),seq_len(320L)))
  expect_true(all(g$edges[,1] <= 1033L & g$edges[,2] > 1033L))
  g <- suitesparse.graph('lock_700__LCC_691_of_700')
  expect_identical(g$provenance$original_graph_id, 'HB/lock_700')
  expect_equal(g$provenance$original_vertices, 700)
  expect_identical(g$vertex_data$source_index, g$vertex_data$original_vertex)
  expect_true(all(diff(g$vertex_data$original_vertex) > 0))
  expect_gt(max(g$vertex_data$original_vertex), g$n)
  expect_error(suitesparse.graph('HB/lock_700'), 'Unknown bundled')
})

test_that('lookup errors provide usable candidates and do not guess ambiguity', {
  for (name in list(NULL, NA_character_, '', 1L, c('494_bus', 'lesmis'))) {
    expect_error(suitesparse.graph(name), 'one non-empty')
  }
  expect_error(suitesparse.graph('494_bu'), 'HB/494_bus')
  catalogue <- data.frame(id=c('A/example','B/example'), name=c('example','example'))
  expect_error(grip:::.suitesparse_resolve('example', catalogue), 'A/example, B/example')
  expect_identical(grip:::.suitesparse_resolve('A/example', catalogue), 'A/example')
  # Changes by a caller must not mutate cached package objects.
  g <- suitesparse.graph('dwt_307'); g$edges[1,1] <- -1L
  expect_gt(suitesparse.graph('dwt_307')$edges[1,1], 0L)
})
