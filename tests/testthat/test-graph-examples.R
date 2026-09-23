test_that('SuiteSparse graph examples preserve source patterns and vertex identities', {
  env <- new.env(); data('zheng.graphs', package = 'grip', envir = env)
  expected <- c(dwt_66 = 66, lesmis = 77, dwt_307 = 307, '494_bus' = 494, dwt_1005 = 1005, '1138_bus' = 1138)
  edge.counts <- c(127,254,1108,586,3808,1458)
  expect_identical(names(env$zheng.graphs), names(expected))
  for (i in seq_along(expected)) {
    id <- names(expected)[i]; g <- env$zheng.graphs[[id]]
    expect_equal(g$n, unname(expected[i])); expect_equal(nrow(g$edges), edge.counts[i])
    expect_identical(g$vertex_data$source_index, seq_len(g$n))
    expect_true(all(g$edge_weights == 1))
    con <- gzfile(system.file('extdata','zheng-graphs',paste0(id,'.mtx.gz'),package='grip'))
    lines <- readLines(con); close(con)
    values <- lines[!grepl('^%',lines)]
    entries <- as.matrix(read.table(text = paste(values[-1],collapse='\n')))
    if (ncol(entries) == 3L) entries <- entries[entries[,3] != 0,,drop=FALSE]
    entries <- entries[entries[,1] != entries[,2],1:2,drop=FALSE]
    E <- unique(cbind(pmin(entries[,1],entries[,2]),pmax(entries[,1],entries[,2])))
    E <- E[order(E[,1],E[,2]),,drop=FALSE]; storage.mode(E) <- 'integer'
    expect_equal(g$edges, unname(E))
    # Independent breadth-first traversal, without allocating all-pairs distances.
    seen <- 1L
    repeat {
      next.seen <- sort(unique(c(seen, g$edges[g$edges[,1] %in% seen,2], g$edges[g$edges[,2] %in% seen,1])))
      if (identical(seen,next.seen)) break
      seen <- next.seen
    }
    expect_identical(seen, seq_len(g$n))
  }
  expect_identical(env$zheng.graphs$lesmis$vertex_data$label[1], 'Myriel')
})

test_that('gallery scores distinguish global distances from edge lengths', {
  env <- new.env(); sys.source(system.file('scripts','graph-examples.R',package='grip'), env)
  g <- list(edges = rbind(c(1,2),c(2,3)), edge_weights = c(1,1))
  D <- as.matrix(dist(0:2))
  Z <- rbind(c(0,0,0),c(1,0,0),c(1,1,0))
  score <- env$graph_example_scores(Z,g,D)
  expect_equal(score[['edge_error']],0)
  expect_equal(score[['graph_error']],sqrt((sqrt(2)-2)^2/6))
  aligned <- env$graph_example_align(2*Z, Z)
  expect_equal(as.vector(dist(aligned)),as.vector(dist(2*Z)),tolerance=1e-12)
})

test_that('Procrustes alignment recovers known rotations, reflections and translations', {
  env <- new.env(); sys.source(system.file('scripts','graph-examples.R',package='grip'), env)
  reference <- rbind(c(0,0,0), c(1,0,0), c(0,2,0), c(0,0,3), c(2,1,-1))
  theta <- .73
  rotation <- rbind(c(cos(theta),-sin(theta),0), c(sin(theta),cos(theta),0), c(0,0,1))
  for (reflection in c(1,-1)) {
    transformed <- sweep(reference %*% rotation %*% diag(c(1,1,reflection)),
                         2, c(7,-4,2), '+')
    aligned <- env$graph_example_align(transformed, reference)
    expect_equal(aligned, reference, tolerance = 1e-12)
    expect_equal(as.vector(dist(aligned)), as.vector(dist(transformed)), tolerance = 1e-12)
  }
})

test_that('gallery recipes preserve the saved graph field names', {
  env <- new.env()
  sys.source(system.file('scripts', 'graph-examples.R', package = 'grip'), env)
  for (g in env$graph_example_cases()) {
    expect_length(g$edge_weights, nrow(g$edges))
    expect_true(all(is.finite(g$edge_weights) & g$edge_weights > 0))
  }
})
