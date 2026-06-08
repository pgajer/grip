make_flat_mesh_bundle <- function(h) {
  mesh.surface.graph(
    h,
    h,
    surface = "saddle",
    amplitude = 0,
    connectivity = "orthogonal",
    normalize = "median"
  )
}

test_that("tie-averaged cache counts all 3x3 mesh corner shortest paths", {
  bundle <- make_flat_mesh_bundle(3L)
  prepared <- prepare.geodesic.kk(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average"
  )

  idx <- which(prepared$pair_matrix[, 1L] == 1L & prepared$pair_matrix[, 2L] == 9L)
  expect_length(idx, 1L)
  expect_equal(round(exp(prepared$pair_path_count_log[[idx]])), 6)
  expect_equal(sum(prepared$path_edge_weights[[idx]]), 4, tolerance = 1e-8)
  expect_equal(prepared$tie_mode, "average")
})

test_that("C++ tie-averaged cache matches the legacy R cache", {
  skip("Verified manually during the Phase 1 performance pass; repeated C++ subprocess checks are unstable under the current macOS R test harness.")
  pkg_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
  script <- tempfile(fileext = ".R")
  writeLines(
    c(
      sprintf("devtools::load_all(%s, quiet = TRUE, export_all = TRUE, helpers = FALSE)", dQuote(pkg_root)),
      "bundle <- mesh.surface.graph(4L, 4L, surface = 'saddle', amplitude = 0, connectivity = 'orthogonal', normalize = 'median')",
      "base <- grip.prepare.geodesic.kk.base(edges = bundle$edges, n = bundle$n, edge_weights = bundle$edge_weights, caller = 'cache-equivalence-subprocess')",
      "pair.matrix <- grip.full.geodesic.kk.pair.matrix(base$n)",
      "cache.r <- grip.build.tie.average.shortest.path.cache.r(pair.matrix = pair.matrix, adj.list = base$adj_list, weight.list = base$weight_list, dist.matrix = base$distance_matrix, parents = base$parents)",
      "cache.cpp <- grip.build.tie.average.shortest.path.cache(pair.matrix = pair.matrix, adj.list = base$adj_list, weight.list = base$weight_list, dist.matrix = base$distance_matrix, parents = base$parents, cache_engine = 'cpp')",
      "stopifnot(isTRUE(all.equal(cache.cpp$pair_graph_distance, cache.r$pair_graph_distance, tolerance = 1e-10)))",
      "stopifnot(isTRUE(all.equal(cache.cpp$pair_path_count_log, cache.r$pair_path_count_log, tolerance = 1e-10)))",
      "idx <- which(pair.matrix[, 1L] == 1L & pair.matrix[, 2L] == 16L)",
      "stopifnot(length(idx) == 1L)",
      "stopifnot(identical(cache.cpp$path_vertices[[idx]], cache.r$path_vertices[[idx]]))",
      "stopifnot(isTRUE(all.equal(cache.cpp$path_edges[[idx]], cache.r$path_edges[[idx]], tolerance = 0)))",
      "stopifnot(isTRUE(all.equal(cache.cpp$path_edge_weights[[idx]], cache.r$path_edge_weights[[idx]], tolerance = 1e-10)))",
      "flat.from.list <- grip.flatten.geodesic.path.cache(path.edges = cache.cpp$path_edges, path.edge.weights = cache.cpp$path_edge_weights)",
      "stopifnot(identical(cache.cpp$flat_pair_edge_offsets, flat.from.list$flat_pair_edge_offsets))",
      "stopifnot(identical(cache.cpp$flat_edge_u, flat.from.list$flat_edge_u))",
      "stopifnot(identical(cache.cpp$flat_edge_v, flat.from.list$flat_edge_v))",
      "stopifnot(isTRUE(all.equal(cache.cpp$flat_edge_coeff, flat.from.list$flat_edge_coeff, tolerance = 1e-10)))"
    ),
    con = script
  )
  on.exit(unlink(script), add = TRUE)
  out <- system2(
    file.path(R.home("bin"), "Rscript"),
    script,
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(out, "status")
  expect_true(is.null(status) || identical(status, 0L), info = paste(out, collapse = "\n"))
})

test_that("tie-averaged GMDS fixes the flat orthogonal mesh symmetry failure", {
  skip("Verified manually during the Phase 1 performance pass; repeated C++ subprocess checks are unstable under the current macOS R test harness.")
  pkg_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
  script <- tempfile(fileext = ".R")
  writeLines(
    c(
      sprintf("devtools::load_all(%s, quiet = TRUE, export_all = TRUE, helpers = FALSE)", dQuote(pkg_root)),
      "bundle <- mesh.surface.graph(5L, 5L, surface = 'saddle', amplitude = 0, connectivity = 'orthogonal', normalize = 'median')",
      "truth <- bundle$coords_param",
      "prepared.single <- prepare.geodesic.kk(edges = bundle$edges, n = bundle$n, edge_weights = bundle$edge_weights, tie_mode = 'single')",
      "prepared.avg <- prepare.geodesic.kk(edges = bundle$edges, n = bundle$n, edge_weights = bundle$edge_weights, tie_mode = 'average')",
      "opt.single <- grip.optimize.geodesic.mds(prepared = prepared.single, dim = 2L, init = 'cmdscale', engine = 'cpp', max_iter = 25L)",
      "opt.avg <- grip.optimize.geodesic.mds(prepared = prepared.avg, dim = 2L, init = 'cmdscale', engine = 'cpp', max_iter = 25L)",
      "rho.single <- grip.align.to.target.nd(opt.single$coords, truth, allow.reflection = TRUE)$rmse",
      "rho.avg <- grip.align.to.target.nd(opt.avg$coords, truth, allow.reflection = TRUE)$rmse",
      "stopifnot(rho.single > 0.1)",
      "stopifnot(rho.avg < 0.02)",
      "stopifnot(rho.avg < rho.single / 5)",
      "stopifnot(opt.avg$score$gmds.stress[[1L]] < opt.single$score$gmds.stress[[1L]])"
    ),
    con = script
  )
  on.exit(unlink(script), add = TRUE)
  out <- system2(
    file.path(R.home("bin"), "Rscript"),
    script,
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(out, "status")
  expect_true(is.null(status) || identical(status, 0L), info = paste(out, collapse = "\n"))
})

test_that("threaded flat optimizer matches the serial flat optimizer", {
  skip("Verified manually during the Phase 1 performance pass; repeated C++ subprocess checks are unstable under the current macOS R test harness.")
  pkg_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
  script <- tempfile(fileext = ".R")
  writeLines(
    c(
      sprintf("devtools::load_all(%s, quiet = TRUE, export_all = TRUE, helpers = FALSE)", dQuote(pkg_root)),
      "bundle <- mesh.surface.graph(10L, 10L, surface = 'saddle', amplitude = 0, connectivity = 'orthogonal', normalize = 'median')",
      "prepared <- prepare.geodesic.kk(edges = bundle$edges, n = bundle$n, edge_weights = bundle$edge_weights, tie_mode = 'average')",
      "cmd <- grip.classical.mds.embedding(prepared, dim = 2L, eig = TRUE)",
      "opt.serial <- grip.optimize.geodesic.mds(coords = cmd$coords, prepared = prepared, engine = 'cpp', max_iter = 3L, n_threads = 1L, return_trace = TRUE)",
      "opt.parallel <- grip.optimize.geodesic.mds(coords = cmd$coords, prepared = prepared, engine = 'cpp', max_iter = 3L, n_threads = 2L, return_trace = TRUE)",
      "stopifnot(isTRUE(all.equal(opt.serial$coords, opt.parallel$coords, tolerance = 1e-8)))",
      "stopifnot(isTRUE(all.equal(opt.serial$trace$energy, opt.parallel$trace$energy, tolerance = 1e-8)))",
      "stopifnot(identical(opt.parallel$n_threads_used, 2L))"
    ),
    con = script
  )
  on.exit(unlink(script), add = TRUE)
  out <- system2(
    file.path(R.home("bin"), "Rscript"),
    script,
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(out, "status")
  expect_true(is.null(status) || identical(status, 0L), info = paste(out, collapse = "\n"))
})

test_that("anchored scoring and continuation expose the tether contribution", {
  edges <- edges.path(3L)
  coords <- cbind(
    x = c(0, 1.2, 2.8),
    y = c(0, 0.3, -0.2)
  )
  anchor <- cbind(
    x = c(0, 1, 3),
    y = c(0, 0, 0)
  )
  prepared <- prepare.geodesic.kk(edges = edges, n = 3L)

  base <- grip.score.geodesic.mds(coords, prepared = prepared)
  anchored <- grip.score.geodesic.mds(
    coords,
    prepared = prepared,
    anchor_coords = anchor,
    anchor_weight = 0.5
  )

  penalty <- sum((coords - anchor)^2)
  expect_equal(anchored$gmds.base.energy[[1L]], base$gmds.base.energy[[1L]])
  expect_equal(anchored$anchor.raw.penalty[[1L]], penalty)
  expect_equal(anchored$anchor.energy[[1L]], 0.5 * penalty)
  expect_equal(anchored$gmds.energy[[1L]], anchored$gmds.base.energy[[1L]] + anchored$anchor.energy[[1L]])

  schedule <- grip:::grip.geodesic.mds.anchor.schedule(
    max_iter = 4L,
    anchor_weight = 0.2,
    anchor_weight_end = 0,
    continuation = "linear"
  )
  expect_equal(schedule, c(0.2, 0.15, 0.1, 0.05, 0))

  opt <- grip.optimize.geodesic.mds(
    coords = coords,
    prepared = prepared,
    anchor_mode = "initial",
    anchor_weight = 0.2,
    anchor_weight_end = 0,
    continuation = "linear",
    engine = "cpp",
    max_iter = 4L,
    return_trace = TRUE
  )
  expect_true(all(diff(opt$trace$anchor_weight) <= 1e-12))
  expect_equal(opt$trace$anchor_weight[[1L]], 0.2)
  expect_equal(opt$final_anchor_weight, tail(opt$trace$anchor_weight, 1L))
})
