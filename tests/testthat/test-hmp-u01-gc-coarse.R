test_that("hmp.u01.gc.coarse bundled object is internally consistent", {
  env <- new.env(parent = emptyenv())
  data("hmp.u01.gc.coarse", package = "grip", envir = env)
  obj <- get("hmp.u01.gc.coarse", envir = env, inherits = FALSE)

  expect_type(obj, "list")
  expect_true(all(c("adj_list", "weight_list", "vertex_data", "graph_info") %in% names(obj)))

  n <- length(obj$adj_list)
  expect_equal(n, 1828L)
  expect_equal(length(obj$weight_list), n)
  expect_equal(nrow(obj$vertex_data), n)

  degs <- vapply(obj$adj_list, length, integer(1))
  wlen <- vapply(obj$weight_list, length, integer(1))
  expect_equal(degs, wlen)

  info <- obj$graph_info
  expect_equal(info$selected_k, 3L)
  expect_equal(info$original_vertices, 6474L)
  expect_equal(info$coarse_vertices, 1828L)
  expect_equal(info$edge_count, 4656L)
})

test_that("bundled HMP/U01 vignette results are available and well formed", {
  path <- system.file(
    "extdata", "hmp_u01_gc_coarse", "vignette_results.rds",
    package = "grip"
  )
  expect_true(nzchar(path))
  expect_true(file.exists(path))

  res <- readRDS(path)
  expect_true(all(c("preset_summary", "local_search_summary", "layouts") %in% names(res)))
  expect_true(all(c("default", "tree", "torus") %in% names(res$layouts$preset)))
  expect_true(nrow(res$preset_summary) >= 3L)
  expect_true(nrow(res$local_search_summary) >= 2L)
})
