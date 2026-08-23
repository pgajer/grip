test_that("hmp.gc contains only HMP Illumina 16S samples", {
  env <- new.env(parent = emptyenv())
  data("hmp.gc", package = "grip", envir = env)
  obj <- get("hmp.gc", envir = env, inherits = FALSE)

  expect_type(obj, "list")
  expect_true(all(c("adj_list", "weight_list", "vertex_data", "graph_info") %in% names(obj)))

  n <- length(obj$adj_list)
  expect_equal(n, 4391L)
  expect_equal(length(obj$weight_list), n)
  expect_equal(nrow(obj$vertex_data), n)
  expect_true(all(grepl("^HMP", obj$vertex_data$project)))
  expect_true(all(obj$vertex_data$platform == "Illumina"))

  expect_equal(lengths(obj$adj_list), lengths(obj$weight_list))
  expect_true(all(vapply(seq_len(n), function(i) {
    all(vapply(obj$adj_list[[i]], function(j) i %in% obj$adj_list[[j]], logical(1L)))
  }, logical(1L))))

  info <- obj$graph_info
  expect_equal(info$selected_k, 3L)
  expect_equal(info$original_vertices, 4391L)
  expect_equal(info$edge_count, 9067L)
  expect_identical(info$assay, "Illumina 16S rRNA amplicon sequencing")
})
