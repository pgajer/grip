adj_from_edges <- function(edges, n) {
  out <- vector("list", n)
  for (i in seq_len(nrow(edges))) {
    u <- edges[i, 1L]
    v <- edges[i, 2L]
    out[[u]] <- c(out[[u]], v)
    out[[v]] <- c(out[[v]], u)
  }
  lapply(out, as.integer)
}

write_tsv <- function(x, path) {
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = TRUE,
    row.names = FALSE,
    col.names = TRUE
  )
}

local_test_dir <- function(prefix) {
  path <- file.path(
    tempdir(),
    paste0(prefix, "-", Sys.getpid(), "-", sample.int(1e6, 1L))
  )
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

test_that("gripui_project normalizes minimal layout catalogs", {
  graph <- list(
    adj_list = list(2L, c(1L, 3L), 2L),
    vertex_data = data.frame(vertex_id = c("v1", "v2", "v3"), stringsAsFactors = FALSE)
  )
  layouts <- data.frame(
    candidate = "toy.layout",
    stage = "layout",
    seed = 1L,
    status = "ok",
    stringsAsFactors = FALSE
  )

  project <- gripui_project(graph = graph, layouts = layouts, title = "Toy")

  expect_s3_class(project, "gripui_project")
  expect_equal(project$meta$title, "Toy")
  expect_true(all(c(
    "layout_id",
    "candidate",
    "stage",
    "seed",
    "status",
    "viewable",
    "availability"
  ) %in% names(project$layouts)))
  expect_false(project$layouts$viewable[[1L]])
  expect_equal(project$layouts$availability[[1L]], "summary-only")
})

test_that("gripui_project_from_compare attaches in-memory coordinates", {
  edges <- edges.path(6)
  cmp <- grip.compare.layouts(
    edges = edges,
    n = 6,
    candidates = c("default"),
    seeds = 1:2,
    return.layouts = TRUE
  )
  graph <- list(adj_list = adj_from_edges(edges, n = 6))

  project <- gripui_project_from_compare(cmp, graph = graph, title = "Path compare")

  expect_s3_class(project, "gripui_project")
  expect_equal(nrow(project$layouts), 2L)
  expect_true(all(project$layouts$viewable))
  expect_true(all(vapply(project$layouts$coords, is.matrix, logical(1L))))

  coords <- gripui_load_layout_coords(project, project$layouts$layout_id[[1L]])
  expect_true(is.matrix(coords))
  expect_equal(nrow(coords), 6L)
})

test_that("gripui_project_from_dir loads a normalized bundle", {
  root <- local_test_dir("gripui-bundle")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  dir.create(file.path(root, "catalog"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "graph"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "artifacts", "layouts"), recursive = TRUE, showWarnings = FALSE)

  graph <- list(
    adj_list = list(2L, c(1L, 3L), 2L),
    vertex_data = data.frame(vertex_id = c("v1", "v2", "v3"), stringsAsFactors = FALSE)
  )
  saveRDS(graph, file.path(root, "graph", "graph.rds"))

  coords.path <- file.path(root, "artifacts", "layouts", "toy_embedding.tsv")
  write_tsv(
    data.frame(
      vertex_id = c("v2", "v1", "v3"),
      x = c(2, 1, 3),
      y = c(0, 0, 0),
      stringsAsFactors = FALSE
    ),
    coords.path
  )

  write_tsv(
    data.frame(
      candidate = "toy.layout",
      stage = "bundle",
      seed = 1L,
      status = "ok",
      coords_path = file.path("artifacts", "layouts", "toy_embedding.tsv"),
      stringsAsFactors = FALSE
    ),
    file.path(root, "catalog", "layout_catalog.tsv")
  )

  project <- gripui_project_from_dir(root)

  expect_s3_class(project, "gripui_project")
  expect_true(project$layouts$viewable[[1L]])
  expect_match(project$layouts$coords_path[[1L]], "toy_embedding.tsv$")

  coords <- gripui_load_layout_coords(project, 1L)
  expect_equal(coords[1L, 1L], 1)
  expect_equal(coords[2L, 1L], 2)
  expect_equal(coords[3L, 1L], 3)
})

test_that("gripui_project_from_dir loads HMP-style run tables and manifests", {
  parent <- local_test_dir("gripui-hmp-parent")
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)

  root <- file.path(parent, "hmp_u01_gc_coarsened_layout_selection_2026-03-23")
  repulsion.root <- file.path(parent, "hmp_u01_gc_repulsion_sweep_from_088_2026-03-23")
  dir.create(file.path(root, "coarse_stage"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "full_stage2"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "visuals", "coarse.search.088", "layouts"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "visuals", "coarse.search.088", "html"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(repulsion.root, "visuals_seed1_cst", "repulsion_1.25"), recursive = TRUE, showWarnings = FALSE)

  write_tsv(
    data.frame(
      candidate = "coarse.search.001",
      seed = 1L,
      status = "ok",
      `error.message` = "",
      stringsAsFactors = FALSE
    ),
    file.path(root, "coarse_stage", "coarse_stage_runs.tsv")
  )

  write_tsv(
    data.frame(
      candidate = "coarse.search.088",
      seed = 1L,
      status = "ok",
      `error.message` = "",
      `repulsion.factor` = 1.25,
      stringsAsFactors = FALSE
    ),
    file.path(root, "full_stage2", "full_stage2_runs.tsv")
  )

  full.embedding <- file.path(
    root, "visuals", "coarse.search.088", "layouts", "coarse.search.088_seed1_embedding.tsv"
  )
  write_tsv(
    data.frame(
      sample_id = c("s3", "s1", "s4", "s2"),
      x = c(30, 10, 40, 20),
      y = c(3, 1, 4, 2),
      z = c(300, 100, 400, 200),
      stringsAsFactors = FALSE
    ),
    full.embedding
  )
  full.html <- file.path(root, "visuals", "coarse.search.088", "html", "coarse.search.088_seed1_plain.html")
  writeLines("<html><body>full</body></html>", full.html)

  write_tsv(
    data.frame(
      candidate = "coarse.search.088",
      seed = 1L,
      plain_html = full.html,
      embedding_tsv = full.embedding,
      stringsAsFactors = FALSE
    ),
    file.path(root, "top_layout_visual_manifest.tsv")
  )

  write_tsv(
    data.frame(
      candidate = "repulsion.1.25",
      seed = 1L,
      status = "ok",
      `error.message` = "",
      `repulsion.factor` = 1.25,
      stringsAsFactors = FALSE
    ),
    file.path(repulsion.root, "repulsion_sweep_runs.tsv")
  )

  repulsion.embedding <- file.path(
    repulsion.root, "visuals_seed1_cst", "repulsion_1.25", "repulsion_1.25_seed1_embedding.tsv"
  )
  write_tsv(
    data.frame(
      sample_id = c("s1", "s2", "s3", "s4"),
      x = c(1, 2, 3, 4),
      y = c(11, 12, 13, 14),
      z = c(21, 22, 23, 24),
      stringsAsFactors = FALSE
    ),
    repulsion.embedding
  )
  repulsion.html <- file.path(
    repulsion.root, "visuals_seed1_cst", "repulsion_1.25", "repulsion_1.25_seed1_cst.html"
  )
  writeLines("<html><body>repulsion</body></html>", repulsion.html)

  write_tsv(
    data.frame(
      repulsion_factor = 1.25,
      seed = 1L,
      cst_html = repulsion.html,
      embedding_tsv = repulsion.embedding,
      stringsAsFactors = FALSE
    ),
    file.path(repulsion.root, "repulsion_sweep_visual_manifest.tsv")
  )

  graphs <- list(
    coarse = list(
      adj_list = list(2L, 1L),
      vertex_data = data.frame(vertex_id = c("c1", "c2"), stringsAsFactors = FALSE)
    ),
    full = list(
      adj_list = list(c(2L, 3L), c(1L, 4L), c(1L, 4L), c(2L, 3L)),
      vertex_data = data.frame(sample_id = c("s1", "s2", "s3", "s4"), stringsAsFactors = FALSE)
    )
  )

  project <- gripui_project_from_dir(root, graph = graphs, title = "HMP")

  expect_s3_class(project, "gripui_project")
  expect_equal(sort(unique(project$layouts$stage)), c("coarse_stage", "full_stage2", "repulsion_sweep"))
  expect_equal(sum(project$layouts$viewable), 2L)
  expect_equal(project$layouts$color_view_default[project$layouts$stage == "repulsion_sweep"], "cst")

  full.row <- project$layouts[project$layouts$stage == "full_stage2", , drop = FALSE]
  coords <- gripui_load_layout_coords(project, full.row$layout_id[[1L]])
  expect_equal(coords[1L, ], c(10, 1, 100))
  expect_equal(coords[3L, ], c(30, 3, 300))
})
