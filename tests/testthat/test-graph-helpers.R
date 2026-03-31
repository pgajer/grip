test_that("basic graph helpers return two-column integer matrices", {
  helpers <- list(
    edges.path(5),
    edges.cycle(6),
    edges.mesh(3, 4),
    edges.cylinder(3, 4),
    edges.torus(3, 4),
    edges.sphere(4, 5),
    edges.cube(3),
    edges.kary.tree(2, 3),
    edges.sierpinski.triangle(3),
    edges.sierpinski.tetrahedron(2),
    edges.sierpinski.carpet(2)
  )

  for (edges in helpers) {
    expect_type(edges, "integer")
    expect_equal(ncol(edges), 2)
    expect_true(all(edges[, 1] != edges[, 2]))
  }
})

test_that("sierpinski triangle counts match the paper depth convention", {
  expect_equal(max(edges.sierpinski.triangle(5)), 366L)
  expect_equal(max(edges.sierpinski.triangle(6)), 1095L)
})

test_that("sierpinski tetrahedron counts match the legacy generator", {
  expect_equal(max(edges.sierpinski.tetrahedron(5)), 2050L)
  expect_equal(max(edges.sierpinski.tetrahedron(6)), 8194L)
})

test_that("sierpinski carpet labels all occupied cells consecutively", {
  edges <- edges.sierpinski.carpet(4)
  expect_equal(max(edges), 4096L)
  expect_true(all(sort(unique(c(edges))) == seq_len(4096L)))
})

test_that("sphere graph labels intermediate latitude rings consecutively", {
  edges <- edges.sphere(5, 8)
  expect_equal(max(edges), 26L)
  expect_true(all(sort(unique(c(edges))) == seq_len(26L)))
})

test_that("cube graph matches the cube-surface vertex count", {
  edges <- edges.cube(4)
  expect_equal(max(edges), 56L)
  expect_true(all(sort(unique(c(edges))) == seq_len(56L)))
})

test_that("mesh surface embedding returns finite 3D coordinates", {
  coords <- mesh.surface.embedding(4, 5, surface = "saddle", amplitude = 0.8)

  expect_true(is.matrix(coords))
  expect_true(is.numeric(coords))
  expect_equal(dim(coords), c(20L, 3L))
  expect_true(all(is.finite(coords)))
  expect_equal(colnames(coords), c("x", "y", "z"))
})

test_that("mesh surface graph returns normalized positive edge weights", {
  spec <- mesh.surface.graph(5, 6, surface = "paraboloid", amplitude = 0.9)

  expect_s3_class(spec, "grip_mesh_surface_graph")
  expect_equal(spec$edges, edges.mesh(5, 6))
  expect_equal(spec$n, 30L)
  expect_equal(length(spec$edge_weights), nrow(spec$edges))
  expect_true(all(is.finite(spec$edge_weights)))
  expect_true(all(spec$edge_weights > 0))
  expect_gt(stats::sd(spec$edge_weights), 0)
  expect_equal(stats::median(spec$edge_weights), 1, tolerance = 1e-10)
  expect_equal(dim(spec$coords_surface), c(30L, 3L))
  expect_equal(dim(spec$coords_param), c(30L, 2L))
  expect_equal(spec$family, "mesh")
  expect_equal(spec$surface, "paraboloid")
})

test_that("ripple mesh graph supports alternate weight normalization", {
  spec <- mesh.surface.graph(
    4, 4,
    surface = "ripple",
    amplitude = 0.6,
    freq_u = 1.5,
    freq_v = 0.5,
    normalize = "mean"
  )

  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_gt(max(spec$edge_weights) - min(spec$edge_weights), 1e-6)
})
