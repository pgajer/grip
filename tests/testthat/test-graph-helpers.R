test_that("basic graph helpers return two-column integer matrices", {
  helpers <- list(
    edges.path(5),
    edges.cycle(6),
    edges.mesh(3, 4),
    edges.cylinder(3, 4),
    edges.torus(3, 4),
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
