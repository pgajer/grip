test_that("basic graph helpers return two-column integer matrices", {
  full_mask <- matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE)
  occupied_keep <- matrix(c(1, 1, 0, 1, 1, 1), nrow = 2, byrow = TRUE)
  helpers <- list(
    edges.path(5),
    edges.cycle(6),
    edges.mesh(3, 4),
    edges.occupied.mesh(occupied_keep),
    edges.cylinder(3, 4),
    edges.torus(3, 4),
    edges.sphere(4, 5),
    edges.cube(3),
    edges.kary.tree(2, 3),
    edges.recursive.mask.grid(full_mask, 2),
    edges.recursive.triangle.mask(mask.triangle.classic(), 2),
    edges.vicsek(2),
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

test_that("recursive mask grid reproduces sierpinski carpet topology", {
  carpet_mask <- matrix(
    c(
      1, 1, 1,
      1, 0, 1,
      1, 1, 1
    ),
    nrow = 3,
    byrow = TRUE
  )

  expect_equal(edges.recursive.mask.grid(carpet_mask, 2), edges.sierpinski.carpet(2))
})

test_that("recursive triangle mask reproduces classic sierpinski triangle topology", {
  expect_equal(
    edges.recursive.triangle.mask(mask.triangle.classic(), 2),
    edges.sierpinski.triangle(2)
  )
})

test_that("vicsek graph labels occupied cells consecutively", {
  edges <- edges.vicsek(4)
  expect_equal(max(edges), 625L)
  expect_true(all(sort(unique(c(edges))) == seq_len(625L)))
})

test_that("mask helpers return expected connected motifs", {
  expect_equal(mask.cross(3, 1), matrix(c(
    FALSE, TRUE, FALSE,
    TRUE, TRUE, TRUE,
    FALSE, TRUE, FALSE
  ), nrow = 3, byrow = TRUE))
  expect_equal(sum(mask.border(5, 1)), 16L)
  expect_equal(sum(mask.corner(5, 2, corner = "top_left")), 4L)
  expect_equal(sum(mask.asymmetric.holes(5, 1)), 22L)
  expect_equal(mask.triangle.classic(), c(left = TRUE, right = TRUE, top = TRUE, center = FALSE))
  expect_equal(mask.triangle.bridge("top"), c(left = TRUE, right = TRUE, top = FALSE, center = TRUE))
  expect_equal(edges.recursive.mask.grid(mask.cross(3, 1), 2), edges.vicsek(2))
})

test_that("occupied mesh graph labels occupied cells consecutively", {
  keep <- keep.periodic.holes(6, 7, hole_period = 3, hole_height = 1, hole_width = 1)
  edges <- edges.occupied.mesh(keep)

  expect_equal(max(edges), sum(keep))
  expect_true(all(sort(unique(c(edges))) == seq_len(sum(keep))))
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

test_that("cylinder surface embedding returns finite 3D coordinates", {
  coords <- cylinder.surface.embedding(5, 8, surface = "barrel", amplitude = 0.25)

  expect_true(is.matrix(coords))
  expect_true(is.numeric(coords))
  expect_equal(dim(coords), c(40L, 3L))
  expect_true(all(is.finite(coords)))
  expect_equal(colnames(coords), c("x", "y", "z"))
})

test_that("cylinder surface graph returns normalized positive edge weights", {
  spec <- cylinder.surface.graph(5, 8, surface = "hourglass", amplitude = 0.35)

  expect_s3_class(spec, "grip_cylinder_surface_graph")
  expect_equal(spec$edges, edges.cylinder(5, 8))
  expect_equal(spec$n, 40L)
  expect_equal(length(spec$edge_weights), nrow(spec$edges))
  expect_true(all(is.finite(spec$edge_weights)))
  expect_true(all(spec$edge_weights > 0))
  expect_gt(stats::sd(spec$edge_weights), 0)
  expect_equal(stats::median(spec$edge_weights), 1, tolerance = 1e-10)
  expect_equal(dim(spec$coords_surface), c(40L, 3L))
  expect_equal(dim(spec$coords_param), c(40L, 2L))
  expect_equal(spec$family, "cylinder")
  expect_equal(spec$surface, "hourglass")
})

test_that("wavy cylinder graph supports alternate weight normalization", {
  spec <- cylinder.surface.graph(
    6, 9,
    surface = "wavy",
    amplitude = 0.2,
    freq_theta = 3,
    freq_z = 1.5,
    twist = 0.4,
    normalize = "mean"
  )

  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_gt(max(spec$edge_weights) - min(spec$edge_weights), 1e-6)
})

test_that("torus surface embedding returns finite 3D coordinates", {
  coords <- torus.surface.embedding(6, 9, surface = "standard")

  expect_true(is.matrix(coords))
  expect_true(is.numeric(coords))
  expect_equal(dim(coords), c(54L, 3L))
  expect_true(all(is.finite(coords)))
  expect_equal(colnames(coords), c("x", "y", "z"))
})

test_that("torus surface graph returns normalized positive edge weights", {
  spec <- torus.surface.graph(6, 9, surface = "pinched", amplitude = 0.2)

  expect_s3_class(spec, "grip_torus_surface_graph")
  expect_equal(spec$edges, edges.torus(6, 9))
  expect_equal(spec$n, 54L)
  expect_equal(length(spec$edge_weights), nrow(spec$edges))
  expect_true(all(is.finite(spec$edge_weights)))
  expect_true(all(spec$edge_weights > 0))
  expect_gt(stats::sd(spec$edge_weights), 0)
  expect_equal(stats::median(spec$edge_weights), 1, tolerance = 1e-10)
  expect_equal(dim(spec$coords_surface), c(54L, 3L))
  expect_equal(dim(spec$coords_param), c(54L, 2L))
  expect_equal(spec$family, "torus")
  expect_equal(spec$surface, "pinched")
})

test_that("wavy torus graph supports alternate weight normalization", {
  spec <- torus.surface.graph(
    7, 10,
    surface = "wavy",
    amplitude = 0.15,
    freq_major = 3,
    freq_minor = 2,
    twist = 0.35,
    normalize = "mean"
  )

  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_gt(max(spec$edge_weights) - min(spec$edge_weights), 1e-6)
})

test_that("sphere surface embedding returns finite 3D coordinates", {
  coords <- sphere.surface.embedding(6, 9, surface = "standard")

  expect_true(is.matrix(coords))
  expect_true(is.numeric(coords))
  expect_equal(dim(coords), c(38L, 3L))
  expect_true(all(is.finite(coords)))
  expect_equal(colnames(coords), c("x", "y", "z"))
})

test_that("sphere surface graph returns normalized positive edge weights", {
  spec <- sphere.surface.graph(6, 9, surface = "ellipsoid", amplitude = 0.2)

  expect_s3_class(spec, "grip_sphere_surface_graph")
  expect_equal(spec$edges, edges.sphere(6, 9))
  expect_equal(spec$n, 38L)
  expect_equal(length(spec$edge_weights), nrow(spec$edges))
  expect_true(all(is.finite(spec$edge_weights)))
  expect_true(all(spec$edge_weights > 0))
  expect_gt(stats::sd(spec$edge_weights), 0)
  expect_equal(stats::median(spec$edge_weights), 1, tolerance = 1e-10)
  expect_equal(dim(spec$coords_surface), c(38L, 3L))
  expect_equal(dim(spec$coords_param), c(38L, 2L))
  expect_equal(spec$family, "sphere")
  expect_equal(spec$surface, "ellipsoid")
})

test_that("wavy sphere graph supports alternate weight normalization", {
  spec <- sphere.surface.graph(
    7, 10,
    surface = "wavy",
    amplitude = 0.18,
    freq_theta = 3,
    freq_lat = 2,
    twist = 0.3,
    normalize = "mean"
  )

  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_gt(max(spec$edge_weights) - min(spec$edge_weights), 1e-6)
})

test_that("sierpinski triangle surface embedding returns finite 3D coordinates", {
  coords <- sierpinski.triangle.surface.embedding(2, surface = "folded", amplitude = 0.8)

  expect_true(is.matrix(coords))
  expect_true(is.numeric(coords))
  expect_equal(dim(coords), c(max(edges.sierpinski.triangle(2)), 3L))
  expect_true(all(is.finite(coords)))
  expect_equal(colnames(coords), c("x", "y", "z"))
})

test_that("sierpinski triangle surface graph returns normalized positive edge weights", {
  spec <- sierpinski.triangle.surface.graph(2, surface = "paraboloid", amplitude = 0.9)

  expect_s3_class(spec, "grip_sierpinski_triangle_surface_graph")
  expect_equal(spec$edges, edges.sierpinski.triangle(2))
  expect_equal(spec$n, max(edges.sierpinski.triangle(2)))
  expect_equal(length(spec$edge_weights), nrow(spec$edges))
  expect_true(all(is.finite(spec$edge_weights)))
  expect_true(all(spec$edge_weights > 0))
  expect_gt(stats::sd(spec$edge_weights), 0)
  expect_equal(stats::median(spec$edge_weights), 1, tolerance = 1e-10)
  expect_equal(dim(spec$coords_surface), c(spec$n, 3L))
  expect_equal(dim(spec$coords_param), c(spec$n, 2L))
  expect_equal(spec$family, "sierpinski.triangle")
  expect_equal(spec$surface, "paraboloid")
})

test_that("ripple sierpinski triangle graph supports alternate weight normalization", {
  spec <- sierpinski.triangle.surface.graph(
    2,
    surface = "ripple",
    amplitude = 0.6,
    freq_u = 1.5,
    freq_v = 0.5,
    normalize = "mean"
  )

  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_gt(max(spec$edge_weights) - min(spec$edge_weights), 1e-6)
})

test_that("sierpinski tetrahedron surface embedding returns finite 3D coordinates", {
  coords <- sierpinski.tetrahedron.surface.embedding(2, surface = "twisted", twist = 0.7)

  expect_true(is.matrix(coords))
  expect_true(is.numeric(coords))
  expect_equal(dim(coords), c(max(edges.sierpinski.tetrahedron(2)), 3L))
  expect_true(all(is.finite(coords)))
  expect_equal(colnames(coords), c("x", "y", "z"))
})

test_that("sierpinski tetrahedron surface graph returns normalized positive edge weights", {
  spec <- sierpinski.tetrahedron.surface.graph(2, surface = "wavy", amplitude = 0.25)

  expect_s3_class(spec, "grip_sierpinski_tetrahedron_surface_graph")
  expect_equal(spec$edges, edges.sierpinski.tetrahedron(2))
  expect_equal(spec$n, max(edges.sierpinski.tetrahedron(2)))
  expect_equal(length(spec$edge_weights), nrow(spec$edges))
  expect_true(all(is.finite(spec$edge_weights)))
  expect_true(all(spec$edge_weights > 0))
  expect_gt(stats::sd(spec$edge_weights), 0)
  expect_equal(stats::median(spec$edge_weights), 1, tolerance = 1e-10)
  expect_equal(dim(spec$coords_surface), c(spec$n, 3L))
  expect_equal(dim(spec$coords_param), c(spec$n, 3L))
  expect_equal(spec$family, "sierpinski.tetrahedron")
  expect_equal(spec$surface, "wavy")
})

test_that("recursive triangle mask surface embedding returns finite 3D coordinates", {
  coords <- recursive.triangle.mask.surface.embedding(
    mask = mask.triangle.bridge("right"),
    level = 2,
    surface = "folded",
    amplitude = 0.7
  )

  expect_true(is.matrix(coords))
  expect_true(is.numeric(coords))
  expect_true(all(is.finite(coords)))
  expect_equal(ncol(coords), 3L)
  expect_equal(colnames(coords), c("x", "y", "z"))
})

test_that("recursive triangle mask surface graph returns normalized positive edge weights", {
  bridge_mask <- mask.triangle.bridge("top")
  spec <- recursive.triangle.mask.surface.graph(
    mask = bridge_mask,
    level = 2,
    surface = "ripple",
    amplitude = 0.55,
    freq_u = 1.25,
    freq_v = 0.75
  )

  expect_s3_class(spec, "grip_recursive_triangle_mask_surface_graph")
  expect_equal(spec$edges, edges.recursive.triangle.mask(bridge_mask, 2))
  expect_equal(spec$n, max(spec$edges))
  expect_equal(spec$mask, bridge_mask)
  expect_equal(length(spec$edge_weights), nrow(spec$edges))
  expect_true(all(is.finite(spec$edge_weights)))
  expect_true(all(spec$edge_weights > 0))
  expect_gt(stats::sd(spec$edge_weights), 0)
  expect_equal(stats::median(spec$edge_weights), 1, tolerance = 1e-10)
  expect_equal(dim(spec$coords_surface), c(spec$n, 3L))
  expect_equal(dim(spec$coords_param), c(spec$n, 2L))
  expect_equal(spec$family, "recursive.triangle.mask")
  expect_equal(spec$surface, "ripple")
})

test_that("recursive mask grid surface graph returns normalized positive edge weights", {
  vicsek_mask <- matrix(
    c(
      0, 1, 0,
      1, 1, 1,
      0, 1, 0
    ),
    nrow = 3,
    byrow = TRUE
  )
  spec <- recursive.mask.grid.surface.graph(
    vicsek_mask,
    2,
    surface = "saddle",
    amplitude = 0.8
  )

  expect_s3_class(spec, "grip_recursive_mask_grid_surface_graph")
  expect_equal(spec$edges, edges.vicsek(2))
  expect_equal(spec$n, 25L)
  expect_equal(spec$mask, vicsek_mask != 0)
  expect_equal(spec$mask_size, 3L)
  expect_equal(spec$side, 9L)
  expect_true(all(spec$edge_weights > 0))
  expect_equal(stats::median(spec$edge_weights), 1, tolerance = 1e-10)
})

test_that("sierpinski carpet surface embedding returns finite 3D coordinates", {
  coords <- sierpinski.carpet.surface.embedding(2, surface = "saddle", amplitude = 0.8)

  expect_true(is.matrix(coords))
  expect_true(is.numeric(coords))
  expect_equal(dim(coords), c(64L, 3L))
  expect_true(all(is.finite(coords)))
  expect_equal(colnames(coords), c("x", "y", "z"))
})

test_that("sierpinski carpet surface graph returns normalized positive edge weights", {
  spec <- sierpinski.carpet.surface.graph(2, surface = "paraboloid", amplitude = 0.9)

  expect_s3_class(spec, "grip_sierpinski_carpet_surface_graph")
  expect_equal(spec$edges, edges.sierpinski.carpet(2))
  expect_equal(spec$n, 64L)
  expect_equal(length(spec$edge_weights), nrow(spec$edges))
  expect_true(all(is.finite(spec$edge_weights)))
  expect_true(all(spec$edge_weights > 0))
  expect_gt(stats::sd(spec$edge_weights), 0)
  expect_equal(stats::median(spec$edge_weights), 1, tolerance = 1e-10)
  expect_equal(dim(spec$coords_surface), c(64L, 3L))
  expect_equal(dim(spec$coords_param), c(64L, 2L))
  expect_equal(spec$family, "sierpinski.carpet")
  expect_equal(spec$surface, "paraboloid")
  expect_equal(spec$level, 2L)
})

test_that("ripple sierpinski carpet graph supports alternate weight normalization", {
  spec <- sierpinski.carpet.surface.graph(
    2,
    surface = "ripple",
    amplitude = 0.6,
    freq_u = 1.5,
    freq_v = 0.5,
    normalize = "mean"
  )

  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_gt(max(spec$edge_weights) - min(spec$edge_weights), 1e-6)
})

test_that("vicsek surface graph returns normalized positive edge weights", {
  spec <- vicsek.surface.graph(
    2,
    surface = "ripple",
    amplitude = 0.6,
    freq_u = 1.5,
    freq_v = 0.5,
    normalize = "mean"
  )

  expect_s3_class(spec, "grip_vicsek_surface_graph")
  expect_equal(spec$edges, edges.vicsek(2))
  expect_equal(spec$n, 25L)
  expect_equal(spec$family, "vicsek")
  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_gt(max(spec$edge_weights) - min(spec$edge_weights), 1e-6)
})

test_that("occupied mesh surface graph returns normalized positive edge weights", {
  keep <- keep.staggered.windows(
    8, 9,
    window_height = 1,
    window_width = 2,
    row_period = 3,
    col_period = 4
  )
  spec <- occupied.mesh.surface.graph(
    keep,
    surface = "paraboloid",
    amplitude = 0.9
  )

  expect_s3_class(spec, "grip_occupied_mesh_surface_graph")
  expect_equal(spec$edges, edges.occupied.mesh(keep))
  expect_equal(spec$n, sum(keep))
  expect_equal(length(spec$edge_weights), nrow(spec$edges))
  expect_true(all(spec$edge_weights > 0))
  expect_equal(stats::median(spec$edge_weights), 1, tolerance = 1e-10)
  expect_equal(dim(spec$coords_param), c(sum(keep), 2L))
  expect_equal(spec$family, "occupied.mesh")
})

test_that("deterministic perforated grids produce distinct occupancy patterns", {
  periodic <- keep.periodic.holes(9, 9, hole_period = 4, hole_height = 1, hole_width = 1)
  staggered <- keep.staggered.windows(9, 9, window_height = 1, window_width = 2)
  slits <- keep.slit.channels(9, 9, orientation = "vertical", slit_period = 4)
  notches <- keep.asymmetric.notches(9, 10, notch_depth = 3, notch_width = 2)

  expect_true(is.matrix(periodic) && is.logical(periodic))
  expect_true(is.matrix(staggered) && is.logical(staggered))
  expect_true(is.matrix(slits) && is.logical(slits))
  expect_true(is.matrix(notches) && is.logical(notches))
  expect_false(identical(periodic, staggered))
  expect_false(identical(periodic, slits))
  expect_true(sum(periodic) < length(periodic))
  expect_true(sum(staggered) < length(staggered))
  expect_true(sum(slits) < length(slits))
  expect_true(sum(notches) < length(notches))
})
