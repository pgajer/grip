test_that("basic graph helpers return two-column integer matrices", {
  full_mask <- matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE)
  occupied_keep <- matrix(c(1, 1, 0, 1, 1, 1), nrow = 2, byrow = TRUE)
  helpers <- list(
    edges.path(5),
    edges.cycle(6),
    edges.mesh(3, 4),
    grip:::edges.occupied.mesh(occupied_keep),
    edges.cylinder(3, 4),
    grip:::edges.irregular.ball(),
    grip:::edges.irregular.shell(),
    edges.torus(3, 4),
    grip:::edges.irregular.torus(),
    grip:::edges.sphere(4, 5),
    edges.cube(3),
    grip:::edges.recursive.cube.mask(array(TRUE, dim = c(2, 2, 2)), 2),
    grip:::edges.cube.periodic.tunnels(1),
    grip:::edges.cube.asymmetric.cavities(1),
    grip:::edges.cube.channel.network(1),
    grip:::edges.triangulated.polyhedron("octahedron", 1),
    grip:::edges.triangulated.annulus(7),
    grip:::edges.triangulated.pair.of.pants(7),
    grip:::edges.irregular.annulus(),
    grip:::edges.irregular.pair.of.pants(),
    grip:::edges.irregular.double.torus(),
    grip:::edges.irregular.sphere(),
    edges.kary.tree(2, 3),
    grip:::edges.recursive.mask.grid(full_mask, 2),
    grip:::edges.recursive.triangle.mask(mask.triangle.classic(), 2),
    grip:::edges.recursive.tetrahedron.mask(mask.tetrahedron.classic(), 2),
    grip:::edges.menger.sponge(1),
    grip:::edges.vicsek(2),
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

test_that("triangulated polyhedron counts match closed-surface formulas", {
  expect_equal(max(grip:::edges.triangulated.polyhedron("tetrahedron", 1)), 10L)
  expect_equal(max(grip:::edges.triangulated.polyhedron("octahedron", 1)), 18L)
  expect_equal(max(grip:::edges.triangulated.polyhedron("icosahedron", 1)), 42L)
})

test_that("weighted kary tree graph follows custom depth and branch rules", {
  spec <- kary.tree.weighted.graph(
    k = 2,
    depth = 2,
    base_length = 2,
    depth_rule = "custom",
    depth_factors = c(1, 3),
    branch_rule = "custom",
    branch_factors = c(1, 2),
    normalize = "none"
  )

  expect_s3_class(spec, "grip_kary_tree_weighted_graph")
  expect_equal(spec$edges, edges.kary.tree(2, 2))
  expect_equal(spec$n, 7L)
  expect_equal(spec$edge_weights, c(2, 4, 6, 12, 6, 12))
  expect_equal(spec$weight_scale, 1)
  expect_equal(spec$vertex_depth, c(0L, 1L, 1L, 2L, 2L, 2L, 2L))
  expect_equal(spec$parent, c(0L, 1L, 1L, 2L, 2L, 3L, 3L))
  expect_equal(spec$depth_factors, c(1, 3))
  expect_equal(spec$branch_factors, c(1, 2))
  expect_equal(spec$edge_table$branch_index, c(1L, 2L, 1L, 2L, 1L, 2L))
  expect_equal(spec$edge_table$raw_weight, c(2, 4, 6, 12, 6, 12))
  expect_equal(spec$family, "kary.tree.weighted")
})

test_that("weighted kary tree graph supports built-in rule normalization", {
  spec <- kary.tree.weighted.graph(
    k = 3,
    depth = 2,
    depth_rule = "geometric",
    depth_decay = 0.7,
    branch_rule = "linear",
    branch_spread = 0.4,
    normalize = "mean"
  )

  expect_equal(spec$edges, edges.kary.tree(3, 2))
  expect_equal(spec$n, 13L)
  expect_equal(length(spec$edge_weights), nrow(spec$edges))
  expect_true(all(spec$edge_weights > 0))
  expect_gt(stats::sd(spec$edge_weights), 0)
  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_equal(spec$depth_factors, c(1, 0.7))
  expect_equal(spec$branch_factors, c(0.8, 1.0, 1.2))
  expect_equal(spec$family, "kary.tree.weighted")
})

test_that("weighted kary tree graph handles depth zero", {
  spec <- kary.tree.weighted.graph(k = 4, depth = 0)

  expect_equal(spec$edges, matrix(integer(), ncol = 2L))
  expect_equal(spec$n, 1L)
  expect_equal(spec$edge_weights, numeric())
  expect_equal(spec$vertex_depth, 0L)
  expect_equal(spec$parent, 0L)
  expect_equal(nrow(spec$edge_table), 0L)
  expect_equal(spec$depth_factors, numeric())
  expect_equal(length(spec$branch_factors), 4L)
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

  expect_equal(grip:::edges.recursive.mask.grid(carpet_mask, 2), edges.sierpinski.carpet(2))
})

test_that("recursive cube mask reproduces Menger sponge topology", {
  menger_mask <- array(TRUE, dim = c(3, 3, 3))
  for (row in seq_len(3)) {
    for (col in seq_len(3)) {
      for (layer in seq_len(3)) {
        if (sum(c(row, col, layer) == 2L) >= 2L) {
          menger_mask[row, col, layer] <- FALSE
        }
      }
    }
  }

  expect_equal(grip:::edges.recursive.cube.mask(menger_mask, 2), grip:::edges.menger.sponge(2))
})

test_that("recursive triangle mask reproduces classic sierpinski triangle topology", {
  expect_equal(
    grip:::edges.recursive.triangle.mask(mask.triangle.classic(), 2),
    edges.sierpinski.triangle(2)
  )
})

test_that("recursive tetrahedron mask reproduces classic sierpinski tetrahedron topology", {
  expect_equal(
    grip:::edges.recursive.tetrahedron.mask(mask.tetrahedron.classic(), 2),
    edges.sierpinski.tetrahedron(2)
  )
})

test_that("vicsek graph labels occupied cells consecutively", {
  edges <- grip:::edges.vicsek(4)
  expect_equal(max(edges), 625L)
  expect_true(all(sort(unique(c(edges))) == seq_len(625L)))
})

test_that("menger sponge graph labels occupied cells consecutively", {
  edges <- grip:::edges.menger.sponge(2)
  expect_equal(max(edges), 400L)
  expect_true(all(sort(unique(c(edges))) == seq_len(400L)))
})

test_that("named porous cube families label occupied cells consecutively", {
  periodic <- grip:::edges.cube.periodic.tunnels(level = 1)
  cavities <- grip:::edges.cube.asymmetric.cavities(level = 1)
  channels <- grip:::edges.cube.channel.network(level = 1)

  expect_equal(max(periodic), sum(mask.cube.periodic.tunnels()))
  expect_equal(max(cavities), sum(mask.cube.asymmetric.cavities()))
  expect_equal(max(channels), sum(mask.cube.channel.network()))
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
  expect_equal(mask.tetrahedron.classic(),
               c(base_left = TRUE, base_right = TRUE, base_back = TRUE, apex = TRUE))
  expect_equal(mask.tetrahedron.corner.missing("apex"),
               c(base_left = TRUE, base_right = TRUE, base_back = TRUE, apex = FALSE))
  expect_equal(grip:::edges.recursive.mask.grid(mask.cross(3, 1), 2), grip:::edges.vicsek(2))
})

test_that("cube mask helpers return distinct porous keep-arrays", {
  periodic <- mask.cube.periodic.tunnels()
  cavities <- mask.cube.asymmetric.cavities()
  channels <- mask.cube.channel.network()

  expect_true(is.array(periodic) && is.logical(periodic))
  expect_true(is.array(cavities) && is.logical(cavities))
  expect_true(is.array(channels) && is.logical(channels))
  expect_equal(dim(periodic), c(5L, 5L, 5L))
  expect_equal(dim(cavities), c(5L, 5L, 5L))
  expect_equal(dim(channels), c(5L, 5L, 5L))
  expect_false(identical(periodic, cavities))
  expect_false(identical(periodic, channels))
  expect_false(identical(cavities, channels))
})

test_that("periodic tunnel cube mask contains Menger sponge as the 3x3 special case", {
  expect_equal(
    grip:::edges.cube.periodic.tunnels(
      level = 2,
      side = 3,
      tunnel_width = 1,
      tunnel_period = 2,
      tunnel_offset = 2
    ),
    grip:::edges.menger.sponge(2)
  )
})

test_that("occupied mesh graph labels occupied cells consecutively", {
  keep <- keep.periodic.holes(6, 7, hole_period = 3, hole_height = 1, hole_width = 1)
  edges <- grip:::edges.occupied.mesh(keep)

  expect_equal(max(edges), sum(keep))
  expect_true(all(sort(unique(c(edges))) == seq_len(sum(keep))))
})

test_that("mesh connectivity option adds diagonal edges deterministically", {
  orth <- edges.mesh(3, 4, connectivity = "orthogonal")
  diag <- edges.mesh(3, 4, connectivity = "diagonal")

  expect_equal(nrow(orth), 17L)
  expect_equal(nrow(diag), 29L)
  expect_true(all(apply(orth, 1L, function(e) any(diag[, 1L] == e[[1L]] & diag[, 2L] == e[[2L]]))))

  keep_full <- matrix(TRUE, nrow = 2, ncol = 2)
  keep_missing <- matrix(c(TRUE, TRUE, TRUE, FALSE), nrow = 2, byrow = TRUE)
  expect_equal(nrow(grip:::edges.occupied.mesh(keep_full, connectivity = "diagonal")), 6L)
  expect_equal(
    grip:::edges.occupied.mesh(keep_missing, connectivity = "diagonal"),
    grip:::edges.occupied.mesh(keep_missing, connectivity = "orthogonal")
  )
})

test_that("sphere graph labels intermediate latitude rings consecutively", {
  edges <- grip:::edges.sphere(5, 8)
  expect_equal(max(edges), 26L)
  expect_true(all(sort(unique(c(edges))) == seq_len(26L)))
})

test_that("cube graph matches the cube-surface vertex count", {
  edges <- edges.cube(4)
  expect_equal(max(edges), 56L)
  expect_true(all(sort(unique(c(edges))) == seq_len(56L)))
})

test_that("triangulated polyhedron graph labels vertices consecutively", {
  edges <- grip:::edges.triangulated.polyhedron("icosahedron", 1)
  expect_equal(max(edges), 42L)
  expect_true(all(sort(unique(c(edges))) == seq_len(42L)))
})

test_that("triangulated annulus graph labels vertices consecutively", {
  edges <- grip:::edges.triangulated.annulus(resolution = 7)
  expect_true(all(sort(unique(c(edges))) == seq_len(max(edges))))
  expect_gt(max(edges), 0L)
})

test_that("triangulated pair-of-pants graph labels vertices consecutively", {
  edges <- grip:::edges.triangulated.pair.of.pants(resolution = 7)
  expect_true(all(sort(unique(c(edges))) == seq_len(max(edges))))
  expect_gt(max(edges), 0L)
})

test_that("irregular annulus graph labels vertices consecutively", {
  edges <- grip:::edges.irregular.annulus(rings = 6, outer_count = 24)
  expect_true(all(sort(unique(c(edges))) == seq_len(max(edges))))
  expect_gt(max(edges), 0L)
})

test_that("irregular ball graph labels vertices consecutively", {
  edges <- grip:::edges.irregular.ball(base = "icosahedron", level = 1, layers = 3)
  expect_true(all(sort(unique(c(edges))) == seq_len(max(edges))))
  expect_gt(max(edges), 0L)
})

test_that("irregular shell graph labels vertices consecutively", {
  edges <- grip:::edges.irregular.shell(base = "octahedron", level = 1, layers = 3)
  expect_true(all(sort(unique(c(edges))) == seq_len(max(edges))))
  expect_gt(max(edges), 0L)
})

test_that("irregular pair-of-pants graph labels vertices consecutively", {
  edges <- grip:::edges.irregular.pair.of.pants(slices = 11, outer_count = 24)
  expect_true(all(sort(unique(c(edges))) == seq_len(max(edges))))
  expect_gt(max(edges), 0L)
})

test_that("irregular torus graph labels vertices consecutively", {
  edges <- grip:::edges.irregular.torus(major_rings = 8, tube_count = 16)
  expect_true(all(sort(unique(c(edges))) == seq_len(max(edges))))
  expect_gt(max(edges), 0L)
})

test_that("irregular double torus graph labels vertices consecutively", {
  edges <- grip:::edges.irregular.double.torus(slices = 11, tube_count = 14)
  expect_true(all(sort(unique(c(edges))) == seq_len(max(edges))))
  expect_gt(max(edges), 0L)
})

test_that("irregular sphere graph labels vertices consecutively", {
  edges <- grip:::edges.irregular.sphere(bands = 6, equator_count = 24)
  expect_true(all(sort(unique(c(edges))) == seq_len(max(edges))))
  expect_gt(max(edges), 0L)
})

test_that("mesh surface embedding returns finite 3D coordinates", {
  coords <- grip:::mesh.surface.embedding(4, 5, surface = "saddle", amplitude = 0.8)

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
  expect_equal(spec$connectivity, "orthogonal")
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

test_that("mesh surface graph supports diagonal connectivity", {
  spec <- mesh.surface.graph(
    4, 4,
    surface = "ripple",
    amplitude = 0.6,
    freq_u = 1.5,
    freq_v = 0.5,
    connectivity = "diagonal",
    normalize = "mean"
  )

  expect_equal(spec$edges, edges.mesh(4, 4, connectivity = "diagonal"))
  expect_equal(spec$connectivity, "diagonal")
  expect_equal(length(spec$edge_weights), nrow(spec$edges))
  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
})

test_that("irregular rectangle parameter coordinates are deterministic and ordered", {
  coords1 <- grip:::irregular.rectangle.param.coords(
    5, 6,
    row_irregularity = 0.22,
    col_irregularity = 0.18,
    interior_warp = 0.06
  )
  coords2 <- grip:::irregular.rectangle.param.coords(
    5, 6,
    row_irregularity = 0.22,
    col_irregularity = 0.18,
    interior_warp = 0.06
  )

  expect_equal(coords1, coords2)
  expect_true(is.matrix(coords1))
  expect_true(all(is.finite(coords1)))
  expect_equal(dim(coords1), c(30L, 2L))
  expect_equal(colnames(coords1), c("u", "v"))

  u_mat <- matrix(coords1[, 1L], nrow = 5L, ncol = 6L, byrow = TRUE)
  v_mat <- matrix(coords1[, 2L], nrow = 5L, ncol = 6L, byrow = TRUE)
  expect_true(all(apply(u_mat, 1L, function(x) all(diff(x) > 0))))
  expect_true(all(apply(v_mat, 2L, function(x) all(diff(x) < 0))))
  expect_equal(u_mat[, 1L], rep(-1, 5L), tolerance = 1e-10)
  expect_equal(u_mat[, 6L], rep(1, 5L), tolerance = 1e-10)
  expect_equal(v_mat[1L, ], rep(1, 6L), tolerance = 1e-10)
  expect_equal(v_mat[5L, ], rep(-1, 6L), tolerance = 1e-10)
})

test_that("irregular rectangle surface embedding supports flat and curved lifts", {
  flat <- grip:::irregular.rectangle.surface.embedding(
    4, 5,
    surface = "flat",
    row_irregularity = 0.2,
    col_irregularity = 0.15
  )
  curved <- grip:::irregular.rectangle.surface.embedding(
    4, 5,
    surface = "paraboloid",
    amplitude = 0.7,
    row_irregularity = 0.2,
    col_irregularity = 0.15
  )

  expect_equal(dim(flat), c(20L, 3L))
  expect_equal(colnames(flat), c("x", "y", "z"))
  expect_true(all(flat[, 3L] == 0))
  expect_true(all(is.finite(curved)))
  expect_gt(max(abs(curved[, 3L])), 0)
})

test_that("irregular rectangle surface graph returns normalized positive edge weights", {
  spec <- irregular.rectangle.surface.graph(
    5, 6,
    surface = "paraboloid",
    amplitude = 0.8,
    row_irregularity = 0.2,
    col_irregularity = 0.15,
    interior_warp = 0.07,
    normalize = "mean"
  )

  expect_s3_class(spec, "grip_irregular_rectangle_surface_graph")
  expect_equal(spec$edges, edges.mesh(5, 6))
  expect_equal(spec$n, 30L)
  expect_equal(length(spec$edge_weights), nrow(spec$edges))
  expect_true(all(is.finite(spec$edge_weights)))
  expect_true(all(spec$edge_weights > 0))
  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_gt(stats::sd(spec$edge_weights), 0)
  expect_equal(dim(spec$coords_surface), c(30L, 3L))
  expect_equal(dim(spec$coords_param), c(30L, 2L))
  expect_equal(dim(spec$coords_regular_param), c(30L, 2L))
  expect_equal(length(spec$row_breaks), 5L)
  expect_equal(length(spec$col_breaks), 6L)
  expect_equal(spec$family, "irregular.rectangle")
  expect_equal(spec$surface, "paraboloid")
  expect_equal(spec$connectivity, "orthogonal")
  expect_gt(max(abs(spec$coords_param - spec$coords_regular_param)), 1e-6)
})

test_that("irregular rectangle surface graph supports diagonal connectivity and square cases", {
  spec <- irregular.rectangle.surface.graph(
    4, 4,
    surface = "ripple",
    amplitude = 0.5,
    freq_u = 1.5,
    freq_v = 0.75,
    connectivity = "diagonal",
    row_irregularity = 0.18,
    col_irregularity = 0.18
  )

  expect_equal(spec$edges, edges.mesh(4, 4, connectivity = "diagonal"))
  expect_equal(spec$connectivity, "diagonal")
  expect_equal(length(spec$edge_weights), nrow(spec$edges))
  expect_true(all(spec$edge_weights > 0))
  expect_equal(stats::median(spec$edge_weights), 1, tolerance = 1e-10)
  expect_equal(spec$n, 16L)
})

test_that("sampled rectangle helpers are reproducible and respect bounds", {
  coords1 <- grip:::sampled.rectangle.param.coords(
    n = 12,
    xmin = -2,
    xmax = 3,
    ymin = -1,
    ymax = 4,
    seed = 42
  )
  coords2 <- grip:::sampled.rectangle.param.coords(
    n = 12,
    xmin = -2,
    xmax = 3,
    ymin = -1,
    ymax = 4,
    seed = 42
  )
  embed <- grip:::sampled.rectangle.surface.embedding(
    n = 12,
    xmin = -2,
    xmax = 3,
    ymin = -1,
    ymax = 4,
    seed = 42,
    surface = "folded",
    amplitude = 0.6
  )

  expect_equal(coords1, coords2)
  expect_equal(dim(coords1), c(12L, 2L))
  expect_true(all(coords1[, 1L] >= -2 & coords1[, 1L] <= 3))
  expect_true(all(coords1[, 2L] >= -1 & coords1[, 2L] <= 4))
  expect_equal(dim(embed), c(12L, 3L))
  expect_equal(colnames(embed), c("x", "y", "z"))
  expect_true(all(is.finite(embed)))
  expect_gt(max(embed[, 3L]) - min(embed[, 3L]), 0)
})

test_that("sampled rectangle surface graph returns normalized positive iKNN weights", {
  spec <- sampled.rectangle.surface.graph(
    n = 40,
    k = 5,
    xmin = -1.5,
    xmax = 1.5,
    ymin = -0.5,
    ymax = 2,
    seed = 7,
    surface = "paraboloid",
    amplitude = 0.8,
    graph_space = "surface",
    normalize = "mean"
  )

  expect_s3_class(spec, "grip_sampled_rectangle_surface_graph")
  expect_equal(spec$n, 40L)
  expect_true(is.matrix(spec$edges))
  expect_equal(ncol(spec$edges), 2L)
  expect_equal(length(spec$edge_weights), nrow(spec$edges))
  expect_true(all(is.finite(spec$edge_weights)))
  expect_true(all(spec$edge_weights > 0))
  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_equal(length(spec$raw_edge_weights), nrow(spec$edges))
  expect_equal(length(spec$iknn_witness_edge_weights), nrow(spec$edges))
  expect_equal(dim(spec$coords_surface), c(40L, 3L))
  expect_equal(dim(spec$coords_param), c(40L, 2L))
  expect_equal(dim(spec$coords_param_unit), c(40L, 2L))
  expect_equal(spec$family, "sampled.rectangle")
  expect_equal(spec$surface, "paraboloid")
  expect_equal(spec$graph_space, "surface")
  expect_equal(spec$k, 5L)
  expect_false(isTRUE(all.equal(spec$raw_edge_weights, spec$iknn_witness_edge_weights)))
})

test_that("sampled rectangle surface graphs reuse one sample across k values", {
  seq_spec <- sampled.rectangle.surface.graphs(
    n = 36,
    k = c(3, 6),
    xmin = -1,
    xmax = 2,
    ymin = -2,
    ymax = 1,
    seed = 11,
    surface = "ripple",
    amplitude = 0.5,
    freq_u = 2,
    freq_v = 3,
    graph_space = "param"
  )

  expect_s3_class(seq_spec, "grip_sampled_rectangle_surface_graphs")
  expect_equal(seq_spec$k, c(3L, 6L))
  expect_equal(names(seq_spec$graphs), c("k3", "k6"))
  expect_equal(dim(seq_spec$coords_surface), c(36L, 3L))
  expect_equal(dim(seq_spec$coords_param), c(36L, 2L))
  expect_equal(dim(seq_spec$k_statistics), c(2L, 5L))
  expect_equal(seq_spec$family, "sampled.rectangle")
  expect_equal(seq_spec$graph_space, "param")
  expect_equal(seq_spec$graphs$k3$coords_param, seq_spec$coords_param)
  expect_equal(seq_spec$graphs$k6$coords_surface, seq_spec$coords_surface)
  expect_equal(length(seq_spec$graphs$k3$iknn_witness_edge_weights), nrow(seq_spec$graphs$k3$edges))
  expect_true(all(vapply(seq_spec$graphs, function(g) all(g$edge_weights > 0), logical(1L))))
  expect_true(nrow(seq_spec$graphs$k6$edges) >= nrow(seq_spec$graphs$k3$edges))
})

test_that("sampled rectangle graphs can be rebuilt from saved parameter coordinates", {
  spec <- sampled.rectangle.surface.graph(
    n = 32,
    k = 4,
    xmin = -1,
    xmax = 2,
    ymin = -2,
    ymax = 1,
    seed = 7,
    surface = "paraboloid",
    amplitude = 0.6,
    graph_space = "surface",
    normalize = "median"
  )

  rebuilt <- grip:::.sampled.rectangle.surface.graph.from.coords(
    coords_param = spec$coords_param,
    k = spec$k,
    xmin = spec$xmin,
    xmax = spec$xmax,
    ymin = spec$ymin,
    ymax = spec$ymax,
    seed = spec$seed,
    surface = spec$surface,
    amplitude = 0.6,
    graph_space = spec$graph_space,
    normalize = spec$normalize
  )

  expect_equal(rebuilt$edges, spec$edges)
  expect_equal(rebuilt$raw_edge_weights, spec$raw_edge_weights)
  expect_equal(rebuilt$coords_surface, spec$coords_surface)
  expect_equal(rebuilt$coords_param_unit, spec$coords_param_unit)
})

test_that("sampled rectangle saved topology can be reweighted on a new surface", {
  spec <- sampled.rectangle.surface.graph(
    n = 28,
    k = 5,
    xmin = -1,
    xmax = 1,
    ymin = -1,
    ymax = 1,
    seed = 13,
    surface = "paraboloid",
    amplitude = 0.7,
    graph_space = "surface",
    normalize = "mean"
  )

  reweighted <- grip:::.sampled.rectangle.surface.graph.reweight.saved.topology(
    graph = spec,
    surface = "saddle",
    amplitude = 0.7,
    normalize = "mean"
  )

  expect_equal(reweighted$edges, spec$edges)
  expect_equal(reweighted$coords_param, spec$coords_param)
  expect_false(isTRUE(all.equal(reweighted$coords_surface, spec$coords_surface)))
  expect_false(isTRUE(all.equal(reweighted$raw_edge_weights, spec$raw_edge_weights)))
  expect_equal(mean(reweighted$edge_weights), 1, tolerance = 1e-10)
})

test_that("cylinder surface embedding returns finite 3D coordinates", {
  coords <- grip:::cylinder.surface.embedding(5, 8, surface = "barrel", amplitude = 0.25)

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
  coords <- grip:::torus.surface.embedding(6, 9, surface = "standard")

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
  coords <- grip:::sphere.surface.embedding(6, 9, surface = "standard")

  expect_true(is.matrix(coords))
  expect_true(is.numeric(coords))
  expect_equal(dim(coords), c(38L, 3L))
  expect_true(all(is.finite(coords)))
  expect_equal(colnames(coords), c("x", "y", "z"))
})

test_that("sphere surface graph returns normalized positive edge weights", {
  spec <- sphere.surface.graph(6, 9, surface = "ellipsoid", amplitude = 0.2)

  expect_s3_class(spec, "grip_sphere_surface_graph")
  expect_equal(spec$edges, grip:::edges.sphere(6, 9))
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

test_that("triangulated polyhedron surface embedding returns finite 3D coordinates", {
  coords <- grip:::triangulated.polyhedron.surface.embedding(
    base = "icosahedron",
    level = 1,
    surface = "twisted",
    twist = 0.5
  )

  expect_true(is.matrix(coords))
  expect_true(is.numeric(coords))
  expect_equal(dim(coords), c(42L, 3L))
  expect_true(all(is.finite(coords)))
  expect_equal(colnames(coords), c("x", "y", "z"))
})

test_that("triangulated polyhedron surface graph returns normalized positive edge weights", {
  spec <- triangulated.polyhedron.surface.graph(
    base = "octahedron",
    level = 2,
    surface = "inflated",
    amplitude = 0.6
  )

  expect_s3_class(spec, "grip_triangulated_polyhedron_surface_graph")
  expect_equal(spec$edges, grip:::edges.triangulated.polyhedron("octahedron", 2))
  expect_equal(spec$n, max(spec$edges))
  expect_equal(spec$base, "octahedron")
  expect_equal(spec$level, 2L)
  expect_equal(spec$subdivision, 4L)
  expect_true(all(is.finite(spec$edge_weights)))
  expect_true(all(spec$edge_weights > 0))
  expect_gt(stats::sd(spec$edge_weights), 0)
  expect_equal(stats::median(spec$edge_weights), 1, tolerance = 1e-10)
  expect_equal(dim(spec$coords_surface), c(spec$n, 3L))
  expect_equal(dim(spec$coords_param), c(spec$n, 3L))
  expect_equal(spec$family, "triangulated.polyhedron")
  expect_equal(spec$surface, "inflated")
})

test_that("wavy triangulated polyhedron graph supports alternate normalization", {
  spec <- triangulated.polyhedron.surface.graph(
    base = "tetrahedron",
    level = 1,
    surface = "wavy",
    amplitude = 0.18,
    freq = 1.5,
    normalize = "mean"
  )

  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_gt(max(spec$edge_weights) - min(spec$edge_weights), 1e-6)
})

test_that("triangulated annulus surface embedding returns finite 3D coordinates", {
  coords <- grip:::triangulated.annulus.surface.embedding(
    resolution = 7,
    surface = "folded",
    amplitude = 0.4
  )

  expect_true(is.matrix(coords))
  expect_true(is.numeric(coords))
  expect_true(all(is.finite(coords)))
  expect_equal(ncol(coords), 3L)
  expect_equal(colnames(coords), c("x", "y", "z"))
})

test_that("triangulated annulus surface graph returns normalized positive edge weights", {
  spec <- triangulated.annulus.surface.graph(
    resolution = 7,
    surface = "ripple",
    amplitude = 0.35,
    freq_u = 1.2,
    freq_v = 0.8,
    normalize = "mean"
  )

  expect_s3_class(spec, "grip_triangulated_annulus_surface_graph")
  expect_equal(spec$edges, grip:::edges.triangulated.annulus(resolution = 7))
  expect_equal(spec$n, max(spec$edges))
  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_gt(max(spec$edge_weights) - min(spec$edge_weights), 1e-6)
  expect_equal(spec$family, "triangulated.annulus")
  expect_equal(spec$surface, "ripple")
  expect_equal(spec$resolution, 7L)
  expect_equal(dim(spec$coords_surface), c(spec$n, 3L))
  expect_equal(dim(spec$coords_param), c(spec$n, 2L))
})

test_that("triangulated pair-of-pants surface graph returns normalized positive edge weights", {
  spec <- triangulated.pair.of.pants.surface.graph(
    resolution = 7,
    surface = "saddle",
    amplitude = 0.45
  )

  expect_s3_class(spec, "grip_triangulated_pair_of_pants_surface_graph")
  expect_equal(spec$edges, grip:::edges.triangulated.pair.of.pants(resolution = 7))
  expect_equal(spec$n, max(spec$edges))
  expect_true(all(spec$edge_weights > 0))
  expect_equal(stats::median(spec$edge_weights), 1, tolerance = 1e-10)
  expect_gt(stats::sd(spec$edge_weights), 0)
  expect_equal(spec$family, "triangulated.pair.of.pants")
  expect_equal(spec$surface, "saddle")
  expect_equal(spec$resolution, 7L)
  expect_equal(dim(spec$coords_surface), c(spec$n, 3L))
  expect_equal(dim(spec$coords_param), c(spec$n, 2L))
})

test_that("irregular annulus surface graph returns normalized positive edge weights", {
  spec <- irregular.annulus.surface.graph(
    rings = 6,
    outer_count = 24,
    surface = "ripple",
    amplitude = 0.3,
    freq_u = 1.3,
    freq_v = 0.7,
    normalize = "mean"
  )

  expect_s3_class(spec, "grip_irregular_annulus_surface_graph")
  expect_equal(spec$edges, grip:::edges.irregular.annulus(rings = 6, outer_count = 24))
  expect_equal(spec$n, max(spec$edges))
  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_gt(max(spec$edge_weights) - min(spec$edge_weights), 1e-6)
  expect_equal(spec$family, "irregular.annulus")
  expect_equal(spec$surface, "ripple")
  expect_equal(spec$rings, 6L)
  expect_equal(length(spec$ring_sizes), 6L)
  expect_equal(dim(spec$coords_surface), c(spec$n, 3L))
  expect_equal(dim(spec$coords_param), c(spec$n, 2L))
})

test_that("irregular ball solid graph returns normalized positive edge weights", {
  spec <- irregular.ball.solid.graph(
    base = "icosahedron",
    level = 1,
    layers = 3,
    surface = "wavy",
    amplitude = 0.12,
    freq_theta = 2,
    freq_phi = 3,
    normalize = "mean"
  )

  expect_s3_class(spec, "grip_irregular_ball_solid_graph")
  expect_equal(spec$edges, grip:::edges.irregular.ball(base = "icosahedron", level = 1, layers = 3))
  expect_equal(spec$n, max(spec$edges))
  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_gt(max(spec$edge_weights) - min(spec$edge_weights), 1e-6)
  expect_equal(spec$family, "irregular.ball")
  expect_equal(spec$surface, "wavy")
  expect_equal(spec$base, "icosahedron")
  expect_equal(spec$level, 1L)
  expect_equal(spec$layers, 3L)
  expect_equal(length(spec$radii), 3L)
  expect_equal(dim(spec$coords_surface), c(spec$n, 3L))
  expect_equal(dim(spec$coords_param), c(spec$n, 3L))
})

test_that("irregular shell solid graph returns normalized positive edge weights", {
  spec <- irregular.shell.solid.graph(
    base = "octahedron",
    level = 1,
    layers = 3,
    inner_radius = 0.42,
    surface = "bulged",
    amplitude = 0.1,
    normalize = "mean"
  )

  expect_s3_class(spec, "grip_irregular_shell_solid_graph")
  expect_equal(spec$edges, grip:::edges.irregular.shell(base = "octahedron", level = 1, layers = 3,
                                                 inner_radius = 0.42))
  expect_equal(spec$n, max(spec$edges))
  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_gt(max(spec$edge_weights) - min(spec$edge_weights), 1e-6)
  expect_equal(spec$family, "irregular.shell")
  expect_equal(spec$surface, "bulged")
  expect_equal(spec$base, "octahedron")
  expect_equal(spec$level, 1L)
  expect_equal(spec$layers, 3L)
  expect_equal(length(spec$radii), 3L)
  expect_equal(dim(spec$coords_surface), c(spec$n, 3L))
  expect_equal(dim(spec$coords_param), c(spec$n, 3L))
})

test_that("irregular torus surface graph returns normalized positive edge weights", {
  spec <- irregular.torus.surface.graph(
    major_rings = 8,
    tube_count = 16,
    surface = "wavy",
    amplitude = 0.18,
    freq_major = 2,
    freq_minor = 1.5,
    normalize = "mean"
  )

  expect_s3_class(spec, "grip_irregular_torus_surface_graph")
  expect_equal(spec$edges, grip:::edges.irregular.torus(major_rings = 8, tube_count = 16))
  expect_equal(spec$n, max(spec$edges))
  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_gt(max(spec$edge_weights) - min(spec$edge_weights), 1e-6)
  expect_equal(spec$family, "irregular.torus")
  expect_equal(spec$surface, "wavy")
  expect_equal(spec$major_rings, 8L)
  expect_equal(length(spec$ring_sizes), 8L)
  expect_equal(dim(spec$coords_surface), c(spec$n, 3L))
  expect_equal(dim(spec$coords_param), c(spec$n, 2L))
})

test_that("irregular pair-of-pants surface graph returns normalized positive edge weights", {
  spec <- irregular.pair.of.pants.surface.graph(
    slices = 11,
    outer_count = 24,
    surface = "ripple",
    amplitude = 0.28,
    freq_u = 1.1,
    freq_v = 0.8,
    normalize = "mean"
  )

  expect_s3_class(spec, "grip_irregular_pair_of_pants_surface_graph")
  expect_equal(spec$edges, grip:::edges.irregular.pair.of.pants(slices = 11, outer_count = 24))
  expect_equal(spec$n, max(spec$edges))
  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_gt(max(spec$edge_weights) - min(spec$edge_weights), 1e-6)
  expect_equal(spec$family, "irregular.pair.of.pants")
  expect_equal(spec$surface, "ripple")
  expect_equal(spec$slices, 11L)
  expect_equal(length(spec$slice_sizes), 11L)
  expect_equal(length(spec$slice_components), 11L)
  expect_true(all(spec$slice_components >= 1L))
  expect_equal(dim(spec$coords_surface), c(spec$n, 3L))
  expect_equal(dim(spec$coords_param), c(spec$n, 2L))
})

test_that("irregular double torus surface graph returns normalized positive edge weights", {
  spec <- irregular.double.torus.surface.graph(
    slices = 11,
    tube_count = 14,
    surface = "wavy",
    amplitude = 0.14,
    freq_x = 2,
    freq_theta = 2.5,
    normalize = "mean"
  )

  expect_s3_class(spec, "grip_irregular_double_torus_surface_graph")
  expect_equal(spec$edges, grip:::edges.irregular.double.torus(slices = 11, tube_count = 14))
  expect_equal(spec$n, max(spec$edges))
  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_gt(max(spec$edge_weights) - min(spec$edge_weights), 1e-6)
  expect_equal(spec$family, "irregular.double.torus")
  expect_equal(spec$surface, "wavy")
  expect_equal(spec$slices, 11L)
  expect_equal(length(spec$slice_sizes), 11L)
  expect_equal(length(spec$slice_components), 11L)
  expect_equal(dim(spec$coords_surface), c(spec$n, 3L))
  expect_equal(dim(spec$coords_param), c(spec$n, 3L))
})

test_that("irregular sphere surface graph returns normalized positive edge weights", {
  spec <- irregular.sphere.surface.graph(
    bands = 6,
    equator_count = 24,
    surface = "wavy",
    amplitude = 0.18,
    freq_theta = 3,
    freq_lat = 2,
    normalize = "mean"
  )

  expect_s3_class(spec, "grip_irregular_sphere_surface_graph")
  expect_equal(spec$edges, grip:::edges.irregular.sphere(bands = 6, equator_count = 24))
  expect_equal(spec$n, max(spec$edges))
  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_gt(max(spec$edge_weights) - min(spec$edge_weights), 1e-6)
  expect_equal(spec$family, "irregular.sphere")
  expect_equal(spec$surface, "wavy")
  expect_equal(spec$bands, 6L)
  expect_equal(length(spec$band_sizes), 6L)
  expect_equal(dim(spec$coords_surface), c(spec$n, 3L))
  expect_equal(dim(spec$coords_param), c(spec$n, 2L))
})

test_that("sierpinski triangle surface embedding returns finite 3D coordinates", {
  coords <- grip:::sierpinski.triangle.surface.embedding(2, surface = "folded", amplitude = 0.8)

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
  coords <- grip:::sierpinski.tetrahedron.surface.embedding(2, surface = "twisted", twist = 0.7)

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

test_that("recursive cube mask surface embedding returns finite 3D coordinates", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  mask[1, 1, 2] <- FALSE
  coords <- grip:::recursive.cube.mask.surface.embedding(
    mask = mask,
    level = 2,
    surface = "twisted",
    twist = 0.5
  )

  expect_true(is.matrix(coords))
  expect_true(is.numeric(coords))
  expect_true(all(is.finite(coords)))
  expect_equal(ncol(coords), 3L)
  expect_equal(colnames(coords), c("x", "y", "z"))
})

test_that("recursive cube mask surface graph returns normalized positive edge weights", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  mask[2, 1, 2] <- FALSE
  spec <- recursive.cube.mask.surface.graph(
    mask = mask,
    level = 2,
    surface = "bulged",
    amplitude = 0.18
  )

  expect_s3_class(spec, "grip_recursive_cube_mask_surface_graph")
  expect_equal(spec$edges, grip:::edges.recursive.cube.mask(mask, 2))
  expect_equal(spec$n, max(spec$edges))
  expect_equal(spec$mask, mask != 0)
  expect_equal(spec$mask_side, 2L)
  expect_equal(spec$side, 4L)
  expect_true(all(is.finite(spec$edge_weights)))
  expect_true(all(spec$edge_weights > 0))
  expect_gt(stats::sd(spec$edge_weights), 0)
  expect_equal(stats::median(spec$edge_weights), 1, tolerance = 1e-10)
  expect_equal(dim(spec$coords_surface), c(spec$n, 3L))
  expect_equal(dim(spec$coords_param), c(spec$n, 3L))
  expect_equal(spec$family, "recursive.cube.mask")
  expect_equal(spec$surface, "bulged")
})

test_that("menger sponge surface graph returns normalized positive edge weights", {
  spec <- menger.sponge.surface.graph(
    2,
    surface = "wavy",
    amplitude = 0.14,
    freq = 1.5,
    normalize = "mean"
  )

  expect_s3_class(spec, "grip_menger_sponge_surface_graph")
  expect_equal(spec$edges, grip:::edges.menger.sponge(2))
  expect_equal(spec$n, 400L)
  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_gt(max(spec$edge_weights) - min(spec$edge_weights), 1e-6)
  expect_equal(spec$family, "menger.sponge")
  expect_equal(spec$surface, "wavy")
})

test_that("periodic tunnel cube surface graph returns normalized positive edge weights", {
  spec <- cube.periodic.tunnels.surface.graph(
    level = 1,
    side = 5,
    surface = "wavy",
    amplitude = 0.16,
    freq = 1.5,
    normalize = "mean"
  )

  expect_s3_class(spec, "grip_cube_periodic_tunnels_surface_graph")
  expect_equal(spec$edges, grip:::edges.cube.periodic.tunnels(level = 1, side = 5))
  expect_equal(spec$n, 81L)
  expect_equal(mean(spec$edge_weights), 1, tolerance = 1e-10)
  expect_gt(max(spec$edge_weights) - min(spec$edge_weights), 1e-6)
  expect_equal(spec$family, "cube.periodic.tunnels")
  expect_equal(spec$surface, "wavy")
})

test_that("asymmetric cavity cube surface graph returns normalized positive edge weights", {
  spec <- cube.asymmetric.cavities.surface.graph(
    level = 1,
    side = 5,
    surface = "bulged",
    amplitude = 0.18
  )

  expect_s3_class(spec, "grip_cube_asymmetric_cavities_surface_graph")
  expect_equal(spec$edges, grip:::edges.cube.asymmetric.cavities(level = 1, side = 5))
  expect_equal(spec$n, 116L)
  expect_true(all(spec$edge_weights > 0))
  expect_equal(stats::median(spec$edge_weights), 1, tolerance = 1e-10)
  expect_equal(spec$family, "cube.asymmetric.cavities")
  expect_equal(spec$surface, "bulged")
})

test_that("channel-network cube surface graph returns normalized positive edge weights", {
  spec <- cube.channel.network.surface.graph(
    level = 1,
    side = 5,
    surface = "twisted",
    twist = 0.45
  )

  expect_s3_class(spec, "grip_cube_channel_network_surface_graph")
  expect_equal(spec$edges, grip:::edges.cube.channel.network(level = 1, side = 5))
  expect_equal(spec$n, 104L)
  expect_true(all(spec$edge_weights > 0))
  expect_equal(stats::median(spec$edge_weights), 1, tolerance = 1e-10)
  expect_equal(spec$family, "cube.channel.network")
  expect_equal(spec$surface, "twisted")
})

test_that("recursive tetrahedron mask surface embedding returns finite 3D coordinates", {
  coords <- grip:::recursive.tetrahedron.mask.surface.embedding(
    mask = mask.tetrahedron.corner.missing("base_right"),
    level = 2,
    surface = "twisted",
    twist = 0.65
  )

  expect_true(is.matrix(coords))
  expect_true(is.numeric(coords))
  expect_true(all(is.finite(coords)))
  expect_equal(ncol(coords), 3L)
  expect_equal(colnames(coords), c("x", "y", "z"))
})

test_that("recursive tetrahedron mask surface graph returns normalized positive edge weights", {
  corner_mask <- mask.tetrahedron.corner.missing("apex")
  spec <- recursive.tetrahedron.mask.surface.graph(
    mask = corner_mask,
    level = 2,
    surface = "wavy",
    amplitude = 0.22,
    freq = 2.5
  )

  expect_s3_class(spec, "grip_recursive_tetrahedron_mask_surface_graph")
  expect_equal(spec$edges, grip:::edges.recursive.tetrahedron.mask(corner_mask, 2))
  expect_equal(spec$n, max(spec$edges))
  expect_equal(spec$mask, corner_mask)
  expect_equal(length(spec$edge_weights), nrow(spec$edges))
  expect_true(all(is.finite(spec$edge_weights)))
  expect_true(all(spec$edge_weights > 0))
  expect_gt(stats::sd(spec$edge_weights), 0)
  expect_equal(stats::median(spec$edge_weights), 1, tolerance = 1e-10)
  expect_equal(dim(spec$coords_surface), c(spec$n, 3L))
  expect_equal(dim(spec$coords_param), c(spec$n, 3L))
  expect_equal(spec$family, "recursive.tetrahedron.mask")
  expect_equal(spec$surface, "wavy")
})

test_that("recursive triangle mask surface embedding returns finite 3D coordinates", {
  coords <- grip:::recursive.triangle.mask.surface.embedding(
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
  expect_equal(spec$edges, grip:::edges.recursive.triangle.mask(bridge_mask, 2))
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
  expect_equal(spec$edges, grip:::edges.vicsek(2))
  expect_equal(spec$n, 25L)
  expect_equal(spec$mask, vicsek_mask != 0)
  expect_equal(spec$mask_size, 3L)
  expect_equal(spec$side, 9L)
  expect_true(all(spec$edge_weights > 0))
  expect_equal(stats::median(spec$edge_weights), 1, tolerance = 1e-10)
})

test_that("sierpinski carpet surface embedding returns finite 3D coordinates", {
  coords <- grip:::sierpinski.carpet.surface.embedding(2, surface = "saddle", amplitude = 0.8)

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
  expect_equal(spec$edges, grip:::edges.vicsek(2))
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
  expect_equal(spec$edges, grip:::edges.occupied.mesh(keep))
  expect_equal(spec$n, sum(keep))
  expect_equal(length(spec$edge_weights), nrow(spec$edges))
  expect_true(all(spec$edge_weights > 0))
  expect_equal(stats::median(spec$edge_weights), 1, tolerance = 1e-10)
  expect_equal(dim(spec$coords_param), c(sum(keep), 2L))
  expect_equal(spec$family, "occupied.mesh")
  expect_equal(spec$connectivity, "orthogonal")
})

test_that("occupied mesh surface graph supports diagonal connectivity", {
  keep <- matrix(
    c(
      TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE,
      TRUE, TRUE, FALSE
    ),
    nrow = 3,
    byrow = TRUE
  )
  spec <- occupied.mesh.surface.graph(
    keep,
    surface = "saddle",
    amplitude = 0.5,
    connectivity = "diagonal"
  )

  expect_equal(spec$edges, grip:::edges.occupied.mesh(keep, connectivity = "diagonal"))
  expect_equal(spec$connectivity, "diagonal")
  expect_true(nrow(spec$edges) > nrow(grip:::edges.occupied.mesh(keep, connectivity = "orthogonal")))
  expect_true(all(spec$edge_weights > 0))
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
