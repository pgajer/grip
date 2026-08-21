# Examples for exported functions whose primary roxygen blocks are distributed
# across the implementation files. Keeping these additions together makes the
# public-API example audit straightforward without moving function docs away
# from their implementations.

#' @name build.weighted.misf
#' @rdname build.weighted.misf
#' @examples
#' cycle_edges <- edges.cycle(6)
#' cycle_lengths <- seq(0.8, 1.3, length.out = nrow(cycle_edges))
#' levels <- build.weighted.misf(
#'   edges = cycle_edges, n = 6, edge_weights = cycle_lengths,
#'   num_init = 3, num_nbrs = 2, seed = 1
#' )
NULL

#' @name porous_cube_surface_helpers
#' @rdname porous_cube_surface_helpers
#' @examples
#' asymmetric <- cube.asymmetric.cavities.surface.graph(
#'   level = 1, side = 5, cavity_size = 1, pocket_size = 1
#' )
#' channels <- cube.channel.network.surface.graph(
#'   level = 1, side = 5, channel_width = 1, branch_offset = 2
#' )
#' tunnels <- cube.periodic.tunnels.surface.graph(
#'   level = 1, side = 5, tunnel_width = 1,
#'   tunnel_period = 2, tunnel_offset = 1
#' )
NULL

#' @name cylinder_surface_helpers
#' @rdname cylinder_surface_helpers
#' @examples
#' cylinder <- cylinder.surface.graph(4, 5, surface = "hourglass")
NULL

#' @name edge.kk
#' @rdname edge.kk
#' @examples
#' cycle_edges <- edges.cycle(6)
#' theta <- seq(0, 2 * pi, length.out = 7)[-7]
#' initial <- cbind(cos(theta), sin(theta))
#' prepared <- prepare.edge.kk(cycle_edges, n = 6)
#' fit <- edge.kk(
#'   coords = initial, prepared = prepared,
#'   max_iter = 1, density_n = 32, engine = "R"
#' )
NULL

#' @name edge.length.density.stiffness
#' @rdname edge.length.density.stiffness
#' @examples
#' stiffness <- edge.length.density.stiffness(
#'   c(0.8, 0.9, 1, 1.1, 1.2, 1.3), density_n = 32
#' )
NULL

#' @name edge.repulsive.stage
#' @rdname edge.repulsive.stage
#' @examples
#' cycle_edges <- edges.cycle(6)
#' theta <- seq(0, 2 * pi, length.out = 7)[-7]
#' initial <- cbind(cos(theta), sin(theta))
#' fit <- edge.repulsive.stage(
#'   initial, cycle_edges, rep(1, nrow(cycle_edges)),
#'   lambda = 0.1, max.iter = 1, engine = "R"
#' )
NULL

#' @name edge.repulsive.state
#' @rdname edge.repulsive.state
#' @examples
#' cycle_edges <- edges.cycle(6)
#' theta <- seq(0, 2 * pi, length.out = 7)[-7]
#' coords <- cbind(cos(theta), sin(theta))
#' state <- edge.repulsive.state(
#'   coords, cycle_edges, rep(1, nrow(cycle_edges)),
#'   lambda = 0.1, engine = "R"
#' )
NULL

#' @name graph_generators
#' @rdname graph_generators
#' @examples
#' cube_edges <- edges.cube(3)
#' cycle_edges <- edges.cycle(6)
#' cylinder_edges <- edges.cylinder(3, 4)
#' tree_edges <- edges.kary.tree(k = 2, depth = 2)
#' mesh_edges <- edges.mesh(3, 4)
#' carpet_edges <- edges.sierpinski.carpet(level = 1)
#' tetrahedron_edges <- edges.sierpinski.tetrahedron(level = 1)
#' torus_edges <- edges.torus(3, 4)
NULL

#' @name geodesic.kk
#' @rdname geodesic.kk
#' @examples
#' cycle_edges <- edges.cycle(6)
#' theta <- seq(0, 2 * pi, length.out = 7)[-7]
#' initial <- cbind(cos(theta), sin(theta))
#' prepared <- prepare.geodesic.kk(cycle_edges, n = 6)
#' fit <- geodesic.kk(initial, prepared = prepared, max_iter = 1)
NULL

#' @name geometry.diagnostics
#' @rdname geometry.diagnostics
#' @examples
#' cycle_edges <- edges.cycle(6)
#' theta <- seq(0, 2 * pi, length.out = 7)[-7]
#' target <- cbind(cos(theta), sin(theta))
#' diagnostics <- geometry.diagnostics(
#'   target + 0.01, target, cycle_edges,
#'   sample.size.symmetry = 10, sample.size.wedges = 10
#' )
NULL

#' @name globalrep.weighted.grip
#' @rdname globalrep.weighted.grip
#' @examples
#' cycle_edges <- edges.cycle(6)
#' cycle_lengths <- seq(0.8, 1.3, length.out = nrow(cycle_edges))
#' layout <- globalrep.weighted.grip(
#'   edges = cycle_edges, n = 6, edge_weights = cycle_lengths,
#'   dim = 2, rounds = 1, final_rounds = 1,
#'   num_init = 3, num_nbrs = 2, seed = 1
#' )
NULL

#' @name graph.riemannian.star.structure
#' @rdname graph.riemannian.star.structure
#' @examples
#' cycle_edges <- edges.cycle(6)
#' theta <- seq(0, 2 * pi, length.out = 7)[-7]
#' coords <- cbind(cos(theta), sin(theta))
#' prepared <- prepare.graph.geodesic.mds(cycle_edges, n = 6)
#' star <- graph.riemannian.star.structure(X = coords, prepared = prepared)
NULL

#' @name irregular_annulus_surface_helpers
#' @rdname irregular_annulus_surface_helpers
#' @examples
#' annulus <- irregular.annulus.surface.graph(rings = 3, outer_count = 12)
NULL

#' @name irregular_ball_solid_helpers
#' @rdname irregular_ball_solid_helpers
#' @examples
#' ball <- irregular.ball.solid.graph(
#'   base = "tetrahedron", level = 0, layers = 2
#' )
NULL

#' @name irregular_double_torus_surface_helpers
#' @rdname irregular_double_torus_surface_helpers
#' @examples
#' double_torus <- irregular.double.torus.surface.graph(
#'   slices = 7, tube_count = 8
#' )
NULL

#' @name irregular_pair_of_pants_surface_helpers
#' @rdname irregular_pair_of_pants_surface_helpers
#' @examples
#' pair_of_pants <- irregular.pair.of.pants.surface.graph(
#'   slices = 5, outer_count = 12
#' )
NULL

#' @name irregular_rectangle_surface_helpers
#' @rdname irregular_rectangle_surface_helpers
#' @examples
#' rectangle <- irregular.rectangle.surface.graph(4, 5, surface = "ripple")
NULL

#' @name irregular_shell_solid_helpers
#' @rdname irregular_shell_solid_helpers
#' @examples
#' shell <- irregular.shell.solid.graph(
#'   base = "tetrahedron", level = 0, layers = 2
#' )
NULL

#' @name irregular_sphere_surface_helpers
#' @rdname irregular_sphere_surface_helpers
#' @examples
#' sphere <- irregular.sphere.surface.graph(bands = 3, equator_count = 12)
NULL

#' @name irregular_torus_surface_helpers
#' @rdname irregular_torus_surface_helpers
#' @examples
#' torus <- irregular.torus.surface.graph(major_rings = 4, tube_count = 8)
NULL

#' @name kary_tree_weighted_graph_helpers
#' @rdname kary_tree_weighted_graph_helpers
#' @examples
#' tree <- kary.tree.weighted.graph(k = 2, depth = 2)
NULL

#' @name perforated_grid_helpers
#' @rdname perforated_grid_helpers
#' @examples
#' notches <- keep.asymmetric.notches(7, 8, notch_depth = 2, notch_width = 2)
#' holes <- keep.periodic.holes(7, 8, hole_period = 3)
#' slits <- keep.slit.channels(7, 8, slit_period = 3)
#' windows <- keep.staggered.windows(7, 8, row_period = 3, col_period = 4)
NULL

#' @name kernel.gram.gkk
#' @rdname kernel.gram.gkk
#' @examples
#' cycle_edges <- edges.cycle(6)
#' theta <- seq(0, 2 * pi, length.out = 7)[-7]
#' initial <- cbind(cos(theta), sin(theta))
#' prepared <- prepare.graph.geodesic.mds(cycle_edges, n = 6)
#' fit <- kernel.gram.gkk(
#'   coords = initial, prepared = prepared, X = initial,
#'   max_iter = 1, density_n = 32, engine = "R"
#' )
NULL

#' @name landmark.geodesic.kk
#' @rdname landmark.geodesic.kk
#' @examples
#' cycle_edges <- edges.cycle(6)
#' theta <- seq(0, 2 * pi, length.out = 7)[-7]
#' initial <- cbind(cos(theta), sin(theta))
#' prepared <- prepare.landmark.geodesic.kk(
#'   cycle_edges, n = 6, local_nbrs = 2, landmark_count = 2
#' )
#' fit <- landmark.geodesic.kk(initial, prepared = prepared, max_iter = 1)
NULL

#' @name mask_pattern_helpers
#' @rdname mask_pattern_helpers
#' @examples
#' asymmetric <- mask.asymmetric.holes(k = 5, hole_size = 1)
#' border <- mask.border(k = 5, thickness = 1)
#' corner <- mask.corner(k = 5, width = 2, corner = "top_left")
#' cross <- mask.cross(k = 5, arm_width = 1)
NULL

#' @name cube_mask_pattern_helpers
#' @rdname cube_mask_pattern_helpers
#' @examples
#' cavities <- mask.cube.asymmetric.cavities(
#'   side = 5, cavity_size = 1, pocket_size = 1
#' )
#' channels <- mask.cube.channel.network(
#'   side = 5, channel_width = 1, branch_offset = 2
#' )
#' tunnels <- mask.cube.periodic.tunnels(
#'   side = 5, tunnel_width = 1, tunnel_period = 2, tunnel_offset = 1
#' )
NULL

#' @name tetrahedron_mask_helpers
#' @rdname tetrahedron_mask_helpers
#' @examples
#' classic <- mask.tetrahedron.classic()
#' missing_apex <- mask.tetrahedron.corner.missing("apex")
NULL

#' @name triangle_mask_helpers
#' @rdname triangle_mask_helpers
#' @examples
#' classic <- mask.triangle.classic()
#' bridge <- mask.triangle.bridge("top")
NULL

#' @name menger_sponge_surface_helpers
#' @rdname menger_sponge_surface_helpers
#' @examples
#' sponge <- menger.sponge.surface.graph(level = 1)
NULL

#' @name mesh_surface_helpers
#' @rdname mesh_surface_helpers
#' @examples
#' mesh <- mesh.surface.graph(4, 5, surface = "saddle")
NULL

#' @name metric.mds
#' @rdname metric.mds
#' @examples
#' cycle_edges <- edges.cycle(6)
#' prepared <- prepare.graph.geodesic.mds(cycle_edges, n = 6)
#' fit <- metric.mds(prepared = prepared, dim = 2)
NULL

#' @name misf.geodesic.kk
#' @rdname misf.geodesic.kk
#' @examples
#' mesh_edges <- edges.mesh(4, 4)
#' fit <- misf.geodesic.kk(
#'   edges = mesh_edges, n = 16, num_init = 4, dim = 2,
#'   top_level_restarts = 2, top_level_max_iter = 3,
#'   insertion_mode = "geodesic", insertion_max_iter = 6,
#'   refinement_pair_mode = "auto", refinement_full_limit = 4,
#'   refinement_max_iter = 3, final_pair_mode = "landmark",
#'   final_full_limit = 4, final_max_iter = 3,
#'   return_trace = TRUE, seed = 1
#' )
NULL

#' @name occupied_mesh_surface_helpers
#' @rdname occupied_mesh_surface_helpers
#' @examples
#' keep <- matrix(TRUE, nrow = 4, ncol = 5)
#' keep[2, 3] <- FALSE
#' occupied <- occupied.mesh.surface.graph(keep, surface = "paraboloid")
NULL

#' @name prepare.edge.kk
#' @rdname prepare.edge.kk
#' @examples
#' cycle_edges <- edges.cycle(6)
#' prepared <- prepare.edge.kk(cycle_edges, n = 6)
NULL

#' @name prepare.geodesic.kk
#' @rdname prepare.geodesic.kk
#' @examples
#' cycle_edges <- edges.cycle(6)
#' prepared <- prepare.geodesic.kk(cycle_edges, n = 6)
NULL

#' @name prepare.graph.geodesic.mds
#' @rdname prepare.graph.geodesic.mds
#' @examples
#' cycle_edges <- edges.cycle(6)
#' prepared <- prepare.graph.geodesic.mds(cycle_edges, n = 6)
NULL

#' @name prepare.landmark.geodesic.kk
#' @rdname prepare.landmark.geodesic.kk
#' @examples
#' cycle_edges <- edges.cycle(6)
#' prepared <- prepare.landmark.geodesic.kk(
#'   cycle_edges, n = 6, local_nbrs = 2, landmark_count = 2
#' )
NULL

#' @name recursive_cube_mask_surface_helpers
#' @rdname recursive_cube_mask_surface_helpers
#' @examples
#' mask <- mask.cube.channel.network(
#'   side = 5, channel_width = 1, branch_offset = 2
#' )
#' cube <- recursive.cube.mask.surface.graph(mask, level = 1)
NULL

#' @name recursive_mask_grid_surface_helpers
#' @rdname recursive_mask_grid_surface_helpers
#' @examples
#' grid <- recursive.mask.grid.surface.graph(mask.cross(3), level = 1)
NULL

#' @name recursive_tetrahedron_mask_surface_helpers
#' @rdname recursive_tetrahedron_mask_surface_helpers
#' @examples
#' tetrahedron <- recursive.tetrahedron.mask.surface.graph(level = 1)
NULL

#' @name recursive_triangle_mask_surface_helpers
#' @rdname recursive_triangle_mask_surface_helpers
#' @examples
#' triangle <- recursive.triangle.mask.surface.graph(level = 1)
NULL

#' @name repulsive.stage
#' @rdname repulsive.stage
#' @examples
#' theta <- seq(0, 2 * pi, length.out = 7)[-7]
#' initial <- cbind(cos(theta), sin(theta))
#' fit <- repulsive.stage(
#'   initial, lambda = 0.1, max.iter = 1, engine = "R"
#' )
NULL

#' @name repulsive.state
#' @rdname repulsive.state
#' @examples
#' theta <- seq(0, 2 * pi, length.out = 7)[-7]
#' coords <- cbind(cos(theta), sin(theta))
#' state <- repulsive.state(coords, lambda = 0.1, engine = "R")
NULL

#' @name sampled_rectangle_surface_helpers
#' @rdname sampled_rectangle_surface_helpers
#' @examples
#' rectangle <- sampled.rectangle.surface.graph(n = 20, k = 4, seed = 1)
#' sequence <- sampled.rectangle.surface.graphs(
#'   n = 20, k = c(3, 4), seed = 1
#' )
NULL

#' @name score.geodesic.kk
#' @rdname score.geodesic.kk
#' @examples
#' cycle_edges <- edges.cycle(6)
#' theta <- seq(0, 2 * pi, length.out = 7)[-7]
#' coords <- cbind(cos(theta), sin(theta))
#' prepared <- prepare.geodesic.kk(cycle_edges, n = 6)
#' score <- score.geodesic.kk(coords, prepared = prepared)
NULL

#' @name score.gmds
#' @rdname score.gmds
#' @examples
#' cycle_edges <- edges.cycle(6)
#' theta <- seq(0, 2 * pi, length.out = 7)[-7]
#' coords <- cbind(cos(theta), sin(theta))
#' prepared <- prepare.graph.geodesic.mds(cycle_edges, n = 6)
#' score <- score.gmds(coords, prepared = prepared)
NULL

#' @name score.landmark.geodesic.kk
#' @rdname score.landmark.geodesic.kk
#' @examples
#' cycle_edges <- edges.cycle(6)
#' theta <- seq(0, 2 * pi, length.out = 7)[-7]
#' coords <- cbind(cos(theta), sin(theta))
#' prepared <- prepare.landmark.geodesic.kk(
#'   cycle_edges, n = 6, local_nbrs = 2, landmark_count = 2
#' )
#' score <- score.landmark.geodesic.kk(coords, prepared = prepared)
NULL

#' @name score.misf.geodesic.kk
#' @rdname score.misf.geodesic.kk
#' @examples
#' mesh_edges <- edges.mesh(4, 4)
#' prepared <- prepare.misf.geodesic.kk(
#'   mesh_edges, n = 16, num_init = 4,
#'   top_level_mode = "skip", seed = 1
#' )
#' coords <- matrix(seq_len(32), nrow = 16, ncol = 2)
#' score <- score.misf.geodesic.kk(coords = coords, prepared = prepared)
NULL

#' @name sierpinski_carpet_surface_helpers
#' @rdname sierpinski_carpet_surface_helpers
#' @examples
#' carpet <- sierpinski.carpet.surface.graph(level = 1)
NULL

#' @name sierpinski_tetrahedron_surface_helpers
#' @rdname sierpinski_tetrahedron_surface_helpers
#' @examples
#' tetrahedron <- sierpinski.tetrahedron.surface.graph(level = 1)
NULL

#' @name sierpinski_triangle_surface_helpers
#' @rdname sierpinski_triangle_surface_helpers
#' @examples
#' triangle <- sierpinski.triangle.surface.graph(level = 1)
NULL

#' @name sphere_surface_helpers
#' @rdname sphere_surface_helpers
#' @examples
#' sphere <- sphere.surface.graph(4, 6, surface = "ellipsoid")
NULL

#' @name torus_surface_helpers
#' @rdname torus_surface_helpers
#' @examples
#' torus <- torus.surface.graph(4, 6, surface = "pinched")
NULL

#' @name triangulated_annulus_surface_helpers
#' @rdname triangulated_annulus_surface_helpers
#' @examples
#' annulus <- triangulated.annulus.surface.graph(resolution = 6)
NULL

#' @name triangulated_pair_of_pants_surface_helpers
#' @rdname triangulated_pair_of_pants_surface_helpers
#' @examples
#' pair_of_pants <- triangulated.pair.of.pants.surface.graph(resolution = 8)
NULL

#' @name triangulated_polyhedron_surface_helpers
#' @rdname triangulated_polyhedron_surface_helpers
#' @examples
#' tetrahedron <- triangulated.polyhedron.surface.graph(
#'   base = "tetrahedron", level = 1
#' )
NULL

#' @name vicsek_surface_helpers
#' @rdname vicsek_surface_helpers
#' @examples
#' vicsek <- vicsek.surface.graph(level = 1)
NULL

#' @name weighted.grip.nd
#' @rdname weighted.grip.nd
#' @examples
#' cycle_edges <- edges.cycle(6)
#' cycle_lengths <- seq(0.8, 1.3, length.out = nrow(cycle_edges))
#' layout <- weighted.grip.nd(
#'   edges = cycle_edges, n = 6, edge_weights = cycle_lengths,
#'   dim = 4, rounds = 1, final_rounds = 1,
#'   num_init = 5, num_nbrs = 2, seed = 1
#' )
NULL
