#' Migrate code to the grip 0.2 API
#'
#' Version 0.2.0 removes the long-form compatibility aliases that were
#' deprecated in the 0.1 series and narrows the exported graph-family API.
#' The canonical replacements are listed below.
#'
#' @section Removed aliases:
#' - `grip.layout()` -> [grip()]
#' - `grip.layout.globalrep()` -> [globalrep.grip()]
#' - `grip.layout.globalrep.weighted()` -> [globalrep.weighted.grip()]
#' - `grip.layout.legacy()` -> [legacy.grip()]
#' - `grip.layout.trace()` -> [trace.grip()]
#' - `grip.layout.trace.legacy()` -> [trace.legacy.grip()]
#' - `grip.layout.trace.weighted()` -> `trace.grip(metric = "edge_length")`
#' - `grip.layout.weighted()` -> `grip(metric = "edge_length")`
#' - `grip.layout.weighted.nd()` -> [weighted.grip.nd()]
#' - `grip.plot()` -> [plot.layout()]
#' - `grip.project.3d()` -> [project.3d()]
#' - `grip.metric.mds.layout()` -> [metric.mds()]
#' - `grip.optimize.kernel.gram.gkk.layout()` -> [kernel.gram.gkk()]
#' - `grip.score.gmds.layout()` -> [score.gmds()]
#' - `grip.gmds.layout.result()` -> `gmds.result()` (now internal)
#' - `grip.edge.length.density.stiffness()` ->
#'   [edge.length.density.stiffness()]
#' - `grip.prepare.landmark.geodesic.kk()` ->
#'   [prepare.landmark.geodesic.kk()]
#' - `grip.prepare.geodesic.kk()` -> [prepare.geodesic.kk()]
#' - `grip.prepare.graph.geodesic.mds()` ->
#'   [prepare.graph.geodesic.mds()]
#' - `grip.prepare.misf.geodesic.kk()` -> [prepare.misf.geodesic.kk()]
#' - `grip.prepare.edge.kk()` -> [prepare.edge.kk()]
#' - `grip.score.landmark.geodesic.kk()` -> [score.landmark.geodesic.kk()]
#' - `grip.score.geodesic.kk()` -> [score.geodesic.kk()]
#' - `grip.score.misf.geodesic.kk()` -> [score.misf.geodesic.kk()]
#' - `grip.optimize.landmark.geodesic.kk()` -> [landmark.geodesic.kk()]
#' - `grip.optimize.geodesic.kk()` -> [geodesic.kk()]
#' - `grip.optimize.misf.geodesic.kk()` -> [misf.geodesic.kk()]
#' - `grip.optimize.edge.kk.layout()` -> [edge.kk()]
#' - `grip.optimize.edge.gkk.layout()` -> [edge.kk()]
#' - `grip.optimize.edge.isometric.layout()` -> [edge.kk()]
#' - `grip.optimize.edge.repulsive.stage()` -> [edge.repulsive.stage()]
#' - `grip.optimize.repulsive.stage()` -> [repulsive.stage()]
#' - `grip.edge.repulsive.state()` -> [edge.repulsive.state()]
#' - `grip.repulsive.state()` -> [repulsive.state()]
#' - `grip.build.misf.weighted()` -> [build.weighted.misf()]
#' - `grip.build.misf()` -> [build.misf()]
#' - `grip.geometry.diagnostics()` -> [geometry.diagnostics()]
#' - `grip.params.from.summary()` -> [params.from.summary()]
#' - `grip.score.layout()` -> [score.layout()]
#' - `grip.compare.layouts()` -> [compare.layouts()]
#'
#' @section Graph-family helpers:
#' Complete `*.surface.graph()` and `*.solid.graph()` constructors remain
#' public and return reusable bundles. Use the bundle's `edges`,
#' `coords_surface`, and `coords_param` components instead of the former
#' specialized edge-list, standalone embedding, and parameter-coordinate
#' exports. Primitive generators such as [edges.path()], [edges.cycle()],
#' [edges.mesh()], and [edges.cube()] remain public.
#'
#' @name grip-0.2-migration
#' @aliases grip.layout
#' @aliases grip.layout.globalrep
#' @aliases grip.layout.globalrep.weighted
#' @aliases grip.layout.legacy
#' @aliases grip.layout.trace
#' @aliases grip.layout.trace.legacy
#' @aliases grip.layout.trace.weighted
#' @aliases grip.layout.weighted
#' @aliases grip.layout.weighted.nd
#' @aliases grip.plot
#' @aliases grip.project.3d
#' @aliases grip.metric.mds.layout
#' @aliases grip.optimize.kernel.gram.gkk.layout
#' @aliases grip.score.gmds.layout
#' @aliases grip.gmds.layout.result
#' @aliases grip.edge.length.density.stiffness
#' @aliases grip.prepare.landmark.geodesic.kk
#' @aliases grip.prepare.geodesic.kk
#' @aliases grip.prepare.graph.geodesic.mds
#' @aliases grip.prepare.misf.geodesic.kk
#' @aliases grip.prepare.edge.kk
#' @aliases grip.score.landmark.geodesic.kk
#' @aliases grip.score.geodesic.kk
#' @aliases grip.score.misf.geodesic.kk
#' @aliases grip.optimize.landmark.geodesic.kk
#' @aliases grip.optimize.geodesic.kk
#' @aliases grip.optimize.misf.geodesic.kk
#' @aliases grip.optimize.edge.kk.layout
#' @aliases grip.optimize.edge.gkk.layout
#' @aliases grip.optimize.edge.isometric.layout
#' @aliases grip.optimize.edge.repulsive.stage
#' @aliases grip.optimize.repulsive.stage
#' @aliases grip.edge.repulsive.state
#' @aliases grip.repulsive.state
#' @aliases grip.build.misf.weighted
#' @aliases grip.build.misf
#' @aliases grip.geometry.diagnostics
#' @aliases grip.params.from.summary
#' @aliases grip.score.layout
#' @aliases grip.compare.layouts
NULL
