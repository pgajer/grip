#' Deprecated long-form API names
#'
#' These long names are deprecated in favor of the shorter canonical API.
#'
#' @param ... Arguments passed to the replacement function.
#' @return The replacement function's return value.
#' @name deprecated-grip-api
NULL

#' @rdname deprecated-grip-api
#' @export
grip.layout <- function(...) {
  .Deprecated("grip")
  grip(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.layout.globalrep <- function(...) {
  .Deprecated("globalrep.grip")
  globalrep.grip(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.layout.globalrep.weighted <- function(...) {
  .Deprecated("globalrep.weighted.grip")
  globalrep.weighted.grip(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.layout.legacy <- function(...) {
  .Deprecated("legacy.grip")
  legacy.grip(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.layout.trace <- function(...) {
  .Deprecated("trace.grip")
  trace.grip(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.layout.trace.legacy <- function(...) {
  .Deprecated("trace.legacy.grip")
  trace.legacy.grip(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.layout.trace.weighted <- function(...) {
  .Deprecated("trace.weighted.grip")
  trace.weighted.grip(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.layout.weighted <- function(...) {
  .Deprecated("weighted.grip")
  weighted.grip(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.layout.weighted.nd <- function(...) {
  .Deprecated("weighted.grip.nd")
  weighted.grip.nd(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.plot <- function(...) {
  .Deprecated("plot.layout")
  plot.layout(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.project.3d <- function(...) {
  .Deprecated("project.3d")
  project.3d(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.metric.mds.layout <- function(...) {
  .Deprecated("metric.mds")
  metric.mds(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.optimize.kernel.gram.gkk.layout <- function(...) {
  .Deprecated("kernel.gram.gkk")
  kernel.gram.gkk(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.score.gmds.layout <- function(...) {
  .Deprecated("score.gmds")
  score.gmds(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.gmds.layout.result <- function(...) {
  .Deprecated("gmds.result")
  gmds.result(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.edge.length.density.stiffness <- function(...) {
  .Deprecated("edge.length.density.stiffness")
  edge.length.density.stiffness(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.prepare.landmark.geodesic.kk <- function(...) {
  .Deprecated("prepare.landmark.geodesic.kk")
  prepare.landmark.geodesic.kk(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.prepare.geodesic.kk <- function(...) {
  .Deprecated("prepare.geodesic.kk")
  prepare.geodesic.kk(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.prepare.graph.geodesic.mds <- function(...) {
  .Deprecated("prepare.graph.geodesic.mds")
  prepare.graph.geodesic.mds(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.prepare.misf.geodesic.kk <- function(...) {
  .Deprecated("prepare.misf.geodesic.kk")
  prepare.misf.geodesic.kk(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.score.landmark.geodesic.kk <- function(...) {
  .Deprecated("score.landmark.geodesic.kk")
  score.landmark.geodesic.kk(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.score.geodesic.kk <- function(...) {
  .Deprecated("score.geodesic.kk")
  score.geodesic.kk(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.score.misf.geodesic.kk <- function(...) {
  .Deprecated("score.misf.geodesic.kk")
  score.misf.geodesic.kk(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.optimize.landmark.geodesic.kk <- function(...) {
  .Deprecated("landmark.geodesic.kk")
  landmark.geodesic.kk(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.optimize.geodesic.kk <- function(...) {
  .Deprecated("geodesic.kk")
  geodesic.kk(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.optimize.misf.geodesic.kk <- function(...) {
  .Deprecated("misf.geodesic.kk")
  misf.geodesic.kk(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.optimize.edge.repulsive.stage <- function(...) {
  .Deprecated("edge.repulsive.stage")
  edge.repulsive.stage(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.optimize.repulsive.stage <- function(...) {
  .Deprecated("repulsive.stage")
  repulsive.stage(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.edge.repulsive.state <- function(...) {
  .Deprecated("edge.repulsive.state")
  edge.repulsive.state(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.repulsive.state <- function(...) {
  .Deprecated("repulsive.state")
  repulsive.state(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.build.misf.weighted <- function(...) {
  .Deprecated("build.weighted.misf")
  build.weighted.misf(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.build.misf <- function(...) {
  .Deprecated("build.misf")
  build.misf(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.geometry.diagnostics <- function(...) {
  .Deprecated("geometry.diagnostics")
  geometry.diagnostics(...)
}

#' @rdname deprecated-grip-api
#' @export
grip.params.from.summary <- function(...) {
  .Deprecated("params.from.summary")
  params.from.summary(...)
}
