knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 5.2
)

if (requireNamespace("rgl", quietly = TRUE)) {
  options(rgl.useNULL = TRUE)
  rgl::setupKnitr(autoprint = TRUE)
}

library(grip)

edges <- edges.cycle(18)

cycle.cmp <- grip.compare.layouts(
  edges = edges,
  n = 18,
  dim = 2,
  candidates = list(
    default = list(),
    circle = list(placement = "circle"),
    circle_refined = list(
      placement = "circle",
      rounds = 64L,
      final_rounds = 128L,
      num_init = 6L,
      num_nbrs = 8L,
      r = 0.05,
      s = 4.5,
      repulsion_factor = 1.25
    ),
    bary_refined = list(
      placement = "barycenter",
      rounds = 64L,
      final_rounds = 128L,
      num_init = 6L,
      num_nbrs = 8L,
      r = 0.05,
      s = 4.5,
      repulsion_factor = 1.25
    )
  ),
  seeds = 1:3,
  return.layouts = TRUE
)

knitr::kable(cycle.cmp$summary[, c(
  "candidate",
  "sampled.stress.mean",
  "edge.length.cv.mean",
  "sampled.nonedge.sep.ratio.mean",
  "score.composite"
)])

cycle.best <- cycle.cmp$summary$candidate[[1L]]
cycle.2d.default <- cycle.cmp$layouts$default[[1L]]
cycle.2d.best <- cycle.cmp$layouts[[cycle.best]][[1L]]

op <- par(mfrow = c(1, 2), mar = c(1.2, 1.2, 3.2, 1.2))
on.exit(par(op), add = TRUE)

grip.plot(cycle.2d.default, edges, main = "Cycle default", pch = 16, cex = 0.7)
grip.plot(cycle.2d.best, edges, main = paste("Cycle tuned:", cycle.best), pch = 16, cex = 0.7)

make.rglwidget.graph <- function(coords,
                                 edges,
                                 theta = 35,
                                 phi = 22,
                                 zoom = 0.8,
                                 point.col = "#1F3B73",
                                 edge.col = "gray80",
                                 point.size = 5,
                                 width = 800,
                                 height = 520) {
  if (!requireNamespace("rgl", quietly = TRUE)) {
    return("Install the optional 'rgl' package to see the interactive 3D widget here.")
  }

  rgl::open3d(useNULL = TRUE)
  rgl::bg3d(color = "white")

  apply(edges, 1, function(e) {
    rgl::segments3d(
      coords[e, 1],
      coords[e, 2],
      coords[e, 3],
      color = edge.col
    )
  })

  rgl::points3d(
    coords[, 1],
    coords[, 2],
    coords[, 3],
    color = point.col,
    size = point.size
  )

  rgl::view3d(theta = theta, phi = phi, zoom = zoom)
  widget <- rgl::rglwidget(width = width, height = height)
  rgl::close3d()
  widget
}

plot.trace.frame <- function(coords, edges, main = "", pch = 16, cex = 0.55) {
  active <- stats::complete.cases(coords[, 1:2, drop = FALSE])
  xy <- coords[active, 1:2, drop = FALSE]

  xlim <- range(xy[, 1])
  ylim <- range(xy[, 2])
  xpad <- 0.08 * diff(xlim)
  ypad <- 0.08 * diff(ylim)
  if (!is.finite(xpad) || xpad == 0) xpad <- 0.2
  if (!is.finite(ypad) || ypad == 0) ypad <- 0.2

  plot(
    xy[, 1], xy[, 2],
    type = "n",
    asp = 1,
    axes = FALSE,
    xlab = "",
    ylab = "",
    xlim = xlim + c(-xpad, xpad),
    ylim = ylim + c(-ypad, ypad),
    main = main
  )

  active.edges <- edges[active[edges[, 1]] & active[edges[, 2]], , drop = FALSE]
  if (nrow(active.edges) > 0) {
    apply(active.edges, 1, function(e) {
      graphics::segments(
        coords[e[1], 1], coords[e[1], 2],
        coords[e[2], 1], coords[e[2], 2],
        col = "gray80"
      )
    })
  }

  points(xy[, 1], xy[, 2], pch = pch, cex = cex, col = "black")
}

edges <- edges.mesh(10, 10)

mesh.2d.default <- grip.layout(edges, n = 100, dim = 2, seed = 11)
mesh.2d.preset <- grip.layout(edges, n = 100, dim = 2, preset = "mesh", seed = 11)

op <- par(mfrow = c(1, 2), mar = c(1.2, 1.2, 3.2, 1.2))
on.exit(par(op), add = TRUE)

grip.plot(mesh.2d.default, edges, main = "Mesh default", pch = 16, cex = 0.45)
grip.plot(mesh.2d.preset, edges, main = "Mesh preset", pch = 16, cex = 0.45)

edges <- edges.sierpinski.carpet(3)
n <- max(edges)

carpet.2d.default <- grip.layout(edges, n = n, dim = 2, seed = 12)
carpet.2d.preset <- grip.layout(edges, n = n, dim = 2, preset = "carpet", seed = 12)

op <- par(mfrow = c(1, 2), mar = c(1.2, 1.2, 3.2, 1.2))
on.exit(par(op), add = TRUE)

grip.plot(carpet.2d.default, edges, main = "Carpet default", pch = 16, cex = 0.28)
grip.plot(carpet.2d.preset, edges, main = "Carpet preset", pch = 16, cex = 0.28)

edges <- edges.kary.tree(2, 6)
n <- max(edges)

tree.2d.default <- grip.layout(edges, n = n, dim = 2, seed = 13)
tree.2d.preset <- grip.layout(edges, n = n, dim = 2, preset = "tree", seed = 13)

op <- par(mfrow = c(1, 2), mar = c(1.2, 1.2, 3.2, 1.2))
on.exit(par(op), add = TRUE)

grip.plot(tree.2d.default, edges, main = "Tree default", pch = 16, cex = 0.45)
grip.plot(tree.2d.preset, edges, main = "Tree preset", pch = 16, cex = 0.45)

edges <- edges.cylinder(8, 12)
n <- max(edges)

cylinder.3d.default <- grip.layout(edges, n = n, dim = 3, seed = 21)
cylinder.3d.preset <- grip.layout(edges, n = n, dim = 3, preset = "torus", seed = 21)

make.rglwidget.graph(cylinder.3d.default, edges)

make.rglwidget.graph(cylinder.3d.preset, edges)

edges <- edges.sierpinski.tetrahedron(4)
n <- max(edges)

tetrahedron.3d.default <- grip.layout(edges, n = n, dim = 3, seed = 22)
tetrahedron.3d.tuned <- grip.layout(
  edges, n = n, dim = 3,
  placement = "barycenter",
  rounds = 128,
  final_rounds = 192,
  num_init = 12,
  num_nbrs = 16,
  r = 0.07,
  s = 9.0,
  repulsion_factor = 1.5,
  seed = 22
)

make.rglwidget.graph(tetrahedron.3d.default, edges)

make.rglwidget.graph(tetrahedron.3d.tuned, edges)

edges <- edges.mesh(6, 6)
mesh.2d.trace <- grip.layout.trace(
  edges,
  n = max(edges),
  dim = 2,
  preset = "mesh",
  seed = 31,
  trace = "level"
)

head(mesh.2d.trace$meta)

frame.idx <- unique(c(
  1L,
  max(2L, floor(length(mesh.2d.trace$frames) / 3)),
  max(3L, floor(2 * length(mesh.2d.trace$frames) / 3)),
  length(mesh.2d.trace$frames)
))

frame.labels <- c(
  "Coarsest start",
  "Early refinement",
  "Late refinement",
  "Final layout"
)

op <- par(mfrow = c(2, 2), mar = c(1.2, 1.2, 3.0, 1.2))
on.exit(par(op), add = TRUE)

for (i in seq_along(frame.idx)) {
  plot.trace.frame(
    mesh.2d.trace$frames[[frame.idx[[i]]]],
    # In early trace frames, only active vertices are shown.
    edges,
    main = frame.labels[[i]],
    cex = 0.45
  )
}
