## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dev = "jpeg",
  dev.args = list(quality = 90),
  fig.retina = 1,
  dpi = 96,
  fig.width = 11,
  fig.height = 4.2
)

## -----------------------------------------------------------------------------
library(grip)

## ----view-helpers, include=FALSE----------------------------------------------
source(system.file('scripts','vignette-views.R',package='grip'),local=TRUE)

## ----mesh-family--------------------------------------------------------------
mesh.graph <- mesh.surface.graph(
  5, 5,
  surface = "saddle",
  amplitude = 0.8
)

mesh.coords <- grip(metric = "edge_length",
  mesh.graph$edges,
  n = mesh.graph$n,
  edge.weights = mesh.graph$edge_weights,
  dim = 3,
  preset = "mesh",
  seed = 1
)

mesh.summary <- data.frame(
  family = mesh.graph$family,
  surface = mesh.graph$surface,
  n = mesh.graph$n,
  m = nrow(mesh.graph$edges),
  weight_scale = mesh.graph$weight_scale
)

knitr::kable(mesh.summary, digits = 3)

## ----mesh-interactive---------------------------------------------------------
documentation.view(list(`Weighted GRIP`=mesh.coords, `Reference geometry`=mesh.graph$coords_surface, `Parameter grid`=cbind(mesh.graph$coords_param,z=0)), edges=mesh.graph$edges,
  aligned=TRUE, description="Reference and fitted geometries are centered without rotation or scale fitting.")

## ----torus-family-------------------------------------------------------------
torus.graph <- torus.surface.graph(
  5, 6,
  surface = "pinched",
  minor.radius = 0.33,
  amplitude = 0.18,
  twist = 0.18
)

torus.coords <- grip(metric = "edge_length",
  torus.graph$edges,
  n = torus.graph$n,
  edge.weights = torus.graph$edge_weights,
  dim = 3,
  preset = "torus",
  seed = 2
)

knitr::kable(
  data.frame(
    family = torus.graph$family,
    surface = torus.graph$surface,
    n = torus.graph$n,
    m = nrow(torus.graph$edges)
  )
)

## ----torus-interactive--------------------------------------------------------
documentation.view(list(`Weighted GRIP`=torus.coords, `Reference geometry`=torus.graph$coords_surface), edges=torus.graph$edges,
  aligned=TRUE, description="Reference and fitted geometries are centered without rotation or scale fitting.")

## ----vicsek-family------------------------------------------------------------
vicsek.graph <- vicsek.surface.graph(
  level = 2,
  surface = "ripple",
  amplitude = 0.55
)

vicsek.coords <- grip(metric = "edge_length",
  vicsek.graph$edges,
  n = vicsek.graph$n,
  edge.weights = vicsek.graph$edge_weights,
  dim = 3,
  preset = "carpet",
  seed = 3
)

## ----vicsek-interactive-------------------------------------------------------
documentation.view(list(`Weighted GRIP`=vicsek.coords, `Reference geometry`=vicsek.graph$coords_surface), edges=vicsek.graph$edges,
  aligned=TRUE, description="Reference and fitted geometries are centered without rotation or scale fitting.")

## ----annulus-family-----------------------------------------------------------
annulus.graph <- irregular.annulus.surface.graph(
  rings = 5,
  outer.count = 24,
  inner.radius = 0.42,
  surface = "folded",
  amplitude = 0.2
)

annulus.coords <- grip(metric = "edge_length",
  annulus.graph$edges,
  n = annulus.graph$n,
  edge.weights = annulus.graph$edge_weights,
  dim = 3,
  preset = "irregular",
  seed = 4
)

## ----annulus-interactive------------------------------------------------------
documentation.view(list(`Weighted GRIP`=annulus.coords, `Reference geometry`=annulus.graph$coords_surface), edges=annulus.graph$edges,
  aligned=TRUE, description="Reference and fitted geometries are centered without rotation or scale fitting.")

## ----cube-family--------------------------------------------------------------
cube.graph <- cube.channel.network.surface.graph(
  side = 5,
  level = 1,
  surface = "twisted",
  twist = 0.2
)

cube.coords <- grip(metric = "edge_length",
  cube.graph$edges,
  n = cube.graph$n,
  edge.weights = cube.graph$edge_weights,
  dim = 3,
  preset = "irregular",
  seed = 5
)

knitr::kable(
  data.frame(
    family = cube.graph$family,
    n = cube.graph$n,
    m = nrow(cube.graph$edges)
  )
)

## ----cube-interactive---------------------------------------------------------
documentation.view(list(`Weighted GRIP`=cube.coords, `Reference geometry`=cube.graph$coords_surface), edges=cube.graph$edges,
  aligned=TRUE, description="Reference and fitted geometries are centered without rotation or scale fitting.")

## ----tree-family--------------------------------------------------------------
tree.graph <- kary.tree.weighted.graph(
  k = 2,
  depth = 4,
  depth.rule = "geometric",
  depth.decay = 0.82,
  branch.rule = "linear",
  branch.spread = 0.25
)

tree.coords <- grip(metric = "edge_length",
  tree.graph$edges,
  n = tree.graph$n,
  edge.weights = tree.graph$edge_weights,
  dim = 2,
  preset = "tree",
  seed = 6
)

knitr::kable(
  head(tree.graph$edge_table[, c(
    "parent",
    "child",
    "child_depth",
    "branch_index",
    "edge_weight"
  )]),
  digits = 3
)

## ----tree-interactive---------------------------------------------------------
documentation.view(list(`Weighted GRIP (2D fit)`=cbind(tree.coords,z=0)), edges=tree.graph$edges,
  aligned=TRUE, description="Reference and fitted geometries are centered without rotation or scale fitting.")

