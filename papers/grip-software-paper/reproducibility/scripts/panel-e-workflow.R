# Run from papers/grip-software-paper (or the extracted submission directory).
# source("reproducibility/scripts/panel-e-workflow.R") recomputes the fit,
# scores it, and creates the interactive panel-E view; nothing is overwritten.
# Routine manuscript builds instead use the supplied saved fit and scores.

## ---- panel-e-input ----
example <- readRDS("reproducibility/precomputed/two-fidelity-saddle.rds")
source("reproducibility/scripts/panel-e-display.R")

## ---- panel-e-fit ----
prepared <- grip::prepare.geodesic.kk(
  example$edges, n = nrow(example$coords),
  edge_weights = example$weights, tie_mode = "single")
mds <- grip::metric.mds(prepared = prepared, dim = 3, diagnostics = FALSE)
fit <- grip::edge.kk(
  coords = mds$coords, prepared = prepared, dim = 3,
  max_iter = 200L, stiffness_method = "density",
  density_mix_schedule = c(0, .25, .5, .75, 1),
  scale_mode = "profiled", edge_length_epsilon = 0,
  diagnostics = FALSE, return_trace = TRUE, seed = 5211005L)
before <- grip::score.gmds(mds$coords, prepared = prepared,
  scale_mode = "profiled", edge_length_epsilon = 0)
after <- grip::score.gmds(fit$coords, prepared = prepared,
  scale_mode = "profiled", edge_length_epsilon = 0)

## ---- panel-e-view ----
display <- panel.e.display(example, fit$coords)
mesh <- ivue::layer3D.mesh(display$triangles,
  col = "gray75", alpha = .23, edge.col = "gray45",
  edge.alpha = .23, edge.width = .65)
panel.e <- ivue::plot3D.plain(display$coords, col = display$colors,
  point.type = "sphere", sphere.radius = .009,
  axes = FALSE, aspect = "equal",
  layers = list(mesh, display$bounds),
  camera = ivue::camera.zup(elevation = 20, turn = -135, zoom = .7),
  width = 720L, height = 640L)
panel.e
