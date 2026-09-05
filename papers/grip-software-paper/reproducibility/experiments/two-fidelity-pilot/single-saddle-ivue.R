# Run this script section by section in R. All coordinates stay in 3D.
# Requires grip >= 0.2.0.9000, dgraphs >= 0.2.0, and the current ivue API.
stopifnot(packageVersion("grip") >= "0.2.0.9000",
          packageVersion("dgraphs") >= "0.2.0")

# Settings ---------------------------------------------------------------
n <- 1000
seed <- 1
C <- 0.8
half.width <- 1                         # Figure 7 uses [-1, 1]^2, not [-2, 2]^2.
ks <- 3:20
max.iter <- 200                         # Per density-mixing stage of edge-KK.

# 1. Sample one saddle ---------------------------------------------------
# Rejection sampling makes the sample uniform in surface area. Accepting
# every proposal instead would give uniform (x, y), not uniform surface area.
set.seed(seed)
xy <- matrix(numeric(), ncol = 2)
while (nrow(xy) < n) {
  proposal <- matrix(runif(4000, -half.width, half.width), ncol = 2)
  area <- sqrt(1 + 4 * C^2 * rowSums(proposal^2))
  accept <- runif(nrow(proposal)) < area / sqrt(1 + 8 * C^2 * half.width^2)
  xy <- rbind(xy, proposal[accept, , drop = FALSE])
}
xy <- xy[seq_len(n), , drop = FALSE]
X <- cbind(x = xy[, 1], y = xy[, 2], z = C * (xy[, 1]^2 - xy[, 2]^2))

# Color always identifies the ORIGINAL x coordinate, including after fitting.
colors <- ivue::color.scale.cont(X[, "x"], limits = c(-half.width, half.width))
ivue::plot3D.cont(X, values = X[, "x"], scale = colors, point.size = 4)
ivue::plot3D.plain(X, point.type = "sphere", sphere.radius = 0.015)

# 2. Symmetric kNN graphs -------------------------------------------------
# Euclidean edge lengths; component-MST repair connects a graph if necessary.
graphs <- setNames(lapply(ks, function(k) {
  dgraphs::create.sknn.graph(
    X, k = k, neighbor.method = "exact", edge.weight = "distance",
    connect.components = TRUE, connect.method = "component.mst",
    prune.method = "none", prune.edges = FALSE
  )
}), as.character(ks))

k <- 10                                # Change k and rerun this view only.
g <- graphs[[as.character(k)]]
ivue::plot3D.graph(list(adj.list = g$adj_list, weight.list = g$weight_list),
                  X = X, vertices = seq_len(n),
                  values = X[, "x"], scale = colors,
                  edge.col = "gray80", edge.width = 1)

# 3. Metric MDS for every k -----------------------------------------------
# diagnostics = FALSE avoids storing all fixed paths: this script needs
# coordinates, not the reference-distance calibration or fidelity audit.
mds <- lapply(graphs, function(g) {
  grip::classical.mds(edges = g$edge_matrix, edge_weights = g$edge_weight,
                   n = n, dim = 3, diagnostics = FALSE)$coords
})

k <- 10
g <- graphs[[as.character(k)]]
ivue::plot3D.graph(list(adj.list = g$adj_list, weight.list = g$weight_list),
                  X = mds[[as.character(k)]],
                  vertices = seq_len(n), values = X[, "x"], scale = colors,
                  edge.col = "gray80", edge.width = 1)
# Point-only alternative:
ivue::plot3D.cont(mds[[as.character(k)]], values = X[, "x"], scale = colors)

# 4. Metric MDS followed by edge-KK for every k ----------------------------
# Same primary refinement settings as Figure 7; no extra audit continuation.
mds.edge.kk <- setNames(lapply(ks, function(k) {
  message("edge-KK: k = ", k)
  g <- graphs[[as.character(k)]]
  grip::edge.kk(
    coords = mds[[as.character(k)]], edges = g$edge_matrix,
    edge_weights = g$edge_weight, n = n, dim = 3, max_iter = max.iter,
    stiffness_method = "density", density_mix_schedule = c(0, .25, .5, .75, 1),
    scale_mode = "profiled", edge_length_epsilon = 0,
    diagnostics = FALSE, return_trace = FALSE, seed = seed
  )$coords
}), as.character(ks))

k <- 10
g <- graphs[[as.character(k)]]
ivue::plot3D.graph(list(adj.list = g$adj_list, weight.list = g$weight_list),
                  X = mds.edge.kk[[as.character(k)]],
                  vertices = seq_len(n), values = X[, "x"], scale = colors,
                  edge.col = "gray80", edge.width = 1)
ivue::plot3D.cont(mds.edge.kk[[as.character(k)]], values = X[, "x"], scale = colors)

# 5. Graph fidelity and reference-saddle recovery --------------------------
# Run from the grip repository root after fitting the layouts above.
# Requires the development grip version with score.coordinates/score.surface.
source("papers/grip-software-paper/reproducibility/experiments/two-fidelity-pilot/score-saddle-reference.R")
reference.results <- score.saddle.reference(
  X, graphs, mds, mds.edge.kk, C,
  sample_sizes = c(2000L, 8000L)
)
# Full results retain both alignments, both sampling resolutions, and all k.
reference.scores <- reference.results$scores
subset(reference.scores, k == 10 & sample_size == 8000,
       select = c(method, alignment, path_rel, edge_rel, stress1,
                  coordinate_rmse, coordinate_relative_rmse, alignment_scale,
                  surface_rms, surface_mc_se))
