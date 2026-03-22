## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")


## -----------------------------------------------------------------------------
library(grip)


## -----------------------------------------------------------------------------
edges <- edges.path(12)
coords <- grip.layout(edges, n = 12, dim = 2,
                      placement = "barycenter",
                      rounds = 25, final_rounds = 25,
                      num_init = 5, num_nbrs = 6,
                      seed = 1)
grip.plot(coords, edges, main = "Path (2D)", pch = 16, cex = 0.6)


## -----------------------------------------------------------------------------
edges <- edges.cycle(16)
coords <- grip.layout(edges, n = 16, dim = 2,
                      placement = "barycenter",
                      rounds = 25, final_rounds = 25,
                      num_init = 5, num_nbrs = 6,
                      seed = 2)
grip.plot(coords, edges, main = "Cycle (2D)", pch = 16, cex = 0.6)


## -----------------------------------------------------------------------------
edges <- edges.mesh(5, 5)
coords <- grip.layout(edges, n = 25, dim = 2,
                      placement = "barycenter",
                      rounds = 25, final_rounds = 25,
                      num_init = 6, num_nbrs = 8,
                      seed = 3)
grip.plot(coords, edges, main = "Mesh (2D)", pch = 16, cex = 0.6)


## -----------------------------------------------------------------------------
edges <- edges.sierpinski.triangle(2)
n <- max(edges)
coords <- grip.layout(edges, n = n, dim = 2,
                      placement = "circle",
                      rounds = 25, final_rounds = 25,
                      num_init = 5, num_nbrs = 7,
                      seed = 4)
grip.plot(coords, edges, main = "Sierpinski (2D, circle placement)", pch = 16, cex = 0.6)


## -----------------------------------------------------------------------------
edges <- edges.mesh(4, 4)
coords <- grip.layout(edges, n = 16, dim = 3,
                      placement = "barycenter",
                      rounds = 25, final_rounds = 25,
                      num_init = 5, num_nbrs = 7,
                      seed = 5)
head(coords)


## -----------------------------------------------------------------------------
adj_list <- list(c(2), c(1, 3), c(2, 4), c(3))
weight_list <- list(c(1.0), c(1.0, 2.0), c(2.0, 1.5), c(1.5))
coords <- grip.layout(adj_list = adj_list,
                      weight_list = weight_list,
                      n = 4,
                      dim = 2,
                      placement = "barycenter",
                      rounds = 25, final_rounds = 25,
                      num_init = 3, num_nbrs = 3,
                      seed = 12)
grip.plot(coords, main = "Adjacency list input (2D)")


## -----------------------------------------------------------------------------
edges <- edges.cylinder(4, 6)
coords <- grip.layout(edges, n = 24, dim = 3,
                      placement = "barycenter",
                      rounds = 25, final_rounds = 25,
                      num_init = 6, num_nbrs = 8,
                      seed = 6)
head(coords)


## -----------------------------------------------------------------------------
edges <- edges.torus(4, 4)
coords <- grip.layout(edges, n = 16, dim = 3,
                      placement = "barycenter",
                      rounds = 25, final_rounds = 25,
                      num_init = 5, num_nbrs = 7,
                      seed = 7)
head(coords)


## ----eval=FALSE---------------------------------------------------------------
# # Optional 3D plotting (requires rgl)
# if (requireNamespace("rgl", quietly = TRUE)) {
#   rgl::plot3d(coords[, 1], coords[, 2], coords[, 3],
#               size = 4, col = "steelblue")
# }
