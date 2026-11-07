## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(fig.retina=1, dpi=72, collapse = TRUE, comment = "#>",
                     fig.width = 6, fig.height = 4.5)
# Website-only articles are not installed with the package.
article_url <- function(name) {
  prefix <- if (identical(Sys.getenv("IN_PKGDOWN"), "true")) "" else
    if (isTRUE(getOption("grip.local.preview"))) "../articles/" else
      "https://pgajer.github.io/grip/articles/"
  paste0(prefix, name, ".html")
}
# Installed help and pkgdown have different reference-page locations.
help_url <- function(topic) {
  prefix <- if (identical(Sys.getenv("IN_PKGDOWN"), "true")) {
    "../reference/"
  } else "../html/"
  paste0(prefix, topic, ".html")
}

## ----first-layout, fig.alt="A 20-vertex rectangular mesh drawn with the mesh preset, showing its rows and columns."----
library(grip)
edges <- edges.mesh(4, 5)
coords <- grip(edges, n = 20, dim = 2, preset = "mesh", seed = 11)
stopifnot(identical(dim(coords), c(20L, 2L)), all(is.finite(coords)))
source(system.file('scripts','vignette-views.R',package='grip'),local=TRUE)
documentation.view(list(`GRIP mesh (2D fit)`=cbind(coords,z=0)),edges,
  aligned=TRUE,description='A planar mesh fit shown at z = 0.')
quality <- score.layout(coords, edges = edges, n = 20,
                        sample.size.stress = 200, stress.seed = 11,
                        edge.crossings = "never")
knitr::kable(quality[, c("sampled.stress", "edge.length.cv")], digits = 3)

## ----isolated-vertex----------------------------------------------------------
path.with.isolate <- grip(edges.path(5), n = 6, dim = 2, seed = 11)
stopifnot(nrow(path.with.isolate) == 6L, all(is.finite(path.with.isolate)))

## ----common-extraction--------------------------------------------------------
small.edges <- edges.path(6)
fits <- list(
  ordinary = grip(small.edges, n = 6, dim = 2, seed = 1),
  trace = trace.grip(small.edges, n = 6, dim = 2, seed = 1),
  classical = classical.mds(edges = small.edges, n = 6)
)
matrices <- lapply(fits, layout.coords)
vapply(matrices, function(z) score.layout(z, edges = small.edges, n = 6,
  sample.size.stress = 15, stress.seed = 1,
  edge.crossings = "never")$sampled.stress, numeric(1))
# Plot any extracted matrix with the same graph:
# plot.layout(matrices[[1]], edges = small.edges)

## ----preparation-estimate-----------------------------------------------------
estimate.preparation(10000, n.edges = 20000)

## ----navigation, eval=FALSE---------------------------------------------------
# help(package = "grip")
# help("grip-package", package = "grip")
# help("grip", package = "grip")
# vignette(package = "grip")
# vignette("synthetic-graph-families", package = "grip")

