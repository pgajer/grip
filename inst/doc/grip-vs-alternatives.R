## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dev = "jpeg",
  dev.args = list(quality = 90),
  fig.retina = 1,
  dpi = 96,
  fig.width = 10,
  fig.height = 4.5
)

## ----score-helper, eval=FALSE-------------------------------------------------
# score_layout <- function(coords, edges, n,
#                          clusters = NULL,
#                          edge.crossings = "auto",
#                          sample.size.stress = 2000L,
#                          sample.size.nonedge = 5000L) {
#   score.layout(
#     coords = coords,
#     edges = edges,
#     n = n,
#     clusters = clusters,
#     sample.size.stress = sample.size.stress,
#     sample.size.nonedge = sample.size.nonedge,
#     stress.seed = 42L,
#     nonedge.seed = 42L,
#     edge.crossings = edge.crossings
#   )
# }

## ----interactive-helpers, include=FALSE---------------------------------------
source(system.file('scripts','vignette-views.R',package='grip'),local=TRUE)
comparison.views <- function(layouts, edges) {
  htmltools::tagList(lapply(names(layouts), function(name) {
    z <- layouts[[name]]
    if (ncol(z)==2L) z <- cbind(z,0)
    documentation.view(setNames(list(z),name),edges,aligned=TRUE,
      description=paste(name, ': saved coordinates; 2D fits remain planar.'))
  }))
}

## ----load-results-------------------------------------------------------------
res.path <- system.file(
  "extdata", "vs_alternatives", "benchmark_results.rds",
  package = "grip"
)
if (!nzchar(res.path)) {
  stop(
    "Precomputed benchmark results not found. ",
    "Run inst/scripts/precompute-vs-alternatives.R first."
  )
}
res <- readRDS(res.path)

## ----karate-table-------------------------------------------------------------
karate.display <- res$karate$scores[, c(
  "method",
  "sampled.stress",
  "edge.length.cv",
  "edge.crossings",
  "sampled.nonedge.sep.ratio",
  "cluster.separation"
)]

knitr::kable(
  karate.display,
  digits = 4,
  row.names = FALSE,
  caption = "Karate club: precomputed comparison; one seeded layout per method."
)

## ----karate-plot, fig.width=12, fig.height=7.5--------------------------------
k <- res$karate
comparison.views(k$layouts, k$edges)

## ----mesh-table---------------------------------------------------------------
mesh.display <- res$mesh$scores[, c(
  "method",
  "sampled.stress",
  "edge.length.cv",
  "edge.crossings",
  "sampled.nonedge.sep.ratio"
)]

knitr::kable(
  mesh.display,
  digits = 4,
  row.names = FALSE,
  caption = "12x12 mesh: precomputed comparison; one seeded layout per method."
)

## ----mesh-plot, fig.width=12, fig.height=7.5----------------------------------
m <- res$mesh
comparison.views(m$layouts, m$edges)

## ----hmp-saved-results--------------------------------------------------------
knitr::kable(res$hmp$timing, digits = 3, row.names = FALSE,
  caption = "HMP graph: elapsed seconds for the precomputed layout runs.")
knitr::kable(res$hmp$scores[, c("method", "sampled.stress", "edge.length.cv")],
  digits = 4, row.names = FALSE,
  caption = "HMP graph: hop-distance quality of the first timed layout per method.")

## ----reproduce, eval=FALSE----------------------------------------------------
# Sys.setenv(GRIP_VS_ALTERNATIVES_OUTPUT = "benchmark_results.rds")
# source(system.file(
#   "scripts", "precompute-vs-alternatives.R",
#   package = "grip"
# ))

