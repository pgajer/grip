## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 5)


## -----------------------------------------------------------------------------
library(grip)


## -----------------------------------------------------------------------------
read.edge.csv <- function(file.name) {
  candidates <- c(
    system.file("extdata", file.name, package = "grip"),
    file.path("inst", "extdata", file.name),
    file.path("..", "inst", "extdata", file.name)
  )
  path <- candidates[file.exists(candidates)][1L]
  if (!length(path) || is.na(path) || !nzchar(path)) {
    stop("could not locate ", file.name)
  }
  as.matrix(utils::read.csv(path))
}

compact.summary <- function(x) {
  keep <- c(
    "candidate",
    "score.composite",
    "sampled.stress.mean",
    "edge.length.cv.mean",
    "sampled.nonedge.sep.ratio.mean",
    "stability.procrustes.mean",
    "edge.crossings.mean"
  )
  x[, keep]
}


## -----------------------------------------------------------------------------
kite.edges <- read.edge.csv("krackhardt-kite-edges.csv")
kite.n <- max(kite.edges)

kite.coords <- grip.layout(
  kite.edges,
  n = kite.n,
  dim = 2,
  preset = "tree",
  seed = 1
)

kite.score <- grip.score.layout(
  coords = kite.coords,
  edges = kite.edges,
  n = kite.n,
  edge.crossings = "always"
)

knitr::kable(kite.score, digits = 3)


## -----------------------------------------------------------------------------
grip.plot(kite.coords, kite.edges, main = "Krackhardt kite with preset = 'tree'",
          pch = 16, cex = 1)


## -----------------------------------------------------------------------------
karate.edges <- read.edge.csv("karate-club-edges.csv")
karate.n <- max(karate.edges)

karate.compare <- grip.compare.layouts(
  edges = karate.edges,
  n = karate.n,
  dim = 2,
  candidates = c("default", "tree", "mesh"),
  seeds = 1:3,
  sample.size.stress = 1000L,
  sample.size.nonedge = 2000L,
  edge.crossings = "always",
  return.layouts = TRUE
)

knitr::kable(compact.summary(karate.compare$summary), digits = 3)


## -----------------------------------------------------------------------------
best.karate <- karate.compare$summary$candidate[[1L]]
best.seed <- karate.compare$runs$seed[
  karate.compare$runs$candidate == best.karate &
    karate.compare$runs$status == "ok"
][[1L]]

best.karate.coords <- karate.compare$layouts[[best.karate]][[as.character(best.seed)]]

grip.plot(best.karate.coords, karate.edges,
          main = paste("Karate club best candidate:", best.karate),
          pch = 16, cex = 0.8)


## -----------------------------------------------------------------------------
kite.compare <- grip.compare.layouts(
  edges = kite.edges,
  n = kite.n,
  dim = 2,
  candidates = list(
    default = NULL,
    tree = "tree",
    tree.long = list(preset = "tree", final_rounds = 224)
  ),
  seeds = 1:3,
  sample.size.stress = 800L,
  sample.size.nonedge = 1500L,
  edge.crossings = "always"
)

knitr::kable(compact.summary(kite.compare$summary), digits = 3)

