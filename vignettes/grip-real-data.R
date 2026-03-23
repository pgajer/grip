## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 5.5
)


## -----------------------------------------------------------------------------
library(grip)


## -----------------------------------------------------------------------------
read.extdata.csv <- function(file.name) {
  candidates <- c(
    system.file("extdata", file.name, package = "grip"),
    file.path("inst", "extdata", file.name),
    file.path("..", "inst", "extdata", file.name)
  )
  path <- candidates[file.exists(candidates)][1L]
  if (!length(path) || is.na(path) || !nzchar(path)) {
    stop("could not locate ", file.name)
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}

read.edge.csv <- function(file.name) {
  as.matrix(read.extdata.csv(file.name))
}

read.karate.club <- function() {
  labels <- read.extdata.csv("karate-club-membership.csv")
  labels <- labels[order(labels$vertex), , drop = FALSE]
  labels$club
}

compact.summary <- function(x) {
  keep <- intersect(c(
    "candidate",
    "preset",
    "rounds",
    "final.rounds",
    "num.nbrs",
    "repulsion.factor",
    "score.composite",
    "sampled.stress.mean",
    "edge.length.cv.mean",
    "sampled.nonedge.sep.ratio.mean",
    "stability.procrustes.mean",
    "cluster.separation.mean"
  ), names(x))
  x[, keep, drop = FALSE]
}

candidate.from.summary <- function(row) {
  row <- row[1, , drop = FALSE]
  out <- list(
    placement = row$placement[[1L]],
    rounds = row$rounds[[1L]],
    final_rounds = row$final.rounds[[1L]],
    num_init = row$num.init[[1L]],
    num_nbrs = row$num.nbrs[[1L]],
    r = row$r[[1L]],
    s = row$s[[1L]],
    repulsion_factor = row$repulsion.factor[[1L]],
    tinit_factor = row$tinit.factor[[1L]]
  )
  if (!is.na(row$preset[[1L]]) && nzchar(row$preset[[1L]])) {
    out$preset <- row$preset[[1L]]
  }
  out
}

project.coords <- function(coords, azimuth = 35, elevation = 22) {
  if (ncol(coords) == 2L) {
    return(coords)
  }
  az <- azimuth * pi / 180
  el <- elevation * pi / 180
  rz <- matrix(c(
    cos(az), -sin(az), 0,
    sin(az),  cos(az), 0,
    0,        0,       1
  ), nrow = 3, byrow = TRUE)
  rx <- matrix(c(
    1, 0, 0,
    0, cos(el), -sin(el),
    0, sin(el),  cos(el)
  ), nrow = 3, byrow = TRUE)
  rotated <- coords %*% rz %*% rx
  rotated[, 1:2, drop = FALSE]
}

plot.projected.graph <- function(coords,
                                 edges,
                                 main = "",
                                 azimuth = 35,
                                 elevation = 22,
                                 vertex.col = "black",
                                 edge.col = "gray80",
                                 pch = 16,
                                 cex = 0.8) {
  xy <- project.coords(coords, azimuth = azimuth, elevation = elevation)
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
  if (!is.null(edges) && nrow(edges) > 0) {
    apply(edges, 1, function(e) {
      graphics::segments(
        xy[e[1], 1], xy[e[1], 2],
        xy[e[2], 1], xy[e[2], 2],
        col = edge.col
      )
    })
  }
  points(xy[, 1], xy[, 2], pch = pch, cex = cex, col = vertex.col)
}

plot.layout.triptych <- function(layouts, edges, titles, vertex.cols = NULL) {
  if (is.null(vertex.cols)) {
    vertex.cols <- rep(list("black"), length(layouts))
  }
  op <- par(mfrow = c(1, length(layouts)), mar = c(1.2, 1.2, 3.2, 1.2))
  on.exit(par(op), add = TRUE)
  for (i in seq_along(layouts)) {
    plot.projected.graph(
      layouts[[i]],
      edges = edges,
      main = titles[[i]],
      vertex.col = vertex.cols[[i]]
    )
  }
}

plot.layout.views <- function(coords, edges, title, vertex.col = "black") {
  op <- par(mfrow = c(1, 2), mar = c(1.2, 1.2, 3.2, 1.2))
  on.exit(par(op), add = TRUE)
  plot.projected.graph(coords, edges, main = paste(title, "\nview 1"),
                       azimuth = 35, elevation = 22, vertex.col = vertex.col)
  plot.projected.graph(coords, edges, main = paste(title, "\nview 2"),
                       azimuth = -45, elevation = 28, vertex.col = vertex.col)
}

plot.search.slices <- function(summary.df, title.prefix = "") {
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
  on.exit(par(op), add = TRUE)

  palette.cols <- grDevices::hcl.colors(length(unique(summary.df$num.nbrs)), "Dark 3")
  color.idx <- match(summary.df$num.nbrs, sort(unique(summary.df$num.nbrs)))

  plot(
    summary.df$repulsion.factor,
    summary.df$sampled.stress.mean,
    pch = 19,
    col = palette.cols[color.idx],
    xlab = "repulsion.factor",
    ylab = "sampled.stress.mean",
    main = paste(title.prefix, "stress")
  )
  legend(
    "topright",
    legend = paste("num.nbrs =", sort(unique(summary.df$num.nbrs))),
    col = palette.cols,
    pch = 19,
    bty = "n",
    cex = 0.8
  )

  plot(
    summary.df$final.rounds,
    summary.df$sampled.nonedge.sep.ratio.mean,
    pch = 19,
    col = palette.cols[color.idx],
    xlab = "final.rounds",
    ylab = "sampled.nonedge.sep.ratio.mean",
    main = paste(title.prefix, "separation")
  )
}


## -----------------------------------------------------------------------------
karate.edges <- read.edge.csv("karate-club-edges.csv")
karate.n <- max(karate.edges)
karate.club <- read.karate.club()

kite.edges <- read.edge.csv("krackhardt-kite-edges.csv")
kite.n <- max(kite.edges)

benchmark.info <- data.frame(
  graph = c("Zachary karate club", "Krackhardt kite"),
  n.vertices = c(karate.n, kite.n),
  n.edges = c(nrow(karate.edges), nrow(kite.edges)),
  stringsAsFactors = FALSE
)

knitr::kable(benchmark.info)


## -----------------------------------------------------------------------------
karate.presets <- grip.compare.layouts(
  edges = karate.edges,
  n = karate.n,
  dim = 3,
  candidates = c("default", "tree", "mesh"),
  seeds = 1:3,
  sample.size.stress = 1000L,
  sample.size.nonedge = 2000L,
  edge.crossings = "never",
  return.layouts = TRUE
)

knitr::kable(compact.summary(karate.presets$summary), digits = 3)


## ----fig.width=10, fig.height=3.8---------------------------------------------
plot.layout.triptych(
  layouts = list(
    karate.presets$layouts$default[["1"]],
    karate.presets$layouts$tree[["1"]],
    karate.presets$layouts$mesh[["1"]]
  ),
  edges = karate.edges,
  titles = c("default", "tree", "mesh")
)


## -----------------------------------------------------------------------------
karate.base <- karate.presets$summary[1, , drop = FALSE]

karate.search <- list(
  candidate.prefix = "karate.local",
  placement = karate.base$placement[[1L]],
  rounds = karate.base$rounds[[1L]],
  final_rounds = sort(unique(as.integer(c(
    max(32L, karate.base$final.rounds[[1L]] - 32L),
    karate.base$final.rounds[[1L]],
    karate.base$final.rounds[[1L]] + 32L
  )))),
  num_init = karate.base$num.init[[1L]],
  num_nbrs = sort(unique(as.integer(c(
    max(4L, karate.base$num.nbrs[[1L]] - 4L),
    karate.base$num.nbrs[[1L]],
    karate.base$num.nbrs[[1L]] + 4L
  )))),
  r = karate.base$r[[1L]],
  s = karate.base$s[[1L]],
  repulsion_factor = sort(unique(c(
    max(0, karate.base$repulsion.factor[[1L]] - 0.5),
    karate.base$repulsion.factor[[1L]],
    karate.base$repulsion.factor[[1L]] + 0.5
  )))
)
if (!is.na(karate.base$preset[[1L]]) && nzchar(karate.base$preset[[1L]])) {
  karate.search$preset <- karate.base$preset[[1L]]
}

karate.local <- grip.compare.layouts(
  edges = karate.edges,
  n = karate.n,
  dim = 3,
  search = karate.search,
  seeds = 1:2,
  sample.size.stress = 1000L,
  sample.size.nonedge = 2000L,
  edge.crossings = "never",
  return.layouts = TRUE
)

knitr::kable(head(compact.summary(karate.local$summary), 8L), digits = 3)


## -----------------------------------------------------------------------------
plot.search.slices(karate.local$summary, title.prefix = "Karate local search")


## -----------------------------------------------------------------------------
karate.local.best <- karate.local$summary[1, , drop = FALSE]
karate.local.best.seed <- karate.local$runs$seed[
  karate.local$runs$candidate == karate.local.best$candidate[[1L]] &
    karate.local$runs$status == "ok"
][[1L]]
karate.local.best.coords <- karate.local$layouts[[karate.local.best$candidate[[1L]]]][[as.character(karate.local.best.seed)]]


## ----fig.width=8, fig.height=4.2----------------------------------------------
plot.layout.views(
  karate.local.best.coords,
  karate.edges,
  title = paste("Karate local-search winner:", karate.local.best$candidate[[1L]])
)


## -----------------------------------------------------------------------------
karate.metadata <- grip.compare.layouts(
  edges = karate.edges,
  n = karate.n,
  dim = 3,
  candidates = list(
    default = NULL,
    tuned = candidate.from.summary(karate.local.best)
  ),
  clusters = karate.club,
  seeds = 1:3,
  sample.size.stress = 1000L,
  sample.size.nonedge = 2000L,
  edge.crossings = "never",
  return.layouts = TRUE
)

knitr::kable(
  karate.metadata$summary[, c(
    "candidate",
    "score.composite",
    "sampled.stress.mean",
    "sampled.nonedge.sep.ratio.mean",
    "cluster.separation.mean"
  )],
  digits = 3
)


## -----------------------------------------------------------------------------
karate.metadata.best <- karate.metadata$summary[1, , drop = FALSE]
karate.metadata.best.seed <- karate.metadata$runs$seed[
  karate.metadata$runs$candidate == karate.metadata.best$candidate[[1L]] &
    karate.metadata$runs$status == "ok"
][[1L]]
karate.metadata.best.coords <- karate.metadata$layouts[[karate.metadata.best$candidate[[1L]]]][[as.character(karate.metadata.best.seed)]]
karate.club.cols <- ifelse(karate.club == "Mr. Hi", "#1b9e77", "#d95f02")


## ----fig.width=8, fig.height=4.2----------------------------------------------
plot.layout.views(
  karate.metadata.best.coords,
  karate.edges,
  title = paste("Karate best layout with metadata:", karate.metadata.best$candidate[[1L]]),
  vertex.col = karate.club.cols
)
legend(
  "bottomleft",
  legend = c("Mr. Hi", "Officer"),
  col = c("#1b9e77", "#d95f02"),
  pch = 16,
  bty = "n"
)


## -----------------------------------------------------------------------------
kite.compare <- grip.compare.layouts(
  edges = kite.edges,
  n = kite.n,
  dim = 3,
  candidates = c("default", "tree"),
  search = list(
    candidate.prefix = "kite.tree.local",
    preset = "tree",
    final_rounds = c(128L, 160L, 192L),
    num_nbrs = c(6L, 8L),
    repulsion_factor = c(0.0, 0.5)
  ),
  seeds = 1:2,
  sample.size.stress = 400L,
  sample.size.nonedge = 800L,
  edge.crossings = "never",
  return.layouts = TRUE
)

knitr::kable(head(compact.summary(kite.compare$summary), 8L), digits = 3)


## -----------------------------------------------------------------------------
kite.best <- kite.compare$summary[1, , drop = FALSE]
kite.best.seed <- kite.compare$runs$seed[
  kite.compare$runs$candidate == kite.best$candidate[[1L]] &
    kite.compare$runs$status == "ok"
][[1L]]
kite.best.coords <- kite.compare$layouts[[kite.best$candidate[[1L]]]][[as.character(kite.best.seed)]]


## ----fig.width=8, fig.height=4.2----------------------------------------------
plot.layout.views(
  kite.best.coords,
  kite.edges,
  title = paste("Krackhardt kite winner:", kite.best$candidate[[1L]])
)


## -----------------------------------------------------------------------------
knitr::kable(
  grip.score.layout(
    coords = karate.metadata.best.coords,
    edges = karate.edges,
    n = karate.n,
    clusters = karate.club,
    edge.crossings = "never"
  ),
  digits = 3
)

