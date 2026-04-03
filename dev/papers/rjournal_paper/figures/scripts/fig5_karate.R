#!/usr/bin/env Rscript
## Figure 5: Zachary karate club, tuned tree-based 3D layout colored by faction.
## Two orthographic views of a static, publication-friendly orthographic render.

library(grip)

figure_output_path <- function(file.name) {
  args <- commandArgs(trailingOnly = FALSE)
  file.arg <- grep("^--file=", args, value = TRUE)
  script.dir <- if (length(file.arg)) {
    dirname(normalizePath(sub("^--file=", "", file.arg[1L])))
  } else {
    getwd()
  }
  out.dir <- file.path(normalizePath(file.path(script.dir, "..")), "figures")
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
  file.path(out.dir, file.name)
}

## --- helpers ---
read.extdata.csv <- function(file.name) {
  candidates <- c(
    system.file("extdata", file.name, package = "grip"),
    file.path("inst", "extdata", file.name),
    file.path("..", "inst", "extdata", file.name)
  )
  path <- candidates[file.exists(candidates)][1L]
  if (!length(path) || is.na(path) || !nzchar(path))
    stop("could not locate ", file.name)
  utils::read.csv(path, stringsAsFactors = FALSE)
}

## --- data ---
karate.edges <- as.matrix(read.extdata.csv("karate-club-edges.csv"))
karate.n <- max(karate.edges)
labels.df <- read.extdata.csv("karate-club-membership.csv")
karate.club <- labels.df$club[order(labels.df$vertex)]

karate.figure.params <- list(
  preset = "tree",
  repulsion_factor = 4,
  rounds = 160L,
  final_rounds = 224L,
  num_nbrs = 16L,
  r = 0.08,
  s = 6
)

karate.figure <- grip.compare.layouts(
  edges = karate.edges, n = karate.n, dim = 3,
  candidates = list(tuned = karate.figure.params),
  clusters = karate.club,
  seeds = 1:3,
  sample.size.stress = 1000L,
  sample.size.nonedge = 2000L,
  edge.crossings = "never",
  return.layouts = TRUE
)

ok.runs <- karate.figure$runs[karate.figure$runs$status == "ok", , drop = FALSE]
best.run <- ok.runs[order(
  -ok.runs$cluster.separation,
  -ok.runs$sampled.nonedge.sep.ratio,
  ok.runs$sampled.stress,
  ok.runs$seed
), , drop = FALSE][1, , drop = FALSE]

best.coords <- karate.figure$layouts$tuned[[as.character(best.run$seed[[1L]])]]
best.coords <- prcomp(best.coords, center = TRUE, scale. = FALSE)$x
club.cols <- ifelse(karate.club == "Mr. Hi", "#1b9e77", "#d95f02")

## --- figure ---
out.file <- figure_output_path("fig5_karate.png")
png(out.file, width = 2000, height = 900, res = 200)
op <- par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
grip.plot(best.coords, karate.edges, projection = "ortho",
          azimuth = 25, elevation = 20,
          vertex.col = club.cols, main = "View 1")
grip.plot(best.coords, karate.edges, projection = "ortho",
          azimuth = -55, elevation = 30,
          vertex.col = club.cols, main = "View 2")
par(op)
dev.off()
cat("Wrote ", out.file, "\n", sep = "")
