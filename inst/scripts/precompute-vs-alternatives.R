#!/usr/bin/env Rscript
## ------------------------------------------------------------------
## Precompute benchmark results for the grip-vs-alternatives article.
##
## Run this script once in an R session that has grip, igraph, and
## graphlayouts installed.  It saves an RDS file that the article
## loads at render time, so the pkgdown build does not need igraph or
## graphlayouts.
##
## Usage:
##   Rscript inst/scripts/precompute-vs-alternatives.R
##
## Set GRIP_VS_ALTERNATIVES_OUTPUT to choose the output file. In a grip
## source checkout the default replaces the bundled package artifact. When
## run from an installed package, the default is benchmark_results.rds in
## the current working directory.
## ------------------------------------------------------------------

script_path <- if (!is.null(sys.frames()[[1]]$ofile)) {
  normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE)
} else {
  normalizePath(
    "inst/scripts/precompute-vs-alternatives.R",
    winslash = "/",
    mustWork = FALSE
  )
}

source_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  winslash = "/",
  mustWork = FALSE
)
description_path <- file.path(source_root, "DESCRIPTION")
is_source_checkout <- file.exists(description_path) &&
  identical(
    unname(read.dcf(description_path, fields = "Package")[[1L]]),
    "grip"
  )

default_output <- if (is_source_checkout) {
  file.path(
    source_root,
    "inst",
    "extdata",
    "vs_alternatives",
    "benchmark_results.rds"
  )
} else {
  file.path(getwd(), "benchmark_results.rds")
}
out.path <- Sys.getenv(
  "GRIP_VS_ALTERNATIVES_OUTPUT",
  unset = default_output
)

library(grip)
library(igraph)
library(graphlayouts)

cat("Precomputing benchmark results for grip-vs-alternatives article...\n")

## ---- helpers -------------------------------------------------------

as_coord_matrix <- function(layout_result) {
  matrix(as.numeric(layout_result), ncol = ncol(layout_result))
}

score_layout <- function(coords, edges, n,
                         clusters = NULL,
                         edge.crossings = "auto",
                         sample.size.stress = 2000L,
                         sample.size.nonedge = 5000L) {
  score.layout(
    coords              = coords,
    edges               = edges,
    n                   = n,
    clusters            = clusters,
    sample.size.stress  = sample.size.stress,
    sample.size.nonedge = sample.size.nonedge,
    stress.seed         = 42L,
    nonedge.seed        = 42L,
    edge.crossings      = edge.crossings
  )
}

edge.matrix.from.adj <- function(adj.list) {
  edges <- list()
  idx <- 0L
  for (u in seq_along(adj.list)) {
    nbrs <- adj.list[[u]]
    nbrs <- nbrs[nbrs > u]
    if (!length(nbrs)) next
    for (v in nbrs) {
      idx <- idx + 1L
      edges[[idx]] <- c(u, v)
    }
  }
  do.call(rbind, edges)
}

results <- list()

## ---- Benchmark 1: Zachary karate club ------------------------------
cat("  Karate club...\n")

karate.edges <- as.matrix(utils::read.csv(
  system.file("extdata", "karate-club-edges.csv", package = "grip"),
  stringsAsFactors = FALSE
))
karate.n <- max(karate.edges)
karate.m <- nrow(karate.edges)

karate.labels <- utils::read.csv(
  system.file("extdata", "karate-club-membership.csv", package = "grip"),
  stringsAsFactors = FALSE
)
karate.labels <- karate.labels[order(karate.labels$vertex), ]
karate.clubs  <- karate.labels$club

karate.ig <- graph_from_edgelist(karate.edges, directed = FALSE)

set.seed(1); karate.fr      <- as_coord_matrix(layout_with_fr(karate.ig))
set.seed(1); karate.kk      <- as_coord_matrix(layout_with_kk(karate.ig))
set.seed(1); karate.drl     <- as_coord_matrix(layout_with_drl(karate.ig))
set.seed(1); karate.stress  <- as_coord_matrix(layout_with_stress(karate.ig))
karate.grip.default <- grip(karate.edges, n = karate.n, dim = 2, seed = 1)
## No preset matches this dense social network, so we use manually tuned
## parameters: more rounds for convergence, moderate repulsion to spread
## the dense core, and barycenter placement.
karate.grip.tuned   <- grip(karate.edges, n = karate.n, dim = 2,
                                    rounds = 128, final_rounds = 192,
                                    num_nbrs = 16, repulsion_factor = 2.0,
                                    r = 0.05, s = 5.0, seed = 1)

karate.scores <- rbind(
  cbind(method = "FR (igraph)",
        score_layout(karate.fr, karate.edges, karate.n,
                     clusters = karate.clubs, edge.crossings = "always")),
  cbind(method = "KK (igraph)",
        score_layout(karate.kk, karate.edges, karate.n,
                     clusters = karate.clubs, edge.crossings = "always")),
  cbind(method = "DrL (igraph)",
        score_layout(karate.drl, karate.edges, karate.n,
                     clusters = karate.clubs, edge.crossings = "always")),
  cbind(method = "Stress (graphlayouts)",
        score_layout(karate.stress, karate.edges, karate.n,
                     clusters = karate.clubs, edge.crossings = "always")),
  cbind(method = "grip default",
        score_layout(karate.grip.default, karate.edges, karate.n,
                     clusters = karate.clubs, edge.crossings = "always")),
  cbind(method = "grip tuned",
        score_layout(karate.grip.tuned, karate.edges, karate.n,
                     clusters = karate.clubs, edge.crossings = "always"))
)

results$karate <- list(
  n = karate.n, m = karate.m,
  clubs = karate.clubs,
  edges = karate.edges,
  scores = karate.scores,
  layouts = list(
    fr = karate.fr, kk = karate.kk, drl = karate.drl,
    stress = karate.stress,
    grip.default = karate.grip.default,
    grip.tuned = karate.grip.tuned
  )
)

## ---- Benchmark 2: 12x12 mesh --------------------------------------
cat("  12x12 mesh...\n")

mesh.edges <- edges.mesh(12, 12)
mesh.n <- 12 * 12
mesh.m <- nrow(mesh.edges)
mesh.ig <- graph_from_edgelist(mesh.edges, directed = FALSE)

set.seed(1); mesh.fr      <- as_coord_matrix(layout_with_fr(mesh.ig))
set.seed(1); mesh.kk      <- as_coord_matrix(layout_with_kk(mesh.ig))
set.seed(1); mesh.drl     <- as_coord_matrix(layout_with_drl(mesh.ig))
set.seed(1); mesh.stress  <- as_coord_matrix(layout_with_stress(mesh.ig))
mesh.grip.default <- grip(mesh.edges, n = mesh.n, dim = 2, seed = 1)
mesh.grip.mesh    <- grip(mesh.edges, n = mesh.n, dim = 2,
                                  preset = "mesh", seed = 1)

mesh.scores <- rbind(
  cbind(method = "FR (igraph)",
        score_layout(mesh.fr, mesh.edges, mesh.n, edge.crossings = "always")),
  cbind(method = "KK (igraph)",
        score_layout(mesh.kk, mesh.edges, mesh.n, edge.crossings = "always")),
  cbind(method = "DrL (igraph)",
        score_layout(mesh.drl, mesh.edges, mesh.n, edge.crossings = "always")),
  cbind(method = "Stress (graphlayouts)",
        score_layout(mesh.stress, mesh.edges, mesh.n, edge.crossings = "always")),
  cbind(method = "grip default",
        score_layout(mesh.grip.default, mesh.edges, mesh.n, edge.crossings = "always")),
  cbind(method = "grip mesh preset",
        score_layout(mesh.grip.mesh, mesh.edges, mesh.n, edge.crossings = "always"))
)

results$mesh <- list(
  n = mesh.n, m = mesh.m,
  edges = mesh.edges,
  scores = mesh.scores,
  layouts = list(
    fr = mesh.fr, kk = mesh.kk, drl = mesh.drl,
    stress = mesh.stress,
    grip.default = mesh.grip.default,
    grip.mesh = mesh.grip.mesh
  )
)

## ---- Benchmark 3: Level-4 Sierpinski carpet ------------------------
cat("  Sierpinski carpet level 4...\n")

carpet.edges <- edges.sierpinski.carpet(4)
carpet.n <- max(carpet.edges)
carpet.m <- nrow(carpet.edges)
carpet.ig <- graph_from_edgelist(carpet.edges, directed = FALSE)

set.seed(1); carpet.time.fr      <- system.time(carpet.fr      <- as_coord_matrix(layout_with_fr(carpet.ig)))
set.seed(1); carpet.time.kk      <- system.time(carpet.kk      <- as_coord_matrix(layout_with_kk(carpet.ig)))
set.seed(1); carpet.time.drl     <- system.time(carpet.drl     <- as_coord_matrix(layout_with_drl(carpet.ig)))
set.seed(1); carpet.time.stress  <- system.time(carpet.stress  <- as_coord_matrix(layout_with_stress(carpet.ig)))
carpet.time.grip.default <- system.time(
  carpet.grip.default <- grip(carpet.edges, n = carpet.n, dim = 2, seed = 1))
carpet.time.grip.carpet  <- system.time(
  carpet.grip.carpet  <- grip(carpet.edges, n = carpet.n, dim = 2, preset = "carpet", seed = 1))

carpet.scores <- rbind(
  cbind(method = "FR (igraph)",
        score_layout(carpet.fr, carpet.edges, carpet.n, edge.crossings = "never")),
  cbind(method = "KK (igraph)",
        score_layout(carpet.kk, carpet.edges, carpet.n, edge.crossings = "never")),
  cbind(method = "DrL (igraph)",
        score_layout(carpet.drl, carpet.edges, carpet.n, edge.crossings = "never")),
  cbind(method = "Stress (graphlayouts)",
        score_layout(carpet.stress, carpet.edges, carpet.n, edge.crossings = "never")),
  cbind(method = "grip default",
        score_layout(carpet.grip.default, carpet.edges, carpet.n, edge.crossings = "never")),
  cbind(method = "grip carpet preset",
        score_layout(carpet.grip.carpet, carpet.edges, carpet.n, edge.crossings = "never"))
)

carpet.timing <- data.frame(
  method = c("FR (igraph)", "KK (igraph)", "DrL (igraph)",
             "Stress (graphlayouts)", "grip default", "grip carpet preset"),
  elapsed.sec = c(carpet.time.fr[["elapsed"]],
                  carpet.time.kk[["elapsed"]],
                  carpet.time.drl[["elapsed"]],
                  carpet.time.stress[["elapsed"]],
                  carpet.time.grip.default[["elapsed"]],
                  carpet.time.grip.carpet[["elapsed"]])
)

results$carpet <- list(
  n = carpet.n, m = carpet.m,
  edges = carpet.edges,
  scores = carpet.scores,
  timing = carpet.timing,
  layouts = list(
    fr = carpet.fr, kk = carpet.kk, drl = carpet.drl,
    stress = carpet.stress,
    grip.default = carpet.grip.default,
    grip.carpet = carpet.grip.carpet
  )
)

## ---- Benchmark 4: HMP microbial network ----------------------------
cat("  HMP microbial network...\n")

data(hmp.u01.gc.coarse)
hmp.adj   <- hmp.u01.gc.coarse$adj_list
hmp.wt    <- hmp.u01.gc.coarse$weight_list
hmp.n     <- length(hmp.adj)
hmp.cst   <- hmp.u01.gc.coarse$vertex_data$cst
hmp.edges <- edge.matrix.from.adj(hmp.adj)
hmp.m     <- nrow(hmp.edges)
hmp.ig    <- graph_from_edgelist(hmp.edges, directed = FALSE)

set.seed(1); hmp.time.fr      <- system.time(hmp.fr     <- as_coord_matrix(layout_with_fr(hmp.ig)))
set.seed(1); hmp.time.drl     <- system.time(hmp.drl    <- as_coord_matrix(layout_with_drl(hmp.ig)))
set.seed(1); hmp.time.stress  <- system.time(hmp.stress <- as_coord_matrix(layout_with_stress(hmp.ig)))
hmp.time.grip.default <- system.time(
  hmp.grip.default <- grip(
    adj_list = hmp.adj, weight_list = hmp.wt, n = hmp.n, dim = 2, seed = 1))

hmp.scores <- rbind(
  cbind(method = "FR (igraph)",
        score_layout(hmp.fr, hmp.edges, hmp.n,
                     clusters = hmp.cst, edge.crossings = "never")),
  cbind(method = "DrL (igraph)",
        score_layout(hmp.drl, hmp.edges, hmp.n,
                     clusters = hmp.cst, edge.crossings = "never")),
  cbind(method = "Stress (graphlayouts)",
        score_layout(hmp.stress, hmp.edges, hmp.n,
                     clusters = hmp.cst, edge.crossings = "never")),
  cbind(method = "grip default",
        score_layout(hmp.grip.default, hmp.edges, hmp.n,
                     clusters = hmp.cst, edge.crossings = "never"))
)

hmp.timing <- data.frame(
  method = c("FR (igraph)", "DrL (igraph)",
             "Stress (graphlayouts)", "grip default"),
  elapsed.sec = c(hmp.time.fr[["elapsed"]],
                  hmp.time.drl[["elapsed"]],
                  hmp.time.stress[["elapsed"]],
                  hmp.time.grip.default[["elapsed"]])
)

results$hmp <- list(
  n = hmp.n, m = hmp.m,
  cst = hmp.cst,
  edges = hmp.edges,
  scores = hmp.scores,
  timing = hmp.timing,
  layouts = list(
    fr = hmp.fr, drl = hmp.drl,
    stress = hmp.stress,
    grip.default = hmp.grip.default
  )
)

## ---- save ----------------------------------------------------------
results$session_info <- utils::sessionInfo()

out.dir <- dirname(out.path)
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
saveRDS(results, out.path, compress = "xz")

cat(sprintf("Saved to %s (%.1f MB)\n", out.path,
            file.size(out.path) / 1e6))
cat("Done.\n")
