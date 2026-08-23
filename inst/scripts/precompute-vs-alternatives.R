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

if (is_source_checkout) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Package 'pkgload' is required when running from a source checkout.")
  }
  pkgload::load_all(source_root, quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else {
  library(grip)
}
library(igraph)
library(graphlayouts)

cat("Precomputing benchmark results for grip-vs-alternatives article...\n")

## ---- helpers -------------------------------------------------------

as_coord_matrix <- function(layout_result) {
  matrix(as.numeric(layout_result), ncol = ncol(layout_result))
}

benchmark_repeats <- 5L

time_layout <- function(layout_function, repeats = benchmark_repeats, seed = 1L) {
  elapsed <- numeric(repeats)
  first_coords <- NULL
  for (i in seq_len(repeats)) {
    set.seed(seed)
    invisible(gc(FALSE))
    started <- proc.time()[["elapsed"]]
    coords <- as_coord_matrix(layout_function())
    elapsed[[i]] <- proc.time()[["elapsed"]] - started
    if (is.null(first_coords)) first_coords <- coords
  }
  list(
    coords = first_coords,
    elapsed.sec = elapsed,
    elapsed.sec.median = stats::median(elapsed),
    elapsed.sec.iqr = stats::IQR(elapsed)
  )
}

system_value <- function(command, args) {
  value <- tryCatch(
    system2(command, args, stdout = TRUE, stderr = FALSE),
    error = function(e) character(0L)
  )
  if (!length(value)) NA_character_ else paste(trimws(value), collapse = " ")
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

## ---- Benchmark 4: weighted saddle mesh ------------------------------
cat("  Weighted saddle mesh...\n")

saddle.graph <- mesh.surface.graph(5, 5, surface = "saddle", amplitude = 0.8)
saddle.ig <- graph_from_edgelist(saddle.graph$edges, directed = FALSE)
E(saddle.ig)$weight <- saddle.graph$edge_weights

saddle.grip.combinatorial <- grip(
  saddle.graph$edges, n = saddle.graph$n, dim = 3, preset = "mesh", seed = 1
)
saddle.grip.weighted <- grip(
  saddle.graph$edges, n = saddle.graph$n,
  edge_weights = saddle.graph$edge_weights,
  dim = 3, preset = "mesh", metric = "edge_length", seed = 1
)
saddle.grip.lgkk <- grip(
  saddle.graph$edges, n = saddle.graph$n,
  edge_weights = saddle.graph$edge_weights,
  dim = 3, preset = "mesh", metric = "edge_length",
  lgkk_polish_rounds = 6L, seed = 1
)
set.seed(1)
saddle.kk.weighted <- as_coord_matrix(layout_with_kk(
  saddle.ig,
  weights = saddle.graph$edge_weights,
  dim = 3
))

saddle.prepared <- prepare.geodesic.kk(
  saddle.graph$edges,
  n = saddle.graph$n,
  edge_weights = saddle.graph$edge_weights
)
saddle.layouts <- list(
  "Combinatorial GRIP" = saddle.grip.combinatorial,
  "Weighted GRIP" = saddle.grip.weighted,
  "Weighted GRIP + LGKK polish" = saddle.grip.lgkk,
  "Weighted KK (igraph)" = saddle.kk.weighted
)
saddle.scores <- do.call(rbind, lapply(names(saddle.layouts), function(method) {
  cbind(
    method = method,
    score.geodesic.kk(saddle.layouts[[method]], prepared = saddle.prepared)
  )
}))

results$weighted_saddle <- list(
  n = saddle.graph$n,
  m = nrow(saddle.graph$edges),
  edges = saddle.graph$edges,
  edge_weights = saddle.graph$edge_weights,
  target_coords = saddle.graph$coords_surface,
  scores = saddle.scores,
  layouts = saddle.layouts
)

## ---- Benchmark 5: HMP microbial network ----------------------------
cat("  HMP microbial network...\n")

data(hmp.gc)
hmp.adj   <- hmp.gc$adj_list
hmp.wt    <- hmp.gc$weight_list
hmp.n     <- length(hmp.adj)
hmp.edges <- edge.matrix.from.adj(hmp.adj)
hmp.m     <- nrow(hmp.edges)
## This benchmark isolates topology-oriented layout. The igraph object has no
## edge-weight attribute, all four layout calls therefore ignore the PCA edge
## lengths, and score.layout() evaluates hop-distance stress.
hmp.ig    <- graph_from_edgelist(hmp.edges, directed = FALSE)

hmp.timed <- list(
  "FR (igraph)" = time_layout(function() layout_with_fr(hmp.ig)),
  "DrL (igraph)" = time_layout(function() layout_with_drl(hmp.ig)),
  "Stress (graphlayouts)" = time_layout(function() layout_with_stress(hmp.ig)),
  "grip default (hop)" = time_layout(function() {
    grip(adj_list = hmp.adj, n = hmp.n, dim = 2, seed = 1)
  })
)
hmp.layouts <- lapply(hmp.timed, `[[`, "coords")

hmp.scores <- do.call(rbind, lapply(names(hmp.layouts), function(method) {
  cbind(
    method = method,
    score_layout(
      hmp.layouts[[method]], hmp.edges, hmp.n,
      edge.crossings = "never"
    )
  )
}))

hmp.timing <- data.frame(
  method = names(hmp.timed),
  n.runs = benchmark_repeats,
  elapsed.sec.median = vapply(hmp.timed, `[[`, numeric(1L), "elapsed.sec.median"),
  elapsed.sec.iqr = vapply(hmp.timed, `[[`, numeric(1L), "elapsed.sec.iqr"),
  stringsAsFactors = FALSE
)

results$hmp <- list(
  n = hmp.n, m = hmp.m,
  edges = hmp.edges,
  graph_metric = "hop distance; edge weights ignored by all benchmark methods",
  scores = hmp.scores,
  timing = hmp.timing,
  layouts = hmp.layouts
)

## ---- save ----------------------------------------------------------
results$session_info <- utils::sessionInfo()
blas_path <- extSoftVersion()[["BLAS"]]
if (!is.na(blas_path) && nzchar(blas_path)) {
  r_home <- normalizePath(R.home(), winslash = "/", mustWork = TRUE)
  blas_path <- sub(r_home, "<R_HOME>", blas_path, fixed = TRUE)
}
results$benchmark_metadata <- list(
  repeats = benchmark_repeats,
  timing_statistic = "median elapsed seconds; IQR also retained",
  seed_policy = "each repeat resets the layout RNG seed to 1",
  cpu = system_value("sysctl", c("-n", "machdep.cpu.brand_string")),
  hardware_model = system_value("sysctl", c("-n", "hw.model")),
  memory_bytes = system_value("sysctl", c("-n", "hw.memsize")),
  system = Sys.info(),
  r_version = R.version.string,
  platform = R.version$platform,
  blas = blas_path,
  package_versions = c(
    grip = if (is_source_checkout) {
      unname(read.dcf(description_path, fields = "Version")[[1L]])
    } else {
      as.character(utils::packageVersion("grip"))
    },
    igraph = as.character(utils::packageVersion("igraph")),
    graphlayouts = as.character(utils::packageVersion("graphlayouts"))
  )
)

out.dir <- dirname(out.path)
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
saveRDS(results, out.path, compress = "xz")

cat(sprintf("Saved to %s (%.1f MB)\n", out.path,
            file.size(out.path) / 1e6))
cat("Done.\n")
