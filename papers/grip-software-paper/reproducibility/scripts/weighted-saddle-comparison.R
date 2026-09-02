#!/usr/bin/env Rscript
# Shared weighted-saddle experiment and independent checks for the paper.
# Source this file from the full benchmark, or update only its saddle component:
# Rscript scripts/weighted-saddle-comparison.R INPUT.rds OUTPUT.rds

check_weighted_saddle <- function(saddle, tolerance = 1e-10) {
  prepared <- saddle$prepared
  stopifnot(saddle$n == 25L, saddle$m == 40L,
            nrow(prepared$pair_matrix) == choose(saddle$n, 2L),
            identical(prepared$tie_mode, "single"))
  edge_lengths <- function(z, edges, epsilon = 0) {
    sqrt(rowSums((z[edges[, 1L], , drop = FALSE] -
                    z[edges[, 2L], , drop = FALSE])^2) + epsilon^2)
  }
  profiled_rmse <- function(observed, input) {
    scale <- sum(observed * input) / sum(input^2)
    sqrt(sum((observed - scale * input)^2) / sum((scale * input)^2))
  }
  # The stored paths, not shortest paths recomputed in the embedding, define
  # every path statistic. Check their original graph lengths as well.
  edge_key <- function(e) paste(pmin(e[, 1L], e[, 2L]), pmax(e[, 1L], e[, 2L]))
  input_lengths <- vapply(prepared$path_edges, function(e) {
    sum(prepared$edge_targets[match(edge_key(e), edge_key(prepared$edges))])
  }, numeric(1L))
  stopifnot(max(abs(input_lengths - prepared$pair_graph_distance)) < tolerance)
  checks <- lapply(names(saddle$layouts), function(method) {
    z <- saddle$layouts[[method]]
    stopifnot(is.matrix(z), identical(dim(z), c(saddle$n, 3L)), all(is.finite(z)))
    edge <- edge_lengths(z, prepared$edges)
    path <- vapply(prepared$path_edges, function(e) {
      sum(edge_lengths(z, e, epsilon = 1e-8))
    }, numeric(1L))
    g <- prepared$pair_graph_distance
    k <- 1 / g^2
    scale <- sum(k * g * path) / sum(k * g^2)
    relative <- (path - scale * g) / (scale * g)
    independent <- c(
      edge.rel.rmse = profiled_rmse(edge, prepared$edge_targets),
      gmds.stress = profiled_rmse(path, g),
      gkk.weighted.rel.rmse = sqrt(sum(k * relative^2) / sum(k)),
      gkk.mean.rel.path.error = mean(abs(relative))
    )
    reported <- unlist(saddle$scores[saddle$scores$method == method,
                                    names(independent)], use.names = FALSE)
    error <- max(abs(independent - reported))
    stopifnot(is.finite(error), error < tolerance)
    data.frame(method = method, max.absolute.score.difference = error)
  })
  # The generating coordinates realize all input edges and retained paths.
  target <- saddle$target_coords
  stopifnot(profiled_rmse(edge_lengths(target, prepared$edges),
                         prepared$edge_targets) < tolerance)
  do.call(rbind, checks)
}

weighted_saddle_comparison <- function(saddle = NULL) {
  for (package in c("grip", "igraph")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Install required package: ", package)
    }
  }
  if (utils::packageVersion("grip") < "0.2.0") stop("grip >= 0.2.0 is required")

  if (is.null(saddle)) {
    graph <- grip::mesh.surface.graph(5, 5, surface = "saddle", amplitude = 0.8)
    ig <- igraph::graph_from_edgelist(graph$edges, directed = FALSE)
    set.seed(1)
    kk <- igraph::layout_with_kk(ig, weights = graph$edge_weights, dim = 3)
    saddle <- list(
      n = graph$n, m = nrow(graph$edges), edges = graph$edges,
      edge_weights = graph$edge_weights, target_coords = graph$coords_surface,
      layouts = list(
        "Combinatorial GRIP" = grip::grip(
          graph$edges, n = graph$n, dim = 3, preset = "mesh", seed = 1),
        "Weighted GRIP" = grip::grip(
          graph$edges, n = graph$n, edge_weights = graph$edge_weights,
          dim = 3, preset = "mesh", metric = "edge_length", seed = 1),
        "Weighted KK (igraph)" = matrix(as.numeric(kk), ncol = 3)
      )
    )
  }
  prepared <- grip::prepare.geodesic.kk(
    saddle$edges, n = saddle$n, edge_weights = saddle$edge_weights)
  weighted <- saddle$layouts[["Weighted GRIP"]]
  mds <- grip::metric.mds(prepared = prepared, dim = 3L)

  # Set the same defaults explicitly for both edge-KK initializers. Do not tune
  # the schedule or iteration budget separately after inspecting the results.
  edge_settings <- list(
    dim = 3L, stiffness_method = "density", stiffness_transform = "identity",
    density_mix_schedule = c(0, 0.25, 0.5, 0.75, 1), density_n = 512L,
    scale_mode = "profiled", max_iter = 50L, initial_step = 1,
    step_shrink = 0.5, armijo_factor = 1e-4, grad_tol = 1e-8,
    min_step = 1e-8, edge_length_epsilon = 1e-8, distance_floor = 1e-8,
    recenter = TRUE, return_trace = TRUE, diagnostics = TRUE,
    seed = 1L, engine = "cpp"
  )
  refine <- function(z) {
    do.call(grip::edge.kk, c(list(coords = z, prepared = prepared), edge_settings))
  }
  weighted_edge <- refine(weighted)
  mds_edge <- refine(mds$coords)

  # Match grip(..., lgkk_polish_rounds = 6), explicitly reusing the displayed
  # weighted-GRIP coordinates rather than generating another initial layout.
  landmark_prepared <- grip::prepare.landmark.geodesic.kk(
    saddle$edges, n = saddle$n, edge_weights = saddle$edge_weights,
    local_nbrs = 20L, landmark_count = 8L)
  lgkk <- grip::landmark.geodesic.kk(
    coords = weighted, prepared = landmark_prepared, max_iter = 6L,
    local_nbrs = 20L, landmark_count = 8L, return_trace = TRUE)
  old_lgkk <- saddle$layouts[["Weighted GRIP + LGKK polish"]]
  old_lgkk_difference <- if (is.null(old_lgkk)) NA_real_ else {
    max(abs(old_lgkk - lgkk$coords))
  }
  saddle$layouts <- list(
    "Combinatorial GRIP" = saddle$layouts[["Combinatorial GRIP"]],
    "Weighted GRIP" = weighted,
    "Weighted GRIP + edge-KK" = weighted_edge$coords,
    "Metric MDS" = mds$coords,
    "Metric MDS + edge-KK" = mds_edge$coords,
    "Weighted GRIP + LGKK polish" = lgkk$coords,
    "Weighted KK (igraph)" = saddle$layouts[["Weighted KK (igraph)"]]
  )
  saddle$prepared <- prepared
  saddle$scores <- do.call(rbind, lapply(names(saddle$layouts), function(method) {
    z <- saddle$layouts[[method]]
    gmds <- grip::score.gmds(z, prepared = prepared)
    cbind(method = method,
          grip::score.geodesic.kk(z, prepared = prepared),
          gmds[, c("edge.rel.rmse", "gmds.stress", "metric.chord.stress")])
  }))
  saddle$refinement_traces <- list(
    weighted_edge = weighted_edge$trace, mds_edge = mds_edge$trace,
    weighted_lgkk = lgkk$trace)
  saddle$comparison_metadata <- list(
    generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    r_version = R.version.string, platform = R.version$platform,
    package_versions = vapply(c("grip", "igraph"), function(p) {
      as.character(utils::packageVersion(p))
    }, character(1L)),
    dimension = 3L, seed = 1L, edge_kk_settings = edge_settings,
    lgkk_settings = list(max_iter = 6L, local_nbrs = 20L, landmark_count = 8L),
    previous_lgkk_max_coordinate_difference = old_lgkk_difference,
    path_policy = "one retained input shortest path for each of 300 unordered pairs",
    baseline_policy = "reuse supplied baseline coordinates when an artifact is supplied"
  )
  saddle$validation <- check_weighted_saddle(saddle)
  saddle
}

if (sys.nframe() == 0L) {
  paper_library <- Sys.getenv("GRIP_RJOURNAL_PACKAGE_LIBRARY")
  if (nzchar(paper_library)) .libPaths(c(paper_library, .libPaths()))
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) stop("Usage: weighted-saddle-comparison.R INPUT.rds OUTPUT.rds")
  original <- readRDS(args[[1L]])
  updated <- original
  updated$weighted_saddle <- weighted_saddle_comparison(original$weighted_saddle)
  untouched <- setdiff(names(original), "weighted_saddle")
  stopifnot(identical(original[untouched], updated[untouched]))
  saveRDS(updated, args[[2L]], compress = "xz")
  restored <- readRDS(args[[2L]])
  stopifnot(identical(updated, restored))
  print(restored$weighted_saddle$scores[, c("method", "edge.rel.rmse",
        "gmds.stress", "gkk.weighted.rel.rmse", "gkk.mean.rel.path.error")])
  print(restored$weighted_saddle$validation)
  print(restored$weighted_saddle$comparison_metadata)
  cat("Verified: all non-saddle components are unchanged.\n")
}
