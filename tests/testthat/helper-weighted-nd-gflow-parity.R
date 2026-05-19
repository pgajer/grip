`%||%` <- function(x, y) if (is.null(x)) y else x

grip_weighted_nd_gflow_parity_default_repo <- function() {
  Sys.getenv("GFLOW_REPO_DIR", "/Users/pgajer/current_projects/gflow")
}

grip_weighted_nd_gflow_parity_load <- function(
    gflow_repo = grip_weighted_nd_gflow_parity_default_repo()) {
  if (!dir.exists(gflow_repo)) {
    stop("gflow repository not found at GFLOW_REPO_DIR: ", gflow_repo)
  }
  if (!requireNamespace("gflow", quietly = TRUE)) {
    stop("Package `gflow` is required for gflow-backed parity tests.")
  }

  required <- c("FNN", "geometry")
  missing <- required[
    !vapply(required, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing) > 0L) {
    stop("Missing packages required by gflow fixtures: ",
         paste(missing, collapse = ", "))
  }

  env <- new.env(parent = globalenv())
  sys.source(
    file.path(gflow_repo, "dev/graph-trend-filtering/geodesic_annihilator_experiments.R"),
    envir = env
  )
  sys.source(file.path(gflow_repo, "R/grids.R"), envir = env)

  validation_file <- file.path(
    gflow_repo,
    "dev/graph-trend-filtering/local-dimension-validation/compute_local_dimension_validation.R"
  )
  lines <- readLines(validation_file, warn = FALSE)
  start <- which(grepl("^`%\\|\\|%`", lines))[1L]
  eval(parse(text = paste(lines[start:length(lines)], collapse = "\n")), envir = env)
  env
}

grip_weighted_nd_gflow_parity_to_graph <- function(graph) {
  adj <- graph$adj.list %||% graph$adj_list
  weights <- graph$weight.list %||% graph$weight_list
  if (is.null(adj) || is.null(weights)) {
    stop("gflow fixture did not contain adjacency and weight lists.")
  }
  list(
    adj_list = adj,
    weight_list = weights,
    edges = graph$edges,
    edge_weights = graph$edge.weights %||% graph$edge_weights,
    coords = graph$coords %||% NULL,
    param = graph$param %||% NULL
  )
}

grip_weighted_nd_gflow_parity_quadform2d_graph <- function(
    env,
    n = 48L,
    seed = 20260518L,
    index_k = 1L,
    coefficients = c(1.0, 0.65),
    graph_type = c("delaunay", "adaptive_radius", "cknn"),
    k = 8L,
    radius_factor = 1.35,
    cknn_delta = 1.35) {
  graph_type <- match.arg(graph_type)
  U <- env$sample.square.points(n, seed)
  U <- 2 * U - 1
  X <- gflow::quadform.embed(
    U,
    index.k = index_k,
    coefficients = coefficients
  )
  edges <- grip_weighted_nd_gflow_quadform_edges(
    U,
    graph_type = graph_type,
    k = k,
    radius_factor = radius_factor,
    cknn_delta = cknn_delta
  )
  weights <- gflow::quadform.edge.lengths(
    U[edges[, 1L], , drop = FALSE],
    U[edges[, 2L], , drop = FALSE],
    index.k = index_k,
    coefficients = coefficients
  )
  graph <- env$make.graph.from.edges(nrow(U), edges, weights)
  graph$coords <- X
  graph$param <- U
  graph$quadform.index.k <- index_k
  graph$quadform.coefficients <- coefficients
  graph$graph.type <- graph_type
  graph$fixture.note <- sprintf(
    "gflow quadform.embed() plus exact quadform.edge.lengths() on a 2D %s graph",
    graph_type
  )
  graph
}

grip_weighted_nd_gflow_quadform_edges <- function(U,
                                                  graph_type = c("delaunay", "adaptive_radius", "cknn"),
                                                  k = 8L,
                                                  radius_factor = 1.35,
                                                  cknn_delta = 1.35) {
  graph_type <- match.arg(graph_type)
  U <- as.matrix(U)
  if (identical(graph_type, "delaunay")) {
    tri <- geometry::delaunayn(U, options = "QJ")
    edge_mat <- do.call(rbind, lapply(seq_len(nrow(tri)), function(i) {
      cmb <- utils::combn(tri[i, ], 2L)
      t(apply(cmb, 2L, sort))
    }))
    edge_mat <- unique(edge_mat)
    return(edge_mat[order(edge_mat[, 1L], edge_mat[, 2L]), , drop = FALSE])
  }

  D <- as.matrix(stats::dist(U))
  k <- as.integer(k)
  if (is.na(k) || k < 1L || k >= nrow(U)) {
    stop("k must be a positive integer smaller than the number of rows in U")
  }
  kth <- apply(D, 1L, function(z) sort(z)[k + 1L])
  threshold <- if (identical(graph_type, "adaptive_radius")) {
    radius_factor * outer(kth, kth, pmax)
  } else {
    cknn_delta * sqrt(outer(kth, kth))
  }
  edges <- which(upper.tri(D) & D <= threshold, arr.ind = TRUE)
  if (nrow(edges) < 1L) {
    stop(sprintf("%s graph has no edges", graph_type))
  }
  edges[order(edges[, 1L], edges[, 2L]), , drop = FALSE]
}

grip_weighted_nd_gflow_parity_quadform3d_graph <- function(
    env,
    n = 80L,
    seed = 20260518L,
    index_k = 3L,
    coefficients = c(1, 0.75, 0.5),
    graph_type = c("delaunay", "adaptive_radius", "cknn"),
    k = 8L,
    radius_factor = 1.35,
    cknn_delta = 1.35) {
  graph_type <- match.arg(graph_type)
  X.param <- env$sample.ball.validation(n, radius = 0.85, seed = seed)
  X.embed <- gflow::quadform.embed(
    X.param,
    index.k = index_k,
    coefficients = coefficients
  )
  edges <- grip_weighted_nd_gflow_quadform_edges(
    X.param,
    graph_type = graph_type,
    k = k,
    radius_factor = radius_factor,
    cknn_delta = cknn_delta
  )
  weights <- gflow::quadform.edge.lengths(
    X.param[edges[, 1L], , drop = FALSE],
    X.param[edges[, 2L], , drop = FALSE],
    index.k = index_k,
    coefficients = coefficients
  )
  graph <- env$make.graph.from.edges(nrow(X.param), edges, weights)
  graph$coords <- X.embed
  graph$param <- X.param
  graph$quadform.index.k <- index_k
  graph$quadform.coefficients <- coefficients
  graph$domain.shape <- "ball"
  graph$domain.radius <- 0.85
  graph$graph.type <- graph_type
  graph$fixture.note <- sprintf(
    "gflow quadform.embed() plus exact quadform.edge.lengths() on a 3D %s graph",
    graph_type
  )
  graph
}

grip_weighted_nd_gflow_parity_nonuniform3d_graph <- function(
    env,
    n_axis = 4L,
    seed = 20260518L,
    k = 7L) {
  set.seed(seed)
  X <- as.matrix(env$create.3D.rect.grid(
    n = n_axis,
    x1.range = c(0, 1),
    x2.range = c(0, 1),
    x3.range = c(0, 1),
    type = "runif",
    f = 0
  ))
  graph <- env$make.symmetric.knn.graph(X, k = k, metric.coords = X)
  graph$coords <- X
  graph$fixture.note <- "gflow create.3D.rect.grid(type = 'runif') plus symmetric kNN graph"
  graph
}

grip_weighted_nd_gflow_coeff_tag <- function(coefficients) {
  paste(
    gsub("-", "m", gsub("\\.", "p", sprintf("%.2f", coefficients))),
    collapse = "_"
  )
}

grip_weighted_nd_gflow_quadform_case <- function(case_id,
                                                 family,
                                                 intrinsic_dim,
                                                 graph) {
  list(
    case_id = case_id,
    family = family,
    source = "gflow quadform utilities",
    intrinsic_dim = intrinsic_dim,
    graph = grip_weighted_nd_gflow_parity_to_graph(graph)
  )
}

grip_weighted_nd_gflow_stress_quadform_cases <- function(env,
                                                         quad2_n,
                                                         quad3_n) {
  quad_graph_types <- c("delaunay", "adaptive_radius", "cknn")
  quad2_coefficients <- list(
    low = c(1.0, 0.35),
    mid = c(1.0, 0.65),
    high = c(1.0, 1.25)
  )
  quad3_coefficients <- list(
    low = c(1.0, 0.5, 0.25),
    mid = c(1.0, 0.75, 0.5),
    high = c(1.0, 1.25, 0.8)
  )

  quad2 <- unlist(lapply(c(0L, 1L), function(index_k) {
    unlist(lapply(names(quad2_coefficients), function(curvature) {
      lapply(quad_graph_types, function(graph_type) {
        coefficients <- quad2_coefficients[[curvature]]
        grip_weighted_nd_gflow_quadform_case(
          case_id = sprintf(
            "gflow_quadform_2d_%s_index_%d_curv_%s_%s_%d",
            graph_type,
            index_k,
            curvature,
            grip_weighted_nd_gflow_coeff_tag(coefficients),
            quad2_n
          ),
          family = "2D quadform stress",
          intrinsic_dim = 2L,
          graph = grip_weighted_nd_gflow_parity_quadform2d_graph(
            env,
            n = quad2_n,
            index_k = index_k,
            coefficients = coefficients,
            graph_type = graph_type,
            k = 8L
          )
        )
      })
    }), recursive = FALSE)
  }), recursive = FALSE)

  quad3 <- unlist(lapply(c(0L, 1L), function(index_k) {
    unlist(lapply(names(quad3_coefficients), function(curvature) {
      lapply(quad_graph_types, function(graph_type) {
        coefficients <- quad3_coefficients[[curvature]]
        grip_weighted_nd_gflow_quadform_case(
          case_id = sprintf(
            "gflow_quadform_3d_%s_index_%d_curv_%s_%s_%d",
            graph_type,
            index_k,
            curvature,
            grip_weighted_nd_gflow_coeff_tag(coefficients),
            quad3_n
          ),
          family = "3D quadform stress",
          intrinsic_dim = 3L,
          graph = grip_weighted_nd_gflow_parity_quadform3d_graph(
            env,
            n = quad3_n,
            index_k = index_k,
            coefficients = coefficients,
            graph_type = graph_type,
            k = 8L
          )
        )
      })
    }), recursive = FALSE)
  }), recursive = FALSE)

  c(quad2, quad3)
}

grip_weighted_nd_gflow_parity_cases <- function(
    mode = c("smoke", "full", "stress"),
    gflow_repo = grip_weighted_nd_gflow_parity_default_repo()) {
  mode <- match.arg(mode)
  env <- grip_weighted_nd_gflow_parity_load(gflow_repo)

  grid2_n <- if (mode == "smoke") 6L else 9L
  grid3_n <- if (mode == "smoke") 4L else 5L
  quad2_n <- if (mode == "smoke") 44L else 90L
  quad3_n <- if (mode == "smoke") 42L else 80L

  base_cases <- list(
    list(
      case_id = sprintf("gflow_grid_2d_flat_%dx%d", grid2_n, grid2_n),
      family = "flat 2D grid",
      source = "gflow local-dimension validation",
      intrinsic_dim = 2L,
      graph = grip_weighted_nd_gflow_parity_to_graph(
        env$make.grid.validation.graph(grid2_n, grid2_n, "none")
      )
    ),
    list(
      case_id = sprintf("gflow_grid_2d_nonuniform_%dx%d", grid2_n + 1L, grid2_n + 1L),
      family = "nonuniform 2D grid",
      source = "gflow local-dimension validation",
      intrinsic_dim = 2L,
      graph = grip_weighted_nd_gflow_parity_to_graph(
        env$make.grid.validation.graph(grid2_n + 1L, grid2_n + 1L, "nonuniform")
      )
    ),
    list(
      case_id = sprintf("gflow_cube_lattice_%dx%dx%d", grid3_n, grid3_n, grid3_n),
      family = "flat 3D lattice",
      source = "gflow local-dimension validation",
      intrinsic_dim = 3L,
      graph = grip_weighted_nd_gflow_parity_to_graph(
        env$make.cube.lattice.validation.graph(grid3_n, grid3_n, grid3_n)
      )
    ),
    list(
      case_id = sprintf("gflow_grid_3d_nonuniform_runif_%d", grid3_n^3),
      family = "nonuniform 3D grid",
      source = "gflow grid utilities",
      intrinsic_dim = 3L,
      graph = grip_weighted_nd_gflow_parity_to_graph(
        grip_weighted_nd_gflow_parity_nonuniform3d_graph(env, n_axis = grid3_n)
      )
    ),
    list(
      case_id = sprintf("gflow_quadform_2d_delaunay_%d", quad2_n),
      family = "2D quadform",
      source = "gflow quadform utilities",
      intrinsic_dim = 2L,
      graph = grip_weighted_nd_gflow_parity_to_graph(
        grip_weighted_nd_gflow_parity_quadform2d_graph(env, n = quad2_n)
      )
    ),
    list(
      case_id = sprintf("gflow_quadform_3d_delaunay_%d", quad3_n),
      family = "3D quadform",
      source = "gflow local-dimension validation",
      intrinsic_dim = 3L,
      graph = grip_weighted_nd_gflow_parity_to_graph(
        env$make.quadform3d.validation.graph(
          n = quad3_n,
          index.k = 3L,
          coefficients = c(1, 0.75, 0.5),
          domain.shape = "ball",
          domain.radius = 0.85,
          seed = 20260518L,
          graph.type = "delaunay.1skeleton"
        )
      )
    )
  )

  if (identical(mode, "smoke")) {
    return(base_cases)
  }

  if (identical(mode, "stress")) {
    return(grip_weighted_nd_gflow_stress_quadform_cases(
      env = env,
      quad2_n = quad2_n,
      quad3_n = quad3_n
    ))
  }

  quad_graph_types <- c("delaunay", "adaptive_radius", "cknn")
  quad2 <- lapply(c(0L, 1L), function(index_k) {
    lapply(quad_graph_types, function(graph_type) {
      list(
        case_id = sprintf("gflow_quadform_2d_%s_index_%d_%d",
                          graph_type, index_k, quad2_n),
        family = "2D quadform",
        source = "gflow quadform utilities",
        intrinsic_dim = 2L,
        graph = grip_weighted_nd_gflow_parity_to_graph(
          grip_weighted_nd_gflow_parity_quadform2d_graph(
            env,
            n = quad2_n,
            index_k = index_k,
            coefficients = c(1.0, 0.65),
            graph_type = graph_type,
            k = 8L
          )
        )
      )
    })
  })
  quad2 <- unlist(quad2, recursive = FALSE)

  quad3 <- lapply(c(1L, 3L), function(index_k) {
    lapply(quad_graph_types, function(graph_type) {
      list(
        case_id = sprintf("gflow_quadform_3d_%s_index_%d_%d",
                          graph_type, index_k, quad3_n),
        family = "3D quadform",
        source = "gflow quadform utilities",
        intrinsic_dim = 3L,
        graph = grip_weighted_nd_gflow_parity_to_graph(
          grip_weighted_nd_gflow_parity_quadform3d_graph(
            env,
            n = quad3_n,
            index_k = index_k,
            coefficients = c(1, 0.75, 0.5),
            graph_type = graph_type,
            k = 8L
          )
        )
      )
    })
  })
  quad3 <- unlist(quad3, recursive = FALSE)

  c(base_cases[seq_len(4L)], quad2, quad3)
}

grip_weighted_nd_gflow_parity_tuning <- function(mode = c("smoke", "full", "stress")) {
  mode <- match.arg(mode)
  list(
    placement = "barycenter",
    rounds = if (mode == "smoke") 8L else 24L,
    final_rounds = if (mode == "smoke") 12L else 48L,
    num_init = if (mode == "smoke") 10L else 16L,
    num_nbrs = if (mode == "smoke") 10L else 18L,
    r = 0.03,
    s = 6.0,
    repulsion_factor = 1.5,
    length_normalization = "median",
    tinit_factor = 6L,
    seed = 6L,
    disconnected = "components"
  )
}

grip_weighted_nd_gflow_parity_align <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  x_center <- sweep(x, 2L, colMeans(x), "-")
  y_center <- sweep(y, 2L, colMeans(y), "-")
  decomp <- svd(t(y_center) %*% x_center)
  rotation <- decomp$u %*% t(decomp$v)
  y_rot <- y_center %*% rotation
  scale <- sum(x_center * y_rot) / sum(y_rot^2)
  y_aligned <- y_rot * scale
  list(
    x_center = x_center,
    y_center = y_center,
    y_aligned = y_aligned,
    rotation = rotation,
    scale = scale
  )
}

grip_weighted_nd_gflow_parity_edge_rmse <- function(coords, edges, weights) {
  if (is.null(edges) || is.null(weights) || nrow(edges) < 1L) {
    return(NA_real_)
  }
  coords <- as.matrix(coords)
  lengths <- sqrt(rowSums(
    (coords[edges[, 1L], , drop = FALSE] -
       coords[edges[, 2L], , drop = FALSE])^2
  ))
  lengths <- lengths / stats::median(lengths)
  weights <- weights / stats::median(weights)
  sqrt(mean((lengths - weights)^2))
}

grip_weighted_nd_gflow_parity_compare_one <- function(case,
                                                      dim,
                                                      tuning) {
  graph <- case$graph
  common <- c(
    list(
      adj_list = graph$adj_list,
      weight_list = graph$weight_list,
      dim = dim
    ),
    tuning
  )
  legacy_time <- system.time({
    legacy <- do.call(grip.layout.weighted, common)
  })[["elapsed"]]
  nd_time <- system.time({
    nd <- do.call(grip.layout.weighted.nd, common)
  })[["elapsed"]]

  aligned <- grip_weighted_nd_gflow_parity_align(legacy, nd)
  direct_delta <- legacy - nd
  centered_delta <- aligned$x_center - aligned$y_center
  aligned_delta <- aligned$x_center - aligned$y_aligned

  data.frame(
    case_id = case$case_id,
    family = case$family,
    source = case$source,
    intrinsic_dim = case$intrinsic_dim,
    layout_dim = dim,
    n = nrow(legacy),
    m = if (is.null(graph$edges)) NA_integer_ else nrow(graph$edges),
    legacy_elapsed_sec = legacy_time,
    nd_elapsed_sec = nd_time,
    direct_rmse = sqrt(mean(direct_delta^2)),
    direct_max_abs = max(abs(direct_delta)),
    centered_rmse = sqrt(mean(centered_delta^2)),
    centered_max_abs = max(abs(centered_delta)),
    procrustes_rmse = sqrt(mean(aligned_delta^2)),
    procrustes_max_abs = max(abs(aligned_delta)),
    procrustes_scale = aligned$scale,
    legacy_edge_rmse = grip_weighted_nd_gflow_parity_edge_rmse(
      legacy, graph$edges, graph$edge_weights
    ),
    nd_edge_rmse = grip_weighted_nd_gflow_parity_edge_rmse(
      nd, graph$edges, graph$edge_weights
    ),
    finite = all(is.finite(legacy)) && all(is.finite(nd)),
    stringsAsFactors = FALSE
  )
}

grip_weighted_nd_gflow_parity_run <- function(
    mode = c("smoke", "full", "stress"),
    dims = c(2L, 3L),
    gflow_repo = grip_weighted_nd_gflow_parity_default_repo(),
    tuning = NULL) {
  mode <- match.arg(mode)
  cases <- grip_weighted_nd_gflow_parity_cases(mode = mode, gflow_repo = gflow_repo)
  if (is.null(tuning)) {
    tuning <- grip_weighted_nd_gflow_parity_tuning(mode = mode)
  }
  rows <- list()
  i <- 1L
  for (case in cases) {
    for (dim in dims) {
      rows[[i]] <- grip_weighted_nd_gflow_parity_compare_one(case, dim, tuning)
      i <- i + 1L
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

grip_weighted_nd_gflow_parallel_workers <- function() {
  value <- suppressWarnings(as.integer(Sys.getenv("GRIP_GFLOW_PARITY_WORKERS", "")))
  if (length(value) != 1L || is.na(value) || value < 1L) {
    value <- max(1L, min(8L, parallel::detectCores(logical = TRUE) - 1L))
  }
  value
}

grip_weighted_nd_gflow_safe_file_id <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

grip_weighted_nd_gflow_parity_run_parallel <- function(
    mode = c("stress", "full", "smoke"),
    dims = c(2L, 3L),
    gflow_repo = grip_weighted_nd_gflow_parity_default_repo(),
    tuning = NULL,
    workers = grip_weighted_nd_gflow_parallel_workers(),
    output_dir = NULL) {
  mode <- match.arg(mode)
  cases <- grip_weighted_nd_gflow_parity_cases(mode = mode, gflow_repo = gflow_repo)
  if (is.null(tuning)) {
    tuning <- grip_weighted_nd_gflow_parity_tuning(mode = mode)
  }
  tasks <- do.call(rbind, lapply(seq_along(cases), function(case_index) {
    data.frame(
      case_index = case_index,
      dim = as.integer(dims),
      stringsAsFactors = FALSE
    )
  }))
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  run_one <- function(task_index) {
    task <- tasks[task_index, , drop = FALSE]
    case <- cases[[task$case_index[[1L]]]]
    layout_dim <- task$dim[[1L]]
    result <- grip_weighted_nd_gflow_parity_compare_one(
      case = case,
      dim = layout_dim,
      tuning = tuning
    )
    if (!is.null(output_dir)) {
      file_id <- sprintf(
        "%03d_%s_dim_%d.csv",
        task_index,
        grip_weighted_nd_gflow_safe_file_id(case$case_id),
        layout_dim
      )
      utils::write.csv(result, file.path(output_dir, file_id), row.names = FALSE)
    }
    result
  }

  workers <- max(1L, as.integer(workers))
  rows <- if (workers == 1L || nrow(tasks) <= 1L) {
    lapply(seq_len(nrow(tasks)), run_one)
  } else {
    parallel::mclapply(
      seq_len(nrow(tasks)),
      run_one,
      mc.cores = min(workers, nrow(tasks)),
      mc.preschedule = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

grip_weighted_nd_gflow_parity_thresholds <- function() {
  list(
    direct_max_abs = 1e-10,
    procrustes_rmse = 1e-8,
    procrustes_max_abs = 1e-7
  )
}

grip_weighted_nd_gflow_trace_cases <- function(
    mode = c("smoke", "full"),
    gflow_repo = grip_weighted_nd_gflow_parity_default_repo()) {
  mode <- match.arg(mode)
  cases <- grip_weighted_nd_gflow_parity_cases(mode = mode, gflow_repo = gflow_repo)
  keep <- vapply(cases, function(case) {
    case$family %in% c(
      "flat 2D grid",
      "nonuniform 2D grid",
      "flat 3D lattice",
      "nonuniform 3D grid"
    )
  }, logical(1))
  cases[keep]
}

grip_weighted_nd_gflow_trace_tuning <- function(mode = c("smoke", "full")) {
  tuning <- grip_weighted_nd_gflow_parity_tuning(mode = mode)
  tuning$disconnected <- NULL
  tuning$trace.every <- 1L
  tuning
}

grip_weighted_nd_gflow_trace_compare_one <- function(case,
                                                     dim = case$intrinsic_dim,
                                                     tuning = grip_weighted_nd_gflow_trace_tuning("smoke")) {
  graph <- case$graph
  common <- c(
    list(
      adj_list = graph$adj_list,
      weight_list = graph$weight_list,
      dim = as.integer(dim)
    ),
    tuning
  )
  final_anchor_factor <- tuning$final_anchor_factor %||% 0
  final_move_scale_after_first <- tuning$final_move_scale_after_first %||% 1
  legacy_common <- common[
    !names(common) %in% c("final_anchor_factor", "final_move_scale_after_first")
  ]
  legacy_args <- c(
    legacy_common,
    list(
      trace = "round",
      coarse_repulsion_factor = tuning$repulsion_factor,
      coarse_repulsion_sample = 100000L,
      coarse_repulsion_exact_below = 100000L,
      final_anchor_factor = final_anchor_factor,
      final_move_scale_after_first = final_move_scale_after_first,
      final_mode = "fr",
      lgkk_polish_rounds = 0L,
      lgkk_multiscale_rounds = 0L,
      metric_neighbor_cap = NULL,
      diagnostics = "none"
    )
  )
  nd_trace_fn <- get("grip.layout.weighted.nd.trace", asNamespace("grip"))
  legacy <- do.call(grip.layout.trace.weighted, legacy_args)
  nd <- do.call(nd_trace_fn, common)

  frame_count <- min(length(legacy$frames), length(nd$frames))
  rows <- vector("list", frame_count)
  for (idx in seq_len(frame_count)) {
    legacy_frame <- as.matrix(legacy$frames[[idx]])
    nd_frame <- as.matrix(nd$frames[[idx]])
    active <- stats::complete.cases(legacy_frame) &
      stats::complete.cases(nd_frame)
    legacy_active <- legacy_frame[active, , drop = FALSE]
    nd_active <- nd_frame[active, , drop = FALSE]
    direct_delta <- legacy_active - nd_active
    centered_rmse <- NA_real_
    centered_max_abs <- NA_real_
    procrustes_rmse <- NA_real_
    procrustes_max_abs <- NA_real_
    procrustes_scale <- NA_real_
    if (nrow(legacy_active) >= 2L) {
      aligned <- grip_weighted_nd_gflow_parity_align(legacy_active, nd_active)
      centered_delta <- aligned$x_center - aligned$y_center
      aligned_delta <- aligned$x_center - aligned$y_aligned
      centered_rmse <- sqrt(mean(centered_delta^2))
      centered_max_abs <- max(abs(centered_delta))
      procrustes_rmse <- sqrt(mean(aligned_delta^2))
      procrustes_max_abs <- max(abs(aligned_delta))
      procrustes_scale <- aligned$scale
    }
    metadata_cols <- c(
      "phase",
      "level_index",
      "misf_level",
      "round_in_level",
      "active_vertices"
    )
    rows[[idx]] <- data.frame(
      case_id = case$case_id,
      family = case$family,
      source = case$source,
      intrinsic_dim = case$intrinsic_dim,
      layout_dim = as.integer(dim),
      n = nrow(legacy$final),
      frame = idx,
      phase = as.character(legacy$meta$phase[[idx]]),
      nd_phase = as.character(nd$meta$phase[[idx]]),
      level_index = as.integer(legacy$meta$level_index[[idx]]),
      nd_level_index = as.integer(nd$meta$level_index[[idx]]),
      misf_level = as.integer(legacy$meta$misf_level[[idx]]),
      nd_misf_level = as.integer(nd$meta$misf_level[[idx]]),
      round_in_level = as.integer(legacy$meta$round_in_level[[idx]]),
      nd_round_in_level = as.integer(nd$meta$round_in_level[[idx]]),
      active_vertices = sum(active),
      direct_rmse = sqrt(mean(direct_delta^2)),
      direct_max_abs = max(abs(direct_delta)),
      centered_rmse = centered_rmse,
      centered_max_abs = centered_max_abs,
      procrustes_rmse = procrustes_rmse,
      procrustes_max_abs = procrustes_max_abs,
      procrustes_scale = procrustes_scale,
      metadata_match = identical(
        legacy$meta[idx, metadata_cols],
        nd$meta[idx, metadata_cols]
      ),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  out$frame_count_match <- length(legacy$frames) == length(nd$frames)
  rownames(out) <- NULL
  out
}

grip_weighted_nd_gflow_trace_run <- function(
    mode = c("smoke", "full"),
    dims = c("intrinsic", "all"),
    gflow_repo = grip_weighted_nd_gflow_parity_default_repo(),
    tuning = NULL) {
  mode <- match.arg(mode)
  dims <- match.arg(dims)
  cases <- grip_weighted_nd_gflow_trace_cases(mode = mode, gflow_repo = gflow_repo)
  if (is.null(tuning)) {
    tuning <- grip_weighted_nd_gflow_trace_tuning(mode = mode)
  }
  rows <- list()
  i <- 1L
  for (case in cases) {
    case_dims <- if (identical(dims, "all")) c(2L, 3L) else case$intrinsic_dim
    for (dim in case_dims) {
      rows[[i]] <- grip_weighted_nd_gflow_trace_compare_one(
        case = case,
        dim = dim,
        tuning = tuning
      )
      i <- i + 1L
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

grip_weighted_nd_gflow_trace_first_divergence <- function(metrics,
                                                          tolerance = 1e-8) {
  split_metrics <- split(
    metrics,
    list(metrics$case_id, metrics$layout_dim),
    drop = TRUE
  )
  rows <- lapply(split_metrics, function(df) {
    finite <- is.finite(df$procrustes_rmse)
    divergent <- finite & df$procrustes_rmse > tolerance
    if (!any(divergent)) {
      idx <- nrow(df)
      status <- "within_tolerance"
    } else {
      idx <- which(divergent)[[1L]]
      status <- "diverged"
    }
    data.frame(
      case_id = df$case_id[[idx]],
      family = df$family[[idx]],
      layout_dim = df$layout_dim[[idx]],
      first_divergent_frame = if (identical(status, "diverged")) df$frame[[idx]] else NA_integer_,
      first_divergent_phase = if (identical(status, "diverged")) df$phase[[idx]] else NA_character_,
      first_divergent_level_index = if (identical(status, "diverged")) df$level_index[[idx]] else NA_integer_,
      first_divergent_misf_level = if (identical(status, "diverged")) df$misf_level[[idx]] else NA_integer_,
      first_divergent_round = if (identical(status, "diverged")) df$round_in_level[[idx]] else NA_integer_,
      first_divergent_active_vertices = if (identical(status, "diverged")) df$active_vertices[[idx]] else NA_integer_,
      first_divergent_procrustes_rmse = if (identical(status, "diverged")) df$procrustes_rmse[[idx]] else NA_real_,
      first_divergent_direct_max_abs = if (identical(status, "diverged")) df$direct_max_abs[[idx]] else NA_real_,
      final_procrustes_rmse = utils::tail(df$procrustes_rmse, 1L),
      metadata_all_match = all(df$metadata_match) && all(df$frame_count_match),
      status = status,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$case_id, out$layout_dim), , drop = FALSE]
}

grip_weighted_nd_gflow_refinement_step_trace_target <- function(
    mode = c("smoke", "full"),
    gflow_repo = grip_weighted_nd_gflow_parity_default_repo(),
    tolerance = 1e-8) {
  mode <- match.arg(mode)
  metrics <- grip_weighted_nd_gflow_trace_run(
    mode = mode,
    dims = "all",
    gflow_repo = gflow_repo
  )
  first <- grip_weighted_nd_gflow_trace_first_divergence(
    metrics,
    tolerance = tolerance
  )
  divergent <- first[first$status == "diverged", , drop = FALSE]
  if (nrow(divergent) == 0L) {
    stop("No divergent gflow trace frame found at tolerance ", tolerance)
  }
  divergent <- divergent[order(
    divergent$first_divergent_procrustes_rmse,
    decreasing = TRUE
  ), , drop = FALSE]
  target <- divergent[1L, , drop = FALSE]
  target$round_start <- pmax(1L, target$first_divergent_round - 2L)
  target$round_end <- target$first_divergent_round
  target
}

grip_weighted_nd_gflow_refinement_step_matrix <- function(trace_list, name) {
  mat <- as.matrix(trace_list[[name]])
  colnames(mat) <- paste0(name, "_Dim", seq_len(ncol(mat)))
  mat
}

grip_weighted_nd_gflow_refinement_step_terms_flatten <- function(trace_list,
                                                                 backend,
                                                                 case_id,
                                                                 family,
                                                                 layout_dim) {
  terms <- as.data.frame(trace_list$attraction_terms)
  if (nrow(terms) == 0L) {
    return(data.frame())
  }
  parent <- as.data.frame(trace_list$trace)[terms$parent_row, , drop = FALSE]
  cbind(
    data.frame(
      backend = backend,
      case_id = case_id,
      family = family,
      layout_dim = layout_dim,
      level_index = parent$level_index,
      misf_level = parent$misf_level,
      round_in_level = parent$round_in_level,
      active_vertices = parent$active_vertices,
      order_index = parent$order_index,
      stringsAsFactors = FALSE
    ),
    terms,
    grip_weighted_nd_gflow_refinement_step_matrix(trace_list, "attraction_term_delta"),
    grip_weighted_nd_gflow_refinement_step_matrix(trace_list, "attraction_term_step"),
    grip_weighted_nd_gflow_refinement_step_matrix(trace_list, "attraction_term_cumulative")
  )
}

grip_weighted_nd_gflow_refinement_step_flatten <- function(trace_list,
                                                           backend,
                                                           case_id,
                                                           family,
                                                           layout_dim) {
  trace <- as.data.frame(trace_list$trace)
  if (nrow(trace) == 0L) {
    return(data.frame())
  }
  cbind(
    data.frame(
      backend = backend,
      case_id = case_id,
      family = family,
      layout_dim = layout_dim,
      stringsAsFactors = FALSE
    ),
    trace,
    grip_weighted_nd_gflow_refinement_step_matrix(trace_list, "coords_before"),
    grip_weighted_nd_gflow_refinement_step_matrix(trace_list, "pre_temp_disp"),
    grip_weighted_nd_gflow_refinement_step_matrix(trace_list, "attraction_disp"),
    grip_weighted_nd_gflow_refinement_step_matrix(trace_list, "repulsion_disp"),
    grip_weighted_nd_gflow_refinement_step_matrix(trace_list, "applied_disp"),
    grip_weighted_nd_gflow_refinement_step_matrix(trace_list, "coords_after")
  )
}

grip_weighted_nd_gflow_refinement_step_compare <- function(legacy,
                                                           nd,
                                                           case_id,
                                                           family,
                                                           layout_dim) {
  legacy_flat <- grip_weighted_nd_gflow_refinement_step_flatten(
    legacy,
    backend = "legacy",
    case_id = case_id,
    family = family,
    layout_dim = layout_dim
  )
  nd_flat <- grip_weighted_nd_gflow_refinement_step_flatten(
    nd,
    backend = "nd",
    case_id = case_id,
    family = family,
    layout_dim = layout_dim
  )
  if (nrow(legacy_flat) == 0L || nrow(nd_flat) == 0L) {
    stop("Refinement step trace did not capture any rows.")
  }

  key <- c("level_index", "misf_level", "round_in_level", "order_index", "vertex")
  legacy_flat <- legacy_flat[do.call(order, legacy_flat[key]), , drop = FALSE]
  nd_flat <- nd_flat[do.call(order, nd_flat[key]), , drop = FALSE]
  if (!identical(legacy_flat[key], nd_flat[key])) {
    stop("Legacy and ND refinement step trace keys do not match.")
  }

  vector_groups <- list(
    coords_before = grep("^coords_before_Dim", names(legacy_flat), value = TRUE),
    pre_temp_disp = grep("^pre_temp_disp_Dim", names(legacy_flat), value = TRUE),
    attraction_disp = grep("^attraction_disp_Dim", names(legacy_flat), value = TRUE),
    repulsion_disp = grep("^repulsion_disp_Dim", names(legacy_flat), value = TRUE),
    applied_disp = grep("^applied_disp_Dim", names(legacy_flat), value = TRUE),
    coords_after = grep("^coords_after_Dim", names(legacy_flat), value = TRUE)
  )
  scalar_cols <- c(
    "heat_before",
    "heat_after",
    "old_cos_before",
    "old_cos_after",
    "old_disp_norm_before",
    "pre_temp_disp_norm"
  )
  out <- legacy_flat[c(
    "case_id",
    "family",
    "layout_dim",
    "level_index",
    "misf_level",
    "round_in_level",
    "active_vertices",
    "order_index",
    "vertex"
  )]
  for (col in scalar_cols) {
    out[[paste0(col, "_legacy")]] <- legacy_flat[[col]]
    out[[paste0(col, "_nd")]] <- nd_flat[[col]]
    out[[paste0(col, "_delta")]] <- legacy_flat[[col]] - nd_flat[[col]]
  }
  out$attraction_edges_legacy <- legacy_flat$attraction_edges
  out$attraction_edges_nd <- nd_flat$attraction_edges
  out$attraction_edges_match <- legacy_flat$attraction_edges == nd_flat$attraction_edges
  out$repulsion_neighbors_legacy <- legacy_flat$repulsion_neighbors
  out$repulsion_neighbors_nd <- nd_flat$repulsion_neighbors
  out$repulsion_neighbors_match <- legacy_flat$repulsion_neighbors == nd_flat$repulsion_neighbors
  for (group in names(vector_groups)) {
    cols <- vector_groups[[group]]
    delta <- as.matrix(legacy_flat[cols]) - as.matrix(nd_flat[cols])
    out[[paste0(group, "_rmse")]] <- sqrt(rowMeans(delta^2))
    out[[paste0(group, "_max_abs")]] <- apply(abs(delta), 1L, max)
  }
  rownames(out) <- NULL
  out
}

grip_weighted_nd_gflow_refinement_step_term_compare <- function(legacy,
                                                                nd,
                                                                case_id,
                                                                family,
                                                                layout_dim) {
  legacy_terms <- grip_weighted_nd_gflow_refinement_step_terms_flatten(
    legacy,
    backend = "legacy",
    case_id = case_id,
    family = family,
    layout_dim = layout_dim
  )
  nd_terms <- grip_weighted_nd_gflow_refinement_step_terms_flatten(
    nd,
    backend = "nd",
    case_id = case_id,
    family = family,
    layout_dim = layout_dim
  )
  if (nrow(legacy_terms) == 0L || nrow(nd_terms) == 0L) {
    stop("Attraction term trace did not capture any rows.")
  }

  key <- c("parent_row", "term_index", "vertex", "neighbor")
  legacy_terms <- legacy_terms[do.call(order, legacy_terms[key]), , drop = FALSE]
  nd_terms <- nd_terms[do.call(order, nd_terms[key]), , drop = FALSE]
  if (!identical(legacy_terms[key], nd_terms[key])) {
    stop("Legacy and ND attraction term trace keys do not match.")
  }

  vector_groups <- list(
    delta = grep("^attraction_term_delta_Dim", names(legacy_terms), value = TRUE),
    step = grep("^attraction_term_step_Dim", names(legacy_terms), value = TRUE),
    cumulative = grep("^attraction_term_cumulative_Dim", names(legacy_terms), value = TRUE)
  )
  scalar_cols <- c("weight", "norm2", "desired", "desired2", "scale")
  out <- legacy_terms[c(
    "case_id",
    "family",
    "layout_dim",
    "level_index",
    "misf_level",
    "round_in_level",
    "active_vertices",
    "order_index",
    "parent_row",
    "term_index",
    "vertex",
    "neighbor"
  )]
  for (col in scalar_cols) {
    out[[paste0(col, "_legacy")]] <- legacy_terms[[col]]
    out[[paste0(col, "_nd")]] <- nd_terms[[col]]
    out[[paste0(col, "_delta")]] <- legacy_terms[[col]] - nd_terms[[col]]
  }
  for (group in names(vector_groups)) {
    cols <- vector_groups[[group]]
    delta <- as.matrix(legacy_terms[cols]) - as.matrix(nd_terms[cols])
    out[[paste0(group, "_rmse")]] <- sqrt(rowMeans(delta^2))
    out[[paste0(group, "_max_abs")]] <- apply(abs(delta), 1L, max)
  }
  rownames(out) <- NULL
  out
}

grip_weighted_nd_gflow_refinement_step_trace_run <- function(
    mode = c("smoke", "full"),
    gflow_repo = grip_weighted_nd_gflow_parity_default_repo(),
    target = NULL,
    tolerance = 1e-8) {
  mode <- match.arg(mode)
  if (is.null(target)) {
    target <- grip_weighted_nd_gflow_refinement_step_trace_target(
      mode = mode,
      gflow_repo = gflow_repo,
      tolerance = tolerance
    )
  }
  cases <- grip_weighted_nd_gflow_trace_cases(mode = mode, gflow_repo = gflow_repo)
  if (!any(vapply(cases, function(x) {
    identical(x$case_id, target$case_id[[1L]])
  }, logical(1)))) {
    cases <- grip_weighted_nd_gflow_parity_cases(mode = mode, gflow_repo = gflow_repo)
  }
  case <- cases[vapply(cases, function(x) {
    identical(x$case_id, target$case_id[[1L]])
  }, logical(1))][[1L]]
  tuning <- grip_weighted_nd_gflow_trace_tuning(mode = mode)
  graph <- case$graph
  layout_dim <- as.integer(target$layout_dim[[1L]])

  legacy_validate <- get("grip.validate.weighted.layout.inputs", asNamespace("grip"))
  nd_validate <- get("grip.validate.weighted.nd.layout.inputs", asNamespace("grip"))
  legacy_cpp <- get("grip_layout_globalrep_weighted_trace_adj_cpp", asNamespace("grip"))
  nd_cpp <- get("grip_layout_weighted_nd_trace_adj_cpp", asNamespace("grip"))

  legacy_input <- legacy_validate(
    adj_list = graph$adj_list,
    weight_list = graph$weight_list,
    dim = layout_dim,
    placement = "barycenter",
    seed = tuning$seed,
    length_normalization = tuning$length_normalization,
    caller = "grip_weighted_nd_gflow_refinement_step_trace_run"
  )
  nd_input <- nd_validate(
    adj_list = graph$adj_list,
    weight_list = graph$weight_list,
    dim = layout_dim,
    seed = tuning$seed,
    length_normalization = tuning$length_normalization,
    caller = "grip_weighted_nd_gflow_refinement_step_trace_run"
  )

  level_index <- as.integer(target$first_divergent_level_index[[1L]])
  misf_level <- as.integer(target$first_divergent_misf_level[[1L]])
  round_start <- as.integer(target$round_start[[1L]])
  round_end <- as.integer(target$round_end[[1L]])

  legacy <- legacy_cpp(
    adj_list = legacy_input$adj_list,
    weight_list = legacy_input$weight_list,
    n = legacy_input$n,
    dim = legacy_input$dim,
    placement = "barycenter",
    rounds = as.integer(tuning$rounds),
    final_rounds = as.integer(tuning$final_rounds),
    num_init = as.integer(tuning$num_init),
    num_nbrs = as.integer(tuning$num_nbrs),
    r = as.double(tuning$r),
    s = as.double(tuning$s),
    repulsion_factor = as.double(tuning$repulsion_factor),
    coarse_repulsion_factor = as.double(tuning$repulsion_factor),
    coarse_repulsion_sample = 100000L,
    coarse_repulsion_exact_below = 100000L,
    final_anchor_factor = 0,
    final_move_scale_after_first = 1,
    insertion_anchor_count = 3L,
    insertion_anchor_scope = "any_higher",
    insertion_anchor_strategy = "first",
    level0_insertion_mode = "inherit",
    level0_anchor_count = 3L,
    level0_local_kk_steps = 3L,
    lgkk_multiscale_rounds = 0L,
    lgkk_rounds_coarse = 0L,
    lgkk_rounds_pre_final = 0L,
    lgkk_rounds_final = 0L,
    lgkk_local_nbrs = 20L,
    lgkk_landmark_count = 8L,
    lgkk_multiscale_scope = "all",
    lgkk_active_limit = 4096L,
    final_mode = "fr",
    tinit_factor = as.integer(tuning$tinit_factor),
    seed = legacy_input$seed,
    trace = "round",
    trace_every = 1L,
    metric_neighbor_cap = 0L,
    refinement_step_trace = TRUE,
    refinement_step_level_index = level_index,
    refinement_step_misf_level = misf_level,
    refinement_step_round_start = round_start,
    refinement_step_round_end = round_end
  )
  nd <- nd_cpp(
    adj_list = nd_input$adj_list,
    weight_list = nd_input$weight_list,
    n = nd_input$n,
    dim = nd_input$dim,
    placement = "barycenter",
    rounds = as.integer(tuning$rounds),
    final_rounds = as.integer(tuning$final_rounds),
    num_init = as.integer(tuning$num_init),
    num_nbrs = as.integer(tuning$num_nbrs),
    r = as.double(tuning$r),
    s = as.double(tuning$s),
    repulsion_factor = as.double(tuning$repulsion_factor),
    tinit_factor = as.integer(tuning$tinit_factor),
    final_anchor_factor = 0,
    final_move_scale_after_first = 1,
    final_mode = "fr",
    metric_neighbor_cap = 0L,
    insertion_anchor_count = 3L,
    insertion_anchor_scope = "any_higher",
    insertion_anchor_strategy = "first",
    level0_insertion_mode = "inherit",
    level0_anchor_count = 3L,
    level0_local_kk_steps = 3L,
    seed = nd_input$seed,
    trace_every = 1L,
    refinement_step_trace = TRUE,
    refinement_step_level_index = level_index,
    refinement_step_misf_level = misf_level,
    refinement_step_round_start = round_start,
    refinement_step_round_end = round_end
  )

  comparison <- grip_weighted_nd_gflow_refinement_step_compare(
    legacy$refinement_step_trace,
    nd$refinement_step_trace,
    case_id = case$case_id,
    family = case$family,
    layout_dim = layout_dim
  )
  terms <- grip_weighted_nd_gflow_refinement_step_term_compare(
    legacy$refinement_step_trace,
    nd$refinement_step_trace,
    case_id = case$case_id,
    family = case$family,
    layout_dim = layout_dim
  )
  attr(comparison, "target") <- target
  attr(comparison, "attraction_terms") <- terms
  comparison
}

grip_weighted_nd_gflow_insertion_trace_compare <- function(legacy,
                                                           nd,
                                                           case_id,
                                                           family,
                                                           layout_dim) {
  legacy_trace <- legacy$insertion_trace$trace
  nd_trace <- nd$insertion_trace$trace
  n <- min(nrow(legacy_trace), nrow(nd_trace))
  if (n < 1L) {
    return(data.frame())
  }
  legacy_trace <- legacy_trace[seq_len(n), , drop = FALSE]
  nd_trace <- nd_trace[seq_len(n), , drop = FALSE]
  legacy_initial <- as.matrix(legacy$insertion_trace$coords_initial)[seq_len(n), , drop = FALSE]
  nd_initial <- as.matrix(nd$insertion_trace$coords_initial)[seq_len(n), , drop = FALSE]
  legacy_after <- as.matrix(legacy$insertion_trace$coords_after)[seq_len(n), , drop = FALSE]
  nd_after <- as.matrix(nd$insertion_trace$coords_after)[seq_len(n), , drop = FALSE]
  legacy_old_initial <- as.matrix(legacy$insertion_trace$old_disp_initial)[seq_len(n), , drop = FALSE]
  nd_old_initial <- as.matrix(nd$insertion_trace$old_disp_initial)[seq_len(n), , drop = FALSE]
  legacy_old_after <- as.matrix(legacy$insertion_trace$old_disp_after)[seq_len(n), , drop = FALSE]
  nd_old_after <- as.matrix(nd$insertion_trace$old_disp_after)[seq_len(n), , drop = FALSE]

  vec_rmse <- function(x, y) sqrt(rowMeans((x - y)^2))
  vec_max <- function(x, y) apply(abs(x - y), 1L, max)
  metadata_cols <- c(
    "level_index",
    "misf_level",
    "previous_active_vertices",
    "active_vertices",
    "order_index",
    "vertex",
    "root_depth",
    "anchor_count_requested",
    "anchor_count_used",
    "insertion_mode",
    "local_kk_steps"
  )

  data.frame(
    case_id = case_id,
    family = family,
    layout_dim = as.integer(layout_dim),
    insertion_row = seq_len(n),
    level_index = legacy_trace$level_index,
    misf_level = legacy_trace$misf_level,
    active_vertices = legacy_trace$active_vertices,
    order_index = legacy_trace$order_index,
    vertex = legacy_trace$vertex,
    root_depth = legacy_trace$root_depth,
    anchor_count_requested = legacy_trace$anchor_count_requested,
    anchor_count_used = legacy_trace$anchor_count_used,
    insertion_mode = legacy_trace$insertion_mode,
    local_kk_steps = legacy_trace$local_kk_steps,
    metadata_match = vapply(seq_len(n), function(i) {
      identical(legacy_trace[i, metadata_cols], nd_trace[i, metadata_cols])
    }, logical(1)),
    anchor_match = legacy_trace$anchors == nd_trace$anchors,
    coords_initial_rmse = vec_rmse(legacy_initial, nd_initial),
    coords_initial_max_abs = vec_max(legacy_initial, nd_initial),
    coords_after_rmse = vec_rmse(legacy_after, nd_after),
    coords_after_max_abs = vec_max(legacy_after, nd_after),
    old_disp_initial_rmse = vec_rmse(legacy_old_initial, nd_old_initial),
    old_disp_initial_max_abs = vec_max(legacy_old_initial, nd_old_initial),
    old_disp_after_rmse = vec_rmse(legacy_old_after, nd_old_after),
    old_disp_after_max_abs = vec_max(legacy_old_after, nd_old_after),
    old_disp_norm_initial_delta =
      legacy_trace$old_disp_norm_initial - nd_trace$old_disp_norm_initial,
    old_disp_norm_after_delta =
      legacy_trace$old_disp_norm_after - nd_trace$old_disp_norm_after,
    legacy_anchors = legacy_trace$anchors,
    nd_anchors = nd_trace$anchors,
    stringsAsFactors = FALSE
  )
}

grip_weighted_nd_gflow_insertion_trace_run <- function(
    mode = c("smoke", "full"),
    case_id = NULL,
    dim = 2L,
    gflow_repo = grip_weighted_nd_gflow_parity_default_repo()) {
  mode <- match.arg(mode)
  cases <- grip_weighted_nd_gflow_parity_cases(mode = mode, gflow_repo = gflow_repo)
  if (is.null(case_id)) {
    case_id <- cases[[1L]]$case_id
  }
  keep <- vapply(cases, function(x) identical(x$case_id, case_id), logical(1))
  if (!any(keep)) {
    stop("case_id not found: ", case_id)
  }
  case <- cases[keep][[1L]]
  tuning <- grip_weighted_nd_gflow_parity_tuning(mode = mode)
  graph <- case$graph

  legacy_validate <- get("grip.validate.weighted.layout.inputs", asNamespace("grip"))
  nd_validate <- get("grip.validate.weighted.nd.layout.inputs", asNamespace("grip"))
  legacy_cpp <- get("grip_layout_globalrep_weighted_trace_adj_cpp", asNamespace("grip"))
  nd_cpp <- get("grip_layout_weighted_nd_trace_adj_cpp", asNamespace("grip"))

  legacy_input <- legacy_validate(
    adj_list = graph$adj_list,
    weight_list = graph$weight_list,
    dim = as.integer(dim),
    placement = "barycenter",
    seed = tuning$seed,
    length_normalization = tuning$length_normalization,
    caller = "grip_weighted_nd_gflow_insertion_trace_run"
  )
  nd_input <- nd_validate(
    adj_list = graph$adj_list,
    weight_list = graph$weight_list,
    dim = as.integer(dim),
    seed = tuning$seed,
    length_normalization = tuning$length_normalization,
    caller = "grip_weighted_nd_gflow_insertion_trace_run"
  )

  legacy <- legacy_cpp(
    adj_list = legacy_input$adj_list,
    weight_list = legacy_input$weight_list,
    n = legacy_input$n,
    dim = legacy_input$dim,
    placement = "barycenter",
    rounds = as.integer(tuning$rounds),
    final_rounds = as.integer(tuning$final_rounds),
    num_init = as.integer(tuning$num_init),
    num_nbrs = as.integer(tuning$num_nbrs),
    r = as.double(tuning$r),
    s = as.double(tuning$s),
    repulsion_factor = as.double(tuning$repulsion_factor),
    coarse_repulsion_factor = as.double(tuning$repulsion_factor),
    coarse_repulsion_sample = 100000L,
    coarse_repulsion_exact_below = 100000L,
    final_anchor_factor = 0,
    final_move_scale_after_first = 1,
    insertion_anchor_count = 3L,
    insertion_anchor_scope = "any_higher",
    insertion_anchor_strategy = "first",
    level0_insertion_mode = "inherit",
    level0_anchor_count = 3L,
    level0_local_kk_steps = 3L,
    lgkk_multiscale_rounds = 0L,
    lgkk_rounds_coarse = 0L,
    lgkk_rounds_pre_final = 0L,
    lgkk_rounds_final = 0L,
    lgkk_local_nbrs = 20L,
    lgkk_landmark_count = 8L,
    lgkk_multiscale_scope = "all",
    lgkk_active_limit = 4096L,
    final_mode = "fr",
    tinit_factor = as.integer(tuning$tinit_factor),
    seed = legacy_input$seed,
    trace = "level",
    trace_every = 1L,
    metric_neighbor_cap = 0L,
    refinement_step_trace = TRUE,
    refinement_step_level_index = -1L,
    refinement_step_misf_level = -1L,
    refinement_step_round_start = -1L,
    refinement_step_round_end = -1L
  )
  nd <- nd_cpp(
    adj_list = nd_input$adj_list,
    weight_list = nd_input$weight_list,
    n = nd_input$n,
    dim = nd_input$dim,
    placement = "barycenter",
    rounds = as.integer(tuning$rounds),
    final_rounds = as.integer(tuning$final_rounds),
    num_init = as.integer(tuning$num_init),
    num_nbrs = as.integer(tuning$num_nbrs),
    r = as.double(tuning$r),
    s = as.double(tuning$s),
    repulsion_factor = as.double(tuning$repulsion_factor),
    tinit_factor = as.integer(tuning$tinit_factor),
    final_anchor_factor = 0,
    final_move_scale_after_first = 1,
    final_mode = "fr",
    metric_neighbor_cap = 0L,
    insertion_anchor_count = 3L,
    insertion_anchor_scope = "any_higher",
    insertion_anchor_strategy = "first",
    level0_insertion_mode = "inherit",
    level0_anchor_count = 3L,
    level0_local_kk_steps = 3L,
    seed = nd_input$seed,
    trace_every = 1L,
    refinement_step_trace = TRUE,
    refinement_step_level_index = -1L,
    refinement_step_misf_level = -1L,
    refinement_step_round_start = -1L,
    refinement_step_round_end = -1L
  )

  grip_weighted_nd_gflow_insertion_trace_compare(
    legacy = legacy,
    nd = nd,
    case_id = case$case_id,
    family = case$family,
    layout_dim = as.integer(dim)
  )
}

grip_weighted_nd_gflow_refinement_step_summary <- function(comparison) {
  metric_cols <- grep("(_delta|_rmse|_max_abs)$", names(comparison), value = TRUE)
  rows <- lapply(metric_cols, function(col) {
    values <- abs(comparison[[col]])
    idx <- which.max(values)
    data.frame(
      metric = col,
      max_abs = values[[idx]],
      round_in_level = comparison$round_in_level[[idx]],
      vertex = comparison$vertex[[idx]],
      order_index = comparison$order_index[[idx]],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$max_abs, decreasing = TRUE), , drop = FALSE]
}

grip_weighted_nd_gflow_refinement_step_term_summary <- function(comparison) {
  terms <- attr(comparison, "attraction_terms")
  if (is.null(terms) || nrow(terms) == 0L) {
    stop("comparison does not include attraction term diagnostics")
  }
  metric_cols <- grep("(_delta|_rmse|_max_abs)$", names(terms), value = TRUE)
  rows <- lapply(metric_cols, function(col) {
    values <- abs(terms[[col]])
    idx <- which.max(values)
    data.frame(
      metric = col,
      max_abs = values[[idx]],
      round_in_level = terms$round_in_level[[idx]],
      vertex = terms$vertex[[idx]],
      order_index = terms$order_index[[idx]],
      term_index = terms$term_index[[idx]],
      neighbor = terms$neighbor[[idx]],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$max_abs, decreasing = TRUE), , drop = FALSE]
}

grip_weighted_nd_gflow_refinement_source_target <- function(comparison) {
  target <- attr(comparison, "target")
  terms <- attr(comparison, "attraction_terms")
  if (is.null(target) || is.null(terms) || nrow(terms) == 0L) {
    stop("comparison must include target and attraction term diagnostics")
  }
  failing_round <- as.integer(target$first_divergent_round[[1L]])
  step_rows <- comparison[
    comparison$round_in_level == failing_round,
    ,
    drop = FALSE
  ]
  if (nrow(step_rows) == 0L) {
    stop("comparison does not contain the first divergent round")
  }
  step_idx <- which.max(step_rows$attraction_disp_max_abs)
  focal_vertex <- step_rows$vertex[[step_idx]]
  focal_order <- step_rows$order_index[[step_idx]]
  term_rows <- terms[
    terms$round_in_level == failing_round & terms$vertex == focal_vertex,
    ,
    drop = FALSE
  ]
  if (nrow(term_rows) == 0L) {
    stop("comparison does not contain attraction terms for the focal vertex")
  }
  term_idx <- which.max(term_rows$step_max_abs)
  data.frame(
    case_id = target$case_id[[1L]],
    family = target$family[[1L]],
    layout_dim = target$layout_dim[[1L]],
    first_divergent_round = failing_round,
    focal_vertex = focal_vertex,
    focal_order_index = focal_order,
    term_index = term_rows$term_index[[term_idx]],
    neighbor = term_rows$neighbor[[term_idx]],
    term_step_max_abs = term_rows$step_max_abs[[term_idx]],
    term_delta_max_abs = term_rows$delta_max_abs[[term_idx]],
    term_norm2_delta = term_rows$norm2_delta[[term_idx]],
    term_scale_delta = term_rows$scale_delta[[term_idx]],
    stringsAsFactors = FALSE
  )
}

grip_weighted_nd_gflow_refinement_source_vertex_trace <- function(comparison,
                                                                  vertices) {
  rows <- comparison[comparison$vertex %in% vertices, , drop = FALSE]
  rows <- rows[order(rows$round_in_level, rows$vertex), , drop = FALSE]
  role <- ifelse(rows$vertex == vertices[[1L]], "focal", "neighbor")
  data.frame(
    round_in_level = rows$round_in_level,
    vertex = rows$vertex,
    role = role,
    order_index = rows$order_index,
    coords_before_max_abs = rows$coords_before_max_abs,
    attraction_disp_max_abs = rows$attraction_disp_max_abs,
    repulsion_disp_max_abs = rows$repulsion_disp_max_abs,
    pre_temp_disp_max_abs = rows$pre_temp_disp_max_abs,
    applied_disp_max_abs = rows$applied_disp_max_abs,
    coords_after_max_abs = rows$coords_after_max_abs,
    heat_before_delta = rows$heat_before_delta,
    heat_after_delta = rows$heat_after_delta,
    stringsAsFactors = FALSE
  )
}

grip_weighted_nd_gflow_refinement_source_term_trace <- function(comparison,
                                                                focal_vertex,
                                                                neighbor) {
  terms <- attr(comparison, "attraction_terms")
  rows <- terms[
    terms$vertex == focal_vertex & terms$neighbor == neighbor,
    ,
    drop = FALSE
  ]
  rows <- rows[order(rows$round_in_level, rows$term_index), , drop = FALSE]
  data.frame(
    round_in_level = rows$round_in_level,
    term_index = rows$term_index,
    vertex = rows$vertex,
    neighbor = rows$neighbor,
    weight_delta = rows$weight_delta,
    norm2_delta = rows$norm2_delta,
    desired_delta = rows$desired_delta,
    desired2_delta = rows$desired2_delta,
    scale_delta = rows$scale_delta,
    delta_max_abs = rows$delta_max_abs,
    step_max_abs = rows$step_max_abs,
    cumulative_max_abs = rows$cumulative_max_abs,
    stringsAsFactors = FALSE
  )
}

grip_weighted_nd_gflow_refinement_source_summary <- function(source,
                                                             tolerance = 1e-12) {
  vertex_trace <- source$vertex_trace
  term_trace <- source$term_trace
  first_vertex <- vertex_trace[
    vertex_trace$coords_before_max_abs > tolerance |
      vertex_trace$applied_disp_max_abs > tolerance |
      vertex_trace$attraction_disp_max_abs > tolerance |
      vertex_trace$repulsion_disp_max_abs > tolerance,
    ,
    drop = FALSE
  ]
  first_term <- term_trace[
    term_trace$delta_max_abs > tolerance |
      abs(term_trace$norm2_delta) > tolerance |
      abs(term_trace$scale_delta) > tolerance |
      term_trace$step_max_abs > tolerance,
    ,
    drop = FALSE
  ]

  data.frame(
    diagnostic = c(
      "first_vertex_state_difference",
      "first_focal_edge_term_difference",
      "first_focal_edge_step_difference"
    ),
    round_in_level = c(
      if (nrow(first_vertex)) first_vertex$round_in_level[[1L]] else NA_integer_,
      if (nrow(first_term)) first_term$round_in_level[[1L]] else NA_integer_,
      if (nrow(first_term[first_term$step_max_abs > tolerance, , drop = FALSE])) {
        first_term[first_term$step_max_abs > tolerance, , drop = FALSE]$round_in_level[[1L]]
      } else {
        NA_integer_
      }
    ),
    vertex = c(
      if (nrow(first_vertex)) first_vertex$vertex[[1L]] else NA_integer_,
      if (nrow(first_term)) first_term$vertex[[1L]] else NA_integer_,
      if (nrow(first_term[first_term$step_max_abs > tolerance, , drop = FALSE])) {
        first_term[first_term$step_max_abs > tolerance, , drop = FALSE]$vertex[[1L]]
      } else {
        NA_integer_
      }
    ),
    neighbor = c(
      NA_integer_,
      if (nrow(first_term)) first_term$neighbor[[1L]] else NA_integer_,
      if (nrow(first_term[first_term$step_max_abs > tolerance, , drop = FALSE])) {
        first_term[first_term$step_max_abs > tolerance, , drop = FALSE]$neighbor[[1L]]
      } else {
        NA_integer_
      }
    ),
    max_abs = c(
      if (nrow(first_vertex)) {
        max(abs(unlist(first_vertex[1L, c(
          "coords_before_max_abs",
          "applied_disp_max_abs",
          "attraction_disp_max_abs",
          "repulsion_disp_max_abs"
        )])), na.rm = TRUE)
      } else {
        NA_real_
      },
      if (nrow(first_term)) {
        max(abs(unlist(first_term[1L, c(
          "delta_max_abs",
          "norm2_delta",
          "scale_delta",
          "step_max_abs"
        )])), na.rm = TRUE)
      } else {
        NA_real_
      },
      if (nrow(first_term[first_term$step_max_abs > tolerance, , drop = FALSE])) {
        first_term[first_term$step_max_abs > tolerance, , drop = FALSE]$step_max_abs[[1L]]
      } else {
        NA_real_
      }
    ),
    stringsAsFactors = FALSE
  )
}

grip_weighted_nd_gflow_refinement_source_run <- function(
    mode = c("smoke", "full"),
    gflow_repo = grip_weighted_nd_gflow_parity_default_repo(),
    tolerance = 1e-8,
    source_tolerance = 1e-12) {
  mode <- match.arg(mode)
  target <- grip_weighted_nd_gflow_refinement_step_trace_target(
    mode = mode,
    gflow_repo = gflow_repo,
    tolerance = tolerance
  )
  target$round_start <- 1L
  target$round_end <- target$first_divergent_round
  comparison <- grip_weighted_nd_gflow_refinement_step_trace_run(
    mode = mode,
    gflow_repo = gflow_repo,
    target = target,
    tolerance = tolerance
  )
  source_target <- grip_weighted_nd_gflow_refinement_source_target(comparison)
  vertices <- c(source_target$focal_vertex[[1L]], source_target$neighbor[[1L]])
  vertex_trace <- grip_weighted_nd_gflow_refinement_source_vertex_trace(
    comparison,
    vertices = vertices
  )
  term_trace <- grip_weighted_nd_gflow_refinement_source_term_trace(
    comparison,
    focal_vertex = source_target$focal_vertex[[1L]],
    neighbor = source_target$neighbor[[1L]]
  )
  out <- list(
    target = source_target,
    vertex_trace = vertex_trace,
    term_trace = term_trace
  )
  out$summary <- grip_weighted_nd_gflow_refinement_source_summary(
    out,
    tolerance = source_tolerance
  )
  out
}
