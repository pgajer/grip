`%||%` <- function(x, y) if (is.null(x)) y else x

gisomap.dimension.resolve.args <- function(args, initializer, dim) {
  if (is.null(args)) {
    return(list())
  }
  if (!is.list(args)) {
    stop("argument collections must be lists")
  }
  specific <- args[[initializer]]
  out <- if (is.list(specific)) specific else args
  out <- out[setdiff(names(out), c("metric_mds", "weighted_grip"))]
  lapply(out, function(value) {
    if (is.function(value)) {
      value(dim = dim, initializer = initializer)
    } else {
      value
    }
  })
}

gisomap.dimension.scalar <- function(x, name, default = NA_real_) {
  if (is.null(x) || !is.data.frame(x) || !name %in% names(x) || nrow(x) < 1L) {
    return(default)
  }
  value <- x[[name]][[nrow(x)]]
  if (length(value) == 0L) default else value
}

gisomap.dimension.trace.scalar <- function(trace, name, default = NA_real_) {
  if (is.null(trace) || !is.data.frame(trace) || !name %in% names(trace) || nrow(trace) < 1L) {
    return(default)
  }
  value <- trace[[name]][[nrow(trace)]]
  if (length(value) == 0L) default else value
}

gisomap.dimension.row <- function(graph_name,
                                  initializer,
                                  dim,
                                  fit,
                                  elapsed_sec,
                                  status = "ok",
                                  error = NA_character_) {
  diag <- fit$diagnostics
  initial <- fit$metadata$initial_diagnostics
  trace <- fit$trace
  edge_rel_rmse <- gisomap.dimension.scalar(diag, "edge.rel.rmse")
  initial_edge_rel_rmse <- gisomap.dimension.scalar(initial, "edge.rel.rmse")
  data.frame(
    graph = graph_name,
    initializer = initializer,
    dim = as.integer(dim),
    status = status,
    elapsed_sec = elapsed_sec,
    edge_rel_rmse = edge_rel_rmse,
    edge_rmse = gisomap.dimension.scalar(diag, "edge.rmse"),
    edge_scale = gisomap.dimension.scalar(diag, "edge.scale"),
    gmds_stress = gisomap.dimension.scalar(diag, "gmds.stress"),
    gmds_short_stress = gisomap.dimension.scalar(diag, "gmds.short.stress"),
    gmds_mid_stress = gisomap.dimension.scalar(diag, "gmds.mid.stress"),
    gmds_long_stress = gisomap.dimension.scalar(diag, "gmds.long.stress"),
    initial_edge_rel_rmse = initial_edge_rel_rmse,
    initial_edge_rmse = gisomap.dimension.scalar(initial, "edge.rmse"),
    edge_rel_rmse_delta = edge_rel_rmse - initial_edge_rel_rmse,
    final_energy = gisomap.dimension.trace.scalar(trace, "energy"),
    final_gradient_norm = gisomap.dimension.trace.scalar(trace, "gradient_norm"),
    final_step = gisomap.dimension.trace.scalar(trace, "step"),
    accepted_steps = if (is.data.frame(trace) && "accepted" %in% names(trace)) {
      sum(as.logical(trace$accepted), na.rm = TRUE)
    } else {
      NA_integer_
    },
    trace_rows = if (is.data.frame(trace)) nrow(trace) else NA_integer_,
    error = error,
    stringsAsFactors = FALSE
  )
}

gisomap.dimension.add.metadata <- function(row, metadata = list()) {
  if (is.null(metadata) || length(metadata) == 0L) {
    return(row)
  }
  for (name in names(metadata)) {
    value <- metadata[[name]]
    if (length(value) != 1L) {
      value <- paste(value, collapse = ", ")
    }
    row[[name]] <- value
  }
  if ("expected_dim" %in% names(row)) {
    row$dim_gap <- row$dim - row$expected_dim
    row$dim_relation <- ifelse(
      row$dim < row$expected_dim,
      "below",
      ifelse(row$dim == row$expected_dim, "expected", "above")
    )
  }
  row
}

gisomap.dimension.error.row <- function(graph_name,
                                        initializer,
                                        dim,
                                        elapsed_sec,
                                        error) {
  data.frame(
    graph = graph_name,
    initializer = initializer,
    dim = as.integer(dim),
    status = "error",
    elapsed_sec = elapsed_sec,
    edge_rel_rmse = NA_real_,
    edge_rmse = NA_real_,
    edge_scale = NA_real_,
    gmds_stress = NA_real_,
    gmds_short_stress = NA_real_,
    gmds_mid_stress = NA_real_,
    gmds_long_stress = NA_real_,
    initial_edge_rel_rmse = NA_real_,
    initial_edge_rmse = NA_real_,
    edge_rel_rmse_delta = NA_real_,
    final_energy = NA_real_,
    final_gradient_norm = NA_real_,
    final_step = NA_real_,
    accepted_steps = NA_integer_,
    trace_rows = NA_integer_,
    error = conditionMessage(error),
    stringsAsFactors = FALSE
  )
}

gisomap.dimension.sweep <- function(graph_name = "graph",
                                    prepared = NULL,
                                    edges = NULL,
                                    n = NULL,
                                    adj_list = NULL,
                                    weight_list = NULL,
                                    edge_weights = NULL,
                                    dims = 2:6,
                                    initializers = c("metric_mds", "weighted_grip"),
                                    init_args = list(),
                                    edge_kk_args = list(),
                                    seed = 1L,
                                    diagnostics = TRUE,
                                    keep_fits = TRUE,
                                    metadata = list(),
                                    verbose = TRUE) {
  initializers <- match.arg(initializers, c("metric_mds", "weighted_grip"), several.ok = TRUE)
  dims <- as.integer(dims)
  if (length(dims) < 1L || any(is.na(dims)) || any(dims < 2L)) {
    stop("dims must contain integer dimensions >= 2")
  }
  prepared <- grip:::grip.gmds.require.prepared(
    prepared = prepared,
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights
  )

  rows <- list()
  fits <- list()
  idx <- 1L
  for (initializer in initializers) {
    for (dim in dims) {
      if (isTRUE(verbose)) {
        message(sprintf("[%s] init=%s dim=%d", graph_name, initializer, dim))
      }
      start <- proc.time()[["elapsed"]]
      result <- tryCatch({
        fit <- grip.layout.gisomap(
          prepared = prepared,
          dim = dim,
          init = initializer,
          init_args = gisomap.dimension.resolve.args(init_args, initializer, dim),
          edge_kk_args = gisomap.dimension.resolve.args(edge_kk_args, initializer, dim),
          diagnostics = diagnostics,
          seed = seed + dim
        )
        elapsed <- proc.time()[["elapsed"]] - start
        rows[[idx]] <- gisomap.dimension.row(
          graph_name = graph_name,
          initializer = initializer,
          dim = dim,
          fit = fit,
          elapsed_sec = elapsed
        )
        rows[[idx]] <- gisomap.dimension.add.metadata(rows[[idx]], metadata)
        if (isTRUE(keep_fits)) {
          fits[[paste(graph_name, initializer, dim, sep = "__")]] <- fit
        }
        TRUE
      }, error = function(e) {
        elapsed <- proc.time()[["elapsed"]] - start
        rows[[idx]] <<- gisomap.dimension.error.row(
          graph_name = graph_name,
          initializer = initializer,
          dim = dim,
          elapsed_sec = elapsed,
          error = e
        )
        rows[[idx]] <<- gisomap.dimension.add.metadata(rows[[idx]], metadata)
        FALSE
      })
      idx <- idx + 1L
    }
  }
  out <- list(
    graph_name = graph_name,
    dims = dims,
    initializers = initializers,
    metadata = metadata,
    results = do.call(rbind, rows),
    fits = fits
  )
  class(out) <- c("gisomap_dimension_sweep", "list")
  out
}

gisomap.dimension.example.graphs <- function() {
  mesh_edges <- edges.mesh(4L, 4L)
  cube_edges <- edges.cube(3L)
  saddle <- mesh.surface.graph(5L, 5L, surface = "saddle", amplitude = 0.6)
  list(
    grid_4x4 = list(
      name = "grid_4x4",
      edges = mesh_edges,
      n = 16L,
      edge_weights = rep(1, nrow(mesh_edges)),
      description = "Uniform 4 x 4 orthogonal grid with unit edge lengths."
    ),
    cube_lattice_3 = list(
      name = "cube_lattice_3",
      edges = cube_edges,
      n = max(cube_edges),
      edge_weights = rep(1, nrow(cube_edges)),
      description = "Small cubic fixture from edges.cube(3) with unit edge lengths."
    ),
    saddle_5x5 = list(
      name = "saddle_5x5",
      edges = saddle$edges,
      n = saddle$n,
      edge_weights = saddle$edge_weights,
      description = "Weighted 5 x 5 saddle surface graph from grip fixtures."
    )
  )
}

gisomap.dimension.sweep.examples <- function(dims = 2:5,
                                             initializers = c("metric_mds", "weighted_grip"),
                                             edge_kk_max_iter = 25L,
                                             weighted_rounds = 8L,
                                             seed = 101L,
                                             verbose = TRUE) {
  graphs <- gisomap.dimension.example.graphs()
  init_args <- list(
    weighted_grip = list(
      placement = "barycenter",
      rounds = weighted_rounds,
      final_rounds = weighted_rounds,
      num_init = function(dim, initializer) max(6L, dim + 1L),
      num_nbrs = 8L,
      repulsion_factor = 0.5
    )
  )
  edge_kk_args <- list(
    stiffness_method = "uniform",
    density_mix_schedule = 1,
    scale_mode = "profiled",
    max_iter = edge_kk_max_iter,
    return_trace = TRUE,
    engine = "cpp"
  )

  lapply(graphs, function(graph) {
    gisomap.dimension.sweep(
      graph_name = graph$name,
      edges = graph$edges,
      n = graph$n,
      edge_weights = graph$edge_weights,
      dims = dims,
      initializers = initializers,
      init_args = init_args,
      edge_kk_args = edge_kk_args,
      seed = seed,
      verbose = verbose
    )
  })
}

gisomap.dimension.bind.results <- function(sweeps) {
  do.call(rbind, lapply(sweeps, `[[`, "results"))
}

gisomap.dimension.load.gflow <- function(
    gflow_repo = "/Users/pgajer/current_projects/gflow") {
  if (dir.exists(gflow_repo) && file.exists(file.path(gflow_repo, "DESCRIPTION")) &&
      requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(gflow_repo, quiet = TRUE)
    return(asNamespace("gflow"))
  }
  if (!requireNamespace("gflow", quietly = TRUE)) {
    stop("gflow is not installed and local pkgload::load_all() is unavailable")
  }
  asNamespace("gflow")
}

gisomap.dimension.gflow_graph_to_fixture <- function(graph,
                                                     name,
                                                     expected_dim,
                                                     family,
                                                     description,
                                                     edge_weights = graph$edge_weight,
                                                     graph_construction = "gflow::create.sknn.graph",
                                                     source_dim = NA_integer_,
                                                     k = NA_integer_) {
  edges <- graph$edge_matrix
  if (is.null(edges)) {
    stop("gflow graph must contain edge_matrix")
  }
  edge_weights <- as.numeric(edge_weights)
  if (length(edge_weights) != nrow(edges)) {
    stop("edge_weights must be parallel to graph$edge_matrix")
  }
  list(
    name = name,
    edges = matrix(as.integer(edges), ncol = 2L),
    n = if (!is.null(graph$n_vertices)) graph$n_vertices else max(edges),
    edge_weights = edge_weights,
    dims = 2:(expected_dim + 3L),
    metadata = list(
      family = family,
      expected_dim = as.integer(expected_dim),
      source_dim = as.integer(source_dim),
      graph_construction = graph_construction,
      k = as.integer(k),
      n_vertices = if (!is.null(graph$n_vertices)) graph$n_vertices else max(edges),
      n_edges = nrow(edges),
      n_components_before = graph$n_components_before %||% NA_integer_,
      n_components_after = graph$n_components_after %||% NA_integer_,
      n_mst_edges_added = graph$n_mst_edges_added %||% NA_integer_,
      description = description
    )
  )
}

gisomap.dimension.quadform_fixture <- function(dim,
                                               n = 70L,
                                               k = 8L,
                                               index.k = 1L,
                                               coefficients = NULL,
                                               seed = 1000L + dim,
                                               graph_space = c("embed", "param"),
                                               gflow_repo = "/Users/pgajer/current_projects/gflow") {
  graph_space <- match.arg(graph_space)
  gflow <- gisomap.dimension.load.gflow(gflow_repo)
  sample_fn <- get("quadform.sample.dataset", envir = gflow)
  sknn_fn <- get("create.sknn.graph", envir = gflow)
  edge_lengths_fn <- get("quadform.edge.lengths", envir = gflow)
  ds <- sample_fn(
    n = n,
    dim = dim,
    index.k = index.k,
    coefficients = coefficients,
    grid.size = if (dim == 2L) 41L else 7L,
    sample.connection.k = k,
    seed = seed
  )
  graph_X <- if (identical(graph_space, "embed")) ds$X_embed else ds$X_param
  graph <- sknn_fn(
    graph_X,
    k = k,
    connect.components = TRUE,
    connect.method = "component.mst",
    edge.weight = "distance",
    neighbor.method = "exact"
  )
  edges <- graph$edge_matrix
  edge_weights <- edge_lengths_fn(
    ds$X_param[edges[, 1L], , drop = FALSE],
    ds$X_param[edges[, 2L], , drop = FALSE],
    index.k = index.k,
    coefficients = coefficients
  )
  coeff_tag <- if (is.null(coefficients)) "default" else paste(coefficients, collapse = "_")
  gisomap.dimension.gflow_graph_to_fixture(
    graph = graph,
    name = sprintf("quadform_%dd_index%d_k%d", dim, index.k, k),
    expected_dim = dim,
    family = "quadform",
    description = sprintf(
      "%dD gflow quadform sample, index %d, coefficients %s; sKNN support with component-MST repair and quadform edge lengths.",
      dim, index.k, coeff_tag
    ),
    edge_weights = edge_weights,
    graph_construction = sprintf("gflow::quadform.sample.dataset(dim=%d) + gflow::create.sknn.graph(k=%d, connect.components=TRUE)", dim, k),
    source_dim = dim,
    k = k
  )
}

gisomap.dimension.euclidean_fixture <- function(expected_dim,
                                                n = 72L,
                                                k = max(8L, 2L * expected_dim + 6L),
                                                seed = 2000L + expected_dim,
                                                gflow_repo = "/Users/pgajer/current_projects/gflow") {
  gflow <- gisomap.dimension.load.gflow(gflow_repo)
  sknn_fn <- get("create.sknn.graph", envir = gflow)
  set.seed(seed)
  X <- matrix(stats::runif(n * expected_dim, min = -1, max = 1), ncol = expected_dim)
  X <- scale(X, center = TRUE, scale = FALSE)
  graph <- sknn_fn(
    X,
    k = k,
    connect.components = TRUE,
    connect.method = "component.mst",
    edge.weight = "distance",
    neighbor.method = "exact"
  )
  fixture <- gisomap.dimension.gflow_graph_to_fixture(
    graph = graph,
    name = sprintf("euclidean_%dd_sknn_k%d", expected_dim, k),
    expected_dim = expected_dim,
    family = "euclidean_random_geometric",
    description = sprintf(
      "Uniform random geometric graph in R^%d; sKNN support with component-MST repair and true Euclidean edge lengths.",
      expected_dim
    ),
    graph_construction = sprintf("gflow::create.sknn.graph(k=%d, connect.components=TRUE)", k),
    source_dim = expected_dim,
    k = k
  )
  fixture$dims <- 2:(expected_dim + 3L)
  fixture
}

gisomap.dimension.known.fixtures <- function(
    gflow_repo = "/Users/pgajer/current_projects/gflow") {
  list(
    quadform_2d = gisomap.dimension.quadform_fixture(
      dim = 2L,
      n = 70L,
      k = 8L,
      index.k = 1L,
      coefficients = c(1.0, 0.45),
      seed = 3102L,
      gflow_repo = gflow_repo
    ),
    quadform_3d = gisomap.dimension.quadform_fixture(
      dim = 3L,
      n = 70L,
      k = 10L,
      index.k = 1L,
      coefficients = c(1.0, 0.45, 0.2),
      seed = 3103L,
      gflow_repo = gflow_repo
    ),
    euclidean_4d = gisomap.dimension.euclidean_fixture(
      expected_dim = 4L,
      n = 72L,
      seed = 3104L,
      gflow_repo = gflow_repo
    ),
    euclidean_5d = gisomap.dimension.euclidean_fixture(
      expected_dim = 5L,
      n = 72L,
      seed = 3105L,
      gflow_repo = gflow_repo
    ),
    euclidean_6d = gisomap.dimension.euclidean_fixture(
      expected_dim = 6L,
      n = 72L,
      seed = 3106L,
      gflow_repo = gflow_repo
    )
  )
}

gisomap.dimension.expected_summary <- function(results,
                                               threshold = 1e-3,
                                               drop_factor = 10) {
  ok <- results[results$status == "ok" & !is.na(results$expected_dim), , drop = FALSE]
  if (nrow(ok) == 0L) {
    return(data.frame())
  }
  keys <- unique(ok[, c("graph", "initializer"), drop = FALSE])
  rows <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    d <- ok[ok$graph == keys$graph[[i]] & ok$initializer == keys$initializer[[i]], , drop = FALSE]
    expected_dim <- unique(d$expected_dim)[[1L]]
    below <- d[d$dim < expected_dim, , drop = FALSE]
    at_or_above <- d[d$dim >= expected_dim, , drop = FALSE]
    at_expected <- d[d$dim == expected_dim, , drop = FALSE]
    best <- d[which.min(d$edge_rel_rmse), , drop = FALSE]
    first_below_threshold <- at_or_above[at_or_above$edge_rel_rmse <= threshold, , drop = FALSE]
    rows[[i]] <- data.frame(
      graph = keys$graph[[i]],
      initializer = keys$initializer[[i]],
      expected_dim = expected_dim,
      best_dim = best$dim[[1L]],
      best_edge_rel_rmse = best$edge_rel_rmse[[1L]],
      expected_edge_rel_rmse = if (nrow(at_expected)) at_expected$edge_rel_rmse[[1L]] else NA_real_,
      best_below_edge_rel_rmse = if (nrow(below)) min(below$edge_rel_rmse, na.rm = TRUE) else NA_real_,
      below_to_expected_ratio = if (nrow(below) && nrow(at_expected)) {
        min(below$edge_rel_rmse, na.rm = TRUE) / at_expected$edge_rel_rmse[[1L]]
      } else {
        NA_real_
      },
      first_dim_below_threshold = if (nrow(first_below_threshold)) first_below_threshold$dim[[1L]] else NA_integer_,
      expected_pass_threshold = if (nrow(at_expected)) at_expected$edge_rel_rmse[[1L]] <= threshold else NA,
      expected_pass_drop = if (nrow(below) && nrow(at_expected)) {
        min(below$edge_rel_rmse, na.rm = TRUE) >= drop_factor * at_expected$edge_rel_rmse[[1L]]
      } else {
        NA
      },
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

gisomap.dimension.sweep.fixtures <- function(fixtures,
                                             initializers = c("metric_mds", "weighted_grip"),
                                             edge_kk_max_iter = 80L,
                                             weighted_rounds = 12L,
                                             seed = 401L,
                                             verbose = TRUE) {
  init_args <- list(
    weighted_grip = list(
      placement = "barycenter",
      rounds = weighted_rounds,
      final_rounds = weighted_rounds,
      num_init = function(dim, initializer) max(8L, dim + 1L),
      num_nbrs = 12L,
      repulsion_factor = 0.5
    )
  )
  edge_kk_args <- list(
    stiffness_method = "uniform",
    density_mix_schedule = 1,
    scale_mode = "profiled",
    max_iter = edge_kk_max_iter,
    return_trace = TRUE,
    engine = "cpp"
  )
  lapply(fixtures, function(fixture) {
    gisomap.dimension.sweep(
      graph_name = fixture$name,
      edges = fixture$edges,
      n = fixture$n,
      edge_weights = fixture$edge_weights,
      dims = fixture$dims,
      initializers = initializers,
      init_args = init_args,
      edge_kk_args = edge_kk_args,
      seed = seed,
      metadata = fixture$metadata,
      verbose = verbose
    )
  })
}

gisomap.dimension.html.escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

gisomap.dimension.format.value <- function(x, digits = 4L) {
  if (length(x) == 0L || is.na(x)) {
    return("")
  }
  if (is.numeric(x)) {
    return(format(signif(x, digits), scientific = TRUE, trim = TRUE))
  }
  gisomap.dimension.html.escape(x)
}

gisomap.dimension.is.integerish <- function(x) {
  finite <- is.finite(x)
  any(finite) && all(abs(x[finite] - round(x[finite])) < 1e-8) &&
    max(abs(x[finite])) < 1e6
}

gisomap.dimension.html.table <- function(x, columns = names(x), digits = 4L) {
  columns <- intersect(columns, names(x))
  x <- x[, columns, drop = FALSE]
  header <- paste0("<th>", gisomap.dimension.html.escape(names(x)), "</th>", collapse = "")
  rows <- apply(x, 1L, function(row) {
    cells <- vapply(seq_along(row), function(i) {
      value <- row[[i]]
      if (is.numeric(x[[i]])) {
        value <- as.numeric(value)
        value <- if (gisomap.dimension.is.integerish(x[[i]])) {
          if (is.na(value)) "" else format(as.integer(round(value)), trim = TRUE)
        } else {
          gisomap.dimension.format.value(value, digits = digits)
        }
      } else {
        value <- if (is.na(value)) "" else gisomap.dimension.html.escape(value)
      }
      paste0("<td>", value, "</td>")
    }, character(1L))
    paste0("<tr>", paste0(cells, collapse = ""), "</tr>")
  })
  paste0(
    "<table><thead><tr>", header, "</tr></thead><tbody>",
    paste0(rows, collapse = "\n"),
    "</tbody></table>"
  )
}

gisomap.dimension.uniquify.svg <- function(svg, prefix) {
  prefix <- gsub("[^A-Za-z0-9_\\-]", "_", prefix)
  ids <- regmatches(svg, gregexpr("id=\"[^\"]+\"", svg, perl = TRUE))[[1L]]
  if (!length(ids) || identical(ids, character(0L)) || ids[[1L]] == "-1") {
    return(svg)
  }
  ids <- unique(sub("^id=\"([^\"]+)\"$", "\\1", ids))
  ids <- ids[order(nchar(ids), decreasing = TRUE)]
  for (id in ids) {
    new_id <- paste0(prefix, "__", id)
    svg <- gsub(paste0("id=\"", id, "\""), paste0("id=\"", new_id, "\""), svg, fixed = TRUE)
    svg <- gsub(paste0("href=\"#", id, "\""), paste0("href=\"#", new_id, "\""), svg, fixed = TRUE)
    svg <- gsub(paste0("xlink:href=\"#", id, "\""), paste0("xlink:href=\"#", new_id, "\""), svg, fixed = TRUE)
    svg <- gsub(paste0("url(#", id, ")"), paste0("url(#", new_id, ")"), svg, fixed = TRUE)
  }
  svg
}

gisomap.dimension.curve.svg <- function(results, graph_name) {
  data <- results[results$graph == graph_name & results$status == "ok", , drop = FALSE]
  if (nrow(data) == 0L) {
    return("<p>No successful fits for this graph.</p>")
  }
  tmp <- tempfile(fileext = ".svg")
  grDevices::svg(tmp, width = 9.6, height = 5.4)
  old <- par(mar = c(5.8, 5.0, 2.6, 1.2), las = 1)
  tryCatch({
    y <- data$edge_rel_rmse
    y <- y[is.finite(y) & y > 0]
    ylim <- if (length(y) > 0L) range(y) else c(1e-6, 1)
    if (diff(ylim) == 0) {
      ylim <- ylim * c(0.8, 1.2)
    }
    plot(
      NA,
      xlim = range(data$dim),
      ylim = ylim,
      log = "y",
      xlab = "Embedding dimension",
      ylab = "Final edge relative RMSE",
      main = graph_name,
      cex.lab = 1.05,
      cex.axis = 0.95
    )
    if ("expected_dim" %in% names(data) && any(!is.na(data$expected_dim))) {
      expected_dim <- unique(data$expected_dim[!is.na(data$expected_dim)])[[1L]]
      graphics::abline(v = expected_dim, col = "#4a5568", lty = 2, lwd = 1.4)
      graphics::mtext(
        sprintf("expected d = %d", expected_dim),
        side = 3,
        line = 0.1,
        adj = 0.02,
        col = "#4a5568",
        cex = 0.85
      )
    }
    colors <- c(metric_mds = "#2b6cb0", weighted_grip = "#b83280")
    pchs <- c(metric_mds = 16, weighted_grip = 17)
    for (initializer in unique(data$initializer)) {
      d <- data[data$initializer == initializer, , drop = FALSE]
      d <- d[order(d$dim), , drop = FALSE]
      lines(d$dim, d$edge_rel_rmse, col = colors[[initializer]], lwd = 2)
      points(d$dim, d$edge_rel_rmse, col = colors[[initializer]], pch = pchs[[initializer]], cex = 1.1)
    }
    legend(
      "bottom",
      legend = unique(data$initializer),
      col = colors[unique(data$initializer)],
      pch = pchs[unique(data$initializer)],
      lwd = 2,
      horiz = TRUE,
      bty = "n"
    )
  }, finally = {
    par(old)
    grDevices::dev.off()
  })
  svg <- paste(readLines(tmp, warn = FALSE), collapse = "\n")
  unlink(tmp)
  gisomap.dimension.uniquify.svg(svg, paste0("gisomap_", graph_name))
}

gisomap.dimension.graph.interpretation <- function(graph_results) {
  ok <- graph_results[graph_results$status == "ok", , drop = FALSE]
  if (nrow(ok) == 0L) {
    return("No successful fits were available for interpretation.")
  }
  expected <- if ("expected_dim" %in% names(ok) && any(!is.na(ok$expected_dim))) {
    unique(ok$expected_dim[!is.na(ok$expected_dim)])[[1L]]
  } else {
    NA_integer_
  }
  lines <- vapply(unique(ok$initializer), function(initializer) {
    d <- ok[ok$initializer == initializer, , drop = FALSE]
    d <- d[order(d$dim), , drop = FALSE]
    best <- d[which.min(d$edge_rel_rmse), , drop = FALSE]
    if (!is.na(expected) && any(d$dim == expected)) {
      at_expected <- d[d$dim == expected, , drop = FALSE]
      below <- d[d$dim < expected, , drop = FALSE]
      drop_text <- if (nrow(below) > 0L) {
        ratio <- min(below$edge_rel_rmse, na.rm = TRUE) / at_expected$edge_rel_rmse[[1L]]
        sprintf(
          "At the expected dimension the residual is %.3g, a %.1fx drop from the best lower-dimensional residual.",
          at_expected$edge_rel_rmse[[1L]],
          ratio
        )
      } else {
        sprintf(
          "At the expected dimension the residual is %.3g.",
          at_expected$edge_rel_rmse[[1L]]
        )
      }
      sprintf(
        "%s: best residual %.3g occurs at dimension %d. %s",
        initializer,
        best$edge_rel_rmse[[1L]],
        best$dim[[1L]],
        drop_text
      )
    } else {
      sprintf(
        "%s: best residual %.3g occurs at dimension %d.",
        initializer,
        best$edge_rel_rmse[[1L]],
        best$dim[[1L]]
      )
    }
  }, character(1L))
  paste(lines, collapse = " ")
}

gisomap.dimension.what.learned <- function(results, expected_summary) {
  if (nrow(expected_summary) == 0L) {
    return("This exploratory run verifies that the sweep machinery executes, but it does not contain expected-dimension fixtures.")
  }
  metric <- expected_summary[expected_summary$initializer == "metric_mds", , drop = FALSE]
  weighted <- expected_summary[expected_summary$initializer == "weighted_grip", , drop = FALSE]
  metric_exact <- sum(metric$best_dim == metric$expected_dim, na.rm = TRUE)
  weighted_exact <- sum(weighted$best_dim == weighted$expected_dim, na.rm = TRUE)
  metric_drop <- sum(metric$expected_pass_drop %in% TRUE, na.rm = TRUE)
  weighted_drop <- sum(weighted$expected_pass_drop %in% TRUE, na.rm = TRUE)
  paste0(
    "The metric-MDS initializer gives the cleaner dimension diagnostic in this run: ",
    metric_exact, " of ", nrow(metric), " metric-MDS curves have their best residual at the expected dimension, and ",
    metric_drop, " pass the lower-to-expected drop heuristic. ",
    "The weighted-GRIP initializer is less stable under the current quick tuning: ",
    weighted_exact, " of ", nrow(weighted), " curves have their best residual at the expected dimension, and ",
    weighted_drop, " pass the drop heuristic. ",
    "The Euclidean high-dimensional fixtures support the proposed diagnostic most strongly; the quadform fixtures show the expected lower-dimensional drop, but still improve above the nominal dimension, which suggests that graph construction and edge-only flexibility remain important failure modes."
  )
}

gisomap.dimension.report <- function(sweeps,
                                     file,
                                     title = "GISOMAP Dimension Sweep Report",
                                     subtitle = NULL,
                                     threshold = 1e-3,
                                     drop_factor = 10) {
  results <- gisomap.dimension.bind.results(sweeps)
  expected_summary <- gisomap.dimension.expected_summary(
    results,
    threshold = threshold,
    drop_factor = drop_factor
  )
  graphs <- unique(results$graph)
  run_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  git_ref <- tryCatch(
    system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_
  )
  run_meta <- data.frame(
    item = c("Build time", "Git commit", "Report file", "Rows", "Errors"),
    value = c(run_time, git_ref, file, nrow(results), sum(results$status != "ok")),
    stringsAsFactors = FALSE
  )
  graph_sections <- vapply(graphs, function(graph) {
    graph_results <- results[results$graph == graph, , drop = FALSE]
    description <- if ("description" %in% names(graph_results)) {
      unique(graph_results$description[!is.na(graph_results$description)])
    } else {
      character()
    }
    description <- if (length(description)) description[[1L]] else "Fixture metadata were not provided."
    paste0(
      "<section><h2>", gisomap.dimension.html.escape(graph), "</h2>\n",
      "<p>", gisomap.dimension.html.escape(description), "</p>\n",
      "<p>This figure plots final edge-KK relative RMSE against target dimension. ",
      "A dashed vertical line marks the expected dimension when the fixture provides one.</p>\n",
      gisomap.dimension.curve.svg(results, graph), "\n",
      "<p class=\"interpretation\">", gisomap.dimension.html.escape(gisomap.dimension.graph.interpretation(graph_results)), "</p>\n",
      "</section>"
    )
  }, character(1L))
  summary_block <- if (nrow(expected_summary) > 0L) {
    paste0(
      "<h2>Expected-Dimension Summary</h2>",
      "<p>This curated table is the main pass/fail diagnostic. The residual at the expected dimension should be small, ",
      "and the best lower-dimensional residual should be substantially larger.</p>",
      gisomap.dimension.html.table(
        expected_summary,
        columns = c(
          "graph", "initializer", "expected_dim", "best_dim",
          "expected_edge_rel_rmse", "best_below_edge_rel_rmse",
          "below_to_expected_ratio", "expected_pass_threshold",
          "expected_pass_drop"
        )
      )
    )
  } else {
    ""
  }
  learned_detail <- if (nrow(expected_summary) > 0L) {
    paste0(
      "The current evidence is mixed. The high-dimensional Euclidean fixtures support the basic dimension-sweep hypothesis, especially with metric MDS initialization. ",
      "The quadform fixtures and the weighted-GRIP initializer show that graph construction, curvature, local minima, and edge-only flexibility can shift the best residual above the nominal dimension. ",
      "The next methodological step is to tune graph density and polish budgets, then add held-out or non-edge geodesic diagnostics so the estimator is not driven by local edge flexibility alone."
    )
  } else {
    paste0(
      "This exploratory report is a smoke test for the GISOMAP sweep machinery. ",
      "Because these fixtures do not all carry expected intrinsic dimensions, the report should be read as an execution and diagnostic-format check rather than a calibration result. ",
      "Known-dimension fixtures should be used before drawing conclusions about dimension recovery."
    )
  }
  html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\"><title>",
    gisomap.dimension.html.escape(title),
    "</title><style>",
    "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:0;line-height:1.55;color:#1f2933;background:#fbfbf8;}",
    "main{max-width:1120px;margin:0 auto;padding:36px 42px 56px;}h1{font-size:32px;margin:0 0 6px;}h2{font-size:23px;margin-top:36px;border-top:1px solid #d9d2c3;padding-top:24px;}h3{font-size:17px;margin-top:22px;}",
    ".meta{color:#5f6b73;margin-bottom:22px;}.note{max-width:980px;background:#f1f5f4;border-left:4px solid #2b6cb0;padding:12px 16px;margin:18px 0;}",
    ".formula{font-family:Georgia,serif;background:#fff;border:1px solid #d9dee2;padding:10px 14px;margin:12px 0;display:inline-block;}.interpretation{max-width:980px;}",
    "table{border-collapse:collapse;width:100%;font-size:13px;margin:10px 0 22px;background:white;}th,td{border:1px solid #d9dee2;padding:6px 8px;text-align:right;vertical-align:top;}th{background:#edf2f7;color:#243746;}td:first-child,th:first-child,td:nth-child(2),th:nth-child(2),td:nth-child(4),th:nth-child(4),td:last-child,th:last-child{text-align:left;}",
    "svg{width:100%;height:auto;background:white;border:1px solid #d9dee2;margin:8px 0 12px;}code{background:#eef2f1;padding:1px 4px;border-radius:3px;}details{margin-top:18px;}",
    "</style></head><body>",
    "<main>",
    "<h1>", gisomap.dimension.html.escape(title), "</h1>",
    if (!is.null(subtitle)) paste0("<p class=\"meta\">", gisomap.dimension.html.escape(subtitle), "</p>") else "",
    "<h2>Purpose</h2>",
    "<p>The question is whether GISOMAP residual curves can diagnose the global embedding dimension of a weighted graph. ",
    "For a fixture with known intrinsic dimension <em>d</em>, a useful diagnostic should keep edge-KK residuals high below <em>d</em>, drop sharply at <em>d</em>, and remain approximately flat above <em>d</em>.</p>",
    "<h2>Diagnostic</h2>",
    "<p>Each run calls <code>grip.layout.gisomap()</code> with either metric MDS or weighted-GRIP initialization, followed by edge-KK polish. ",
    "The plotted quantity is the final relative edge RMSE from the existing edge-KK diagnostics.</p>",
    "<div class=\"formula\">E<sub>edge</sub>(Z) = 1/2 &sum;<sub>(i,j)&isin;E</sub> k<sub>ij</sub>(&#8741;z<sub>i</sub>-z<sub>j</sub>&#8741; - s w<sub>ij</sub>)<sup>2</sup></div>",
    summary_block,
    "<h2>Residual Curves</h2>",
    paste0(graph_sections, collapse = "\n"),
    "<h2>Results Summary And Discussion</h2>",
    "<p>", gisomap.dimension.html.escape(gisomap.dimension.what.learned(results, expected_summary)), "</p>",
    "<h2>What We Learned</h2>",
    if (nrow(expected_summary) > 0L) gisomap.dimension.html.table(expected_summary) else "",
    "<p>", gisomap.dimension.html.escape(learned_detail), "</p>",
    "<h2>Appendix</h2>",
    "<h3>Run Metadata</h3>", gisomap.dimension.html.table(run_meta),
    "<h3>Full Result Table</h3>",
    gisomap.dimension.html.table(
      results,
      columns = c(
        "graph", "family", "expected_dim", "initializer", "dim",
        "dim_relation", "status", "edge_rel_rmse",
        "edge_rmse", "gmds_stress", "initial_edge_rel_rmse",
        "edge_rel_rmse_delta", "elapsed_sec", "error"
      )
    ),
    "<h3>Reproducibility</h3>",
    "<p>Run the corresponding script in <code>dev/gisomap</code> from the grip repository root. ",
    "Large tables are also written as CSV artifacts in <code>output/gisomap</code>.</p>",
    "</main></body></html>"
  )
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, file)
  invisible(file)
}
