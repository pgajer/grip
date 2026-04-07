gripui.gmds.default.catalog <- function(catalog = gripui_graph_family_catalog()) {
  gripui.family.validate.catalog(catalog)
  catalog
}

gripui.gmds.initial.state <- function(catalog) {
  family_id <- names(catalog)[[1L]]
  desc <- catalog[[family_id]]
  category <- desc$category
  choices <- gripui.family.choices(catalog, category)
  list(
    category = category,
    choices = choices,
    family_id = family_id,
    desc = desc
  )
}

gripui.gmds.default.seed <- function(values, offset = 0L, fallback = 1L) {
  base <- fallback
  if ("seed" %in% names(values)) {
    seed <- suppressWarnings(as.integer(round(as.numeric(values$seed[[1L]]))))
    if (length(seed) == 1L && is.finite(seed)) {
      base <- seed
    }
  }
  as.integer(base + offset)
}

gripui.gmds.filter.edge.matrix <- function(edges, vertex_ids) {
  edges <- as.matrix(edges)
  storage.mode(edges) <- "integer"
  vertex_ids <- as.integer(vertex_ids)
  if (!length(vertex_ids) || !nrow(edges)) {
    return(matrix(integer(), ncol = 2L))
  }
  keep <- edges[, 1L] %in% vertex_ids & edges[, 2L] %in% vertex_ids
  edges[keep, , drop = FALSE]
}

gripui.gmds.layer_partition <- function(prepared) {
  levels <- prepared$misf$levels
  partition <- integer(prepared$n)
  for (level in seq.int(from = 0L, to = length(levels) - 1L, by = 1L)) {
    partition[as.integer(levels[[level + 1L]])] <- level
  }
  partition
}

gripui.gmds.layer.colors <- function(prepared) {
  max.level <- length(prepared$misf$levels) - 1L
  palette <- grDevices::hcl.colors(max(max.level + 1L, 3L), palette = "Dark 3")
  palette <- palette[seq_len(max.level + 1L)]
  stats::setNames(palette, as.character(seq.int(0L, max.level, by = 1L)))
}

gripui.gmds.level.choices <- function(prepared) {
  levels <- seq.int(from = length(prepared$misf$levels) - 1L, to = 0L, by = -1L)
  stats::setNames(as.character(levels), sprintf("V_%d", levels))
}

gripui.gmds.method.catalog <- function() {
  list(
    gmds = list(
      id = "gmds",
      label = "GMDS",
      objective_label = "Geodesic MDS",
      refinement_phrase = "pure-GMDS",
      prepare_fn = "grip.prepare.misf.geodesic.mds",
      optimize_fn = "grip.optimize.misf.geodesic.mds",
      fit_class = "grip_misf_gmds_fit"
    ),
    gkk = list(
      id = "gkk",
      label = "GKK",
      objective_label = "Full geodesic-KK",
      refinement_phrase = "full geodesic-KK",
      pair_mode = "full",
      prepare_fn = "grip.prepare.misf.geodesic.kk",
      optimize_fn = "grip.optimize.misf.geodesic.kk",
      fit_class = "grip_misf_gkk_fit"
    ),
    lgkk = list(
      id = "lgkk",
      label = "LGKK",
      objective_label = "Landmark geodesic-KK",
      refinement_phrase = "landmark geodesic-KK",
      pair_mode = "landmark",
      prepare_fn = "grip.prepare.misf.geodesic.kk",
      optimize_fn = "grip.optimize.misf.geodesic.kk",
      fit_class = "grip_misf_gkk_fit"
    )
  )
}

gripui.gmds.method.spec <- function(method = "gmds") {
  method <- as.character(method[[1L]])
  catalog <- gripui.gmds.method.catalog()
  if (!method %in% names(catalog)) {
    stop("unknown GMDS app method: ", method)
  }
  catalog[[method]]
}

gripui.gmds.method.choices <- function() {
  catalog <- gripui.gmds.method.catalog()
  stats::setNames(
    vapply(catalog, `[[`, character(1L), "id"),
    vapply(catalog, `[[`, character(1L), "label")
  )
}

gripui.gmds.min.hop.separation <- function(level) {
  level <- as.integer(level)
  if (level <= 0L) {
    return(1L)
  }
  as.integer(2^(level - 1L) + 1L)
}

gripui.gmds.level.table <- function(prepared) {
  levels <- seq.int(from = length(prepared$misf$levels) - 1L, to = 0L, by = -1L)
  rows <- lapply(levels, function(level) {
    vertices <- as.integer(prepared$misf$levels[[level + 1L]])
    next.vertices <- if (level < length(prepared$misf$levels) - 1L) {
      as.integer(prepared$misf$levels[[level + 2L]])
    } else {
      integer(0L)
    }
    role <- character(0L)
    if (identical(level, prepared$coarsest_level_level)) {
      role <- c(role, "Coarsest filtration level")
    }
    if (identical(level, prepared$top_level_level)) {
      role <- c(role, "Top solve level")
    }
    if (identical(level, 0L)) {
      role <- c(role, "Full graph")
    }
    if (!length(role)) {
      role <- "Intermediate level"
    }
    data.frame(
      level = sprintf("V_%d", level),
      size = length(vertices),
      new_vertices = if (length(next.vertices)) length(setdiff(vertices, next.vertices)) else length(vertices),
      min_hop_sep = gripui.gmds.min.hop.separation(level),
      role = paste(role, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

gripui.gmds.bundle.summary <- function(bundle) {
  prepared <- bundle$prepared
  payload <- bundle$payload
  method <- bundle$method
  coarsest.vertices <- prepared$misf$levels[[prepared$coarsest_level_level + 1L]]
  gripui.family.summary.table(list(
    method = method$label,
    objective = method$objective_label,
    family = payload$family_label,
    category = payload$category,
    vertices = payload$n,
    edges = gripui.family.count.edges(payload$edges),
    coarsest_level = sprintf("V_%d (|V| = %d)", prepared$coarsest_level_level, length(coarsest.vertices)),
    top_solve_level = sprintf("V_%d (|V| = %d)", prepared$top_level_level, length(prepared$top_level_vertices)),
    top_level_reason = prepared$top_level_selection_reason,
    trace_rows = nrow(bundle$stage_trace),
    elapsed_sec = sprintf("%.3f", as.double(bundle$fit$timing$total))
  ))
}

gripui.gmds.stage.trace.table <- function(bundle) {
  trace <- bundle$stage_trace
  if (is.null(trace) || !nrow(trace)) {
    return(data.frame())
  }
  keep <- intersect(
    c("stage_index", "label", "level", "active_n", "inserted_n", "weighted_rel_rmse", "elapsed_sec"),
    names(trace)
  )
  out <- trace[, keep, drop = FALSE]
  names(out) <- c("Step", "Label", "Level", "Active", "Inserted", "Weighted_RMSE", "Elapsed_sec")[seq_along(keep)]
  if ("Weighted_RMSE" %in% names(out)) {
    out$Weighted_RMSE <- ifelse(is.finite(out$Weighted_RMSE), sprintf("%.4f", out$Weighted_RMSE), "")
  }
  if ("Elapsed_sec" %in% names(out)) {
    out$Elapsed_sec <- ifelse(is.finite(out$Elapsed_sec), sprintf("%.3f", out$Elapsed_sec), "")
  }
  if ("Inserted" %in% names(out)) {
    out$Inserted[is.na(out$Inserted)] <- ""
  }
  out
}

gripui.gmds.selected.level.table <- function(bundle, level) {
  prepared <- bundle$prepared
  level <- as.integer(level)
  vertices <- as.integer(prepared$misf$levels[[level + 1L]])
  next.vertices <- if (level < length(prepared$misf$levels) - 1L) {
    as.integer(prepared$misf$levels[[level + 2L]])
  } else {
    integer(0L)
  }
  role <- character(0L)
  if (identical(level, prepared$coarsest_level_level)) {
    role <- c(role, "coarsest filtration level")
  }
  if (identical(level, prepared$top_level_level)) {
    role <- c(role, "top solve level")
  }
  if (identical(level, 0L)) {
    role <- c(role, "full graph")
  }
  gripui.family.summary.table(list(
    level = sprintf("V_%d", level),
    size = length(vertices),
    new_vertices = if (length(next.vertices)) length(setdiff(vertices, next.vertices)) else length(vertices),
    min_hop_sep = gripui.gmds.min.hop.separation(level),
    role = if (length(role)) paste(role, collapse = "; ") else "intermediate level"
  ))
}

gripui.gmds.stage.payload <- function(bundle, state) {
  if (is.null(bundle) || is.null(bundle$stage_payloads)) {
    return(NULL)
  }
  bundle$stage_payloads[[state]]
}

gripui.gmds.level.stage.payload <- function(bundle, stage, level) {
  if (is.null(bundle) || is.null(bundle$fit)) {
    return(NULL)
  }
  record <- grip.geodesic.misf.trace.stage.lookup(bundle$fit, stage = stage, level = as.integer(level))
  if (is.null(record)) {
    return(NULL)
  }
  coords.full <- if (is.null(record$coords_full)) NULL else as.matrix(record$coords_full)
  if (is.null(coords.full)) {
    display.coords <- NULL
    rmse <- NA_real_
  } else {
    aligned <- grip.geodesic.misf.align.partial.to.target(
      coords = coords.full,
      target = bundle$payload$coords_display,
      allow_reflection = TRUE
    )
    display.coords <- aligned$aligned
    rmse <- aligned$rmse
  }
  list(
    state = stage,
    stage = record$stage,
    stage_key = record$stage_key,
    label = record$label,
    level = as.integer(record$level),
    active_vertices = as.integer(record$active_vertices),
    inserted_vertices = as.integer(record$inserted_vertices),
    active_edges = grip.geodesic.misf.filter.edge.matrix(bundle$prepared$edges, record$active_vertices),
    coords = coords.full,
    display_coords = display.coords,
    pair_mode = record$pair_mode,
    rmse = rmse,
    summary = record$summary
  )
}

gripui.gmds.expansion.levels <- function(bundle) {
  if (is.null(bundle) || is.null(bundle$stage_data) || !length(bundle$stage_data)) {
    return(integer())
  }
  insertion.levels <- unique(vapply(Filter(function(record) {
    identical(record$stage, "insertion")
  }, bundle$stage_data), function(record) as.integer(record$level), integer(1L)))
  refinement.levels <- unique(vapply(Filter(function(record) {
    identical(record$stage, "refinement")
  }, bundle$stage_data), function(record) as.integer(record$level), integer(1L)))
  levels <- intersect(insertion.levels, refinement.levels)
  sort(levels, decreasing = TRUE)
}

gripui.gmds.expansion.level.choices <- function(bundle) {
  levels <- gripui.gmds.expansion.levels(bundle)
  if (!length(levels)) {
    return(stats::setNames(character(), character()))
  }
  stats::setNames(as.character(levels), sprintf("V_%d", levels))
}

gripui.gmds.stage.metric.value <- function(x, digits = 4L) {
  if (length(x) != 1L || !is.finite(x)) {
    return("")
  }
  formatC(x, format = "f", digits = digits)
}

gripui.gmds.stage.display.coords <- function(bundle, payload, fill_reference = TRUE) {
  coords <- payload$display_coords
  if (is.null(coords)) {
    return(NULL)
  }
  coords <- gripui.family.pad_coords(coords, target_cols = 3L)
  if (!isTRUE(fill_reference) || is.null(bundle$payload$coords_display)) {
    return(coords)
  }
  reference <- gripui.family.pad_coords(bundle$payload$coords_display, target_cols = 3L)
  keep <- gripui.valid.coord.rows(coords)
  reference[keep, ] <- coords[keep, , drop = FALSE]
  reference
}

gripui.gmds.stage.note <- function(bundle, state) {
  prepared <- bundle$prepared
  refinement.phrase <- bundle$method$refinement_phrase
  switch(
    state,
    seed = if (prepared$top_level_level < prepared$coarsest_level_level) {
      sprintf(
        "The highlighted d+1 seed is chosen on V_%d rather than on the coarsest filtration level V_%d because V_%d is undersized for a %dD start.",
        prepared$top_level_level,
        prepared$coarsest_level_level,
        prepared$coarsest_level_level,
        bundle$settings$dim
      )
    } else {
      sprintf(
        "The highlighted d+1 seed is chosen directly on the coarsest filtration level V_%d.",
        prepared$top_level_level
      )
    },
    top_level = sprintf(
      "The selected top solve level V_%d after the restartable %s solve that starts from the initial placement shown at left.",
      prepared$top_level_level,
      refinement.phrase
    ),
    initial_placement = sprintf(
      "All vertices of V_%d after geometric seed placement and anchor-based insertion, before the top-level %s refinement.",
      prepared$top_level_level,
      refinement.phrase
    ),
    ""
  )
}

gripui.gmds.level.stage.note <- function(bundle, stage, level, inserted_n = NULL) {
  stage <- match.arg(stage, c("insertion", "refinement"))
  level <- as.integer(level)
  scaffold.level <- level + 1L
  refinement.phrase <- bundle$method$refinement_phrase
  switch(
    stage,
    insertion = if (!is.null(inserted_n) && is.finite(inserted_n)) {
      sprintf(
        "This panel shows the expansion from V_%d to V_%d. %d newly activated vertices were inserted into the scaffold inherited from V_%d.",
        scaffold.level,
        level,
        as.integer(inserted_n),
        scaffold.level
      )
    } else {
      sprintf(
        "This panel shows the expansion from V_%d to V_%d using the scaffold inherited from V_%d.",
        scaffold.level,
        level,
        scaffold.level
      )
    },
    refinement = sprintf(
      "This panel shows the active-level %s refinement on V_%d after insertion of that level.",
      refinement.phrase,
      level
    )
  )
}

gripui.gmds.stage.summary.table <- function(payload) {
  if (is.null(payload)) {
    return(data.frame())
  }
  summary <- payload$summary
  if (is.null(summary)) {
    summary <- list()
  }
  gripui.family.summary.table(list(
    label = payload$label,
    level = sprintf("V_%d", payload$level),
    active_vertices = length(payload$active_vertices),
    inserted_vertices = if (length(payload$inserted_vertices)) length(payload$inserted_vertices) else if (!is.null(summary$inserted_n) && is.finite(summary$inserted_n)) as.integer(summary$inserted_n) else "",
    pair_mode = gripui.family.value_or_default(payload$pair_mode, ""),
    rmse = gripui.gmds.stage.metric.value(payload$rmse, 4L),
    pair_constraints = if (!is.null(summary$pair_n) && is.finite(summary$pair_n)) as.integer(summary$pair_n) else "",
    energy = gripui.gmds.stage.metric.value(summary$energy, 4L),
    stress = gripui.gmds.stage.metric.value(summary$stress, 4L),
    weighted_rel_rmse = gripui.gmds.stage.metric.value(summary$weighted_rel_rmse, 4L),
    elapsed_sec = gripui.gmds.stage.metric.value(summary$elapsed_sec, 3L)
  ))
}

gripui.gmds.repro.code <- function(bundle) {
  payload <- bundle$payload
  prepared <- bundle$prepared
  settings <- bundle$settings
  if (identical(bundle$method$id, "gmds")) {
    return(paste(
      payload$code,
      "",
      sprintf(
        "prepared <- grip.prepare.misf.geodesic.mds(edges = graph$edges, n = graph$n, edge_weights = graph$edge_weights, tie_mode = \"average\", num_init = %dL, dim = %dL, top_level_mode = \"skip\", seed = %dL)",
        settings$num_init,
        settings$dim,
        settings$prepare_seed
      ),
      sprintf(
        paste0(
          "fit <- grip.optimize.misf.geodesic.mds(prepared = prepared, dim = %dL, ",
          "top_level_restarts = 1L, top_level_max_iter = %dL, top_level_engine = \"cpp\", ",
          "insertion_anchor_policy = \"prev_level_spread\", insertion_max_iter = %dL, ",
          "refinement_local_nbrs = %dL, refinement_landmark_count = %dL, ",
          "refinement_pair_mode = \"sparse\", refinement_anchor_weight = 0.05, ",
          "refinement_anchor_weight_end = 0.01, refinement_continuation = \"linear\", ",
          "refinement_max_iter = %dL, refinement_engine = \"cpp\", ",
          "final_polish_max_iter = %dL, final_polish_engine = \"cpp\", ",
          "n_threads = %dL, return_trace = TRUE, return_frames = FALSE, seed = %dL)"
        ),
        settings$dim,
        settings$top_level_max_iter,
        settings$insertion_max_iter,
        settings$refinement_local_nbrs,
        settings$refinement_landmark_count,
        settings$refinement_max_iter,
        settings$final_polish_max_iter,
        settings$n_threads,
        settings$optimizer_seed
      ),
      sprintf(
        "# selected top solve level: V_%d (reason: %s)",
        prepared$top_level_level,
        prepared$top_level_selection_reason
      ),
      sep = "\n"
    ))
  }

  paste(
    payload$code,
    "",
    sprintf(
      "prepared <- grip.prepare.misf.geodesic.kk(edges = graph$edges, n = graph$n, edge_weights = graph$edge_weights, tie_mode = \"average\", num_init = %dL, dim = %dL, top_level_mode = \"skip\", top_level_pair_mode = \"%s\", top_level_local_nbrs = %dL, top_level_landmark_count = %dL, seed = %dL)",
      settings$num_init,
      settings$dim,
      settings$pair_mode,
      settings$refinement_local_nbrs,
      settings$refinement_landmark_count,
      settings$prepare_seed
    ),
    sprintf(
      paste0(
        "fit <- grip.optimize.misf.geodesic.kk(prepared = prepared, dim = %dL, ",
        "top_level_pair_mode = \"%s\", top_level_full_limit = %dL, top_level_local_nbrs = %dL, top_level_landmark_count = %dL, ",
        "top_level_restarts = 1L, top_level_max_iter = %dL, ",
        "insertion_anchor_policy = \"prev_level_spread\", insertion_max_iter = %dL, ",
        "refinement_pair_mode = \"%s\", refinement_full_limit = %dL, refinement_local_nbrs = %dL, refinement_landmark_count = %dL, ",
        "refinement_anchor_weight = 0.05, refinement_anchor_weight_end = 0.01, refinement_continuation = \"linear\", ",
        "refinement_max_iter = %dL, final_pair_mode = \"%s\", final_full_limit = %dL, final_local_nbrs = %dL, final_landmark_count = %dL, final_max_iter = %dL, ",
        "return_trace = TRUE, return_frames = FALSE, seed = %dL)"
      ),
      settings$dim,
      settings$pair_mode,
      settings$pair_full_limit,
      settings$refinement_local_nbrs,
      settings$refinement_landmark_count,
      settings$top_level_max_iter,
      settings$insertion_max_iter,
      settings$pair_mode,
      settings$pair_full_limit,
      settings$refinement_local_nbrs,
      settings$refinement_landmark_count,
      settings$refinement_max_iter,
      settings$pair_mode,
      settings$pair_full_limit,
      settings$refinement_local_nbrs,
      settings$refinement_landmark_count,
      settings$final_polish_max_iter,
      settings$optimizer_seed
    ),
    sprintf(
      "# selected top solve level: V_%d (reason: %s)",
      prepared$top_level_level,
      prepared$top_level_selection_reason
    ),
    sep = "\n"
  )
}

gripui.gmds.render.widget <- function(coords,
                                      edges,
                                      vertex_colors,
                                      vertex_alpha = 0.95,
                                      edge_alpha = 0.25,
                                      highlight_vertices = integer(0L),
                                      highlight_edge_color = "#8a5a44",
                                      highlight_size = 8,
                                      base_size = 5,
                                      background_alpha = NULL) {
  old <- gripui.enable.rgl.null.device()
  on.exit(options(rgl.useNULL = old), add = TRUE)

  coords <- gripui.family.pad_coords(coords, target_cols = 3L)
  edges <- as.matrix(edges)
  storage.mode(edges) <- "integer"
  n <- nrow(coords)
  keep <- gripui.valid.coord.rows(coords)
  if (is.null(background_alpha)) {
    background_alpha <- vertex_alpha
  }

  base_cols <- grDevices::adjustcolor(vertex_colors, alpha.f = background_alpha)
  highlight_cols <- grDevices::adjustcolor(vertex_colors, alpha.f = vertex_alpha)
  highlight_vertices <- as.integer(highlight_vertices)
  highlight_vertices <- highlight_vertices[highlight_vertices >= 1L & highlight_vertices <= n]

  rgl::open3d(useNULL = TRUE)
  on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)
  rgl::clear3d(type = "shapes")
  rgl::clear3d(type = "lights")
  rgl::bg3d(color = "white")
  rgl::light3d()

  if (nrow(edges) > 0L) {
    good.edges <- keep[edges[, 1L]] & keep[edges[, 2L]]
    edges <- edges[good.edges, , drop = FALSE]
    if (nrow(edges) > 0L) {
      for (i in seq_len(nrow(edges))) {
        e <- edges[i, ]
        rgl::segments3d(
          x = coords[e, 1L],
          y = coords[e, 2L],
          z = coords[e, 3L],
          col = grDevices::adjustcolor("#9ca3af", alpha.f = edge_alpha),
          lwd = 1
        )
      }
    }
  }

  if (any(keep)) {
    rgl::points3d(
      coords[keep, 1L],
      coords[keep, 2L],
      coords[keep, 3L],
      col = base_cols[keep],
      size = base_size
    )
  }

  if (length(highlight_vertices)) {
    selected.edges <- gripui.gmds.filter.edge.matrix(edges, highlight_vertices)
    if (nrow(selected.edges) > 0L) {
      for (i in seq_len(nrow(selected.edges))) {
        e <- selected.edges[i, ]
        rgl::segments3d(
          x = coords[e, 1L],
          y = coords[e, 2L],
          z = coords[e, 3L],
          col = grDevices::adjustcolor(highlight_edge_color, alpha.f = min(1, edge_alpha + 0.35)),
          lwd = 3
        )
      }
    }
    selected.keep <- keep[highlight_vertices]
    selected.vertices <- highlight_vertices[selected.keep]
    if (length(selected.vertices)) {
      rgl::points3d(
        coords[selected.vertices, 1L],
        coords[selected.vertices, 2L],
        coords[selected.vertices, 3L],
        col = highlight_cols[selected.vertices],
        size = highlight_size
      )
    }
  }

  rgl::rglwidget(height = 520)
}

gripui.gmds.compute.bundle <- function(desc,
                                       values,
                                       method = "gmds",
                                       dim = 3L,
                                       num_init = 24L,
                                       prepare_seed = NULL,
                                       optimizer_seed = NULL,
                                       top_level_max_iter = 4L,
                                       insertion_max_iter = 24L,
                                       refinement_max_iter = 2L,
                                       final_polish_max_iter = 2L,
                                       n_threads = 0L) {
  payload <- gripui.family.build.payload(desc, values)
  method.spec <- gripui.gmds.method.spec(method)
  dim <- as.integer(round(dim))
  num_init <- as.integer(round(num_init))
  top_level_max_iter <- as.integer(round(top_level_max_iter))
  insertion_max_iter <- as.integer(round(insertion_max_iter))
  refinement_max_iter <- as.integer(round(refinement_max_iter))
  final_polish_max_iter <- as.integer(round(final_polish_max_iter))
  n_threads <- as.integer(round(n_threads))
  if (is.null(prepare_seed) || !length(prepare_seed) || !is.finite(prepare_seed)) {
    prepare_seed <- gripui.gmds.default.seed(values, offset = 1000L)
  }
  if (is.null(optimizer_seed) || !length(optimizer_seed) || !is.finite(optimizer_seed)) {
    optimizer_seed <- gripui.gmds.default.seed(values, offset = 2000L)
  }
  prepare_seed <- as.integer(round(prepare_seed))
  optimizer_seed <- as.integer(round(optimizer_seed))
  refinement_local_nbrs <- max(4L, dim + 1L)
  top_level_landmark_count <- max(2L, min(8L, num_init))

  if (identical(method.spec$id, "gmds")) {
    prepared <- grip.prepare.misf.geodesic.mds(
      edges = payload$edges,
      n = payload$n,
      edge_weights = payload$edge_weights,
      tie_mode = "average",
      num_init = num_init,
      dim = dim,
      top_level_mode = "skip",
      seed = prepare_seed
    )

    refinement_landmark_count <- max(2L, min(8L, length(prepared$top_level_vertices)))

    fit <- grip.optimize.misf.geodesic.mds(
      prepared = prepared,
      dim = dim,
      top_level_restarts = 1L,
      top_level_max_iter = top_level_max_iter,
      top_level_engine = "cpp",
      insertion_anchor_policy = "prev_level_spread",
      insertion_max_iter = insertion_max_iter,
      refinement_local_nbrs = refinement_local_nbrs,
      refinement_landmark_count = refinement_landmark_count,
      refinement_pair_mode = "sparse",
      refinement_anchor_weight = 0.05,
      refinement_anchor_weight_end = 0.01,
      refinement_continuation = "linear",
      refinement_max_iter = refinement_max_iter,
      refinement_engine = "cpp",
      final_polish_max_iter = final_polish_max_iter,
      final_polish_engine = "cpp",
      n_threads = n_threads,
      return_trace = TRUE,
      return_frames = FALSE,
      seed = optimizer_seed
    )
    pair.mode <- NA_character_
    pair.full.limit <- NA_integer_
  } else {
    prepared <- grip.prepare.misf.geodesic.kk(
      edges = payload$edges,
      n = payload$n,
      edge_weights = payload$edge_weights,
      tie_mode = "average",
      num_init = num_init,
      dim = dim,
      top_level_mode = "skip",
      top_level_pair_mode = method.spec$pair_mode,
      top_level_full_limit = max(payload$n, num_init),
      top_level_local_nbrs = refinement_local_nbrs,
      top_level_landmark_count = top_level_landmark_count,
      seed = prepare_seed
    )

    refinement_landmark_count <- max(2L, min(8L, length(prepared$top_level_vertices)))
    pair.mode <- method.spec$pair_mode
    pair.full.limit <- max(payload$n, length(prepared$top_level_vertices))

    fit <- grip.optimize.misf.geodesic.kk(
      prepared = prepared,
      dim = dim,
      top_level_pair_mode = pair.mode,
      top_level_full_limit = pair.full.limit,
      top_level_local_nbrs = refinement_local_nbrs,
      top_level_landmark_count = refinement_landmark_count,
      top_level_restarts = 1L,
      top_level_max_iter = top_level_max_iter,
      insertion_anchor_policy = "prev_level_spread",
      insertion_max_iter = insertion_max_iter,
      refinement_pair_mode = pair.mode,
      refinement_full_limit = pair.full.limit,
      refinement_local_nbrs = refinement_local_nbrs,
      refinement_landmark_count = refinement_landmark_count,
      refinement_anchor_weight = 0.05,
      refinement_anchor_weight_end = 0.01,
      refinement_continuation = "linear",
      refinement_max_iter = refinement_max_iter,
      final_pair_mode = pair.mode,
      final_full_limit = pair.full.limit,
      final_local_nbrs = refinement_local_nbrs,
      final_landmark_count = refinement_landmark_count,
      final_max_iter = final_polish_max_iter,
      return_trace = TRUE,
      return_frames = FALSE,
      seed = optimizer_seed
    )
  }

  list(
    payload = payload,
    method = method.spec,
    prepared = prepared,
    fit = fit,
    stage_trace = fit$stage_trace,
    stage_data = grip.geodesic.misf.trace.stage.data(fit),
    stage_payloads = grip.geodesic.misf.trace.stage.payloads(
      fit,
      target = payload$coords_display,
      states = c("seed", "initial_placement", "top_level", "after_insertion", "after_refinement", "final_polish")
    ),
    settings = list(
      dim = dim,
      method = method.spec$id,
      num_init = num_init,
      prepare_seed = prepare_seed,
      optimizer_seed = optimizer_seed,
      top_level_max_iter = top_level_max_iter,
      insertion_max_iter = insertion_max_iter,
      refinement_max_iter = refinement_max_iter,
      final_polish_max_iter = final_polish_max_iter,
      pair_mode = pair.mode,
      pair_full_limit = pair.full.limit,
      refinement_local_nbrs = refinement_local_nbrs,
      refinement_landmark_count = refinement_landmark_count,
      n_threads = n_threads
    )
  )
}

gripui.gmds.ui <- function(catalog, title, subtitle = NULL) {
  initial_state <- gripui.gmds.initial.state(catalog)
  categories <- c("All families", gripui.family.categories(catalog))
  css.path <- system.file("app/www/gripui.css", package = "grip")

  bslib::page_sidebar(
    title = shiny::div(
      class = "gripui-titlebar",
      shiny::div(class = "gripui-brand", title)
    ),
    theme = gripui.family.theme(),
    sidebar = bslib::sidebar(
      width = 390,
      if (!is.null(subtitle) && nzchar(subtitle)) {
        shiny::tags$p(class = "gripui-project-summary", subtitle)
      },
      shiny::uiOutput("gmds_family_meta"),
      shiny::selectInput("gmds_family_category", "Category", choices = categories, selected = initial_state$category),
      shiny::selectInput("gmds_family_id", "Family", choices = initial_state$choices, selected = initial_state$family_id),
      shiny::selectInput("gmds_family_preset", "Preset", choices = gripui.family.preset.choices(initial_state$desc), selected = "default"),
      shiny::uiOutput("gmds_family_param_panel"),
      shiny::tags$hr(),
      shiny::tags$h5("Method Controls"),
      shiny::selectInput("gmds_method", "Method", choices = gripui.gmds.method.choices(), selected = "gmds"),
      shiny::numericInput("gmds_dim", "Embedding dimension", value = 3L, min = 2L, max = 3L, step = 1L),
      shiny::numericInput("gmds_num_init", "MIS filtration num_init", value = 24L, min = 3L, step = 1L),
      shiny::numericInput("gmds_prepare_seed", "MIS filtration seed", value = 1001L, min = 0L, step = 1L),
      shiny::numericInput("gmds_optimizer_seed", "Optimizer seed", value = 2001L, min = 0L, step = 1L),
      shiny::numericInput("gmds_top_level_max_iter", "Top-level max iter", value = 4L, min = 0L, step = 1L),
      shiny::numericInput("gmds_insertion_max_iter", "Insertion max iter", value = 24L, min = 0L, step = 1L),
      shiny::numericInput("gmds_refinement_max_iter", "Refinement max iter", value = 2L, min = 0L, step = 1L),
      shiny::numericInput("gmds_final_polish_max_iter", "Final polish max iter", value = 2L, min = 0L, step = 1L),
      shiny::numericInput("gmds_n_threads", "Threads (0 = auto)", value = 0L, min = 0L, step = 1L),
      shiny::tags$hr(),
      shiny::tags$h5("Display"),
      shiny::sliderInput("gmds_vertex_alpha", "Selected vertex opacity", min = 0.10, max = 1, value = 0.95, step = 0.05),
      shiny::sliderInput("gmds_background_alpha", "Background vertex opacity", min = 0.05, max = 1, value = 0.35, step = 0.05),
      shiny::sliderInput("gmds_edge_alpha", "Edge opacity", min = 0.05, max = 1, value = 0.20, step = 0.05),
      shiny::uiOutput("gmds_level_control"),
      shiny::uiOutput("gmds_expansion_level_control"),
      shiny::actionButton("render_gmds_case", "Render MISF stages", class = "btn-primary"),
      shiny::uiOutput("gmds_render_status")
    ),
    if (nzchar(css.path)) shiny::tags$head(shiny::includeCSS(css.path)),
    bslib::card(
      full_screen = TRUE,
      class = "gripui-card",
      bslib::card_header("Bundle Summary"),
      shiny::tableOutput("gmds_bundle_summary")
    ),
    bslib::accordion(
      bslib::accordion_panel(
        "Stage 0: Graph and Geometry",
        bslib::layout_columns(
          col_widths = c(7, 5),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header("Reference graph geometry"),
            shiny::uiOutput("gmds_reference_view")
          ),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header("Reference summary"),
            shiny::tableOutput("gmds_reference_summary")
          )
        )
      ),
      bslib::accordion_panel(
        "Stage 1: MIS Filtration",
        bslib::layout_columns(
          col_widths = c(8, 4),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header(shiny::uiOutput("gmds_misf_title")),
            shiny::uiOutput("gmds_misf_view")
          ),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header("MIS filtration levels"),
            shiny::tableOutput("gmds_misf_table")
          )
        ),
        bslib::card(
          full_screen = TRUE,
          class = "gripui-card",
          bslib::card_header("Selected level detail"),
          shiny::tableOutput("gmds_selected_level_summary")
        )
      ),
      bslib::accordion_panel(
        "Stage 2: Seed Selection and Initial Placement",
        bslib::layout_columns(
          col_widths = c(4, 4, 4),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header(shiny::uiOutput("gmds_seed_title")),
            shiny::uiOutput("gmds_seed_view")
          ),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header(shiny::uiOutput("gmds_initial_title")),
            shiny::uiOutput("gmds_initial_view")
          ),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header(shiny::uiOutput("gmds_top_level_title")),
            shiny::uiOutput("gmds_top_level_view")
          )
        ),
        bslib::layout_columns(
          col_widths = c(4, 4, 4),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header("Seed summary"),
            shiny::tableOutput("gmds_seed_summary")
          ),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header("Initial placement summary"),
            shiny::tableOutput("gmds_initial_summary")
          ),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header("Top-level solve summary"),
            shiny::tableOutput("gmds_top_level_summary")
          )
        )
      ),
      bslib::accordion_panel(
        "Stage 3: Per-Level Insertion and Refinement",
        bslib::layout_columns(
          col_widths = c(6, 6),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header(shiny::uiOutput("gmds_insertion_title")),
            shiny::uiOutput("gmds_insertion_view")
          ),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header(shiny::uiOutput("gmds_refinement_title")),
            shiny::uiOutput("gmds_refinement_view")
          )
        ),
        bslib::layout_columns(
          col_widths = c(6, 6),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header("Insertion summary"),
            shiny::tableOutput("gmds_insertion_summary")
          ),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header("Refinement summary"),
            shiny::tableOutput("gmds_refinement_summary")
          )
        )
      ),
      bslib::accordion_panel(
        "Canonical Trace Summary",
        bslib::layout_columns(
          col_widths = c(7, 5),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header("Stage trace"),
            shiny::tableOutput("gmds_stage_trace")
          ),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header("Reproducible R call"),
            shiny::verbatimTextOutput("gmds_code")
          )
        )
      )
    )
  )
}

gripui.gmds.server <- function(catalog) {
  function(input, output, session) {
    bundle_state <- shiny::reactiveVal(NULL)
    error_state <- shiny::reactiveVal(NULL)

    current_desc <- shiny::reactive({
      catalog[[input$gmds_family_id]]
    })

    current_level <- shiny::reactive({
      bundle <- bundle_state()
      if (is.null(bundle)) {
        return(0L)
      }
      selected <- suppressWarnings(as.integer(input$gmds_focus_level))
      if (!length(selected) || !is.finite(selected)) {
        return(bundle$prepared$top_level_level)
      }
      selected
    })

    current_expansion_level <- shiny::reactive({
      bundle <- bundle_state()
      if (is.null(bundle)) {
        return(0L)
      }
      available <- gripui.gmds.expansion.levels(bundle)
      if (!length(available)) {
        return(bundle$prepared$top_level_level)
      }
      selected <- suppressWarnings(as.integer(input$gmds_expansion_level))
      if (!length(selected) || !is.finite(selected) || !selected %in% available) {
        return(available[[1L]])
      }
      selected
    })

    build_bundle <- function(values) {
      desc <- current_desc()
      result <- tryCatch(
        gripui.gmds.compute.bundle(
          desc = desc,
          values = values,
          method = gripui.family.value_or_default(input$gmds_method, "gmds"),
          dim = gripui.family.value_or_default(input$gmds_dim, 3L),
          num_init = gripui.family.value_or_default(input$gmds_num_init, 24L),
          prepare_seed = gripui.family.value_or_default(input$gmds_prepare_seed, gripui.gmds.default.seed(values, offset = 1000L)),
          optimizer_seed = gripui.family.value_or_default(input$gmds_optimizer_seed, gripui.gmds.default.seed(values, offset = 2000L)),
          top_level_max_iter = gripui.family.value_or_default(input$gmds_top_level_max_iter, 4L),
          insertion_max_iter = gripui.family.value_or_default(input$gmds_insertion_max_iter, 24L),
          refinement_max_iter = gripui.family.value_or_default(input$gmds_refinement_max_iter, 2L),
          final_polish_max_iter = gripui.family.value_or_default(input$gmds_final_polish_max_iter, 2L),
          n_threads = gripui.family.value_or_default(input$gmds_n_threads, 0L)
        ),
        error = function(e) e
      )
      if (inherits(result, "error")) {
        error_state(conditionMessage(result))
        shiny::showNotification(conditionMessage(result), type = "error")
        return(invisible(NULL))
      }
      bundle_state(result)
      error_state(NULL)
      shiny::updateSelectInput(
        session,
        "gmds_focus_level",
        choices = gripui.gmds.level.choices(result$prepared),
        selected = as.character(result$prepared$top_level_level)
      )
      expansion.choices <- gripui.gmds.expansion.level.choices(result)
      if (length(expansion.choices)) {
        shiny::updateSelectInput(
          session,
          "gmds_expansion_level",
          choices = expansion.choices,
          selected = unname(expansion.choices[[1L]])
        )
      }
      invisible(result)
    }

    shiny::observe({
      choices <- gripui.family.choices(catalog, input$gmds_family_category)
      current <- shiny::isolate(input$gmds_family_id)
      selected <- if (!is.null(current) && current %in% unname(choices)) current else unname(choices[[1L]])
      shiny::updateSelectInput(session, "gmds_family_id", choices = choices, selected = selected)
    })

    shiny::observe({
      desc <- current_desc()
      shiny::updateSelectInput(
        session,
        "gmds_family_preset",
        choices = gripui.family.preset.choices(desc),
        selected = "default"
      )
    })

    shiny::observeEvent(list(input$gmds_family_id, input$gmds_family_preset, input$gmds_method), {
      desc <- current_desc()
      values <- .gripui.family.merge.values(desc, preset_id = input$gmds_family_preset)
      shiny::updateNumericInput(session, "gmds_prepare_seed", value = gripui.gmds.default.seed(values, offset = 1000L))
      shiny::updateNumericInput(session, "gmds_optimizer_seed", value = gripui.gmds.default.seed(values, offset = 2000L))
      build_bundle(values)
    }, ignoreInit = FALSE)

    shiny::observeEvent(input$render_gmds_case, {
      desc <- current_desc()
      values <- gripui.family.collect.values(desc, input, preset_id = input$gmds_family_preset)
      build_bundle(values)
    })

    output$gmds_family_meta <- shiny::renderUI({
      desc <- current_desc()
      shiny::tagList(
        shiny::tags$div(class = "gripui-project-title", desc$label),
        shiny::tags$p(class = "gripui-project-summary", desc$summary),
        shiny::tags$p(
          class = "gripui-project-summary",
          shiny::tags$strong("Source: "),
          shiny::tags$code(desc$function_name),
          shiny::tags$br(),
          shiny::tags$code(desc$implementation)
        )
      )
    })

    output$gmds_family_param_panel <- shiny::renderUI({
      desc <- current_desc()
      values <- .gripui.family.merge.values(desc, preset_id = input$gmds_family_preset)
      gripui.family.param.ui(desc, values)
    })

    output$gmds_level_control <- shiny::renderUI({
      bundle <- bundle_state()
      if (is.null(bundle)) {
        return(NULL)
      }
      shiny::selectInput(
        "gmds_focus_level",
        "Focus filtration level",
        choices = gripui.gmds.level.choices(bundle$prepared),
        selected = as.character(bundle$prepared$top_level_level)
      )
    })

    output$gmds_expansion_level_control <- shiny::renderUI({
      bundle <- bundle_state()
      if (is.null(bundle)) {
        return(NULL)
      }
      choices <- gripui.gmds.expansion.level.choices(bundle)
      if (!length(choices)) {
        return(NULL)
      }
      shiny::selectInput(
        "gmds_expansion_level",
        "Expansion/refinement level",
        choices = choices,
        selected = unname(choices[[1L]])
      )
    })

    output$gmds_render_status <- shiny::renderUI({
      err <- error_state()
      bundle <- bundle_state()
      if (!is.null(err)) {
        return(shiny::tags$p(style = "color:#8a1c1c;", err))
      }
      if (is.null(bundle)) {
        return(shiny::tags$p(class = "gripui-muted", "Render a family to inspect its canonical MISF stage traces."))
      }
      prepared <- bundle$prepared
      method.label <- bundle$method$label
      lead <- if (prepared$top_level_level < prepared$coarsest_level_level) {
        sprintf(
          "True MIS filtration starts at V_%d, but the %s solve starts from the coarsest admissible level V_%d.",
          prepared$coarsest_level_level,
          method.label,
          prepared$top_level_level
        )
      } else {
        sprintf("The %s solve starts from the coarsest filtration level V_%d.", method.label, prepared$top_level_level)
      }
      shiny::tags$p(
        class = "gripui-selection-status",
        sprintf(
          "%s Rendered with %d vertices, %d edges, and %d canonical stages.",
          lead,
          bundle$payload$n,
          gripui.family.count.edges(bundle$payload$edges),
          nrow(bundle$stage_trace)
        )
      )
    })

    output$gmds_bundle_summary <- shiny::renderTable({
      bundle <- bundle_state()
      shiny::req(bundle)
      gripui.gmds.bundle.summary(bundle)
    }, striped = TRUE, bordered = FALSE, spacing = "s")

    output$gmds_reference_summary <- shiny::renderTable({
      bundle <- bundle_state()
      shiny::req(bundle)
      summary <- gripui.family.graph.summary(bundle$payload)
      note <- gripui.family.geometry.summary(bundle$payload)
      rbind(summary, note)
    }, striped = TRUE, bordered = FALSE, spacing = "s")

    output$gmds_reference_view <- shiny::renderUI({
      bundle <- bundle_state()
      shiny::req(bundle)
      widget <- gripui.gmds.render.widget(
        coords = bundle$payload$coords_display,
        edges = bundle$payload$edges,
        vertex_colors = rep("#8a5a44", bundle$payload$n),
        vertex_alpha = 0.95,
        edge_alpha = as.numeric(input$gmds_edge_alpha),
        base_size = 6
      )
      shiny::tagList(
        widget,
        shiny::tags$p(
          style = "padding:0 1rem 1rem;color:#5f5445;line-height:1.45;margin-bottom:0;",
          if (is.null(bundle$payload$note) || !nzchar(bundle$payload$note)) {
            "Reference graph geometry supplied by the selected synthetic family."
          } else {
            bundle$payload$note
          }
        )
      )
    })

    output$gmds_misf_title <- shiny::renderUI({
      bundle <- bundle_state()
      shiny::req(bundle)
      shiny::tags$span(sprintf("MIS filtration centered on V_%d", current_level()))
    })

    output$gmds_misf_view <- shiny::renderUI({
      bundle <- bundle_state()
      shiny::req(bundle)
      level <- current_level()
      prepared <- bundle$prepared
      partition <- gripui.gmds.layer_partition(prepared)
      palette <- gripui.gmds.layer.colors(prepared)
      colors <- unname(palette[as.character(partition)])
      active.vertices <- as.integer(prepared$misf$levels[[level + 1L]])
      widget <- gripui.gmds.render.widget(
        coords = bundle$payload$coords_display,
        edges = bundle$payload$edges,
        vertex_colors = colors,
        vertex_alpha = as.numeric(input$gmds_vertex_alpha),
        edge_alpha = as.numeric(input$gmds_edge_alpha),
        highlight_vertices = active.vertices,
        base_size = 5,
        highlight_size = 9,
        background_alpha = as.numeric(input$gmds_background_alpha)
      )
      note <- if (level < prepared$coarsest_level_level && identical(level, prepared$top_level_level) && prepared$top_level_level < prepared$coarsest_level_level) {
        sprintf(
          "Vertices are colored by the MIS filtration partition, and the highlighted set is V_%d. The solve starts here because the coarsest level V_%d has fewer than d + 1 vertices.",
          level,
          prepared$coarsest_level_level
        )
      } else {
        sprintf(
          "Vertices are colored by the MIS filtration partition, and the highlighted set is V_%d.",
          level
        )
      }
      shiny::tagList(
        widget,
        shiny::tags$p(
          style = "padding:0 1rem 1rem;color:#5f5445;line-height:1.45;margin-bottom:0;",
          note
        )
      )
    })

    output$gmds_misf_table <- shiny::renderTable({
      bundle <- bundle_state()
      shiny::req(bundle)
      gripui.gmds.level.table(bundle$prepared)
    }, striped = TRUE, bordered = FALSE, spacing = "s")

    output$gmds_selected_level_summary <- shiny::renderTable({
      bundle <- bundle_state()
      shiny::req(bundle)
      gripui.gmds.selected.level.table(bundle, current_level())
    }, striped = TRUE, bordered = FALSE, spacing = "s")

    output$gmds_seed_title <- shiny::renderUI({
      bundle <- bundle_state()
      shiny::req(bundle)
      payload <- gripui.gmds.stage.payload(bundle, "seed")
      shiny::tags$span(gripui.family.value_or_default(if (!is.null(payload)) payload$label else NULL, "Geometric seed"))
    })

    output$gmds_initial_title <- shiny::renderUI({
      bundle <- bundle_state()
      shiny::req(bundle)
      payload <- gripui.gmds.stage.payload(bundle, "initial_placement")
      shiny::tags$span(gripui.family.value_or_default(if (!is.null(payload)) payload$label else NULL, "Initial placement"))
    })

    output$gmds_top_level_title <- shiny::renderUI({
      bundle <- bundle_state()
      shiny::req(bundle)
      payload <- gripui.gmds.stage.payload(bundle, "top_level")
      shiny::tags$span(gripui.family.value_or_default(if (!is.null(payload)) payload$label else NULL, "Top-level GMDS solve"))
    })

    output$gmds_seed_view <- shiny::renderUI({
      bundle <- bundle_state()
      shiny::req(bundle)
      payload <- gripui.gmds.stage.payload(bundle, "seed")
      shiny::req(payload)
      coords <- gripui.gmds.stage.display.coords(bundle, payload, fill_reference = TRUE)
      widget <- gripui.gmds.render.widget(
        coords = coords,
        edges = bundle$payload$edges,
        vertex_colors = rep("#8a5a44", nrow(coords)),
        vertex_alpha = as.numeric(input$gmds_vertex_alpha),
        edge_alpha = as.numeric(input$gmds_edge_alpha),
        highlight_vertices = payload$active_vertices,
        base_size = 5,
        highlight_size = 10,
        background_alpha = as.numeric(input$gmds_background_alpha)
      )
      shiny::tagList(
        widget,
        shiny::tags$p(
          style = "padding:0 1rem 1rem;color:#5f5445;line-height:1.45;margin-bottom:0;",
          paste(
            gripui.gmds.stage.note(bundle, "seed"),
            "Inactive vertices stay at their aligned reference positions for context."
          )
        )
      )
    })

    output$gmds_initial_view <- shiny::renderUI({
      bundle <- bundle_state()
      shiny::req(bundle)
      payload <- gripui.gmds.stage.payload(bundle, "initial_placement")
      shiny::req(payload)
      coords <- gripui.gmds.stage.display.coords(bundle, payload, fill_reference = TRUE)
      widget <- gripui.gmds.render.widget(
        coords = coords,
        edges = bundle$payload$edges,
        vertex_colors = rep("#1f3b73", nrow(coords)),
        vertex_alpha = as.numeric(input$gmds_vertex_alpha),
        edge_alpha = as.numeric(input$gmds_edge_alpha),
        highlight_vertices = payload$active_vertices,
        base_size = 5,
        highlight_size = 9,
        background_alpha = as.numeric(input$gmds_background_alpha)
      )
      shiny::tagList(
        widget,
        shiny::tags$p(
          style = "padding:0 1rem 1rem;color:#5f5445;line-height:1.45;margin-bottom:0;",
          paste(
            gripui.gmds.stage.note(bundle, "initial_placement"),
            sprintf(
              "This is the warm-start state used by the top-level %s refinement.",
              bundle$method$refinement_phrase
            )
          )
        )
      )
    })

    output$gmds_top_level_view <- shiny::renderUI({
      bundle <- bundle_state()
      shiny::req(bundle)
      payload <- gripui.gmds.stage.payload(bundle, "top_level")
      shiny::req(payload)
      coords <- gripui.gmds.stage.display.coords(bundle, payload, fill_reference = TRUE)
      widget <- gripui.gmds.render.widget(
        coords = coords,
        edges = bundle$payload$edges,
        vertex_colors = rep("#206a5d", nrow(coords)),
        vertex_alpha = as.numeric(input$gmds_vertex_alpha),
        edge_alpha = as.numeric(input$gmds_edge_alpha),
        highlight_vertices = payload$active_vertices,
        base_size = 5,
        highlight_size = 9,
        background_alpha = as.numeric(input$gmds_background_alpha)
      )
      shiny::tagList(
        widget,
        shiny::tags$p(
          style = "padding:0 1rem 1rem;color:#5f5445;line-height:1.45;margin-bottom:0;",
          paste(
            gripui.gmds.stage.note(bundle, "top_level"),
            "The same active set is shown after objective-specific top-level refinement."
          )
        )
      )
    })

    output$gmds_seed_summary <- shiny::renderTable({
      bundle <- bundle_state()
      shiny::req(bundle)
      gripui.gmds.stage.summary.table(gripui.gmds.stage.payload(bundle, "seed"))
    }, striped = TRUE, bordered = FALSE, spacing = "s")

    output$gmds_initial_summary <- shiny::renderTable({
      bundle <- bundle_state()
      shiny::req(bundle)
      gripui.gmds.stage.summary.table(gripui.gmds.stage.payload(bundle, "initial_placement"))
    }, striped = TRUE, bordered = FALSE, spacing = "s")

    output$gmds_top_level_summary <- shiny::renderTable({
      bundle <- bundle_state()
      shiny::req(bundle)
      gripui.gmds.stage.summary.table(gripui.gmds.stage.payload(bundle, "top_level"))
    }, striped = TRUE, bordered = FALSE, spacing = "s")

    output$gmds_insertion_title <- shiny::renderUI({
      bundle <- bundle_state()
      shiny::req(bundle)
      payload <- gripui.gmds.level.stage.payload(bundle, "insertion", current_expansion_level())
      shiny::tags$span(gripui.family.value_or_default(if (!is.null(payload)) payload$label else NULL, sprintf("Insertion of V_%d", current_expansion_level())))
    })

    output$gmds_refinement_title <- shiny::renderUI({
      bundle <- bundle_state()
      shiny::req(bundle)
      payload <- gripui.gmds.level.stage.payload(bundle, "refinement", current_expansion_level())
      shiny::tags$span(gripui.family.value_or_default(if (!is.null(payload)) payload$label else NULL, sprintf("Refinement of V_%d", current_expansion_level())))
    })

    output$gmds_insertion_view <- shiny::renderUI({
      bundle <- bundle_state()
      shiny::req(bundle)
      level <- current_expansion_level()
      payload <- gripui.gmds.level.stage.payload(bundle, "insertion", level)
      shiny::req(payload)
      coords <- gripui.gmds.stage.display.coords(bundle, payload, fill_reference = TRUE)
      widget <- gripui.gmds.render.widget(
        coords = coords,
        edges = bundle$payload$edges,
        vertex_colors = rep("#1f3b73", nrow(coords)),
        vertex_alpha = as.numeric(input$gmds_vertex_alpha),
        edge_alpha = as.numeric(input$gmds_edge_alpha),
        highlight_vertices = payload$active_vertices,
        base_size = 5,
        highlight_size = 9,
        background_alpha = as.numeric(input$gmds_background_alpha)
      )
      inserted.n <- if (length(payload$inserted_vertices)) length(payload$inserted_vertices) else payload$summary$inserted_n
      shiny::tagList(
        widget,
        shiny::tags$p(
          style = "padding:0 1rem 1rem;color:#5f5445;line-height:1.45;margin-bottom:0;",
          paste(
            gripui.gmds.level.stage.note(bundle, "insertion", level, inserted_n = inserted.n),
            "The highlighted active set is the whole current level after insertion."
          )
        )
      )
    })

    output$gmds_refinement_view <- shiny::renderUI({
      bundle <- bundle_state()
      shiny::req(bundle)
      level <- current_expansion_level()
      payload <- gripui.gmds.level.stage.payload(bundle, "refinement", level)
      shiny::req(payload)
      coords <- gripui.gmds.stage.display.coords(bundle, payload, fill_reference = TRUE)
      widget <- gripui.gmds.render.widget(
        coords = coords,
        edges = bundle$payload$edges,
        vertex_colors = rep("#206a5d", nrow(coords)),
        vertex_alpha = as.numeric(input$gmds_vertex_alpha),
        edge_alpha = as.numeric(input$gmds_edge_alpha),
        highlight_vertices = payload$active_vertices,
        base_size = 5,
        highlight_size = 9,
        background_alpha = as.numeric(input$gmds_background_alpha)
      )
      shiny::tagList(
        widget,
        shiny::tags$p(
          style = "padding:0 1rem 1rem;color:#5f5445;line-height:1.45;margin-bottom:0;",
          paste(
            gripui.gmds.level.stage.note(bundle, "refinement", level),
            "The same active set is shown after objective-specific refinement of that level."
          )
        )
      )
    })

    output$gmds_insertion_summary <- shiny::renderTable({
      bundle <- bundle_state()
      shiny::req(bundle)
      gripui.gmds.stage.summary.table(gripui.gmds.level.stage.payload(bundle, "insertion", current_expansion_level()))
    }, striped = TRUE, bordered = FALSE, spacing = "s")

    output$gmds_refinement_summary <- shiny::renderTable({
      bundle <- bundle_state()
      shiny::req(bundle)
      gripui.gmds.stage.summary.table(gripui.gmds.level.stage.payload(bundle, "refinement", current_expansion_level()))
    }, striped = TRUE, bordered = FALSE, spacing = "s")

    output$gmds_stage_trace <- shiny::renderTable({
      bundle <- bundle_state()
      shiny::req(bundle)
      gripui.gmds.stage.trace.table(bundle)
    }, striped = TRUE, bordered = FALSE, spacing = "s")

    output$gmds_code <- shiny::renderText({
      bundle <- bundle_state()
      shiny::req(bundle)
      gripui.gmds.repro.code(bundle)
    })
  }
}

#' Build the GMDS stage explorer Shiny application
#'
#' The GMDS stage explorer reuses the synthetic family catalog from
#' [gripui_family_app()] and computes a canonical MISF-GMDS bundle for the
#' selected graph. The current implementation covers Milestones 1 through 4:
#' graph and geometry inspection, visualization of the MIS filtration,
#' explicit seed and initial-placement panels, per-level insertion and
#' refinement panels, method switching across canonical GMDS/GKK/LGKK traces,
#' and a canonical stage-trace summary that later panels will extend.
#'
#' @param catalog Family catalog, usually [gripui_graph_family_catalog()].
#' @param title Application title.
#' @param subtitle Optional subtitle shown in the sidebar header.
#'
#' @return A `shiny.appobj`.
#' @export
#'
#' @examplesIf local({ old <- getOption("rgl.useNULL"); options(rgl.useNULL = TRUE); on.exit(options(rgl.useNULL = old), add = TRUE); requireNamespace("shiny", quietly = TRUE) && requireNamespace("bslib", quietly = TRUE) && requireNamespace("rgl", quietly = TRUE) })
#' app <- gripui_gmds_app()
#' inherits(app, "shiny.appobj")
gripui_gmds_app <- function(catalog = gripui.gmds.default.catalog(),
                            title = "GMDS Stage Explorer",
                            subtitle = "Milestones 1-4: graph and geometry selection, MIS filtration, seed selection, insertion/refinement, and method switching across canonical GMDS/GKK/LGKK stage traces.") {
  catalog <- gripui.gmds.default.catalog(catalog)
  old <- gripui.require.family.app.packages()
  on.exit(options(rgl.useNULL = old), add = TRUE)

  shiny::shinyApp(
    ui = gripui.gmds.ui(catalog = catalog, title = title, subtitle = subtitle),
    server = gripui.gmds.server(catalog = catalog)
  )
}

#' Run the GMDS stage explorer Shiny application
#'
#' @param catalog Family catalog, usually [gripui_graph_family_catalog()].
#' @param title Application title.
#' @param subtitle Optional subtitle shown in the sidebar header.
#' @param host Host passed to [shiny::runApp()].
#' @param port Port passed to [shiny::runApp()].
#' @param launch.browser Whether to launch a browser.
#' @param auto.stop.after Optional delay, in seconds, after which the app stops
#'   itself. This is mainly useful for automated examples and tests.
#' @param ... Additional arguments passed to [shiny::runApp()].
#'
#' @return Invisibly returns the result of [shiny::runApp()].
#' @export
#'
#' @examplesIf local({ old <- getOption("rgl.useNULL"); options(rgl.useNULL = TRUE); on.exit(options(rgl.useNULL = old), add = TRUE); requireNamespace("shiny", quietly = TRUE) && requireNamespace("bslib", quietly = TRUE) && requireNamespace("rgl", quietly = TRUE) && requireNamespace("later", quietly = TRUE) })
#' run_gripui_gmds(launch.browser = FALSE, quiet = TRUE, auto.stop.after = 0.1)
run_gripui_gmds <- function(catalog = gripui.gmds.default.catalog(),
                            title = "GMDS Stage Explorer",
                            subtitle = "Milestones 1-4: graph and geometry selection, MIS filtration, seed selection, insertion/refinement, and method switching across canonical GMDS/GKK/LGKK stage traces.",
                            host = "127.0.0.1",
                            port = getOption("shiny.port"),
                            launch.browser = interactive(),
                            auto.stop.after = NULL,
                            ...) {
  app <- gripui_gmds_app(catalog = catalog, title = title, subtitle = subtitle)

  if (!is.null(auto.stop.after)) {
    if (!requireNamespace("later", quietly = TRUE)) {
      stop("Package 'later' is required when auto.stop.after is used.", call. = FALSE)
    }
    if (!is.numeric(auto.stop.after) ||
        length(auto.stop.after) != 1L ||
        is.na(auto.stop.after) ||
        auto.stop.after < 0) {
      stop("auto.stop.after must be a single non-negative number of seconds.", call. = FALSE)
    }

    previous.on.start <- app$onStart
    app$onStart <- function() {
      on.stop <- if (is.function(previous.on.start)) previous.on.start() else NULL
      later::later(function() shiny::stopApp(invisible(NULL)), delay = auto.stop.after)
      on.stop
    }
  }

  shiny::runApp(
    app,
    host = host,
    port = port,
    launch.browser = launch.browser,
    ...
  )
}
