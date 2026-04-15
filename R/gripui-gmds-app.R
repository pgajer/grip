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

gripui.gmds.valid.levels <- function(prepared) {
  seq.int(from = 0L, to = length(prepared$misf$levels) - 1L, by = 1L)
}

gripui.gmds.normalize.level <- function(prepared,
                                        level = NULL,
                                        default = prepared$top_level_level) {
  valid <- gripui.gmds.valid.levels(prepared)
  selected <- suppressWarnings(as.integer(level[[1L]]))
  fallback <- suppressWarnings(as.integer(default[[1L]]))
  if (!length(fallback) || !is.finite(fallback) || !fallback %in% valid) {
    fallback <- valid[[length(valid)]]
  }
  if (!length(selected) || !is.finite(selected) || !selected %in% valid) {
    return(fallback)
  }
  selected
}

gripui.gmds.method.catalog <- function() {
  list(
    grip = list(
      id = "grip",
      label = "GRIP",
      objective_label = "Classical GRIP multiscale layout",
      refinement_phrase = "GRIP local-KK",
      prepare_fn = "grip.build.misf",
      optimize_fn = "grip.layout.trace",
      fit_class = "grip_misf_grip_fit"
    ),
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
  level <- gripui.gmds.normalize.level(prepared, level = level)
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
  if (identical(bundle$method$id, "grip")) {
    return(switch(
      state,
      seed = sprintf(
        "GRIP does not expose a separate geometric d+1 seed. This panel mirrors the first recorded coarsest-level state on V_%d.",
        prepared$top_level_level
      ),
      top_level = sprintf(
        "This panel shows the coarsest active set V_%d after the recorded GRIP multiscale rounds on that level.",
        prepared$top_level_level
      ),
      initial_placement = sprintf(
        "The compiled GRIP trace starts at the coarsest active set V_%d, so the app reuses that first recorded state as the initial-placement panel.",
        prepared$top_level_level
      ),
      ""
    ))
  }
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
    level = if (!is.null(payload$level) && length(payload$level) == 1L && is.finite(payload$level)) {
      sprintf("V_%d", payload$level)
    } else {
      ""
    },
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

gripui.gmds.paper.sync.table <- function() {
  data.frame(
    app_panel = c(
      "Stage 0: Graph and Geometry",
      "Stage 1: MIS Filtration",
      "Stage 2: Seed Selection and Initial Placement",
      "Stage 3: Per-Level Insertion and Refinement",
      "Canonical Trace Summary"
    ),
    manuscript_section = c(
      "Synthetic graph family setup and experiment descriptions",
      "Filtration",
      "Coarsest seed, expansion, and refinement",
      "Initial placement and refinement of V_l",
      "Multiscale case studies and recent experiments"
    ),
    figures_tables = c(
      "Case-study setup tables and geometry figures",
      "Figure 7 (MIS filtration cases)",
      "Top-level stage figures and interactive HTML panels",
      "Per-level stage figures and tables",
      "Stage trace tables and experiment discussion"
    ),
    stringsAsFactors = FALSE
  )
}

gripui.gmds.manuscript.source.path <- function() {
  "dev/papers/geodesic_mds_paper/manuscript/geodesic_mds.tex"
}

gripui.gmds.export.preset.choices <- function() {
  stats::setNames(
    c("paper_figure_bundle", "audit_bundle", "tables_only"),
    c("Paper figure bundle", "Interactive audit bundle", "Tables only")
  )
}

gripui.gmds.figure.preset.catalog <- function() {
  list(
    paper_panel = list(
      id = "paper_panel",
      label = "Paper panel",
      width_in = 6.6,
      height_in = 4.8,
      png_width = 1980L,
      png_height = 1440L,
      res = 300L,
      edge_alpha = 0.42,
      highlight_edge_alpha = 0.92,
      background_alpha = 0.24,
      vertex_alpha = 0.90,
      highlight_alpha = 0.98,
      base_cex = 0.78,
      highlight_cex = 1.38,
      edge_lwd = 1.10,
      highlight_edge_lwd = 2.35
    ),
    paper_wide = list(
      id = "paper_wide",
      label = "Paper wide",
      width_in = 7.2,
      height_in = 4.2,
      png_width = 2160L,
      png_height = 1260L,
      res = 300L,
      edge_alpha = 0.38,
      highlight_edge_alpha = 0.88,
      background_alpha = 0.22,
      vertex_alpha = 0.88,
      highlight_alpha = 0.98,
      base_cex = 0.74,
      highlight_cex = 1.30,
      edge_lwd = 1.05,
      highlight_edge_lwd = 2.20
    )
  )
}

gripui.gmds.figure.preset.spec <- function(preset = "paper_panel") {
  preset <- as.character(preset[[1L]])
  catalog <- gripui.gmds.figure.preset.catalog()
  if (!preset %in% names(catalog)) {
    stop("unknown GMDS figure preset: ", preset)
  }
  catalog[[preset]]
}

gripui.gmds.figure.preset.choices <- function() {
  catalog <- gripui.gmds.figure.preset.catalog()
  stats::setNames(
    vapply(catalog, `[[`, character(1L), "id"),
    vapply(catalog, `[[`, character(1L), "label")
  )
}

gripui.gmds.static.project <- function(coords) {
  coords <- gripui.family.pad_coords(coords, target_cols = 3L)
  keep <- gripui.valid.coord.rows(coords)
  if (!any(keep)) {
    return(list(
      projected = matrix(NA_real_, nrow = nrow(coords), ncol = 2L),
      depth = rep(NA_real_, nrow(coords)),
      keep = keep
    ))
  }
  centered <- coords
  centered[keep, ] <- scale(centered[keep, , drop = FALSE], center = TRUE, scale = FALSE)
  theta <- 35 * pi / 180
  phi <- 22 * pi / 180
  rot.z <- matrix(c(
    cos(theta), -sin(theta), 0,
    sin(theta), cos(theta), 0,
    0, 0, 1
  ), nrow = 3L, byrow = TRUE)
  rot.x <- matrix(c(
    1, 0, 0,
    0, cos(phi), -sin(phi),
    0, sin(phi), cos(phi)
  ), nrow = 3L, byrow = TRUE)
  rotated <- centered %*% rot.z %*% rot.x
  out <- matrix(NA_real_, nrow = nrow(coords), ncol = 2L)
  out[keep, ] <- rotated[keep, 1:2, drop = FALSE]
  list(
    projected = out,
    depth = rotated[, 3L],
    keep = keep
  )
}

gripui.gmds.static.edge.order <- function(edges, depth, keep) {
  if (!nrow(edges)) {
    return(integer())
  }
  valid <- keep[edges[, 1L]] & keep[edges[, 2L]]
  if (!any(valid)) {
    return(integer())
  }
  edges <- edges[valid, , drop = FALSE]
  edge.depth <- rowMeans(cbind(depth[edges[, 1L]], depth[edges[, 2L]]))
  order(edge.depth, decreasing = FALSE)
}

gripui.gmds.draw.static.figure <- function(export_payload,
                                           figure_preset = "paper_panel") {
  spec <- gripui.gmds.figure.preset.spec(figure_preset)
  coords <- if (!is.null(export_payload$display_coords)) {
    export_payload$display_coords
  } else {
    export_payload$coords
  }
  coords <- gripui.family.pad_coords(coords, target_cols = 3L)
  projection <- gripui.gmds.static.project(coords)
  proj <- projection$projected
  keep <- projection$keep
  depth <- projection$depth

  x.range <- range(proj[keep, 1L], finite = TRUE)
  y.range <- range(proj[keep, 2L], finite = TRUE)
  x.pad <- max(0.05, 0.06 * diff(x.range))
  y.pad <- max(0.05, 0.06 * diff(y.range))

  graphics::par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  graphics::plot.new()
  graphics::plot.window(
    xlim = c(x.range[[1L]] - x.pad, x.range[[2L]] + x.pad),
    ylim = c(y.range[[1L]] - y.pad, y.range[[2L]] + y.pad),
    asp = 1
  )

  edges <- as.matrix(export_payload$edges)
  storage.mode(edges) <- "integer"
  if (nrow(edges)) {
    edge.order <- gripui.gmds.static.edge.order(edges, depth, keep)
    if (length(edge.order)) {
      edge.keep <- keep[edges[, 1L]] & keep[edges[, 2L]]
      edges.valid <- edges[edge.keep, , drop = FALSE]
      edges.valid <- edges.valid[edge.order, , drop = FALSE]
      highlight.vertices <- if (!is.null(export_payload$highlight_vertices)) {
        as.integer(export_payload$highlight_vertices)
      } else if (!is.null(export_payload$active_vertices)) {
        as.integer(export_payload$active_vertices)
      } else {
        integer(0L)
      }
      edge.base.col <- grDevices::adjustcolor("#9ca3af", alpha.f = spec$edge_alpha)
      edge.highlight.col <- grDevices::adjustcolor("#8a5a44", alpha.f = spec$highlight_edge_alpha)
      edge.highlight <- edges.valid[, 1L] %in% highlight.vertices & edges.valid[, 2L] %in% highlight.vertices
      graphics::segments(
        x0 = proj[edges.valid[, 1L], 1L],
        y0 = proj[edges.valid[, 1L], 2L],
        x1 = proj[edges.valid[, 2L], 1L],
        y1 = proj[edges.valid[, 2L], 2L],
        col = ifelse(edge.highlight, edge.highlight.col, edge.base.col),
        lwd = ifelse(edge.highlight, spec$highlight_edge_lwd, spec$edge_lwd)
      )
    }
  }

  point.order <- order(depth[keep], decreasing = FALSE)
  vertices <- which(keep)[point.order]
  point.cols <- grDevices::adjustcolor(export_payload$vertex_colors[vertices], alpha.f = spec$background_alpha)
  highlight.vertices <- if (!is.null(export_payload$highlight_vertices)) {
    as.integer(export_payload$highlight_vertices)
  } else if (!is.null(export_payload$active_vertices)) {
    as.integer(export_payload$active_vertices)
  } else {
    integer(0L)
  }
  is.highlight <- vertices %in% highlight.vertices
  point.cols[is.highlight] <- grDevices::adjustcolor(export_payload$vertex_colors[vertices[is.highlight]], alpha.f = spec$highlight_alpha)
  point.cex <- rep(spec$base_cex, length(vertices))
  point.cex[is.highlight] <- spec$highlight_cex
  graphics::points(
    proj[vertices, 1L],
    proj[vertices, 2L],
    pch = 16,
    col = point.cols,
    cex = point.cex
  )
}

gripui.gmds.write.static.figure <- function(export_payload,
                                            png_path = NULL,
                                            pdf_path = NULL,
                                            figure_preset = "paper_panel") {
  spec <- gripui.gmds.figure.preset.spec(figure_preset)
  render_device <- function(open_device) {
    open_device()
    tryCatch(
      gripui.gmds.draw.static.figure(export_payload, figure_preset = figure_preset),
      finally = grDevices::dev.off()
    )
  }
  if (!is.null(pdf_path)) {
    render_device(function() {
      grDevices::pdf(pdf_path, width = spec$width_in, height = spec$height_in, useDingbats = FALSE)
    })
  }
  if (!is.null(png_path)) {
    render_device(function() {
      grDevices::png(
        filename = png_path,
        width = spec$png_width,
        height = spec$png_height,
        res = spec$res,
        bg = "#ffffff"
      )
    })
  }
  invisible(list(
    png_path = png_path,
    pdf_path = pdf_path,
    figure_preset = figure_preset
  ))
}

gripui.gmds.paper.panel.catalog <- function() {
  data.frame(
    stage_id = c(
      "reference",
      "misf",
      "seed",
      "initial_placement",
      "top_level",
      "insertion",
      "refinement",
      "final_polish",
      "trace_summary"
    ),
    app_panel = c(
      "Stage 0: Graph and Geometry",
      "Stage 1: MIS Filtration",
      "Stage 2: Seed Selection and Initial Placement",
      "Stage 2: Seed Selection and Initial Placement",
      "Stage 2: Seed Selection and Initial Placement",
      "Stage 3: Per-Level Insertion and Refinement",
      "Stage 3: Per-Level Insertion and Refinement",
      "Canonical Trace Summary",
      "Canonical Trace Summary"
    ),
    manuscript_section = c(
      "Synthetic graph family setup and experiment descriptions",
      "Filtration",
      "Coarsest seed, expansion, and refinement",
      "Coarsest seed, expansion, and refinement",
      "Coarsest seed, expansion, and refinement",
      "Initial placement and refinement of V_l",
      "Initial placement and refinement of V_l",
      "Multiscale case studies and recent experiments",
      "Multiscale case studies and recent experiments"
    ),
    figure_focus = c(
      "Case-study setup tables and geometry figures",
      "Figure 7 (MIS filtration cases)",
      "Top-level stage figures and interactive HTML panels",
      "Top-level stage figures and interactive HTML panels",
      "Top-level stage figures and interactive HTML panels",
      "Per-level stage figures and tables",
      "Per-level stage figures and tables",
      "Final comparison figures and result tables",
      "Stage trace tables and experiment discussion"
    ),
    narrative_focus = c(
      "Describe the synthetic graph family, geometry, and experiment controls.",
      "Define the MIS filtration and interpret the highlighted level sizes and separations.",
      "Explain how the admissible top level and the seed are selected.",
      "Describe how the seed expands into the initial placement before refinement.",
      "Discuss the objective-specific top-level solve and its qualitative outcome.",
      "Describe how new vertices are inserted into the inherited scaffold at the selected level.",
      "Discuss the level-specific refinement and how the objective changes the embedding.",
      "Summarize final embedding quality and how the selected method behaves on the full graph.",
      "Use the canonical stage table to connect the app trace to manuscript figures and tables."
    ),
    stringsAsFactors = FALSE
  )
}

gripui.gmds.paper.context <- function(bundle,
                                      stage_id,
                                      focus_level = NULL,
                                      expansion_level = NULL) {
  catalog <- gripui.gmds.paper.panel.catalog()
  stage_id <- as.character(stage_id[[1L]])
  row <- catalog[catalog$stage_id == stage_id, , drop = FALSE]
  if (!nrow(row)) {
    row <- data.frame(
      stage_id = stage_id,
      app_panel = "Unmapped panel",
      manuscript_section = "",
      figure_focus = "",
      narrative_focus = "",
      stringsAsFactors = FALSE
    )
  }
  level <- switch(
    stage_id,
    reference = NA_integer_,
    trace_summary = NA_integer_,
    misf = as.integer(focus_level[[1L]]),
    seed = as.integer(bundle$prepared$top_level_level),
    initial_placement = as.integer(bundle$prepared$top_level_level),
    top_level = as.integer(bundle$prepared$top_level_level),
    insertion = as.integer(expansion_level[[1L]]),
    refinement = as.integer(expansion_level[[1L]]),
    final_polish = 0L,
    as.integer(focus_level[[1L]])
  )
  data.frame(
    method = bundle$method$label,
    family = bundle$payload$family_label,
    export_stage = stage_id,
    level = level,
    app_panel = row$app_panel[[1L]],
    manuscript_section = row$manuscript_section[[1L]],
    figure_focus = row$figure_focus[[1L]],
    narrative_focus = row$narrative_focus[[1L]],
    manuscript_source = gripui.gmds.manuscript.source.path(),
    stringsAsFactors = FALSE
  )
}

gripui.gmds.paper.context.table <- function(bundle,
                                            stage_id,
                                            focus_level = NULL,
                                            expansion_level = NULL) {
  context <- gripui.gmds.paper.context(
    bundle = bundle,
    stage_id = stage_id,
    focus_level = focus_level,
    expansion_level = expansion_level
  )
  gripui.family.summary.table(list(
    method = context$method,
    family = context$family,
    stage = context$export_stage,
    level = if (is.finite(context$level)) sprintf("V_%d", context$level) else "",
    manuscript_section = context$manuscript_section,
    figure_focus = context$figure_focus,
    manuscript_source = context$manuscript_source
  ))
}

gripui.gmds.paper.inline.note <- function(bundle,
                                          stage_id,
                                          focus_level = NULL,
                                          expansion_level = NULL) {
  context <- gripui.gmds.paper.context(
    bundle = bundle,
    stage_id = stage_id,
    focus_level = focus_level,
    expansion_level = expansion_level
  )
  sprintf(
    "Paper link: %s. Figure focus: %s.",
    context$manuscript_section[[1L]],
    context$figure_focus[[1L]]
  )
}

gripui.gmds.export.stage.slug <- function(label) {
  label <- tolower(as.character(label[[1L]]))
  label <- gsub("[^A-Za-z0-9]+", "_", label)
  gsub("^_+|_+$", "", label)
}

gripui.gmds.paper.caption.template <- function(bundle,
                                               export_payload,
                                               stage_id,
                                               focus_level = NULL,
                                               expansion_level = NULL) {
  context <- gripui.gmds.paper.context(
    bundle = bundle,
    stage_id = stage_id,
    focus_level = focus_level,
    expansion_level = expansion_level
  )
  level.text <- if (!is.null(export_payload$level) &&
    length(export_payload$level) == 1L &&
    is.finite(export_payload$level)) {
    sprintf(" at V_%d", export_payload$level)
  } else {
    ""
  }
  sprintf(
    "%s on the %s family using %s%s. Suggested manuscript placement: %s.",
    export_payload$label,
    bundle$payload$family_label,
    bundle$method$label,
    level.text,
    context$manuscript_section[[1L]]
  )
}

gripui.gmds.paper.note.lines <- function(bundle,
                                         export_payload,
                                         stage_id,
                                         focus_level = NULL,
                                         expansion_level = NULL,
                                         preset = "paper_figure_bundle",
                                         figure_preset = "paper_panel") {
  context <- gripui.gmds.paper.context(
    bundle = bundle,
    stage_id = stage_id,
    focus_level = focus_level,
    expansion_level = expansion_level
  )
  c(
    sprintf("# GMDS %s export", preset),
    "",
    sprintf("- Method: %s", bundle$method$label),
    sprintf("- Family: %s", bundle$payload$family_label),
    sprintf("- Stage: %s", export_payload$label),
    sprintf("- Figure preset: %s", figure_preset),
    sprintf("- Manuscript source: `%s`", context$manuscript_source[[1L]]),
    sprintf("- Suggested section: %s", context$manuscript_section[[1L]]),
    sprintf("- Figure focus: %s", context$figure_focus[[1L]]),
    sprintf("- Narrative focus: %s", context$narrative_focus[[1L]]),
    "",
    "## Suggested caption",
    "",
    gripui.gmds.paper.caption.template(
      bundle = bundle,
      export_payload = export_payload,
      stage_id = stage_id,
      focus_level = focus_level,
      expansion_level = expansion_level
    )
  )
}

gripui.gmds.write.export.bundle <- function(bundle,
                                            export_payload,
                                            stage_id,
                                            focus_level = NULL,
                                            expansion_level = NULL,
                                            preset = "paper_figure_bundle",
                                            figure_preset = "paper_panel",
                                            dir) {
  preset <- as.character(preset[[1L]])
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  stage.trace.path <- file.path(dir, "stage_trace.csv")
  utils::write.csv(bundle$stage_trace, stage.trace.path, row.names = FALSE)

  vertices.path <- file.path(dir, "selected_stage_vertices.csv")
  utils::write.csv(
    gripui.gmds.export.vertex.table(bundle, export_payload),
    vertices.path,
    row.names = FALSE
  )

  summary.path <- file.path(dir, "selected_stage_summary.csv")
  utils::write.csv(
    gripui.gmds.stage.summary.table(export_payload),
    summary.path,
    row.names = FALSE
  )

  context.path <- file.path(dir, "paper_context.csv")
  utils::write.csv(
    gripui.gmds.paper.context(
      bundle = bundle,
      stage_id = stage_id,
      focus_level = focus_level,
      expansion_level = expansion_level
    ),
    context.path,
    row.names = FALSE
  )

  note.path <- file.path(dir, "paper_note.md")
  writeLines(
    gripui.gmds.paper.note.lines(
      bundle = bundle,
      export_payload = export_payload,
      stage_id = stage_id,
      focus_level = focus_level,
      expansion_level = expansion_level,
      preset = preset,
      figure_preset = figure_preset
    ),
    con = note.path
  )

  caption.path <- file.path(dir, "suggested_caption.txt")
  writeLines(
    gripui.gmds.paper.caption.template(
      bundle = bundle,
      export_payload = export_payload,
      stage_id = stage_id,
      focus_level = focus_level,
      expansion_level = expansion_level
    ),
    con = caption.path
  )

  repro.path <- file.path(dir, "repro.R")
  writeLines(gripui.gmds.repro.code(bundle), con = repro.path)

  manifest.path <- file.path(dir, "bundle_manifest.txt")
  writeLines(c(
    sprintf("preset=%s", preset),
    sprintf("method=%s", bundle$method$id),
    sprintf("family=%s", bundle$payload$family_id),
    sprintf("stage_id=%s", stage_id),
    sprintf("stage_label=%s", export_payload$label),
    sprintf("figure_preset=%s", figure_preset)
  ), con = manifest.path)

  written <- c(
    stage.trace.path,
    vertices.path,
    summary.path,
    context.path,
    note.path,
    caption.path,
    repro.path,
    manifest.path
  )

  if (preset %in% c("paper_figure_bundle", "audit_bundle")) {
    bundle.summary.path <- file.path(dir, "bundle_summary.csv")
    utils::write.csv(gripui.gmds.bundle.summary(bundle), bundle.summary.path, row.names = FALSE)
    paper.sync.path <- file.path(dir, "paper_sync_map.csv")
    utils::write.csv(gripui.gmds.paper.sync.table(), paper.sync.path, row.names = FALSE)
    written <- c(written, bundle.summary.path, paper.sync.path)
  }

  if (identical(preset, "audit_bundle")) {
    misf.levels.path <- file.path(dir, "misf_levels.csv")
    utils::write.csv(gripui.gmds.level.table(bundle$prepared), misf.levels.path, row.names = FALSE)
    written <- c(written, misf.levels.path)
  }

  if (!identical(preset, "tables_only")) {
    static.png.path <- file.path(dir, sprintf("selected_stage_%s.png", figure_preset))
    static.pdf.path <- file.path(dir, sprintf("selected_stage_%s.pdf", figure_preset))
    gripui.gmds.write.static.figure(
      export_payload = export_payload,
      png_path = static.png.path,
      pdf_path = static.pdf.path,
      figure_preset = figure_preset
    )
    written <- c(written, static.png.path, static.pdf.path)
  }

  if (!identical(preset, "tables_only") && requireNamespace("htmlwidgets", quietly = TRUE)) {
    snapshot.path <- file.path(dir, "selected_stage_snapshot.html")
    widget <- gripui.gmds.export.widget(
      bundle = bundle,
      export_payload = export_payload
    )
    htmlwidgets::saveWidget(widget, file = snapshot.path, selfcontained = TRUE)
    written <- c(written, snapshot.path)
  }

  figure.spec.path <- file.path(dir, "figure_preset.csv")
  utils::write.csv(as.data.frame(gripui.gmds.figure.preset.spec(figure_preset), stringsAsFactors = FALSE), figure.spec.path, row.names = FALSE)
  written <- c(written, figure.spec.path)

  invisible(normalizePath(written, winslash = "/", mustWork = FALSE))
}

gripui.gmds.export.stage.choices <- function(bundle) {
  base <- c(
    reference = "Reference graph geometry",
    misf = "Selected MIS filtration level",
    seed = "Seed stage",
    initial_placement = "Initial placement",
    top_level = "Top-level solve",
    final_polish = "Final full-graph polish"
  )
  if (!is.null(bundle) && length(gripui.gmds.expansion.levels(bundle))) {
    base <- c(
      base[seq_len(5L)],
      insertion = "Selected insertion level",
      refinement = "Selected refinement level",
      final_polish = base[["final_polish"]]
    )
  }
  stats::setNames(names(base), unname(base))
}

gripui.gmds.export.stage.payload <- function(bundle,
                                             stage_id,
                                             focus_level = NULL,
                                             expansion_level = NULL) {
  if (is.null(bundle)) {
    return(NULL)
  }
  stage_id <- as.character(stage_id[[1L]])
  if (identical(stage_id, "reference")) {
    coords <- gripui.family.pad_coords(bundle$payload$coords_display, target_cols = 3L)
    return(list(
      kind = "reference",
      stage = "reference",
      label = "Reference graph geometry",
      level = NA_integer_,
      coords = coords,
      display_coords = coords,
      edges = bundle$payload$edges,
      vertex_colors = rep("#8a5a44", nrow(coords)),
      active_vertices = seq_len(nrow(coords)),
      highlight_vertices = integer(0L),
      pair_mode = NA_character_,
      summary = list(
        active_n = bundle$payload$n,
        pair_n = gripui.family.count.edges(bundle$payload$edges)
      )
    ))
  }
  if (identical(stage_id, "misf")) {
    prepared <- bundle$prepared
    level <- gripui.gmds.normalize.level(prepared, level = focus_level)
    partition <- gripui.gmds.layer_partition(prepared)
    palette <- gripui.gmds.layer.colors(prepared)
    colors <- unname(palette[as.character(partition)])
    vertices <- as.integer(prepared$misf$levels[[level + 1L]])
    coords <- gripui.family.pad_coords(bundle$payload$coords_display, target_cols = 3L)
    return(list(
      kind = "misf",
      stage = "misf",
      label = sprintf("MIS filtration level V_%d", level),
      level = level,
      coords = coords,
      display_coords = coords,
      edges = bundle$payload$edges,
      vertex_colors = colors,
      active_vertices = vertices,
      highlight_vertices = vertices,
      pair_mode = NA_character_,
      summary = list(
        active_n = length(vertices),
        pair_n = nrow(gripui.gmds.filter.edge.matrix(bundle$payload$edges, vertices))
      )
    ))
  }
  if (stage_id %in% c("seed", "initial_placement", "top_level", "final_polish")) {
    payload <- gripui.gmds.stage.payload(bundle, stage_id)
  } else if (stage_id %in% c("insertion", "refinement")) {
    payload <- gripui.gmds.level.stage.payload(bundle, stage_id, as.integer(expansion_level[[1L]]))
  } else {
    payload <- NULL
  }
  if (is.null(payload)) {
    return(NULL)
  }
  colors <- switch(
    stage_id,
    seed = rep("#8a5a44", nrow(payload$display_coords)),
    initial_placement = rep("#1f3b73", nrow(payload$display_coords)),
    top_level = rep("#206a5d", nrow(payload$display_coords)),
    insertion = rep("#1f3b73", nrow(payload$display_coords)),
    refinement = rep("#206a5d", nrow(payload$display_coords)),
    final_polish = rep("#206a5d", nrow(payload$display_coords)),
    rep("#8a5a44", nrow(payload$display_coords))
  )
  payload$kind <- "stage"
  payload$edges <- if (!is.null(payload$edge_matrix)) {
    payload$edge_matrix
  } else {
    bundle$payload$edges
  }
  payload$full_edges <- bundle$payload$edges
  payload$vertex_colors <- colors
  payload
}

gripui.gmds.export.vertex.table <- function(bundle, export_payload) {
  coords <- if (!is.null(export_payload$display_coords)) {
    export_payload$display_coords
  } else {
    export_payload$coords
  }
  coords <- gripui.family.pad_coords(coords, target_cols = 3L)
  active <- rep(FALSE, nrow(coords))
  highlight <- rep(FALSE, nrow(coords))
  if (!is.null(export_payload$active_vertices)) {
    active[as.integer(export_payload$active_vertices)] <- TRUE
  }
  if (!is.null(export_payload$highlight_vertices)) {
    highlight[as.integer(export_payload$highlight_vertices)] <- TRUE
  } else if (!is.null(export_payload$active_vertices)) {
    highlight[as.integer(export_payload$active_vertices)] <- TRUE
  }
  data.frame(
    vertex_id = seq_len(nrow(coords)),
    x = coords[, 1L],
    y = coords[, 2L],
    z = coords[, 3L],
    active = active,
    highlighted = highlight,
    stringsAsFactors = FALSE
  )
}

gripui.gmds.export.widget <- function(bundle,
                                      export_payload,
                                      vertex_alpha = 0.95,
                                      edge_alpha = 0.20,
                                      background_alpha = 0.35) {
  coords <- if (!is.null(export_payload$display_coords)) {
    export_payload$display_coords
  } else {
    export_payload$coords
  }
  gripui.gmds.render.widget(
    coords = coords,
    edges = export_payload$edges,
    vertex_colors = export_payload$vertex_colors,
    vertex_alpha = vertex_alpha,
    edge_alpha = edge_alpha,
    highlight_vertices = if (!is.null(export_payload$highlight_vertices)) export_payload$highlight_vertices else export_payload$active_vertices,
    base_size = 5,
    highlight_size = 9,
    background_alpha = background_alpha
  )
}

gripui.gmds.repro.code <- function(bundle) {
  payload <- bundle$payload
  prepared <- bundle$prepared
  settings <- bundle$settings
  if (identical(bundle$method$id, "grip")) {
    return(paste(
      payload$code,
      "",
      sprintf(
        "misf <- grip.build.misf(edges = graph$edges, n = graph$n, edge_weights = graph$edge_weights, num_init = %dL, num_nbrs = %dL, seed = %dL)",
        settings$num_init,
        settings$grip_num_nbrs,
        settings$grip_seed
      ),
      sprintf(
        paste0(
          "tr <- grip.layout.trace(edges = graph$edges, n = graph$n, dim = %dL, ",
          "placement = \"barycenter\", rounds = %dL, final_rounds = %dL, ",
          "num_init = %dL, num_nbrs = %dL, trace = \"round\", trace.every = 1L, seed = %dL)"
        ),
        settings$dim,
        settings$grip_rounds,
        settings$grip_final_rounds,
        settings$num_init,
        settings$grip_num_nbrs,
        settings$grip_seed
      ),
      sprintf(
        "# top solve level: V_%d (reason: %s)",
        prepared$top_level_level,
        prepared$top_level_selection_reason
      ),
      "# note: the app uses one GRIP seed for both MIS filtration and the traced layout.",
      sep = "\n"
    ))
  }
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

  if (identical(method.spec$id, "grip")) {
    grip.seed <- prepare_seed
    grip.rounds <- max(1L, top_level_max_iter, refinement_max_iter)
    grip.final.rounds <- max(1L, final_polish_max_iter)
    grip.num.nbrs <- max(20L, refinement_local_nbrs)
    refinement_landmark_count <- NA_integer_

    prepared <- grip.new.misf.grip.prepared(list(
      n = payload$n,
      edges = payload$edges,
      edge_weights = payload$edge_weights,
      misf = grip.build.misf(
        edges = payload$edges,
        n = payload$n,
        edge_weights = payload$edge_weights,
        num_init = num_init,
        num_nbrs = grip.num.nbrs,
        seed = grip.seed
      )
    ))
    prepared$level_vertices <- prepared$misf$levels
    prepared$active_levels <- prepared$misf$levels
    prepared$insertion_order <- prepared$misf$mish_order
    prepared$coarsest_level_index <- length(prepared$misf$levels)
    prepared$coarsest_level_level <- as.integer(length(prepared$misf$levels) - 1L)
    prepared$top_level_index <- prepared$coarsest_level_index
    prepared$top_level_level <- prepared$coarsest_level_level
    prepared$top_level_vertices <- as.integer(prepared$misf$levels[[prepared$top_level_index]])
    prepared$top_level_selection_reason <- "GRIP starts from the coarsest filtration level."
    prepared$top_level_dim <- dim

    elapsed <- system.time({
      raw.trace <- grip.layout.trace(
        edges = payload$edges,
        n = payload$n,
        edge_weights = payload$edge_weights,
        dim = dim,
        placement = "barycenter",
        rounds = grip.rounds,
        final_rounds = grip.final.rounds,
        num_init = num_init,
        num_nbrs = grip.num.nbrs,
        trace = "round",
        trace.every = 1L,
        seed = grip.seed
      )
    })[["elapsed"]]

    fit <- grip.new.misf.grip.fit(list(
      coords = raw.trace$final,
      prepared = prepared,
      raw_trace = raw.trace,
      stage_trace = raw.trace$stage_trace,
      stage_data = raw.trace$stage_data,
      trace = list(
        trace_schema_version = grip.geodesic.misf.stage.trace.schema.version(),
        stage_data = raw.trace$stage_data,
        trace_mode = raw.trace$trace,
        meta = raw.trace$meta,
        diagnostics = raw.trace$diagnostics
      ),
      frames = list(
        stage_data = raw.trace$stage_data,
        raw_frames = raw.trace$frames,
        final = raw.trace$final
      ),
      timing = list(
        top_level = NA_real_,
        insertion = NA_real_,
        refinement = NA_real_,
        final_polish = NA_real_,
        total = as.double(elapsed)
      )
    ))
    pair.mode <- NA_character_
    pair.full.limit <- NA_integer_
  } else if (identical(method.spec$id, "gmds")) {
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
      grip_seed = if (identical(method.spec$id, "grip")) grip.seed else NA_integer_,
      grip_rounds = if (identical(method.spec$id, "grip")) grip.rounds else NA_integer_,
      grip_final_rounds = if (identical(method.spec$id, "grip")) grip.final.rounds else NA_integer_,
      grip_num_nbrs = if (identical(method.spec$id, "grip")) grip.num.nbrs else NA_integer_,
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
      ),
      bslib::accordion_panel(
        "Stage 4: Exports and Paper Sync",
        bslib::layout_columns(
          col_widths = c(6, 6),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header("Paper synchronization map"),
            shiny::tags$p(
              style = "padding:0 1rem;color:#5f5445;line-height:1.45;margin-bottom:0.75rem;",
              shiny::tags$strong("Active manuscript source: "),
              shiny::tags$code(gripui.gmds.manuscript.source.path())
            ),
            shiny::tableOutput("gmds_paper_sync_table")
          ),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header("Exports"),
            shiny::uiOutput("gmds_export_stage_control"),
            shiny::uiOutput("gmds_export_preset_control"),
            shiny::uiOutput("gmds_figure_preset_control"),
            shiny::tableOutput("gmds_export_context"),
            shiny::tableOutput("gmds_export_summary"),
            shiny::verbatimTextOutput("gmds_export_note"),
            shiny::tags$div(
              style = "padding:0 1rem 1rem;",
              shiny::downloadButton("gmds_download_stage_trace_csv", "Download stage trace CSV"),
              shiny::downloadButton("gmds_download_stage_vertices_csv", "Download selected stage vertices CSV"),
              shiny::downloadButton("gmds_download_stage_snapshot_html", "Download selected stage snapshot HTML"),
              shiny::downloadButton("gmds_download_stage_figure_png", "Download selected stage PNG"),
              shiny::downloadButton("gmds_download_stage_figure_pdf", "Download selected stage PDF"),
              shiny::downloadButton("gmds_download_stage_bundle_zip", "Download paper-ready bundle ZIP")
            )
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
      gripui.gmds.normalize.level(
        bundle$prepared,
        level = input$gmds_focus_level,
        default = bundle$prepared$top_level_level
      )
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

    current_export_stage <- shiny::reactive({
      bundle <- bundle_state()
      if (is.null(bundle)) {
        return("top_level")
      }
      choices <- gripui.gmds.export.stage.choices(bundle)
      selected <- gripui.family.value_or_default(input$gmds_export_stage, "")
      if (!nzchar(selected) || !selected %in% unname(choices)) {
        return(unname(choices[[1L]]))
      }
      selected
    })

    current_export_payload <- shiny::reactive({
      bundle <- bundle_state()
      shiny::req(bundle)
      gripui.gmds.export.stage.payload(
        bundle = bundle,
        stage_id = current_export_stage(),
        focus_level = current_level(),
        expansion_level = current_expansion_level()
      )
    })

    current_export_preset <- shiny::reactive({
      preset <- gripui.family.value_or_default(input$gmds_export_preset, "paper_figure_bundle")
      choices <- unname(gripui.gmds.export.preset.choices())
      if (!preset %in% choices) {
        preset <- choices[[1L]]
      }
      preset
    })

    current_figure_preset <- shiny::reactive({
      preset <- gripui.family.value_or_default(input$gmds_figure_preset, "paper_panel")
      choices <- unname(gripui.gmds.figure.preset.choices())
      if (!preset %in% choices) {
        preset <- choices[[1L]]
      }
      preset
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

    output$gmds_export_stage_control <- shiny::renderUI({
      bundle <- bundle_state()
      if (is.null(bundle)) {
        return(NULL)
      }
      choices <- gripui.gmds.export.stage.choices(bundle)
      shiny::selectInput(
        "gmds_export_stage",
        "Export target",
        choices = choices,
        selected = unname(choices[[1L]])
      )
    })

    output$gmds_export_preset_control <- shiny::renderUI({
      bundle <- bundle_state()
      if (is.null(bundle)) {
        return(NULL)
      }
      choices <- gripui.gmds.export.preset.choices()
      shiny::selectInput(
        "gmds_export_preset",
        "Export preset",
        choices = choices,
        selected = unname(choices[[1L]])
      )
    })

    output$gmds_figure_preset_control <- shiny::renderUI({
      bundle <- bundle_state()
      if (is.null(bundle)) {
        return(NULL)
      }
      choices <- gripui.gmds.figure.preset.choices()
      shiny::selectInput(
        "gmds_figure_preset",
        "Figure preset",
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
          paste(
            if (is.null(bundle$payload$note) || !nzchar(bundle$payload$note)) {
              "Reference graph geometry supplied by the selected synthetic family."
            } else {
              bundle$payload$note
            },
            gripui.gmds.paper.inline.note(bundle, "reference")
          )
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
          paste(
            note,
            gripui.gmds.paper.inline.note(bundle, "misf", focus_level = level)
          )
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
            "Inactive vertices stay at their aligned reference positions for context.",
            gripui.gmds.paper.inline.note(bundle, "seed", focus_level = payload$level)
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
            ),
            gripui.gmds.paper.inline.note(bundle, "initial_placement", focus_level = payload$level)
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
            "The same active set is shown after objective-specific top-level refinement.",
            gripui.gmds.paper.inline.note(bundle, "top_level", focus_level = payload$level)
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
            "The highlighted active set is the whole current level after insertion.",
            gripui.gmds.paper.inline.note(bundle, "insertion", expansion_level = level)
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
            "The same active set is shown after objective-specific refinement of that level.",
            gripui.gmds.paper.inline.note(bundle, "refinement", expansion_level = level)
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

    output$gmds_paper_sync_table <- shiny::renderTable({
      gripui.gmds.paper.sync.table()
    }, striped = TRUE, bordered = FALSE, spacing = "s")

    output$gmds_export_context <- shiny::renderTable({
      bundle <- bundle_state()
      shiny::req(bundle)
      gripui.gmds.paper.context.table(
        bundle = bundle,
        stage_id = current_export_stage(),
        focus_level = current_level(),
        expansion_level = current_expansion_level()
      )
    }, striped = TRUE, bordered = FALSE, spacing = "s")

    output$gmds_export_summary <- shiny::renderTable({
      payload <- current_export_payload()
      shiny::req(payload)
      gripui.gmds.stage.summary.table(payload)
    }, striped = TRUE, bordered = FALSE, spacing = "s")

    output$gmds_export_note <- shiny::renderText({
      bundle <- bundle_state()
      payload <- current_export_payload()
      shiny::req(bundle, payload)
      paste(
        gripui.gmds.paper.note.lines(
          bundle = bundle,
          export_payload = payload,
          stage_id = current_export_stage(),
          focus_level = current_level(),
          expansion_level = current_expansion_level(),
          preset = current_export_preset(),
          figure_preset = current_figure_preset()
        ),
        collapse = "\n"
      )
    })

    output$gmds_download_stage_trace_csv <- shiny::downloadHandler(
      filename = function() {
        sprintf(
          "gmds_stage_trace_%s_%s.csv",
          gripui.family.value_or_default(input$gmds_family_id, "family"),
          gripui.family.value_or_default(input$gmds_method, "method")
        )
      },
      content = function(file) {
        bundle <- bundle_state()
        shiny::req(bundle)
        utils::write.csv(bundle$stage_trace, file, row.names = FALSE)
      }
    )

    output$gmds_download_stage_vertices_csv <- shiny::downloadHandler(
      filename = function() {
        payload <- current_export_payload()
        shiny::req(payload)
        label <- gsub("[^A-Za-z0-9]+", "_", tolower(payload$label))
        sprintf(
          "gmds_stage_vertices_%s_%s.csv",
          gripui.family.value_or_default(label, "stage"),
          gripui.family.value_or_default(input$gmds_method, "method")
        )
      },
      content = function(file) {
        bundle <- bundle_state()
        payload <- current_export_payload()
        shiny::req(bundle, payload)
        utils::write.csv(
          gripui.gmds.export.vertex.table(bundle, payload),
          file,
          row.names = FALSE
        )
      }
    )

    output$gmds_download_stage_snapshot_html <- shiny::downloadHandler(
      filename = function() {
        payload <- current_export_payload()
        shiny::req(payload)
        label <- gsub("[^A-Za-z0-9]+", "_", tolower(payload$label))
        sprintf(
          "gmds_stage_snapshot_%s_%s.html",
          gripui.family.value_or_default(label, "stage"),
          gripui.family.value_or_default(input$gmds_method, "method")
        )
      },
      content = function(file) {
        bundle <- bundle_state()
        payload <- current_export_payload()
        shiny::req(bundle, payload)
        if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
          stop("Package 'htmlwidgets' is required for HTML snapshot export.", call. = FALSE)
        }
        widget <- gripui.gmds.export.widget(
          bundle = bundle,
          export_payload = payload,
          vertex_alpha = as.numeric(input$gmds_vertex_alpha),
          edge_alpha = as.numeric(input$gmds_edge_alpha),
          background_alpha = as.numeric(input$gmds_background_alpha)
        )
        htmlwidgets::saveWidget(widget, file = file, selfcontained = TRUE)
      }
    )

    output$gmds_download_stage_figure_png <- shiny::downloadHandler(
      filename = function() {
        payload <- current_export_payload()
        shiny::req(payload)
        sprintf(
          "gmds_stage_figure_%s_%s_%s.png",
          gripui.gmds.export.stage.slug(payload$label),
          gripui.family.value_or_default(input$gmds_method, "method"),
          current_figure_preset()
        )
      },
      content = function(file) {
        payload <- current_export_payload()
        shiny::req(payload)
        gripui.gmds.write.static.figure(
          export_payload = payload,
          png_path = file,
          figure_preset = current_figure_preset()
        )
      }
    )

    output$gmds_download_stage_figure_pdf <- shiny::downloadHandler(
      filename = function() {
        payload <- current_export_payload()
        shiny::req(payload)
        sprintf(
          "gmds_stage_figure_%s_%s_%s.pdf",
          gripui.gmds.export.stage.slug(payload$label),
          gripui.family.value_or_default(input$gmds_method, "method"),
          current_figure_preset()
        )
      },
      content = function(file) {
        payload <- current_export_payload()
        shiny::req(payload)
        gripui.gmds.write.static.figure(
          export_payload = payload,
          pdf_path = file,
          figure_preset = current_figure_preset()
        )
      }
    )

    output$gmds_download_stage_bundle_zip <- shiny::downloadHandler(
      filename = function() {
        payload <- current_export_payload()
        shiny::req(payload)
        sprintf(
          "gmds_%s_%s_%s.zip",
          gripui.gmds.export.stage.slug(payload$label),
          gripui.family.value_or_default(input$gmds_method, "method"),
          current_export_preset()
        )
      },
      content = function(file) {
        bundle <- bundle_state()
        payload <- current_export_payload()
        shiny::req(bundle, payload)
        tmp.dir <- tempfile("gmds-export-")
        dir.create(tmp.dir, recursive = TRUE, showWarnings = FALSE)
        on.exit(unlink(tmp.dir, recursive = TRUE, force = TRUE), add = TRUE)
        files <- gripui.gmds.write.export.bundle(
          bundle = bundle,
          export_payload = payload,
          stage_id = current_export_stage(),
          focus_level = current_level(),
          expansion_level = current_expansion_level(),
          preset = current_export_preset(),
          figure_preset = current_figure_preset(),
          dir = tmp.dir
        )
        old.wd <- getwd()
        on.exit(setwd(old.wd), add = TRUE)
        setwd(tmp.dir)
        utils::zip(
          zipfile = file,
          files = basename(files),
          flags = "-r9Xq"
        )
      }
    )
  }
}

gripui.gmds.default.subtitle <- function() {
  paste(
    "Milestones 1-5: graph and geometry selection, MIS filtration, seed selection,",
    "insertion/refinement, method switching, and paper-ready export/paper sync",
    "with static PNG/PDF figure presets across canonical GRIP/GMDS/GKK/LGKK",
    "stage traces."
  )
}

#' Build the GMDS stage explorer Shiny application
#'
#' The GMDS stage explorer reuses the synthetic family catalog from
#' [gripui_family_app()] and computes a canonical MISF stage bundle for the
#' selected graph and method. The current implementation covers Milestones 1
#' through 5: graph and geometry inspection, visualization of the MIS
#' filtration, explicit seed and initial-placement panels, per-level insertion
#' and refinement panels, method switching across canonical
#' GRIP/GMDS/GKK/LGKK traces, a canonical stage-trace summary, and paper-ready
#' export/paper-synchronization helpers including static PNG/PDF figure presets.
#'
#' @param catalog Family catalog, usually [gripui_graph_family_catalog()].
#' @param title Application title.
#' @param subtitle Optional subtitle shown in the sidebar header.
#'
#' @return A `shiny.appobj`.
#' @noRd
#'
#' @examplesIf local({ old <- getOption("rgl.useNULL"); options(rgl.useNULL = TRUE); on.exit(options(rgl.useNULL = old), add = TRUE); requireNamespace("shiny", quietly = TRUE) && requireNamespace("bslib", quietly = TRUE) && requireNamespace("rgl", quietly = TRUE) })
#' app <- gripui_gmds_app()
#' inherits(app, "shiny.appobj")
gripui_gmds_app <- function(catalog = gripui.gmds.default.catalog(),
                            title = "GMDS Stage Explorer",
                            subtitle = gripui.gmds.default.subtitle()) {
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
#' @noRd
#'
#' @examplesIf local({ old <- getOption("rgl.useNULL"); options(rgl.useNULL = TRUE); on.exit(options(rgl.useNULL = old), add = TRUE); requireNamespace("shiny", quietly = TRUE) && requireNamespace("bslib", quietly = TRUE) && requireNamespace("rgl", quietly = TRUE) && requireNamespace("later", quietly = TRUE) })
#' run_gripui_gmds(launch.browser = FALSE, quiet = TRUE, auto.stop.after = 0.1)
run_gripui_gmds <- function(catalog = gripui.gmds.default.catalog(),
                            title = "GMDS Stage Explorer",
                            subtitle = gripui.gmds.default.subtitle(),
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
