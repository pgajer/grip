gripui.enable.rgl.null.device <- function() {
  old <- getOption("rgl.useNULL")
  if (!isTRUE(old)) {
    options(rgl.useNULL = TRUE)
  }
  old
}

gripui.require.app.packages <- function() {
  old <- gripui.enable.rgl.null.device()
  pkgs <- c("shiny", "bslib", "DT", "htmltools", "rgl")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing) > 0L) {
    options(rgl.useNULL = old)
    stop(
      "gripui requires optional packages that are not installed: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  old
}

gripui.catalog.numeric.columns <- function(catalog) {
  if (!is.data.frame(catalog) || nrow(catalog) == 0L) {
    return(character(0))
  }
  keep <- vapply(catalog, function(x) is.numeric(x) || is.integer(x), logical(1L))
  setdiff(names(catalog)[keep], c("viewable"))
}

gripui.default.axis <- function(choices, preferred, fallback = NULL) {
  hit <- preferred[preferred %in% choices]
  if (length(hit) > 0L) {
    return(hit[[1L]])
  }
  if (!is.null(fallback) && fallback %in% choices) {
    return(fallback)
  }
  if (length(choices) > 0L) {
    return(choices[[1L]])
  }
  ""
}

gripui.display.columns <- function(catalog) {
  preferred <- c(
    "candidate",
    "stage",
    "seed",
    "availability",
    "status",
    "score_composite_extended",
    "score_composite",
    "sampled_stress",
    "cst_cluster_separation",
    "subcst_cluster_separation",
    "sampled_nonedge_sep_ratio",
    "repulsion_factor",
    "rounds",
    "final_rounds",
    "num_nbrs",
    "r",
    "s"
  )
  cols <- preferred[preferred %in% names(catalog)]
  if (length(cols) == 0L) {
    cols <- names(catalog)
  }
  cols
}

gripui.selected.info <- function(row) {
  keep <- c(
    "layout_id",
    "candidate",
    "stage",
    "graph_id",
    "seed",
    "availability",
    "status",
    "score_composite_extended",
    "score_composite",
    "sampled_stress",
    "cst_cluster_separation",
    "subcst_cluster_separation",
    "sampled_nonedge_sep_ratio",
    "repulsion_factor",
    "placement",
    "rounds",
    "final_rounds",
    "num_init",
    "num_nbrs",
    "r",
    "s",
    "elapsed_sec"
  )
  keep <- keep[keep %in% names(row)]
  if (length(keep) == 0L) {
    return(data.frame(field = character(0), value = character(0), stringsAsFactors = FALSE))
  }
  vals <- lapply(row[1, keep, drop = FALSE], function(x) {
    if (length(x) == 0L || is.null(x) || is.na(x)) {
      return(NA_character_)
    }
    if (is.numeric(x)) {
      return(format(x, digits = 6))
    }
    as.character(x)
  })
  data.frame(
    field = keep,
    value = unlist(vals, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

gripui.artifact.links <- function(row) {
  make.link <- function(path, label) {
    if (length(path) != 1L || is.na(path) || !nzchar(path)) {
      return(NULL)
    }
    shiny::tags$li(
      shiny::tags$span(
        class = "gripui-artifact-entry",
        shiny::tags$span(class = "gripui-artifact-label", label),
        shiny::tags$code(class = "gripui-artifact-path", path)
      )
    )
  }

  links <- list(
    make.link(row$coords_path[[1L]], "Open embedding TSV"),
    make.link(row$html_path[[1L]], "Open saved HTML"),
    make.link(row$gif_path[[1L]], "Open saved GIF"),
    make.link(row$thumbnail_path[[1L]], "Open thumbnail")
  )
  links <- Filter(Negate(is.null), links)
  if (length(links) == 0L) {
    return(shiny::tags$p(class = "gripui-muted", "No external artifact files are attached to this layout."))
  }
  shiny::tags$ul(class = "gripui-links", links)
}

gripui.preview.image.path <- function(row) {
  candidates <- c("gif_path", "thumbnail_path")
  for (name in candidates) {
    if (!(name %in% names(row))) {
      next
    }
    path <- row[[name]][[1L]]
    if (is.character(path) && length(path) == 1L && !is.na(path) && nzchar(path)) {
      return(path)
    }
  }
  NULL
}

gripui.image.content.type <- function(path) {
  ext <- tolower(tools::file_ext(path))
  switch(
    ext,
    gif = "image/gif",
    png = "image/png",
    jpg = "image/jpeg",
    jpeg = "image/jpeg",
    svg = "image/svg+xml",
    NULL
  )
}

gripui.graph.edges <- function(graph) {
  if (is.null(graph) || is.null(graph$adj_list)) {
    return(NULL)
  }
  grip.edges.from.adj.list(graph$adj_list)
}

gripui.vertex.color.values <- function(graph, n, color_by) {
  if (is.null(graph) || is.null(graph$vertex_data) ||
      !is.character(color_by) || length(color_by) != 1L ||
      !nzchar(color_by) || identical(color_by, "plain") ||
      !(color_by %in% names(graph$vertex_data))) {
    return(rep("#1f3b73", n))
  }

  vals <- graph$vertex_data[[color_by]]
  if (length(vals) != n) {
    return(rep("#1f3b73", n))
  }
  if (is.numeric(vals)) {
    if (!any(is.finite(vals))) {
      return(rep("#bdbdbd", length(vals)))
    }
    rng <- range(vals, finite = TRUE, na.rm = TRUE)
    pal <- grDevices::colorRampPalette(c("#d8e2dc", "#1f3b73"))(100)
    out <- rep("#bdbdbd", length(vals))
    ok <- is.finite(vals)
    if (all(ok) && diff(rng) == 0) {
      out[ok] <- pal[[75L]]
      return(out)
    }
    if (any(ok)) {
      scaled <- (vals[ok] - rng[[1L]]) / max(diff(rng), .Machine$double.eps)
      idx <- pmax(1L, pmin(100L, floor(scaled * 99L) + 1L))
      out[ok] <- pal[idx]
    }
    return(out)
  }

  fac <- as.factor(vals)
  lev <- levels(fac)
  pal <- stats::setNames(grDevices::hcl.colors(max(length(lev), 1L), "Dark 3"), lev)
  out <- pal[as.character(fac)]
  out[is.na(out)] <- "#bdbdbd"
  unname(out)
}

gripui.valid.coord.rows <- function(coords) {
  apply(coords, 1L, function(x) all(is.finite(x)))
}

gripui.render.layout.plot2d <- function(coords,
                                        graph = NULL,
                                        color_by = "plain",
                                        show_edges = TRUE) {
  keep <- gripui.valid.coord.rows(coords)
  if (!any(keep)) {
    graphics::plot.new()
    graphics::title("No finite coordinates available")
    return(invisible(NULL))
  }
  cols <- gripui.vertex.color.values(graph, nrow(coords), color_by = color_by)
  xx <- coords[, 1L]
  yy <- coords[, 2L]

  xlim <- range(xx[keep], finite = TRUE)
  ylim <- range(yy[keep], finite = TRUE)
  xpad <- max(diff(xlim) * 0.06, 1e-8)
  ypad <- max(diff(ylim) * 0.06, 1e-8)

  graphics::plot(
    xx, yy,
    type = "n",
    asp = 1,
    xlab = "x",
    ylab = "y",
    xlim = xlim + c(-xpad, xpad),
    ylim = ylim + c(-ypad, ypad)
  )

  edges <- gripui.graph.edges(graph)
  if (isTRUE(show_edges) && !is.null(edges) && nrow(edges) > 0L) {
    good.edges <- keep[edges[, 1L]] & keep[edges[, 2L]]
    edges <- edges[good.edges, , drop = FALSE]
    if (nrow(edges) > 0L) {
      apply(edges, 1L, function(e) {
        graphics::segments(
          xx[e[1L]], yy[e[1L]],
          xx[e[2L]], yy[e[2L]],
          col = grDevices::adjustcolor("#9ca3af", alpha.f = 0.45)
        )
      })
    }
  }

  graphics::points(xx[keep], yy[keep], pch = 16, cex = 0.7, col = cols[keep])
  invisible(NULL)
}

gripui.render.rglwidget <- function(coords,
                                    graph = NULL,
                                    color_by = "plain",
                                    show_edges = TRUE) {
  old <- gripui.enable.rgl.null.device()
  on.exit(options(rgl.useNULL = old), add = TRUE)

  rgl::open3d(useNULL = TRUE)
  on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)
  rgl::clear3d(type = "shapes")
  rgl::clear3d(type = "lights")
  rgl::bg3d(color = "white")
  rgl::light3d()

  keep <- gripui.valid.coord.rows(coords)
  cols <- gripui.vertex.color.values(graph, nrow(coords), color_by = color_by)
  if (any(keep)) {
    rgl::points3d(
      coords[keep, 1L],
      coords[keep, 2L],
      coords[keep, 3L],
      col = cols[keep],
      size = 5
    )
  }

  edges <- gripui.graph.edges(graph)
  if (isTRUE(show_edges) && !is.null(edges) && nrow(edges) > 0L) {
    good.edges <- keep[edges[, 1L]] & keep[edges[, 2L]]
    edges <- edges[good.edges, , drop = FALSE]
    if (nrow(edges) > 0L) {
      for (i in seq_len(nrow(edges))) {
        e <- edges[i, ]
        rgl::segments3d(
          x = coords[e, 1L],
          y = coords[e, 2L],
          z = coords[e, 3L],
          col = grDevices::adjustcolor("#9ca3af", alpha.f = 0.30),
          lwd = 1
        )
      }
    }
  }

  rgl::rglwidget(width = "100%", height = 520)
}

gripui.landscape.colors <- function(values) {
  if (is.numeric(values)) {
    if (!any(is.finite(values))) {
      return(rep("#9ca3af", length(values)))
    }
    rng <- range(values, finite = TRUE, na.rm = TRUE)
    pal <- grDevices::colorRampPalette(c("#d8e2dc", "#1f3b73"))(100)
    out <- rep("#9ca3af", length(values))
    ok <- is.finite(values)
    if (all(ok) && diff(rng) == 0) {
      out[ok] <- pal[[75L]]
      return(out)
    }
    if (any(ok)) {
      scaled <- (values[ok] - rng[[1L]]) / max(diff(rng), .Machine$double.eps)
      idx <- pmax(1L, pmin(100L, floor(scaled * 99L) + 1L))
      out[ok] <- pal[idx]
    }
    return(out)
  }
  fac <- as.factor(values)
  lev <- levels(fac)
  pal <- stats::setNames(grDevices::hcl.colors(max(length(lev), 1L), "Dark 3"), lev)
  out <- pal[as.character(fac)]
  out[is.na(out)] <- "#9ca3af"
  unname(out)
}
