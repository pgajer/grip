gripui.is.empty.string <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && !nzchar(trimws(x))
}

gripui.normalize.catalog.names <- function(x) {
  if (!is.data.frame(x)) {
    stop("x must be a data.frame")
  }
  nm <- names(x)
  nm <- tolower(nm)
  nm <- gsub("[^[:alnum:]]+", "_", nm)
  nm <- gsub("^_+|_+$", "", nm)
  nm <- gsub("_+", "_", nm)
  names(x) <- nm
  x
}

gripui.make.layout.id <- function(stage, candidate, seed = NA_integer_) {
  parts <- c(stage, candidate)
  if (length(seed) == 1L && is.finite(seed) && !is.na(seed)) {
    parts <- c(parts, sprintf("seed_%03d", as.integer(seed)))
  }
  out <- paste(parts, collapse = "__")
  out <- tolower(out)
  out <- gsub("[^[:alnum:]]+", "_", out)
  out <- gsub("^_+|_+$", "", out)
  gsub("_+", "_", out)
}

gripui.normalize.path.value <- function(path, root) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    return(NA_character_)
  }
  path <- as.character(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    path <- file.path(root, path)
  }
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (!file.exists(path)) {
    return(NA_character_)
  }
  path
}

gripui.resolve.catalog.paths <- function(catalog, root) {
  if (!is.data.frame(catalog) || nrow(catalog) == 0L) {
    return(catalog)
  }
  path.cols <- grep("(_path|_html|_gif|_png|_tsv)$", names(catalog), value = TRUE)
  if (length(path.cols) == 0L) {
    return(catalog)
  }
  for (col in path.cols) {
    catalog[[col]] <- vapply(
      catalog[[col]],
      gripui.normalize.path.value,
      character(1L),
      root = root
    )
  }
  catalog
}

gripui.row.bind.fill <- function(items) {
  items <- items[!vapply(items, is.null, logical(1L))]
  if (length(items) == 0L) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  cols <- unique(unlist(lapply(items, names), use.names = FALSE))
  out <- lapply(items, function(df) {
    miss <- setdiff(cols, names(df))
    for (col in miss) {
      df[[col]] <- NA
    }
    df[, cols, drop = FALSE]
  })
  rownames(out[[1L]]) <- NULL
  do.call(rbind, out)
}

gripui.first.nonempty.column <- function(df, choices) {
  for (name in choices) {
    if (!(name %in% names(df))) {
      next
    }
    vals <- as.character(df[[name]])
    vals[!nzchar(vals)] <- NA_character_
    if (any(!is.na(vals))) {
      return(vals)
    }
  }
  rep(NA_character_, nrow(df))
}

gripui.mark.layout.availability <- function(catalog) {
  coords.in.memory <- "coords" %in% names(catalog)
  coords.available <- if (coords.in.memory) {
    vapply(
      catalog$coords,
      function(x) is.matrix(x) && nrow(x) > 0L,
      logical(1L)
    )
  } else {
    rep(FALSE, nrow(catalog))
  }

  coords.path.available <- if ("coords_path" %in% names(catalog)) {
    is.character(catalog$coords_path) & !is.na(catalog$coords_path) & nzchar(catalog$coords_path)
  } else {
    rep(FALSE, nrow(catalog))
  }

  html.available <- if ("html_path" %in% names(catalog)) {
    is.character(catalog$html_path) & !is.na(catalog$html_path) & nzchar(catalog$html_path)
  } else {
    rep(FALSE, nrow(catalog))
  }

  gif.available <- if ("gif_path" %in% names(catalog)) {
    is.character(catalog$gif_path) & !is.na(catalog$gif_path) & nzchar(catalog$gif_path)
  } else {
    rep(FALSE, nrow(catalog))
  }

  thumb.available <- if ("thumbnail_path" %in% names(catalog)) {
    is.character(catalog$thumbnail_path) & !is.na(catalog$thumbnail_path) & nzchar(catalog$thumbnail_path)
  } else {
    rep(FALSE, nrow(catalog))
  }

  interactive <- coords.available | coords.path.available
  artifact.only <- !interactive & (html.available | gif.available | thumb.available)
  availability <- rep("summary-only", nrow(catalog))
  availability[artifact.only] <- "artifact-only"
  availability[interactive] <- "interactive"
  catalog$availability <- availability
  catalog$viewable <- interactive | artifact.only
  catalog
}

gripui.attach.artifact.manifest <- function(catalog,
                                            manifest,
                                            by.cols,
                                            root,
                                            html.preference = c("plain_html", "cst_html"),
                                            gif.preference = c("plain_gif", "cst_gif"),
                                            coords.preference = c("coords_path", "embedding_tsv")) {
  if (is.null(manifest) || !is.data.frame(manifest) || nrow(manifest) == 0L) {
    return(catalog)
  }
  manifest <- gripui.normalize.catalog.names(manifest)
  manifest <- gripui.resolve.catalog.paths(manifest, root = root)

  missing.by <- setdiff(by.cols, intersect(names(catalog), names(manifest)))
  if (length(missing.by) > 0L) {
    return(catalog)
  }

  keep.cols <- unique(c(
    by.cols,
    intersect(c(html.preference, gif.preference, coords.preference,
                "plain_html", "cst_html", "subcst_html",
                "dcst_depth1_absorb_html", "dcst_depth2_absorb_html",
                "continuous_html", "plain_gif", "cst_gif",
                "embedding_tsv", "thumbnail_path"),
              names(manifest))
  ))
  manifest <- manifest[, keep.cols, drop = FALSE]
  merged <- merge(catalog, manifest, by = by.cols, all.x = TRUE, sort = FALSE)
  ord <- match(catalog$layout_id, merged$layout_id)
  if (all(!is.na(ord))) {
    merged <- merged[ord, , drop = FALSE]
  }

  if (!("coords_path" %in% names(merged))) {
    merged$coords_path <- gripui.first.nonempty.column(merged, coords.preference)
  } else {
    fill <- gripui.first.nonempty.column(merged, coords.preference)
    use.fill <- is.na(merged$coords_path) | !nzchar(merged$coords_path)
    merged$coords_path[use.fill] <- fill[use.fill]
  }
  if (!("html_path" %in% names(merged))) {
    merged$html_path <- gripui.first.nonempty.column(merged, html.preference)
  } else {
    fill <- gripui.first.nonempty.column(merged, html.preference)
    use.fill <- is.na(merged$html_path) | !nzchar(merged$html_path)
    merged$html_path[use.fill] <- fill[use.fill]
  }
  if (!("gif_path" %in% names(merged))) {
    merged$gif_path <- gripui.first.nonempty.column(merged, gif.preference)
  } else {
    fill <- gripui.first.nonempty.column(merged, gif.preference)
    use.fill <- is.na(merged$gif_path) | !nzchar(merged$gif_path)
    merged$gif_path[use.fill] <- fill[use.fill]
  }

  merged
}

gripui.default.graph.id <- function(stage) {
  if (!is.character(stage) || length(stage) != 1L || is.na(stage)) {
    return("graph")
  }
  if (grepl("^coarse", stage)) {
    return("coarse")
  }
  if (grepl("^full|^repulsion", stage)) {
    return("full")
  }
  "graph"
}

gripui.build.catalog.from.runs <- function(runs,
                                           stage,
                                           root = NULL,
                                           graph_id = NULL,
                                           manifest = NULL,
                                           manifest_by = c("candidate", "seed"),
                                           html.preference = c("plain_html", "cst_html"),
                                           gif.preference = c("plain_gif", "cst_gif"),
                                           default_color_view = "plain") {
  if (!is.data.frame(runs)) {
    stop("runs must be a data.frame")
  }
  runs <- gripui.normalize.catalog.names(runs)
  if (!all(c("candidate", "seed", "status") %in% names(runs))) {
    stop("runs must contain candidate, seed, and status columns")
  }

  out <- runs
  out$stage <- rep(as.character(stage), nrow(out))
  if (!("graph_id" %in% names(out))) {
    gid <- if (is.null(graph_id)) gripui.default.graph.id(stage) else as.character(graph_id)
    out$graph_id <- rep(gid, nrow(out))
  }
  out$layout_id <- mapply(
    gripui.make.layout.id,
    stage = out$stage,
    candidate = out$candidate,
    seed = out$seed,
    USE.NAMES = FALSE
  )
  out$color_view_default <- rep(default_color_view, nrow(out))
  out <- gripui.resolve.catalog.paths(out, root = if (is.null(root)) getwd() else root)

  if (!is.null(manifest)) {
    out <- gripui.attach.artifact.manifest(
      catalog = out,
      manifest = manifest,
      by.cols = manifest_by,
      root = if (is.null(root)) getwd() else root,
      html.preference = html.preference,
      gif.preference = gif.preference
    )
  }

  out <- gripui.mark.layout.availability(out)
  out
}

gripui.attach.compare.layouts <- function(catalog, layouts) {
  if (is.null(layouts) || !is.list(layouts) || nrow(catalog) == 0L) {
    return(gripui.mark.layout.availability(catalog))
  }
  coords <- vector("list", nrow(catalog))
  for (i in seq_len(nrow(catalog))) {
    cand <- as.character(catalog$candidate[[i]])
    seed <- as.character(catalog$seed[[i]])
    if (!(cand %in% names(layouts))) {
      coords[[i]] <- NULL
      next
    }
    cand.layouts <- layouts[[cand]]
    if (!is.list(cand.layouts) || !(seed %in% names(cand.layouts))) {
      coords[[i]] <- NULL
      next
    }
    coords[[i]] <- as.matrix(cand.layouts[[seed]])
  }
  catalog$coords <- coords
  gripui.mark.layout.availability(catalog)
}
