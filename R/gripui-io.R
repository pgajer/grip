gripui.open.text.connection <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("path must be a single non-empty character value")
  }
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    return(gzfile(path, open = "rt"))
  }
  file(path, open = "rt")
}

gripui.read.tsv <- function(path) {
  con <- gripui.open.text.connection(path)
  on.exit(close(con), add = TRUE)
  utils::read.delim(
    con,
    sep = "\t",
    header = TRUE,
    quote = "\"",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

gripui_read_embedding_tsv <- function(path) {
  tbl <- gripui.read.tsv(path)
  if (ncol(tbl) < 2L) {
    stop("embedding TSV must contain an id column and at least one coordinate column")
  }

  coord.cols <- intersect(c("x", "y", "z"), names(tbl))
  if (length(coord.cols) == 0L) {
    numeric.cols <- names(tbl)[vapply(tbl, is.numeric, logical(1L))]
    coord.cols <- numeric.cols
  }
  if (length(coord.cols) < 2L) {
    stop("embedding TSV must contain at least two numeric coordinate columns")
  }

  id.cols <- setdiff(names(tbl), coord.cols)
  id.col <- if (length(id.cols) > 0L) id.cols[[1L]] else NULL
  list(
    data = tbl,
    id_col = id.col,
    coords = as.matrix(tbl[, coord.cols, drop = FALSE])
  )
}

gripui_read_graph_rds <- function(path) {
  obj <- readRDS(path)
  if (!is.list(obj)) {
    stop("graph RDS must contain a list")
  }
  obj
}

gripui_resolve_project_path <- function(root, path) {
  if (!is.character(root) || length(root) != 1L || is.na(root) || !nzchar(root)) {
    stop("root must be a single non-empty character value")
  }
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  gripui.normalize.path.value(path, root = root)
}

gripui.is.graph.object <- function(x) {
  is.list(x) && any(c("adj_list", "weight_list", "vertex_data", "graph_info") %in% names(x))
}

gripui.get.graph <- function(project, graph_id = NULL) {
  graph <- project$graph
  if (is.null(graph)) {
    return(NULL)
  }
  if (gripui.is.graph.object(graph)) {
    return(graph)
  }
  if (!is.list(graph) || is.null(names(graph))) {
    return(NULL)
  }
  if (is.null(graph_id) || is.na(graph_id) || !nzchar(graph_id)) {
    return(graph[[1L]])
  }
  graph[[graph_id]]
}

gripui.best.vertex.key <- function(vertex_data, ids) {
  if (!is.data.frame(vertex_data) || length(ids) == 0L) {
    return(NULL)
  }
  scores <- vapply(vertex_data, function(col) {
    sum(ids %in% as.character(col))
  }, integer(1L))
  if (length(scores) == 0L || max(scores) == 0L) {
    return(NULL)
  }
  names(scores)[[which.max(scores)]]
}

gripui_load_layout_coords <- function(project, layout) {
  if (!inherits(project, "gripui_project")) {
    stop("project must inherit from 'gripui_project'")
  }
  catalog <- project$layouts

  row <- if (is.numeric(layout) && length(layout) == 1L && !is.na(layout)) {
    catalog[as.integer(layout), , drop = FALSE]
  } else if (is.character(layout) && length(layout) == 1L) {
    catalog[catalog$layout_id == layout, , drop = FALSE]
  } else if (is.data.frame(layout) && nrow(layout) == 1L) {
    layout
  } else {
    stop("layout must be a row index, layout_id, or one-row data.frame")
  }

  if (nrow(row) != 1L) {
    stop("layout selection must resolve to exactly one row")
  }

  if ("coords" %in% names(row) && is.list(row$coords) &&
      length(row$coords) == 1L && is.matrix(row$coords[[1L]])) {
    return(row$coords[[1L]])
  }

  if (!("coords_path" %in% names(row)) || is.na(row$coords_path[[1L]]) || !nzchar(row$coords_path[[1L]])) {
    stop("selected layout does not have in-memory or on-disk coordinates")
  }

  emb <- gripui_read_embedding_tsv(row$coords_path[[1L]])
  graph <- gripui.get.graph(project, graph_id = if ("graph_id" %in% names(row)) row$graph_id[[1L]] else NULL)
  if (is.null(graph) || is.null(emb$id_col) || is.null(graph$vertex_data)) {
    return(emb$coords)
  }

  key <- gripui.best.vertex.key(graph$vertex_data, emb$data[[emb$id_col]])
  if (is.null(key)) {
    return(emb$coords)
  }

  ids <- as.character(emb$data[[emb$id_col]])
  graph.ids <- as.character(graph$vertex_data[[key]])
  match.idx <- match(ids, graph.ids)
  if (all(is.na(match.idx))) {
    return(emb$coords)
  }

  n <- if (!is.null(graph$adj_list)) length(graph$adj_list) else nrow(graph$vertex_data)
  out <- matrix(NA_real_, nrow = n, ncol = ncol(emb$coords))
  ok <- !is.na(match.idx)
  out[match.idx[ok], ] <- emb$coords[ok, , drop = FALSE]
  out
}
