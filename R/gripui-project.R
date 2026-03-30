gripui.normalize.graph.input <- function(graph) {
  if (is.null(graph)) {
    return(NULL)
  }
  if (is.character(graph) && length(graph) == 1L && !is.na(graph)) {
    return(gripui_read_graph_rds(graph))
  }
  if (gripui.is.graph.object(graph)) {
    return(graph)
  }
  if (is.list(graph) && length(graph) > 0L &&
      !is.null(names(graph)) &&
      all(vapply(graph, gripui.is.graph.object, logical(1L)))) {
    return(graph)
  }
  stop("graph must be NULL, a graph list, a named list of graph lists, or a path to an RDS file")
}

gripui.validate.graph.object <- function(graph, what = "graph") {
  if (!is.list(graph)) {
    stop(sprintf("%s must be a list", what))
  }
  if (!is.null(graph$adj_list) && !is.list(graph$adj_list)) {
    stop(sprintf("%s$adj_list must be a list when provided", what))
  }
  if (!is.null(graph$weight_list) && !is.list(graph$weight_list)) {
    stop(sprintf("%s$weight_list must be a list when provided", what))
  }
  if (!is.null(graph$adj_list) && !is.null(graph$weight_list) &&
      length(graph$adj_list) != length(graph$weight_list)) {
    stop(sprintf("%s$weight_list must match %s$adj_list in length", what, what))
  }
  if (!is.null(graph$vertex_data) && !is.data.frame(graph$vertex_data)) {
    stop(sprintf("%s$vertex_data must be a data.frame when provided", what))
  }
  if (!is.null(graph$vertex_data) && !is.null(graph$adj_list) &&
      nrow(graph$vertex_data) != length(graph$adj_list)) {
    stop(sprintf("%s$vertex_data must have one row per vertex", what))
  }
  invisible(TRUE)
}

gripui.normalize.layouts.df <- function(layouts) {
  if (!is.data.frame(layouts)) {
    stop("layouts must be a data.frame")
  }
  layouts <- gripui.normalize.catalog.names(layouts)
  if (!("candidate" %in% names(layouts))) {
    layouts$candidate <- sprintf("layout_%03d", seq_len(nrow(layouts)))
  }
  if (!("stage" %in% names(layouts))) {
    layouts$stage <- rep("layout", nrow(layouts))
  }
  if (!("seed" %in% names(layouts))) {
    layouts$seed <- NA_integer_
  }
  if (!("status" %in% names(layouts))) {
    layouts$status <- rep("ok", nrow(layouts))
  }
  if (!("graph_id" %in% names(layouts))) {
    layouts$graph_id <- rep(NA_character_, nrow(layouts))
  }
  if (!("layout_id" %in% names(layouts))) {
    layouts$layout_id <- mapply(
      gripui.make.layout.id,
      stage = layouts$stage,
      candidate = layouts$candidate,
      seed = layouts$seed,
      USE.NAMES = FALSE
    )
  }
  if (!("color_view_default" %in% names(layouts))) {
    layouts$color_view_default <- rep("plain", nrow(layouts))
  }
  if (!("coords_path" %in% names(layouts))) {
    layouts$coords_path <- rep(NA_character_, nrow(layouts))
  }
  if (!("html_path" %in% names(layouts))) {
    layouts$html_path <- rep(NA_character_, nrow(layouts))
  }
  if (!("gif_path" %in% names(layouts))) {
    layouts$gif_path <- rep(NA_character_, nrow(layouts))
  }
  if (!("thumbnail_path" %in% names(layouts))) {
    layouts$thumbnail_path <- rep(NA_character_, nrow(layouts))
  }
  gripui.mark.layout.availability(layouts)
}

#' Create a normalized `gripui` project object
#'
#' @param graph Optional graph object used by layout viewers. This may be a
#'   single list containing fields such as `adj_list`, `weight_list`,
#'   `vertex_data`, and `graph_info`, a named list of such graph objects for
#'   multi-stage projects, or `NULL` when only catalog exploration is needed.
#' @param layouts Data frame with one row per realized layout.
#' @param title Optional project title.
#' @param subtitle Optional subtitle.
#' @param notes Optional notes string.
#'
#' @return An object of class `gripui_project`.
#' @export
#'
#' @examples
#' graph <- list(adj_list = list(2L, c(1L, 3L), c(2L, 4L), 3L))
#' layouts <- data.frame(
#'   candidate = "path.default",
#'   stage = "layout",
#'   seed = 1L,
#'   status = "ok",
#'   stringsAsFactors = FALSE
#' )
#' project <- gripui_project(graph = graph, layouts = layouts, title = "Path graph")
#' project$meta$title
gripui_project <- function(graph = NULL,
                           layouts,
                           title = NULL,
                           subtitle = NULL,
                           notes = NULL) {
  graph <- gripui.normalize.graph.input(graph)
  if (!is.null(graph)) {
    if (gripui.is.graph.object(graph)) {
      gripui.validate.graph.object(graph, what = "graph")
    } else {
      for (name in names(graph)) {
        gripui.validate.graph.object(graph[[name]], what = sprintf("graph[['%s']]", name))
      }
    }
  }

  layouts <- gripui.normalize.layouts.df(layouts)
  out <- list(
    graph = graph,
    layouts = layouts,
    meta = list(
      title = title,
      subtitle = subtitle,
      notes = notes
    )
  )
  class(out) <- "gripui_project"
  gripui_validate_project(out)
  out
}

#' Validate a `gripui_project`
#'
#' @param project Object to validate.
#'
#' @return Invisibly returns `TRUE` when validation succeeds.
#' @export
#'
#' @examples
#' layouts <- data.frame(
#'   candidate = "test.layout",
#'   stage = "layout",
#'   seed = 1L,
#'   status = "ok",
#'   stringsAsFactors = FALSE
#' )
#' project <- gripui_project(graph = NULL, layouts = layouts, title = "Validation example")
#' gripui_validate_project(project)
gripui_validate_project <- function(project) {
  if (!is.list(project)) {
    stop("project must be a list")
  }
  if (!("layouts" %in% names(project))) {
    stop("project must contain a 'layouts' field")
  }
  if (!is.data.frame(project$layouts)) {
    stop("project$layouts must be a data.frame")
  }

  req <- c("layout_id", "candidate", "stage", "seed", "status", "viewable")
  miss <- setdiff(req, names(project$layouts))
  if (length(miss) > 0L) {
    stop("project$layouts is missing required columns: ", paste(miss, collapse = ", "))
  }
  if (anyDuplicated(project$layouts$layout_id)) {
    stop("project$layouts$layout_id values must be unique")
  }
  if ("coords" %in% names(project$layouts) && !is.list(project$layouts$coords)) {
    stop("project$layouts$coords must be a list-column when provided")
  }

  graph <- gripui.normalize.graph.input(project$graph)
  if (!is.null(graph)) {
    if (gripui.is.graph.object(graph)) {
      gripui.validate.graph.object(graph, what = "project$graph")
    } else {
      for (name in names(graph)) {
        gripui.validate.graph.object(graph[[name]], what = sprintf("project$graph[['%s']]", name))
      }
      if ("graph_id" %in% names(project$layouts)) {
        used <- unique(as.character(project$layouts$graph_id))
        used <- used[!is.na(used) & nzchar(used)]
        bad <- setdiff(used, names(graph))
        if (length(bad) > 0L) {
          stop("project$layouts$graph_id references unknown graph ids: ", paste(bad, collapse = ", "))
        }
      }
    }
  }

  invisible(TRUE)
}

#' Convert `grip.compare.layouts()` output into a `gripui_project`
#'
#' @param compare_obj Result of `grip.compare.layouts()`.
#' @param graph Optional graph object, named graph list, or path to a graph RDS.
#'   Defaults to `NULL` when only catalog exploration is needed.
#' @param vertex_data Optional vertex metadata added to `graph`.
#' @param graph_info Optional graph-level metadata added to `graph`.
#' @param title Optional project title.
#'
#' @return A `gripui_project`.
#' @export
#'
#' @examples
#' edges <- edges.path(5)
#' cmp <- grip.compare.layouts(
#'   edges = edges,
#'   n = 5,
#'   candidates = "default",
#'   seeds = 1L,
#'   return.layouts = TRUE
#' )
#' graph <- list(
#'   adj_list = list(2L, c(1L, 3L), c(2L, 4L), c(3L, 5L), 4L)
#' )
#' project <- gripui_project_from_compare(cmp, graph = graph, title = "Path compare")
#' nrow(project$layouts)
gripui_project_from_compare <- function(compare_obj,
                                        graph = NULL,
                                        vertex_data = NULL,
                                        graph_info = NULL,
                                        title = NULL) {
  if (!is.list(compare_obj) || !is.data.frame(compare_obj$runs)) {
    stop("compare_obj must be a result from grip.compare.layouts()")
  }
  graph <- gripui.normalize.graph.input(graph)
  if (is.null(graph)) {
    graph <- list()
  }
  if (gripui.is.graph.object(graph)) {
    if (!is.null(vertex_data)) {
      graph$vertex_data <- vertex_data
    }
    if (!is.null(graph_info)) {
      graph$graph_info <- graph_info
    }
  }

  layouts <- gripui.build.catalog.from.runs(
    runs = compare_obj$runs,
    stage = "compare",
    graph_id = "graph"
  )
  layouts <- gripui.attach.compare.layouts(layouts, compare_obj$layouts)
  gripui_project(
    graph = graph,
    layouts = layouts,
    title = if (is.null(title)) "Layout comparison" else title
  )
}

gripui.read.optional.tsv <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  gripui.read.tsv(path)
}

gripui.read.optional.text <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

gripui.build.bundle.project <- function(root, graph, title, subtitle) {
  catalog.path <- file.path(root, "catalog", "layout_catalog.tsv")
  layouts <- gripui.read.tsv(catalog.path)
  layouts <- gripui.normalize.catalog.names(layouts)
  layouts <- gripui.resolve.catalog.paths(layouts, root = root)
  layouts <- gripui.normalize.layouts.df(layouts)

  if (is.null(graph)) {
    graph.path <- file.path(root, "graph", "graph.rds")
    if (file.exists(graph.path)) {
      graph <- gripui_read_graph_rds(graph.path)
    }
  }

  gripui_project(
    graph = graph,
    layouts = layouts,
    title = if (is.null(title)) basename(root) else title,
    subtitle = subtitle,
    notes = gripui.read.optional.text(file.path(root, "metadata", "notes.md"))
  )
}

gripui.find.hmp.repulsion.root <- function(root) {
  if (file.exists(file.path(root, "repulsion_sweep_runs.tsv"))) {
    return(root)
  }
  parent <- dirname(root)
  kids <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
  kids <- kids[grepl("^hmp_u01_gc_repulsion_sweep_from_088_", basename(kids))]
  kids <- kids[file.exists(file.path(kids, "repulsion_sweep_runs.tsv"))]
  if (length(kids) == 0L) {
    return(NULL)
  }
  sort(kids, decreasing = TRUE)[[1L]]
}

gripui.is.hmp.root <- function(root) {
  any(file.exists(c(
    file.path(root, "coarse_stage", "coarse_stage_runs.tsv"),
    file.path(root, "full_stage1", "full_stage1_runs.tsv"),
    file.path(root, "full_stage2", "full_stage2_runs.tsv"),
    file.path(root, "repulsion_sweep_runs.tsv")
  )))
}

gripui.build.hmp.project <- function(root, graph, title, subtitle) {
  pieces <- list()

  coarse.path <- file.path(root, "coarse_stage", "coarse_stage_runs.tsv")
  if (file.exists(coarse.path)) {
    pieces[[length(pieces) + 1L]] <- gripui.build.catalog.from.runs(
      runs = gripui.read.tsv(coarse.path),
      stage = "coarse_stage",
      root = root,
      graph_id = "coarse"
    )
  }

  full1.path <- file.path(root, "full_stage1", "full_stage1_runs.tsv")
  if (file.exists(full1.path)) {
    pieces[[length(pieces) + 1L]] <- gripui.build.catalog.from.runs(
      runs = gripui.read.tsv(full1.path),
      stage = "full_stage1",
      root = root,
      graph_id = "full"
    )
  }

  full2.path <- file.path(root, "full_stage2", "full_stage2_runs.tsv")
  if (file.exists(full2.path)) {
    pieces[[length(pieces) + 1L]] <- gripui.build.catalog.from.runs(
      runs = gripui.read.tsv(full2.path),
      stage = "full_stage2",
      root = root,
      graph_id = "full",
      manifest = gripui.read.optional.tsv(file.path(root, "top_layout_visual_manifest.tsv")),
      manifest_by = c("candidate", "seed"),
      html.preference = c("plain_html", "cst_html")
    )
  }

  repulsion.root <- gripui.find.hmp.repulsion.root(root)
  if (!is.null(repulsion.root)) {
    pieces[[length(pieces) + 1L]] <- gripui.build.catalog.from.runs(
      runs = gripui.read.tsv(file.path(repulsion.root, "repulsion_sweep_runs.tsv")),
      stage = "repulsion_sweep",
      root = repulsion.root,
      graph_id = "full",
      manifest = gripui.read.optional.tsv(file.path(repulsion.root, "repulsion_sweep_visual_manifest.tsv")),
      manifest_by = c("repulsion_factor", "seed"),
      html.preference = c("cst_html"),
      gif.preference = character(0L),
      default_color_view = "cst"
    )
  }

  layouts <- gripui.row.bind.fill(pieces)
  if (nrow(layouts) == 0L) {
    stop("no HMP/U01 run tables were found under root")
  }
  layouts <- gripui.normalize.layouts.df(layouts)

  notes <- gripui.read.optional.text(
    file.path(root, "hmp_u01_gc_coarsened_layout_selection_summary_2026-03-23.md")
  )
  if (is.null(notes) && !is.null(repulsion.root)) {
    notes <- gripui.read.optional.text(file.path(repulsion.root, "repulsion_sweep_summary.md"))
  }

  gripui_project(
    graph = graph,
    layouts = layouts,
    title = if (is.null(title)) basename(root) else title,
    subtitle = subtitle,
    notes = notes
  )
}

#' Create a `gripui_project` from a directory of saved artifacts
#'
#' @param root Directory containing either a normalized `gripui` bundle or a
#'   saved search output tree such as the current HMP/U01 layout-selection
#'   outputs. The normalized bundle format is the preferred long-term producer
#'   contract; support for HMP-style run tables and manifests is a compatibility
#'   adapter for older outputs.
#' @param graph Optional graph object, named graph list, or path to a graph RDS.
#' @param title Optional project title.
#' @param subtitle Optional subtitle.
#'
#' @return A `gripui_project`.
#' @export
#'
#' @examples
#' root <- tempfile("gripui-bundle-")
#' dir.create(root)
#' on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
#' dir.create(file.path(root, "catalog"), recursive = TRUE, showWarnings = FALSE)
#' dir.create(file.path(root, "graph"), recursive = TRUE, showWarnings = FALSE)
#' dir.create(file.path(root, "artifacts", "layouts"), recursive = TRUE, showWarnings = FALSE)
#'
#' graph <- list(adj_list = list(2L, c(1L, 3L), 2L))
#' saveRDS(graph, file.path(root, "graph", "graph.rds"))
#'
#' coords.path <- file.path(root, "artifacts", "layouts", "toy_embedding.tsv")
#' utils::write.table(
#'   data.frame(
#'     vertex_id = c("v1", "v2", "v3"),
#'     x = c(1, 2, 3),
#'     y = c(0, 0, 0),
#'     stringsAsFactors = FALSE
#'   ),
#'   file = coords.path,
#'   sep = "\t",
#'   quote = TRUE,
#'   row.names = FALSE,
#'   col.names = TRUE
#' )
#' utils::write.table(
#'   data.frame(
#'     candidate = "toy.layout",
#'     stage = "bundle",
#'     seed = 1L,
#'     status = "ok",
#'     coords_path = file.path("artifacts", "layouts", "toy_embedding.tsv"),
#'     stringsAsFactors = FALSE
#'   ),
#'   file = file.path(root, "catalog", "layout_catalog.tsv"),
#'   sep = "\t",
#'   quote = TRUE,
#'   row.names = FALSE,
#'   col.names = TRUE
#' )
#'
#' project <- gripui_project_from_dir(root, title = "My project")
#' project$layouts$availability
gripui_project_from_dir <- function(root,
                                    graph = NULL,
                                    title = NULL,
                                    subtitle = NULL) {
  if (!is.character(root) || length(root) != 1L || is.na(root) || !nzchar(root)) {
    stop("root must be a single non-empty character value")
  }
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  graph <- gripui.normalize.graph.input(graph)

  if (file.exists(file.path(root, "catalog", "layout_catalog.tsv"))) {
    return(gripui.build.bundle.project(
      root = root,
      graph = graph,
      title = title,
      subtitle = subtitle
    ))
  }

  if (gripui.is.hmp.root(root)) {
    return(gripui.build.hmp.project(
      root = root,
      graph = graph,
      title = title,
      subtitle = subtitle
    ))
  }

  stop("root does not look like a supported gripui project directory")
}
