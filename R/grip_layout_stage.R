grip.new.misf.grip.prepared <- function(x = list()) {
  structure(x, class = c("grip_misf_grip_prepared", "list"))
}

grip.new.misf.grip.fit <- function(x = list()) {
  structure(x, class = c("grip_misf_grip_fit", "list"))
}

grip.layout.trace.active.vertices <- function(frame) {
  frame <- as.matrix(frame)
  if (!nrow(frame) || !ncol(frame)) {
    return(integer())
  }
  keep <- rowSums(is.finite(frame)) == ncol(frame)
  as.integer(which(keep))
}

grip.layout.trace.frame.at <- function(trace, row_index) {
  row_index <- as.integer(row_index[[1L]])
  if (!is.finite(row_index) || row_index < 1L || row_index > nrow(trace$meta)) {
    stop("row_index must identify a traced frame")
  }
  frame_id <- as.integer(trace$meta$frame[[row_index]])
  as.matrix(trace$frames[[frame_id]])
}

grip.layout.trace.frames.at <- function(trace, row_indices) {
  row_indices <- unique(as.integer(row_indices))
  row_indices <- row_indices[is.finite(row_indices) & row_indices >= 1L & row_indices <= nrow(trace$meta)]
  if (!length(row_indices)) {
    return(list())
  }
  lapply(row_indices, function(row_index) {
    grip.layout.trace.frame.at(trace, row_index)
  })
}

grip.layout.trace.induced.edge.count <- function(edges, active_vertices) {
  if (is.null(edges)) {
    return(NA_integer_)
  }
  as.integer(nrow(grip.geodesic.misf.filter.edge.matrix(edges, active_vertices)))
}

grip.layout.trace.refinement.row <- function(meta, level_rows) {
  if (!length(level_rows)) {
    return(NA_integer_)
  }
  phase <- as.character(meta$phase[level_rows])
  preferred <- level_rows[phase %in% c("round", "level_start", "init")]
  if (length(preferred)) {
    return(as.integer(utils::tail(preferred, 1L)))
  }
  as.integer(utils::tail(level_rows, 1L))
}

grip.layout.trace.final.rows <- function(meta) {
  lgkk.rows <- which(meta$phase == "lgkk")
  if (length(lgkk.rows)) {
    final.rows <- which(meta$phase == "final")
    return(unique(as.integer(c(utils::tail(final.rows, 1L), lgkk.rows))))
  }
  final.rows <- which(meta$phase == "final")
  if (length(final.rows)) {
    return(as.integer(utils::tail(final.rows, 1L)))
  }
  as.integer(nrow(meta))
}

grip.layout.trace.as.stage.bundle <- function(trace, edges = NULL) {
  meta <- if (is.null(trace$meta)) data.frame() else {
    as.data.frame(trace$meta, stringsAsFactors = FALSE)
  }
  if (!nrow(meta) || is.null(trace$frames) || !length(trace$frames)) {
    return(list(
      stage_trace = grip.geodesic.misf.empty.stage.trace(),
      stage_data = list()
    ))
  }

  if (!is.null(edges)) {
    edges <- as.matrix(edges)
    storage.mode(edges) <- "integer"
  }

  levels <- sort(unique(as.integer(meta$misf_level)), decreasing = TRUE)
  top.level <- as.integer(levels[[1L]])
  top.rows <- which(as.integer(meta$misf_level) == top.level)
  top.init.row <- as.integer(top.rows[[1L]])
  top.init.coords <- grip.layout.trace.frame.at(trace, top.init.row)
  top.active.vertices <- grip.layout.trace.active.vertices(top.init.coords)
  top.trace.rows <- top.rows[meta$phase[top.rows] %in% c("init", "round")]
  if (!length(top.trace.rows)) {
    top.trace.rows <- top.rows
  }
  top.refinement.row <- grip.layout.trace.refinement.row(meta, top.rows)
  top.refinement.coords <- grip.layout.trace.frame.at(trace, top.refinement.row)
  top.refinement.vertices <- grip.layout.trace.active.vertices(top.refinement.coords)

  records <- list()
  records[[length(records) + 1L]] <- grip.geodesic.misf.new.stage.record(
    stage = "seed",
    level = top.level,
    method_family = "grip",
    coords_full = top.init.coords,
    active_vertices = top.active.vertices,
    trace = meta[top.init.row, , drop = FALSE],
    frames = list(top.init.coords),
    summary = list(
      active_n = length(top.active.vertices),
      pair_n = grip.layout.trace.induced.edge.count(edges, top.active.vertices),
      trace_rows = 1L,
      frame_count = 1L
    )
  )
  records[[length(records) + 1L]] <- grip.geodesic.misf.new.stage.record(
    stage = "initial_placement",
    level = top.level,
    method_family = "grip",
    coords_full = top.init.coords,
    active_vertices = top.active.vertices,
    trace = meta[top.init.row, , drop = FALSE],
    frames = list(top.init.coords),
    summary = list(
      active_n = length(top.active.vertices),
      pair_n = grip.layout.trace.induced.edge.count(edges, top.active.vertices),
      trace_rows = 1L,
      frame_count = 1L
    )
  )
  records[[length(records) + 1L]] <- grip.geodesic.misf.new.stage.record(
    stage = "top_level",
    level = top.level,
    method_family = "grip",
    coords_full = top.refinement.coords,
    active_vertices = top.refinement.vertices,
    trace = meta[top.trace.rows, , drop = FALSE],
    frames = grip.layout.trace.frames.at(trace, top.trace.rows),
    summary = list(
      active_n = length(top.refinement.vertices),
      pair_n = grip.layout.trace.induced.edge.count(edges, top.refinement.vertices),
      trace_rows = length(top.trace.rows),
      frame_count = length(top.trace.rows)
    )
  )

  prev.active.vertices <- top.refinement.vertices
  finer.levels <- levels[-1L]
  if (length(finer.levels)) {
    for (level in finer.levels) {
      level <- as.integer(level)
      level.rows <- which(as.integer(meta$misf_level) == level)
      if (!length(level.rows)) {
        next
      }
      insertion.rows <- level.rows[meta$phase[level.rows] == "level_start"]
      insertion.row <- if (length(insertion.rows)) {
        as.integer(insertion.rows[[1L]])
      } else {
        as.integer(level.rows[[1L]])
      }
      insertion.coords <- grip.layout.trace.frame.at(trace, insertion.row)
      insertion.vertices <- grip.layout.trace.active.vertices(insertion.coords)
      inserted.vertices <- setdiff(insertion.vertices, prev.active.vertices)

      refinement.row <- grip.layout.trace.refinement.row(meta, level.rows)
      refinement.rows <- level.rows[level.rows >= insertion.row & level.rows <= refinement.row]
      if (!length(refinement.rows)) {
        refinement.rows <- insertion.row
      }
      refinement.coords <- grip.layout.trace.frame.at(trace, refinement.row)
      refinement.vertices <- grip.layout.trace.active.vertices(refinement.coords)

      records[[length(records) + 1L]] <- grip.geodesic.misf.new.stage.record(
        stage = "insertion",
        level = level,
        method_family = "grip",
        coords_full = insertion.coords,
        active_vertices = insertion.vertices,
        inserted_vertices = inserted.vertices,
        trace = meta[insertion.row, , drop = FALSE],
        frames = list(insertion.coords),
        summary = list(
          active_n = length(insertion.vertices),
          inserted_n = length(inserted.vertices),
          pair_n = grip.layout.trace.induced.edge.count(edges, insertion.vertices),
          trace_rows = 1L,
          frame_count = 1L
        )
      )
      records[[length(records) + 1L]] <- grip.geodesic.misf.new.stage.record(
        stage = "refinement",
        level = level,
        method_family = "grip",
        coords_full = refinement.coords,
        active_vertices = refinement.vertices,
        trace = meta[refinement.rows, , drop = FALSE],
        frames = grip.layout.trace.frames.at(trace, refinement.rows),
        summary = list(
          active_n = length(refinement.vertices),
          pair_n = grip.layout.trace.induced.edge.count(edges, refinement.vertices),
          trace_rows = length(refinement.rows),
          frame_count = length(refinement.rows)
        )
      )

      prev.active.vertices <- refinement.vertices
    }
  }

  final.rows <- grip.layout.trace.final.rows(meta)
  final.row <- as.integer(utils::tail(final.rows, 1L))
  final.coords <- grip.layout.trace.frame.at(trace, final.row)
  final.vertices <- grip.layout.trace.active.vertices(final.coords)
  records[[length(records) + 1L]] <- grip.geodesic.misf.new.stage.record(
    stage = "final_polish",
    level = 0L,
    method_family = "grip",
    coords_full = final.coords,
    active_vertices = final.vertices,
    pair_mode = if (any(meta$phase[final.rows] == "lgkk")) "landmark" else NA_character_,
    trace = meta[final.rows, , drop = FALSE],
    frames = grip.layout.trace.frames.at(trace, final.rows),
    summary = list(
      active_n = length(final.vertices),
      pair_n = grip.layout.trace.induced.edge.count(edges, final.vertices),
      trace_rows = length(final.rows),
      frame_count = length(final.rows)
    )
  )

  names(records) <- vapply(records, `[[`, character(1L), "stage_key")
  list(
    stage_trace = grip.geodesic.misf.stage.trace.from.records(records),
    stage_data = records
  )
}
