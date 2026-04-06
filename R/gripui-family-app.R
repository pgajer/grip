gripui.family.theme <- function() {
  bslib::bs_theme(
    version = 5,
    base_font = bslib::font_google("Space Grotesk"),
    heading_font = bslib::font_google("Fraunces"),
    code_font = bslib::font_google("IBM Plex Mono"),
    bg = "#f4f1ea",
    fg = "#17323a",
    primary = "#1f3b73",
    secondary = "#8a5a44",
    "border-radius" = "0.9rem",
    "card-border-radius" = "0.9rem",
    "btn-border-radius" = "999px"
  )
}

gripui.require.family.app.packages <- function() {
  old <- gripui.enable.rgl.null.device()
  pkgs <- c("shiny", "bslib", "rgl")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing) > 0L) {
    options(rgl.useNULL = old)
    stop(
      "gripui_family_app requires optional packages that are not installed: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  old
}

gripui.family.validate.catalog <- function(catalog) {
  if (!is.list(catalog) || length(catalog) == 0L || is.null(names(catalog))) {
    stop("catalog must be a named non-empty list of family descriptors.", call. = FALSE)
  }
  required <- c("id", "label", "category", "function_name", "summary", "params", "builder", "presets", "implementation")
  for (nm in names(catalog)) {
    desc <- catalog[[nm]]
    missing <- setdiff(required, names(desc))
    if (length(missing) > 0L) {
      stop(sprintf("catalog entry '%s' is missing fields: %s", nm, paste(missing, collapse = ", ")), call. = FALSE)
    }
  }
  invisible(catalog)
}

gripui.family.categories <- function(catalog) {
  unique(vapply(catalog, `[[`, character(1L), "category"))
}

gripui.family.choices <- function(catalog, category = NULL) {
  keep <- rep(TRUE, length(catalog))
  if (!is.null(category) && nzchar(category) && !identical(category, "All families")) {
    keep <- vapply(catalog, function(x) identical(x$category, category), logical(1L))
  }
  items <- catalog[keep]
  stats::setNames(names(items), vapply(items, `[[`, character(1L), "label"))
}

gripui.family.preset.choices <- function(desc) {
  choices <- c(Default = "default")
  if (length(desc$presets) == 0L) {
    return(choices)
  }
  labels <- gsub("_", " ", names(desc$presets), fixed = TRUE)
  labels <- tools::toTitleCase(labels)
  c(choices, stats::setNames(names(desc$presets), labels))
}

gripui.family.input_id <- function(id) {
  paste0("family_param_", id)
}

gripui.family.numeric.vector.text <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return("")
  }
  paste(format(as.numeric(x), trim = TRUE, digits = 6), collapse = ", ")
}

gripui.family.value_or_default <- function(x, default) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) {
    return(default)
  }
  x
}

gripui.family.is.stochastic <- function(desc) {
  is.list(desc) && isTRUE(desc$stochastic)
}

gripui.family.param.spec <- function(desc, id) {
  idx <- which(vapply(desc$params, function(spec) identical(spec$id, id), logical(1L)))
  if (length(idx) == 0L) {
    return(NULL)
  }
  desc$params[[idx[[1L]]]]
}

gripui.family.seed.spec <- function(desc) {
  gripui.family.param.spec(desc, "seed")
}

gripui.family.resample.seed <- function(current_seed = NULL, max_seed = 1000000L) {
  max_seed <- max(as.integer(gripui.family.value_or_default(max_seed, 1000000L)), 0L)
  modulus <- as.numeric(max_seed) + 1
  now_ms <- floor(as.numeric(Sys.time()) * 1000)
  seed <- as.integer(now_ms %% modulus)
  if (!is.null(current_seed) && !is.na(current_seed) && identical(seed, as.integer(current_seed)) && max_seed > 0L) {
    seed <- as.integer((as.numeric(seed) + 1) %% modulus)
  }
  seed
}

gripui.family.filename.token <- function(x) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) {
    return("na")
  }
  if (is.logical(x)) {
    return(if (isTRUE(x[[1L]])) "true" else "false")
  }
  if (is.numeric(x)) {
    x <- format(as.numeric(x), scientific = FALSE, trim = TRUE)
  } else {
    x <- as.character(x)
  }
  x <- paste(x, collapse = "-")
  x <- tolower(gsub("[^A-Za-z0-9]+", "-", x))
  x <- gsub("^-+|-+$", "", x)
  if (nzchar(x)) x else "na"
}

gripui.family.save.key.ids <- function(values) {
  preferred <- c("n", "k", "seed", "surface", "graph_space", "normalize", "xmin", "xmax", "ymin", "ymax")
  ids <- intersect(preferred, names(values))
  if (length(ids) > 0L) {
    return(ids)
  }
  names(values)[seq_len(min(4L, length(values)))]
}

gripui.family.default.save.dir <- function(payload, root = getwd()) {
  file.path(root, "tmp", "gripui-family-graphs", payload$family_id)
}

gripui.family.save.stub <- function(payload, timestamp = Sys.time()) {
  ids <- gripui.family.save.key.ids(payload$values)
  key_bits <- vapply(ids, function(id) {
    paste0(id, "-", gripui.family.filename.token(payload$values[[id]]))
  }, character(1L))
  stamp <- format(as.POSIXct(timestamp, tz = Sys.timezone()), "%Y%m%d-%H%M%S")
  paste(c(payload$family_id, key_bits, stamp), collapse = "__")
}

gripui.family.save.path <- function(payload, root = getwd(), timestamp = Sys.time()) {
  file.path(
    gripui.family.default.save.dir(payload, root = root),
    paste0(gripui.family.save.stub(payload, timestamp = timestamp), ".rds")
  )
}

gripui.family.unique.save.path <- function(path) {
  if (!file.exists(path)) {
    return(path)
  }
  dir <- dirname(path)
  stem <- sub("\\.[Rr][Dd][Ss]$", "", basename(path))
  ext <- tools::file_ext(path)
  counter <- 1L
  repeat {
    candidate <- file.path(dir, sprintf("%s__%02d.%s", stem, counter, ext))
    if (!file.exists(candidate)) {
      return(candidate)
    }
    counter <- counter + 1L
  }
}

gripui.family.save.bundle <- function(payload, path, saved_at = Sys.time()) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  bundle <- list(
    family_id = payload$family_id,
    family_label = payload$family_label,
    category = payload$category,
    values = payload$values,
    rendered_at = payload$rendered_at,
    saved_at = as.POSIXct(saved_at, tz = Sys.timezone()),
    code = payload$code,
    payload = payload,
    session = list(
      grip_version = tryCatch(as.character(utils::packageVersion("grip")), error = function(e) NA_character_),
      working_directory = getwd()
    )
  )
  saveRDS(bundle, path)
  invisible(path)
}

gripui.family.saved.bundle.files <- function(family_id, root = getwd()) {
  dir <- file.path(root, "tmp", "gripui-family-graphs", family_id)
  if (!dir.exists(dir)) {
    return(character(0L))
  }
  paths <- list.files(dir, pattern = "\\.[Rr][Dd][Ss]$", full.names = TRUE)
  if (length(paths) == 0L) {
    return(character(0L))
  }
  info <- file.info(paths)
  paths[order(info$mtime, decreasing = TRUE, na.last = TRUE)]
}

gripui.family.saved.bundle.choices <- function(desc, root = getwd()) {
  paths <- gripui.family.saved.bundle.files(desc$id, root = root)
  if (length(paths) == 0L) {
    return(character(0L))
  }
  labels <- sub("\\.[Rr][Dd][Ss]$", "", basename(paths))
  stats::setNames(paths, labels)
}

gripui.family.saved.load.mode.choices <- function(desc) {
  if (identical(desc$id, "sampled_rectangle")) {
    return(c(
      "Exact saved graph" = "exact",
      "Reuse sample, keep saved topology" = "sample_topology",
      "Reuse sample and rebuild iKNN" = "sample_rebuild"
    ))
  }
  c("Exact saved graph" = "exact")
}

gripui.family.saved.load.mode.help <- function(desc) {
  if (!identical(desc$id, "sampled_rectangle")) {
    return("Exact saved graph restores the saved topology, geometry, and weights.")
  }
  paste(
    "Exact restores the saved graph bundle as rendered.",
    "Keep saved topology reuses the saved sample points and edges, then recomputes 3D weights under the current surface controls.",
    "Rebuild iKNN reuses the saved sample points and rebuilds topology and weights from the current controls."
  )
}

gripui.family.read.saved.bundle <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop("Choose an existing saved graph bundle before loading.", call. = FALSE)
  }
  bundle <- readRDS(path)
  if (!is.list(bundle)) {
    stop("Saved graph bundle is not a list.", call. = FALSE)
  }
  bundle
}

gripui.family.saved.bundle.values <- function(bundle) {
  values <- bundle$values
  if (is.null(values) && is.list(bundle$payload)) {
    values <- bundle$payload$values
  }
  if (is.null(values)) {
    stop("Saved graph bundle does not contain parameter values.", call. = FALSE)
  }
  values
}

gripui.family.saved.bundle.raw <- function(bundle) {
  raw <- NULL
  if (is.list(bundle$payload)) {
    raw <- bundle$payload$raw
  }
  if (is.null(raw)) {
    stop("Saved graph bundle does not contain raw graph data.", call. = FALSE)
  }
  raw
}

gripui.family.sample.reuse.values <- function(current_values, saved_values) {
  out <- current_values
  fixed_ids <- c("n", "xmin", "xmax", "ymin", "ymax", "seed")
  for (id in intersect(fixed_ids, names(saved_values))) {
    out[[id]] <- saved_values[[id]]
  }
  out
}

gripui.family.payload.from.raw <- function(desc, values, raw) {
  display <- gripui.family.display.coords(raw)
  display_coords <- display$coords
  plot_coords <- gripui.family.plot.coords(raw, display_coords)
  graph_obj <- grip.build.adj.from.edges(raw$edges, n = raw$n, edge_weights = raw$edge_weights)
  graph_obj$vertex_data <- gripui.family.vertex.data(raw, graph_obj, display_coords, plot_coords)
  code <- if (is.function(desc$code)) {
    desc$code(values)
  } else {
    .gripui.family.call.code(desc$function_name, values)
  }
  list(
    family_id = desc$id,
    family_label = desc$label,
    category = desc$category,
    values = values,
    raw = raw,
    graph = graph_obj,
    edges = raw$edges,
    n = raw$n,
    edge_weights = raw$edge_weights,
    coords_display = display_coords,
    coords_plot = plot_coords,
    code = code,
    implementation = desc$implementation,
    rendered_at = Sys.time(),
    note = display$note
  )
}

gripui.family.update.param.inputs <- function(session, desc, values) {
  for (spec in desc$params) {
    if (!spec$id %in% names(values)) {
      next
    }
    input_id <- gripui.family.input_id(spec$id)
    value <- values[[spec$id]]
    switch(
      spec$type,
      integer = shiny::updateNumericInput(session, input_id, value = value),
      double = shiny::updateNumericInput(session, input_id, value = value),
      choice = shiny::updateSelectInput(session, input_id, selected = value),
      logical = shiny::updateCheckboxInput(session, input_id, value = isTRUE(value)),
      numeric_vector = shiny::updateTextInput(session, input_id, value = gripui.family.numeric.vector.text(value)),
      stop("Unsupported parameter type: ", spec$type, call. = FALSE)
    )
  }
  invisible(values)
}

gripui.family.sampled.rectangle.raw.from.bundle <- function(desc, current_values, bundle, mode) {
  saved_values <- gripui.family.saved.bundle.values(bundle)
  saved_raw <- gripui.family.saved.bundle.raw(bundle)
  if (!identical(desc$id, "sampled_rectangle")) {
    stop("Sample reuse is currently supported only for the sampled rectangle family.", call. = FALSE)
  }
  values <- gripui.family.sample.reuse.values(current_values, saved_values)
  if (identical(mode, "sample_topology")) {
    topology_ids <- c(
      "k",
      "graph_space",
      "max.path.edge.ratio.deviation.thld",
      "path.edge.ratio.percentile",
      "threshold.percentile"
    )
    for (id in intersect(topology_ids, names(saved_values))) {
      values[[id]] <- saved_values[[id]]
    }
  }
  raw <- switch(
    mode,
    sample_topology = .sampled.rectangle.surface.graph.reweight.saved.topology(
      graph = saved_raw,
      surface = values$surface,
      amplitude = values$amplitude,
      freq_u = values$freq_u,
      freq_v = values$freq_v,
      normalize = values$normalize
    ),
    sample_rebuild = .sampled.rectangle.surface.graph.from.coords(
      coords_param = saved_raw$coords_param,
      k = values$k,
      xmin = values$xmin,
      xmax = values$xmax,
      ymin = values$ymin,
      ymax = values$ymax,
      seed = values$seed,
      surface = values$surface,
      amplitude = values$amplitude,
      freq_u = values$freq_u,
      freq_v = values$freq_v,
      graph_space = values$graph_space,
      max.path.edge.ratio.deviation.thld = values$max.path.edge.ratio.deviation.thld,
      path.edge.ratio.percentile = values$path.edge.ratio.percentile,
      threshold.percentile = values$threshold.percentile,
      normalize = values$normalize
    ),
    stop("Unsupported load mode: ", mode, call. = FALSE)
  )
  list(values = values, raw = raw)
}

gripui.family.js.value <- function(x) {
  if (is.logical(x)) {
    if (isTRUE(x)) "true" else "false"
  } else if (is.numeric(x)) {
    format(x[[1L]], scientific = FALSE, trim = TRUE)
  } else {
    paste0("'", gsub("'", "\\\\'", as.character(x[[1L]]), fixed = TRUE), "'")
  }
}

gripui.family.visible.condition <- function(spec) {
  vis <- spec$visible_if
  if (is.null(vis) || !length(vis)) {
    return(NULL)
  }
  pieces <- vapply(names(vis), function(id) {
    values <- vis[[id]]
    input_name <- gripui.family.input_id(id)
    options <- vapply(values, function(val) {
      sprintf("input['%s'] == %s", input_name, gripui.family.js.value(val))
    }, character(1L))
    sprintf("(%s)", paste(options, collapse = " || "))
  }, character(1L))
  paste(pieces, collapse = " && ")
}

gripui.family.param.control <- function(spec, value) {
  input_id <- gripui.family.input_id(spec$id)
  widget <- switch(
    spec$type,
    integer = shiny::numericInput(
      inputId = input_id,
      label = spec$label,
      value = value,
      min = spec$min,
      max = spec$max,
      step = spec$step
    ),
    double = shiny::numericInput(
      inputId = input_id,
      label = spec$label,
      value = value,
      min = spec$min,
      max = spec$max,
      step = spec$step
    ),
    choice = shiny::selectInput(
      inputId = input_id,
      label = spec$label,
      choices = stats::setNames(spec$choices, tools::toTitleCase(gsub("_", " ", spec$choices, fixed = TRUE))),
      selected = value
    ),
    logical = shiny::checkboxInput(
      inputId = input_id,
      label = spec$label,
      value = isTRUE(value)
    ),
    numeric_vector = shiny::textInput(
      inputId = input_id,
      label = spec$label,
      value = gripui.family.numeric.vector.text(value),
      placeholder = "e.g. 1, 0.8, 0.6"
    ),
    stop("Unsupported parameter type: ", spec$type, call. = FALSE)
  )

  block <- shiny::tagList(
    widget,
    if (!is.null(spec$help) && nzchar(spec$help)) {
      shiny::tags$div(style = "margin-top:-0.5rem;color:#6b7280;font-size:0.9rem;", spec$help)
    }
  )
  cond <- gripui.family.visible.condition(spec)
  if (!is.null(cond)) {
    shiny::conditionalPanel(condition = cond, block)
  } else {
    block
  }
}

gripui.family.param.ui <- function(desc, values) {
  groups <- unique(vapply(desc$params, `[[`, character(1L), "group"))
  shiny::tagList(lapply(groups, function(group) {
    items <- Filter(function(x) identical(x$group, group), desc$params)
    shiny::tags$div(
      class = "gripui-family-group",
      shiny::tags$h5(style = "margin-top:0.8rem;", group),
      shiny::tagList(lapply(items, function(spec) {
        gripui.family.param.control(spec, values[[spec$id]])
      }))
    )
  }))
}

gripui.family.collect.values <- function(desc, input, preset_id = "default") {
  defaults <- .gripui.family.merge.values(desc, preset_id = preset_id)
  out <- defaults
  for (spec in desc$params) {
    raw <- input[[gripui.family.input_id(spec$id)]]
    value <- if (is.null(raw)) spec$default else spec$coerce(raw)
    out[[spec$id]] <- value
  }
  out
}

gripui.family.as_matrix <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  out <- as.matrix(x)
  storage.mode(out) <- "double"
  out
}

gripui.family.pad_coords <- function(coords, target_cols = 3L, colnames_out = c("x", "y", "z")) {
  coords <- gripui.family.as_matrix(coords)
  if (is.null(coords)) {
    return(NULL)
  }
  if (ncol(coords) < target_cols) {
    coords <- cbind(coords, matrix(0, nrow(coords), target_cols - ncol(coords)))
  }
  coords <- coords[, seq_len(target_cols), drop = FALSE]
  colnames(coords) <- colnames_out[seq_len(target_cols)]
  coords
}

gripui.family.tree.leaf_order <- function(parent) {
  children <- split(seq_along(parent), parent)
  children <- children[names(children) != "0"]

  order <- numeric(length(parent))
  next_leaf <- 1

  assign_order <- function(node) {
    kids <- children[[as.character(node)]]
    if (is.null(kids) || length(kids) == 0L) {
      order[[node]] <<- next_leaf
      next_leaf <<- next_leaf + 1
      return(order[[node]])
    }
    child_vals <- vapply(kids, assign_order, numeric(1L))
    order[[node]] <<- mean(child_vals)
    order[[node]]
  }

  assign_order(1L)
  order
}

gripui.family.tree.display.coords <- function(raw) {
  n <- as.integer(raw$n)
  parent <- as.integer(raw$parent)
  if (length(parent) != n) {
    return(cbind(seq_len(n), rep(0, n), rep(0, n)))
  }
  depth <- as.integer(raw$vertex_depth)
  leaf_order <- gripui.family.tree.leaf_order(parent)
  edge_table <- raw$edge_table

  root_dist <- numeric(n)
  branch_z <- numeric(n)
  if (is.data.frame(edge_table) && nrow(edge_table) > 0L) {
    for (i in seq_len(nrow(edge_table))) {
      child <- edge_table$child[[i]]
      par <- edge_table$parent[[i]]
      w <- edge_table$edge_weight[[i]]
      idx <- edge_table$branch_index[[i]]
      k_local <- max(edge_table$branch_index[edge_table$parent == par])
      offset <- idx - (k_local + 1) / 2
      root_dist[[child]] <- root_dist[[par]] + w
      branch_z[[child]] <- branch_z[[par]] * 0.45 + offset * w * 0.35
    }
  }

  y <- leaf_order - mean(leaf_order)
  if (length(y) > 1L && max(abs(y)) > 0) {
    y <- y / max(abs(y))
  }

  coords <- cbind(root_dist, y, branch_z)
  colnames(coords) <- c("x", "y", "z")
  coords
}

gripui.family.display.coords <- function(raw) {
  coords <- gripui.family.pad_coords(raw$coords_surface, target_cols = 3L)
  note <- NULL
  if (!is.null(coords)) {
    return(list(coords = coords, note = note))
  }
  coords <- gripui.family.tree.display.coords(raw)
  note <- paste(
    "This family has intrinsic weighted-tree geometry but no canonical ambient embedding.",
    "The viewer uses a deterministic display layout for inspection."
  )
  list(coords = coords, note = note)
}

gripui.family.plot.coords <- function(raw, display_coords) {
  coords <- raw$coords_param
  if (is.null(coords)) {
    coords <- display_coords
  }
  gripui.family.pad_coords(coords, target_cols = 2L, colnames_out = c("x", "y"))
}

gripui.family.vertex.data <- function(raw, graph_obj, display_coords, plot_coords) {
  degree <- vapply(graph_obj$adj_list, length, integer(1L))
  out <- data.frame(
    x = display_coords[, 1L],
    y = display_coords[, 2L],
    z = display_coords[, 3L],
    plot_x = plot_coords[, 1L],
    plot_y = plot_coords[, 2L],
    degree = degree,
    stringsAsFactors = FALSE
  )
  if (!is.null(raw$coords_surface)) {
    surface <- gripui.family.pad_coords(raw$coords_surface, target_cols = 3L)
    out$surface_x <- surface[, 1L]
    out$surface_y <- surface[, 2L]
    out$surface_z <- surface[, 3L]
  }
  if (!is.null(raw$coords_param)) {
    param <- gripui.family.as_matrix(raw$coords_param)
    for (j in seq_len(ncol(param))) {
      out[[paste0("param_", j)]] <- param[, j]
    }
  }
  if (!is.null(raw$vertex_depth) && length(raw$vertex_depth) == nrow(out)) {
    out$depth <- raw$vertex_depth
  }
  out
}

gripui.family.count.edges <- function(edges) {
  if (is.null(edges)) 0L else nrow(as.matrix(edges))
}

gripui.family.build.payload <- function(desc, values) {
  raw <- desc$builder(values)
  gripui.family.payload.from.raw(desc, values, raw)
}

gripui.family.summary.table <- function(items) {
  data.frame(
    field = names(items),
    value = unname(vapply(items, as.character, character(1L))),
    stringsAsFactors = FALSE
  )
}

gripui.family.graph.summary <- function(payload) {
  degrees <- vapply(payload$graph$adj_list, length, integer(1L))
  gripui.family.summary.table(list(
    family = payload$family_label,
    vertices = payload$n,
    edges = gripui.family.count.edges(payload$edges),
    average_degree = sprintf("%.3f", mean(degrees)),
    max_degree = max(degrees),
    implementation = basename(payload$implementation)
  ))
}

gripui.family.weight.summary <- function(payload) {
  w <- as.numeric(payload$edge_weights)
  cv <- if (length(w) > 1L && mean(w) > 0) stats::sd(w) / mean(w) else 0
  gripui.family.summary.table(list(
    min = sprintf("%.4f", min(w)),
    median = sprintf("%.4f", stats::median(w)),
    mean = sprintf("%.4f", mean(w)),
    max = sprintf("%.4f", max(w)),
    cv = sprintf("%.4f", cv),
    normalize = payload$raw$normalize[[1L]]
  ))
}

gripui.family.geometry.summary <- function(payload) {
  coords <- payload$coords_display
  spans <- apply(coords, 2L, function(x) diff(range(x, finite = TRUE)))
  gripui.family.summary.table(list(
    display_x_span = sprintf("%.4f", spans[[1L]]),
    display_y_span = sprintf("%.4f", spans[[2L]]),
    display_z_span = sprintf("%.4f", spans[[3L]]),
    plot_dimensions = ncol(payload$coords_plot),
    surface = if ("surface" %in% names(payload$raw)) payload$raw$surface[[1L]] else "intrinsic",
    note = if (is.null(payload$note)) "" else payload$note
  ))
}

gripui.family.color.choices <- function(payload) {
  if (is.null(payload) || is.null(payload$graph$vertex_data)) {
    return(c(Plain = "plain"))
  }
  fields <- setdiff(names(payload$graph$vertex_data), c("x", "y", "z"))
  c(Plain = "plain", stats::setNames(fields, tools::toTitleCase(gsub("_", " ", fields, fixed = TRUE))))
}

gripui.family.default.color <- function(choices) {
  preferred <- c("surface_z", "depth", "degree", "plain")
  for (nm in preferred) {
    if (nm %in% unname(choices)) {
      return(nm)
    }
  }
  unname(choices[[1L]])
}

gripui.family.compare.family.input_id <- function(idx) {
  sprintf("compare_family_%d", idx)
}

gripui.family.compare.preset.input_id <- function(idx) {
  sprintf("compare_preset_%d", idx)
}

gripui.family.compare.plot.output_id <- function(idx) {
  sprintf("compare_viewer_2d_%d", idx)
}

gripui.family.compare.preferred.ids <- function(catalog) {
  preferred <- c(
    "mesh",
    "sierpinski_carpet",
    "sphere",
    "irregular_torus",
    "menger_sponge",
    "kary_tree"
  )
  unique(c(preferred[preferred %in% names(catalog)], names(catalog)))
}

gripui.family.compare.default.family_id <- function(catalog, idx = 1L) {
  ids <- gripui.family.compare.preferred.ids(catalog)
  ids[[((idx - 1L) %% length(ids)) + 1L]]
}

gripui.family.compare.default.preset <- function(desc, rank = 1L) {
  ids <- c("default", names(desc$presets))
  ids[[((rank - 1L) %% length(ids)) + 1L]]
}

gripui.family.compare.slot.selection <- function(idx,
                                                 input,
                                                 catalog,
                                                 current_family_id,
                                                 lock_family = FALSE,
                                                 include_current = TRUE) {
  if (isTRUE(include_current) && idx == 1L) {
    return(list(source = "current", slot = idx))
  }

  family_id <- if (isTRUE(lock_family)) {
    current_family_id
  } else {
    raw <- input[[gripui.family.compare.family.input_id(idx)]]
    if (is.null(raw) || !raw %in% names(catalog)) {
      gripui.family.compare.default.family_id(catalog, idx = idx)
    } else {
      raw
    }
  }

  desc <- catalog[[family_id]]
  preset_raw <- input[[gripui.family.compare.preset.input_id(idx)]]
  preset_id <- if (is.null(preset_raw) || !preset_raw %in% c("default", names(desc$presets))) {
    gripui.family.compare.default.preset(desc, rank = idx)
  } else {
    preset_raw
  }

  list(
    source = "preset",
    slot = idx,
    family_id = family_id,
    preset_id = preset_id
  )
}

gripui.family.compare.slot.ui <- function(idx,
                                          input,
                                          catalog,
                                          current_family_id,
                                          include_current = TRUE,
                                          lock_family = FALSE,
                                          current_desc = NULL) {
  selection <- gripui.family.compare.slot.selection(
    idx = idx,
    input = input,
    catalog = catalog,
    current_family_id = current_family_id,
    lock_family = lock_family,
    include_current = include_current
  )

  if (identical(selection$source, "current")) {
    label <- if (is.null(current_desc)) "Current explore controls" else current_desc$label
    return(
      shiny::tags$div(
        class = "gripui-family-group",
        shiny::tags$h5(style = "margin-top:0.8rem;", "Current Explore Selection"),
        shiny::tags$p(
          class = "gripui-project-summary",
          sprintf("Slot %d uses the current Explore family and the current sidebar parameter values.", idx)
        ),
        shiny::tags$p(class = "gripui-selection-status", label)
      )
    )
  }

  desc <- catalog[[selection$family_id]]
  family_choices <- gripui.family.choices(catalog)

  shiny::tags$div(
    class = "gripui-family-group",
    shiny::tags$h5(style = "margin-top:0.8rem;", sprintf("Variant %d", idx)),
    if (isTRUE(lock_family)) {
      shiny::tags$p(
        class = "gripui-selection-status",
        sprintf("Family locked to %s.", catalog[[current_family_id]]$label)
      )
    } else {
      shiny::selectInput(
        inputId = gripui.family.compare.family.input_id(idx),
        label = "Family",
        choices = family_choices,
        selected = selection$family_id
      )
    },
    shiny::selectInput(
      inputId = gripui.family.compare.preset.input_id(idx),
      label = "Preset",
      choices = gripui.family.preset.choices(desc),
      selected = selection$preset_id
    )
  )
}

gripui.family.compare.color.choices <- function(payloads) {
  if (is.null(payloads) || length(payloads) == 0L) {
    return(c(Plain = "plain"))
  }
  fields <- unique(unlist(lapply(payloads, function(payload) {
    if (is.null(payload$graph$vertex_data)) {
      return(character(0))
    }
    setdiff(names(payload$graph$vertex_data), c("x", "y", "z"))
  }), use.names = FALSE))
  c(Plain = "plain", stats::setNames(fields, tools::toTitleCase(gsub("_", " ", fields, fixed = TRUE))))
}

gripui.family.compare.summary <- function(payloads) {
  if (is.null(payloads) || length(payloads) == 0L) {
    return(data.frame())
  }
  rows <- lapply(payloads, function(payload) {
    weights <- as.numeric(payload$edge_weights)
    coords <- payload$coords_display
    spans <- apply(coords, 2L, function(x) diff(range(x, finite = TRUE)))
    data.frame(
      slot = payload$compare_slot,
      label = payload$compare_label,
      family = payload$family_label,
      preset = payload$compare_preset,
      vertices = payload$n,
      edges = gripui.family.count.edges(payload$edges),
      mean_weight = sprintf("%.4f", mean(weights)),
      weight_cv = sprintf("%.4f", if (length(weights) > 1L && mean(weights) > 0) stats::sd(weights) / mean(weights) else 0),
      z_span = sprintf("%.4f", spans[[3L]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

gripui.family.initial.explore.state <- function(catalog) {
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

gripui.family.ui <- function(catalog, title, subtitle = NULL) {
  categories <- c("All families", gripui.family.categories(catalog))
  initial_state <- gripui.family.initial.explore.state(catalog)
  css.path <- system.file("app/www/gripui.css", package = "grip")

  bslib::page_sidebar(
    title = shiny::div(
      class = "gripui-titlebar",
      shiny::div(class = "gripui-brand", title)
    ),
    theme = gripui.family.theme(),
    sidebar = bslib::sidebar(
      width = 380,
      if (!is.null(subtitle) && nzchar(subtitle)) {
        shiny::tags$p(class = "gripui-project-summary", subtitle)
      },
      shiny::conditionalPanel(
        condition = "input['family_app_mode'] == 'Explore'",
        shiny::uiOutput("family_meta"),
        shiny::selectInput("family_category", "Category", choices = categories, selected = initial_state$category),
        shiny::selectInput("family_id", "Family", choices = initial_state$choices, selected = initial_state$family_id),
        shiny::selectInput("family_preset", "Preset", choices = gripui.family.preset.choices(initial_state$desc), selected = "default"),
        shiny::uiOutput("family_param_panel"),
        shiny::selectInput("viewer_color_by", "Color by", choices = c(Plain = "plain"), selected = "plain"),
        shiny::checkboxInput("show_edges", "Show edges", value = TRUE),
        shiny::actionButton("render_family_geometry", "Render geometry", class = "btn-primary"),
        shiny::uiOutput("family_stochastic_actions"),
        shiny::uiOutput("render_status")
      ),
      shiny::conditionalPanel(
        condition = "input['family_app_mode'] == 'Compare'",
        shiny::tags$div(class = "gripui-project-title", "Compare"),
        shiny::tags$p(
          class = "gripui-project-summary",
          "Compare current controls and family presets side by side."
        ),
        shiny::numericInput("compare_panel_count", "Panels", value = 3L, min = 2L, max = 4L, step = 1L),
        shiny::checkboxInput("compare_include_current", "Include current Explore selection", value = TRUE),
        shiny::checkboxInput("compare_lock_family", "Lock family to current Explore family", value = FALSE),
        shiny::uiOutput("compare_control_panel"),
        shiny::selectInput("compare_color_by", "Color by", choices = c(Plain = "plain"), selected = "plain"),
        shiny::checkboxInput("compare_show_edges", "Show edges", value = TRUE),
        shiny::actionButton("render_family_compare", "Render comparison", class = "btn-primary"),
        shiny::uiOutput("compare_render_status")
      )
    ),
    if (nzchar(css.path)) shiny::tags$head(shiny::includeCSS(css.path)),
    bslib::navset_card_tab(
      id = "family_app_mode",
      bslib::nav_panel(
        "Explore",
        bslib::layout_columns(
          col_widths = c(7, 5),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header("3D Geometry"),
            shiny::uiOutput("family_viewer_3d")
          ),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header("2D Auxiliary View"),
            shiny::plotOutput("family_viewer_2d", height = 520)
          )
        ),
        bslib::layout_columns(
          col_widths = c(4, 4, 4),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header("Graph Summary"),
            shiny::tableOutput("family_graph_summary")
          ),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header("Weight Summary"),
            shiny::tableOutput("family_weight_summary")
          ),
          bslib::card(
            full_screen = TRUE,
            class = "gripui-card",
            bslib::card_header("Geometry Summary"),
            shiny::tableOutput("family_geometry_summary")
          )
        ),
        bslib::card(
          full_screen = TRUE,
          class = "gripui-card",
          bslib::card_header("Reproducible R Call"),
          shiny::verbatimTextOutput("family_code")
        )
      ),
      bslib::nav_panel(
        "Compare",
        shiny::uiOutput("compare_grid"),
        bslib::card(
          full_screen = TRUE,
          class = "gripui-card",
          bslib::card_header("Shared Summary"),
          shiny::tableOutput("compare_summary")
        )
      )
    )
  )
}

gripui.family.server <- function(catalog) {
  function(input, output, session) {
    payload_state <- shiny::reactiveVal(NULL)
    error_state <- shiny::reactiveVal(NULL)
    save_state <- shiny::reactiveVal(NULL)
    compare_payloads_state <- shiny::reactiveVal(NULL)
    compare_error_state <- shiny::reactiveVal(NULL)

    current_desc <- shiny::reactive({
      catalog[[input$family_id]]
    })

    shiny::observe({
      choices <- gripui.family.choices(catalog, input$family_category)
      current <- shiny::isolate(input$family_id)
      selected <- if (!is.null(current) && current %in% unname(choices)) current else unname(choices[[1L]])
      shiny::updateSelectInput(session, "family_id", choices = choices, selected = selected)
    })

    shiny::observe({
      desc <- current_desc()
      shiny::updateSelectInput(
        session,
        "family_preset",
        choices = gripui.family.preset.choices(desc),
        selected = "default"
      )
    })

    output$family_meta <- shiny::renderUI({
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

    output$family_param_panel <- shiny::renderUI({
      desc <- current_desc()
      values <- .gripui.family.merge.values(desc, preset_id = input$family_preset)
      gripui.family.param.ui(desc, values)
    })

    build_payload <- function(values) {
      desc <- current_desc()
      result <- tryCatch(
        gripui.family.build.payload(desc, values),
        error = function(e) e
      )
      if (inherits(result, "error")) {
        error_state(conditionMessage(result))
        shiny::showNotification(conditionMessage(result), type = "error")
        return(invisible(NULL))
      }
      payload_state(result)
      error_state(NULL)
      save_state(NULL)
      invisible(result)
    }

    shiny::observeEvent(list(input$family_id, input$family_preset), {
      desc <- current_desc()
      values <- .gripui.family.merge.values(desc, preset_id = input$family_preset)
      build_payload(values)
    }, ignoreInit = FALSE)

    shiny::observeEvent(input$render_family_geometry, {
      desc <- current_desc()
      values <- gripui.family.collect.values(desc, input, preset_id = input$family_preset)
      build_payload(values)
    })

    output$family_stochastic_actions <- shiny::renderUI({
      desc <- current_desc()
      if (!gripui.family.is.stochastic(desc)) {
        return(NULL)
      }
      save_state()
      save_dir <- file.path(getwd(), "tmp", "gripui-family-graphs", desc$id)
      saved_choices <- gripui.family.saved.bundle.choices(desc)
      selected_bundle <- if (!is.null(input$family_saved_bundle) &&
                             input$family_saved_bundle %in% unname(saved_choices)) {
        input$family_saved_bundle
      } else if (length(saved_choices) > 0L) {
        unname(saved_choices[[1L]])
      } else {
        NULL
      }
      mode_choices <- gripui.family.saved.load.mode.choices(desc)
      selected_mode <- if (!is.null(input$family_saved_load_mode) &&
                           input$family_saved_load_mode %in% unname(mode_choices)) {
        input$family_saved_load_mode
      } else {
        unname(mode_choices[[1L]])
      }
      shiny::tagList(
        shiny::tags$div(
          style = "display:flex;gap:0.5rem;flex-wrap:wrap;margin-top:0.5rem;",
          shiny::actionButton("resample_family_graph", "Resample graph"),
          shiny::actionButton("save_family_graph", "Save graph (.rds)")
        ),
        if (length(saved_choices) > 0L) {
          shiny::tagList(
            shiny::selectInput(
              "family_saved_bundle",
              "Saved graph",
              choices = saved_choices,
              selected = selected_bundle
            ),
            shiny::selectInput(
              "family_saved_load_mode",
              "Load mode",
              choices = mode_choices,
              selected = selected_mode
            ),
            shiny::actionButton("load_family_graph", "Load saved graph"),
            shiny::tags$div(
              style = "margin-top:0.5rem;color:#6b7280;font-size:0.9rem;",
              gripui.family.saved.load.mode.help(desc)
            )
          )
        } else {
          shiny::tags$div(
            style = "margin-top:0.5rem;color:#6b7280;font-size:0.9rem;",
            "No saved graph bundles yet for this family."
          )
        },
        shiny::tags$div(
          style = "margin-top:0.5rem;color:#6b7280;font-size:0.9rem;",
          sprintf("Resample assigns a fresh seed. Save writes a complete graph bundle under %s.", save_dir)
        )
      )
    })

    shiny::observeEvent(input$resample_family_graph, {
      desc <- current_desc()
      values <- gripui.family.collect.values(desc, input, preset_id = input$family_preset)
      seed_spec <- gripui.family.seed.spec(desc)
      seed_label <- NULL
      if (!is.null(seed_spec)) {
        values$seed <- seed_spec$coerce(gripui.family.resample.seed(
          current_seed = values$seed,
          max_seed = seed_spec$max
        ))
        shiny::updateNumericInput(session, gripui.family.input_id(seed_spec$id), value = values$seed)
        seed_label <- as.character(values$seed)
      }
      result <- build_payload(values)
      if (!is.null(result)) {
        shiny::showNotification(
          if (!is.null(seed_label)) {
            sprintf("Rendered a new %s sample with seed %s.", desc$label, seed_label)
          } else {
            sprintf("Rendered a new %s sample.", desc$label)
          },
          type = "message"
        )
      }
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$save_family_graph, {
      payload <- payload_state()
      if (is.null(payload)) {
        shiny::showNotification("Render a graph before saving it.", type = "error")
        return(invisible(NULL))
      }
      path <- gripui.family.unique.save.path(gripui.family.save.path(payload))
      result <- tryCatch(
        gripui.family.save.bundle(payload, path),
        error = function(e) e
      )
      if (inherits(result, "error")) {
        shiny::showNotification(conditionMessage(result), type = "error")
        return(invisible(NULL))
      }
      save_state(path)
      shiny::showNotification(sprintf("Saved graph bundle to %s", path), type = "message", duration = 6)
      invisible(path)
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$load_family_graph, {
      desc <- current_desc()
      path <- input$family_saved_bundle
      mode <- gripui.family.value_or_default(input$family_saved_load_mode, "exact")
      result <- tryCatch({
        bundle <- gripui.family.read.saved.bundle(path)
        bundle_values <- gripui.family.saved.bundle.values(bundle)
        if (!identical(gripui.family.value_or_default(bundle$family_id, bundle$payload$family_id), desc$id)) {
          stop("The selected saved graph belongs to a different family.", call. = FALSE)
        }
        built <- if (identical(mode, "exact")) {
          list(values = bundle_values, raw = gripui.family.saved.bundle.raw(bundle))
        } else {
          current_values <- gripui.family.collect.values(desc, input, preset_id = input$family_preset)
          gripui.family.sampled.rectangle.raw.from.bundle(
            desc = desc,
            current_values = current_values,
            bundle = bundle,
            mode = mode
          )
        }
        payload <- gripui.family.payload.from.raw(desc, built$values, built$raw)
        list(
          values = built$values,
          payload = payload,
          path = path,
          mode = mode
        )
      }, error = function(e) e)
      if (inherits(result, "error")) {
        shiny::showNotification(conditionMessage(result), type = "error")
        return(invisible(NULL))
      }
      gripui.family.update.param.inputs(session, desc, result$values)
      payload_state(result$payload)
      error_state(NULL)
      save_state(NULL)
      shiny::showNotification(
        sprintf("Loaded %s using %s.", basename(result$path), gsub("_", " ", result$mode, fixed = TRUE)),
        type = "message",
        duration = 6
      )
      invisible(result$payload)
    }, ignoreInit = TRUE)

    shiny::observe({
      payload <- payload_state()
      choices <- gripui.family.color.choices(payload)
      current <- shiny::isolate(input$viewer_color_by)
      selected <- if (!is.null(current) && current %in% unname(choices)) {
        current
      } else {
        gripui.family.default.color(choices)
      }
      shiny::updateSelectInput(session, "viewer_color_by", choices = choices, selected = selected)
    })

    output$compare_control_panel <- shiny::renderUI({
      count <- min(max(as.integer(gripui.family.value_or_default(input$compare_panel_count, 3L)), 2L), 4L)
      include_current <- isTRUE(input$compare_include_current)
      lock_family <- isTRUE(input$compare_lock_family)
      current_family_id <- if (!is.null(input$family_id) && input$family_id %in% names(catalog)) {
        input$family_id
      } else {
        names(catalog)[[1L]]
      }

      shiny::tagList(lapply(seq_len(count), function(idx) {
        gripui.family.compare.slot.ui(
          idx = idx,
          input = input,
          catalog = catalog,
          current_family_id = current_family_id,
          include_current = include_current,
          lock_family = lock_family,
          current_desc = current_desc()
        )
      }))
    })

    build_current_compare_payload <- function() {
      desc <- current_desc()
      values <- gripui.family.collect.values(desc, input, preset_id = input$family_preset)
      payload <- gripui.family.build.payload(desc, values)
      payload$compare_slot <- 1L
      payload$compare_label <- "Current"
      payload$compare_preset <- "current_controls"
      payload
    }

    build_compare_payloads <- function() {
      count <- min(max(as.integer(gripui.family.value_or_default(input$compare_panel_count, 3L)), 2L), 4L)
      include_current <- isTRUE(input$compare_include_current)
      lock_family <- isTRUE(input$compare_lock_family)
      current_family_id <- if (!is.null(input$family_id) && input$family_id %in% names(catalog)) {
        input$family_id
      } else {
        names(catalog)[[1L]]
      }

      payloads <- lapply(seq_len(count), function(idx) {
        selection <- gripui.family.compare.slot.selection(
          idx = idx,
          input = input,
          catalog = catalog,
          current_family_id = current_family_id,
          lock_family = lock_family,
          include_current = include_current
        )

        if (identical(selection$source, "current")) {
          payload <- build_current_compare_payload()
          payload$compare_slot <- idx
          return(payload)
        }

        desc <- catalog[[selection$family_id]]
        values <- .gripui.family.merge.values(desc, preset_id = selection$preset_id)
        payload <- gripui.family.build.payload(desc, values)
        payload$compare_slot <- idx
        payload$compare_label <- sprintf("Variant %d", idx)
        payload$compare_preset <- selection$preset_id
        payload
      })

      compare_payloads_state(payloads)
      compare_error_state(NULL)
      invisible(payloads)
    }

    shiny::observe({
      payloads <- compare_payloads_state()
      choices <- gripui.family.compare.color.choices(payloads)
      current <- shiny::isolate(input$compare_color_by)
      selected <- if (!is.null(current) && current %in% unname(choices)) {
        current
      } else {
        gripui.family.default.color(choices)
      }
      shiny::updateSelectInput(session, "compare_color_by", choices = choices, selected = selected)
    })

    shiny::observeEvent(
      list(
        input$family_id,
        input$family_preset,
        input$compare_panel_count,
        input$compare_include_current,
        input$compare_lock_family
      ),
      {
        result <- tryCatch(build_compare_payloads(), error = function(e) e)
        if (inherits(result, "error")) {
          compare_error_state(conditionMessage(result))
        }
      },
      ignoreInit = FALSE
    )

    shiny::observeEvent(input$render_family_compare, {
      result <- tryCatch(build_compare_payloads(), error = function(e) e)
      if (inherits(result, "error")) {
        compare_error_state(conditionMessage(result))
        shiny::showNotification(conditionMessage(result), type = "error")
      }
    })

    output$render_status <- shiny::renderUI({
      payload <- payload_state()
      err <- error_state()
      saved_path <- save_state()
      if (!is.null(err)) {
        return(shiny::tags$p(style = "color:#8a1c1c;", err))
      }
      if (is.null(payload)) {
        return(shiny::tags$p(class = "gripui-muted", "Render a family to inspect its geometry."))
      }
      shiny::tagList(
        shiny::tags$p(
          class = "gripui-selection-status",
          sprintf("%s rendered with %d vertices and %d edges.", payload$family_label, payload$n, gripui.family.count.edges(payload$edges))
        ),
        if (!is.null(saved_path) && nzchar(saved_path)) {
          shiny::tags$p(
            class = "gripui-muted",
            style = "word-break:break-word;",
            sprintf("Last saved bundle: %s", saved_path)
          )
        }
      )
    })

    output$compare_render_status <- shiny::renderUI({
      payloads <- compare_payloads_state()
      err <- compare_error_state()
      if (!is.null(err)) {
        return(shiny::tags$p(style = "color:#8a1c1c;", err))
      }
      if (is.null(payloads) || length(payloads) == 0L) {
        return(shiny::tags$p(class = "gripui-muted", "Render a comparison to inspect multiple families side by side."))
      }
      shiny::tags$p(
        class = "gripui-selection-status",
        sprintf("Comparison rendered for %d panel%s.", length(payloads), if (length(payloads) == 1L) "" else "s")
      )
    })

    output$family_viewer_3d <- shiny::renderUI({
      payload <- payload_state()
      if (is.null(payload)) {
        return(shiny::tags$p(class = "gripui-muted", "No geometry rendered yet."))
      }
      shiny::tagList(
        gripui.render.rglwidget(
          payload$coords_display,
          graph = payload$graph,
          color_by = input$viewer_color_by,
          show_edges = isTRUE(input$show_edges)
        ),
        if (!is.null(payload$note) && nzchar(payload$note)) {
          shiny::tags$p(style = "margin-top:0.9rem;color:#6b7280;", payload$note)
        }
      )
    })

    output$family_viewer_2d <- shiny::renderPlot({
      payload <- payload_state()
      if (is.null(payload)) {
        graphics::plot.new()
        graphics::title("No geometry rendered yet")
        return(invisible(NULL))
      }
      gripui.render.layout.plot2d(
        payload$coords_plot,
        graph = payload$graph,
        color_by = input$viewer_color_by,
        show_edges = isTRUE(input$show_edges)
      )
    })

    for (idx in seq_len(4L)) {
      local({
        i <- idx
        output[[gripui.family.compare.plot.output_id(i)]] <- shiny::renderPlot({
          payloads <- compare_payloads_state()
          if (is.null(payloads) || length(payloads) < i || is.null(payloads[[i]])) {
            graphics::plot.new()
            graphics::title("No geometry rendered")
            return(invisible(NULL))
          }

          payload <- payloads[[i]]
          gripui.render.layout.plot2d(
            payload$coords_plot,
            graph = payload$graph,
            color_by = input$compare_color_by,
            show_edges = isTRUE(input$compare_show_edges)
          )
        })
      })
    }

    output$compare_grid <- shiny::renderUI({
      payloads <- compare_payloads_state()
      if (is.null(payloads) || length(payloads) == 0L) {
        return(shiny::tags$p(class = "gripui-muted", "No comparison rendered yet."))
      }

      cards <- lapply(seq_along(payloads), function(idx) {
        payload <- payloads[[idx]]
        preset_label <- if (identical(payload$compare_preset, "current_controls")) {
          "Current Explore controls"
        } else {
          tools::toTitleCase(gsub("_", " ", payload$compare_preset, fixed = TRUE))
        }

        bslib::card(
          full_screen = TRUE,
          class = "gripui-card",
          bslib::card_header(sprintf("%s: %s", payload$compare_label, payload$family_label)),
          shiny::tags$p(class = "gripui-project-summary", sprintf("Preset: %s", preset_label)),
          gripui.render.rglwidget(
            payload$coords_display,
            graph = payload$graph,
            color_by = input$compare_color_by,
            show_edges = isTRUE(input$compare_show_edges)
          ),
          shiny::plotOutput(gripui.family.compare.plot.output_id(idx), height = 220),
          if (!is.null(payload$note) && nzchar(payload$note)) {
            shiny::tags$p(style = "margin-top:0.9rem;color:#6b7280;", payload$note)
          }
        )
      })

      rows <- split(cards, ceiling(seq_along(cards) / 2L))
      shiny::tagList(lapply(rows, function(row_cards) {
        do.call(
          bslib::layout_columns,
          c(list(col_widths = if (length(row_cards) == 1L) 12 else rep(6, length(row_cards))), row_cards)
        )
      }))
    })

    output$family_graph_summary <- shiny::renderTable({
      payload <- payload_state()
      if (is.null(payload)) {
        return(NULL)
      }
      gripui.family.graph.summary(payload)
    }, striped = TRUE, bordered = FALSE, spacing = "s")

    output$family_weight_summary <- shiny::renderTable({
      payload <- payload_state()
      if (is.null(payload)) {
        return(NULL)
      }
      gripui.family.weight.summary(payload)
    }, striped = TRUE, bordered = FALSE, spacing = "s")

    output$family_geometry_summary <- shiny::renderTable({
      payload <- payload_state()
      if (is.null(payload)) {
        return(NULL)
      }
      gripui.family.geometry.summary(payload)
    }, striped = TRUE, bordered = FALSE, spacing = "s")

    output$family_code <- shiny::renderText({
      payload <- payload_state()
      if (is.null(payload)) {
        return("")
      }
      payload$code
    })

    output$compare_summary <- shiny::renderTable({
      payloads <- compare_payloads_state()
      if (is.null(payloads) || length(payloads) == 0L) {
        return(NULL)
      }
      gripui.family.compare.summary(payloads)
    }, striped = TRUE, bordered = FALSE, spacing = "s")
  }
}

#' Build the graph-family geometry explorer Shiny application
#'
#' @param catalog Family catalog, usually `gripui_graph_family_catalog()`.
#' @param title Application title.
#' @param subtitle Optional subtitle shown in the sidebar.
#'
#' @return A `shiny.appobj`.
#' @export
#'
#' @examplesIf local({ old <- getOption("rgl.useNULL"); options(rgl.useNULL = TRUE); on.exit(options(rgl.useNULL = old), add = TRUE); requireNamespace("shiny", quietly = TRUE) && requireNamespace("bslib", quietly = TRUE) && requireNamespace("rgl", quietly = TRUE) })
#' app <- gripui_family_app()
#' inherits(app, "shiny.appobj")
gripui_family_app <- function(catalog = gripui_graph_family_catalog(),
                              title = "Graph Family Geometry Explorer",
                              subtitle = "Interactive geometry browser for synthetic benchmark families.") {
  gripui.family.validate.catalog(catalog)
  old <- gripui.require.family.app.packages()
  on.exit(options(rgl.useNULL = old), add = TRUE)

  shiny::shinyApp(
    ui = gripui.family.ui(catalog = catalog, title = title, subtitle = subtitle),
    server = gripui.family.server(catalog = catalog)
  )
}

#' Run the graph-family geometry explorer Shiny application
#'
#' @param catalog Family catalog, usually `gripui_graph_family_catalog()`.
#' @param title Application title.
#' @param subtitle Optional subtitle shown in the sidebar.
#' @param host Host passed to `shiny::runApp()`.
#' @param port Port passed to `shiny::runApp()`.
#' @param launch.browser Whether to launch a browser.
#' @param auto.stop.after Optional delay, in seconds, after which the app stops
#'   itself. This is mainly useful for automated examples and tests.
#' @param ... Additional arguments passed to `shiny::runApp()`.
#'
#' @return Invisibly returns the result of `shiny::runApp()`.
#' @export
#'
#' @examplesIf local({ old <- getOption("rgl.useNULL"); options(rgl.useNULL = TRUE); on.exit(options(rgl.useNULL = old), add = TRUE); requireNamespace("shiny", quietly = TRUE) && requireNamespace("bslib", quietly = TRUE) && requireNamespace("rgl", quietly = TRUE) && requireNamespace("later", quietly = TRUE) })
#' run_gripui_family(launch.browser = FALSE, quiet = TRUE, auto.stop.after = 0.1)
run_gripui_family <- function(catalog = gripui_graph_family_catalog(),
                              title = "Graph Family Geometry Explorer",
                              subtitle = "Interactive geometry browser for synthetic benchmark families.",
                              host = "127.0.0.1",
                              port = getOption("shiny.port"),
                              launch.browser = interactive(),
                              auto.stop.after = NULL,
                              ...) {
  app <- gripui_family_app(catalog = catalog, title = title, subtitle = subtitle)

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
