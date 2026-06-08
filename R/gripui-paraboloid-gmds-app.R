gripui.paraboloid.gmds.metric.value <- function(x, digits = 4L) {
  if (!is.finite(x)) {
    return("NA")
  }
  formatC(x, format = "f", digits = digits)
}

gripui.paraboloid.gmds.grid.triangles <- function(side) {
  index <- matrix(seq_len(side * side), nrow = side, ncol = side, byrow = TRUE)
  triangles <- vector("list", 2L * (side - 1L) * (side - 1L))
  k <- 1L
  for (r in seq_len(side - 1L)) {
    for (c in seq_len(side - 1L)) {
      triangles[[k]] <- c(index[r, c], index[r + 1L, c], index[r, c + 1L])
      k <- k + 1L
      triangles[[k]] <- c(index[r + 1L, c], index[r + 1L, c + 1L], index[r, c + 1L])
      k <- k + 1L
    }
  }
  do.call(rbind, triangles)
}

gripui.paraboloid.gmds.triangle.areas <- function(coords, triangles) {
  coords <- as.matrix(coords)
  if (ncol(coords) == 2L) {
    coords <- cbind(coords, 0)
  }
  v1 <- coords[triangles[, 2L], , drop = FALSE] - coords[triangles[, 1L], , drop = FALSE]
  v2 <- coords[triangles[, 3L], , drop = FALSE] - coords[triangles[, 1L], , drop = FALSE]
  cross <- cbind(
    v1[, 2L] * v2[, 3L] - v1[, 3L] * v2[, 2L],
    v1[, 3L] * v2[, 1L] - v1[, 1L] * v2[, 3L],
    v1[, 1L] * v2[, 2L] - v1[, 2L] * v2[, 1L]
  )
  0.5 * sqrt(rowSums(cross^2))
}

gripui.paraboloid.gmds.area.floor.ratio <- function(coords, triangles) {
  areas <- gripui.paraboloid.gmds.triangle.areas(coords, triangles)
  med <- stats::median(areas)
  if (!is.finite(med) || med <= 0) {
    return(NA_real_)
  }
  as.double(stats::quantile(areas, probs = 0.05, names = FALSE)) / med
}

gripui.paraboloid.gmds.mesh.roughness <- function(coords, adj_list, edges) {
  centered <- sweep(coords, 2L, colMeans(coords), FUN = "-", check.margin = FALSE)
  median.edge <- stats::median(sqrt(rowSums(
    (centered[edges[, 1L], , drop = FALSE] - centered[edges[, 2L], , drop = FALSE])^2
  )))
  if (!is.finite(median.edge) || median.edge <= 0) {
    return(NA_real_)
  }
  residuals <- vapply(seq_len(nrow(centered)), function(i) {
    nbrs <- adj_list[[i]]
    if (length(nbrs) == 0L) {
      return(0)
    }
    delta <- centered[i, ] - colMeans(centered[nbrs, , drop = FALSE])
    sum(delta^2)
  }, numeric(1L))
  sqrt(mean(residuals)) / median.edge
}

gripui.paraboloid.gmds.make.case <- function(side = 12L, amplitude = 0.35) {
  side <- grip.validate.count(side, "side")
  if (side < 4L) {
    stop("side must be at least 4", call. = FALSE)
  }
  grip.validate.scalar(amplitude, "amplitude", lower = 0, open.lower = TRUE)
  bundle <- mesh.surface.graph(
    side,
    side,
    surface = "paraboloid",
    amplitude = amplitude,
    connectivity = "orthogonal",
    normalize = "median"
  )
  graph <- grip.build.adj.from.edges(bundle$edges, n = bundle$n)
  prepared <- prepare.graph.geodesic.mds(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average"
  )
  cmd <- grip.classical.mds.embedding(prepared, dim = 3L, eig = TRUE)
  list(
    side = side,
    amplitude = amplitude,
    label = sprintf("Orthogonal paraboloid mesh %dx%d", side, side),
    n = bundle$n,
    m = nrow(bundle$edges),
    truth = bundle$coords_surface,
    edges = bundle$edges,
    adj_list = graph$adj_list,
    graph = list(adj_list = graph$adj_list),
    triangles = gripui.paraboloid.gmds.grid.triangles(side),
    prepared = prepared,
    cmd = cmd
  )
}

gripui.paraboloid.gmds.metric.row <- function(case,
                                              coords,
                                              method,
                                              elapsed_sec = NA_real_) {
  score <- grip.score.geodesic.mds(coords = coords, prepared = case$prepared)
  aligned <- grip.align.to.target.nd(coords, case$truth, allow.reflection = TRUE)
  list(
    display_coords = aligned$aligned,
    metrics = data.frame(
      Method = method,
      sigma = score$gmds.stress[[1L]],
      rho = aligned$rmse,
      eta = gripui.paraboloid.gmds.mesh.roughness(coords, case$adj_list, case$edges),
      alpha_0.05 = gripui.paraboloid.gmds.area.floor.ratio(coords, case$triangles),
      runtime_sec = as.double(elapsed_sec),
      stringsAsFactors = FALSE
    )
  )
}

gripui.paraboloid.gmds.constant.schedule <- function(lambda) {
  if (!is.numeric(lambda) || length(lambda) != 1L || is.na(lambda) || lambda < 0) {
    stop("All lambda values must be single non-negative numbers.", call. = FALSE)
  }
  as.double(lambda)
}

gripui.paraboloid.gmds.compute.payload <- function(side = 12L,
                                                   amplitude = 0.35,
                                                   lambda_anchor = 0.10,
                                                   lambda_edge = 0.25,
                                                   lambda_repulsion = 0.20,
                                                   max_iter = 15L,
                                                   n_threads = 0L) {
  case <- gripui.paraboloid.gmds.make.case(side = side, amplitude = amplitude)
  max_iter <- grip.validate.count(max_iter, "max_iter")
  n_threads <- grip.validate.count(n_threads, "n_threads")
  lambda_anchor <- gripui.paraboloid.gmds.constant.schedule(lambda_anchor)
  lambda_edge <- gripui.paraboloid.gmds.constant.schedule(lambda_edge)
  lambda_repulsion <- gripui.paraboloid.gmds.constant.schedule(lambda_repulsion)

  anchor_mode <- if (lambda_anchor > 0) "cmdscale" else "none"

  reference <- gripui.paraboloid.gmds.metric.row(
    case = case,
    coords = case$truth,
    method = "Reference paraboloid",
    elapsed_sec = 0
  )
  reference$display_coords <- case$truth

  rep_started <- proc.time()[["elapsed"]]
  rep_fit <- grip.optimize.geodesic.mds(
    coords = case$cmd$coords,
    prepared = case$prepared,
    anchor_mode = anchor_mode,
    anchor_weight = lambda_anchor,
    anchor_weight_end = lambda_anchor,
    continuation = "constant",
    repulsion_weight = lambda_repulsion,
    repulsion_weight_end = lambda_repulsion,
    repulsion_continuation = "constant",
    engine = "cpp",
    max_iter = max_iter,
    n_threads = n_threads,
    return_trace = FALSE
  )
  rep_elapsed <- proc.time()[["elapsed"]] - rep_started
  repulsion <- gripui.paraboloid.gmds.metric.row(
    case = case,
    coords = rep_fit$coords,
    method = "Anchor + Repulsion",
    elapsed_sec = rep_elapsed
  )

  edge_started <- proc.time()[["elapsed"]]
  edge_fit <- grip.optimize.geodesic.mds(
    coords = case$cmd$coords,
    prepared = case$prepared,
    anchor_mode = anchor_mode,
    anchor_weight = lambda_anchor,
    anchor_weight_end = lambda_anchor,
    continuation = "constant",
    edge_spring_weight = lambda_edge,
    edge_spring_weight_end = lambda_edge,
    edge_spring_continuation = "constant",
    engine = "cpp",
    max_iter = max_iter,
    n_threads = n_threads,
    return_trace = FALSE
  )
  edge_elapsed <- proc.time()[["elapsed"]] - edge_started
  edge <- gripui.paraboloid.gmds.metric.row(
    case = case,
    coords = edge_fit$coords,
    method = "Anchor + Edge Spring",
    elapsed_sec = edge_elapsed
  )

  metric_table <- rbind(
    reference$metrics,
    repulsion$metrics,
    edge$metrics
  )
  metric_table$sigma <- round(metric_table$sigma, 4L)
  metric_table$rho <- round(metric_table$rho, 4L)
  metric_table$eta <- round(metric_table$eta, 4L)
  metric_table$alpha_0.05 <- round(metric_table$alpha_0.05, 4L)
  metric_table$runtime_sec <- round(metric_table$runtime_sec, 3L)

  list(
    case = case,
    reference = reference,
    anchor_repulsion = repulsion,
    anchor_edge_spring = edge,
    metric_table = metric_table,
    settings = list(
      lambda_anchor = lambda_anchor,
      lambda_edge = lambda_edge,
      lambda_repulsion = lambda_repulsion,
      max_iter = max_iter,
      n_threads = n_threads
    )
  )
}

gripui.paraboloid.gmds.panel.note <- function(payload, which = c("reference", "anchor_repulsion", "anchor_edge_spring")) {
  which <- match.arg(which)
  settings <- payload$settings
  switch(
    which,
    reference = sprintf(
      "Reference surface for %s (%d vertices, %d edges).",
      payload$case$label,
      payload$case$n,
      payload$case$m
    ),
    anchor_repulsion = sprintf(
      "Computed from the classical-MDS start with constant lambda_A = %s and lambda_R = %s for %d iterations.",
      gripui.paraboloid.gmds.metric.value(settings$lambda_anchor, 3L),
      gripui.paraboloid.gmds.metric.value(settings$lambda_repulsion, 3L),
      settings$max_iter
    ),
    anchor_edge_spring = sprintf(
      "Computed from the classical-MDS start with constant lambda_A = %s and lambda_E = %s for %d iterations.",
      gripui.paraboloid.gmds.metric.value(settings$lambda_anchor, 3L),
      gripui.paraboloid.gmds.metric.value(settings$lambda_edge, 3L),
      settings$max_iter
    )
  )
}

gripui.paraboloid.gmds.panel.card <- function(title, coords, graph, note) {
  bslib::card(
    full_screen = TRUE,
    bslib::card_header(title),
    shiny::tags$div(
      style = "padding:0.6rem 0.6rem 0.2rem;",
      gripui.render.rglwidget(coords = coords, graph = graph, show_edges = TRUE)
    ),
    shiny::tags$p(
      style = "padding:0 1rem 1rem;color:#5f5445;line-height:1.45;margin-bottom:0;",
      note
    )
  )
}

gripui.paraboloid.gmds.ui <- function(title, subtitle) {
  bslib::page_sidebar(
    theme = gripui.family.theme(),
    title = shiny::div(
      style = "display:flex;flex-direction:column;gap:0.25rem;",
      shiny::div(class = "gripui-brand", title),
      shiny::div(style = "font-size:0.98rem;color:#5f5445;font-weight:400;", subtitle)
    ),
    sidebar = bslib::sidebar(
      width = 340,
      shiny::numericInput("mesh_side", "Mesh side length", value = 12L, min = 4L, step = 1L),
      shiny::numericInput("mesh_amplitude", "Paraboloid amplitude", value = 0.35, min = 0.01, step = 0.05),
      shiny::numericInput("lambda_anchor", "lambda_A", value = 0.10, min = 0, step = 0.05),
      shiny::numericInput("lambda_repulsion", "lambda_R", value = 0.20, min = 0, step = 0.05),
      shiny::numericInput("lambda_edge", "lambda_E", value = 0.25, min = 0, step = 0.05),
      shiny::numericInput("max_iter", "Max iterations", value = 15L, min = 0L, step = 1L),
      shiny::numericInput("n_threads", "Threads (0 = auto)", value = 0L, min = 0L, step = 1L),
      shiny::actionButton("render_layouts", "Compute layouts", class = "btn-primary"),
      shiny::uiOutput("render_status"),
      shiny::tags$hr(),
      shiny::tags$p(
        style = "color:#5f5445;line-height:1.5;",
        "The app keeps the user-entered lambda values constant across all iterations. The two computed panels always use the same reference paraboloid and the same classical-MDS initialization."
      )
    ),
    bslib::card(
      full_screen = TRUE,
      bslib::card_header("Metric Summary"),
      shiny::tableOutput("metric_table")
    ),
    shiny::uiOutput("layout_panels")
  )
}

gripui.paraboloid.gmds.server <- function() {
  function(input, output, session) {
    payload_state <- shiny::reactiveVal(NULL)
    error_state <- shiny::reactiveVal(NULL)

    compute_payload <- function() {
      shiny::withProgress(message = "Computing paraboloid layouts", value = 0, {
        shiny::incProgress(0.15, detail = "Building reference mesh")
        side <- as.integer(round(input$mesh_side))
        amplitude <- as.double(input$mesh_amplitude)
        lambda_anchor <- as.double(input$lambda_anchor)
        lambda_edge <- as.double(input$lambda_edge)
        lambda_repulsion <- as.double(input$lambda_repulsion)
        max_iter <- as.integer(round(input$max_iter))
        n_threads <- as.integer(round(input$n_threads))
        shiny::incProgress(0.30, detail = "Optimizing anchor + repulsion layout")
        payload <- gripui.paraboloid.gmds.compute.payload(
          side = side,
          amplitude = amplitude,
          lambda_anchor = lambda_anchor,
          lambda_edge = lambda_edge,
          lambda_repulsion = lambda_repulsion,
          max_iter = max_iter,
          n_threads = n_threads
        )
        shiny::incProgress(1, detail = "Done")
        payload
      })
    }

    shiny::observeEvent(input$render_layouts, {
      result <- tryCatch(
        compute_payload(),
        error = function(e) e
      )
      if (inherits(result, "error")) {
        error_state(conditionMessage(result))
        shiny::showNotification(conditionMessage(result), type = "error")
        return(invisible(NULL))
      }
      payload_state(result)
      error_state(NULL)
      invisible(NULL)
    }, ignoreNULL = FALSE)

    output$render_status <- shiny::renderUI({
      err <- error_state()
      payload <- payload_state()
      if (!is.null(err)) {
        return(shiny::tags$p(style = "color:#8a1c1c;", err))
      }
      if (is.null(payload)) {
        return(shiny::tags$p(style = "color:#5f5445;", "Computing the initial view..."))
      }
      shiny::tags$p(
        style = "color:#5f5445;line-height:1.45;",
        sprintf(
          "Rendered %s with %d vertices and %d edges. User weights: lambda_A = %s, lambda_R = %s, lambda_E = %s.",
          payload$case$label,
          payload$case$n,
          payload$case$m,
          gripui.paraboloid.gmds.metric.value(payload$settings$lambda_anchor, 3L),
          gripui.paraboloid.gmds.metric.value(payload$settings$lambda_repulsion, 3L),
          gripui.paraboloid.gmds.metric.value(payload$settings$lambda_edge, 3L)
        )
      )
    })

    output$metric_table <- shiny::renderTable({
      payload <- payload_state()
      shiny::req(payload)
      payload$metric_table
    }, striped = TRUE, bordered = FALSE, digits = 4L)

    output$layout_panels <- shiny::renderUI({
      payload <- payload_state()
      shiny::req(payload)
      cards <- list(
        gripui.paraboloid.gmds.panel.card(
          title = "Reference Paraboloid",
          coords = payload$reference$display_coords,
          graph = payload$case$graph,
          note = gripui.paraboloid.gmds.panel.note(payload, "reference")
        ),
        gripui.paraboloid.gmds.panel.card(
          title = sprintf(
            "Anchor + Repulsion (lambda_A = %s, lambda_R = %s)",
            gripui.paraboloid.gmds.metric.value(payload$settings$lambda_anchor, 3L),
            gripui.paraboloid.gmds.metric.value(payload$settings$lambda_repulsion, 3L)
          ),
          coords = payload$anchor_repulsion$display_coords,
          graph = payload$case$graph,
          note = gripui.paraboloid.gmds.panel.note(payload, "anchor_repulsion")
        ),
        gripui.paraboloid.gmds.panel.card(
          title = sprintf(
            "Anchor + Edge Spring (lambda_A = %s, lambda_E = %s)",
            gripui.paraboloid.gmds.metric.value(payload$settings$lambda_anchor, 3L),
            gripui.paraboloid.gmds.metric.value(payload$settings$lambda_edge, 3L)
          ),
          coords = payload$anchor_edge_spring$display_coords,
          graph = payload$case$graph,
          note = gripui.paraboloid.gmds.panel.note(payload, "anchor_edge_spring")
        )
      )
      do.call(bslib::layout_column_wrap, c(list(width = 1 / 3), cards))
    })
  }
}

#' Build the paraboloid GMDS explorer Shiny application
#'
#' This app lets you choose a paraboloid mesh size, set constant non-negative
#' values for \eqn{\lambda_A}, \eqn{\lambda_R}, and \eqn{\lambda_E}, and compare
#' the reference paraboloid with the anchor-plus-repulsion and
#' anchor-plus-edge-spring regularized GMDS embeddings computed on the fly.
#'
#' @param title Application title.
#' @param subtitle Optional subtitle shown in the sidebar header.
#'
#' @return A `shiny.appobj`.
#' @noRd
#'
#' @examplesIf local({ old <- getOption("rgl.useNULL"); options(rgl.useNULL = TRUE); on.exit(options(rgl.useNULL = old), add = TRUE); requireNamespace("shiny", quietly = TRUE) && requireNamespace("bslib", quietly = TRUE) && requireNamespace("rgl", quietly = TRUE) })
#' app <- gripui_paraboloid_gmds_app()
#' inherits(app, "shiny.appobj")
gripui_paraboloid_gmds_app <- function(title = "Paraboloid GMDS Explorer",
                                       subtitle = "Reference geometry plus two on-the-fly regularized GMDS variants.") {
  old <- gripui.require.family.app.packages()
  on.exit(options(rgl.useNULL = old), add = TRUE)

  shiny::shinyApp(
    ui = gripui.paraboloid.gmds.ui(title = title, subtitle = subtitle),
    server = gripui.paraboloid.gmds.server()
  )
}

#' Run the paraboloid GMDS explorer Shiny application
#'
#' @param title Application title.
#' @param subtitle Optional subtitle shown in the sidebar header.
#' @param host Host passed to `shiny::runApp()`.
#' @param port Port passed to `shiny::runApp()`.
#' @param launch.browser Whether to launch a browser.
#' @param auto.stop.after Optional delay, in seconds, after which the app stops
#'   itself. This is mainly useful for automated examples and tests.
#' @param ... Additional arguments passed to `shiny::runApp()`.
#'
#' @return Invisibly returns the result of `shiny::runApp()`.
#' @noRd
#'
#' @examplesIf local({ old <- getOption("rgl.useNULL"); options(rgl.useNULL = TRUE); on.exit(options(rgl.useNULL = old), add = TRUE); requireNamespace("shiny", quietly = TRUE) && requireNamespace("bslib", quietly = TRUE) && requireNamespace("rgl", quietly = TRUE) && requireNamespace("later", quietly = TRUE) })
#' run_gripui_paraboloid_gmds(launch.browser = FALSE, quiet = TRUE, auto.stop.after = 0.1)
run_gripui_paraboloid_gmds <- function(title = "Paraboloid GMDS Explorer",
                                       subtitle = "Reference geometry plus two on-the-fly regularized GMDS variants.",
                                       host = "127.0.0.1",
                                       port = getOption("shiny.port"),
                                       launch.browser = interactive(),
                                       auto.stop.after = NULL,
                                       ...) {
  app <- gripui_paraboloid_gmds_app(title = title, subtitle = subtitle)

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
