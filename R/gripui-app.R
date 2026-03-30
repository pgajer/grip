#' Build the `gripui` Shiny application
#'
#' @param project A `gripui_project`.
#'
#' @return A `shiny.appobj`.
#' @export
#'
#' @examplesIf local({ old <- getOption("rgl.useNULL"); options(rgl.useNULL = TRUE); on.exit(options(rgl.useNULL = old), add = TRUE); requireNamespace("shiny", quietly = TRUE) && requireNamespace("bslib", quietly = TRUE) && requireNamespace("DT", quietly = TRUE) && requireNamespace("htmltools", quietly = TRUE) && requireNamespace("rgl", quietly = TRUE) })
#' graph <- list(adj_list = list(2L, c(1L, 3L), 2L))
#' layouts <- data.frame(
#'   candidate = "toy.layout",
#'   stage = "layout",
#'   seed = 1L,
#'   status = "ok",
#'   stringsAsFactors = FALSE
#' )
#' project <- gripui_project(graph = graph, layouts = layouts, title = "Toy project")
#' app <- gripui_app(project)
#' inherits(app, "shiny.appobj")
gripui_app <- function(project) {
  gripui_validate_project(project)
  old <- gripui.require.app.packages()
  on.exit(options(rgl.useNULL = old), add = TRUE)

  shiny::shinyApp(
    ui = gripui_ui(),
    server = gripui_server(project)
  )
}

#' Run the `gripui` Shiny application
#'
#' @param project A `gripui_project`.
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
#' @examplesIf local({ old <- getOption("rgl.useNULL"); options(rgl.useNULL = TRUE); on.exit(options(rgl.useNULL = old), add = TRUE); requireNamespace("shiny", quietly = TRUE) && requireNamespace("bslib", quietly = TRUE) && requireNamespace("DT", quietly = TRUE) && requireNamespace("htmltools", quietly = TRUE) && requireNamespace("rgl", quietly = TRUE) && requireNamespace("later", quietly = TRUE) })
#' graph <- list(adj_list = list(2L, c(1L, 3L), 2L))
#' layouts <- data.frame(
#'   candidate = "toy.layout",
#'   stage = "layout",
#'   seed = 1L,
#'   status = "ok",
#'   stringsAsFactors = FALSE
#' )
#' project <- gripui_project(graph = graph, layouts = layouts, title = "Toy project")
#' run_gripui(
#'   project,
#'   launch.browser = FALSE,
#'   quiet = TRUE,
#'   auto.stop.after = 0.1
#' )
run_gripui <- function(project,
                       host = "127.0.0.1",
                       port = getOption("shiny.port"),
                       launch.browser = interactive(),
                       auto.stop.after = NULL,
                       ...) {
  app <- gripui_app(project)

  if (!is.null(auto.stop.after)) {
    if (!requireNamespace("later", quietly = TRUE)) {
      stop("Package 'later' is required when auto.stop.after is used.")
    }
    if (!is.numeric(auto.stop.after) ||
        length(auto.stop.after) != 1L ||
        is.na(auto.stop.after) ||
        auto.stop.after < 0) {
      stop("auto.stop.after must be a single non-negative number of seconds.")
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
