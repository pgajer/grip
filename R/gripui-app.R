#' Build the `gripui` Shiny application
#'
#' @param project A `gripui_project`.
#'
#' @return A `shiny.appobj`.
#' @export
#'
#' @examplesIf requireNamespace("shiny", quietly = TRUE) && requireNamespace("bslib", quietly = TRUE) && requireNamespace("DT", quietly = TRUE) && requireNamespace("htmltools", quietly = TRUE) && requireNamespace("rgl", quietly = TRUE)
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
#'
#' @return Invisibly returns the result of `shiny::runApp()`.
#' @export
#'
#' @examples
#' \dontrun{
#' graph <- list(adj_list = list(2L, c(1L, 3L), 2L))
#' layouts <- data.frame(
#'   candidate = "toy.layout",
#'   stage = "layout",
#'   seed = 1L,
#'   status = "ok",
#'   stringsAsFactors = FALSE
#' )
#' project <- gripui_project(graph = graph, layouts = layouts, title = "Toy project")
#' run_gripui(project)
#' }
run_gripui <- function(project,
                       host = "127.0.0.1",
                       port = getOption("shiny.port"),
                       launch.browser = interactive()) {
  app <- gripui_app(project)
  shiny::runApp(
    app,
    host = host,
    port = port,
    launch.browser = launch.browser
  )
}
