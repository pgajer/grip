gripui_ui <- function() {
  css.path <- system.file("app/www/gripui.css", package = "grip")
  theme <- bslib::bs_theme(
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

  bslib::page_sidebar(
    title = shiny::div(
      class = "gripui-titlebar",
      shiny::div(class = "gripui-brand", "gripui"),
      shiny::uiOutput("project_badge")
    ),
    theme = theme,
    sidebar = bslib::sidebar(
      width = 360,
      shiny::uiOutput("project_summary"),
      shiny::selectizeInput("stage_filter", "Stages", choices = NULL, multiple = TRUE),
      shiny::textInput("candidate_search", "Candidate search", placeholder = "candidate or layout id"),
      shiny::checkboxInput("viewable_only", "Only show viewable layouts", value = FALSE),
      shiny::selectInput("x_axis", "Landscape x", choices = NULL),
      shiny::selectInput("y_axis", "Landscape y", choices = NULL),
      shiny::selectInput("landscape_color_by", "Landscape color", choices = NULL),
      shiny::selectInput("viewer_color_by", "Viewer color", choices = "plain", selected = "plain"),
      shiny::checkboxInput("show_edges", "Show edges when graph is available", value = TRUE),
      shiny::uiOutput("selection_status")
    ),
    if (nzchar(css.path)) shiny::tags$head(shiny::includeCSS(css.path)),
    bslib::layout_columns(
      col_widths = c(5, 7),
      bslib::card(
        full_screen = TRUE,
        class = "gripui-card",
        bslib::card_header("Layout Landscape"),
        shiny::plotOutput("landscape_plot", height = 380, click = "landscape_click")
      ),
      bslib::card(
        full_screen = TRUE,
        class = "gripui-card",
        bslib::card_header(shiny::uiOutput("selected_title")),
        shiny::uiOutput("viewer_panel")
      )
    ),
    bslib::layout_columns(
      col_widths = c(8, 4),
      bslib::card(
        full_screen = TRUE,
        class = "gripui-card",
        bslib::card_header("Catalog"),
        DT::DTOutput("catalog_table")
      ),
      bslib::card(
        full_screen = TRUE,
        class = "gripui-card",
        bslib::card_header("Selected Layout"),
        shiny::tableOutput("selected_info"),
        shiny::uiOutput("selected_links")
      )
    )
  )
}
