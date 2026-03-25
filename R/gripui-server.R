gripui_server <- function(project) {
  function(input, output, session) {
    catalog <- project$layouts
    selected.layout.id <- shiny::reactiveVal(
      if (nrow(catalog) > 0L) catalog$layout_id[[1L]] else NULL
    )

    output$project_badge <- shiny::renderUI({
      shiny::span(
        class = "gripui-badge",
        sprintf("%d layouts", nrow(catalog))
      )
    })

    output$project_summary <- shiny::renderUI({
      title <- if (!is.null(project$meta$title) && nzchar(project$meta$title)) {
        project$meta$title
      } else {
        "Layout explorer"
      }
      subtitle <- if (!is.null(project$meta$subtitle) && nzchar(project$meta$subtitle)) {
        project$meta$subtitle
      } else {
        NULL
      }
      graphs <- if (is.null(project$graph)) {
        "no graph attached"
      } else if (gripui.is.graph.object(project$graph)) {
        "1 graph attached"
      } else {
        sprintf("%d graphs attached", length(project$graph))
      }
      shiny::tagList(
        shiny::tags$div(class = "gripui-project-title", title),
        if (!is.null(subtitle)) shiny::tags$div(class = "gripui-project-subtitle", subtitle),
        shiny::tags$p(
          class = "gripui-project-summary",
          sprintf(
            "%d total rows, %d viewable rows, %s",
            nrow(catalog),
            sum(catalog$viewable, na.rm = TRUE),
            graphs
          )
        )
      )
    })

    shiny::observe({
      stage.choices <- sort(unique(as.character(catalog$stage)))
      shiny::updateSelectizeInput(
        session = session,
        inputId = "stage_filter",
        choices = stage.choices,
        selected = stage.choices,
        server = TRUE
      )

      numeric.cols <- gripui.catalog.numeric.columns(catalog)
      x.default <- gripui.default.axis(
        numeric.cols,
        preferred = c("score_composite_extended", "score_composite", "sampled_stress", "elapsed_sec")
      )
      y.default <- gripui.default.axis(
        setdiff(numeric.cols, x.default),
        preferred = c("cst_cluster_separation", "subcst_cluster_separation", "sampled_nonedge_sep_ratio", "edge_length_cv"),
        fallback = x.default
      )
      shiny::updateSelectInput(session, "x_axis", choices = numeric.cols, selected = x.default)
      shiny::updateSelectInput(session, "y_axis", choices = numeric.cols, selected = y.default)

      color.choices <- unique(c("availability", "stage", "status", numeric.cols))
      shiny::updateSelectInput(
        session,
        "landscape_color_by",
        choices = color.choices,
        selected = if ("availability" %in% color.choices) "availability" else color.choices[[1L]]
      )
    })

    filtered.catalog <- shiny::reactive({
      rows <- catalog

      stages <- input$stage_filter
      if (!is.null(stages) && length(stages) > 0L) {
        rows <- rows[rows$stage %in% stages, , drop = FALSE]
      }
      if (isTRUE(input$viewable_only)) {
        rows <- rows[!is.na(rows$viewable) & rows$viewable, , drop = FALSE]
      }
      search <- trimws(as.character(if (is.null(input$candidate_search)) "" else input$candidate_search))
      if (nzchar(search)) {
        keep <- (!is.na(rows$candidate) & grepl(search, rows$candidate, ignore.case = TRUE)) |
          (!is.na(rows$layout_id) & grepl(search, rows$layout_id, ignore.case = TRUE))
        rows <- rows[keep, , drop = FALSE]
      }
      rows
    })

    shiny::observe({
      rows <- filtered.catalog()
      current <- selected.layout.id()
      if (nrow(rows) == 0L) {
        selected.layout.id(NULL)
        return()
      }
      if (is.null(current) || !(current %in% rows$layout_id)) {
        selected.layout.id(rows$layout_id[[1L]])
      }
    })

    selected.row <- shiny::reactive({
      id <- selected.layout.id()
      if (is.null(id) || !nzchar(id)) {
        return(NULL)
      }
      row <- filtered.catalog()[filtered.catalog()$layout_id == id, , drop = FALSE]
      if (nrow(row) == 1L) {
        return(row)
      }
      row <- catalog[catalog$layout_id == id, , drop = FALSE]
      if (nrow(row) == 1L) {
        return(row)
      }
      NULL
    })

    selected.graph <- shiny::reactive({
      row <- selected.row()
      if (is.null(row)) {
        return(NULL)
      }
      gripui.get.graph(project, if ("graph_id" %in% names(row)) row$graph_id[[1L]] else NULL)
    })

    selected.coords <- shiny::reactive({
      row <- selected.row()
      if (is.null(row) || !isTRUE(row$viewable[[1L]])) {
        return(NULL)
      }
      tryCatch(
        gripui_load_layout_coords(project, row),
        error = function(e) NULL
      )
    })

    shiny::observe({
      row <- selected.row()
      graph <- selected.graph()
      choices <- "plain"
      selected <- "plain"
      if (!is.null(graph) && !is.null(graph$vertex_data)) {
        meta.cols <- setdiff(names(graph$vertex_data), c("x", "y", "z"))
        if (length(meta.cols) > 0L) {
          choices <- c("plain", meta.cols)
        }
      }
      if (!is.null(row) && "color_view_default" %in% names(row)) {
        preferred <- as.character(row$color_view_default[[1L]])
        if (preferred %in% choices) {
          selected <- preferred
        }
      }
      current <- shiny::isolate(input$viewer_color_by)
      if (is.character(current) && length(current) == 1L && current %in% choices) {
        selected <- current
      }
      shiny::updateSelectInput(
        session,
        "viewer_color_by",
        choices = choices,
        selected = selected
      )
    })

    output$selection_status <- shiny::renderUI({
      row <- selected.row()
      if (is.null(row)) {
        return(shiny::tags$p(class = "gripui-muted", "No layout selected."))
      }
      shiny::tags$p(
        class = "gripui-selection-status",
        sprintf(
          "Selected: %s (%s, %s)",
          row$candidate[[1L]],
          row$stage[[1L]],
          row$availability[[1L]]
        )
      )
    })

    output$landscape_plot <- shiny::renderPlot({
      rows <- filtered.catalog()
      xvar <- input$x_axis
      yvar <- input$y_axis
      if (nrow(rows) == 0L ||
          !is.character(xvar) || length(xvar) != 1L || !nzchar(xvar) ||
          !is.character(yvar) || length(yvar) != 1L || !nzchar(yvar) ||
          !(xvar %in% names(rows)) || !(yvar %in% names(rows))) {
        graphics::plot.new()
        graphics::title("No layouts available for the current filters")
        return(invisible(NULL))
      }

      x <- rows[[xvar]]
      y <- rows[[yvar]]
      keep <- is.finite(x) & is.finite(y)
      if (!any(keep)) {
        graphics::plot.new()
        graphics::title("Selected axes do not have finite values")
        return(invisible(NULL))
      }

      color.by <- input$landscape_color_by
      color.values <- if (is.character(color.by) && length(color.by) == 1L &&
                          nzchar(color.by) && color.by %in% names(rows)) {
        rows[[color.by]]
      } else {
        rows$availability
      }
      cols <- gripui.landscape.colors(color.values)

      graphics::plot(
        x[keep],
        y[keep],
        pch = 19,
        cex = 1.1,
        col = cols[keep],
        xlab = xvar,
        ylab = yvar
      )

      selected.id <- selected.layout.id()
      if (!is.null(selected.id) && selected.id %in% rows$layout_id) {
        i <- match(selected.id, rows$layout_id)
        if (!is.na(i) && keep[[i]]) {
          graphics::points(
            x[[i]],
            y[[i]],
            pch = 21,
            cex = 1.8,
            lwd = 2,
            bg = NA,
            col = "#111827"
          )
        }
      }
      invisible(NULL)
    })

    shiny::observeEvent(input$landscape_click, {
      rows <- filtered.catalog()
      xvar <- input$x_axis
      yvar <- input$y_axis
      if (nrow(rows) == 0L || !(xvar %in% names(rows)) || !(yvar %in% names(rows))) {
        return()
      }
      hit <- shiny::nearPoints(rows, input$landscape_click, xvar = xvar, yvar = yvar, maxpoints = 1L)
      if (nrow(hit) == 1L) {
        selected.layout.id(hit$layout_id[[1L]])
      }
    })

    output$catalog_table <- DT::renderDT({
      rows <- filtered.catalog()
      show.cols <- gripui.display.columns(rows)
      DT::datatable(
        rows[, show.cols, drop = FALSE],
        selection = "single",
        rownames = FALSE,
        options = list(pageLength = 12, scrollX = TRUE)
      )
    })

    shiny::observeEvent(input$catalog_table_rows_selected, {
      idx <- input$catalog_table_rows_selected
      rows <- filtered.catalog()
      if (length(idx) == 1L && nrow(rows) >= idx) {
        selected.layout.id(rows$layout_id[[idx]])
      }
    }, ignoreNULL = TRUE)

    shiny::observe({
      rows <- filtered.catalog()
      id <- selected.layout.id()
      proxy <- DT::dataTableProxy("catalog_table")
      if (nrow(rows) == 0L || is.null(id) || !(id %in% rows$layout_id)) {
        DT::selectRows(proxy, NULL)
      } else {
        DT::selectRows(proxy, match(id, rows$layout_id))
      }
    })

    output$selected_title <- shiny::renderUI({
      row <- selected.row()
      if (is.null(row)) {
        return("Selected Layout")
      }
      shiny::tags$div(
        class = "gripui-selected-title",
        shiny::tags$span(class = "gripui-selected-candidate", row$candidate[[1L]]),
        shiny::tags$span(class = "gripui-selected-meta", paste(row$stage[[1L]], row$availability[[1L]], sep = " | "))
      )
    })

    output$viewer_panel <- shiny::renderUI({
      row <- selected.row()
      coords <- selected.coords()
      image.path <- if (is.null(row)) NULL else gripui.preview.image.path(row)
      if (is.null(row)) {
        return(shiny::tags$p(class = "gripui-muted", "No layout selected."))
      }
      if (is.matrix(coords) && ncol(coords) >= 3L) {
        return(rgl::rglwidgetOutput("layout_view_3d", width = "100%", height = "520px"))
      }
      if (is.matrix(coords) && ncol(coords) >= 2L) {
        return(shiny::plotOutput("layout_view_2d", height = "520px"))
      }
      if (!is.null(image.path)) {
        return(shiny::tagList(
          shiny::imageOutput("layout_artifact_image", height = "520px"),
          shiny::tags$p(
            class = "gripui-muted",
            "Showing the saved artifact preview because raw coordinates are not attached for this row."
          )
        ))
      }
      if (!is.na(row$html_path[[1L]]) && nzchar(row$html_path[[1L]])) {
        return(shiny::tagList(
          shiny::tags$p(
            class = "gripui-muted",
            "This layout does not have raw coordinates attached in the project catalog, but a saved HTML artifact is available below."
          ),
          gripui.artifact.links(row)
        ))
      }
      shiny::tags$p(
        class = "gripui-muted",
        "This row is summary-only. Filter to a viewable row to open an interactive layout."
      )
    })

    output$layout_view_2d <- shiny::renderPlot({
      coords <- selected.coords()
      shiny::validate(shiny::need(is.matrix(coords) && ncol(coords) >= 2L, "No 2D coordinates available"))
      graph <- if (isTRUE(input$show_edges)) selected.graph() else NULL
      gripui.render.layout.plot2d(coords, graph = graph, color_by = input$viewer_color_by)
    })

    output$layout_view_3d <- rgl::renderRglwidget({
      coords <- selected.coords()
      shiny::validate(shiny::need(is.matrix(coords) && ncol(coords) >= 3L, "No 3D coordinates available"))
      graph <- if (isTRUE(input$show_edges)) selected.graph() else NULL
      gripui.render.rglwidget(coords, graph = graph, color_by = input$viewer_color_by)
    })

    output$layout_artifact_image <- shiny::renderImage({
      row <- selected.row()
      shiny::validate(shiny::need(!is.null(row), "No layout selected"))
      path <- gripui.preview.image.path(row)
      shiny::validate(shiny::need(!is.null(path), "No preview image available"))
      list(
        src = path,
        contentType = gripui.image.content.type(path),
        alt = sprintf("Saved preview for %s", row$candidate[[1L]])
      )
    }, deleteFile = FALSE)

    output$selected_info <- shiny::renderTable({
      row <- selected.row()
      if (is.null(row)) {
        return(NULL)
      }
      gripui.selected.info(row)
    }, striped = TRUE, bordered = FALSE, spacing = "xs")

    output$selected_links <- shiny::renderUI({
      row <- selected.row()
      if (is.null(row)) {
        return(NULL)
      }
      gripui.artifact.links(row)
    })
  }
}
