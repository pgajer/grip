#!/usr/bin/env Rscript
# Rscript tools/pkg/render-graph-examples.R <layouts.rds> <output-dir>
args <- commandArgs(TRUE); stopifnot(length(args) == 2L)
source('inst/scripts/graph-examples.R')
bundle <- readRDS(args[1])
dir.create(args[2], recursive = TRUE, showWarnings = FALSE); out <- normalizePath(args[2])
for (id in names(bundle$cases)) {
  graph <- bundle$cases[[id]]; fits <- bundle$fits[[id]]
  if (!is.null(fits$error)) next
  cat('Rendering', id, '\n'); flush.console()
  htmlwidgets::saveWidget(graph_example_view(graph, fits), file.path(out, paste0(id, '.html')), selfcontained = TRUE)
  modes <- c(if (!is.null(graph$coords_surface)) 'reference', 'mds', 'refined')
  width <- 1200L / length(modes)
  panels <- lapply(modes, function(mode) {
    label <- c(reference = 'Generating coordinates', mds = 'Metric MDS', refined = 'MDS + edge-KK')[[mode]]
    detail <- if (mode == 'reference') graph$lengths else if (is.null(fits[[mode]]$scores)) 'Unavailable' else
      sprintf('Graph error %.3g; edge error %.3g', fits[[mode]]$scores[['graph_error']], fits[[mode]]$scores[['edge_error']])
    htmltools::tags$section(style = paste0('width:', width, 'px;'), htmltools::tags$h3(label), htmltools::tags$p(detail),
      if (mode != 'reference' && is.null(fits[[mode]]$coords)) htmltools::tags$p('Fit unavailable') else
        graph_example_view(graph, fits, show = mode, width = width, height = 350L, controls = FALSE))
  })
  htmltools::save_html(htmltools::tagList(htmltools::tags$style(htmltools::HTML(
    'body {margin:8px; font:15px system-ui; width:1200px;} .panels {display:flex;} h3,p {text-align:center; margin:4px;} .panels .ivue-description {display:none;}')),
    htmltools::tags$div(class = 'panels', panels)), file.path(out, paste0(id, '-panels.html')), libdir = 'panel-libs')
}
writeLines(c(paste('ivue', packageVersion('ivue')), paste('rgl', packageVersion('rgl')),
  paste('Bundle MD5', unname(tools::md5sum(args[1])))), file.path(out, 'render-provenance.txt'))
