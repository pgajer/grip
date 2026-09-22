#!/usr/bin/env Rscript
# Render the saved, predeclared example fits through ivue. No fitting occurs.
# Run from the package root: Rscript tools/pkg/render-metric-mds-comparison.R <bundle.rds> <output-dir>
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
source('inst/scripts/metric-mds-comparison.R')
bundle <- readRDS(args[[1]])
dir.create(args[[2]], recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(args[[2]])
for (case in bundle$cases) {
  if (case$dimension == 2L) next
  fits <- comparison_representative(bundle, case$id)
  w <- comparison_view(case, fits)
  w$elementId <- paste0('comparison-', case$id)
  htmlwidgets::saveWidget(w, file.path(out, paste0(case$id, '.html')), selfcontained = TRUE)
  panels <- lapply(c('reference','sgd','smacof'), function(mode) {
    label <- if (mode == 'reference') 'Generating reference' else toupper(mode)
    detail <- if (mode == 'reference') case$target else sprintf(
      'Error %.3g; fit %.3f s', fits[[mode]]$row$error, fits[[mode]]$row$seconds)
    htmltools::tags$section(htmltools::tags$h3(label), htmltools::tags$p(detail),
      if (mode != 'reference' && is.null(fits[[mode]]$coords)) htmltools::tags$p('Fit unavailable') else
        comparison_view(case, fits, show = mode, width = 400L, height = 350L,
                        legend = FALSE, controls = FALSE))
  })
  page <- htmltools::tagList(htmltools::tags$style(htmltools::HTML(
    'body {margin:8px; font:15px system-ui; width:1200px;} .panels {display:flex;} h3,p {text-align:center; margin:4px;} section {width:400px;} .panels .ivue-description {display:none;}')),
    htmltools::tags$div(class = 'panels', panels))
  htmltools::save_html(page, file.path(out, paste0(case$id, '-panels.html')), libdir = 'panel-libs')
}
writeLines(c('Views use fixed seed 1 and allowance 100, with rigid alignment and no display rescaling.',
  paste('ivue:', as.character(packageVersion('ivue'))),
  paste('rgl:', as.character(packageVersion('rgl'))),
  paste('Bundle MD5:', unname(tools::md5sum(args[[1]])))), file.path(out, 'render-provenance.txt'))
