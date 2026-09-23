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
}
writeLines(c(paste('ivue', packageVersion('ivue')), paste('rgl', packageVersion('rgl')),
  paste('Bundle MD5', unname(tools::md5sum(args[1])))), file.path(out, 'render-provenance.txt'))
