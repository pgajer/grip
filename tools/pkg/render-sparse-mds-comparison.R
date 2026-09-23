#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
source('inst/scripts/metric-mds-comparison.R')
source('inst/scripts/metric-mds-surface-alignment.R')
source('inst/scripts/sparse-mds-comparison.R')
bundle <- readRDS(args[1]); dir.create(args[2],recursive=TRUE,showWarnings=FALSE)
for (id in names(bundle$cases)) {
  selected <- Filter(function(f) f$row$case==id && f$row$seed==1,bundle$fits)
  fits <- stats::setNames(selected,vapply(selected,function(f) f$row$method,''))
  htmlwidgets::saveWidget(sparse_comparison_view(bundle$cases[[id]],fits),
    file.path(normalizePath(args[2]),paste0(id,'.html')),selfcontained=TRUE)
}
