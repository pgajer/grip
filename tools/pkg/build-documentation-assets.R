#!/usr/bin/env Rscript
# Rebuild presentation transforms only; benchmark fits are never rerun here.
pkgload::load_all('.',quiet=TRUE,export_all=FALSE,helpers=FALSE)
source('inst/scripts/metric-mds-comparison.R')
source('inst/scripts/metric-mds-surface-alignment.R')
source('inst/scripts/sparse-mds-comparison.R')
options(grip.documentation.rebuild=TRUE,
        grip.documentation.display.dir=normalizePath('inst/extdata/documentation-display'))
full <- readRDS('inst/extdata/metric-mds-comparison/benchmark.rds')
for(id in c('saddle_graph-64','saddle_graph-144')) {
  message('Registering ',id)
  invisible(comparison.surface.saved(full$cases[[id]],comparison_representative(full,id)))
}
weighted <- readRDS('inst/extdata/metric-mds-weighting/benchmark.rds')
id <- 'saddle_graph-64'
invisible(comparison.surface.saved(weighted$cases[[id]],comparison_weighting_fits(weighted,id)))
sparse <- readRDS('inst/extdata/sparse-mds-comparison/benchmark.rds')
for(id in c('saddle_graph-64','saddle-4096')) {
  message('Registering sparse ',id)
  case <- sparse$cases[[id]]
  selected <- Filter(function(f) f$row$case==id && f$row$seed==1,sparse$fits)
  fits <- setNames(selected,vapply(selected,function(f) f$row$method,''))
  fits <- lapply(fits,function(f) {f$coords <- f$coords*f$row$relative_scale;f})
  case$triangles <- sparse.comparison.surface.triangles(case)
  invisible(comparison.surface.saved(case,fits,kind='sparse'))
}
message('Saved all display registrations without changing fitted layouts.')
# Navigation thumbnails show the saved MDS layouts, without screenshots of controls.
bundle <- readRDS('inst/extdata/graph-examples/layouts.rds')
for(id in c('dwt_1005','1138_bus','carpet','porous')) {
  png(file.path('inst','extdata','documentation-display',paste0(id,'.png')),
      width=600,height=330,res=110)
  par(mar=c(0,0,0,0))
  z <- grip::project.3d(bundle$fits[[id]]$mds$coords)
  edges <- bundle$cases[[id]]$edges
  plot(z,type='n',asp=1,axes=FALSE,xlab='',ylab='')
  segments(z[edges[,1],1],z[edges[,1],2],z[edges[,2],1],z[edges[,2],2],
           col='#1769AA80',lwd=.6)
  points(z,pch=16,cex=.3,col='#1769AA')
  dev.off()
}
