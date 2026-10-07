#!/usr/bin/env Rscript
options(rgl.useNULL=TRUE)
pkgload::load_all('.',quiet=TRUE,export_all=FALSE,helpers=FALSE)
output <- file.path(getwd(),'output','vignette-previews','doc')
dir.create(output,recursive=TRUE,showWarnings=FALSE)
rmarkdown::render('vignettes/high-dimensional-mds.Rmd',output_dir=output,
  envir=new.env(parent=globalenv()),quiet=FALSE)
