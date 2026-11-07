## ----setup, include=FALSE-----------------------------------------------------
library(grip)
knitr::opts_chunk$set(dev="jpeg", dev.args=list(quality=90), fig.retina=1, dpi=96, collapse=TRUE, comment='#>', fig.width=9, fig.height=5,
                     out.width='100%')
source(system.file('scripts','vignette-views.R',package='grip'),local=TRUE)

## ----candidates---------------------------------------------------------------
g <- suitesparse.graph('494_bus')
cmp <- compare.layouts(edges=g$edges,n=g$n,dim=3,
  candidates=c('default','tree','mesh'),seeds=1:3,
  sample.size.stress=2000L,sample.size.nonedge=3000L,
  edge.crossings='never',return.layouts=TRUE)
knitr::kable(cmp$summary[,c('candidate','sampled.stress.mean',
  'edge.length.cv.mean','sampled.nonedge.sep.ratio.mean',
  'stability.procrustes.mean')],digits=3)

## ----candidate-view, echo=FALSE-----------------------------------------------
layouts <- lapply(cmp$layouts,function(z) z[['1']])
names(layouts) <- paste('GRIP:',names(layouts))
documentation.view(layouts,g$edges,
  description='494_bus: three GRIP presets, seed 1 in every view. Rigid alignment; no resizing.')

## ----runs---------------------------------------------------------------------
knitr::kable(cmp$runs[,intersect(c('candidate','seed','status','sampled.stress'),
                               names(cmp$runs))],digits=3)

