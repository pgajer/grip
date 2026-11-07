## ----setup, include=FALSE-----------------------------------------------------
library(grip)
knitr::opts_chunk$set(dev="jpeg", dev.args=list(quality=90), fig.retina=1, dpi=96, collapse=TRUE, comment='#>', fig.width=9, fig.height=5,
                     out.width='100%')
source(system.file('scripts','vignette-views.R',package='grip'),local=TRUE)

## ----first-layout-------------------------------------------------------------
g <- suitesparse.graph('dwt_307')
fit <- metric.mds(edges=g$edges, n=g$n, dim=3, backend="sgd",
                  seed=1, diagnostics=FALSE)
coords <- fit$coords

## ----first-view, echo=FALSE---------------------------------------------------
documentation.view(list(`SGD metric MDS`=coords),g$edges,
  description='dwt_307: 307 vertices. SGD metric MDS with seed 1.')

## ----plot-recipe, eval=FALSE--------------------------------------------------
# ivue::plot3D.plain(coords, layers=list(ivue::layer3D.edges(g$edges)), controls = FALSE)

## ----score--------------------------------------------------------------------
quality <- score.layout(coords, edges=g$edges, n=g$n,
                        edge.crossings='never')
knitr::kable(quality[,c('sampled.stress','edge.length.cv',
                       'sampled.nonedge.sep.ratio')],digits=3)

## ----metric-mds-default-------------------------------------------------------
mds <- suppressWarnings(metric.mds(
  edges=g$edges, n=g$n, dim=3, init='random',
  max.iter=80, seed=1, diagnostics=FALSE))
mds$metadata$termination

## ----methods-view, echo=FALSE-------------------------------------------------
documentation.view(list('Random start, 80 passes'=mds$coords,'Default initialization, 1000 passes'=coords),g$edges,common.size=TRUE,
  description='dwt_307: two SGD metric-MDS runs with different initialization and pass budgets. Aligned at a common overall size for display. Seed 1.')

