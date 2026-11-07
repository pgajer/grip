## ----setup, include=FALSE-----------------------------------------------------
library(grip)
knitr::opts_chunk$set(dev="jpeg", dev.args=list(quality=90), fig.retina=1, dpi=96, collapse=TRUE,comment='#>',fig.width=9,fig.height=5,out.width='100%')
source(system.file('scripts','vignette-views.R',package='grip'),local=TRUE)

## ----compare-example----------------------------------------------------------
g <- suitesparse.graph('dwt_307')
cmp <- compare.layouts(edges=g$edges,n=g$n,dim=3,
  candidates=c('default','tree'),seeds=1:2,return.layouts=TRUE,
  edge.crossings='never')
adj <- vector('list',g$n)
for(i in seq_len(nrow(g$edges))) {
  a <- g$edges[i,1];b <- g$edges[i,2]
  adj[[a]] <- c(adj[[a]],b);adj[[b]] <- c(adj[[b]],a)
}
project <- gripui.project.from.compare(cmp,
  graph=list(adj_list=adj),title='dwt_307: presets and seeds')
gripui.validate.project(project)

## ----launch, eval=FALSE-------------------------------------------------------
# run.gripui(project)

## ----project-view, echo=FALSE-------------------------------------------------
documentation.view(list('Default: seed 1'=cmp$layouts$default[['1']],
  'Tree: seed 1'=cmp$layouts$tree[['1']]),g$edges,
  description='Saved project layouts: default and tree GRIP presets. Rigid alignment; no resizing.')

## ----project-table------------------------------------------------------------
knitr::kable(project$layouts[,c('layout_id','candidate','seed','status','availability')])

## ----families, eval=FALSE-----------------------------------------------------
# run.gripui.family()

## ----from-directory, eval=FALSE-----------------------------------------------
# project <- gripui.project.from.dir(root='/path/to/saved/outputs',
#   title='My saved search')
# gripui.validate.project(project)
# run.gripui(project)

