## ----setup, include=FALSE-----------------------------------------------------
library(grip)
knitr::opts_chunk$set(dev="jpeg", dev.args=list(quality=90), fig.retina=1, dpi=96, collapse=TRUE,comment='#>',fig.width=9,fig.height=5,out.width='100%')
source(system.file('scripts','vignette-views.R',package='grip'),local=TRUE)

## ----recipe, eval=FALSE-------------------------------------------------------
# g <- suitesparse.graph('dwt_307')
# sgd <- metric.mds(edges=g$edges,n=g$n,dim=3,seed=1,
#                   init='random',max.iter=100,diagnostics=FALSE)
# smacof.fit <- metric.mds(edges=g$edges,n=g$n,dim=3,seed=1,
#                   init='random',max.iter=100,diagnostics=FALSE,
#                   backend='smacof',eps=1e-8)

## ----data, include=FALSE------------------------------------------------------
bundle <- readRDS(system.file('extdata','metric-mds-comparison','benchmark.rds',package='grip'))

## ----selected-saddle-144, echo=FALSE------------------------------------------
documentation.full.view(case,fits)

## ----selected-saddle_graph-144, echo=FALSE------------------------------------
documentation.full.view(case,fits)

## ----selected-saddles, results='asis', echo=FALSE-----------------------------
for(id in c('saddle-144','saddle_graph-144')) {
  case <- bundle$cases[[id]]; fits <- comparison_representative(bundle,id)
  cat('\n\n### ',case$target,'\n\n',sep='')
  child <- c(sprintf('```{r selected-%s, echo=FALSE}',id),
    'documentation.full.view(case,fits)','```')
  cat(knitr::knit_child(text=child,envir=environment(),quiet=TRUE),sep='\n')
  tab <- do.call(rbind,lapply(fits,function(f) f$row[,c('backend','error','seconds','iterations','termination')]))
  cat(knitr::kable(tab,digits=4,row.names=FALSE,
    col.names=c('Backend','Normalized error','Fit seconds','Passes/iterations','Termination')),sep='\n')
}

## ----tradeoff, echo=FALSE-----------------------------------------------------
comparison_tradeoff(bundle,c('saddle_graph-64','saddle_graph-144'))

## ----weighting-recipe, eval=FALSE---------------------------------------------
# fit <- metric.mds(edges=g$edges,n=g$n,dim=3,
#   pair.weights='inverse_squared',backend='sgd',seed=1)

## ----weighting, echo=FALSE----------------------------------------------------
w <- readRDS(system.file('extdata','metric-mds-weighting','benchmark.rds',package='grip'))
case <- w$cases[['saddle_graph-64']]
fits <- comparison_weighting_fits(w,case$id)
documentation.full.view(case,fits,weighting=TRUE)

## ----weighting-table, echo=FALSE----------------------------------------------
knitr::kable(do.call(rbind,lapply(fits,function(f)
  f$row[,intersect(c('backend','pair_weights','uniform_error','relative_error','seconds'),names(f$row))])),digits=4)

