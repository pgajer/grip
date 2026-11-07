## ----setup, include=FALSE-----------------------------------------------------
library(grip)
knitr::opts_chunk$set(dev="jpeg", dev.args=list(quality=90), fig.retina=1, dpi=96, collapse=TRUE,comment='#>',fig.width=9,fig.height=5,out.width='100%')
source(system.file('scripts','vignette-views.R',package='grip'),local=TRUE)

## ----sparse-recipe, eval=FALSE------------------------------------------------
# g <- suitesparse.graph('dwt_307')
# fit <- metric.mds(edges=g$edges,n=g$n,edge.weights=g$edge_weights,
#   dim=3,approximation='sparse',backend='sgd',seed=1,
#   sparse.control=list(n.pivots=32),max.iter=100)
# ivue::plot3D.plain(fit$coords, layers=list(ivue::layer3D.edges(g$edges)), controls = FALSE)

## ----observation-sparse-recipe, eval=FALSE------------------------------------
# prep <- prepare.sparse.mds(X, method="euclidean", n.pivots=200,
#   q=10, neighbor.variance.target=0.90, seed=1)
# fit <- metric.mds(prepared=prep, approximation="sparse", backend="sgd",
#   dim=3, max.iter=100, seed=1,
#   sgd.control=list(scheduler="zheng", schedule.epsilon=0.1))

## ----linf-sparse-recipe, eval=FALSE-------------------------------------------
# stopifnot(all(is.finite(X)), all(X >= 0), all(apply(X, 1, max) > 0))
# Xmax <- X / apply(X, 1, max)
# # Consolidate any duplicate normalized rows and retain their sample mapping.
# prep <- prepare.sparse.mds(Xmax, method="linf", n.pivots=200, q=10, seed=1)

## ----sparse-data, include=FALSE-----------------------------------------------
source(system.file('scripts','sparse-mds-comparison.R',package='grip'),local=TRUE)
bundle <- readRDS(system.file('extdata','sparse-mds-comparison','benchmark.rds',package='grip'))
summary <- sparse_comparison_summary(bundle)
method_label <- function(x) ifelse(x=='full','Full SGD',ifelse(x=='grip','Weighted GRIP',
  paste('Sparse SGD:',sub('sparse','',x),'pivots')))

## ----sparse-view-dwt_307, echo=FALSE------------------------------------------
sparse_comparison_view(case,fits,saved.display=TRUE)

## ----sparse-view-saddle-4096, echo=FALSE--------------------------------------
sparse_comparison_view(case,fits,saved.display=TRUE)

## ----selected-sparse, results='asis', echo=FALSE------------------------------
for(id in c('dwt_307','saddle-4096')) {
  case <- bundle$cases[[id]]
  cat('\n\n### ',case$label,' (',case$n,' vertices)\n\n',sep='')
  med <- summary[summary$case==id,];med$method <- method_label(med$method)
  cat(knitr::kable(med[,c('method','returned','attempted','relative_error',
    'profiled_relative_error','seconds','peak_mib')],digits=3,row.names=FALSE,
    col.names=c('Method','Returned','Attempted','Original error','Scale-adjusted error','Total seconds','Peak MiB')),sep='\n')
  selected <- Filter(function(f) f$row$case==id && f$row$seed==1,bundle$fits)
  fits <- setNames(selected,vapply(selected,function(f) f$row$method,''))
  child <- c(sprintf('```{r sparse-view-%s, echo=FALSE}',id),
    'sparse_comparison_view(case,fits,saved.display=TRUE)','```')
  cat(knitr::knit_child(text=child,envir=environment(),quiet=TRUE),sep='\n')
  cat('\n\nSeed-1 display scale multipliers: ',paste(vapply(fits,function(f)
    paste(method_label(f$row$method),format(f$row$relative_scale,digits=3)),''),collapse='; '),'.\n\n',sep='')
}

## ----explicit-constraints-----------------------------------------------------
constraints <- list(pairs = rbind(c(1, 2), c(2, 3), c(1, 3)),
                    targets = c(1, 1, 2))
fit <- metric.mds(n = 3, approximation = "sparse", constraints = constraints,
                  seed = 1, max.iter = 30)
fit$metadata$constraint_source

