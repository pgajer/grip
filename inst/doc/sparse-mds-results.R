## ----setup, include=FALSE-----------------------------------------------------
library(grip)
knitr::opts_chunk$set(dev="jpeg", dev.args=list(quality=90), fig.retina=1, dpi=96, collapse=TRUE,comment='#>',out.width='100%')
source(system.file('scripts','vignette-views.R',package='grip'),local=TRUE)
source(system.file('scripts','metric-mds-surface-alignment.R',package='grip'),local=TRUE)
source(system.file('scripts','sparse-mds-comparison.R',package='grip'),local=TRUE)
bundle <- readRDS(system.file('extdata','sparse-mds-comparison','benchmark.rds',package='grip'))
rows <- sparse_comparison_rows(bundle)
summary <- sparse_comparison_summary(bundle)
have.ivue <- isTRUE(getOption('grip.vignette.interactive',FALSE)) &&
  requireNamespace('ivue',quietly=TRUE) && requireNamespace('rgl',quietly=TRUE) &&
  'layer3D.callback' %in% getNamespaceExports('ivue') &&
  'limits' %in% names(formals(ivue::plot3D.plain))
method_label <- function(x) ifelse(x=='full','Full SGD',ifelse(x=='grip','Weighted GRIP',
  paste('Sparse SGD:',sub('sparse','',x),'pivots')))

## ----sparse-view-saddle_graph-64, echo=FALSE----------------------------------
sparse_comparison_view(case,fits,saved.display=TRUE)

## ----sparse-view-paraboloid_graph-64, echo=FALSE------------------------------
sparse_comparison_view(case,fits,saved.display=TRUE)

## ----sparse-view-helix_graph-64, echo=FALSE-----------------------------------
sparse_comparison_view(case,fits,saved.display=TRUE)

## ----sparse-gallery, results='asis', echo=FALSE-------------------------------
for (id in names(bundle$cases)) {
  case <- bundle$cases[[id]]
  cat('\n\n### ',case$label,' (',case$n,' vertices)\n\n',sep='')
  cat(if(bundle$evaluation[[id]]$mode=='full') 'Error uses all vertex pairs.\n\n' else
    'Error is estimated from 32 independently sampled source vertices. Full SGD was not attempted.\n\n')
  med <- summary[summary$case==id,]; med$method <- method_label(med$method)
  cat(knitr::kable(med[,c('method','returned','attempted','relative_error','profiled_relative_error','seconds','peak_mib')],
    row.names=FALSE,digits=3,col.names=c('Method','Returned','Attempted','Original error','Scale-adjusted error','Total seconds','Peak MiB')),sep='\n')
  selected <- Filter(function(f) f$row$case==id && f$row$seed==1,bundle$fits)
  fits <- stats::setNames(selected,vapply(selected,function(f) f$row$method,''))
  cat('\n\nDisplayed seed-1 scale multipliers: ',paste(vapply(fits,function(f) paste(method_label(f$row$method),format(f$row$relative_scale,digits=3)),''),collapse='; '),'.\n\n',sep='')
  child <- c(sprintf('```{r sparse-view-%s, echo=FALSE}',id),
    'sparse_comparison_view(case,fits,saved.display=TRUE)','```')
  if (id %in% c('dwt_307','saddle-4096')) {
    cat('\n\n[Interactive figure in the main guide](sparse-metric-mds.html).\n\n')
    next
  }
  cat(knitr::knit_child(text=child,envir=environment(),quiet=TRUE),sep='\n')
}
failed <- rows[rows$status!='ok',c('case','seed','method','failure')]
if(nrow(failed)) { cat('\n\n### Unsuccessful attempts\n\n'); cat(knitr::kable(failed),sep='\n') }

