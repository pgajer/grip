## ----setup, include=FALSE-----------------------------------------------------
library(grip)
knitr::opts_chunk$set(dev="jpeg", dev.args=list(quality=90), fig.retina=1, dpi=96, collapse=TRUE,comment='#>',fig.width=9,fig.height=5,out.width='100%')
source(system.file('scripts','vignette-views.R',package='grip'),local=TRUE)

## ----graph--------------------------------------------------------------------
data(hmp.u01.gc.coarse)
g <- hmp.u01.gc.coarse
knitr::kable(data.frame(Quantity=c('Original giant-component vertices',
  'Coarsened vertices','Undirected edges','Selected k'),
  Value=unlist(g$graph_info[c('original_vertices','coarse_vertices','edge_count','selected_k')],use.names=FALSE)))

## ----saved, include=FALSE-----------------------------------------------------
results <- readRDS(system.file('extdata','hmp_u01_gc_coarse','vignette_results.rds',package='grip'))
edges <- do.call(rbind,lapply(seq_along(g$adj_list),function(i) {
  j <- g$adj_list[[i]];j <- j[j>i];cbind(rep(i,length(j)),j)
}))

## ----hmp-view, echo=FALSE-----------------------------------------------------
fresh <- readRDS(system.file('extdata','hmp_u01_gc_coarse','weighted_layouts.rds',package='grip'))
documentation.group.layouts(list(
  `SGD metric MDS`=fresh$sgd$value$coords,
  `Weighted GRIP`=fresh$weighted_grip$value,
  `Hop-based GRIP`=results$layouts$preset$default),
  edges,g$vertex_data$cst,
  description='HMP/U01: CST colors are identical across layouts.')

## ----new-status---------------------------------------------------------------
knitr::kable(fresh$sgd$value$metadata$starts,digits=4)
knitr::kable(data.frame(Method=c("SGD metric MDS","Weighted GRIP"),
  Seconds=c(fresh$sgd$elapsed_seconds,fresh$weighted_grip$elapsed_seconds)),digits=2)

## ----preset-scores------------------------------------------------------------
knitr::kable(results$preset_summary[,c('candidate','sampled.stress.mean',
  'edge.length.cv.mean','cluster.separation.mean','score.composite')],digits=3)

## ----local-scores-------------------------------------------------------------
knitr::kable(results$local_search_summary[,c('candidate','final.rounds',
  'repulsion.factor','sampled.stress.mean','cluster.separation.mean')],digits=3)

## ----new-weighted, eval=FALSE-------------------------------------------------
# mds <- metric.mds(adj.list=g$adj_list,weight.list=g$weight_list,
#   n=length(g$adj_list),dim=3,backend='sgd',init='random',seed=1,
#   max.iter=200,pair.weights='uniform',diagnostics=FALSE)
# coords <- grip(adj.list=g$adj_list,weight.list=g$weight_list,
#   n=length(g$adj_list),dim=3,metric='edge_length',seed=1)
# score.layout(coords,adj.list=g$adj_list,weight.list=g$weight_list,
#   n=length(g$adj_list),clusters=g$vertex_data$cst,edge.crossings='never')

## ----coarsening---------------------------------------------------------------
knitr::kable(g$graph_info$coarsening_rounds)

## ----reproduce, eval=FALSE----------------------------------------------------
# Sys.setenv(GRIP_HMP_RESULTS_OUTPUT='fresh-hmp-results.rds')
# source(system.file('scripts','precompute-hmp-results.R',package='grip'))

