## ----setup, include=FALSE-----------------------------------------------------
library(grip)
knitr::opts_chunk$set(dev="jpeg", dev.args=list(quality=90), fig.retina=1, dpi=96, collapse=TRUE, comment='#>', fig.width=9, fig.height=5,
                     out.width='100%')
source(system.file('scripts','vignette-views.R',package='grip'),local=TRUE)

## ----torus-fits---------------------------------------------------------------
g <- torus.surface.graph(12,8,normalize='none')
hop <- grip(edges=g$edges,n=g$n,edge.weights=g$edge_weights,
  dim=3,metric='hop',seed=1)
weighted <- grip(edges=g$edges,n=g$n,edge.weights=g$edge_weights,
  dim=3,metric='edge_length',length.normalization='none',seed=1)

## ----torus-view, echo=FALSE---------------------------------------------------
documentation.view(list('GRIP: hop'=hop,'GRIP: edge lengths'=weighted),
  g$edges,reference=g$coords_surface,common.size=TRUE,
  description='Torus: identical input lengths in both fits. Gray: generating coordinates. Aligned at a common overall size for display.')

## ----weighted-scores----------------------------------------------------------
fits <- list(Hop=hop,'Edge lengths'=weighted)
scores <- do.call(rbind,lapply(fits,function(z)
  score.layout(z,edges=g$edges,n=g$n,edge.weights=g$edge_weights,
    sample.size.stress=2000L,edge.crossings='never')))
knitr::kable(scores[,c('sampled.stress','edge.length.cv',
                     'sampled.nonedge.sep.ratio')],digits=3)

## ----tree---------------------------------------------------------------------
tree <- kary.tree.weighted.graph(k=2,depth=5,depth.decay=.85,
  branch.spread=.2,normalize='none')
tree.coords <- grip(edges=tree$edges,n=tree$n,edge.weights=tree$edge_weights,
  dim=3,metric='edge_length',length.normalization='none',seed=1)

## ----tree-view, echo=FALSE----------------------------------------------------
documentation.view(list('Weighted GRIP'=tree.coords),tree$edges,
  description='An intrinsic weighted tree: shorter branches at greater depth. No generating reference.')

