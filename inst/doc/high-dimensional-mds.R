## ----setup, include=FALSE-----------------------------------------------------
library(grip)
knitr::opts_chunk$set(dev="jpeg", dev.args=list(quality=90), fig.retina=1, dpi=96, collapse=TRUE,comment='#>',fig.width=10,fig.height=6,out.width='100%')
source(system.file('scripts','high-dimensional-mds.R',package='grip'),local=TRUE)
bundle <- readRDS(system.file('extdata','high-dimensional-mds','examples.rds',package='grip'))

## ----recipe-input-------------------------------------------------------------
library(grip)
bundle <- readRDS(system.file('extdata','high-dimensional-mds','examples.rds',package='grip'))
example <- bundle$cases[['paraboloid_17_6']]
prepared <- prepare.graph.geodesic.mds(
  edges=example$edges, n=example$n, edge.weights=example$weights)
D <- prepared$distance_matrix
stopifnot(isTRUE(all.equal(unname(D),unname(example$distances),tolerance=1e-10)))

## ----recipe-fit---------------------------------------------------------------
set.seed(17)
start10 <- matrix(rnorm(example$n * 10),example$n,10)
fit <- function(dimension, initial) {
  # Record schedule termination below; suppress its expected warning here.
  withCallingHandlers(
    metric.mds(prepared=prepared, dim=dimension, init=initial,
      backend='sgd', approximation='full', pair.weights='uniform',
      n.init=1, max.iter=20, seed=17, diagnostics=FALSE,
      sgd.control=list(checkpoint.every=1L)),
    warning=function(w) {
      if (grepl('iteration_limit',conditionMessage(w),fixed=TRUE))
        invokeRestart('muffleWarning')
    })
}
direct <- fit(3,start10[,1:3])
high <- fit(10,start10)
pca3 <- prcomp(high$coords,center=TRUE,scale.=FALSE,rank.=3)$x[,1:3,drop=FALSE]
refined <- fit(3,pca3)
classical.start <- fit(3,'classical')

## ----recipe-score-------------------------------------------------------------
normalized.stress <- function(z,D) {
  target <- as.vector(as.dist(D))
  keep <- is.finite(target)
  sum((as.vector(dist(z))[keep]-target[keep])^2)/sum(target[keep]^2)
}
layouts <- list('Direct 3D'=direct$coords, '10D before projection'=high$coords,
  '10D + PCA'=pca3, '10D + PCA + 3D refinement'=refined$coords,
  'Classical + 3D refinement'=classical.start$coords)
knitr::kable(data.frame(Route=names(layouts),
  Dimensions=vapply(layouts,ncol,0L),
  `Normalized stress`=vapply(layouts,normalized.stress,0,D=D),
  check.names=FALSE),digits=7,row.names=FALSE)

## ----recipe-termination, echo=FALSE-------------------------------------------
knitr::kable(data.frame(Stage=c('Direct 3D','10D','3D refinement','Classical-start 3D'),
  Termination=vapply(list(direct,high,refined,classical.start),function(f) f$metadata$termination,'')),
  row.names=FALSE)

## ----ambient-view, echo=FALSE-------------------------------------------------
highdim.view(bundle$cases[['paraboloid_17_ambient']])

## ----ambient-scores, echo=FALSE-----------------------------------------------
highdim.table(highdim.routes(bundle$cases[['paraboloid_17_ambient']]))

## ----graph-view, echo=FALSE---------------------------------------------------
highdim.view(bundle$cases[['paraboloid_17_6']])

## ----graph-scores, echo=FALSE-------------------------------------------------
highdim.table(highdim.routes(bundle$cases[['paraboloid_17_6']]))

## ----graph-dimensions, echo=FALSE---------------------------------------------
highdim.table(highdim.dimensions(bundle$cases[['paraboloid_17_6']]))

## ----illc-view, echo=FALSE----------------------------------------------------
highdim.view(bundle$cases[['HB__illc1033']])

## ----illc-scores, echo=FALSE--------------------------------------------------
highdim.table(highdim.routes(bundle$cases[['HB__illc1033']]))
highdim.table(highdim.dimensions(bundle$cases[['HB__illc1033']]))

## ----lock-view, echo=FALSE----------------------------------------------------
highdim.view(bundle$cases[['HB__lock_700']])

## ----lock-scores, echo=FALSE--------------------------------------------------
highdim.table(highdim.routes(bundle$cases[['HB__lock_700']]))
highdim.table(highdim.dimensions(bundle$cases[['HB__lock_700']]))

## ----contrast, echo=FALSE, results='asis'-------------------------------------
a <- highdim.dimensions(bundle$cases[['HB__illc1033']])
b <- highdim.dimensions(bundle$cases[['HB__lock_700']])
cat(sprintf('\nAt intermediate dimension 10, mean paired time-matched gain is **%.2f%%** for HB/illc1033 and **%.2f%%** for HB/lock_700. These examples show why an additional fitting stage should be evaluated at its total cost, rather than compared only with one short direct stage.\n',
  a$`Mean paired gain (%)`[a$Dimension==10],b$`Mean paired gain (%)`[b$Dimension==10]))

