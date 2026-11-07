## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(fig.retina=1, dpi=72, collapse=TRUE,comment='#>',fig.width=8,fig.height=5,warning=FALSE,message=FALSE)
available <- all(vapply(c('dgraphs','igraph','plotly'),requireNamespace,logical(1),quietly=TRUE))
available <- available && file.exists(system.file('examples','fermat-surfaces.R',package='dgraphs')) &&
  'return.graph' %in% names(formals(dgraphs::fermat.distances)) &&
  'distance.matrix' %in% names(formals(grip::metric.mds))
run.experiment <- identical(Sys.getenv('GRIP_RUN_FERMAT_EXPERIMENT'), 'true') ||
  nzchar(Sys.getenv('FERMAT_SURFACE_RESULTS'))
knitr::opts_chunk$set(eval=available && run.experiment)

## ----dependencies,echo=FALSE,eval=TRUE,results='asis'-------------------------
if(!available) cat('This executable example requires development versions of dgraphs and grip with the direct distance and supporting-graph APIs, plus igraph and plotly. Install these dependencies to render the calculated examples.\n')
if(available && !run.experiment) cat('The reference-geodesic calculations are an explicit reproduction step. Set GRIP_RUN_FERMAT_EXPERIMENT=true to run them, or FERMAT_SURFACE_RESULTS to a directory containing the saved results. Routine package builds show the code without rerunning this experiment.\n')

## ----dense-example------------------------------------------------------------
# library(grip)
# s <- dgraphs::sample.synthetic.geometry(
#   dgraphs::synthetic.quadform(2,3,forms=list(diag(c(1,-1)))),
#   dgraphs::synthetic.sampling.quadform.lab('disk',mode='area'),
#   n=150,seed=4101)
# X <- s$predictors
# D <- dgraphs::fermat.distances(points=X,p=2,backend='implicit',rooted=FALSE)
# fit <- metric.mds(distance.matrix=D,dim=3,backend='sgd',init='random',
#   seed=3051,max.iter=300,diagnostics=FALSE,
#   sgd.control=list(final.rate=1e-5))
# Y <- fit$coords

## ----shared, include=FALSE----------------------------------------------------
# source(Sys.getenv('FERMAT_SURFACE_HELPER',
#   unset=system.file('examples','fermat-surfaces.R',package='dgraphs')))
# cache.root <- Sys.getenv('FERMAT_SURFACE_RESULTS')
# results <- lapply(c('paraboloid','saddle','helix'),function(shape) {
#   path <- file.path(cache.root,paste0(shape,'_4101.rds'))
#   if(nzchar(cache.root) && file.exists(path)) {
#     cached<-readRDS(path)
#     stopifnot(nrow(cached$reference$pairs)==1000L,
#       shape=='helix' || cached$reference$method=='best_of_six')
#     return(cached)
#   }
#   s <- fermat.surface.sample(shape,150L,4101L)
#   fermat.surface.experiment(s,fermat.surface.reference(s))
# })
# metrics <- do.call(rbind,lapply(results,`[[`,'records'))

## ----comparison-table, echo=FALSE---------------------------------------------
# knitr::kable(metrics[metrics$k==0,c('shape','p','distance_nrmse','embedding_nrmse',
#  'embedded_geodesic_nrmse')],digits=3,col.names=c('Example','p',
#  'Input vs geodesic','Embedding vs input','Embedding vs geodesic'))

## ----k-comparison, echo=FALSE-------------------------------------------------
# old<-par(mfrow=c(1,3),mar=c(4,4,2,1))
# for(shape in unique(metrics$shape)) {
#  z<-metrics[metrics$shape==shape,]
#  plot(c(0,10),c(0,1.3*max(c(0,z$embedded_geodesic_nrmse),na.rm=TRUE)),type='n',
#  xlab='k (0 = complete targets)',ylab='Embedding vs geodesic error',main=shape)
#  for(p in 1:4){q<-z[z$p==p,];lines(q$k,q$embedded_geodesic_nrmse,col=p,pch=p,type='b')}
#  legend('topright',legend=paste('p =',1:4),col=1:4,lty=1,pch=1:4,bty='n',cex=.7)
# };par(old)

## ----coverage,echo=FALSE------------------------------------------------------
# z<-metrics[metrics$components>1,]
# knitr::kable(aggregate(cbind(retained,components,reference_pairs)~shape+p,z,range))

## ----paraboloid,echo=FALSE----------------------------------------------------
# fermat.surface.overlay(results[[1]],p=2,k=0)

## ----saddle,echo=FALSE--------------------------------------------------------
# fermat.surface.overlay(results[[2]],p=2,k=0)

## ----restricted,echo=FALSE----------------------------------------------------
# fermat.surface.overlay(results[[2]],p=2,k=3)

## ----helix,echo=FALSE---------------------------------------------------------
# fermat.surface.overlay(results[[3]],p=2,k=0)

## ----helix-baseline,echo=FALSE------------------------------------------------
# h<-results[[3]]; b<-h$geodesic_baseline
# knitr::kable(data.frame(target='Exact helix arc length',
#  embedding_error=fermat.surface.score(as.matrix(dist(b$coords))[lower.tri(b$target)],
#  b$target[lower.tri(b$target)],rescale=FALSE)['nrmse']),digits=6,row.names=FALSE)
# plotly::plot_ly(x=b$coords[,1],y=b$coords[,2],z=b$coords[,3],type='scatter3d',
#  mode='markers',marker=list(size=3,color=h$sample$latent[,1],colorscale='Viridis'))

## ----chord-audit,echo=FALSE---------------------------------------------------
# chord.audit <- do.call(rbind,lapply(results,function(r) {
#  ij<-r$reference$pairs; g<-r$reference$distances
#  c<-sqrt(rowSums((r$sample$predictors[ij[,1],]-r$sample$predictors[ij[,2],])^2))
#  ok<-is.finite(g)&g>0&c>0; g<-g[ok];c<-c[ok];a<-sum(c*g)/sum(c^2)
#  data.frame(Example=r$sample$shape,Pairs=length(g),Scale=a,
#  Raw.error=sqrt(sum((c-g)^2)/sum(g^2)),Scaled.error=sqrt(sum((a*c-g)^2)/sum(g^2)),
#  Median.ratio=median(g/c),Maximum.ratio=max(g/c))
# }))
# knitr::kable(chord.audit,digits=4)

## ----reference-success,echo=FALSE---------------------------------------------
# knitr::kable(do.call(rbind,lapply(results,function(r) {
#  v<-r$reference$check.distances
#  data.frame(Example=r$sample$shape,Attempted=length(v),Finite=sum(is.finite(v)),
#  Failed=sum(!is.finite(v)),Check=if(r$sample$shape=='helix')'Analytic identity' else if(!is.null(r$reference$check.method))r$reference$check.method else 'None')
# })))

## ----clairaut-comparison,echo=FALSE-------------------------------------------
# r<-results[[1]]; b<-r$reference$distances; h<-r$reference$check.distances
# ok<-is.finite(b)&is.finite(h)&h>0; delta<-(b[ok]-h[ok])/h[ok]
# knitr::kable(data.frame(Pairs=sum(ok),Median.signed=median(delta),
#  Maximum.absolute=max(abs(delta))),digits=6)
# knitr::kable(do.call(rbind,lapply(1:4,function(p){
#  d<-r$matrices[[as.character(p)]][r$reference$pairs]
#  data.frame(p=p,Best.of.six=unname(fermat.surface.score(d,b)['nrmse']),
#  Clairaut=unname(fermat.surface.score(d,h)['nrmse']))})),digits=5)

