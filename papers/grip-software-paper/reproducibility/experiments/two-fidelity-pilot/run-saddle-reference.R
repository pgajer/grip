#!/usr/bin/env Rscript
# From the repository root:
# Rscript papers/grip-software-paper/reproducibility/experiments/two-fidelity-pilot/run-saddle-reference.R output/single-saddle-reference-2026-09-04
# Uses the local package source, not a replacement of the installed CRAN package.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)>=1L)
out <- args[1]; dir.create(out,recursive=TRUE,showWarnings=FALSE)
pkgload::load_all(".",quiet=TRUE)
home <- "papers/grip-software-paper/reproducibility/experiments/two-fidelity-pilot"
source(file.path(home,"score-saddle-reference.R"))
n <- if(length(args)>1L) as.integer(args[2]) else 500L
seed <- 1L; C <- .8; half.width <- 1; ks <- 3:20; max.iter <- 200L
cache <- file.path(out,"layouts.rds")
settings <- list(n=n,seed=seed,C=C,half.width=half.width,ks=ks,max.iter=max.iter)
if(file.exists(cache)) {
  fit <- readRDS(cache); stopifnot(identical(fit$settings,settings))
  X <- fit$X; graphs <- fit$graphs; mds <- fit$mds; mds.edge.kk <- fit$mds.edge.kk
} else {
  set.seed(seed); xy <- matrix(numeric(),ncol=2)
  while(nrow(xy)<n) {
    proposal <- matrix(runif(4000,-half.width,half.width),ncol=2)
    area <- sqrt(1+4*C^2*rowSums(proposal^2))
    accept <- runif(nrow(proposal))<area/sqrt(1+8*C^2*half.width^2)
    xy <- rbind(xy,proposal[accept,,drop=FALSE])
  }
  xy <- xy[seq_len(n),,drop=FALSE]
  X <- cbind(x=xy[,1],y=xy[,2],z=C*(xy[,1]^2-xy[,2]^2))
  graphs <- setNames(lapply(ks,function(k) dgraphs::create.sknn.graph(X,k=k,
    neighbor.method="exact",edge.weight="distance",connect.components=TRUE,
    connect.method="component.mst",prune.method="none",prune.edges=FALSE)),as.character(ks))
  mds <- lapply(graphs,function(g) grip::metric.mds(edges=g$edge_matrix,
    edge_weights=g$edge_weight,n=n,dim=3,diagnostics=FALSE)$coords)
  mds.edge.kk <- setNames(lapply(ks,function(k) {
    message("Fitting k = ",k)
    g <- graphs[[as.character(k)]]
    grip::edge.kk(coords=mds[[as.character(k)]],edges=g$edge_matrix,
      edge_weights=g$edge_weight,n=n,dim=3,max_iter=max.iter,stiffness_method="density",
      density_mix_schedule=c(0,.25,.5,.75,1),scale_mode="profiled",edge_length_epsilon=0,
      diagnostics=FALSE,return_trace=FALSE,seed=seed)$coords
  }),as.character(ks))
  saveRDS(list(settings=settings,X=X,graphs=graphs,mds=mds,mds.edge.kk=mds.edge.kk),cache)
}
result <- score.saddle.reference(X,graphs,mds,mds.edge.kk,C)
write.csv(result$scores,file.path(out,"scores.csv"),row.names=FALSE)
saveRDS(result,file.path(out,"reference-scores.rds"))
# Denser surface discretization and more samples at the user's selected k=10.
convergence <- score.saddle.reference(X,graphs["10"],mds["10"],mds.edge.kk["10"],C,
                                     sample_sizes=c(8000L,32000L),subdivisions=3L)
write.csv(convergence$scores,file.path(out,"convergence-k10.csv"),row.names=FALSE)
writeLines(capture.output(list(settings=settings,R=R.version.string,
  packages=vapply(c("grip","dgraphs","geometry"),function(p) as.character(packageVersion(p)),""),
  BLAS=extSoftVersion()["BLAS"],RNG=RNGkind())),file.path(out,"environment.txt"))
source(file.path(home,"plot-saddle-reference.R"))
plot.saddle.reference(result$scores,file.path(out,"reference-comparison.png"))
print(subset(result$scores,k==10 & sample_size==8000,
  select=c(method,alignment,path_rel,edge_rel,stress1,coordinate_rmse,coordinate_relative_rmse,
           alignment_scale,surface_rms,surface_mc_se)),row.names=FALSE)
