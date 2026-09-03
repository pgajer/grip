#!/usr/bin/env Rscript
# Separate from the manuscript and from the earlier 2D flattening study.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
out <- normalizePath(args[1],mustWork=FALSE)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
RNGkind("Mersenne-Twister","Inversion","Rejection")
for(r in 1:5) {
  set.seed(2211000L+r)
  uv <- matrix(numeric(),ncol=2)
  while(nrow(uv)<1000L) {
    a <- matrix(runif(4000,-1,1),ncol=2)
    accept <- runif(nrow(a)) < sqrt(1+2.56*rowSums(a*a))/sqrt(6.12)
    uv <- rbind(uv,a[accept,,drop=FALSE])
  }
  uv <- uv[1:1000,,drop=FALSE]
  x <- cbind(x=uv[,1],y=uv[,2],z=.8*(uv[,1]^2-uv[,2]^2))
  # Fixed random source subset available if reference timing requires it.
  set.seed(3211000L+r); sources <- sort(sample.int(1000,128))
  saveRDS(list(coords=x,seed=2211000L+r,sources=sources,session=sessionInfo()),
          file.path(out,sprintf("cloud-%02d.rds",r)))
  write.csv(x,file.path(out,sprintf("cloud-%02d.csv",r)),row.names=FALSE)
  write.table(sources-1L,file.path(out,sprintf("sources-%02d.txt",r)),row.names=FALSE,col.names=FALSE)
}
