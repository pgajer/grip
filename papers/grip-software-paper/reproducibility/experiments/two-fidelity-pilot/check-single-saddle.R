#!/usr/bin/env Rscript
# Fast checks; --integration additionally runs all five stages on 80 points.
file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
source(file.path(dirname(normalizePath(file)), "single-saddle-ivue.R"))
for (args in list(list(n=7), list(n=20, ks=20), list(reference.sources=1),
                 list(reference.grids=c(81,41)), list(edge.alpha=2), list(max.iter=0))) {
  stopifnot(inherits(try(do.call(single.saddle.config,args),silent=TRUE),"try-error"))
}
stopifnot(identical(single.saddle.config(n=20)$ks,3:19))
dir <- tempfile("single-saddle-check-"); dir.create(dir)
cfg <- single.saddle.config(out=dir)
set.seed(19); prior <- .Random.seed; kind <- RNGkind()
cloud <- sample.single.saddle(cfg)
stopifnot(identical(prior,.Random.seed), identical(kind,RNGkind()),
  nrow(cloud$coords)==1000L, length(cloud$sources)==128L,
  max(abs(cloud$coords[,3]-.8*(cloud$coords[,1]^2-cloud$coords[,2]^2)))==0)
paper <- readRDS(file.path(.single.saddle.source,"../../precomputed/two-fidelity-saddle.rds"))
stopifnot(paper$representative==5, identical(cloud$coords,paper$coords))
for (n in c(8L,37L,100L)) for (s in c(2L,n)) {
  D <- abs(outer(seq_len(n),seq_len(n),"-"))
  pairs <- .ss.reference.pairs(list(sources=seq_len(s),D=D[seq_len(s),,drop=FALSE]))
  stopifnot(nrow(pairs$pairs)==s*(n-1)-choose(s,2), !anyDuplicated(as.data.frame(pairs$pairs)),
    all(pairs$pairs[,1]<pairs$pairs[,2]), identical(as.numeric(D[pairs$pairs]),as.numeric(pairs$d)))
}
stopifnot(.ss.profile(3*(1:10),1:10)==0)
cat("Input validation, RNG preservation, exact published-cloud reproduction, and pair bookkeeping passed.\n")

if ("--integration" %in% commandArgs(TRUE)) {
  result <- run.single.saddle(single.saddle.config(n=80, ks=c(1,3,6,10,20),
    reference.sources=8,reference.grids=c(11,21),fine.grid=31,fine.sources=4,
    max.iter=8,audit.iter=12,fit.ks=c(3,20),out=tempfile("single-saddle-integration-")),open=FALSE)
  stopifnot(result$calibration$boundary, result$calibration$curve$bridges[1]>0,
    length(result$fits)==2L, length(result$views$widgets)>6L)
  for (fit in result$fits) stopifnot(fit$pairs==choose(80,2),
    fit$scores$path_rel[1]<1e-7,fit$scores$edge_rel[1]<1e-12,
    max(fit$scores$package_score_difference)<1e-10)
  scene <- attr(result$views$widgets[["k-020-mds-edge-kk-graph"]],"ivue")$scene
  lines <- Filter(function(o) identical(o$type,"lines"),scene$objects)
  # rgl stores material alpha at 8-bit precision.
  stopifnot(any(vapply(lines,function(o) length(o$material$alpha)==1L &&
    abs(o$material$alpha-.03)<=1/255,logical(1))))
  cat("All five stages and explicit edge opacity passed; artifacts: ",result$config$out,"\n",sep="")
}
