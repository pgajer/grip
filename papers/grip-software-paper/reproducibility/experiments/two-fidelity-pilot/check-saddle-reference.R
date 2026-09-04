#!/usr/bin/env Rscript
# Run from the repository root after run-saddle-reference.R.
args <- commandArgs(trailingOnly=TRUE); stopifnot(length(args)==1L)
out <- args[1]
d <- read.csv(file.path(out,"scores.csv"))
c <- read.csv(file.path(out,"convergence-k10.csv"))
fit <- readRDS(file.path(out,"layouts.rds"))
stopifnot(nrow(d)==18L*3L*2L*2L, !anyNA(d), all(d$zero_area_faces==0),
          all(is.finite(as.matrix(d[vapply(d,is.numeric,logical(1))]))),
          all(d$coordinate_rmse>=0),all(d$surface_rms>=0))
base <- subset(d,method=="Original saddle")
stopifnot(max(base$coordinate_rmse)<1e-12,max(base$path_rel)<1e-7,max(base$edge_rel)<1e-12)
# Independently reproduce the sampled-saddle experiment's similarity Procrustes formula.
old.procrustes <- function(z,x) {
  z <- sweep(z,2,colMeans(z),"-"); xx <- sweep(x,2,colMeans(x),"-")
  f <- svd(crossprod(z,xx))
  aligned <- sum(f$d)/sum(z*z)*z%*%f$u%*%t(f$v)
  sqrt(sum((aligned-xx)^2)/sum(xx^2))
}
errors <- numeric()
for (key in names(fit$mds)) for(method in c("Metric MDS","MDS + edge-KK")) {
  Z <- if(method=="Metric MDS") fit$mds[[key]] else fit$mds.edge.kk[[key]]
  reported <- d[d$k==as.integer(key) & d$method==method & d$alignment=="similarity",]
  errors <- c(errors,abs(reported$coordinate_relative_rmse-old.procrustes(Z,fit$X)))
}
stopifnot(max(errors)<1e-12)
sampling <- merge(subset(d,sample_size==2000),subset(d,sample_size==8000),
                  by=c("k","method","alignment"),suffixes=c(".small",".large"))
sampling$absolute_change <- abs(sampling$surface_rms.small-sampling$surface_rms.large)
sampling$relative_change <- sampling$absolute_change/pmax(sampling$surface_rms.large,1e-12)
refinement <- merge(subset(d,k==10 & sample_size==8000),subset(c,sample_size==32000),
                    by=c("k","method","alignment"),suffixes=c(".main",".fine"))
refinement$absolute_change <- abs(refinement$surface_rms.main-refinement$surface_rms.fine)
refinement$relative_change <- refinement$absolute_change/pmax(refinement$surface_rms.fine,1e-12)
write.csv(sampling,file.path(out,"sampling-comparison.csv"),row.names=FALSE)
write.csv(refinement,file.path(out,"refinement-comparison.csv"),row.names=FALSE)
report <- capture.output({
  cat("Complete single-cloud grid: 216 rows; no nonfinite values or zero-area faces.\n")
  cat("Maximum discrepancy from historical similarity Procrustes formula:",max(errors),"\n")
  cat("Surface sampling relative changes (2,000 to 8,000 per direction):\n")
  print(quantile(sampling$relative_change,c(0,.5,.95,1)))
  cat("k=10: main mesh/8,000 samples versus denser mesh/32,000 samples:\n")
  print(refinement[,c("method","alignment","surface_rms.main","surface_rms.fine","absolute_change","relative_change")],row.names=FALSE)
  cat("Monte Carlo variation and mesh discretization are distinct from between-cloud uncertainty.\n")
})
writeLines(report,file.path(out,"validation.txt")); cat(report,sep="\n")
