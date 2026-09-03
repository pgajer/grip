#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1)
out <- normalizePath(args[1]); curves <- list(); validation <- list(); selection <- list()
read_ref <- function(r,m,s) {
  x<-as.matrix(read.csv(file.path(out,sprintf("reference-r%02d-m%d-s%d.csv",r,m,s)),header=FALSE))
  list(sources=as.integer(x[,1]),D=x[,-1,drop=FALSE])
}
ref_pairs <- function(ref) {
  grid<-expand.grid(row=seq_along(ref$sources),j=1:1000)
  i<-ref$sources[grid$row];j<-grid$j
  keep<-i!=j & (!(j %in% ref$sources) | i<j)
  list(pairs=cbind(pmin(i,j),pmax(i,j))[keep,,drop=FALSE],
       d=ref$D[as.matrix(grid)][keep])
}
rel <- function(y,d) sqrt(sum((y-d)^2)/sum(d*d))
for(r in 1:5) {
  if(!file.exists(file.path(out,sprintf("reference-r%02d-m81-s128.csv",r)))) next
  graph_data<-readRDS(file.path(out,sprintf("graphs-%02d.rds",r)))
  ref<-read_ref(r,81,128); coarse<-read_ref(r,41,128); pairs<-ref_pairs(ref)
  stopifnot(length(pairs$d)==119744L,all(pairs$d>0))
  checks<-read.csv(file.path(out,sprintf("smooth-checks-%02d.csv",r)))
  observed<-ref$D[cbind(match(checks$i,ref$sources),checks$j)]
  v<-data.frame(replicate=r,reference_pairs=length(pairs$d),
    coarse_fine_relative_difference=rel(coarse$D,ref$D),
    bvp_relative_difference=rel(observed,checks$distance),
    bvp_max_absolute_difference=max(abs(observed-checks$distance)),
    bvp_max_pair_relative_difference=max(abs(observed/checks$distance-1)))
  finer_file<-file.path(out,sprintf("reference-r%02d-m161-s16.csv",r))
  finer<-if(file.exists(finer_file)) read_ref(r,161,16) else NULL
  full_finer_file<-file.path(out,sprintf("reference-r%02d-m161-s128.csv",r))
  full_finer<-if(file.exists(full_finer_file)) read_ref(r,161,128) else NULL
  v$finer_relative_difference<-if(!is.null(finer)) rel(ref$D[match(finer$sources,ref$sources),],finer$D) else NA_real_
  rows<-lapply(graph_data$graphs,function(g) {
    y<-g$distances[pairs$pairs];d<-pairs$d;a<-sum(y*d)/sum(d*d)
    fine<-if(!is.null(finer)) ref_pairs(finer) else NULL
    sub<-function(s) {rr<-list(sources=ref$sources[1:s],D=ref$D[1:s,,drop=FALSE]);pp<-ref_pairs(rr);rel(g$distances[pp$pairs],pp$d)}
    data.frame(replicate=r,k=g$k,edges=nrow(g$edges),bridges=g$bridges,
      components_before=g$components_before,graph_seconds=g$elapsed,
      xg_error=rel(y,d),xg_profiled_error=rel(y,a*d),fitted_scale=a,
      mean_pair_relative_bias=mean(y/d-1),fraction_shorter=mean(y<d),
      coarse_error=rel(g$distances[pairs$pairs],ref_pairs(coarse)$d),
      sources32_error=sub(32),sources64_error=sub(64),
      fine16_error=if(!is.null(fine)) rel(g$distances[fine$pairs],fine$d) else NA_real_,
      fine128_error=if(!is.null(full_finer)) rel(g$distances[pairs$pairs],ref_pairs(full_finer)$d) else NA_real_,
      medium16_error=sub(16))
  })
  d<-do.call(rbind,rows);d<-d[order(d$k),];best<-d$k[which.min(d$xg_error)]
  v$best_k<-best;v$coarse_best_k<-d$k[which.min(d$coarse_error)]
  v$medium16_best_k<-d$k[which.min(d$medium16_error)]
  v$fine16_best_k<-if(!is.null(finer)) d$k[which.min(d$fine16_error)] else NA_integer_
  v$sources32_best_k<-d$k[which.min(d$sources32_error)]
  v$sources64_best_k<-d$k[which.min(d$sources64_error)]
  v$second_best_gap<-sort(d$xg_error)[2]-min(d$xg_error)
  near<-d$k[d$xg_error<=1.01*min(d$xg_error)]
  v$within_one_percent_k_min<-min(near);v$within_one_percent_k_max<-max(near)
  v$fine128_best_k<-if(!is.null(full_finer)) d$k[which.min(d$fine128_error)] else NA_integer_
  v$fine128_relative_difference<-if(!is.null(full_finer)) rel(ref$D,full_finer$D) else NA_real_
  v$search_boundary<-best %in% range(d$k)
  saveRDS(list(replicate=r,selected_k=best,reference=pairs,curve=d,checks=v,
    source_ids=ref$sources,scope="128 random sources, all targets; unordered duplicates removed"),
    file.path(out,sprintf("calibration-%02d.rds",r)))
  curves[[length(curves)+1L]]<-d;validation[[length(validation)+1L]]<-v
  cat(sprintf("cloud=%d best_k=%d XG=%.6f boundary=%s BVP_Rel=%.3g fine_Rel=%.3g\n",r,best,min(d$xg_error),v$search_boundary,v$bvp_relative_difference,v$finer_relative_difference));flush.console()
}
if(length(curves)) {
  write.csv(do.call(rbind,curves),file.path(out,"calibration-curves.csv"),row.names=FALSE)
  write.csv(do.call(rbind,validation),file.path(out,"reference-validation.csv"),row.names=FALSE)
}
