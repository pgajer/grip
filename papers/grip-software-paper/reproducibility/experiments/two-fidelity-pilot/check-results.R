#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)>=1)
reps<-if(length(args)>1) as.integer(strsplit(args[2],",")[[1]]) else 1:5
out<-normalizePath(args[1]);curves<-read.csv(file.path(out,"calibration-curves.csv"))
validation<-read.csv(file.path(out,"reference-validation.csv"));rows<-list()
for(r in reps) {
  cloud<-readRDS(file.path(out,sprintf("cloud-%02d.rds",r)))
  x<-cloud$coords;stopifnot(nrow(x)==1000,ncol(x)==3,anyDuplicated(x)==0L,
    max(abs(x[,1:2]))<=1,max(abs(x[,3]-.8*(x[,1]^2-x[,2]^2)))<1e-14)
  cal<-readRDS(file.path(out,sprintf("calibration-%02d.rds",r)))
  pairs<-cal$reference$pairs;dx<-cal$reference$d
  stopifnot(nrow(pairs)==119744,all(pairs[,1]<pairs[,2]),anyDuplicated(pairs)==0L)
  chords<-sqrt(rowSums((x[pairs[,1],]-x[pairs[,2],])^2))
  stopifnot(min(dx-chords)>-1e-10)
  d<-curves[curves$replicate==r,];stopifnot(identical(as.integer(d$k),3:80))
  f<-readRDS(file.path(out,sprintf("fit-%02d-k%d.rds",r,cal$selected_k)))
  stopifnot(f$k==d$k[which.min(d$xg_error)],nrow(f$prepared$pair_matrix)==499500,
    length(f$prepared$flat_pair_edge_offsets)==499501,
    tail(f$prepared$flat_pair_edge_offsets,1)==length(f$prepared$flat_edge_u),
    all(f$prepared$flat_edge_coeff==1),
    all(vapply(f$candidates,ncol,integer(1))==3L),
    all(is.finite(as.matrix(f$scores[,-3]))),f$scores$path_rel[1]<1e-7)
  g<-f$graph;dg<-g$distances[pairs]
  xg<-sqrt(sum((dg-dx)^2)/sum(dx^2))
  stopifnot(abs(xg-min(d$xg_error))<1e-14,abs(xg-f$scores$xg_error[1])<1e-14)
  # A second implementation of all fixed-path sums, from the saved offsets.
  p<-f$prepared;path_agreement<-0;edge_agreement<-0;stress_agreement<-0
  for(m in 1:3) {
    z<-f$candidates[[m]]
    ell<-sqrt(rowSums((z[p$flat_edge_u+1L,]-z[p$flat_edge_v+1L,])^2))
    cumul<-c(0,cumsum(ell));y<-diff(cumul[p$flat_pair_edge_offsets+1L])
    target<-p$pair_graph_distance;a<-sum(y*target)/sum(target^2)
    path<-sqrt(sum((y-a*target)^2)/sum((a*target)^2))
    edge<-sqrt(rowSums((z[p$edges[,1],]-z[p$edges[,2],])^2));a<-sum(edge*p$edge_targets)/sum(p$edge_targets^2)
    edge<-sqrt(sum((edge-a*p$edge_targets)^2)/sum((a*p$edge_targets)^2))
    chord<-sqrt(rowSums((z[p$pair_matrix[,1],]-z[p$pair_matrix[,2],])^2));a<-sum(chord*target)/sum(target^2)
    stress<-sqrt(sum((chord-a*target)^2)/sum(chord^2))
    path_agreement<-max(path_agreement,abs(path-f$scores$path_rel[m]))
    edge_agreement<-max(edge_agreement,abs(edge-f$scores$edge_rel[m]))
    stress_agreement<-max(stress_agreement,abs(stress-f$scores$stress1[m]))
  }
  stopifnot(path_agreement<1e-8,edge_agreement<1e-12,stress_agreement<1e-12)
  rows[[r]]<-data.frame(replicate=r,reference_pairs=nrow(pairs),layout_pairs=499500,
    minimum_surface_minus_chord=min(dx-chords),path_score_difference=path_agreement,
    edge_score_difference=edge_agreement,stress_score_difference=stress_agreement,
    strict_distance_max_difference=f$strict_distance_max_difference)
  cat("Checked cloud",r,"\n");flush.console()
}
write.csv(do.call(rbind,rows),file.path(out,"independent-validation.csv"),row.names=FALSE)
print(do.call(rbind,rows))
