# Independent all-pair identity versus both the R evaluator and native trace.
out<-"output/mds-edge-kk-radius";files<-list.files(file.path(out,"fits-v2"),pattern="k(239|479)-result[.]rds$",recursive=TRUE,full.names=TRUE)
stopifnot(length(files)>0)
rows<-lapply(files,function(file) {
 r<-readRDS(file)$value;p<-readRDS(sub("-result.rds","-prepared.rds",file,fixed=TRUE))$value$p
 z<-r$candidates$stress/r$scale_unit
 state<-grip:::grip.edge.isometric.energy.gradient(z,p$edges,p$edge_targets,rep(1,nrow(p$edges)),scale=1,edge_length_epsilon=0)
 chord<-sqrt(rowSums((z[p$pair_matrix[,1],]-z[p$pair_matrix[,2],])^2))
 direct<-matrix(0,p$n,p$n);direct[p$edges]<-p$edge_targets;direct<-direct+t(direct)
 raw<-sum((chord-direct[p$pair_matrix])^2)
 native<-readRDS(sub("-result.rds","-stress_fixed_uniform.rds",file,fixed=TRUE))$value$trace$energy[1]
 discrepancy<-abs(2*state$energy-raw);native_discrepancy<-abs(2*native-raw)
 stopifnot(discrepancy<1e-8,native_discrepancy<1e-8)
 strict_raw<-sum((chord-p$distance_matrix[p$pair_matrix])^2)
 data.frame(case=r$scores$case[1],regime=r$scores$regime[1],twice_edge_energy=2*state$energy,
  raw_stress_direct_targets=raw,absolute_difference=discrepancy,native_trace_difference=native_discrepancy,
  strict_target_stress_difference=abs(strict_raw-raw),max_strict_direct_target_difference=max(abs(p$distance_matrix-direct)))
})
write.csv(do.call(rbind,rows),file.path(out,"objective-validation.csv"),row.names=FALSE)
cat("Identity scale, uniform stiffness:",length(rows),"complete-graph R/native controls passed (2E = raw stress on the same direct targets).\n")
