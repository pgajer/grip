# Independent all-pair identity versus the package's fixed-scale edge evaluator.
out<-"output/mds-edge-kk-radius";files<-list.files(file.path(out,"fits-v2"),pattern="k239-result[.]rds$",recursive=TRUE,full.names=TRUE)
stopifnot(length(files)>0)
rows<-lapply(files,function(file) {
 r<-readRDS(file)$value;p<-readRDS(sub("-result.rds","-prepared.rds",file,fixed=TRUE))$value$p
 z<-r$candidates$stress/r$scale_unit
 state<-grip:::grip.edge.isometric.energy.gradient(z,p$edges,p$edge_targets,rep(1,nrow(p$edges)),scale=1,edge_length_epsilon=0)
 chord<-sqrt(rowSums((z[p$pair_matrix[,1],]-z[p$pair_matrix[,2],])^2))
 target<-p$distance_matrix[p$pair_matrix];raw<-sum((chord-target)^2)
 # Near-tie numerical geodesics may violate triangles by a few ulps/tolerances.
 discrepancy<-abs(2*state$energy-raw)
 stopifnot(discrepancy<1e-5)
 data.frame(case=r$scores$case[1],regime=r$scores$regime[1],twice_edge_energy=2*state$energy,raw_stress_normalized_units=raw,absolute_difference=discrepancy)
})
write.csv(do.call(rbind,rows),file.path(out,"objective-validation.csv"),row.names=FALSE)
cat("Identity scale, uniform stiffness:",length(rows),"complete-graph controls passed (2E = raw stress).\n")
