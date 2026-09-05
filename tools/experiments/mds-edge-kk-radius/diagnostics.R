# Additional locality and timing diagnostics from frozen graphs; does not refit.
src<-"tools/experiments/mds-edge-kk-radius";out<-"output/mds-edge-kk-radius"
files<-list.files(file.path(out,"fits-v2"),pattern="-result[.]rds$",recursive=TRUE,full.names=TRUE);files<-files[!grepl("-extra-result",files,fixed=TRUE)]
q<-function(v) as.numeric(quantile(v,c(0,.05,.5,.95,1)))
rows<-lapply(files,function(file) {
 v<-readRDS(file)$value;id<-v$scores[1,c("case","surface","sampling","replicate","n","radius","regime","k")]
 p<-readRDS(sub("-result.rds","-prepared.rds",file,fixed=TRUE))$value$p
 x<-as.matrix(read.csv(file.path(out,"inputs",paste0(id$case,"-truth.csv")),header=FALSE))
 theta<-atan2(x[,2],x[,1]);rho<-sqrt(rowSums(x[,1:2]^2))/id$radius
 angle<-abs(theta[p$edges[,1]]-theta[p$edges[,2]]);angle<-pmin(angle,2*pi-angle)
 radial<-abs(rho[p$edges[,1]]-rho[p$edges[,2]])
 cbind(id,t(c(setNames(q(angle),paste0("angular_gap_",c("min","q05","median","q95","max"))),
  setNames(q(radial),paste0("radial_gap_",c("min","q05","median","q95","max"))))))
})
write.csv(do.call(rbind,rows),file.path(src,"summary/locality.csv"),row.names=FALSE)
timing<-do.call(rbind,lapply(files,function(file) {v<-readRDS(file)$value;cbind(v$scores[1,c("case","regime","k")],t(v$timing))}))
write.csv(timing,file.path(src,"summary/timings.csv"),row.names=FALSE)
# Deterministic representative: first sample, base-disk measure, k32 and complete.
exports<-list()
for(file in files) {
 v<-readRDS(file)$value;id<-v$scores[1,]
 if(id$replicate!=1 || id$n!=240 || id$sampling!="disk" || !(id$radius %in% c(1,8,64)) || !(id$k %in% c(32,239))) next
 for(nm in c("original","classical","stress","stress_primary","stress_fixed_primary")) {
  z<-v$candidates[[nm]];score<-v$scores[v$scores$method==nm,]
  exports[[length(exports)+1L]]<-cbind(score[rep(1,nrow(z)),c("case","surface","radius","regime","k","method","edge_scale")],vertex=seq_len(nrow(z)),x=z[,1],y=z[,2],z=z[,3],x_r2=z[,1]/id$radius^2,y_r2=z[,2]/id$radius^2,z_r2=z[,3]/id$radius^2)
 }
}
write.csv(do.call(rbind,exports),file.path(src,"summary/snapshot-coordinates.csv"),row.names=FALSE)
cat("Locality, timing, and snapshot exports complete:",length(files),"graphs\n")
