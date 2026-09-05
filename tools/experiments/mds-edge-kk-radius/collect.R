src<-"tools/experiments/mds-edge-kk-radius";out<-"output/mds-edge-kk-radius"
mode<-if(length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "partial"
files<-list.files(file.path(out,"fits-v2"),pattern="-result[.]rds$",recursive=TRUE,full.names=TRUE)
results<-lapply(files,readRDS);vals<-lapply(results,`[[`,"value")
summary<-if(mode=="final") file.path(src,"summary") else file.path(out,"summary")
dir.create(summary,recursive=TRUE,showWarnings=FALSE)
write<-function(df,name) write.csv(df,file.path(summary,paste0(name,".csv")),row.names=FALSE)
isextra<-grepl("-extra-result.rds",files,fixed=TRUE)
write(do.call(rbind,lapply(vals[!isextra],`[[`,"scores")),"scores")
if(any(isextra)) write(do.call(rbind,lapply(vals[isextra],`[[`,"scores")),"optimizer-sensitivity")
write(do.call(rbind,lapply(vals[!isextra],`[[`,"graph_stats")),"graphs")
starts<-do.call(rbind,lapply(vals[!isextra],function(v) {a<-v$starts;a$selected<-a$start==v$selected_start;a$initialization_description<-ifelse(a$start==1,"classical","Gaussian");a$raw_stress<-a$raw_stress*v$scale_unit^2;a$initial_raw_stress<-a$initial_raw_stress*v$scale_unit^2;id<-v$scores[1,c("case","surface","sampling","replicate","n","radius","regime","k")];cbind(id[rep(1,nrow(a)),],a)}));write(starts,"starts")
status<-do.call(rbind,lapply(vals,function(v) do.call(rbind,lapply(names(v$fits),function(nm) {
 f<-v$fits[[nm]];t<-f$trace;id<-v$scores[1,c("case","surface","sampling","replicate","n","radius","regime","k")]
 do.call(rbind,lapply(split(t,t$stage),function(a) cbind(id,method=nm,stage=a$stage[1],mix=a$mix[1],iterations=max(a$iteration),last_energy=tail(a$energy,1),gradient_norm=tail(a$gradient_norm,1),accepted=sum(a$accepted[a$iteration>0]),elapsed=f$elapsed)))
}))));write(status,"optimizer-status")
if(mode=="final") {
 stopifnot(sum(!isextra)==1216,sum(isextra)==32,nrow(starts)==1216*3)
 for(s in c("paraboloid","saddle")) for(m in c("disk","surface_area")) for(rp in 1:3) for(r in c(1,2,4,8,16,32,64)) for(g in c("ambient","geodesic")) {
  key<-sprintf("%s-%s-rep%d-n240-r%d",s,m,rp,r)
  stopifnot(sum(vapply(vals[!isextra],function(v) v$scores$case[1]==key && v$scores$regime[1]==g,logical(1)))==7)
 }
 saveRDS(list(files=files,identities=lapply(results,`[[`,"identity")),file.path(summary,"fitting-provenance.rds"),compress=TRUE)
 compact<-lapply(vals,function(v) v[c("scores","candidates","scale_unit","graph_stats","starts","selected_start","timing")])
 dir.create(file.path(summary,"coordinates"),showWarnings=FALSE)
 ids<-vapply(compact,function(v) v$scores$case[1],character(1))
 for(id in unique(ids)) saveRDS(compact[ids==id],file.path(summary,"coordinates",paste0(id,".rds")),compress=TRUE)
}
cat("Collected",sum(!isextra),"primary graphs;",sum(isextra),"optimizer-control graphs\n")
