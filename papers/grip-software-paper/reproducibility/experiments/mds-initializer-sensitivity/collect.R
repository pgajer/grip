source("papers/grip-software-paper/reproducibility/experiments/mds-initializer-sensitivity/common.R")
summary_dir<-file.path(source_dir,"summary");dir.create(summary_dir,showWarnings=FALSE)
files<-list.files(out,"^result-r[0-9]+-k[0-9]+[.]rds$",full.names=TRUE)
stopifnot(length(files)==25L)
fits<-lapply(files,readRDS);freeze<-readRDS(file.path(out,"freeze.rds"))
rows<-function(field) do.call(rbind,lapply(fits,`[[`,field))
write_table<-function(x,name) write.csv(x,file.path(summary_dir,paste0(name,".csv")),row.names=FALSE,na="")
scores<-rows("scores");scores<-merge(scores,freeze$selection[,c("replicate","k","selected","role")],by=c("replicate","k"),sort=FALSE)
write_table(scores,"scores");write_table(rows("starts"),"starts");write_table(rows("extra_scores"),"additional-budget-scores")
write_table(rows("graph_stats"),"graph-statistics");write_table(freeze$selection,"graph-selection")
write_table(freeze$inputs,"input-checksums")
surface_files<-sub("result-","surface-",files,fixed=TRUE);stopifnot(all(file.exists(surface_files)))
write_table(do.call(rbind,lapply(surface_files,function(p) readRDS(p)$scores)),"surface-scores")
status<-list();timing<-list();coords<-list();extra_starts<-list()
trace_rows<-function(trace,f,method,budget,run) {
 do.call(rbind,lapply(split(trace,trace$stage),function(t) data.frame(replicate=f$replicate,k=f$k,
  method=method,run=run,stage=t$stage[1],mix=t$mix[1],iterations=max(t$iteration),budget=budget,
  hit_budget=max(t$iteration)>=budget,final_gradient=tail(t$gradient_norm,1),
  final_energy=tail(t$energy,1),final_edge_scale=tail(t$edge.scale,1),
  final_edge_rel=tail(t$edge.rel.rmse,1),accepted_steps=sum(t$accepted),
  max_energy_increase=max(c(0,diff(t$energy))))))
}
for(f in fits) {
 for(nm in names(f$traces)) status[[length(status)+1L]]<-trace_rows(f$traces[[nm]],f,nm,200L,"primary")
 if(f$selected) {
  for(nm in c("classical_kk","stress_kk")) status[[length(status)+1L]]<-trace_rows(f$extras[[nm]]$trace,f,nm,1000L,"additional")
  a<-f$extras$stress_mds$metadata$starts;a$replicate<-f$replicate;a$k<-f$k
  extra_starts[[length(extra_starts)+1L]]<-a
 }
 timing[[length(timing)+1L]]<-data.frame(replicate=f$replicate,k=f$k,
  component=names(f$timing),elapsed=unname(f$timing))
 if(f$selected) timing[[length(timing)+1L]]<-data.frame(replicate=f$replicate,k=f$k,
  component=paste0("extra_",names(f$extras)),elapsed=vapply(f$extras,`[[`,numeric(1),"elapsed"))
 # All layouts are small; retain coordinates for independent score and figure checks.
 for(nm in names(f$candidates)) {
  z<-f$candidates[[nm]]
  coords[[length(coords)+1L]]<-data.frame(replicate=f$replicate,k=f$k,method=nm,vertex=1:nrow(z),x=z[,1],y=z[,2],z=z[,3])
 }
}
write_table(do.call(rbind,status),"edge-optimizer-status");write_table(do.call(rbind,timing),"timings")
write_table(do.call(rbind,extra_starts),"additional-mds-starts")
con<-gzfile(file.path(summary_dir,"coordinates.csv.gz"),"wt");write.csv(do.call(rbind,coords),con,row.names=FALSE);close(con)
write_table(data.frame(replicate=vapply(fits,`[[`,integer(1),"replicate"),k=vapply(fits,`[[`,integer(1),"k"),
 selected=vapply(fits,`[[`,logical(1),"selected"),
 baseline_max_difference=vapply(fits,`[[`,numeric(1),"baseline_max_difference")),"baseline-parity")
write_table(data.frame(path=c(files,surface_files),md5=unname(tools::md5sum(c(files,surface_files)))),"result-checksums")
jsonlite::write_json(list(protocol=protocol,frozen_at=freeze$frozen_at,baseline_commit=freeze$git,
  source_md5=as.list(source_id()),R=R.version.string,platform=R.version$platform,
  geometry=as.character(packageVersion("geometry")),Rcpp=as.character(packageVersion("Rcpp")),
  BLAS=extSoftVersion()["BLAS"],fit_count=length(fits),candidate_rows=nrow(scores)),
  file.path(summary_dir,"manifest.json"),pretty=TRUE,auto_unbox=TRUE)
writeLines(trimws(capture.output(sessionInfo()),which="right"),file.path(summary_dir,"session-info.txt"))
logmsg("Collected",length(fits),"graphs and",nrow(scores),"primary scores")
