out<-"output/mds-edge-kk-radius";src<-"tools/experiments/mds-edge-kk-radius"
files<-list.files(file.path(out,"fits"),pattern="(primary|uniform)[.]rds$",recursive=TRUE,full.names=TRUE)
rows<-lapply(files,function(file) {
 f<-readRDS(file)$value;name<-basename(file);prefix<-sub("-(classical|stress|full_classical)_(primary|uniform)[.]rds$","",name)
 p<-readRDS(file.path(dirname(file),paste0(prefix,"-prepared.rds")))$value$p;z<-f$coords
 l<-sqrt(rowSums((z[p$edges[,1],]-z[p$edges[,2],])^2));w<-p$edge_targets;b<-sum(l*w)/sum(w*w)
 e<-sum((l-b*w)^2)/2;eh<-sum((l/2-b*w/2)^2)/2
 rel<-sqrt(sum((l-b*w)^2)/sum((b*w)^2))
 stopifnot(abs(eh-e/4)<1e-12)
 data.frame(case=basename(dirname(file)),fit=name,edge_scale=b,profiled_uniform_energy=e,energy_after_halving=eh,edge_relative_rmse=rel,
  final_gradient=tail(f$trace$gradient_norm,1),elapsed=f$elapsed)
})
dir.create(file.path(src,"summary"),showWarnings=FALSE)
write.csv(do.call(rbind,rows),file.path(src,"summary/initial-pilot-scale-audit.csv"),row.names=FALSE)
cat(length(rows),"archived pilot fits; minimum edge scale",min(vapply(rows,`[[`,numeric(1),"edge_scale")),"\n")
