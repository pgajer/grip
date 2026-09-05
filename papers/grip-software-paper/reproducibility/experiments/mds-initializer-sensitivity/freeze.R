source("papers/grip-software-paper/reproducibility/experiments/mds-initializer-sensitivity/common.R")
rows<-list(); manifests<-list()
for(r in 1:5) {
  path<-file.path(input_dir,sprintf("calibration-%02d.rds",r));cal<-readRDS(path)
  ks<-sort(unique(c(32,cal$selected_k+c(-5,0,5),80)))
  q<-cal$curve[match(ks,cal$curve$k),];q$selected<-q$k==cal$selected_k
  q$role<-ifelse(q$selected,"selected",ifelse(q$k==32,"sparse",ifelse(q$k==80,"dense","nearby")))
  q$within_one_percent<-q$xg_error<=1.01*min(cal$curve$xg_error)
  q$plateau_min<-cal$checks$within_one_percent_k_min;q$plateau_max<-cal$checks$within_one_percent_k_max
  q$reference_fine_difference<-cal$checks$finer_relative_difference
  rows[[r]]<-q
  files<-file.path(input_dir,c(sprintf("cloud-%02d.rds",r),sprintf("graphs-%02d.rds",r),
    sprintf("calibration-%02d.rds",r),sprintf("fit-%02d-k%d.rds",r,cal$selected_k)))
  manifests[[r]]<-data.frame(path=files,md5=unname(tools::md5sum(files)))
}
selection<-do.call(rbind,rows);manifest<-do.call(rbind,manifests)
if(file.exists(file.path(out,"freeze.rds"))) {
 old<-readRDS(file.path(out,"freeze.rds"));stopifnot(identical(old$protocol,protocol),identical(old$selection,selection),identical(old$inputs,manifest))
} else atomic_save(list(protocol=protocol,selection=selection,inputs=manifest,
  frozen_at=format(Sys.time(),tz="America/New_York",usetz=TRUE),session=sessionInfo(),
  git=system("git rev-parse HEAD",intern=TRUE)),file.path(out,"freeze.rds"))
write.csv(selection,file.path(out,"graph-selection.csv"),row.names=FALSE)
write.csv(manifest,file.path(out,"input-checksums.csv"),row.names=FALSE)
logmsg("Frozen",nrow(selection),"graphs before new fits")
