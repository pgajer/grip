source("papers/grip-software-paper/reproducibility/experiments/mds-initializer-sensitivity/common.R")
source(file.path(pilot_dir,"score-saddle-reference.R"))
args<-commandArgs(TRUE); reps<-if(length(args)) as.integer(strsplit(args[1],",")[[1]]) else 1:5
for(r in reps) {
 x<-readRDS(file.path(input_dir,sprintf("cloud-%02d.rds",r)))$coords
 tri<-geometry::delaunayn(x[,1:2]);ref<-saddle.reference.mesh(x,tri,.8,2L)
 files<-list.files(out,sprintf("^result-r%02d-k[0-9]+[.]rds$",r),full.names=TRUE)
 for(file in files) {
  f<-readRDS(file);dest<-sub("result-","surface-",file,fixed=TRUE)
  identity<-list(result_md5=unname(tools::md5sum(file)),script=unname(tools::md5sum(c(
    file.path(source_dir,"surface.R"),file.path(pilot_dir,"score-saddle-reference.R")))),
    protocol=protocol,geometry=as.character(packageVersion("geometry")))
  if(file.exists(dest)) {stopifnot(identical(readRDS(dest)$identity,identity));next}
  logmsg("Surface",r,f$k);rows<-list()
  candidates<-f$candidates
  if(f$selected) candidates<-c(candidates,lapply(f$extras,`[[`,"coords"))
  for(nm in names(candidates)) for(alignment in c("rigid","similarity")) {
   fit<-grip::score.coordinates(candidates[[nm]],x,alignment=alignment)
   s<-grip::score.surface(fit$coords,tri,ref$coords,ref$triangles,sample_size=8000L,seed=1901L+r)
   rows[[length(rows)+1L]]<-data.frame(replicate=r,k=f$k,method=nm,alignment=alignment,
    sample_size=8000L,seed=1901L+r,subdivisions=2L,coordinate_relative_rmse=fit$relative_rmse,
    coordinate_rmse=fit$rmse,alignment_scale=fit$scale,surface_rms=s$rms,surface_mc_se=s$rms_mc_se,
    surface_forward_rms=s$forward_rms,surface_reverse_rms=s$reverse_rms,area=s$area,
    reference_area=s$reference_area,zero_area_faces=s$zero_area_faces)
  }
  atomic_save(list(identity=identity,scores=do.call(rbind,rows)),dest)
 }
}
