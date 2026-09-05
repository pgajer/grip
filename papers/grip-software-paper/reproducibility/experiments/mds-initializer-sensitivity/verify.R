source("papers/grip-software-paper/reproducibility/experiments/mds-initializer-sensitivity/common.R")
sd<-file.path(source_dir,"summary");read<-function(n) read.csv(file.path(sd,paste0(n,".csv")))
scores<-read("scores");starts<-read("starts");selection<-read("graph-selection");status<-read("edge-optimizer-status")
surface<-read("surface-scores");extra<-read("additional-budget-scores")
stopifnot(nrow(scores)==125L,nrow(starts)==150L,nrow(selection)==25L,nrow(extra)==15L,nrow(surface)==280L,
  !anyDuplicated(scores[,c("replicate","k","method")]),all(table(starts$replicate,starts$k)[table(starts$replicate,starts$k)>0]==6L),
  all(starts$raw_stress<=starts$initial_raw_stress+1e-6),all(is.finite(starts$raw_stress)),
  all(starts$termination %in% c("stress_tolerance","iteration_limit")),
  all(status$max_energy_increase<1e-7),all(is.finite(surface$surface_rms)),
  all(surface$zero_area_faces==0),all(scores$package_score_difference<1e-10),
  all(scores$independent_path_difference<1e-9))
controls<-subset(scores,method=="Original saddle")
stopifnot(max(controls$path_rel)<1e-7,max(controls$edge_rel)<1e-12,
  max(abs(controls$xz_path_error-controls$xg_error))<1e-7)
mds<-subset(scores,method=="Stress MDS")
stopifnot(max(abs(mds$raw_optimal_scale-1))<1e-7,
  max(abs(mds$raw_target_rmse-mds$stress1))<1e-7)
for(r in 1:5) {
 sel<-subset(selection,replicate==r);k0<-sel$k[sel$selected]
 stopifnot(identical(sort(sel$k),sort(c(32L,k0+c(-5L,0L,5L),80L))))
 for(kk in sel$k) {
  s<-subset(scores,replicate==r & scores$k==kk)
  a<-subset(starts,replicate==r & starts$k==kk)
  stopifnot(sum(a$selected)==1L,abs(a$raw_stress[a$selected]-min(a$raw_stress))<1e-8,
    abs(s$raw_stress[s$method=="Stress MDS"]-min(a$raw_stress))<1e-6,
    s$stress1[s$method=="Stress MDS"]<=s$stress1[s$method=="Classical MDS"]+1e-8)
 }
}
# Recompute every candidate against its saved routes and unchanged references.
Rcpp::sourceCpp(file.path(pilot_dir,"path-lengths.cpp"),cacheDir=file.path(out,"cpp-verify"))
max_recomputed<-0;max_raw_identity<-0
for(file in list.files(out,"^result-r[0-9]+-k[0-9]+[.]rds$",full.names=TRUE)) {
 f<-readRDS(file);stopifnot(identical(f$identity$protocol,protocol),identical(f$identity$source,source_id()))
 r<-f$replicate;k<-f$k;logmsg("Verify",r,k)
 p<-readRDS(file.path(out,sprintf("r%02d-k%02d-prepared.rds",r,k)))$value$p
 cal<-readRDS(file.path(input_dir,sprintf("calibration-%02d.rds",r)))
 x<-f$candidates[["Original saddle"]]
 recomputed<-do.call(rbind,lapply(names(f$candidates),function(nm) score_layout(f$candidates[[nm]],nm,r,k,p,x,cal$reference)))
 cols<-names(f$scores)[vapply(f$scores,is.numeric,logical(1))]
 max_recomputed<-max(max_recomputed,max(abs(as.matrix(recomputed[,cols])-as.matrix(f$scores[,cols]))))
 # Independent scale-profile identity on full target matrix, avoiding sqrt cancellation.
 delta<-as.vector(as.dist(p$distance_matrix));C<-sum(delta^2)
 for(z in f$candidates) {
  d<-as.vector(dist(z));A<-sum(d*delta);B<-sum(d^2)
  sraw<-sum((A/B*d-delta)^2)/C;sp<-sum((d-A/C*delta)^2)/B
  max_raw_identity<-max(max_raw_identity,abs(sraw-sp))
 }
}
stopifnot(max_recomputed<1e-9,max_raw_identity<1e-12)
for(name in c("input-checksums","result-checksums")) {
 h<-read(name);stopifnot(identical(unname(tools::md5sum(h$path)),h$md5))
}
parity<-read("baseline-parity");stopifnot(max(parity$baseline_max_difference,na.rm=TRUE)<1e-7)
# Ensure surface coordinate diagnostics use the same alignment convention.
sim<-subset(surface,alignment=="similarity" & method %in% methods)
pair<-merge(scores,sim,by=c("replicate","k","method"))
stopifnot(nrow(pair)==125L,max(abs(pair$procrustes-pair$coordinate_relative_rmse))<1e-10)
validation<-list(status="passed",graphs=25L,primary_candidates=125L,smacof_starts=150L,
 surface_rows=280L,additional_fits=15L,max_score_recomputation_error=max_recomputed,
 max_scale_identity_error=max_raw_identity,max_package_score_difference=max(scores$package_score_difference),
 max_independent_path_difference=max(scores$independent_path_difference),
 max_baseline_score_difference=max(parity$baseline_max_difference,na.rm=TRUE),
 original_inputs_unchanged=TRUE,verified_at=format(Sys.time(),tz="America/New_York",usetz=TRUE))
jsonlite::write_json(validation,file.path(sd,"validation.json"),pretty=TRUE,auto_unbox=TRUE)
print(validation)
