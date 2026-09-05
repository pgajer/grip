source("papers/grip-software-paper/reproducibility/experiments/mds-initializer-sensitivity/common.R")
r<-as.integer(commandArgs(TRUE)[1]);stopifnot(r %in% 1:5)
mode<-if(length(commandArgs(TRUE))>1) commandArgs(TRUE)[2] else "selected"
freeze<-readRDS(file.path(out,"freeze.rds"));stopifnot(identical(freeze$protocol,protocol))
Rcpp::sourceCpp(file.path(pilot_dir,"path-lengths.cpp"),cacheDir=file.path(out,paste0("cpp-",r)))
selection<-subset(freeze$selection,replicate==r)
selection<-selection[order(!selection$selected,selection$k),]
if(mode=="selected") selection<-subset(selection,selected)
if(mode=="sensitivity") selection<-subset(selection,!selected)
x<-readRDS(file.path(input_dir,sprintf("cloud-%02d.rds",r)))$coords
cal<-readRDS(file.path(input_dir,sprintf("calibration-%02d.rds",r)))
input_hash<-subset(freeze$inputs,grepl(sprintf("-%02d",r),path))
identity<-list(protocol=protocol,source=source_id(),inputs=input_hash)
cache <- function(label,expr) {
 path<-file.path(out,sprintf("r%02d-k%02d-%s.rds",r,k,label))
 if(file.exists(path)) {a<-readRDS(path);stopifnot(identical(a$identity,identity));return(a$value)}
 logmsg("BEGIN",r,k,label);tick<-proc.time()[[3]];a<-force(expr)
 a$elapsed<-proc.time()[[3]]-tick;atomic_save(list(identity=identity,value=a),path)
 logmsg("DONE",r,k,label,sprintf("%.2f seconds",a$elapsed));a
}
allg<-NULL
if(any(!selection$selected)) {logmsg("Loading graphs",r);allg<-readRDS(file.path(input_dir,sprintf("graphs-%02d.rds",r)))$graphs}
for(k in selection$k) {
 selected<-k==cal$selected_k
 final_path<-file.path(out,sprintf("result-r%02d-k%02d.rds",r,k))
 if(file.exists(final_path)) {stopifnot(identical(readRDS(final_path)$identity,identity));next}
 old<-if(selected) readRDS(file.path(input_dir,sprintf("fit-%02d-k%d.rds",r,k))) else NULL
 g<-if(selected) old$graph else allg[[as.character(k)]]
 prepared<-cache("prepared",{
  if(selected) list(p=restore(old$prepared,g)) else {
   p<-grip::prepare.geodesic.kk(g$edges,n=nrow(x),edge_weights=g$weights,tie_mode="single")
   p$path_edges<-p$path_vertices<-p$path_edge_weights<-NULL
   list(p=restore(p,g))
  }
 })
 p<-prepared$p
 stopifnot(nrow(p$pair_matrix)==499500L,all(p$flat_edge_coeff==1),
   max(abs(p$pair_graph_distance-g$distances[p$pair_matrix]))<1e-7)
 cl<-cache("classical",slim(grip::classical.mds(prepared=p,dim=3,diagnostics=FALSE)))
 if(selected) stopifnot(max(abs(as.vector(dist(cl$coords))-as.vector(dist(old$candidates[[2]]))))<1e-7)
 starts<-lapply(1:protocol$starts,function(j) cache(paste0("smacof-start",j),slim(grip::metric.mds(
   prepared=p,dim=3,init=if(j==1) cl$coords else "random",n_init=1L,
   max_iter=protocol$max_iter,eps=protocol$eps,seed=7300000L+10000L*r+j,diagnostics=FALSE))))
 losses<-vapply(starts,function(s) s$metadata$raw_stress,numeric(1));best<-which.min(losses);mds<-starts[[best]]
 start_rows<-do.call(rbind,lapply(seq_along(starts),function(j) {
   a<-starts[[j]]$metadata$starts;a$start<-j;a$initialization<-if(j==1) "classical" else "random"
   a$replicate<-r;a$k<-k;a$selected<-j==best;a$elapsed<-starts[[j]]$elapsed;a
 }))
 runkk<-function(z,extra=FALSE) slim(grip::edge.kk(coords=z,prepared=p,dim=3,
   max_iter=if(extra) protocol$extra_edge_steps else protocol$edge_max_iter,
   stiffness_method="density",density_mix_schedule=if(extra) 1 else protocol$mix,
   scale_mode="profiled",edge_length_epsilon=0,diagnostics=FALSE,return_trace=TRUE,
   seed=if(extra) 6211000L+r else 5211000L+r))
 # Rerun selected classical refinement once to check package parity; use matched timings everywhere.
 ck<-cache("classical-kk",runkk(cl$coords));sk<-cache("stress-kk",runkk(mds$coords))
 if(selected) stopifnot(max(abs(as.vector(dist(ck$coords))-as.vector(dist(old$candidates[[3]]))))<1e-7)
 candidates<-setNames(list(x,cl$coords,mds$coords,ck$coords,sk$coords),methods)
 scores<-do.call(rbind,lapply(names(candidates),function(nm) score_layout(candidates[[nm]],nm,r,k,p,x,cal$reference)))
 stopifnot(scores$edge_rel[1]<1e-12,scores$path_rel[1]<1e-7,
  abs(scores$xz_path_error[1]-scores$xg_error[1])<1e-7,
  abs(scores$raw_stress[3]-min(losses))<1e-6,
  scores$stress1[3]<=scores$stress1[2]+1e-8)
 baseline_error<-NA_real_
 if(selected) {
  cols<-c("path_rel","edge_rel","stress1","xz_path_error","procrustes")
  baseline_error<-max(abs(as.matrix(scores[c(1,2,4),cols])-as.matrix(old$scores[,cols])))
  stopifnot(baseline_error<1e-7)
 }
 extras<-list();extra_scores<-NULL
 if(selected) {
  extras$stress_mds<-cache("extra-smacof",slim(grip::metric.mds(prepared=p,dim=3,init=mds$coords,
    max_iter=protocol$extra_mds_steps,eps=protocol$extra_mds_eps,diagnostics=FALSE)))
  extras$classical_kk<-cache("extra-classical-kk",runkk(ck$coords,TRUE))
  extras$stress_kk<-cache("extra-stress-kk",runkk(sk$coords,TRUE))
  extra_scores<-do.call(rbind,lapply(names(extras),function(nm) score_layout(extras[[nm]]$coords,nm,r,k,p,x,cal$reference)))
 }
 deg<-tabulate(as.vector(g$edges),nbins=nrow(x));q<-function(v) as.numeric(quantile(v,c(0,.05,.5,.95,1)))
 graph_stats<-data.frame(replicate=r,k=k,selected=selected,edges=nrow(g$edges),bridges=g$bridges,
  components_before=g$components_before,t(c(setNames(q(deg),paste0("degree_",c("min","q05","median","q95","max"))),
    setNames(q(g$weights),paste0("length_",c("min","q05","median","q95","max"))))),check.names=FALSE)
 result<-list(identity=identity,replicate=r,k=k,selected=selected,candidates=candidates,
  scores=scores,starts=start_rows,selected_start=best,extras=extras,extra_scores=extra_scores,
  traces=list(classical=ck$trace,stress=sk$trace),
  stage_summaries=list(classical=ck$metadata$stage_summaries,stress=sk$metadata$stage_summaries),
  graph_stats=graph_stats,baseline_max_difference=baseline_error,
  timing=c(preparation=prepared$elapsed,classical=cl$elapsed,smacof=sum(start_rows$elapsed),
   classical_kk=ck$elapsed,stress_kk=sk$elapsed),session=sessionInfo())
 atomic_save(result,final_path)
 logmsg("RESULT",r,k,"best start",best,"path errors",paste(signif(scores$path_rel,4),collapse=","))
 rm(p,prepared,old,g,starts,ck,sk,result,extras);gc(FALSE)
}
