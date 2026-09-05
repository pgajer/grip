#!/usr/bin/env Rscript
lib<-Sys.getenv("GRIP_RJOURNAL_PACKAGE_LIBRARY");if(nzchar(lib)) .libPaths(c(lib,.libPaths()))
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)>=1)
out<-normalizePath(args[1]); reps<-if(length(args)>1) as.integer(strsplit(args[2],",")[[1]]) else 1:5
script<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
Rcpp::sourceCpp(file.path(dirname(normalizePath(script)),"path-lengths.cpp"),cacheDir=file.path(out,"cpp-cache"))
lengths_of<-function(z,e) sqrt(rowSums((z[e[,1],,drop=FALSE]-z[e[,2],,drop=FALSE])^2))
relative<-function(y,d) sqrt(sum((y-d)^2)/sum(d*d))
profile<-function(y,d) {a<-sum(y*d)/sum(d*d);relative(y,a*d)}
align<-function(z,x) {
  z<-sweep(z,2,colMeans(z),"-");xx<-sweep(x,2,colMeans(x),"-")
  f<-svd(crossprod(z,xx));sweep(sum(f$d)/sum(z*z)*z%*%f$u%*%t(f$v),2,colMeans(x),"+")
}
for(r in reps) {
  cal<-readRDS(file.path(out,sprintf("calibration-%02d.rds",r)))
  if(length(args)>2L) cal$selected_k<-as.integer(args[3])
  else if(cal$checks$search_boundary) stop("Calibration optimum is at the search boundary")
  dest<-file.path(out,sprintf("fit-%02d-k%d.rds",r,cal$selected_k))
  if(file.exists(dest)) {cat("Existing",dest,"\n");next}
  started<-proc.time()[[3]];x<-readRDS(file.path(out,sprintf("cloud-%02d.rds",r)))$coords
  allg<-readRDS(file.path(out,sprintf("graphs-%02d.rds",r)))
  g<-allg$graphs[[as.character(cal$selected_k)]];rm(allg);gc(FALSE)
  tick<-proc.time()[[3]]
  p<-grip::prepare.geodesic.kk(g$edges,n=1000,edge_weights=g$weights,tie_mode="single")
  prepare_seconds<-proc.time()[[3]]-tick
  dg<-p$pair_graph_distance;strict<-g$distances[p$pair_matrix]
  stopifnot(nrow(p$pair_matrix)==499500L,all(p$flat_edge_coeff==1),
    max(abs(dg-strict)/pmax(1,dg,strict))<=sqrt(.Machine$double.eps))
  key<-function(e) e[,1]*1000+e[,2]
  ref_index<-match(key(cal$reference$pairs),key(p$pair_matrix));stopifnot(!anyNA(ref_index))
  dx<-cal$reference$d
  tick<-proc.time()[[3]]
  mds<-grip::classical.mds(prepared=p,dim=3,diagnostics=FALSE)
  mds_seconds<-proc.time()[[3]]-tick
  tick<-proc.time()[[3]]
  kk<-grip::edge.kk(coords=mds$coords,prepared=p,dim=3,max_iter=200L,
    stiffness_method="density",density_mix_schedule=c(0,.25,.5,.75,1),
    scale_mode="profiled",edge_length_epsilon=0,diagnostics=FALSE,return_trace=TRUE,seed=5211000+r)
  kk_seconds<-proc.time()[[3]]-tick
  candidates<-list("Original saddle"=x,"Metric MDS"=mds$coords,"MDS + edge-KK"=kk$coords)
  score<-function(z,name) {
    stopifnot(ncol(z)==3L,all(is.finite(z)))
    edge<-lengths_of(z,p$edges);chord<-lengths_of(z,p$pair_matrix)
    path<-pilot_path_lengths(z,p$flat_pair_edge_offsets,p$flat_edge_u,p$flat_edge_v)
    a<-sum(chord*dg)/sum(dg*dg);edge_scale<-sum(edge*p$edge_targets)/sum(p$edge_targets^2)
    check<-grip::score.gmds(z,prepared=p,edge_length_epsilon=0)
    err<-max(abs(c(profile(path,dg)-check$gmds.stress,profile(edge,p$edge_targets)-check$edge.rel.rmse)))
    stopifnot(err<1e-10)
    # Scalar conversion determined only from graph edges, never from D_X.
    physical_path<-path/edge_scale
    independent_ids<-as.integer(seq.int(1L,length(path),length.out=101L))
    independent<-vapply(independent_ids,function(i) sum(lengths_of(z,p$path_edges[[i]])),numeric(1))
    independent_error<-max(abs(path[independent_ids]-independent))
    stopifnot(independent_error<1e-9)
    data.frame(replicate=r,k=g$k,method=name,path_rel=profile(path,dg),
      edge_rel=profile(edge,p$edge_targets),stress1=sqrt(sum((chord-a*dg)^2)/sum(chord^2)),
      xg_error=relative(strict[ref_index],dx),
      xz_path_error=relative(physical_path[ref_index],dx),
      graph_path_error_edge_scale=relative(physical_path[ref_index],strict[ref_index]),
      procrustes=sqrt(sum((align(z,x)-x)^2)/sum(sweep(x,2,colMeans(x),"-")^2)),
      edge_scale=edge_scale,package_score_difference=err,independent_path_difference=independent_error)
  }
  scores<-do.call(rbind,lapply(names(candidates),function(name) score(candidates[[name]],name)))
  stopifnot(scores$edge_rel[1]<1e-12,scores$path_rel[1]<1e-7,
    abs(scores$xz_path_error[1]-scores$xg_error[1])<1e-7)
  # Keep an additional-budget sensitivity run separate from the primary scores.
  tick<-proc.time()[[3]]
  audit<-grip::edge.kk(coords=kk$coords,prepared=p,dim=3,max_iter=1000L,
    stiffness_method="density",density_mix_schedule=1,scale_mode="profiled",edge_length_epsilon=0,
    diagnostics=FALSE,return_trace=TRUE,seed=6211000+r)
  audit_seconds<-proc.time()[[3]]-tick;audit_score<-score(audit$coords,"Additional 1000 steps")
  cat(sprintf("cloud=%d k=%d MDS path=%.5f edge-KK path=%.5f audit=%.5f; MDS %.2fs KK %.2fs\n",r,g$k,
    scores$path_rel[2],scores$path_rel[3],audit_score$path_rel,mds_seconds,kk_seconds));flush.console()
  # Store flat routes and pairs, not the duplicated list representation.
  saved_p<-p[c("n","edges","edge_targets","pair_matrix","pair_graph_distance",
    "flat_pair_edge_offsets","flat_edge_u","flat_edge_v","flat_edge_coeff")]
  result<-list(replicate=r,k=g$k,coords=x,graph=g,prepared=saved_p,reference=cal$reference,
    candidates=candidates,scores=scores,edge_trace=kk$trace,edge_metadata=kk$metadata,
    audit=list(coords=audit$coords,scores=audit_score,trace=audit$trace,metadata=audit$metadata),
    timing=c(preparation=prepare_seconds,mds=mds_seconds,edge_kk=kk_seconds,audit=audit_seconds,
      total=proc.time()[[3]]-started),
    strict_distance_max_difference=max(abs(dg-strict)),session=sessionInfo(),
    protocol=list(dim=3L,edge_max_iter_per_stage=200L,mix=c(0,.25,.5,.75,1),
      audit_max_iter=1000L,scoring_epsilon=0,reference_mesh=81L,reference_sources=128L))
  saveRDS(result,dest)
  rm(result,p,saved_p,kk,audit,g);gc(FALSE)
}
