# A single case worker. Run from repository root; all checkpoint identities checked.
src <- "tools/experiments/mds-edge-kk-radius"
out <- "output/mds-edge-kk-radius"
fit_root <- file.path(out,"fits-v2")
args <- commandArgs(TRUE); stopifnot(length(args)==2)
key <- args[1]; mode <- args[2]; stopifnot(mode %in% c("pilot","main","extra"))
a <- strsplit(key,"-",fixed=TRUE)[[1]]
surface<-a[1];sampling<-a[2];replicate<-as.integer(sub("rep","",a[3]));n<-as.integer(sub("n","",a[4]));radius<-as.numeric(sub("r","",a[5]))
stopifnot(as.character(packageVersion("grip"))=="0.2.0.9001",as.character(packageVersion("smacof"))=="2.1.7")
logmsg <- function(...) {cat(format(Sys.time(),"%H:%M:%S"),key,...,"\n");flush.console()}
files <- c(file.path(src,c("run.R","path-lengths.cpp","PROTOCOL.md")),"R/metric_mds.R","R/grip_quality.R","R/gmds_layout_interface.R",list.files("src",pattern="[.]cpp$",full.names=TRUE))
input <- file.path(out,"inputs",paste0(key,c("-distance.csv","-truth.csv")))
identity<-list(source=tools::md5sum(files),input=tools::md5sum(input),grip=as.character(packageVersion("grip")),smacof=as.character(packageVersion("smacof")))
d<-as.matrix(read.csv(input[1],header=FALSE));x<-as.matrix(read.csv(input[2],header=FALSE))
scale_unit<-sqrt(mean(d[upper.tri(d)]^2));d<-d/scale_unit;x<-x/scale_unit
Rcpp::sourceCpp(file.path(src,"path-lengths.cpp"),cacheDir=file.path(out,"cpp",key))
dir.create(file.path(fit_root,key),recursive=TRUE,showWarnings=FALSE)
cache <- function(label,expr) {
 path<-file.path(fit_root,key,paste0(label,".rds"))
 if(file.exists(path)) {a<-readRDS(path);stopifnot(identical(a$identity,identity));return(a$value)}
 tick<-proc.time()[[3]];a<-force(expr);a$elapsed<-proc.time()[[3]]-tick
 tmp<-paste0(path,".tmp");saveRDS(list(identity=identity,value=a),tmp,compress=TRUE);stopifnot(file.rename(tmp,path));a
}
slim<-function(f) {f$prepared<-NULL;f$metadata$frames<-NULL;f}
lengths_of<-function(z,e) sqrt(rowSums((z[e[,1],,drop=FALSE]-z[e[,2],,drop=FALSE])^2))
relative<-function(y,d) sqrt(sum((y-d)^2)/sum(d*d))
profile<-function(y,d) {a<-sum(y*d)/sum(d*d);relative(y,a*d)}
score<-function(z,name,p) {
 edge<-lengths_of(z,p$edges);chord<-lengths_of(z,p$pair_matrix);strict<-p$distance_matrix[p$pair_matrix];ref<-d[p$pair_matrix]
 path<-pilot_path_lengths(z,p$flat_pair_edge_offsets,p$flat_edge_u,p$flat_edge_v)
 a<-sum(chord*strict)/sum(strict^2);b<-sum(edge*p$edge_targets)/sum(p$edge_targets^2)
 check<-grip::score.gmds(z/b,prepared=p,edge_length_epsilon=0)
 err<-max(abs(c(profile(path,p$pair_graph_distance)-check$gmds.stress,profile(edge,p$edge_targets)-check$edge.rel.rmse)))
 ids<-unique(as.integer(seq(1,length(path),length.out=31)))
 explicit<-vapply(ids,function(i) {j<-seq.int(p$flat_pair_edge_offsets[i]+1L,p$flat_pair_edge_offsets[i+1L]);sum(lengths_of(z,cbind(p$flat_edge_u[j]+1L,p$flat_edge_v[j]+1L)))},numeric(1))
 stopifnot(err<1e-10,max(abs(explicit-path[ids]))<1e-9,all(is.finite(z)))
 sv<-svd(scale(z,scale=FALSE),nu=0,nv=0)$d
 data.frame(case=key,surface=surface,sampling=sampling,replicate=replicate,n=n,radius=radius,regime=regime,k=k,method=name,
  path_rel=profile(path,p$pair_graph_distance),edge_rel=profile(edge,p$edge_targets),
  raw_stress=sum((chord-strict)^2)*scale_unit^2,raw_target_rmse=relative(chord,strict),
  stress1=sqrt(sum((chord-a*strict)^2)/sum(chord^2)),graph_reference=relative(strict,ref),
  path_reference=relative(path/b,ref),chord_reference=relative(chord/b,ref),
  procrustes=grip::score.coordinates(z,x,alignment="similarity")$relative_rmse,edge_scale=b,
  sigma1=sv[1]*scale_unit,sigma2=sv[2]*scale_unit,sigma3=sv[3]*scale_unit,
  sigma2_sigma1=sv[2]/sv[1],sigma3_sigma2=sv[3]/sv[2],sigma3_sigma1=sv[3]/sv[1],
  max_edge_relative=max(abs(edge/b/p$edge_targets-1)),package_score_difference=err,
  independent_path_difference=max(abs(explicit-path[ids])))
}
full<-cache("full-classical",list(coords=unname(stats::cmdscale(as.dist(d),k=3)),engine="cmdscale on full smooth distances"))
ks<-if(n==240) c(4,8,16,32,64,128,239) else c(8,16,64,128,479)
if(mode=="pilot") ks<-intersect(ks,c(4,32,239))
if(mode=="extra") ks<-intersect(ks,c(8,32,128,239))
for(regime in c("geodesic","ambient")) {
 base_dist<-if(regime=="geodesic") d else as.matrix(dist(x))
 neighbors<-lapply(seq_len(n),function(i) setdiff(order(base_dist[i,],seq_len(n)),i))
 for(k in ks) {
  prefix<-paste0(regime,"-k",k);final<-file.path(fit_root,key,paste0(prefix,if(mode=="extra") "-extra-result.rds" else "-result.rds"))
  if(file.exists(final)) {stopifnot(identical(readRDS(final)$identity,identity));next}
  if(mode=="extra" && !file.exists(file.path(fit_root,key,paste0(prefix,"-result.rds")))) next
  logmsg("BEGIN",prefix,mode)
  prep<-cache(paste0(prefix,"-prepared"),{
   e<-unique(t(apply(cbind(rep(seq_len(n),each=k),unlist(lapply(neighbors,head,k))),1,sort)))
   e<-e[order(e[,1],e[,2]),,drop=FALSE];g<-igraph::graph_from_edgelist(e,directed=FALSE);components<-igraph::components(g)$no
   bridges<-0L
   if(components>1) {
    dense<-igraph::graph_from_adjacency_matrix(base_dist,mode="undirected",weighted=TRUE,diag=FALSE)
    mst<-igraph::as_edgelist(igraph::mst(dense),names=FALSE);before<-nrow(e)
    e<-unique(rbind(e,t(apply(mst,1,sort))));bridges<-nrow(e)-before;e<-e[order(e[,1],e[,2]),,drop=FALSE]
   }
   w<-base_dist[e];p<-grip::prepare.geodesic.kk(edges=e,n=n,edge_weights=w,tie_mode="single")
   p$path_edges<-p$path_vertices<-p$path_edge_weights<-NULL
   strict<-igraph::distances(igraph::graph_from_edgelist(e,directed=FALSE),weights=w)
   stopifnot(max(abs(strict-p$distance_matrix))<1e-10)
   degree<-tabulate(as.vector(e),nbins=n);q<-function(v) setNames(as.numeric(quantile(v,c(0,.05,.5,.95,1))),c("min","q05","median","q95","max"))
   stats<-data.frame(case=key,regime=regime,k=k,components_before=components,bridges=bridges,edges=nrow(e),
    route_strict_max=max(abs(p$pair_graph_distance-p$distance_matrix[p$pair_matrix])),
    t(c(setNames(q(degree),paste0("degree_",names(q(degree)))),setNames(q(w*scale_unit),paste0("length_",names(q(w)))),
    setNames(q(lengths_of(x[,1:2,drop=FALSE]*scale_unit/radius,e)),paste0("base_separation_",names(q(w)))))),check.names=FALSE)
   list(p=p,stats=stats)
  });p<-prep$p
  cl<-cache(paste0(prefix,"-classical"),slim(grip::classical.mds(prepared=p,dim=3,diagnostics=FALSE)))
  sm<-cache(paste0(prefix,"-stress"),slim(grip::metric.mds(prepared=p,dim=3,init=cl$coords,n_init=3L,max_iter=1000L,eps=1e-8,seed=8300000L+1000L*replicate,diagnostics=FALSE)))
  kk<-function(z,uniform=FALSE,extra=FALSE,fixed=FALSE) slim(grip::edge.kk(coords=z,prepared=p,dim=3,max_iter=if(extra) 2000L else if(uniform) 1000L else 200L,
   stiffness_method="density",density_mix_schedule=if(uniform||extra) 1 else c(0,.25,.5,.75,1),scale_mode=if(fixed) "identity" else "profiled",edge_length_epsilon=0,diagnostics=FALSE,return_trace=TRUE,seed=9300000L+replicate))
  initial<-list(classical=cl$coords,stress=sm$coords)
  if(regime=="geodesic") initial$full_classical<-full$coords
  fits<-list();candidates<-c(list(original=x),initial)
  if(mode!="extra") {
   for(nm in names(initial)) for(control in c("primary","uniform","fixed_primary","fixed_uniform")) {
    label<-paste(nm,control,sep="_");fits[[label]]<-cache(paste0(prefix,"-",label),kk(initial[[nm]],grepl("uniform",control),fixed=grepl("fixed",control)));candidates[[label]]<-fits[[label]]$coords
   }
  } else {
   candidates<-list()
   for(nm in c("classical","stress")) {
    set.seed(10300000L+replicate);perturb<-matrix(rnorm(n*3),n,3);perturb<-perturb/sqrt(mean(perturb^2))*1e-4
    for(control in c("perturbed","random","original","extended","fixed_perturbed","fixed_random","fixed_original","fixed_extended")) {
     label<-paste(nm,control,sep="_")
     z<-switch(sub("fixed_","",control,fixed=TRUE),perturbed=initial[[nm]]+perturb,random=perturb*1e4,original=x,
      extended=readRDS(file.path(fit_root,key,paste0(prefix,"-",nm,if(grepl("fixed",control)) "_fixed_primary.rds" else "_primary.rds")))$value$coords)
     fits[[label]]<-cache(paste0(prefix,"-",label),kk(z,extra=grepl("extended",control),fixed=grepl("fixed",control)));candidates[[label]]<-fits[[label]]$coords
    }
   }
  }
  scores<-do.call(rbind,lapply(names(candidates),function(nm) score(candidates[[nm]],nm,p)))
  if(mode!="extra") {
   stopifnot(scores$raw_stress[scores$method=="stress"]<=scores$raw_stress[scores$method=="classical"]*(1+1e-7)+1e-8)
   if(regime=="ambient") stopifnot(scores$edge_rel[1]<1e-12,scores$path_rel[1]<1e-7)
   if(k==n-1) {
    # Uniform stiffness and fixed scale: edge energy equals all-pair raw stress.
    echeck<-sum((lengths_of(sm$coords,p$edges)-p$edge_targets)^2)
    pcheck<-sum((lengths_of(sm$coords,p$pair_matrix)-base_dist[p$pair_matrix])^2)
    stopifnot(abs(echeck-pcheck)<1e-8)
   }
  }
  result<-list(identity=identity,value=list(scores=scores,candidates=lapply(candidates,function(z) z*scale_unit),
   scale_unit=scale_unit,graph_stats=prep$stats,starts=sm$metadata$starts,selected_start=sm$metadata$selected_start,
   fits=lapply(fits,function(f) {f$coords<-NULL;f}),timing=c(preparation=prep$elapsed,classical=cl$elapsed,stress=sm$elapsed),session=sessionInfo()))
  tmp<-paste0(final,".tmp");saveRDS(result,tmp,compress=TRUE);stopifnot(file.rename(tmp,final))
  logmsg("DONE",prefix,"path",paste(signif(scores$path_rel,3),collapse=","));rm(p,prep,fits,result);gc(FALSE)
 }
}
