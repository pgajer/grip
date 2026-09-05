# Run entry points from the repository root. No private-tree dependencies.
base <- "papers/grip-software-paper"
source_dir <- file.path(base,"reproducibility/experiments/mds-initializer-sensitivity")
pilot_dir <- file.path(base,"reproducibility/experiments/two-fidelity-pilot")
input_dir <- file.path(base,"build/two-fidelity-pilot")
out <- Sys.getenv("GRIP_MDS_SENSITIVITY_OUTPUT",file.path(base,"build/mds-initializer-sensitivity-strict"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
methods <- c("Original saddle","Classical MDS","Stress MDS","Classical MDS + edge-KK","Stress MDS + edge-KK")
protocol <- list(id="saddle-initializers-k-v2",mds_targets="symmetrized saved strict igraph distances",dim=3L,starts=6L,max_iter=1000L,eps=1e-8,
  seed_rule="7300000 + 10000*replicate + start; same random starts across k",
  edge_max_iter=200L,mix=c(0,.25,.5,.75,1),edge_scale_mode="profiled",
  edge_epsilon=0,extra_edge_steps=1000L,extra_mds_steps=2000L,extra_mds_eps=1e-10,
  graph_rule="unique(32, selected_k-5, selected_k, selected_k+5, 80)",
  surface_samples=8000L,surface_subdivisions=2L,surface_seed_rule="1901 + replicate",
  grip=as.character(packageVersion("grip")),smacof=as.character(packageVersion("smacof")))
stopifnot(protocol$grip=="0.2.0.9000",protocol$smacof=="2.1.7")
logmsg <- function(...) {cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z"),...,"\n");flush.console()}
atomic_save <- function(x,path) {tmp<-paste0(path,".tmp");saveRDS(x,tmp,compress=FALSE);stopifnot(file.rename(tmp,path))}
# Cache identity checks the protocol, graph input, and experiment/package sources.
source_id <- function() {
 files<-c(file.path(source_dir,c("common.R","run.R")),
   "R/metric_mds.R","R/gmds_layout_interface.R","R/grip_quality.R",
   list.files("src",pattern="[.]cpp$",full.names=TRUE))
 tools::md5sum(files)
}
lengths_of <- function(z,e) sqrt(rowSums((z[e[,1],,drop=FALSE]-z[e[,2],,drop=FALSE])^2))
relative <- function(y,d) sqrt(sum((y-d)^2)/sum(d*d))
profile <- function(y,d) {a<-sum(y*d)/sum(d*d);relative(y,a*d)}
restore <- function(p,g) {
  p$distance_matrix<-(g$distances+t(g$distances))/2;p$graph_diameter<-max(g$distances)
  p$pair_mode<-"all_pairs";p$tie_mode<-"single"
  class(p)<-c("grip_gmds_prepared","grip_gkk_prepared","grip_geodesic_kk_prepared","list");p
}
slim <- function(f) {f$prepared<-NULL;f$metadata$frames<-NULL;f}
score_layout <- function(z,name,r,k,p,x,ref) {
  dg<-p$pair_graph_distance; strict<-p$distance_matrix[p$pair_matrix]
  idx<-match(ref$pairs[,1]*p$n+ref$pairs[,2],p$pair_matrix[,1]*p$n+p$pair_matrix[,2]);stopifnot(!anyNA(idx))
  edge<-lengths_of(z,p$edges);chord<-lengths_of(z,p$pair_matrix)
  path<-pilot_path_lengths(z,p$flat_pair_edge_offsets,p$flat_edge_u,p$flat_edge_v)
  a<-sum(chord*strict)/sum(strict^2);b<-sum(edge*p$edge_targets)/sum(p$edge_targets^2)
  check<-grip::score.gmds(z,prepared=p,edge_length_epsilon=0)
  err<-max(abs(c(profile(path,dg)-check$gmds.stress,profile(edge,p$edge_targets)-check$edge.rel.rmse)))
  ids<-unique(as.integer(seq(1,length(path),length.out=101)))
  independent<-vapply(ids,function(i) {
    j<-seq.int(p$flat_pair_edge_offsets[i]+1L,p$flat_pair_edge_offsets[i+1L])
    sum(lengths_of(z,cbind(p$flat_edge_u[j]+1L,p$flat_edge_v[j]+1L)))
  },numeric(1))
  independent_error<-max(abs(path[ids]-independent))
  stopifnot(err<1e-10,independent_error<1e-9,all(is.finite(z)))
  sim<-grip::score.coordinates(z,x,alignment="similarity")
  rigid<-grip::score.coordinates(z,x,alignment="rigid")
  sv<-svd(scale(z,scale=FALSE),nu=0,nv=0)$d
  data.frame(replicate=r,k=k,method=name,path_rel=profile(path,dg),edge_rel=profile(edge,p$edge_targets),
    raw_stress=sum((chord-strict)^2),raw_target_rmse=relative(chord,strict),
    raw_optimal_scale=sum(chord*strict)/sum(chord^2),
    stress1=sqrt(sum((chord-a*strict)^2)/sum(chord^2)),
    stress1_identity=sqrt(sum((chord-strict)^2)/sum(chord^2)),
    xg_error=relative(strict[idx],ref$d),xz_path_error=relative(path[idx]/b,ref$d),
    graph_path_error_edge_scale=relative(path[idx]/b,strict[idx]),
    procrustes=sim$relative_rmse,rigid_relative_rmse=rigid$relative_rmse,
    similarity_scale=sim$scale,edge_scale=b,sigma2_sigma1=sv[2]/sv[1],sigma3_sigma1=sv[3]/sv[1],
    package_score_difference=err,independent_path_difference=independent_error)
}
