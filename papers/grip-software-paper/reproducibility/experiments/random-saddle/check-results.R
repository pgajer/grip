#!/usr/bin/env Rscript
# Post-run independent checks and a clearly separate extended-budget diagnostic.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1) stop("Usage: check-results.R EXPERIMENT_OUTPUT")
script<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)),"experiment.R"))
output<-normalizePath(args[1]);Rcpp::sourceCpp(file.path(experiment_dir,"objective.cpp"),cacheDir=file.path(output,"cpp-cache"))
data<-read.csv(file.path(output,"candidate-results.csv"),check.names=FALSE)
files<-sort(list.files(output,pattern="^sample-n[0-9]+-r[0-9]+[.]rds$",full.names=TRUE))
max_path_difference<-0; max_distance_difference<-0; count<-0L
max_relative_distance_difference<-0; near_tie_graphs<-0L; near_tie_pairs<-0L
max_strict_distance_score_difference<-0
indirect_edge_pairs<-0L
for(file in files) {
  s<-readRDS(file)
  for(g in s$graphs) {
    cache<-g$cache;edges<-cache$edges
    # Near-tie handling can replace an almost collinear edge with an indirect
    # route. Count these cases; verify all six witness edges independently.
    key<-function(e) paste(pmin(e[,1],e[,2]),pmax(e[,1],e[,2]),sep=":")
    positions<-match(key(edges),key(cache$pairs))
    stopifnot(!anyNA(positions))
    single<-diff(cache$offsets)[positions]==1L
    indirect_edge_pairs<-indirect_edge_pairs+sum(!single)
    stopifnot(all(cache$index[cache$offsets[positions[single]]+1L]+1L==which(single)))
    if(g$certificate$certified) {
      witness_pairs<-t(combn(g$certificate$witness$vertices,2))
      witness_positions<-match(key(witness_pairs),key(cache$pairs))
      stopifnot(!anyNA(witness_positions),all(diff(cache$offsets)[witness_positions]==1L))
    }
    # Independent shortest-path computation, with a separate library.
    ig<-igraph::graph_from_edgelist(edges,directed=FALSE)
    D<-igraph::distances(ig,weights=cache$weights)
    reference_distances<-D[cache$pairs]
    distance_errors<-abs(reference_distances-cache$distances)
    distance_difference<-max(distance_errors)
    relative_distance_difference<-max(distance_errors/pmax(1,reference_distances,cache$distances))
    # grip 0.2.0 uses sqrt(machine epsilon) relative near-tie comparisons in
    # grip.dijkstra.tree(). Report the discrepancy and its score impact instead
    # of claiming exact equality with igraph's stricter distance calculation.
    stopifnot(relative_distance_difference<=sqrt(.Machine$double.eps))
    max_relative_distance_difference<-max(max_relative_distance_difference,relative_distance_difference)
    near_tie_graphs<-near_tie_graphs+(distance_difference>1e-8)
    near_tie_pairs<-near_tie_pairs+sum(distance_errors>1e-8)
    max_distance_difference<-max(max_distance_difference,distance_difference)
    # Independent R path summation for the selected winner.
    z<-g$candidates[[g$best_path]];ell<-edge_lengths(z,edges)
    sums<-c(0,cumsum(ell[cache$index+1L]));y<-diff(sums[cache$offsets+1L])
    reference<-relative_rmse(y,cache$distances)
    strict_score<-relative_rmse(y,reference_distances)
    max_strict_distance_score_difference<-max(max_strict_distance_score_difference,abs(strict_score-reference))
    difference<-abs(reference-min(g$scores$path_rel));stopifnot(difference<1e-8)
    max_path_difference<-max(max_path_difference,difference)
    stopifnot(all(vapply(g$candidates,ncol,integer(1))==2L),identical(s$protocol,protocol))
    count<-count+1L
  }
}
validation<-list(samples=length(files),graphs=count,
  max_independent_distance_difference=max_distance_difference,
  max_relative_distance_difference=max_relative_distance_difference,
  near_tie_graphs=near_tie_graphs,near_tie_pairs=near_tie_pairs,
  max_strict_distance_score_difference=max_strict_distance_score_difference,
  max_independent_R_path_score_difference=max_path_difference,
  indirect_edge_pairs=indirect_edge_pairs,
  every_witness_edge_is_a_retained_single_edge_path=TRUE,all_candidates_strictly_2D=TRUE)
saveRDS(validation,file.path(output,"postrun-validation.rds"));print(validation)

# Median-selected sample within each n/k-rule stratum, using the fixed-budget
# results. Extend only its winning candidate; never silently replace main results.
audit<-list()
for(n in sort(unique(data$n))) for(rule in protocol$k_rules) {
  d<-data[data$n==n & data$rule==rule & data$method=="Best path candidate",]
  distance<-abs(d$path_rel-median(d$path_rel))
  eligible<-d[distance<=min(distance)+1e-14,,drop=FALSE]
  row<-eligible[which.min(eligible$replicate),]
  s<-readRDS(file.path(output,sprintf("sample-n%d-r%03d.rds",n,row$replicate)))
  g<-s$graphs[[as.character(row$k)]]
  started<-proc.time()[[3]];fit<-refine_path(g$candidates[[g$best_path]],g$cache,maxit=2000L)
  after<-scores(fit$coords,g$cache);before<-scores(g$candidates[[g$best_path]],g$cache)
  audit[[length(audit)+1]]<-list(n=n,rule=rule,replicate=row$replicate,k=row$k,
    before=before,after=after,fit=fit,elapsed=proc.time()[[3]]-started)
  cat(sprintf("Extended budget n=%d rule=%s: path %.6g -> %.6g; convergence=%d\n",
    n,rule,before[["path_rel"]],after[["path_rel"]],fit$convergence));flush.console()
}
saveRDS(audit,file.path(output,"extended-budget-audit.rds"))
table<-do.call(rbind,lapply(audit,function(a) data.frame(n=a$n,rule=a$rule,replicate=a$replicate,k=a$k,
  before_path=a$before[["path_rel"]],after_path=a$after[["path_rel"]],before_edge=a$before[["edge_rel"]],
  after_edge=a$after[["edge_rel"]],additional_maxit=2000L,convergence=a$fit$convergence,elapsed=a$elapsed)))
write.csv(table,file.path(output,"extended-budget-audit.csv"),row.names=FALSE)
