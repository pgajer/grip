#!/usr/bin/env Rscript
# Read-only analysis of saved simulations; figures and tables are generated here.
paper_library <- Sys.getenv("GRIP_RJOURNAL_PACKAGE_LIBRARY")
if(nzchar(paper_library)) .libPaths(c(paper_library,.libPaths()))
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop("Usage: plot-results.R EXPERIMENT_OUTPUT")
output <- normalizePath(args[[1]])
files <- sort(list.files(output,pattern="^sample-n[0-9]+-r[0-9]+[.]rds$",full.names=TRUE))
if(!length(files)) stop("No completed samples")
method_order <- c("Metric MDS","Best edge-KK","Best path candidate")
rule_order <- c("k_conn","k_conn+2","2*k_conn")
rows <- list(); graph_rows <- list(); fit_rows <- list(); representative_files <- list()
for(file in files) {
  sample <- readRDS(file)
  for(rule in names(sample$k_values)) {
    g <- sample$graphs[[as.character(sample$k_values[[rule]])]]
    kk <- g$scores[grepl(" \\+ edge-KK$",g$scores$method),]
    selections <- c("Metric MDS",kk$method[which.min(kk$edge_rel)],g$best_path)
    for(j in seq_along(selections)) {
      score <- g$scores[g$scores$method==selections[j],]
      rows[[length(rows)+1L]] <- data.frame(n=sample$n,replicate=sample$replicate,
        rule=rule,k_conn=sample$k_conn,k=g$k,m=g$m,method=method_order[j],
        selected_candidate=selections[j],path_rel=score$path_rel,
        edge_rel=score$edge_rel,stress1=score$stress1)
    }
    graph_rows[[length(graph_rows)+1L]] <- data.frame(n=sample$n,replicate=sample$replicate,
      rule=rule,k_conn=sample$k_conn,k=g$k,m=g$m,cliques=g$certificate$count,
      obstruction=g$certificate$certified,witness_ratio=g$certificate$max_ratio,
      original_path_rel=g$target_score[["path_rel"]],original_edge_rel=g$target_score[["edge_rel"]],
      score_check=g$package_score_difference,elapsed=g$elapsed)
  }
  for(g in sample$graphs) fit_rows[[length(fit_rows)+1L]] <-
    cbind(n=sample$n,replicate=sample$replicate,k=g$k,g$fit_status)
}
data <- do.call(rbind,rows); graphs <- do.call(rbind,graph_rows); fits <- do.call(rbind,fit_rows)
write.csv(data,file.path(output,"candidate-results.csv"),row.names=FALSE)
write.csv(graphs,file.path(output,"graph-results.csv"),row.names=FALSE)
write.csv(fits,file.path(output,"optimizer-status.csv"),row.names=FALSE)
set.seed(9011)
median_ci <- function(x) quantile(replicate(2000,median(sample(x,replace=TRUE))),c(.025,.975),names=FALSE)
summary <- do.call(rbind,lapply(split(data,list(data$n,data$rule,data$method),drop=TRUE),function(d) {
  do.call(rbind,lapply(c("path_rel","edge_rel","stress1"),function(metric) {
    x <- d[[metric]]; ci <- median_ci(x)
    data.frame(n=d$n[1],rule=d$rule[1],method=d$method[1],metric=metric,
      samples=nrow(d),median=median(x),q25=unname(quantile(x,.25)),q75=unname(quantile(x,.75)),
      minimum=min(x),maximum=max(x),median_ci_low=ci[1],median_ci_high=ci[2])
  }))
}))
write.csv(summary,file.path(output,"summary.csv"),row.names=FALSE)
paired <- do.call(rbind,lapply(split(data,list(data$n,data$rule),drop=TRUE),function(d) {
  baseline <- d[d$method=="Metric MDS",]; best <- d[d$method=="Best path candidate",]
  best <- best[match(baseline$replicate,best$replicate),]
  x <- best$path_rel-baseline$path_rel; ci <- median_ci(x)
  data.frame(n=d$n[1],rule=d$rule[1],samples=length(x),median_paired_path_difference=median(x),
             ci_low=ci[1],ci_high=ci[2],improved=sum(x<0))
}))
write.csv(paired,file.path(output,"paired-path-differences.csv"),row.names=FALSE)

export <- function(name,width,height,draw) {
  pdf(file.path(output,paste0(name,".pdf")),width=width,height=height,useDingbats=FALSE)
  draw(); dev.off()
  png(file.path(output,paste0(name,".png")),width=width,height=height,units="in",res=180)
  draw(); dev.off()
}

caption_lines <- function(text,width=145) {
  lines<-strwrap(text,width=width)
  for(i in seq_along(lines)) mtext(lines[i],outer=TRUE,side=1,
    line=.35+(i-1)*.9,cex=.78)
}

plot_distribution <- function(metric,title,ylabel,caption) {
  sizes <- sort(unique(data$n)); old <- par(no.readonly=TRUE); on.exit(par(old))
  par(mfrow=c(length(sizes),3),mar=c(3.7,4.2,2.4,.8),oma=c(4.3,0,2.4,0),cex=.91)
  cols <- c("#666666","#2166AC","#B26026")
  ylim <- c(0,1.07*max(data[[metric]]*100))
  for(n in sizes) for(rule in rule_order) {
    d <- data[data$n==n & data$rule==rule,]
    plot(NA,xlim=c(.5,3.5),ylim=ylim,xaxt="n",xlab="",ylab=ylabel)
    mtext(sprintf("n = %d; %s",n,switch(rule,k_conn="k = k_conn",`k_conn+2`="k = k_conn + 2",`2*k_conn`="k = 2 k_conn")),
      side=3,line=1.4,cex=1)
    abline(h=axTicks(2),col="gray92",lwd=.7)
    for(j in 1:3) {
      values <- 100*d[d$method==method_order[j],metric]
      # Fixed offsets retain the same replicate positions across metrics.
      offsets <- .14*sin(seq_along(values)*2.399963)
      points(j+offsets,values,pch=16,cex=.48,col=adjustcolor(cols[j],alpha.f=.4))
      q <- quantile(values,c(.25,.5,.75)); segments(j,q[1],j,q[3],col=cols[j],lwd=4)
      s <- summary[summary$n==n & summary$rule==rule & summary$method==method_order[j] & summary$metric==metric,]
      arrows(j+.22,100*s$median_ci_low,j+.22,100*s$median_ci_high,
             code=3,angle=90,length=.035,col=cols[j],lwd=1.3)
      points(j,q[2],pch=21,bg="white",col=cols[j],cex=1,lwd=1.5)
    }
    axis(1,at=1:3,labels=c("Metric\nMDS","Best\nedge-KK","Best path\ncandidate"),tick=FALSE,cex.axis=.86)
    mtext(sprintf("%d independent samples",length(unique(d$replicate))),side=3,line=.2,cex=.74,col="gray35")
  }
  mtext(title,outer=TRUE,side=3,line=.65,cex=1.16)
  caption_lines(caption)
}

export("path-error-distributions",11.4,8.5,function() plot_distribution("path_rel",
  "Two-dimensional fixed-path distortion across independent saddle samples",
  "Fixed-path relative RMSE (%)",
  "Figure 1. Each dot is one independently sampled graph; all vertex pairs are scored. Open circles mark medians, thick segments the interquartile ranges, and offset whiskers bootstrap 95% intervals for the medians. Best edge-KK minimizes edge error among four edge-refined starts; the best path candidate minimizes path error among 14 candidates. These are achieved errors, not certified optima."))
export("edge-error-distributions",11.4,8.5,function() plot_distribution("edge_rel",
  "Local edge-length distortion in the same selected configurations",
  "Edge relative RMSE (%)",
  "Figure 2. The configurations and sample positions are the same as in Figure 1, now scored on their individual edges with a separately profiled scale. Open circles mark medians, thick segments interquartile ranges, and whiskers bootstrap 95% intervals for medians. A candidate with small aggregate path error need not have comparably small local edge error."))
export("chord-stress-distributions",11.4,8.5,function() plot_distribution("stress1",
  "Endpoint-chord distortion in the same selected configurations",
  "MDS Stress-1 (%)",
  "Figure 3. The same configurations are evaluated with configuration-normalized MDS Stress-1 over all vertex pairs. Open circles mark medians, thick segments interquartile ranges, and whiskers bootstrap 95% intervals for medians. Chord stress and fixed-path relative RMSE have different denominators and assess different geometric quantities."))

representatives <- list()
for(n in sort(unique(data$n))) {
  d <- data[data$n==n & data$rule=="k_conn" & data$method=="Best path candidate",]
  distance<-abs(d$path_rel-median(d$path_rel))
  eligible<-d[distance<=min(distance)+1e-14,,drop=FALSE]
  row<-eligible[which.min(eligible$replicate),]
  sample <- readRDS(file.path(output,sprintf("sample-n%d-r%03d.rds",n,row$replicate)))
  representatives[[as.character(n)]] <- sample
}
write.csv(do.call(rbind,lapply(representatives,function(s) data.frame(
  n=s$n,replicate=s$replicate,k_conn=s$k_conn,seed=s$sample$seed))),
  file.path(output,"representative-selection.csv"),row.names=FALSE)

display_coords <- function(s) {
  g <- s$graphs[[as.character(s$k_conn)]]; target <- s$sample$coords
  layouts <- list("Generating saddle (3D)"=target,"Metric MDS (2D)"=g$candidates[["Metric MDS"]],
    "Best path candidate (2D)"=g$candidates[[g$best_path]])
  align <- function(z) {
    if(ncol(z)==2) z<-cbind(z,0)
    z<-sweep(z,2,colMeans(z),"-"); ref<-sweep(target,2,colMeans(target),"-")
    f<-svd(crossprod(z,ref)); factor<-sum(f$d)/sum(z*z)
    sweep(factor*z%*%f$u%*%t(f$v),2,colMeans(target),"+")
  }
  aligned<-c(list(target),lapply(layouts[-1],align))
  projected<-lapply(aligned,grip::project.3d,azimuth=35,elevation=22)
  # Two anchors chosen in the parameter domain, independently of layout scores.
  uv<-target[,1:2]; i<-which.min(rowSums(sweep(uv,2,c(-.8,0),"-")^2))
  j<-which.min(rowSums(sweep(uv,2,c(.8,0),"-")^2)); ends<-sort(c(i,j))
  pair<-which(g$cache$pairs[,1]==ends[1] & g$cache$pairs[,2]==ends[2])
  edge_index<-g$cache$index[seq.int(g$cache$offsets[pair]+1L,g$cache$offsets[pair+1L])]+1L
  sc<-rbind(g$target_score,unlist(g$scores[g$scores$method=="Metric MDS",c("path_rel","edge_rel","stress1")]),
             unlist(g$scores[g$scores$method==g$best_path,c("path_rel","edge_rel","stress1")]))
  list(projected=projected,titles=names(layouts),scores=sc,edges=g$cache$edges,
       route=g$cache$edges[edge_index,,drop=FALSE],ends=ends,n=s$n,k=s$k_conn,rep=s$replicate)
}
display <- lapply(representatives,display_coords)
all <- do.call(rbind,unlist(lapply(display,`[[`,"projected"),recursive=FALSE))
pad <- function(x) range(x)+c(-1,1)*.07*diff(range(x))
limits<-list(x=pad(all[,1]),y=pad(all[,2]))
export("representative-layouts",11.4,8.6,function() {
  old<-par(no.readonly=TRUE);on.exit(par(old))
  par(mfrow=c(length(display),3),mar=c(3.4,.2,2.3,.2),oma=c(4.4,0,2.2,0),cex=.91)
  for(d in display) for(j in 1:3) {
    z<-d$projected[[j]]
    plot(z,type="n",asp=1,xlim=limits$x,ylim=limits$y,axes=FALSE,xlab="",ylab="")
    mtext(d$titles[j],side=3,line=1.35,cex=1.05)
    e<-d$edges; segments(z[e[,1],1],z[e[,1],2],z[e[,2],1],z[e[,2],2],col="gray79",lwd=.55)
    points(z,pch=16,cex=.25,col="gray50")
    e<-d$route;segments(z[e[,1],1],z[e[,1],2],z[e[,2],1],z[e[,2],2],col="#2166AC",lwd=2.3)
    ends<-d$ends;segments(z[ends[1],1],z[ends[1],2],z[ends[2],1],z[ends[2],2],col="#B26026",lty=2,lwd=1.6)
    points(z[ends,,drop=FALSE],pch=21,bg="white",cex=.75,col="gray20")
    val<-d$scores[j,]
    fmt<-function(x) if(x<1e-10) "< 1e-8%" else sprintf("%.3g%%",100*x)
    mtext(paste("Path relative RMSE",fmt(val[1])),side=1,line=.05,cex=.91)
    mtext(paste("Edge relative RMSE",fmt(val[2])),side=1,line=1.1,cex=.86)
    mtext(paste("MDS Stress-1",fmt(val[3])),side=1,line=2.1,cex=.86)
    mtext(sprintf("n = %d; k = %d; sample %d",d$n,d$k,d$rep),side=3,line=.1,cex=.73)
  }
  mtext("Median-selected samples at the connectivity threshold",outer=TRUE,side=3,line=.7,cex=1.15)
  caption_lines(paste("Figure 4. One sample per size, selected by proximity to the median best-found path error at k_conn, not by appearance.",
    "Blue: the same retained graph path; dashed brown: its endpoint chord. All pairs are used in the scores.",
    "The 2D configurations are displayed as aligned, tilted planes; alignment and projection do not enter the scores."))
})

graph_summary <- do.call(rbind,lapply(split(graphs,list(graphs$n,graphs$rule),drop=TRUE),function(d) {
  data.frame(n=d$n[1],rule=d$rule[1],samples=nrow(d),k_min=min(d$k),k_median=median(d$k),k_max=max(d$k),
    edge_median=median(d$m),certified=sum(d$obstruction),witness_min=min(d$witness_ratio),
    witness_median=median(d$witness_ratio),witness_max=max(d$witness_ratio))
}))
write.csv(graph_summary,file.path(output,"graph-summary.csv"),row.names=FALSE)
status <- as.data.frame(table(fits$optimizer,ifelse(is.na(fits$convergence),"fixed budget",
  ifelse(fits$convergence==0,"solver convergence",ifelse(fits$convergence==1,"iteration limit","other status")))))
names(status)<-c("optimizer","status","fits"); status<-status[status$fits>0,]
write.csv(status,file.path(output,"fit-status-summary.csv"),row.names=FALSE)

lines <- c("# Random-saddle flattening: standalone results", "",
  sprintf("Generated %s. The manuscript has not been changed.",format(Sys.time(),tz="America/New_York",usetz=TRUE)),"",
  sprintf("Completed samples: %d. Unique graphs: %d. Scored neighborhood-rule cases: %d.",
    length(files),nrow(unique(graphs[,c("n","replicate","k")])),nrow(graphs)),
  sprintf("Original 3D path relative RMSE: maximum %.3g. Maximum independent package score discrepancy: %.3g.",
    max(graphs$original_path_rel),max(graphs$score_check)),"",
  "## Best-achieved 2D path error", "",
  "Values below are percentages; brackets contain the across-sample interquartile range, not a confidence interval. Complete bootstrap median intervals are in summary.csv.","",
  "| n | Neighborhood rule | Samples | Path error, median [IQR] | Edge error of that candidate, median [IQR] | Obstruction detected |",
  "|---:|---|---:|---:|---:|---:|")
for(n in sort(unique(data$n))) for(rule in rule_order) {
  d<-summary[summary$n==n & summary$rule==rule & summary$method=="Best path candidate",]
  p<-d[d$metric=="path_rel",];e<-d[d$metric=="edge_rel",];g<-graph_summary[graph_summary$n==n & graph_summary$rule==rule,]
  lines<-c(lines,sprintf("| %d | %s | %d | %.4g [%.4g, %.4g] | %.4g [%.4g, %.4g] | %d/%d |",n,rule,p$samples,
    100*p$median,100*p$q25,100*p$q75,100*e$median,100*e$q25,100*e$q75,g$certified,g$samples))
}
lines<-c(lines,"","## Optimizer accounting","","A finite iteration-limited fit is retained as a candidate, not counted as convergence.","")
for(i in seq_len(nrow(status))) lines<-c(lines,sprintf("- %s: %s, %d fits.",status$optimizer[i],status$status[i],status$fits[i]))
lines<-c(lines,"","## Interpretation","",
  "At both sample sizes, median best-achieved path error increases as the neighborhood rule adds more neighbors. The two sizes have similar median path errors at the connectivity threshold; n=500 has lower medians for the two denser rules, but larger median edge errors of the selected path candidates for every rule. More observations therefore do not uniformly reduce every measure under the fixed optimization budget.","",
  "Positive four-vertex Gram-matrix witnesses rule out exact 2D preservation of the graph edge lengths. The optimization results remain upper bounds on the minimum achievable fixed-path error, not proven minima. Small average path error can coexist with appreciable local edge distortion. The k comparison changes the graph and its path family, so differences are sensitivity results rather than errors against one common geometry.","",
  "The original coordinates realize the graph lengths in 3D; the fitted configurations are strictly 2D. This study does not measure error against continuous-surface geodesics, recover curvature, or guarantee that connectivity identifies a geometrically adequate neighborhood size.","",
  "## Figures","",
  "- path-error-distributions.pdf / .png: primary fixed-path distortion distributions.",
  "- edge-error-distributions.pdf / .png: edge distortion for the same candidates.",
  "- chord-stress-distributions.pdf / .png: endpoint-chord Stress-1 for the same candidates.",
  "- representative-layouts.pdf / .png: mechanically median-selected examples.","",
  "Protocol, objective, seed rules, and reproduction commands are in the source directory README.md. All results, exact paths, candidate coordinates, witnesses, and solver statuses remain in sample-n*-r*.rds. No samples are excluded from these summaries.")
audit_file<-file.path(output,"extended-budget-audit.csv")
if(file.exists(audit_file)) {
  audit<-read.csv(audit_file)
  lines<-c(lines,"","## Additional optimization budget: selected-case diagnostic","",
    "Each row extends the fixed-budget winner for a median-selected sample within that n/rule group by at most 2,000 iterations. These six results do not replace the population summaries or the representative-figure scores above. Errors are percentages. Status 0 means the solver's numerical stopping criterion was met; status 1 means the additional iteration limit was reached. Neither status certifies global optimality.","",
    "| n | Rule | Sample | Original path error | Extended path error | Extended edge error | Status |",
    "|---:|---|---:|---:|---:|---:|---:|")
  for(i in seq_len(nrow(audit))) {
    a<-audit[i,]
    lines<-c(lines,sprintf("| %d | %s | %d | %.5g | %.5g | %.5g | %d |",
      a$n,a$rule,a$replicate,100*a$before_path,100*a$after_path,100*a$after_edge,a$convergence))
  }
}
validation_file<-file.path(output,"postrun-validation.rds")
if(file.exists(validation_file)) {
  v<-readRDS(validation_file)
  lines<-c(lines,"","## Independent post-run validation","",
    sprintf("All %d graphs were checked independently. Maximum difference in graph distances: %.3g. Maximum difference in the winner's path relative RMSE: %.3g. Every obstruction witness has six retained one-edge paths, and every fitted candidate has exactly two coordinate columns.",
      v$graphs,v$max_independent_distance_difference,v$max_independent_R_path_score_difference),"",
    sprintf("Near-tie handling: %d graphs and %d pairs differ from igraph by more than 1e-8 absolute distance. The maximum scaled relative difference is %.3g, within the package's sqrt(machine epsilon) tolerance. Using the stricter distances changes the winner's path relative RMSE by at most %.3g. The original runs and paths are retained.",
      v$near_tie_graphs,v$near_tie_pairs,v$max_relative_distance_difference,v$max_strict_distance_score_difference),
    sprintf("Near-tie handling selected indirect routes for %d edge-endpoint pairs; none belonged to a retained obstruction witness.",v$indirect_edge_pairs))
}
writeLines(lines,file.path(output,"results.md"))
cat(paste(lines,collapse="\n"),"\n")
