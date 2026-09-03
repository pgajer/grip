#!/usr/bin/env Rscript
lib<-Sys.getenv("GRIP_RJOURNAL_PACKAGE_LIBRARY");if(nzchar(lib)) .libPaths(c(lib,.libPaths()))
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L);out<-normalizePath(args[1])
curves<-read.csv(file.path(out,"calibration-curves.csv"))
validation<-read.csv(file.path(out,"reference-validation.csv"))
fits<-lapply(1:5,function(r) {
  k<-validation$best_k[validation$replicate==r]
  readRDS(file.path(out,sprintf("fit-%02d-k%d.rds",r,k)))
})
scores<-do.call(rbind,lapply(fits,`[[`,"scores"))
audit<-do.call(rbind,lapply(fits,function(f) f$audit$scores))
times<-do.call(rbind,lapply(fits,function(f) data.frame(replicate=f$replicate,k=f$k,t(f$timing))))
write.csv(scores,file.path(out,"layout-scores.csv"),row.names=FALSE)
write.csv(audit,file.path(out,"additional-budget-scores.csv"),row.names=FALSE)
write.csv(times,file.path(out,"layout-timings.csv"),row.names=FALSE)
status<-do.call(rbind,lapply(fits,function(f) do.call(rbind,lapply(c("primary","audit"),function(stage) {
  trace<-if(stage=="primary") f$edge_trace else f$audit$trace
  cap<-if(stage=="primary") 200L else 1000L
  do.call(rbind,lapply(split(trace,trace$stage),function(t) data.frame(replicate=f$replicate,
    run=stage,stage=t$stage[1],iterations=max(t$iteration),budget=cap,
    hit_budget=max(t$iteration)>=cap,final_gradient=tail(t$gradient_norm,1))))
}))))
write.csv(status,file.path(out,"fit-status.csv"),row.names=FALSE)
export<-function(name,w,h,draw) {
  pdf(file.path(out,paste0(name,".pdf")),width=w,height=h,useDingbats=FALSE);draw();dev.off()
  png(file.path(out,paste0(name,".png")),width=w,height=h,units="in",res=180);draw();dev.off()
}
caption<-function(s,width=148) {
  lines<-strwrap(s,width=width)
  for(i in seq_along(lines)) mtext(lines[i],side=1,outer=TRUE,line=.3+.9*(i-1),cex=.78)
}
colors<-c("#1E5C89","#B26026","#777777","#5085AD","#927457")
export("graph-calibration",11.5,5.7,function() {
  par(mfrow=c(1,2),mar=c(4,4.5,2.6,1),oma=c(4.5,0,2.1,0),cex=.92)
  for(zoom in c(FALSE,TRUE)) {
    d<-if(zoom) curves[curves$k>=40,] else curves
    plot(NA,xlim=range(d$k),ylim=c(0,max(d$xg_error)*106),xlab="Number of neighbors k",
      ylab="Surface-to-graph relative RMSE (%)")
    abline(h=axTicks(2),col="gray92",lwd=.6)
    for(r in 1:5) {
      a<-d[d$replicate==r,];lines(a$k,100*a$xg_error,col=colors[r],lty=r,lwd=1.5)
      b<-a[which.min(a$xg_error),];points(b$k,100*b$xg_error,pch=21,bg="white",col=colors[r],cex=1.1,lwd=1.4)
    }
    title(if(zoom) "Neighborhood of the minima" else "Full calibration range",cex.main=1.05,font.main=1)
    if(!zoom) legend("topright",paste("Cloud",1:5),col=colors,lty=1:5,lwd=1.5,bty="n",cex=.87)
  }
  mtext("Calibrating graph geometry against the saddle surface",outer=TRUE,side=3,line=.5,cex=1.15)
  caption("Figure 1. Five independent surface-area-uniform clouds, each with 1,000 observations. Distances from 128 random sources to all other vertices give 119,744 distinct reference pairs per cloud. Each curve compares the repaired symmetric-kNN graph with a numerical surface-distance reference without rescaling. Open circles mark within-cloud minima. The right panel enlarges the denser-graph range; its vertical scale differs.")
})
align<-function(z,x) {
  zz<-sweep(z,2,colMeans(z),"-");xx<-sweep(x,2,colMeans(x),"-");f<-svd(crossprod(zz,xx))
  sweep(sum(f$d)/sum(zz*zz)*zz%*%f$u%*%t(f$v),2,colMeans(x),"+")
}
# Mechanically select the middle cloud by its minimum X->G loss, before viewing layouts.
loss<-vapply(fits,function(f) f$scores$xg_error[1],numeric(1))
representative<-order(loss,seq_along(loss))[3];f<-fits[[representative]]
write.csv(data.frame(replicate=f$replicate,k=f$k,xg_error=loss[representative]),
  file.path(out,"representative-selection.csv"),row.names=FALSE)
plot_layouts<-function(f) {
  z<-lapply(f$candidates,function(z) grip::project.3d(align(z,f$coords),azimuth=35,elevation=22))
  all<-do.call(rbind,z);xr<-range(all[,1]);yr<-range(all[,2]);xr<-xr+c(-1,1)*.035*diff(xr);yr<-yr+c(-1,1)*.035*diff(yr)
  par(mfrow=c(1,3),mar=c(4.8,.15,1.6,.15),oma=c(4.5,0,3.7,0),cex=.94)
  uv<-f$coords[,1:2];i<-which.min(rowSums(sweep(uv,2,c(-.8,0),"-")^2));j<-which.min(rowSums(sweep(uv,2,c(.8,0),"-")^2))
  ends<-sort(c(i,j));p<-f$prepared;pid<-which(p$pair_matrix[,1]==ends[1] & p$pair_matrix[,2]==ends[2])
  idx<-seq.int(p$flat_pair_edge_offsets[pid]+1L,p$flat_pair_edge_offsets[pid+1L])
  route<-cbind(p$flat_edge_u[idx]+1L,p$flat_edge_v[idx]+1L)
  vertex_cols<-colorRampPalette(c("#173D65","#86AFC4","#D9B18B","#8E4921"))(100)[
    pmin(100,pmax(1,1+floor((f$coords[,1]+1)*49.5)))]
  for(m in seq_along(z)) {
    a<-z[[m]];e<-f$graph$edges
    plot(a,type="n",xlim=xr,ylim=yr,asp=1,axes=FALSE,xlab="",ylab="")
    segments(a[e[,1],1],a[e[,1],2],a[e[,2],1],a[e[,2],2],col=adjustcolor("gray45",alpha.f=.12),lwd=.35)
    points(a,pch=16,cex=.34,col=vertex_cols)
    e<-route;segments(a[e[,1],1],a[e[,1],2],a[e[,2],1],a[e[,2],2],col="#1E5C89",lwd=2.4)
    segments(a[ends[1],1],a[ends[1],2],a[ends[2],1],a[ends[2],2],col="#B26026",lwd=1.6,lty=2)
    points(a[ends,],pch=21,cex=.85,col="gray20",bg="white")
    mtext(c("Original saddle sample","Metric MDS","Metric MDS + edge-KK")[m],side=3,line=.2,cex=1.03)
    sc<-f$scores[m,];fmt<-function(v) if(v<1e-10) "< 1e-8%" else sprintf("%.3f%%",100*v)
    for(v in 1:3) mtext(paste(c("Path relative RMSE","Edge relative RMSE","MDS Stress-1")[v],
      fmt(sc[[c("path_rel","edge_rel","stress1")[v]]])),side=1,line=.6+(v-1)*1.08,cex=.9)
    mtext(sprintf("Shape discrepancy: %.2f%%",100*sc$procrustes),side=1,line=3.8,cex=.8,col="gray30")
  }
  mtext("Graph-to-embedding fidelity: three-dimensional configurations",outer=TRUE,side=3,line=1.65,cex=1.12)
  mtext(sprintf("Cloud %d; n = 1,000; selected k = %d; surface-to-graph error = %.3f%%",f$replicate,f$k,100*f$scores$xg_error[1]),
    outer=TRUE,side=3,line=.4,cex=.88)
  caption("Figure 2. The cloud is selected by median minimum surface-to-graph error, not by appearance. All layouts are 3D, similarity-aligned to the original observations and shown from one view. Colors retain original x-coordinate identity. Blue marks one fixed graph path; dashed brown is its endpoint chord. Scores use unprojected coordinates and all 499,500 pairs; shape discrepancy is normalized Procrustes error.")
}
export("representative-layouts",12.2,6.8,function() plot_layouts(f))
export("paired-layout-scores",11.5,5.8,function() {
  par(mfrow=c(1,3),mar=c(4.2,4.2,2.5,.6),oma=c(4.4,0,2.2,0),cex=.92)
  metrics<-c("path_rel","edge_rel","stress1");titles<-c("Fixed-path relative RMSE","Edge relative RMSE","MDS Stress-1")
  for(m in seq_along(metrics)) {
    metric<-metrics[m];ylim<-c(0,max(c(scores[[metric]],audit[[metric]]))*106)
    plot(NA,xlim=c(.8,3.25),ylim=ylim,xaxt="n",xlab="",ylab="Error (%)")
    abline(h=axTicks(2),col="gray92",lwd=.6)
    for(r in 1:5) {
      sc<-scores[scores$replicate==r,];lines(1:3,100*sc[[metric]],type="b",pch=r,col=colors[r],lty=r,cex=.85)
      points(3.13,100*audit[audit$replicate==r,metric],pch=4,col=colors[r],cex=.9)
    }
    axis(1,at=1:3,labels=c("Original\nsaddle","Metric\nMDS","MDS +\nedge-KK"),tick=FALSE,cex.axis=.9)
    title(titles[m],font.main=1,cex.main=1)
  }
  mtext("Paired layout comparisons on the five selected graphs",outer=TRUE,side=3,line=.55,cex=1.13)
  caption("Figure 3. Each line joins configurations of the same cloud on its selected graph. All five cloud-level results are shown; these pilot results are not population confidence intervals. Primary edge-KK uses five stages of at most 200 steps each. Crosses to its right show a separate additional 1,000-step sensitivity run. The metrics have different normalizations, and panel vertical scales differ.")
})
export("two-stage-errors",9.3,5.9,function() {
  par(mar=c(4.5,4.5,2,1),oma=c(4.4,0,1.6,0),cex=.94)
  ys<-cbind(vapply(fits,function(f) f$scores$xg_error[1],numeric(1)),
    vapply(fits,function(f) f$scores$xz_path_error[2],numeric(1)),
    vapply(fits,function(f) f$scores$xz_path_error[3],numeric(1)))
  plot(NA,xlim=c(.85,3.15),ylim=c(0,106*max(ys)),xaxt="n",xlab="",ylab="Relative RMSE against surface distances (%)")
  abline(h=axTicks(2),col="gray92",lwd=.6)
  for(r in 1:5) lines(1:3,100*ys[r,],col=colors[r],lty=r,pch=r,type="b",lwd=1.4,cex=.9)
  axis(1,at=1:3,labels=c("Input graph","MDS embedded paths","MDS + edge-KK paths"),tick=FALSE)
  legend("topright",paste("Cloud",1:5),col=colors,lty=1:5,pch=1:5,bty="n",cex=.85)
  mtext("From surface distances to embedded graph-path lengths",outer=TRUE,side=3,line=.3,cex=1.13)
  caption("Figure 4. All three quantities are evaluated against the same surface-reference pairs. Each layout is converted to input units using its fitted edge scale, determined from graph edges without access to surface distances. No additional scale is fitted against the surface. This end-to-end comparison supplements, rather than replaces, the separate graph-to-embedding diagnostics.",115)
})
graph_data<-readRDS(file.path(out,sprintf("graphs-%02d.rds",f$replicate)))
export("why-both-fidelities-matter",10.8,6.8,function() {
  par(mfrow=c(1,2),mar=c(3.8,.4,2.3,.4),oma=c(4.6,0,2.6,0),cex=.96)
  z<-grip::project.3d(f$coords,azimuth=35,elevation=22)
  reference_sources<-readRDS(file.path(out,sprintf("cloud-%02d.rds",f$replicate)))$sources
  i<-reference_sources[which.min(rowSums(sweep(f$coords[reference_sources,1:2],2,c(-.8,0),"-")^2))]
  j<-which.min(rowSums(sweep(f$coords[,1:2],2,c(.8,0),"-")^2))
  for(k in c(3,f$k)) {
    g<-graph_data$graphs[[as.character(k)]];e<-g$edges
    plot(z,type="n",asp=1,axes=FALSE,xlab="",ylab="")
    segments(z[e[,1],1],z[e[,1],2],z[e[,2],1],z[e[,2],2],col=adjustcolor("gray40",alpha.f=.14),lwd=.4)
    points(z,pch=16,cex=.32,col="gray45")
    ig<-igraph::graph_from_edgelist(e,directed=FALSE)
    route<-as.integer(igraph::shortest_paths(ig,from=i,to=j,weights=g$weights)$vpath[[1]])
    lines(z[route,],col="#1E5C89",lwd=2.2);points(z[c(i,j),],pch=21,bg="white",col="gray15",cex=.8)
    if(g$bridges>0) {
      e<-g$mst_edges;segments(z[e[,1],1],z[e[,1],2],z[e[,2],1],z[e[,2],2],col="#B26026",lwd=2,lty=2)
    }
    xg<-curves$xg_error[curves$replicate==f$replicate & curves$k==k]
    mtext(sprintf("k = %d; %d MST bridges",k,g$bridges),side=3,line=.6,cex=1.03)
    mtext(sprintf("Surface-to-graph error: %.3f%%",100*xg),side=1,line=.3,cex=.96)
    mtext("Original-coordinate path error: 0",side=1,line=1.5,cex=.93)
    mtext(sprintf("Highlighted graph-path length: %.3f",g$distances[i,j]),side=1,line=2.6,cex=.87)
  }
  mtext(sprintf("The same saddle coordinates can preserve a poor graph perfectly (cloud %d)",f$replicate),
    outer=TRUE,side=3,line=.7,cex=1.1)
  caption("Figure 5. The same original 3D observations are shown with a sparse repaired graph and the calibrated graph. Blue highlights a shortest path between the same endpoints; dashed brown marks MST bridges. Original coordinates preserve every Euclidean graph edge and hence every retained graph path exactly for either graph. That zero graph-to-embedding loss does not certify graph-to-surface fidelity. Surface errors use the full reference-pair subset.")
})
rm(graph_data)
lines<-c("# Two-fidelity saddle pilot: results","",sprintf("Generated %s.",format(Sys.time(),tz="America/New_York",usetz=TRUE)),"",
  "Five independent surface-area-uniform clouds, each n=1000. All fitted layouts are 3D. No manuscript files were modified.","",
  "## Graph calibration","","Loss values below are percentages. The numerical surface reference uses 128 random sources and 119,744 distinct pairs per cloud. X->G uses physical units, without scale profiling.","",
  "| Cloud | Selected k | X->G error | MDS path error | MDS + edge-KK path error | MDS edge error | Refined edge error | MDS Stress-1 | Refined Stress-1 |",
  "|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
for(f in fits) {
  a<-f$scores[2,];b<-f$scores[3,]
  lines<-c(lines,sprintf("| %d | %d | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f |",f$replicate,f$k,
    100*a$xg_error,100*a$path_rel,100*b$path_rel,100*a$edge_rel,100*b$edge_rel,100*a$stress1,100*b$stress1))
}
lines<-c(lines,"","## Accounting and numerical checks","",
  sprintf("- %d graph cases; five MDS fits, five primary edge-KK fits, and five separate additional-budget fits. No failed or excluded clouds.",nrow(curves)),
  sprintf("- Primary edge-KK: %d/%d stages reached their iteration cap; additional-budget runs: %d/%d. Earlier termination is not assumed to mean certified convergence.",
    sum(status$hit_budget[status$run=="primary"]),sum(status$run=="primary"),
    sum(status$hit_budget[status$run=="audit"]),sum(status$run=="audit")),
  sprintf("- Original-coordinate edge error <= %.3g; path error <= %.3g (dimensionless).",max(scores$edge_rel[scores$method=="Original saddle"]),max(scores$path_rel[scores$method=="Original saddle"])),
  sprintf("- Maximum disagreement with package scoring: %.3g. Independent R route-sum check: %.3g.",max(scores$package_score_difference),max(scores$independent_path_difference)),
  sprintf("- Smooth-surface BVP spot-check error: %.4g%% to %.4g%% relative RMSE across clouds; each cloud has 128 controls.",100*min(validation$bvp_relative_difference),100*max(validation$bvp_relative_difference)),
  sprintf("- Mesh 81 versus mesh 161 on 16 sources per cloud: %.4g%% to %.4g%% relative difference.",100*min(validation$finer_relative_difference),100*max(validation$finer_relative_difference)),"",
  sprintf("- The extra 128-source fine-mesh check on cloud 5 retains k=%d (primary k=%d); the full-source reference difference is %.4g%%.",validation$fine128_best_k[5],validation$best_k[5],100*validation$fine128_relative_difference[5]),"",
  "## Pilot findings","",
  sprintf("The initial k=3:20 and expanded k=3:40 sweeps had boundary minima. The final k=3:80 sweep gives interior reference minima at %s. Their surface-to-graph errors span %.3f%% to %.3f%%. Nearby k values can be nearly indistinguishable, and smaller reference-source subsets sometimes select a different integer.",
    paste(validation$best_k,collapse=", "),100*min(loss),100*max(loss)),"",
  sprintf("Primary edge-KK reduces fixed-path error in %d/5 clouds and edge error in %d/5; it increases MDS Stress-1 in %d/5. This is a change in which geometric quantity is preserved, not uniform improvement under all criteria.",
    sum(scores$path_rel[scores$method=="MDS + edge-KK"]<scores$path_rel[scores$method=="Metric MDS"]),
    sum(scores$edge_rel[scores$method=="MDS + edge-KK"]<scores$edge_rel[scores$method=="Metric MDS"]),
    sum(scores$stress1[scores$method=="MDS + edge-KK"]>scores$stress1[scores$method=="Metric MDS"])),"",
  sprintf("Procrustes shape discrepancy changes from a median %.1f%% under MDS to %.1f%% after edge-KK. Small graph-path error therefore should not be presented as recovery of the original saddle coordinates. These finite-budget results do not show that exact recovery is impossible: the original coordinates provide an exact zero-loss realization.",
    100*median(scores$procrustes[scores$method=="Metric MDS"]),100*median(scores$procrustes[scores$method=="MDS + edge-KK"])),"",
  sprintf("At k=3, repaired-graph surface errors span %.1f%% to %.1f%% even though the original coordinates preserve all graph paths exactly. Figure 5 makes that separation explicit.",
    100*min(curves$xg_error[curves$k==3]),100*max(curves$xg_error[curves$k==3])),"",
  "## Interpretation limits","",
  "The per-cloud k is an oracle choice using a numerical reference for the known surface. It is not an automatic rule available for unknown real-data geometry. Calibration uses a random source subset, while graph-to-embedding scoring uses all unordered pairs. Exact mesh distances are approximations to smooth-surface distances; mesh comparisons, smooth BVP checks, and source-subset sensitivity are recorded separately. Positive fitted errors are achieved values, not proven minima. Five clouds establish pilot behavior, not population precision.","",
  "The original 3D observations realize all Euclidean graph edges and every retained path exactly in mathematics, including MST bridges. This does not imply zero chord stress, exact surface-geodesic approximation, or shape recovery by the fitted methods. Endpoint chords and surface geodesics are distinct quantities.","",
  "## Files","",
  "- graph-calibration.pdf/png: complete error-versus-k curves and enlarged minima.",
  "- representative-layouts.pdf/png: median-selected three-panel 3D comparison.",
  "- paired-layout-scores.pdf/png: all five clouds and additional-budget diagnostics.",
  "- two-stage-errors.pdf/png: surface-to-graph and end-to-end embedded-path errors.",
  "- why-both-fidelities-matter.pdf/png: the same original coordinates preserve a poor graph and a calibrated graph exactly.",
  "- calibration-curves.csv, reference-validation.csv, layout-scores.csv, additional-budget-scores.csv, layout-timings.csv, fit-status.csv: complete numeric summaries.",
  "- Source README.md: protocol, equations, reproduction commands, dependency versions, and interpretation.")
writeLines(lines,file.path(out,"results.md"))
print(scores[,c("replicate","k","method","path_rel","edge_rel","stress1","xz_path_error","procrustes")])
