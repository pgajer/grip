# Chart contract: single cloud, ordered k; compare graph fidelity and reference
# recovery, without treating k values as independent replicates. Static R PNG.
# Two method colors plus neutral baseline, with redundant line styles; zero axes.
plot.saddle.reference <- function(scores, file, columns=4L) {
  d <- subset(scores,sample_size==max(sample_size))
  methods <- c("Original saddle","Metric MDS","MDS + edge-KK")
  colors <- c("#555555","#2166AC","#B35806"); types <- c(3,1,2)
  stopifnot(columns %in% c(2L,4L))
  grDevices::png(file,width=if(columns==4L)2000 else 1400,
                height=if(columns==4L)1120 else 1800,res=160)
  on.exit(grDevices::dev.off())
  par(mfrow=c(8/columns,columns),mar=c(3.6,4.2,2.5,.8),oma=c(6,0,3.8,0),mgp=c(2.6,.8,0),las=1)
  specs <- list(c("path_rel","rigid","Fixed-path relative RMSE","%"),
    c("edge_rel","rigid","Edge relative RMSE","%"),
    c("stress1","rigid","MDS Stress-1","%"),
    c("alignment_scale","similarity","Fitted coordinate scale","Scale"),
    c("coordinate_relative_rmse","rigid","Coordinate error: rigid","% of reference RMS radius"),
    c("coordinate_relative_rmse","similarity","Coordinate error: similarity","% of reference RMS radius"),
    c("surface_rms","rigid","Surface distance: rigid","Coordinate units"),
    c("surface_rms","similarity","Surface distance: similarity","Coordinate units"))
  for(spec in specs) {
    z <- d[d$alignment==spec[2],]; factor <- if(startsWith(spec[4],"%")) 100 else 1
    value <- z[[spec[1]]]*factor
    plot(range(z$k),c(0,max(value)*1.08),type="n",xlab="Neighborhood size k",ylab=spec[4],
         main=spec[3],cex.main=.92,bty="l")
    abline(h=axTicks(2),col="gray92",lwd=.7)
    for(i in seq_along(methods)) {
      m <- z[z$method==methods[i],]; m <- m[order(m$k),]
      lines(m$k,m[[spec[1]]]*factor,col=colors[i],lty=types[i],lwd=2)
    }
  }
  mtext(if(columns==4L) "All graph pairs scored; graph losses use profiled scales. Surface scores: 8,000 samples per direction; matched parameter footprint." else
        "All graph pairs; profiled graph scales. Surface RMS: 8,000 samples per direction; matched footprint.",
        outer=TRUE,side=1,line=2.3,cex=.76)
  mtext(if(columns==4L) "Rigid alignment preserves size; similarity alignment fits scale. Curves describe one sample, not sampling-population estimates." else
        "Rigid alignment preserves size; similarity fits scale. One cloud, not independent replicates across k.",
        outer=TRUE,side=1,line=3.3,cex=.76)
  par(fig=c(0,1,0,1),mar=rep(0,4),oma=rep(0,4),new=TRUE)
  plot.new(); plot.window(xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i")
  text(.5,.98,sprintf("Graph fidelity and reference-saddle recovery | one cloud, n = %d",d$n[1]),
       adj=c(.5,1),cex=1.4)
  legend(.5,if(columns==4L).09 else .055,xjust=.5,yjust=.5,legend=c("Original saddle","metric-MDS","metric-MDS + edge-KK"),col=colors,lty=types,lwd=2,
         bty="n",horiz=TRUE,cex=.9)
}
