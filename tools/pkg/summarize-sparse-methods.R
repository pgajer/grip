args <- commandArgs(trailingOnly=TRUE)
input <- args[1]; output <- args[2]
z <- read.csv(file.path(input,"fits.csv"),stringsAsFactors=FALSE)
dir.create(output,recursive=TRUE,showWarnings=FALSE)
keys <- c("case","output_dimension","method","target_metric")
summ <- do.call(rbind,lapply(split(z,interaction(z[keys],drop=TRUE)),function(a) {
  data.frame(a[1,keys],n=nrow(a),median=median(a$rms_relative_percent),
    minimum=min(a$rms_relative_percent),maximum=max(a$rms_relative_percent),
    median_fit_seconds=median(a$fit_seconds))
}))
write.csv(summ,file.path(output,"summary.csv"),row.names=FALSE)
ref <- function(name) if(grepl("geodesic",name)) paste0("dense_",name) else "dense_euclidean"
sparse <- z[!startsWith(z$method,"dense_"),]
sparse$reference <- vapply(sparse$method,ref,"")
sparse$reference_error <- mapply(function(cs,k,seed,m) z$rms_relative_percent[
  z$case==cs & z$output_dimension==k & z$seed==seed & z$method==m],
  sparse$case,sparse$output_dimension,sparse$seed,sparse$reference)
sparse$delta_pp <- sparse$rms_relative_percent-sparse$reference_error
write.csv(sparse,file.path(output,"paired_comparisons.csv"),row.names=FALSE)
case_title <- function(id) {
  d <- sub("quadform_d([23])_index([01])","\\1",id)
  j <- sub("quadform_d([23])_index([01])","\\2",id)
  paste0(d,"D quadratic form, index ",j," (120 profiles)")
}
pdf(file.path(output,"error-comparison.pdf"),width=10,height=8,onefile=TRUE)
for(cs in unique(z$case)) {
  par(mfrow=c(2,2),mar=c(5,4.5,3.5,1),oma=c(2,0,2,0))
  for(k in 2:3) for(metric in c("Euclidean","Graph geodesic")) {
    methods <- if(metric=="Euclidean") c("dense_euclidean","uniform","euclidean_q10","euclidean_farthest_q10") else c("dense_geodesic_q10","geodesic_q10")
    labels <- if(metric=="Euclidean") c("Dense\nreference","Uniform\nlandmarks","Regional\nrandomized","Regional\nfarthest-first") else c("Dense\nreference","Regional\nrandomized")
    a <- z[z$case==cs & z$output_dimension==k & z$method %in% methods,]
    plot(NA,xlim=c(.6,length(methods)+.4),ylim=c(0,max(a$rms_relative_percent)*1.15+.1),
      xaxt="n",xlab="",ylab="RMS relative distance error (%)",main=paste(k,"output dimensions -",metric))
    axis(1,seq_along(methods),labels,cex.axis=.8)
    for(j in seq_along(methods)) {
      v <- a$rms_relative_percent[a$method==methods[j]]
      points(j+c(-.08,0,.08),v,pch=16,col="#777777")
      segments(j,min(v),j,max(v),col="#C74440",lwd=2)
      points(j,median(v),pch=18,col="#C74440",cex=1.4)
    }
  }
  mtext(paste0("Figure ",match(cs,unique(z$case)),". ",case_title(cs)),outer=TRUE,side=3,cex=1.1)
  mtext("q=10, 16 pivots, 100 epochs; gray: 3 seeds; red: median and range (descriptive, not a confidence interval)",outer=TRUE,side=1,cex=.75)
}
dev.off()
pdf(file.path(output,"neighborhood-sensitivity.pdf"),width=10,height=8,onefile=TRUE)
for(cs in unique(z$case)) {
  par(mfrow=c(2,2),mar=c(4.5,4.5,3.5,1),oma=c(3,0,2,0))
  for(k in 2:3) for(metric in c("euclidean","geodesic")) {
    a <- sparse[sparse$case==cs & sparse$output_dimension==k &
      sparse$method %in% paste0(metric,"_q",c(5,10,20)),]
    limits <- range(c(0,a$delta_pp)); span<-max(.1,diff(limits))
    limits<-limits+c(-.12,.35)*span
    plot(NA,xlim=c(.6,3.4),ylim=limits,xaxt="n",xlab="Number of neighbors q",
      ylab="Sparse minus dense error (percentage points)",
      main=paste(k,"output dimensions -",metric))
    axis(1,1:3,c(5,10,20)); abline(h=0,col="#888888",lty=2)
    for(j in 1:3) {
      q<-c(5,10,20)[j]; v<-a$delta_pp[a$method==paste0(metric,"_q",q)]
      points(j+c(-.08,0,.08),v,pch=16,col="#777777")
      segments(j,min(v),j,max(v),col="#C74440",lwd=2)
      points(j,median(v),pch=18,col="#C74440",cex=1.4)
      text(j,max(v)+.13*span,paste0(sum(v<0),"/",length(v)," better"),cex=.8)
    }
  }
  mtext(paste0("Figure ",4+match(cs,unique(z$case)),". ",case_title(cs)),outer=TRUE,side=3,cex=1.1)
  mtext("Lower is better; positive values favor dense. Graph comparisons use the same q-specific target.",outer=TRUE,side=1,line=1,cex=.8)
  mtext("Gray: 3 paired seeds; red: median and range. Descriptive; insufficient replication for population uncertainty claims.",outer=TRUE,side=1,line=0,cex=.75)
}
dev.off()
# Scores for the PCA baseline and neighbor agreement, computed from saved inputs.
bundle<-readRDS(file.path(input,"validation.rds"))
audit<-list(); pca_scores<-list()
for(cs in names(bundle$cases)) {
  X<-bundle$cases[[cs]]$X; D<-as.matrix(dist(X)); P<-prcomp(X)$x
  for(q in c(5,10,20)) {
    prep<-bundle$preparations[[paste(cs,1,paste0("euclidean_q",q),sep="/")]]
    nn<-prep$graph$metadata$knn_index
    boundary<-apply(D,1,function(v)sort(v[v>0])[q])
    agreement<-vapply(seq_len(nrow(X)),function(i) mean(D[i,nn[i,]] <= boundary[i]+1e-12),0)
    audit[[length(audit)+1L]]<-data.frame(case=cs,q=q,search_dimension=prep$pca$dimension,
      retained_variance=prep$pca$achieved,mean_agreement=mean(agreement),minimum_agreement=min(agreement))
  }
  for(k in 2:3) {
    dd<-as.numeric(dist(P[,seq_len(k),drop=FALSE])); target<-as.numeric(as.dist(D))
    pca_scores[[length(pca_scores)+1L]]<-data.frame(case=cs,output_dimension=k,
      rms_relative_percent=100*sqrt(mean((dd/target-1)^2)))
  }
}
write.csv(do.call(rbind,audit),file.path(output,"neighbor-agreement.csv"),row.names=FALSE)
write.csv(do.call(rbind,pca_scores),file.path(output,"pca-baselines.csv"),row.names=FALSE)
