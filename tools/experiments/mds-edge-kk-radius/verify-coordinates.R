# Recompute dimensional and chord/edge diagnostics from portable coordinates.
src<-"tools/experiments/mds-edge-kk-radius";out<-"output/mds-edge-kk-radius"
smoke<-"--smoke" %in% commandArgs(TRUE)
home<-if(smoke) file.path(out,"verify-smoke") else file.path(src,"summary")
files<-list.files(file.path(home,"coordinates"),full.names=TRUE,pattern="[.]rds$");if(!smoke) stopifnot(length(files)==88)
stopifnot(length(files)>0)
lengths_of<-function(z,e) sqrt(rowSums((z[e[,1],,drop=FALSE]-z[e[,2],,drop=FALSE])^2))
maxdiff<-0;count<-0
for(file in files) {
 cases<-readRDS(file);key<-sub("[.]rds$","",basename(file));d<-as.matrix(read.csv(file.path(out,"inputs",paste0(key,"-distance.csv")),header=FALSE));truth<-as.matrix(read.csv(file.path(out,"inputs",paste0(key,"-truth.csv")),header=FALSE));n<-nrow(d)
 unit<-sqrt(mean(d[upper.tri(d)]^2));pair<-which(upper.tri(d),arr.ind=TRUE);graphcache<-list()
 for(v in cases) {
  id<-v$scores[1,];regime<-id$regime;k<-id$k;gid<-paste(regime,k)
  if(is.null(graphcache[[gid]])) {
   base<-if(regime=="geodesic") d/unit else as.matrix(dist(truth/unit));nn<-lapply(seq_len(n),function(i) head(setdiff(order(base[i,],seq_len(n)),i),k))
   ii<-rep(seq_len(n),each=k);jj<-unlist(nn);e<-unique(cbind(pmin(ii,jj),pmax(ii,jj)));e<-e[order(e[,1],e[,2]),,drop=FALSE]
   if(igraph::components(igraph::graph_from_edgelist(e,directed=FALSE))$no>1) {
    dense<-igraph::graph_from_adjacency_matrix(unname(base),mode="undirected",weighted=TRUE,diag=FALSE);mst<-igraph::as_edgelist(igraph::mst(dense),names=FALSE)
    e<-unique(rbind(e,t(apply(mst,1,sort))));e<-e[order(e[,1],e[,2]),,drop=FALSE]
   }
   w<-base[e]*unit;dg<-igraph::distances(igraph::graph_from_edgelist(e,directed=FALSE),weights=base[e])*unit
   graphcache[[gid]]<-list(e=e,w=w,dg=dg)
  }
  g<-graphcache[[gid]];target<-g$dg[pair]
  for(nm in names(v$candidates)) {
   z<-v$candidates[[nm]];row<-v$scores[v$scores$method==nm,];sv<-svd(scale(z,scale=FALSE),nu=0,nv=0)$d
   chord<-lengths_of(z,pair);edge<-lengths_of(z,g$e);b<-sum(edge*g$w)/sum(g$w*g$w)
   raw<-sum((chord-target)^2);edge_rel<-sqrt(sum((edge-b*g$w)^2)/sum((b*g$w)^2))
   zc<-scale(z,scale=FALSE);xc<-scale(truth,scale=FALSE)
   zn<-zc/sqrt(sum(zc^2));xn<-xc/sqrt(sum(xc^2));rotation<-svd(crossprod(zn,xn))
   coordinate_error<-sqrt(sum((sum(rotation$d)*zn %*% (rotation$u %*% t(rotation$v))-xn)^2))
   a<-sum(chord*target)/sum(target*target);stress1<-sqrt(sum((chord-a*target)^2)/sum(chord*chord))
   differences<-abs(c(sv[1]/max(1,row$sigma1)-row$sigma1/max(1,row$sigma1),sv[2]/sv[1]-row$sigma2_sigma1,
    sv[3]/sv[2]-row$sigma3_sigma2,edge_rel-row$edge_rel,b-row$edge_scale,stress1-row$stress1,coordinate_error-row$procrustes,
    (raw-row$raw_stress)/max(1,row$raw_stress),sqrt(sum((g$dg[pair]-d[pair])^2)/sum(d[pair]^2))-row$graph_reference))
   stopifnot(max(differences)<1e-7,all(is.finite(z)),nrow(z)==n,ncol(z)==3)
   maxdiff<-max(maxdiff,differences);count<-count+1L
  }
 }
}
if(!smoke) stopifnot(count==16928)
if("--record" %in% commandArgs(TRUE)) write.csv(data.frame(candidates=count,max_scaled_discrepancy=maxdiff),file.path(src,"summary/coordinate-validation.csv"),row.names=FALSE)
cat("Portable-coordinate verification passed:",count,"candidates; max scaled discrepancy",maxdiff,"\n")
