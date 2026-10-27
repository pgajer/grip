# Shared, source-visible calculations for the two representation articles.
# Small full-distance examples; no sparse approximation or graph shortest paths.
study_hash <- function(x) {
  f <- tempfile(); on.exit(unlink(f)); saveRDS(x, f, version=3)
  unname(tools::md5sum(f))
}
study_data <- function() {
  e <- new.env(); utils::data('valencia2k', package='linf', envir=e)
  X <- e$valencia2k$rel
  stopifnot(all(is.finite(X)), all(X>=0), all(rowSums(X)>0))
  P <- X/rowSums(X)
  # One composition per fit row; preserve multiplicities separately.
  keys <- apply(P,1,paste,collapse='|'); unique.rows <- which(!duplicated(keys))
  set.seed(20261027)
  ix <- sort(sample(unique.rows, min(300L,length(unique.rows))))
  P <- P[ix,,drop=FALSE]
  li <- match('Lactobacillus_iners',colnames(P)); stopifnot(!is.na(li))
  groups <- colnames(P)[max.col(P,ties.method='first')]
  common <- names(sort(table(groups),decreasing=TRUE))[1:4]
  groups[!groups %in% common] <- 'Other dominant taxa'
  # A small four-part composition: L1 Euclidean targets admit an exact 3D fit.
  set.seed(20261028)
  toy <- matrix(rexp(48*4,rate=1),48,4)
  toy <- toy/rowSums(toy); colnames(toy)<-c('Li','A','B','C')
  rownames(toy)<-paste0('toy',seq_len(nrow(toy)))
  list(P=P, input.rows=ix, source.n=nrow(X), source.distinct=length(unique.rows),
    multiplicity=as.integer(table(factor(keys,levels=keys[ix]))),
    li=li, groups=groups, toy=toy,
    selection=P[,li]>=apply(P,1,max))
}
study_targets <- function(P) {
  Q <- P/apply(P,1,max)
  list('L1 / Euclidean'=as.matrix(dist(P)),
       'Max / Euclidean'=as.matrix(dist(Q)),
       'Max / L-infinity'=as.matrix(dist(Q,method='maximum')))
}
study_score <- function(coords,D) {
  d <- as.numeric(dist(coords)); target <- as.numeric(as.dist(D))
  stopifnot(length(d)==length(target),all(target>0),all(is.finite(d)))
  error <- (d-target)/target
  q <- quantile(target,c(.1,.9),names=FALSE)
  c(RMS_percent=100*sqrt(mean(error^2)),
    Shortest_10pct_RMS=100*sqrt(mean(error[target<=q[1]]^2)),
    Longest_10pct_RMS=100*sqrt(mean(error[target>=q[2]]^2)),
    Median_absolute_percent=100*median(abs(error)))
}
study_fit <- function(D,seed,init='random') {
  cfg <- list(D=D,seed=seed,init=init,epochs=300L,dimensions=3L,
    weights='inverse_squared',scheduler='zheng',schedule.epsilon=.01,
    fingerprint=getOption('grip.article.fingerprint',NULL),
    version=as.character(utils::packageVersion('grip')))
  root <- getOption('grip.article.cache',file.path(tempdir(),'grip-article-cache'))
  dir.create(root,recursive=TRUE,showWarnings=FALSE)
  key <- study_hash(cfg); path <- file.path(root,paste0(key,'.rds'))
  if(file.exists(path)) {
    z<-readRDS(path); stopifnot(identical(z$config,cfg)); return(z$fit)
  }
  warning.text<-character()
  fit <- withCallingHandlers(grip::metric.mds(distance.matrix=D, dim=3,
    init=init,n.init=1,max.iter=300,seed=seed,backend='sgd',
    pair.weights='inverse_squared',diagnostics=FALSE,
    sgd.control=list(scheduler='zheng',schedule.epsilon=.01,
      checkpoint.every=10,max.workspace.bytes=256*1024^2)),warning=function(w) {
        if(!grepl('budget exhaustion does not certify convergence',conditionMessage(w),fixed=TRUE)) return()
        warning.text<<-c(warning.text,conditionMessage(w)); invokeRestart('muffleWarning')
      })
  fit$article.warnings<-warning.text
  score<-study_score(fit$coords,D)
  stopifnot(abs(score[1]/100-fit$metadata$target_normalized_rmse)<1e-8)
  saveRDS(list(config=cfg,fit=fit),path)
  fit
}
study_runs <- function(targets,seeds=c(11L,29L,47L)) {
  fits<-lapply(targets,function(D) lapply(seeds,function(s) study_fit(D,s)))
  rows<-do.call(rbind,lapply(names(fits),function(m) do.call(rbind,
    lapply(seq_along(seeds),function(i) data.frame(Method=m,Seed=seeds[i],
      t(study_score(fits[[m]][[i]]$coords,targets[[m]])),
      Termination=fits[[m]][[i]]$metadata$termination,check.names=FALSE)))))
  best<-lapply(names(fits),function(m) {
    losses<-vapply(fits[[m]],function(f) study_score(f$coords,targets[[m]])[1],numeric(1))
    fits[[m]][[which.min(losses)]]
  }); names(best)<-names(fits)
  list(fits=fits,scores=rows,best=best)
}
study_align <- function(X,reference) {
  A<-scale(X,scale=FALSE); B<-scale(reference,scale=FALSE)
  s<-svd(crossprod(A,B)); Y<-A%*%(s$u%*%t(s$v))
  stopifnot(max(abs(as.numeric(dist(X))-as.numeric(dist(Y))))<1e-10)
  Y
}
study_widgets <- function(coords,groups,descriptions) {
  ref<-coords[[1]]; aligned<-lapply(coords,study_align,reference=ref)
  lim<-t(apply(do.call(rbind,aligned),2,range)); span<-max(lim[,2]-lim[,1])
  lim[,1]<-lim[,1]-.08*span;lim[,2]<-lim[,2]+.08*span
  lev<-sort(unique(groups)); colors<-setNames(c('#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9')[seq_along(lev)],lev)
  groups<-factor(groups,levels=lev)
  sc<-ivue::color.scale.groups(groups,colors=colors)
  lapply(seq_along(aligned),function(i) ivue::plot3D.groups(aligned[[i]],groups,
    scale=sc,point.type='sphere',sphere.radius=.007*span,controls=FALSE,
    legend.title='Dominant taxon',legend.width=210,height=450L,
    aspect='equal',limits=lim,axes=FALSE,xlab='',ylab='',zlab='',
    camera=ivue::camera.zup(elevation=20,turn=-125,zoom=.85),
    description=descriptions[i]))
}
