#!/usr/bin/env Rscript
# One numerical thread; run from the grip source root.
args<-commandArgs(trailingOnly=TRUE)
article<-if(length(args)) args[1] else 'embedding-fidelity'
stopifnot(article %in% c('embedding-fidelity','local-global-embeddings'))
options(rgl.useNULL=TRUE)
root<-normalizePath('.')
# Authoritative sibling sources when present, installed packages otherwise.
for(p in c('linf','ivue')) {
  path<-Sys.getenv(paste0(toupper(p),'_SOURCE'),file.path('..',p))
  if(file.exists(file.path(path,'DESCRIPTION'))) pkgload::load_all(path,quiet=TRUE)
}
if(pkgbuild::needs_compile('.')) pkgbuild::compile_dll('.',debug=FALSE)
pkgload::load_all('.',recompile=FALSE,quiet=TRUE)
out<-if(length(args)>1L) normalizePath(args[2],mustWork=FALSE) else
  file.path(root,'build','representation-guides')
cache.root<-file.path(root,'build','representation-guides','cache')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source.files<-c(list.files('R',full.names=TRUE,pattern='[.]R$'),
  list.files('src',full.names=TRUE,pattern='[.](cpp|h|hpp|so|dll|dylib)$'),
  'vignettes/articles/embedding-study.R')
fingerprint<-tools::md5sum(source.files)
options(grip.article.fingerprint=fingerprint,grip.article.cache=cache.root)
e<-new.env(parent=globalenv())
rmarkdown::render(file.path('vignettes','articles',paste0(article,'.Rmd')),
  output_dir=out,envir=e,quiet=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,paste0(article,'-session.txt')))
write.csv(data.frame(file=names(fingerprint),md5=unname(fingerprint)),
  file.path(out,paste0(article,'-source-fingerprints.csv')),row.names=FALSE)
write.csv(data.frame(row=e$z$input.rows,id=rownames(e$z$P),multiplicity=e$z$multiplicity),
  file.path(out,paste0(article,'-observations.csv')),row.names=FALSE)
if(article=='embedding-fidelity') {
  check<-rbind(c(0,0,0),c(1,0,0),c(0,2,0),c(0,0,3))
  stopifnot(abs(e$study_score(1.1*check,as.matrix(dist(check)))[1]-10)<1e-10)
  stopifnot(e$study_score(check,as.matrix(dist(check)))[1]<1e-10)
  write.csv(e$real$scores,file.path(out,'fidelity-scores.csv'),row.names=FALSE)
  stopifnot(nrow(e$real$scores)==9L,all(is.finite(e$real$scores$RMS_percent)))
}
if(article=='local-global-embeddings') {
  write.csv(e$comparison,file.path(out,'local-global-scores.csv'),row.names=FALSE)
  write.csv(e$local.seed.table,file.path(out,'local-seed-scores.csv'),row.names=FALSE)
  stopifnot(length(e$selected)==sum(e$z$selection),
    all(e$comparison$Local<=e$comparison$Restricted+1e-7))
  for(m in names(e$targets)) stopifnot(identical(e$local.targets[[m]],
    e$targets[[m]][e$selected,e$selected,drop=FALSE]))
}
cat('PASS:',article,'rendered; embedded numerical assertions passed.\n')
