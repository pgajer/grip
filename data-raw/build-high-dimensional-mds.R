#!/usr/bin/env Rscript
# Extract a small, verified teaching bundle; never fit or overwrite study results.
# Usage: Rscript data-raw/build-high-dimensional-mds.R /path/to/high_dimensional_mds_to_3d
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L, file.exists('DESCRIPTION'))
study <- normalizePath(args[[1]])
manifest <- jsonlite::fromJSON(file.path(study,'results/evidence_manifest.json'))$files
used <- character()
input <- function(relative) {
  file <- file.path(study,relative)
  row <- manifest[manifest$path==relative,,drop=FALSE]
  stopifnot(nrow(row)==1L, file.exists(file),
    identical(digest::digest(file=file,algo='sha256'),row$sha256))
  used <<- union(used,relative)
  file
}
scores <- read.csv(input('results/all_scores.csv'),stringsAsFactors=FALSE)
times <- read.csv(input('results/all_time_comparisons.csv'),stringsAsFactors=FALSE)
ids <- c('paraboloid_17_ambient','paraboloid_17_6','HB__illc1033','HB__lock_700')
labels <- c('Paraboloid: ambient distances','Paraboloid: graph-path distances',
            'HB/illc1033','HB/lock_700')
dimensions <- c(4L,6L,8L,10L,15L,20L)
cases <- setNames(vector('list',length(ids)),ids)
validation <- list()
for (i in seq_along(ids)) {
  id <- ids[[i]]; prefix <- paste0('output/',id,'/')
  f <- readRDS(input(paste0(prefix,'fixture.rds')))
  D <- matrix(readBin(input(paste0(prefix,'targets.bin')),'double',n=f$n^2,
                      size=8L,endian='little'),f$n,f$n)
  selected <- scores[scores$dataset==id &
    (scores$route %in% c('classical','direct','classical_refined') |
       (scores$dimension %in% dimensions & scores$route %in% c('high','pca','refined'))),]
  stopifnot(nrow(selected)==61L,all(selected$status=='completed'))
  shown <- selected[selected$start %in% c(0L,17L) & selected$route!='high',]
  coords <- list()
  target <- as.vector(as.dist(D)); eligible <- is.finite(target)
  for (j in seq_len(nrow(shown))) {
    row <- shown[j,]
    z <- as.matrix(read.csv(input(paste0(prefix,row$coordinate_file)),header=FALSE))
    stopifnot(identical(dim(z),c(as.integer(f$n),3L)),all(is.finite(z)))
    stress <- sum((as.vector(dist(z))[eligible]-target[eligible])^2)/sum(target[eligible]^2)
    stopifnot(abs(stress-row$normalized_stress)<1e-10*max(1,row$normalized_stress))
    key <- paste(row$route,row$dimension,sep='_')
    coords[[key]] <- unname(z)
    validation[[paste(id,key,sep='/')]] <- abs(stress-row$normalized_stress)
  }
  cases[[id]] <- list(id=id,label=labels[[i]],n=f$n,
    edges=unname(f$edges),weights=unname(f$weights),vertex_ids=f$ids,
    components=as.integer(f$components),reference=unname(f$X),
    # Only the small synthetic target matrices are needed for the runnable recipe.
    distances=if(i<=2L) D else NULL,coords=coords,
    scores=selected[,c('start','route','dimension','normalized_stress','seconds','epochs','converged')],
    matched_time=times[times$dataset==id & times$dimension %in% dimensions,
      c('seed','dimension','route_seconds','route_normalized_stress','direct_stress',
        'direct_best_checkpoint_seconds','status')])
}
bundle <- list(study='GRIP-EXP-039',date='2026-09-29',display_start=17L,
  starts=c(17L,29L,43L),dimensions=dimensions,cases=cases)
out <- 'inst/extdata/high-dimensional-mds'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
saveRDS(bundle,file.path(out,'examples.rds'),compress='xz',version=2)
stopifnot(identical(bundle,readRDS(file.path(out,'examples.rds'))))
write.csv(manifest[match(used,manifest$path),],file.path(out,'source-hashes.csv'),row.names=FALSE)
cat(length(cases),'cases;',length(validation),'3D layouts independently rescored;',
    'maximum absolute stress difference',max(unlist(validation)),'\n')
