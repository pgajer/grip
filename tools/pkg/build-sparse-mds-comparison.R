#!/usr/bin/env Rscript
# Explicit bounded experiment; not run by package installation or vignettes.
# Usage: Rscript tools/pkg/build-sparse-mds-comparison.R <fresh-dir> <library> <install-log>
args <- commandArgs(TRUE)
if (length(args)!=3L) stop('Supply fresh output directory, candidate library, compilation log')
if (dir.exists(args[1])) stop('Use a fresh directory; previous attempts must be retained')
lib <- normalizePath(args[2]); log <- readLines(args[3])
compile <- grep('-c (sparse_metric_mds|metric_mds_sgd)[.]cpp',log,value=TRUE)
optimization <- regmatches(compile,gregexpr('-O[0-3sgfast]+',compile))
if (length(compile)!=2L || any(vapply(optimization,function(x) !length(x) || tail(x,1L)=='-O0',logical(1)))) stop('Need fresh optimized native compilation')
dir.create(args[1],recursive=TRUE); out <- normalizePath(args[1])
.libPaths(c(lib,.libPaths())); library(grip,lib.loc=lib)
source('inst/scripts/metric-mds-comparison.R'); source('inst/scripts/sparse-mds-comparison.R')
file.copy(c('inst/scripts/metric-mds-comparison.R','inst/scripts/sparse-mds-comparison.R'),out)
protocol <- list(seeds=1:3,epochs=100L,small_pivots=c(8L,32L),large_pivots=c(32L,128L,200L),
  dense_vertex_limit=512L,per_fit_seconds=90,experiment_seconds=1200,
  assessment_seed=4901L,assessment_sources=32L,display_seed=1L,
  status='Provisional calibration choices; independent failures retained',
  timing='Total layout call includes preparation; process peak RSS also includes namespace loading and serialization, not independent evaluation. Phase times exposed only by sparse SGD.',
  starts='Same supplied random matrix for both SGD methods; full MDS profiles scale and selects best checkpoint, sparse returns terminal without rescaling; GRIP uses its own multiscale initialization.',
  settings='SGD default hybrid controls, 100 epochs. Weighted GRIP defaults, length_normalization=none. Equal epoch counts do not imply equal work.',
  source_commit=system('git rev-parse HEAD',intern=TRUE),
  source_md5=as.list(tools::md5sum(c('R/sparse_metric_mds.R','src/sparse_metric_mds.cpp',
    'src/metric_mds_sgd.cpp','inst/scripts/sparse-mds-comparison.R',
    'tools/pkg/build-sparse-mds-comparison.R'))),
  compile=compile,R=R.version.string,system=Sys.info()[c("sysname","release","version","machine")],
  versions=vapply(c('grip','igraph','ivue'),function(p) as.character(packageVersion(p)),''))
saveRDS(protocol,file.path(out,'protocol.rds'))
cases <- sparse_comparison_cases(); bundle <- list(protocol=protocol,cases=cases,fits=list(),evaluation=list())
worker <- file.path(out,'worker.R')
writeLines(c('a <- commandArgs(TRUE); cfg <- readRDS(a[1]); .libPaths(c(cfg$lib,.libPaths()))',
  'library(grip,lib.loc=cfg$lib); source(cfg$recipe)',
  'answer <- tryCatch(sparse_comparison_fit(cfg$case,cfg$method,cfg$seed,cfg$epochs),error=function(e) list(failure=conditionMessage(e)))',
  'saveRDS(answer,a[2])'),worker)
started <- proc.time()[['elapsed']]
for (id in names(cases)) {
  case <- cases[[id]]
  methods <- if (case$n<=protocol$dense_vertex_limit) c('full','sparse8','sparse32','grip') else c('sparse32','sparse128','sparse200','grip')
  graph <- igraph::graph_from_edgelist(case$edges,directed=FALSE)
  sources <- if (case$n<=protocol$dense_vertex_limit) seq_len(case$n) else
    withr::with_seed(protocol$assessment_seed,sample.int(case$n,protocol$assessment_sources))
  # Independent evaluator is outside timed/peak-memory workers, shared across methods.
  D <- igraph::distances(graph,v=sources,weights=case$weights,algorithm='dijkstra')
  bundle$evaluation[[id]] <- list(sources=sources,targets=D,
    mode=if (length(sources)==case$n) 'full' else 'sampled sources')
  for (seed in protocol$seeds) for (method in if (seed %% 2) methods else rev(methods)) {
    key <- paste(id,seed,method,sep='_'); cfg <- file.path(out,paste0(key,'-input.rds'))
    result <- file.path(out,paste0(key,'-result.rds'))
    saveRDS(list(case=case,method=method,seed=seed,epochs=protocol$epochs,lib=lib,
      recipe=file.path(out,'sparse-mds-comparison.R')),cfg)
    remaining <- protocol$experiment_seconds-(proc.time()[['elapsed']]-started)
    is.mac <- Sys.info()[['sysname']]=='Darwin'
    res <- if (remaining<=0) simpleError('Experiment time allowance exhausted') else tryCatch(
      processx::run('/usr/bin/time',c(if(is.mac) '-l' else '-v',file.path(R.home('bin'),'Rscript'),worker,cfg,result),
        timeout=1000*min(remaining,protocol$per_fit_seconds),error_on_status=FALSE),error=function(e)e)
    fail <- if (inherits(res,'error')) conditionMessage(res) else if (res$status!=0) paste('Worker exit',res$status) else ''
    stderr <- if (inherits(res,'error')) '' else res$stderr
    writeLines(stderr,file.path(out,paste0(key,'-resources.txt')))
    answer <- if (file.exists(result)) tryCatch(readRDS(result),error=function(e)
      list(failure=paste('Unreadable worker result:',conditionMessage(e)))) else list()
    if (!is.null(answer$failure)) fail <- answer$failure
    if (is.null(answer$coords) && !nzchar(fail)) fail <- 'Worker returned no coordinates'
    peak.line <- grep('maximum resident set size|Maximum resident set size',strsplit(stderr,'\n')[[1]],value=TRUE)
    peak <- if(length(peak.line)) {
      value <- if(is.mac) as.double(sub('^ *([0-9]+).*','\\1',peak.line[1])) else as.double(sub('.*: *','',peak.line[1]))
      value/(if(is.mac) 1024^2 else 1024)
    } else NA_real_
    error <- NA_real_
    if (!nzchar(fail)) {
      residuals <- lapply(seq_along(sources),function(k) {
        d <- sqrt(rowSums(sweep(answer$coords,2,answer$coords[sources[k],])^2))
        keep <- seq_len(case$n)!=sources[k]
        ((d[keep]-D[k,keep])/D[k,keep])^2
      })
      error <- sqrt(mean(unlist(residuals)))
    }
    row <- data.frame(case=id,seed=seed,method=method,status=if(nzchar(fail)) 'failed' else 'ok',failure=fail,
      relative_error=error,seconds=if(is.null(answer$seconds)) NA_real_ else answer$seconds,
      preparation_seconds=if(is.null(answer$preparation_seconds)) NA_real_ else answer$preparation_seconds,
      fitting_seconds=if(is.null(answer$fitting_seconds)) NA_real_ else answer$fitting_seconds,peak_mib=peak)
    # Full worker results remain in separate private files. The publication
    # bundle needs coordinates/provenance, not every sparse constraint per seed.
    if (!is.null(answer$metadata$sparse)) {
      sp <- answer$metadata$sparse
      answer$metadata$constraint_count <- nrow(sp$pairs)
      answer$metadata$sparse <- sp[c('pivots','region','workspace_estimate_bytes')]
    }
    answer$row <- row; bundle$fits[[key]] <- answer
    saveRDS(bundle,file.path(out,'progress.rds'))
    cat(key,row$status,'seconds',row$seconds,'RSS MiB',peak,'\n'); flush.console()
  }
}
bundle$protocol$total_seconds <- proc.time()[['elapsed']]-started
bundle <- sparse_comparison_profile(bundle)
bundle$protocol$scale_diagnostic <- 'Original and optimally scale-adjusted relative errors are both reported; figures show the adjusted scale.'
saveRDS(bundle,file.path(out,'benchmark.rds'),compress='xz')
write.csv(sparse_comparison_rows(bundle),file.path(out,'results.csv'),row.names=FALSE)
write.csv(sparse_comparison_summary(bundle),file.path(out,'summary.csv'),row.names=FALSE)
