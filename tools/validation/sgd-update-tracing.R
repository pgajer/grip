#!/usr/bin/env Rscript
# Reproducible comparison with the pre-tracing native solver. Run from repo root.
# One R worker; optimized C++17 builds. Outputs are local evidence, not package data.
args<-commandArgs(trailingOnly=TRUE)
revision<-if(length(args)) args[[1]] else "69c533bc"
out<-if(length(args)>1) args[[2]] else "output/sgd-tracing"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1)
flags<-system2("R",c("CMD","config","CXX17FLAGS"),stdout=TRUE)
stopifnot(any(grepl("-O[23]",flags)),!any(grepl("-O0",flags)))
old<-system2("git",c("show",paste0(revision,":src/metric_mds_sgd.cpp")),stdout=TRUE)
stopifnot(is.null(attr(old,"status")))
writeLines(c("// [[Rcpp::plugins(cpp17)]]",old),file.path(out,"baseline-kernel.cpp"))
writeLines(c("// [[Rcpp::plugins(cpp17)]]",readLines("src/metric_mds_sgd.cpp")),file.path(out,"tracing-kernel.cpp"))
baseline<-new.env();updated<-new.env()
Rcpp::sourceCpp(file.path(out,"baseline-kernel.cpp"),env=baseline,rebuild=TRUE,
 cacheDir=file.path(out,"baseline-build"),verbose=TRUE)
Rcpp::sourceCpp(file.path(out,"tracing-kernel.cpp"),env=updated,rebuild=TRUE,
 cacheDir=file.path(out,"tracing-build"),verbose=TRUE)
oldfit<-baseline$grip_sgd_mds_cpp;newfit<-updated$grip_sgd_mds_cpp
set.seed(20261020)
make_case<-function(n,dim,epochs,sparse=FALSE) {
 truth<-matrix(rnorm(n*5),n,5);start<-matrix(rnorm(n*dim),n,dim)
 pairs<-t(utils::combn(n,2));if(sparse)pairs<-pairs[sample(nrow(pairs),min(4000,nrow(pairs))),,drop=FALSE]
 targets<-sqrt(rowSums((truth[pairs[,1],,drop=FALSE]-truth[pairs[,2],,drop=FALSE])^2))
 weights<-1/targets^2
 list(start=start,targets=targets,rates=exp(seq(log(.5),log(.02),length.out=epochs)),
 seed=7123L,checkpoint_every=5L,max_workspace_bytes=512*1024^2,shuffle=TRUE,
 weights=weights,endpoints=if(sparse)pairs else NULL,
 reverse_weights=if(sparse)weights*rep(c(0,1,3),length.out=length(weights)) else NULL,
 retain_best=!sparse,keep_state=TRUE)
}
cases<-list(dense_250_3d=make_case(250,3,10),dense_500_5d=make_case(500,5,10),
 sparse_600_3d=make_case(600,3,30,TRUE))
strip<-function(x){x$update_trace<-NULL;x$workspace_bytes<-NULL;x}
# Equivalence includes default-disabled output, state, RNG, stopping and collisions.
validation<-list()
small<-make_case(10,3,7)
collision<-small;collision$start[2,]<-collision$start[1,];collision$targets[1]<-0
rule<-list(relative.tol=1e10,motion.tol=1e10,residual.tol=1e10,patience=1L,min.epochs=2L)
stopped<-small;stopped$stopping<-rule
for(name in names(c(cases,list(collision=collision,stopped=stopped)))) {
 a<-c(cases,list(collision=collision,stopped=stopped))[[name]]
 original<-do.call(oldfit,a);off<-do.call(newfit,a)
 stopifnot(identical(original,off))
 on<-do.call(newfit,c(a,list(update_trace=list(every=1000,max.bytes=128*1024^2))))
 stopifnot(identical(strip(original),strip(on)))
 # Continue both versions from their identical state, using identical rates.
 a$start<-original$state$coords;a$state<-original$state;a$rates<-c(.015,.01)
 continued.old<-do.call(oldfit,a);continued.new<-do.call(newfit,a)
 stopifnot(identical(continued.old,continued.new))
 validation[[name]]<-list(default_bit_identical=TRUE,traced_numerics_and_rng_identical=TRUE,
   continuation_bit_identical=TRUE,epochs=off$niter,achieved_stress=tail(off$trace$raw_stress,1),
   termination=if(off$stabilized)"practical_stabilization" else "epoch_budget")
}
jsonlite::write_json(validation,file.path(out,"baseline-equivalence.json"),pretty=TRUE,auto_unbox=TRUE,digits=17)
timed<-function(f,a,repeats) {
 begin<-proc.time()[["elapsed"]]
 for(i in seq_len(repeats))invisible(do.call(f,a))
 (proc.time()[["elapsed"]]-begin)/repeats
}
rows<-list();summary<-list();set.seed(61427)
for(name in names(cases)) {
 a<-cases[[name]]
 # Warm both implementations, then time repeated randomized-order paired blocks.
 invisible(do.call(oldfit,a));invisible(do.call(newfit,a))
 calibration<-timed(oldfit,a,1)
 repeats<-min(50L,max(1L,ceiling(.25/max(.001,calibration))))
 for(block in 1:25) {
   order<-sample(c("baseline","untraced"));times<-numeric(2);names(times)<-order
   for(method in order){gc();times[method]<-timed(if(method=="baseline")oldfit else newfit,a,repeats)}
   rows[[length(rows)+1]]<-data.frame(case=name,block=block,repeats=repeats,
    baseline_seconds=times["baseline"],untraced_seconds=times["untraced"])
 }
 tab<-do.call(rbind,rows);tab<-tab[tab$case==name,]
 ratios<-tab$untraced_seconds/tab$baseline_seconds
 boots<-replicate(5000,median(sample(ratios,replace=TRUE)))
 summary[[name]]<-list(baseline_median_seconds=median(tab$baseline_seconds),
   untraced_median_seconds=median(tab$untraced_seconds),paired_median_ratio=median(ratios),
   paired_median_ratio_bootstrap_95=unname(quantile(boots,c(.025,.975))),blocks=25,repeats=repeats)
 write.table(do.call(rbind,rows),file.path(out,"timings.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
 cat(name,": median untraced/baseline =",median(ratios),"\n")
}
jsonlite::write_json(list(baseline_revision=revision,compiler_flags=flags,
 kernel_md5=tools::md5sum(file.path(out,c("baseline-kernel.cpp","tracing-kernel.cpp"))),
 cases=summary,session=capture.output(sessionInfo())),file.path(out,"performance.json"),pretty=TRUE,auto_unbox=TRUE,digits=17)
