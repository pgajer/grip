# Reproducible, deliberately small dense-reference validation. Source defines
# the runner; Rscript executes it. Sparse production code never calls dist().
run_sparse_methods_validation <- function(root, output) {
  Sys.setenv(TZ="America/New_York",OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",
             VECLIB_MAXIMUM_THREADS="1")
  pkgload::load_all(file.path(dirname(root),"dgraphs"),recompile=FALSE,quiet=TRUE)
  pkgload::load_all(root,recompile=FALSE,quiet=TRUE)
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  if (file.exists(file.path(output,"fits.csv"))) stop("Use a new output directory; preserve prior runs")
  code <- c(list.files(file.path(root,"R"),full.names=TRUE,pattern="\\.R$"),
    list.files(file.path(root,"src"),full.names=TRUE,pattern="\\.cpp$"),
    file.path(dirname(root),"dgraphs/R/sknn_graphs.R"),
    file.path(dirname(root),"dgraphs/src/sknn_graphs.cpp"),
    file.path(dirname(root),"dgraphs/src/kNN.cpp"),
    file.path(root,"src/grip.so"),file.path(dirname(root),"dgraphs/src/dgraphs.so"),
    file.path(root,"inst/scripts/sparse-mds-methods-validation.R"))
  fingerprint <- tools::md5sum(code)
  write.csv(data.frame(path=names(fingerprint),md5=unname(fingerprint)),
    file.path(output,"source_manifest.csv"),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(output,"session-info.txt"))
  n <- 120L; epochs <- 100L; h <- 16L; seeds <- 1:3
  pairs <- t(utils::combn(n,2)); storage.mode(pairs)<-"integer"
  rows <- list(); fits <- list(); preparations <- list(); cases <- list(); graph_rows <- list()
  assess <- function(Y,D) {
    target <- D[pairs]
    delta <- sqrt(rowSums((Y[pairs[,1],,drop=FALSE]-Y[pairs[,2],,drop=FALSE])^2))
    100*sqrt(mean((delta/target-1)^2))
  }
  # All-pair constraints use the SAME terminal-iterate SGD interface as sparse
  # fits. This is a full objective, deliberately restricted to 120-point tests.
  dense_constraints <- function(D) list(pairs=pairs,targets=D[pairs],
    count_i=rep(1,nrow(pairs)),count_j=rep(1,nrow(pairs)))
  for (d in 2:3) for (index in 0:1) {
    id <- paste0("quadform_d",d,"_index",index)
    set.seed(20261020+10*d+index)
    U <- matrix(runif(n*d,-1,1),n,d)
    A <- diag(d); if(index==1L) A[1,1]<- -1
    geometry <- dgraphs::synthetic.quadform(d,d+1L,list(.8*A))
    X <- dgraphs::embed.synthetic.geometry(geometry,U)
    cases[[id]] <- list(X=X,U=U,A=.8*A,data_seed=20261020+10*d+index)
    DE <- as.matrix(stats::dist(X))
    for (seed in seeds) {
      specs <- list(uniform=list(method="uniform",q=10L),
        euclidean_farthest_q10=list(method="euclidean",q=10L,selection="farthest"))
      for(q in c(5L,10L,20L)) {
        specs[[paste0("euclidean_q",q)]]<-list(method="euclidean",q=q)
        specs[[paste0("geodesic_q",q)]]<-list(method="geodesic",q=q)
      }
      ps <- lapply(specs,function(z) grip::prepare.sparse.mds(X,method=z$method,
        n.pivots=h,q=z$q,pivot.selection=if(is.null(z$selection)) "randomized" else z$selection,
        neighbor.variance.target=.9,seed=20261030+seed,max.workspace.bytes=256*1024^2))
      targets <- list(euclidean=DE)
      for(q in c(5L,10L,20L)) {
        pp<-ps[[paste0("geodesic_q",q)]]; ed<-dgraphs::graph.edges(pp$graph)
        gg<-igraph::graph_from_data_frame(ed[,1:2],directed=FALSE,vertices=seq_len(n))
        DG<-unname(igraph::distances(gg,weights=ed$length))
        stopifnot(all(is.finite(DG)),all(DG >= DE-1e-10))
        targets[[paste0("geodesic_q",q)]]<-DG
        graph_rows[[length(graph_rows)+1L]]<-data.frame(case=id,seed=seed,q=q,
          components=pp$graph$metadata$n_components_before,
          bridges=pp$graph$metadata$n_mst_edges_added,
          bridge_distance_evaluations=pp$graph$metadata$bridge_distance_evaluations,
          dense_allocated=pp$graph$metadata$dense_distance_matrix_allocated,
          pca_dimension=pp$pca$dimension,pca_variance=pp$pca$achieved)
      }
      for(name in names(ps)) preparations[[paste(id,seed,name,sep="/")]]<-ps[[name]]
      for(k in 2:3) {
        set.seed(20261100+100*d+10*index+seed)
        start<-matrix(rnorm(n*k),n,k)
        start<-start*sqrt(mean(DE[lower.tri(DE)]^2)/(2*k))
        jobs<-names(ps)
        jobs<-c(jobs,"dense_euclidean",paste0("dense_geodesic_q",c(5,10,20)))
        for(name in jobs) {
          target_name<-if(grepl("geodesic",name)) sub("dense_","",name) else "euclidean"
          D<-targets[[target_name]]
          began<-proc.time()[["elapsed"]]
          ans<-tryCatch({
            if(startsWith(name,"dense_")) grip::metric.mds(n=n,
              constraints=dense_constraints(D),approximation="sparse",backend="sgd",
              dim=k,init=start,max.iter=epochs,seed=seed,
              sgd.control=list(scheduler="zheng",schedule.epsilon=.1)) else
              grip::metric.mds(prepared=ps[[name]],approximation="sparse",backend="sgd",
                dim=k,init=start,max.iter=epochs,seed=seed,
                sgd.control=list(scheduler="zheng",schedule.epsilon=.1))
          },error=function(e)e)
          elapsed<-proc.time()[["elapsed"]]-began
          ok<-!inherits(ans,"error")
          row<-data.frame(case=id,intrinsic_dimension=d,index=index,n=n,output_dimension=k,
            seed=seed,method=name,target_metric=target_name,epochs=epochs,status=if(ok) "ok" else "failed",
            error=if(ok) "" else conditionMessage(ans),
            rms_relative_percent=if(ok) assess(ans$coords,D) else NA_real_,
            euclidean_error_percent=if(ok) assess(ans$coords,DE) else NA_real_,
            proxy=if(ok) ans$metadata$sparse_proxy_stress else NA_real_,
            pair_count=if(ok) nrow(ans$metadata$sparse$pairs) else NA_integer_,
            fit_seconds=elapsed,
            preparation_seconds=if(startsWith(name,"dense_")) NA_real_ else ps[[name]]$metadata$preparation_seconds,
            backend="sgd",schedule="zheng",schedule_epsilon=.1,starts=1L)
          key<-paste(id,k,seed,name,sep="/")
          fits[[key]]<-list(coords=if(ok) ans$coords else NULL,row=row,
            effective_control=if(ok) ans$metadata$settings$sgd_control else NULL,
            fitting_seed=if(ok) ans$metadata$fitting_seed else NULL)
          rows[[length(rows)+1L]]<-row
          if(length(rows) %% 24L==0L) {
            utils::write.csv(do.call(rbind,rows),file.path(output,"fits.partial.csv"),row.names=FALSE)
            cat(length(rows),"fits recorded\n"); flush.console()
          }
        }
      }
    }
  }
  results<-do.call(rbind,rows)
  utils::write.csv(results,file.path(output,"fits.csv"),row.names=FALSE)
  utils::write.csv(do.call(rbind,graph_rows),file.path(output,"graphs.csv"),row.names=FALSE)
  saveRDS(list(cases=cases,preparations=preparations,fits=fits,
    generated_at=format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z"),
    specification=list(n=n,epochs=epochs,pivots=h,seeds=seeds,PCA_retained=.9)),
    file.path(output,"validation.rds"),compress="xz")
  stopifnot(nrow(results)==288L,identical(fingerprint,tools::md5sum(code)))
  cat("Completed:",sum(results$status=="ok"),"of",nrow(results),"fits.\n")
  if(any(results$status!="ok")) stop("Some validation fits failed; outputs preserved")
  invisible(results)
}
if (sys.nframe()==0L) {
  args<-commandArgs(trailingOnly=TRUE)
  if(length(args)!=2L) stop("Usage: Rscript sparse-mds-methods-validation.R ROOT OUTPUT")
  run_sparse_methods_validation(normalizePath(args[1]),args[2])
}
