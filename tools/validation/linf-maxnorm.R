# Controlled comparison against preserved experiment constraints. Called by Python.
a <- commandArgs(TRUE)
stopifnot(length(a)==3L)
repo <- a[1]; original <- a[2]; out <- a[3]
pkgload::load_all(repo,recompile=FALSE,quiet=TRUE)
read <- function(path,n,integer=FALSE) readBin(path,if(integer)integer() else numeric(),
  n=n,size=if(integer)4L else 8L,endian="little")
write <- function(x,path) writeBin(as.vector(x),path,size=if(is.integer(x))4L else 8L,endian="little")
reports <- list()
for (scope in c("pilot","full")) {
  task <- jsonlite::fromJSON(file.path(out,scope,"task.json"))
  folder <- file.path(out,scope)
  X <- matrix(read(file.path(folder,"input.bin"),task$n*task$p),task$n,task$p)
  started <- proc.time()[[3]]
  prep <- prepare.sparse.mds(X,method="linf",n.pivots=task$h,q=10,seed=20261025,
    max.workspace.bytes=8*1024^3)
  elapsed <- proc.time()[[3]]-started
  cat(scope,"public preparation:",elapsed,"seconds;",nrow(prep$constraints$pairs),"pairs\n")
  e <- prep$graph$edges
  write(as.matrix(e[,1:2]),file.path(folder,"edges.bin"))
  write(e$length,file.path(folder,"lengths.bin"))
  write(prep$pivots,file.path(folder,"pivots.bin"))
  if(scope=="pilot") for(name in names(prep$constraints))
    write(prep$constraints[[name]],file.path(folder,paste0("public_",name,".bin")))
  # Supplied pivots are used only here, to compare the same historical experiment.
  m <- jsonlite::fromJSON(file.path(original,scope,"linf",paste0("h",task$h),"manifest.json"))
  fixed <- grip:::grip_sparse_prepare_cpp(task$n,as.matrix(e[,1:2]),e$length,task$h,
    as.integer(m$pivot_rows+1L),20261025L,8*1024^3,observations=X,
    save_distances=FALSE,metric="linf")
  for(name in c("pairs","targets","count_i","count_j"))
    write(fixed[[name]],file.path(folder,paste0("fixed_",name,".bin")))
  report <- list(n=task$n,p=task$p,h=task$h,edges=nrow(e),
    public_pairs=nrow(prep$constraints$pairs),fixed_pairs=nrow(fixed$pairs),
    preparation_seconds=elapsed,metadata=prep$metadata)
  if(scope=="pilot") {
    # Refit the historical 200-pivot, seed-17, 200-epoch sparse pilot.
    init <- matrix(read(file.path(original,"pilot/linf/initial17.bin"),task$n*3L),
      task$n,3L,byrow=TRUE)
    fit <- metric.mds(n=task$n,constraints=fixed[c("pairs","targets","count_i","count_j")],
      approximation="sparse",backend="sgd",dim=3,init=init,n.init=1,max.iter=200,
      seed=17,diagnostics=FALSE,pair.weights="inverse_squared",
      sgd.control=list(scheduler="zheng",schedule.epsilon=.1,checkpoint.every=10,
        max.workspace.bytes=8*1024^3))
    write(fit$coords,file.path(folder,"coordinates.bin"))
    report$fit <- list(backend="sgd",seed=17,starts=1,epochs=200,
      sparse_proxy_stress=fit$metadata$sparse_proxy_stress,converged=fit$metadata$converged)
  }
  reports[[scope]] <- report
  jsonlite::write_json(report,file.path(folder,"report.json"),auto_unbox=TRUE,pretty=TRUE,digits=17)
  rm(X,prep,fixed); gc()
}
writeLines(capture.output(sessionInfo()),file.path(out,"session.txt"))
