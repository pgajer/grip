#!/usr/bin/env Rscript
# Worker for benchmark-metric-sgd.py. Run against an installed candidate build.
args <- commandArgs(TRUE)
candidate.library <- Sys.getenv("GRIP_BENCHMARK_LIBRARY")
if (!nzchar(candidate.library)) stop("Set GRIP_BENCHMARK_LIBRARY to the installed candidate library")
.libPaths(c(normalizePath(candidate.library, mustWork=TRUE), .libPaths()))
library(grip, lib.loc=candidate.library)
mode <- args[[1]]
out <- args[[2]]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
if (mode == "prepare") {
  stopifnot("backend" %in% names(formals(metric.mds)))
  installed.files <- list.files(system.file(package="grip"),recursive=TRUE,full.names=TRUE)
  installed.files <- installed.files[grepl("(grip\\.(so|dll|dylib)|grip\\.rdb)$",installed.files)]
  write.csv(data.frame(file=basename(installed.files),md5=unname(tools::md5sum(installed.files))),
            file.path(out,"installed-candidate-files.csv"),row.names=FALSE)
  specifications <- list(c("euclidean",64), c("euclidean",256),
                         c("weighted_graph",64), c("weighted_graph",256),
                         c("karate_club",34))
  manifest <- list()
  counter <- 0L
  for (spec in specifications) {
    family <- spec[[1]]; n <- as.integer(spec[[2]])
    started <- Sys.time()
    if (family == "karate_club") {
      edges <- as.matrix(read.csv(system.file("extdata", "karate-club-edges.csv", package="grip")))
      p <- grip:::grip.metric.mds.distance.prepared(edges=edges,n=n)
    } else if (family == "weighted_graph") {
      i <- seq_len(n)
      edges <- rbind(cbind(i,i%%n+1L),cbind(i,(i+6L)%%n+1L))
      weights <- 1 + .3 * sin(seq_len(nrow(edges)))^2
      p <- grip:::grip.metric.mds.distance.prepared(edges=edges,n=n,edge_weights=weights)
    } else {
      set.seed(701+n)
      truth <- matrix(rnorm(n*3),n,3)
      p <- grip:::grip.metric.mds.distance.prepared(edges=edges.cycle(n),n=n)
      p$distance_matrix <- as.matrix(dist(truth))
      p$graph_diameter <- max(p$distance_matrix)
      p$graph_build_mode <- "euclidean_distance_fixture"
    }
    preparation_seconds <- as.numeric(difftime(Sys.time(),started,units="secs"))
    target <- as.double(as.dist(p$distance_matrix))
    for (dimension in 2:3) {
      baseline.started <- Sys.time()
      baseline <- classical.mds(prepared=p,dim=dimension,diagnostics=FALSE)$coords
      observed <- as.double(dist(baseline))
      baseline <- baseline * sum(observed*target)/sum(observed^2)
      baseline_rmse <- sqrt(sum((as.double(dist(baseline))-target)^2)/sum(target^2))
      baseline_seconds <- as.numeric(difftime(Sys.time(),baseline.started,units="secs"))
      for (seed in 1:5) {
        set.seed(1000+seed)
        start <- matrix(rnorm(n*dimension),n,dimension)
        counter <- counter+1L
        id <- sprintf("case-%03d",counter)
        saveRDS(list(prepared=p,start=start,dimension=dimension,seed=seed),file.path(out,paste0(id,".rds")))
        manifest[[counter]] <- data.frame(case=id,family=family,n=n,dimension=dimension,seed=seed,
          preparation_seconds=preparation_seconds,classical_seconds=baseline_seconds,
          classical_rmse=baseline_rmse)
      }
    }
  }
  write.csv(do.call(rbind,manifest),file.path(out,"cases.csv"),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(out,"session-info.txt"))
} else if (mode == "fit") {
  # Load the same optional dependency before timing either backend. This keeps
  # first-use namespace loading out of the optimizer comparison.
  loadNamespace("smacof")
  case <- args[[3]]; backend <- args[[4]]; budget <- as.integer(args[[5]])
  destination <- args[[6]]
  input <- readRDS(file.path(out,paste0(case,".rds")))
  warnings <- character()
  started <- Sys.time()
  result <- tryCatch(withCallingHandlers(metric.mds(prepared=input$prepared,
      dim=input$dimension,init=input$start,seed=input$seed,n_init=1,
      backend=backend,max_iter=budget,diagnostics=FALSE),
    warning=function(w) {warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")}),
    error=function(e) e)
  elapsed <- as.numeric(difftime(Sys.time(),started,units="secs"))
  failed <- inherits(result,"error")
  row <- data.frame(status=if (failed) "error" else "ok",fit_seconds=elapsed,
    rmse=if (failed) NA_real_ else result$metadata$target_normalized_rmse,
    raw_stress=if (failed) NA_real_ else result$metadata$raw_stress,
    initial_raw_stress=if (failed) NA_real_ else result$metadata$starts$initial_raw_stress,
    iterations=if (failed) NA_real_ else result$metadata$starts$iterations,
    pair_updates=if (failed || backend != "sgd") NA_real_ else result$metadata$starts$pair_updates,
    termination=if (failed) "error" else result$metadata$termination,
    warning=paste(warnings,collapse=" | "),error=if (failed) conditionMessage(result) else "")
  write.csv(row,destination,row.names=FALSE)
  if (!failed) saveRDS(result,paste0(destination,".rds"))
} else stop("Unknown worker mode")
