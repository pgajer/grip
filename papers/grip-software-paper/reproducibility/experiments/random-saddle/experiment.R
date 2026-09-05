#!/usr/bin/env Rscript
# Standalone experiment: no manuscript or existing benchmark is modified.
paper_library <- Sys.getenv("GRIP_RJOURNAL_PACKAGE_LIBRARY")
if (nzchar(paper_library)) .libPaths(c(paper_library, .libPaths()))
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
experiment_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else getwd()

protocol <- list(version = 1L, amplitude = 0.8, domain = c(-1, 1),
  sampling = "uniform surface area; rejection sampling", noise = 0,
  k_rules = c("k_conn", "k_conn+2", "2*k_conn"),
  starts = c("Metric MDS", "PCA", "Weighted GRIP 1", "Weighted GRIP 2"),
  edge_iterations_per_stage = 50L, edge_density_schedule = c(0, .25, .5, .75, 1),
  path_maxit = 200L, path_factr = 1e5, path_pgtol = 1e-8,
  path_starts = "four edge-KK results plus unrefined MDS and PCA",
  pair_policy = "all unordered pairs; one retained input shortest path",
  scores = "epsilon=0; separate unweighted profiled scales; configuration-normalized Stress-1",
  certificate_eigen_ratio_tolerance = 1e-8,
  sample_seed = "1000000 + 1000*n + replicate", rng = "Mersenne-Twister/Inversion/Rejection")

sample_saddle <- function(n, seed, amplitude = protocol$amplitude) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed)
  accepted <- matrix(numeric(), ncol = 2); proposed <- 0L
  max_jacobian <- sqrt(1 + 8*amplitude^2)
  while (nrow(accepted) < n) {
    batch <- max(64L, 2L*(n-nrow(accepted)))
    uv <- matrix(runif(2L*batch, -1, 1), ncol = 2)
    keep <- runif(batch) < sqrt(1+4*amplitude^2*rowSums(uv^2))/max_jacobian
    accepted <- rbind(accepted, uv[keep,,drop=FALSE]); proposed <- proposed+batch
  }
  uv <- accepted[seq_len(n),,drop=FALSE]
  list(coords = cbind(x=uv[,1], y=uv[,2], z=amplitude*(uv[,1]^2-uv[,2]^2)),
       seed=seed, proposals=proposed)
}

sknn <- function(x, k) dgraphs::create.sknn.graph(x, k=k,
  neighbor.method="exact", connect.components=FALSE, prune.method="none",
  prune.edges=FALSE, edge.weight="distance")

connected_threshold <- function(x) {
  # Sequential search certifies that no smaller k is connected.
  for(k in seq_len(nrow(x)-1L)) {
    g <- sknn(x,k)
    if (g$n_components_before == 1L) return(k)
  }
  stop("No connected graph")
}

edge_lengths <- function(z,e) sqrt(rowSums((z[e[,1],,drop=FALSE]-z[e[,2],,drop=FALSE])^2))
relative_rmse <- function(y,g) {
  a <- sum(y*g)/sum(g*g)
  sqrt(sum((y-a*g)^2)/sum((a*g)^2))
}

path_cache <- function(p) {
  key <- function(e) paste(pmin(e[,1],e[,2]),pmax(e[,1],e[,2]),sep=":")
  edge_key <- key(p$edges)
  index <- lapply(p$path_edges,function(e) match(key(e),edge_key)-1L)
  stopifnot(!anyNA(unlist(index)))
  list(n=p$n, edges=p$edges, weights=p$edge_targets, pairs=p$pair_matrix,
    distances=p$pair_graph_distance, offsets=as.integer(c(0,cumsum(lengths(index)))),
    index=as.integer(unlist(index,use.names=FALSE)))
}

evaluate <- function(z, cache, edge_only=FALSE) saddle_objective(as.numeric(z),nrow(z),
  cache$edges,cache$offsets,cache$index,
  if(edge_only) cache$weights else cache$distances,edge_only)

scores <- function(z,cache) {
  pairs <- cache$pairs; g <- cache$distances
  chords <- edge_lengths(z,pairs); a <- sum(chords*g)/sum(g*g)
  c(path_rel=sqrt(evaluate(z,cache)$value),
    edge_rel=relative_rmse(edge_lengths(z,cache$edges),cache$weights),
    stress1=sqrt(sum((chords-a*g)^2)/sum(chords^2)))
}

normalize_coords <- function(z,cache) {
  z <- sweep(z,2,colMeans(z),"-")
  z/evaluate(z,cache,TRUE)$scale
}

refine_path <- function(z,cache,maxit=protocol$path_maxit) {
  z <- normalize_coords(z,cache); previous <- NULL; result <- NULL
  calc <- function(par) {
    if(!identical(par,previous)) {
      result <<- evaluate(matrix(par,nrow(z)),cache); previous <<- par
    }
    result
  }
  start_value <- calc(as.numeric(z))$value
  fit <- optim(as.numeric(z), function(v) calc(v)$value,
    function(v) calc(v)$gradient, method="L-BFGS-B",
    control=list(maxit=maxit, factr=protocol$path_factr, pgtol=protocol$path_pgtol))
  if(!is.finite(fit$value) || fit$value > start_value+1e-10) stop("Path optimizer increased objective")
  list(coords=matrix(fit$par,nrow(z)), convergence=fit$convergence,
       message=fit$message, counts=fit$counts, initial_value=start_value, value=fit$value)
}

clique_obstruction <- function(x,edges,weights) {
  graph <- igraph::graph_from_edgelist(edges,directed=FALSE)
  cliques <- igraph::cliques(graph,min=4,max=4)
  if(!length(cliques)) return(list(count=0L,certified=FALSE,max_ratio=NA_real_,witness=NULL))
  weight_matrix <- matrix(NA_real_,nrow(x),nrow(x))
  weight_matrix[edges] <- weights; weight_matrix[edges[,2:1]] <- weights
  diag(weight_matrix) <- 0
  H <- diag(4)-matrix(1/4,4,4)
  ratios <- vapply(cliques,function(ids) {
    ids <- as.integer(ids); d <- weight_matrix[ids,ids]
    stopifnot(all(is.finite(d)))
    eig <- eigen(-.5*H%*%(d*d)%*%H,symmetric=TRUE,only.values=TRUE)$values
    eig[3]/eig[1]
  },numeric(1))
  selected <- which.max(ratios); ids <- as.integer(cliques[[selected]])
  list(count=length(cliques),certified=ratios[selected] > protocol$certificate_eigen_ratio_tolerance,
       max_ratio=ratios[selected], certified_count=sum(ratios>protocol$certificate_eigen_ratio_tolerance),
       witness=list(vertices=ids,distances=weight_matrix[ids,ids],coords=x[ids,,drop=FALSE]))
}

run_graph <- function(x,k,seed) {
  started <- proc.time()[[3]]; g <- sknn(x,k)
  stopifnot(g$n_components_before==1L,g$n_mst_edges_added==0L)
  scale <- median(g$edge_weight); weights <- g$edge_weight/scale; n <- nrow(x)
  p <- grip::prepare.geodesic.kk(edges=g$edge_matrix,n=n,edge_weights=weights,tie_mode="single")
  cache <- path_cache(p); cache$n <- n
  stopifnot(nrow(cache$pairs)==choose(n,2),max(abs(edge_lengths(x/scale,cache$edges)-cache$weights))<1e-10)
  target_score <- scores(x,cache)
  stopifnot(target_score[["path_rel"]]<1e-10,target_score[["edge_rel"]]<1e-10)
  # A selected single-edge route must realize every input edge in this generic sample.
  input_paths <- vapply(p$path_edges,function(e) sum(edge_lengths(x/scale,e)),numeric(1))
  stopifnot(max(abs(input_paths-cache$distances))<1e-8)
  certificate <- clique_obstruction(x,cache$edges,cache$weights)
  initials <- list("Metric MDS"=grip::classical.mds(prepared=p,dim=2,diagnostics=FALSE)$coords,
    "PCA"=prcomp(x,center=TRUE,scale.=FALSE)$x[,1:2,drop=FALSE])
  for(i in 1:2) initials[[paste("Weighted GRIP",i)]] <- grip::grip(
    cache$edges,n=n,edge_weights=cache$weights,dim=2,metric="edge_length",seed=seed+i)
  candidates <- initials; fit_rows <- list()
  for(name in names(initials)) {
    tm <- proc.time()[[3]]
    fit <- grip::edge.kk(coords=normalize_coords(initials[[name]],cache),prepared=p,dim=2,
      stiffness_method="density",density_mix_schedule=protocol$edge_density_schedule,
      max_iter=protocol$edge_iterations_per_stage,scale_mode="profiled",seed=seed,
      diagnostics=FALSE,return_trace=TRUE)
    candidates[[paste0(name," + edge-KK")]] <- fit$coords
    tr <- fit$trace
    fit_rows[[length(fit_rows)+1]] <- data.frame(method=paste0(name," + edge-KK"),
      optimizer="edge-KK",convergence=NA_integer_,evaluations=nrow(tr),
      message="Fixed package continuation budget; no convergence code returned",
      elapsed=proc.time()[[3]]-tm)
  }
  path_names <- c(paste0(names(initials)," + edge-KK"),"Metric MDS","PCA")
  for(name in path_names) {
    tm <- proc.time()[[3]]; fit <- refine_path(candidates[[name]],cache)
    label <- paste0(name," + path refinement"); candidates[[label]] <- fit$coords
    fit_rows[[length(fit_rows)+1]] <- data.frame(method=label,optimizer="L-BFGS-B",
      convergence=fit$convergence,evaluations=unname(fit$counts[[1]]),
      message=if(is.null(fit$message)) "" else fit$message,elapsed=proc.time()[[3]]-tm)
  }
  all_scores <- do.call(rbind,lapply(names(candidates),function(name) {
    cbind(method=name,as.data.frame(as.list(scores(candidates[[name]],cache))))
  }))
  best_path <- all_scores$method[which.min(all_scores$path_rel)]
  best_edge <- all_scores$method[which.min(all_scores$edge_rel)]
  # Check the reported winner against the package scorer, using epsilon=0.
  independent <- grip::score.gmds(candidates[[best_path]],prepared=p,edge_length_epsilon=0)
  stopifnot(abs(independent$gmds.stress-min(all_scores$path_rel))<1e-10)
  list(k=k,n=n,m=nrow(cache$edges),edge_median=scale,certificate=certificate,
       cache=cache,candidates=candidates,scores=all_scores,target_score=target_score,
       best_path=best_path,best_edge=best_edge,fit_status=do.call(rbind,fit_rows),
       package_score_difference=abs(independent$gmds.stress-min(all_scores$path_rel)),
       elapsed=proc.time()[[3]]-started)
}

run_replicate <- function(n,replicate,output) {
  filename <- file.path(output,sprintf("sample-n%d-r%03d.rds",n,replicate))
  if(file.exists(filename)) {
    previous <- readRDS(filename)
    stopifnot(identical(previous$protocol,protocol)); return(filename)
  }
  seed <- 1000000L+1000L*n+replicate
  sample <- sample_saddle(n,seed); x <- sample$coords; kc <- connected_threshold(x)
  ks <- c(kc,kc+2L,2L*kc)
  graphs <- lapply(unique(ks),function(k) run_graph(x,k,seed+100L*k))
  names(graphs) <- as.character(unique(ks))
  result <- list(protocol=protocol,n=n,replicate=replicate,sample=sample,k_conn=kc,
    k_values=setNames(ks,protocol$k_rules),graphs=graphs,session=sessionInfo(),
    completed=format(Sys.time(),tz="UTC",usetz=TRUE))
  temp <- paste0(filename,".partial"); saveRDS(result,temp,compress="gzip")
  stopifnot(file.rename(temp,filename))
  cat(sprintf("Completed n=%d replicate=%d k_conn=%d; %.1f graph-seconds\n",
    n,replicate,kc,sum(vapply(graphs,`[[`,numeric(1),"elapsed"))))
  flush.console(); filename
}

run_checks <- function(output) {
  x <- sample_saddle(14,1001)$coords; k <- connected_threshold(x)
  g <- sknn(x,k); p <- grip::prepare.geodesic.kk(g$edge_matrix,n=nrow(x),edge_weights=g$edge_weight)
  cache <- path_cache(p); z <- x[,1:2]; set.seed(1002); z <- z+matrix(rnorm(length(z),sd=.1),ncol=2)
  checks <- list()
  for(edge_only in c(FALSE,TRUE)) {
    analytic <- evaluate(z,cache,edge_only)$gradient; numeric_gradient <- numeric(length(z)); h <- 1e-6
    for(j in seq_along(z)) {
      plus <- minus <- z; plus[j] <- plus[j]+h; minus[j] <- minus[j]-h
      numeric_gradient[j] <- (evaluate(plus,cache,edge_only)$value-evaluate(minus,cache,edge_only)$value)/(2*h)
    }
    error <- max(abs(analytic-numeric_gradient)); stopifnot(error<1e-6)
    checks[[if(edge_only) "edge_gradient_error" else "path_gradient_error"]] <- error
  }
  stopifnot(abs(evaluate(z*7,cache)$value-evaluate(z,cache)$value)<1e-12)
  pscore <- grip::score.gmds(z,prepared=p,edge_length_epsilon=0)
  stopifnot(abs(pscore$gmds.stress-scores(z,cache)[["path_rel"]])<1e-12)
  # Plane control: same Euclidean-weight construction has an exact 2D solution.
  plane <- sample_saddle(60,1003,amplitude=0)$coords
  pg <- sknn(plane,connected_threshold(plane))
  pp <- grip::prepare.geodesic.kk(pg$edge_matrix,n=nrow(plane),edge_weights=pg$edge_weight)
  pc <- path_cache(pp); plane_error <- scores(plane[,1:2],pc)[["path_rel"]]
  witness <- clique_obstruction(plane,pc$edges,pc$weights)
  stopifnot(plane_error<1e-12,!witness$certified)
  checks$plane_path_error <- plane_error
  checks$plane_obstruction <- witness$certified
  # A tetrahedron is impossible in 2D even when every edge is a graph geodesic.
  tetra <- rbind(c(0,0,0),c(1,0,0),c(0,1,0),c(0,0,1)); te <- t(combn(4,2))
  stopifnot(clique_obstruction(tetra,te,edge_lengths(tetra,te))$certified)
  checks$tetrahedron_obstruction <- TRUE
  # Sampling Jacobian check by deterministic quadrature versus a large independent draw.
  nodes <- seq(-1,1,length.out=401); uv <- as.matrix(expand.grid(nodes,nodes))
  w <- sqrt(1+4*.8^2*rowSums(uv^2)); expected <- weighted.mean(rowSums(uv^2),w)
  empirical <- mean(rowSums(sample_saddle(100000,1004)$coords[,1:2]^2))
  stopifnot(abs(expected-empirical)<.01)
  checks$surface_area_expected_radius2 <- expected; checks$sample_mean_radius2 <- empirical
  saveRDS(checks,file.path(output,"validation-checks.rds")); print(checks)
}

if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)<2) stop("Usage: experiment.R checks|run OUTPUT [SIZES=250] [REPS=20] [WORKERS=2]")
  mode <- args[1]; output <- normalizePath(args[2],mustWork=FALSE)
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  Rcpp::sourceCpp(file.path(experiment_dir,"objective.cpp"),cacheDir=file.path(output,"cpp-cache"))
  if(mode=="checks") run_checks(output)
  else if(mode=="run") {
    sizes <- if(length(args)>=3) as.integer(strsplit(args[3],",")[[1]]) else 250L
    reps <- if(length(args)>=4) as.integer(args[4]) else 20L
    workers <- if(length(args)>=5) as.integer(args[5]) else 2L
    jobs <- expand.grid(n=sizes,replicate=seq_len(reps))
    began <- Sys.time()
    status <- parallel::mclapply(seq_len(nrow(jobs)),function(i) tryCatch(
      run_replicate(jobs$n[i],jobs$replicate[i],output),error=function(e) {
        message(sprintf("FAILED n=%d r=%d: %s",jobs$n[i],jobs$replicate[i],conditionMessage(e)))
        structure(conditionMessage(e),class="experiment_failure")
      }),mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
    saveRDS(list(jobs=jobs,status=status,started=began,finished=Sys.time(),workers=workers),
      file.path(output,sprintf("run-status-%s.rds",format(began,"%Y%m%d-%H%M%S"))))
    if(any(vapply(status,inherits,logical(1),"experiment_failure"))) stop("One or more replicates failed; inspect run status")
  } else stop("Unknown mode")
}
