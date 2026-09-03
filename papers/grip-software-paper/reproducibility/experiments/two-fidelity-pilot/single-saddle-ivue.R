#!/usr/bin/env Rscript
# Source this file to test the stages individually; sourcing does not run them.
# See README-single-saddle.md for R-console and command-line examples.
# Defaults reproduce the sampling and primary fitting protocol for cloud 5 of
# Figure 7. Smaller n and smaller reference grids are useful for quick tests.

.single.saddle.source <- local({
  frame.files <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
  file <- if (length(frame.files)) tail(frame.files, 1)[[1]] else {
    sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
  }
  dirname(normalizePath(file, mustWork=TRUE))
})

single.saddle.config <- function(n=1000L, cloud=5L, ks=3:min(80L, n-1L),
    reference.sources=min(128L, n), reference.grids=c(41L, 81L),
    fine.grid=161L, fine.sources=min(16L, reference.sources),
    max.iter=200L, audit.iter=1000L, fit.ks=NULL, plot.ks=NULL,
    edge.alpha=.03, align.display=TRUE, out=NULL, python=NULL,
    package.library=Sys.getenv("GRIP_RJOURNAL_PACKAGE_LIBRARY"),
    ivue.source=path.expand("~/current_projects/ivue")) {
  integer.value <- function(x, lower, upper=Inf) {
    is.numeric(x) && length(x)==1L && is.finite(x) && x==floor(x) && x>=lower && x<=upper
  }
  stopifnot(integer.value(n, 8, 100000), integer.value(cloud, 1, 1000000),
    length(ks)>0, is.numeric(ks), all(is.finite(ks)), all(ks==floor(ks)),
    all(ks>=1 & ks<n), integer.value(reference.sources, 2, n),
    length(reference.grids)==2L, all(is.finite(reference.grids)),
    all(reference.grids==floor(reference.grids)), all(reference.grids>=3),
    reference.grids[1]<reference.grids[2], integer.value(fine.grid, reference.grids[2]+1),
    integer.value(fine.sources, 2, reference.sources), integer.value(max.iter, 1),
    integer.value(audit.iter, 0), length(edge.alpha)==1L, is.finite(edge.alpha),
    edge.alpha>=0, edge.alpha<=1, is.logical(align.display), length(align.display)==1L,
    !is.na(align.display))
  ks <- sort(unique(as.integer(ks)))
  for (v in list(fit.ks, plot.ks)) if (!is.null(v)) stopifnot(length(v)>0, all(v %in% ks))
  paper <- normalizePath(file.path(.single.saddle.source, "../../.."))
  if (is.null(out)) out <- file.path(paper, "build", sprintf("single-saddle-n%d-cloud%d", n, cloud))
  if (is.null(python) || !nzchar(python)) {
    python <- Sys.getenv("GRIP_SADDLE_PYTHON")
    local.python <- file.path(paper, "build/two-fidelity-pilot/venv/bin/python")
    if (!nzchar(python)) python <- if (file.exists(local.python)) local.python else Sys.which("python3")
  }
  list(n=as.integer(n), cloud=as.integer(cloud), ks=ks,
    reference.sources=as.integer(reference.sources), reference.grids=as.integer(reference.grids),
    fine.grid=as.integer(fine.grid), fine.sources=as.integer(fine.sources),
    max.iter=as.integer(max.iter), audit.iter=as.integer(audit.iter),
    fit.ks=if (is.null(fit.ks)) NULL else sort(unique(as.integer(fit.ks))),
    plot.ks=plot.ks, edge.alpha=edge.alpha, align.display=align.display,
    out=path.expand(out), python=python, package.library=package.library, ivue.source=ivue.source)
}

setup.single.saddle <- function(config) {
  if (nzchar(config$package.library)) .libPaths(c(config$package.library, .libPaths()))
  for (p in c("grip", "dgraphs", "igraph", "Rcpp", "htmlwidgets", "htmltools")) {
    if (!requireNamespace(p, quietly=TRUE)) stop("Install required R package: ", p)
  }
  for (p in c("grip", "dgraphs")) if (utils::packageVersion(p)<"0.2.0") {
    stop(p, " >= 0.2.0 is required; set package.library or GRIP_RJOURNAL_PACKAGE_LIBRARY.")
  }
  old <- options(rgl.useNULL=TRUE); on.exit(options(old))
  if (!is.null(config$ivue.source)) {
    if (!file.exists(file.path(config$ivue.source, "DESCRIPTION"))) stop("ivue.source is not a package directory.")
    if (!requireNamespace("pkgload", quietly=TRUE)) stop("Install pkgload to load the local ivue source.")
    pkgload::load_all(config$ivue.source, quiet=TRUE, export_all=FALSE)
  }
  if (!requireNamespace("ivue", quietly=TRUE) || !requireNamespace("rgl", quietly=TRUE))
    stop("ivue and rgl are required for interactive views.")
  if (!nzchar(config$python) || !file.exists(config$python)) stop("Set python to a Python executable with requirements.txt installed.")
  check <- system2(config$python, c("-c", shQuote("import numpy, scipy, pygeodesic")), stdout=TRUE, stderr=TRUE)
  if (!is.null(attr(check, "status"))) stop("Python dependencies unavailable:\n", paste(check, collapse="\n"))
  scientific <- config[setdiff(names(config), c("out", "python", "package.library", "ivue.source",
                                               "plot.ks", "edge.alpha", "align.display"))]
  scientific$packages <- sapply(c("grip", "dgraphs", "igraph"), function(p) as.character(utils::packageVersion(p)))
  scientific$R <- R.version.string
  scientific$python.packages <- system2(config$python, c("-c", shQuote(
    "from importlib.metadata import version; print(';'.join(p+'='+version(p) for p in ('numpy','scipy','pygeodesic')))")), stdout=TRUE)
  scientific$source.md5 <- tools::md5sum(file.path(.single.saddle.source,
    c("single-saddle-ivue.R", "single-saddle-reference.py", "surface-reference.py", "path-lengths.cpp")))
  manifest <- file.path(config$out, "config.rds")
  if (file.exists(manifest)) {
    if (!identical(scientific, readRDS(manifest))) stop("Cached protocol differs. Choose a new output directory; existing results were not changed.")
  } else {
    if (dir.exists(config$out) && length(list.files(config$out, all.files=TRUE, no..=TRUE)))
      stop("Output directory is nonempty but has no protocol record. Choose a new directory.")
    dir.create(config$out, recursive=TRUE, showWarnings=FALSE)
    saveRDS(scientific, manifest)
  }
  config$out <- normalizePath(config$out)
  capture.output(sessionInfo(), file=file.path(config$out, "session-info.txt"))
  config
}

.ss.relative <- function(y, d) sqrt(sum((y-d)^2)/sum(d*d))
.ss.profile <- function(y, d) {
  scale <- sum(y*d)/sum(d*d)
  if (!is.finite(scale) || scale<=0) stop("Degenerate fitted scale.")
  .ss.relative(y, scale*d)
}
.ss.lengths <- function(z, e) sqrt(rowSums((z[e[,1],,drop=FALSE]-z[e[,2],,drop=FALSE])^2))
.ss.align <- function(z, x) {
  zz <- sweep(z, 2, colMeans(z), "-"); xx <- sweep(x, 2, colMeans(x), "-")
  f <- svd(crossprod(zz, xx))
  sweep(sum(f$d)/sum(zz*zz)*zz%*%f$u%*%t(f$v), 2, colMeans(x), "+")
}
.ss.with.seed <- function(seed, expr) {
  had <- exists(".Random.seed", .GlobalEnv, inherits=FALSE)
  if (had) previous <- get(".Random.seed", .GlobalEnv)
  kind <- RNGkind()
  on.exit({ do.call(RNGkind, as.list(kind)); if (had) assign(".Random.seed", previous, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits=FALSE)) rm(".Random.seed", envir=.GlobalEnv) })
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed); force(expr)
}

# Stage 1: surface-area-uniform sampling; identical to cloud 5 when n=1000.
sample.single.saddle <- function(config) {
  file <- file.path(config$out, sprintf("cloud-%02d.rds", config$cloud))
  if (file.exists(file)) return(readRDS(file))
  x <- .ss.with.seed(2211000L+config$cloud, {
    uv <- matrix(numeric(), ncol=2)
    while (nrow(uv)<config$n) {
      a <- matrix(runif(4000, -1, 1), ncol=2)
      accept <- runif(nrow(a)) < sqrt(1+2.56*rowSums(a*a))/sqrt(6.12)
      uv <- rbind(uv, a[accept,,drop=FALSE])
    }
    uv <- uv[seq_len(config$n),,drop=FALSE]
    cbind(x=uv[,1], y=uv[,2], z=.8*(uv[,1]^2-uv[,2]^2))
  })
  sources <- .ss.with.seed(3211000L+config$cloud, sort(sample.int(config$n, config$reference.sources)))
  cloud <- list(coords=x, sources=sources, seed=2211000L+config$cloud)
  saveRDS(cloud, file)
  write.csv(x, file.path(config$out, sprintf("cloud-%02d.csv", config$cloud)), row.names=FALSE)
  write.table(sources-1L, file.path(config$out, sprintf("sources-%02d.txt", config$cloud)), row.names=FALSE, col.names=FALSE)
  cloud
}

.ss.reference.pairs <- function(ref) {
  grid <- expand.grid(row=seq_along(ref$sources), j=seq_len(ncol(ref$D)))
  i <- ref$sources[grid$row]; j <- grid$j
  keep <- i!=j & (!(j %in% ref$sources) | i<j)
  list(pairs=cbind(pmin(i,j), pmax(i,j))[keep,,drop=FALSE], d=ref$D[as.matrix(grid)][keep])
}

# Stage 2: reference on two meshes, a finer source subset, plane and BVP checks.
reference.single.saddle <- function(config) {
  dest <- file.path(config$out, "reference.rds")
  if (file.exists(dest)) return(readRDS(dest))
  status <- system2(config$python, c(shQuote(file.path(.single.saddle.source, "single-saddle-reference.py")),
    shQuote(config$out), "--cloud", config$cloud, "--grids", paste(config$reference.grids, collapse=","),
    "--sources", config$reference.sources, "--fine-grid", config$fine.grid, "--fine-sources", config$fine.sources))
  if (status!=0L) stop("Numerical reference failed; inspect its console output.")
  read.ref <- function(grid, sources) {
    a <- as.matrix(read.csv(file.path(config$out, sprintf("reference-r%02d-m%d-s%d.csv", config$cloud, grid, sources)), header=FALSE))
    list(sources=as.integer(a[,1]), D=a[,-1,drop=FALSE])
  }
  coarse <- read.ref(config$reference.grids[1], config$reference.sources)
  main <- read.ref(config$reference.grids[2], config$reference.sources)
  fine <- read.ref(config$fine.grid, config$fine.sources)
  pairs <- .ss.reference.pairs(main)
  stopifnot(length(pairs$d)==config$reference.sources*(config$n-1)-choose(config$reference.sources, 2),
    all(is.finite(pairs$d)), all(pairs$d>0))
  bvp <- read.csv(file.path(config$out, "smooth-checks.csv"))
  bvp.ref <- main$D[cbind(match(bvp$i, main$sources), bvp$j)]
  validation <- data.frame(coarse_to_main=.ss.relative(coarse$D, main$D),
    main_to_fine=.ss.relative(main$D[match(fine$sources, main$sources),,drop=FALSE], fine$D),
    main_to_smooth=.ss.relative(bvp.ref, bvp$distance), smooth_pairs=nrow(bvp))
  result <- list(main=main, coarse=coarse, fine=fine, pairs=pairs, validation=validation)
  saveRDS(result, dest); write.csv(validation, file.path(config$out, "reference-checks.csv"), row.names=FALSE)
  result
}

# Stage 3: same sample at every k; store one graph per file, not 78 dense D_Gs.
graphs.single.saddle <- function(config, cloud, reference) {
  dir.create(file.path(config$out, "graphs"), showWarnings=FALSE)
  curves <- lapply(config$ks, function(k) {
    file <- file.path(config$out, "graphs", sprintf("k-%03d.rds", k))
    if (file.exists(file)) return(readRDS(file)$summary)
    tick <- proc.time()[[3]]
    g <- dgraphs::create.sknn.graph(cloud$coords, k=k, neighbor.method="exact",
      connect.components=TRUE, connect.method="component.mst", prune.method="none",
      prune.edges=FALSE, edge.weight="distance")
    ig <- igraph::graph_from_edgelist(g$edge_matrix, directed=FALSE)
    stopifnot(igraph::vcount(ig)==config$n, igraph::is_connected(ig),
      g$n_mst_edges_added==g$n_components_before-1L,
      max(abs(.ss.lengths(cloud$coords,g$edge_matrix)-g$edge_weight))<1e-12)
    D <- igraph::distances(ig, v=reference$main$sources, weights=g$edge_weight)
    stopifnot(all(is.finite(D)))
    y <- .ss.reference.pairs(list(sources=reference$main$sources, D=D))$d
    yf <- .ss.reference.pairs(list(sources=reference$fine$sources,
      D=D[match(reference$fine$sources, reference$main$sources),,drop=FALSE]))$d
    d <- reference$pairs$d
    row <- data.frame(k=k, edges=nrow(g$edge_matrix), components_before=g$n_components_before,
      bridges=g$n_mst_edges_added, xg_error=.ss.relative(y,d),
      coarse_error=.ss.relative(y,.ss.reference.pairs(reference$coarse)$d),
      fine_subset_error=.ss.relative(yf,.ss.reference.pairs(reference$fine)$d),
      fitted_scale=sum(y*d)/sum(d*d), mean_relative_bias=mean(y/d-1),
      fraction_shorter=mean(y<d), seconds=proc.time()[[3]]-tick)
    saveRDS(list(k=k, edges=g$edge_matrix, weights=g$edge_weight,
      mst_edges=g$mst_edge_matrix, summary=row), file)
    message(sprintf("k=%d: %d edges, %d MST bridges, X->G RMSE %.4f%%", k, row$edges, row$bridges, 100*row$xg_error))
    row
  })
  curve <- do.call(rbind, curves)
  best <- curve$k[which.min(curve$xg_error)]
  boundary <- best %in% range(config$ks)
  if (boundary) warning("Minimum is on the tested k boundary; expand ks in a NEW output directory before interpreting it as an interior minimum.")
  write.csv(curve, file.path(config$out, "calibration.csv"), row.names=FALSE)
  list(curve=curve, selected.k=best, boundary=boundary)
}

# Stage 4: full fixed routes, 3D MDS, density-continuation edge-KK, independent scores.
fit.single.saddle <- function(config, cloud, reference, calibration) {
  cpp <- new.env(parent=globalenv())
  Rcpp::sourceCpp(file.path(.single.saddle.source, "path-lengths.cpp"),
    cacheDir=file.path(config$out, "cpp-cache"), env=cpp)
  fit.ks <- if (is.null(config$fit.ks)) calibration$selected.k else config$fit.ks
  fits <- lapply(fit.ks, function(k) {
    dest <- file.path(config$out, sprintf("fit-k-%03d.rds", k))
    if (file.exists(dest)) return(readRDS(dest))
    g <- readRDS(file.path(config$out, "graphs", sprintf("k-%03d.rds", k)))
    tick <- proc.time()[[3]]
    p <- grip::prepare.geodesic.kk(g$edges, n=config$n, edge_weights=g$weights, tie_mode="single")
    ig <- igraph::graph_from_edgelist(g$edges, directed=FALSE)
    strict <- igraph::distances(ig, weights=g$weights)[p$pair_matrix]
    dg <- p$pair_graph_distance
    stopifnot(nrow(p$pair_matrix)==choose(config$n,2), all(p$flat_edge_coeff==1),
      max(abs(dg-strict)/pmax(1,dg,strict))<=sqrt(.Machine$double.eps))
    prepare.seconds <- proc.time()[[3]]-tick
    key <- function(e) paste(e[,1],e[,2],sep=":")
    ref.index <- match(key(reference$pairs$pairs),key(p$pair_matrix)); stopifnot(!anyNA(ref.index))
    tick <- proc.time()[[3]]
    mds <- grip::metric.mds(prepared=p, dim=3, diagnostics=FALSE)
    mds.seconds <- proc.time()[[3]]-tick
    tick <- proc.time()[[3]]
    kk <- grip::edge.kk(coords=mds$coords, prepared=p, dim=3, max_iter=config$max.iter,
      stiffness_method="density", density_mix_schedule=c(0,.25,.5,.75,1),
      scale_mode="profiled", edge_length_epsilon=0, diagnostics=FALSE, return_trace=TRUE, seed=5211000+config$cloud)
    kk.seconds <- proc.time()[[3]]-tick
    candidates <- list("Original saddle"=cloud$coords, "Metric MDS"=mds$coords, "MDS + edge-KK"=kk$coords)
    score <- function(z, method) {
      stopifnot(nrow(z)==config$n, ncol(z)==3, all(is.finite(z)))
      edge <- .ss.lengths(z,p$edges); chord <- .ss.lengths(z,p$pair_matrix)
      path <- cpp$pilot_path_lengths(z,p$flat_pair_edge_offsets,p$flat_edge_u,p$flat_edge_v)
      check <- grip::score.gmds(z,prepared=p,edge_length_epsilon=0)
      error <- max(abs(c(.ss.profile(path,dg)-check$gmds.stress,
        .ss.profile(edge,p$edge_targets)-check$edge.rel.rmse)))
      ids <- unique(as.integer(seq(1,length(path),length.out=min(101,length(path)))))
      direct <- vapply(ids, function(i) sum(.ss.lengths(z,p$path_edges[[i]])), numeric(1))
      direct.error <- max(abs(path[ids]-direct))
      stopifnot(error<1e-10, direct.error<1e-9)
      edge.scale <- sum(edge*p$edge_targets)/sum(p$edge_targets^2)
      chord.scale <- sum(chord*dg)/sum(dg*dg)
      data.frame(k=k,method=method,path_rel=.ss.profile(path,dg),edge_rel=.ss.profile(edge,p$edge_targets),
        stress1=sqrt(sum((chord-chord.scale*dg)^2)/sum(chord^2)),
        xg_error=.ss.relative(strict[ref.index],reference$pairs$d),
        xz_path_error=.ss.relative(path[ref.index]/edge.scale,reference$pairs$d),
        procrustes=sqrt(sum((.ss.align(z,cloud$coords)-cloud$coords)^2)/
          sum(sweep(cloud$coords,2,colMeans(cloud$coords),"-")^2)),
        package_score_difference=error,independent_path_difference=direct.error)
    }
    scores <- do.call(rbind,lapply(names(candidates),function(name) score(candidates[[name]],name)))
    stopifnot(scores$edge_rel[1]<1e-12, scores$path_rel[1]<1e-7)
    audit <- NULL
    if (config$audit.iter>0L) {
      extra <- grip::edge.kk(coords=kk$coords,prepared=p,dim=3,max_iter=config$audit.iter,
        stiffness_method="density",density_mix_schedule=1,scale_mode="profiled",edge_length_epsilon=0,
        diagnostics=FALSE,return_trace=TRUE,seed=6211000+config$cloud)
      audit <- list(coords=extra$coords,scores=score(extra$coords,"Additional budget"),trace=extra$trace,metadata=extra$metadata)
    }
    nearest <- function(t) which.min(rowSums(sweep(cloud$coords[,1:2],2,t,"-")^2))
    ends <- sort(c(nearest(c(-.8,0)),nearest(c(.8,0))))
    if (ends[1]==ends[2]) ends <- p$pair_matrix[which.max(dg),]
    id <- which(p$pair_matrix[,1]==ends[1] & p$pair_matrix[,2]==ends[2])
    fit <- list(k=k,candidates=candidates,scores=scores,audit=audit,ends=ends,
      route=p$path_edges[[id]],edge.trace=kk$trace,edge.metadata=kk$metadata,
      timing=c(preparation=prepare.seconds,mds=mds.seconds,edge_kk=kk.seconds),
      pairs=nrow(p$pair_matrix),strict_distance_max_difference=max(abs(dg-strict)))
    saveRDS(fit,dest); message("Fitted and checked k=",k); fit
  })
  names(fits) <- as.character(fit.ks)
  write.csv(do.call(rbind,lapply(fits, `[[`, "scores")),file.path(config$out,"layout-scores.csv"),row.names=FALSE)
  fits
}

# Stage 5: native ivue widgets; coordinates are never projected to 2D here.
view.single.saddle <- function(config, cloud, calibration, fits) {
  old <- options(rgl.useNULL=TRUE); on.exit(options(old))
  output <- file.path(config$out,"views"); dir.create(output,showWarnings=FALSE)
  colors <- ivue::color.scale.cont(cloud$coords[,1], limits=c(-1,1),
    palette=c("#173D65","#86AFC4","#D9B18B","#8E4921"))
  scene <- list(values=cloud$coords[,1],scale=colors,legend.title="Original x coordinate",
    point.size=4,axes=TRUE,xlab="x",ylab="y",zlab="z",aspect="equal",height=650L,
    camera=list(theta=35,phi=22,fov=0,zoom=.75))
  # RGBA colors alone may lose alpha on rgl line primitives. Replace only the
  # built-in graph lines with an explicitly translucent layer on ivue's private
  # device, BEFORE adding the opaque route/chord. No device is opened or closed.
  edge.opacity <- function(edges) ivue::layer3D.callback(function(ctx, edges, opacity) {
    shapes <- rgl::ids3d(type="shapes")
    for (id in shapes$id[shapes$type=="lines"]) rgl::pop3d(id=id)
    if (opacity>0) rgl::segments3d(ctx$X[as.vector(t(edges)),,drop=FALSE],
      col="gray45",alpha=opacity,lwd=1,lit=FALSE)
  }, args=list(edges=edges,opacity=config$edge.alpha))
  files <- character(); widgets <- list()
  save.view <- function(widget, slug, title, description) {
    htmlwidgets::saveWidget(widget,file.path(output,paste0(slug,".html")),selfcontained=TRUE,title=title)
    stopifnot(nrow(attr(widget,"ivue")$X)==config$n)
    files[slug] <<- paste0(slug,".html"); widgets[[slug]] <<- widget
    list(title=title,description=description,file=paste0(slug,".html"))
  }
  entries <- list(save.view(do.call(ivue::plot3D.cont,c(list(X=cloud$coords),scene)),
    "original-saddle","Original saddle — points only","Uniform surface-area sample of z = 0.8(x² − y²)."))
  plot.ks <- if (is.null(config$plot.ks)) sort(unique(c(head(config$ks,1),
    intersect(c(10L,20L,40L),config$ks),calibration$selected.k,tail(config$ks,1)))) else config$plot.ks
  for (k in plot.ks) {
    g <- readRDS(file.path(config$out,"graphs",sprintf("k-%03d.rds",k)))
    layers <- c(list(edge.opacity(g$edges)), if (nrow(g$mst_edges)) list(ivue::layer3D.edges(g$mst_edges,col="#AD4C15",width=3)) else list())
    ig <- igraph::graph_from_edgelist(g$edges,directed=FALSE)
    # The endpoints stay the same as in the first fitted graph, but the route
    # is computed on THIS k graph; it is never recomputed after embedding.
    ends <- fits[[1]]$ends
    route <- as.integer(igraph::shortest_paths(ig,ends[1],ends[2],weights=g$weights)$vpath[[1]])
    layers <- c(layers,list(ivue::layer3D.path(route,col="#1E5C89",width=4),
      ivue::layer3D.edges(matrix(ends,nrow=1),col="#B26026",width=2)))
    graph <- data.frame(from=g$edges[,1],to=g$edges[,2],weight=g$weights)
    w <- do.call(ivue::plot3D.graph,c(list(graph=graph,vertices=seq_len(config$n),X=cloud$coords,
      weight.type="distance",edge.col=grDevices::adjustcolor("gray45",alpha.f=config$edge.alpha),edge.width=1,layers=layers),scene))
    stopifnot(nrow(attr(w,"ivue")$graph$edges)==nrow(g$edges))
    entries[[length(entries)+1]] <- save.view(w,sprintf("graph-k-%03d",k),sprintf("Original coordinates: symmetric kNN, k = %d",k),
      sprintf("%d edges; %d MST bridges (orange); X→G RMSE %.4f%%. Blue: retained route; brown: endpoint chord.",
        nrow(g$edges),g$summary$bridges,100*g$summary$xg_error))
  }
  for (fit in fits) {
    g <- readRDS(file.path(config$out,"graphs",sprintf("k-%03d.rds",fit$k)))
    for (name in names(fit$candidates)) for (edges in c(TRUE,FALSE)) {
      z <- fit$candidates[[name]]
      if (config$align.display) z <- .ss.align(z,cloud$coords)
      layers <- list(ivue::layer3D.edges(fit$route,col="#1E5C89",width=4),
        ivue::layer3D.edges(matrix(fit$ends,nrow=1),col="#B26026",width=2))
      w <- if (edges) do.call(ivue::plot3D.graph,c(list(
        graph=data.frame(from=g$edges[,1],to=g$edges[,2],weight=g$weights),vertices=seq_len(config$n),X=z,
        weight.type="distance",edge.col="gray45",edge.width=1,layers=c(list(edge.opacity(g$edges)),layers)),scene)) else
        do.call(ivue::plot3D.cont,c(list(X=z,layers=layers),scene))
      slug <- sprintf("k-%03d-%s-%s",fit$k,tolower(gsub("[^[:alnum:]]+","-",name)),if(edges) "graph" else "points")
      s <- fit$scores[fit$scores$method==name,]
      description <- sprintf("Path RMSE %.4f%%; edge RMSE %.4f%%; Stress-1 %.4f%%. All %s pairs; scores precede display alignment. %s Background edges: %s.",
        100*s$path_rel,100*s$edge_rel,100*s$stress1,format(fit$pairs,big.mark=","),
        if(config$align.display) "Similarity-aligned for display only." else "Unaligned coordinates.",if(edges) "shown" else "hidden")
      entries[[length(entries)+1]] <- save.view(w,slug,sprintf("k = %d: %s (%s)",fit$k,name,if(edges) "graph" else "points"),description)
    }
  }
  grDevices::pdf(file.path(output,"calibration.pdf"),width=7,height=4.5)
  plot(calibration$curve$k,100*calibration$curve$xg_error,type="b",pch=16,
    xlab="Number of neighbors k",ylab="Surface-to-graph relative RMSE (%)")
  abline(v=calibration$selected.k,lty=2,col="#B26026"); grDevices::dev.off()
  tags <- htmltools::tags
  page <- tags$html(tags$head(tags$title("Single-saddle ivue exploration"),tags$meta(charset="utf-8")),
    tags$body(style="font-family:system-ui;max-width:960px;margin:40px auto;line-height:1.5;padding:0 20px",
      tags$h1(sprintf("One saddle sample: n = %s, cloud %d",format(config$n,big.mark=","),config$cloud)),
      tags$p(sprintf("Calibrated over k = %d…%d; minimum at k = %d.%s",min(config$ks),max(config$ks),calibration$selected.k,
        if(calibration$boundary) " WARNING: the minimum is at the tested boundary." else "")),
      tags$p("Drag a view to rotate; scroll to zoom. All scenes are three-dimensional, with equal data-unit axes. Colors always identify the original x coordinate. Each widget auto-fits its own bounds; compare shape, not apparent size between windows."),
      tags$p("Graph edges use input Euclidean lengths. Blue shows a fixed graph route; brown joins its endpoints. Edge-KK is a finite-budget edge objective, not direct optimization of fixed-path GMDS. Reference and source-subset checks are numerical sensitivity checks, not certified error bounds."),
      tags$p(tags$a(href="calibration.pdf","Calibration curve")," · ",tags$a(href="../calibration.csv","Calibration data")," · ",
        tags$a(href="../layout-scores.csv","Layout scores")," · ",tags$a(href="../reference-checks.csv","Reference checks")),
      lapply(entries,function(e) tags$section(tags$h2(tags$a(href=e$file,e$title)),tags$p(e$description)))))
  htmltools::save_html(page,file.path(output,"index.html"))
  list(index=file.path(output,"index.html"),files=files,widgets=widgets)
}

run.single.saddle <- function(config=single.saddle.config(), open=interactive()) {
  config <- setup.single.saddle(config)
  message("1/5 Sampling"); cloud <- sample.single.saddle(config)
  message("2/5 Numerical surface reference and checks"); reference <- reference.single.saddle(config)
  message("3/5 Graph sweep and calibration"); calibration <- graphs.single.saddle(config,cloud,reference)
  message("4/5 MDS, edge-KK, and independent score checks"); fits <- fit.single.saddle(config,cloud,reference,calibration)
  message("5/5 Interactive ivue views"); views <- view.single.saddle(config,cloud,calibration,fits)
  result <- list(config=config,cloud=cloud,reference=reference,calibration=calibration,fits=fits,views=views)
  saveRDS(result[setdiff(names(result),"views")],file.path(config$out,"experiment.rds"))
  message("Open ",views$index)
  if (open) utils::browseURL(views$index)
  invisible(result)
}

if (sys.nframe()==0L) {
  args <- commandArgs(TRUE)
  if ("--help" %in% args) {
    cat("Usage: Rscript single-saddle-ivue.R [n] [cloud] [output-directory]\n",
        "For all controls, source the script and call single.saddle.config() / run.single.saddle().\n")
  } else {
    if (length(args)>3L) stop("Use --help or the R-console interface for additional controls.")
    config <- single.saddle.config(n=if(length(args)) as.integer(args[1]) else 1000L,
      cloud=if(length(args)>1) as.integer(args[2]) else 5L,
      out=if(length(args)>2) args[3] else NULL)
    run.single.saddle(config,open=FALSE)
  }
}
