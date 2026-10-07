# Presentation and summaries for the higher-dimensional MDS vignette.
# Uses the shared ivue selector; does not fit layouts or change saved scores.
source(system.file('scripts','metric-mds-comparison.R',package='grip'),local=TRUE)

highdim.available <- function() {
  isTRUE(getOption('grip.vignette.interactive',TRUE)) &&
    requireNamespace('ivue',quietly=TRUE) && requireNamespace('rgl',quietly=TRUE) &&
    all(c('plot3D.groups','layer3D.callback','layer3D.axes','camera.zup') %in%
          getNamespaceExports('ivue')) &&
    'limits' %in% names(formals(ivue::plot3D.plain))
}

highdim.view <- function(case) {
  if (!highdim.available()) return(knitr::asis_output(paste(
    '\n\nInteractive view unavailable: install ivue with plot3D.groups(),',
    'layer3D.callback(), layer3D.axes(), camera.zup(), and limits support,',
    'plus rgl, then rebuild with options(grip.vignette.interactive=TRUE).',
    'The numerical tables remain available. No static substitute is drawn.\n\n')))
  z <- case$coords
  labels <- c(direct_3='Direct 3D MDS',classical_3='Classical MDS',
              classical_refined_3='Classical MDS + 3D refinement')
  for (d in c(4L,6L,8L,10L,15L,20L)) {
    labels[paste0('pca_',d)] <- paste0(d,'D MDS + PCA')
    labels[paste0('refined_',d)] <- paste0(d,'D MDS + PCA + 3D refinement')
  }
  z <- z[names(labels)]
  components <- split(seq_len(case$n),case$components)
  nontrivial <- Filter(function(x) length(x)>1L,components)
  stopifnot(length(nontrivial)==1L)
  main <- nontrivial[[1L]]; isolates <- setdiff(seq_len(case$n),main)
  target <- if(is.null(case$reference)) z[['direct_3']] else case$reference
  target <- sweep(target,2,colMeans(target[main,,drop=FALSE]))
  z <- lapply(z,function(x) {
    x[main,] <- comparison_align(x[main,,drop=FALSE],target[main,,drop=FALSE]); x
  })
  # Isolate locations have no target distances to the main component. Put all
  # of them in the same display-only row; no vertex is dropped or scored here.
  if (length(isolates)) {
    bounds <- apply(do.call(rbind,lapply(z,function(x) x[main,,drop=FALSE])),2,range)
    spacing <- max(bounds[2,]-bounds[1,])*.04
    position <- cbind(bounds[2,1]+2*spacing,
                      bounds[1,2]+seq_along(isolates)*spacing,0)
    z <- lapply(z,function(x) {x[isolates,] <- position; x})
  }
  # Dense meshes can hide the vertices. The deterministic display subset has
  # no effect on the full-distance fits or exact stress tables.
  edges <- if(grepl('ambient$',case$id)) matrix(integer(),0,2) else case$edges
  if(nrow(edges)>2500L) edges <- edges[unique(round(seq(1,nrow(edges),length.out=2500L))),,drop=FALSE]
  colors <- setNames(ifelse(grepl('refined_',names(z)),'#009E73',
                    ifelse(grepl('pca_',names(z)),'#D66A19','#1769AA')),unname(labels))
  fits <- lapply(z,function(x) list(coords=x))
  comparison_view(list(X=target,n=case$n,dimension=3L,draw_edges=edges,triangles=NULL),
    fits,series_labels=labels,series_colors=c(Reference='#888888',colors),
    reference=!is.null(case$reference),surface.alignments=fits,legend=FALSE,
    width=900,height=500,camera=ivue::camera.zup(elevation=22,turn=-125,zoom=.45),
    description=paste(case$label,'— all',case$n,'vertices; display start 17.',
      nrow(edges),'edges drawn. Rigid matched-vertex alignment; no resizing.',
      if(length(isolates)) paste(length(isolates),'isolates placed in a separate row.') else '',
      if(!is.null(case$reference)) 'Gray points: generating cloud, always visible.' else ''))
}

highdim.routes <- function(case,dimension=10L) {
  keys <- c('direct','classical','classical_refined','high','pca','refined')
  labels <- c('Direct 3D','Classical 3D','Classical + 3D refinement',
              paste0(dimension,'D before projection'),paste0(dimension,'D + PCA'),
              paste0(dimension,'D + PCA + 3D refinement'))
  rows <- lapply(seq_along(keys),function(i) {
    s <- case$scores[case$scores$route==keys[i] &
      case$scores$dimension==if(keys[i] %in% c('high','pca','refined')) dimension else 3L,]
    stopifnot(nrow(s)==if(keys[i]=='classical') 1L else 3L)
    data.frame(Route=labels[i],Starts=nrow(s),
      `Mean stress`=mean(s$normalized_stress),`Minimum stress`=min(s$normalized_stress),
      `Maximum stress`=max(s$normalized_stress),`Mean route seconds`=mean(s$seconds),check.names=FALSE)
  })
  do.call(rbind,rows)
}

highdim.dimensions <- function(case) {
  do.call(rbind,lapply(sort(unique(case$matched_time$dimension)),function(d) {
    t <- case$matched_time[case$matched_time$dimension==d,]
    p <- case$scores[case$scores$route=='pca' & case$scores$dimension==d,]
    gain <- 100*(t$direct_stress-t$route_normalized_stress)/t$direct_stress
    data.frame(Dimension=d,`PCA stress`=mean(p$normalized_stress),
      `Refined stress`=mean(t$route_normalized_stress),
      `Time-matched direct stress`=mean(t$direct_stress),
      `Mean route seconds`=mean(t$route_seconds),
      `Mean paired gain (%)`=mean(gain),`Gain range (%)`=sprintf('%.4g to %.4g',min(gain),max(gain)),
      check.names=FALSE)
  }))
}

# Retain near-zero nonzero values instead of rounding them to exact zero.
highdim.table <- function(x) {
  for (name in names(x)) {
    if (grepl('stress',name,ignore.case=TRUE)) x[[name]] <- sprintf('%.7g',x[[name]])
    else if (grepl('seconds',name)) x[[name]] <- sprintf('%.3f',x[[name]])
    else if (grepl('Mean paired gain',name)) x[[name]] <- sprintf('%.4g',x[[name]])
  }
  knitr::kable(x,row.names=FALSE)
}
