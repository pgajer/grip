# Shared presentation recipes. These helpers do not fit or score layouts.
source(system.file('scripts', 'metric-mds-comparison.R', package='grip'), local=TRUE)
source(system.file('scripts', 'metric-mds-surface-alignment.R', package='grip'), local=TRUE)

documentation.interactive <- function() {
  isTRUE(getOption('grip.vignette.interactive', TRUE)) &&
    requireNamespace('ivue', quietly=TRUE) && requireNamespace('rgl', quietly=TRUE) &&
    all(c('plot3D.groups','layer3D.callback','layer3D.axes','camera.zup') %in%
          getNamespaceExports('ivue')) &&
    'limits' %in% names(formals(ivue::plot3D.plain))
}

documentation.static <- function(coordinates, edges, colors=NULL) {
  if (is.null(colors)) colors <- setNames(c('#1769AA','#D66A19','#009E73','#AA3377',
                                          '#888888')[seq_along(coordinates)], names(coordinates))
  xy <- lapply(coordinates, function(x) grip::project.3d(x, azimuth=35, elevation=22))
  bounds <- do.call(rbind, xy)
  graphics::plot(bounds, type='n', asp=1, axes=FALSE, xlab='', ylab='')
  for (i in seq_along(xy)) {
    z <- xy[[i]]; color <- colors[[i]]
    if (nrow(edges)) graphics::segments(z[edges[,1],1], z[edges[,1],2],
      z[edges[,2],1], z[edges[,2],2], col=grDevices::adjustcolor(color,.5), lwd=.6)
    graphics::points(z, pch=16, cex=if (nrow(z)>500) .25 else .45, col=color)
  }
  graphics::legend('topleft', legend=names(coordinates), col=colors, pch=16, bty='n', cex=.85)
  invisible(NULL)
}

documentation.view <- function(layouts, edges, reference=NULL, description='', aligned=FALSE,
                               common.size=FALSE) {
  stopifnot(length(layouts)>0, !is.null(names(layouts)))
  target <- if (is.null(reference)) layouts[[1]] else reference
  center <- colMeans(target)
  target <- sweep(target,2,center)
  if (!is.null(reference)) reference <- target
  moved <- if (aligned) lapply(layouts,function(z) sweep(z,2,center)) else
    lapply(layouts, comparison_align, reference=target)
  if (common.size) {
    multipliers <- vapply(moved,function(z) sqrt(sum(target^2)/sum(z^2)),0)
    stopifnot(all(is.finite(multipliers)),all(multipliers>0))
    moved <- Map(function(z,k) z*k,moved,multipliers)
    description <- paste(description,'Display multipliers:',paste(names(multipliers),
      format(multipliers,digits=3),collapse='; '),'. Scores use original coordinates.')
  }
  palette <- setNames(c('#1769AA','#D66A19','#009E73','#AA3377','#CC79A7')[seq_along(moved)], names(moved))
  if (!is.null(reference)) palette <- c(Reference='#888888',palette)
  if (!documentation.interactive()) {
    points <- c(if (!is.null(reference)) list(Reference=reference), moved)
    documentation.static(points, edges, palette)
    return(knitr::asis_output(paste0('\n\n',description,'\n\n')))
  }
  case <- list(X=target, n=nrow(target), dimension=3L, draw_edges=edges, triangles=NULL)
  fits <- lapply(moved, function(z) list(coords=z))
  comparison_view(case, fits, series_labels=setNames(names(moved),names(moved)),
    series_colors=palette, reference=!is.null(reference), description=description,
    surface.alignments=fits)
}

documentation.groups <- function(coords, edges, groups, description='') {
  groups <- as.character(groups); groups[is.na(groups)] <- 'Unclassified'
  levels <- sort(unique(groups)); colors <- setNames(grDevices::hcl.colors(length(levels),'Dark 3'),levels)
  if (!documentation.interactive()) {
    grip::plot.layout(coords, edges=edges, projection='ortho', vertex.col=colors[groups], edge.col='gray65')
    graphics::legend('topleft', legend=levels, col=colors, pch=16, bty='n', cex=.8)
    return(invisible(NULL))
  }
  widget <- ivue::plot3D.groups(coords,groups=groups,
    scale=ivue::color.scale.groups(levels,colors=colors),
    layers=list(ivue::layer3D.edges(edges,col='#80808055',width=.6)),
    axes=FALSE, aspect='equal', camera=ivue::camera.zup(),
    legend.show=TRUE, controls=FALSE, description=description, width=900,height=480)
  htmlwidgets::onRender(widget,"function(el) {
    el.querySelectorAll('.ivue-legend details').forEach(function(x){x.remove();});
    el.querySelectorAll('.ivue-legend').forEach(function(x){x.style.width='max-content';});
  }")
}

documentation.full.view <- function(case,fits,weighting=FALSE) {
  surface <- if (grepl('^saddle_graph-',case$id)) comparison.surface.saved(case,fits) else NULL
  if (documentation.interactive()) return(comparison_view(case,fits,weighting=weighting,
    surface.alignments=surface$registrations))
  # Match the initial interactive selection, with the same reference and transform.
  key <- if (weighting) 'sgd_uniform' else 'sgd'
  z <- if (is.null(surface)) comparison_align(fits[[key]]$coords,case$X) else
    surface$registrations[[key]]$coords
  documentation.static(list(Reference=case$X,SGD=z),case$draw_edges,c('#888888','#1769AA'))
}

documentation.trace <- function(trace, edges) {
  # Planar graph: use an SVG with the same bounds at every recorded stage.
  frames <- trace$frames
  all <- do.call(rbind,frames); all <- all[stats::complete.cases(all),,drop=FALSE]
  center <- (apply(all[,1:2],2,min)+apply(all[,1:2],2,max))/2
  extent <- max(apply(all[,1:2],2,function(x) diff(range(x))))
  transform <- function(z) sweep(sweep(z[,1:2,drop=FALSE],2,center),2,c(1,-1),'*') * 440/max(extent,1e-8) + 250
  groups <- lapply(seq_along(frames),function(i) {
    z <- transform(frames[[i]]); active <- stats::complete.cases(z)
    e <- edges[active[edges[,1]] & active[edges[,2]],,drop=FALSE]
    segments <- if(nrow(e)) sprintf('<line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f"/>',
      z[e[,1],1],z[e[,1],2],z[e[,2],1],z[e[,2],2]) else character()
    dots <- sprintf('<circle cx="%.3f" cy="%.3f" r="1.7"/>',z[active,1],z[active,2])
    sprintf('<g class="trace-frame" style="display:%s"><g stroke="#547a99" stroke-width=".8">%s</g><g fill="#1769aa">%s</g></g>',
      if(i==length(frames)) 'inline' else 'none',paste(segments,collapse=''),paste(dots,collapse=''))
  })
  labels <- vapply(seq_along(frames),function(i) paste('Frame',i,'of',length(frames),
    '—',trace$meta$phase[i],'—',sum(stats::complete.cases(frames[[i]])),'active vertices'),'')
  htmltools::tags$div(class='grip-trace',
    htmltools::HTML(paste0('<svg viewBox="0 0 500 500" role="img" aria-label="Recorded recursive triangle construction" style="width:100%;max-height:500px">',paste(groups,collapse=''),'</svg>')),
    htmltools::tags$label('Recorded stage: ',htmltools::tags$input(type='range',min=1,
      max=length(frames),value=length(frames),step=1,
      oninput="var p=this.closest('.grip-trace');var i=Number(this.value)-1;p.querySelectorAll('.trace-frame').forEach(function(g,j){g.style.display=j===i?'inline':'none';});p.querySelectorAll('.trace-label').forEach(function(g,j){g.style.display=j===i?'inline':'none';});")),
    lapply(seq_along(labels),function(i) htmltools::tags$span(class='trace-label',
      style=paste0('display:',if(i==length(frames)) 'inline' else 'none'),labels[i])))
}
