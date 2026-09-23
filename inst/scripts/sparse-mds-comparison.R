# Public-API recipes for the sparse-SGD comparison. Source the metric comparison
# recipe first for alignment and ivue display; sourcing never runs a benchmark.
sparse_comparison_cases <- function() {
  cases <- comparison_cases(8L)[c('saddle_graph-64','paraboloid_graph-64','helix_graph-64')]
  env <- new.env(); utils::data('zheng.graphs',package='grip',envir=env)
  g <- env$zheng.graphs[['dwt_307']]
  cases[['dwt_307']] <- list(id='dwt_307',label='Structural graph dwt_307',n=g$n,
    dimension=3L,edges=g$edges,weights=g$edge_weights,draw_edges=g$edges,
    X=NULL,triangles=NULL,target='Graph shortest paths')
  side <- 64L; uv <- as.matrix(expand.grid(seq(-1,1,length.out=side),seq(-1,1,length.out=side)))
  X <- cbind(uv,.8*(uv[,1]^2-uv[,2]^2)); idx <- matrix(seq_len(nrow(X)),side,side)
  a <- as.vector(idx[-side,-side]); b <- a+1L; c <- a+side; d <- c+1L
  E <- unique(rbind(cbind(a,b),cbind(a,c),cbind(b,d),cbind(c,d),cbind(a,d),cbind(b,c)))
  w <- sqrt(rowSums((X[E[,1],]-X[E[,2],])^2))
  cases[['saddle-4096']] <- list(id='saddle-4096',label='Larger saddle graph',n=nrow(X),
    dimension=3L,edges=E,weights=w,draw_edges=E,X=X,triangles=NULL,target='Graph shortest paths')
  cases
}

sparse_comparison_fit <- function(case, method, seed, epochs=100L) {
  start <- withr::with_seed(1000L+seed,matrix(runif(case$n*3),case$n,3))
  began <- proc.time()[['elapsed']]
  if (method=='full') {
    fit <- grip::metric.mds(edges=case$edges,n=case$n,edge_weights=case$weights,
      dim=3,init=start,seed=seed,max_iter=epochs,backend='sgd',
      pair_weights='inverse_squared',diagnostics=FALSE)
    prep <- NA_real_; fitting <- NA_real_
  } else if (method=='grip') {
    fit <- grip::grip(edges=case$edges,n=case$n,edge_weights=case$weights,dim=3,
      metric='edge_length',length_normalization='none',seed=seed,disconnected='error')
    prep <- NA_real_; fitting <- NA_real_
  } else {
    h <- as.integer(sub('sparse','',method))
    fit <- grip::metric.mds(approximation="sparse",edges=case$edges,n=case$n,edge_weights=case$weights,
      dim=3,init=start,seed=seed,max_iter=epochs,sparse_control=list(n_pivots=h))
    prep <- fit$metadata$preparation_seconds; fitting <- fit$metadata$fitting_seconds
  }
  seconds <- proc.time()[['elapsed']]-began
  list(coords=grip::layout.coords(fit),seconds=seconds,preparation_seconds=prep,
    fitting_seconds=fitting,metadata=if (method=='grip') NULL else fit$metadata)
}

sparse_comparison_view <- function(case, fits, adjust_scale = TRUE) {
  fits <- fits[vapply(fits,function(f) !is.null(f$coords),logical(1))]
  if (!length(fits)) stop('No successful fits to display')
  if (adjust_scale) fits <- lapply(fits,function(f) { f$coords <- f$coords*f$row$relative_scale; f })
  reference <- !is.null(case$X)
  if (!reference) case$X <- fits[[1]]$coords
  labels <- vapply(names(fits),function(k) if (k=='full') 'Full SGD' else if (k=='grip')
    'Weighted GRIP' else paste('Sparse SGD:',sub('sparse','',k),'pivots'),'')
  palette <- stats::setNames(c('#1769AA','#009E73','#D66A19','#AA3377','#CC79A7')[seq_along(fits)],labels)
  if (reference) palette <- c(Reference='#888888',palette)
  comparison_view(case,fits,series_labels=labels,series_colors=palette,reference=reference,layout_selector=TRUE,
    description=paste(case$label,case$n,'vertices. Rotate to compare layouts.',
      if (reference) 'Gray: generating reference.' else 'No generating reference is available.',
      if (adjust_scale) 'Each method is resized by its independently fitted scalar; Procrustes then aligns orientation.' else 'Input-distance units; Procrustes preserves scale.'))
}

sparse_comparison_rows <- function(bundle) do.call(rbind,lapply(bundle$fits,`[[`,'row'))
sparse_comparison_summary <- function(bundle) {
  rows <- sparse_comparison_rows(bundle)
  keys <- unique(rows[,c('case','method')])
  do.call(rbind,lapply(seq_len(nrow(keys)),function(i) {
    z <- rows[rows$case==keys$case[i] & rows$method==keys$method[i],]
    ok <- z$status=='ok'
    med <- function(name) if (any(ok)) stats::median(z[[name]][ok],na.rm=TRUE) else NA_real_
    data.frame(keys[i,],returned=sum(ok),attempted=nrow(z),
      relative_error=med('relative_error'),profiled_relative_error=med('profiled_relative_error'),
      relative_scale=med('relative_scale'),seconds=med('seconds'),peak_mib=med('peak_mib'))
  }))
}

# One scalar fitted to the independently evaluated distances, not a new layout.
# Retain original error and coordinates; this diagnoses global size differences.
sparse_comparison_profile <- function(bundle) {
  for (key in names(bundle$fits)) {
    f <- bundle$fits[[key]]
    if (f$row$status!='ok') {
      f$row$relative_scale <- f$row$profiled_relative_error <- NA_real_
    } else {
      eval <- bundle$evaluation[[f$row$case]]; n <- nrow(f$coords)
      ratios <- unlist(lapply(seq_along(eval$sources),function(k) {
        i <- eval$sources[k]; keep <- seq_len(n)!=i
        sqrt(rowSums(sweep(f$coords,2,f$coords[i,])^2))[keep]/eval$targets[k,keep]
      }))
      stopifnot(isTRUE(all.equal(sqrt(mean((ratios-1)^2)),f$row$relative_error,tolerance=1e-10)))
      multiplier <- sum(ratios)/sum(ratios^2)
      f$row$relative_scale <- multiplier
      f$row$profiled_relative_error <- sqrt(mean((multiplier*ratios-1)^2))
    }
    bundle$fits[[key]] <- f
  }
  bundle
}
