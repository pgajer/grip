test_that('saved surface displays reject changed inputs and preserve fitted distances', {
  e <- new.env()
  for(file in c('metric-mds-comparison.R','metric-mds-surface-alignment.R'))
    sys.source(system.file('scripts',file,package='grip'),e)
  bundle <- readRDS(system.file('extdata','metric-mds-comparison','benchmark.rds',package='grip'))
  case <- bundle$cases[['saddle_graph-64']]
  fits <- e$comparison_representative(bundle,case$id)
  original <- serialize(fits,NULL)
  display <- e$comparison.surface.saved(case,fits)
  expect_identical(serialize(fits,NULL),original)
  for(method in names(fits)) expect_equal(as.vector(dist(display$registrations[[method]]$coords)),
    as.vector(dist(fits[[method]]$coords)),tolerance=1e-10)
  fits$sgd$coords[1,1] <- fits$sgd$coords[1,1]+.001
  expect_error(e$comparison.surface.saved(case,fits),'outdated display registration')
})

test_that('teaching views center without changing pairwise distances', {
  e <- new.env()
  sys.source(system.file('scripts','vignette-views.R',package='grip'),e)
  e$documentation.interactive <- function() TRUE
  e$comparison_view <- function(case,fits,...) list(case=case,fits=fits,options=list(...))
  x <- rbind(c(10,20,30),c(11,20,30),c(10,22,30),c(10,20,33))
  y <- sweep(2*x,2,c(5,6,7),'+')
  view <- e$documentation.view(list('metric-MDS'=x,GRIP=y),matrix(c(1,2),1,2))
  expect_equal(unname(colMeans(view$case$X)),c(0,0,0),tolerance=1e-12)
  expect_equal(as.vector(dist(view$fits$GRIP$coords)),as.vector(dist(y)),tolerance=1e-10)
  expect_false(view$options$reference)
  same.size <- e$documentation.view(list('metric-MDS'=x,GRIP=y),
    matrix(c(1,2),1,2),common.size=TRUE)
  expect_equal(as.vector(dist(same.size$fits$GRIP$coords)),as.vector(dist(x)),tolerance=1e-10)
  expect_match(same.size$options$description,'Display multipliers',fixed=TRUE)
})
