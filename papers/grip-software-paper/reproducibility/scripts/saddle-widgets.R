# Fixed-connectivity views: triangles are defined in the observed parameter plane.
# Render-only dependencies: ivue (872f9d4), rgl, geometry, htmlwidgets.
saddle_widgets <- function(pilot, reference) {
  cloud <- reference$clouds[[as.character(pilot$representative)]]
  X <- cloud$coords; triangles <- cloud$triangles; Z <- cloud$aligned
  bounds <- apply(rbind(do.call(rbind,Z),c(-1,-1,-.8),c(1,1,.8)),2,range)
  bounds[1,] <- bounds[1,]-.05
  bounds[2,] <- bounds[2,]+.05
  corners <- as.matrix(expand.grid(bounds[,1],bounds[,2],bounds[,3]))
  common <- ivue::layer3D.callback(function(ctx) {
    rgl::points3d(corners,col="white",alpha=0,size=1)
  })
  colors <- colorRampPalette(c("#173D65","#86AFC4","#D9B18B","#8E4921"))(100)[
    pmin(100,pmax(1,1+floor((X[,1]+1)*49.5)))]
  mesh <- ivue::layer3D.mesh(triangles,col="gray75",alpha=.23,
    edge.col="gray45",edge.alpha=.23,edge.width=.65)
  grid <- seq(-1,1,length.out=81)
  analytic <- ivue::layer3D.surface(grid,grid,outer(grid,grid,function(x,y).8*(x*x-y*y)),
    col="lightblue",alpha=.25,edges=FALSE)
  draw <- function(z,layers,col=colors) ivue::plot3D.plain(z,col=col,
    point.type="sphere",sphere.radius=.009,axes=FALSE,aspect="equal",
    layers=c(layers,list(common)),camera=ivue::camera.zup(elevation=20,turn=-135,zoom=.7),
    width=720L,height=640L)
  views <- list()
  for(j in seq_along(Z)) {
    # A route and its chord are diagnostic overlays, not mesh edges.
    route <- ivue::layer3D.edges(pilot$route,col="#1E5C89",width=4)
    chord <- ivue::layer3D.edges(matrix(pilot$ends,1,2),col="#B26026",width=3)
    views[[paste0("mesh",j)]] <- draw(Z[[j]],list(mesh,route,chord))
    views[[paste0("reference",j)]] <- draw(Z[[j]],list(analytic))
  }
  joined <- rbind(Z[[2]],Z[[3]]); n <- nrow(X)
  overlay_layers <- list(
    ivue::layer3D.mesh(triangles,col="#286EAB",alpha=.20,edges=FALSE),
    ivue::layer3D.mesh(triangles+n,col="#B66027",alpha=.20,edges=FALSE))
  overlay_colors <- rep(c("#286EAB","#B66027"),each=n)
  views$overlay <- draw(joined,overlay_layers,overlay_colors)
  # Same aligned fits, triangulation, and camera as A; only the reference is added.
  reference_gray <- ivue::layer3D.surface(grid,grid,
    outer(grid,grid,function(x,y).8*(x*x-y*y)),
    col="gray70",alpha=.22,edges=FALSE)
  views$overlay_reference <- draw(joined,c(overlay_layers,list(reference_gray)),overlay_colors)
  # A fixed systematic subset keeps correspondence segments readable.
  selected <- seq(1L,n,by=20L)
  views$displacement <- draw(joined,list(
    ivue::layer3D.edges(cbind(selected,selected+n),col="gray35",width=2)),
    rep(c("#286EAB","#B66027"),each=n))
  views
}
