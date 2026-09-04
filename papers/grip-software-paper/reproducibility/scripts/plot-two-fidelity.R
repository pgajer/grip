# Shared publication figures. Functions draw only; the Rmd supplies captions.
pilot_cloud_colors <- c("#1E5C89", "#B26026", "#777777", "#5085AD", "#927457")

pilot_align <- function(z, x) {
  zz <- sweep(z, 2, colMeans(z), "-")
  xx <- sweep(x, 2, colMeans(x), "-")
  f <- svd(crossprod(zz, xx))
  sweep(sum(f$d) / sum(zz * zz) * zz %*% f$u %*% t(f$v), 2, colMeans(x), "+")
}

plot_pilot_calibration_layouts <- function(pilot, edge_alpha = .03,
                                         mesh_dir = "reproducibility/figures/saddle") {
  stopifnot(length(edge_alpha) == 1L, is.finite(edge_alpha),
            edge_alpha >= 0, edge_alpha <= 1)
  op <- par(no.readonly = TRUE)
  on.exit({layout(1); par(op)})
  layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5), nrow = 2, byrow = TRUE),
         heights = c(1, 1.05))
  par(oma = c(0, 0, 0, 0), cex = 1, las = 1, mgp = c(2.5, .7, 0))
  for (zoom in c(FALSE, TRUE)) {
    par(mar = c(4.0, 4.5, 2.1, .8))
    d <- if (zoom) pilot$curves[pilot$curves$k >= 40, ] else pilot$curves
    plot(NA, xlim = range(d$k), ylim = c(0, max(d$xg_error) * 106),
         xlab = "Number of neighbors k", ylab = "Surface-to-graph RMSE (%)", bty = "l")
    abline(h = axTicks(2), col = "gray92", lwd = .6)
    for (r in 1:5) {
      a <- d[d$replicate == r, ]
      lines(a$k, 100 * a$xg_error, col = pilot_cloud_colors[r], lty = r, lwd = 1.4)
      b <- a[which.min(a$xg_error), ]
      points(b$k, 100 * b$xg_error, pch = 21, bg = "white",
             col = pilot_cloud_colors[r], cex = 1.05, lwd = 1.3)
    }
    title(if (zoom) "B   Enlargement near the minima" else "A   Graph calibration",
          cex.main = 1.1, font.main = 1)
    if (!zoom) legend("topright", paste("Cloud", 1:5), col = pilot_cloud_colors,
                      lty = 1:5, lwd = 1.4, bty = "n", cex = .9)
  }
  z <- lapply(pilot$candidates, function(a)
    grip::project.3d(pilot_align(a, pilot$coords), azimuth = 35, elevation = 22))
  all <- do.call(rbind, z)
  xr <- range(all[, 1]); yr <- range(all[, 2])
  xr <- xr + c(-1, 1) * .04 * diff(xr); yr <- yr + c(-1, 1) * .04 * diff(yr)
  vc <- colorRampPalette(c("#173D65", "#86AFC4", "#D9B18B", "#8E4921"))(100)[
    pmin(100, pmax(1, 1 + floor((pilot$coords[, 1] + 1) * 49.5)))]
  for (j in 1:3) {
    par(mar = c(.3, .1, 2.2, .1))
    if (!is.null(mesh_dir)) {
      picture <- png::readPNG(file.path(mesh_dir,paste0("mesh",j,".png")))
      plot.new(); plot.window(xlim=c(0,4),ylim=c(0,3),asp=1,xaxs="i",yaxs="i")
      rasterImage(picture,0,0,4,3)
      title(c("C   Original saddle","D   metric-MDS","E   metric-MDS + edge-KK")[j],
            font.main=1,cex.main=1.08,line=.5)
      next
    }
    a <- z[[j]]; e <- pilot$edges
    plot(a, type = "n", xlim = xr, ylim = yr, asp = 1, axes = FALSE, xlab = "", ylab = "")
    if (edge_alpha > 0) {
      segments(a[e[,1],1], a[e[,1],2], a[e[,2],1], a[e[,2],2],
               col = adjustcolor("gray45", alpha.f = edge_alpha), lwd = .35)
    }
    points(a, pch = 16, cex = .34, col = vc)
    e <- pilot$route
    segments(a[e[,1],1], a[e[,1],2], a[e[,2],1], a[e[,2],2], col = "#1E5C89", lwd = 2.3)
    e <- pilot$ends
    segments(a[e[1],1], a[e[1],2], a[e[2],1], a[e[2],2], col = "#B26026", lty = 2, lwd = 1.5)
    points(a[e, ], pch = 21, bg = "white", col = "gray20", cex = .8)
    title(c("C   Original saddle", "D   Metric MDS", "E   MDS + edge-KK")[j],
          font.main = 1, cex.main = 1.08, line = .5)
  }
}

plot_pilot_scores <- function(pilot, audit = FALSE, reference = NULL) {
  op <- par(no.readonly = TRUE); on.exit(par(op))
  methods <- c("Original saddle", "Metric MDS", "MDS + edge-KK")
  metrics <- c("path_rel", "edge_rel", "stress1")
  if(!is.null(reference)) metrics <- c(metrics,"coordinate_relative_rmse","surface_rms")
  colors <- c("#737373", "#286EAB", "#B66027")
  offsets <- seq(-.16, .16, length.out = 5)
  par(mfrow = if(is.null(reference)) c(1,3) else c(2,3), mar = c(4.2,4.2,2.1,.7), oma = c(0,0,0,0),
      cex = 1, mgp = c(2.6,.7,0), las = 1, tcl = -.25)
  for (m in seq_along(metrics)) {
    metric <- metrics[m]
    y <- vapply(methods, function(method) {
      d <- if(m<=3) pilot$scores[pilot$scores$method == method, ] else
        subset(reference$scores,alignment=="similarity" & sample_size==8000)
      if(m>3) d <- d[d$method==method,]
      (if(m==5) 1 else 100) * d[[metric]][match(1:5, d$replicate)]
    }, numeric(5))
    ymax <- max(y); medians <- apply(y, 2, median)
    plot(NA, xlim = c(.65, 3.35), ylim = c(-.045 * ymax, 1.14 * ymax),
         xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
         xlab = "", ylab = c("Relative RMSE (%)","Relative RMSE (%)","Stress-1 (%)",
                              "Reference-radius RMSE (%)","Coordinate units")[m])
    ticks <- pretty(c(0, ymax), n = 4); ticks <- ticks[ticks >= 0 & ticks <= 1.1 * ymax]
    abline(h = ticks, col = "#E9E9E9", lwd = .7)
    axis(2, at = ticks, col = "#555555", lwd = .75)
    axis(1, at = 1:3, labels = c("Original\ncoordinates", "metric-MDS", "metric-MDS\n+ edge-KK"),
         tick = FALSE, line = .35, cex.axis = .78, gap.axis = -1)
    box(bty = "l", col = "#555555", lwd = .75)
    for (r in 1:5) lines(1:3 + offsets[r], y[r, ], col = "#C2C2C2", lwd = .8)
    for (j in 1:3) {
      points(j + offsets, y[,j], pch = 16, cex = .88, col = adjustcolor(colors[j], alpha.f = .7))
      points(j, medians[j], pch = 1, cex = 1.6, lwd = 1.65, col = colors[j])
      text(j, max(y[,j]) + .064 * ymax, sprintf("%.3f", medians[j]), col = colors[j], cex = .85)
    }
    if (audit && m<=3) {
      a <- 100 * pilot$audit[[metric]][match(1:5, pilot$audit$replicate)]
      points(3.24 + offsets * .25, a, pch = 4, col = "#555555", cex = .8)
    }
    title(c("A   Fixed-path RMSE", "B   Edge RMSE", "C   MDS Stress-1",
            "D   Coordinate agreement","E   Symmetric surface distance")[m],
          font.main = 1, cex.main = 1.03, line = .8)
  }
  if(!is.null(reference)) {
    par(mar=c(1,1,2,1))
    plot.new()
    text(.05,.8,"Five independent clouds",adj=0,cex=1)
    text(.05,.65,"Dots: individual clouds\nLines: paired observations\nOpen circles and labels:\nmedians",
         adj=c(0,1),cex=.9)
    text(.05,.3,"D-E: similarity alignment\nE: 8,000 area-uniform samples\nper direction",
         adj=c(0,1),cex=.9)
  }
}

plot_pilot_controls <- function(pilot) {
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mfrow = c(1,2), mar = c(2.7,.2,2.2,.2), cex = 1)
  z <- grip::project.3d(pilot$coords, azimuth = 35, elevation = 22)
  for (g in pilot$controls) {
    e <- g$edges
    plot(z, type = "n", asp = 1, axes = FALSE, xlab = "", ylab = "")
    segments(z[e[,1],1],z[e[,1],2],z[e[,2],1],z[e[,2],2],col = adjustcolor("gray40", alpha.f = .14), lwd = .4)
    points(z, pch = 16, cex = .32, col = "gray45")
    lines(z[g$route, ], col = "#1E5C89", lwd = 2.2)
    points(z[pilot$control_ends, ], pch = 21, bg = "white", col = "gray15", cex = .8)
    if (g$bridges > 0) {
      e <- g$mst_edges
      segments(z[e[,1],1],z[e[,1],2],z[e[,2],1],z[e[,2],2],col = "#B26026", lty = 2, lwd = 2)
    }
    error <- with(pilot$curves, xg_error[replicate == pilot$representative & k == g$k])
    title(sprintf("k = %d; %d MST bridges", g$k, g$bridges), font.main = 1)
    mtext(sprintf("Surface-to-graph error: %.3f%%", 100 * error), 1, line = .2)
    mtext("Original-coordinate fixed-path error: 0", 1, line = 1.4)
  }
}

plot_pilot_end_to_end <- function(pilot) {
  op <- par(no.readonly = TRUE); on.exit(par(op))
  y <- vapply(c("Original saddle", "Metric MDS", "MDS + edge-KK"), function(method) {
    d <- pilot$scores[pilot$scores$method == method, ]
    100 * d$xz_path_error[match(1:5, d$replicate)]
  }, numeric(5))
  par(mar = c(3.5,4.3,1,.8), cex = .95, las = 1)
  plot(NA, xlim = c(.8,3.2), ylim = c(0, max(y) * 1.08), xaxt = "n",
       xlab = "", ylab = "RMSE against surface distances (%)", bty = "l")
  abline(h = axTicks(2), col = "gray92")
  for (r in 1:5) lines(1:3, y[r, ], type = "b", pch = r, lty = r, col = pilot_cloud_colors[r])
  axis(1, at = 1:3, labels = c("Input graph", "MDS embedded paths", "MDS + edge-KK paths"), tick = FALSE)
  legend("topright", paste("Cloud", 1:5), col = pilot_cloud_colors, lty = 1:5, pch = 1:5, bty = "n", cex = .85)
}
