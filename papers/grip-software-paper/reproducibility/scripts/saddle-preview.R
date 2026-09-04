# Repository-root entry point; static preview uses the manuscript's reviewed views.
base <- "papers/grip-software-paper"
png("man/figures/readme-saddle-reference.png",width=1800,height=540,res=160)
par(mfrow=c(1,3),mar=c(0,0,2,0))
for(j in 1:3) {
  im <- png::readPNG(file.path(base,"reproducibility/figures/saddle",paste0("mesh",j,".png")))
  plot.new(); plot.window(xlim=c(0,4),ylim=c(0,3),asp=1,xaxs="i",yaxs="i")
  rasterImage(im,0,0,4,3)
  title(c("Original saddle","metric-MDS","metric-MDS + edge-KK")[j],font.main=1,cex.main=1)
}
dev.off()
