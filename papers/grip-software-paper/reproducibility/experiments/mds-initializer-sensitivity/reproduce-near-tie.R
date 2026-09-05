# Minimal reproducer of the preparation integration limitation found in Phase 2.
# Diagnostic, not a package regression that enshrines the current failure.
e<-rbind(c(1L,2L),c(1L,3L),c(2L,3L));w<-c(1,2+1e-8,1)
p<-grip::prepare.geodesic.kk(e,n=3L,edge_weights=w,tie_mode="single")
asym<-max(abs(p$distance_matrix-t(p$distance_matrix)))
fit<-tryCatch(grip::metric.mds(prepared=p,dim=2L,diagnostics=FALSE),error=function(e) e)
strict<-grip::metric.mds(edges=e,n=3L,edge_weights=w,dim=2L,diagnostics=FALSE)
cat("Prepared matrix asymmetry:",format(asym,digits=12),"\n")
cat("Path-prepared call:",if(inherits(fit,"error")) conditionMessage(fit) else "succeeded","\n")
cat("Strict distance-only raw stress:",strict$metadata$raw_stress,"\n")
stopifnot(strict$metadata$raw_stress<1e-12)
