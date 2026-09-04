# Independently check distributed coordinates, scores, and original-score agreement.
pkgload::load_all(".",quiet=TRUE)
base <- "papers/grip-software-paper/reproducibility"
source(file.path(base,"experiments/two-fidelity-pilot/score-saddle-reference.R"))
r <- readRDS(file.path(base,"precomputed/saddle-reference-diagnostics.rds"))
p <- readRDS(file.path(base,"precomputed/two-fidelity-saddle.rds"))
d <- read.csv(file.path(base,"precomputed/saddle-reference-diagnostics.csv"))
stopifnot(isTRUE(all.equal(d,r$scores,tolerance=1e-12)),nrow(d)==60L,
          !anyDuplicated(d[c("replicate","method","alignment","sample_size")]))
for(id in names(r$clouds)) {
  cloud <- r$clouds[[id]]; X <- cloud$coords; t <- cloud$triangles
  stopifnot(nrow(X)==1000L,max(abs(X[,3]-.8*(X[,1]^2-X[,2]^2)))<1e-12,
            min(t)>=1L,max(t)<=1000L)
  ref <- saddle.reference.mesh(X,t,.8,2L)
  for(method in names(cloud$candidates)) {
    s <- subset(d,replicate==as.integer(id) & alignment=="similarity" & sample_size==8000)
    s <- s[s$method==method,]
    z <- cloud$aligned[[method]]
    independent <- sqrt(sum((z-X)^2)/sum(sweep(X,2,colMeans(X),"-")^2))
    old <- subset(p$scores,replicate==as.integer(id))
    stopifnot(abs(independent-s$coordinate_relative_rmse)<1e-12,
      abs(independent-old$procrustes[match(method,old$method)])<1e-10)
    surface <- grip::score.surface(z,t,ref$coords,ref$triangles,
      sample_size=8000L,seed=1901L+as.integer(id))
    stopifnot(abs(surface$rms-s$surface_rms)<1e-12)
  }
}
single <- read.csv(file.path(base,"precomputed/single-saddle-reference.csv"))
stopifnot(nrow(single)==216L,all(single$n==500L),identical(sort(unique(single$k)),3:20))
message("PASS: five saved clouds, 60 diagnostic records, 15 independent coordinate checks, ",
        "15 deterministic surface recomputations, unchanged historical Procrustes values, separate 500-point data.")
