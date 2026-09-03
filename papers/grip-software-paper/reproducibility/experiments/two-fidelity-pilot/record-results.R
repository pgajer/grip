#!/usr/bin/env Rscript
# Preserve the compact pilot readout in Git; large datasets and figures stay in build/.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
out<-normalizePath(args[1]);script<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
dest<-file.path(dirname(normalizePath(script)),"pilot-summary")
files<-c("results.md","calibration-curves.csv","reference-validation.csv","layout-scores.csv",
  "additional-budget-scores.csv","layout-timings.csv","fit-status.csv","independent-validation.csv",
  "representative-selection.csv","reference-environment.json","reference-run.json")
stopifnot(all(file.exists(file.path(out,files))))
dir.create(dest,showWarnings=FALSE)
stopifnot(all(file.copy(file.path(out,files),file.path(dest,files),overwrite=TRUE)))
inputs<-c(list.files(out,pattern="^cloud-[0-9]+[.]csv$",full.names=TRUE),
  list.files(out,pattern="^reference-r[0-9]+-m81-s128[.]npz$",full.names=TRUE),
  list.files(out,pattern="^fit-[0-9]+-k[0-9]+[.]rds$",full.names=TRUE))
write.csv(data.frame(file=basename(inputs),md5=unname(tools::md5sum(inputs))),
  file.path(dest,"input-checksums.csv"),row.names=FALSE)
