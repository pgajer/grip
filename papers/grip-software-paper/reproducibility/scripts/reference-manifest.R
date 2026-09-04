# Regenerate after exporting data and rendering S4 and its static views.
root <- "papers/grip-software-paper"
files <- c("reproducibility/precomputed/saddle-reference-diagnostics.rds",
 "reproducibility/precomputed/saddle-reference-diagnostics.csv",
 "reproducibility/precomputed/single-saddle-reference.csv",
 paste0("reproducibility/figures/saddle/",c("mesh1","mesh2","mesh3","overlay","displacement"),".png"),
 "reproducibility/figures/saddle/panel-e-workflow.png",
 "supplement/S4-interactive-saddle.html")
stopifnot(all(file.exists(file.path(root,files))))
write.table(data.frame(md5=unname(tools::md5sum(file.path(root,files))),file=files),
 file.path(root,"reproducibility/reference-assets-md5.tsv"),
 sep="\t",row.names=FALSE,quote=FALSE)
