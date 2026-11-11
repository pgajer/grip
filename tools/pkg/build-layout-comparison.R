#!/usr/bin/env Rscript
# Reproduce every fit and score in the fixed nine-layout GRIP demonstration.
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else ".", mustWork=TRUE)
output <- if (length(args)>1L) args[[2L]] else
  file.path(root,"inst","extdata","layout-comparison","comparison.rds")
library(grip)
version <- read.dcf(file.path(root,"DESCRIPTION"))[1L,"Version"]
stopifnot(as.character(packageVersion("grip")) == version)
lines <- readLines(file.path(root,"vignettes","grip-real-data.Rmd"))
chunk <- function(label) {
  first <- grep(paste0("^```\\{r ",label,"([,}])"),lines)
  stopifnot(length(first)==1L)
  last <- first + which(lines[(first+1L):length(lines)]=="```")[[1L]]
  paste(lines[(first+1L):(last-1L)],collapse="\n")
}
recipe <- paste(chunk("candidates-graph"),chunk("candidates-fit"),sep="\n")
env <- new.env(parent=globalenv()); warnings <- character()
withCallingHandlers(eval(parse(text=recipe),envir=env),warning=function(w) {
  warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")
})
stopifnot(nrow(env$cmp$runs)==9L, length(env$cmp$layouts)==3L,
          all(vapply(env$cmp$layouts,length,integer(1))==3L))
source <- sort(c(list.files(file.path(root,"R"),full.names=TRUE),
                 list.files(file.path(root,"src"),pattern="[.](cpp|h|hpp)$",full.names=TRUE)))
hashes <- tools::md5sum(source); names(hashes) <- substring(names(hashes),nchar(root)+2L)
bundle <- list(graph=env$g,comparison=env$cmp,recipe=recipe,warnings=warnings,
  provenance=list(package.version=version,source.md5=hashes,
    installed.dll.md5=unname(tools::md5sum(getLoadedDLLs()[["grip"]][["path"]])),
    session.info=capture.output(sessionInfo())))
dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE)
saveRDS(bundle,output,compress="xz",version=2)
cat("Saved all nine original GRIP layouts and scores:",output,"\n")
