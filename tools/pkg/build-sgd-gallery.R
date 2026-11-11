#!/usr/bin/env Rscript
# Reproduce the gallery's fixed illustrations without changing their fit settings.
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else ".", mustWork=TRUE)
output <- if (length(args)>1L) args[[2L]] else
  file.path(root,"inst","extdata","sgd-gallery","layouts.rds")
library(grip)
version <- read.dcf(file.path(root,"DESCRIPTION"))[1L,"Version"]
stopifnot(as.character(packageVersion("grip")) == version)
lines <- readLines(file.path(root,"vignettes","gallery.Rmd"))
chunk <- function(label) {
  first <- grep(paste0("^```\\{r ",label,"([,}])"),lines)
  stopifnot(length(first)==1L)
  last <- first + which(lines[(first+1L):length(lines)]=="```")[[1L]]
  paste(lines[(first+1L):(last-1L)],collapse="\n")
}
names <- c("mesh","torus","vicsek","annulus","cube","tree")
env <- new.env(parent=globalenv())
graphs <- fits <- recipes <- warnings <- list()
for (name in names) {
  recipe <- paste(chunk(paste0(name,"-family")),chunk(paste0(name,"-fit")),sep="\n")
  recorded <- character()
  withCallingHandlers(eval(parse(text=recipe),envir=env),warning=function(w) {
    recorded <<- c(recorded,conditionMessage(w));invokeRestart("muffleWarning")
  })
  graph <- get(paste0(name,".graph"),env)
  fit <- get(paste0(name,".fit"),env)
  stopifnot(nrow(fit$coords)==graph$n,ncol(fit$coords)==3L,
            all(is.finite(fit$coords)),
            identical(as.integer(fit$metadata$settings$max_iter),1000L))
  graphs[[name]] <- graph
  fits[[name]] <- fit[c("coords","metadata")]
  recipes[[name]] <- recipe
  warnings[[name]] <- recorded
}
source <- sort(c(list.files(file.path(root,"R"),full.names=TRUE),
                 list.files(file.path(root,"src"),pattern="[.](cpp|h|hpp)$",full.names=TRUE)))
hashes <- tools::md5sum(source);names(hashes) <- substring(names(hashes),nchar(root)+2L)
bundle <- list(graphs=graphs,fits=fits,recipes=recipes,warnings=warnings,
  provenance=list(package.version=version,source.md5=hashes,
    installed.dll.md5=unname(tools::md5sum(getLoadedDLLs()[["grip"]][["path"]])),
    session.info=capture.output(sessionInfo()),
    scope="Fixed gallery illustrations; coordinates and metadata only, not complete fit objects"))
dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE)
saveRDS(bundle,output,compress="xz",version=2)
cat("Saved six unchanged 1,000-epoch illustrations:",output,"\n")
