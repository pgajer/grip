args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript build.R ARTICLE_SOURCE OUTPUT_DIRECTORY")
article <- normalizePath(args[1], mustWork = TRUE)
dir.create(path.expand(args[2]), recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(path.expand(args[2]), mustWork = TRUE)
if (identical(article, output) || startsWith(output, paste0(article, "/")))
  stop("Choose an output directory outside the article source directory.")
assets <- list.files(article, pattern = "\\.(tex|pdf|csv)$", full.names = TRUE)
assets <- assets[!basename(assets) %in% c("main.pdf", "build-info.tex")]
stopifnot(all(file.copy(assets, output, overwrite = TRUE)))
stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", tz = "America/New_York")
writeLines(sprintf("\\renewcommand{\\reportbuilddatetime}{%s}", stamp),
           file.path(output, "build-info.tex"))
old <- setwd(output)
on.exit(setwd(old))
for (i in 1:2) {
  code <- system2("pdflatex", c("-interaction=nonstopmode", "-halt-on-error", "main.tex"),
                  stdout = "build-console.log", stderr = "build-errors.log")
  if (code != 0L) stop("LaTeX build failed; see build-console.log")
}
cat(file.path(output, "main.pdf"), "\n")
