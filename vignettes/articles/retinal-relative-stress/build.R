args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
article <- normalizePath(args[1])
stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", tz = "America/New_York")
writeLines(sprintf("\\renewcommand{\\reportbuilddatetime}{%s}", stamp),
           file.path(article, "build-info.tex"))
old <- setwd(article)
on.exit(setwd(old))
for (i in 1:2) {
  code <- system2("pdflatex", c("-interaction=nonstopmode", "-halt-on-error", "main.tex"),
                  stdout = "build-console.log", stderr = "build-errors.log")
  if (code != 0L) stop("LaTeX build failed; see build-console.log")
}
cat(file.path(article, "main.pdf"), "\n")
