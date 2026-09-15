#!/usr/bin/env Rscript
# Run from the package root. The guide, not a generated inventory, is maintained.
ns <- parse("NAMESPACE")
entries <- function(tag) Filter(function(x) is.call(x) && identical(x[[1L]], as.name(tag)), as.list(ns))
exports <- vapply(entries("export"), function(x) as.character(x[[2L]]), "")
methods <- vapply(entries("S3method"), function(x) paste(x[[2L]], x[[3L]], sep = "."), "")
guide <- readLines("vignettes/function-guide.Rmd", warn = FALSE)
catalog_section <- guide[seq.int(match("## Function catalog", guide) + 1L,
                                match("## Find detailed help", guide) - 1L)]
rows <- grep("^\\| `[^`]+[(][)]` \\|", catalog_section, value = TRUE)
if (!all(grepl('help_url\\("', rows))) stop("Every catalog row needs a help target")
catalog <- sub("^\\| `([^`]+)[(][)]`.*$", "\\1", rows)
fail <- function(label, x) if (length(x)) stop(label, ": ", paste(x, collapse = ", "))
fail("Missing exports", setdiff(exports, catalog))
fail("Non-exports in catalog", setdiff(catalog, exports))
fail("Duplicate catalog rows", unique(catalog[duplicated(catalog)]))
rd <- lapply(list.files("man", pattern = "[.]Rd$", full.names = TRUE), tools::parse_Rd)
get_tag <- function(doc, tag) vapply(Filter(function(x) identical(attr(x, "Rd_tag"), tag), doc), function(x) paste(unlist(x), collapse = ""), "")
aliases <- lapply(rd, get_tag, tag = "\\alias")
topics <- vapply(rd, function(x) get_tag(x, "\\name")[[1L]], "")
all_aliases <- unique(unlist(aliases))
fail("Exports/methods without help", setdiff(c(exports, methods), all_aliases))
targets <- sub('.*help_url\\("([^"]+)"\\).*', '\\1', grep('help_url\\("', guide, value = TRUE))
fail("Unknown guide help targets", setdiff(targets, topics))
row_targets <- sub('.*help_url\\("([^"]+)"\\).*', '\\1', rows)
for (i in seq_along(catalog)) {
  if (!catalog[[i]] %in% aliases[[match(row_targets[[i]], topics)]]) {
    stop("Help target does not document ", catalog[[i]])
  }
}
fail("Unexplained registered methods", methods[!vapply(methods, function(x) any(grepl(x, guide, fixed = TRUE)), logical(1))])
count_line <- sprintf("**%d explicit public exports and %s registered S3", length(exports), c("zero", "one", "two", "three")[[length(methods) + 1L]])
if (!any(grepl(count_line, guide, fixed = TRUE))) stop("Update the guide's headline counts")
if (!any(grepl(sprintf("so there are %d", length(union(exports, methods))), guide, fixed = TRUE))) stop("Update unique interface count")
# Every source-level public function exists, without evaluating any body/defaults.
pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
for (name in union(exports, methods)) stopifnot(is.function(get(name, asNamespace("grip"))))
config <- yaml::read_yaml("_pkgdown.yml")
reference <- unlist(lapply(config$reference, `[[`, "contents"))
fail("Duplicate pkgdown topics", unique(reference[duplicated(reference)]))
fail("Unknown pkgdown topics", setdiff(reference, topics))
fail("Missing pkgdown topics", setdiff(unique(c(row_targets, methods)), reference))
for (name in c("function-guide", "synthetic-graph-families")) {
  path <- paste0("vignettes/", name, ".Rmd")
  lines <- readLines(path)
  stopifnot(any(grepl("%\\VignetteEngine{knitr::rmarkdown}", lines, fixed = TRUE)))
  fail("Private/local path in vignette", grep("/Users/|/tmp/|\\.codex/private|geosmooth::|dgraphs::|grip:::", lines, value = TRUE))
  links <- regmatches(lines, gregexpr("\\]\\([a-zA-Z0-9-]+[.]html\\)", lines))
  files <- sub("^]\\((.*)[.]html\\)$", "\\1.Rmd", unlist(links))
  fail("Missing linked vignette", files[!file.exists(file.path("vignettes", files))])
}
cat(sprintf("Verified %d/%d catalog exports; %d S3 registrations; %d overlap; %d unique interfaces.\n", length(catalog), length(exports), length(methods), length(intersect(exports, methods)), length(union(exports, methods))))
cat("All catalog help aliases, method help, pkgdown topics, metadata, and local vignette links resolve.\n")
