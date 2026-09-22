#!/usr/bin/env Rscript
# Repository-only checks: public API documentation and installed-data privacy.
stopifnot(file.exists("DESCRIPTION"), file.exists("NAMESPACE"))
public <- sub("^export\\((.*)\\)$", "\\1",
              grep("^export\\(", readLines("NAMESPACE"), value = TRUE))
rd_files <- list.files("man", pattern = "\\.Rd$", full.names = TRUE)
if (!length(public) || !length(rd_files)) {
  stop("Public exports and package Rd documentation must both be present.")
}
aliases <- unique(unlist(lapply(rd_files, function(file) {
  rd <- tools::parse_Rd(file)
  nodes <- Filter(function(node) identical(attr(node, "Rd_tag"), "\\alias"), rd)
  vapply(nodes, function(node) paste(unlist(node), collapse = ""), character(1))
}), use.names = FALSE))
missing <- setdiff(public, aliases)
if (length(missing)) {
  stop("Public exports lack package Rd aliases: ", paste(missing, collapse = ", "))
}

inspect_paths <- function(x, location) {
  if (is.character(x) && any(grepl("(/Users/|/home/)[^/]+/", x))) {
    stop("Personal filesystem path in serialized data: ", location)
  }
  if (is.list(x)) for (i in seq_along(x)) {
    inspect_paths(x[[i]], paste0(location, "/", i))
  }
  # Inspect attribute values, not the attributes-list wrapper recursively.
  for (name in names(attributes(x))) {
    inspect_paths(attr(x, name), paste0(location, "@", name))
  }
}
tracked <- system2("git", "ls-files", stdout = TRUE)
for (file in tracked[grepl("\\.(rds|rda)$", tracked)]) {
  if (!file.exists(file)) next
  if (grepl("\\.rds$", file)) {
    inspect_paths(readRDS(file), file)
  } else {
    values <- new.env(parent = emptyenv())
    objects <- load(file, envir = values)
    for (name in objects) inspect_paths(values[[name]], paste0(file, "/", name))
  }
}
text_files <- tracked[grepl("\\.(md|Rmd|R|tex|cpp|h|json|yml)$", tracked)]
for (file in setdiff(text_files, c("AGENTS.md", "tools/check-release-content.R"))) {
  if (!file.exists(file)) next
  if (any(grepl("\\.codex/private/", readLines(file, warn = FALSE)))) {
    stop("Public content refers to private working material: ", file)
  }
}
cat(length(public), "public exports have Rd documentation; serialized data and private-reference checks passed.\n")
