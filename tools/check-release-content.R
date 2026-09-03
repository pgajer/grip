#!/usr/bin/env Rscript
# Repository-only checks: installed-data privacy and design-catalog API status.
stopifnot(file.exists("DESCRIPTION"), file.exists("NAMESPACE"))
public <- sub("^export\\((.*)\\)$", "\\1",
              grep("^export\\(", readLines("NAMESPACE"), value = TRUE))
catalog <- readLines("dev/design/graph-families/graph_families_generated_in_thread_2026-03-31.md")
rows <- grep("^\\| `", catalog, value = TRUE)
for (row in rows) {
  name <- sub("^\\| `([^`]+)`.*", "\\1", row)
  expected <- if (name %in% public) "exported" else "internal"
  if (!endsWith(row, paste0("| ", expected, " |"))) {
    stop("Graph-family catalog API status does not match NAMESPACE: ", name)
  }
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
cat(length(rows), "catalog entries match NAMESPACE; serialized data and private-reference checks passed.\n")
