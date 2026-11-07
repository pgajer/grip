library(testthat)
library(grip)

# The designated full-dependency check must exercise the optional methods.
# Ordinary minimal installations may still skip tests of unavailable backends.
full <- identical(Sys.getenv("GRIP_FULL_DEPENDENCIES"), "true")
if (full) {
  cat("Full-dependency check library paths:\n")
  print(.libPaths())
  for (pkg in c("smacof", "shiny", "bslib", "DT", "htmltools", "htmlwidgets",
                "rgl", "later", "httpuv")) {
    cat("Loading", pkg, "from", find.package(pkg), "\n")
    loadNamespace(pkg) # Preserve the actual load error rather than silently skip.
  }
}
results <- test_check("grip")
if (full && any(as.data.frame(results)$skipped)) {
  stop("Full-dependency check skipped tests; inspect test results above")
}
