# Stored parameter records and native result fields retain their existing schema.
# Translate those records only at internal R call boundaries. This is not an
# alias mechanism for public calls: obsolete named arguments fail normally.
.grip.invoke <- function(what, args, quote = FALSE, envir = parent.frame()) {
  if (is.character(what)) what <- get(what, envir = envir, mode = "function")
  declared <- names(formals(what))
  old <- names(args)
  dotted <- gsub("_", ".", old, fixed = TRUE)
  change <- !old %in% declared & dotted %in% declared
  if (any(change)) {
    old[change] <- dotted[change]
    if (anyDuplicated(old[nzchar(old)])) stop("Conflicting internal parameter names")
    names(args) <- old
  }
  do.call(what, args, quote = quote, envir = envir)
}
