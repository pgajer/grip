# Stored parameter records and native result fields retain their existing schema.
# Translate those records at internal R call boundaries. Public compatibility
# is limited to the explicit transition entry points below.
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

# Normalize only exact former underscore spellings of documented arguments.
# The canonical signature stays visible to callers inspecting formals().
.grip.transition.names <- function(nms, declared) {
  dotted <- gsub("_", ".", nms, fixed = TRUE)
  aliases <- grepl("_", nms, fixed = TRUE) & dotted %in% declared
  nms[aliases] <- dotted[aliases]
  if (anyDuplicated(nms[nzchar(nms)]))
    stop("Conflicting argument names: supply each argument once, using either its old or current spelling", call. = FALSE)
  nms
}

# Re-evaluate the normalized call in its caller, without forcing arguments or
# filling defaults. This preserves missing(), preset handling and lazy values.
.grip.transition.call <- function(call, fun) {
  args <- as.list(call)[-1L]
  dots <- as.list(args[["..."]])
  args[["..."]] <- NULL
  nms <- names(dots)
  declared <- setdiff(names(formals(fun)), "...")
  if (is.null(nms) || any(!nzchar(nms)))
    stop("unused unnamed arguments in ...", call. = FALSE)
  translated <- gsub("_", ".", nms, fixed = TRUE)
  valid <- grepl("_", nms, fixed = TRUE) & translated %in% declared
  if (!all(valid))
    stop("unused argument(s): ", paste(nms[!valid], collapse = ", "), call. = FALSE)
  args <- c(args, dots)
  names(args) <- .grip.transition.names(names(args), declared)
  as.call(c(list(fun), args))
}
