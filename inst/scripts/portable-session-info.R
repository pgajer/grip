# Preserve reproducibility metadata without local DESCRIPTION-file attributes.
# This helper has no package dependencies and is also shipped with the paper.
portable_session_info <- function(info = utils::sessionInfo()) {
  for (section in c("otherPkgs", "loadedOnly")) {
    if (!is.null(info[[section]])) {
      info[[section]] <- lapply(info[[section]], function(description) {
        attr(description, "file") <- NULL
        description
      })
    }
  }
  info
}
