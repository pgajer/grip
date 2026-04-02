repo_root <- normalizePath(file.path(getwd(), "..", "..", ".."), winslash = "/", mustWork = FALSE)

if (file.exists(file.path(repo_root, "DESCRIPTION")) &&
    requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo_root, quiet = TRUE, export_all = FALSE, helpers = FALSE)
}

if (!requireNamespace("grip", quietly = TRUE)) {
  stop("The 'grip' package must be installed or loadable to run this app.", call. = FALSE)
}

grip::run_gripui_paraboloid_gmds()
