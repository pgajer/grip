#!/usr/bin/env Rscript

package_library <- Sys.getenv("GRIP_RJOURNAL_PACKAGE_LIBRARY", unset = "")
if (nzchar(package_library)) {
  if (!dir.exists(package_library)) {
    stop("GRIP_RJOURNAL_PACKAGE_LIBRARY does not exist: ", package_library)
  }
  .libPaths(c(normalizePath(package_library, mustWork = TRUE), .libPaths()))
}

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to render grip-paper.Rmd.")
}

if (!requireNamespace("rjtools", quietly = TRUE)) {
  stop("Package 'rjtools' is required to render the R Journal paper.")
}

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (!length(file_arg)) {
  stop("Could not determine the path to this script.")
}

script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1L])))
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)

paper_candidates <- list(
  list(
    dir = file.path(
      repo_root,
      "papers",
      "grip-software-paper"
    ),
    input = "grip-software-paper.Rmd"
  ),
  list(
    dir = file.path(
      repo_root,
      "dev",
      "papers",
      "rjournal_paper",
      "manuscript",
      "legacy_r_journal_drafts"
    ),
    input = "grip-paper.Rmd"
  )
)

selected <- NULL
for (candidate in paper_candidates) {
  if (file.exists(file.path(candidate$dir, candidate$input))) {
    selected <- candidate
    break
  }
}

if (is.null(selected)) {
  stop("Could not locate the GRIP software-paper manuscript source.")
}

paper_dir <- normalizePath(selected$dir, mustWork = TRUE)
input_file <- selected$input
build_dir <- file.path(paper_dir, "build")
build_info_tex <- file.path(build_dir, "manuscript_build_info.tex")

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(paper_dir)

args <- commandArgs(trailingOnly = TRUE)
render_pdf <- !("--html-only" %in% args)
render_html <- "--html" %in% args || "--all" %in% args
timestamped <- !("--no-timestamp" %in% args)

timestamp_suffix <- format(Sys.time(), "%Y%m%d_%H%M%S", tz = "America/New_York")
base_name <- tools::file_path_sans_ext(input_file)

escape_tex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x, perl = TRUE)
  x
}

write_build_info <- function() {
  dir.create(build_dir, recursive = TRUE, showWarnings = FALSE)

  show_build_stamp <- tolower(Sys.getenv(
    "GRIP_RJOURNAL_SHOW_BUILD_STAMP",
    unset = "false"
  )) %in% c("1", "true", "yes", "on")

  git_version <- tryCatch(
    system2("git", c("-C", repo_root, "describe", "--tags", "--always", "--dirty"), stdout = TRUE, stderr = FALSE),
    error = function(e) "unversioned"
  )
  if (!length(git_version)) git_version <- "unversioned"

  git_build_number <- tryCatch(
    system2("git", c("-C", repo_root, "rev-list", "--count", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) "NA"
  )
  if (!length(git_build_number)) git_build_number <- "NA"

  git_commit <- tryCatch(
    system2("git", c("-C", repo_root, "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) "unknown"
  )
  if (!length(git_commit)) git_commit <- "unknown"

  build_datetime <- format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S %Z",
    tz = "America/New_York"
  )

  lines <- c(
    if (show_build_stamp) "\\grippapershowbuildstamptrue" else "\\grippapershowbuildstampfalse",
    sprintf("\\renewcommand{\\manuscriptversion}{%s}", escape_tex(git_version[[1L]])),
    sprintf("\\renewcommand{\\manuscriptbuildnumber}{%s}", escape_tex(git_build_number[[1L]])),
    sprintf("\\renewcommand{\\manuscriptcommit}{%s}", escape_tex(git_commit[[1L]])),
    sprintf("\\renewcommand{\\manuscriptbuilddatetime}{%s}", escape_tex(build_datetime))
  )
  writeLines(lines, build_info_tex, useBytes = TRUE)
}

write_build_info()

output_name <- function(ext) {
  if (timestamped) {
    file.path(build_dir, sprintf("%s_%s.%s", base_name, timestamp_suffix, ext))
  } else {
    file.path(build_dir, sprintf("%s.%s", base_name, ext))
  }
}

stable_name <- function(ext) {
  file.path(build_dir, sprintf("%s.%s", base_name, ext))
}

move_generated_directory <- function(source, target) {
  if (!dir.exists(source)) {
    return(invisible(FALSE))
  }
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  entries <- list.files(source, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  for (entry in entries) {
    destination <- file.path(target, basename(entry))
    if (dir.exists(entry)) {
      if (dir.exists(destination)) {
        move_generated_directory(entry, destination)
      } else if (!file.rename(entry, destination)) {
        if (!file.copy(entry, destination, recursive = TRUE)) {
          stop("Failed to move generated directory ", entry, " to ", destination)
        }
        unlink(entry, recursive = TRUE, force = TRUE)
      }
    } else {
      if (file.exists(destination)) {
        unlink(destination, force = TRUE)
      }
      if (!file.rename(entry, destination)) {
        if (!file.copy(entry, destination, overwrite = TRUE)) {
          stop("Failed to move generated file ", entry, " to ", destination)
        }
        unlink(entry, force = TRUE)
      }
    }
  }
  unlink(source, recursive = TRUE, force = TRUE)
  invisible(TRUE)
}

move_pdf_products_to_build <- function() {
  generated_files <- c(
    sprintf("%s.%s", base_name, c("R", "tex", "pdf", "log", "knit.md")),
    basename(Sys.glob(file.path(paper_dir, "RJwrapper.*")))
  )
  for (generated_file in unique(generated_files)) {
    source <- file.path(paper_dir, generated_file)
    if (!file.exists(source)) {
      next
    }
    destination <- file.path(build_dir, generated_file)
    if (file.exists(destination)) {
      unlink(destination, recursive = TRUE, force = TRUE)
    }
    if (!file.rename(source, destination)) {
      if (!file.copy(source, destination, overwrite = TRUE)) {
        stop("Failed to move generated file ", source, " to ", destination)
      }
      unlink(source, recursive = TRUE, force = TRUE)
    }
  }
  move_generated_directory(
    file.path(paper_dir, paste0(base_name, "_files")),
    file.path(build_dir, paste0(base_name, "_files"))
  )
  invisible(TRUE)
}

strip_trailing_whitespace <- function(path) {
  if (!file.exists(path)) {
    return(invisible(FALSE))
  }

  lines <- readLines(path, warn = FALSE)
  cleaned <- sub("[[:space:]]+$", "", lines)
  if (!identical(lines, cleaned)) {
    writeLines(cleaned, path, useBytes = TRUE)
  }
  invisible(TRUE)
}

add_table_aria_labels <- function(path) {
  if (!file.exists(path)) {
    stop("Cannot add table accessibility labels; missing HTML: ", path)
  }
  lines <- readLines(path, warn = FALSE)
  table_lines <- which(trimws(lines) == "<table>")
  if (!length(table_lines)) {
    stop("No HTML tables found while adding accessibility labels: ", path)
  }
  for (table_line in table_lines) {
    following <- seq.int(table_line + 1L, min(length(lines), table_line + 5L))
    caption_line <- following[grepl("<caption>", lines[following], fixed = TRUE)][1L]
    if (is.na(caption_line)) {
      stop("HTML table has no nearby caption at line ", table_line)
    }
    label <- gsub("<[^>]+>", "", lines[[caption_line]])
    label <- sub("^Table [0-9]+: *", "", label)
    label <- gsub("&", "&amp;", label, fixed = TRUE)
    label <- gsub('"', "&quot;", label, fixed = TRUE)
    label <- gsub("<", "&lt;", label, fixed = TRUE)
    label <- gsub(">", "&gt;", label, fixed = TRUE)
    lines[[table_line]] <- sprintf('<table aria-label="%s">', label)
  }
  writeLines(lines, path, useBytes = TRUE)
  invisible(length(table_lines))
}

style_package_name_html <- function(path) {
  if (!file.exists(path)) {
    stop("Cannot style package names; missing HTML: ", path)
  }

  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  literal_title <- "grip: Multiscale Graph Layout and Geodesic Embedding Tools in R"
  marked_title <- "`grip`: Multiscale Graph Layout and Geodesic Embedding Tools in R"
  display_title <- paste0(
    '<code class="grip-title">grip</code>',
    ': Multiscale Graph Layout and Geodesic Embedding Tools in R'
  )
  title_css <- paste0(
    '<style id="grip-package-typography">',
    '.grip-title{font-family:ui-monospace,SFMono-Regular,Menlo,Monaco,',
    'Consolas,"Liberation Mono",monospace;font-size:.9em;font-weight:700;',
    'color:inherit;background:transparent;border:0;padding:0;white-space:nowrap}',
    '.pkg-name{font-weight:700}',
    '</style>'
  )

  # Keep bibliographic and social metadata literal: formatting marks belong
  # only in the displayed title.
  html <- gsub(marked_title, literal_title, html, fixed = TRUE)
  html <- gsub("Kamada--Kawai", "Kamada\u2013Kawai", html, fixed = TRUE)
  html <- gsub("\\\\pkg{grip}", "grip", html, fixed = TRUE)
  html <- gsub("\\pkg{grip}", "grip", html, fixed = TRUE)

  # rjtools maps \pkg{} to a placeholder link in HTML. Replace that
  # representation with the bold package-name convention used by the journal.
  for (package_name in c("grip", "dgraphs", "igraph", "ggraph", "graphlayouts")) {
    html <- gsub(
      sprintf('<a href="#">%s</a>', package_name),
      sprintf('<strong class="pkg-name">%s</strong>', package_name),
      html,
      fixed = TRUE
    )
  }
  html <- sub(
    paste0("<h1>", literal_title, "</h1>"),
    paste0("<h1>", display_title, "</h1>"),
    html,
    fixed = TRUE
  )
  html <- sub("</head>", paste0(title_css, "\n</head>"), html, fixed = TRUE)

  writeLines(html, path, useBytes = TRUE)
  invisible(TRUE)
}

render_one <- function(output_format, ext) {
  stable_out <- stable_name(ext)
  render_dir <- if (identical(ext, "pdf")) paper_dir else build_dir
  rmarkdown::render(
    input = input_file,
    output_format = output_format,
    output_dir = render_dir,
    intermediates_dir = render_dir,
    knit_root_dir = paper_dir,
    quiet = FALSE,
    clean = FALSE,
    envir = new.env(parent = globalenv())
  )
  if (identical(ext, "pdf")) {
    move_pdf_products_to_build()
    strip_trailing_whitespace(stable_name("tex"))
    # rjtools compiles RJwrapper.tex. A direct-LaTeX attempt can leave a
    # same-basename log containing a fatal \maketitle error even though the
    # wrapper build succeeded; do not retain that misleading transient log.
    direct_tex_log <- stable_name("log")
    if (file.exists(direct_tex_log)) unlink(direct_tex_log)
  } else if (identical(ext, "html")) {
    add_table_aria_labels(stable_out)
    style_package_name_html(stable_out)
  }
  if (!timestamped) {
    return(stable_out)
  }

  stamped_out <- output_name(ext)
  if (file.exists(stamped_out)) {
    unlink(stamped_out)
  }

  ok <- file.rename(stable_out, stamped_out)
  if (!isTRUE(ok)) {
    ok <- file.copy(stable_out, stamped_out, overwrite = TRUE)
    if (!isTRUE(ok)) {
      stop("Failed to move ", stable_out, " to ", stamped_out)
    }
    unlink(stable_out)
  }

  stamped_out
}

outputs <- character()

# Render HTML first because rmarkdown reuses and refreshes the common
# <article>_files directory. Rendering PDF last preserves the LaTeX figure
# sidecars required to compile the submitted .tex source.
if (isTRUE(render_html)) {
  outputs <- c(outputs, render_one("rjtools::rjournal_article", "html"))
}

if (isTRUE(render_pdf)) {
  outputs <- c(outputs, render_one("rjtools::rjournal_pdf_article", "pdf"))
}

message("Rendered: ", paste(outputs, collapse = ", "))
