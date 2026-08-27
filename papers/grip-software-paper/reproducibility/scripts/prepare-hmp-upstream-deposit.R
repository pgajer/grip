#!/usr/bin/env Rscript

## Prepare the minimal UMB-HMP-only inputs needed to rebuild the graph used by
## the grip R Journal paper. The source tables may contain additional cohorts
## and extensive clinical metadata; none of that material is copied.

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (!length(file_arg)) {
  stop("Could not determine the path to this script.")
}

script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
reproducibility_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

metadata_path <- Sys.getenv("GRIP_HMP_COMBINED_METADATA_TSV", unset = "")
feature_path <- Sys.getenv("GRIP_HMP_COMBINED_FEATURE_MATRIX_TSV", unset = "")
output_dir <- Sys.getenv(
  "GRIP_HMP_DEPOSIT_DIR",
  unset = file.path(reproducibility_dir, "hmp_gc", "upstream")
)

if (!nzchar(metadata_path) || !nzchar(feature_path)) {
  stop(
    "Set GRIP_HMP_COMBINED_METADATA_TSV and ",
    "GRIP_HMP_COMBINED_FEATURE_MATRIX_TSV to the source working tables."
  )
}
for (path in c(metadata_path, feature_path)) {
  if (!file.exists(path)) stop("Missing source file: ", path)
}

metadata <- utils::read.delim(
  metadata_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
features <- utils::read.delim(
  feature_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

metadata_columns <- c(
  "canonical_sample_id", "Project", "16S_Platform", "16S_Phase"
)
missing_metadata <- setdiff(metadata_columns, names(metadata))
if (length(missing_metadata)) {
  stop("Metadata is missing columns: ", paste(missing_metadata, collapse = ", "))
}
if (!"sample_id" %in% names(features)) {
  stop("Feature matrix must contain a sample_id column.")
}

keep <- !is.na(metadata$Project) &
  !is.na(metadata$`16S_Platform`) &
  grepl("^HMP", metadata$Project) &
  metadata$`16S_Platform` == "Illumina" &
  !is.na(metadata$canonical_sample_id) &
  nzchar(metadata$canonical_sample_id)

hmp_metadata <- metadata[keep, metadata_columns, drop = FALSE]
if (anyDuplicated(hmp_metadata$canonical_sample_id)) {
  stop("The selected HMP metadata contains duplicate canonical sample IDs.")
}

feature_index <- match(hmp_metadata$canonical_sample_id, features$sample_id)
if (anyNA(feature_index)) {
  stop("Some selected HMP samples are absent from the feature matrix.")
}
hmp_features <- features[feature_index, , drop = FALSE]
if (!identical(hmp_metadata$canonical_sample_id, hmp_features$sample_id)) {
  stop("Metadata and feature rows are not in identical sample order.")
}

count_columns <- setdiff(names(hmp_features), "sample_id")
counts <- as.matrix(hmp_features[, count_columns, drop = FALSE])
storage.mode(counts) <- "double"
if (any(!is.finite(counts)) || any(counts < 0)) {
  stop("Feature counts must be finite and non-negative.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
metadata_output <- file.path(output_dir, "hmp_illumina_metadata.tsv.gz")
feature_output <- file.path(output_dir, "hmp_illumina_feature_counts.tsv.gz")

write.table(
  hmp_metadata,
  gzfile(metadata_output),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  hmp_features,
  gzfile(feature_output),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

sha256 <- function(path) {
  result <- system2(
    "shasum",
    c("-a", "256", normalizePath(path, mustWork = TRUE)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L) {
    stop("Could not calculate SHA-256 for ", path, ": ", paste(result, collapse = "\n"))
  }
  sub("[[:space:]].*$", "", result[[1L]])
}

data_files <- c(metadata_output, feature_output)
manifest <- data.frame(
  file = basename(data_files),
  description = c(
    "UMB-HMP-only Illumina sample identifiers and technical cohort fields",
    "UMB-HMP-only Illumina microbial feature-count matrix"
  ),
  rows = c(nrow(hmp_metadata), nrow(hmp_features)),
  columns = c(ncol(hmp_metadata), ncol(hmp_features)),
  bytes = as.numeric(file.info(data_files)$size),
  sha256 = vapply(data_files, sha256, character(1L)),
  stringsAsFactors = FALSE
)

write.table(
  manifest,
  file.path(output_dir, "MANIFEST.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
writeLines(
  sprintf("%s  %s", manifest$sha256, manifest$file),
  file.path(output_dir, "SHA256SUMS"),
  useBytes = TRUE
)

message(
  "Prepared ", nrow(hmp_metadata), " HMP Illumina samples and ",
  length(count_columns), " microbial features in ", output_dir, "."
)
