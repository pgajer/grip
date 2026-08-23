#!/usr/bin/env Rscript

## Build the HMP-only Illumina 16S graph used by the R Journal paper.
##
## The source table also contains rows from a related U01 cohort and a small
## number of legacy 454 runs. This script retains only rows explicitly labeled
## as HMP and Illumina before feature screening, PCA, and graph construction.

script_path <- if (!is.null(sys.frames()[[1]]$ofile)) {
  normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE)
} else {
  normalizePath("data-raw/hmp_gc.R", winslash = "/", mustWork = FALSE)
}
repo_root <- normalizePath(
  file.path(dirname(script_path), ".."),
  winslash = "/",
  mustWork = FALSE
)

metadata_path <- Sys.getenv(
  "GRIP_HMP_METADATA_TSV",
  unset = ""
)
feature_path <- Sys.getenv(
  "GRIP_HMP_FEATURE_MATRIX_TSV",
  unset = ""
)

if (!nzchar(metadata_path) || !nzchar(feature_path)) {
  stop(
    "Set GRIP_HMP_METADATA_TSV and GRIP_HMP_FEATURE_MATRIX_TSV to the ",
    "combined upstream metadata and feature tables. These source tables are ",
    "not distributed with grip or the paper supplement."
  )
}
for (path in c(metadata_path, feature_path)) {
  if (!file.exists(path)) stop("Missing HMP source file: ", path)
}
if (!requireNamespace("dgraphs", quietly = TRUE)) {
  stop("Package 'dgraphs' is required to build the HMP graph.")
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

required_metadata <- c(
  "canonical_sample_id", "Project", "16S_Platform", "16S_Phase", "CST", "subCST"
)
missing_metadata <- setdiff(required_metadata, names(metadata))
if (length(missing_metadata)) {
  stop("Metadata is missing columns: ", paste(missing_metadata, collapse = ", "))
}
if (!"sample_id" %in% names(features)) {
  stop("Feature matrix must contain a sample_id column.")
}

selected_metadata <- metadata[
  grepl("^HMP", metadata$Project) & metadata$`16S_Platform` == "Illumina",
  ,
  drop = FALSE
]
selected_metadata <- selected_metadata[
  !is.na(selected_metadata$canonical_sample_id) &
    nzchar(selected_metadata$canonical_sample_id),
  ,
  drop = FALSE
]

feature_index <- match(selected_metadata$canonical_sample_id, features$sample_id)
if (anyNA(feature_index)) {
  stop("Some selected HMP samples are absent from the feature matrix.")
}

counts <- as.matrix(features[feature_index, setdiff(names(features), "sample_id"), drop = FALSE])
storage.mode(counts) <- "double"
rownames(counts) <- selected_metadata$canonical_sample_id
if (any(!is.finite(counts)) || any(counts < 0)) {
  stop("HMP feature counts must be finite and non-negative.")
}

prevalence_cutoff <- 0.01
retained_reads_min <- 500
retained_fraction_min <- 0.25
pca_dimensions <- 10L
selected_k <- 3L

detection_count <- colSums(counts > 0)
detection_prevalence <- detection_count / nrow(counts)
keep_features <- detection_prevalence >= prevalence_cutoff
filtered_counts <- counts[, keep_features, drop = FALSE]

original_reads <- rowSums(counts)
retained_reads <- rowSums(filtered_counts)
retained_fraction <- ifelse(original_reads > 0, retained_reads / original_reads, 0)
usable <- retained_reads >= retained_reads_min &
  retained_fraction >= retained_fraction_min

filtered_counts <- filtered_counts[usable, , drop = FALSE]
selected_metadata <- selected_metadata[usable, , drop = FALSE]
original_reads <- original_reads[usable]
retained_reads <- retained_reads[usable]
retained_fraction <- retained_fraction[usable]

relative_abundance <- filtered_counts / rowSums(filtered_counts)
pca_fit <- stats::prcomp(
  relative_abundance,
  center = TRUE,
  scale. = FALSE,
  rank. = min(pca_dimensions, ncol(relative_abundance), nrow(relative_abundance) - 1L)
)
pca_coordinates <- pca_fit$x[, seq_len(min(pca_dimensions, ncol(pca_fit$x))), drop = FALSE]

graph <- dgraphs::create.sknn.graph(
  pca_coordinates,
  k = selected_k,
  connect.components = FALSE,
  neighbor.method = "ann"
)
adj_list <- graph$adj_list
weight_list <- graph$weight_list

connected_components <- function(adj) {
  component <- integer(length(adj))
  component_id <- 0L
  for (start in seq_along(adj)) {
    if (component[[start]] != 0L) next
    component_id <- component_id + 1L
    queue <- start
    component[[start]] <- component_id
    head <- 1L
    while (head <= length(queue)) {
      vertex <- queue[[head]]
      head <- head + 1L
      neighbors <- adj[[vertex]]
      unseen <- neighbors[component[neighbors] == 0L]
      if (length(unseen)) {
        component[unseen] <- component_id
        queue <- c(queue, unseen)
      }
    }
  }
  component
}

component <- connected_components(adj_list)
component_sizes <- table(component)
giant_component_id <- as.integer(names(component_sizes)[which.max(component_sizes)])
giant_vertices <- which(component == giant_component_id)
old_to_new <- integer(length(adj_list))
old_to_new[giant_vertices] <- seq_along(giant_vertices)

gc_adj_list <- vector("list", length(giant_vertices))
gc_weight_list <- vector("list", length(giant_vertices))
for (i in seq_along(giant_vertices)) {
  old_vertex <- giant_vertices[[i]]
  keep <- component[adj_list[[old_vertex]]] == giant_component_id
  gc_adj_list[[i]] <- as.integer(old_to_new[adj_list[[old_vertex]][keep]])
  gc_weight_list[[i]] <- as.numeric(weight_list[[old_vertex]][keep])
}

selected_metadata <- selected_metadata[giant_vertices, , drop = FALSE]
original_reads <- original_reads[giant_vertices]
retained_reads <- retained_reads[giant_vertices]
retained_fraction <- retained_fraction[giant_vertices]

edge_rows <- list()
edge_index <- 0L
for (u in seq_along(gc_adj_list)) {
  keep <- which(gc_adj_list[[u]] > u)
  for (j in keep) {
    edge_index <- edge_index + 1L
    edge_rows[[edge_index]] <- data.frame(
      source = selected_metadata$canonical_sample_id[[u]],
      target = selected_metadata$canonical_sample_id[[gc_adj_list[[u]][[j]]]],
      weight = gc_weight_list[[u]][[j]],
      stringsAsFactors = FALSE
    )
  }
}
edges <- do.call(rbind, edge_rows)

vertex_data <- data.frame(
  vertex_id = selected_metadata$canonical_sample_id,
  project = selected_metadata$Project,
  platform = selected_metadata$`16S_Platform`,
  phase = selected_metadata$`16S_Phase`,
  cst = selected_metadata$CST,
  subcst = selected_metadata$subCST,
  original_reads = as.numeric(original_reads),
  retained_reads = as.numeric(retained_reads),
  retained_fraction = as.numeric(retained_fraction),
  stringsAsFactors = FALSE
)

hmp.gc <- list(
  adj_list = lapply(gc_adj_list, as.integer),
  weight_list = lapply(gc_weight_list, as.numeric),
  vertex_data = vertex_data,
  graph_info = list(
    source_dataset = "Human Microbiome Project",
    assay = "Illumina 16S rRNA amplicon sequencing",
    source_rows = nrow(metadata),
    hmp_illumina_rows = sum(
      grepl("^HMP", metadata$Project) & metadata$`16S_Platform` == "Illumina",
      na.rm = TRUE
    ),
    usable_vertices = sum(usable),
    original_vertices = length(gc_adj_list),
    edge_count = nrow(edges),
    retained_features = sum(keep_features),
    discarded_components = length(component_sizes) - 1L,
    discarded_vertices = sum(component_sizes) - max(component_sizes),
    representation = ">=1% prevalence features; relative abundance; 10-component PCA",
    graph_constructor = "symmetric k-nearest-neighbor graph",
    selected_k = selected_k,
    prevalence_cutoff = prevalence_cutoff,
    retained_reads_min = retained_reads_min,
    retained_fraction_min = retained_fraction_min
  )
)

extdata_dir <- file.path(repo_root, "inst", "extdata", "hmp_gc")
data_dir <- file.path(repo_root, "data")
dir.create(extdata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

write.table(
  edges,
  gzfile(file.path(extdata_dir, "graph_edges.tsv.gz")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  vertex_data,
  gzfile(file.path(extdata_dir, "vertex_metadata.tsv.gz")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

feature_manifest <- data.frame(
  feature = colnames(counts),
  detection_count = as.integer(detection_count),
  detection_prevalence = as.numeric(detection_prevalence),
  retained = as.logical(keep_features),
  stringsAsFactors = FALSE
)
write.table(
  feature_manifest,
  gzfile(file.path(extdata_dir, "feature_manifest.tsv.gz")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

graph_summary <- data.frame(
  quantity = c(
    "source_rows", "hmp_illumina_rows", "usable_vertices",
    "giant_component_vertices", "graph_edges", "discarded_components",
    "retained_features", "selected_k"
  ),
  value = c(
    hmp.gc$graph_info$source_rows,
    hmp.gc$graph_info$hmp_illumina_rows,
    hmp.gc$graph_info$usable_vertices,
    hmp.gc$graph_info$original_vertices,
    hmp.gc$graph_info$edge_count,
    hmp.gc$graph_info$discarded_components,
    hmp.gc$graph_info$retained_features,
    selected_k
  ),
  stringsAsFactors = FALSE
)
write.table(
  graph_summary,
  file.path(extdata_dir, "graph_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

save(hmp.gc, file = file.path(data_dir, "hmp.gc.rda"), compress = "xz")

message(
  "Saved HMP-only graph with ", length(hmp.gc$adj_list), " vertices and ",
  hmp.gc$graph_info$edge_count, " edges."
)
