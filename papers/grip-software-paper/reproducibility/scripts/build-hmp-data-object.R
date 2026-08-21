#!/usr/bin/env Rscript

## Build the exact coarsened HMP/U01 graph object used by the paper from the
## tabular inputs distributed in this reproduction supplement.

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (!length(file_arg)) {
  stop("Could not determine the path to this script.")
}

script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
reproduction_root <- normalizePath(
  file.path(script_dir, ".."),
  mustWork = TRUE
)
input_dir <- Sys.getenv(
  "GRIP_HMP_INPUT_DIR",
  unset = file.path(reproduction_root, "hmp_u01_gc_coarse")
)
out.path <- Sys.getenv(
  "GRIP_HMP_DATA_OUTPUT",
  unset = file.path(reproduction_root, "generated", "hmp_u01_gc_coarse.rds")
)

input_paths <- c(
  edges = file.path(input_dir, "coarse_graph_edges.tsv.gz"),
  labels = file.path(input_dir, "coarse_vertex_labels.tsv.gz"),
  membership = file.path(input_dir, "coarse_membership.tsv.gz"),
  rounds = file.path(input_dir, "coarsening_rounds.tsv")
)
missing <- input_paths[!file.exists(input_paths)]
if (length(missing)) {
  stop("Missing HMP reproduction input(s): ", paste(missing, collapse = ", "))
}

read_table <- function(path) {
  utils::read.delim(
    path,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

edges <- read_table(input_paths[["edges"]])
labels <- read_table(input_paths[["labels"]])
members <- read_table(input_paths[["membership"]])
rounds <- read_table(input_paths[["rounds"]])

vertex_ids <- labels$coarse_sample_id
vertex_index <- seq_along(vertex_ids)
names(vertex_index) <- vertex_ids

adj_list <- vector("list", length(vertex_ids))
weight_list <- vector("list", length(vertex_ids))
for (i in seq_len(nrow(edges))) {
  a <- vertex_index[[edges$source[[i]]]]
  b <- vertex_index[[edges$target[[i]]]]
  w <- as.numeric(edges$weight[[i]])
  adj_list[[a]] <- c(adj_list[[a]], b)
  adj_list[[b]] <- c(adj_list[[b]], a)
  weight_list[[a]] <- c(weight_list[[a]], w)
  weight_list[[b]] <- c(weight_list[[b]], w)
}

vertex_data <- labels
names(vertex_data)[names(vertex_data) == "coarse_sample_id"] <- "vertex_id"

hmp.u01.gc.coarse <- list(
  adj_list = lapply(adj_list, as.integer),
  weight_list = lapply(weight_list, as.numeric),
  vertex_data = vertex_data,
  graph_info = list(
    source_dataset = "HMP+U01 16S amplicon",
    representation = ">=1% relative abundance + PCA",
    selected_k = 3L,
    original_vertices = 6474L,
    coarse_vertices = nrow(vertex_data),
    edge_count = nrow(edges),
    discarded_components = 7L,
    coarsening_rounds = rounds
  )
)

dir.create(dirname(out.path), recursive = TRUE, showWarnings = FALSE)
saveRDS(hmp.u01.gc.coarse, out.path, compress = "xz")
message("Saved HMP graph object to ", normalizePath(out.path, mustWork = TRUE))
