script_path <- if (!is.null(sys.frames()[[1]]$ofile)) {
  normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE)
} else {
  normalizePath("data-raw/hmp_u01_gc_coarse.R", winslash = "/", mustWork = FALSE)
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
default_source_root <- normalizePath(
  file.path(
    repo_root,
    "..",
    "ZB",
    "chm_paper",
    "analysis_output",
    "hmp_u01_gc_coarsened_layout_selection_2026-03-23",
    "coarsening"
  ),
  winslash = "/",
  mustWork = FALSE
)
source_root <- Sys.getenv("GRIP_HMP_U01_COARSE_SOURCE_DIR", unset = default_source_root)
extdata_dir <- file.path(repo_root, "inst", "extdata", "hmp_u01_gc_coarse")
data_dir <- file.path(repo_root, "data")

dir.create(extdata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

source_files <- c(
  "coarse_graph_edges.tsv",
  "coarse_vertex_labels.tsv",
  "coarse_membership.tsv",
  "coarsening_rounds.tsv"
)

for (nm in source_files) {
  src <- file.path(source_root, nm)
  if (!file.exists(src)) stop("Missing source file: ", src)
}

edges <- read.delim(file.path(source_root, "coarse_graph_edges.tsv"), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
labels <- read.delim(file.path(source_root, "coarse_vertex_labels.tsv"), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
members <- read.delim(file.path(source_root, "coarse_membership.tsv"), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
rounds <- read.delim(file.path(source_root, "coarsening_rounds.tsv"), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

gz_path <- function(name) file.path(extdata_dir, paste0(name, ".gz"))
write.table(edges, gzfile(gz_path("coarse_graph_edges.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(labels, gzfile(gz_path("coarse_vertex_labels.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(members, gzfile(gz_path("coarse_membership.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(rounds, file.path(extdata_dir, "coarsening_rounds.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

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

save(hmp.u01.gc.coarse, file = file.path(data_dir, "hmp.u01.gc.coarse.rda"), compress = "xz")
