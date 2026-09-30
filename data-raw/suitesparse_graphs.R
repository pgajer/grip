# Invoked by build-suitesparse.py after source and graph validation.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
records <- jsonlite::fromJSON(args[1], simplifyVector = FALSE)
int <- function(x) as.integer(unlist(x))
chr <- function(x) as.character(unlist(x))
graphs <- setNames(lapply(records, function(x) {
  edges <- matrix(int(x$edges), ncol = 2L, byrow = TRUE)
  list(n = as.integer(x$n), edges = edges, edge_weights = rep(1, nrow(edges)),
       vertex_data = data.frame(vertex = seq_len(x$n), source_index = int(x$source_indices),
         label = chr(x$vertex_labels), source_partition = chr(x$partitions),
         original_vertex = int(x$original_indices),
         vertex_id = chr(x$vertex_ids), source_vertex_id = chr(x$source_vertex_ids), component = int(x$components)),
       provenance = x$provenance)
}), vapply(records, `[[`, '', 'id'))
# Preserve the original six-example metadata, reconstructed from archived matrices.
source('data-raw/original_six.R', local = TRUE)
for (name in names(original.six)) {
  old <- original.six[[name]]; id <- old$graph_info$source_name
  stopifnot(identical(graphs[[id]]$edges, old$edges), identical(graphs[[id]]$n, old$n),
            identical(graphs[[id]]$vertex_data$label, old$vertex_data$label))
  graphs[[id]]$provenance$original_example <- old$graph_info
}
catalogue <- do.call(rbind, lapply(records, function(x) data.frame(
  id=x$id, name=sub('^[^/]+/', '', x$id), group=sub('/.*', '', x$id),
  vertices=as.integer(x$n), edges=length(x$edges), conversion_type=x$provenance$conversion_type,
  components=as.integer(x$n_components), original_vertices=as.integer(x$provenance$original_vertices),
  retained_fraction=x$provenance$retained_fraction, gallery=x$gallery, original_six=x$original_six,
  stringsAsFactors=FALSE)))
rownames(catalogue) <- NULL
bundle <- list(schema_version=1L, catalogue=catalogue, graphs=graphs)
file <- 'inst/extdata/suitesparse/graphs.rds'
saveRDS(bundle, file, compress='xz', version=2)
stopifnot(identical(bundle, readRDS(file)))
cat(nrow(catalogue), 'graphs;', file.info(file)$size, 'xz bytes; exact round-trip verified\n')
