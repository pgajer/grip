# Rebuild the six derived graph datasets from the archived Matrix Market files.
# Run from the grip source directory; no download is needed.
root <- 'inst/extdata/zheng-graphs'
names <- c('dwt_66', 'lesmis', 'dwt_307', '494_bus', 'dwt_1005', '1138_bus')
groups <- c('HB', 'Newman', 'HB', 'HB', 'HB', 'HB')
zheng.graphs <- setNames(lapply(seq_along(names), function(k) {
  name <- names[k]
  file <- file.path(root, paste0(name, '.mtx.gz'))
  con <- gzfile(file, 'rt'); lines <- readLines(con); close(con)
  stopifnot(grepl('coordinate (real|integer|pattern) symmetric$', lines[1]))
  content <- lines[!grepl('^%', lines)]
  shape <- scan(text = content[1], quiet = TRUE)
  entries <- as.matrix(read.table(text = paste(content[-1], collapse = '\n')))
  if (grepl(' pattern ', lines[1])) entries <- cbind(entries, 1)
  stopifnot(shape[1] == shape[2], nrow(entries) == shape[3], ncol(entries) == 3L)
  keep <- entries[,1] != entries[,2] & entries[,3] != 0
  edges <- cbind(pmin(entries[keep,1], entries[keep,2]), pmax(entries[keep,1], entries[keep,2]))
  edges <- unique(edges); edges <- edges[order(edges[,1], edges[,2]), , drop = FALSE]
  storage.mode(edges) <- 'integer'; dimnames(edges) <- NULL
  labels <- if (name == 'lesmis') readLines(file.path(root, 'lesmis_nodename.txt')) else as.character(seq_len(shape[1]))
  stopifnot(length(labels) == shape[1], !anyNA(edges))
  list(n = as.integer(shape[1]), edges = edges, edge_weights = rep(1, nrow(edges)),
    vertex_data = data.frame(vertex = seq_len(shape[1]), source_index = seq_len(shape[1]), label = labels),
    graph_info = list(name = paste0(name, '_unweighted_graph'), source_name = paste(groups[k], name, sep = '/'),
      source_url = paste0('https://sparse.tamu.edu/', groups[k], '/', name),
      download_url = paste0('https://www.cise.ufl.edu/research/sparse/MM/', groups[k], '/', name, '.tar.gz'),
      license = 'CC BY 4.0', license_url = 'https://creativecommons.org/licenses/by/4.0/',
      source_header = lines[grepl('^%', lines)], source_compressed_md5 = unname(tools::md5sum(file)),
      original_vertices = as.integer(shape[1]), stored_matrix_entries = nrow(entries),
      diagonal_entries_removed = sum(entries[,1] == entries[,2]),
      modifications = 'Nonzero off-diagonal pattern converted to unique undirected edges with unit traversal lengths. All vertices and source indices retained; no component selection or added edges. Original values and metadata retained in the archived matrix.'))
}), names)
save(zheng.graphs, file = 'data/zheng.graphs.rda', compress = 'xz', version = 2)
