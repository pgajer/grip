#!/usr/bin/env Rscript

manual_root <- file.path("dev", "manual")
pdf_dir <- file.path(manual_root, "pdf", "sierpinski-baseline")
tmp_dir <- file.path(manual_root, "tmp", "sierpinski-baseline")
preview_dir <- file.path(tmp_dir, "pdf-previews")
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'devtools' or the 'grip' package to run the Sierpinski baseline.")
}

baseline_profile <- list(
  triangle = list(
    placement = "circle",
    rounds = 64L,
    final_rounds = 128L,
    num_init = 7L,
    num_nbrs = 9L,
    r = 0.15,
    s = 3.0,
    repulsion_factor = 1.0,
    seeds = 1:3
  ),
  carpet = list(
    placement = "barycenter",
    rounds = 64L,
    final_rounds = 128L,
    num_init = 24L,
    num_nbrs = 12L,
    r = 0.15,
    s = 3.0,
    repulsion_factor = 1.0,
    seeds = 1:3
  )
)

build_sierpinski_carpet <- function(level) {
  stopifnot(level >= 0)
  side <- 3L^level
  grid <- expand.grid(x = 0:(side - 1L), y = 0:(side - 1L))

  keep_cell <- function(x, y) {
    while (x > 0L || y > 0L) {
      if ((x %% 3L) == 1L && (y %% 3L) == 1L) {
        return(FALSE)
      }
      x <- x %/% 3L
      y <- y %/% 3L
    }
    TRUE
  }

  keep <- mapply(keep_cell, grid$x, grid$y)
  cells <- grid[keep, , drop = FALSE]
  coords <- cbind(
    x = cells$x + 0.5,
    y = (side - 1L - cells$y) + 0.5
  )

  id_map <- matrix(0L, nrow = side, ncol = side)
  for (i in seq_len(nrow(cells))) {
    id_map[cells$x[i] + 1L, cells$y[i] + 1L] <- i
  }

  edges <- list()
  for (i in seq_len(nrow(cells))) {
    x <- cells$x[i]
    y <- cells$y[i]
    if (x + 1L < side) {
      nbr <- id_map[x + 2L, y + 1L]
      if (nbr > 0L) edges[[length(edges) + 1L]] <- c(i, nbr)
    }
    if (y + 1L < side) {
      nbr <- id_map[x + 1L, y + 2L]
      if (nbr > 0L) edges[[length(edges) + 1L]] <- c(i, nbr)
    }
  }

  list(edges = do.call(rbind, edges), coords = coords)
}

build_sierpinski_triangle <- function(level) {
  stopifnot(level >= 0)

  merge_nodes <- function(edges, from, to) {
    edges[edges == from] <- to
    edges
  }

  build <- function(k) {
    if (k == 0L) {
      coords <- rbind(
        c(0, 0),
        c(1, 0),
        c(0.5, sqrt(3) / 2)
      )
      edges <- rbind(c(1L, 2L), c(2L, 3L), c(3L, 1L))
      return(list(edges = edges, coords = coords, corners = c(1L, 2L, 3L), n = 3L))
    }

    left <- build(k - 1L)
    right <- build(k - 1L)
    top <- build(k - 1L)

    left_coords <- left$coords / 2
    right_coords <- right$coords / 2 + matrix(c(0.5, 0), nrow(right$coords), 2, byrow = TRUE)
    top_coords <- top$coords / 2 + matrix(c(0.25, sqrt(3) / 4), nrow(top$coords), 2, byrow = TRUE)

    off1 <- left$n
    off2 <- left$n + right$n

    edges <- rbind(left$edges,
                   right$edges + off1,
                   top$edges + off2)
    coords <- rbind(left_coords, right_coords, top_coords)

    L <- left$corners
    R <- right$corners + off1
    T <- top$corners + off2

    edges <- merge_nodes(edges, R[1L], L[2L])
    edges <- merge_nodes(edges, T[1L], L[3L])
    edges <- merge_nodes(edges, T[2L], R[3L])

    ids <- sort(unique(c(edges)))
    map <- seq_along(ids)
    names(map) <- ids
    edges <- cbind(map[as.character(edges[, 1L])],
                   map[as.character(edges[, 2L])])
    coords <- coords[ids, , drop = FALSE]
    edges <- t(apply(edges, 1, sort))
    edges <- unique(edges)

    corners <- c(map[as.character(L[1L])],
                 map[as.character(R[2L])],
                 map[as.character(T[3L])])

    list(edges = edges, coords = coords, corners = corners, n = length(ids))
  }

  out <- build(level)
  list(edges = out$edges, coords = out$coords)
}

normalize_coords <- function(coords) {
  centered <- scale(coords, center = TRUE, scale = FALSE)
  radius <- max(sqrt(rowSums(centered^2)))
  if (!is.finite(radius) || radius <= 0) {
    return(centered)
  }
  centered / radius
}

align_to_target <- function(source, target) {
  src <- normalize_coords(source)
  dst <- normalize_coords(target)
  cross <- t(src) %*% dst
  sv <- svd(cross)
  rot <- sv$u %*% t(sv$v)
  aligned <- src %*% rot
  list(
    aligned = aligned,
    target = dst,
    rmse = sqrt(mean(rowSums((aligned - dst)^2)))
  )
}

make_adj_list <- function(edges, n) {
  adj <- vector("list", n)
  for (i in seq_len(nrow(edges))) {
    u <- edges[i, 1L]
    v <- edges[i, 2L]
    adj[[u]] <- c(adj[[u]], v)
    adj[[v]] <- c(adj[[v]], u)
  }
  lapply(adj, as.integer)
}

bfs_distances <- function(adj, source) {
  n <- length(adj)
  dist <- rep.int(-1L, n)
  q <- integer(n)
  head <- 1L
  tail <- 1L
  q[[tail]] <- source
  dist[[source]] <- 0L
  while (head <= tail) {
    v <- q[[head]]
    head <- head + 1L
    for (u in adj[[v]]) {
      if (dist[[u]] == -1L) {
        tail <- tail + 1L
        q[[tail]] <- u
        dist[[u]] <- dist[[v]] + 1L
      }
    }
  }
  dist
}

sample_vertex_pairs <- function(n, sample_size, rng_seed) {
  set.seed(rng_seed)
  out <- matrix(0L, nrow = sample_size, ncol = 2L)
  for (i in seq_len(sample_size)) {
    pair <- sample.int(n, 2L, replace = FALSE)
    out[i, ] <- sort(pair)
  }
  unique(out)
}

edge_length_stats <- function(coords, edges) {
  diffs <- coords[edges[, 1L], , drop = FALSE] - coords[edges[, 2L], , drop = FALSE]
  lengths <- sqrt(rowSums(diffs^2))
  list(
    lengths = lengths,
    median = stats::median(lengths),
    cv = stats::sd(lengths) / mean(lengths)
  )
}

sampled_stress <- function(coords, adj, sample_size = 2000L, rng_seed = 1L) {
  pairs <- sample_vertex_pairs(nrow(coords), sample_size, rng_seed)
  if (nrow(pairs) == 0L) {
    return(NA_real_)
  }

  src_ids <- unique(pairs[, 1L])
  dist_map <- vector("list", length(src_ids))
  names(dist_map) <- as.character(src_ids)
  for (src in src_ids) {
    dist_map[[as.character(src)]] <- bfs_distances(adj, src)
  }

  gd <- integer(nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    gd[[i]] <- dist_map[[as.character(pairs[i, 1L])]][[pairs[i, 2L]]]
  }
  keep <- gd > 0L
  if (!any(keep)) {
    return(NA_real_)
  }

  pairs <- pairs[keep, , drop = FALSE]
  gd <- as.double(gd[keep])
  ed <- sqrt(rowSums((coords[pairs[, 1L], , drop = FALSE] -
                      coords[pairs[, 2L], , drop = FALSE])^2))
  scale_factor <- sum(ed * gd) / sum(gd * gd)
  sqrt(mean((ed - scale_factor * gd)^2))
}

sampled_nonedge_separation_ratio <- function(coords, edges, sample_size = 5000L, rng_seed = 1L) {
  n <- nrow(coords)
  edge_stats <- edge_length_stats(coords, edges)
  edge_keys <- paste(edges[, 1L], edges[, 2L], sep = "-")
  edge_set <- unique(edge_keys)

  set.seed(rng_seed)
  samples <- numeric(0L)
  attempts <- 0L
  max_attempts <- sample_size * 20L
  while (length(samples) < sample_size && attempts < max_attempts) {
    pair <- sort(sample.int(n, 2L, replace = FALSE))
    key <- paste(pair[[1L]], pair[[2L]], sep = "-")
    attempts <- attempts + 1L
    if (key %in% edge_set) {
      next
    }
    diff <- coords[pair[[1L]], ] - coords[pair[[2L]], ]
    samples <- c(samples, sqrt(sum(diff * diff)))
  }

  if (length(samples) == 0L || !is.finite(edge_stats$median) || edge_stats$median <= 0) {
    return(NA_real_)
  }
  min(samples) / edge_stats$median
}

plot_layout_panel <- function(coords, edges, title_text, subtitle_text,
                              xlim = NULL, ylim = NULL) {
  coords <- normalize_coords(coords)
  if (is.null(xlim) || is.null(ylim)) {
    xr <- range(coords[, 1L])
    yr <- range(coords[, 2L])
    xpad <- max(diff(xr) * 0.12, 0.12)
    ypad <- max(diff(yr) * 0.12, 0.12)
    xlim <- xr + c(-xpad, xpad)
    ylim <- yr + c(-ypad, ypad)
  }

  n <- nrow(coords)
  node_cex <- if (n <= 150L) {
    1.0
  } else if (n <= 800L) {
    0.5
  } else {
    0.22
  }
  edge_lwd <- if (n <= 150L) {
    1.4
  } else if (n <= 800L) {
    0.65
  } else {
    0.2
  }

  plot(NA,
       xlim = xlim,
       ylim = ylim,
       asp = 1,
       axes = FALSE,
       xlab = "",
       ylab = "")
  apply(edges, 1, function(e) {
    graphics::segments(coords[e[1L], 1L], coords[e[1L], 2L],
                       coords[e[2L], 1L], coords[e[2L], 2L],
                       col = grDevices::adjustcolor("#0f3b5f", alpha.f = 0.18),
                       lwd = edge_lwd)
  })
  graphics::points(coords[, 1L], coords[, 2L],
                   pch = 21,
                   bg = "#f05a28",
                   col = "#16324f",
                   cex = node_cex,
                   lwd = 0.6)
  graphics::title(main = title_text, sub = subtitle_text,
                  col.main = "#16324f", cex.main = 1.15,
                  col.sub = "#466074", cex.sub = 0.9)
}

write_compare_pdf <- function(path, canonical_coords, layout_coords, edges,
                              title_text, subtitle_text) {
  left_norm <- normalize_coords(canonical_coords)
  right_norm <- normalize_coords(layout_coords)
  xs <- c(left_norm[, 1L], right_norm[, 1L])
  ys <- c(left_norm[, 2L], right_norm[, 2L])
  xr <- range(xs)
  yr <- range(ys)
  xpad <- max(diff(xr) * 0.12, 0.12)
  ypad <- max(diff(yr) * 0.12, 0.12)
  xlim <- xr + c(-xpad, xpad)
  ylim <- yr + c(-ypad, ypad)

  grDevices::pdf(path, width = 15.5, height = 8.5,
                 paper = "special", bg = "#f7f3ea", useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(0, 0, 2.8, 0), xaxs = "i", yaxs = "i")
  plot_layout_panel(canonical_coords, edges, paste(title_text, "- canonical"), subtitle_text,
                    xlim = xlim, ylim = ylim)
  plot_layout_panel(layout_coords, edges, paste(title_text, "- GRIP"), subtitle_text,
                    xlim = xlim, ylim = ylim)
}

render_pdf_previews <- function(pdf_paths, out_dir) {
  pdftoppm <- Sys.which("pdftoppm")
  if (nzchar(pdftoppm)) {
    for (pdf in pdf_paths) {
      prefix <- file.path(out_dir, tools::file_path_sans_ext(basename(pdf)))
      system2(pdftoppm, c("-png", "-f", "1", "-singlefile", pdf, prefix))
    }
    return(invisible(TRUE))
  }

  python_candidates <- unique(c(
    Sys.which("python3"),
    "/Library/Developer/CommandLineTools/usr/bin/python3",
    "/usr/bin/python3"
  ))
  python_candidates <- python_candidates[file.exists(python_candidates)]
  if (length(python_candidates) == 0L) {
    message("Skipping PDF preview rendering: neither pdftoppm nor python3 is available.")
    return(invisible(FALSE))
  }

  script <- paste(
    "import os, sys",
    "try:",
    "    import fitz",
    "except Exception:",
    "    sys.exit(2)",
    "out_dir = sys.argv[1]",
    "for pdf in sys.argv[2:]:",
    "    doc = fitz.open(pdf)",
    "    page = doc.load_page(0)",
    "    pix = page.get_pixmap(matrix=fitz.Matrix(2, 2), alpha=False)",
    "    out = os.path.join(out_dir, os.path.splitext(os.path.basename(pdf))[0] + '.png')",
    "    pix.save(out)",
    sep = "\n"
  )
  script_path <- tempfile(fileext = ".py")
  on.exit(unlink(script_path), add = TRUE)
  writeLines(script, con = script_path)

  for (python in python_candidates) {
    result <- suppressWarnings(system2(python, c(script_path, out_dir, pdf_paths)))
    if (identical(result, 0L)) {
      return(invisible(TRUE))
    }
  }
  message("Skipping PDF preview rendering: no available python3 interpreter could import PyMuPDF.")
  invisible(FALSE)
}

build_graph_specs <- function() {
  triangle_levels <- 4:6
  carpet_levels <- 3:4

  triangle_specs <- lapply(triangle_levels, function(level) {
    built <- build_sierpinski_triangle(level)
    package_edges <- edges.sierpinski.triangle(level)
    if (!identical(unname(built$edges), unname(package_edges))) {
      stop(sprintf("Canonical triangle builder does not match edges.sierpinski.triangle(%d)", level))
    }
    list(family = "triangle", level = level, edges = package_edges, canonical = built$coords)
  })

  carpet_specs <- lapply(carpet_levels, function(level) {
    built <- build_sierpinski_carpet(level)
    package_edges <- edges.sierpinski.carpet(level)
    if (!identical(unname(built$edges), unname(package_edges))) {
      stop(sprintf("Canonical carpet builder does not match edges.sierpinski.carpet(%d)", level))
    }
    list(family = "carpet", level = level, edges = package_edges, canonical = built$coords)
  })

  c(triangle_specs, carpet_specs)
}

run_one_layout <- function(spec, cfg, seed) {
  n <- max(spec$edges)
  adj <- make_adj_list(spec$edges, n)
  started <- proc.time()[["elapsed"]]
  coords <- grip.layout(
    edges = spec$edges,
    n = n,
    dim = 2,
    placement = cfg$placement,
    rounds = cfg$rounds,
    final_rounds = cfg$final_rounds,
    num_init = cfg$num_init,
    num_nbrs = cfg$num_nbrs,
    r = cfg$r,
    s = cfg$s,
    repulsion_factor = cfg$repulsion_factor,
    seed = seed
  )
  elapsed <- proc.time()[["elapsed"]] - started
  aligned <- align_to_target(coords, spec$canonical)
  edge_stats <- edge_length_stats(coords, spec$edges)

  data.frame(
    family = spec$family,
    level = spec$level,
    seed = seed,
    vertices = n,
    edges = nrow(spec$edges),
    placement = cfg$placement,
    rounds = cfg$rounds,
    final_rounds = cfg$final_rounds,
    num_init = cfg$num_init,
    num_nbrs = cfg$num_nbrs,
    r = cfg$r,
    s = cfg$s,
    repulsion_factor = cfg$repulsion_factor,
    procrustes_rmse = aligned$rmse,
    edge_length_cv = edge_stats$cv,
    median_edge_length = edge_stats$median,
    sampled_stress = sampled_stress(coords, adj, sample_size = 2000L, rng_seed = 1000L + seed),
    sampled_nonedge_sep_ratio = sampled_nonedge_separation_ratio(
      coords, spec$edges, sample_size = 5000L, rng_seed = 2000L + seed
    ),
    elapsed_sec = elapsed,
    stringsAsFactors = FALSE
  )
}

format_num <- function(x, digits = 4) {
  ifelse(is.finite(x), format(round(x, digits), nsmall = digits, trim = TRUE), "NA")
}

write_summary_markdown <- function(raw_metrics, summary_metrics, pdf_paths, path) {
  lines <- c(
    "# Sierpinski Baseline",
    "",
    "Baseline profile:",
    "- Triangle: placement=`circle`, rounds=`64`, final_rounds=`128`, num_init=`7`, num_nbrs=`9`, r=`0.15`, s=`3.0`, repulsion_factor=`1.0`, seeds=`1:3`",
    "- Carpet: placement=`barycenter`, rounds=`64`, final_rounds=`128`, num_init=`24`, num_nbrs=`12`, r=`0.15`, s=`3.0`, repulsion_factor=`1.0`, seeds=`1:3`",
    "",
    "Summary by family/level:",
    "",
    "| Family | Level | Vertices | Seeds | Procrustes RMSE | Edge-length CV | Sampled stress | Non-edge sep ratio | Best seed |",
    "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"
  )

  for (i in seq_len(nrow(summary_metrics))) {
    row <- summary_metrics[i, ]
    lines <- c(lines, sprintf(
      "| %s | %d | %d | %d | %s +/- %s | %s +/- %s | %s +/- %s | %s +/- %s | %d |",
      row$family,
      row$level,
      row$vertices,
      row$n_seeds,
      format_num(row$procrustes_rmse_mean),
      format_num(row$procrustes_rmse_sd),
      format_num(row$edge_length_cv_mean),
      format_num(row$edge_length_cv_sd),
      format_num(row$sampled_stress_mean),
      format_num(row$sampled_stress_sd),
      format_num(row$sampled_nonedge_sep_ratio_mean),
      format_num(row$sampled_nonedge_sep_ratio_sd),
      row$best_seed
    ))
  }

  lines <- c(lines, "", "Selected comparison PDFs:", "")
  for (pdf in pdf_paths) {
    lines <- c(lines, sprintf("- `%s`", pdf))
  }

  utils::write.table(lines, file = path, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

run_sierpinski_baseline <- function() {
  graphs <- build_graph_specs()
  raw_metrics <- do.call(
    rbind,
    lapply(graphs, function(spec) {
      cfg <- baseline_profile[[spec$family]]
      do.call(
        rbind,
        lapply(cfg$seeds, function(seed) run_one_layout(spec, cfg, seed))
      )
    })
  )

  raw_metrics <- raw_metrics[order(raw_metrics$family, raw_metrics$level, raw_metrics$seed), , drop = FALSE]

  summary_metrics <- do.call(
    rbind,
    lapply(split(raw_metrics, paste(raw_metrics$family, raw_metrics$level, sep = "-")), function(df) {
      best_idx <- which.min(df$procrustes_rmse)
      data.frame(
        family = df$family[[1L]],
        level = df$level[[1L]],
        vertices = df$vertices[[1L]],
        n_seeds = nrow(df),
        procrustes_rmse_mean = mean(df$procrustes_rmse),
        procrustes_rmse_sd = stats::sd(df$procrustes_rmse),
        edge_length_cv_mean = mean(df$edge_length_cv),
        edge_length_cv_sd = stats::sd(df$edge_length_cv),
        sampled_stress_mean = mean(df$sampled_stress),
        sampled_stress_sd = stats::sd(df$sampled_stress),
        sampled_nonedge_sep_ratio_mean = mean(df$sampled_nonedge_sep_ratio),
        sampled_nonedge_sep_ratio_sd = stats::sd(df$sampled_nonedge_sep_ratio),
        best_seed = df$seed[[best_idx]],
        stringsAsFactors = FALSE
      )
    })
  )

  pdf_paths <- character()
  for (spec in graphs) {
    rows <- raw_metrics$family == spec$family & raw_metrics$level == spec$level
    group <- raw_metrics[rows, , drop = FALSE]
    best_seed <- group$seed[[which.min(group$procrustes_rmse)]]
    cfg <- baseline_profile[[spec$family]]
    coords <- grip.layout(
      edges = spec$edges,
      n = max(spec$edges),
      dim = 2,
      placement = cfg$placement,
      rounds = cfg$rounds,
      final_rounds = cfg$final_rounds,
      num_init = cfg$num_init,
      num_nbrs = cfg$num_nbrs,
      r = cfg$r,
      s = cfg$s,
      repulsion_factor = cfg$repulsion_factor,
      seed = best_seed
    )
    pdf_path <- file.path(
      pdf_dir,
      sprintf("sierpinski-%s-level-%d-baseline-best-seed-%d.pdf",
              spec$family, spec$level, best_seed)
    )
    subtitle <- sprintf(
      "best seed=%d, rounds=%d, final_rounds=%d, num_nbrs=%d, r=%.2f, s=%.1f, repulsion_factor=%.2f",
      best_seed, cfg$rounds, cfg$final_rounds, cfg$num_nbrs, cfg$r, cfg$s,
      cfg$repulsion_factor
    )
    write_compare_pdf(
      pdf_path,
      canonical_coords = spec$canonical,
      layout_coords = coords,
      edges = spec$edges,
      title_text = sprintf("Sierpinski %s level %d", spec$family, spec$level),
      subtitle_text = subtitle
    )
    pdf_paths <- c(pdf_paths, pdf_path)
  }

  csv_path <- file.path(tmp_dir, "sierpinski-baseline-raw-metrics.csv")
  summary_csv_path <- file.path(tmp_dir, "sierpinski-baseline-summary-metrics.csv")
  summary_md_path <- file.path(tmp_dir, "sierpinski-baseline-summary.md")

  utils::write.csv(raw_metrics, csv_path, row.names = FALSE)
  utils::write.csv(summary_metrics, summary_csv_path, row.names = FALSE)
  write_summary_markdown(raw_metrics, summary_metrics, pdf_paths, summary_md_path)
  render_pdf_previews(pdf_paths, preview_dir)

  message(sprintf("Raw metrics written to %s", csv_path))
  message(sprintf("Summary metrics written to %s", summary_csv_path))
  message(sprintf("Markdown summary written to %s", summary_md_path))
  message(sprintf("Comparison PDFs written to %s", pdf_dir))

  invisible(list(
    raw_metrics = raw_metrics,
    summary_metrics = summary_metrics,
    pdf_paths = pdf_paths,
    csv_path = csv_path,
    summary_csv_path = summary_csv_path,
    summary_md_path = summary_md_path
  ))
}

if (sys.nframe() == 0L) {
  run_sierpinski_baseline()
}
