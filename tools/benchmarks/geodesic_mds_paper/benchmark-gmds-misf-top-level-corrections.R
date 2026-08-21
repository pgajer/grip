#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) == 0L) {
  stop("Unable to determine script path from commandArgs().")
}
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

args <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% args

run_tag <- if (smoke) {
  sprintf("gmds-misf-top-level-corrections-smoke-%s", format(Sys.Date(), "%Y-%m-%d"))
} else {
  "gmds-misf-top-level-corrections-2026-04-02"
}

design_root <- file.path(repo_root, "output", "geodesic_mds_paper")
tmp_dir <- file.path(design_root, "tmp", run_tag)
pdf_dir <- file.path(design_root, "reports", run_tag)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

tex_path <- file.path(design_root, "reports", "gmds_misf_top_level_corrections_report_2026-04-02.tex")
pdf_path <- file.path(design_root, "reports", "gmds_misf_top_level_corrections_report_2026-04-02.pdf")
rds_path <- file.path(tmp_dir, "gmds_misf_top_level_corrections_results.rds")
metrics_csv <- file.path(tmp_dir, "gmds_misf_top_level_corrections_metrics.csv")
phase_a_rds <- file.path(
  design_root,
  "tmp",
  "gmds-misf-top-level-initializers-2026-04-02",
  "gmds_misf_top_level_initializer_results.rds"
)

if (!file.exists(phase_a_rds)) {
  stop("Expected Phase A bundle at: ", phase_a_rds)
}

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo_root, quiet = TRUE, helpers = FALSE)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(repo_root, quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'pkgload', 'devtools', or the 'grip' package to run this benchmark.")
}

ns <- asNamespace("grip")
align_to_target_nd <- get("grip.align.to.target.nd", envir = ns)

cfg <- list(
  dim = 3L,
  anchor_weight = 0.10,
  anchor_weight_end = 0.02,
  continuation = "linear",
  regularized_max_iter = if (smoke) 4L else 8L,
  pure_max_iter = if (smoke) 4L else 8L,
  trace_frames = if (smoke) 3L else 4L
)

fmt_num <- function(x, digits = 4L) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}

fmt_time <- function(x) {
  ifelse(
    is.finite(x),
    ifelse(x < 1, formatC(x, format = "f", digits = 3L), formatC(x, format = "f", digits = 2L)),
    "--"
  )
}

tex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x, perl = TRUE)
  x
}

center_coords <- function(coords) {
  coords <- as.matrix(coords)
  sweep(coords, 2L, colMeans(coords), "-", check.margin = FALSE)
}

select_trace_indices <- function(n_frames, max_frames = 4L) {
  if (n_frames <= 0L) {
    return(integer(0L))
  }
  max_frames <- max(1L, as.integer(max_frames))
  idx <- unique(round(seq(1, n_frames, length.out = min(max_frames, n_frames))))
  as.integer(idx)
}

extract_gmds_trace_frames <- function(fit, target, max_frames = 4L) {
  if (is.null(fit$frames) || !length(fit$frames)) {
    return(list())
  }
  idx <- select_trace_indices(length(fit$frames), max_frames = max_frames)
  trace_df <- if (is.data.frame(fit$trace)) fit$trace else NULL
  lapply(idx, function(i) {
    frame.coords <- as.matrix(fit$frames[[i]])
    aligned <- align_to_target_nd(frame.coords, target, allow.reflection = TRUE)
    trace_row <- if (!is.null(trace_df) && nrow(trace_df) >= i) trace_df[i, , drop = FALSE] else NULL
    list(
      frame_index = as.integer(i),
      total_frames = as.integer(length(fit$frames)),
      iteration = if (is.null(trace_row) || !"iteration" %in% names(trace_row)) NA_integer_ else as.integer(trace_row$iteration[[1L]]),
      energy = if (is.null(trace_row) || !"energy" %in% names(trace_row)) NA_real_ else as.double(trace_row$energy[[1L]]),
      gmds_energy = if (is.null(trace_row) || !"gmds_energy" %in% names(trace_row)) NA_real_ else as.double(trace_row$gmds_energy[[1L]]),
      anchor_weight = if (is.null(trace_row) || !"anchor_weight" %in% names(trace_row)) NA_real_ else as.double(trace_row$anchor_weight[[1L]]),
      display_coords = aligned$aligned
    )
  })
}

compute_stage_metrics <- function(case,
                                  initializer_id,
                                  initializer_label,
                                  stage_id,
                                  stage_label,
                                  coords,
                                  elapsed_sec = NA_real_,
                                  note = NULL) {
  if (is.null(coords) || !is.matrix(coords) || any(!is.finite(coords))) {
    return(data.frame(
      case_id = case$id,
      case_label = case$label,
      family = case$family,
      top_level = case$top_level,
      top_n = case$top_n,
      initializer_id = initializer_id,
      initializer_label = initializer_label,
      stage_id = stage_id,
      stage_label = stage_label,
      elapsed_sec = as.double(elapsed_sec),
      sigma_geo = NA_real_,
      rho = NA_real_,
      eta = NA_real_,
      alpha_0_05 = NA_real_,
      note = if (is.null(note)) "non-finite coordinates" else as.character(note),
      stringsAsFactors = FALSE
    ))
  }
  score <- grip.score.geodesic.mds(coords = coords, prepared = case$top_prepared)
  aligned <- align_to_target_nd(coords, case$truth, allow.reflection = TRUE)
  display.coords <- aligned$aligned
  data.frame(
    case_id = case$id,
    case_label = case$label,
    family = case$family,
    top_level = case$top_level,
    top_n = case$top_n,
    initializer_id = initializer_id,
    initializer_label = initializer_label,
    stage_id = stage_id,
    stage_label = stage_label,
    elapsed_sec = as.double(elapsed_sec),
    sigma_geo = as.double(score$gmds.stress[[1L]]),
    rho = as.double(aligned$rmse),
    eta = as.double(sample_roughness(display.coords, case$display_adj, case$display_edges)),
    alpha_0_05 = as.double(area_floor_ratio(display.coords, case$display_triangles)),
    note = if (is.null(note)) "" else as.character(note),
    stringsAsFactors = FALSE
  )
}

triangle_areas <- function(coords, triangles) {
  coords <- as.matrix(coords)
  triangles <- as.matrix(triangles)
  if (!nrow(triangles)) {
    return(numeric(0L))
  }
  if (ncol(coords) == 2L) {
    coords <- cbind(coords, 0)
  }
  v1 <- coords[triangles[, 2L], , drop = FALSE] - coords[triangles[, 1L], , drop = FALSE]
  v2 <- coords[triangles[, 3L], , drop = FALSE] - coords[triangles[, 1L], , drop = FALSE]
  cross <- cbind(
    v1[, 2L] * v2[, 3L] - v1[, 3L] * v2[, 2L],
    v1[, 3L] * v2[, 1L] - v1[, 1L] * v2[, 3L],
    v1[, 1L] * v2[, 2L] - v1[, 2L] * v2[, 1L]
  )
  0.5 * sqrt(rowSums(cross^2))
}

area_floor_ratio <- function(coords, triangles) {
  areas <- triangle_areas(coords, triangles)
  if (!length(areas)) {
    return(NA_real_)
  }
  med <- stats::median(areas)
  if (!is.finite(med) || med <= 0) {
    return(NA_real_)
  }
  as.double(stats::quantile(areas, probs = 0.05, names = FALSE)) / med
}

sample_roughness <- function(coords, adj_list, edges) {
  coords <- as.matrix(coords)
  centered <- sweep(coords, 2L, colMeans(coords), FUN = "-", check.margin = FALSE)
  if (!nrow(edges)) {
    return(NA_real_)
  }
  median.edge <- stats::median(sqrt(rowSums(
    (centered[edges[, 1L], , drop = FALSE] - centered[edges[, 2L], , drop = FALSE])^2
  )))
  if (!is.finite(median.edge) || median.edge <= 0) {
    return(NA_real_)
  }
  residuals <- vapply(seq_len(nrow(centered)), function(i) {
    nbrs <- adj_list[[i]]
    if (length(nbrs) == 0L) {
      return(0)
    }
    delta <- centered[i, ] - colMeans(centered[nbrs, , drop = FALSE])
    sum(delta^2)
  }, numeric(1L))
  sqrt(mean(residuals)) / median.edge
}

make_reference_metrics <- function(case) {
  data.frame(
    stage_label = "Reference sample",
    elapsed_sec = NA_real_,
    sigma_geo = 0,
    rho = 0,
    eta = as.double(sample_roughness(case$truth, case$display_adj, case$display_edges)),
    alpha_0_05 = as.double(area_floor_ratio(case$truth, case$display_triangles)),
    note = "Reference active-set sample",
    stringsAsFactors = FALSE
  )
}

run_correction <- function(case, method_result) {
  initializer.coords <- as.matrix(method_result$coords)
  if (any(!is.finite(initializer.coords))) {
    coords <- matrix(NA_real_, nrow = case$top_n, ncol = cfg$dim)
    failure.metrics <- list(
      initializer = compute_stage_metrics(
        case = case,
        initializer_id = method_result$method_id,
        initializer_label = method_result$method_label,
        stage_id = "initializer",
        stage_label = "Initializer",
        coords = coords,
        elapsed_sec = NA_real_,
        note = "initializer unavailable"
      ),
      anchor = compute_stage_metrics(
        case = case,
        initializer_id = method_result$method_id,
        initializer_label = method_result$method_label,
        stage_id = "anchor_regularized",
        stage_label = "Anchor-GMDS",
        coords = coords,
        elapsed_sec = NA_real_,
        note = "initializer unavailable"
      ),
      pure = compute_stage_metrics(
        case = case,
        initializer_id = method_result$method_id,
        initializer_label = method_result$method_label,
        stage_id = "pure_followup",
        stage_label = "Pure follow-up",
        coords = coords,
        elapsed_sec = NA_real_,
        note = "initializer unavailable"
      )
    )
    return(list(
      initializer = list(
        coords = coords,
        display_coords = coords,
        metrics = failure.metrics$initializer
      ),
      anchor = list(
        coords = coords,
        display_coords = coords,
        metrics = failure.metrics$anchor,
        trace_selected = list()
      ),
      pure = list(
        coords = coords,
        display_coords = coords,
        metrics = failure.metrics$pure,
        trace_selected = list()
      )
    ))
  }

  initializer.centered <- center_coords(initializer.coords)

  initializer.metrics <- compute_stage_metrics(
    case = case,
    initializer_id = method_result$method_id,
    initializer_label = method_result$method_label,
    stage_id = "initializer",
    stage_label = "Initializer",
    coords = initializer.centered,
    elapsed_sec = method_result$metrics$elapsed_sec[[1L]],
    note = method_result$metrics$note[[1L]]
  )

  started <- proc.time()[["elapsed"]]
  anchor.fit <- grip.optimize.geodesic.mds(
    coords = initializer.centered,
    prepared = case$top_prepared,
    init = "user",
    anchor_mode = "user",
    anchor_coords = initializer.centered,
    anchor_weight = cfg$anchor_weight,
    anchor_weight_end = cfg$anchor_weight_end,
    continuation = cfg$continuation,
    engine = "cpp",
    n_threads = 0L,
    max_iter = cfg$regularized_max_iter,
    return_trace = TRUE,
    recenter = TRUE
  )
  anchor.elapsed <- proc.time()[["elapsed"]] - started

  started <- proc.time()[["elapsed"]]
  pure.fit <- grip.optimize.geodesic.mds(
    coords = anchor.fit$coords,
    prepared = case$top_prepared,
    init = "user",
    anchor_mode = "none",
    engine = "cpp",
    n_threads = 0L,
    max_iter = cfg$pure_max_iter,
    return_trace = TRUE,
    recenter = TRUE
  )
  pure.elapsed <- proc.time()[["elapsed"]] - started

  anchor.aligned <- align_to_target_nd(anchor.fit$coords, case$truth, allow.reflection = TRUE)
  pure.aligned <- align_to_target_nd(pure.fit$coords, case$truth, allow.reflection = TRUE)

  list(
    initializer = list(
      coords = initializer.centered,
      display_coords = align_to_target_nd(initializer.centered, case$truth, allow.reflection = TRUE)$aligned,
      metrics = initializer.metrics
    ),
    anchor = list(
      coords = anchor.fit$coords,
      display_coords = anchor.aligned$aligned,
      metrics = compute_stage_metrics(
        case = case,
        initializer_id = method_result$method_id,
        initializer_label = method_result$method_label,
        stage_id = "anchor_regularized",
        stage_label = "Anchor-GMDS",
        coords = anchor.fit$coords,
        elapsed_sec = anchor.elapsed,
        note = sprintf(
          "Anchor schedule %.2f -> %.2f over %d iterations",
          cfg$anchor_weight,
          cfg$anchor_weight_end,
          cfg$regularized_max_iter
        )
      ),
      trace_selected = extract_gmds_trace_frames(anchor.fit, target = case$truth, max_frames = cfg$trace_frames)
    ),
    pure = list(
      coords = pure.fit$coords,
      display_coords = pure.aligned$aligned,
      metrics = compute_stage_metrics(
        case = case,
        initializer_id = method_result$method_id,
        initializer_label = method_result$method_label,
        stage_id = "pure_followup",
        stage_label = "Pure follow-up",
        coords = pure.fit$coords,
        elapsed_sec = pure.elapsed,
        note = sprintf(
          "Pure GMDS from anchor-regularized result, %d iterations",
          cfg$pure_max_iter
        )
      ),
      trace_selected = extract_gmds_trace_frames(pure.fit, target = case$truth, max_frames = cfg$trace_frames)
    )
  )
}

safe_run_correction <- function(case, method_result) {
  message("  correction: ", method_result$method_label)
  tryCatch(
    run_correction(case, method_result),
    error = function(e) {
      msg <- conditionMessage(e)
      coords <- matrix(NA_real_, nrow = case$top_n, ncol = cfg$dim)
      list(
        initializer = list(
          coords = coords,
          display_coords = coords,
          metrics = compute_stage_metrics(
            case = case,
            initializer_id = method_result$method_id,
            initializer_label = method_result$method_label,
            stage_id = "initializer",
            stage_label = "Initializer",
            coords = coords,
            elapsed_sec = method_result$metrics$elapsed_sec[[1L]],
            note = paste("failed:", msg)
          )
        ),
        anchor = list(
          coords = coords,
          display_coords = coords,
          metrics = compute_stage_metrics(
            case = case,
            initializer_id = method_result$method_id,
            initializer_label = method_result$method_label,
            stage_id = "anchor_regularized",
            stage_label = "Anchor-GMDS",
            coords = coords,
            elapsed_sec = NA_real_,
            note = paste("failed:", msg)
          ),
          trace_selected = list()
        ),
        pure = list(
          coords = coords,
          display_coords = coords,
          metrics = compute_stage_metrics(
            case = case,
            initializer_id = method_result$method_id,
            initializer_label = method_result$method_label,
            stage_id = "pure_followup",
            stage_label = "Pure follow-up",
            coords = coords,
            elapsed_sec = NA_real_,
            note = paste("failed:", msg)
          ),
          trace_selected = list()
        )
      )
    }
  )
}

save_stage_panel <- function(case_result, stage, file_path, title_line) {
  reference.metrics <- make_reference_metrics(case_result$case)
  entries <- c(
    list(list(
      display_coords = case_result$case$truth,
      metrics = reference.metrics,
      label = "Reference sample"
    )),
    lapply(case_result$methods, function(method) {
      entry <- method[[stage]]
      list(
        display_coords = entry$display_coords,
        metrics = entry$metrics,
        label = entry$metrics$initializer_label[[1L]]
      )
    })
  )

  grDevices::png(file_path, width = 2200, height = 2200, res = 220, bg = "#ffffff")
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old.par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(3L, 3L), mar = c(1.1, 1.1, 3.2, 0.4), oma = c(0, 0, 1.0, 0))

  for (entry in entries) {
    row <- entry$metrics[1L, , drop = FALSE]
    if (all(is.finite(entry$display_coords))) {
      plot.layout(
        coords = entry$display_coords,
        edges = case_result$case$display_edges,
        projection = "ortho",
        azimuth = 35,
        elevation = 24,
        vertex.col = if (identical(entry$label, "Reference sample")) "#bc6c25" else "#3a5a40",
        edge.col = "#adb5bd",
        main = ""
      )
    } else {
      graphics::plot.new()
      graphics::text(0.5, 0.55, labels = entry$label, cex = 1.0, font = 2L)
      graphics::text(0.5, 0.40, labels = "stage failed", cex = 0.95, col = "#8d0801")
    }
    ttl <- if (identical(entry$label, "Reference sample")) {
      "Reference sample"
    } else {
      sprintf(
        "%s\nsigma %s, rho %s",
        row$initializer_label[[1L]],
        fmt_num(row$sigma_geo[[1L]], 4L),
        fmt_num(row$rho[[1L]], 4L)
      )
    }
    graphics::mtext(ttl, side = 3L, line = 0.3, cex = 0.82)
  }
  graphics::mtext(title_line, side = 3L, outer = TRUE, line = -0.3, cex = 1.15, font = 2L)
}

make_case_summary_table <- function(df) {
  rows <- vapply(seq_len(nrow(df)), function(i) {
    paste(
      tex_escape(df$initializer_label[[i]]),
      tex_escape(df$stage_label[[i]]),
      fmt_time(df$elapsed_sec[[i]]),
      fmt_num(df$sigma_geo[[i]], 4L),
      fmt_num(df$rho[[i]], 4L),
      fmt_num(df$eta[[i]], 4L),
      fmt_num(df$alpha_0_05[[i]], 4L),
      sep = " & "
    )
  }, character(1L))
  paste(rows, collapse = " \\\\\n")
}

phase_a_bundle <- readRDS(phase_a_rds)

case_results <- lapply(phase_a_bundle$case_results, function(phase_a_case_result) {
  case <- phase_a_case_result$case
  message("Running Phase B top-level correction panel for: ", case$label)
  methods <- lapply(phase_a_case_result$methods, function(method_result) {
    correction <- safe_run_correction(case, method_result)
    list(
      method_id = method_result$method_id,
      method_label = method_result$method_label,
      initializer = correction$initializer,
      anchor = correction$anchor,
      pure = correction$pure
    )
  })
  metrics <- do.call(rbind, unlist(lapply(methods, function(method) {
    list(method$initializer$metrics, method$anchor$metrics, method$pure$metrics)
  }), recursive = FALSE))
  initializer_figure <- file.path(pdf_dir, sprintf("%s_initializer_grid.png", case$id))
  anchor_figure <- file.path(pdf_dir, sprintf("%s_anchor_regularized_grid.png", case$id))
  pure_figure <- file.path(pdf_dir, sprintf("%s_pure_followup_grid.png", case$id))

  case_result <- list(
    case = case,
    methods = methods,
    metrics = metrics,
    initializer_figure = initializer_figure,
    anchor_figure = anchor_figure,
    pure_figure = pure_figure
  )

  save_stage_panel(
    case_result,
    stage = "initializer",
    file_path = initializer_figure,
    title_line = sprintf("%s: initializer basins", case$label)
  )
  save_stage_panel(
    case_result,
    stage = "anchor",
    file_path = anchor_figure,
    title_line = sprintf("%s: after short anchor-regularized GMDS", case$label)
  )
  save_stage_panel(
    case_result,
    stage = "pure",
    file_path = pure_figure,
    title_line = sprintf("%s: after pure-GMDS follow-up", case$label)
  )

  case_result
})

metrics_df <- do.call(rbind, lapply(case_results, `[[`, "metrics"))
utils::write.csv(metrics_df, metrics_csv, row.names = FALSE)

summary_df <- do.call(rbind, lapply(case_results, function(case_result) {
  anchor.df <- subset(case_result$metrics, stage_id == "anchor_regularized")
  pure.df <- subset(case_result$metrics, stage_id == "pure_followup")
  data.frame(
    case_id = case_result$case$id,
    case_label = case_result$case$label,
    top_level = case_result$case$top_level,
    top_n = case_result$case$top_n,
    best_anchor_sigma = anchor.df$initializer_label[[which.min(anchor.df$sigma_geo)]],
    best_anchor_rho = anchor.df$initializer_label[[which.min(anchor.df$rho)]],
    best_pure_sigma = pure.df$initializer_label[[which.min(pure.df$sigma_geo)]],
    best_pure_rho = pure.df$initializer_label[[which.min(pure.df$rho)]],
    stringsAsFactors = FALSE
  )
}))

bundle <- list(
  run_tag = run_tag,
  generated_at = as.character(Sys.time()),
  cfg = cfg,
  phase_a_rds = phase_a_rds,
  case_results = case_results,
  metrics = metrics_df,
  summary = summary_df,
  output = list(
    tex = tex_path,
    pdf = pdf_path,
    rds = rds_path,
    metrics_csv = metrics_csv
  )
)
saveRDS(bundle, rds_path)

overall_lines <- c(
  "\\documentclass[11pt]{article}",
  "\\usepackage[margin=1in]{geometry}",
  "\\usepackage{graphicx}",
  "\\usepackage{booktabs}",
  "\\usepackage{longtable}",
  "\\usepackage{xcolor}",
  "\\usepackage{hyperref}",
  "\\usepackage{float}",
  "\\hypersetup{colorlinks=true,linkcolor=blue,urlcolor=blue,citecolor=blue}",
  "\\setlength{\\parskip}{0.6em}",
  "\\setlength{\\parindent}{0pt}",
  "\\begin{document}",
  "\\begin{center}",
  "{\\LARGE Phase B: MISF-GMDS Top-Level Initializer Corrections}\\\\[0.5em]",
  "{\\large 2026-04-02}",
  "\\end{center}",
  paste(
    "This report implements Group~B of the MISF-GMDS top-level initializer plan.",
    "Each Phase~A coarse initializer is reused as a fixed top-level basin, corrected by a short anchor-only regularized GMDS stage, and then polished by a short pure-GMDS follow-up on the same complete coarse geodesic graph."
  ),
  sprintf(
    "The correction schedule is held constant across initializers: anchor regularization uses $\\lambda_A = %.2f \\rightarrow %.2f$ with linear continuation over %d iterations, followed by %d pure-GMDS iterations without anchoring.",
    cfg$anchor_weight,
    cfg$anchor_weight_end,
    cfg$regularized_max_iter,
    cfg$pure_max_iter
  ),
  paste(
    "As in Phase~A, all displayed lines come from the fixed active-set parameter-space triangulation rather than the complete coarse graph, so the snapshots emphasize coarse shape recovery instead of graph density."
  ),
  "\\section*{Overall Summary}",
  "\\small",
  "\\begin{longtable}{p{4.4cm}rrp{2.8cm}p{2.8cm}p{2.8cm}p{2.8cm}}",
  "\\toprule",
  "Case & $L$ & $|V_L|$ & Best anchor $\\sigma$ & Best anchor $\\rho$ & Best pure $\\sigma$ & Best pure $\\rho$ \\\\",
  "\\midrule",
  "\\endhead"
)

overall_rows <- vapply(seq_len(nrow(summary_df)), function(i) {
  paste(
    tex_escape(summary_df$case_label[[i]]),
    summary_df$top_level[[i]],
    summary_df$top_n[[i]],
    tex_escape(summary_df$best_anchor_sigma[[i]]),
    tex_escape(summary_df$best_anchor_rho[[i]]),
    tex_escape(summary_df$best_pure_sigma[[i]]),
    tex_escape(summary_df$best_pure_rho[[i]]),
    sep = " & "
  )
}, character(1L))

overall_lines <- c(
  overall_lines,
  paste(overall_rows, collapse = " \\\\\n"),
  " \\\\",
  "\\bottomrule",
  "\\end{longtable}",
  "\\normalsize"
)

for (case_result in case_results) {
  case <- case_result$case
  df <- case_result$metrics
  pure.df <- subset(df, stage_id == "pure_followup")
  df <- df[
    order(
      match(df$initializer_label, unique(df$initializer_label)),
      match(df$stage_id, c("initializer", "anchor_regularized", "pure_followup"))
    ),
    ,
    drop = FALSE
  ]
  best.sigma.idx <- which.min(pure.df$sigma_geo)
  best.rho.idx <- which.min(pure.df$rho)

  overall_lines <- c(
    overall_lines,
    sprintf("\\section*{%s}", tex_escape(case$label)),
    sprintf(
      "Top-level MISF level $L=%d$ with $|V_L|=%d$ active vertices. Under the common Group~B correction schedule, the best final pure-follow-up $\\sigma_{\\mathrm{geo}}$ is achieved by \\textbf{%s}, while the best final Procrustes error $\\rho$ is achieved by \\textbf{%s}. The three figure grids below show the initializer basins, the anchor-regularized intermediate results, and the final pure-GMDS follow-up layouts.",
      case$top_level,
      case$top_n,
      tex_escape(pure.df$initializer_label[[best.sigma.idx]]),
      tex_escape(pure.df$initializer_label[[best.rho.idx]])
    ),
    "\\begin{table}[H]",
    "\\centering",
    "\\begin{tabular}{llrrrrr}",
    "\\toprule",
    "Initializer & Stage & $t$ (s) & $\\sigma$ & $\\rho$ & $\\eta$ & $\\alpha_{0.05}$ \\\\",
    "\\midrule",
    make_case_summary_table(df),
    " \\\\",
    "\\bottomrule",
    "\\end{tabular}",
    sprintf("\\caption{Phase~B top-level correction metrics for %s. Stage runtimes are reported separately for the initializer, anchor-regularized GMDS stage, and pure-GMDS follow-up. Lower $\\sigma$, $\\rho$, and $\\eta$ are better; larger $\\alpha_{0.05}$ is better.}", tex_escape(case$label)),
    "\\end{table}",
    "\\begin{figure}[H]",
    "\\centering",
    sprintf("\\includegraphics[width=0.94\\textwidth]{%s}", normalizePath(case_result$initializer_figure, winslash = "/", mustWork = TRUE)),
    sprintf("\\caption{Reference coarse sample and reused Phase~A initializers for %s.}", tex_escape(case$label)),
    "\\end{figure}",
    "\\begin{figure}[H]",
    "\\centering",
    sprintf("\\includegraphics[width=0.94\\textwidth]{%s}", normalizePath(case_result$anchor_figure, winslash = "/", mustWork = TRUE)),
    sprintf("\\caption{Same initializer order after the short anchor-regularized GMDS correction on %s.}", tex_escape(case$label)),
    "\\end{figure}",
    "\\begin{figure}[H]",
    "\\centering",
    sprintf("\\includegraphics[width=0.94\\textwidth]{%s}", normalizePath(case_result$pure_figure, winslash = "/", mustWork = TRUE)),
    sprintf("\\caption{Same initializer order after the pure-GMDS follow-up on %s.}", tex_escape(case$label)),
    "\\end{figure}"
  )
}

overall_lines <- c(
  overall_lines,
  "\\section*{Interactive Companion}",
  paste(
    "The companion HTML gallery is generated by",
    "\\texttt{tools/reports/geodesic_mds_paper/render-gmds-misf-top-level-corrections-html.R}.",
    "It shows all top-level layouts from this experiment as interactive \\texttt{rglwidget} panels: the reused initializer basins, the anchor-regularized corrections, the pure-GMDS follow-up layouts, and selected correction snapshots for both GMDS stages."
  ),
  "\\end{document}"
)

writeLines(overall_lines, tex_path)

message("Wrote Phase B metrics: ", metrics_csv)
message("Wrote Phase B bundle: ", rds_path)
message("Wrote Phase B report TeX: ", tex_path)
