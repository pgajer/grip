test_that("GMDS app bundle builds canonical trace bundles for GRIP, GMDS, GKK, and LGKK", {
  compute_bundle <- getFromNamespace("gripui.gmds.compute.bundle", "grip")
  merge_values <- getFromNamespace(".gripui.family.merge.values", "grip")

  catalog <- gripui_graph_family_catalog()
  desc <- catalog$mesh
  values <- merge_values(desc, preset_id = "default")
  values$h <- 4L
  values$w <- 4L
  values$surface <- "paraboloid"
  values$amplitude <- 0.25

  cases <- list(
    grip = list(fit_class = "grip_misf_grip_fit", top_label = "Top-level GRIP solve"),
    gmds = list(fit_class = "grip_misf_gmds_fit", top_label = "Top-level GMDS solve"),
    gkk = list(fit_class = "grip_misf_gkk_fit", top_label = "Top-level GKK solve"),
    lgkk = list(fit_class = "grip_misf_gkk_fit", top_label = "Top-level LGKK solve")
  )

  level_stage_payload <- getFromNamespace("gripui.gmds.level.stage.payload", "grip")
  export_stage_payload <- getFromNamespace("gripui.gmds.export.stage.payload", "grip")
  export_stage_choices <- getFromNamespace("gripui.gmds.export.stage.choices", "grip")
  export_preset_choices <- getFromNamespace("gripui.gmds.export.preset.choices", "grip")
  figure_preset_choices <- getFromNamespace("gripui.gmds.figure.preset.choices", "grip")
  paper_context <- getFromNamespace("gripui.gmds.paper.context", "grip")
  paper_note <- getFromNamespace("gripui.gmds.paper.inline.note", "grip")
  write_static_figure <- getFromNamespace("gripui.gmds.write.static.figure", "grip")
  write_export_bundle <- getFromNamespace("gripui.gmds.write.export.bundle", "grip")
  paper_sync_table <- getFromNamespace("gripui.gmds.paper.sync.table", "grip")

  for (method in names(cases)) {
    bundle <- compute_bundle(
      desc = desc,
      values = values,
      method = method,
      dim = 3L,
      num_init = 6L,
      prepare_seed = 1101L,
      optimizer_seed = 2101L,
      top_level_max_iter = 1L,
      insertion_max_iter = 4L,
      refinement_max_iter = 1L,
      final_polish_max_iter = 1L,
      n_threads = 0L
    )

    expect_true(is.list(bundle))
    expect_equal(bundle$payload$family_id, "mesh")
    expect_equal(bundle$method$id, method)
    expect_s3_class(bundle$fit, cases[[method]]$fit_class)
    expect_true(nrow(bundle$stage_trace) >= 3L)
    expect_true(length(bundle$stage_data) >= 3L)
    expect_true(bundle$prepared$top_level_level >= 0L)
    expect_true(all(c("seed", "top_level", "final_polish") %in% bundle$stage_trace$stage))
    expect_true(is.list(bundle$stage_payloads))
    expect_true(all(c("seed", "initial_placement", "top_level") %in% names(bundle$stage_payloads)))
    expect_true(length(bundle$stage_payloads$seed$active_vertices) >= 3L)
    expect_equal(bundle$stage_payloads$initial_placement$level, bundle$prepared$top_level_level)
    expect_equal(bundle$stage_payloads$top_level$level, bundle$prepared$top_level_level)
    expect_equal(bundle$stage_payloads$top_level$label, cases[[method]]$top_label)

    expansion_levels <- getFromNamespace("gripui.gmds.expansion.levels", "grip")(bundle)
    expect_true(length(expansion_levels) >= 1L)

    insertion_payload <- level_stage_payload(bundle, "insertion", expansion_levels[[1L]])
    refinement_payload <- level_stage_payload(bundle, "refinement", expansion_levels[[1L]])
    expect_false(is.null(insertion_payload))
    expect_false(is.null(refinement_payload))
    expect_equal(insertion_payload$level, expansion_levels[[1L]])
    expect_equal(refinement_payload$level, expansion_levels[[1L]])

    export.choices <- export_stage_choices(bundle)
    expect_true(all(c("reference", "misf", "top_level", "final_polish") %in% unname(export.choices)))
    expect_true(all(c("paper_figure_bundle", "audit_bundle", "tables_only") %in% unname(export_preset_choices())))
    expect_true(all(c("paper_panel", "paper_wide") %in% unname(figure_preset_choices())))
    export.payload <- export_stage_payload(
      bundle = bundle,
      stage_id = "top_level",
      focus_level = bundle$prepared$top_level_level,
      expansion_level = expansion_levels[[1L]]
    )
    expect_true(is.list(export.payload))
    expect_equal(export.payload$label, cases[[method]]$top_label)
    expect_true(is.matrix(export.payload$display_coords))

    context <- paper_context(
      bundle = bundle,
      stage_id = "top_level",
      focus_level = bundle$prepared$top_level_level,
      expansion_level = expansion_levels[[1L]]
    )
    expect_true(is.data.frame(context))
    expect_equal(context$manuscript_section[[1L]], "Coarsest seed, expansion, and refinement")
    expect_match(
      paper_note(bundle, "top_level", focus_level = bundle$prepared$top_level_level),
      "Paper link:"
    )

    png.path <- tempfile(fileext = ".png")
    pdf.path <- tempfile(fileext = ".pdf")
    write_static_figure(
      export_payload = export.payload,
      png_path = png.path,
      pdf_path = pdf.path,
      figure_preset = "paper_panel"
    )
    expect_true(file.exists(png.path))
    expect_true(file.exists(pdf.path))
    expect_gt(file.info(png.path)$size, 0)
    expect_gt(file.info(pdf.path)$size, 0)
    unlink(c(png.path, pdf.path), force = TRUE)

    dir <- tempfile("gmds-export-test-")
    files <- write_export_bundle(
      bundle = bundle,
      export_payload = export.payload,
      stage_id = "top_level",
      focus_level = bundle$prepared$top_level_level,
      expansion_level = expansion_levels[[1L]],
      preset = "paper_figure_bundle",
      figure_preset = "paper_panel",
      dir = dir
    )
    expect_true(all(file.exists(files)))
    expect_true(any(grepl("\\.png$", files)))
    expect_true(any(grepl("\\.pdf$", files)))
    unlink(dir, recursive = TRUE, force = TRUE)
  }

  paper.map <- paper_sync_table()
  expect_true(is.data.frame(paper.map))
  expect_true(nrow(paper.map) >= 5L)
})

test_that("GMDS stage explorer app builds", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("rgl")

  old <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old), add = TRUE)

  app <- gripui_gmds_app(catalog = gripui_graph_family_catalog()[c("mesh", "sampled_rectangle")])
  expect_s3_class(app, "shiny.appobj")
})
