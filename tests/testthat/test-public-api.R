test_that("the public API matches the documented development contract", {
  expected <- sort(c(
    "build.misf",
    "build.weighted.misf",
    "compare.layouts",
    "cube.asymmetric.cavities.surface.graph",
    "cube.channel.network.surface.graph",
    "cube.periodic.tunnels.surface.graph",
    "cylinder.surface.graph",
    "edge.kk",
    "edge.length.density.stiffness",
    "edge.repulsive.stage",
    "edge.repulsive.state",
    "edges.cube",
    "edges.cycle",
    "edges.cylinder",
    "edges.kary.tree",
    "edges.mesh",
    "edges.path",
    "edges.sierpinski.carpet",
    "edges.sierpinski.tetrahedron",
    "edges.sierpinski.triangle",
    "edges.torus",
    "geodesic.kk",
    "geometry.diagnostics",
    "globalrep.grip",
    "globalrep.weighted.grip",
    "graph.riemannian.star.structure",
    "grip",
    "gripui_app",
    "gripui_family_app",
    "gripui_graph_family_catalog",
    "gripui_project",
    "gripui_project_from_compare",
    "gripui_project_from_dir",
    "gripui_validate_project",
    "irregular.annulus.surface.graph",
    "irregular.ball.solid.graph",
    "irregular.double.torus.surface.graph",
    "irregular.pair.of.pants.surface.graph",
    "irregular.rectangle.surface.graph",
    "irregular.shell.solid.graph",
    "irregular.sphere.surface.graph",
    "irregular.torus.surface.graph",
    "kary.tree.weighted.graph",
    "keep.asymmetric.notches",
    "keep.periodic.holes",
    "keep.slit.channels",
    "keep.staggered.windows",
    "kernel.gram.gkk",
    "landmark.geodesic.kk",
    "legacy.grip",
    "mask.asymmetric.holes",
    "mask.border",
    "mask.corner",
    "mask.cross",
    "mask.cube.asymmetric.cavities",
    "mask.cube.channel.network",
    "mask.cube.periodic.tunnels",
    "mask.tetrahedron.classic",
    "mask.tetrahedron.corner.missing",
    "mask.triangle.bridge",
    "mask.triangle.classic",
    "menger.sponge.surface.graph",
    "mesh.surface.graph",
    "classical.mds",
    "metric.mds",
    "misf.geodesic.kk",
    "occupied.mesh.surface.graph",
    "params.from.summary",
    "plot.layout",
    "prepare.edge.kk",
    "prepare.geodesic.kk",
    "prepare.graph.geodesic.mds",
    "prepare.landmark.geodesic.kk",
    "prepare.misf.geodesic.kk",
    "project.3d",
    "recursive.cube.mask.surface.graph",
    "recursive.mask.grid.surface.graph",
    "recursive.tetrahedron.mask.surface.graph",
    "recursive.triangle.mask.surface.graph",
    "repulsive.stage",
    "repulsive.state",
    "run_gripui",
    "run_gripui_family",
    "sampled.rectangle.surface.graph",
    "sampled.rectangle.surface.graphs",
    "score.geodesic.kk",
    "score.coordinates",
    "score.gmds",
    "score.landmark.geodesic.kk",
    "score.layout",
    "score.surface",
    "score.misf.geodesic.kk",
    "sierpinski.carpet.surface.graph",
    "sierpinski.tetrahedron.surface.graph",
    "sierpinski.triangle.surface.graph",
    "sphere.surface.graph",
    "torus.surface.graph",
    "trace.grip",
    "trace.legacy.grip",
    "triangulated.annulus.surface.graph",
    "triangulated.pair.of.pants.surface.graph",
    "triangulated.polyhedron.surface.graph",
    "vicsek.surface.graph",
    "weighted.grip.nd"
  ))

  expect_length(expected, 104L)
  expect_identical(sort(getNamespaceExports("grip")), expected)
})

test_that("removed compatibility aliases are absent from the namespace", {
  removed <- c(
    "grip.layout",
    "grip.layout.globalrep",
    "grip.layout.globalrep.weighted",
    "grip.layout.legacy",
    "grip.layout.trace",
    "grip.layout.trace.legacy",
    "grip.layout.trace.weighted",
    "grip.layout.weighted",
    "grip.layout.weighted.nd",
    "grip.plot",
    "grip.project.3d",
    "grip.metric.mds.layout",
    "grip.optimize.kernel.gram.gkk.layout",
    "grip.score.gmds.layout",
    "grip.gmds.layout.result",
    "grip.edge.length.density.stiffness",
    "grip.prepare.landmark.geodesic.kk",
    "grip.prepare.geodesic.kk",
    "grip.prepare.graph.geodesic.mds",
    "grip.prepare.misf.geodesic.kk",
    "grip.prepare.edge.kk",
    "grip.score.landmark.geodesic.kk",
    "grip.score.geodesic.kk",
    "grip.score.misf.geodesic.kk",
    "grip.optimize.landmark.geodesic.kk",
    "grip.optimize.geodesic.kk",
    "grip.optimize.misf.geodesic.kk",
    "grip.optimize.edge.kk.layout",
    "grip.optimize.edge.gkk.layout",
    "grip.optimize.edge.isometric.layout",
    "grip.optimize.edge.repulsive.stage",
    "grip.optimize.repulsive.stage",
    "grip.edge.repulsive.state",
    "grip.repulsive.state",
    "grip.build.misf.weighted",
    "grip.build.misf",
    "grip.geometry.diagnostics",
    "grip.params.from.summary",
    "grip.score.layout",
    "grip.compare.layouts"
  )
  ns <- asNamespace("grip")

  expect_length(removed, 40L)
  expect_false(any(vapply(
    removed,
    exists,
    logical(1L),
    envir = ns,
    inherits = FALSE
  )))
})

test_that("internalized graph builders remain available to graph bundles", {
  internal <- c(
    "gmds.result",
    "mesh.surface.embedding",
    "sampled.rectangle.param.coords",
    "edges.occupied.mesh",
    "edges.triangulated.polyhedron"
  )
  ns <- asNamespace("grip")

  expect_true(all(vapply(
    internal,
    exists,
    logical(1L),
    envir = ns,
    inherits = FALSE
  )))
  expect_false(any(internal %in% getNamespaceExports("grip")))

  graph <- mesh.surface.graph(3L, 4L, normalize = "none")
  expect_identical(graph$edges, grip:::edges.mesh(3L, 4L))
  expect_identical(
    graph$coords_surface,
    grip:::mesh.surface.embedding(3L, 4L)
  )
})

test_that("every exported function is called by an associated Rd example", {
  collect_rd_tag <- function(x, tag) {
    found <- list()
    walk <- function(node) {
      if (inherits(node, "Rd") || is.list(node)) {
        if (identical(attr(node, "Rd_tag"), tag)) {
          found[[length(found) + 1L]] <<- node
        } else {
          for (child in node) {
            walk(child)
          }
        }
      }
    }
    walk(x)
    found
  }

  source_man <- testthat::test_path("..", "..", "man")
  if (dir.exists(source_man)) {
    rd <- lapply(
      list.files(source_man, pattern = "[.]Rd$", full.names = TRUE),
      tools::parse_Rd
    )
  } else {
    rd <- unname(tools::Rd_db("grip"))
  }

  topic_info <- lapply(rd, function(topic) {
    aliases <- vapply(
      collect_rd_tag(topic, "\\alias"),
      function(tag) paste0(unlist(tag), collapse = ""),
      character(1L)
    )
    examples <- collect_rd_tag(topic, "\\examples")
    example_text <- paste(
      vapply(
        examples,
        function(tag) paste0(unlist(tag), collapse = ""),
        character(1L)
      ),
      collapse = "\n"
    )
    list(aliases = aliases, example_text = example_text)
  })

  uncovered <- Filter(
    function(export) {
      topics <- Filter(function(topic) export %in% topic$aliases, topic_info)
      example_text <- paste(
        vapply(topics, function(topic) topic$example_text, character(1L)),
        collapse = "\n"
      )
      call_pattern <- paste0(
        "(?<![[:alnum:]_.])",
        gsub("[.]", "\\\\.", export),
        "[[:space:]]*[(]"
      )
      !nzchar(example_text) || !grepl(call_pattern, example_text, perl = TRUE)
    },
    getNamespaceExports("grip")
  )

  if (length(uncovered)) {
    testthat::fail(paste(
      "Exports without explicit Rd example calls:",
      paste(uncovered, collapse = ", ")
    ))
  }
  succeed()
})
