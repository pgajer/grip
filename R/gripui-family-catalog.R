.gripui.family.param <- function(id,
                                 label,
                                 type,
                                 default,
                                 group = "Topology",
                                 min = NULL,
                                 max = NULL,
                                 step = NULL,
                                 choices = NULL,
                                 help = NULL,
                                 visible_if = NULL,
                                 advanced = FALSE,
                                 coerce = NULL,
                                 validate = NULL) {
  list(
    id = id,
    label = label,
    type = type,
    default = default,
    group = group,
    min = min,
    max = max,
    step = step,
    choices = choices,
    help = help,
    visible_if = visible_if,
    advanced = advanced,
    coerce = coerce,
    validate = validate
  )
}

.gripui.family.param.int <- function(id, label, default, min, max, step = 1L, ...) {
  .gripui.family.param(
    id = id,
    label = label,
    type = "integer",
    default = as.integer(default),
    min = as.integer(min),
    max = as.integer(max),
    step = as.integer(step),
    coerce = function(x) as.integer(round(as.numeric(x))),
    ...
  )
}

.gripui.family.param.double <- function(id, label, default, min, max, step = 0.1, ...) {
  .gripui.family.param(
    id = id,
    label = label,
    type = "double",
    default = as.numeric(default),
    min = as.numeric(min),
    max = as.numeric(max),
    step = as.numeric(step),
    coerce = function(x) as.numeric(x),
    ...
  )
}

.gripui.family.param.choice <- function(id, label, default, choices, ...) {
  .gripui.family.param(
    id = id,
    label = label,
    type = "choice",
    default = as.character(default),
    choices = as.character(choices),
    coerce = function(x) as.character(x),
    ...
  )
}

.gripui.family.param.logical <- function(id, label, default = FALSE, ...) {
  .gripui.family.param(
    id = id,
    label = label,
    type = "logical",
    default = isTRUE(default),
    coerce = function(x) isTRUE(x),
    ...
  )
}

.gripui.family.parse.numeric.vector <- function(x) {
  if (is.null(x)) {
    return(numeric())
  }
  if (is.numeric(x)) {
    return(as.numeric(x))
  }
  if (length(x) != 1L || is.na(x) || !nzchar(trimws(as.character(x)))) {
    return(numeric())
  }
  pieces <- trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1L]])
  pieces <- pieces[nzchar(pieces)]
  if (length(pieces) == 0L) {
    return(numeric())
  }
  as.numeric(pieces)
}

.gripui.family.param.numeric_vector <- function(id, label, default = numeric(), ...) {
  .gripui.family.param(
    id = id,
    label = label,
    type = "numeric_vector",
    default = as.numeric(default),
    coerce = .gripui.family.parse.numeric.vector,
    ...
  )
}

.gripui.family.desc <- function(id,
                                label,
                                category,
                                function_name,
                                summary,
                                params,
                                builder,
                                presets = NULL,
                                code = NULL,
                                implementation = "/Users/pgajer/current_projects/grip/R/graph_helpers.R",
                                stochastic = FALSE) {
  if (is.null(presets)) {
    presets <- list()
  }
  list(
    id = id,
    label = label,
    category = category,
    function_name = function_name,
    summary = summary,
    params = params,
    builder = builder,
    presets = presets,
    code = code,
    implementation = implementation,
    stochastic = isTRUE(stochastic)
  )
}

.gripui.family.simple.desc <- function(id,
                                       label,
                                       category,
                                       function_name,
                                       summary,
                                       params,
                                       arg_ids = NULL,
                                       presets = NULL,
                                       implementation = "/Users/pgajer/current_projects/grip/R/graph_helpers.R",
                                       stochastic = FALSE) {
  .gripui.family.desc(
    id = id,
    label = label,
    category = category,
    function_name = function_name,
    summary = summary,
    params = params,
    presets = presets,
    implementation = implementation,
    stochastic = stochastic,
    builder = function(p) {
      args <- if (is.null(arg_ids)) p else p[arg_ids]
      builder_fun <- get(function_name, mode = "function")
      do.call(builder_fun, args)
    },
    code = function(p) {
      args <- if (is.null(arg_ids)) p else p[arg_ids]
      .gripui.family.call.code(function_name, args)
    }
  )
}

.gripui.family.param.defaults <- function(desc) {
  out <- lapply(desc$params, `[[`, "default")
  names(out) <- vapply(desc$params, `[[`, character(1L), "id")
  out
}

.gripui.family.preset.values <- function(desc, preset_id = "default") {
  if (is.null(preset_id) || !nzchar(preset_id) || identical(preset_id, "default")) {
    return(list())
  }
  if (!preset_id %in% names(desc$presets)) {
    return(list())
  }
  desc$presets[[preset_id]]
}

.gripui.family.merge.values <- function(desc, preset_id = "default", current = list()) {
  defaults <- .gripui.family.param.defaults(desc)
  preset <- .gripui.family.preset.values(desc, preset_id = preset_id)
  utils::modifyList(utils::modifyList(defaults, preset), current)
}

.gripui.family.call.code <- function(fun_name, args) {
  if (length(args) == 0L) {
    return(sprintf("%s()", fun_name))
  }
  pieces <- mapply(
    function(nm, val) sprintf("%s = %s", nm, paste(deparse(val), collapse = " ")),
    names(args),
    args,
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )
  sprintf("%s(%s)", fun_name, paste(pieces, collapse = ", "))
}

.gripui.family.call.code.parts <- function(fun_name, parts) {
  sprintf(
    "%s(%s)",
    fun_name,
    paste(sprintf("%s = %s", names(parts), unname(parts)), collapse = ", ")
  )
}

.gripui.family.normalize.param <- function() {
  .gripui.family.param.choice(
    "normalize",
    "Normalize weights",
    "median",
    c("median", "mean", "none"),
    group = "Weights"
  )
}

.gripui.family.square.mask <- function(p) {
  if (identical(p$mask_kind, "asymmetric_holes") && (p$k %% 2L) == 0L) {
    stop("Square asymmetric-hole masks require odd k values.", call. = FALSE)
  }
  switch(
    p$mask_kind,
    cross = mask.cross(k = p$k, arm_width = p$arm_width),
    border = mask.border(k = p$k, thickness = p$thickness),
    corner = mask.corner(k = p$k, width = p$width, corner = p$corner),
    asymmetric_holes = mask.asymmetric.holes(k = p$k, hole_size = p$hole_size),
    stop("Unsupported square mask kind: ", p$mask_kind, call. = FALSE)
  )
}

.gripui.family.square.mask.code <- function(p) {
  switch(
    p$mask_kind,
    cross = sprintf("mask.cross(k = %sL, arm_width = %sL)", p$k, p$arm_width),
    border = sprintf("mask.border(k = %sL, thickness = %sL)", p$k, p$thickness),
    corner = sprintf(
      "mask.corner(k = %sL, width = %sL, corner = %s)",
      p$k,
      p$width,
      paste(deparse(p$corner), collapse = " ")
    ),
    asymmetric_holes = sprintf("mask.asymmetric.holes(k = %sL, hole_size = %sL)", p$k, p$hole_size)
  )
}

.gripui.family.keep.matrix <- function(p) {
  switch(
    p$pattern,
    periodic_holes = keep.periodic.holes(
      h = p$h,
      w = p$w,
      hole_period = p$hole_period,
      hole_height = p$hole_height,
      hole_width = p$hole_width,
      row_offset = p$row_offset,
      col_offset = p$col_offset
    ),
    staggered_windows = keep.staggered.windows(
      h = p$h,
      w = p$w,
      window_height = p$window_height,
      window_width = p$window_width,
      row_period = p$row_period,
      col_period = p$col_period,
      row_offset = p$row_offset,
      col_offset = p$col_offset
    ),
    slit_channels = keep.slit.channels(
      h = p$h,
      w = p$w,
      orientation = p$orientation,
      slit_period = p$slit_period,
      slit_width = p$slit_width,
      bridge_spacing = p$bridge_spacing,
      bridge_size = p$bridge_size,
      offset = p$offset
    ),
    asymmetric_notches = keep.asymmetric.notches(
      h = p$h,
      w = p$w,
      notch_depth = p$notch_depth,
      notch_width = p$notch_width
    ),
    stop("Unsupported occupied mesh pattern: ", p$pattern, call. = FALSE)
  )
}

.gripui.family.keep.matrix.code <- function(p) {
  switch(
    p$pattern,
    periodic_holes = sprintf(
      paste0(
        "keep.periodic.holes(h = %sL, w = %sL, hole_period = %sL, ",
        "hole_height = %sL, hole_width = %sL, row_offset = %sL, col_offset = %sL)"
      ),
      p$h, p$w, p$hole_period, p$hole_height, p$hole_width, p$row_offset, p$col_offset
    ),
    staggered_windows = sprintf(
      paste0(
        "keep.staggered.windows(h = %sL, w = %sL, window_height = %sL, ",
        "window_width = %sL, row_period = %sL, col_period = %sL, ",
        "row_offset = %sL, col_offset = %sL)"
      ),
      p$h, p$w, p$window_height, p$window_width, p$row_period, p$col_period, p$row_offset, p$col_offset
    ),
    slit_channels = sprintf(
      paste0(
        "keep.slit.channels(h = %sL, w = %sL, orientation = %s, slit_period = %sL, ",
        "slit_width = %sL, bridge_spacing = %sL, bridge_size = %sL, offset = %sL)"
      ),
      p$h, p$w, paste(deparse(p$orientation), collapse = " "), p$slit_period,
      p$slit_width, p$bridge_spacing, p$bridge_size, p$offset
    ),
    asymmetric_notches = sprintf(
      "keep.asymmetric.notches(h = %sL, w = %sL, notch_depth = %sL, notch_width = %sL)",
      p$h, p$w, p$notch_depth, p$notch_width
    )
  )
}

.gripui.family.triangle.mask <- function(p) {
  switch(
    p$mask_kind,
    classic = mask.triangle.classic(),
    bridge = mask.triangle.bridge(missing = p$missing),
    stop("Unsupported triangle mask kind: ", p$mask_kind, call. = FALSE)
  )
}

.gripui.family.triangle.mask.code <- function(p) {
  switch(
    p$mask_kind,
    classic = "mask.triangle.classic()",
    bridge = sprintf("mask.triangle.bridge(missing = %s)", paste(deparse(p$missing), collapse = " "))
  )
}

.gripui.family.tetrahedron.mask <- function(p) {
  switch(
    p$mask_kind,
    classic = mask.tetrahedron.classic(),
    corner_missing = mask.tetrahedron.corner.missing(omit = p$omit),
    stop("Unsupported tetrahedron mask kind: ", p$mask_kind, call. = FALSE)
  )
}

.gripui.family.tetrahedron.mask.code <- function(p) {
  switch(
    p$mask_kind,
    classic = "mask.tetrahedron.classic()",
    corner_missing = sprintf(
      "mask.tetrahedron.corner.missing(omit = %s)",
      paste(deparse(p$omit), collapse = " ")
    )
  )
}

.gripui.family.cube.mask <- function(p) {
  switch(
    p$mask_kind,
    menger = .menger.sponge.mask(),
    periodic_tunnels = mask.cube.periodic.tunnels(
      side = p$side,
      tunnel_width = p$tunnel_width,
      tunnel_period = p$tunnel_period,
      tunnel_offset = p$tunnel_offset
    ),
    asymmetric_cavities = mask.cube.asymmetric.cavities(
      side = p$side,
      cavity_size = p$cavity_size,
      pocket_size = p$pocket_size
    ),
    channel_network = mask.cube.channel.network(
      side = p$side,
      channel_width = p$channel_width,
      branch_offset = p$branch_offset
    ),
    stop("Unsupported cube mask kind: ", p$mask_kind, call. = FALSE)
  )
}

.gripui.family.cube.mask.code <- function(p) {
  switch(
    p$mask_kind,
    menger = ".menger.sponge.mask()",
    periodic_tunnels = sprintf(
      paste0(
        "mask.cube.periodic.tunnels(side = %sL, tunnel_width = %sL, ",
        "tunnel_period = %sL, tunnel_offset = %sL)"
      ),
      p$side, p$tunnel_width, p$tunnel_period, p$tunnel_offset
    ),
    asymmetric_cavities = sprintf(
      "mask.cube.asymmetric.cavities(side = %sL, cavity_size = %sL, pocket_size = %sL)",
      p$side, p$cavity_size, p$pocket_size
    ),
    channel_network = sprintf(
      "mask.cube.channel.network(side = %sL, channel_width = %sL, branch_offset = %sL)",
      p$side, p$channel_width, p$branch_offset
    )
  )
}

.gripui.family.surface2d.params <- function(surface_choices,
                                            amplitude_default = 0.75,
                                            freq_u_default = 1L,
                                            freq_v_default = 1L,
                                            include_scales = TRUE) {
  out <- list(
    .gripui.family.param.choice("surface", "Surface", surface_choices[[1L]], surface_choices, group = "Geometry"),
    .gripui.family.param.double("amplitude", "Amplitude", amplitude_default, 0, 2, 0.05, group = "Geometry"),
    .gripui.family.param.int("freq_u", "Frequency u", freq_u_default, 1L, 8L, group = "Geometry"),
    .gripui.family.param.int("freq_v", "Frequency v", freq_v_default, 1L, 8L, group = "Geometry")
  )
  if (isTRUE(include_scales)) {
    out <- c(
      out,
      list(
        .gripui.family.param.double("x_scale", "x scale", 1, 0.25, 4, 0.05, group = "Geometry"),
        .gripui.family.param.double("y_scale", "y scale", 1, 0.25, 4, 0.05, group = "Geometry")
      )
    )
  }
  c(out, list(.gripui.family.normalize.param()))
}

.gripui.family.surface3d.params <- function(surface_choices,
                                            amplitude_default = 0.2,
                                            freq_default = 2L,
                                            twist_default = 0.6,
                                            include_xyz_scales = FALSE) {
  out <- list(
    .gripui.family.param.choice("surface", "Surface", surface_choices[[1L]], surface_choices, group = "Geometry"),
    .gripui.family.param.double("amplitude", "Amplitude", amplitude_default, 0, 1.5, 0.05, group = "Geometry"),
    .gripui.family.param.int("freq", "Frequency", freq_default, 1L, 8L, group = "Geometry"),
    .gripui.family.param.double("twist", "Twist", twist_default, 0, 2, 0.05, group = "Geometry")
  )
  if (isTRUE(include_xyz_scales)) {
    out <- c(
      out,
      list(
        .gripui.family.param.double("x_scale", "x scale", 1, 0.25, 4, 0.05, group = "Geometry"),
        .gripui.family.param.double("y_scale", "y scale", 1, 0.25, 4, 0.05, group = "Geometry"),
        .gripui.family.param.double("z_scale", "z scale", 1, 0.25, 4, 0.05, group = "Geometry")
      )
    )
  }
  c(out, list(.gripui.family.normalize.param()))
}

#' Catalog of graph families for the geometry explorer app
#'
#' The catalog is a registry used by `gripui_family_app()` to populate its
#' family selector, presets, parameter controls, builder calls, and source
#' references.
#'
#' @return A named list of family descriptors.
#' @export
#'
#' @examples
#' catalog <- gripui_graph_family_catalog()
#' names(catalog)
#' catalog$mesh$function_name
gripui_graph_family_catalog <- function() {
  graph_impl <- "/Users/pgajer/current_projects/grip/R/graph_helpers.R"

  list(
    mesh = .gripui.family.simple.desc(
      id = "mesh",
      label = "Mesh surface",
      category = "Regular lifted lattices",
      function_name = "mesh.surface.graph",
      summary = "Regular rectangular mesh lifted into a smooth 3D surface.",
      implementation = graph_impl,
      arg_ids = c("h", "w", "surface", "amplitude", "freq_u", "freq_v", "x_scale", "y_scale", "normalize"),
      params = c(
        list(
          .gripui.family.param.int("h", "Height", 6L, 2L, 24L),
          .gripui.family.param.int("w", "Width", 6L, 2L, 24L)
        ),
        .gripui.family.surface2d.params(c("saddle", "paraboloid", "ripple"))
      ),
      presets = list(
        paraboloid = list(surface = "paraboloid", amplitude = 0.55),
        ripple = list(surface = "ripple", amplitude = 0.55, freq_u = 2L, freq_v = 2L)
      )
    ),
    irregular_rectangle = .gripui.family.simple.desc(
      id = "irregular_rectangle",
      label = "Irregular rectangle surface",
      category = "Irregular lifted lattices",
      function_name = "irregular.rectangle.surface.graph",
      summary = "Simply connected rectangular mesh with deterministic irregular planar spacing and smooth lifted-surface geometry.",
      implementation = graph_impl,
      arg_ids = c(
        "h", "w", "surface", "amplitude", "freq_u", "freq_v", "x_scale", "y_scale",
        "row_irregularity", "col_irregularity", "row_phase", "col_phase",
        "interior_warp", "shear", "min_step_ratio", "connectivity", "normalize"
      ),
      params = c(
        list(
          .gripui.family.param.int("h", "Height", 8L, 2L, 32L),
          .gripui.family.param.int("w", "Width", 8L, 2L, 32L),
          .gripui.family.param.choice(
            "connectivity",
            "Connectivity",
            "orthogonal",
            c("orthogonal", "diagonal"),
            group = "Topology"
          ),
          .gripui.family.param.double(
            "row_irregularity",
            "Row irregularity",
            0.20,
            0,
            0.9,
            0.05,
            group = "Irregularity",
            help = "Controls deterministic nonuniform spacing along rows."
          ),
          .gripui.family.param.double(
            "col_irregularity",
            "Column irregularity",
            0.20,
            0,
            0.9,
            0.05,
            group = "Irregularity",
            help = "Controls deterministic nonuniform spacing along columns."
          ),
          .gripui.family.param.double("row_phase", "Row phase", 0.35, 0, 2, 0.05, group = "Irregularity"),
          .gripui.family.param.double("col_phase", "Column phase", 0.65, 0, 2, 0.05, group = "Irregularity"),
          .gripui.family.param.double(
            "interior_warp",
            "Interior warp",
            0.08,
            0,
            0.5,
            0.02,
            group = "Irregularity",
            help = "Smooth boundary-vanishing warp applied before lifting."
          ),
          .gripui.family.param.double(
            "shear",
            "Shear",
            0,
            -0.75,
            0.75,
            0.05,
            group = "Irregularity"
          ),
          .gripui.family.param.double(
            "min_step_ratio",
            "Min step ratio",
            0.30,
            0.05,
            0.95,
            0.05,
            group = "Irregularity",
            help = "Lower bound on local row/column spacing after irregularization."
          )
        ),
        .gripui.family.surface2d.params(
          c("paraboloid", "saddle", "ripple", "flat"),
          amplitude_default = 0.75
        )
      ),
      presets = list(
        square = list(h = 8L, w = 8L),
        flat = list(surface = "flat", amplitude = 0),
        more_irregular = list(
          row_irregularity = 0.35,
          col_irregularity = 0.35,
          interior_warp = 0.14
        ),
        diagonal = list(connectivity = "diagonal")
      )
    ),
    sampled_rectangle = .gripui.family.simple.desc(
      id = "sampled_rectangle",
      label = "Sampled rectangle iKNN surface",
      category = "Sampled lifted surfaces",
      function_name = "sampled.rectangle.surface.graph",
      summary = "Uniformly sampled rectangle lifted into 3D, with topology built from an exact iKNN graph on the sampled points.",
      implementation = graph_impl,
      stochastic = TRUE,
      arg_ids = c(
        "n", "k", "xmin", "xmax", "ymin", "ymax", "seed",
        "surface", "amplitude", "freq_u", "freq_v",
        "graph_space", "normalize"
      ),
      params = c(
        list(
          .gripui.family.param.int("n", "Samples", 80L, 8L, 400L, group = "Topology"),
          .gripui.family.param.int("k", "k", 6L, 1L, 40L, group = "Topology"),
          .gripui.family.param.choice(
            "graph_space",
            "Graph space",
            "surface",
            c("surface", "param"),
            group = "Topology",
            help = "Build the iKNN graph in the 3D embedding or in the planar rectangle."
          ),
          .gripui.family.param.double("xmin", "xmin", -1, -4, 4, 0.1, group = "Rectangle"),
          .gripui.family.param.double("xmax", "xmax", 1, -4, 4, 0.1, group = "Rectangle"),
          .gripui.family.param.double("ymin", "ymin", -1, -4, 4, 0.1, group = "Rectangle"),
          .gripui.family.param.double("ymax", "ymax", 1, -4, 4, 0.1, group = "Rectangle"),
          .gripui.family.param.int("seed", "Seed", 1L, 0L, 1000000L, group = "Rectangle")
        ),
        .gripui.family.surface2d.params(
          c("flat", "saddle", "paraboloid", "ripple", "folded"),
          amplitude_default = 0.75,
          include_scales = FALSE
        )
      ),
      presets = list(
        paraboloid = list(surface = "paraboloid", amplitude = 0.7),
        ripple = list(surface = "ripple", amplitude = 0.55, freq_u = 2L, freq_v = 2L),
        planar_graph = list(surface = "flat", amplitude = 0, graph_space = "param"),
        wide_rectangle = list(xmin = -2, xmax = 2, ymin = -1, ymax = 1, n = 120L)
      )
    ),
    cylinder = .gripui.family.simple.desc(
      id = "cylinder",
      label = "Cylinder surface",
      category = "Regular lifted lattices",
      function_name = "cylinder.surface.graph",
      summary = "Wrapped cylindrical grid with barrel, hourglass, or wavy radial profiles.",
      implementation = graph_impl,
      arg_ids = c("h", "w", "surface", "radius", "height", "amplitude", "freq_theta", "freq_z", "twist", "normalize"),
      params = list(
        .gripui.family.param.int("h", "Axial samples", 7L, 2L, 24L),
        .gripui.family.param.int("w", "Wrapped samples", 14L, 3L, 48L),
        .gripui.family.param.choice("surface", "Surface", "barrel", c("barrel", "hourglass", "wavy"), group = "Geometry"),
        .gripui.family.param.double("radius", "Radius", 1, 0.2, 4, 0.05, group = "Geometry"),
        .gripui.family.param.double("height", "Height", 2.4, 0.5, 6, 0.05, group = "Geometry"),
        .gripui.family.param.double("amplitude", "Amplitude", 0.28, 0, 1.5, 0.02, group = "Geometry"),
        .gripui.family.param.int("freq_theta", "Angular frequency", 2L, 1L, 8L, group = "Geometry"),
        .gripui.family.param.int("freq_z", "Axial frequency", 1L, 1L, 8L, group = "Geometry"),
        .gripui.family.param.double("twist", "Twist", 0.25, 0, 2, 0.05, group = "Geometry"),
        .gripui.family.normalize.param()
      ),
      presets = list(
        hourglass = list(surface = "hourglass"),
        wavy = list(surface = "wavy", amplitude = 0.24, freq_theta = 3L, freq_z = 2L)
      )
    ),
    torus = .gripui.family.simple.desc(
      id = "torus",
      label = "Torus surface",
      category = "Regular lifted lattices",
      function_name = "torus.surface.graph",
      summary = "Wrapped toroidal grid with standard, pinched, or wavy tube geometry.",
      implementation = graph_impl,
      arg_ids = c("h", "w", "surface", "major_radius", "minor_radius", "amplitude", "freq_major", "freq_minor", "twist", "normalize"),
      params = list(
        .gripui.family.param.int("h", "Major-ring samples", 8L, 3L, 24L),
        .gripui.family.param.int("w", "Tube samples", 14L, 3L, 48L),
        .gripui.family.param.choice("surface", "Surface", "standard", c("standard", "pinched", "wavy"), group = "Geometry"),
        .gripui.family.param.double("major_radius", "Major radius", 2, 0.5, 6, 0.05, group = "Geometry"),
        .gripui.family.param.double("minor_radius", "Minor radius", 0.75, 0.1, 3, 0.05, group = "Geometry"),
        .gripui.family.param.double("amplitude", "Amplitude", 0.2, 0, 1, 0.02, group = "Geometry"),
        .gripui.family.param.int("freq_major", "Major frequency", 2L, 1L, 8L, group = "Geometry"),
        .gripui.family.param.int("freq_minor", "Minor frequency", 1L, 1L, 8L, group = "Geometry"),
        .gripui.family.param.double("twist", "Twist", 0.25, 0, 2, 0.05, group = "Geometry"),
        .gripui.family.normalize.param()
      ),
      presets = list(
        pinched = list(surface = "pinched"),
        wavy = list(surface = "wavy", amplitude = 0.18, freq_major = 3L, freq_minor = 2L)
      )
    ),
    sphere = .gripui.family.simple.desc(
      id = "sphere",
      label = "Sphere surface",
      category = "Regular lifted lattices",
      function_name = "sphere.surface.graph",
      summary = "Regular spherical surface family with ellipsoid and wavy variants.",
      implementation = graph_impl,
      arg_ids = c("h", "w", "surface", "radius", "amplitude", "freq_theta", "freq_lat", "twist", "normalize"),
      params = list(
        .gripui.family.param.int("h", "Latitude bands", 7L, 3L, 24L),
        .gripui.family.param.int("w", "Band samples", 14L, 3L, 48L),
        .gripui.family.param.choice("surface", "Surface", "standard", c("standard", "ellipsoid", "wavy"), group = "Geometry"),
        .gripui.family.param.double("radius", "Radius", 1, 0.2, 4, 0.05, group = "Geometry"),
        .gripui.family.param.double("amplitude", "Amplitude", 0.2, 0, 1, 0.02, group = "Geometry"),
        .gripui.family.param.int("freq_theta", "Longitudinal frequency", 3L, 1L, 8L, group = "Geometry"),
        .gripui.family.param.int("freq_lat", "Latitudinal frequency", 2L, 1L, 8L, group = "Geometry"),
        .gripui.family.param.double("twist", "Twist", 0.25, 0, 2, 0.05, group = "Geometry"),
        .gripui.family.normalize.param()
      ),
      presets = list(
        ellipsoid = list(surface = "ellipsoid"),
        wavy = list(surface = "wavy", amplitude = 0.18, freq_theta = 4L, freq_lat = 3L)
      )
    ),
    recursive_mask_grid = .gripui.family.desc(
      id = "recursive_mask_grid",
      label = "Recursive square-mask grid",
      category = "Recursive square masks",
      function_name = "recursive.mask.grid.surface.graph",
      summary = "Generic recursive square-mask family for carpet-like and asymmetric grid fractals.",
      implementation = graph_impl,
      params = c(
        list(
          .gripui.family.param.choice("mask_kind", "Mask pattern", "cross", c("cross", "border", "corner", "asymmetric_holes")),
          .gripui.family.param.int("k", "Mask side", 5L, 3L, 11L),
          .gripui.family.param.int("arm_width", "Arm width", 1L, 1L, 4L, visible_if = list(mask_kind = "cross")),
          .gripui.family.param.int("thickness", "Border thickness", 1L, 1L, 4L, visible_if = list(mask_kind = "border")),
          .gripui.family.param.int("width", "Corner width", 2L, 1L, 5L, visible_if = list(mask_kind = "corner")),
          .gripui.family.param.choice(
            "corner",
            "Corner",
            "top_left",
            c("top_left", "top_right", "bottom_left", "bottom_right"),
            visible_if = list(mask_kind = "corner")
          ),
          .gripui.family.param.int("hole_size", "Hole size", 1L, 1L, 3L, visible_if = list(mask_kind = "asymmetric_holes")),
          .gripui.family.param.int("level", "Recursion level", 2L, 1L, 4L)
        ),
        .gripui.family.surface2d.params(c("saddle", "paraboloid", "ripple"))
      ),
      presets = list(
        border = list(mask_kind = "border"),
        corner = list(mask_kind = "corner", corner = "bottom_right"),
        asymmetric = list(mask_kind = "asymmetric_holes", k = 5L, hole_size = 1L)
      ),
      builder = function(p) {
        args <- list(
          mask = .gripui.family.square.mask(p),
          level = p$level,
          surface = p$surface,
          amplitude = p$amplitude,
          freq_u = p$freq_u,
          freq_v = p$freq_v,
          x_scale = p$x_scale,
          y_scale = p$y_scale,
          normalize = p$normalize
        )
        do.call(recursive.mask.grid.surface.graph, args)
      },
      code = function(p) {
        .gripui.family.call.code.parts("recursive.mask.grid.surface.graph", c(
          mask = .gripui.family.square.mask.code(p),
          level = paste(deparse(p$level), collapse = " "),
          surface = paste(deparse(p$surface), collapse = " "),
          amplitude = paste(deparse(p$amplitude), collapse = " "),
          freq_u = paste(deparse(p$freq_u), collapse = " "),
          freq_v = paste(deparse(p$freq_v), collapse = " "),
          x_scale = paste(deparse(p$x_scale), collapse = " "),
          y_scale = paste(deparse(p$y_scale), collapse = " "),
          normalize = paste(deparse(p$normalize), collapse = " ")
        ))
      }
    ),
    sierpinski_carpet = .gripui.family.simple.desc(
      id = "sierpinski_carpet",
      label = "Sierpinski carpet",
      category = "Recursive square masks",
      function_name = "sierpinski.carpet.surface.graph",
      summary = "Classic recursive carpet derived from the square-mask grid family.",
      implementation = graph_impl,
      arg_ids = c("level", "surface", "amplitude", "freq_u", "freq_v", "x_scale", "y_scale", "normalize"),
      params = c(
        list(.gripui.family.param.int("level", "Recursion level", 2L, 1L, 4L)),
        .gripui.family.surface2d.params(c("saddle", "paraboloid", "ripple"))
      ),
      presets = list(
        paraboloid = list(surface = "paraboloid"),
        ripple = list(surface = "ripple", amplitude = 0.55, freq_u = 2L, freq_v = 2L)
      )
    ),
    vicsek = .gripui.family.simple.desc(
      id = "vicsek",
      label = "Vicsek fractal",
      category = "Recursive square masks",
      function_name = "vicsek.surface.graph",
      summary = "Recursive cross-mask grid family with strong axial bottlenecks.",
      implementation = graph_impl,
      arg_ids = c("level", "surface", "amplitude", "freq_u", "freq_v", "x_scale", "y_scale", "normalize"),
      params = c(
        list(.gripui.family.param.int("level", "Recursion level", 2L, 1L, 4L)),
        .gripui.family.surface2d.params(c("saddle", "paraboloid", "ripple"))
      ),
      presets = list(
        ripple = list(surface = "ripple", amplitude = 0.6, freq_u = 2L, freq_v = 2L)
      )
    ),
    occupied_mesh = .gripui.family.desc(
      id = "occupied_mesh",
      label = "Perforated occupied mesh",
      category = "Recursive square masks",
      function_name = "occupied.mesh.surface.graph",
      summary = "Finite perforated mesh family built from deterministic keep-pattern constructors.",
      implementation = graph_impl,
      params = c(
        list(
          .gripui.family.param.choice(
            "pattern",
            "Keep pattern",
            "periodic_holes",
            c("periodic_holes", "staggered_windows", "slit_channels", "asymmetric_notches")
          ),
          .gripui.family.param.int("h", "Height", 12L, 4L, 48L),
          .gripui.family.param.int("w", "Width", 14L, 4L, 48L),
          .gripui.family.param.int("hole_period", "Hole period", 4L, 2L, 10L, visible_if = list(pattern = "periodic_holes")),
          .gripui.family.param.int("hole_height", "Hole height", 1L, 1L, 6L, visible_if = list(pattern = "periodic_holes")),
          .gripui.family.param.int("hole_width", "Hole width", 1L, 1L, 6L, visible_if = list(pattern = "periodic_holes")),
          .gripui.family.param.int("row_offset", "Row offset", 2L, 1L, 10L, visible_if = list(pattern = c("periodic_holes", "staggered_windows"))),
          .gripui.family.param.int("col_offset", "Column offset", 2L, 1L, 10L, visible_if = list(pattern = c("periodic_holes", "staggered_windows"))),
          .gripui.family.param.int("window_height", "Window height", 1L, 1L, 6L, visible_if = list(pattern = "staggered_windows")),
          .gripui.family.param.int("window_width", "Window width", 2L, 1L, 8L, visible_if = list(pattern = "staggered_windows")),
          .gripui.family.param.int("row_period", "Row period", 4L, 2L, 10L, visible_if = list(pattern = "staggered_windows")),
          .gripui.family.param.int("col_period", "Column period", 5L, 2L, 10L, visible_if = list(pattern = "staggered_windows")),
          .gripui.family.param.choice("orientation", "Slit orientation", "vertical", c("vertical", "horizontal"), visible_if = list(pattern = "slit_channels")),
          .gripui.family.param.int("slit_period", "Slit period", 5L, 2L, 10L, visible_if = list(pattern = "slit_channels")),
          .gripui.family.param.int("slit_width", "Slit width", 1L, 1L, 4L, visible_if = list(pattern = "slit_channels")),
          .gripui.family.param.int("bridge_spacing", "Bridge spacing", 4L, 2L, 10L, visible_if = list(pattern = "slit_channels")),
          .gripui.family.param.int("bridge_size", "Bridge size", 1L, 1L, 4L, visible_if = list(pattern = "slit_channels")),
          .gripui.family.param.int("offset", "Slit offset", 2L, 1L, 10L, visible_if = list(pattern = "slit_channels")),
          .gripui.family.param.int("notch_depth", "Notch depth", 3L, 1L, 10L, visible_if = list(pattern = "asymmetric_notches")),
          .gripui.family.param.int("notch_width", "Notch width", 2L, 1L, 10L, visible_if = list(pattern = "asymmetric_notches"))
        ),
        .gripui.family.surface2d.params(c("saddle", "paraboloid", "ripple"))
      ),
      presets = list(
        windows = list(pattern = "staggered_windows"),
        slits = list(pattern = "slit_channels"),
        notches = list(pattern = "asymmetric_notches")
      ),
      builder = function(p) {
        args <- list(
          keep = .gripui.family.keep.matrix(p),
          surface = p$surface,
          amplitude = p$amplitude,
          freq_u = p$freq_u,
          freq_v = p$freq_v,
          x_scale = p$x_scale,
          y_scale = p$y_scale,
          normalize = p$normalize
        )
        do.call(occupied.mesh.surface.graph, args)
      },
      code = function(p) {
        .gripui.family.call.code.parts("occupied.mesh.surface.graph", c(
          keep = .gripui.family.keep.matrix.code(p),
          surface = paste(deparse(p$surface), collapse = " "),
          amplitude = paste(deparse(p$amplitude), collapse = " "),
          freq_u = paste(deparse(p$freq_u), collapse = " "),
          freq_v = paste(deparse(p$freq_v), collapse = " "),
          x_scale = paste(deparse(p$x_scale), collapse = " "),
          y_scale = paste(deparse(p$y_scale), collapse = " "),
          normalize = paste(deparse(p$normalize), collapse = " ")
        ))
      }
    ),
    recursive_triangle_mask = .gripui.family.desc(
      id = "recursive_triangle_mask",
      label = "Recursive triangle mask",
      category = "Recursive triangle masks",
      function_name = "recursive.triangle.mask.surface.graph",
      summary = "Generic recursive triangle-mask family covering classic and bridge gaskets.",
      implementation = graph_impl,
      params = c(
        list(
          .gripui.family.param.choice("mask_kind", "Mask pattern", "classic", c("classic", "bridge")),
          .gripui.family.param.choice("missing", "Bridge omission", "top", c("top", "left", "right"), visible_if = list(mask_kind = "bridge")),
          .gripui.family.param.int("level", "Recursion level", 2L, 1L, 5L)
        ),
        .gripui.family.surface2d.params(c("flat", "saddle", "paraboloid", "ripple", "folded"))
      ),
      presets = list(
        bridge_top = list(mask_kind = "bridge", missing = "top"),
        bridge_left = list(mask_kind = "bridge", missing = "left"),
        folded = list(surface = "folded", amplitude = 0.65)
      ),
      builder = function(p) {
        args <- list(
          mask = .gripui.family.triangle.mask(p),
          level = p$level,
          surface = p$surface,
          amplitude = p$amplitude,
          freq_u = p$freq_u,
          freq_v = p$freq_v,
          x_scale = p$x_scale,
          y_scale = p$y_scale,
          normalize = p$normalize
        )
        do.call(recursive.triangle.mask.surface.graph, args)
      },
      code = function(p) {
        .gripui.family.call.code.parts("recursive.triangle.mask.surface.graph", c(
          mask = .gripui.family.triangle.mask.code(p),
          level = paste(deparse(p$level), collapse = " "),
          surface = paste(deparse(p$surface), collapse = " "),
          amplitude = paste(deparse(p$amplitude), collapse = " "),
          freq_u = paste(deparse(p$freq_u), collapse = " "),
          freq_v = paste(deparse(p$freq_v), collapse = " "),
          x_scale = paste(deparse(p$x_scale), collapse = " "),
          y_scale = paste(deparse(p$y_scale), collapse = " "),
          normalize = paste(deparse(p$normalize), collapse = " ")
        ))
      }
    ),
    sierpinski_triangle = .gripui.family.simple.desc(
      id = "sierpinski_triangle",
      label = "Sierpinski triangle",
      category = "Recursive triangle masks",
      function_name = "sierpinski.triangle.surface.graph",
      summary = "Classic Sierpinski triangle with flat and lifted surface variants.",
      implementation = graph_impl,
      arg_ids = c("level", "surface", "amplitude", "freq_u", "freq_v", "x_scale", "y_scale", "normalize"),
      params = c(
        list(.gripui.family.param.int("level", "Recursion level", 2L, 1L, 5L)),
        .gripui.family.surface2d.params(c("flat", "saddle", "paraboloid", "ripple", "folded"))
      ),
      presets = list(
        folded = list(surface = "folded"),
        ripple = list(surface = "ripple", amplitude = 0.55, freq_u = 2L, freq_v = 2L)
      )
    ),
    recursive_tetrahedron_mask = .gripui.family.desc(
      id = "recursive_tetrahedron_mask",
      label = "Recursive tetrahedron mask",
      category = "Recursive tetrahedron masks",
      function_name = "recursive.tetrahedron.mask.surface.graph",
      summary = "Corner-mask tetrahedral gasket family, including asymmetric omissions.",
      implementation = graph_impl,
      params = c(
        list(
          .gripui.family.param.choice("mask_kind", "Mask pattern", "classic", c("classic", "corner_missing")),
          .gripui.family.param.choice(
            "omit",
            "Omit corner",
            "apex",
            c("apex", "base_left", "base_right", "base_back"),
            visible_if = list(mask_kind = "corner_missing")
          ),
          .gripui.family.param.int("level", "Recursion level", 2L, 1L, 4L)
        ),
        .gripui.family.surface3d.params(c("standard", "squashed", "twisted", "wavy"), amplitude_default = 0.3)
      ),
      presets = list(
        apex_missing = list(mask_kind = "corner_missing", omit = "apex"),
        wavy = list(surface = "wavy")
      ),
      builder = function(p) {
        args <- list(
          mask = .gripui.family.tetrahedron.mask(p),
          level = p$level,
          surface = p$surface,
          amplitude = p$amplitude,
          freq = p$freq,
          twist = p$twist,
          normalize = p$normalize
        )
        do.call(recursive.tetrahedron.mask.surface.graph, args)
      },
      code = function(p) {
        .gripui.family.call.code.parts("recursive.tetrahedron.mask.surface.graph", c(
          mask = .gripui.family.tetrahedron.mask.code(p),
          level = paste(deparse(p$level), collapse = " "),
          surface = paste(deparse(p$surface), collapse = " "),
          amplitude = paste(deparse(p$amplitude), collapse = " "),
          freq = paste(deparse(p$freq), collapse = " "),
          twist = paste(deparse(p$twist), collapse = " "),
          normalize = paste(deparse(p$normalize), collapse = " ")
        ))
      }
    ),
    sierpinski_tetrahedron = .gripui.family.simple.desc(
      id = "sierpinski_tetrahedron",
      label = "Sierpinski tetrahedron",
      category = "Recursive tetrahedron masks",
      function_name = "sierpinski.tetrahedron.surface.graph",
      summary = "Classic tetrahedral gasket with squashed, twisted, and wavy variants.",
      implementation = graph_impl,
      arg_ids = c("level", "surface", "amplitude", "freq", "twist", "normalize"),
      params = c(
        list(.gripui.family.param.int("level", "Recursion level", 2L, 1L, 4L)),
        .gripui.family.surface3d.params(c("standard", "squashed", "twisted", "wavy"), amplitude_default = 0.3)
      ),
      presets = list(
        twisted = list(surface = "twisted"),
        wavy = list(surface = "wavy")
      )
    ),
    recursive_cube_mask = .gripui.family.desc(
      id = "recursive_cube_mask",
      label = "Recursive cube mask",
      category = "Recursive cube masks",
      function_name = "recursive.cube.mask.surface.graph",
      summary = "Generic recursive cube-mask family spanning Menger and porous-cube style masks.",
      implementation = graph_impl,
      params = c(
        list(
          .gripui.family.param.choice(
            "mask_kind",
            "Mask pattern",
            "menger",
            c("menger", "periodic_tunnels", "asymmetric_cavities", "channel_network")
          ),
          .gripui.family.param.int("level", "Recursion level", 1L, 1L, 3L),
          .gripui.family.param.int("side", "Mask side", 5L, 3L, 9L, visible_if = list(mask_kind = c("periodic_tunnels", "asymmetric_cavities", "channel_network"))),
          .gripui.family.param.int("tunnel_width", "Tunnel width", 1L, 1L, 3L, visible_if = list(mask_kind = "periodic_tunnels")),
          .gripui.family.param.int("tunnel_period", "Tunnel period", 2L, 2L, 6L, visible_if = list(mask_kind = "periodic_tunnels")),
          .gripui.family.param.int("tunnel_offset", "Tunnel offset", 2L, 1L, 6L, visible_if = list(mask_kind = "periodic_tunnels")),
          .gripui.family.param.int("cavity_size", "Cavity size", 2L, 1L, 4L, visible_if = list(mask_kind = "asymmetric_cavities")),
          .gripui.family.param.int("pocket_size", "Pocket size", 1L, 1L, 4L, visible_if = list(mask_kind = "asymmetric_cavities")),
          .gripui.family.param.int("channel_width", "Channel width", 1L, 1L, 3L, visible_if = list(mask_kind = "channel_network")),
          .gripui.family.param.int("branch_offset", "Branch offset", 2L, 1L, 6L, visible_if = list(mask_kind = "channel_network"))
        ),
        .gripui.family.surface3d.params(c("standard", "bulged", "twisted", "wavy"), include_xyz_scales = TRUE)
      ),
      presets = list(
        tunnels = list(mask_kind = "periodic_tunnels"),
        cavities = list(mask_kind = "asymmetric_cavities"),
        channels = list(mask_kind = "channel_network"),
        wavy = list(surface = "wavy")
      ),
      builder = function(p) {
        args <- list(
          mask = .gripui.family.cube.mask(p),
          level = p$level,
          surface = p$surface,
          amplitude = p$amplitude,
          freq = p$freq,
          twist = p$twist,
          x_scale = p$x_scale,
          y_scale = p$y_scale,
          z_scale = p$z_scale,
          normalize = p$normalize
        )
        do.call(recursive.cube.mask.surface.graph, args)
      },
      code = function(p) {
        .gripui.family.call.code.parts("recursive.cube.mask.surface.graph", c(
          mask = .gripui.family.cube.mask.code(p),
          level = paste(deparse(p$level), collapse = " "),
          surface = paste(deparse(p$surface), collapse = " "),
          amplitude = paste(deparse(p$amplitude), collapse = " "),
          freq = paste(deparse(p$freq), collapse = " "),
          twist = paste(deparse(p$twist), collapse = " "),
          x_scale = paste(deparse(p$x_scale), collapse = " "),
          y_scale = paste(deparse(p$y_scale), collapse = " "),
          z_scale = paste(deparse(p$z_scale), collapse = " "),
          normalize = paste(deparse(p$normalize), collapse = " ")
        ))
      }
    ),
    menger_sponge = .gripui.family.simple.desc(
      id = "menger_sponge",
      label = "Menger sponge",
      category = "Recursive cube masks",
      function_name = "menger.sponge.surface.graph",
      summary = "Classic cubical recursive porous family.",
      implementation = graph_impl,
      arg_ids = c("level", "surface", "amplitude", "freq", "twist", "x_scale", "y_scale", "z_scale", "normalize"),
      params = c(
        list(.gripui.family.param.int("level", "Recursion level", 2L, 1L, 3L)),
        .gripui.family.surface3d.params(c("standard", "bulged", "twisted", "wavy"), include_xyz_scales = TRUE)
      ),
      presets = list(
        bulged = list(surface = "bulged"),
        wavy = list(surface = "wavy")
      )
    ),
    cube_periodic_tunnels = .gripui.family.simple.desc(
      id = "cube_periodic_tunnels",
      label = "Cube periodic tunnels",
      category = "Recursive cube masks",
      function_name = "cube.periodic.tunnels.surface.graph",
      summary = "Recursive cube family with periodic tunnel drilling.",
      implementation = graph_impl,
      arg_ids = c("level", "side", "tunnel_width", "tunnel_period", "tunnel_offset", "surface", "amplitude", "freq", "twist", "x_scale", "y_scale", "z_scale", "normalize"),
      params = c(
        list(
          .gripui.family.param.int("level", "Recursion level", 1L, 1L, 3L),
          .gripui.family.param.int("side", "Mask side", 5L, 3L, 9L),
          .gripui.family.param.int("tunnel_width", "Tunnel width", 1L, 1L, 3L),
          .gripui.family.param.int("tunnel_period", "Tunnel period", 2L, 2L, 6L),
          .gripui.family.param.int("tunnel_offset", "Tunnel offset", 2L, 1L, 6L)
        ),
        .gripui.family.surface3d.params(c("standard", "bulged", "twisted", "wavy"), include_xyz_scales = TRUE)
      ),
      presets = list(
        bulged = list(surface = "bulged"),
        dense = list(side = 7L, level = 1L, tunnel_period = 3L)
      )
    ),
    cube_asymmetric_cavities = .gripui.family.simple.desc(
      id = "cube_asymmetric_cavities",
      label = "Cube asymmetric cavities",
      category = "Recursive cube masks",
      function_name = "cube.asymmetric.cavities.surface.graph",
      summary = "Recursive cubical family with asymmetric interior cavities and pockets.",
      implementation = graph_impl,
      arg_ids = c("level", "side", "cavity_size", "pocket_size", "surface", "amplitude", "freq", "twist", "x_scale", "y_scale", "z_scale", "normalize"),
      params = c(
        list(
          .gripui.family.param.int("level", "Recursion level", 1L, 1L, 3L),
          .gripui.family.param.int("side", "Mask side", 5L, 3L, 9L),
          .gripui.family.param.int("cavity_size", "Cavity size", 2L, 1L, 4L),
          .gripui.family.param.int("pocket_size", "Pocket size", 1L, 1L, 4L)
        ),
        .gripui.family.surface3d.params(c("standard", "bulged", "twisted", "wavy"), include_xyz_scales = TRUE)
      ),
      presets = list(
        twisted = list(surface = "twisted"),
        deeper = list(side = 7L, level = 1L, cavity_size = 3L, pocket_size = 2L)
      )
    ),
    cube_channel_network = .gripui.family.simple.desc(
      id = "cube_channel_network",
      label = "Cube channel network",
      category = "Recursive cube masks",
      function_name = "cube.channel.network.surface.graph",
      summary = "Recursive cubical family with branching channel networks.",
      implementation = graph_impl,
      arg_ids = c("level", "side", "channel_width", "branch_offset", "surface", "amplitude", "freq", "twist", "x_scale", "y_scale", "z_scale", "normalize"),
      params = c(
        list(
          .gripui.family.param.int("level", "Recursion level", 1L, 1L, 3L),
          .gripui.family.param.int("side", "Mask side", 5L, 3L, 9L),
          .gripui.family.param.int("channel_width", "Channel width", 1L, 1L, 3L),
          .gripui.family.param.int("branch_offset", "Branch offset", 2L, 1L, 6L)
        ),
        .gripui.family.surface3d.params(c("standard", "bulged", "twisted", "wavy"), include_xyz_scales = TRUE)
      ),
      presets = list(
        wavy = list(surface = "wavy"),
        wide = list(side = 7L, level = 1L, channel_width = 2L)
      )
    ),
    triangulated_polyhedron = .gripui.family.simple.desc(
      id = "triangulated_polyhedron",
      label = "Triangulated polyhedron",
      category = "Triangulated manifolds",
      function_name = "triangulated.polyhedron.surface.graph",
      summary = "Closed triangulated manifold from subdivided tetrahedron, octahedron, or icosahedron bases.",
      implementation = graph_impl,
      arg_ids = c("base", "level", "surface", "amplitude", "freq", "twist", "normalize"),
      params = c(
        list(
          .gripui.family.param.choice("base", "Base polyhedron", "icosahedron", c("tetrahedron", "octahedron", "icosahedron")),
          .gripui.family.param.int("level", "Subdivision level", 1L, 0L, 3L)
        ),
        .gripui.family.surface3d.params(c("standard", "inflated", "twisted", "wavy"), amplitude_default = 0.25)
      ),
      presets = list(
        octahedron = list(base = "octahedron"),
        inflated = list(surface = "inflated")
      )
    ),
    triangulated_annulus = .gripui.family.simple.desc(
      id = "triangulated_annulus",
      label = "Triangulated annulus",
      category = "Triangulated manifolds",
      function_name = "triangulated.annulus.surface.graph",
      summary = "Boundary triangulated annulus built from a clipped triangular lattice.",
      implementation = graph_impl,
      arg_ids = c("resolution", "outer_radius", "inner_radius", "surface", "amplitude", "freq_u", "freq_v", "normalize"),
      params = c(
        list(
          .gripui.family.param.int("resolution", "Resolution", 12L, 4L, 36L),
          .gripui.family.param.double("outer_radius", "Outer radius", 1, 0.2, 4, 0.05),
          .gripui.family.param.double("inner_radius", "Inner radius", 0.45, 0.05, 2.5, 0.05)
        ),
        .gripui.family.surface2d.params(c("flat", "saddle", "paraboloid", "ripple", "folded"), amplitude_default = 0.6, include_scales = FALSE)
      ),
      presets = list(
        folded = list(surface = "folded"),
        ripple = list(surface = "ripple", amplitude = 0.5, freq_u = 2L, freq_v = 2L)
      )
    ),
    triangulated_pair_of_pants = .gripui.family.simple.desc(
      id = "triangulated_pair_of_pants",
      label = "Triangulated pair of pants",
      category = "Triangulated manifolds",
      function_name = "triangulated.pair.of.pants.surface.graph",
      summary = "Boundary triangulated pair-of-pants family built from a clipped triangular lattice.",
      implementation = graph_impl,
      arg_ids = c("resolution", "outer_radius", "hole_radius", "hole_offset", "hole_height", "surface", "amplitude", "freq_u", "freq_v", "normalize"),
      params = c(
        list(
          .gripui.family.param.int("resolution", "Resolution", 12L, 4L, 36L),
          .gripui.family.param.double("outer_radius", "Outer radius", 1.1, 0.2, 4, 0.05),
          .gripui.family.param.double("hole_radius", "Hole radius", 0.24, 0.05, 1.5, 0.02),
          .gripui.family.param.double("hole_offset", "Hole offset", 0.38, 0.05, 1.5, 0.02),
          .gripui.family.param.double("hole_height", "Hole height", 0.18, 0.01, 1, 0.02)
        ),
        .gripui.family.surface2d.params(c("flat", "saddle", "paraboloid", "ripple", "folded"), amplitude_default = 0.6, include_scales = FALSE)
      ),
      presets = list(
        folded = list(surface = "folded"),
        ripple = list(surface = "ripple", amplitude = 0.5, freq_u = 2L, freq_v = 2L)
      )
    ),
    irregular_annulus = .gripui.family.simple.desc(
      id = "irregular_annulus",
      label = "Irregular annulus",
      category = "Irregular manifolds",
      function_name = "irregular.annulus.surface.graph",
      summary = "Point-sampled annulus with irregular counts, spacing, and phases.",
      implementation = graph_impl,
      arg_ids = c("rings", "outer_count", "outer_radius", "inner_radius", "count_irregularity", "radial_irregularity", "phase_twist", "surface", "amplitude", "freq_u", "freq_v", "normalize"),
      params = c(
        list(
          .gripui.family.param.int("rings", "Rings", 6L, 3L, 18L),
          .gripui.family.param.int("outer_count", "Outer count", 28L, 8L, 80L),
          .gripui.family.param.double("outer_radius", "Outer radius", 1, 0.2, 4, 0.05),
          .gripui.family.param.double("inner_radius", "Inner radius", 0.45, 0.05, 2.5, 0.05),
          .gripui.family.param.double("count_irregularity", "Count irregularity", 0.2, 0, 0.8, 0.05),
          .gripui.family.param.double("radial_irregularity", "Radial irregularity", 0.35, 0, 1, 0.05),
          .gripui.family.param.double("phase_twist", "Phase twist", 0.35, 0, 1.5, 0.05)
        ),
        .gripui.family.surface2d.params(c("flat", "saddle", "paraboloid", "ripple", "folded"), amplitude_default = 0.6, include_scales = FALSE)
      ),
      presets = list(
        folded = list(surface = "folded"),
        more_irregular = list(count_irregularity = 0.35, radial_irregularity = 0.5)
      )
    ),
    irregular_sphere = .gripui.family.simple.desc(
      id = "irregular_sphere",
      label = "Irregular sphere",
      category = "Irregular manifolds",
      function_name = "irregular.sphere.surface.graph",
      summary = "Point-sampled sphere with irregular band densities and phases.",
      implementation = graph_impl,
      arg_ids = c("bands", "equator_count", "count_irregularity", "lat_irregularity", "phase_twist", "surface", "radius", "amplitude", "freq_theta", "freq_lat", "twist", "normalize"),
      params = list(
        .gripui.family.param.int("bands", "Bands", 6L, 3L, 18L),
        .gripui.family.param.int("equator_count", "Equator count", 28L, 8L, 80L),
        .gripui.family.param.double("count_irregularity", "Count irregularity", 0.2, 0, 0.8, 0.05),
        .gripui.family.param.double("lat_irregularity", "Latitudinal irregularity", 0.35, 0, 1, 0.05),
        .gripui.family.param.double("phase_twist", "Phase twist", 0.35, 0, 1.5, 0.05),
        .gripui.family.param.choice("surface", "Surface", "standard", c("standard", "ellipsoid", "wavy"), group = "Geometry"),
        .gripui.family.param.double("radius", "Radius", 1, 0.2, 4, 0.05, group = "Geometry"),
        .gripui.family.param.double("amplitude", "Amplitude", 0.2, 0, 1, 0.02, group = "Geometry"),
        .gripui.family.param.int("freq_theta", "Longitudinal frequency", 3L, 1L, 8L, group = "Geometry"),
        .gripui.family.param.int("freq_lat", "Latitudinal frequency", 2L, 1L, 8L, group = "Geometry"),
        .gripui.family.param.double("twist", "Twist", 0.25, 0, 2, 0.05, group = "Geometry"),
        .gripui.family.normalize.param()
      ),
      presets = list(
        ellipsoid = list(surface = "ellipsoid"),
        more_irregular = list(count_irregularity = 0.35, lat_irregularity = 0.5)
      )
    ),
    irregular_pair_of_pants = .gripui.family.simple.desc(
      id = "irregular_pair_of_pants",
      label = "Irregular pair of pants",
      category = "Irregular manifolds",
      function_name = "irregular.pair.of.pants.surface.graph",
      summary = "Point-sampled pair-of-pants family with slice irregularity and phase twist.",
      implementation = graph_impl,
      arg_ids = c("slices", "outer_count", "outer_radius", "hole_radius", "hole_offset", "hole_height", "count_irregularity", "vertical_irregularity", "phase_twist", "surface", "amplitude", "freq_u", "freq_v", "normalize"),
      params = c(
        list(
          .gripui.family.param.int("slices", "Slices", 11L, 5L, 24L),
          .gripui.family.param.int("outer_count", "Outer count", 28L, 8L, 80L),
          .gripui.family.param.double("outer_radius", "Outer radius", 1.1, 0.2, 4, 0.05),
          .gripui.family.param.double("hole_radius", "Hole radius", 0.24, 0.05, 1.5, 0.02),
          .gripui.family.param.double("hole_offset", "Hole offset", 0.38, 0.05, 1.5, 0.02),
          .gripui.family.param.double("hole_height", "Hole height", 0.18, 0.01, 1, 0.02),
          .gripui.family.param.double("count_irregularity", "Count irregularity", 0.2, 0, 0.8, 0.05),
          .gripui.family.param.double("vertical_irregularity", "Vertical irregularity", 0.35, 0, 1, 0.05),
          .gripui.family.param.double("phase_twist", "Phase twist", 0.35, 0, 1.5, 0.05)
        ),
        .gripui.family.surface2d.params(c("flat", "saddle", "paraboloid", "ripple", "folded"), amplitude_default = 0.6, include_scales = FALSE)
      ),
      presets = list(
        folded = list(surface = "folded"),
        more_irregular = list(count_irregularity = 0.35, vertical_irregularity = 0.5)
      )
    ),
    irregular_torus = .gripui.family.simple.desc(
      id = "irregular_torus",
      label = "Irregular torus",
      category = "Irregular manifolds",
      function_name = "irregular.torus.surface.graph",
      summary = "Point-sampled torus with irregular major-ring counts and phases.",
      implementation = graph_impl,
      arg_ids = c("major_rings", "tube_count", "count_irregularity", "major_irregularity", "phase_twist", "surface", "major_radius", "minor_radius", "amplitude", "freq_major", "freq_minor", "twist", "normalize"),
      params = list(
        .gripui.family.param.int("major_rings", "Major rings", 8L, 4L, 24L),
        .gripui.family.param.int("tube_count", "Tube count", 16L, 6L, 64L),
        .gripui.family.param.double("count_irregularity", "Count irregularity", 0.2, 0, 0.8, 0.05),
        .gripui.family.param.double("major_irregularity", "Major irregularity", 0.25, 0, 1, 0.05),
        .gripui.family.param.double("phase_twist", "Phase twist", 0.35, 0, 1.5, 0.05),
        .gripui.family.param.choice("surface", "Surface", "standard", c("standard", "pinched", "wavy"), group = "Geometry"),
        .gripui.family.param.double("major_radius", "Major radius", 2, 0.5, 6, 0.05, group = "Geometry"),
        .gripui.family.param.double("minor_radius", "Minor radius", 0.75, 0.1, 3, 0.05, group = "Geometry"),
        .gripui.family.param.double("amplitude", "Amplitude", 0.2, 0, 1, 0.02, group = "Geometry"),
        .gripui.family.param.int("freq_major", "Major frequency", 2L, 1L, 8L, group = "Geometry"),
        .gripui.family.param.int("freq_minor", "Minor frequency", 1L, 1L, 8L, group = "Geometry"),
        .gripui.family.param.double("twist", "Twist", 0.25, 0, 2, 0.05, group = "Geometry"),
        .gripui.family.normalize.param()
      ),
      presets = list(
        pinched = list(surface = "pinched"),
        more_irregular = list(count_irregularity = 0.35, major_irregularity = 0.45)
      )
    ),
    irregular_double_torus = .gripui.family.simple.desc(
      id = "irregular_double_torus",
      label = "Irregular double torus",
      category = "Irregular manifolds",
      function_name = "irregular.double.torus.surface.graph",
      summary = "Closed genus-2 point-sampled surface with irregular slice and tube structure.",
      implementation = graph_impl,
      arg_ids = c("slices", "tube_count", "branch_length", "branch_offset", "tube_radius", "transition_width", "count_irregularity", "axial_irregularity", "phase_twist", "surface", "amplitude", "freq_x", "freq_theta", "twist", "normalize"),
      params = list(
        .gripui.family.param.int("slices", "Slices", 11L, 5L, 24L),
        .gripui.family.param.int("tube_count", "Tube count", 14L, 6L, 64L),
        .gripui.family.param.double("branch_length", "Branch length", 0.85, 0.1, 2.5, 0.05),
        .gripui.family.param.double("branch_offset", "Branch offset", 0.72, 0.1, 2, 0.05),
        .gripui.family.param.double("tube_radius", "Tube radius", 0.28, 0.05, 1.5, 0.02),
        .gripui.family.param.double("transition_width", "Transition width", 0.42, 0.05, 1.5, 0.02),
        .gripui.family.param.double("count_irregularity", "Count irregularity", 0.2, 0, 0.8, 0.05),
        .gripui.family.param.double("axial_irregularity", "Axial irregularity", 0.3, 0, 1, 0.05),
        .gripui.family.param.double("phase_twist", "Phase twist", 0.35, 0, 1.5, 0.05),
        .gripui.family.param.choice("surface", "Surface", "standard", c("standard", "bulged", "twisted", "wavy"), group = "Geometry"),
        .gripui.family.param.double("amplitude", "Amplitude", 0.25, 0, 1.5, 0.05, group = "Geometry"),
        .gripui.family.param.int("freq_x", "Axial frequency", 2L, 1L, 8L, group = "Geometry"),
        .gripui.family.param.int("freq_theta", "Angular frequency", 2L, 1L, 8L, group = "Geometry"),
        .gripui.family.param.double("twist", "Twist", 0.6, 0, 2, 0.05, group = "Geometry"),
        .gripui.family.normalize.param()
      ),
      presets = list(
        bulged = list(surface = "bulged"),
        twisted = list(surface = "twisted")
      )
    ),
    irregular_ball = .gripui.family.simple.desc(
      id = "irregular_ball",
      label = "Irregular ball solid",
      category = "Volumetric solids",
      function_name = "irregular.ball.solid.graph",
      summary = "Tetrahedralized volumetric ball built from nested subdivided polyhedral shells.",
      implementation = graph_impl,
      arg_ids = c("base", "level", "layers", "outer_radius", "radial_irregularity", "layer_twist", "surface", "amplitude", "freq_theta", "freq_phi", "twist", "normalize"),
      params = list(
        .gripui.family.param.choice("base", "Base polyhedron", "icosahedron", c("tetrahedron", "octahedron", "icosahedron")),
        .gripui.family.param.int("level", "Subdivision level", 1L, 0L, 2L),
        .gripui.family.param.int("layers", "Radial layers", 3L, 2L, 6L),
        .gripui.family.param.double("outer_radius", "Outer radius", 1, 0.2, 4, 0.05),
        .gripui.family.param.double("radial_irregularity", "Radial irregularity", 0.25, 0, 1, 0.05),
        .gripui.family.param.double("layer_twist", "Layer twist", 0.35, 0, 1.5, 0.05),
        .gripui.family.param.choice("surface", "Surface", "standard", c("standard", "bulged", "twisted", "wavy"), group = "Geometry"),
        .gripui.family.param.double("amplitude", "Amplitude", 0.2, 0, 1, 0.02, group = "Geometry"),
        .gripui.family.param.int("freq_theta", "Theta frequency", 2L, 1L, 8L, group = "Geometry"),
        .gripui.family.param.int("freq_phi", "Phi frequency", 2L, 1L, 8L, group = "Geometry"),
        .gripui.family.param.double("twist", "Twist", 0.6, 0, 2, 0.05, group = "Geometry"),
        .gripui.family.normalize.param()
      ),
      presets = list(
        octahedron = list(base = "octahedron"),
        wavy = list(surface = "wavy")
      )
    ),
    irregular_shell = .gripui.family.simple.desc(
      id = "irregular_shell",
      label = "Irregular shell solid",
      category = "Volumetric solids",
      function_name = "irregular.shell.solid.graph",
      summary = "Tetrahedralized hollow shell built from nested subdivided polyhedral shells.",
      implementation = graph_impl,
      arg_ids = c("base", "level", "layers", "inner_radius", "outer_radius", "radial_irregularity", "layer_twist", "surface", "amplitude", "freq_theta", "freq_phi", "twist", "normalize"),
      params = list(
        .gripui.family.param.choice("base", "Base polyhedron", "icosahedron", c("tetrahedron", "octahedron", "icosahedron")),
        .gripui.family.param.int("level", "Subdivision level", 1L, 0L, 2L),
        .gripui.family.param.int("layers", "Radial layers", 3L, 2L, 6L),
        .gripui.family.param.double("inner_radius", "Inner radius", 0.45, 0.05, 2, 0.05),
        .gripui.family.param.double("outer_radius", "Outer radius", 1, 0.2, 4, 0.05),
        .gripui.family.param.double("radial_irregularity", "Radial irregularity", 0.25, 0, 1, 0.05),
        .gripui.family.param.double("layer_twist", "Layer twist", 0.35, 0, 1.5, 0.05),
        .gripui.family.param.choice("surface", "Surface", "standard", c("standard", "bulged", "twisted", "wavy"), group = "Geometry"),
        .gripui.family.param.double("amplitude", "Amplitude", 0.2, 0, 1, 0.02, group = "Geometry"),
        .gripui.family.param.int("freq_theta", "Theta frequency", 2L, 1L, 8L, group = "Geometry"),
        .gripui.family.param.int("freq_phi", "Phi frequency", 2L, 1L, 8L, group = "Geometry"),
        .gripui.family.param.double("twist", "Twist", 0.6, 0, 2, 0.05, group = "Geometry"),
        .gripui.family.normalize.param()
      ),
      presets = list(
        octahedron = list(base = "octahedron"),
        twisted = list(surface = "twisted")
      )
    ),
    kary_tree = .gripui.family.desc(
      id = "kary_tree",
      label = "Intrinsic weighted k-ary tree",
      category = "Intrinsic trees",
      function_name = "kary.tree.weighted.graph",
      summary = "Tree family with intrinsic edge lengths controlled directly by depth and branch rules.",
      implementation = graph_impl,
      params = list(
        .gripui.family.param.int("k", "Branching factor", 3L, 2L, 6L),
        .gripui.family.param.int("depth", "Depth", 3L, 1L, 7L),
        .gripui.family.param.double("base_length", "Base length", 1, 0.1, 4, 0.05),
        .gripui.family.param.choice("depth_rule", "Depth rule", "geometric", c("geometric", "constant", "custom")),
        .gripui.family.param.double("depth_decay", "Depth decay", 0.85, 0.1, 1.5, 0.05, visible_if = list(depth_rule = "geometric")),
        .gripui.family.param.numeric_vector(
          "depth_factors",
          "Depth factors",
          default = numeric(),
          visible_if = list(depth_rule = "custom"),
          help = "Comma-separated values, one per edge depth."
        ),
        .gripui.family.param.choice("branch_rule", "Branch rule", "linear", c("linear", "uniform", "custom")),
        .gripui.family.param.double("branch_spread", "Branch spread", 0.3, 0, 2, 0.05, visible_if = list(branch_rule = "linear")),
        .gripui.family.param.numeric_vector(
          "branch_factors",
          "Branch factors",
          default = numeric(),
          visible_if = list(branch_rule = "custom"),
          help = "Comma-separated values, one per child slot."
        ),
        .gripui.family.normalize.param()
      ),
      presets = list(
        binary = list(k = 2L),
        uniform = list(branch_rule = "uniform", depth_rule = "constant"),
        custom = list(depth_rule = "custom", depth_factors = c(1, 0.8, 0.6), branch_rule = "custom", branch_factors = c(0.85, 1, 1.15))
      ),
      builder = function(p) {
        p$depth_factors <- if (length(p$depth_factors) == 0L) NULL else p$depth_factors
        p$branch_factors <- if (length(p$branch_factors) == 0L) NULL else p$branch_factors
        args <- list(
          k = p$k,
          depth = p$depth,
          base_length = p$base_length,
          depth_rule = p$depth_rule,
          depth_decay = p$depth_decay,
          depth_factors = p$depth_factors,
          branch_rule = p$branch_rule,
          branch_spread = p$branch_spread,
          branch_factors = p$branch_factors,
          normalize = p$normalize
        )
        do.call(kary.tree.weighted.graph, args)
      },
      code = function(p) {
        p$depth_factors <- if (length(p$depth_factors) == 0L) NULL else p$depth_factors
        p$branch_factors <- if (length(p$branch_factors) == 0L) NULL else p$branch_factors
        .gripui.family.call.code("kary.tree.weighted.graph", list(
          k = p$k,
          depth = p$depth,
          base_length = p$base_length,
          depth_rule = p$depth_rule,
          depth_decay = p$depth_decay,
          depth_factors = p$depth_factors,
          branch_rule = p$branch_rule,
          branch_spread = p$branch_spread,
          branch_factors = p$branch_factors,
          normalize = p$normalize
        ))
      }
    )
  )
}
