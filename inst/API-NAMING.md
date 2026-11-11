# Migrating to dot-delimited API names

Public functions and parameters use lowercase words separated by dots. This is a
breaking change: replace underscores in old argument names with dots.
Positional calls retain their argument order.

For one transition release, `grip()` and `edge.kk()` also accept their former
underscore argument names, including names inside `edge.kk()`'s
`weighted.grip.args` list. These calls do not emit deprecation warnings.
Supplying both spellings of the same argument is an error, as are unknown
arguments. This lets existing downstream packages migrate before the old
spellings are removed in the following release. It does not restore the old
app entry points, other functions' argument names, or earlier MDS behavior.

For example:

```r
fit <- metric.mds(edges = graph$edges, edge.weights = graph$edge_weights,
                  dim = 3, max.iter = 100, pair.weights = "inverse_squared",
                  sgd.control = list(learning.rate = 0.5, checkpoint.every = 1))
fit <- metric.mds(edges = graph$edges, approximation = "sparse",
                  sparse.control = list(n.pivots = 20))
```

Apply the same change to named argument lists used with `do.call()`, to
`weighted.grip.args`, and to the layout settings in `compare.layouts()`'s
`candidates` and `search` lists. `params.from.summary()` returns the new names.

The renamed function entry points are:

| Previous name | Current name |
| --- | --- |
| `gripui_app` | `gripui.app` |
| `gripui_family_app` | `gripui.family.app` |
| `gripui_graph_family_catalog` | `gripui.graph.family.catalog` |
| `gripui_project` | `gripui.project` |
| `gripui_project_from_compare` | `gripui.project.from.compare` |
| `gripui_project_from_dir` | `gripui.project.from.dir` |
| `gripui_validate_project` | `gripui.validate.project` |
| `run_gripui` | `run.gripui` |
| `run_gripui_family` | `run.gripui.family` |

Mathematical notation such as `X` and `scale.L0` is retained. Argument names
required by external interfaces, including native bindings and S3 methods, retain
those interfaces' spelling. Internal helpers are not a supported public API.

This change concerns functions, parameters, and public control lists. It does not
rename enum values such as `"inverse_squared"` or `"edge_length"`, graph fields
such as `graph$edge_weights`, result columns, or S3 classes. Existing saved graph
and layout objects keep their field names. Newly returned optimizer control
settings use the new dot-delimited keys.
