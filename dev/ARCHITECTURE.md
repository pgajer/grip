# Developer map

## Follow one computation

| Stage | Maintained implementation | Contract |
|---|---|---|
| Public graph input | `R/grip_layout.R`, `R/grip_layout_weighted_nd.R`, `R/gmds_layout_interface.R`, `R/metric_mds.R` | Preserve established argument order and matrix/list returns. |
| Common graph validation | `R/grip_graph_inputs.R` | Check integer identity, representation exclusivity, reciprocal adjacency/lengths, and self-loops before conversion or dispatch. Ordinary and higher-dimensional layout validation share the same graph validator. |
| Graph construction | `R/graph_helpers.R` | Primitive and surface/recursive graph builders; the task guide and grouped reference index identify public families. File size alone is not a reason to split these coherent families. |
| Preparation | `R/grip_quality.R`, `R/geodesic_mds.R`, `R/preparation_resources.R` | Graph-specific caches; full and landmark preparation use dense distances. Estimate/warn before allocation. Never reuse caches after graph identity changes. |
| Ordinary native layout | `src/DrawGraph.cpp`, `MishEngine.cpp`, `MishWeighted.cpp`, `DrawGraphND.cpp` and wrappers | Vertex order and seed semantics remain stable. Interrupts occur on the main thread. Legacy array ownership uses shared cleanup; MISF queue temporaries use automatic ownership. |
| MDS/refinement engines | `R/metric_mds.R`, `R/geodesic_mds.R`, `R/geodesic_mds_bending*.R`, `src/geodesic_mds*.cpp` | Classical strain, metric raw stress, and retained-path energy are different objectives. Bending extensions explicitly call `.base` implementations. Each function has one definition; no load-time capture/rebinding or `zz_` load-order requirement. |
| Results and scoring | `R/gmds_layout_interface.R`, `R/grip_quality.R`, `R/reference_scores.R`, `R/layout_coords.R` | Numerical fields remain available independently of printing. Coordinate extraction does not change ordering, dimension, units, or alignment. |
| Display | `R/grip_plot.R`, `R/gripui-*.R`, `pkgdown/extra.*` | Static drawing is headless. Viewers remain optional. Website animations load only on request. |

The non-bending `.base` functions retain their existing bodies. The final bending
functions call those helpers by name; source filenames no longer decide which
implementation was captured. Keep native R API calls outside worker threads and
outside live-worker scopes. Cancellation returns no partial fitted object.

The installed **Finding your way around grip** guide contains the public input,
result, scale, initializer, cache-reuse, and experimental-interface contracts.
Do not bury user-relevant contracts in this developer map. Initializer spellings
remain method-specific as documented there; do not add aliases with merely similar
meanings.

## Maintained sources and generated output

- Roxygen comments in `R/` → `make document` → `man/` and `NAMESPACE`.
- `README.Rmd` and `tools/pkg/` media recipes → `make readme` → `README.md` and figures.
- `vignettes/*.Rmd` are installed guides; `vignettes/articles/*.Rmd` are website-only.
- `_pkgdown.yml` and `pkgdown/extra.*` control the generated website. Do not hand-edit `docs/`.
- Frozen benchmark results have `inst/extdata/vs_alternatives/BENCHMARK_PROVENANCE.md`.
  Saddle display provenance is in `tools/pkg/fixtures/saddle-preview.json`;
  unknown fitting settings must not be invented.
- `make audit-api-guide` checks every export, registered method, and help target.
- Public source and package archive checks remain in `make repo-hygiene`.
  Manuscript work and investigation logs stay in the existing private workspace.

## Verification

Use focused regression tests while editing. Changes to native or refinement
code need the existing finite-difference, bending, R/C++ parity, tied-path,
threading, and graph-identity tests, followed by a clean package check.

```sh
make document
make audit-api-guide
make build
Rscript tools/pkg/check-full.R grip_0.2.0.9001.tar.gz
```

Use the archive version produced by the build. `check-full.R` stages installed
optional packages in a temporary library and runs `R CMD check --as-cran` with
`GRIP_FULL_DEPENDENCIES=true`. It does not install packages or change the user's
libraries. This avoids a split-library failure where a system-library optional
package has a dependency installed only in a user library. Missing dependencies
still fail with the actual load error and library paths. Full-dependency runs
must have zero skipped tests; ordinary checks can intentionally lack optional
backends. Tests also exercise actionable failures for absent smacof and viewers.

CI runs the full-dependency contract on Ubuntu (R development, release, previous
release), Windows (release), and macOS (release), using headless rgl. A green
build does not certify all numerical cases; retain failures and supplementary
scientific evidence alongside the final candidate's results.
