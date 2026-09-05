# MDS followed by edge-KK on expanding paraboloids and saddles

This study compares classical scaling, metric stress MDS, and edge-KK across
radius, symmetric-kNN graph construction, sampling measure, and sample size.
The geodesic graph regime and the manuscript's ambient graph regime answer
different questions. The experiment also exposed contraction of the profiled,
unnormalized edge objective; fixed-target-scale controls were added before the
full sweep. Read [the pilot decision](PILOT.md) and [amended protocol](PROTOCOL.md).

The final deliverables are [the PDF report](report.pdf),
[the complete results](summary/RESULTS.md), and the full stratified tables under
`summary/`. Integration into the main software manuscript is a separate phase.

## Scope and controls

- Two graph surfaces: z=x²+y² and z=x²−y² over disks of radius r.
- Uniform base-disk and uniform surface-area samples, three independent n=240
  samples, radii 1,2,4,8,16,32,64; quantiles and angles coupled across radii.
- Geodesic ranking/weights and ambient ranking/weights kept separate.
- Symmetric union kNN, k=4,8,16,32,64,128,239; deterministic ties and labeled
  minimum-spanning-tree augmentation of disconnected graphs.
- Graph-classical and graph-stress MDS; fixed full-geodesic classical starts
  additionally isolate the refinement graph from initialization changes.
- Density-to-uniform edge-KK and equal-budget uniform-only controls, each at
  profiled and identity target scale. Three MDS starts, 1,000 steps per start.
- Nested n=480 check at r=64, k=8,16,64,128,479, both measures and surfaces.
- Perturbation, random/original starts, and additional budgets at r=64,
  k=8,32,128,239 on sample 1. These do not replace primary fits.

The complete design has **1,216 primary graphs** and **32 optimizer-sensitivity
graphs**. The pilot first used 72 graphs, then checked 16 of those graphs with
additional initializations/budgets before expansion. All achieved configurations
are local, finite-budget results. A small objective, graph-to-layout error,
or singular-value ratio is not a theorem of geometric recovery.

## Reproduction

Run commands from this repository root unless using `make -C` as below.
Use grip 0.2.0.9001 and smacof 2.1-7, plus igraph, Rcpp, Python NumPy/SciPy,
Numba, pandas, Matplotlib, and a LaTeX installation. The package and backend
versions are checked by the worker. Python and R versions are recorded in the
artifact manifest and saved session records. The manifest records HEAD at build
time and exact artifact hashes; the hashes identify the render inputs even when
the final result commit follows the build. No private notes are inputs.

To verify and render the compact, tracked results without refitting:

```sh
make -C tools/experiments/mds-edge-kk-radius verify
make -C tools/experiments/mds-edge-kk-radius report
```

To restore the exact fitting inputs and rerun (substantial computation):

```sh
python3 tools/experiments/mds-edge-kk-radius/restore-inputs.py
python3 tools/experiments/mds-edge-kk-radius/run-grid.py pilot --workers 4
python3 tools/experiments/mds-edge-kk-radius/run-grid.py extra --workers 4
python3 tools/experiments/mds-edge-kk-radius/run-grid.py main --workers 6
python3 tools/experiments/mds-edge-kk-radius/run-grid.py size --workers 2
python3 tools/experiments/mds-edge-kk-radius/run-grid.py extra --workers 4
make -C tools/experiments/mds-edge-kk-radius collect
```

The first `extra` invocation uses only completed pilot graphs; the final one
fills the remaining declared settings. The read-only `progress.py` reports checkpoint counts and worker errors.
Do not run overlapping workers on the
same case. Each fit checks protocol, input, package-source, and driver hashes
before reusing a checkpoint. Edits to fitting sources require a separate cache
or an intentional rerun, not bypassing identity checks.

To regenerate smooth geodesics rather than use the portable matrices:

```sh
python3 tools/experiments/mds-edge-kk-radius/prepare-inputs.py original
python3 tools/experiments/mds-edge-kk-radius/prepare-inputs.py replicate
python3 tools/experiments/mds-edge-kk-radius/prepare-inputs.py size
NUMBA_NUM_THREADS=2 python3 tools/experiments/mds-edge-kk-radius/validate-geodesics.py
```

Original audited caches are reused when available; otherwise the same numerical
solvers generate the inputs. Warm radii 0.1,0.25,0.5 are solver continuation
checkpoints, not extra fitted study cases. Portable matrices are under
`summary/inputs/`; full velocity/route/optimizer caches are under the ignored
`output/mds-edge-kk-radius/` directory. Regenerated numerical geodesics must be
validated; cross-platform bitwise equality is not promised.

## Data definitions

`scores.csv` stores each primary candidate, with exact surface, measure, sample,
radius, n, regime, k, and method. `_primary` means density continuation with
profiled scale; `_uniform` means uniform-only at profiled scale;
`_fixed_primary` and `_fixed_uniform` fix target scale to one. Unqualified
`classical` and `stress` are MDS before refinement. `full_classical` uses smooth
Δ directly and only appears in the geodesic regime. `original` is the generating
coordinate control. `optimizer-sensitivity.csv` holds the separate controls. Random and original
starts are shared by the classical and stress branch comparisons; duplicated
control rows are not independent replications.

- `path_rel` / `edge_rel`: target-profiled relative RMSE, each with its own
  diagnostic scale; targets are retained route lengths / graph edges.
- `raw_stress`: chord stress against strict graph distances at returned physical
  scale. This is the MDS criterion, not profiled edge-KK's objective.
- `raw_target_rmse`: raw chord residual norm divided by strict target norm.
- `stress1`: chord residual after target profiling, divided by chord norm.
- `graph_reference`: strict graph-distance error versus smooth Δ, physical units.
- `path_reference`: embedded retained-path error versus smooth Δ, calibrated
  using graph edges only. `chord_reference` uses endpoint chords at that scale.
- `procrustes`: relative coordinate error after similarity alignment with known
  correspondence; it is not a closest-surface or topology diagnostic.
- `sigma1`, `sigma2`, `sigma3`: returned physical singular values;
  the three ratio columns are invariant to uniform scale.
- `edge_scale`: unweighted graph-edge calibration b. Divide returned coordinates
  by b for the edge-calibrated physical geometry used in the figures.
- `max_edge_relative`: the maximum edge error at that common calibrated scale;
  unlike a mean residual, it can support the stated uniform edge-to-path bound.

`graphs.csv`, `locality.csv`, and `input-validation.csv` expose connectivity,
edge locality, and radial/angular sample coverage. `starts.csv` includes all
MDS starts and their termination; raw stresses are in physical input units.
Gaussian draws are matched across k, with each initial scale fitted to its own
target matrix by the adapter.
`optimizer-status.csv` records each edge-KK stage and gradient. `timings.csv`
records preparation and MDS times; edge-KK elapsed times are in the status table.
Timings were measured with concurrent workers and are not controlled benchmarks.
`summary/coordinates/` holds compact per-case coordinates, scores, and start
records. Uniform r²-normalized coordinates are obtained by dividing the saved physical
coordinates by the recorded radius squared; the spatial exports also retain
these values explicitly. Graph routes and full per-iteration traces remain in the restartable
local fitting cache and are regenerated by the fitting scripts.

Numerical geodesic validators report sampled discrepancies, not certified
all-pair accuracy bounds. Both shape ratios must be considered: s2/s1 diagnoses
elongation toward a line; s3/s2 diagnoses planarity only when s2 remains nonzero.
All spatial figures use a single uniform scale and equal physical axis units.
