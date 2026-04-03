# Weighted GRIP Phase 3

## Scope

Phase 3 adds three things on top of the Phase 1 and Phase 2 weighted GRIP API:

- weighted tuning presets exposed through
  `grip.layout.globalrep.weighted()`,
  `grip.layout.weighted()`, and
  `grip.layout.trace.weighted()`
- a representative benchmark suite that compares weighted GRIP against
  `KK`, `KK->GKK`, and `KK->LGKK`
- a cross-family weighted benchmark panel covering the family regimes called
  out in the phase plan

The existing combinatorial GRIP functions remain unchanged.

## Weighted Presets

The weighted wrappers now accept:

- `carpet`
- `mesh`
- `cylinder`
- `torus`
- `sphere`
- `irregular`
- `tree`

These presets only affect the shared high-level tuning fields:

- `placement`
- `rounds`
- `final_rounds`
- `num_init`
- `num_nbrs`
- `r`
- `s`
- `repulsion_factor`

The more specialized weighted-core knobs such as
`coarse_repulsion_factor`, `coarse_repulsion_sample`,
`coarse_repulsion_exact_below`, and LGKK post-polish settings still use the
explicit function arguments. Presets are applied field by field, so explicit
arguments override the preset cleanly.

### Preset Intent

- `mesh`: rectangular and near-mesh weighted surfaces
- `cylinder`: cylindrical wrapped grids
- `torus`: toroidal wrapped grids
- `sphere`: closed near-spherical weighted meshes
- `irregular`: irregular manifold-like weighted families
- `tree`: intrinsic weighted trees
- `carpet`: high-neighborhood recursive lattice families

## Benchmark Driver

The new Phase 3 benchmark driver is:

- `tools/benchmarks/gkk_lgkk_paper/benchmark-weighted-grip-family-panel.R`

It compares these methods:

- `GRIP`
- `Weighted GRIP`
- `Weighted GRIP + LGKK`
- `KK`
- `KK->GKK`
- `KK->LGKK`

### Family Panel

The current representative panel includes:

- mesh saddle
- cylinder hourglass
- torus pinched
- sphere wavy
- irregular annulus folded
- irregular sphere ellipsoid

Each family is run in both `d = 2` and `d = 3`.

### Metrics

The script records:

- runtime
- full-GKK energy and relative RMSE
- sparse-LGKK energy and relative RMSE
- classical KK relative RMSE
- sampled stress
- edge-length CV
- Procrustes RMSE to the family target geometry

### Running It

Full run from the repo root:

```sh
Rscript tools/benchmarks/gkk_lgkk_paper/benchmark-weighted-grip-family-panel.R
```

Smoke run:

```sh
Rscript tools/benchmarks/gkk_lgkk_paper/benchmark-weighted-grip-family-panel.R --smoke
```

The smoke run uses fewer seeds and fewer optimizer iterations and writes to a
temporary output directory unless `--out=...` is supplied.

## Validation

Phase 3 adds preset-equivalence and preset-rejection tests for:

- `grip.layout.weighted()`
- `grip.layout.globalrep.weighted()`
- `grip.layout.trace.weighted()`

The benchmark script is also designed to be smoke-testable with the `--smoke`
flag so it can be validated without producing a full report-sized run.
