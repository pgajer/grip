# EXP-6 Concrete Benchmark Spec: MISF-GKK Multiscale Benchmark

Date: 2026-04-05

This note turns `EXP-6` from
[experiments_plan_2026-04-05.md](/Users/pgajer/current_projects/grip/dev/papers/gkk_lgkk_paper/experiments_plan_2026-04-05.md)
into a concrete runnable benchmark specification.

Current driver:

- `tools/benchmarks/gkk_lgkk_paper/benchmark-misf-gkk-panel.R`

Current smoke artifact:

- `output/gkk_lgkk_paper/tmp/misf-gkk-panel-smoke-2026-04-05/`

## Purpose

This benchmark should justify the main-text multiscale contribution of Paper 2:

- `MISF-GKK` is a real standalone multiscale optimizer in the same method family as `GKK` and `LGKK`
- `MISF-GKK` retains most of the geodesic-quality gain of direct geodesic polishing in the overlap regime
- `MISF-GKK` becomes the practical option once graph size grows enough that direct full `GKK` is expensive

At the moment, that is still the target claim rather than the observed result:

- the smoke benchmark runs successfully
- but the current default `MISF-GKK` configuration is not yet competitive with direct `KK -> LGKK` or `KK -> GKK` on the reduced overlap panel
- this benchmark spec therefore doubles as a tuning target before the full paper-grade run

This benchmark is **not** the place to make weighted `GRIP` the main scaling story. That is a separate comparison track.

## Primary Evaluation Setting

- embedding dimension: `3D` only
- shortest-path tie handling: `tie_mode = "average"`
- seeds: `1, 2, 3`
- primary target metric: exact full-`GKK` relative RMSE scored on the final coordinates
- primary cost metric: wall-clock runtime in seconds

Why `3D` only:

- the multiscale question is about scalable geodesic fidelity on intrinsically geometric weighted graphs
- the earlier mesh and cross-family sections already carry the broader 2D versus 3D story
- keeping `EXP-6` 3D-only prevents the multiscale section from becoming too diffuse

## Family Panel

Use four families only. This is deliberate: `EXP-6` should stay disciplined.

### Family F1: Mesh saddle

Constructor:

```r
mesh.surface.graph(h, w, surface = "saddle", amplitude = 0.75)
```

Exact size ladder:

- `12 x 12` with `n = 144`, `m = 264`
- `16 x 16` with `n = 256`, `m = 480`
- `20 x 20` with `n = 400`, `m = 760`

### Family F2: Torus pinched

Constructor:

```r
torus.surface.graph(h, w, surface = "pinched", amplitude = 0.22)
```

Exact size ladder:

- `12 x 12` with `n = 144`, `m = 288`
- `16 x 16` with `n = 256`, `m = 512`
- `20 x 20` with `n = 400`, `m = 800`

### Family F3: Irregular annulus folded

Constructor:

```r
irregular.annulus.surface.graph(
  rings = rings,
  outer_count = outer_count,
  surface = "folded",
  amplitude = 0.45
)
```

Exact size ladder:

- `rings = 8, outer_count = 32` with `n = 192`, `m = 544`
- `rings = 10, outer_count = 40` with `n = 297`, `m = 851`
- `rings = 12, outer_count = 48` with `n = 430`, `m = 1240`

### Family F4: Irregular torus pinched

Constructor:

```r
irregular.torus.surface.graph(
  major_rings = major_rings,
  tube_count = tube_count,
  surface = "pinched",
  amplitude = 0.18
)
```

Exact size ladder:

- `major_rings = 9, tube_count = 18` with `n = 162`, `m = 502`
- `major_rings = 12, tube_count = 24` with `n = 288`, `m = 884`
- `major_rings = 14, tube_count = 28` with `n = 392`, `m = 1200`

## Regime Split

The benchmark should be analyzed as two regimes.

### Regime A: overlap regime

Goal:

- compare `MISF-GKK` directly against exact full `GKK` and direct `LGKK` where full-geodesic optimization is still feasible

Instances:

- all families at the first two size tiers

That is:

- mesh `12 x 12`, `16 x 16`
- torus `12 x 12`, `16 x 16`
- annulus `8 x 32`, `10 x 40`
- irregular torus `9 x 18`, `12 x 24`

### Regime B: scale regime

Goal:

- show where `MISF-GKK` becomes the practical multiscale option

Instances:

- all families at the largest size tier

That is:

- mesh `20 x 20`
- torus `20 x 20`
- annulus `12 x 48`
- irregular torus `14 x 28`

## Seed Policy

Use exactly:

- `seed = 1`
- `seed = 2`
- `seed = 3`

Seed usage should be deterministic and consistent:

- for direct methods: shared-start coordinate jitter
- for `MISF-GKK`: MISF extraction, top-level restarts, and any restart-driven initialization

## Initialization Policy

### Direct-method shared start

Use the same full-graph start construction for all direct baselines:

1. take the family’s parameter coordinates
2. normalize them as in the current `KK/GKK/LGKK` mesh benchmark
3. add seed-specific Gaussian jitter
4. center the layout

This should match the current style already used in
[benchmark-kk-gkk-lgkk-mesh-suite.R](/Users/pgajer/current_projects/grip/tools/benchmarks/gkk_lgkk_paper/benchmark-kk-gkk-lgkk-mesh-suite.R).

### MISF-GKK start

Use:

- `top_level_init = "cmdscale"`
- `top_level_restarts = 8`

This keeps the multiscale method deterministic enough for comparison while still allowing a reasonable top-level search.

## Comparison Methods

The benchmark should use two method panels, one per regime.

### Regime A method panel

Run exactly these methods:

1. `KK`
2. `KK -> GKK`
3. `KK -> LGKK`
4. `MISF-GKK-auto`
5. `MISF-GKK-full`

### Regime B method panel

Run exactly these methods:

1. `KK`
2. `KK -> LGKK`
3. `MISF-GKK-auto`

Do **not** require `KK -> GKK` on the largest tier in the core benchmark.
If you later find it affordable and helpful, it can be added as an appendix or supplementary result, but it should not block `EXP-6`.

## Exact Method Definitions

### M1: `KK`

Definition:

- classical `igraph::layout_with_kk()` on the full graph
- weighted by graph edge lengths
- initialized from the shared full-graph start

Label:

- `KK`

### M2: `KK -> GKK`

Definition:

- run `KK` first
- polish with exact full `GKK`

Recommended settings:

- `max_iter = 12`
- `scale_mode = "profiled"`
- `return_trace = TRUE`

Label:

- `KK->GKK`

### M3: `KK -> LGKK`

Definition:

- run `KK` first
- polish with sparse `LGKK`

Recommended settings:

- `local_nbrs = 10`
- `landmark_count = 10`
- `max_iter = 12`
- `return_trace = TRUE`

Label:

- `KK->LGKK`

### M4: `MISF-GKK-auto`

Definition:

- full public multiscale method using the new family
- exact `GKK` on sufficiently small active sets
- sparse `LGKK` when the active set grows beyond the configured limit

Use exactly these settings:

```r
grip.optimize.misf.geodesic.kk(
  edges = spec$edges,
  n = spec$n,
  edge_weights = spec$edge_weights,
  tie_mode = "average",
  num_init = 24L,
  num_nbrs = 20L,
  dim = 3L,
  top_level_pair_mode = "auto",
  top_level_full_limit = 256L,
  top_level_local_nbrs = 12L,
  top_level_landmark_count = 8L,
  top_level_restarts = 8L,
  top_level_max_iter = 16L,
  top_level_init = "cmdscale",
  insertion_mode = "weighted_kk",
  insertion_max_iter = 32L,
  refinement_pair_mode = "auto",
  refinement_full_limit = 192L,
  refinement_local_nbrs = 8L,
  refinement_landmark_count = 4L,
  refinement_anchor_weight = 0.05,
  refinement_anchor_weight_end = 0.05,
  refinement_continuation = "constant",
  refinement_max_iter = 8L,
  final_pair_mode = "auto",
  final_full_limit = 256L,
  final_local_nbrs = 12L,
  final_landmark_count = 8L,
  final_max_iter = 8L,
  return_trace = TRUE,
  seed = seed
)
```

Label:

- `MISF-GKK-auto`

### M5: `MISF-GKK-full`

Definition:

- exact-heavy multiscale reference
- used only in Regime A to quantify the approximation gap introduced by `auto`

Use exactly the same settings as `MISF-GKK-auto`, except:

```r
top_level_pair_mode = "full"
refinement_pair_mode = "full"
final_pair_mode = "full"
```

Label:

- `MISF-GKK-full`

## Metrics To Record

### Primary metrics

- exact full-`GKK` energy
- exact full-`GKK` relative RMSE
- runtime in seconds

### Secondary metrics

- classical `KK` relative RMSE
- sparse `LGKK` relative RMSE
- Procrustes RMSE against the target surface
- sampled stress

### MISF-GKK-specific diagnostics

For every `MISF-GKK` run, also record:

- top-level active-set size
- effective pair mode at the top level
- number of multiscale refinement stages
- per-stage pair mode
- timing breakdown:
  - top level
  - insertion
  - refinement
  - final polish
- final scoring mode

## Output Files

Recommended benchmark script:

- `tools/benchmarks/gkk_lgkk_paper/benchmark-misf-gkk-panel.R`

Recommended run tag:

- `misf-gkk-panel-<date>`

Recommended output root:

- `output/gkk_lgkk_paper/benchmarks/misf-gkk-panel-<date>/`

Required data outputs:

- `data/case_manifest.csv`
- `data/raw_results.csv`
- `data/method_summary.csv`
- `data/family_size_summary.csv`
- `data/misf_stage_summary.csv`
- `data/pair_mode_usage_summary.csv`

Required figure outputs:

- `figures/misf_gkk_quality_by_size_dim3.png`
- `figures/misf_gkk_runtime_by_size_dim3.png`
- `figures/misf_gkk_quality_runtime_pareto_dim3.png`
- `figures/misf_gkk_representative_gallery_dim3.png`
- `figures/misf_gkk_stage_breakdown_dim3.png`

Required table outputs:

- `tables/misf_gkk_method_summary_dim3.tex`
- `tables/misf_gkk_overlap_summary_dim3.tex`
- `tables/misf_gkk_scale_summary_dim3.tex`

Recommended rendered report:

- `output/gkk_lgkk_paper/reports/benchmarks/misf_gkk_panel_report_<date>/`

## Figures and Tables For The Paper

### Main-text figure

- `Figure 8`
  - runtime-quality or quality-by-size comparison for `MISF-GKK-auto`
  - should show Regime A and Regime B clearly

### Main-text table

- `Table 5`
  - compact multiscale benchmark summary
  - include:
    - method
    - regime
    - mean full-`GKK` relative RMSE
    - mean runtime

### Optional appendix material

- stage breakdown figure
- pair-mode usage table
- representative gallery panel

## Interpretation Rules

The benchmark should support these exact paper questions:

1. In Regime A, how close is `MISF-GKK-auto` to `KK->GKK` and `KK->LGKK` on the primary exact full-`GKK` metric?
2. In Regime A, how much approximation gap is introduced by `auto` relative to `MISF-GKK-full`?
3. In Regime B, does `MISF-GKK-auto` remain materially better than plain `KK` while keeping runtime practical?
4. Across both regimes, does the stage trace show sensible switching between full and landmark pair modes?

## Exclusions

These should **not** be in the core `EXP-6` benchmark:

- weighted `GRIP`
- weighted `GRIP + core LGKK`
- weighted `GRIP + polish LGKK`
- Sierpinski carpet
- cube-channel network
- 2D runs

They can appear elsewhere:

- weighted `GRIP` in `EXP-8`
- additional families in `EXP-2` or appendices

## Bottom Line

`EXP-6` should be a focused 3D benchmark over four families, three seeds, and two regimes:

- an overlap regime where direct `GKK` is still comparable
- a scale regime where `MISF-GKK` should become the practical multiscale method

That is the cleanest way to make `MISF-GKK` a real main-text contribution in Paper 2.
