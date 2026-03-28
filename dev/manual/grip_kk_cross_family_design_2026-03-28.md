# GRIP vs KK Cross-Family Benchmark Design

## Goal

Evaluate how well the current best generic GRIP setting transfers across the full implemented graph
suite, excluding sphere graphs, and compare it with `igraph::layout_with_kk()`. The benchmark should
answer three questions:

1. Does the current best generic setting remain strong outside the level-4 carpet?
2. Which GRIP setting is best for each family or graph?
3. Is there a single global GRIP setting that is best overall across the non-sphere suite?

## Working Assumption

For this experiment, "the current best setting" means the best generic universal setting discovered
 in the recent carpet round-2 study:

- sampled coarse repulsion
- GRIP multiscale base defaults
- `final_rounds = 32`

This is the `grip_f32` candidate below.

## Graph Suite

Use every implemented graph family except sphere:

- Sierpinski carpets, levels 2 through 6
- Sierpinski triangles, levels 2 through 6
- Sierpinski tetrahedrons, levels 2 through 6
- meshes: `8x8`, `10x10`, `12x12`
- cylinders: `8x8`, `10x10`, `12x12`
- tori: `8x8`, `12x12`, `16x16`
- paths: `32`, `64`, `128`
- cycles: `32`, `64`, `128`
- cubes: sides `4`, `6`, `8`
- trees: `k=2 depth=4`, `k=2 depth=5`, `k=2 depth=6`, `k=3 depth=4`

The only graph that will be pre-skipped by scale is the level-6 carpet if a global vertex cap is
needed for runtime safety.

## Benchmark Structure

### Stage 1: Broad GRIP candidate screen

Purpose:
- prune the universal candidate list before the expensive full benchmark

Representative graph rule:
- for each non-sphere family, take the largest graph with at most `10000` vertices

Universal GRIP candidates screened:
- `grip_default_adaptive`
- `grip_f32`
- `grip_f64`
- `grip_f96`
- `grip_f128`
- `grip_f192`
- `grip_f384`
- `grip_exact_f64`
- `grip_exact_f96`

Selection rule:
- rank candidates by mean quality score across the Stage-1 representatives
- carry the top four universal candidates into Stage 2
- always carry `grip_f32` and `grip_default_adaptive` even if they miss the top four

### Stage 2: Full non-sphere benchmark

Candidates:
- selected universal GRIP candidates from Stage 1
- `grip_f32`
- `grip_default_adaptive`
- family-specific preset references:
  - `grip_mesh_preset_reference` on mesh graphs
  - `grip_torus_preset_reference` on torus graphs
  - `grip_tree_preset_reference` on tree graphs
- `igraph_kk_default`

Special handling for `igraph_kk_default`:
- attempt on all Stage-2 graphs up to a dedicated KK scale cap
- record explicit `skipped` or `timeout` rows outside that cap
- use `dim = 2` or `3` to match the graph specification

## Metrics

Per-run metrics:
- Procrustes RMSE to the canonical graph embedding
- edge-length coefficient of variation
- sampled stress
- sampled non-edge separation ratio
- elapsed time
- run status: `ok`, `error`, `timeout`, or `skipped`

Derived scores:
- `quality_score`
  - same weighted rank score used in earlier GRIP family benchmarks
  - weights:
    - RMSE `0.45`
    - edge CV `0.20`
    - sampled stress `0.20`
    - non-edge separation `0.15`
- `value_score`
  - `quality_score` plus a runtime-rank component
  - runtime is treated as lower-is-better
  - this is included because KK quality may be excellent on some graphs while remaining far slower

## Seeding Policy

GRIP:
- `n <= 500`: seeds `1:3`
- `500 < n <= 5000`: seeds `1:2`
- `n > 5000`: seed `1`

KK:
- seed `1` in Stage 2
- reason:
  - earlier work showed deterministic behavior on the carpet
  - the cross-family benchmark is primarily about method tradeoffs, not KK seed distributions

## Runtime Policy

- Stage-2 global graph cap: `100000` vertices
- Stage-2 KK cap: `10000` vertices
- timeouts:
  - GRIP candidates: adaptive elapsed time limit by size
  - KK: subprocess timeout via `callr`

## Deliverables

- experiment action plan
- benchmark script
- Stage-1 summary CSV, markdown, and plot
- Stage-2 raw CSV, graph summary CSV, family summary CSV, candidate summary CSV
- Stage-2 markdown summary
- quality-vs-runtime scatter
- family-winner heatmap or equivalent matrix plot
- tracked high-level summary noting:
  - best universal setting
  - best per-family settings
  - how `grip_f32` compares with KK on the comparable subset
