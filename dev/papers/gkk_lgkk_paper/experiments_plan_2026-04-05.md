# GKK LGKK Paper Experiments, Tests, Figures, and Tables Plan

Date: 2026-04-05

This document is intentionally coupled to:

- `dev/papers/gkk_lgkk_paper/outline_2026-04-05.md`

The idea is simple:

- the outline states what the paper wants to claim
- this document states what evidence is already in hand
- and what experiments, tests, figures, and tables still need to be produced before the paper is ready

## Status Legend

- `DONE` means the evidence or artifact already exists in a usable form
- `PARTIAL` means something real exists, but it is not yet at paper quality
- `TODO` means the evidence is still missing
- `OPTIONAL` means it would strengthen the paper but is not required for a minimal strong version
- `CONDITIONAL` means include it only if the result is crisp enough to justify main-text space

## Current Evidence Already In Hand

### Existing code and tests

- `TEST-1` `DONE`
  - targeted validation run passed `212` checks across:
    - `geodesic-kk`
    - `landmark-geodesic-kk`
    - `layout-globalrep`
    - `layout-trace`
    - `layout-weighted-globalrep`
    - `layout-weighted-trace`
  - role in paper:
    - supports the implementation/reproducibility section
    - should appear in a short validation paragraph or appendix table

- `TEST-2` `DONE`
  - dedicated `MISF-GKK` scaffold tests currently pass (`32` checks)
  - current coverage includes:
    - prepared-object construction
    - immediate top-level solve
    - scoring against external coordinates
    - landmark scoring mode
    - multistage optimizer output with traces
  - role in paper:
    - supports the claim that `MISF-GKK` is a complete public family, not just an internal prototype

### Existing benchmark evidence

- `EXP-1` mesh benchmark `PARTIAL`
  - already run and documented
  - current scope:
    - 3 surface families
    - 3 sizes
    - 2 embedding dimensions
    - 3 seeds
    - 54 shared-start benchmark cases
  - strongest current result:
    - in 3D, mean full-`GKK` relative RMSE is about `0.0254` for `KK`
    - about `0.0067` for `KK -> GKK`
    - about `0.0069` for `KK -> LGKK`
    - `LGKK` is about seven times faster than full `GKK`
  - current artifacts:
    - benchmark data under `output/gkk_lgkk_paper/benchmarks/kk-gkk-lgkk-mesh-suite-2026-03-31/`
    - rendered report under `output/gkk_lgkk_paper/reports/benchmarks/kk_gkk_lgkk_mesh_benchmark_report_2026-03-31/`

- `EXP-2` cross-family weighted smoke panel `PARTIAL`
  - already run at smoke scale
  - current 3D aggregate over five families:
    - `KK`: mean full-`GKK` relative RMSE about `0.0448`
    - `KK -> GKK`: about `0.0280`
    - `KK -> LGKK`: about `0.0228`
  - this is encouraging, but it is still smoke-level evidence
  - current artifacts:
    - benchmark data under `output/gkk_lgkk_paper/benchmarks/weighted-grip-phase5-smoke-check/`
    - summary report under `output/gkk_lgkk_paper/reports/benchmarks/weighted_grip_phase5_report_2026-04-02/`

### Existing multiscale evidence

- `EXP-6` `MISF-GKK` method-family implementation `PARTIAL`
  - the public `prepare / optimize / score` family now exists
  - public documentation exists
  - basic scaffold tests pass
  - a concrete benchmark driver now exists:
    - `tools/benchmarks/gkk_lgkk_paper/benchmark-misf-gkk-panel.R`
  - a persistent smoke artifact now exists:
    - `output/gkk_lgkk_paper/tmp/misf-gkk-panel-smoke-2026-04-05/`
  - current smoke status:
    - the driver runs end to end
    - but the current default `MISF-GKK` configuration is not yet competitive with direct `KK -> LGKK` or exact `KK -> GKK` on the reduced overlap panel
  - what is still missing is paper-grade benchmark evidence and likely one tuning / ablation pass rather than implementation existence

- `EXP-8` weighted `GRIP` comparison `PARTIAL`
  - weighted `GRIP` plus in-core and post-polish `LGKK` already exists
  - current evidence is useful for internal direction
  - current evidence is not yet clean enough to guarantee a main-text section in the journal version

## Minimum Evidence Bundle For A Strong Paper

These are the pieces I think we need before calling the paper submission-ready.

1. `EXP-0` toy motivation example
2. `EXP-1` final mesh benchmark section
3. `EXP-2` final cross-family benchmark section
4. `EXP-6` a paper-grade `MISF-GKK` multiscale benchmark
5. `EXP-3` tie-handling illustration or appendix note
6. `EXP-5` convergence/runtime quality tradeoff evidence
7. `TEST-1` and `TEST-2` implementation validation summary

`EXP-4` and `EXP-8` would strengthen the paper further, but the paper can still work without them if the core experiments are strong.

## Experiment Inventory

### EXP-0 Toy chord-versus-path motivation figure

Status:

- `PARTIAL`

Question:

- what is the smallest example that makes the `KK` chord-versus-path mismatch obvious?

What to build:

- one small weighted graph with a clear geometric interpretation
- one drawing where two vertices are close by chord but far by embedded path, or vice versa

Deliverables:

- `Figure 1`
- one paragraph in the introduction

Notes:

- this does not need to be a large experiment
- it needs to be visually and conceptually clean

### EXP-1 Final mesh benchmark

Status:

- `PARTIAL`

Paper question:

- on the cleanest weighted surface families, how much do `GKK` and `LGKK` improve geodesic fidelity over plain `KK`, and what runtime do they cost?

Comparisons:

- `KK`
- `KK -> GKK`
- `KK -> LGKK`

Recommended setup:

- keep the current three mesh families
- keep the current three sizes
- keep both 2D and 3D
- keep at least three seeds
- shared-start initialization is acceptable, but state clearly that `GKK` and `LGKK` are warm-started polishers in these experiments

Already in place:

- benchmark run
- summary CSVs
- internal report
- gallery figure

Still needed:

- convert the report outputs into final paper figures/tables
- rerun only if implementation changed materially
- produce one concise case gallery figure suitable for the manuscript
- add convergence traces if feasible

Primary metrics:

- full-`GKK` relative RMSE
- runtime
- secondary:
  - classical `KK` relative RMSE
  - Procrustes RMSE
  - sampled stress

Figures:

- `Figure 3` mesh quality comparison
- `Figure 4` mesh runtime-quality tradeoff
- `Figure 5` representative mesh gallery

Tables:

- `Table 1` mesh benchmark method summary
- `Table 2` per-case best-method summary or compact case summary

### EXP-2 Final cross-family weighted geometric benchmark

Status:

- `TODO`

Paper question:

- do the mesh conclusions generalize to a broader panel of weighted geometric families?

Core comparisons:

- `KK`
- `KK -> GKK`
- `KK -> LGKK`

Recommended family panel:

- mesh
- torus
- sphere
- Sierpinski carpet
- cube channel network
- irregular annulus
- irregular torus

Recommended dimensions:

- 3D is primary
- 2D is secondary and can be summarized more compactly

Recommended seeds:

- at least `3`
- more if runtime is acceptable

Recommended outputs:

- aggregate summary CSV
- per-family summary CSV
- representative layout gallery
- one compact final report

Still needed:

- rerun the family panel in non-smoke mode
- finalize the exact family list
- decide whether all families use the same size scale or matched complexity bands

Primary metrics:

- full-`GKK` relative RMSE
- runtime

Secondary metrics:

- Procrustes RMSE
- sampled stress

Figures:

- `Figure 6` cross-family aggregate quality comparison
- `Figure 7` cross-family runtime-quality plot
- optional family gallery panel

Tables:

- `Table 3` family panel definition
- `Table 4` cross-family method summary

### EXP-3 Tie-handling illustration

Status:

- `TODO`

Paper question:

- when shortest paths are highly tied, does `single` versus `average` tie handling materially affect the realized path metric?

Scope:

- one or two small carefully chosen examples
- not a huge benchmark section

Preferred outcome:

- a compact figure or mini-table showing a symmetric graph where tie averaging removes an artifact or changes the realized path computation in an interpretable way

Editorial use:

- if strong and clear, keep in the main text
- otherwise move to an appendix

Figure:

- `Figure 2` tie-handling or sparse-pair illustration, depending on which is cleaner

Table:

- optional appendix mini-table only

### EXP-4 LGKK sparsity and parameter sensitivity

Status:

- `OPTIONAL`

Paper question:

- how sensitive is `LGKK` to `local_nbrs` and `landmark_count`, and where is the quality/runtime knee?

Scope:

- one mesh family plus one irregular family is enough
- no need for a full benchmark sweep across every family

Why it matters:

- strengthens the claim that `LGKK` is a practical sparse approximation rather than a narrowly tuned construction

Figures:

- optional `Figure 8` quality/runtime versus landmark budget

Tables:

- optional appendix parameter summary

### EXP-5 Convergence and runtime profile

Status:

- `TODO`

Paper question:

- how quickly do `GKK` and `LGKK` improve over their shared starts, and what is the runtime tradeoff?

Scope:

- use the mesh benchmark first
- optionally reuse one or two cross-family cases

Desired outputs:

- objective or error versus iteration
- runtime versus final quality

Why this matters:

- strengthens the practical story of `LGKK`
- makes the full `GKK` versus sparse `LGKK` tradeoff more explicit than endpoint tables alone

Figures:

- can be folded into `Figure 4`
- or split into a dedicated convergence figure if space allows

### EXP-6 MISF-GKK multiscale benchmark

Concrete spec:

- see `dev/papers/gkk_lgkk_paper/exp6_misf_gkk_benchmark_spec_2026-04-05.md`

Status:

- `TODO`

Paper question:

- does the new standalone `MISF-GKK` family provide a practical multiscale path from coarse to fine while preserving most of the geodesic-quality gain of direct `GKK/LGKK`?

Recommended comparisons:

- `MISF-GKK` with `pair_mode = auto`
- `MISF-GKK` with stronger exact settings on smaller graphs
- `KK -> GKK`
- `KK -> LGKK`
- optionally plain `KK`

Recommended editorial rule:

- this experiment should stay small and focused
- it exists to justify the multiscale contribution as a real method-family result

Recommended graph sizes:

- overlap regime:
  - mesh `12 x 12`, `16 x 16`
  - torus `12 x 12`, `16 x 16`
  - irregular annulus `8 x 32`, `10 x 40`
  - irregular torus `9 x 18`, `12 x 24`
- scale regime:
  - mesh `20 x 20`
  - torus `20 x 20`
  - irregular annulus `12 x 48`
  - irregular torus `14 x 28`

Exact seeds:

- `1, 2, 3`

Exact main comparison methods:

- overlap regime:
  - `KK`
  - `KK->GKK`
  - `KK->LGKK`
  - `MISF-GKK-auto`
  - `MISF-GKK-full`
- scale regime:
  - `KK`
  - `KK->LGKK`
  - `MISF-GKK-auto`

What this section should establish:

- quality retention relative to direct `GKK/LGKK`
- runtime advantage as graph size grows
- sensible behavior of `auto` pair-mode switching

Current evidence:

- implementation and tests exist
- a concrete benchmark driver now exists
- a smoke benchmark has been run successfully
- current smoke results are cautionary rather than manuscript-ready:
  - the method runs cleanly
  - the present default settings underperform direct `KK -> LGKK` and `KK -> GKK` on the reduced panel

Still needed:

- run one targeted tuning pass before the full benchmark
- run a disciplined benchmark with final settings
- summarize `auto` versus explicit pair-mode choices

Figures:

- one main-text scalability figure
- one runtime-quality figure if space allows

Tables:

- one compact multiscale summary table

### EXP-7 MISF-GKK ablation study

Status:

- `TODO`

Paper question:

- which design choices in `MISF-GKK` matter most?

Recommended ablations:

- top-level pair mode: `full` vs `landmark` vs `auto`
- refinement pair mode: `full` vs `landmark` vs `auto`
- final polish pair mode
- insertion mode:
  - `geodesic`
  - `weighted_kk`
  - optionally `weighted_grip`

Scope:

- one mesh family plus one irregular family is enough
- keep graph sizes moderate

Why it matters:

- this is the cleanest way to make `MISF-GKK` scientifically legible
- otherwise the new family risks reading like a bag of engineering choices

Figures:

- optional ablation figure or appendix panel

Tables:

- appendix ablation table strongly recommended

### EXP-8 Weighted GRIP comparison section

Status:

- `CONDITIONAL`

Paper question:

- once `MISF-GKK` exists as the cleaner multiscale geodesic family, does weighted `GRIP` still deserve space in the main paper?

Recommended comparisons:

- weighted `GRIP`
- weighted `GRIP + core LGKK`
- weighted `GRIP + polish LGKK`
- `MISF-GKK`
- optionally `KK -> LGKK` as a quality reference

Editorial rule:

- keep this section only if it adds clarity rather than diluting the paper
- it should be shorter than the `MISF-GKK` section

Use:

- comparison section or appendix
- not the primary scaling section

### TEST-1 General implementation validation

Status:

- `DONE`

What it supports:

- reproducibility and implementation confidence

What should appear in the paper:

- one short paragraph in the implementation section
- or one appendix table summarizing covered components

Current summary:

- `212` targeted checks currently pass across the core `GKK`, `LGKK`, and weighted/combinatorial layout pathways

### TEST-2 MISF-GKK family validation

Status:

- `DONE`

What it supports:

- the new public `MISF-GKK` family is implemented end-to-end

What should appear in the paper:

- one sentence in the implementation section
- or one appendix row in the validation table

Current summary:

- `32` dedicated scaffold checks currently pass for `prepare`, `optimize`, `score`, trace output, and landmark scoring mode

## Figure Plan

### Main-text figure list

- `Figure 1` toy chord-versus-path mismatch example
  - source: `EXP-0`
  - status: `TODO`

- `Figure 2` either:
  - shortest-path tie-handling illustration
  - or `LGKK` sparse pair construction figure
  - source: `EXP-3` or `EXP-4`
  - status: `TODO`

- `Figure 3` mesh benchmark quality figure
  - source: `EXP-1`
  - status: `PARTIAL`

- `Figure 4` mesh benchmark runtime-quality or convergence figure
  - source: `EXP-1` plus `EXP-5`
  - status: `TODO`

- `Figure 5` representative mesh gallery
  - source: `EXP-1`
  - status: `PARTIAL`

- `Figure 6` cross-family benchmark aggregate comparison
  - source: `EXP-2`
  - status: `TODO`

- `Figure 7` cross-family runtime-quality tradeoff or 2D versus 3D summary
  - source: `EXP-2`
  - status: `TODO`

- `Figure 8` `MISF-GKK` scalability or multiscale comparison figure
  - source: `EXP-6`
  - status: `TODO`

### Optional figure list

- `Figure 9` `LGKK` sensitivity or landmark-budget tradeoff
  - source: `EXP-4`
  - status: `OPTIONAL`

- `Figure 10` weighted-`GRIP` comparison figure
  - source: `EXP-8`
  - status: `CONDITIONAL`

## Table Plan

### Main-text tables

- `Table 1` mesh benchmark method summary
  - source: `EXP-1`
  - status: `PARTIAL`

- `Table 2` mesh benchmark per-case or per-family best-method summary
  - source: `EXP-1`
  - status: `PARTIAL`

- `Table 3` cross-family benchmark panel definition
  - source: `EXP-2`
  - status: `TODO`

- `Table 4` cross-family benchmark aggregate summary
  - source: `EXP-2`
  - status: `TODO`

- `Table 5` `MISF-GKK` multiscale benchmark summary
  - source: `EXP-6`
  - status: `TODO`

### Appendix tables

- `Appendix Table A1` validation/test coverage summary
  - source: `TEST-1` and `TEST-2`
  - status: `DONE`

- `Appendix Table A2` optional `LGKK` sensitivity summary
  - source: `EXP-4`
  - status: `OPTIONAL`

- `Appendix Table A3` `MISF-GKK` ablation summary
  - source: `EXP-7`
  - status: `TODO`

## Practical Run Order

If we want the fastest route to a paper-ready evidence base, I would run the missing work in this order:

1. `EXP-0`
   - quick conceptual win for the introduction
2. `EXP-1`
   - finalize the strongest section first
3. `EXP-2`
   - this is the most important missing paper-grade benchmark
4. `EXP-6`
   - this determines whether the paper can honestly claim a full multiscale method-family contribution
5. `EXP-5`
   - makes the `LGKK` practical story sharper
6. `EXP-3`
   - decide main text versus appendix
7. `EXP-7`
   - helps make `MISF-GKK` scientifically legible
8. `EXP-8`
   - only if we still want weighted `GRIP` in the main paper
9. `EXP-4`
   - useful, but easiest to trim if schedule gets tight

## Decision Gates

### Gate A: core paper readiness

The paper is ready for a strong core draft once these are complete:

- `EXP-0`
- `EXP-1`
- `EXP-2`
- `EXP-6`
- `EXP-5`
- `TEST-1`
- `TEST-2`

### Gate B: whether weighted GRIP stays in the main text

Keep weighted `GRIP` in the main paper only if:

- `EXP-8` produces a clean and short story
- and the section does not blur the main `GKK/LGKK/MISF-GKK` identity

Otherwise:

- move weighted `GRIP` comparison to an appendix
- or mention it briefly in the discussion and reserve a fuller treatment for later work

## Bottom Line

The paper already has enough substance to justify a strong `GKK + LGKK` direction.

What it still needs is not a new conceptual center, but a tighter evidence package:

1. finalize the mesh section
2. run the full cross-family benchmark
3. run a real `MISF-GKK` multiscale benchmark
4. add one or two compact supporting studies
5. keep weighted `GRIP` only if the evidence stays crisp
