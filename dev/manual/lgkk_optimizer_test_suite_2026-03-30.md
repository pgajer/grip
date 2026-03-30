# LGKK Optimizer Test Suite

## Goal

Evaluate the first experimental LGKK optimizer as an opt-in post-polish for
`grip.layout()`, with emphasis on three questions:

1. Does the optimizer improve layout quality on the structured graphs that have
   been driving the GRIP/KK investigation?
2. Can one bounded parameter setting work well across both unweighted and
   weighted graphs?
3. How expensive is the polish relative to the current GRIP default and to
   `igraph::layout_with_kk()`?

## Why This Suite

The new optimizer introduces two real algorithmic knobs:

- `lgkk_local_nbrs = q`
- `lgkk_landmark_count = M`

and one budget knob:

- `lgkk_polish_rounds`

The suite therefore separates:

- a bounded search stage used to identify promising `(rounds, q, M)` settings;
- a confirmation stage used to test those settings on larger or different
  graphs.

This keeps the experiment informative without exploding runtime on the level-4
carpet.

## Graph Set

### Stage 1: Training / Bounded Search

- Sierpinski carpet level 3
- Weighted mesh 6x6 with horizontal edge weight `1` and vertical edge weight
  `2`, and a stretched rectangular target embedding

These two graphs were chosen deliberately:

- the level-3 carpet is the fastest graph that still exposes the symmetry and
  rectilinearity issues motivating the KK/LGKK work;
- the weighted mesh is the smallest nontrivial graph where the weighted-graph
  motivation of LGKK can be tested directly.

### Stage 2: Confirmation

- Sierpinski carpet level 3
- Sierpinski carpet level 4
- Sierpinski triangle level 4
- Unweighted mesh 8x8
- Weighted mesh 6x6

## Methods Compared

### Stage 1 Candidate Grid

- Baseline `grip.layout()` with `lgkk_polish_rounds = 0`
- LGKK grid:
  - `lgkk_polish_rounds in {2, 4, 8}`
  - `lgkk_local_nbrs in {6, 12, 20}`
  - `lgkk_landmark_count in {4, 8}`
- `igraph::layout_with_kk()` reference

### Stage 2 Confirmation Methods

- Baseline `grip.layout()`
- `igraph::layout_with_kk()`
- Best LGKK quality candidate from Stage 1
- Best LGKK value candidate from Stage 1

## Metrics

All methods are scored with:

- runtime
- Procrustes RMSE to the canonical/target embedding
- LGKK relative weighted RMSE
- edge-axis deviation
- edge-length CV

Carpet-specific metrics:

- central-hole skew
- central-hole aspect error
- central-hole center error

## Selection Rule

Stage 1 produces two recommended LGKK candidates:

- `best_quality`: lowest mean rank across
  `procrustes.rmse`, `lgkk.weighted.rel.rmse`, and `edge.axis.deviation`
- `best_value`: lowest mean rank across the same three quality metrics plus
  runtime

This separation matters because an LGKK setting may be geometrically better but
 too expensive to be a practical default.

## Deliverables

- raw and summarized CSV files for Stage 1 and Stage 2
- Stage 1 heatmaps and runtime/quality frontier plot
- Stage 2 per-graph comparison PNGs
- a LaTeX report consolidating the benchmark design, results, and
  recommendations
