# GRIP Coarse Global-Repulsion Layout Manifest

Date: 2026-03-27

## Summary

This document specifies a new experimental GRIP layout variant that preserves the current `grip.layout()` behavior and introduces a sibling API with additional coarse-level anti-folding repulsion.

The new algorithm keeps the existing multiscale GRIP structure:

- MISF construction is unchanged.
- Initial placement of newly introduced vertices is unchanged.
- Coarse and intermediate refinement still use local Kamada-Kawai (KK).
- Final full-graph refinement still uses the current local Fruchterman-Reingold (FR) force.

The only algorithmic change is:

- on coarse MISF levels (`misfLevel > 0`), add an extra repulsive term over the currently active vertices
- use exact global repulsion when the active set is small
- switch to sampled global repulsion when the active set is larger

## Motivation

Local KK refinement can preserve neighborhood structure well, but because it only sees the retained local graph-distance neighborhood, distant active regions may drift into each other and create foldovers or self-overlaps before the finest level is reached.

Replacing coarse KK with pure FR is not the preferred first intervention because:

- the current FR implementation is also local, not globally repulsive
- coarse MISF levels do not yet contain all original vertices, so naive FR-on-original-edges is poorly anchored

This design addresses the observed failure mode more directly by adding a weak, geometry-aware, active-set-wide repulsion term while keeping the existing KK scaffold intact.

## Public API

Add a new exported R function:

- `grip.layout.globalrep()`

This function mirrors the existing `grip.layout()` interface and disconnected-component handling, but adds the following tuning arguments:

- `coarse_repulsion_factor`
- `coarse_repulsion_sample`
- `coarse_repulsion_exact_below`

Initial defaults:

- `coarse_repulsion_factor = 0.2`
- `coarse_repulsion_sample = 16`
- `coarse_repulsion_exact_below = 128`

The original `grip.layout()` and `grip.layout.trace()` APIs remain unchanged.

## Algorithm Definition

Let `Vi` denote the current active MISF level and `csize = |Vi|`.

At each coarse refinement step (`misfLevel > 0`):

1. Compute the existing local KK displacement using the retained graph-distance neighborhood `Ni(v)`.
2. Compute an additional repulsive displacement using only currently active vertices in `Vi`.
3. Add the repulsive displacement to the KK displacement before temperature scaling and normalization.

At the final level (`misfLevel == 0`):

- preserve the current `FR_spring()` behavior exactly

## Global Repulsion Policy

For a vertex `v` at a coarse level:

- if `csize <= coarse_repulsion_exact_below`, compute repulsion against every other active vertex in `Vi`
- otherwise, sample up to `coarse_repulsion_sample` distinct active vertices uniformly without replacement and scale the resulting sum by `(csize - 1) / sample_size`

Repulsion is never computed against inactive vertices.

The repulsive force uses the same inverse-square geometric form already used by the current FR-style repulsion:

- direction is `pos[v] - pos[u]`
- magnitude scales with `1 / ||pos[v] - pos[u]||^2`

Its strength is controlled by:

- `coarse_repulsion_factor * 0.05 * edge^2`

which mirrors the existing finest-level `repulsion_factor` scaling convention.

## Determinism

All sampling must use the existing graph RNG (`fast_Rand()` seeded through the current `seed` flow) so that fixed seeds remain deterministic.

## Scope Boundaries

Included in this slice:

- new experimental layout API
- new C++ engine path for coarse global repulsion
- validation and tests for the new API
- Rd/NAMESPACE/Rcpp export generation

Explicitly not included in this slice:

- any change to `grip.layout()` behavior
- a trace variant of the new algorithm
- new presets tuned for the new algorithm
- documentation or vignette updates beyond function-level docs
- coarse-graph reconstruction or Barnes-Hut style acceleration

## Compatibility Requirements

- `grip.layout()` must return byte-for-byte identical results to its current behavior for the same seed and arguments
- disconnected-component handling for the new API should match `grip.layout()`
- setting `coarse_repulsion_factor = 0` in the new API should recover the current GRIP layout behavior

## Validation Requirements

The implementation should be verified with focused tests for:

- finite output matrix with expected dimensions
- deterministic seeded runs
- parameter validation for new arguments
- changed output when coarse repulsion is enabled
- equality with `grip.layout()` when `coarse_repulsion_factor = 0`
- disconnected-component packing parity with the existing layout wrapper

## Follow-Up Work

After this slice lands, the next most useful follow-ups are:

1. Add `grip.layout.globalrep.trace()`
2. Build tuning scripts on likely foldover-prone families
3. Decide whether this should eventually become a preset or remain a sibling API
