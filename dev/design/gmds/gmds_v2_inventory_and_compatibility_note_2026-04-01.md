# GMDS v2 Inventory and Compatibility Note

Date: 2026-04-01

This note records the Step 1 compatibility pass for the GMDS v2 cleanup.

## Goal

Before changing semantics, identify:

- the current public GMDS entry points,
- the main internal assumptions baked into tests and reports,
- the compatibility strategy for introducing the manuscript-aligned graph-first
  API.

## Current Public Surface

The currently exported GMDS-facing functions are:

- [`grip.prepare.geodesic.mds()`](/Users/pgajer/current_projects/grip/R/grip_quality.R#L3363)
- [`grip.score.geodesic.mds()`](/Users/pgajer/current_projects/grip/R/grip_quality.R#L3472)
- [`grip.optimize.geodesic.mds()`](/Users/pgajer/current_projects/grip/R/grip_quality.R#L3649)

Their current behavior is:

- `grip.prepare.geodesic.mds()` is data-native and builds a symmetric weighted
  `k`-NN graph before preparing the all-pairs path cache.
- `grip.score.geodesic.mds()` and `grip.optimize.geodesic.mds()` accept a
  prepared object or, if omitted, construct one through the data-native path.

## Existing Graph-Native Building Blocks

The graph-native low-level machinery already exists:

- graph input normalization:
  [`grip.prepare.geodesic.kk.base()`](/Users/pgajer/current_projects/grip/R/grip_quality.R#L1038)
- all-pairs chosen-path cache:
  [`grip.prepare.geodesic.kk()`](/Users/pgajer/current_projects/grip/R/grip_quality.R#L1820)

That means the new manuscript-aligned graph-first GMDS API can be added without
rewriting the core path-cache logic.

## Compatibility-Sensitive Call Sites

At the time of this pass, the directly relevant call sites were:

### Tests

- [`test-geodesic-mds.R`](/Users/pgajer/current_projects/grip/tests/testthat/test-geodesic-mds.R)

This file assumed:

- `grip.prepare.geodesic.mds()` is deterministic,
- `grip.optimize.geodesic.mds(data = ..., k = ...)` uses the data-native path,
- prepared objects from the data-native path carry class
  `grip_gmds_prepared`.

### Documentation

- [`grip.prepare.geodesic.mds.Rd`](/Users/pgajer/current_projects/grip/man/grip.prepare.geodesic.mds.Rd)
- [`grip.score.geodesic.mds.Rd`](/Users/pgajer/current_projects/grip/man/grip.score.geodesic.mds.Rd)
- [`grip.optimize.geodesic.mds.Rd`](/Users/pgajer/current_projects/grip/man/grip.optimize.geodesic.mds.Rd)

These pages still described the data-native entry point as if it were the only
public preparation path.

### Reports and design notes

The current design and experiment notes already point to a graph-first
understanding of GMDS, but the code-level naming had not caught up fully. The
main relevant notes are:

- [geodesic_mds_infrastructure_design_2026-03-31.md](/Users/pgajer/current_projects/grip/dev/design/geodesic_mds_infrastructure_design_2026-03-31.md)
- [gmds_algorithm_implementation_and_pathology_overview_2026-03-31.md](/Users/pgajer/current_projects/grip/dev/design/gmds_algorithm_implementation_and_pathology_overview_2026-03-31.md)
- [gmds_v2_four_phase_cleanup_design_2026-04-01.md](/Users/pgajer/current_projects/grip/dev/design/gmds_v2_four_phase_cleanup_design_2026-04-01.md)

## Compatibility Strategy

The chosen compatibility strategy is:

1. add a new graph-first public function rather than renaming the existing one,
2. preserve `grip.prepare.geodesic.mds(data, k, ...)` exactly as a convenience
   wrapper,
3. route the data-native function through the new graph-first implementation,
4. keep the returned class hierarchy unchanged so scoring and optimization code
   continue to work without modification.

Concretely:

- new function:
  `grip.prepare.graph.geodesic.mds(...)`
- existing function retained:
  `grip.prepare.geodesic.mds(data, k, ...)`

This minimizes breakage while aligning the code with the manuscript's
graph-first formulation.

## What Should Remain Stable

The following behavior should remain stable after Step 2:

- existing user code that calls `grip.prepare.geodesic.mds(data, k, ...)`
- existing user code that calls `grip.optimize.geodesic.mds(data = ..., k = ...)`
- prepared objects continuing to inherit `grip_gmds_prepared`
- the internal all-pairs cache contents for a given graph and `tie_mode`

## Step 2 Verification Criterion

The key compatibility test for Step 2 is:

`If we build the same weighted graph explicitly and prepare it through the new graph-first API, the common GMDS cache fields should match the current data-native wrapper exactly.`

That criterion is now the basis of the new regression test added alongside the
Step 2 code changes.
