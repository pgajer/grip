# Weighted GRIP Phase 5

## Goal

Phase 5 completes two pieces that were intentionally left for after the core
weighted GRIP engine stabilized:

- bring weighted multiscale LGKK rounds inside the weighted GRIP core, and
- broaden the benchmark rollout beyond the original mesh/cylinder/torus/sphere
  panel.

## What Changed

### 1. In-core weighted multiscale LGKK refinement

The weighted engine now runs the same compiled multiscale LGKK level-refinement
hook that the combinatorial globalrep engine already uses.

This means:

- `grip.layout.globalrep.weighted()`
- `grip.layout.weighted()`
- `grip.layout.trace.weighted()`

now support the full multiscale LGKK control set:

- `lgkk_multiscale_rounds`
- `lgkk_rounds_coarse`
- `lgkk_rounds_pre_final`
- `lgkk_rounds_final`
- `lgkk_local_nbrs`
- `lgkk_landmark_count`
- `lgkk_multiscale_scope`
- `lgkk_active_limit`

The weighted wrappers no longer expose only post-layout LGKK polish. They can
now refine layouts during the weighted GRIP hierarchy itself.

### 2. Weighted trace coverage for multiscale LGKK

The weighted trace API now records the in-core LGKK phase as part of the traced
trajectory and is tested against the final weighted layout result for both:

- shared multiscale LGKK budgets, and
- staged coarse / pre-final / final budgets.

### 3. Broader family benchmark rollout

The family benchmark driver now covers a wider Phase 5 panel:

- mesh
- cylinder
- torus
- sphere
- Sierpinski carpet
- porous cube channel network
- irregular annulus
- irregular torus

It also distinguishes two weighted-LGKK variants:

- `Weighted GRIP + core LGKK`
- `Weighted GRIP + polish LGKK`

This makes it possible to compare:

- combinatorial GRIP,
- weighted GRIP,
- weighted GRIP with in-core multiscale LGKK,
- weighted GRIP with post-layout LGKK polish,
- and the KK / GKK / LGKK family.

## Validation

Phase 5 adds and checks:

- weighted globalrep tests for multiscale LGKK knob effects,
- weighted globalrep tests for staged LGKK round budgets,
- weighted trace tests that final traced coordinates match weighted layout
  coordinates under multiscale LGKK,
- weighted trace checks that the trace metadata contains an `lgkk` phase,
- and a refreshed benchmark smoke run for the broadened family panel.
