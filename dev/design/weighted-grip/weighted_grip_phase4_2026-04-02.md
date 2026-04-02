# Weighted GRIP Phase 4

## Goal

Phase 4 focuses on the weighted shortest-path hot path inside the weighted
GRIP core.

The main target is the per-vertex weighted neighborhood construction used by:

- weighted initial MISF seeding
- weighted insertion
- weighted local refinement

## What Changed

### 1. Settled-order weighted search

The weighted insertion path no longer computes a full all-vertex distance
vector and then sorts it just to discover the nearest useful weighted
neighbors.

Instead, it now consumes vertices directly in Dijkstra settled order and fills:

- weighted neighborhood caches
- insertion-anchor candidates

as soon as they become available.

This removes one of the largest avoidable costs in the weighted core:

- full distance-vector initialization
- full reachable-vertex sorting for every inserted vertex

### 2. Scratch-state reuse

The weighted Dijkstra traversal now reuses internal scratch arrays instead of
reallocating and reinitializing them on every weighted search.

This keeps the exact weighted path semantics while cutting repeated setup cost.

### 3. Exact early stopping

For inserted vertices, the exact weighted search now stops as soon as:

- the required weighted neighborhood caches are filled, and
- the insertion-anchor selection has enough information to finalize.

That is still exact for the weighted insertion/cache-building task because the
search frontier is processed in nondecreasing weighted shortest-path order.

### 4. Optional approximate weighted neighborhoods

The weighted public API now accepts:

- `metric_neighbor_cap`

When left as `NULL`, the weighted search uses the exact Phase 4 path above.

When set to a positive integer, the weighted search stops after that many
settled Dijkstra vertices while building weighted neighborhood caches for
inserted vertices. This gives an explicit approximation mode for larger graphs
without changing the default exact semantics.

This cap currently applies only to weighted neighborhood construction during
vertex insertion. It does not change:

- weighted MISF radius searches
- the initial exact diameter-style pass used to size the layout box

## Public API Impact

Phase 4 extends:

- `grip.layout.globalrep.weighted()`
- `grip.layout.weighted()`
- `grip.layout.trace.weighted()`

with the optional argument:

- `metric_neighbor_cap = NULL`

## Validation

Phase 4 adds tests that verify:

- sufficiently large `metric_neighbor_cap` reproduces the exact default result
- approximate weighted runs remain deterministic for fixed seeds
- weighted trace stays consistent with weighted layout under the cap
- invalid cap values are rejected cleanly
