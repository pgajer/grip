# LGKK Weighted-Graph Fix Results

## Root Cause

The compiled weighted shortest-path cache in
`src/MishSupport.cpp::compute_active_shortest_paths()` used a scaled tolerance
test

`alt + tol * scale < best`

with

`scale = max(1, abs(alt), abs(best))`.

When `best` was still `Inf`, `scale` also became `Inf`, so the improve test was
never true. As a result:

- no weighted shortest-path distances were ever recorded,
- the compiled LGKK pair cache was empty on weighted graphs,
- the compiled multiscale LGKK stage accepted no steps and became a no-op.

## Fix

The weighted Dijkstra branch now special-cases non-finite `best` values:

- `scale` ignores `best` when `best` is not finite,
- `improve` is immediately true when `best` is not finite.

I also kept the earlier deterministic queue-order improvement so equal-distance
frontier ties are resolved by actual vertex id rather than active-set index.

## Validation

- `devtools::test(filter='layout-(globalrep|trace)|landmark-geodesic-kk')`
  passed with `121` tests.
- Focused benchmark rerun:
  [summary](/Users/pgajer/current_projects/grip/dev/manual/tmp/lgkk-multiscale-linesearch-2026-03-30/lgkk-multiscale-integration-summary.md)

## Key Results

### Weighted mesh `6x6`

- Baseline: RMSE `0.0782`, LGKK rel. RMSE `0.1367`, `0.024s`
- Compiled LGKK x4 after fix: RMSE `0.0295`, LGKK rel. RMSE `0.0503`, `0.030s`
- R LGKK polish x4: RMSE `0.0292`, LGKK rel. RMSE `0.0518`, `0.082s`

So the compiled stage now improves the weighted graph strongly and lands very
close to the R post-polish result, while running about `2.7x` faster on this
benchmark.

### Level-4 carpet

- Baseline: RMSE `0.0419`, axis deviation `0.0246`, skew `0.0776`, `2.701s`
- Compiled LGKK x4 after fix: RMSE `0.0409`, axis deviation `0.0174`, skew
  `0.0770`, `16.820s`
- R LGKK polish x4: RMSE `0.0409`, axis deviation `0.0174`, skew `0.0770`,
  `60.214s`

So the carpet improvement was preserved.

## Recommendation

- Keep the weighted-cache fix.
- Keep the compiled LGKK multiscale stage experimental, but it is now viable on
  both the carpet and weighted-mesh benchmarks.
- Drop `balanced_band` as the current insertion fix; it still worsens the
  carpet skew defect at the first fully active frame.
