# GRIP vs KK Cross-Family Benchmark Summary

Run tag: `grip-kk-cross-family-full-2026-03-28`

## Stage 1

- screened universal candidates: `grip_default_adaptive`, `grip_f32`, `grip_f64`, `grip_f96`, `grip_f128`, `grip_f192`, `grip_f384`, `grip_exact_f64`, `grip_exact_f96`
- selected universal candidates for Stage 2: `grip_f96`, `grip_f128`, `grip_exact_f96`, `grip_f64`, `grip_f32`, `grip_default_adaptive`
- Stage-1 winner: `grip_f96` (mean quality score 0.3831)

## Stage 2

- global quality winner: `grip_torus_preset_reference` (mean quality score 0.2214)
- global value winner: `igraph_kk_default` (mean value score 0.2998)
- note: the torus/tree/mesh presets are family-specific and should not be treated as universal default candidates

Family quality winners:
- `cube`: `igraph_kk_default`
- `cycle`: `igraph_kk_default`
- `cylinder`: `igraph_kk_default`
- `kary.tree`: `igraph_kk_default`
- `mesh`: `igraph_kk_default`
- `path`: `grip_f32`
- `sierpinski.carpet`: `igraph_kk_default`
- `sierpinski.tetrahedron`: `grip_default_adaptive`
- `sierpinski.triangle`: `grip_default_adaptive`
- `torus`: `grip_torus_preset_reference`

Comparable subset (`grip_f32` vs `igraph_kk_default`):
- graphs where both completed: `35`
- mean RMSE: grip_f32 `0.1923`, igraph KK `0.1956`
- mean elapsed sec: grip_f32 `0.113`, igraph KK `3.906`

Universal comparable subset (`grip_*` universal candidates plus `igraph_kk_default`):
- comparable graphs: `35`
- universal quality winner: `igraph_kk_default` (mean quality score `0.3100`)
- universal value winner: `igraph_kk_default` (mean value score `0.2756`)
- best GRIP universal quality candidate: `grip_f128` (mean quality score `0.4899`)
- best GRIP universal value candidate: `grip_f96` (mean value score `0.4996`)

Best GRIP setting by family:
- `cube`: `grip_exact_f96`
- `cycle`: `grip_default_adaptive`
- `cylinder`: `grip_exact_f96`
- `kary.tree`: `grip_tree_preset_reference`
- `mesh`: `grip_mesh_preset_reference`
- `path`: `grip_f32`
- `sierpinski.carpet`: `grip_default_adaptive`
- `sierpinski.tetrahedron`: `grip_default_adaptive`
- `sierpinski.triangle`: `grip_default_adaptive`
- `torus`: `grip_torus_preset_reference`

Primary outputs:
- Stage-1 summary: `dev/manual/tmp/grip-kk-cross-family-full-2026-03-28/stage1-summary.md`
- Stage-2 summary: `dev/manual/tmp/grip-kk-cross-family-full-2026-03-28/stage2-summary.md`
- Stage-2 candidate CSV: `dev/manual/tmp/grip-kk-cross-family-full-2026-03-28/stage2-candidate-summary.csv`
- Stage-2 family CSV: `dev/manual/tmp/grip-kk-cross-family-full-2026-03-28/stage2-family-summary.csv`
- Universal common-subset summary: `dev/manual/tmp/grip-kk-cross-family-full-2026-03-28/stage2-universal-common-summary.md`
- Best GRIP-by-family summary: `dev/manual/tmp/grip-kk-cross-family-full-2026-03-28/stage2-best-grip-family-summary.md`
