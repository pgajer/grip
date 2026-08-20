# Carpet Layout PNGs

This report renders the seed-1 layouts corresponding to the timing comparison for
`grip.layout()`, `igraph::layout_with_fr()`, and `igraph::layout_with_kk()` on
Sierpinski carpet levels 4 and 5.

- contact sheet: `dev/design/sierpinski-carpet-level4-5-layouts-seed1-grid.png`

| Level | Method | Status | Elapsed sec | PNG |
| ---: | --- | --- | ---: | --- |
| 4 | grip.layout() | ok | 3.255 | `output/gkk_lgkk_paper/tmp/carpet-layout-timing-level4-5/png/sierpinski-carpet-level-4-grip-layout-seed-1.png` |
| 4 | igraph::layout_with_fr() | ok | 0.134 | `output/gkk_lgkk_paper/tmp/carpet-layout-timing-level4-5/png/sierpinski-carpet-level-4-igraph-fr-seed-1.png` |
| 4 | igraph::layout_with_kk() | ok | 23.546 | `output/gkk_lgkk_paper/tmp/carpet-layout-timing-level4-5/png/sierpinski-carpet-level-4-igraph-kk-seed-1.png` |
| 5 | grip.layout() | ok | 22.542 | `output/gkk_lgkk_paper/tmp/carpet-layout-timing-level4-5/png/sierpinski-carpet-level-5-grip-layout-seed-1.png` |
| 5 | igraph::layout_with_fr() | ok | 1.309 | `output/gkk_lgkk_paper/tmp/carpet-layout-timing-level4-5/png/sierpinski-carpet-level-5-igraph-fr-seed-1.png` |
| 5 | igraph::layout_with_kk() | timeout | NA | `output/gkk_lgkk_paper/tmp/carpet-layout-timing-level4-5/png/sierpinski-carpet-level-5-igraph-kk-seed-1.png` |
