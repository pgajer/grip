# GRIP / KK Round-2 Summary

## Scope

Round 2 followed directly from the round-1 trace finding that the level-4 carpet reaches its best
global symmetry before the final full-graph phase completes. The second round therefore tested two
questions:

1. How much of the tradeoff is controlled by `final_rounds` alone?
2. Does exact coarse repulsion materially improve the result enough to justify its cost?

## Design

- graph: Sierpinski carpet level 4
- size: `4096` vertices, `6424` edges
- seeds per configuration: `6`
- Stage A:
  sampled coarse repulsion with `final_rounds = 0, 32, 64, 96, 128, 160, 192, 256, 384`
- Stage B:
  exact coarse repulsion with `final_rounds = 0, 64, 96, 128, 160, 192`

Metrics:
- Procrustes RMSE to canonical carpet
- edge-length CV
- sampled non-edge separation ratio
- sampled stress
- elapsed seconds
- edge axis-deviation proxy after Procrustes alignment

## Main Result

The best overall balance on the level-4 carpet came from the sampled-coarse configuration with
`final_rounds = 32`.

Top configurations by quality-rank sum:

| Config | Coarse mode | `final_rounds` | Mean RMSE | Mean edge CV | Mean axis dev | Mean sep ratio | Mean sec |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `sampled_f32` | sampled | `32` | `0.0305` | `0.1469` | `0.1652` | `1.1862` | `0.927` |
| `sampled_f64` | sampled | `64` | `0.0341` | `0.1536` | `0.1688` | `1.1747` | `1.273` |
| `exact_f64` | exact | `64` | `0.0368` | `0.1570` | `0.1709` | `1.1768` | `18.146` |
| `sampled_f96` | sampled | `96` | `0.0365` | `0.1574` | `0.1713` | `1.1759` | `1.601` |

Reference rows:

- current default `sampled_f384`:
  RMSE `0.0446`, edge CV `0.1623`, axis deviation `0.1748`, sep ratio `1.2104`, mean runtime `3.722s`
- no final stage `sampled_f0`:
  RMSE `0.0220`, edge CV `0.1708`, axis deviation `0.1926`, sep ratio `1.1798`, mean runtime `0.400s`

## Interpretation

The second round sharpens the round-1 conclusion:

- Cutting the finest-level budget is the high-leverage knob.
- `final_rounds = 0` gives the best symmetry but is visibly too wavy and too non-rectilinear.
- A small but nonzero finest-level budget fixes a lot of that waviness without giving away all of
  the early global symmetry.
- On this graph, `final_rounds = 32` is the best balance tested.
- Exact coarse repulsion was not worth its cost here.
  At matched `final_rounds`, it was consistently much slower and did not beat the sampled-coarse
  variants on the quality balance that matters for this carpet.

So the main practical conclusion from round 2 is:

- if the target is the level-4 carpet, reducing `final_rounds` is much more promising than making
  coarse repulsion exact
- the most promising next GRIP modification is a new or adaptive finest-level schedule, not a more
  expensive coarse-level repulsion mode

## Artifacts

Tracked planning document:
- `dev/design/grip_kk_round2_action_plan_2026-03-28.md`

Generated experiment outputs:
- `output/gkk_lgkk_paper/tmp/carpet-level4-finalstage-round2/carpet-level4-finalstage-round2-summary.md`
- `output/gkk_lgkk_paper/tmp/carpet-level4-finalstage-round2/carpet-level4-finalstage-round2-summary.csv`
- `output/gkk_lgkk_paper/tmp/carpet-level4-finalstage-round2/carpet-level4-finalstage-round2-contact-sheet.png`
- `output/gkk_lgkk_paper/tmp/carpet-level4-finalstage-round2/carpet-level4-finalstage-round2-final-rounds-lines.png`
- `output/gkk_lgkk_paper/tmp/carpet-level4-finalstage-round2/carpet-level4-finalstage-round2-tradeoff-scatter.png`
