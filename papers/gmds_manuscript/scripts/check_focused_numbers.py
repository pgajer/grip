#!/usr/bin/env python3
"""Recompute manuscript numerical claims from the frozen study tables.

Writes both LaTeX macros and a machine-readable claim ledger. Aggregate
validation summaries are explicitly distinguished from row-level calculations.
This does not rerun the original geodesic solvers or embedding optimizers.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

P = Path(__file__).resolve().parents[1]
E = P / "evidence"
ledger = {}


def read(name):
    return pd.read_csv(E / name)


def record(key, value, tex, source, rule):
    ledger[key] = dict(value=value, latex=tex, source=source, selection=rule)


def integer(key, value, source, rule):
    value = int(value)
    record(key, value, f"{value:,}".replace(",", r"{,}"), source, rule)


def decimal(key, value, places, source, rule):
    value = float(value)
    record(key, value, f"{value:.{places}f}", source, rule)


def scientific(key, value, source, rule):
    value = float(value)
    coefficient, exponent = f"{value:.2e}".split("e")
    record(key, value, coefficient + r"\times10^{" + str(int(exponent)) + "}", source, rule)


asrc = "study-initializers/summary/scores.csv.gz"
bsrc = "study-radius/summary/scores.csv.gz"
a, b = read(asrc), read(bsrc)
main = b[b.n.eq(240)]
v = a.pivot(index=["replicate", "k"], columns="method", values="path_rel")
c = a.pivot(index=["replicate", "k"], columns="method", values="stress1")
for init in ["Classical MDS", "Stress MDS"]:
    assert (v[init + " + edge-KK"] < v[init]).all()
    assert (c[init + " + edge-KK"] > c[init]).all()
integer("BoundedGraphs", len(v), asrc, "Distinct replicate/k; both initializers improve path and worsen chord on every row")
integer("StressWins", (v["Stress MDS + edge-KK"] < v["Classical MDS + edge-KK"]).sum(), asrc, "Paired final retained-path error on all graphs")
sel = a[a.selected]
for key, method in [("ClassicalBefore", "Classical MDS"), ("ClassicalAfter", "Classical MDS + edge-KK"),
                    ("StressBefore", "Stress MDS"), ("StressAfter", "Stress MDS + edge-KK")]:
    decimal(key, 100*sel[sel.method.eq(method)].path_rel.median(), 3, asrc, f"Selected graphs; median percent path error; {method}")
q = sel[sel.method.eq("Stress MDS + edge-KK")].procrustes*100
decimal("CoordLow", q.min(), 2, asrc, "Selected stress + edge-KK; minimum coordinate percent error")
decimal("CoordHigh", q.max(), 2, asrc, "Selected stress + edge-KK; maximum coordinate percent error")
assert sel.groupby("replicate").k.first().tolist() == [71, 70, 67, 72, 73]

gsrc = "study-radius/summary/graphs.csv.gz"
g = read(gsrc)
integer("RadiusGraphs", len(g), gsrc, "All primary n=240 and n=480 graph rows")
integer("RadiusCandidates", len(b), bsrc, "All primary configurations, before and after refinement")
integer("RepairedGraphs", (g.bridges > 0).sum(), gsrc, "Graphs with at least one added MST edge")
integer("InputMatrices", b['case'].nunique(), bsrc, "Distinct surface/sample/radius input cases")
ssrc = "study-radius/summary/starts.csv.gz"
s = read(ssrc)
integer("RadiusStarts", len(s), ssrc, "All primary MDS starts")
integer("CappedStarts", (s.iterations >= 1000).sum(), ssrc, "Starts reaching declared 1000-iteration cap")
assert len(s) == 3*len(g)
vsrc = "study-radius/summary/validation.json"
validation = json.loads((E / vsrc).read_text())
integer("ValidationPairs", validation["geodesic_validation_pairs"], vsrc, "Recorded independent high-radius validation; not a new solver run")
scientific("GeodesicDiscrepancy", validation["max_sampled_geodesic_discrepancy"], vsrc, "Recorded sampled maximum; not an all-pairs bound")
csrc = "study-radius/summary/coordinate-validation.csv.gz"
cv = read(csrc)
integer("ValidatedCoordinates", cv.candidates.iloc[0], csrc, "Recorded primary plus control coordinate reconstructions")
scientific("CoordinateDiscrepancy", cv.max_scaled_discrepancy.iloc[0], csrc, "Recorded independently reconstructed diagnostic maximum")

q = main[main.surface.eq("saddle") & main.sampling.eq("disk") & main.radius.eq(64) & main.k.eq(32) & main.method.eq("stress_fixed_primary")]
for reg, prefix in [("geodesic", "GeodesicSaddle"), ("ambient", "AmbientSaddle")]:
    for metric, suffix in [("path_reference", "Path"), ("procrustes", "Coord")]:
        vals = q[q.regime.eq(reg)][metric]
        assert len(vals) == 3
        decimal(prefix+suffix, 100*vals.median(), 2, bsrc, f"n240 r64 k32 base-disk saddle, {reg}, stress_fixed_primary; median percent {metric}")

q = main[main.method.eq("full_classical_fixed_primary")]
ranges = q.groupby("case").path_reference.agg(["min", "max"])
integer("FullInitializerCases", ((ranges['max']-ranges['min']) > 1e-8).sum(), bsrc, "Distinct n240 cases; full classical initializer, fixed-scale continuation; variation across k exceeds 1e-8")
assert len(ranges) == 84
spread = 100*(ranges['max']-ranges['min'])
for suffix, value in [("Min", spread.min()), ("Median", spread.median()), ("Max", spread.max())]:
    decimal("NeighborhoodSpread"+suffix, value, 2, bsrc,
            "Within-case max minus min across all seven k, 84 n240 full-geodesic classical fixed-scale cases; percentage points")
for key, n, k in [("NestedBase", 240, 8), ("NestedSameK", 480, 8), ("NestedSameFraction", 480, 16)]:
    vals = b[b.surface.eq("paraboloid") & b.sampling.eq("disk") & b.replicate.eq(1) & b.radius.eq(64) & b.n.eq(n) & b.k.eq(k) & b.regime.eq("geodesic") & b.method.eq("stress_fixed_primary")].path_reference
    assert len(vals) == 1
    decimal(key, 100*vals.iloc[0], 2, bsrc, f"Nested rep1 disk paraboloid, r64, geodesic, fixed-scale stress continuation; n={n}, k={k}; percent path/reference error")

osrc = "study-radius/summary/optimizer-sensitivity.csv.gz"
o = read(osrc)
for init in ["classical", "stress"]:
    for control in ["perturbed", "random"]:
        maxima = []
        for policy in ["", "fixed_"]:
            z = o[o.method.eq(f"{init}_{policy}{control}")].merge(
                b[b.method.eq(f"{init}_{policy}primary")], on=["case", "regime", "k"], suffixes=("_control", "_base"), validate="one_to_one")
            maxima.append(float(abs(z.sigma3_sigma2_control-z.sigma3_sigma2_base).max()))
        if control == "perturbed":
            decimal(init.capitalize()+"Perturbation", max(maxima), 5, osrc+"; "+bsrc,
                    f"Maximum absolute change of s3/s2 from matched primary fit, both scales, {init} starts")
        else:
            decimal(init.capitalize()+"RandomFixed", maxima[1], 2, osrc+"; "+bsrc,
                    f"Maximum absolute change of s3/s2 from matched fixed-scale primary fit, {init} starts")

profiled = b[b.method.str.endswith(("_primary", "_uniform")) & ~b.method.str.contains("_fixed_")]
fixed = b[b.method.str.contains("_fixed_")]
integer("ProfiledCount", len(profiled), bsrc, "All primary profiled refinements, both sizes, all initializers and schedules")
for suffix, threshold in [("Half", .5), ("Tenth", .1), ("Thousandth", .001)]:
    integer("ProfiledBelow"+suffix, (profiled.edge_scale < threshold).sum(), bsrc,
            f"All primary profiled refinements; edge calibration strictly below {threshold}")
decimal("ProfiledMedian", profiled.edge_scale.median(), 3, bsrc, "Median calibration across all primary profiled refinements")
scientific("ProfiledMinimum", profiled.edge_scale.min(), bsrc, "All primary refinements, all initializers/schedules/sample sizes; minimum edge calibration")
decimal("FixedMinimum", fixed.edge_scale.min(), 5, bsrc, "All primary fixed-scale refinements, all initializers/schedules/sample sizes; minimum edge calibration")
decimal("FixedMaximum", fixed.edge_scale.max(), 5, bsrc, "All primary fixed-scale refinements, all initializers/schedules/sample sizes; maximum edge calibration")
v = main.pivot(index=["case", "regime", "k"], columns="method", values="path_rel")
integer("PrimaryGraphs", len(v), bsrc, "Distinct n240 primary graph settings")
for prefix, method in [("Profiled", "stress_primary"), ("Fixed", "stress_fixed_primary")]:
    d = v[method]-v.stress
    for suffix, count in [("Improved", (d < -1e-8).sum()), ("Worsened", (d > 1e-8).sum()), ("Tied", (abs(d) <= 1e-8).sum())]:
        integer(prefix+suffix, count, bsrc, f"Paired n240 stress MDS versus {method}; absolute tie tolerance 1e-8")
assert [ledger[k]['value'] for k in ['ProfiledImproved','ProfiledWorsened','ProfiledTied','FixedImproved','FixedWorsened','FixedTied']] == [1016,76,84,1075,17,84]

# Paired continuation contrasts, not unpaired differences of medians.
continuation = []
for label, policy in [("Profiled", ""), ("Fixed", "fixed_")]:
    primary = main[main.method.eq("stress_"+policy+"primary")].set_index(["case", "regime", "k"])
    uniform = main[main.method.eq("stress_"+policy+"uniform")].set_index(["case", "regime", "k"])
    assert len(primary) == len(uniform) == 1176 and primary.index.equals(uniform.index)
    for metric in ["path_rel", "procrustes"]:
        d = primary[metric]-uniform[metric]
        continuation.append(dict(scale=label, metric=metric, pairs=len(d),
            lower=int((d < -1e-8).sum()), higher=int((d > 1e-8).sum()), tied=int((abs(d)<=1e-8).sum()),
            q25_pp=100*d.quantile(.25), median_pp=100*d.median(), q75_pp=100*d.quantile(.75),
            minimum_pp=100*d.min(), maximum_pp=100*d.max()))
fixed_path = next(row for row in continuation if row["scale"] == "Fixed" and row["metric"] == "path_rel")
assert [fixed_path[k] for k in ["pairs", "lower", "higher", "tied"]] == [1176, 389, 541, 246]
pd.DataFrame(continuation).to_csv(E/"continuation-comparison.csv", index=False, float_format="%.12g")
components = pd.read_csv(E/"coordinate-components/component-medians.csv")
for surface, regime, prefix in [("saddle", "ambient", "AmbientSaddle"),
                                ("paraboloid", "ambient", "AmbientParaboloid"),
                                ("paraboloid", "geodesic", "GeodesicParaboloid")]:
    row = components[(components.surface==surface)&(components.regime==regime)].iloc[0]
    for metric, suffix in [("horizontal_error", "Horizontal"), ("vertical_error", "Vertical"), ("height_only_error", "HeightOnly")]:
        decimal(prefix+suffix, 100*row[metric], 2, "coordinate-components/coordinates.csv.gz",
            f"Median over three n240 r64 k32 base-disk {surface} {regime} fixed-scale stress continuation fits; one global similarity alignment; percent {metric}")
paths = pd.read_csv(E/"retained-route-validation/independent-check.csv")
assert len(paths) == 24
frozen = b[b.n.eq(240)&b.radius.eq(64)&b.k.eq(32)&b.sampling.eq("disk")&b.method.isin(["stress", "stress_fixed_primary"])]
matched = paths.merge(frozen, on=["case", "regime", "method"], suffixes=("_check", "_saved"), validate="one_to_one")
assert len(matched)==24
discrepancy = max(abs(matched.path_rel_check-matched.path_rel_saved).max(),
                  abs(matched.path_reference_check-matched.path_reference_saved).max())
assert discrepancy < 8e-15
scientific("RetainedRouteDiscrepancy", discrepancy, "retained-route-validation/independent-check.csv; "+bsrc,
           "Maximum absolute disagreement of all-route recalculations for 24 configurations with frozen path and path-to-surface scores; caches required only to rerun optional R check")
(P / "tables/focused/numbers.tex").write_text("% Generated by check_focused_numbers.py; do not edit.\n" + "".join(r"\newcommand{\num"+k+"}{"+v["latex"]+"}\n" for k,v in ledger.items()))
(E / "numerical-claims.json").write_text(json.dumps(ledger, indent=2)+"\n")
print(f"Recomputed {len(ledger)} numerical claims; all paired and selection checks passed.")
