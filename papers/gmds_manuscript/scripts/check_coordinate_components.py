#!/usr/bin/env python3
"""Reconstruct component errors after ONE similarity fit, and local rigidity ranks."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

P = Path(__file__).resolve().parents[1]
E = P/"evidence/coordinate-components"
manifest = json.loads((E/"source-manifest.json").read_text())
assert hashlib.sha256((E/"coordinates.csv.gz").read_bytes()).hexdigest() == manifest["export_sha256"]
data = pd.read_csv(E/"coordinates.csv.gz", float_precision="round_trip")
scores = pd.read_csv(P/"evidence/study-radius/summary/scores.csv.gz")
assert len(data) == manifest["rows"] == 7200
records, rigidity = [], []
for case, group in data.groupby("case"):
    truth = group[group.method.eq("original")].sort_values("vertex")
    x = truth[["x", "y", "z"]].to_numpy()
    xc = x-x.mean(axis=0)
    for (regime, method), fit in group[~group.method.eq("original")].groupby(["regime", "method"]):
        fit = fit.sort_values("vertex")
        assert np.array_equal(fit.vertex, truth.vertex)
        z = fit[["x", "y", "z"]].to_numpy()
        z -= z.mean(axis=0)
        u, s, vt = np.linalg.svd(z.T@xc)
        aligned = s.sum()/np.sum(z*z) * z@(u@vt)
        residual = aligned-xc
        global_error = np.linalg.norm(residual)/np.linalg.norm(xc)
        horizontal = np.linalg.norm(residual[:, :2])/np.linalg.norm(xc[:, :2])
        vertical = np.linalg.norm(residual[:, 2])/np.linalg.norm(xc[:, 2])
        fraction = np.sum(xc[:, :2]**2)/np.sum(xc**2)
        assert np.isclose(global_error**2, fraction*horizontal**2+(1-fraction)*vertical**2, atol=1e-14)
        saved = scores[(scores['case']==case)&(scores.regime==regime)&(scores.k==32)&(scores.method==method)]
        assert len(saved) == 1 and abs(saved.procrustes.iloc[0]-global_error) < 1e-10
        records.append(dict(case=case, surface=case.split('-')[0], regime=regime, method=method,
            coordinate_error=global_error, horizontal_error=horizontal, vertical_error=vertical,
            height_only_error=np.sqrt(fraction), saved_discrepancy=abs(saved.procrustes.iloc[0]-global_error)))
    if case.startswith("saddle"):
        n = len(x)
        delta = x[:, None, :]-x[None, :, :]
        distance = np.linalg.norm(delta, axis=2)
        edges = set()
        for i in range(n):
            order = np.lexsort((np.arange(n), distance[i]))
            edges.update(tuple(sorted((i, int(j)))) for j in order[order != i][:32])
        edges = sorted(edges)
        seen = {0}
        while True:
            expanded = seen | {j for i,j in edges if i in seen} | {i for i,j in edges if j in seen}
            if expanded == seen: break
            seen = expanded
        assert len(seen) == n  # No MST repair in these highlighted graphs.
        jacobian = np.zeros((len(edges), 3*n))
        for row, (i,j) in enumerate(edges):
            direction = delta[i,j]/distance[i,j]
            jacobian[row,3*i:3*i+3] = direction
            jacobian[row,3*j:3*j+3] = -direction
        s = np.linalg.svd(jacobian, compute_uv=False)
        tolerance = max(jacobian.shape)*np.finfo(float).eps*s[0]
        rank = int(np.sum(s > tolerance))
        assert rank == 3*n-6
        rigidity.append(dict(case=case, edges=len(edges), rank=rank, tolerance=tolerance,
            smallest_nongauge_singular_value=s[3*n-7], largest_gauge_singular_value=s[3*n-6]))
result = pd.DataFrame(records)
medians = result[result.method.eq("stress_fixed_primary")].groupby(["surface", "regime"])[
    ["coordinate_error", "horizontal_error", "vertical_error", "height_only_error"]].median().reset_index()
frames = {"component-errors.csv": result, "component-medians.csv": medians,
          "rigidity.csv": pd.DataFrame(rigidity)}
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--refresh', action='store_true', help='Explicitly replace the frozen diagnostic exports')
args = parser.parse_args()
for name, actual in frames.items():
    saved = pd.read_csv(E/name)
    assert list(actual.columns) == list(saved.columns) and actual.shape == saved.shape, name
    for col in actual:
        if pd.api.types.is_numeric_dtype(actual[col]):
            # Export rounding and BLAS roundoff differ across environments; rank and
            # scientific errors are checked, without requiring identical zero residuals.
            assert np.allclose(actual[col], saved[col], rtol=1e-9, atol=1e-10), (name, col)
        else:
            assert actual[col].equals(saved[col]), (name, col)
    out = E if args.refresh else P/'build/validation/coordinate-components'
    out.mkdir(parents=True, exist_ok=True)
    actual.to_csv(out/name, index=False, float_format="%.12g")
print("Verified 24 globally aligned coordinate diagnostics and three rank-714 ambient saddle frameworks.")
