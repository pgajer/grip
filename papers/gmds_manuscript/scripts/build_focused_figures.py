#!/usr/bin/env python3
"""Generate focused-paper figures, tables, and numeric macros from frozen data."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

P = Path(__file__).resolve().parents[1]
F = P / "figures/focused"
T = P / "tables/focused"
F.mkdir(parents=True, exist_ok=True)
T.mkdir(parents=True, exist_ok=True)
a = pd.read_csv(P / "evidence/study-initializers/summary/scores.csv.gz")
b = pd.read_csv(P / "evidence/study-radius/summary/scores.csv.gz")
main = b[b.n.eq(240)]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42})
blue, orange, green, purple = "#24679b", "#be681e", "#227a68", "#89508f"


def save(fig, name):
    fig.savefig(F / (name + ".pdf"), bbox_inches="tight",
                metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(F / (name + ".png"), dpi=160, bbox_inches="tight")
    plt.close(fig)


fig, axes = plt.subplots(1, 2, figsize=(8.2, 2.7), layout="constrained")
for ax, theta, title in zip(axes, [2*np.pi/3, np.pi/6], ["Wide angle", "Narrow angle"]):
    z = np.array([[1, 0], [0, 0], [np.cos(theta), np.sin(theta)]])
    ax.plot(*z.T, "o-", color=blue, lw=3)
    ax.plot(*z[[0, 2]].T, "--", color=orange, lw=2)
    for i, (x, y) in enumerate(z):
        ax.text(x+.05, y+.05, str(i+1), fontsize=12)
    ax.set_aspect("equal")
    ax.set(xlim=(-.7, 1.3), ylim=(-.25, 1.2), title=title)
    ax.text(-.65, -.16, f"Path: 2     Chord: {np.linalg.norm(z[0]-z[2]):.2f}")
    ax.axis("off")
save(fig, "paths")

fig, axes = plt.subplots(2, 2, figsize=(8.4, 6.0), layout="constrained")
for row, surface in enumerate(["paraboloid", "saddle"]):
    for col, metric in enumerate(["sigma2_sigma1", "sigma3_sigma2"]):
        ax = axes[row, col]
        for method, color in [("classical", blue), ("stress", orange)]:
            q = main[(main.surface==surface)&(main.sampling=="disk")&
                     (main.regime=="geodesic")&(main.k==239)&(main.method==method)]
            v = q.groupby("radius")[metric].agg(["median", "min", "max"])
            ax.plot(v.index, v["median"], "o-", color=color, label=method.capitalize()+" MDS")
            ax.fill_between(v.index, v["min"], v["max"], color=color, alpha=.14)
        ax.set(xscale="log", ylim=(0, 1.02), xlabel="Disk radius r",
               ylabel=["Second / first singular value", "Third / second singular value"][col],
               title=surface.capitalize())
        ax.grid(alpha=.2)
axes[0, 0].legend(fontsize=9)
save(fig, "shape-ratios")

fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.4), layout="constrained")
methods = ["Classical MDS", "Classical MDS + edge-KK", "Stress MDS", "Stress MDS + edge-KK"]
selected = a[a.selected]
for ax, metric, title in zip(axes, ["path_rel", "stress1", "procrustes"],
                            ["Retained-path error (%)", "Chord Stress-1 (%)", "Coordinate error (%)"]):
    for start, color in [(0, blue), (2, orange)]:
        subset = selected[selected.method.isin(methods[start:start+2])]
        v = subset.pivot(index="replicate", columns="method", values=metric)[methods[start:start+2]]*100
        for row in v.to_numpy():
            ax.plot([start, start+1], row, "o-", color=color, alpha=.4, ms=4)
        ax.plot([start, start+1], v.median(), "s-", color=color, lw=2.4, ms=7)
    ax.set_xticks(range(4), ["Classical", "+ edge-KK", "Stress", "+ edge-KK"], rotation=30, ha="right")
    ax.set_ylabel(title)
    ax.grid(axis="y", alpha=.2)
save(fig, "initializer-comparison")

fig, axes = plt.subplots(2, 2, figsize=(8.4, 6.2), layout="constrained")
for col, surface in enumerate(["paraboloid", "saddle"]):
    for row, metric in enumerate(["graph_reference", "path_reference"]):
        ax = axes[row, col]
        for regime, color in [("geodesic", green), ("ambient", purple)]:
            q = main[(main.surface==surface)&(main.sampling=="disk")&(main.radius==64)&
                     (main.regime==regime)&(main.method=="stress_fixed_primary")]
            v = q.groupby("k")[metric].agg(["median", "min", "max"])*100
            # Exact complete-graph reference matches can lie below a log display floor.
            low = 1e-7
            ax.plot(v.index, np.maximum(v["median"], low), "o-", color=color, label=regime.capitalize())
            ax.fill_between(v.index, np.maximum(v["min"], low), np.maximum(v["max"], low), color=color, alpha=.14)
        ax.set(xscale="log", yscale="log", xlabel="Neighbors k", title=surface.capitalize(),
               ylabel=["Graph / surface distance error (%)", "Embedded path / surface error (%)"][row])
        ax.grid(alpha=.2)
axes[0,0].legend(fontsize=9)
save(fig, "graph-sensitivity")

fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5), layout="constrained")
for method, color, label in [("stress_primary", green, "Profiled scale"),
                              ("stress_fixed_primary", purple, "Fixed scale")]:
    q = main[main.method.eq(method)]
    v = q.groupby("radius").edge_scale.agg(["median", "min", "max"])
    axes[0].plot(v.index, v["median"], "o-", color=color, label=label)
    axes[0].fill_between(v.index, v["min"], v["max"], color=color, alpha=.14)
    vv = main.pivot(index=["case", "regime", "k"], columns="method", values="path_rel")
    mask = (vv.stress>1e-8)&(vv[method]>1e-8)
    axes[1].scatter(vv.loc[mask, "stress"]*100, vv.loc[mask, method]*100, s=9,
                    alpha=.3, color=color, label=label, rasterized=True)
axes[0].set(xscale="log", yscale="log", xlabel="Disk radius r", ylabel="Edge calibration factor b")
axes[0].axhline(1, color="#555", ls="--", lw=1)
axes[0].legend(fontsize=9)
axes[1].plot([1e-6, 100], [1e-6, 100], "--", color="#555", lw=1)
axes[1].set(xscale="log", yscale="log", xlabel="Path error before edge-KK (%)",
            ylabel="Path error after edge-KK (%)")
for ax in axes: ax.grid(alpha=.15)
save(fig, "scale")

snap = pd.read_csv(P / "evidence/study-radius/summary/snapshot-coordinates.csv.gz")
branch_colors = np.array(["#24679b", "#be681e", "#227a68", "#89508f"])
fig = plt.figure(figsize=(10.2, 7.0), layout="constrained")
spatial_checks = []
for row, radius in enumerate([1, 64]):
    case = f"saddle-disk-rep1-n240-r{radius}"
    q = snap[snap.case.eq(case) & snap.k.eq(32)]
    truth = q[q.regime.eq("geodesic") & q.method.eq("original")].sort_values("vertex")
    x = truth[["x", "y", "z"]].to_numpy()
    assert len(x) == 240
    center = x.mean(axis=0)
    xc = (x-center)/radius**2
    branch = np.where(abs(x[:,0]) >= abs(x[:,1]), np.where(x[:,0] >= 0, 0, 1),
                      np.where(x[:,1] >= 0, 2, 3))
    panels = [(xc, "Original surface")]
    for regime in ["geodesic", "ambient"]:
        zq = q[q.regime.eq(regime) & q.method.eq("stress_fixed_primary")].sort_values("vertex")
        assert np.array_equal(truth.vertex, zq.vertex)
        z = zq[["x", "y", "z"]].to_numpy()
        z -= z.mean(axis=0)
        u, _, vt = np.linalg.svd(z.T @ xc)
        rot = u @ vt
        scale = zq.edge_scale.iloc[0]
        assert np.allclose(rot.T@rot, np.eye(3), atol=1e-12)
        display = z @ rot / scale / radius**2
        assert np.allclose(np.linalg.norm(display, axis=1), np.linalg.norm(z, axis=1)/scale/radius**2)
        panels.append((display, regime.capitalize()+" graph + edge-KK"))
        spatial_checks.append(dict(case=case, regime=regime, vertices=len(z), edge_scale=scale,
                                   rotation=rot.tolist(), uniform_divisor=scale*radius**2))
    rho, theta = np.meshgrid(np.linspace(0, radius, 12), np.linspace(0, 2*np.pi, 49))
    xx, yy = rho*np.cos(theta), rho*np.sin(theta)
    mesh = np.stack([xx, yy, xx*xx-yy*yy], axis=-1)
    mesh = (mesh-center)/radius**2
    limit = 1.08*max(max(abs(z).max() for z, _ in panels), abs(mesh).max())
    for col, (z, title) in enumerate(panels):
        ax = fig.add_subplot(2, 3, 1+row*3+col, projection="3d")
        ax.set_proj_type("ortho")
        if col == 0:
            ax.plot_wireframe(mesh[:,:,0], mesh[:,:,1], mesh[:,:,2], rstride=4, cstride=1,
                              color="#777777", alpha=.35, linewidth=.45)
            ax.plot(*mesh[:,-1,:].T, color="#333333", lw=.8, alpha=.65)
        ax.scatter(*z.T, c=branch_colors[branch], s=9, alpha=.75, depthshade=False)
        ax.set(xlim=(-limit, limit), ylim=(-limit, limit), zlim=(-limit, limit),
               title=f"{title}\nr = {radius}", xlabel="X", ylabel="Y", zlabel="Z")
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=22, azim=-58)
        ax.locator_params(nbins=3)
        ax.tick_params(labelsize=7, pad=0)
        ax.xaxis.labelpad = ax.yaxis.labelpad = ax.zaxis.labelpad = -3
from matplotlib.lines import Line2D
fig.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=c, label=l)
                    for c,l in zip(branch_colors, ["+x branch", "−x branch", "+y branch", "−y branch"])],
           loc="outside lower center", ncol=4, frameon=False)
save(fig, "saddle-spatial")
(P/"evidence/spatial-display-checks.json").write_text(json.dumps(spatial_checks, indent=2)+"\n")

med = selected.groupby("method")[["path_rel", "edge_rel", "stress1", "procrustes"]].median()*100
rows = [r"\begin{tabular}{lrrrr}", r"\toprule", r"Method & Retained path (\%) & Edge (\%) & Chord (\%) & Coordinate (\%)\\", r"\midrule"]
for method in ["Classical MDS", "Stress MDS", "Classical MDS + edge-KK", "Stress MDS + edge-KK"]:
    vals = med.loc[method]
    rows.append(method + " & " + " & ".join(f"{v:.3f}" if i<3 else f"{v:.2f}" for i,v in enumerate(vals)) + r"\\")
rows += [r"\bottomrule", r"\end{tabular}"]
(T/"initializers.tex").write_text("\n".join(rows)+"\n")
q = main[(main.sampling=="disk")&(main.radius==64)&(main.k==32)&(main.method=="stress_fixed_primary")]
med2=q.groupby(["surface", "regime"])[["sigma2_sigma1", "sigma3_sigma2", "procrustes", "path_reference"]].median()
rows=[r"\begin{tabular}{llrrrr}",r"\toprule",r"Surface & Graph & $s_2/s_1$ & $s_3/s_2$ & Coordinate (\%) & \shortstack{Retained path/\\surface (\%)}\\",r"\midrule"]
for surface in ["paraboloid", "saddle"]:
    for regime in ["geodesic", "ambient"]:
        v=med2.loc[surface,regime]
        rows.append(f"{surface.capitalize()} & {regime} & {v.iloc[0]:.3f} & {v.iloc[1]:.3f} & {100*v.iloc[2]:.2f} & {100*v.iloc[3]:.2f}"+r"\\")
rows += [r"\bottomrule",r"\end{tabular}"]
(T/"radius.tex").write_text("\n".join(rows)+"\n")

components = pd.read_csv(P/"evidence/coordinate-components/component-medians.csv")
rows = [r"\begin{tabular}{llrrr}", r"\toprule",
        r"Surface & Graph & Global (\%) & Horizontal (\%) & Vertical (\%)\\", r"\midrule"]
for v in components.itertuples():
    rows.append(f"{v.surface.capitalize()} & {v.regime} & {100*v.coordinate_error:.2f} & {100*v.horizontal_error:.2f} & {100*v.vertical_error:.3f}"+r"\\")
rows += [r"\bottomrule", r"\end{tabular}"]
(T/"components.tex").write_text("\n".join(rows)+"\n")

continuation = pd.read_csv(P/"evidence/continuation-comparison.csv")
rows = [r"\begin{tabular}{llrrr}", r"\toprule",
        r"Scale & Error & \shortstack{Lower/higher/\\tied} & \shortstack{Median [quartiles]\\(pp)} & \shortstack{Minimum, maximum\\(pp)}\\", r"\midrule"]
for v in continuation.itertuples():
    label = "Retained path" if v.metric == "path_rel" else "Coordinate"
    rows.append(f"{v.scale} & {label} & {v.lower}/{v.higher}/{v.tied} & {v.median_pp:.4f} [{v.q25_pp:.4f}, {v.q75_pp:.4f}] & {v.minimum_pp:.3f}, {v.maximum_pp:.3f}"+r"\\")
rows += [r"\bottomrule", r"\end{tabular}"]
(T/"continuation.tex").write_text("\n".join(rows)+"\n")

v=main.pivot(index=["case", "regime", "k"],columns="method",values="path_rel")
counts={}
for method in ["stress_primary","stress_fixed_primary"]:
    delta=v[method]-v.stress
    counts[method]=dict(improved=int((delta < -1e-8).sum()),worsened=int((delta > 1e-8).sum()),ties=int((abs(delta)<=1e-8).sum()))
assert counts["stress_primary"]==dict(improved=1016,worsened=76,ties=84)
assert counts["stress_fixed_primary"]==dict(improved=1075,worsened=17,ties=84)
counts["initializers"]={}
for init in ["Classical MDS", "Stress MDS"]:
    v=a.pivot(index=["replicate","k"],columns="method",values="path_rel")
    assert (v[init+" + edge-KK"] < v[init]).all()
    counts["initializers"][init]=int((v[init+" + edge-KK"] < v[init]).sum())
(P/"evidence/focused-number-checks.json").write_text(json.dumps(counts,indent=2)+"\n")
print("Generated six figures, four tables, paired-count checks, and spatial transform records.")
