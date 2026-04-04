#!/usr/bin/env python3
"""Figure 7: 2D vs. 3D Dimensional Comparison.

Side-by-side grouped bars showing Full-GKK Relative RMSE for a
selected set of methods in 2D (left bars) vs. 3D (right bars),
for each benchmark family.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[5]
SUMMARY_PATH = (
    REPO_ROOT
    / "output"
    / "gkk_lgkk_paper"
    / "benchmarks"
    / "weighted-grip-phase5-smoke-check"
    / "data"
    / "summary_metrics.csv"
)

# --- Load data ---
df = pd.read_csv(SUMMARY_PATH)
df = df[df["method"] != "start"].copy()

# --- Focus on key methods for clarity ---
methods_show = ["wgrip_polish_lgkk", "lgkk", "gkk"]
method_labels = {
    "wgrip_polish_lgkk": "W-GRIP+pLGKK",
    "lgkk":              "KK→LGKK",
    "gkk":               "KK→GKK",
}

family_order = [
    "mesh", "torus", "sierpinski_carpet", "cube_channel_network", "irregular_torus"
]
family_labels = {
    "mesh":                 "Mesh",
    "torus":                "Torus",
    "sierpinski_carpet":    "Carpet",
    "cube_channel_network": "Cube ch.",
    "irregular_torus":      "Irr. torus",
}

colors_2d = {"wgrip_polish_lgkk": "#1abc9c", "lgkk": "#8e44ad", "gkk": "#2980b9"}
colors_3d = {"wgrip_polish_lgkk": "#0e6655", "lgkk": "#6c3483", "gkk": "#1a5276"}

# --- Plot ---
mpl.rcParams.update({
    "font.family": "serif",
    "font.size": 9,
    "axes.linewidth": 0.5,
})

fig, ax = plt.subplots(figsize=(7.0, 3.5))

n_families = len(family_order)
n_methods = len(methods_show)
n_bars = n_methods * 2  # 2D + 3D per method
x = np.arange(n_families)
bar_width = 0.11
gap = 0.03  # small gap between 2D and 3D within same method

for i, m in enumerate(methods_show):
    vals_2d = []
    vals_3d = []
    for f in family_order:
        r2 = df[(df["method"] == m) & (df["family_id"] == f) & (df["dim"] == 2)]
        r3 = df[(df["method"] == m) & (df["family_id"] == f) & (df["dim"] == 3)]
        vals_2d.append(r2["gkk_rel_rmse_mean"].values[0] if len(r2) > 0 else 0)
        vals_3d.append(r3["gkk_rel_rmse_mean"].values[0] if len(r3) > 0 else 0)

    base_offset = (i - (n_methods - 1) / 2) * (2 * bar_width + gap)
    ax.bar(
        x + base_offset - bar_width / 2 - gap / 4,
        vals_2d, bar_width * 0.92,
        label=f"{method_labels[m]} 2D",
        color=colors_2d[m], edgecolor="white", linewidth=0.3,
        hatch="//", alpha=0.85,
    )
    ax.bar(
        x + base_offset + bar_width / 2 + gap / 4,
        vals_3d, bar_width * 0.92,
        label=f"{method_labels[m]} 3D",
        color=colors_3d[m], edgecolor="white", linewidth=0.3,
    )

ax.set_xticks(x)
ax.set_xticklabels([family_labels[f] for f in family_order])
ax.set_ylabel("Full-GKK Relative RMSE")
ax.set_title("2D vs. 3D Geodesic Fidelity Comparison", fontsize=10, fontweight="bold")
ax.legend(
    fontsize=6.5, ncol=3, loc="upper left",
    framealpha=0.9, edgecolor="gray", handlelength=1.4,
)
ax.set_ylim(0, None)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="y", alpha=0.3, linewidth=0.4)

fig.tight_layout()
fig.savefig("fig7_2d_vs_3d.pdf", bbox_inches="tight", dpi=300)
fig.savefig("fig7_2d_vs_3d.png", bbox_inches="tight", dpi=200)
print("Figure 7 saved.")
