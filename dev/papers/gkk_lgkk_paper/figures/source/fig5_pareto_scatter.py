#!/usr/bin/env python3
"""Figure 5: Runtime vs. Quality Pareto Scatter (3D).

Each point is a (method, family) pair. X = runtime (log scale),
Y = Full-GKK Relative RMSE. Methods are colored; families are
distinguished by marker shape.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# --- Load data ---
df = pd.read_csv(
    "/sessions/keen-magical-bell/mnt/grip/output/benchmarks/weighted-grip-phase5-smoke-check/data/summary_metrics.csv"
)
df3 = df[(df["dim"] == 3) & (df["method"] != "start")].copy()

# --- Method styling ---
method_labels = {
    "grip":              "GRIP",
    "wgrip":             "W-GRIP",
    "wgrip_core_lgkk":   "W-GRIP+cLGKK",
    "wgrip_polish_lgkk": "W-GRIP+pLGKK",
    "kk":                "KK",
    "gkk":               "KK→GKK",
    "lgkk":              "KK→LGKK",
}
colors = {
    "grip":              "#7f8c8d",
    "wgrip":             "#27ae60",
    "wgrip_core_lgkk":   "#2ecc71",
    "wgrip_polish_lgkk": "#1abc9c",
    "kk":                "#3498db",
    "gkk":               "#2980b9",
    "lgkk":              "#8e44ad",
}

# --- Family markers ---
family_markers = {
    "mesh":                 "o",
    "torus":                "s",
    "sierpinski_carpet":    "D",
    "cube_channel_network": "^",
    "irregular_torus":      "v",
}
family_labels = {
    "mesh":                 "Mesh",
    "torus":                "Torus",
    "sierpinski_carpet":    "Carpet",
    "cube_channel_network": "Cube ch.",
    "irregular_torus":      "Irr. torus",
}

# --- Plot ---
mpl.rcParams.update({
    "font.family": "serif",
    "font.size": 9,
    "axes.linewidth": 0.5,
})

fig, ax = plt.subplots(figsize=(6.5, 4.0))

# Plot each method × family
for _, row in df3.iterrows():
    m = row["method"]
    f = row["family_id"]
    if m not in method_labels:
        continue
    rt = row["runtime_sec_mean"]
    if pd.isna(rt) or rt <= 0:
        rt = 0.001  # floor for log scale
    rmse = row["gkk_rel_rmse_mean"]
    ax.scatter(
        rt, rmse,
        c=colors[m],
        marker=family_markers[f],
        s=60, edgecolors="white", linewidth=0.4, zorder=3,
    )

# Method legend (color)
from matplotlib.lines import Line2D
method_handles = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor=colors[m],
           markersize=7, label=method_labels[m])
    for m in method_labels
]
leg1 = ax.legend(
    handles=method_handles, title="Method", fontsize=7, title_fontsize=7.5,
    loc="upper left", framealpha=0.9, edgecolor="gray",
)
ax.add_artist(leg1)

# Family legend (shape)
family_handles = [
    Line2D([0], [0], marker=family_markers[f], color="w",
           markerfacecolor="#555", markersize=7, label=family_labels[f])
    for f in family_markers
]
ax.legend(
    handles=family_handles, title="Family", fontsize=7, title_fontsize=7.5,
    loc="lower right", framealpha=0.9, edgecolor="gray",
)

ax.set_xscale("log")
ax.set_xlabel("Runtime (seconds, log scale)")
ax.set_ylabel("Full-GKK Relative RMSE")
ax.set_title("Runtime–Quality Trade-off (3D)", fontsize=10, fontweight="bold")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(True, alpha=0.25, linewidth=0.4)

fig.tight_layout()
fig.savefig("fig5_pareto_scatter.pdf", bbox_inches="tight", dpi=300)
fig.savefig("fig5_pareto_scatter.png", bbox_inches="tight", dpi=200)
print("Figure 5 saved.")
