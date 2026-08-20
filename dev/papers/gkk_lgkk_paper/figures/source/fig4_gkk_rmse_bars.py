#!/usr/bin/env python3
"""Figure 4: Full-GKK Relative RMSE by Graph Family (3D).

Grouped bar chart showing geodesic fidelity across all seven methods
for each of the five benchmark families, in 3D embedding.
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
FIGURE_DIR = REPO_ROOT / "output" / "gkk_lgkk_paper" / "figures"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

# --- Load data ---
df = pd.read_csv(SUMMARY_PATH)
df3 = df[(df["dim"] == 3) & (df["method"] != "start")].copy()

# --- Method order and display names ---
method_order = [
    "grip", "wgrip", "wgrip_core_lgkk", "wgrip_polish_lgkk",
    "kk", "gkk", "lgkk"
]
method_labels = {
    "grip":              "GRIP",
    "wgrip":             "W-GRIP",
    "wgrip_core_lgkk":   "W-GRIP+cLGKK",
    "wgrip_polish_lgkk": "W-GRIP+pLGKK",
    "kk":                "KK",
    "gkk":               "KK→GKK",
    "lgkk":              "KK→LGKK",
}

# --- Family order ---
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

# --- Colors ---
colors = {
    "grip":              "#7f8c8d",
    "wgrip":             "#27ae60",
    "wgrip_core_lgkk":   "#2ecc71",
    "wgrip_polish_lgkk": "#1abc9c",
    "kk":                "#3498db",
    "gkk":               "#2980b9",
    "lgkk":              "#8e44ad",
}

# --- Build matrix ---
vals = {}
for m in method_order:
    vals[m] = []
    for f in family_order:
        row = df3[(df3["method"] == m) & (df3["family_id"] == f)]
        if len(row) > 0:
            vals[m].append(row["gkk_rel_rmse_mean"].values[0])
        else:
            vals[m].append(0)

# --- Plot ---
mpl.rcParams.update({
    "font.family": "serif",
    "font.size": 9,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
})

fig, ax = plt.subplots(figsize=(7.0, 3.2))

n_families = len(family_order)
n_methods = len(method_order)
x = np.arange(n_families)
bar_width = 0.11
offsets = np.arange(n_methods) - (n_methods - 1) / 2

for i, m in enumerate(method_order):
    ax.bar(
        x + offsets[i] * bar_width,
        vals[m],
        bar_width * 0.92,
        label=method_labels[m],
        color=colors[m],
        edgecolor="white",
        linewidth=0.3,
    )

ax.set_xticks(x)
ax.set_xticklabels([family_labels[f] for f in family_order])
ax.set_ylabel("Full-GKK Relative RMSE")
ax.set_title("Geodesic Fidelity by Graph Family (3D)", fontsize=10, fontweight="bold")
ax.legend(
    fontsize=7, ncol=4, loc="upper left",
    framealpha=0.9, edgecolor="gray", handlelength=1.2
)
ax.set_ylim(0, None)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="y", alpha=0.3, linewidth=0.4)

fig.tight_layout()
fig.savefig(FIGURE_DIR / "fig4_gkk_rmse_bars.pdf", bbox_inches="tight", dpi=300)
fig.savefig(FIGURE_DIR / "fig4_gkk_rmse_bars.png", bbox_inches="tight", dpi=200)
print("Figure 4 saved.")
