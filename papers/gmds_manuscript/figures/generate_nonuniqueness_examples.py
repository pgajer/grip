from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "fig_nonuniqueness_examples.pdf"


def draw_path_family(ax):
    v0 = np.array([0.0, 0.0])
    v1 = np.array([1.5, 0.0])
    angles = np.deg2rad([35, 75, 115])
    radius = 1.4
    colors = ["#7aa6d8", "#1f4e79", "#c45b3c"]

    ax.plot([v0[0], v1[0]], [v0[1], v1[1]], color="black", lw=2.2)
    for angle, color in zip(angles, colors):
        v2 = v1 + radius * np.array([np.cos(angle), np.sin(angle)])
        ax.plot([v1[0], v2[0]], [v1[1], v2[1]], color=color, lw=2.2)
        ax.scatter([v2[0]], [v2[1]], s=40, color=color, zorder=5)

    theta = np.linspace(angles.min(), angles.max(), 100)
    arc = v1 + radius * np.column_stack([np.cos(theta), np.sin(theta)])
    ax.plot(arc[:, 0], arc[:, 1], ls="--", lw=1.2, color="#666666")

    ax.scatter([v0[0], v1[0]], [v0[1], v1[1]], s=45, color="black", zorder=5)
    ax.text(v0[0] - 0.12, v0[1] - 0.18, r"$v_0$", fontsize=11)
    ax.text(v1[0] - 0.03, v1[1] - 0.18, r"$v_1$", fontsize=11)
    ax.text(arc[35, 0] + 0.04, arc[35, 1] + 0.05, r"$v_2$", fontsize=11, color="#1f4e79")
    ax.text(0.05, 0.95, "Path graph:\ncontinuous family", transform=ax.transAxes,
            ha="left", va="top", fontsize=11)
    ax.set_xlim(-0.35, 3.2)
    ax.set_ylim(-0.35, 1.75)
    ax.set_aspect("equal")
    ax.axis("off")


def draw_triangle_realizations(ax):
    a = np.array([0.0, 0.0])
    b = np.array([1.4, 0.0])
    top = np.array([0.7, 1.1])
    bottom = np.array([0.7, -1.1])
    overlap = np.array([0.7, 1.1])

    # Non-overlapping realization on the left.
    shift_left = np.array([-1.25, 0.0])
    pts_left = {
        "v1": a + shift_left,
        "v2": b + shift_left,
        "v0": top + shift_left,
        "v3": bottom + shift_left,
    }
    # Overlapping realization on the right.
    shift_right = np.array([1.35, 0.0])
    pts_right = {
        "v1": a + shift_right,
        "v2": b + shift_right,
        "v0": top + shift_right,
        "v3": overlap + shift_right,
    }

    def draw_realization(points, color_top, color_bottom, title):
        for p, q in [("v0", "v1"), ("v0", "v2"), ("v1", "v2"), ("v3", "v1"), ("v3", "v2")]:
            color = color_top if "v0" in (p, q) else color_bottom
            if {"v1", "v2"} == {p, q}:
                color = "black"
            ax.plot([points[p][0], points[q][0]], [points[p][1], points[q][1]], color=color, lw=2.0)
        for label, point in points.items():
            ax.scatter([point[0]], [point[1]], s=42, color="black", zorder=5)
            dx = -0.08 if label in {"v1", "v0"} else 0.04
            dy = 0.07 if label in {"v0", "v3"} else -0.18
            ax.text(point[0] + dx, point[1] + dy, rf"${label}$", fontsize=11)
        centroid = sum(points.values()) / 4.0
        ax.text(centroid[0], -1.45, title, ha="center", va="top", fontsize=10)

    draw_realization(pts_left, "#1f4e79", "#c45b3c", "non-overlapping")
    draw_realization(pts_right, "#1f4e79", "#c45b3c", "overlapping")
    ax.text(0.03, 0.95, "Shared-edge triangles:\ndiscrete minima", transform=ax.transAxes,
            ha="left", va="top", fontsize=11)
    ax.set_xlim(-2.2, 3.2)
    ax.set_ylim(-1.7, 1.6)
    ax.set_aspect("equal")
    ax.axis("off")


def main():
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2))
    draw_path_family(axes[0])
    draw_triangle_realizations(axes[1])
    fig.tight_layout(w_pad=1.5)
    fig.savefig(OUT, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
