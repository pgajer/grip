from pathlib import Path

import matplotlib.pyplot as plt


EDGE = "#b1bac8"
TEXT = "#243447"
HILITE = "#d95f02"
ALT = "#2c7fb8"
NODE = "#ffffff"
CENTER = "#ffe7bf"


def draw_edges(ax, positions, edges, color=EDGE, lw=2.2, zorder=1):
    for u, v in edges:
        x1, y1 = positions[u]
        x2, y2 = positions[v]
        ax.plot(
            [x1, x2],
            [y1, y2],
            color=color,
            linewidth=lw,
            solid_capstyle="round",
            zorder=zorder,
        )


def draw_nodes(ax, positions, labels, special=None):
    special = special or set()
    for node, (x, y) in positions.items():
        face = CENTER if node in special else NODE
        ax.scatter([x], [y], s=430, facecolor=face, edgecolor="#26313f", linewidth=1.9, zorder=3)
        dx, dy = labels[node][1]
        ax.text(
            x + dx,
            y + dy,
            labels[node][0],
            fontsize=16,
            ha="center",
            va="center",
            color="#111111",
        )


def paths_panel(out_path: Path):
    fig, ax = plt.subplots(figsize=(5.1, 3.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    pos = {
        "a": (0.10, 0.75),
        "b": (0.18, 0.38),
        "u": (0.34, 0.61),
        "l": (0.51, 0.50),
        "c": (0.70, 0.51),
        "d": (0.88, 0.74),
        "e": (0.89, 0.31),
        "r": (0.55, 0.21),
    }
    labels = {
        "a": (r"$i_1$", (-0.05, 0.09)),
        "b": (r"$i_2$", (-0.05, -0.02)),
        "u": (r"$u$", (0.00, 0.09)),
        "l": (r"$\ell$", (0.00, 0.11)),
        "c": (r"$v$", (0.00, -0.11)),
        "d": (r"$j_1$", (0.05, 0.08)),
        "e": (r"$j_2$", (0.06, -0.02)),
        "r": ("", (0.00, -0.10)),
    }
    base_edges = [("a", "u"), ("b", "l"), ("u", "l"), ("l", "c"), ("c", "d"), ("c", "e"), ("l", "r")]
    active_edges = [("a", "u"), ("u", "l"), ("b", "l"), ("l", "c"), ("c", "d"), ("c", "e")]

    draw_edges(ax, pos, base_edges, color=EDGE, lw=2.4, zorder=1)
    draw_edges(ax, pos, active_edges, color=HILITE, lw=4.8, zorder=2)
    draw_nodes(ax, pos, labels, special={"l"})

    fig.savefig(out_path, dpi=300, bbox_inches="tight")


def switch_panel(out_path: Path):
    fig, ax = plt.subplots(figsize=(5.1, 4.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    pos = {
        "i": (0.10, 0.50),
        "u": (0.40, 0.70),
        "v": (0.40, 0.28),
        "j": (0.90, 0.50),
    }
    labels = {
        "i": (r"$i$", (-0.05, 0.08)),
        "u": (r"$u$", (0.00, 0.11)),
        "v": (r"$v$", (0.00, -0.11)),
        "j": (r"$j$", (0.05, 0.08)),
    }
    top = [("i", "u"), ("u", "j")]
    bottom = [("i", "v"), ("v", "j")]

    ax.text(0.5, 0.97, r"(b) Path switching / tie point", ha="center", va="top", fontsize=15, fontweight="bold")
    ax.text(
        0.5,
        0.90,
        r"equal-length routes create a nonsmooth point in $d_k^{\mathrm{emb}}(i,j;\mathbf{Z})$",
        ha="center",
        va="top",
        fontsize=9.8,
        color=TEXT,
    )

    draw_edges(ax, pos, top, color=HILITE, lw=4.8, zorder=1)
    draw_edges(ax, pos, bottom, color=ALT, lw=4.8, zorder=1)
    draw_nodes(ax, pos, labels)

    ax.text(0.53, 0.69, r"$d_{\mathrm{top}}(i,j)$", fontsize=12, color=HILITE)
    ax.text(0.53, 0.19, r"$d_{\mathrm{bot}}(i,j)$", fontsize=12, color=ALT)
    ax.text(
        0.5,
        0.07,
        r"$d_{\mathrm{top}}=d_{\mathrm{bot}}\ \Rightarrow\ \partial d_k^{\mathrm{emb}}(i,j;\mathbf{Z})$ uses convex combinations",
        ha="center",
        va="center",
        fontsize=9.8,
        color=TEXT,
    )

    fig.savefig(out_path, dpi=300, bbox_inches="tight")


def main():
    here = Path(__file__).resolve().parent
    paths_panel(here / "fig_gradient_nonlocality_paths.pdf")
    switch_panel(here / "fig_gradient_nonlocality_switch.pdf")


if __name__ == "__main__":
    main()
