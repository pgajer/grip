from pathlib import Path

import matplotlib.pyplot as plt


EDGE_COLOR = "#a9b2c2"
NODE_EDGE = "#1f2937"
HILITE_LEFT = "#d95f02"
HILITE_RIGHT = "#1f77b4"
NODE_FILL = "#fff5e6"


def draw_panel(
    ax,
    positions,
    label_offsets,
    panel_title,
    subtitle_1,
    subtitle_2,
    highlighted_edges,
    highlight_color,
    bottom_line_1,
    bottom_line_2,
    bottom_color,
):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(
        0.5,
        0.97,
        panel_title,
        ha="center",
        va="top",
        fontsize=16,
        fontweight="bold",
    )
    ax.text(0.5, 0.89, subtitle_1, ha="center", va="top", fontsize=11)
    ax.text(0.5, 0.79, subtitle_2, ha="center", va="top", fontsize=11)

    all_edges = [
        ("i", "a"),
        ("i", "c"),
        ("a", "c"),
        ("a", "b"),
        ("b", "d"),
        ("c", "d"),
        ("b", "j"),
        ("d", "j"),
    ]

    highlighted = {frozenset(edge) for edge in highlighted_edges}
    for u, v in all_edges:
        x1, y1 = positions[u]
        x2, y2 = positions[v]
        is_hilited = frozenset((u, v)) in highlighted
        ax.plot(
            [x1, x2],
            [y1, y2],
            color=highlight_color if is_hilited else EDGE_COLOR,
            linewidth=4.6 if is_hilited else 2.2,
            solid_capstyle="round",
            zorder=1,
        )

    for node, (x, y) in positions.items():
        ax.scatter(
            [x],
            [y],
            s=430,
            facecolor=NODE_FILL if node in {"i", "j"} else "white",
            edgecolor=NODE_EDGE,
            linewidth=1.8,
            zorder=3,
        )

        dx, dy = label_offsets[node]
        ax.text(
            x + dx,
            y + dy,
            rf"${node}$",
            fontsize=15,
            ha="center",
            va="center",
        )

    ax.text(
        0.5,
        0.10,
        bottom_line_1,
        ha="center",
        va="center",
        fontsize=12,
        color=bottom_color,
    )
    ax.text(
        0.5,
        0.03,
        bottom_line_2,
        ha="center",
        va="center",
        fontsize=13,
        color=bottom_color,
    )


def main():
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.2))
    fig.subplots_adjust(left=0.035, right=0.985, top=0.96, bottom=0.06, wspace=0.13)

    left_positions = {
        "i": (0.10, 0.48),
        "a": (0.34, 0.60),
        "b": (0.64, 0.54),
        "j": (0.91, 0.40),
        "c": (0.40, 0.25),
        "d": (0.77, 0.20),
    }
    left_offsets = {
        "i": (-0.06, 0.08),
        "a": (-0.02, 0.06),
        "b": (0.00, 0.07),
        "j": (0.05, 0.08),
        "c": (0.00, -0.10),
        "d": (0.03, -0.09),
    }
    draw_panel(
        axes[0],
        left_positions,
        left_offsets,
        r"(a) Input graph geodesic",
        r"$\delta^{(k)}_{ij}$ on the fixed $k$-NN graph $G_k$",
        "edge weights are measured in the input space",
        [("i", "a"), ("a", "b"), ("b", "j")],
        HILITE_LEFT,
        "",
        r"$\delta^{(k)}_{ij}=w_{ia}+w_{ab}+w_{bj}$",
        HILITE_LEFT,
    )

    right_positions = {
        "i": (0.12, 0.52),
        "a": (0.38, 0.59),
        "b": (0.66, 0.52),
        "j": (0.93, 0.36),
        "c": (0.46, 0.30),
        "d": (0.82, 0.21),
    }
    right_offsets = {
        "i": (-0.06, 0.09),
        "a": (0.00, 0.07),
        "b": (0.00, 0.09),
        "j": (0.05, 0.08),
        "c": (0.00, -0.11),
        "d": (0.06, -0.02),
    }
    draw_panel(
        axes[1],
        right_positions,
        right_offsets,
        r"(b) Embedded input geodesic",
        r"$d_k^{\mathrm{emb}}(i,j;\mathbf{Z})$ along $\gamma_{ij}^{(k)}$",
        "same route as in (a), but edge lengths come from the embedding",
        [("i", "a"), ("a", "b"), ("b", "j")],
        HILITE_RIGHT,
        "route fixed by the input geometry",
        r"$d_k^{\mathrm{emb}}(i,j;\mathbf{Z})=\Vert z_i-z_a\Vert+\Vert z_a-z_b\Vert+\Vert z_b-z_j\Vert$",
        HILITE_RIGHT,
    )

    out_path = Path(__file__).with_name("fig_geodesic_stress_definition.pdf")
    fig.savefig(out_path, dpi=300)


if __name__ == "__main__":
    main()
