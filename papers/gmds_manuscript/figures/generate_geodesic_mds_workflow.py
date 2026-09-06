from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


BOX_IN = "#f4e4c1"
BOX_GRAPH = "#f7cfb4"
BOX_INIT = "#d7eadf"
BOX_LOOP = "#d7e8f7"
BOX_OUT = "#d8f0d2"
TEXT = "#233142"
ARROW = "#5b6b80"
LOOP = "#2c7fb8"


def add_box(
    ax,
    x,
    y,
    w,
    h,
    title,
    lines,
    facecolor,
    edgecolor="#51606f",
    title_fs=12,
    body_fs=10,
):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.02",
        linewidth=1.7,
        edgecolor=edgecolor,
        facecolor=facecolor,
    )
    ax.add_patch(patch)
    ax.text(
        x + w / 2,
        y + h * 0.72,
        title,
        ha="center",
        va="center",
        fontsize=title_fs,
        fontweight="bold",
        color=TEXT,
    )
    ax.text(
        x + w / 2,
        y + h * 0.29,
        "\n".join(lines),
        ha="center",
        va="center",
        fontsize=body_fs,
        color=TEXT,
        linespacing=1.35,
    )


def add_arrow(ax, start, end, color=ARROW, rad=0.0, lw=2.0):
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=14,
        linewidth=lw,
        color=color,
        connectionstyle=f"arc3,rad={rad}",
    )
    ax.add_patch(arrow)


def main():
    fig, ax = plt.subplots(figsize=(6.8, 9.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.11, 0.95, "Input metric", fontsize=12, fontweight="bold", color=TEXT)

    box_x = 0.12
    box_w = 0.68
    box_h = 0.10

    y_input = [0.84, 0.70, 0.56, 0.42]
    input_specs = [
        (r"1. Input data $\mathbf{X}$", [r"$n$ points in $\mathbb{R}^m$", "ambient coordinates"], BOX_IN),
        (r"2. Build symmetric $k$-NN graph $G_k$", ["connect local neighbors", "Euclidean edge weights"], BOX_GRAPH),
        (r"3. Compute $\Delta^{(k)}$ and $\Gamma_k$", ["all-pairs shortest paths", "store canonical input geodesics"], BOX_GRAPH),
        (r"4. Initialize $\mathbf{Z}^{(0)}$", ["classical MDS on", r"$\Delta^{(k)}$"], BOX_INIT),
    ]

    for y, (title, lines, color) in zip(y_input, input_specs):
        add_box(ax, box_x, y, box_w, box_h, title, lines, color, title_fs=12, body_fs=10.5)

    for y1, y2 in zip(y_input, y_input[1:]):
        add_arrow(ax, (box_x + box_w / 2, y1), (box_x + box_w / 2, y2 + box_h), lw=2.2)

    ax.text(0.11, 0.37, "Iterative geodesic SMACOF", fontsize=12, fontweight="bold", color=LOOP)

    loop_bg = FancyBboxPatch(
        (0.07, 0.02),
        0.86,
        0.30,
        boxstyle="round,pad=0.015,rounding_size=0.02",
        linewidth=1.6,
        edgecolor=LOOP,
        facecolor="#f7fbff",
        linestyle="--",
        alpha=0.95,
    )
    ax.add_patch(loop_bg)

    add_arrow(ax, (box_x + box_w / 2, y_input[-1]), (box_x + box_w / 2, 0.33), color=LOOP, lw=2.3)

    box_h_loop = 0.07
    y_loop = [0.22, 0.125]
    loop_specs = [
        (r"5. Evaluate fixed geodesics", [r"reuse stored $\gamma_{ij}^{(k)}$", r"compute $d_k^{\mathrm{emb}}$"]),
        (r"6. Accumulate path weights", [r"$\tilde{\mathbf{V}}^{(t)},\tilde{\mathbf{W}}^{(t)}$", "from the fixed geodesic family"]),
    ]

    for y, (title, lines) in zip(y_loop, loop_specs):
        add_box(ax, box_x, y, box_w, box_h_loop, title, lines, BOX_LOOP, title_fs=10.1, body_fs=9.0)

    for y1, y2 in zip(y_loop, y_loop[1:]):
        add_arrow(ax, (box_x + box_w / 2, y1), (box_x + box_w / 2, y2 + box_h_loop), color=LOOP, lw=2.2)

    update_x = 0.12
    update_y = 0.03
    update_w = 0.60
    update_h = 0.075
    add_box(
        ax,
        update_x,
        update_y,
        update_w,
        update_h,
        r"7. Update $\mathbf{Z}^{(t+1)}$",
        [r"$(\tilde{\mathbf{W}}^{(t)})^+\tilde{\mathbf{V}}^{(t)}\mathbf{Z}^{(t)}$", "evaluate geodesic stress"],
        BOX_LOOP,
        title_fs=10.6,
        body_fs=9.0,
    )
    add_arrow(ax, (box_x + box_w / 2, y_loop[-1]), (update_x + update_w / 2, update_y + update_h), color=LOOP, lw=2.2)

    add_arrow(
        ax,
        (update_x + update_w * 0.86, update_y + update_h * 0.12),
        (box_x + box_w * 0.86, y_loop[0] + box_h_loop * 0.88),
        color=LOOP,
        rad=0.50,
        lw=2.2,
    )
    ax.text(0.84, 0.18, "if not converged", fontsize=9.8, color=LOOP, rotation=90, ha="center", va="center")

    out_path = Path(__file__).with_name("fig_geodesic_mds_workflow.pdf")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
