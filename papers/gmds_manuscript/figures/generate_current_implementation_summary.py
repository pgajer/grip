from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

FIG_DIR = Path(__file__).resolve().parent
GRIP_TMP = Path("/Users/pgajer/current_projects/grip/dev/design/tmp")

PHASE_E_FINAL = GRIP_TMP / "gmds-misf-phase-e-integrated-2026-04-02" / "gmds_misf_phase_e_integrated_final_metrics.csv"
BLOCK1_SUMMARY = GRIP_TMP / "gmds-misf-phase-e-mesh-geometries-2026-04-02" / "gmds_misf_phase_e_mesh_geometries_summary_metrics.csv"
BLOCK1_FINAL = GRIP_TMP / "gmds-misf-phase-e-mesh-geometries-2026-04-02" / "gmds_misf_phase_e_mesh_geometries_final_metrics.csv"

OUT_PDF = FIG_DIR / "fig_current_implementation_summary.pdf"
OUT_PNG = FIG_DIR / "fig_current_implementation_summary.png"


CASE_LABELS = {
    "paraboloid_regular_12x12": "Regular\n12x12",
    "paraboloid_regular_15x15": "Regular\n15x15",
    "paraboloid_irregular_rectangle_15x15": "Irregular\n15x15",
}

METHOD_ORDER = [
    ("cmd_pure_gmds", "Direct\ncMDS"),
    ("phase_c_best_seeded", "Phase C\nseeded"),
    ("phase_d_fixed_top_weighted_kk", "Phase D\nfixed+WKK"),
    ("phase_e_proxy", "Phase E\nproxy"),
]

METHOD_COLORS = {
    "cmd_pure_gmds": "#8d99ae",
    "phase_c_best_seeded": "#457b9d",
    "phase_d_fixed_top_weighted_kk": "#2a9d8f",
    "phase_e_proxy": "#e76f51",
}


def load_progression():
    df_phase_e = pd.read_csv(PHASE_E_FINAL)
    df_block1 = pd.read_csv(BLOCK1_FINAL)

    direct = (
        df_block1[df_block1["baseline"] == "direct_cmd_pure_gmds"]
        .copy()
        .assign(
            case_id=lambda x: x["case_id"].map(
                {
                    "regular_paraboloid_12x12": "paraboloid_regular_12x12",
                    "regular_paraboloid_15x15": "paraboloid_regular_15x15",
                    "irregular_rectangle_paraboloid_15x15": "paraboloid_irregular_rectangle_15x15",
                }
            )
        )
    )
    direct = direct[["case_id", "gmds_stress", "procrustes_rmse"]].copy()
    direct["baseline"] = "cmd_pure_gmds"

    phase_e = df_phase_e[["case_id", "baseline", "gmds_stress", "procrustes_rmse"]].copy()

    all_rows = pd.concat([direct, phase_e], ignore_index=True)
    all_rows = all_rows[all_rows["case_id"].isin(CASE_LABELS)]
    return all_rows


def plot_progression(ax, df):
    case_ids = list(CASE_LABELS)
    x = np.arange(len(case_ids))
    width = 0.18

    for idx, (baseline, label) in enumerate(METHOD_ORDER):
        sub = (
            df[df["baseline"] == baseline]
            .set_index("case_id")
            .reindex(case_ids)
        )
        xpos = x + (idx - 1.5) * width
        bars = ax.bar(
            xpos,
            sub["procrustes_rmse"].values,
            width=width,
            color=METHOD_COLORS[baseline],
            edgecolor="black",
            linewidth=0.5,
            label=label,
        )
        for bar, sigma in zip(bars, sub["gmds_stress"].values):
            if pd.notna(sigma):
                ax.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    bar.get_height() + 0.005,
                    f"$\\sigma$={sigma:.3f}",
                    ha="center",
                    va="bottom",
                    fontsize=7,
                    rotation=90,
                )

    ax.set_xticks(x)
    ax.set_xticklabels([CASE_LABELS[c] for c in case_ids], fontsize=9)
    ax.set_ylabel("Final Procrustes RMSE $\\rho$")
    ax.set_title("Current best pipeline progression on paraboloid families", fontsize=11)
    ax.set_ylim(0, max(df["procrustes_rmse"]) * 1.35)
    ax.grid(axis="y", alpha=0.25, linestyle=":")
    ax.legend(frameon=False, fontsize=8, ncol=2, loc="upper right")


def load_selector_counts():
    df = pd.read_csv(BLOCK1_SUMMARY)

    methods = [
        "cmdscale",
        "kk",
        "grip",
        "weighted_grip",
        "weighted_grip_polish_lgkk",
    ]
    rows = []
    for source_col, source_label in [
        ("proxy_method_id", "Proxy pick"),
        ("oracle_rho_method_id", "Oracle best $\\rho$"),
        ("oracle_sigma_method_id", "Oracle best $\\sigma$"),
    ]:
        counts = df[source_col].value_counts()
        for method in methods:
            rows.append(
                {
                    "source": source_label,
                    "method": method,
                    "count": int(counts.get(method, 0)),
                }
            )
    return pd.DataFrame(rows), methods


def plot_selector_counts(ax, counts_df, methods):
    source_order = ["Proxy pick", "Oracle best $\\rho$", "Oracle best $\\sigma$"]
    colors = {
        "Proxy pick": "#e76f51",
        "Oracle best $\\rho$": "#457b9d",
        "Oracle best $\\sigma$": "#2a9d8f",
    }
    x = np.arange(len(methods))
    width = 0.25

    for idx, source in enumerate(source_order):
        sub = (
            counts_df[counts_df["source"] == source]
            .set_index("method")
            .reindex(methods)
        )
        ax.bar(
            x + (idx - 1) * width,
            sub["count"].values,
            width=width,
            color=colors[source],
            edgecolor="black",
            linewidth=0.5,
            label=source,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(
        ["cMDS", "KK", "GRIP", "W-GRIP", "W-GRIP\n+LGKK"],
        fontsize=9,
    )
    ax.set_ylabel("Count across 12 Block-1 cases")
    ax.set_title("Phase E Block 1 seed-selection mismatch", fontsize=11)
    ax.grid(axis="y", alpha=0.25, linestyle=":")
    ax.legend(frameon=False, fontsize=8)


def main():
    progression = load_progression()
    selector_counts, methods = load_selector_counts()

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
        }
    )

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.8))
    plot_progression(axes[0], progression)
    plot_selector_counts(axes[1], selector_counts, methods)

    fig.subplots_adjust(left=0.055, right=0.99, bottom=0.19, top=0.86, wspace=0.26)

    fig.savefig(OUT_PDF, bbox_inches="tight")
    fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
    print(f"Wrote {OUT_PDF}")
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
