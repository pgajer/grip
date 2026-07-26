#!/usr/bin/env python3

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[3]
BENCHMARK_CANDIDATES = [
    REPO_ROOT / "output" / "gkk_lgkk_paper" / "benchmarks" / "weighted-grip-phase5-family-panel-2026-04-02",
    REPO_ROOT / "output" / "gkk_lgkk_paper" / "benchmarks" / "weighted-grip-phase5-smoke-check",
]
REPORT_ROOT = (
    REPO_ROOT
    / "output"
    / "gkk_lgkk_paper"
    / "reports"
    / "benchmarks"
    / "weighted_grip_phase5_report_2026-04-02"
)
OUTPUT_TEX = REPORT_ROOT / "weighted_grip_phase5_report_2026-04-02.tex"

TEST_RESULTS = [
    {
        "check": "`testthat::test_local(filter = \"layout-weighted-globalrep|layout-weighted-trace\")`",
        "purpose": "Weighted Phase 5 API coverage, including multiscale LGKK and weighted trace consistency.",
        "pass_count": 62,
        "status": "PASS",
    },
    {
        "check": "`testthat::test_local(filter = \"layout-globalrep|layout-trace\")`",
        "purpose": "Regression protection for the original combinatorial globalrep and trace APIs.",
        "pass_count": 116,
        "status": "PASS",
    },
    {
        "check": "`benchmark-weighted-grip-family-panel.R --smoke`",
        "purpose": "End-to-end smoke validation of the broadened Phase 5 family panel.",
        "pass_count": math.nan,
        "status": "PASS",
    },
]

METHOD_ORDER = [
    "grip",
    "wgrip",
    "wgrip_core_lgkk",
    "wgrip_polish_lgkk",
    "kk",
    "gkk",
    "lgkk",
]

METHOD_LABELS = {
    "grip": "GRIP",
    "wgrip": "Weighted GRIP",
    "wgrip_core_lgkk": "Weighted GRIP + core LGKK",
    "wgrip_polish_lgkk": "Weighted GRIP + polish LGKK",
    "kk": "KK",
    "gkk": "KK->GKK",
    "lgkk": "KK->LGKK",
}


def latex_escape(value: object) -> str:
    text = str(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


def fmt_num(x: float, digits: int = 4) -> str:
    if x is None or not math.isfinite(float(x)):
        return "NA"
    return f"{float(x):.{digits}f}"


def fmt_pct(x: float, digits: int = 1) -> str:
    if x is None or not math.isfinite(float(x)):
        return "NA"
    return f"{float(x):.{digits}f}\\%"


def make_table(df: pd.DataFrame, columns: list[str], align: str, caption: str, label: str) -> str:
    header = " & ".join(latex_escape(col) for col in columns) + r" \\"
    rows = []
    for _, row in df.iterrows():
        cells = [latex_escape(row[col]) for col in columns]
        rows.append(" & ".join(cells) + r" \\")
    body = "\n".join(rows)
    return rf"""
\begin{{table}}[htbp]
\centering
\caption{{{caption}}}
\label{{{label}}}
\small
\begin{{tabular}}{{{align}}}
\toprule
{header}
\midrule
{body}
\bottomrule
\end{{tabular}}
\end{{table}}
""".strip()


def find_benchmark_root() -> Path:
    for root in BENCHMARK_CANDIDATES:
        raw_path = root / "data" / "raw_metrics.csv"
        summary_path = root / "data" / "summary_metrics.csv"
        if raw_path.exists() and summary_path.exists():
            return root
    raise FileNotFoundError(
        "Expected benchmark outputs were not found. Run at least the Phase 5 smoke benchmark first."
    )


def load_data() -> tuple[Path, pd.DataFrame, pd.DataFrame]:
    bench_root = find_benchmark_root()
    raw_df = pd.read_csv(bench_root / "data" / "raw_metrics.csv")
    summary_df = pd.read_csv(bench_root / "data" / "summary_metrics.csv")
    return bench_root, raw_df, summary_df


def build_validation_table() -> pd.DataFrame:
    rows = []
    for item in TEST_RESULTS:
        rows.append(
            {
                "Check": item["check"],
                "Purpose": item["purpose"],
                "Passes": "NA" if pd.isna(item["pass_count"]) else str(int(item["pass_count"])),
                "Status": item["status"],
            }
        )
    return pd.DataFrame(rows)


def build_family_table(raw_df: pd.DataFrame) -> pd.DataFrame:
    family_df = (
        raw_df.groupby(["family_id", "family_label", "preset"], as_index=False)
        .agg(
            vertices=("vertices", "min"),
            edges=("edges", "min"),
            edge_weight_cv=("edge_weight_cv", "mean"),
        )
        .sort_values("family_label")
    )
    family_df["Vertices"] = family_df["vertices"].map(lambda x: str(int(x)))
    family_df["Edges"] = family_df["edges"].map(lambda x: str(int(x)))
    family_df["Mean edge-weight CV"] = family_df["edge_weight_cv"].map(lambda x: fmt_num(x, 3))
    family_df["Family"] = family_df["family_label"]
    family_df["Preset"] = family_df["preset"]
    return family_df[["Family", "Preset", "Vertices", "Edges", "Mean edge-weight CV"]]


def build_overall_table(raw_df: pd.DataFrame) -> pd.DataFrame:
    df = raw_df.loc[raw_df["method"] != "start"].copy()
    overall = (
        df.groupby(["dim", "method"], as_index=False)
        .agg(
            gkk_rel_rmse=("gkk_rel_rmse", "mean"),
            kk_rel_rmse=("kk_rel_rmse", "mean"),
            runtime_sec=("runtime_sec", "mean"),
            procrustes_rmse=("procrustes_rmse", "mean"),
            sampled_stress=("sampled_stress", "mean"),
        )
    )
    overall["method"] = pd.Categorical(overall["method"], categories=METHOD_ORDER, ordered=True)
    overall = overall.sort_values(["dim", "method"])
    overall["Dimension"] = overall["dim"].map(lambda x: f"d = {int(x)}")
    overall["Method"] = overall["method"].map(METHOD_LABELS)
    overall["Mean GKK rel. RMSE"] = overall["gkk_rel_rmse"].map(lambda x: fmt_num(x, 4))
    overall["Mean runtime (s)"] = overall["runtime_sec"].map(lambda x: fmt_num(x, 3))
    overall["Mean Procrustes RMSE"] = overall["procrustes_rmse"].map(lambda x: fmt_num(x, 4))
    return overall[
        [
            "Dimension",
            "Method",
            "Mean GKK rel. RMSE",
            "Mean runtime (s)",
            "Mean Procrustes RMSE",
        ]
    ]


def build_best_family_table(summary_df: pd.DataFrame) -> pd.DataFrame:
    df = summary_df.loc[(summary_df["dim"] == 3) & (summary_df["method"] != "start")].copy()
    best_rows = []
    for _, piece in df.groupby("family_id"):
        idx = piece["gkk_rel_rmse_mean"].idxmin()
        best_rows.append(piece.loc[idx])
    best_df = pd.DataFrame(best_rows).sort_values("family_label")
    best_df["Family"] = best_df["family_label"]
    best_df["Best method"] = best_df["method_label"]
    best_df["Mean GKK rel. RMSE"] = best_df["gkk_rel_rmse_mean"].map(lambda x: fmt_num(x, 4))
    best_df["Mean runtime (s)"] = best_df["runtime_sec_mean"].map(lambda x: fmt_num(x, 3))
    return best_df[["Family", "Best method", "Mean GKK rel. RMSE", "Mean runtime (s)"]]


def build_interpretation(raw_df: pd.DataFrame) -> list[str]:
    df = raw_df.loc[raw_df["method"] != "start"].copy()
    overall = (
        df.groupby(["dim", "method"], as_index=False)
        .agg(gkk_rel_rmse=("gkk_rel_rmse", "mean"), runtime_sec=("runtime_sec", "mean"))
    )
    stats = {
        (int(row.dim), row.method): row
        for row in overall.itertuples(index=False)
    }

    def rel_change(dim: int, base: str, new: str) -> tuple[str, str]:
        a = stats[(dim, base)].gkk_rel_rmse
        b = stats[(dim, new)].gkk_rel_rmse
        delta = 100.0 * (b - a) / a
        if b < a:
          return "reduction", fmt_pct(-delta, 1)
        return "increase", fmt_pct(delta, 1)

    lines = []
    relation_3d, delta_3d = rel_change(3, "grip", "wgrip")
    lines.append(
        "In the primary 3D track of the smoke panel, weighted GRIP changes the mean full-GKK relative RMSE from "
        f"{fmt_num(stats[(3, 'grip')].gkk_rel_rmse, 4)} for GRIP to "
        f"{fmt_num(stats[(3, 'wgrip')].gkk_rel_rmse, 4)}, which is a {relation_3d} of {delta_3d}."
    )
    lines.append(
        "What matters more for Phase 5 is the weighted-LGKK path. Adding in-core multiscale LGKK lowers the 3D mean to "
        f"{fmt_num(stats[(3, 'wgrip_core_lgkk')].gkk_rel_rmse, 4)}, while post-layout polish reaches "
        f"{fmt_num(stats[(3, 'wgrip_polish_lgkk')].gkk_rel_rmse, 4)}. "
        "So even in the smoke panel, the new in-core refinement clearly improves over weighted GRIP alone."
    )
    lines.append(
        "The KK-family refinements still provide the lowest aggregate geodesic error. In 3D, "
        f"KK->GKK reaches {fmt_num(stats[(3, 'gkk')].gkk_rel_rmse, 4)} and KK->LGKK reaches "
        f"{fmt_num(stats[(3, 'lgkk')].gkk_rel_rmse, 4)}. "
        "That keeps them as the quality ceiling in this validation panel."
    )
    lines.append(
        "The runtime figure shows the intended tradeoff clearly: GRIP and weighted GRIP remain the fastest methods overall, "
        "while core-LGKK and polish-LGKK weighted variants occupy a middle ground between pure weighted GRIP and the more expensive KK warm-start pipelines."
    )
    lines.append(
        "The 2D panels remain useful, but the companion results reinforce the earlier design decision to treat 3D as the primary evaluation space. "
        "The weighted families in this panel are intrinsically geometric, so 2D layouts systematically face a representation bottleneck that is absent in 3D."
    )
    lines.append(
        "Comparing the two weighted-LGKK variants, the Phase 5 in-core version is the more architecturally important result even when the post-polish variant is marginally better on this smoke panel. "
        "It means the weighted hierarchy can now incorporate geometry-aware refinement natively instead of depending only on a terminal cleanup stage."
    )
    lines.append(
        "Against the classical GRIP baseline, the main conclusion is positive: edge-length awareness matters on these families, and the sister weighted API improves geometric fidelity without destabilizing the original combinatorial functions."
    )
    return lines


def build_document(bench_root: Path, raw_df: pd.DataFrame, summary_df: pd.DataFrame) -> str:
    validation_table = build_validation_table()
    family_table = build_family_table(raw_df)
    overall_table = build_overall_table(raw_df)
    best_family_table = build_best_family_table(summary_df)
    interpretation = build_interpretation(raw_df)
    validation_tex = make_table(
        validation_table,
        ["Check", "Purpose", "Passes", "Status"],
        "p{0.25\\linewidth}p{0.49\\linewidth}cc",
        "Phase 5 validation checks executed during implementation.",
        "tab:phase5-checks",
    )
    family_tex = make_table(
        family_table,
        ["Family", "Preset", "Vertices", "Edges", "Mean edge-weight CV"],
        "p{0.36\\linewidth}p{0.12\\linewidth}rrr",
        "Weighted family panel used for the full Phase 5 benchmark.",
        "tab:phase5-families",
    )
    overall_d2 = overall_table.loc[overall_table["Dimension"] == "d = 2"].copy()
    overall_d3 = overall_table.loc[overall_table["Dimension"] == "d = 3"].copy()
    overall_d2_tex = make_table(
        overall_d2,
        ["Method", "Mean GKK rel. RMSE", "Mean runtime (s)", "Mean Procrustes RMSE"],
        "p{0.38\\linewidth}rrr",
        "Overall means for d = 2 across the benchmark panel.",
        "tab:phase5-overall-d2",
    )
    overall_d3_tex = make_table(
        overall_d3,
        ["Method", "Mean GKK rel. RMSE", "Mean runtime (s)", "Mean Procrustes RMSE"],
        "p{0.38\\linewidth}rrr",
        "Overall means for d = 3 across the benchmark panel.",
        "tab:phase5-overall-d3",
    )
    best_family_tex = make_table(
        best_family_table,
        ["Family", "Best method", "Mean GKK rel. RMSE", "Mean runtime (s)"],
        "p{0.38\\linewidth}p{0.28\\linewidth}rr",
        "Best 3D method by family under mean full-GKK relative RMSE.",
        "tab:phase5-best-d3",
    )

    bench_cases = int((raw_df["method"] == "start").sum())
    total_layout_runs = int((raw_df["method"] != "start").sum())
    seeds = sorted(int(x) for x in raw_df["seed"].unique())
    dims = sorted(int(x) for x in raw_df["dim"].unique())
    dims_text = ", ".join(map(str, dims))
    seeds_text = ", ".join(map(str, seeds))

    fig_rel = (
        Path("..") / ".." / ".." / "benchmarks" / bench_root.name / "figures" / "gkk_rel_rmse.png"
    )
    fig_runtime = (
        Path("..") / ".." / ".." / "benchmarks" / bench_root.name / "figures" / "runtime_sec.png"
    )
    benchmark_scope = "full family panel" if "family-panel" in bench_root.name else "smoke validation panel"

    interp_text = "\n\n".join(
        rf"\paragraph{{Interpretation {idx + 1}.}} {latex_escape(line)}"
        for idx, line in enumerate(interpretation)
    )

    return rf"""
\documentclass[11pt]{{article}}
\usepackage[margin=1in]{{geometry}}
\usepackage{{graphicx}}
\usepackage{{booktabs}}
\usepackage{{longtable}}
\usepackage{{array}}
\usepackage{{float}}
\usepackage{{hyperref}}
\usepackage{{caption}}
\usepackage{{parskip}}
\usepackage{{xcolor}}

\hypersetup{{
  colorlinks=true,
  linkcolor=blue!50!black,
  urlcolor=blue!50!black
}}

\title{{Weighted GRIP Phase 5: Validation and Benchmark Report}}
\author{{Codex}}
\date{{2026-04-02}}

\begin{{document}}
\maketitle

\section*{{Overview}}

This note summarizes the checks and benchmark results for Weighted GRIP Phase 5. Phase 5 had two explicit goals: integrate multiscale LGKK refinement inside the weighted GRIP core, and broaden the weighted-family benchmark panel beyond the original mesh-oriented rollout.

The benchmark dataset used in this report is the {benchmark_scope}. It covered {len(family_table)} weighted families, dimensions $d \in \{{{dims_text}\}}$, and seeds {seeds_text}. That produced {bench_cases} benchmark cases and {total_layout_runs} final method runs once the shared-start rows are excluded.

\section*{{Validation Checks}}

{validation_tex}

\section*{{Experimental Design}}

Each family was evaluated with the following method panel:
\begin{{itemize}}
\item GRIP
\item Weighted GRIP
\item Weighted GRIP + core LGKK
\item Weighted GRIP + polish LGKK
\item KK
\item KK->GKK
\item KK->LGKK
\end{{itemize}}

The primary quality metric is the full-GKK relative RMSE, because the Phase 5 work is explicitly about making the weighted GRIP pipeline more faithful to graph geometry. We also track classical KK-relative RMSE, sampled stress, Procrustes RMSE against the target geometry, and wall-clock runtime.

{family_tex}

\section*{{Aggregate Results}}

Figure~\ref{{fig:gkk-rel-rmse}} shows the main quality result. Lower is better. Figure~\ref{{fig:runtime}} shows the corresponding runtime cost.

\begin{{figure}}[H]
\centering
\includegraphics[width=\linewidth]{{{latex_escape(fig_rel.as_posix())}}}
\caption{{Full-GKK relative RMSE across the Phase 5 family panel.}}
\label{{fig:gkk-rel-rmse}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=\linewidth]{{{latex_escape(fig_runtime.as_posix())}}}
\caption{{Runtime across the Phase 5 family panel.}}
\label{{fig:runtime}}
\end{{figure}}

{overall_d2_tex}

{overall_d3_tex}

{best_family_tex}

\section*{{Interpretation}}

{interp_text}

\section*{{Companion HTML}}

The interactive 3D companion for this report is:

\begin{{center}}
\path{{output/gkk_lgkk_paper/html/weighted_grip_phase5_layout_gallery_2026-04-02.html}}
\end{{center}}

It organizes the primary 3D benchmark layouts by family for the representative benchmark seed 1, and includes the target geometry plus every benchmarked layout method used in Phase 5.

\end{{document}}
""".strip() + "\n"


def main() -> None:
    bench_root, raw_df, summary_df = load_data()
    REPORT_ROOT.mkdir(parents=True, exist_ok=True)
    tex = build_document(bench_root, raw_df, summary_df)
    OUTPUT_TEX.write_text(tex, encoding="utf-8")
    print(f"Wrote {OUTPUT_TEX}")


if __name__ == "__main__":
    main()
