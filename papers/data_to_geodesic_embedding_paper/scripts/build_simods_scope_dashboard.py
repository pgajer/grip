#!/usr/bin/env python3
"""Build the SIMODS paper planning dashboard.

The dashboard is a self-contained local HTML artifact generated from the
SIMODS paper scope note plus a small structured summary of the current paper
plan.  It is intended to support decisions: what belongs in the paper, what is
missing, what evidence exists, and what should happen next.
"""

from __future__ import annotations

import datetime as _dt
import html
import pathlib
import re
from dataclasses import dataclass


ROOT = pathlib.Path(__file__).resolve().parents[1]
SOURCE = ROOT / "notes" / "simods_data_geodesic_mds_paper_scope.md"
OUTPUT = ROOT / "notes" / "simods_paper_planning_dashboard.html"


@dataclass(frozen=True)
class Task:
    key: str
    title: str
    priority: str
    status: str
    owner: str
    source: str
    summary: str
    next_action: str


@dataclass(frozen=True)
class Evidence:
    claim: str
    status: str
    sources: tuple[int, ...]
    note: str


@dataclass(frozen=True)
class Experiment:
    name: str
    purpose: str
    status: str
    needed_for_paper: str
    assets: tuple[str, ...]


def read_source() -> str:
    if not SOURCE.exists():
        raise FileNotFoundError(f"Missing source note: {SOURCE}")
    return SOURCE.read_text(encoding="utf-8")


def extract_build_timestamp(markdown: str) -> str:
    match = re.search(r"^Build timestamp:\s*(.+)$", markdown, re.MULTILINE)
    return match.group(1).strip() if match else "unknown"


def extract_outline_titles(markdown: str) -> list[str]:
    lines = markdown.splitlines()
    in_outline = False
    titles: list[str] = []
    for line in lines:
        if line.strip() == "## Proposed Paper Outline":
            in_outline = True
            continue
        if in_outline and line.startswith("## "):
            break
        match = re.match(r"^###\s+(.+)$", line)
        if in_outline and match:
            titles.append(match.group(1).strip())
    return titles


TASKS: tuple[Task, ...] = (
    Task(
        "A",
        "Decide manuscript asset structure",
        "P1",
        "started",
        "Manuscript / Orchestrator",
        "SIMODS scope note, Section A",
        "The SIMODS workspace exists, but the actual LaTeX skeleton and "
        "manuscript-local references/figure conventions are still pending.",
        "Create the first LaTeX skeleton, references.bib, and a policy for "
        "copying or regenerating publication figures into data_to_geodesic_embedding_paper.",
    ),
    Task(
        "B",
        "Consolidate graph reconstruction results",
        "P1",
        "needs synthesis",
        "Geodesic data geometry orchestrator",
        "geodesic_data_geometry reports and benchmark CSVs",
        "The paper needs a small decisive result package, not the full history "
        "of exploratory graph-construction reports.",
        "Ask the GDG lane to produce a paper-ready table/figure set showing "
        "adaptive-radius and cKNN as current leading graph families.",
    ),
    Task(
        "C",
        "Generalize non-oracle graph parameter selection",
        "P1",
        "not yet standardized",
        "Implementation + Experiment agents",
        "gflow graph-selection machinery",
        "This is the main missing bridge from oracle benchmarks to usable real "
        "data analysis.",
        "Audit existing gflow k-selection functions, define a common interface, "
        "and test GCV, JS-degree, graph edit/stability, and connectivity burden "
        "criteria across adaptive-radius and cKNN.",
    ),
    Task(
        "D",
        "Show oracle recovery of non-oracle choices",
        "P1",
        "not yet run",
        "Experiment agent",
        "quadratic-surface oracle benchmarks",
        "The paper needs to show that non-oracle selectors recover near-oracle "
        "graph parameters on controlled geometry.",
        "Run selector-vs-oracle comparisons on the 2D curvature suite, then "
        "extend selectively to 3D quadratic hypersurface examples.",
    ),
    Task(
        "E",
        "Define final edge-KK stopping criteria",
        "P1",
        "partly formulated",
        "GMDS implementation lane",
        "edge-KK timing and diagnostic reports",
        "The paper needs one coherent stopping rule rather than fixed budgets "
        "scattered across experiments.",
        "Choose edge rRMSE, q95 residual, improvement window, max-iteration, and "
        "failure/reporting thresholds; retrofit main scripts to record them.",
    ),
    Task(
        "F",
        "Decide whether repulsive unfolding is in the main paper",
        "P2",
        "open decision",
        "Theory + Experiment agents",
        "edge-isometric repulsive unfolding reports",
        "Initial visual evidence suggests weighted GRIP -> edge-KK and metric "
        "MDS -> edge-KK may already be spread enough for the paper.",
        "Keep repulsive unfolding as future/discussion unless it clearly fixes a "
        "failure mode in a compact experiment.",
    ),
    Task(
        "G",
        "Choose real data case studies",
        "P1",
        "not decided",
        "Human + Experiment agent",
        "AGP, cell-cycle, Valencia 13k, VIRGO2 candidates",
        "The SIMODS paper should include one or two interpretable real-data "
        "examples, but they must not drown the mathematical story.",
        "Choose a first low-friction dataset, likely cell cycle or VIRGO2, and "
        "define preprocessing, distance, metadata, and comparison baselines.",
    ),
    Task(
        "H",
        "Clarify software ownership",
        "P2",
        "needs concise map",
        "Provenance / Librarian",
        "gflow, grip, geodesicMDS, geodesic_data_geometry, gmdsui",
        "The code can remain cross-package, but the paper needs a clean "
        "reproducibility story.",
        "Write a software map: graph construction in gflow, layout in grip, "
        "experiments/manuscript in geodesicMDS and geodesic_data_geometry, "
        "visual inspection in gmdsui.",
    ),
    Task(
        "I",
        "Set SIMODS length and supplement strategy",
        "P2",
        "open decision",
        "Manuscript lane",
        "SIMODS author guidelines",
        "The project has more experiments than can fit in one paper unless the "
        "main figures are selected aggressively.",
        "Decide what stays in the main text, what moves to supplement, and what "
        "remains project documentation.",
    ),
)


EVIDENCE: tuple[Evidence, ...] = (
    Evidence(
        "Adaptive-radius and cKNN are current leading graph families on 2D quadratic surfaces.",
        "tested / needs distilled paper table",
        (
            1,
            2,
        ),
        "Use a distilled figure/table, not the full exploratory report sequence.",
    ),
    Evidence(
        "Fixed-radius and mKNN are weaker and can be sparse baselines.",
        "tested / needs concise citation",
        (
            2,
        ),
        "Good candidate for supplement or a short negative-result paragraph.",
    ),
    Evidence(
        "Pruning usually preserves graph-geodesic geometry while reducing edge count.",
        "tested / probably secondary",
        (
            3,
        ),
        "Useful as sparsification context; likely not central to the SIMODS story.",
    ),
    Evidence(
        "Surface-target errors are cleaner than sample-oracle errors.",
        "tested / interpretation settled for now",
        (
            4,
        ),
        "Future benchmark summaries should focus on surface-target metrics.",
    ),
    Evidence(
        "Metric MDS -> edge-KK and weighted GRIP -> edge-KK are strong GMDS baselines.",
        "tested / needs main-paper figure",
        (
            5,
            6,
        ),
        "This should be one of the central layout-method comparisons.",
    ),
    Evidence(
        "edge-KK is a quadratic edge-barrier special case with a_ij = k_ij ell_ij^2.",
        "tested / documented",
        (
            7,
        ),
        "Use to connect the older edge-KK language to the newer edge-isometric framework.",
    ),
    Evidence(
        "Repulsive unfolding is conceptually rich but not yet necessary for the headline method.",
        "inferred from current visual tests",
        (
            6,
        ),
        "Keep optional until it earns a clear role.",
    ),
)


SOURCE_REFERENCES: dict[int, str] = {
    1: "/Users/pgajer/current_projects/geodesic_data_geometry",
    2: "/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction",
    3: "/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/quadform-pruning-method-comparison",
    4: "/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/quadform-first-benchmark",
    5: "/Users/pgajer/current_projects/geodesicMDS/experiments",
    6: "/Users/pgajer/current_projects/geodesicMDS/notes/edge_isometric_repulsive_unfolding_gmds.pdf",
    7: "/Users/pgajer/current_projects/geodesicMDS/notes/edge_gkk_quadratic_barrier_equivalence_report.pdf",
}


EXPERIMENTS: tuple[Experiment, ...] = (
    Experiment(
        "Oracle graph reconstruction on quadratic surfaces",
        "Determine which graph constructions best recover intrinsic surface geodesics.",
        "substantial results exist",
        "Main evidence for data-to-graph layer.",
        (
            "2D curvature suite",
            "3D Delaunay oracle stress tests",
            "surface-target relative RMS tables/figures",
        ),
    ),
    Experiment(
        "Non-oracle graph parameter selection",
        "Select graph parameters when intrinsic geodesic truth is unavailable.",
        "major missing piece",
        "Essential bridge to real data.",
        (
            "GCV smoothing",
            "JS degree distribution",
            "graph edit/stability",
            "connectivity/MST burden",
        ),
    ),
    Experiment(
        "GMDS layout comparison",
        "Compare metric MDS, weighted GRIP, and edge-KK repaired layouts.",
        "results exist, needs distillation",
        "Main evidence for graph-to-layout layer.",
        (
            "metric MDS -> edge-KK",
            "weighted GRIP -> edge-KK",
            "timing through n up to 3200 where available",
        ),
    ),
    Experiment(
        "Stopping-rule and scaling study",
        "Turn edge-KK from an experimental repair step into a reproducible algorithm.",
        "partly run",
        "Needed for clean methods section.",
        (
            "edge rRMSE",
            "q95 residual",
            "improvement window",
            "elapsed time",
        ),
    ),
    Experiment(
        "Real data case study",
        "Demonstrate the full non-oracle pipeline on data where geometry matters.",
        "not selected",
        "Needed for broad data-science relevance.",
        (
            "cell cycle",
            "VIRGO2",
            "Valencia 13k",
            "AGP",
        ),
    ),
)


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")


def render_task_cards() -> str:
    cards = []
    for task in TASKS:
        cards.append(
            f"""
            <article class="task-card" data-priority="{esc(task.priority)}" data-status="{esc(task.status)}" data-owner="{esc(task.owner)}">
              <div class="task-topline">
                <span class="task-key">{esc(task.key)}</span>
                <span class="pill priority">{esc(task.priority)}</span>
                <span class="pill status">{esc(task.status)}</span>
              </div>
              <h3>{esc(task.title)}</h3>
              <p>{esc(task.summary)}</p>
              <dl>
                <dt>Owner</dt><dd>{esc(task.owner)}</dd>
                <dt>Source</dt><dd>{esc(task.source)}</dd>
                <dt>Next action</dt><dd>{esc(task.next_action)}</dd>
              </dl>
            </article>
            """
        )
    return "\n".join(cards)


def render_outline(outline_titles: list[str]) -> str:
    items = []
    for i, title in enumerate(outline_titles, start=1):
        section_slug = slug(title)
        items.append(
            f"""
            <article class="outline-card" id="outline-{section_slug}">
              <span class="outline-number">{i}</span>
              <h3>{esc(title)}</h3>
              <p>{outline_hint(title)}</p>
            </article>
            """
        )
    return "\n".join(items)


def outline_hint(title: str) -> str:
    hints = {
        "1. Introduction": "Frame the paper as a full data-to-geodesic-embedding pipeline.",
        "2. Related Work": "Connect MDS, Isomap, graph construction, PHATE, graph drawing, GRIP, and graph signal smoothing.",
        "3. Data Geodesic Geometry Reconstruction": "Formalize X -> G(X) and define graph families, MST repair, and pruning.",
        "4. Oracle Benchmarks On Quadratic Surfaces": "Show which graph constructions recover known surface geodesics.",
        "5. Non-Oracle Graph Parameter Selection": "Turn oracle benchmark insight into practical graph selection criteria.",
        "6. Geodesic MDS": "Define graph-to-layout isometry through embedded graph path lengths.",
        "7. GMDS Diagnostics And Stopping Rules": "Make edge fidelity, spread, runtime, and convergence measurable.",
        "8. Experiments": "Assemble oracle, non-oracle, layout, timing, and real-data evidence.",
        "9. Discussion": "Position the work as a bridge between data analysis and graph drawing.",
    }
    return esc(hints.get(title, "Paper section from the scope document."))


def render_evidence() -> str:
    rows = []
    for item in EVIDENCE:
        sources = ", ".join(f"[{source_id}]" for source_id in item.sources)
        rows.append(
            f"""
            <tr>
              <td>{esc(item.claim)}</td>
              <td><span class="pill evidence-status">{esc(item.status)}</span></td>
              <td>{esc(item.note)}</td>
              <td class="source-refs">{esc(sources)}</td>
            </tr>
            """
        )
    return "\n".join(rows)


def render_source_references() -> str:
    items = []
    for source_id, path in SOURCE_REFERENCES.items():
        items.append(f"<li><span>[{source_id}]</span><code>{esc(path)}</code></li>")
    return "\n".join(items)


def render_experiments() -> str:
    rows = []
    for exp in EXPERIMENTS:
        assets = "".join(f"<li>{esc(asset)}</li>" for asset in exp.assets)
        rows.append(
            f"""
            <article class="experiment-card">
              <h3>{esc(exp.name)}</h3>
              <p>{esc(exp.purpose)}</p>
              <div class="experiment-meta">
                <span><strong>Status:</strong> {esc(exp.status)}</span>
                <span><strong>Paper role:</strong> {esc(exp.needed_for_paper)}</span>
              </div>
              <ul>{assets}</ul>
            </article>
            """
        )
    return "\n".join(rows)


def render_tasks_markdown() -> str:
    lines = ["# SIMODS Paper Planning Dashboard Export", ""]
    lines.append("## Priority Tasks")
    for task in TASKS:
        lines.extend(
            [
                "",
                f"### {task.key}. {task.title}",
                f"- Priority: {task.priority}",
                f"- Status: {task.status}",
                f"- Owner: {task.owner}",
                f"- Source: {task.source}",
                f"- Summary: {task.summary}",
                f"- Next action: {task.next_action}",
            ]
        )
    lines.append("")
    lines.append("## Evidence Items")
    for item in EVIDENCE:
        lines.extend(
            [
                "",
                f"- Claim: {item.claim}",
                f"  - Status: {item.status}",
                f"  - Note: {item.note}",
                f"  - Sources: {', '.join(f'[{source_id}]' for source_id in item.sources)}",
            ]
        )
    lines.append("")
    lines.append("## Source References")
    for source_id, path in SOURCE_REFERENCES.items():
        lines.append(f"- [{source_id}] {path}")
    return "\n".join(lines)


def build_html() -> str:
    markdown = read_source()
    source_timestamp = extract_build_timestamp(markdown)
    outline_titles = extract_outline_titles(markdown)
    now = _dt.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")
    export_markdown = html.escape(render_tasks_markdown())

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>SIMODS Paper Planning Dashboard</title>
  <style>
    :root {{
      --bg: #f7f7f3;
      --ink: #1d252c;
      --muted: #5e6871;
      --panel: #ffffff;
      --line: #d8d6cd;
      --accent: #0f6d7a;
      --accent-2: #8b3f1f;
      --good: #24745b;
      --warn: #9a6a00;
      --risk: #9a2f2f;
      --shadow: 0 10px 30px rgba(20, 30, 40, 0.08);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      color: var(--ink);
      background: var(--bg);
      line-height: 1.55;
    }}
    a {{ color: var(--accent); }}
    code {{
      padding: 0.1rem 0.28rem;
      border-radius: 4px;
      background: #ece9df;
      font-size: 0.92em;
    }}
    header {{
      padding: 28px clamp(20px, 5vw, 56px) 22px;
      background: #14232b;
      color: #fff;
    }}
    header h1 {{
      margin: 0 0 10px;
      font-size: clamp(2rem, 4vw, 4rem);
      line-height: 1.05;
      letter-spacing: 0;
    }}
    header p {{
      max-width: 980px;
      margin: 0;
      color: #d8e2e5;
      font-size: 1.06rem;
    }}
    .meta-strip {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      margin-top: 18px;
    }}
    .meta-strip span {{
      border: 1px solid rgba(255,255,255,0.24);
      color: #eef6f7;
      border-radius: 999px;
      padding: 6px 10px;
      font-size: 0.88rem;
    }}
    nav {{
      position: sticky;
      top: 0;
      z-index: 10;
      display: flex;
      gap: 8px;
      overflow-x: auto;
      padding: 10px clamp(20px, 5vw, 56px);
      background: rgba(247, 247, 243, 0.95);
      border-bottom: 1px solid var(--line);
      backdrop-filter: blur(10px);
    }}
    nav a {{
      flex: 0 0 auto;
      text-decoration: none;
      color: var(--ink);
      border: 1px solid var(--line);
      border-radius: 999px;
      padding: 7px 12px;
      background: #fff;
      font-size: 0.92rem;
    }}
    main {{
      width: min(1320px, calc(100% - 40px));
      margin: 24px auto 56px;
    }}
    section {{
      margin: 24px 0;
    }}
    .section-title {{
      display: flex;
      align-items: baseline;
      justify-content: space-between;
      gap: 20px;
      margin-bottom: 12px;
    }}
    h2 {{
      margin: 0;
      font-size: 1.55rem;
      letter-spacing: 0;
    }}
    .section-title p {{
      max-width: 720px;
      margin: 0;
      color: var(--muted);
    }}
    .grid {{
      display: grid;
      gap: 14px;
    }}
    .grid.two {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }}
    .grid.three {{ grid-template-columns: repeat(3, minmax(0, 1fr)); }}
    .panel, .outline-card, .task-card, .experiment-card {{
      border: 1px solid var(--line);
      border-radius: 8px;
      background: var(--panel);
      box-shadow: var(--shadow);
    }}
    .panel {{
      padding: 18px;
    }}
    .pipeline {{
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 10px;
      align-items: stretch;
    }}
    .pipe-step {{
      position: relative;
      min-height: 130px;
      padding: 16px;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: #fff;
    }}
    .pipe-step:not(:last-child)::after {{
      content: ">";
      position: absolute;
      right: -12px;
      top: 50%;
      transform: translateY(-50%);
      color: var(--accent);
      font-weight: 700;
    }}
    .pipe-step h3, .task-card h3, .outline-card h3, .experiment-card h3 {{
      margin: 0 0 8px;
      font-size: 1.05rem;
    }}
    .pipe-step p, .task-card p, .outline-card p, .experiment-card p {{
      margin: 0;
      color: var(--muted);
    }}
    .outline-grid {{
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 12px;
    }}
    .outline-card {{
      padding: 16px;
    }}
    .outline-number {{
      display: inline-grid;
      place-items: center;
      width: 30px;
      height: 30px;
      margin-bottom: 10px;
      border-radius: 50%;
      color: #fff;
      background: var(--accent);
      font-weight: 700;
    }}
    .controls {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      margin-bottom: 12px;
    }}
    .controls button, .controls select, .export-button {{
      appearance: none;
      border: 1px solid var(--line);
      border-radius: 6px;
      background: #fff;
      color: var(--ink);
      padding: 8px 10px;
      font: inherit;
    }}
    .controls button.active {{
      border-color: var(--accent);
      background: #e6f4f5;
      color: #073f47;
    }}
    .task-grid {{
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 12px;
    }}
    .task-card {{
      padding: 15px;
    }}
    .task-topline {{
      display: flex;
      align-items: center;
      gap: 6px;
      margin-bottom: 10px;
    }}
    .task-key {{
      display: inline-grid;
      place-items: center;
      width: 26px;
      height: 26px;
      border-radius: 50%;
      color: #fff;
      background: var(--ink);
      font-size: 0.85rem;
      font-weight: 700;
    }}
    .pill {{
      display: inline-block;
      border-radius: 999px;
      padding: 4px 8px;
      background: #ece9df;
      color: var(--ink);
      font-size: 0.78rem;
      font-weight: 650;
      white-space: nowrap;
    }}
    .priority {{
      background: #fff0d1;
      color: var(--warn);
    }}
    .status {{
      background: #e8f1ef;
      color: var(--good);
    }}
    dl {{
      display: grid;
      grid-template-columns: 82px minmax(0, 1fr);
      gap: 4px 10px;
      margin: 12px 0 0;
      font-size: 0.9rem;
    }}
    dt {{
      color: var(--muted);
      font-weight: 700;
    }}
    dd {{
      margin: 0;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      background: #fff;
      border: 1px solid var(--line);
      border-radius: 8px;
      overflow: hidden;
      box-shadow: var(--shadow);
    }}
    th, td {{
      padding: 11px 12px;
      border-bottom: 1px solid var(--line);
      vertical-align: top;
      text-align: left;
    }}
    th:nth-child(2), td:nth-child(2) {{
      width: 170px;
    }}
    th:last-child, td:last-child {{
      width: 110px;
      white-space: nowrap;
    }}
    th {{
      background: #ebe8de;
      font-size: 0.88rem;
    }}
    td ul {{
      margin: 0;
      padding-left: 18px;
    }}
    .source-refs {{
      font-weight: 700;
      color: var(--accent);
    }}
    .source-list {{
      display: grid;
      gap: 8px;
      margin: 12px 0 0;
      padding: 14px 16px;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: #fff;
      list-style: none;
      box-shadow: var(--shadow);
    }}
    .source-list li {{
      display: grid;
      grid-template-columns: 42px minmax(0, 1fr);
      gap: 8px;
      align-items: baseline;
    }}
    .source-list span {{
      font-weight: 800;
      color: var(--accent);
    }}
    .experiment-grid {{
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 12px;
    }}
    .experiment-card {{
      padding: 16px;
    }}
    .experiment-meta {{
      display: grid;
      gap: 4px;
      margin: 10px 0;
      color: var(--muted);
      font-size: 0.92rem;
    }}
    .risk-list {{
      display: grid;
      gap: 10px;
      margin: 0;
      padding: 0;
      list-style: none;
    }}
    .risk-list li {{
      padding: 12px;
      border-left: 4px solid var(--accent-2);
      background: #fff;
      border-radius: 6px;
    }}
    textarea {{
      width: 100%;
      min-height: 220px;
      padding: 12px;
      border: 1px solid var(--line);
      border-radius: 8px;
      font: 0.9rem ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      background: #fff;
      color: var(--ink);
    }}
    footer {{
      width: min(1320px, calc(100% - 40px));
      margin: 0 auto 42px;
      color: var(--muted);
      font-size: 0.92rem;
    }}
    @media (max-width: 980px) {{
      .grid.two, .grid.three, .pipeline, .outline-grid, .task-grid, .experiment-grid {{
        grid-template-columns: 1fr;
      }}
      .pipe-step:not(:last-child)::after {{ display: none; }}
    }}
  </style>
</head>
<body>
  <header>
    <h1>SIMODS Paper Planning Dashboard</h1>
    <p>A local decision artifact for the integrated paper on data geodesic geometry and geodesic multidimensional scaling.</p>
    <div class="meta-strip">
      <span>Generated: {esc(now)}</span>
      <span>Source scope timestamp: {esc(source_timestamp)}</span>
      <span>Target: SIAM Journal on Mathematics of Data Science</span>
      <span>Mode: planning, evidence, and next actions</span>
    </div>
  </header>
  <nav>
    <a href="#thesis">Thesis</a>
    <a href="#outline">Outline</a>
    <a href="#tasks">Tasks</a>
    <a href="#evidence">Evidence</a>
    <a href="#experiments">Experiments</a>
    <a href="#risks">Risks</a>
    <a href="#export">Export</a>
  </nav>
  <main>
    <section id="thesis">
      <div class="section-title">
        <h2>Central Thesis</h2>
        <p>The paper is strongest as a complete route from data to geodesic quasi-isometric embedding.</p>
      </div>
      <div class="pipeline">
        <div class="pipe-step">
          <h3>Data</h3>
          <p>Finite sample X in an ambient feature space, often from a structured geometric or biological process.</p>
        </div>
        <div class="pipe-step">
          <h3>Data Graph</h3>
          <p>Construct G(X) so graph geodesics approximate intrinsic geometry. Main candidates: adaptive-radius and cKNN.</p>
        </div>
        <div class="pipe-step">
          <h3>Graph Metric</h3>
          <p>Use edge lengths and shortest paths to obtain d^G, with non-oracle parameter selection for real data.</p>
        </div>
        <div class="pipe-step">
          <h3>GMDS Layout</h3>
          <p>Represent graph geodesics through embedded graph path lengths using weighted GRIP -> edge-KK or metric MDS -> edge-KK.</p>
        </div>
      </div>
    </section>

    <section id="outline">
      <div class="section-title">
        <h2>Paper Outline</h2>
        <p>These sections come from the scope note. Each card names the role the section must play in the final SIMODS manuscript.</p>
      </div>
      <div class="outline-grid">
        {render_outline(outline_titles)}
      </div>
    </section>

    <section id="tasks">
      <div class="section-title">
        <h2>What Still Needs To Be Added</h2>
        <p>Filter the current task board by priority or status. This is the working version of the open-items section.</p>
      </div>
      <div class="controls" aria-label="Task filters">
        <button class="active" data-filter="all">All</button>
        <button data-filter="P1">P1</button>
        <button data-filter="P2">P2</button>
        <button data-filter="not yet">Not yet</button>
        <button data-filter="needs">Needs synthesis/standardization</button>
        <button data-filter="open">Open decisions</button>
      </div>
      <div class="task-grid" id="task-grid">
        {render_task_cards()}
      </div>
    </section>

    <section id="evidence">
      <div class="section-title">
        <h2>Claim To Evidence Map</h2>
        <p>Use this table to keep manuscript claims traceable to reports, scripts, or benchmark outputs.</p>
      </div>
      <table>
        <thead>
          <tr>
            <th>Claim</th>
            <th>Status</th>
            <th>Paper note</th>
            <th>Sources</th>
          </tr>
        </thead>
        <tbody>
          {render_evidence()}
        </tbody>
      </table>
      <ul class="source-list">
        {render_source_references()}
      </ul>
    </section>

    <section id="experiments">
      <div class="section-title">
        <h2>Experiment Roadmap</h2>
        <p>The paper needs enough experiments to support the whole pipeline without reproducing every exploratory report.</p>
      </div>
      <div class="experiment-grid">
        {render_experiments()}
      </div>
    </section>

    <section id="risks">
      <div class="section-title">
        <h2>Risks And Unknowns</h2>
        <p>These are the places where the manuscript can lose coherence if they remain implicit.</p>
      </div>
      <ul class="risk-list">
        <li><strong>Non-oracle selection is the main missing scientific bridge.</strong> Without it, the paper risks becoming an oracle benchmark plus a layout method.</li>
        <li><strong>Real-data scope needs restraint.</strong> One clean biological or structured-data example may be stronger than several half-digested examples.</li>
        <li><strong>Repulsive unfolding should earn its place.</strong> If edge-KK baselines already give useful spread, repulsive methods may belong in discussion or future work.</li>
        <li><strong>Software ownership must be explicit.</strong> gflow, grip, geodesicMDS, geodesic_data_geometry, and gmdsui each have a role; the manuscript should make that map simple.</li>
      </ul>
    </section>

    <section id="export">
      <div class="section-title">
        <h2>Export Next-Step Summary</h2>
        <p>Copy this into Codex when you want to turn the dashboard into implementation or writing tasks.</p>
      </div>
      <button class="export-button" id="copy-export">Copy summary</button>
      <textarea id="export-text" spellcheck="false">{export_markdown}</textarea>
    </section>
  </main>
  <footer>
    Source note: <code>{esc(SOURCE)}</code><br>
    Generated by: <code>{esc(pathlib.Path(__file__).resolve())}</code>
  </footer>
  <script>
    const buttons = document.querySelectorAll('[data-filter]');
    const cards = document.querySelectorAll('.task-card');
    buttons.forEach((button) => {{
      button.addEventListener('click', () => {{
        buttons.forEach((b) => b.classList.remove('active'));
        button.classList.add('active');
        const filter = button.getAttribute('data-filter');
        cards.forEach((card) => {{
          const haystack = [
            card.dataset.priority || '',
            card.dataset.status || '',
            card.dataset.owner || '',
            card.textContent || ''
          ].join(' ').toLowerCase();
          card.style.display = (filter === 'all' || haystack.includes(filter.toLowerCase())) ? '' : 'none';
        }});
      }});
    }});

    document.getElementById('copy-export').addEventListener('click', async () => {{
      const text = document.getElementById('export-text').value;
      try {{
        await navigator.clipboard.writeText(text);
        document.getElementById('copy-export').textContent = 'Copied';
        setTimeout(() => document.getElementById('copy-export').textContent = 'Copy summary', 1200);
      }} catch (error) {{
        document.getElementById('export-text').focus();
        document.getElementById('export-text').select();
      }}
    }});
  </script>
</body>
</html>
"""


def main() -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    html_text = build_html()
    OUTPUT.write_text(html_text, encoding="utf-8")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
