#!/usr/bin/env python3
"""Build the SIMODS figure/table readiness board from a YAML registry."""

from __future__ import annotations

import datetime as _dt
import html
import json
import pathlib
from collections import Counter, defaultdict
from typing import Any

import yaml


ROOT = pathlib.Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "evidence" / "figures_tables.yml"
RESOURCES = ROOT / "evidence" / "resources.yml"
CLAIMS = ROOT / "evidence" / "claims.yml"
OUTPUT = ROOT / "notes" / "figure_table_readiness_board.html"


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def load_yaml(path: pathlib.Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Missing YAML file: {path}")
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a YAML mapping.")
    return data


def load_registry() -> tuple[dict[str, Any], set[str], set[str]]:
    data = load_yaml(REGISTRY)
    resource_data = load_yaml(RESOURCES)
    claim_data = load_yaml(CLAIMS)
    resource_ids = {r["id"] for r in resource_data.get("resources", [])}
    claim_ids = {c["id"] for c in claim_data.get("claims", [])}
    if "records" not in data:
        raise ValueError("Figure/table registry must contain a 'records' list.")
    required = [
        "id",
        "title",
        "kind",
        "status",
        "manuscript_section",
        "main_or_supplement",
        "source_resource_ids",
        "supported_claim_ids",
        "rendered_path",
        "source_path",
        "build_command",
        "readiness_note",
        "action_needed",
        "owner_role",
        "caveats",
    ]
    schema = data.get("schema", {})
    allowed = {
        "kind": set(schema.get("allowed_kinds", [])),
        "status": set(schema.get("allowed_statuses", [])),
        "main_or_supplement": set(schema.get("allowed_main_or_supplement", [])),
    }
    seen: set[str] = set()
    for record in data["records"]:
        missing = [key for key in required if key not in record]
        if missing:
            raise ValueError(f"Record {record.get('id', '<unknown>')} missing fields: {missing}")
        if record["id"] in seen:
            raise ValueError(f"Duplicate record id: {record['id']}")
        seen.add(record["id"])
        for field, choices in allowed.items():
            if choices and record[field] not in choices:
                raise ValueError(f"Record {record['id']} has unsupported {field}: {record[field]}")
        if not isinstance(record["source_resource_ids"], list):
            raise ValueError(f"Record {record['id']} source_resource_ids must be a list.")
        if not isinstance(record["supported_claim_ids"], list):
            raise ValueError(f"Record {record['id']} supported_claim_ids must be a list.")
        for resource_id in record["source_resource_ids"]:
            if resource_id not in resource_ids:
                raise ValueError(f"Record {record['id']} references unknown resource {resource_id}")
        for claim_id in record["supported_claim_ids"]:
            if claim_id not in claim_ids:
                raise ValueError(f"Record {record['id']} references unknown claim {claim_id}")
    return data, resource_ids, claim_ids


def class_token(value: str) -> str:
    return value.lower().replace(" ", "-").replace("/", "-")


def path_exists(path_value: str | None) -> bool:
    return bool(path_value) and pathlib.Path(path_value).exists()


def missing_reasons(record: dict[str, Any]) -> list[str]:
    reasons: list[str] = []
    rendered_path = record.get("rendered_path")
    source_path = record.get("source_path")
    if not rendered_path:
        reasons.append("no rendered path")
    elif not pathlib.Path(rendered_path).exists():
        reasons.append("rendered path missing")
    if not source_path:
        reasons.append("no source path")
    elif not pathlib.Path(source_path).exists():
        reasons.append("source path missing")
    if not record.get("supported_claim_ids"):
        reasons.append("no claim IDs")
    return reasons


def render_options(values: list[str]) -> str:
    return "\n".join(f"<option value=\"{esc(value)}\">{esc(value)}</option>" for value in values)


def render_summary(records: list[dict[str, Any]]) -> str:
    kind_counts = Counter(r["kind"] for r in records)
    ready_main = sum(
        1 for r in records
        if r["main_or_supplement"] == "main" and r["status"] in {"ready", "usable_with_edits"}
    )
    needs_work = sum(1 for r in records if r["status"] in {"needs_redraw", "needs_distillation"})
    tiles = [
        ("Figures", kind_counts.get("figure", 0)),
        ("Tables", kind_counts.get("table", 0)),
        ("Interactives", kind_counts.get("interactive", 0)),
        ("Planned/Missing", kind_counts.get("planned", 0)),
        ("Ready Main", ready_main),
        ("Needs Redraw/Distill", needs_work),
    ]
    return f"""
    <div class="summary-grid">
      {''.join(f'<div class="summary-tile"><strong>{value}</strong><span>{esc(label)}</span></div>' for label, value in tiles)}
    </div>
    """


def list_items(items: list[str]) -> str:
    return "".join(f"<li>{esc(item)}</li>" for item in items)


def inline_refs(values: list[str]) -> str:
    return ", ".join(f"<span class=\"ref-pill\">{esc(value)}</span>" for value in values)


def render_cards(records: list[dict[str, Any]]) -> str:
    cards = []
    for record in records:
        miss = missing_reasons(record)
        miss_text = "; ".join(miss)
        cards.append(
            f"""
            <article class="asset-card"
              data-kind="{esc(record['kind'])}"
              data-status="{esc(record['status'])}"
              data-section="{esc(record['manuscript_section'])}"
              data-placement="{esc(record['main_or_supplement'])}"
              data-owner="{esc(record['owner_role'])}"
              data-resources="{esc(','.join(record['source_resource_ids']))}"
              data-claims="{esc(','.join(record['supported_claim_ids']))}"
              data-missing="{esc(miss_text)}"
              data-search="{esc(json.dumps(record, sort_keys=True))}">
              <div class="asset-topline">
                <span class="asset-id">{esc(record['id'])}</span>
                <span class="pill kind">{esc(record['kind'])}</span>
                <span class="pill status status-{class_token(record['status'])}">{esc(record['status'])}</span>
                <span class="pill placement">{esc(record['main_or_supplement'])}</span>
              </div>
              <h3>{esc(record['title'])}</h3>
              <dl>
                <dt>Section</dt><dd>{esc(record['manuscript_section'])}</dd>
                <dt>Owner</dt><dd>{esc(record['owner_role'])}</dd>
                <dt>Resources</dt><dd>{inline_refs(record['source_resource_ids'])}</dd>
                <dt>Claims</dt><dd>{inline_refs(record['supported_claim_ids'])}</dd>
                <dt>Rendered</dt><dd><code>{esc(record.get('rendered_path') or '')}</code></dd>
                <dt>Source</dt><dd><code>{esc(record.get('source_path') or '')}</code></dd>
                <dt>Build</dt><dd><code>{esc(record.get('build_command') or '')}</code></dd>
                <dt>Readiness</dt><dd>{esc(record['readiness_note'])}</dd>
                <dt>Action</dt><dd>{esc(record['action_needed'])}</dd>
                <dt>Missing</dt><dd>{esc(miss_text or 'none')}</dd>
              </dl>
              <div>
                <h4>Caveats</h4>
                <ul>{list_items(record['caveats'])}</ul>
              </div>
            </article>
            """
        )
    return "\n".join(cards)


def render_section_matrix(records: list[dict[str, Any]]) -> str:
    matrix: dict[str, Counter] = defaultdict(Counter)
    buckets = ["ready", "usable_with_edits", "needs_redraw", "needs_distillation", "needed", "supplement_only", "background_only", "deprecated"]
    for record in records:
        matrix[record["manuscript_section"]][record["status"]] += 1
    rows = []
    for section, counts in sorted(matrix.items()):
        cells = "".join(f"<td>{counts.get(bucket, 0)}</td>" for bucket in buckets)
        rows.append(f"<tr><td>{esc(section)}</td>{cells}<td>{sum(counts.values())}</td></tr>")
    headers = "".join(f"<th>{esc(bucket)}</th>" for bucket in buckets)
    return f"""
    <table>
      <thead><tr><th>Manuscript section</th>{headers}<th>Total</th></tr></thead>
      <tbody>{''.join(rows)}</tbody>
    </table>
    """


def render_missing(records: list[dict[str, Any]]) -> str:
    rows = []
    for record in records:
        reasons = missing_reasons(record)
        if not reasons:
            continue
        rows.append(
            f"""
            <tr>
              <td>{esc(record['id'])}</td>
              <td>{esc(record['title'])}</td>
              <td>{esc(record['kind'])}</td>
              <td>{esc('; '.join(reasons))}</td>
              <td>{esc(record['action_needed'])}</td>
            </tr>
            """
        )
    if not rows:
        rows.append("<tr><td colspan=\"5\">No missing asset issues detected.</td></tr>")
    return f"""
    <table>
      <thead><tr><th>ID</th><th>Title</th><th>Kind</th><th>Issue</th><th>Action</th></tr></thead>
      <tbody>{''.join(rows)}</tbody>
    </table>
    """


def export_markdown(records: list[dict[str, Any]]) -> str:
    lines = ["# SIMODS Figure/Table Readiness Export", ""]
    for record in records:
        lines.extend(
            [
                f"## {record['id']}: {record['title']}",
                f"- Kind: {record['kind']}",
                f"- Status: {record['status']}",
                f"- Manuscript section: {record['manuscript_section']}",
                f"- Main/supplement: {record['main_or_supplement']}",
                f"- Source resources: {', '.join(record['source_resource_ids'])}",
                f"- Supported claims: {', '.join(record['supported_claim_ids'])}",
                f"- Rendered path: {record.get('rendered_path') or ''}",
                f"- Source path: {record.get('source_path') or ''}",
                f"- Build command: {record.get('build_command') or ''}",
                f"- Readiness note: {record['readiness_note']}",
                f"- Action needed: {record['action_needed']}",
                f"- Owner role: {record['owner_role']}",
                f"- Missing flags: {', '.join(missing_reasons(record)) or 'none'}",
                "- Caveats:",
            ]
        )
        lines.extend(f"  - {item}" for item in record["caveats"])
        lines.append("")
    return "\n".join(lines)


def build_html() -> str:
    data, _resource_ids, _claim_ids = load_registry()
    metadata = data.get("metadata", {})
    records = data["records"]
    now = _dt.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")
    kinds = sorted({r["kind"] for r in records})
    statuses = sorted({r["status"] for r in records})
    sections = sorted({r["manuscript_section"] for r in records})
    placements = sorted({r["main_or_supplement"] for r in records})
    owners = sorted({r["owner_role"] for r in records})
    resource_ids = sorted({rid for r in records for rid in r["source_resource_ids"]})
    claim_ids = sorted({cid for r in records for cid in r["supported_claim_ids"]})
    export_text = esc(export_markdown(records))

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>SIMODS Figure/Table Readiness Board</title>
  <style>
    :root {{
      --bg: #f7f6f1;
      --ink: #1d252c;
      --muted: #5c6670;
      --panel: #fff;
      --line: #d7d1c3;
      --accent: #0f6d7a;
      --accent-soft: #e2f2f3;
      --ready: #1f7356;
      --work: #9a6a00;
      --need: #a33f2f;
      --quiet: #516a7a;
      --shadow: 0 10px 30px rgba(20, 30, 40, 0.08);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      background: var(--bg);
      color: var(--ink);
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      line-height: 1.5;
    }}
    header {{
      padding: 28px clamp(20px, 5vw, 56px);
      color: #fff;
      background: #17232b;
    }}
    header h1 {{ margin: 0 0 8px; font-size: clamp(2rem, 4vw, 3.4rem); line-height: 1.05; letter-spacing: 0; }}
    header p {{ max-width: 980px; margin: 0; color: #dbe5e7; }}
    .meta-strip {{ display: flex; flex-wrap: wrap; gap: 10px; margin-top: 16px; }}
    .meta-strip span {{ border: 1px solid rgba(255,255,255,0.26); border-radius: 999px; padding: 6px 10px; font-size: 0.88rem; }}
    nav {{
      position: sticky;
      top: 0;
      z-index: 10;
      display: flex;
      gap: 8px;
      overflow-x: auto;
      padding: 10px clamp(20px, 5vw, 56px);
      border-bottom: 1px solid var(--line);
      background: rgba(247, 246, 241, 0.96);
      backdrop-filter: blur(10px);
    }}
    nav a {{ flex: 0 0 auto; text-decoration: none; color: var(--ink); background: #fff; border: 1px solid var(--line); border-radius: 999px; padding: 7px 12px; font-size: 0.92rem; }}
    main {{ width: min(1400px, calc(100% - 40px)); margin: 24px auto 56px; }}
    section {{ margin: 24px 0; }}
    .section-title {{ display: flex; justify-content: space-between; align-items: baseline; gap: 20px; margin-bottom: 12px; }}
    h2 {{ margin: 0; font-size: 1.55rem; }}
    .section-title p {{ max-width: 780px; margin: 0; color: var(--muted); }}
    .summary-grid {{
      display: grid;
      grid-template-columns: repeat(6, minmax(0, 1fr));
      gap: 12px;
      margin-bottom: 12px;
    }}
    .summary-tile {{ padding: 18px; border: 1px solid var(--line); border-radius: 8px; background: var(--panel); box-shadow: var(--shadow); }}
    .summary-tile strong {{ display: block; font-size: 2rem; line-height: 1; color: var(--accent); }}
    .summary-tile span {{ color: var(--muted); }}
    .filters {{
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 10px;
      padding: 14px;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: var(--panel);
      box-shadow: var(--shadow);
      margin-bottom: 14px;
    }}
    .filters label {{ display: grid; gap: 4px; font-size: 0.84rem; color: var(--muted); font-weight: 700; }}
    select, input, button {{
      border: 1px solid var(--line);
      border-radius: 6px;
      background: #fff;
      padding: 8px 10px;
      font: inherit;
      color: var(--ink);
    }}
    .asset-grid {{ display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 12px; }}
    .asset-card {{ padding: 16px; border: 1px solid var(--line); border-radius: 8px; background: var(--panel); box-shadow: var(--shadow); }}
    .asset-topline {{ display: flex; flex-wrap: wrap; gap: 6px; align-items: center; margin-bottom: 10px; }}
    .asset-id {{ display: inline-grid; place-items: center; min-width: 44px; padding: 4px 8px; border-radius: 999px; color: #fff; background: #17232b; font-weight: 800; font-size: 0.8rem; }}
    .pill, .ref-pill {{ display: inline-block; border-radius: 999px; padding: 4px 8px; background: #ece8dc; color: var(--ink); font-size: 0.78rem; font-weight: 700; white-space: nowrap; margin: 0 3px 3px 0; }}
    .kind {{ background: var(--accent-soft); color: #064951; }}
    .placement {{ background: #e8edf8; color: #24415f; }}
    .status-ready, .status-usable_with_edits {{ background: var(--ready); color: #fff; }}
    .status-needs_redraw, .status-needs_distillation {{ background: var(--work); color: #fff; }}
    .status-needed {{ background: var(--need); color: #fff; }}
    .status-supplement_only, .status-background_only, .status-deprecated {{ background: var(--quiet); color: #fff; }}
    .asset-card h3 {{ margin: 0 0 10px; font-size: 1.08rem; line-height: 1.35; }}
    dl {{ display: grid; grid-template-columns: 96px minmax(0, 1fr); gap: 5px 10px; margin: 0 0 12px; font-size: 0.92rem; }}
    dt {{ color: var(--muted); font-weight: 800; }}
    dd {{ margin: 0; }}
    code {{ display: block; padding: 6px 8px; border-radius: 5px; background: #eeeae0; overflow-wrap: anywhere; font-size: 0.84rem; }}
    h4 {{ margin: 0 0 5px; font-size: 0.9rem; color: var(--accent); }}
    ul {{ margin: 0; padding-left: 18px; }}
    table {{ width: 100%; border-collapse: collapse; border: 1px solid var(--line); border-radius: 8px; overflow: hidden; background: #fff; box-shadow: var(--shadow); }}
    th, td {{ padding: 10px 12px; border-bottom: 1px solid var(--line); text-align: left; vertical-align: top; }}
    th {{ background: #ebe6d9; }}
    textarea {{
      width: 100%;
      min-height: 300px;
      padding: 12px;
      border: 1px solid var(--line);
      border-radius: 8px;
      font: 0.9rem ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      background: #fff;
      color: var(--ink);
    }}
    footer {{ width: min(1400px, calc(100% - 40px)); margin: 0 auto 42px; color: var(--muted); font-size: 0.92rem; }}
    @media (max-width: 1100px) {{ .summary-grid, .filters, .asset-grid {{ grid-template-columns: 1fr; }} }}
  </style>
</head>
<body>
  <header>
    <h1>SIMODS Figure/Table Readiness Board</h1>
    <p>Source-driven inventory of existing, candidate, supplemental, background, and missing manuscript visuals and tables.</p>
    <div class="meta-strip">
      <span>Generated: {esc(now)}</span>
      <span>Registry: {esc(REGISTRY)}</span>
      <span>Registry updated: {esc(metadata.get('updated', 'unknown'))}</span>
      <span>Records: {len(records)}</span>
    </div>
  </header>
  <nav>
    <a href="#summary">Summary</a>
    <a href="#assets">Assets</a>
    <a href="#matrix">Section Matrix</a>
    <a href="#missing">Missing Assets</a>
    <a href="#export">Export</a>
  </nav>
  <main>
    <section id="summary">
      <div class="section-title">
        <h2>Summary</h2>
        <p>Counts separate existing figures/tables from planned missing assets and highlight what can move toward main text.</p>
      </div>
      {render_summary(records)}
    </section>
    <section id="assets">
      <div class="section-title">
        <h2>Assets</h2>
        <p>Filter by kind, status, section, placement, owner, resource ID, claim ID, or free text.</p>
      </div>
      <div class="filters">
        <label>Kind
          <select id="filter-kind"><option value="">All</option>{render_options(kinds)}</select>
        </label>
        <label>Status
          <select id="filter-status"><option value="">All</option>{render_options(statuses)}</select>
        </label>
        <label>Section
          <select id="filter-section"><option value="">All</option>{render_options(sections)}</select>
        </label>
        <label>Main/Supplement
          <select id="filter-placement"><option value="">All</option>{render_options(placements)}</select>
        </label>
        <label>Owner
          <select id="filter-owner"><option value="">All</option>{render_options(owners)}</select>
        </label>
        <label>Resource ID
          <select id="filter-resource"><option value="">All</option>{render_options(resource_ids)}</select>
        </label>
        <label>Claim ID
          <select id="filter-claim"><option value="">All</option>{render_options(claim_ids)}</select>
        </label>
        <label>Search
          <input id="filter-search" type="search" placeholder="title, action, path">
        </label>
      </div>
      <div class="asset-grid" id="asset-grid">
        {render_cards(records)}
      </div>
    </section>
    <section id="matrix">
      <div class="section-title">
        <h2>Section Coverage Matrix</h2>
        <p>Shows visual/table readiness by manuscript section.</p>
      </div>
      {render_section_matrix(records)}
    </section>
    <section id="missing">
      <div class="section-title">
        <h2>Missing Assets</h2>
        <p>Records with no rendered path, no source path, missing files, or no claim IDs.</p>
      </div>
      {render_missing(records)}
    </section>
    <section id="export">
      <div class="section-title">
        <h2>Export For Report Agent Or Auditor</h2>
        <p>Markdown summary suitable for a drafting or audit handoff.</p>
      </div>
      <button id="copy-export">Copy export</button>
      <textarea id="export-text" spellcheck="false">{export_text}</textarea>
    </section>
  </main>
  <footer>
    Generated from <code>{esc(REGISTRY)}</code>. Edit the YAML registry, not the generated HTML.
  </footer>
  <script>
    const cards = Array.from(document.querySelectorAll('.asset-card'));
    const filters = {{
      kind: document.getElementById('filter-kind'),
      status: document.getElementById('filter-status'),
      section: document.getElementById('filter-section'),
      placement: document.getElementById('filter-placement'),
      owner: document.getElementById('filter-owner'),
      resource: document.getElementById('filter-resource'),
      claim: document.getElementById('filter-claim'),
      search: document.getElementById('filter-search')
    }};

    function exactOrEmpty(card, field, value) {{
      if (!value) return true;
      return (card.dataset[field] || '') === value;
    }}

    function listContains(card, field, value) {{
      if (!value) return true;
      return (card.dataset[field] || '').split(',').map(x => x.trim()).includes(value);
    }}

    function applyFilters() {{
      const text = filters.search.value.trim().toLowerCase();
      cards.forEach((card) => {{
        const ok =
          exactOrEmpty(card, 'kind', filters.kind.value) &&
          exactOrEmpty(card, 'status', filters.status.value) &&
          exactOrEmpty(card, 'section', filters.section.value) &&
          exactOrEmpty(card, 'placement', filters.placement.value) &&
          exactOrEmpty(card, 'owner', filters.owner.value) &&
          listContains(card, 'resources', filters.resource.value) &&
          listContains(card, 'claims', filters.claim.value) &&
          (!text || (card.dataset.search || '').toLowerCase().includes(text) || card.textContent.toLowerCase().includes(text));
        card.style.display = ok ? '' : 'none';
      }});
    }}

    Object.values(filters).forEach((el) => {{
      el.addEventListener('input', applyFilters);
      el.addEventListener('change', applyFilters);
    }});

    document.getElementById('copy-export').addEventListener('click', async () => {{
      const text = document.getElementById('export-text').value;
      try {{
        await navigator.clipboard.writeText(text);
        document.getElementById('copy-export').textContent = 'Copied';
        setTimeout(() => document.getElementById('copy-export').textContent = 'Copy export', 1200);
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
    OUTPUT.write_text(build_html(), encoding="utf-8")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
