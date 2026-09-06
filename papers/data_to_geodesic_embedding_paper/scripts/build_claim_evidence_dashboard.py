#!/usr/bin/env python3
"""Build the SIMODS claim-to-evidence dashboard from a YAML registry."""

from __future__ import annotations

import datetime as _dt
import html
import json
import pathlib
from collections import Counter, defaultdict
from typing import Any

import yaml


ROOT = pathlib.Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "evidence" / "claims.yml"
OUTPUT = ROOT / "notes" / "claim_evidence_dashboard.html"


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def load_registry() -> dict[str, Any]:
    if not REGISTRY.exists():
        raise FileNotFoundError(f"Missing evidence registry: {REGISTRY}")
    with REGISTRY.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError("Evidence registry must be a YAML mapping.")
    if "claims" not in data or "sources" not in data:
        raise ValueError("Evidence registry requires 'claims' and 'sources'.")
    for claim in data["claims"]:
        missing = [
            key
            for key in [
                "id",
                "claim",
                "status",
                "evidence_strength",
                "manuscript_section",
                "needed_for",
                "source_ids",
                "risk",
                "next_action",
            ]
            if key not in claim
        ]
        if missing:
            raise ValueError(f"Claim {claim.get('id', '<unknown>')} missing fields: {missing}")
        for source_id in claim["source_ids"]:
            if source_id not in data["sources"]:
                raise ValueError(f"Claim {claim['id']} references unknown source {source_id}")
    return data


def status_class(status: str) -> str:
    return "status-" + status.lower().replace(" ", "-").replace("/", "-")


def render_summary(claims: list[dict[str, Any]]) -> str:
    status_counts = Counter(claim["status"] for claim in claims)
    needed_counts = Counter(claim["needed_for"] for claim in claims)
    audit_count = sum(1 for claim in claims if "audit" in claim["status"] or "audit" in claim["risk"].lower())
    source_count = len({sid for claim in claims for sid in claim["source_ids"]})
    tiles = [
        ("Claims", len(claims)),
        ("Sources Used", source_count),
        ("Needs Audit Signals", audit_count),
        ("Main Text Claims", needed_counts.get("main_text", 0)),
    ]
    status_html = "".join(
        f"<span class=\"mini-pill {status_class(status)}\">{esc(status)}: {count}</span>"
        for status, count in sorted(status_counts.items())
    )
    return f"""
    <div class="summary-grid">
      {''.join(f'<div class="summary-tile"><strong>{esc(value)}</strong><span>{esc(label)}</span></div>' for label, value in tiles)}
    </div>
    <div class="status-row">{status_html}</div>
    """


def render_filter_options(values: list[str]) -> str:
    return "\n".join(f"<option value=\"{esc(value)}\">{esc(value)}</option>" for value in values)


def render_claims(claims: list[dict[str, Any]], sources: dict[str, Any]) -> str:
    rows = []
    for claim in claims:
        source_refs = ", ".join(f"[{sid}]" for sid in claim["source_ids"])
        source_titles = "; ".join(f"{sid}: {sources[sid]['label']}" for sid in claim["source_ids"])
        rows.append(
            f"""
            <article class="claim-card"
              data-status="{esc(claim['status'])}"
              data-section="{esc(claim['manuscript_section'])}"
              data-needed="{esc(claim['needed_for'])}"
              data-strength="{esc(claim['evidence_strength'])}"
              data-search="{esc(json.dumps(claim, sort_keys=True))}">
              <div class="claim-topline">
                <span class="claim-id">{esc(claim['id'])}</span>
                <span class="pill {status_class(claim['status'])}">{esc(claim['status'])}</span>
                <span class="pill strength">{esc(claim['evidence_strength'])}</span>
                <span class="pill needed">{esc(claim['needed_for'])}</span>
              </div>
              <h3>{esc(claim['claim'])}</h3>
              <dl>
                <dt>Section</dt><dd>{esc(claim['manuscript_section'])}</dd>
                <dt>Sources</dt><dd class="source-ref" title="{esc(source_titles)}">{esc(source_refs)}</dd>
                <dt>Risk</dt><dd>{esc(claim['risk'])}</dd>
                <dt>Next action</dt><dd>{esc(claim['next_action'])}</dd>
                <dt>Notes</dt><dd>{esc(claim.get('notes', ''))}</dd>
              </dl>
            </article>
            """
        )
    return "\n".join(rows)


def render_sources(sources: dict[str, Any]) -> str:
    items = []
    for source_id, source in sources.items():
        path = source.get("path", "")
        items.append(
            f"""
            <li>
              <span class="source-id">[{esc(source_id)}]</span>
              <div>
                <strong>{esc(source.get('label', source_id))}</strong>
                <em>{esc(source.get('type', 'source'))}</em>
                <code>{esc(path)}</code>
              </div>
            </li>
            """
        )
    return "\n".join(items)


def render_section_matrix(claims: list[dict[str, Any]]) -> str:
    matrix: dict[str, Counter] = defaultdict(Counter)
    for claim in claims:
        matrix[claim["manuscript_section"]][claim["status"]] += 1
    statuses = sorted({claim["status"] for claim in claims})
    rows = []
    for section, counts in sorted(matrix.items()):
        cells = "".join(f"<td>{counts.get(status, 0)}</td>" for status in statuses)
        rows.append(f"<tr><td>{esc(section)}</td>{cells}<td>{sum(counts.values())}</td></tr>")
    headers = "".join(f"<th>{esc(status)}</th>" for status in statuses)
    return f"""
    <table class="matrix">
      <thead><tr><th>Manuscript section</th>{headers}<th>Total</th></tr></thead>
      <tbody>{''.join(rows)}</tbody>
    </table>
    """


def export_markdown(claims: list[dict[str, Any]], sources: dict[str, Any]) -> str:
    lines = ["# SIMODS Claim-to-Evidence Export", ""]
    lines.append("## Claims")
    for claim in claims:
        lines.extend(
            [
                "",
                f"### {claim['id']}: {claim['claim']}",
                f"- Status: {claim['status']}",
                f"- Evidence strength: {claim['evidence_strength']}",
                f"- Manuscript section: {claim['manuscript_section']}",
                f"- Needed for: {claim['needed_for']}",
                f"- Sources: {', '.join(f'[{sid}]' for sid in claim['source_ids'])}",
                f"- Risk: {claim['risk']}",
                f"- Next action: {claim['next_action']}",
                f"- Notes: {claim.get('notes', '')}",
            ]
        )
    lines.extend(["", "## Source References"])
    for source_id, source in sources.items():
        lines.append(f"- [{source_id}] {source.get('label', source_id)} ({source.get('type', 'source')}): {source.get('path', '')}")
    return "\n".join(lines)


def build_html() -> str:
    data = load_registry()
    metadata = data.get("metadata", {})
    claims = data["claims"]
    sources = data["sources"]
    statuses = sorted({claim["status"] for claim in claims})
    sections = sorted({claim["manuscript_section"] for claim in claims})
    needed = sorted({claim["needed_for"] for claim in claims})
    strengths = sorted({claim["evidence_strength"] for claim in claims})
    now = _dt.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")
    export_text = esc(export_markdown(claims, sources))

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>SIMODS Claim-to-Evidence Dashboard</title>
  <style>
    :root {{
      --bg: #f6f6f1;
      --ink: #1d252c;
      --muted: #5c6670;
      --panel: #fff;
      --line: #d6d1c4;
      --accent: #0f6d7a;
      --accent-soft: #e2f2f3;
      --tested: #1f7356;
      --inferred: #6b5e9c;
      --proposed: #9a6a00;
      --audit: #a43a35;
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
    header h1 {{
      margin: 0 0 8px;
      font-size: clamp(2rem, 4vw, 3.6rem);
      line-height: 1.05;
      letter-spacing: 0;
    }}
    header p {{
      max-width: 980px;
      margin: 0;
      color: #dbe5e7;
    }}
    .meta-strip {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      margin-top: 16px;
    }}
    .meta-strip span {{
      border: 1px solid rgba(255,255,255,0.26);
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
      border-bottom: 1px solid var(--line);
      background: rgba(246, 246, 241, 0.96);
      backdrop-filter: blur(10px);
    }}
    nav a {{
      flex: 0 0 auto;
      text-decoration: none;
      color: var(--ink);
      background: #fff;
      border: 1px solid var(--line);
      border-radius: 999px;
      padding: 7px 12px;
      font-size: 0.92rem;
    }}
    main {{
      width: min(1340px, calc(100% - 40px));
      margin: 24px auto 56px;
    }}
    section {{ margin: 24px 0; }}
    .section-title {{
      display: flex;
      justify-content: space-between;
      align-items: baseline;
      gap: 20px;
      margin-bottom: 12px;
    }}
    h2 {{ margin: 0; font-size: 1.55rem; }}
    .section-title p {{
      max-width: 720px;
      margin: 0;
      color: var(--muted);
    }}
    .summary-grid {{
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 12px;
      margin-bottom: 12px;
    }}
    .summary-tile {{
      padding: 18px;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: var(--panel);
      box-shadow: var(--shadow);
    }}
    .summary-tile strong {{
      display: block;
      font-size: 2rem;
      line-height: 1;
      color: var(--accent);
    }}
    .summary-tile span {{ color: var(--muted); }}
    .status-row {{
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
    }}
    .mini-pill, .pill {{
      display: inline-block;
      border-radius: 999px;
      padding: 4px 8px;
      background: #ece8dc;
      color: var(--ink);
      font-size: 0.78rem;
      font-weight: 700;
      white-space: nowrap;
    }}
    .status-tested {{ color: #fff; background: var(--tested); }}
    .status-inferred {{ color: #fff; background: var(--inferred); }}
    .status-proposed {{ color: #fff; background: var(--proposed); }}
    .status-needs-audit {{ color: #fff; background: var(--audit); }}
    .strength {{ background: var(--accent-soft); color: #064951; }}
    .needed {{ background: #f0eadc; color: #60440b; }}
    .filters {{
      display: grid;
      grid-template-columns: repeat(5, minmax(0, 1fr));
      gap: 10px;
      padding: 14px;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: var(--panel);
      box-shadow: var(--shadow);
      margin-bottom: 14px;
    }}
    .filters label {{
      display: grid;
      gap: 4px;
      font-size: 0.84rem;
      color: var(--muted);
      font-weight: 700;
    }}
    select, input, button {{
      border: 1px solid var(--line);
      border-radius: 6px;
      background: #fff;
      padding: 8px 10px;
      font: inherit;
      color: var(--ink);
    }}
    .claim-grid {{
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 12px;
    }}
    .claim-card {{
      padding: 16px;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: var(--panel);
      box-shadow: var(--shadow);
    }}
    .claim-topline {{
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      align-items: center;
      margin-bottom: 10px;
    }}
    .claim-id {{
      display: inline-grid;
      place-items: center;
      min-width: 44px;
      padding: 4px 8px;
      border-radius: 999px;
      color: #fff;
      background: #17232b;
      font-weight: 800;
      font-size: 0.8rem;
    }}
    .claim-card h3 {{
      margin: 0 0 10px;
      font-size: 1.05rem;
      line-height: 1.35;
    }}
    dl {{
      display: grid;
      grid-template-columns: 92px minmax(0, 1fr);
      gap: 5px 10px;
      margin: 0;
      font-size: 0.92rem;
    }}
    dt {{ color: var(--muted); font-weight: 800; }}
    dd {{ margin: 0; }}
    .source-ref {{
      color: var(--accent);
      font-weight: 800;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      border: 1px solid var(--line);
      border-radius: 8px;
      overflow: hidden;
      background: #fff;
      box-shadow: var(--shadow);
    }}
    th, td {{
      padding: 10px 12px;
      border-bottom: 1px solid var(--line);
      text-align: left;
      vertical-align: top;
    }}
    th {{ background: #ebe6d9; }}
    .source-list {{
      display: grid;
      gap: 8px;
      margin: 0;
      padding: 0;
      list-style: none;
    }}
    .source-list li {{
      display: grid;
      grid-template-columns: 72px minmax(0, 1fr);
      gap: 10px;
      padding: 12px;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: #fff;
      box-shadow: var(--shadow);
    }}
    .source-id {{
      color: var(--accent);
      font-weight: 900;
    }}
    .source-list em {{
      display: inline-block;
      margin-left: 8px;
      color: var(--muted);
      font-style: normal;
      font-size: 0.88rem;
    }}
    code {{
      display: block;
      margin-top: 5px;
      padding: 6px 8px;
      border-radius: 5px;
      background: #eeeae0;
      overflow-wrap: anywhere;
      font-size: 0.86rem;
    }}
    textarea {{
      width: 100%;
      min-height: 250px;
      padding: 12px;
      border: 1px solid var(--line);
      border-radius: 8px;
      font: 0.9rem ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      background: #fff;
      color: var(--ink);
    }}
    footer {{
      width: min(1340px, calc(100% - 40px));
      margin: 0 auto 42px;
      color: var(--muted);
      font-size: 0.92rem;
    }}
    @media (max-width: 1100px) {{
      .summary-grid, .filters, .claim-grid {{
        grid-template-columns: 1fr;
      }}
    }}
  </style>
</head>
<body>
  <header>
    <h1>SIMODS Claim-to-Evidence Dashboard</h1>
    <p>Truth-control surface for manuscript claims: what is tested, inferred, proposed, missing, or ready for audit.</p>
    <div class="meta-strip">
      <span>Generated: {esc(now)}</span>
      <span>Registry: {esc(REGISTRY)}</span>
      <span>Registry updated: {esc(metadata.get('updated', 'unknown'))}</span>
      <span>Claims: {len(claims)}</span>
    </div>
  </header>
  <nav>
    <a href="#summary">Summary</a>
    <a href="#claims">Claims</a>
    <a href="#matrix">Section Matrix</a>
    <a href="#sources">Sources</a>
    <a href="#export">Export</a>
  </nav>
  <main>
    <section id="summary">
      <div class="section-title">
        <h2>Summary</h2>
        <p>Use this dashboard before drafting or auditing: unsupported claims should become tasks, not prose.</p>
      </div>
      {render_summary(claims)}
    </section>

    <section id="claims">
      <div class="section-title">
        <h2>Claims</h2>
        <p>Filter by status, manuscript section, evidence strength, paper role, or free text.</p>
      </div>
      <div class="filters">
        <label>Status
          <select id="filter-status"><option value="">All</option>{render_filter_options(statuses)}</select>
        </label>
        <label>Section
          <select id="filter-section"><option value="">All</option>{render_filter_options(sections)}</select>
        </label>
        <label>Evidence strength
          <select id="filter-strength"><option value="">All</option>{render_filter_options(strengths)}</select>
        </label>
        <label>Needed for
          <select id="filter-needed"><option value="">All</option>{render_filter_options(needed)}</select>
        </label>
        <label>Search
          <input id="filter-search" type="search" placeholder="claim, source, next action">
        </label>
      </div>
      <div class="claim-grid" id="claim-grid">
        {render_claims(claims, sources)}
      </div>
    </section>

    <section id="matrix">
      <div class="section-title">
        <h2>Section Matrix</h2>
        <p>Quick check for sections that lean too heavily on proposed or inferred claims.</p>
      </div>
      {render_section_matrix(claims)}
    </section>

    <section id="sources">
      <div class="section-title">
        <h2>Source References</h2>
        <p>Long paths stay here so claim cards remain readable.</p>
      </div>
      <ul class="source-list">
        {render_sources(sources)}
      </ul>
    </section>

    <section id="export">
      <div class="section-title">
        <h2>Export For Auditor Or Writer</h2>
        <p>Copy this summary into a handoff when a manuscript section needs drafting or audit.</p>
      </div>
      <button id="copy-export">Copy export</button>
      <textarea id="export-text" spellcheck="false">{export_text}</textarea>
    </section>
  </main>
  <footer>
    Generated from <code>{esc(REGISTRY)}</code>. Edit the YAML registry, not the generated HTML.
  </footer>
  <script>
    const cards = Array.from(document.querySelectorAll('.claim-card'));
    const filters = {{
      status: document.getElementById('filter-status'),
      section: document.getElementById('filter-section'),
      strength: document.getElementById('filter-strength'),
      needed: document.getElementById('filter-needed'),
      search: document.getElementById('filter-search')
    }};

    function matches(card, field, value) {{
      if (!value) return true;
      return (card.dataset[field] || '') === value;
    }}

    function applyFilters() {{
      const text = filters.search.value.trim().toLowerCase();
      cards.forEach((card) => {{
        const ok =
          matches(card, 'status', filters.status.value) &&
          matches(card, 'section', filters.section.value) &&
          matches(card, 'strength', filters.strength.value) &&
          matches(card, 'needed', filters.needed.value) &&
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
