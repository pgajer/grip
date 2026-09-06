#!/usr/bin/env python3
"""Build the SIMODS manuscript resource dashboard from a YAML registry."""

from __future__ import annotations

import datetime as _dt
import html
import json
import pathlib
import urllib.parse
from collections import Counter, defaultdict
from typing import Any

import yaml


ROOT = pathlib.Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "evidence" / "resources.yml"
OUTPUT = ROOT / "notes" / "manuscript_resource_dashboard.html"


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def path_href(value: object) -> str:
    path = str(value)
    parsed = urllib.parse.urlparse(path)
    if parsed.scheme in {"http", "https", "file"}:
        return path
    candidate = pathlib.Path(path)
    if not candidate.is_absolute():
        candidate = ROOT / candidate
    return candidate.resolve().as_uri()


def linked_path(value: object) -> str:
    return f"<a class=\"path-link\" href=\"{esc(path_href(value))}\" target=\"_blank\" rel=\"noopener\"><code>{esc(value)}</code></a>"


def load_registry() -> dict[str, Any]:
    if not REGISTRY.exists():
        raise FileNotFoundError(f"Missing resource registry: {REGISTRY}")
    data = yaml.safe_load(REGISTRY.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or "resources" not in data:
        raise ValueError("Resource registry must contain a 'resources' list.")
    required = [
        "id",
        "title",
        "project",
        "type",
        "readiness",
        "priority",
        "manuscript_relevance",
        "resource_role",
        "path",
        "sections_supported",
        "reusable_items",
        "edits_needed",
        "caveats",
    ]
    for resource in data["resources"]:
        missing = [key for key in required if key not in resource]
        if missing:
            raise ValueError(f"Resource {resource.get('id', '<unknown>')} missing fields: {missing}")
    schema = data.get("schema", {})
    allowed_fields = {
        "readiness": schema.get("allowed_readiness_values"),
        "manuscript_relevance": schema.get("allowed_manuscript_relevance"),
        "resource_role": schema.get("allowed_resource_roles"),
    }
    for field, allowed in allowed_fields.items():
        if not allowed:
            continue
        allowed_set = set(allowed)
        for resource in data["resources"]:
            value = resource[field]
            if value not in allowed_set:
                raise ValueError(
                    f"Resource {resource.get('id', '<unknown>')} has unsupported {field}: {value}"
                )
    return data


def class_token(value: str) -> str:
    return value.lower().replace(" ", "-").replace("/", "-")


def render_options(values: list[str]) -> str:
    return "\n".join(f"<option value=\"{esc(value)}\">{esc(value)}</option>" for value in values)


def render_summary(resources: list[dict[str, Any]]) -> str:
    projects = len({r["project"] for r in resources})
    p1 = sum(1 for r in resources if r["priority"] == "P1")
    reusable = sum(1 for r in resources if r["readiness"] in {"usable with edits", "ready to cite internally", "source of planning truth", "active planning artifact", "active audit artifact"})
    main_text = sum(1 for r in resources if r["manuscript_relevance"] == "reusable for main text")
    canonical = sum(1 for r in resources if r["resource_role"] == "canonical")
    types = Counter(r["type"] for r in resources)
    type_pills = "".join(
        f"<span class=\"mini-pill\">{esc(kind)}: {count}</span>"
        for kind, count in sorted(types.items())
    )
    tiles = [
        ("Resources", len(resources)),
        ("Projects", projects),
        ("P1 Resources", p1),
        ("Reusable Now", reusable),
        ("Main Text", main_text),
        ("Canonical", canonical),
    ]
    return f"""
    <div class="summary-grid">
      {''.join(f'<div class="summary-tile"><strong>{value}</strong><span>{esc(label)}</span></div>' for label, value in tiles)}
    </div>
    <div class="type-row">{type_pills}</div>
    """


def list_items(items: list[str]) -> str:
    return "".join(f"<li>{esc(item)}</li>" for item in items)


def render_cards(resources: list[dict[str, Any]]) -> str:
    cards = []
    for resource in resources:
        sections = ", ".join(resource["sections_supported"])
        cards.append(
            f"""
            <article class="resource-card"
              data-project="{esc(resource['project'])}"
              data-type="{esc(resource['type'])}"
              data-readiness="{esc(resource['readiness'])}"
              data-relevance="{esc(resource['manuscript_relevance'])}"
              data-role="{esc(resource['resource_role'])}"
              data-priority="{esc(resource['priority'])}"
              data-sections="{esc(sections)}"
              data-search="{esc(json.dumps(resource, sort_keys=True))}">
              <div class="resource-topline">
                <span class="resource-id">{esc(resource['id'])}</span>
                <span class="pill priority-{class_token(resource['priority'])}">{esc(resource['priority'])}</span>
                <span class="pill readiness">{esc(resource['readiness'])}</span>
                <span class="pill relevance">{esc(resource['manuscript_relevance'])}</span>
                <span class="pill role">{esc(resource['resource_role'])}</span>
                <span class="pill type">{esc(resource['type'])}</span>
              </div>
              <h3>{esc(resource['title'])}</h3>
              <dl>
                <dt>Project</dt><dd>{esc(resource['project'])}</dd>
                <dt>Relevance</dt><dd>{esc(resource['manuscript_relevance'])}</dd>
                <dt>Role</dt><dd>{esc(resource['resource_role'])}</dd>
                <dt>Sections</dt><dd>{esc(sections)}</dd>
                <dt>Rendered</dt><dd>{linked_path(resource['path'])}</dd>
                <dt>Source</dt><dd><code>{esc(resource.get('source_path', ''))}</code></dd>
              </dl>
              <div class="resource-lists">
                <div>
                  <h4>Reusable Items</h4>
                  <ul>{list_items(resource['reusable_items'])}</ul>
                </div>
                <div>
                  <h4>Edits Needed</h4>
                  <ul>{list_items(resource['edits_needed'])}</ul>
                </div>
                <div>
                  <h4>Caveats</h4>
                  <ul>{list_items(resource['caveats'])}</ul>
                </div>
              </div>
            </article>
            """
        )
    return "\n".join(cards)


def render_section_matrix(resources: list[dict[str, Any]]) -> str:
    matrix: dict[str, Counter] = defaultdict(Counter)
    for resource in resources:
        for section in resource["sections_supported"]:
            matrix[section][resource["priority"]] += 1
    priorities = sorted({r["priority"] for r in resources})
    rows = []
    for section, counts in sorted(matrix.items()):
        cells = "".join(f"<td>{counts.get(priority, 0)}</td>" for priority in priorities)
        rows.append(f"<tr><td>{esc(section)}</td>{cells}<td>{sum(counts.values())}</td></tr>")
    headers = "".join(f"<th>{esc(priority)}</th>" for priority in priorities)
    return f"""
    <table>
      <thead><tr><th>Manuscript section</th>{headers}<th>Total resources</th></tr></thead>
      <tbody>{''.join(rows)}</tbody>
    </table>
    """


def export_markdown(resources: list[dict[str, Any]]) -> str:
    lines = ["# SIMODS Manuscript Resource Export", ""]
    for resource in resources:
        lines.extend(
            [
                f"## {resource['id']}: {resource['title']}",
                f"- Project: {resource['project']}",
                f"- Type: {resource['type']}",
                f"- Readiness: {resource['readiness']}",
                f"- Manuscript relevance: {resource['manuscript_relevance']}",
                f"- Resource role: {resource['resource_role']}",
                f"- Priority: {resource['priority']}",
                f"- Rendered path: {resource['path']}",
                f"- Source path: {resource.get('source_path', '')}",
                f"- Sections supported: {', '.join(resource['sections_supported'])}",
                "- Reusable items:",
            ]
        )
        lines.extend(f"  - {item}" for item in resource["reusable_items"])
        lines.append("- Edits needed:")
        lines.extend(f"  - {item}" for item in resource["edits_needed"])
        lines.append("- Caveats:")
        lines.extend(f"  - {item}" for item in resource["caveats"])
        lines.append("")
    return "\n".join(lines)


def build_html() -> str:
    data = load_registry()
    metadata = data.get("metadata", {})
    resources = data["resources"]
    now = _dt.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")
    projects = sorted({r["project"] for r in resources})
    types = sorted({r["type"] for r in resources})
    readiness = sorted({r["readiness"] for r in resources})
    relevance = sorted({r["manuscript_relevance"] for r in resources})
    roles = sorted({r["resource_role"] for r in resources})
    priorities = sorted({r["priority"] for r in resources})
    sections = sorted({section for r in resources for section in r["sections_supported"]})
    export_text = esc(export_markdown(resources))

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>SIMODS Manuscript Resource Dashboard</title>
  <style>
    :root {{
      --bg: #f7f6f1;
      --ink: #1d252c;
      --muted: #5c6670;
      --panel: #fff;
      --line: #d7d1c3;
      --accent: #0f6d7a;
      --accent-soft: #e2f2f3;
      --p1: #a33f2f;
      --p2: #9a6a00;
      --p3: #516a7a;
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
      background: rgba(247, 246, 241, 0.96);
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
      width: min(1360px, calc(100% - 40px));
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
    .type-row {{
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
    .priority-p1 {{ background: var(--p1); color: #fff; }}
    .priority-p2 {{ background: var(--p2); color: #fff; }}
    .priority-p3 {{ background: var(--p3); color: #fff; }}
    .readiness {{ background: var(--accent-soft); color: #064951; }}
    .relevance {{ background: #e8edf8; color: #24415f; }}
    .role {{ background: #e8f1e7; color: #284f2d; }}
    .type {{ background: #f0eadc; color: #60440b; }}
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
    .resource-grid {{
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 12px;
    }}
    .resource-card {{
      padding: 16px;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: var(--panel);
      box-shadow: var(--shadow);
    }}
    .resource-topline {{
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      align-items: center;
      margin-bottom: 10px;
    }}
    .resource-id {{
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
    .resource-card h3 {{
      margin: 0 0 10px;
      font-size: 1.08rem;
      line-height: 1.35;
    }}
    dl {{
      display: grid;
      grid-template-columns: 82px minmax(0, 1fr);
      gap: 5px 10px;
      margin: 0 0 12px;
      font-size: 0.92rem;
    }}
    dt {{ color: var(--muted); font-weight: 800; }}
    dd {{ margin: 0; }}
    code {{
      display: block;
      padding: 6px 8px;
      border-radius: 5px;
      background: #eeeae0;
      overflow-wrap: anywhere;
      font-size: 0.84rem;
    }}
    .path-link {{
      color: var(--accent);
      text-decoration: none;
    }}
    .path-link code {{
      color: inherit;
      border: 1px solid transparent;
    }}
    .path-link:hover code {{
      border-color: var(--accent);
      background: #e6f3f1;
    }}
    .resource-lists {{
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 10px;
    }}
    h4 {{
      margin: 0 0 5px;
      font-size: 0.9rem;
      color: var(--accent);
    }}
    ul {{
      margin: 0;
      padding-left: 18px;
    }}
    li {{ margin-bottom: 3px; }}
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
    textarea {{
      width: 100%;
      min-height: 280px;
      padding: 12px;
      border: 1px solid var(--line);
      border-radius: 8px;
      font: 0.9rem ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      background: #fff;
      color: var(--ink);
    }}
    footer {{
      width: min(1360px, calc(100% - 40px));
      margin: 0 auto 42px;
      color: var(--muted);
      font-size: 0.92rem;
    }}
    @media (max-width: 1100px) {{
      .summary-grid, .filters, .resource-grid, .resource-lists {{
        grid-template-columns: 1fr;
      }}
    }}
  </style>
</head>
<body>
  <header>
    <h1>SIMODS Manuscript Resource Dashboard</h1>
    <p>Inventory of reusable manuscripts, reports, dashboards, figures, notes, and experiment outputs for the SIMODS paper.</p>
    <div class="meta-strip">
      <span>Generated: {esc(now)}</span>
      <span>Registry: {esc(REGISTRY)}</span>
      <span>Registry updated: {esc(metadata.get('updated', 'unknown'))}</span>
      <span>Resources: {len(resources)}</span>
    </div>
  </header>
  <nav>
    <a href="#summary">Summary</a>
    <a href="#resources">Resources</a>
    <a href="#matrix">Section Matrix</a>
    <a href="#export">Export</a>
  </nav>
  <main>
    <section id="summary">
      <div class="section-title">
        <h2>Summary</h2>
        <p>Use this to decide what can be lifted, adapted, cited, distilled, or left as background.</p>
      </div>
      {render_summary(resources)}
    </section>

    <section id="resources">
      <div class="section-title">
        <h2>Resources</h2>
        <p>Filter by project, type, readiness, manuscript relevance, resource role, priority, manuscript section, or free text.</p>
      </div>
      <div class="filters">
        <label>Project
          <select id="filter-project"><option value="">All</option>{render_options(projects)}</select>
        </label>
        <label>Type
          <select id="filter-type"><option value="">All</option>{render_options(types)}</select>
        </label>
        <label>Readiness
          <select id="filter-readiness"><option value="">All</option>{render_options(readiness)}</select>
        </label>
        <label>Manuscript Relevance
          <select id="filter-relevance"><option value="">All</option>{render_options(relevance)}</select>
        </label>
        <label>Resource Role
          <select id="filter-role"><option value="">All</option>{render_options(roles)}</select>
        </label>
        <label>Priority
          <select id="filter-priority"><option value="">All</option>{render_options(priorities)}</select>
        </label>
        <label>Section
          <select id="filter-section"><option value="">All</option>{render_options(sections)}</select>
        </label>
        <label>Search
          <input id="filter-search" type="search" placeholder="title, reusable item, path">
        </label>
      </div>
      <div class="resource-grid" id="resource-grid">
        {render_cards(resources)}
      </div>
    </section>

    <section id="matrix">
      <div class="section-title">
        <h2>Section Coverage Matrix</h2>
        <p>Shows which manuscript sections already have candidate resources and their priority levels.</p>
      </div>
      {render_section_matrix(resources)}
    </section>

    <section id="export">
      <div class="section-title">
        <h2>Export For Writer Or Provenance/Librarian</h2>
        <p>Copy this into a handoff when asking an agent to draft a section or audit resource coverage.</p>
      </div>
      <button id="copy-export">Copy export</button>
      <textarea id="export-text" spellcheck="false">{esc(export_markdown(resources))}</textarea>
    </section>
  </main>
  <footer>
    Generated from <code>{esc(REGISTRY)}</code>. Edit the YAML registry, not the generated HTML.
  </footer>
  <script>
    const cards = Array.from(document.querySelectorAll('.resource-card'));
    const filters = {{
      project: document.getElementById('filter-project'),
      type: document.getElementById('filter-type'),
      readiness: document.getElementById('filter-readiness'),
      relevance: document.getElementById('filter-relevance'),
      role: document.getElementById('filter-role'),
      priority: document.getElementById('filter-priority'),
      section: document.getElementById('filter-section'),
      search: document.getElementById('filter-search')
    }};

    function exactOrEmpty(card, field, value) {{
      if (!value) return true;
      return (card.dataset[field] || '') === value;
    }}

    function sectionMatches(card, value) {{
      if (!value) return true;
      return (card.dataset.sections || '').split(',').map(x => x.trim()).includes(value);
    }}

    function applyFilters() {{
      const text = filters.search.value.trim().toLowerCase();
      cards.forEach((card) => {{
        const ok =
          exactOrEmpty(card, 'project', filters.project.value) &&
          exactOrEmpty(card, 'type', filters.type.value) &&
          exactOrEmpty(card, 'readiness', filters.readiness.value) &&
          exactOrEmpty(card, 'relevance', filters.relevance.value) &&
          exactOrEmpty(card, 'role', filters.role.value) &&
          exactOrEmpty(card, 'priority', filters.priority.value) &&
          sectionMatches(card, filters.section.value) &&
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
