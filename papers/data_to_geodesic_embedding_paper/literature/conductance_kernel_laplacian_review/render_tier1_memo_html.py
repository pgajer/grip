#!/usr/bin/env python3
"""Render Tier 1 Markdown memos to HTML with formula-like code as LaTeX math."""

from __future__ import annotations

import re
import subprocess
import tempfile
from pathlib import Path

import yaml


BASE = Path(__file__).resolve().parent
MEMO_DIR = BASE / "paper_memos"
HTML_DIR = MEMO_DIR / "html"
CSS = "../../review_html.css"
FIGURE_MANIFEST = BASE / "paper_figure_screenshots.yml"

GREEK = {
    "alpha": r"\alpha",
    "beta": r"\beta",
    "delta": r"\delta",
    "epsilon": r"\epsilon",
    "lambda": r"\lambda",
    "mu": r"\mu",
    "pi": r"\pi",
    "rho": r"\rho",
    "sigma": r"\sigma",
    "theta": r"\theta",
}

FUNCTIONS = {
    "argmin": r"\argmin",
    "cos": r"\cos",
    "dim": r"\dim",
    "div": r"\operatorname{div}",
    "exp": r"\exp",
    "grad": r"\nabla",
    "int": r"\int",
    "log": r"\log",
    "max": r"\max",
    "min": r"\min",
    "sin": r"\sin",
    "sqrt": r"\sqrt",
    "sum": r"\sum",
    "tr": r"\operatorname{tr}",
}


def load_figure_manifest() -> dict[str, list[dict]]:
    if not FIGURE_MANIFEST.exists():
        return {}
    data = yaml.safe_load(FIGURE_MANIFEST.read_text(encoding="utf-8"))
    by_paper: dict[str, list[dict]] = {}
    for fig in data.get("figures", []):
        by_paper.setdefault(fig["paper_id"], []).append(fig)
    return by_paper


FIGURES_BY_PAPER = load_figure_manifest()


def looks_like_path_or_code(text: str) -> bool:
    lowered = text.lower()
    path_markers = [
        "/users/",
        "../",
        "./",
        "literature/",
        "paper_memos/",
        "figures/",
        "sources/",
    ]
    if "\\" in text or any(marker in lowered for marker in path_markers):
        return True
    if any(lowered.endswith(ext) for ext in [".md", ".pdf", ".yml", ".yaml", ".html", ".png", ".svg", ".tex"]):
        return True
    if any(token in text for token in ["fit.rdgraph", "refit.rdgraph", "source_manifest", "review_html", "gflow", "SIMODS"]):
        return True
    if re.fullmatch(r"P\d{2}|R\d{3}|H\d{3}|Eq\.?\s*\(?\d+\)?|Figure\s+\d+", text):
        return True
    return False


def paper_id_from_path(path: Path) -> str:
    match = re.match(r"(P\d{2})_", path.name)
    if not match:
        raise ValueError(f"Cannot determine paper ID from {path}")
    return match.group(1)


def figure_block(fig: dict) -> str:
    image_rel = "../../" + fig["image_path"]
    alt = f"{fig['paper_id']} {fig['label']} cropped figure panel"
    page_label = fig.get("display_page", fig["pdf_page"])
    caption = (
        f"Reproduced from {fig['citation']}, {fig['label']}, PDF p. {page_label} "
        f"(cropped from canonical reading copy), internal review only. {fig['relevance']}"
    )
    return (
        f'\n<figure class="paper-figure" id="{fig["id"]}">\n'
        f'  <img src="{image_rel}" alt="{alt}">\n'
        f"  <figcaption>{caption}</figcaption>\n"
        f"</figure>\n"
    )


def mentions_figure_label(line: str, label: str) -> bool:
    if re.search(rf"\b{re.escape(label)}\b", line):
        return True
    if label.startswith("Figure "):
        tail = label.removeprefix("Figure ")
        if re.search(rf"\bFigures\s+{re.escape(tail)}\s+and\b", line):
            return True
        if re.search(rf"\bFigures\b[^\n]*?\band\s+{re.escape(tail)}\b", line):
            return True
    return False


def inject_figure_screenshots(text: str, paper_id: str) -> str:
    figures = FIGURES_BY_PAPER.get(paper_id, [])
    if not figures:
        return text
    inserted: set[str] = set()
    in_figures_section = False
    output: list[str] = []
    for line in text.splitlines():
        stripped = line.strip()
        if stripped == "### Figures And Experiments":
            in_figures_section = True
        elif stripped.startswith("### ") and stripped != "### Figures And Experiments":
            in_figures_section = False
        output.append(line)
        if not in_figures_section:
            continue
        for fig in figures:
            if mentions_figure_label(line, fig["label"]) and fig["id"] not in inserted:
                output.append(figure_block(fig))
                inserted.add(fig["id"])
    return "\n".join(output) + "\n"


def looks_like_formula(text: str) -> bool:
    if looks_like_path_or_code(text):
        return False
    stripped = text.strip()
    if len(stripped) <= 1:
        return False
    strong_markers = ["=", "^", "_", "||", "<", ">", "+", "/", "sum", "int", "sqrt", "exp", "lambda", "alpha", "epsilon", "sigma", "rho", "Delta", "grad"]
    if any(marker in stripped for marker in strong_markers):
        return True
    if re.search(r"\b[A-Z]\s*-\s*[A-Z]\b", stripped):
        return True
    return False


def convert_subscripts(text: str) -> str:
    text = re.sub(r"\b([A-Za-z][A-Za-z0-9]*)_([A-Za-z0-9]+)\b", r"\1_{\2}", text)
    text = re.sub(r"\b([A-Za-z][A-Za-z0-9]*)_\{([^}]+)\}", r"\1_{\2}", text)
    return text


def convert_latex(text: str) -> str:
    s = text.strip()
    s = s.replace("~=", r"\approx")
    s = s.replace(">=", r"\ge ")
    s = s.replace("<=", r"\le ")
    s = s.replace("->", r"\to")
    s = s.replace("<-", r"\leftarrow")
    s = s.replace("...", r"\ldots")
    s = s.replace("*", r"\,")
    s = s.replace("||", r"\Vert ")
    s = re.sub(r"\bR\^\(([^)]+)\)", r"\\mathbb{R}^{\1}", s)
    s = re.sub(r"\bR\^([A-Za-z0-9]+)", r"\\mathbb{R}^{\1}", s)
    s = s.replace("mathcal L", r"\mathcal{L}")
    s = s.replace("Delta", r"\Delta")
    s = s.replace("dmu", r"d\mu")
    s = s.replace("dot{x}", r"\dot{x}")
    s = s.replace("dot{w}", r"\dot{w}")
    s = convert_subscripts(s)
    for name, repl in sorted(FUNCTIONS.items(), key=lambda item: -len(item[0])):
        s = re.sub(rf"\b{name}\b", lambda _match, replacement=repl: replacement, s)
    for name, repl in sorted(GREEK.items(), key=lambda item: -len(item[0])):
        s = re.sub(rf"\b{name}\b", lambda _match, replacement=repl: replacement, s)
    s = re.sub(r"\bif\b", r"\\text{if}", s)
    s = re.sub(r"\botherwise\b", r"\\text{otherwise}", s)
    return s


def align_line(line: str) -> str:
    s = convert_latex(line)
    if "&" not in s and "=" in s:
        s = s.replace("=", "&=", 1)
    return s


def convert_block(code: str) -> str:
    lines = [line.rstrip() for line in code.strip("\n").splitlines()]
    lines = [line for line in lines if line.strip()]
    if not lines:
        return ""
    body = " \\\\\n".join(align_line(line) for line in lines)
    env = "aligned" if any("=" in line for line in lines) else "gathered"
    return f"\\[\n\\begin{{{env}}}\n{body}\n\\end{{{env}}}\n\\]\n"


def preprocess_markdown(text: str) -> str:
    placeholders: list[str] = []

    def fence_repl(match: re.Match[str]) -> str:
        lang = (match.group(1) or "").strip().lower()
        code = match.group(2)
        token = f"@@FENCE_{len(placeholders)}@@"
        if lang == "mermaid" or not looks_like_formula(code):
            placeholders.append(match.group(0))
        else:
            placeholders.append(convert_block(code))
        return token

    text = re.sub(r"```([A-Za-z0-9_-]*)\n(.*?)```", fence_repl, text, flags=re.DOTALL)

    def inline_repl(match: re.Match[str]) -> str:
        code = match.group(1)
        if looks_like_formula(code):
            return r"\(" + convert_latex(code) + r"\)"
        return match.group(0)

    text = re.sub(r"`([^`\n]+)`", inline_repl, text)
    for i, replacement in enumerate(placeholders):
        text = text.replace(f"@@FENCE_{i}@@", replacement)
    return text


def title_for(path: Path) -> str:
    first = path.read_text(encoding="utf-8").splitlines()[0]
    return first.removeprefix("# ").strip()


def strip_duplicate_top_heading(text: str) -> str:
    lines = text.splitlines()
    if lines and lines[0].startswith("# "):
        return "\n".join(lines[1:]).lstrip() + "\n"
    return text


def render_one(path: Path) -> Path:
    HTML_DIR.mkdir(exist_ok=True)
    out = HTML_DIR / f"{path.stem}.html"
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp) / path.name
        source = path.read_text(encoding="utf-8")
        source = strip_duplicate_top_heading(source)
        source = inject_figure_screenshots(source, paper_id_from_path(path))
        source = source.replace("](../figures/", "](../../figures/")
        tmp_path.write_text(preprocess_markdown(source), encoding="utf-8")
        subprocess.run(
            [
                "pandoc",
                str(tmp_path),
                "--from",
                "markdown+tex_math_dollars+tex_math_single_backslash",
                "--to",
                "html5",
                "--standalone",
                "--mathjax",
                "--metadata",
                f"title={title_for(path)}",
                "--css",
                CSS,
                "-o",
                str(out),
            ],
            check=True,
            cwd=BASE,
        )
    return out


def write_index(outputs: list[Path]) -> None:
    rows = "\n".join(
        f'    <li><a href="{out.name}">{title_for(MEMO_DIR / (out.stem + ".md"))}</a></li>'
        for out in outputs
    )
    (HTML_DIR / "index.html").write_text(
        f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>H005 Tier 1 Review Memo HTML Index</title>
  <link rel="stylesheet" href="../../review_html.css">
</head>
<body>
  <h1>H005 Tier 1 Review Memo HTML Index</h1>
  <p>Standalone HTML versions of the audited Tier 1 review memos. Formula-like memo notation is rendered through MathJax.</p>
  <ul>
{rows}
  </ul>
</body>
</html>
""",
        encoding="utf-8",
    )


def main() -> None:
    outputs = [render_one(path) for path in sorted(MEMO_DIR.glob("P*.md"))]
    write_index(outputs)
    print(f"Rendered {len(outputs)} Tier 1 memo HTML pages to {HTML_DIR}")


if __name__ == "__main__":
    main()
