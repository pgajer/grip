#!/usr/bin/env python3
"""Render internal-review screenshots for cited paper figures."""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import yaml


BASE = Path(__file__).resolve().parent
MANIFEST = BASE / "paper_figure_screenshots.yml"


def load_figures() -> list[dict]:
    data = yaml.safe_load(MANIFEST.read_text(encoding="utf-8"))
    return data["figures"]


def render_page(source_pdf: Path, page: int, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        prefix = Path(tmp) / "page"
        subprocess.run(
            [
                "pdftoppm",
                "-png",
                "-r",
                "140",
                "-f",
                str(page),
                "-l",
                str(page),
                str(source_pdf),
                str(prefix),
            ],
            check=True,
        )
        rendered = next(Path(tmp).glob("page-*.png"))
        output.write_bytes(rendered.read_bytes())


def main() -> None:
    rendered: set[Path] = set()
    for fig in load_figures():
        source_pdf = BASE / fig["source_pdf"]
        output = BASE / fig["image_path"]
        if output in rendered:
            continue
        if not source_pdf.exists():
            raise FileNotFoundError(f"Missing source PDF for {fig['id']}: {source_pdf}")
        render_page(source_pdf, int(fig["pdf_page"]), output)
        rendered.add(output)
    print(f"Rendered {len(rendered)} unique paper-figure page screenshots")


if __name__ == "__main__":
    main()
