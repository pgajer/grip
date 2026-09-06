#!/usr/bin/env python3
"""Build polished H005 paper reviews."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


BASE = Path(__file__).resolve().parent
REVIEW_DIRS = sorted(path for path in BASE.glob("P??_*") if path.is_dir())
LOG_DIR = BASE / "build_logs"
HTML_DIR = BASE / "html"


def review_tex(review_dir: Path) -> Path | None:
    candidates = sorted(review_dir.glob("P??_review.tex"))
    return candidates[0] if candidates else None


def build_one(review_dir: Path) -> bool:
    tex = review_tex(review_dir)
    if tex is None:
        print(f"skip {review_dir.name}: no P##_review.tex")
        return True
    LOG_DIR.mkdir(exist_ok=True)
    HTML_DIR.mkdir(exist_ok=True)
    log_path = LOG_DIR / f"{review_dir.name}.log"
    with log_path.open("w", encoding="utf-8") as log:
        pdf = subprocess.run(
            ["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", tex.name],
            cwd=review_dir,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    if pdf.returncode != 0:
        print(f"FAIL pdf {review_dir.name}; see {log_path}")
        return False
    html_out = HTML_DIR / f"{tex.stem}.html"
    html = subprocess.run(
        [
            "pandoc",
            str(tex),
            "--from",
            "latex",
            "--to",
            "html5",
            "--standalone",
            "--mathjax",
            "--citeproc",
            "--bibliography",
            str((BASE / "../references.bib").resolve()),
            "--css",
            str((BASE / "../review_html.css").resolve()),
            "-o",
            str(html_out),
        ],
        cwd=review_dir,
    )
    if html.returncode != 0:
        print(f"FAIL html {review_dir.name}")
        return False
    print(f"built {review_dir.name}")
    return True


def main() -> int:
    wanted = set(sys.argv[1:])
    dirs = REVIEW_DIRS
    if wanted:
        dirs = [path for path in REVIEW_DIRS if path.name[:3] in wanted]
    ok = True
    for review_dir in dirs:
        ok = build_one(review_dir) and ok
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
