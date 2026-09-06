#!/bin/zsh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PAPER_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_DIR="$(git -C "$PAPER_DIR" rev-parse --show-toplevel 2>/dev/null || (cd "$PAPER_DIR/.." && pwd))"
MANUSCRIPT_DIR="$PAPER_DIR"
BUILD_DIR="$PAPER_DIR/build"
LATEXMK="${LATEXMK:-$(command -v latexmk || print /Library/TeX/texbin/latexmk)}"
INPUT_TEX="geodesic_mds.tex"
OUTPUT_PDF="$BUILD_DIR/geodesic_mds.pdf"
BUILD_INFO_TEX="$BUILD_DIR/manuscript_build_info.tex"
# The paper is canonical by repository-relative path, including clean worktrees.
if [[ ! -f "$PAPER_DIR/$INPUT_TEX" || ! -d "$PAPER_DIR/evidence" ]]; then
  print -u2 "Missing focused manuscript source or portable evidence at $PAPER_DIR"
  exit 1
fi

escape_for_tex() {
  python3 - "$1" <<'PY'
import sys

s = sys.argv[1]
repl = {
    "\\": r"\textbackslash{}",
    "_": r"\_",
    "&": r"\&",
    "%": r"\%",
    "#": r"\#",
    "$": r"\$",
    "{": r"\{",
    "}": r"\}",
}
print("".join(repl.get(ch, ch) for ch in s))
PY
}

mkdir -p "$BUILD_DIR"

GIT_VERSION="$(git -C "$REPO_DIR" describe --tags --always 2>/dev/null || print 'unversioned')"
GIT_BUILD_NUMBER="$(git -C "$REPO_DIR" rev-list --count HEAD 2>/dev/null || print 'NA')"
GIT_COMMIT="$(git -C "$REPO_DIR" rev-parse --short HEAD 2>/dev/null || print 'unknown')"
if [[ -n "$(git -C "$PAPER_DIR" status --porcelain --untracked-files=normal -- . ':!review-bundles' 2>/dev/null)" ]]; then
  GIT_VERSION="$GIT_VERSION-paper-dirty"
fi
BUILD_DATETIME="$(TZ=America/New_York date '+%Y-%m-%d %H:%M:%S %Z')"

cat > "$BUILD_INFO_TEX" <<EOF
\renewcommand{\manuscriptversion}{$(escape_for_tex "$GIT_VERSION")}
\renewcommand{\manuscriptbuildnumber}{$(escape_for_tex "$GIT_BUILD_NUMBER")}
\renewcommand{\manuscriptcommit}{$(escape_for_tex "$GIT_COMMIT")}
\renewcommand{\manuscriptbuilddatetime}{$(escape_for_tex "$BUILD_DATETIME")}
EOF

cd "$MANUSCRIPT_DIR"

"$LATEXMK" \
  -pdf \
  -interaction=nonstopmode \
  -halt-on-error \
  -file-line-error \
  -outdir="$BUILD_DIR" \
  "$INPUT_TEX"

# Record the exact active inputs and output used by the immutable packager.
python3 - "$PAPER_DIR" <<'PYBUILD'
from pathlib import Path
import hashlib, json, subprocess, sys
p = Path(sys.argv[1])
git = subprocess.run(['git', '-C', str(p), 'rev-parse', 'HEAD'], capture_output=True, text=True)
files = [p/'geodesic_mds.tex', p/'geodesic_mds.bib', p/'build/geodesic_mds.bbl',
         p/'build/geodesic_mds.pdf', p/'build/manuscript_build_info.tex']
files += sorted((p/'figures/focused').glob('*.pdf')) + sorted((p/'tables/focused').glob('*.tex'))
(p/'build/build-inputs.json').write_text(json.dumps({'commit': git.stdout.strip() or None,
    'sha256': {str(f.relative_to(p)): hashlib.sha256(f.read_bytes()).hexdigest() for f in files}}, indent=2)+'\n')
PYBUILD
print "$OUTPUT_PDF"
