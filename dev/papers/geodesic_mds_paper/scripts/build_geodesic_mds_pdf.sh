#!/bin/zsh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PAPER_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_DIR="$(git -C "$PAPER_DIR" rev-parse --show-toplevel 2>/dev/null || (cd "$PAPER_DIR/../../.." && pwd))"
MANUSCRIPT_DIR="$PAPER_DIR/manuscript"
BUILD_DIR="$PAPER_DIR/build"
LATEXMK="${LATEXMK:-/Library/TeX/texbin/latexmk}"
INPUT_TEX="geodesic_mds.tex"
OUTPUT_PDF="$BUILD_DIR/geodesic_mds.pdf"
BUILD_INFO_TEX="$BUILD_DIR/manuscript_build_info.tex"

escape_for_tex() {
  python - "$1" <<'PY'
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

GIT_VERSION="$(git -C "$REPO_DIR" describe --tags --always --dirty 2>/dev/null || print 'unversioned')"
GIT_BUILD_NUMBER="$(git -C "$REPO_DIR" rev-list --count HEAD 2>/dev/null || print 'NA')"
GIT_COMMIT="$(git -C "$REPO_DIR" rev-parse --short HEAD 2>/dev/null || print 'unknown')"
BUILD_DATETIME="$(date '+%Y-%m-%d %H:%M:%S %Z')"

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

print "$OUTPUT_PDF"
