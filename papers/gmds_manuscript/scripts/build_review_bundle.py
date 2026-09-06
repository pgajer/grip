#!/usr/bin/env python3
"""Package the active paper and frozen-evidence reproduction, using explicit inputs.

Run make pdf and make verify first. This creates author-review artifacts, not an
approved submission. Full experiment fitting remains in the pinned grip sources.
"""
import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import zipfile
from zoneinfo import ZoneInfo

P = Path(__file__).resolve().parents[1]
FIGURES = ['paths', 'shape-ratios', 'initializer-comparison',
           'graph-sensitivity', 'scale', 'saddle-spatial']
SCRIPTS = ['build_focused_figures.py', 'check_focused_numbers.py',
           'check_coordinate_components.py', 'export_coordinate_components.py',
           'check_retained_routes.R',
           'check_focused_mathematics.py', 'check_focused_citations.py',
           'verify_focused_evidence.py', 'export_focused_evidence.py',
           'build_geodesic_mds_pdf.sh', 'build_review_bundle.py']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git(*args):
    result = subprocess.run(['git', '-C', str(P), *args], text=True,
                            capture_output=True)
    return result.stdout.strip() if result.returncode == 0 else None


def archive(destination, root, files, extra):
    manifest = {}
    with zipfile.ZipFile(destination, 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for relative in sorted(files):
            source = P / relative
            if not source.is_file() or source.is_symlink():
                raise ValueError(f'Missing or symlinked input: {relative}')
            data = source.read_bytes()
            manifest[relative] = hashlib.sha256(data).hexdigest()
            z.writestr(f'{root}/{relative}', data)
        for relative, value in sorted(extra.items()):
            data = value.encode()
            manifest[relative] = hashlib.sha256(data).hexdigest()
            z.writestr(f'{root}/{relative}', data)
        z.writestr(f'{root}/bundle-files.json', json.dumps(manifest, indent=2)+'\n')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='New immutable review-package directory')
    args = parser.parse_args()
    commit = git('rev-parse', 'HEAD')
    dirty = git('status', '--porcelain', '--untracked-files=normal', '.', ':!review-bundles')
    if not commit or dirty:
        raise ValueError('Commit the active paper inputs before packaging (review-bundles excluded).')
    built = json.loads((P/'build/build-inputs.json').read_text())
    if built['commit'] != commit:
        raise ValueError('Rebuild the PDF at the current source commit before packaging.')
    for relative, expected in built['sha256'].items():
        if digest(P/relative) != expected:
            raise ValueError(f'Input or PDF changed since build: {relative}')
    name = datetime.now(ZoneInfo('America/New_York')).strftime('%Y-%m-%d_%H%M%S') + '-' + commit[:7]
    out = (args.output or P/'review-bundles'/name).resolve()
    out.mkdir(parents=True, exist_ok=False)
    core = ['geodesic_mds.tex', 'geodesic_mds.bib',
            'build/manuscript_build_info.tex']
    core += [f'figures/focused/{name}.pdf' for name in FIGURES]
    core += [f'tables/focused/{name}.tex' for name in ['initializers', 'radius', 'numbers', 'components', 'continuation']]
    # Minimal typesetting inputs exclude source reports and historical drafts.
    typesetting = archive(out/'typesetting-source.zip', 'geodesic-mds-typesetting',
        core, {
        'geodesic_mds.bbl': (P/'build/geodesic_mds.bbl').read_text(),
        'Makefile': 'pdf:\n\tlatexmk -pdf -interaction=nonstopmode -halt-on-error geodesic_mds.tex\n',
        'README.md': '''# Author-review typesetting source

Run `make pdf` with latexmk and a standard TeX installation. The build timestamp
is retained from the reviewed version. This archive contains only active
manuscript inputs and the generated bibliography; it is not an approved arXiv
upload. Author review, authorship/disclosure approval, and submission choices
remain pending. See the companion frozen-evidence reproduction archive for
figure generation and numerical checks.
'''
    })
    reproduction_files = core + ['Makefile', 'README.md', 'citation_verification.html']
    reproduction_files += [f'scripts/{name}' for name in SCRIPTS]
    reproduction_files += [str(path.relative_to(P)) for path in (P/'evidence').rglob('*')
                           if path.is_file()]
    reproduction = archive(out/'evidence-reproduction.zip', 'geodesic-mds-reproduction',
        reproduction_files, {})
    shutil.copy2(P/'build/geodesic_mds.pdf', out/'geodesic_mds.pdf')
    provenance = {
        'created_eastern': datetime.now(ZoneInfo('America/New_York')).isoformat(),
        'paper_commit': git('rev-parse', 'HEAD'),
        'paper_branch': git('branch', '--show-current'),
        'tracked_paper_changes': git('status', '--porcelain', '--untracked-files=no', '.', ':!review-bundles'),
        'status': 'author_review_draft_not_approved_for_submission',
        'scope': 'Frozen-evidence verification and manuscript/figure regeneration; no full experiment refit.',
        'files': {name: {'sha256': digest(out/name), 'bytes': (out/name).stat().st_size}
                  for name in ['geodesic_mds.pdf', 'typesetting-source.zip', 'evidence-reproduction.zip']},
        'typesetting_inputs': len(typesetting),
        'reproduction_inputs': len(reproduction)
    }
    (out/'artifact-manifest.json').write_text(json.dumps(provenance, indent=2)+'\n')
    (out/'README.md').write_text('''# Graph-geodesic methods paper: author-review package

Start with `geodesic_mds.pdf`. The focused paper combines the completed MDS
initializer and radius studies, with their experimental designs kept separate.

- `typesetting-source.zip`: minimal active TeX/BibTeX/bibliography, six figure
  PDFs, five generated table/macro files, and a Makefile. Extract and run
  `make pdf` with latexmk on PATH.
- `evidence-reproduction.zip`: active source, figure builders, validation
  scripts, citation-support HTML, frozen numerical exports and provenance.
  Requires Python with NumPy, pandas and Matplotlib, zsh, and latexmk. Extract,
  run `make verify`, then `make pdf`. Set LATEXMK if needed.
- `artifact-manifest.json`: paper commit, branch, input counts and SHA-256 hashes.
  Each source archive also has `bundle-files.json` with per-input hashes.

The reproduction archive regenerates the figures and checks selected claims
from saved evidence. It does not rerun the original MDS/geodesic/edge-KK fits.
Their pinned source paths and commits are in the evidence provenance records.
The original study documents are preserved as source records; their relative
links refer to the source study trees and are not a complete mirror here.

This is an author-review draft. Internal checks are complete; independent
scientific review informed this revision; approval by both authors, disclosure review, overlap review
with the software article, and submission/category/license choices remain.
The canonical source is grip/papers/gmds_manuscript. No submission has been performed.

Distribute this directory as one matched set. Its PDF, source archives and
manifest describe the same revision; older dated review folders are historical
artifacts and must not be substituted for any member of this set.
''')
    print(str(out))
    print(json.dumps(provenance, indent=2))


if __name__ == '__main__':
    main()
