#!/usr/bin/env python3
"""Write or verify the source/summary/figure/PDF checksum manifest."""
from pathlib import Path
import hashlib, sys
p=Path(__file__).resolve().parent
manifest=p/'summary'/'artifact-checksums.sha256'
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
if '--write' in sys.argv:
    files=sorted(f for f in p.rglob('*') if f.is_file() and f!=manifest and
                 'build' not in f.relative_to(p).parts and '__pycache__' not in f.relative_to(p).parts)
    manifest.write_text(''.join(f'{digest(f)}  {f.relative_to(p)}\n' for f in files))
for line in manifest.read_text().splitlines():
    expected,name=line.split('  ',1)
    assert digest(p/name)==expected, f'Artifact changed: {name}'
print('Source, summary, figure and PDF checksums passed.')
