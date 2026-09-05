#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
p=Path(__file__).resolve().parent;m=json.loads((p/'artifacts.json').read_text())
for f in m['files']:
 file=p/f['path'];assert file.is_file(),file
 assert file.stat().st_size==f['bytes'],file
 assert hashlib.sha256(file.read_bytes()).hexdigest()==f['sha256'],file
print('Artifact checksums passed:',len(m['files']))
