#!/usr/bin/env python3
"""Read-only checkpoint census; cache files are atomically published."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import collections,json
p=Path(__file__).resolve().parents[3]/'output/mds-edge-kk-radius'
files=list((p/'fits-v2').rglob('*-result.rds'));main=[f for f in files if 'extra-result' not in f.name];extra=[f for f in files if 'extra-result' in f.name]
counts=collections.Counter(f.parent.name.split('-')[0]+' / '+f.parent.name.split('-')[1] for f in main)
errors=[]
for f in (p/'logs').glob('*.log'):
 if f.name.startswith('pilot-'):continue
 text=f.read_text()
 if 'Execution halted' in text or '\nError' in text:errors.append(str(f.relative_to(p)))
print(json.dumps(dict(time=datetime.now(ZoneInfo('America/New_York')).isoformat(),primary=len(main),primary_target=1216,optimizer_graphs=len(extra),optimizer_target=32,by_surface_measure=dict(counts),errors=errors),indent=2))
