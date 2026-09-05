#!/usr/bin/env python3
"""Hash the completed source/data/figure/PDF bundle, excluding build byproducts."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import hashlib,json,platform,subprocess
import numpy,scipy,pandas,matplotlib,numba
p=Path(__file__).resolve().parent
files=[f for f in p.rglob('*') if f.is_file() and not any(x in f.relative_to(p).parts for x in ['build','__pycache__']) and f.name not in ['artifacts.json','report-partial.pdf']]
manifest=dict(created=datetime.now(ZoneInfo('America/New_York')).isoformat(),git_head_at_build=subprocess.check_output(['git','rev-parse','HEAD'],cwd=p,text=True).strip(),
 git_tree_clean_at_build=not bool(subprocess.check_output(['git','status','--porcelain','--',str(p)],cwd=p,text=True).strip()),
 python=platform.python_version(),numpy=numpy.__version__,scipy=scipy.__version__,pandas=pandas.__version__,matplotlib=matplotlib.__version__,numba=numba.__version__,
 R=subprocess.check_output(['Rscript','-e','cat(R.version.string,"; grip",as.character(packageVersion("grip")),"; smacof",as.character(packageVersion("smacof")),"; igraph",as.character(packageVersion("igraph")))'],text=True).strip(),
 files=[dict(path=str(f.relative_to(p)),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in sorted(files)])
(p/'artifacts.json').write_text(json.dumps(manifest,indent=2)+'\n');print('Frozen',len(files),'artifacts')
