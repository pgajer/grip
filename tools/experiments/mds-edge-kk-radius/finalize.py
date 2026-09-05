#!/usr/bin/env python3
"""Finish aggregation and numerical/report gates after all declared jobs complete."""
from pathlib import Path
import os,subprocess,time
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[2];OUT=REPO/'output/mds-edge-kk-radius'
while True:
 files=list((OUT/'fits-v2').rglob('*-result.rds'));extra=sum('extra-result' in f.name for f in files);main=len(files)-extra
 if main==1216 and extra==32:break
 errors=[str(f) for f in (OUT/'logs').glob('*.log') if not f.name.startswith('pilot-') and 'Execution halted' in f.read_text()]
 if errors:raise RuntimeError(errors)
 print('WAIT',main,'/1216 primary graphs;',extra,'/32 optimizer graphs',flush=True);time.sleep(10)
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',OMP_NUM_THREADS='1')
commands=[['make','-C',str(HERE),'collect'],['Rscript',str(HERE/'verify-coordinates.R'),'--record'],['make','-C',str(HERE),'report'],['python3',str(HERE/'check-artifacts.py')]]
for cmd in commands:
 print('RUN',' '.join(cmd),flush=True);subprocess.run(cmd,cwd=REPO,env=env,check=True)
print('Numerical and build gates passed; visual PDF review remains required.',flush=True)
