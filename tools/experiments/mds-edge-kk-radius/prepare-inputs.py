#!/usr/bin/env python3
import argparse,concurrent.futures,os,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('mode',choices=['original','replicate','size']);args=p.parse_args()
here=Path(__file__).resolve().parent;out=here.parents[2]/'output/mds-edge-kk-radius';out.mkdir(exist_ok=True)
jobs=[(s,m,rep,n) for s in ['paraboloid','saddle'] for m in ['disk','surface_area'] for rep,n in ([(1,240)] if args.mode=='original' else [(2,240),(3,240)] if args.mode=='replicate' else [(1,480)])]
def run(job):
 s,m,rep,n=job;env=dict(os.environ,NUMBA_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',OMP_NUM_THREADS='1')
 with (out/f'inputs-{s}-{m}-{rep}-n{n}.log').open('w') as log:
  subprocess.run([sys.executable,str(here/'inputs.py'),'--surface',s,'--sampling',m,'--rep',str(rep),'--n',str(n)],stdout=log,stderr=subprocess.STDOUT,env=env,check=True)
 print('DONE',job,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(run,jobs))
