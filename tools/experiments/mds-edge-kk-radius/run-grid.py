#!/usr/bin/env python3
"""Bounded process concurrency; individual fits are restartable, identity checked."""
import argparse,concurrent.futures,os,re,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('mode',choices=['pilot','main','extra','size']);p.add_argument('--workers',type=int,default=4);p.add_argument('--skip',default='');p.add_argument('--case-filter',default='');a=p.parse_args()
here=Path(__file__).resolve().parent;out=here.parents[2]/'output/mds-edge-kk-radius';(out/'logs').mkdir(exist_ok=True)
jobs=[f'{s}-{m}-rep{rep}-n{n}-r{r}' for s in ['paraboloid','saddle'] for m in ['disk','surface_area'] for rep in ([1,2,3] if a.mode=='main' else [1]) for n in ([480] if a.mode=='size' else [240]) for r in ([1,8,64] if a.mode=='pilot' else [64] if a.mode in ['extra','size'] else [1,2,4,8,16,32,64])]
jobs=[j for j in jobs if j != a.skip and (not a.case_filter or re.search(a.case_filter,j))]
def run(key):
 file=out/'inputs'/(key+'-distance.csv')
 if not file.exists():raise RuntimeError(f'Input missing: {file}; complete prepare-inputs first')
 env=dict(os.environ,OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',OMP_NUM_THREADS='1')
 start=time.time()
 with (out/'logs'/f'{a.mode}-{key}.log').open('w') as log:
  subprocess.run(['Rscript',str(here/'run.R'),key,'main' if a.mode=='size' else a.mode],stdout=log,stderr=subprocess.STDOUT,env=env,check=True)
 print('DONE',a.mode,key,round(time.time()-start,2),'seconds',flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:list(pool.map(run,jobs))
