#!/usr/bin/env python3
"""Freeze compact, portable fitting inputs and sample coverage diagnostics."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
HERE=Path(__file__).resolve().parent;OUT=HERE.parents[2]/'output/mds-edge-kk-radius';SD=HERE/'summary';(SD/'inputs').mkdir(parents=True,exist_ok=True)
files=[f for f in sorted((OUT/'inputs').glob('*.npz')) if float(f.stem.split('-r')[-1])>=1 and ('n240' in f.stem or f.stem.endswith('-r64'))]
assert len(files)==88,len(files)
rows=[]
for f in files:
 c=np.load(f);n=len(c['truth']);target=SD/'inputs'/f.name
 np.savez_compressed(target,d=c['d'],truth=c['truth'],u=c['u'],uniform=c['uniform'],theta=c['theta'])
 angle=c['theta'];rho=np.linalg.norm(c['u'],axis=1);hist=np.histogram(angle,bins=np.linspace(0,2*np.pi,13))[0]
 diag=json.loads(str(c['diagnostics']))
 rows.append(dict(case=f.stem,seed=int(c['seed']),n=n,radius=float(f.stem.split('-r')[-1]),
  radial_min=float(rho.min()),radial_median=float(np.median(rho)),radial_max=float(rho.max()),
  angular_12bin_min=int(hist.min()),angular_12bin_max=int(hist.max()),
  max_angular_gap=float(np.diff(np.r_[np.sort(angle),np.min(angle)+2*np.pi]).max()),
  solver_sha256=str(c['solver_hash']),input_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
  chord_lower_violation=diag['chord_lower_violation'],path_upper_violation=diag['path_upper_violation'],triangle_violation=diag['triangle_violation']))
pd.DataFrame(rows).to_csv(SD/'input-validation.csv',index=False)
print('Frozen',len(rows),'portable matrices and sample-coverage diagnostics')
