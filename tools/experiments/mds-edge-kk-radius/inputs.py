#!/usr/bin/env python3
"""Checkpoint smooth inputs for the frozen radius protocol; run from repo root."""
import argparse, hashlib, importlib.util, json, sys, time
from pathlib import Path
import numpy as np
from scipy.spatial.distance import pdist, squareform
HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
OUT=REPO/'output/mds-edge-kk-radius/inputs'
OUT.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(HERE.parent/'saddle-mmds-radius'))
from smooth_geodesic import solve_pairs,pair_coordinates
spec=importlib.util.spec_from_file_location('parageo',HERE.parent/'paraboloid-mmds-radius/geodesic.py')
pg=importlib.util.module_from_spec(spec);spec.loader.exec_module(pg)
spec=importlib.util.spec_from_file_location('saddleexp',HERE.parent/'saddle-mmds-radius/experiment.py')
sg=importlib.util.module_from_spec(spec);spec.loader.exec_module(sg)

def main():
 p=argparse.ArgumentParser();p.add_argument('--surface',choices=['paraboloid','saddle'],required=True);p.add_argument('--rep',type=int,default=1);p.add_argument('--n',type=int,default=240);p.add_argument('--sampling',choices=['disk','surface_area'],required=True);a=p.parse_args()
 seed=20260904+a.rep;rng=np.random.default_rng(seed)
 q=rng.uniform(size=240);theta=rng.uniform(0,2*np.pi,240)
 if a.n==480:
  assert a.rep==1
  extra=np.random.default_rng(20261905);q=np.r_[q,extra.uniform(size=240)];theta=np.r_[theta,extra.uniform(0,2*np.pi,240)]
 assert len(q)==a.n
 old=np.load(REPO/f'output/{a.surface}-mmds-geodesic/embeddings.npz')
 solver=HERE.parent/('saddle-mmds-radius/smooth_geodesic.py' if a.surface=='saddle' else 'paraboloid-mmds-radius/geodesic.py')
 solver_hash=hashlib.sha256(solver.read_bytes()).hexdigest()
 previous=None;warm=None
 for r in [.1,.25,.5,1,2,4,8,16,32,64]:
  start=time.time();key=f'{a.surface}-{a.sampling}-rep{a.rep}-n{a.n}-r{r:g}'
  file=OUT/(key+'.npz')
  rho=r*np.sqrt(q) if a.sampling=='disk' else .5*np.sqrt(np.expm1((2/3)*np.log1p(q*np.expm1(1.5*np.log1p(4*r*r)))))
  u=(rho/r)[:,None]*np.c_[np.cos(theta),np.sin(theta)]
  truth=np.c_[r*u, rho*rho if a.surface=='paraboloid' else r*r*(u[:,0]**2-u[:,1]**2)]
  aa,bb=pair_coordinates(u);velocity=np.zeros_like(aa);diag={}
  if file.exists():
   cached=np.load(file);assert str(cached['solver_hash'])==solver_hash and np.array_equal(cached['u'],u)
   d=cached['d'];velocity=cached['velocity'];diag=json.loads(str(cached['diagnostics']))
  else:
   if a.rep==1 and a.n==240:
    d=old[f'{a.sampling}_r{r:g}_geodesic'];assert np.allclose(old[f'{a.sampling}_r{r:g}_truth'],truth,rtol=1e-14,atol=1e-14)
    if a.surface=='saddle':
     c=np.load(REPO/f'output/saddle-mmds-geodesic/distance_cache/{a.sampling}_r{r:g}_n240.npz');assert str(c['solver_hash'])==solver_hash
     velocity=c['velocity']
    diag={'origin':'audited original matrix'}
   elif a.surface=='paraboloid':d=pg.distance_matrix(rho,theta)
   else:
    velocity,distances,diag=solve_pairs(aa,bb,r,bb-aa if warm is None else warm,previous)
    d=squareform(distances)
   rms=np.sqrt(np.mean(squareform(d)**2));lower=pdist(truth)
   if a.surface=='saddle':upper=sg.straight_path_length(aa,bb,r)
   else:
    mer=pg.primitives(np.zeros_like(rho),rho)[1];ii,jj=np.triu_indices(a.n,1);gap=np.abs(theta[ii]-theta[jj]);gap=np.minimum(gap,2*np.pi-gap)
    upper=np.abs(mer[ii]-mer[jj])+np.minimum(rho[ii],rho[jj])*gap
   lower_violation=float(np.max(lower-squareform(d))/rms);upper_violation=float(np.max(squareform(d)-upper)/rms)
   # Exhaustive triangle check is cheap at these sample sizes.
   tri=max(float(np.max(d-d[:,j,None]-d[None,j,:])) for j in range(a.n))/rms
   assert lower_violation<1e-7 and upper_violation<1e-7 and tri<1e-7,(key,lower_violation,upper_violation,tri)
   diag.update(chord_lower_violation=lower_violation,path_upper_violation=upper_violation,triangle_violation=tri)
   np.savez_compressed(file,u=u,truth=truth,d=d,velocity=velocity,uniform=q,theta=theta,solver_hash=solver_hash,diagnostics=json.dumps(diag),seed=seed,elapsed=time.time()-start)
  previous=(aa,bb,r,velocity);warm=velocity
  if r>=1 and (a.n==240 or r==64):
   np.savetxt(OUT/(key+'-distance.csv'),d,delimiter=',',fmt='%.17g');np.savetxt(OUT/(key+'-truth.csv'),truth,delimiter=',',fmt='%.17g')
  print(key,round(time.time()-start,2),diag,flush=True)
if __name__=='__main__':main()
