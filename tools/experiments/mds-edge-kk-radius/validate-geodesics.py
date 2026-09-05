#!/usr/bin/env python3
"""Independent checks on actual high-radius input pairs, not all-pair bounds."""
import csv,importlib.util,json,sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp,solve_bvp,quad
HERE=Path(__file__).resolve().parent;OUT=HERE.parents[2]/'output/mds-edge-kk-radius'
sys.path.insert(0,str(HERE.parent/'saddle-mmds-radius'))
from smooth_geodesic import shoot
spec=importlib.util.spec_from_file_location('pg',HERE.parent/'paraboloid-mmds-radius/geodesic.py');pg=importlib.util.module_from_spec(spec);spec.loader.exec_module(pg)
rows=[]
files=sorted((OUT/'inputs').glob('*-r64.npz'));assert len(files)==16
for file in files:
 c=np.load(file);d=c['d'];n=len(d);i,j=np.triu_indices(n,1);ds=d[i,j]
 rng=np.random.default_rng(6731)
 # Include shortest/farthest pairs plus evenly spaced length ranks and random pairs.
 ids=np.unique(np.r_[np.argsort(ds)[np.linspace(0,len(ds)-1,8).astype(int)],rng.choice(len(ds),8,replace=False)])
 for pair in ids:
  a,b=c['u'][i[pair]],c['u'][j[pair]];target=ds[pair];r=64
  if file.name.startswith('saddle'):
   v,refined,ok,it,steps=shoot(a,b,r,c['velocity'][pair],tolerance=2e-11,rtol=2e-12)
   if not ok:v,refined,ok,it,steps=shoot(a,b,r,v,tolerance=2e-10,rtol=5e-13)
   assert ok,(file.stem,int(pair),'refined shooting did not converge')
   def ode(t,y):
    xx,yy,vv,ww=y;k=4*r*r;acc=k*(vv*vv-ww*ww)/(1+k*(xx*xx+yy*yy));return np.array([vv,ww,-acc*xx,acc*yy])
   sol=solve_ivp(ode,[0,1],np.r_[a,v],method='DOP853',rtol=2e-12,atol=2e-13,dense_output=True);assert sol.success
   endpoint=np.linalg.norm(sol.y[:2,-1]-b)
   grid=np.unique(np.r_[sol.t,np.linspace(0,1,60)])
   bvp=solve_bvp(ode,lambda ya,yb:np.r_[ya[:2]-a,yb[:2]-b],grid,sol.sol(grid),tol=1e-6,max_nodes=30000);assert bvp.success,(file,pair,bvp.message)
   vv=bvp.y[2:,0];independent=r*np.sqrt(np.dot(vv,vv)+4*r*r*(a[0]*vv[0]-a[1]*vv[1])**2)
   # Quadrature along the independent collocation solution checks its speed integral.
   length=quad(lambda t: r*np.sqrt(np.sum(bvp.sol(t)[2:]**2)+4*r*r*(bvp.sol(t)[0]*bvp.sol(t)[2]-bvp.sol(t)[1]*bvp.sol(t)[3])**2),0,1,epsabs=1e-6,epsrel=1e-8,limit=500)[0]
  else:
   ra,rb=r*np.linalg.norm(a),r*np.linalg.norm(b);angle=c['theta'][i[pair]]-c['theta'][j[pair]]
   refined,info=pg.pair_distances(ra,rb,angle,iterations=74,details=True)
   cc=float(info['c']);ta=float(info['ta']);tb=float(info['tb']);intervals=[(0,ta),(0,tb)] if bool(info['turning']) else [(ta,tb)]
   independent=sum(quad(lambda t:np.sqrt(1+4*cc*cc+4*t*t),l,h,epsabs=1e-10,epsrel=1e-12)[0] for l,h in intervals)
   length=independent;endpoint=np.nan
  rec=dict(case=file.stem,i=int(i[pair]+1),j=int(j[pair]+1),target=target,refined=float(refined),independent=independent,
   relative_refinement=abs(target-refined)/target,relative_independent=abs(target-independent)/target,
   relative_integrated_length=abs(length-independent)/independent,endpoint_over_radius=endpoint)
  assert rec['relative_independent']<1e-5 and rec['relative_refinement']<1e-5 and rec['relative_integrated_length']<1e-5,rec
  rows.append(rec)
 print('VALIDATED',file.stem,'max relative',max(z['relative_independent'] for z in rows if z['case']==file.stem),flush=True)
 with (OUT/'geodesic-validation.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print('DONE',len(rows),'sampled pairs; these are not all-pair error bounds',flush=True)
