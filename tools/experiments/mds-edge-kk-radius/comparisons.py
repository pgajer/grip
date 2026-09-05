#!/usr/bin/env python3
"""Paired summaries; counts describe comparisons, never independent trials."""
from pathlib import Path
import pandas as pd
p=Path(__file__).resolve().parent/'summary';s=pd.read_csv(p/'scores.csv');rows=[]
main=s[s.n==240]
for (surface,sampling,regime),group in main.groupby(['surface','sampling','regime']):
 v=group.pivot(index=['replicate','radius','k'],columns='method')
 for policy in ['primary','fixed_primary']:
  path=v.path_rel['stress_'+policy]-v.path_rel['classical_'+policy]
  coord=v.procrustes['stress_'+policy]-v.procrustes['classical_'+policy]
  rows.append(dict(surface=surface,sampling=sampling,regime=regime,policy=policy,comparisons=len(path),stress_lower_path=int((path < -1e-8).sum()),classical_lower_path=int((path>1e-8).sum()),ties=int((abs(path)<=1e-8).sum()),median_path_difference=path.median(),median_coordinate_difference=coord.median()))
pd.DataFrame(rows).to_csv(p/'initializer-comparison.csv',index=False)
rows=[]
for surface in ['paraboloid','saddle']:
 for sampling in ['disk','surface_area']:
  for regime in ['geodesic','ambient']:
   case=s[(s.surface==surface)&(s.sampling==sampling)&(s.regime==regime)&(s.replicate==1)&(s.radius==64)]
   for k240,k480,comparison in [(8,8,'fixed k=8'),(8,16,'fixed k/n from k=8'),(64,64,'fixed k=64'),(64,128,'fixed k/n from k=64'),(239,479,'complete graph')]:
    for method in ['classical','stress','stress_primary','stress_fixed_primary']:
     a=case[(case.n==240)&(case.k==k240)&(case.method==method)].iloc[0];b=case[(case.n==480)&(case.k==k480)&(case.method==method)].iloc[0]
     row=dict(surface=surface,sampling=sampling,regime=regime,method=method,comparison=comparison,k240=k240,k480=k480)
     for metric in ['path_rel','graph_reference','path_reference','procrustes','sigma2_sigma1','sigma3_sigma2']:
      row[metric+'_n240']=a[metric];row[metric+'_n480']=b[metric];row[metric+'_difference']=b[metric]-a[metric]
     rows.append(row)
pd.DataFrame(rows).to_csv(p/'size-comparison.csv',index=False)
# Scale contraction and start variation are stratified rather than pooled away.
main.groupby(['surface','sampling','regime','radius','method']).edge_scale.agg(['min','median','max']).reset_index().to_csv(p/'scale-by-radius.csv',index=False)
print('Initializer, size, and scale comparisons written')
