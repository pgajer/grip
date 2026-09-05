#!/usr/bin/env python3
"""Validate complete coverage and numerical invariants in compact result tables."""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
p=argparse.ArgumentParser();p.add_argument('--pilot',action='store_true');a=p.parse_args()
HERE=Path(__file__).resolve().parent;OUT=HERE.parents[2]/'output/mds-edge-kk-radius';sd=OUT/'summary' if a.pilot else HERE/'summary'
s=pd.read_csv(sd/'scores.csv');g=pd.read_csv(sd/'graphs.csv');st=pd.read_csv(sd/'starts.csv')
assert len(g)==(72 if a.pilot else 1216),len(g)
assert not s.duplicated(['case','regime','k','method']).any()
assert not g.duplicated(['case','regime','k']).any()
assert np.isfinite(s.select_dtypes('number')).all().all()
assert s.package_score_difference.max()<1e-10 and s.independent_path_difference.max()<1e-9
assert len(st)==3*len(g)
assert np.isfinite(st.raw_stress).all()
assert (st.groupby(['case','regime','k']).selected.sum()==1).all()
for (case,regime,k),group in s.groupby(['case','regime','k']):
 initial=['classical','stress']+(['full_classical'] if regime=='geodesic' else [])
 expected_methods=set(['original']+initial+[nm+'_'+policy for nm in initial for policy in ['primary','uniform','fixed_primary','fixed_uniform']])
 assert set(group.method)==expected_methods,(case,regime,k)
if not a.pilot:
 expected={(f'{surf}-{measure}-rep{rep}-n240-r{r}',regime,k) for surf in ['paraboloid','saddle'] for measure in ['disk','surface_area'] for rep in [1,2,3] for r in [1,2,4,8,16,32,64] for regime in ['ambient','geodesic'] for k in [4,8,16,32,64,128,239]}
 expected|={(f'{surf}-{measure}-rep1-n480-r64',regime,k) for surf in ['paraboloid','saddle'] for measure in ['disk','surface_area'] for regime in ['ambient','geodesic'] for k in [8,16,64,128,479]}
 assert set(g[['case','regime','k']].itertuples(index=False,name=None))==expected
 assert len(s)==16416
v=s.pivot(index=['case','regime','k'],columns='method')
assert ((v.raw_stress.stress-v.raw_stress.classical)<1e-7*np.maximum(v.raw_stress.classical,1)+1e-7).all()
original=s[(s.method=='original')&(s.regime=='ambient')]
assert original.edge_rel.max()<1e-12 and original.path_rel.max()<1e-7
complete=s[(s.k==s.n-1)&(s.regime=='ambient')&(s.method.isin(['classical','stress']))]
assert complete.raw_target_rmse.max()<1e-6
fixed=s[s.method.str.contains('fixed')]
assert fixed.edge_scale.between(.05,20).all(),fixed.edge_scale.describe()
# Primary conclusions use all cases; neither better stress nor more dimensions is assumed.
summary=dict(status='passed',graphs=len(g),candidates=len(s),mds_starts=len(st),
 repaired_graphs=int((g.bridges>0).sum()),max_route_strict_discrepancy=float(g.route_strict_max.max()),
 max_package_score_difference=float(s.package_score_difference.max()),
 max_independent_path_difference=float(s.independent_path_difference.max()),
 min_profiled_edge_scale=float(s[s.method.isin(['classical_primary','stress_primary'])].edge_scale.min()),
 fixed_scale_range=[float(fixed.edge_scale.min()),float(fixed.edge_scale.max())],
 mds_iteration_limit=int((st.termination=='iteration_limit').sum()))
if not a.pilot:
 obj=pd.read_csv(sd/'objective-validation.csv');assert len(obj)==176
 assert obj.absolute_difference.max()<1e-8 and obj.native_trace_difference.max()<1e-8
 extra=pd.read_csv(sd/'optimizer-sensitivity.csv');assert len(extra.groupby(['case','regime','k']))==32
 assert len(extra)==32*16
 geo=pd.read_csv(sd/'geodesic-validation.csv');assert len(geo.groupby('case'))==16
 assert geo.relative_independent.max()<1e-5
 summary.update(geodesic_validation_pairs=len(geo),max_sampled_geodesic_discrepancy=float(geo.relative_independent.max()))
(sd/('pilot-validation.json' if a.pilot else 'validation.json')).write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
