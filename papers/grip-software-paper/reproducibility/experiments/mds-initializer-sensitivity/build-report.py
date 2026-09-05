#!/usr/bin/env python3
"""Build the empirical readout and PDF from tracked compact summaries."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import json, os, shutil, subprocess
import numpy as np
import pandas as pd
p=Path(__file__).resolve().parent
sd=p/'summary';bd=p/'build';bd.mkdir(exist_ok=True)
s=pd.read_csv(sd/'scores.csv');a=pd.read_csv(sd/'starts.csv')
u=pd.read_csv(sd/'surface-scores.csv');extra=pd.read_csv(sd/'additional-budget-scores.csv')
sel=pd.read_csv(sd/'graph-selection.csv');gs=pd.read_csv(sd/'graph-statistics.csv')
status=pd.read_csv(sd/'edge-optimizer-status.csv');timing=pd.read_csv(sd/'timings.csv')
validation=json.loads((sd/'validation.json').read_text());assert validation['status']=='passed'
methods=['Classical MDS','Stress MDS','Classical MDS + edge-KK','Stress MDS + edge-KK']
short=['Classical MDS','Stress MDS','Classical + edge-KK','Stress + edge-KK']
prim=s[s.selected].copy()
piv=s.pivot(index=['replicate','k'],columns='method')
path_improve={m:int((piv['path_rel'][m+' + edge-KK']<piv['path_rel'][m]).sum()) for m in methods[:2]}
chord_worse={m:int((piv['stress1'][m+' + edge-KK']>piv['stress1'][m]).sum()) for m in methods[:2]}
win=int((piv['path_rel'][methods[3]]<piv['path_rel'][methods[2]]).sum())
near=s[s.role.isin(['nearby','selected'])].pivot(index=['replicate','k'],columns='method')
nearwin=int((near['path_rel'][methods[3]]<near['path_rel'][methods[2]]).sum())
rows=[]
for m,label in zip(methods,short):
 d=prim[prim.method==m]
 surf=u[(u.method==m)&(u.alignment=='similarity')].merge(sel[sel.selected][['replicate','k']],on=['replicate','k'])
 rows.append([label]+[100*d[v].median() for v in ['path_rel','edge_rel','stress1','xz_path_error','procrustes']]+[surf.surface_rms.median()])
table=pd.DataFrame(rows,columns=['Method','Path (%)','Edge (%)','Chord Stress-1 (%)','X→Z path (%)','Coordinate (%)','Surface RMS'])
def mdtable(df,digits=4):
 lines=['| '+' | '.join(df.columns)+' |','| '+' | '.join(['---']*len(df.columns))+' |']
 for row in df.itertuples(index=False,name=None):
  lines.append('| '+' | '.join(f'{x:.{digits}f}' if isinstance(x,(float,np.floating)) else str(x) for x in row)+' |')
 return '\n'.join(lines)
def rng(values,mult=1,digits=3):return f'{mult*values.min():.{digits}f}–{mult*values.max():.{digits}f}'
def tx(text):
 return text.replace('→',r'$\to$').replace('%',r'\%').replace('–','--').replace('&',r'\&')
# Paired selected-graph budget changes.
sel_starts=a.merge(sel[sel.selected][['replicate','k']],on=['replicate','k'])
spread=sel_starts.groupby('replicate').raw_stress.agg(['min','max'])
spread['excess']=100*(spread['max']/spread['min']-1)
best=prim[prim.method=='Stress MDS'].set_index('replicate')
em=extra[extra.method=='stress_mds'].set_index('replicate')
mds_drop=100*(1-em.raw_stress/best.raw_stress)
base_kk=prim[prim.method==methods[3]]
coord_range=rng(base_kk.procrustes,100,2)
mins=sel.groupby('replicate').xg_error.min()
sparse=sel[sel.role=='sparse'];dense=sel[sel.role=='dense']
edgecontrol=s[s.method=='Original saddle']
worst_near=near['procrustes'][methods[3]].max()*100
findings=(f'Across all 25 evaluated graphs, edge-KK reduces fixed-path error in '
 f'{path_improve[methods[0]]}/25 classical-initialized and {path_improve[methods[1]]}/25 '
 f'stress-initialized comparisons. It increases chord profiled Stress-1 in '
 f'{chord_worse[methods[0]]}/25 and {chord_worse[methods[1]]}/25 comparisons, respectively. '
 f'Stress initialization gives the lower final path error on {win}/25 graphs '
 f'({nearwin}/15 at the selected and nearby choices). This is an empirical preference, '
 f'not uniform superiority over graph choices or global minimizers.')
selected_text=(f'At the selected graphs, median path error is {table.iloc[0,1]:.3f}% for classical MDS '
 f'and {table.iloc[1,1]:.3f}% for stress MDS; after edge-KK it is '
 f'{table.iloc[2,1]:.3f}% and {table.iloc[3,1]:.3f}%, respectively. '
 f'The corresponding median chord scores increase from {table.iloc[0,3]:.3f}% '
 f'to {table.iloc[2,3]:.3f}% and from {table.iloc[1,3]:.3f}% to {table.iloc[3,3]:.3f}%.')
geometry_text=(f'Geometric recovery is less consistent. The stress-initialized edge-KK coordinate '
 f'error at selected graphs ranges from {coord_range}% after similarity alignment. '
 f'Good fixed-path agreement therefore does not establish correspondence or surface recovery. '
 f'Across the selected and nearby choices, its largest coordinate error is {worst_near:.2f}%.')
optimizer_text=(f'At selected graphs, the worst of six achieved raw stresses exceeds the best by '
 f'{rng(spread.excess,1,2)}%. Across the full sweep, '
 f'{int((a.termination=="iteration_limit").sum())}/150 primary starts reach the iteration limit; '
 f'{int(((a.termination=="iteration_limit") & a.selected).sum())}/25 selected starts do so. '
 f'The other primary starts stop at the backend stress-change tolerance. The selected-graph '
 f'2,000-step tighter-tolerance continuation lowers raw stress by {rng(mds_drop,1,3)}%. '
 f'All {len(status)} recorded edge-KK stages '
 f'{"reach their iteration caps" if status.hit_budget.all() else "have stopping and gradient records in the accompanying table"}. '
 f'The achieved finite-budget fits must not be called global optima.')
(bd/'report-build-info.tex').write_text(r'\renewcommand{\reportbuilddatetime}{'+datetime.now(ZoneInfo('America/New_York')).strftime('%Y-%m-%d %H:%M:%S %Z')+'}\n')
(bd/'report-values.tex').write_text('% Numerical text is generated from the tracked summaries.\n')
(bd/'report-findings.tex').write_text(tx(findings)+'\n\n'+tx(selected_text)+'\n\n'+tx(geometry_text)+'\n')
(bd/'report-optimizer.tex').write_text(tx(optimizer_text)+'\n\nElapsed component times were measured with concurrent cloud workers. They are not a controlled performance benchmark.\n')
latex=[r'\begin{center}\small',r'\begin{tabular}{lrrrrrr}\toprule',
 r'Method & Path (\%) & Edge (\%) & Chord (\%) & $X\to Z$ (\%) & Coord. (\%) & Surface RMS\\\midrule']
for row in table.itertuples(index=False,name=None):latex.append(row[0]+' & '+' & '.join(f'{x:.3f}' for x in row[1:])+r'\\')
latex.extend([r'\bottomrule\end{tabular}',r'\end{center}',r'\normalsize'])
(bd/'report-selected-table.tex').write_text('\n'.join(latex)+'\n')
# Keep all exact per-cloud selected values easy to inspect.
selected_detail=prim[prim.method!='Original saddle'][['replicate','k','method','path_rel','edge_rel','stress1','xz_path_error','procrustes']].copy()
for v in selected_detail.columns[3:]:selected_detail[v]*=100
selected_detail.to_csv(sd/'selected-scores-percent.csv',index=False)
# Characterize reversals without selecting them out of the experiment.
reverse=piv['path_rel'][methods[3]]>=piv['path_rel'][methods[2]]
reversal='; '.join(f'cloud {r}, k={k}' for r,k in reverse[reverse].index) or 'none'
nearby=sel[sel.role=='nearby'];selected_g=sel[sel.selected]
range_text=(f'The selected-graph X→G errors are {rng(selected_g.xg_error,100,3)}%. '
 f'They are {rng(sparse.xg_error,100,3)}% at k=32 and {rng(dense.xg_error,100,3)}% at k=80. '
 f'The original-coordinate control has essentially zero graph-to-layout edge/path error at every k, '
 f'while retaining the graph/reference error. This demonstrates why graph-to-layout fidelity '
 f'cannot substitute for X→G validation.')
text=f'''# MDS initialization and neighborhood sensitivity: results

{findings}

{selected_text}

{geometry_text}

[Full PDF report](../report.pdf) · [Protocol and reproduction](../README.md)

## Scope

Five saved area-uniform clouds, n=1,000, on z=0.8(x²−y²) over [-1,1]²; 25 saved
symmetric-union kNN graphs; four fitted 3D candidates and the original-coordinate
control. MDS input is the strict graph shortest-path matrix built from ambient
chord-weighted edges. Numerical smooth-surface geodesics are evaluation references,
not the MDS inputs in this phase. Increasing-radius/geodesic-input experiments
remain a separate phase.

## Selected graphs

All entries below are medians across five paired clouds. Edge and path scores
use separately fitted target scales; chord is profiled Stress-1. Coordinate and
surface scores use similarity alignment. Surface RMS is in coordinate units.

{mdtable(table,3)}

![Selected comparison](../figures/selected-comparison.png)

Exact per-cloud scores, in percent where appropriate: [selected-scores-percent.csv](selected-scores-percent.csv).
All returned-scale raw stresses, raw target-normalized RMSE, literal Stress-1,
rigid coordinate error, and singular-value ratios: [scores.csv](scores.csv).

## Neighborhood sensitivity

{range_text}

The final-path advantage of stress over classical initialization reverses at:
**{reversal}**. These are retained results. Each branch still has its own
initializer-to-refinement comparison in the full tables. The experiment does
not establish monotonic improvement with k or a universally optimal initializer.

All {len(gs)} evaluated graphs have connectivity and degree/length statistics
in [graph-statistics.csv](graph-statistics.csv). Components before repair range
from {int(gs.components_before.min())} to {int(gs.components_before.max())}; the
total number of added bridges over these graphs is {int(gs.bridges.sum())}.
The full calibration plateau and reference-sensitivity fields are retained in
[graph-selection.csv](graph-selection.csv). The five k values per cloud do not
cover every integer in the plausible region.

![Distance sensitivity](../figures/k-distance-sensitivity.png)

## Geometric agreement

![Geometric sensitivity](../figures/k-geometry-sensitivity.png)

Coordinate correspondence and closest-surface agreement answer different
questions from edge/path fidelity. Surface scoring uses the original parameter
triangulation, a twice-subdivided lifted reference, and 8,000 samples per direction.
Original-mesh surface RMS is {rng(u[(u.method=='Original saddle') & (u.alignment=='similarity')].surface_rms,1,4)};
this discretization control is not subtracted from candidate scores. Both rigid
and similarity alignment and surface Monte Carlo SE are recorded in
[surface-scores.csv](surface-scores.csv). Surface RMS is not a fold, topology,
or injectivity certificate.

![First-cloud illustration](../figures/selected-shapes.png)

The visual example uses cloud 1 by a fixed first-cloud rule. All panels have equal
units and the same limits after similarity alignment. Color identifies the
original height of corresponding vertices. No anisotropic display scaling is used.

## Optimization and additional budgets

{optimizer_text}

![Optimizer sensitivity](../figures/optimizer-sensitivity.png)

The additional MDS result is a continuation of the selected fit, not a new
multi-start search, and does not replace the primary initializer. Each extra
edge-KK fit extends its own primary fit. Small raw-stress changes can accompany
noticeable coordinate changes; an objective tolerance is not a recovery guarantee.
See [starts.csv](starts.csv), [additional-budget-scores.csv](additional-budget-scores.csv),
[additional-mds-starts.csv](additional-mds-starts.csv), and
[edge-optimizer-status.csv](edge-optimizer-status.csv).
Elapsed times in [timings.csv](timings.csv) are concurrent-worker component
measurements, not a controlled performance benchmark.

## Numerical validation and input limitation

The complete verifier passed: {validation['graphs']} graphs, {validation['primary_candidates']}
primary score rows, {validation['smacof_starts']} starts, {validation['surface_rows']}
surface rows, and {validation['additional_fits']} additional fits. Independent
recomputation, control realizations, raw/profile scale identities, least-stress
selection, descent within edge-KK stages, and input/result checksums were checked.
The maximum historical score discrepancy is {validation['max_baseline_score_difference']:.3g};
all original pilot inputs remain unchanged. See [validation.json](validation.json).

An initial partial run exposed near-tie asymmetry in grip's retained-path
preparation (2.67e-8 on one graph), rejected by the metric-MDS symmetry validator.
Every final fit was rerun with the saved symmetric strict graph distances under
protocol v2. Retained routes still define path diagnostics; strict distances
define the chord objective. The partial run is excluded. The package validator
was not weakened and no package algorithm was modified. A general package repair
for path-prepared near-tie inputs remains a separate integration issue; the
current experiment supplies its verified symmetric inputs explicitly.

Reference-geodesic numerical validation is inherited from the original pilot.
Its sampled discrepancies are not certified all-pairs bounds or statistical
confidence intervals. These five samples and finite optimization budgets do not
establish population convergence or uniqueness of an embedding.

## Recommendation for the R Journal manuscript

Include a concise initializer comparison and plausible-k sensitivity result in
the main paper. The supported message is that edge-KK improves retained-path
fidelity under both MDS initializers, with a chord-fidelity tradeoff; evaluating
the graph against a reference remains a separate requirement. Describe the
stress-initializer advantage as frequent in these samples, not uniform. Avoid
claiming that low path error recovers the generating saddle or its pointwise
geometry. Put the full k grid, starts, graph statistics, surface/coordinate
checks, and budget sensitivity in supplementary reproducibility material.

The study uses the existing ratio-SMACOF implementation described by
[Mair, Groenen, and de Leeuw (2022)](https://doi.org/10.18637/jss.v102.i10);
its substantive citation evidence is in [citation_verification.html](../citation_verification.html).
The algebraic scale identity and numerical results here are independently checked.

Phase 2 does not edit or replace the main manuscript, old fits, or public package
API. Radius studies and manuscript integration remain later phases.
'''
(sd/'results.md').write_text(text)
env=os.environ.copy();env['TEXINPUTS']=str(bd)+os.pathsep+env.get('TEXINPUTS','')
for cmd in [['pdflatex','-interaction=nonstopmode','-halt-on-error','-output-directory=build','report.tex'],
 ['bibtex','build/report'],['pdflatex','-interaction=nonstopmode','-halt-on-error','-output-directory=build','report.tex'],
 ['pdflatex','-interaction=nonstopmode','-halt-on-error','-output-directory=build','report.tex']]:
 result=subprocess.run(cmd,cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
 with (bd/'build.log').open('a') as f:f.write(result.stdout)
 if result.returncode:raise RuntimeError(result.stdout[-5000:])
shutil.copyfile(bd/'report.pdf',p/'report.pdf')
subprocess.run(['python3',str(p/'check-citations.py')],check=True)
subprocess.run(['python3',str(p/'check-artifacts.py'),'--write'],check=True)
print(p/'report.pdf')
