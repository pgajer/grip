#!/usr/bin/env python3
"""Generate numerical prose, assemble atlas, and render the validated PDF."""
import argparse,json,os,shutil,subprocess
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import numpy as np,pandas as pd
p=argparse.ArgumentParser();p.add_argument('--partial',action='store_true');a=p.parse_args()
HERE=Path(__file__).resolve().parent;OUT=HERE.parents[2]/'output/mds-edge-kk-radius';SD=OUT/'summary' if a.partial else HERE/'summary';BD=HERE/'build';BD.mkdir(exist_ok=True)
s=pd.read_csv(SD/'scores.csv');main=s[s.n==240];g=pd.read_csv(SD/'graphs.csv');starts=pd.read_csv(SD/'starts.csv')
if not a.partial:assert json.loads((SD/'validation.json').read_text())['status']=='passed'
def tex(s):return s.replace('%',r'\%').replace('→',r'$\to$').replace('–','--')
def write(name,text):(BD/(name+'.tex')).write_text(tex(text)+'\n')
def median(q,col):return q[col].median()
def interval(v):return f'{v.min():.4g}--{v.max():.4g}'
v=main.pivot(index=['case','regime','k'],columns='method');graphs=len(v)
counts={}
for policy in ['primary','fixed_primary']:
 for init in ['classical','stress']:
  delta=v.path_rel[init+'_'+policy]-v.path_rel[init]
  counts[init+'_'+policy]=(int((delta < -1e-8).sum()),int((delta>1e-8).sum()))
a64=main[(main.radius==64)&(main.k==32)&(main.method=='stress_fixed_primary')]
shape=[]
for surf in ['paraboloid','saddle']:
 for regime in ['geodesic','ambient']:
  q=a64[(a64.surface==surf)&(a64.regime==regime)&(a64.sampling=='disk')]
  shape.append(f'{surf.capitalize()} & {regime} & {q.sigma2_sigma1.median():.3f} & {q.sigma3_sigma2.median():.3f} & {100*q.procrustes.median():.2f} & {100*q.path_reference.median():.2f}'+r'\\')
findings=(f'The graph regime materially changes the geometric answer. At $r=64$, $k=32$, '
 'the fixed-scale stress-MDS-to-edge-KK branch gives the following base-disk medians across independent samples. Coordinate error uses similarity alignment; path/reference error compares edge-calibrated retained paths with smooth geodesics. '
 +'\n'+r'\begin{center}\begin{tabular}{llrrrr}\toprule Surface & Graph regime & $s_2/s_1$ & $s_3/s_2$ & Coordinate (%) & Path/reference (%)\\\midrule'+'\n'+'\n'.join(shape)+'\n'+r'\bottomrule\end{tabular}\end{center}'+'\n\n'
 f'Across {graphs} primary $n=240$ graphs, profiled-scale edge-KK lowers fixed-path error '
 f'in {counts["stress_primary"][0]} stress-initialized comparisons and increases it in '
 f'{counts["stress_primary"][1]}; fixed-scale edge-KK lowers it in '
 f'{counts["stress_fixed_primary"][0]} and increases it in {counts["stress_fixed_primary"][1]}. '
 'Changes smaller than $10^{-8}$ in absolute relative-error units are counted as ties. '
 'Good path fidelity is therefore neither a universal refinement improvement nor evidence of surface recovery.\n\n'
 'The pilot also exposed a scale issue: the profiled, unnormalized edge objective can decrease by '
 'uniform contraction. Fixed-scale controls prevent this degeneracy and are essential when interpreting '
 'absolute loss, stopping values, and spatial extent.')
write('report-findings',findings)
geo=pd.read_csv(HERE/'summary/geodesic-validation.csv');inp=pd.read_csv(HERE/'summary/input-validation.csv')
validation=(f'The experiment contains {len(g)} primary graphs, {len(s)} scored primary candidates, and '
 f'{len(starts)} MDS starts. {int((g.bridges>0).sum())} graphs require MST augmentation; '
 f'the largest retained-route/strict-distance discrepancy is {g.route_strict_max.max():.3g} in normalized units. '
 f'Independent edge/path calculations agree with grip to {s.package_score_difference.max():.3g}; '
 f'explicit route sums agree to {s.independent_path_difference.max():.3g}.\n\n'
 f'All 88 numerical input matrices pass symmetry, all-pair chord lower bounds, explicit surface-path '
 f'upper bounds, and exhaustive sampled-vertex triangle checks at tolerance $10^{{-7}}$ times RMS distance. '
 f'The largest normalized triangle violation is {inp.triangle_violation.max():.3g}. '
 f'Independent high-radius validation covers {len(geo)} actual pairs in {geo.case.nunique()} cases, '
 f'including short/far and random pairs. The largest sampled relative discrepancy is '
 f'{geo.relative_independent.max():.3g}. Saddle checks use tighter shooting and independent SciPy collocation/integration; '
 'paraboloid checks use independent quadrature. This sampled maximum is not an all-pair accuracy bound.')
if (HERE/'summary/coordinate-validation.csv').exists():
 cv=pd.read_csv(HERE/'summary/coordinate-validation.csv').iloc[0]
 validation+=f' A separate reconstruction of graph distances and diagnostics from {int(cv.candidates)} saved coordinate arrays agrees to {cv.max_scaled_discrepancy:.3g} after scaling dimensional discrepancies.'
profiled=main[main.method.isin(['classical_primary','stress_primary'])]
fixed=s[s.method.str.contains('fixed')]
validation+=f' Primary profiled fits reach edge-calibration factors as small as {profiled.edge_scale.min():.3g}, whereas the fixed-scale controls range from {fixed.edge_scale.min():.5f} to {fixed.edge_scale.max():.5f}. The contraction occurs at intermediate radii too, so it is not an asymptotic property of the surface.'
write('report-validation',validation)
spread=starts.groupby(['case','regime','k']).normalized_stress.agg(['min','max']);spread=np.sqrt(spread['max'])-np.sqrt(spread['min'])
opt=(f'{int((starts.termination=="iteration_limit").sum())}/{len(starts)} MDS starts reach the iteration cap; '
 f'the largest between-start difference in target-normalized RMSE is {100*spread.max():.3f} percentage points. '
 'Iteration limits and small stress changes do not certify global optimality.')
if (SD/'optimizer-sensitivity.csv').exists():
 ex=pd.read_csv(SD/'optimizer-sensitivity.csv');opt+=f' The additional controls contain {len(ex)} scored candidates on {ex.groupby(["case","regime","k"]).ngroups} graphs.'
 for pol,label in [('', 'profiled'),('fixed_', 'fixed')]:
  base=main[main.method=='stress_'+pol+'primary']
  perturb=ex[ex.method=='stress_'+pol+'perturbed'].merge(base,on=['case','regime','k'],suffixes=('_extra','_base'))
  rand=ex[ex.method=='stress_'+pol+'random'].merge(base,on=['case','regime','k'],suffixes=('_extra','_base'))
  ext=ex[ex.method=='stress_'+pol+'extended'].merge(base,on=['case','regime','k'],suffixes=('_extra','_base'))
  opt+=(f' Under {label} scale, the largest $s_3/s_2$ change is '
    f'{abs(perturb.sigma3_sigma2_extra-perturb.sigma3_sigma2_base).max():.4f} for tiny perturbations, '
    f'{abs(ext.sigma3_sigma2_extra-ext.sigma3_sigma2_base).max():.4f} for extra steps, and '
    f'{abs(rand.sigma3_sigma2_extra-rand.sigma3_sigma2_base).max():.4f} for independent random starts.')
 opt+=' These distinguish local stability around MDS starts from dependence on the choice of a different starting configuration.'
write('report-optimizer',opt)
implications=('The results support presenting edge-KK as a graph-edge refinement with separately assessed '
 'path and geometric fidelity. They do not support a general claim that it reverses MDS flattening or '
 'recovers the generating manifold. The graph/reference discrepancy, neighborhood locality, and '
 'scale policy must be visible alongside the final graph-to-layout score. '
 'The complete ambient control is exact Euclidean realization, not evidence of intrinsic-geodesic recovery.')
write('report-implications',implications)

sensitivity=[]
for surf in ['paraboloid','saddle']:
 for regime in ['geodesic','ambient']:
  q=main[(main.surface==surf)&(main.regime==regime)&(main.radius==64)&(main.method=='stress_fixed_primary')]
  sensitivity.append(f'{surf.capitalize()}, {regime} graphs: $s_2/s_1$ ranges from {q.sigma2_sigma1.min():.3f} to {q.sigma2_sigma1.max():.3f}, and $s_3/s_2$ from {q.sigma3_sigma2.min():.3f} to {q.sigma3_sigma2.max():.3f}.')
sensitivity_text='At $r=64$, the fixed-scale stress-initialized branch gives the following ranges across all seven $k$ values, both sampling measures, and three samples. '+' '.join(sensitivity)+' These are observed finite-radius ranges, not limiting-dimension estimates.\n\n'
fixed_init=main[(main.regime=='geodesic')&(main.method=='full_classical_fixed_primary')].groupby('case').path_reference.agg(['min','max'])
sensitivity_text+=f'Holding the full-geodesic classical initializer fixed across $k$ still gives a nonzero path/reference-error range (above $10^{{-8}}$) in {int(((fixed_init["max"]-fixed_init["min"])>1e-8).sum())}/{len(fixed_init)} cases. Thus neighborhood effects are not explained solely by changing the MDS input. The initializer-control figures compare this branch with graph-dependent classical starts.\n\n'
if (SD/'size-comparison.csv').exists():
 size=pd.read_csv(SD/'size-comparison.csv');q=size[(size.method=='stress_fixed_primary')&(size.regime=='geodesic')]
 sensitivity_text+='Increasing sample size does not hold graph locality fixed. In the four nested geodesic cases below, keeping $k=8$ increases path/reference error, whereas doubling $k$ to preserve $k/n$ reduces it relative to the $n=240,k=8$ baseline. Each row is one nested sample pair under fixed-scale stress-MDS-to-edge-KK; values are percentages, not medians across three samples.\n'
 sensitivity_text+=r'\begin{center}\begin{tabular}{llrrr}\toprule Surface & Sampling & $n=240,k=8$ & $n=480,k=8$ & $n=480,k=16$\\\midrule'+'\n'
 for surf in ['paraboloid','saddle']:
  for measure in ['disk','surface_area']:
   one=q[(q.surface==surf)&(q.sampling==measure)&(q.comparison=='fixed k=8')].iloc[0];two=q[(q.surface==surf)&(q.sampling==measure)&(q.comparison=='fixed k/n from k=8')].iloc[0]
   sensitivity_text+=f'{surf.capitalize()} & '+('Base disk' if measure=='disk' else 'Surface area')+f' & {100*one.path_reference_n240:.2f} & {100*one.path_reference_n480:.2f} & {100*two.path_reference_n480:.2f}'+r'\\'+'\n'
 sensitivity_text+=r'\bottomrule\end{tabular}\end{center}'+'\nThe higher-$k$ and ambient comparisons are retained in the size figures and full table. One nested size increase does not establish a general sample-size law.\n\n'
if (SD/'initializer-comparison.csv').exists():
 init=pd.read_csv(SD/'initializer-comparison.csv')
 sensitivity_text+='Stress initialization is not uniformly superior after refinement. '
 for policy,label in [('primary','profiled'),('fixed_primary','fixed')]:
  q=init[init.policy==policy]
  sensitivity_text+=f'With {label} scale, stress initialization gives lower final path error in {q.stress_lower_path.sum()}/1176 comparisons, classical initialization in {q.classical_lower_path.sum()}, with {q.ties.sum()} ties at the declared tolerance. '
 sensitivity_text+='The stratified table retains reversals by surface, sampling measure, and graph regime; these paired counts are not independent statistical trials.'
write('report-sensitivity',sensitivity_text)

(BD/'report-build-info.tex').write_text(r'\renewcommand{\reportbuilddatetime}{'+datetime.now(ZoneInfo('America/New_York')).strftime('%Y-%m-%d %H:%M:%S %Z')+'}\n')
figures=[]
def add(name,title,caption):
 if (HERE/'figures'/(name+'.pdf')).exists():figures.append(r'\figpage{'+name+'}{'+title+'}{'+tex(caption)+'}')
for surf in ['paraboloid','saddle']:
 for regime in ['geodesic','ambient']:
  add(f'spatial-{surf}-{regime}-k32',f'{surf.capitalize()}: {regime} graph, $k=32$',
      'Sample 1, uniform base-disk sampling, n=240. Coordinates are aligned by an orthogonal transformation and translation, edge-calibrated, '
      'and divided uniformly by r². All axes use equal units and common limits. The original wireframe spans the full disk. '
      'A thin equal-unit view does not truncate the original parameter domain. Colors preserve vertex correspondence.')
add('graph-reference','Graph fidelity precedes embedding','The original-coordinate row supplies one graph/reference value per graph. Physical graph distances are compared with smooth geodesics without fitting a scale. Curves show medians over three n=240 samples at representative radii. Values below $10^{-6}$, including complete-geodesic numerical controls, fall outside the displayed range.')
add('paired-refinement','Path fidelity and geometry are different comparisons',
    'Each point is one sample/radius/graph setting. Points below the diagonal have smaller values after refinement. '
    'Logarithmic axes clip values below $10^{-6}$, including exact ambient controls; complete numerical values are in scores.csv. '
    'A lower singular-value ratio describes shape, and is not inherently a better fit.')
add('replication','Variation across independent samples',
    'At r=64, lines show the median and shaded bands the observed range over three independent samples. '
    'These ranges are not confidence intervals. Both sampling measures and graph regimes are shown. '
    'Values below $10^{-5}$ are outside the displayed logarithmic range, including exact complete-ambient controls.')
for metric,label in [('path_reference','Path agreement with smooth geodesics'),('sigma2_sigma1','Elongation'),('sigma3_sigma2','Planarity'),('procrustes','Coordinate agreement')]:
 for surf in ['paraboloid','saddle']:
  for regime in ['geodesic','ambient']:
   add(f'{surf}-{regime}-{metric}',f'{surf.capitalize()}, {regime} graph: {label.lower()}',
    'Every cell is the median over three independent n=240 samples at the indicated radius and neighborhood size. '
    'Rows compare base-disk and surface-area sampling; columns compare initializers and scale controls. '
    'The separated final column is the complete graph, k=239. Color scales are shared across methods and measures. '
    +('Values below $10^{-4}$ use the lowest color; the exact values are retained in the tables.' if metric in ['path_reference','procrustes'] else 'Both singular-value ratios must be read together to distinguish a line from a plane.'))
for metric in ['path_reference','procrustes']:
 add('initializer-control-'+metric,'Isolating the graph from the initializer','Solid curves initialize from graph distances; dashed curves reuse full-geodesic classical MDS across k at each radius and sample. Both use fixed-scale density continuation. Medians over three samples; both measures are shown. Panels use separate logarithmic y ranges. Changing k can still alter the refined geometry even when the starting configuration is fixed.')
add('optimizer-controls','Sensitivity to starts and additional steps','Each point compares a stress-MDS-initialized primary fit with a separate perturbation, random/original start, or extra-budget fit. Both scale policies are shown. Logarithmic panels omit values below $10^{-6}$; all values remain in the tables. These diagnose particular fits, not global optimality.')
add('stiffness-controls','Uniform stiffness versus density continuation','Paired runs use the same 1,000-step allowance. All radii, k values, sampling measures and replicates are included. Values below $10^{-6}$ lie outside logarithmic panels. The diagonal denotes equal values, not an optimality condition.')
for metric in ['path_reference','sigma3_sigma2']:
 add('size-'+metric,'Nested sample-size sensitivity','At r=64, the n=480 sample retains the original 240 points. Solid and dashed lines identify the two sample sizes; colors identify scale policies. Fixed k and fixed k/n comparisons are distinct. Path-error panels use separate logarithmic y ranges; shape-ratio panels share [0,1]. Only one nested sample is available, so these are sensitivity checks, not large-n estimates.')
for surf in ['paraboloid','saddle']:
 for regime in ['geodesic','ambient']:
  add(f'spatial-{surf}-{regime}-k239',f'Complete-graph spatial control: {surf}, {regime}',
   'The complete ambient graph has an exact original-coordinate realization. The geodesic control tests global distance constraints. '
   'Coordinates use edge calibration followed by one uniform r² normalization, equal axis units, and common limits. First base-disk sample only.')
(BD/'report-figures.tex').write_text('\n'.join(figures)+'\n')
if a.partial:(BD/'report-findings.tex').write_text('\\textbf{Partial development render: not the final study.}\n'+(BD/'report-findings.tex').read_text())
shutil.copy(HERE/'report.tex',BD/'report.tex')
for cmd in [['pdflatex','-interaction=nonstopmode','-halt-on-error','report.tex'],['bibtex','report'],['pdflatex','-interaction=nonstopmode','-halt-on-error','report.tex'],['pdflatex','-interaction=nonstopmode','-halt-on-error','report.tex']]:
 subprocess.run(cmd,cwd=BD,check=True,stdout=subprocess.DEVNULL)
shutil.copy(BD/'report.pdf',HERE/('report-partial.pdf' if a.partial else 'report.pdf'))
if not a.partial:
 table_md='| Surface | Graph regime | s₂/s₁ | s₃/s₂ | Coordinate error (%) | Path/reference error (%) |\n|---|---|---:|---:|---:|---:|\n'+'\n'.join('| '+' | '.join(x.strip().replace('\\','') for x in row.split('&'))+' |' for row in shape)
 md_findings=findings[:findings.index(r'\begin{center}')].rstrip()+'\n\n'+table_md+findings[findings.index(r'\end{center}')+len(r'\end{center}'):]
 (SD/'RESULTS.md').write_text('# Radius study: results\n\n[Full PDF](../report.pdf) · [Protocol and reproduction](../README.md)\n\n'+md_findings+'\n\n![Graph fidelity](../figures/graph-reference.png)\n\n'+validation+'\n\n'+opt+'\n\n'+implications+'\n\nAll factor combinations are in [scores.csv](scores.csv); paired initializer and sample-size comparisons are in [initializer-comparison.csv](initializer-comparison.csv) and [size-comparison.csv](size-comparison.csv).\n')
print('Report built with',len(figures),'figures')
