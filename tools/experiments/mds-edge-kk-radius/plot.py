#!/usr/bin/env python3
"""Scientific figures from compact summaries; no fitting or private inputs."""
import argparse
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,LogNorm
from scipy.linalg import orthogonal_procrustes
p=argparse.ArgumentParser();p.add_argument('--partial',action='store_true');a=p.parse_args()
HERE=Path(__file__).resolve().parent;OUT=HERE.parents[2]/'output/mds-edge-kk-radius';SD=OUT/'summary' if a.partial else HERE/'summary'
FIG=HERE/'figures';FIG.mkdir(exist_ok=True)
s=pd.read_csv(SD/'scores.csv');main=s[s.n==240]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
methods=['classical','stress','stress_primary','stress_fixed_primary']
labels=['Classical MDS','Stress MDS','Stress MDS + edge-KK\nprofiled scale','Stress MDS + edge-KK\nfixed scale']
colors=['#3775ad','#d07d23','#237f6f','#8d5198']
def save(fig,name):
 fig.savefig(FIG/(name+'.png'),dpi=150);fig.savefig(FIG/(name+'.pdf'));plt.close(fig)
# Every radius/k cell: median of three independent, paired samples. Same scales by metric.
for surf in ['paraboloid','saddle']:
 for regime in ['geodesic','ambient']:
  for metric,title,norm in [('sigma2_sigma1','Second / first singular value',Normalize(0,1)),('sigma3_sigma2','Third / second singular value',Normalize(0,1)),('path_reference','Edge-calibrated path / smooth-geodesic relative error',LogNorm(1e-4,1)),('procrustes','Coordinate relative error after similarity alignment',LogNorm(1e-4,1))]:
   fig,axes=plt.subplots(2,4,figsize=(13,6.5),layout='constrained')
   for row,measure in enumerate(['disk','surface_area']):
    for col,(method,label) in enumerate(zip(methods,labels)):
     q=main[(main.surface==surf)&(main.regime==regime)&(main.sampling==measure)&(main.method==method)]
     v=q.pivot_table(index='radius',columns='k',values=metric,aggfunc='median').reindex(index=[1,2,4,8,16,32,64],columns=[4,8,16,32,64,128,239])
     ax=axes[row,col];im=ax.imshow(v,origin='lower',aspect='auto',norm=norm,cmap='viridis')
     ax.set_xticks(range(7),[4,8,16,32,64,128,239],fontsize=9);ax.set_yticks(range(7),[1,2,4,8,16,32,64]);ax.set_xlabel('Neighbors k')
     ax.axvline(5.5,color='white',linewidth=1.5)
     if row==0:ax.set_title(label)
     if col==0:ax.set_ylabel(('Uniform base disk' if row==0 else 'Uniform surface area')+'\nRadius r')
   fig.colorbar(im,ax=axes,shrink=.8,label=title,extend='min' if isinstance(norm,LogNorm) else 'neither')
   fig.suptitle(f'{surf.capitalize()} · {regime} neighbors and edge targets\n{title}; median of three samples, n=240',fontsize=14)
   save(fig,f'{surf}-{regime}-{metric}')
# Explicit replication spread at r64; no sampling confidence interval is implied.
fig,axes=plt.subplots(2,4,figsize=(14,7),layout='constrained')
for col,(surf,regime) in enumerate([(s,g) for s in ['paraboloid','saddle'] for g in ['geodesic','ambient']]):
 for row,measure in enumerate(['disk','surface_area']):
  ax=axes[row,col]
  for method,label,color in zip(methods,labels,colors):
   q=main[(main.surface==surf)&(main.regime==regime)&(main.sampling==measure)&(main.radius==64)&(main.method==method)]
   v=q.groupby('k').procrustes.agg(['median','min','max'])
   ax.plot(v.index,v['median'],marker='o',color=color,label=label.replace('\n',' '));ax.fill_between(v.index,v['min'],v['max'],color=color,alpha=.12)
  ax.set(xscale='log',yscale='log',xlabel='Neighbors k',title=f'{surf.capitalize()} · {regime}',ylim=(1e-5,1.1));ax.grid(alpha=.2)
  if col==0:ax.set_ylabel(('Base disk' if row==0 else 'Surface area')+'\nCoordinate relative error')
fig.legend(*axes[0,0].get_legend_handles_labels(),loc='outside lower center' if matplotlib.__version__>='3.7' else 'lower center',ncol=2,fontsize=9)
fig.suptitle('Sampling variation at r=64 · median and observed range across three samples',fontsize=14)
fig.subplots_adjust(bottom=.17) if not fig.get_constrained_layout() else None
save(fig,'replication')
# Differences from initializer and scale controls, preserving all k/r cases.
fig,axes=plt.subplots(2,3,figsize=(12,8),layout='constrained')
for row,policy in enumerate(['primary','fixed_primary']):
 for col,(metric,label) in enumerate([('path_rel','Fixed-path relative error'),('procrustes','Coordinate relative error'),('sigma3_sigma2','Third / second singular value')]):
  v=main.pivot(index=['case','regime','k'],columns='method',values=metric);ax=axes[row,col]
  for surf,c in [('paraboloid',colors[0]),('saddle',colors[1])]:
   use=v.index.get_level_values('case').str.startswith(surf)
   ax.scatter(v.loc[use,'stress'],v.loc[use,'stress_'+policy],s=8,alpha=.35,label=surf,color=c,rasterized=True)
  lo,hi=(1e-6,1) if col<2 else (0,1)
  ax.plot([lo,hi],[lo,hi],color='#666',ls='--',lw=1)
  if col<2:ax.set(xscale='log',yscale='log')
  ax.set(xlim=(lo,hi),ylim=(lo,hi),xlabel='Stress MDS',ylabel='After edge-KK',title=f'{label}\n'+('Profiled scale' if row==0 else 'Fixed scale'));ax.grid(alpha=.15)
axes[0,0].legend();fig.suptitle('Refinement changes path fidelity and shape differently · all n=240 cases',fontsize=14)
save(fig,'paired-refinement')
# Uniform-only versus continuation, same total iteration allowance.
fig,axes=plt.subplots(1,3,figsize=(12,4.4),layout='constrained')
v=main.pivot(index=['case','regime','k'],columns='method')
for ax,metric,label in zip(axes,['path_rel','procrustes','sigma3_sigma2'],['Fixed-path error','Coordinate error','Third / second singular value']):
 for fixed,c in [('',colors[2]),('fixed_',colors[3])]:
  ax.scatter(v[metric]['stress_'+fixed+'primary'],v[metric]['stress_'+fixed+'uniform'],s=8,alpha=.3,color=c,label='Fixed scale' if fixed else 'Profiled scale',rasterized=True)
 lo,hi=(0,1) if metric=='sigma3_sigma2' else (1e-6,1)
 ax.plot([lo,hi],[lo,hi],ls='--',c='#777');ax.set(xlim=(lo,hi),ylim=(lo,hi),xlabel='Density continuation',ylabel='Uniform only',title=label)
 if lo:ax.set(xscale='log',yscale='log')
axes[0].legend();fig.suptitle('Stiffness schedule sensitivity · paired fits on every n=240 graph');save(fig,'stiffness-controls')
if (s.n==480).any():
 for metric,title in [('path_reference','Path / smooth-geodesic error'),('sigma3_sigma2','Third / second singular value')]:
  fig,axes=plt.subplots(2,4,figsize=(14,7),layout='constrained')
  for col,(surf,regime) in enumerate([(s,g) for s in ['paraboloid','saddle'] for g in ['geodesic','ambient']]):
   for row,measure in enumerate(['disk','surface_area']):
    ax=axes[row,col]
    for n,ls,mk in [(240,'-','o'),(480,'--','s')]:
     for method,c in [('stress_primary',colors[2]),('stress_fixed_primary',colors[3])]:
      q=s[(s.surface==surf)&(s.regime==regime)&(s.sampling==measure)&(s.radius==64)&(s.replicate==1)&(s.n==n)&(s.method==method)].sort_values('k')
      ax.plot(q.k,q[metric],ls=ls,marker=mk,c=c,label=f'n={n}, '+('fixed' if 'fixed' in method else 'profiled'))
    ax.set(xscale='log',xlabel='Neighbors k',title=f'{surf.capitalize()} · {regime}');ax.grid(alpha=.2)
    if metric!='sigma3_sigma2':ax.set_yscale('log')
    else:ax.set_ylim(0,1)
    if col==0:ax.set_ylabel(('Base disk' if row==0 else 'Surface area')+'\n'+title)
  axes[0,0].legend(fontsize=8);fig.suptitle('Nested sample-size check at r=64 · fixed k and fixed k/n are distinct',fontsize=14);save(fig,'size-'+metric)
print('Figures written to',FIG)
# Deterministically selected spatial views, with isotropic scale only.
coordfile=SD/'snapshot-coordinates.csv'
if coordfile.exists():
 c=pd.read_csv(coordfile)
 for surf in ['paraboloid','saddle']:
  for regime in ['geodesic','ambient']:
   for k in [32,239]:
    fig=plt.figure(figsize=(13,6.8));gs=fig.add_gridspec(3,4,left=.065,right=.995,top=.85,bottom=.09,wspace=.015,hspace=.03)
    views={};limit=0
    for r in [1,8,64]:
     subset=c[(c.surface==surf)&(c.regime==regime)&(c.k==k)&(c.radius==r)]
     truth=subset[subset.method=='original'].sort_values('vertex')[['x','y','z']].to_numpy();centered=truth-truth.mean(0)
     if surf=='saddle':
      height=truth[:,0]**2-truth[:,1]**2;arms=np.where(height>=0,np.where(truth[:,0]>=0,0,1),np.where(truth[:,1]>=0,2,3));cc=np.array(['#2874a6','#ca8c20','#577b45','#b05789'])[arms]
     else:cc=np.arctan2(truth[:,1],truth[:,0])
     for method in ['original','stress','stress_primary','stress_fixed_primary']:
      q=subset[subset.method==method].sort_values('vertex');z=q[['x','y','z']].to_numpy();z=z-z.mean(0)
      if method!='original':rot,_=orthogonal_procrustes(z,centered);z=z@rot/q.edge_scale.iloc[0]
      z=z/(r*r);views[r,method]=(z,cc,truth.mean(0));limit=max(limit,abs(z).max())
    limit=np.ceil(limit*4)/4
    for row,r in enumerate([1,8,64]):
     for col,method in enumerate(['original','stress','stress_primary','stress_fixed_primary']):
      ax=fig.add_subplot(gs[row,col],projection='3d');z,cc,mean=views[r,method]
      ax.scatter(*z.T,c=cc,s=7,depthshade=False,cmap='twilight' if surf=='paraboloid' else None)
      if method=='original':
       theta=np.linspace(0,2*np.pi,49);rho=np.linspace(0,r,15);rr,tt=np.meshgrid(rho,theta);xx=rr*np.cos(tt);yy=rr*np.sin(tt);zz=xx*xx+(yy*yy if surf=='paraboloid' else -yy*yy)
       ax.plot_wireframe((xx-mean[0])/(r*r),(yy-mean[1])/(r*r),(zz-mean[2])/(r*r),rstride=4,cstride=3,color='#999999',linewidth=.35,alpha=.3)
      ax.set(xlim=(-limit,limit),ylim=(-limit,limit),zlim=(-limit,limit));ax.set_box_aspect((1,1,1));ax.view_init(elev=20,azim=-52)
      ax.set_xticks([-1,0,1]);ax.set_yticks([-1,0,1]);ax.set_zticks([-1,0,1]);ax.tick_params(labelsize=9,pad=-1)
      if row==0:ax.set_title(['Original surface','Stress MDS','Edge-KK: profiled scale','Edge-KK: fixed scale'][col],fontsize=11)
    for label,y in zip(['r = 1','r = 8','r = 64'],[.73,.47,.21]):fig.text(.02,y,label,rotation=90,va='center',fontsize=12)
    fig.suptitle(f'{surf.capitalize()} · {regime} graph · k={k} · first base-disk sample\nEdge-calibrated physical coordinates / r²; equal spatial units and common limits',fontsize=14)
    fig.text(.5,.015,'Original wireframe covers the full parameter disk. Thin equal-unit views do not mean the surface domain was truncated.',ha='center',fontsize=9)
    save(fig,f'spatial-{surf}-{regime}-k{k}')
# Hold full smooth-distance classical initialization fixed while k changes.
for metric,label in [('path_reference','Path / smooth-geodesic error'),('procrustes','Coordinate relative error')]:
 fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
 for col,surf in enumerate(['paraboloid','saddle']):
  for row,measure in enumerate(['disk','surface_area']):
   ax=axes[row,col]
   for r,color in [(1,'#3775ad'),(8,'#d07d23'),(64,'#237f6f')]:
    for method,ls,mk in [('classical_fixed_primary','-','o'),('full_classical_fixed_primary','--','s')]:
     q=main[(main.surface==surf)&(main.regime=='geodesic')&(main.sampling==measure)&(main.radius==r)&(main.method==method)].groupby('k')[metric].median()
     ax.plot(q.index,q,ls=ls,marker=mk,c=color,label=f'r={r}, '+('full Δ start' if method.startswith('full') else 'graph start'))
   ax.set(xscale='log',yscale='log',xlabel='Neighbors k',ylabel=label,title=f'{surf.capitalize()} · '+('base disk' if row==0 else 'surface area'));ax.grid(alpha=.2)
 axes[0,0].legend(fontsize=8,ncol=2);fig.suptitle('Fixed full-geodesic initializer versus graph-dependent initializer\nFixed-scale density continuation; median of three samples, n=240',fontsize=13)
 save(fig,'initializer-control-'+metric)
# The graph itself can introduce error before any embedding is fitted.
fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
for col,surf in enumerate(['paraboloid','saddle']):
 for row,measure in enumerate(['disk','surface_area']):
  ax=axes[row,col]
  for r,color in [(1,'#3775ad'),(8,'#d07d23'),(64,'#237f6f')]:
   for regime,ls,mk in [('geodesic','-','o'),('ambient','--','s')]:
    q=main[(main.surface==surf)&(main.regime==regime)&(main.sampling==measure)&(main.radius==r)&(main.method=='original')].groupby('k').graph_reference.median()
    ax.plot(q.index,q,ls=ls,marker=mk,c=color,label=f'r={r}, {regime}')
  ax.set(xscale='log',yscale='log',ylim=(1e-6,2),xlabel='Neighbors k',ylabel='Graph / smooth-geodesic relative error',title=f'{surf.capitalize()} · '+('base disk' if row==0 else 'surface area'));ax.grid(alpha=.2)
axes[0,0].legend(fontsize=8,ncol=2);fig.suptitle('Graph construction changes the target geometry before embedding\nUnscaled physical distances; median of three samples, n=240',fontsize=13)
save(fig,'graph-reference')
# Optimizer variation around and away from MDS initializations.
if (SD/'optimizer-sensitivity.csv').exists():
 ex=pd.read_csv(SD/'optimizer-sensitivity.csv');fig,axes=plt.subplots(2,3,figsize=(12,8),layout='constrained')
 for row,policy in enumerate(['','fixed_']):
  base=main[main.method=='stress_'+policy+'primary']
  for control,color,marker in [('perturbed','#3775ad','o'),('extended','#237f6f','s'),('random','#d07d23','^'),('original','#8d5198','x')]:
   q=ex[ex.method=='stress_'+policy+control].merge(base,on=['case','regime','k'],suffixes=('_extra','_base'))
   for col,(metric,label) in enumerate([('path_rel','Fixed-path error'),('procrustes','Coordinate error'),('sigma3_sigma2','Third / second singular value')]):
    ax=axes[row,col];ax.scatter(q[metric+'_base'],q[metric+'_extra'],c=color,marker=marker,s=22,alpha=.65,label=control)
    lo,hi=(0,1) if col==2 else (1e-6,1);ax.set(xlim=(lo,hi),ylim=(lo,hi),xlabel='Primary MDS-initialized fit',ylabel='Additional fit',title=label+'\n'+('Fixed scale' if policy else 'Profiled scale'))
    if col<2:ax.set(xscale='log',yscale='log')
  for ax in axes[row]:lo=0 if ax==axes[row,2] else 1e-6;ax.plot([lo,1],[lo,1],c='#777',ls='--',lw=1)
 axes[0,0].legend(fontsize=8);fig.suptitle('Optimizer sensitivity at r=64 · sample 1, both measures and graph regimes',fontsize=13);save(fig,'optimizer-controls')
