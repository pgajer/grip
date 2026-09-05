#!/usr/bin/env python3
"""Static scientific figures from tracked summaries; never fits a layout.

Chart contract: paired comparisons of five clouds, four fitted candidates,
25 graphs; no pooling that hides k. Blue=classical, orange=stress; circles
and solid lines=MDS, squares and dashed lines=refinement. Neutral control.
PNG/PDF export for a LaTeX research report, inspected at final page size.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

root = Path(__file__).resolve().parent
sd = root / 'summary'
figdir = root / 'figures'
figdir.mkdir(exist_ok=True)
s = pd.read_csv(sd / 'scores.csv')
g = pd.read_csv(sd / 'graph-selection.csv')
u = pd.read_csv(sd / 'surface-scores.csv')
a = pd.read_csv(sd / 'starts.csv')
e = pd.read_csv(sd / 'additional-budget-scores.csv')
methods = ['Classical MDS', 'Stress MDS', 'Classical MDS + edge-KK', 'Stress MDS + edge-KK']
short = ['Classical', 'Stress', 'Classical\n+ edge-KK', 'Stress\n+ edge-KK']
colors = ['#2463A0', '#D16B24', '#2463A0', '#D16B24']
markers = ['o', 'o', 's', 's']
styles = ['-', '-', '--', '--']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,
    'axes.labelsize':10,'axes.spines.top':False,'axes.spines.right':False,
    'pdf.fonttype':42,'axes.edgecolor':'#777777','text.color':'#252525',
    'axes.labelcolor':'#252525','xtick.color':'#444444','ytick.color':'#444444'})
def save(fig,name):
    fig.savefig(figdir / (name+'.pdf'),bbox_inches='tight')
    fig.savefig(figdir / (name+'.png'),dpi=170,bbox_inches='tight')
    plt.close(fig)
def legend(fig,y=.96):
    handles = [Line2D([0],[0],color=c,marker=m,linestyle=l,label=n)
               for c,m,l,n in zip(colors,markers,styles,methods)]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,y),ncol=2,frameon=False)
def style(ax):
    ax.grid(axis='y',color='#E5E5E5',lw=.6);ax.set_axisbelow(True)

primary = s[s.selected]
surf = u[(u.alignment=='similarity') & u.method.isin(methods)]
merged = primary.merge(surf[['replicate','k','method','surface_rms']],on=['replicate','k','method'])
metrics = [('path_rel','Fixed-path relative RMSE (%)',100),('edge_rel','Edge relative RMSE (%)',100),
 ('stress1','Chord profiled Stress-1 (%)',100),('xz_path_error','Surface-reference path error (%)',100),
 ('procrustes','Coordinate relative RMSE (%)',100),('surface_rms','Surface RMS (coordinate units)',1)]
fig,axes=plt.subplots(2,3,figsize=(11,5.2))
for ax,(key,title,mult) in zip(axes.flat,metrics):
    for r in range(1,6):
        y=[float(merged[(merged.replicate==r)&(merged.method==m)][key].iloc[0])*mult for m in methods]
        ax.plot(range(4),y,color='#BBBBBB',lw=.7,zorder=1)
        for j in range(4):
            ax.scatter(j+(r-3)*.018,y[j],color=colors[j],marker=markers[j],s=25,zorder=3)
    ax.set_xticks(range(4),short,fontsize=8);ax.set_title(title);ax.set_ylim(bottom=0);style(ax)
fig.suptitle('MDS initialization at the five selected graphs',fontsize=15,y=1.02)
fig.text(.5,.94,'Five paired clouds; n = 1,000 each. Geometry uses similarity alignment.',ha='center',fontsize=10)
fig.tight_layout(rect=[0,0,1,.94],h_pad=2)
save(fig,'selected-comparison')

fig,axes=plt.subplots(3,5,figsize=(11.2,6.4),sharex='col',sharey='row')
for col,r in enumerate(range(1,6)):
    k0=int(g[(g.replicate==r)&g.selected].k.iloc[0])
    for row,(key,label) in enumerate([('path_rel','Fixed-path error\n(%, log scale)'),('stress1','Chord Stress-1 (%)'),('xz_path_error','X→Z path error\n(%, log scale)')]):
        ax=axes[row,col]
        for j,m in enumerate(methods):
            d=s[(s.replicate==r)&(s.method==m)].sort_values('k')
            ax.plot(d.k,100*d[key],color=colors[j],marker=markers[j],ls=styles[j],lw=1.2,ms=4)
        if key=='xz_path_error':
            d=s[(s.replicate==r)&(s.method=='Original saddle')].sort_values('k')
            ax.plot(d.k,100*d.xg_error,color='#555555',ls=':',lw=1.3,label='X→G control')
        ax.axvline(k0,color='#999999',lw=.7,ls=':');style(ax)
        if row==1:ax.set_ylim(bottom=0)
        if row in (0,2):
            ax.set_yscale('log');ax.set_ylim((.07,5) if row==0 else (.15,6))
            ax.set_yticks([.1,.3,1,3],['0.1','0.3','1','3'])
            ax.minorticks_off()
        ax.set_xticks([32,k0,80]);ax.tick_params(axis='x',labelsize=9)
        if col==0:ax.set_ylabel(label)
        if row==0:ax.set_title(f'Cloud {r}; selected k = {k0}')
        if row==2:ax.set_xlabel('Neighbors k')
axes[1,0].set_ylim(0,100*s[s.method.isin(methods)].stress1.max()*1.08)
fig.suptitle('Distance fidelity across the frozen neighborhood choices',y=1.01,fontsize=15)
legend(fig,.98)
fig.text(.5,.005,'Five evaluated k values per cloud; segments are visual guides. Vertical dotted line: selected k. Gray dotted curve: X→G error.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.025,1,.9])
save(fig,'k-distance-sensitivity')

fig,axes=plt.subplots(2,5,figsize=(11.2,5.7),sharex='col',sharey='row')
for col,r in enumerate(range(1,6)):
    k0=int(g[(g.replicate==r)&g.selected].k.iloc[0])
    for row,key in enumerate(['procrustes','surface_rms']):
        ax=axes[row,col]
        for j,m in enumerate(methods):
            d=s[(s.replicate==r)&(s.method==m)].sort_values('k') if row==0 else surf[(surf.replicate==r)&(surf.method==m)].sort_values('k')
            ax.plot(d.k,d[key]*(100 if row==0 else 1),color=colors[j],marker=markers[j],ls=styles[j],lw=1.2,ms=4)
        ax.axvline(k0,color='#999999',lw=.7,ls=':');ax.set_ylim(bottom=0);style(ax)
        ax.set_xticks([32,k0,80])
        if col==0:ax.set_ylabel('Coordinate relative\nRMSE (%)' if row==0 else 'Surface RMS\n(coordinate units)')
        if row==0:ax.set_title(f'Cloud {r}')
        else:ax.set_xlabel('Neighbors k')
axes[0,0].set_ylim(0,100*s[s.method.isin(methods)].procrustes.max()*1.08)
axes[1,0].set_ylim(0,surf.surface_rms.max()*1.08)
fig.suptitle('Geometric agreement across neighborhood choices',y=1.03,fontsize=15)
legend(fig,1.0)
fig.text(.5,0,'Similarity alignment uses known vertex correspondence. Surface RMS compares meshes over the same parameter footprint.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.03,1,.9]);save(fig,'k-geometry-sensitivity')

fig,axes=plt.subplots(1,2,figsize=(10.5,3.8))
ax=axes[0]
for r in range(1,6):
    k0=int(g[(g.replicate==r)&g.selected].k.iloc[0]);d=a[(a.replicate==r)&(a.k==k0)]
    best=d.raw_stress.min()
    for _,row in d.iterrows():
        ax.scatter(r+(row.start-3.5)*.045,100*(row.raw_stress/best-1),
            marker='o' if row.initialization=='classical' else 'x',
            color='#2463A0' if row.initialization=='classical' else '#D16B24',s=32)
ax.set_xticks(range(1,6));ax.set_xlabel('Cloud');ax.set_ylabel('Raw stress above best start (%)');ax.set_ylim(bottom=-.3);style(ax)
ax.set_title('Six starts at each selected graph')
ax.legend(handles=[Line2D([],[],ls='',marker='o',color='#2463A0',label='Classical start'),Line2D([],[],ls='',marker='x',color='#D16B24',label='Random start')],frameon=False,fontsize=8)
ax=axes[1]
for j,(name,label) in enumerate([('classical_kk','Classical MDS + edge-KK'),('stress_kk','Stress MDS + edge-KK')]):
    p=primary[primary.method==label].sort_values('replicate');d=e[e.method==name].sort_values('replicate')
    for i,r in enumerate(range(1,6)):
        ax.plot([r-.12,r+.12],100*np.array([p.path_rel.iloc[i],d.path_rel.iloc[i]]),color=colors[j],lw=1.1)
        ax.scatter(r-.12,100*p.path_rel.iloc[i],marker='s',facecolors='none',edgecolors=colors[j],s=30)
        ax.scatter(r+.12,100*d.path_rel.iloc[i],marker='s',color=colors[j],s=30)
ax.set_xticks(range(1,6));ax.set_xlabel('Cloud');ax.set_ylabel('Fixed-path relative RMSE (%)');ax.set_ylim(bottom=0);style(ax)
ax.set_title('Additional 1,000 uniform edge-KK steps')
ax.text(.98,.04,'Open: primary; filled: extended\nBlue: classical; orange: stress',transform=ax.transAxes,ha='right',va='bottom',fontsize=8)
fig.tight_layout(w_pad=3);save(fig,'optimizer-sensitivity')

# First cloud fixed by protocol, no selection for a favorable recovered shape.
c=pd.read_csv(sd/'coordinates.csv.gz');r=1;k0=int(g[(g.replicate==r)&g.selected].k.iloc[0])
c=c[(c.replicate==r)&(c.k==k0)]
order=['Original saddle']+methods
X=c[c.method==order[0]].sort_values('vertex')[['x','y','z']].to_numpy()
def align(z,x):
    zz=z-z.mean(0);xx=x-x.mean(0);U,D,V=np.linalg.svd(zz.T@xx)
    return (D.sum()/(zz*zz).sum())*zz@U@V+x.mean(0)
fig=plt.figure(figsize=(13,3.4))
for j,m in enumerate(order):
    Z=c[c.method==m].sort_values('vertex')[['x','y','z']].to_numpy();Z=align(Z,X)
    ax=fig.add_subplot(1,5,j+1,projection='3d')
    ax.scatter(*Z.T,c=X[:,2],cmap='cividis',s=2,alpha=.85,rasterized=True)
    ax.set(xlim=(-1.35,1.35),ylim=(-1.35,1.35),zlim=(-1.35,1.35))
    ax.set_box_aspect((1,1,1));ax.view_init(elev=25,azim=-60)
    ax.set_xticks([-1,0,1]);ax.set_yticks([-1,0,1]);ax.set_zticks([-1,0,1]);ax.tick_params(labelsize=7,pad=0)
    ax.set_title(m.replace(' + ','\n+ '),fontsize=10,pad=5)
    ax.set_xlabel('x',labelpad=-6);ax.set_ylabel('y',labelpad=-6);ax.set_zlabel('z',labelpad=-6)
fig.suptitle('Corresponding points after similarity alignment: cloud 1, k = 71',fontsize=14,y=1.05)
fig.text(.5,0,'Identical limits and equal Euclidean units on all axes; color identifies original height. Similarity alignment cannot correct folds or correspondence error.',ha='center',fontsize=9)
fig.subplots_adjust(left=0,right=.97,wspace=.06,bottom=.08,top=.89);save(fig,'selected-shapes')
