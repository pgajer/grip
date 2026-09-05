#!/usr/bin/env python3
"""MDS of continuous smooth-paraboloid geodesic distances; see GEODESIC.md."""
import argparse
from datetime import datetime
from zoneinfo import ZoneInfo
import hashlib
import json
from pathlib import Path
import platform

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import scipy
from scipy.linalg import orthogonal_procrustes
from scipy.spatial.distance import pdist, squareform

from experiment import center, classical, stress_fit, metrics, write_csv, BLUE, GOLD, INK, CMAP
from geodesic import distance_matrix, validate, primitives

ROOT=Path(__file__).resolve().parent
RADII=[.1,.25,.5,1.,2.,4.,8.,16.,32.,64.]


def pca_display(x,q):
    x=center(x)
    _,_,vt=np.linalg.svd(x,full_matrices=False)
    x=x@vt.T
    if np.corrcoef(x[:,0],q)[0,1]<0:
        x[:,0]*=-1
    return x


def quadratic_shape(x,truth):
    """Native-unit z=a(x²+y²)+b x+c y+d after rigid alignment to input."""
    x=center(x)
    rotation,_=orthogonal_procrustes(x,center(truth))
    x=x@rotation
    design=np.c_[np.sum(x[:,:2]**2,axis=1),x[:,:2],np.ones(len(x))]
    coeff=np.linalg.lstsq(design,x[:,2],rcond=None)[0]
    residual=x[:,2]-design@coeff
    return dict(quadratic_coefficient=float(coeff[0]),
                quadratic_r_squared=float(1-np.sum(residual**2)/np.sum(x[:,2]**2)))


def plots(out,rows,saved,n):
    # Scientific static chart contract: exact same-unit spatial panels and
    # ten-point radius sweeps; method styles and radial-position colors.
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,
                         'axes.spines.top':False,'axes.spines.right':False,
                         'text.color':INK,'axes.labelcolor':INK})
    selected=[.25,1.,4.,32.]
    fig=plt.figure(figsize=(13,11))
    gs=fig.add_gridspec(3,4,left=.09,right=.91,top=.85,bottom=.10,wspace=.12,hspace=.10)
    for j,r in enumerate(selected):
        item=saved[('disk',r)]; q=item['q']; rms=item['rms']
        for row,method in enumerate(['truth','classical3','stress3']):
            # Rigid alignment to the input only; no rescaling in the alignment.
            x=center(item[method])
            rot,_=orthogonal_procrustes(x,center(item['truth']))
            x=x@rot/rms
            ax=fig.add_subplot(gs[row,j],projection='3d')
            pts=ax.scatter(*x.T,c=q,cmap=CMAP,vmin=0,vmax=1,s=8,alpha=.9,depthshade=False)
            ax.set(xlim=(-1.5,1.5),ylim=(-1.5,1.5),zlim=(-1.5,1.5))
            if row==0: ax.set_title(f'r = {r:g}',fontsize=13)
            ax.set_box_aspect((1,1,1)); ax.view_init(elev=18,azim=-55)
            for axis in [ax.xaxis,ax.yaxis,ax.zaxis]: axis.set_ticks([-1,0,1])
            ax.tick_params(labelsize=7,pad=-1)
            ax.set_xlabel('x',labelpad=-3);ax.set_ylabel('y',labelpad=-3);ax.set_zlabel('z',labelpad=-3)
    for label,y in [('Original paraboloid',.735),('Classical MDS · 3D',.48),('Raw-stress MDS · 3D',.235)]:
        fig.text(.025,y,label,rotation=90,va='center',fontsize=12)
    fig.suptitle('MDS from smooth geodesic distances on the paraboloid',y=.98,fontsize=16)
    fig.text(.5,.90,f'n = {n}; uniform base-disk sampling; no graph-distance approximation\n'
             'All clouds divided by input RMS geodesic distance; equal spatial units and common limits',ha='center',fontsize=11)
    cax=fig.add_axes([.93,.25,.012,.46]);fig.colorbar(pts,cax=cax,label='Squared base radius / r²')
    fig.text(.5,.035,'Rigid alignment to the original coordinate frame. Color identifies radial position.',ha='center',fontsize=10)
    fig.savefig(out/'geodesic_3d_snapshots.png',dpi=180);fig.savefig(out/'geodesic_3d_snapshots.pdf');plt.close(fig)

    styles=[('classical_3d',GOLD,'s','--','Classical MDS, 3D'),
            ('stress_3d',BLUE,'o','-','Raw-stress MDS, 3D'),
            ('stress_2d',INK,'^',':','Raw-stress MDS, 2D')]
    fig,axes=plt.subplots(1,3,figsize=(15,4.7),layout='constrained')
    for method,color,marker,ls,label in styles:
        subset=[a for a in rows if a['sampling']=='disk' and a['method']==method]
        for ax,key in zip(axes,['relative_distance_rmse','second_over_first','third_over_second']):
            if key=='third_over_second' and method=='stress_2d': continue
            ax.plot([a['radius'] for a in subset],[a[key] for a in subset],color=color,marker=marker,ls=ls,label=label,ms=4)
    for ax in axes:
        ax.set_xscale('log');ax.set_xlabel('Disk radius r');ax.grid(alpha=.18)
    axes[0].set(yscale='log',ylabel='Relative distance RMSE',title='Distance approximation')
    axes[0].legend(fontsize=9,loc='lower right')
    axes[1].set(ylabel='Second / first singular value',title='Elongation',ylim=(0,1.05))
    axes[2].set(ylabel='Third / second singular value',title='Planarity within the smaller axes',ylim=(0,1.05))
    fig.suptitle(f'Smooth-geodesic MDS · uniform base disk · n = {n}\n'
                 'Singular values describe centered coordinates: a plane has s₃/s₂ ≈ 0; a line has s₂/s₁ ≈ 0',fontsize=12)
    fig.savefig(out/'geodesic_diagnostics.png',dpi=180);fig.savefig(out/'geodesic_diagnostics.pdf');plt.close(fig)

    fig,axes=plt.subplots(1,4,figsize=(12,3.8),layout='constrained')
    for ax,r in zip(axes,selected):
        item=saved[('disk',r)]; x=pca_display(item['stress2'],item['q'])/item['rms']
        ax.scatter(*x.T,c=item['q'],cmap=CMAP,vmin=0,vmax=1,s=9)
        ax.set(xlim=(-1.5,1.5),ylim=(-1.5,1.5),title=f'r = {r:g}',xlabel='PC1');ax.set_aspect('equal');ax.grid(alpha=.18)
    axes[0].set_ylabel('PC2')
    fig.suptitle(f'Raw-stress MDS, 2D · smooth geodesic Δ · n = {n}\nCoordinates divided by input RMS geodesic distance',fontsize=12)
    fig.savefig(out/'geodesic_2d_snapshots.png',dpi=180);fig.savefig(out/'geodesic_2d_snapshots.pdf');plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for sampling,color,marker,ls,label in [('disk',BLUE,'o','-','Uniform base disk'),('surface_area',GOLD,'s','--','Uniform surface area')]:
        subset=[a for a in rows if a['sampling']==sampling and a['method']=='stress_3d']
        for ax,key in zip(axes,['relative_distance_rmse','second_over_first']):
            ax.plot([a['radius'] for a in subset],[a[key] for a in subset],color=color,marker=marker,ls=ls,label=label)
            ax.set_xscale('log');ax.set_xlabel('Disk radius r');ax.grid(alpha=.18)
    axes[0].set(yscale='log',ylabel='Relative distance RMSE',title='Distance approximation')
    axes[1].set(ylabel='Second / first singular value',title='Elongation',ylim=(0,1.05));axes[1].legend()
    fig.suptitle(f'Sampling-measure comparison · raw-stress 3D MDS · smooth geodesic Δ · n = {n}',fontsize=12)
    fig.savefig(out/'geodesic_sampling.png',dpi=180);fig.savefig(out/'geodesic_sampling.pdf');plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    for method,color in [('classical_3d',GOLD),('stress_3d',BLUE)]:
        for sampling,marker,ls,fill in [('disk','o','-','full'),('surface_area','s','--','none')]:
            subset=[a for a in rows if a['sampling']==sampling and a['method']==method]
            label=('Classical' if method=='classical_3d' else 'Raw stress')+(' · base disk' if sampling=='disk' else ' · surface area')
            for ax,key in zip(axes,['quadratic_coefficient','quadratic_r_squared']):
                ax.plot([a['radius'] for a in subset],[a[key] for a in subset],color=color,marker=marker,ls=ls,fillstyle=fill,label=label,ms=5)
                ax.set_xscale('log');ax.set_xlabel('Disk radius r');ax.grid(alpha=.18)
    axes[0].axhline(1,color=INK,ls=':',label='Original coefficient = 1')
    axes[0].set(ylabel='Fitted quadratic coefficient a',ylim=(0,1.05),title='Paraboloid shape approximation')
    axes[0].legend(fontsize=8,loc='lower right')
    axes[1].set(ylabel='R² of the quadratic shape fit',ylim=(.9,1.005),title='Approximation quality (focused scale)')
    fig.suptitle('Smooth-geodesic MDS in 3D · quadratic fit after rigid alignment\n'
                 'z = a(x²+y²) + bx + cy + d, in original units; coefficients describe fitted shapes, not exact surfaces',fontsize=12)
    fig.savefig(out/'geodesic_shape.png',dpi=180);fig.savefig(out/'geodesic_shape.pdf');plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n',type=int,default=240)
    parser.add_argument('--seed',type=int,default=20260905)
    parser.add_argument('--out',type=Path,default=ROOT.parents[2]/'output/paraboloid-mmds-geodesic')
    args=parser.parse_args();out=args.out;out.mkdir(parents=True,exist_ok=True)
    tests=validate();write_csv(out/'geodesic_validation.csv',tests)
    rng=np.random.default_rng(args.seed);u=rng.uniform(size=args.n);theta=rng.uniform(0,2*np.pi,args.n)
    rows=[];runs=[];saved={};arrays={};bounds=[]
    for sampling in ['disk','surface_area']:
        for r in RADII:
            rho=r*np.sqrt(u) if sampling=='disk' else .5*np.sqrt(np.expm1((2/3)*np.log1p(u*np.expm1(1.5*np.log1p(4*r*r)))))
            truth=np.c_[rho*np.cos(theta),rho*np.sin(theta),rho*rho]
            target=distance_matrix(rho,theta);target_flat=squareform(target)
            rms=np.sqrt(np.mean(target_flat**2));d=target/rms
            lower=pdist(truth)
            angle=np.abs(theta[:,None]-theta[None,:]);angle=np.minimum(angle,2*np.pi-angle)
            meridian=primitives(np.zeros_like(rho),rho)[1]
            upper=np.abs(meridian[:,None]-meridian[None,:])+np.minimum(rho[:,None],rho[None,:])*angle
            assert np.all(target_flat>=lower-1e-9*rms)
            assert np.all(target<=upper+1e-9*rms)
            triangle_violation=max(float(np.max(target-target[:,k,None]-target[k,None,:])) for k in range(args.n))
            assert triangle_violation<1e-8*rms
            bounds.append(dict(sampling=sampling,radius=r,min_geodesic_minus_chord=float(np.min(target_flat-lower)),
                               max_triangle_violation_over_rms=triangle_violation/rms))
            c3=classical(d);c2=c3[:,:2]
            starts2=[('classical',c2),('xy',truth[:,:2]/rms),('xz',truth[:,[0,2]]/rms),('yz',truth[:,[1,2]]/rms)]
            starts2 += [(f'random_{j}',rng.normal(size=(args.n,2))) for j in range(2)]
            s2,diag2,best2=stress_fit(d,starts2)
            starts3=[('classical',c3),('original',truth/rms),
                     ('2d_perturbed_01',np.c_[s2,.01*rng.normal(size=args.n)]),
                     ('2d_perturbed_1',np.c_[s2,.1*rng.normal(size=args.n)])]
            starts3 += [(f'random_{j}',rng.normal(size=(args.n,3))) for j in range(2)]
            s3,diag3,best3=stress_fit(d,starts3)
            for k,diag,best in [(2,diag2,best2),(3,diag3,best3)]:
                for j,run in enumerate(diag):runs.append(dict(sampling=sampling,radius=r,dimension=k,selected=j==best,**run))
            assert metrics(s3,d)['relative_distance_rmse'] <= metrics(s2,d)['relative_distance_rmse']+1e-6
            assert metrics(s3,d)['relative_distance_rmse'] <= metrics(c3,d)['relative_distance_rmse']+1e-6
            for method,x in [('original',truth),('classical_2d',c2*rms),('classical_3d',c3*rms),
                             ('stress_2d',s2*rms),('stress_3d',s3*rms),('height_line',rho[:,None]**2)]:
                stats=metrics(x,target); sv=np.linalg.svd(center(x),compute_uv=False)
                row=dict(sampling=sampling,radius=r,n=args.n,method=method,**stats,
                         third_over_second=float(sv[2]/sv[1]) if len(sv)>2 else 0.,
                         **(quadratic_shape(x,truth) if x.shape[1]==3 else dict(quadratic_coefficient='',quadratic_r_squared='')))
                rows.append(row)
            saved[(sampling,r)]=dict(q=rho*rho/(r*r),rms=rms,truth=truth,classical3=c3*rms,stress3=s3*rms,stress2=s2*rms)
            key=f'{sampling}_r{r:g}'
            arrays.update({key+'_truth':truth,key+'_geodesic':target,key+'_classical3':c3*rms,key+'_stress3':s3*rms,key+'_stress2':s2*rms})
            m=rows[-2]
            print(f'{sampling:12s} r={r:5g} stress3={m["relative_distance_rmse"]:.7g} s2/s1={m["second_over_first"]:.5g} s3/s2={m["third_over_second"]:.5g} best={diag3[best3]["start"]}',flush=True)
    write_csv(out/'metrics.csv',rows);write_csv(out/'optimizer_runs.csv',runs);write_csv(out/'distance_bounds.csv',bounds)
    np.savez_compressed(out/'embeddings.npz',uniform=u,theta=theta,**arrays)
    plots(out,rows,saved,args.n)
    checks=dict(selected_runs_success=all(a['success'] for a in runs if a['selected']),
                failed_runs=sum(not a['success'] for a in runs),
                max_selected_gradient=max(a['gradient_max'] for a in runs if a['selected']),
                max_quad_length_error=max(a['length_quad_relative_error'] for a in tests),
                max_quad_angle_error=max(a['angle_quad_error'] for a in tests),
                max_ode_endpoint_error_over_radius=max(a['ode_endpoint_error_over_radius'] for a in tests))
    manifest=dict(generated_at=datetime.now(ZoneInfo('America/New_York')).strftime('%Y-%m-%d %H:%M:%S %Z'),
                  seed=args.seed,n=args.n,radii=RADII,input_distances='Continuous shortest geodesics on z=x²+y²; no graph',
                  python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,checks=checks,
                  source_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['geodesic.py','geodesic_experiment.py','experiment.py']},
                  files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir()) if p.is_file() and p.name!='manifest.json'})
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
