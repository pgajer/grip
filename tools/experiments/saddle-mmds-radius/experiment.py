#!/usr/bin/env python3
"""Radius sweep for geodesic MDS of z=x²-y². Run make run; see README.md."""
import argparse
import csv
from datetime import datetime
from zoneinfo import ZoneInfo
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import numba
import scipy
from scipy.linalg import orthogonal_procrustes
from scipy.spatial.distance import pdist,squareform

ROOT=Path(__file__).resolve().parent
SHARED=ROOT.parent/'paraboloid-mmds-radius'/'experiment.py'
spec=importlib.util.spec_from_file_location('mds_shared',SHARED)
shared=importlib.util.module_from_spec(spec);spec.loader.exec_module(shared)
center,classical,stress_fit,metrics,write_csv=shared.center,shared.classical,shared.stress_fit,shared.metrics,shared.write_csv
from smooth_geodesic import solve_pairs,pair_coordinates,validate

RADII=[.1,.25,.5,1.,2.,4.,8.,16.,32.,64.]
PALETTE=['#2463A6','#AD7918','#72783D','#B45586']
LABELS=['Upper +x arm','Upper −x arm','Lower +y arm','Lower −y arm']


def tree_geometry(u):
    height=u[:,0]**2-u[:,1]**2
    arm=np.where(height>=0,np.where(u[:,0]>=0,0,1),np.where(u[:,1]>=0,2,3))
    length=np.abs(height)
    d=np.where(arm[:,None]==arm[None,:],np.abs(length[:,None]-length[None,:]),length[:,None]+length[None,:])
    return d,arm


def straight_path_length(a,b,r):
    delta=b-a
    c=r*np.linalg.norm(delta,axis=1)
    p=2*r*r*(a[:,0]*delta[:,0]-a[:,1]*delta[:,1])
    slope=2*r*r*(delta[:,0]**2-delta[:,1]**2)
    def primitive(v):return .5*(v*np.sqrt(c*c+v*v)+c*c*np.arcsinh(v/c))
    result=np.sqrt(c*c+(p+slope/2)**2)
    regular=np.abs(slope)>1e-8*np.maximum(c,np.abs(p))
    result[regular]=((primitive(p+slope)-primitive(p))/np.where(slope==0,1,slope))[regular]
    return result


def quadratic_fit(x,truth):
    x=center(x);rot,_=orthogonal_procrustes(x,center(truth));x=x@rot
    design=np.c_[x[:,0]**2,x[:,0]*x[:,1],x[:,1]**2,x[:,:2],np.ones(len(x))]
    coeff=np.linalg.lstsq(design,x[:,2],rcond=None)[0]
    fitted=design@coeff
    return dict(quadratic_xx=float(coeff[0]),quadratic_xy=float(coeff[1]),quadratic_yy=float(coeff[2]),
                quadratic_r_squared=float(1-np.sum((fitted-x[:,2])**2)/np.sum(x[:,2]**2)))


def plots(out,rows,saved,n):
    # Scientific static chart contract: paired spatial views, ordered radius
    # sweeps, actual branch membership as color. Equal spatial units in 3D.
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    selected=[.25,1.,4.,32.]
    fig=plt.figure(figsize=(13,11));gs=fig.add_gridspec(3,4,left=.085,right=.98,top=.84,bottom=.13,wspace=.04,hspace=.10)
    for j,r in enumerate(selected):
        item=saved[('disk',r)];colors=np.array(PALETTE)[item['arm']]
        for k,method in enumerate(['truth','classical3','stress3']):
            x=center(item[method]);rot,_=orthogonal_procrustes(x,center(item['truth']));x=x@rot/item['rms']
            ax=fig.add_subplot(gs[k,j],projection='3d');ax.scatter(*x.T,c=colors,s=8,alpha=.9,depthshade=False)
            ax.set(xlim=(-2,2),ylim=(-2,2),zlim=(-2,2));ax.set_box_aspect((1,1,1));ax.view_init(elev=20,azim=-52)
            for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.set_ticks([-1,0,1])
            ax.tick_params(labelsize=7,pad=-1)
            ax.set_xlabel('x',labelpad=-3);ax.set_ylabel('y',labelpad=-3);ax.set_zlabel('z',labelpad=-3)
            if k==0:ax.set_title(f'r = {r:g}',fontsize=13)
    for label,y in [('Original saddle',.735),('Classical MDS · 3D',.49),('Raw-stress MDS · 3D',.26)]:
        fig.text(.025,y,label,rotation=90,va='center',fontsize=12)
    handles=[plt.Line2D([],[],color=c,marker='o',linestyle='',label=l) for c,l in zip(PALETTE,LABELS)]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.53,.05),ncol=4,frameon=False,fontsize=10)
    fig.suptitle('Saddle MDS from smooth-surface geodesic distances',y=.98,fontsize=17)
    fig.text(.5,.895,f'z = x² − y²; n = {n}; uniform base-disk sampling\n'
             'Each case divided by RMS input geodesic distance; equal spatial units and common limits',ha='center',fontsize=11)
    fig.savefig(out/'saddle_3d_snapshots.png',dpi=180);fig.savefig(out/'saddle_3d_snapshots.pdf');plt.close(fig)

    fig,axes=plt.subplots(1,3,figsize=(15,4.8),layout='constrained')
    styles=[('classical_3d',PALETTE[1],'s','--','Classical MDS, 3D'),('stress_3d',PALETTE[0],'o','-','Raw-stress MDS, 3D'),
            ('stress_2d','#30343B','^',':','Raw-stress MDS, 2D')]
    for method,color,marker,ls,label in styles:
        subset=[a for a in rows if a['sampling']=='disk' and a['method']==method]
        for ax,key in zip(axes,['relative_distance_rmse','second_over_first','third_over_second']):
            if key=='third_over_second' and method=='stress_2d':continue
            ax.plot([a['radius'] for a in subset],[a[key] for a in subset],color=color,marker=marker,ls=ls,label=label,ms=4)
    for ax in axes:ax.set_xscale('log');ax.set_xlabel('Disk radius r');ax.grid(alpha=.18)
    axes[0].set(yscale='log',ylabel='Relative distance RMSE',title='Distance approximation');axes[0].legend(fontsize=9)
    axes[1].set(ylabel='Second / first singular value',ylim=(0,1.05),title='Elongation')
    axes[2].set(ylabel='Third / second singular value',ylim=(0,1.05),title='Planarity')
    fig.suptitle(f'Saddle radius sweep · smooth geodesic Δ · uniform base disk · n = {n}\n'
                 'A nondegenerate plane has s₃/s₂ → 0; a line has s₂/s₁ → 0',fontsize=12)
    fig.savefig(out/'saddle_diagnostics.png',dpi=180);fig.savefig(out/'saddle_diagnostics.pdf');plt.close(fig)

    fig,axes=plt.subplots(1,4,figsize=(12,3.8),layout='constrained')
    for ax,r in zip(axes,selected):
        item=saved[('disk',r)];x=center(item['stress2']);_,_,vt=np.linalg.svd(x,full_matrices=False);x=x@vt.T/item['rms']
        ax.scatter(*x.T,c=np.array(PALETTE)[item['arm']],s=9)
        ax.set(xlim=(-2,2),ylim=(-2,2),title=f'r = {r:g}',xlabel='PC1');ax.set_aspect('equal');ax.grid(alpha=.18)
    axes[0].set_ylabel('PC2')
    fig.suptitle(f'Raw-stress saddle MDS in 2D · smooth geodesic Δ · n = {n}\nCoordinates divided by RMS input geodesic distance',fontsize=12)
    fig.savefig(out/'saddle_2d_snapshots.png',dpi=180);fig.savefig(out/'saddle_2d_snapshots.pdf');plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(11,4.7),layout='constrained')
    for sampling,color,marker,ls,label in [('disk',PALETTE[0],'o','-','Uniform base disk'),('surface_area',PALETTE[1],'s','--','Uniform surface area')]:
        subset=[a for a in rows if a['sampling']==sampling and a['method']=='stress_3d']
        for ax,key in zip(axes,['relative_distance_rmse','third_over_second']):
            ax.plot([a['radius'] for a in subset],[a[key] for a in subset],color=color,marker=marker,ls=ls,label=label)
            ax.set_xscale('log');ax.set_xlabel('Disk radius r');ax.grid(alpha=.18)
    axes[0].set(yscale='log',ylabel='Relative distance RMSE',title='Distance approximation')
    axes[1].set(ylabel='Third / second singular value',ylim=(0,1.05),title='Planarity');axes[1].legend()
    fig.suptitle(f'Sampling-measure comparison · saddle raw-stress 3D MDS · n = {n}',fontsize=12)
    fig.savefig(out/'saddle_sampling.png',dpi=180);fig.savefig(out/'saddle_sampling.pdf');plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
    for method,color,marker,ls,label in styles[:2]:
        subset=[a for a in rows if a['sampling']=='disk' and a['method']==method]
        norm=[np.sqrt(a['quadratic_xx']**2+a['quadratic_yy']**2+a['quadratic_xy']**2/2) for a in subset]
        axes[0].plot([a['radius'] for a in subset],norm,color=color,marker=marker,ls=ls,label=label)
        axes[1].plot([a['radius'] for a in subset],[a['quadratic_r_squared'] for a in subset],color=color,marker=marker,ls=ls,label=label)
    axes[0].axhline(np.sqrt(2),color='#30343B',ls=':',label='Original norm = √2')
    for ax in axes:ax.set_xscale('log');ax.set_xlabel('Disk radius r');ax.grid(alpha=.18)
    axes[0].set(ylabel='Quadratic coefficient norm (original units)',title='Quadratic approximation');axes[0].legend(fontsize=9)
    axes[1].set(ylabel='R² of the quadratic fit',title='Approximation quality',ylim=(0,1.05))
    fig.suptitle('Saddle shape diagnostics after rigid alignment\n'
                 'z = ax² + bxy + cy² + dx + ey + f; declining native-unit coefficients alone do not establish planarity',fontsize=11)
    fig.savefig(out/'saddle_shape.png',dpi=180);fig.savefig(out/'saddle_shape.pdf');plt.close(fig)

    previous=ROOT.parents[2]/'output/paraboloid-mmds-geodesic/metrics.csv'
    if previous.exists():
        other=list(csv.DictReader(previous.open()))
        fig,axes=plt.subplots(1,2,figsize=(11,4.6),layout='constrained')
        for data,color,marker,ls,label in [(rows,PALETTE[0],'s','-','Saddle'),(other,PALETTE[1],'o','--','Paraboloid')]:
            subset=[a for a in data if a['sampling']=='disk' and a['method']=='stress_3d']
            for ax,key in zip(axes,['second_over_first','third_over_second']):
                ax.plot([float(a['radius']) for a in subset],[float(a[key]) for a in subset],color=color,marker=marker,ls=ls,label=label)
                ax.set(xscale='log',xlabel='Disk radius r',ylim=(0,1.05));ax.grid(alpha=.18)
        axes[0].set(ylabel='Second / first singular value',title='Elongation')
        axes[1].set(ylabel='Third / second singular value',title='Planarity');axes[1].legend()
        fig.suptitle('Saddle versus paraboloid · raw-stress MDS in 3D\nSmooth geodesic distances; uniform base-disk sampling',fontsize=12)
        fig.savefig(out/'saddle_vs_paraboloid.png',dpi=180);fig.savefig(out/'saddle_vs_paraboloid.pdf');plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--n',type=int,default=240)
    parser.add_argument('--seed',type=int,default=20260905);parser.add_argument('--out',type=Path,default=ROOT.parents[2]/'output/saddle-mmds-geodesic')
    args=parser.parse_args();out=args.out;out.mkdir(parents=True,exist_ok=True);cache=out/'distance_cache';cache.mkdir(exist_ok=True)
    solver_hash=hashlib.sha256((ROOT/'smooth_geodesic.py').read_bytes()).hexdigest()
    validation=validate();write_csv(out/'geodesic_validation.csv',validation)
    rng=np.random.default_rng(args.seed);uniform=rng.uniform(size=args.n);theta=rng.uniform(0,2*np.pi,args.n)
    rows=[];runs=[];checks=[];saved={};arrays={}
    for sampling in ['disk','surface_area']:
        previous=None;warm=None
        for r in RADII:
            start=time.time()
            rho=r*np.sqrt(uniform) if sampling=='disk' else .5*np.sqrt(np.expm1((2/3)*np.log1p(uniform*np.expm1(1.5*np.log1p(4*r*r)))))
            u=(rho/r)[:,None]*np.c_[np.cos(theta),np.sin(theta)]
            a,b=pair_coordinates(u);truth=np.c_[r*u,r*r*(u[:,0]**2-u[:,1]**2)]
            file=cache/f'{sampling}_r{r:g}_n{args.n}.npz';cached=False
            if file.exists():
                stored=np.load(file)
                if str(stored['solver_hash'])==solver_hash and np.array_equal(stored['u'],u):
                    velocity=stored['velocity'];distances=stored['distances'];cached=True;diag=json.loads(str(stored['diagnostics']))
            if not cached:
                velocity,distances,diag=solve_pairs(a,b,r,b-a if warm is None else warm,previous)
                np.savez_compressed(file,u=u,velocity=velocity,distances=distances,solver_hash=solver_hash,diagnostics=json.dumps(diag))
            previous=(a,b,r,velocity);warm=velocity
            target=squareform(distances);rms=np.sqrt(np.mean(distances**2));d=target/rms
            tree,arm=tree_geometry(u);tree_scaled=tree*r*r
            upper=straight_path_length(a,b,r);lower=pdist(truth)
            assert np.all(distances>=lower-1e-7*rms)
            assert np.all(distances<=upper+1e-7*rms)
            assert np.all(target>=tree_scaled-1e-7*rms)
            assert np.all(target<=tree_scaled+4*r+1e-7*rms)
            triangle=max(float(np.max(target-target[:,j,None]-target[j,None,:])) for j in range(args.n))
            assert triangle<1e-7*rms
            checks.append(dict(sampling=sampling,radius=r,tree_relative_rmse=float(np.linalg.norm(pdist(truth)*0+squareform(tree_scaled)-distances)/np.linalg.norm(distances)),
                               max_triangle_violation_over_rms=triangle/rms,distance_seconds=time.time()-start,cache_used=cached,**diag))
            c3=classical(d);c2=c3[:,:2]
            starts2=[('classical',c2),('xy',truth[:,:2]/rms),('xz',truth[:,[0,2]]/rms),('yz',truth[:,[1,2]]/rms)]
            starts2 += [(f'random_{j}',rng.normal(size=(args.n,2))) for j in range(2)]
            s2,diag2,best2=stress_fit(d,starts2)
            starts3=[('classical',c3),('original',truth/rms),('2d_perturbed_01',np.c_[s2,.01*rng.normal(size=args.n)]),
                     ('2d_perturbed_1',np.c_[s2,.1*rng.normal(size=args.n)])]
            starts3 += [(f'random_{j}',rng.normal(size=(args.n,3))) for j in range(2)]
            s3,diag3,best3=stress_fit(d,starts3)
            assert metrics(s3,d)['relative_distance_rmse']<=min(metrics(s2,d)['relative_distance_rmse'],metrics(c3,d)['relative_distance_rmse'])+1e-6
            for dim,diagnostics,best in [(2,diag2,best2),(3,diag3,best3)]:
                for j,run in enumerate(diagnostics):runs.append(dict(sampling=sampling,radius=r,dimension=dim,selected=j==best,**run))
            tree3=classical(tree_scaled/rms)*rms
            for method,x in [('original',truth),('classical_2d',c2*rms),('classical_3d',c3*rms),('stress_2d',s2*rms),('stress_3d',s3*rms),('tree_classical_3d',tree3)]:
                sv=np.linalg.svd(center(x),compute_uv=False)
                shape=quadratic_fit(x,truth) if x.shape[1]==3 else dict(quadratic_xx='',quadratic_xy='',quadratic_yy='',quadratic_r_squared='')
                rows.append(dict(sampling=sampling,radius=r,n=args.n,method=method,**metrics(x,target),third_over_second=float(sv[2]/sv[1]) if len(sv)>2 else 0.,**shape))
            saved[(sampling,r)]=dict(u=u,arm=arm,rms=rms,truth=truth,classical3=c3*rms,stress3=s3*rms,stress2=s2*rms)
            key=f'{sampling}_r{r:g}'
            arrays.update({key+'_truth':truth,key+'_geodesic':target,key+'_tree':tree_scaled,key+'_classical3':c3*rms,key+'_stress3':s3*rms,key+'_stress2':s2*rms})
            m=rows[-2]
            print(f'{sampling:12s} r={r:5g} raw3={m["relative_distance_rmse"]:.6g} s2/s1={m["second_over_first"]:.5g} s3/s2={m["third_over_second"]:.5g} tree_error={checks[-1]["tree_relative_rmse"]:.5g} retries={diag["retry_solves"]} seconds={time.time()-start:.1f}',flush=True)
            write_csv(out/'metrics.csv',rows);write_csv(out/'optimizer_runs.csv',runs);write_csv(out/'distance_checks.csv',checks)
    np.savez_compressed(out/'embeddings.npz',uniform=uniform,theta=theta,**arrays)
    plots(out,rows,saved,args.n)
    provenance={str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'experiment.py',ROOT/'smooth_geodesic.py',SHARED]}
    prior=ROOT.parents[2]/'output/paraboloid-mmds-geodesic/metrics.csv'
    manifest=dict(generated_at=datetime.now(ZoneInfo('America/New_York')).strftime('%Y-%m-%d %H:%M:%S %Z'),seed=args.seed,n=args.n,
                  input_distances='Unique smooth saddle geodesics; adaptive shooting; no graph approximation',
                  python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,numba=numba.__version__,matplotlib=matplotlib.__version__,
                  max_ode_endpoint_error=max(a['endpoint_error_over_radius'] for a in validation),max_bvp_relative_error=max(a['relative_distance_error_bvp'] for a in validation),
                  selected_runs_success=all(a['success'] for a in runs if a['selected']),failed_optimizer_runs=sum(not a['success'] for a in runs),
                  max_selected_gradient=max(a['gradient_max'] for a in runs if a['selected']),source_hashes=provenance,
                  paraboloid_comparison_sha256=hashlib.sha256(prior.read_bytes()).hexdigest() if prior.exists() else None,
                  files={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob('*')) if p.is_file() and p.name!='manifest.json'})
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps({k:v for k,v in manifest.items() if 'error' in k or 'success' in k or 'failed' in k},indent=2))


if __name__=='__main__':main()
