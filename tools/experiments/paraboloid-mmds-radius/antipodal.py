#!/usr/bin/env python3
"""Check antipodal-boundary geodesic and classical-MDS asymptotics.

See ANTIPODAL.md for derivations and the separate fixed-sample raw-stress theorem.
"""
from datetime import datetime
from zoneinfo import ZoneInfo
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh
from scipy.spatial.distance import pdist, squareform

from experiment import classical, write_csv
from geodesic import distance_matrix

ROOT=Path(__file__).resolve().parent


def logmean(q):
    difference=q[:,None]-q[None,:]
    logs=np.log(q[:,None])-np.log(q[None,:])
    return np.divide(difference,logs,out=np.broadcast_to(q[:,None],difference.shape).copy(),where=logs!=0)


def fixed_sample():
    rng=np.random.default_rng(20260905)
    q=np.r_[rng.uniform(size=240),1,1]
    theta=np.r_[rng.uniform(0,2*np.pi,240),0,np.pi]
    alpha=np.abs(theta[:,None]-theta[None,:]);alpha=np.minimum(alpha,2*np.pi-alpha)
    difference=q[:,None]-q[None,:];logs=np.log(q[:,None])-np.log(q[None,:])
    c=.25*difference*logs+logmean(q)*alpha**2
    qc=q-q.mean()
    h=np.eye(len(q))-np.ones((len(q),len(q)))/len(q)-np.outer(qc,qc)/np.dot(qc,qc)
    t=-.5*h@c@h
    eigenvalues,eigenvectors=eigh(t,subset_by_index=[len(q)-2,len(q)-1])
    assert np.all(eigenvalues>0)
    kappa=np.sqrt(np.sum(eigenvalues*(eigenvectors[-2]-eigenvectors[-1])**2))
    rows=[]
    for r in [4.,16.,64.,256.,1024.]:
        target=distance_matrix(r*np.sqrt(q),theta)
        x=classical(target/r**2,3)*r**2
        ratio=np.linalg.norm(x[-2]-x[-1])/r
        # This feasible raw-stress candidate matches the antipodal pair exactly.
        candidate=np.c_[np.zeros(len(q)),r*r*q]
        candidate[-2,0]=target[-2,-1]/2;candidate[-1,0]=-target[-2,-1]/2
        raw_loss=float(np.sum((pdist(candidate)-squareform(target))**2))
        rows.append(dict(radius=r,geodesic_over_r=target[-2,-1]/r,classical_3d_over_r=ratio,
                         classical_limit=kappa,classical_limit_error=abs(ratio-kappa),
                         geodesic_correction_times_r=(np.pi*r-target[-2,-1])*r,
                         feasible_raw_stress=raw_loss))
    assert rows[-1]['classical_limit_error']<1e-6
    assert abs(rows[-2]['geodesic_correction_times_r']-np.pi**3/96)<1e-6
    return rows,dict(q=q,theta=theta,subleading_squared_distances=c,transverse_gram=t)


def population():
    rows=[]
    for measure in ['disk','surface_area']:
        for n in [64,128,256,512]:
            x,w=leggauss(n);q=(x+1)/2;w=w/2
            if measure=='surface_area':w=w*1.5*np.sqrt(q)
            operator=logmean(q)*np.sqrt(w[:,None]*w[None,:])
            values,vectors=eigh(operator)
            lam=values[-1];f=vectors[:,-1]/np.sqrt(w)
            f1=np.dot((1-q)/(-np.log(q)),w*f)/lam
            kappa=abs(2*np.sqrt(2*lam)*f1)
            # Check competing angular and radial sectors after projecting away
            # the constant and leading height coordinates.
            a=np.sqrt(w);a/=np.linalg.norm(a)
            b=np.sqrt(w)*(q-np.dot(w,q)/np.sum(w));b/=np.linalg.norm(b)
            h=np.eye(n)-np.outer(a,a)-np.outer(b,b)
            radial_max=float(eigh(-np.pi**2/6*h@operator@h,subset_by_index=[n-1,n-1],eigvals_only=True)[0])
            even_max=max(0,-values[0]/4)
            odd_max=lam/9
            assert lam>max(radial_max,even_max,odd_max)
            rows.append(dict(sampling=measure,quadrature_nodes=n,classical_3d_limit=kappa,
                             leading_radial_kernel_eigenvalue=lam,
                             next_radial_sector_eigenvalue=radial_max,
                             max_even_angular_eigenvalue=even_max,
                             max_remaining_odd_angular_eigenvalue=odd_max))
    for measure in ['disk','surface_area']:
        selected=[a for a in rows if a['sampling']==measure]
        assert abs(selected[-1]['classical_3d_limit']-selected[-2]['classical_3d_limit'])<1e-7
    return rows


def main():
    out=ROOT.parents[2]/'output/paraboloid-mmds-antipodal';out.mkdir(parents=True,exist_ok=True)
    finite,arrays=fixed_sample();pop=population()
    write_csv(out/'fixed_sample.csv',finite);write_csv(out/'population_quadrature.csv',pop)
    np.savez_compressed(out/'limiting_matrices.npz',**arrays)
    manifest=dict(generated_at=datetime.now(ZoneInfo('America/New_York')).strftime('%Y-%m-%d %H:%M:%S %Z'),
                  seed=20260905,interior_points=240,added_boundary_points=2,
                  sources={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['antipodal.py','geodesic.py','experiment.py']},
                  files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir()) if p.is_file() and p.name!='manifest.json'})
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Fixed-sample classical 3D coefficient:',finite[-1]['classical_limit'])
    for row in pop:
        if row['quadrature_nodes']==512:print(row)


if __name__=='__main__':main()
