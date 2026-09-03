#!/usr/bin/env python3
"""Exact polyhedral distances on successively refined saddle meshes.

Sample points are mesh vertices, so there is no endpoint snapping. Delaunay
triangulation is in the parameter square; paths cross triangle interiors.
Independent smooth-surface boundary-value solves provide spot checks.
"""
import argparse
import csv
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from importlib.metadata import version
import numpy as np
from scipy.spatial import Delaunay
from scipy.integrate import solve_bvp, simpson
from pygeodesic import geodesic


def mesh(x, resolution, amplitude=.8):
    grid = np.linspace(-1, 1, resolution)
    u, v = np.meshgrid(grid, grid)
    uv = np.vstack((x[:, :2], np.column_stack((u.ravel(), v.ravel()))))
    faces = Delaunay(uv).simplices.astype(np.int32)
    vertices = np.column_stack((uv, amplitude*(uv[:, 0]**2-uv[:, 1]**2)))
    return vertices, faces


def smooth_distance(a, b, tolerance=1e-8):
    t = np.linspace(0, 1, 25)
    guess = np.vstack((a[:, None]+(b-a)[:, None]*t, np.repeat((b-a)[:, None], len(t), axis=1)))
    def ode(t, y):
        xy, velocity = y[:2], y[2:]
        grad = np.vstack((1.6*xy[0], -1.6*xy[1]))
        hessian_term = 1.6*(velocity[0]**2-velocity[1]**2)
        return np.vstack((velocity, -grad*hessian_term/(1+np.sum(grad*grad,axis=0))))
    def boundary(ya, yb):
        return np.r_[ya[:2]-a, yb[:2]-b]
    sol = solve_bvp(ode, boundary, t, guess, tol=tolerance,max_nodes=10000)
    tt = np.linspace(0,1,1001); y = sol.sol(tt)
    speed = np.sqrt(np.sum(y[2:]**2,axis=0)+(1.6*y[0]*y[2]-1.6*y[1]*y[3])**2)
    return float(simpson(speed,x=tt)), int(sol.status), float(np.max(np.abs(y[:2])))


def run(out, replicate, resolutions, source_count):
    x = np.loadtxt(out/f"cloud-{replicate:02d}.csv",delimiter=",",skiprows=1)
    if source_count == 1000:
        sources = np.arange(len(x),dtype=np.int32)
    else:
        sources = np.loadtxt(out/f"sources-{replicate:02d}.txt",dtype=np.int32)[:source_count]
    targets = np.arange(len(x),dtype=np.int32)
    for resolution in resolutions:
        dest = out/f"reference-r{replicate:02d}-m{resolution}-s{len(sources)}.npz"
        if dest.exists():
            print(f"Existing {dest.name}",flush=True); continue
        started = time.perf_counter(); V,F = mesh(x,resolution)
        solver = geodesic.PyGeodesicAlgorithmExact(V,F)
        values = np.empty((len(sources),len(x)))
        for row,source in enumerate(sources):
            values[row],_ = solver.geodesicDistances(np.array([source],dtype=np.int32),targets)
            if (row+1)%16==0 or row==0:
                print(f"cloud={replicate} mesh={resolution} sources={row+1}/{len(sources)} elapsed={time.perf_counter()-started:.1f}s",flush=True)
        assert np.isfinite(values).all() and (values>=-1e-12).all()
        assert np.max(np.abs(values[np.arange(len(sources)),sources])) < 1e-12
        elapsed = time.perf_counter()-started
        np.savez_compressed(dest,distances=values,sources=sources,mesh_resolution=resolution,
                            mesh_vertices=len(V),mesh_faces=len(F),elapsed=elapsed)
        # R ingestion uses a wide table: first column is the 1-based source id.
        np.savetxt(dest.with_suffix(".csv"),np.column_stack((sources+1,values)),delimiter=",",fmt="%.17g")
        print(f"SAVED {dest.name}: {elapsed:.2f}s",flush=True)


def checks(out):
    x = np.loadtxt(out/"cloud-01.csv",delimiter=",",skiprows=1)
    V,F = mesh(x[:40],25,amplitude=0)
    solver = geodesic.PyGeodesicAlgorithmExact(V,F)
    d,_ = solver.geodesicDistances(np.array([0],dtype=np.int32),np.arange(40,dtype=np.int32))
    plane_error = np.max(np.abs(d-np.linalg.norm(V[:40]-V[0],axis=1)))
    assert plane_error<1e-9,plane_error
    metadata = {"numpy":version("numpy"),"scipy":version("scipy"),"pygeodesic":version("pygeodesic"),
                "plane_max_absolute_error":float(plane_error)}
    (out/"reference-environment.json").write_text(json.dumps(metadata,indent=2)+"\n")
    print(metadata,flush=True)
    for r in range(1,6):
        x = np.loadtxt(out/f"cloud-{r:02d}.csv",delimiter=",",skiprows=1)
        sources=np.loadtxt(out/f"sources-{r:02d}.txt",dtype=np.int32)
        rng=np.random.default_rng(4211000+r)
        # 128 independently selected target partners, using every reference source.
        pairs=[(int(s),int(rng.choice(np.delete(np.arange(1000),s)))) for s in sources]
        rows=[]
        for s,t in pairs:
            d,status,extent=smooth_distance(x[s,:2],x[t,:2])
            if status!=0 or extent>1+1e-8:
                raise RuntimeError(f"BVP control outside patch or unconverged: {r} {s} {t} {status} {extent}")
            rows.append([s+1,t+1,d,status,extent])
        with (out/f"smooth-checks-{r:02d}.csv").open("w") as f:
            writer=csv.writer(f);writer.writerow(["i","j","distance","status","max_abs_parameter"]);writer.writerows(rows)
        print(f"Smooth BVP controls cloud {r}: {len(rows)} converged and inside patch",flush=True)


def pilot_cloud(job):
    out,r=job
    run(out,r,[41,81],128)
    run(out,r,[161],16)
    return r


if __name__ == "__main__":
    ap=argparse.ArgumentParser();ap.add_argument("out",type=Path)
    ap.add_argument("--cloud",type=int,default=1);ap.add_argument("--resolutions",default="41,81,161")
    ap.add_argument("--sources",type=int,default=128);ap.add_argument("--checks",action="store_true")
    ap.add_argument("--pilot",action="store_true")
    args=ap.parse_args()
    if args.pilot:
        started=time.perf_counter()
        with ProcessPoolExecutor(max_workers=3) as pool:
            completed=list(pool.map(pilot_cloud,[(args.out,r) for r in range(1,6)]))
        (args.out/"reference-run.json").write_text(json.dumps({"completed":completed,"wall_seconds":time.perf_counter()-started,"workers":3},indent=2)+"\n")
    elif args.checks: checks(args.out)
    else: run(args.out,args.cloud,list(map(int,args.resolutions.split(","))),args.sources)
