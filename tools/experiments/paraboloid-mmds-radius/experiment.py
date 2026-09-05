#!/usr/bin/env python3
"""Reproducible Euclidean-distance MDS experiment on expanding paraboloid disks.

Run from any directory. Outputs default to output/paraboloid-mmds-radius.
Requires Python >=3.9, numpy, scipy, matplotlib. No grip package files are modified.
"""
import argparse
import csv
import hashlib
import json
import platform
from datetime import datetime
from zoneinfo import ZoneInfo
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import scipy
from scipy.linalg import eigh, orthogonal_procrustes
from scipy.optimize import minimize, check_grad
from scipy.spatial.distance import pdist, squareform

BLUE, GOLD, INK = "#2463A6", "#AD7918", "#30343B"
CMAP = LinearSegmentedColormap.from_list("radius_blue", ["#C7DDF0", BLUE, "#12385E"])


def center(x):
    return x - x.mean(axis=0)


def classical(d, k=3):
    """Positive spectral truncation of B = -J D^2 J / 2, as in cmdscale."""
    a = d * d
    b = -0.5 * (a - a.mean(0)[None, :] - a.mean(1)[:, None] + a.mean())
    ev, vec = eigh(b, subset_by_index=[len(d)-k, len(d)-1])
    return vec[:, ::-1] * np.sqrt(np.maximum(ev[::-1], 0))


def objective(flat, d, k):
    """Mean squared distance residual and analytic coordinate gradient."""
    n = len(d)
    x = flat.reshape(n, k)
    distances = squareform(pdist(x))
    residual = distances - d
    pairs = n * (n - 1) / 2
    ratio = residual / np.maximum(distances, 1e-14)
    np.fill_diagonal(ratio, 0)
    grad = 2 * (ratio.sum(1)[:, None] * x - ratio @ x) / pairs
    return np.sum(residual * residual) / (2 * pairs), grad.ravel()


def stress_fit(d, starts):
    candidates, diagnostics = [], []
    for label, start in starts:
        start = center(start)
        # Profile the initial overall scale analytically.
        distances = pdist(start)
        start *= np.dot(distances, squareform(d)) / np.dot(distances, distances)
        fit = minimize(objective, start.ravel(), args=(d, start.shape[1]), jac=True,
                       method="L-BFGS-B", options=dict(ftol=1e-14, gtol=1e-10,
                                                       maxiter=4000, maxls=40))
        x = center(fit.x.reshape(start.shape))
        loss, grad = objective(x.ravel(), d, start.shape[1])
        candidates.append((loss, x))
        diagnostics.append(dict(start=label, loss=loss, success=bool(fit.success),
                                message=str(fit.message), iterations=int(fit.nit),
                                gradient_max=float(np.max(np.abs(grad)))))
    best = min(range(len(candidates)), key=lambda j: candidates[j][0])
    return candidates[best][1], diagnostics, best


def metrics(x, target):
    residual = pdist(x) - squareform(target)
    sv = np.linalg.svd(center(x), compute_uv=False)
    return dict(relative_distance_rmse=float(np.linalg.norm(residual) / np.linalg.norm(squareform(target))),
                absolute_distance_rmse=float(np.sqrt(np.mean(residual**2))),
                leading_variance_fraction=float(sv[0]**2 / np.sum(sv**2)),
                second_over_first=float(sv[1]/sv[0]) if len(sv)>1 else 0.,
                third_over_first=float(sv[2]/sv[0]) if len(sv)>2 else 0.)


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def figures(out, rows, saved, n):
    # Chart contract: scientific standalone PNG/PDF; radius sweep (ordered line)
    # and matched point clouds (scatter); one point per sampled observation.
    # Compare objective, dimension, scale, radial ordering. Two palette roots;
    # method also encoded by line style/marker. Equal spatial units in snapshots.
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.labelcolor": INK, "text.color": INK,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout="constrained")
    for method, color, marker, ls in [("classical_2d", GOLD, "s", "--"),
                                       ("stress_2d", BLUE, "o", "-"),
                                       ("height_line", INK, "^", ":")]:
        subset = [a for a in rows if a["sampling"] == "disk" and a["method"] == method]
        x = [a["radius"] for a in subset]
        label = {"classical_2d":"Classical MDS, 2D", "stress_2d":"Raw-stress MDS, 2D",
                 "height_line":"Height-only line (candidate)"}[method]
        axes[0].plot(x, [a["relative_distance_rmse"] for a in subset], color=color,
                     marker=marker, ls=ls, label=label, ms=5)
    axes[0].set(xscale="log", yscale="log", xlabel="Disk radius r",
                ylabel="Relative distance RMSE", title="Distance approximation")
    axes[0].legend(fontsize=9, loc="lower center")
    for method, color, marker, ls in [("classical_3d", INK, "s", "--"), ("stress_2d", BLUE, "o", "-")]:
        subset = [a for a in rows if a["sampling"] == "disk" and a["method"] == method]
        axes[1].plot([a["radius"] for a in subset], [a["leading_variance_fraction"] for a in subset],
                     color=color, marker=marker, ls=ls,
                     label="Original / exact 3D MDS" if method == "classical_3d" else "Raw-stress MDS, 2D")
    rr = np.geomspace(.1, 64, 400)
    axes[1].plot(rr, np.maximum(3, rr*rr)/(6+rr*rr), color=GOLD, ls=":", label="Original: population formula")
    axes[1].set(xscale="log", xlabel="Disk radius r", ylabel="Variance in first principal component",
                ylim=(.3, 1.03), title="Concentration along one dimension")
    axes[1].legend(fontsize=9, loc="lower right")
    for ax in axes:
        ax.grid(alpha=.18)
    fig.suptitle(f"Paraboloid radius experiment · Euclidean input distances · n = {n}\n"
                 "Uniform sampling in the base disk; the same random quantiles at every radius", fontsize=13)
    fig.savefig(out/"radius_diagnostics.png", dpi=180)
    fig.savefig(out/"radius_diagnostics.pdf")
    plt.close(fig)

    selected = [.25, 1., 4., 32.]
    fig = plt.figure(figsize=(13, 8.6))
    gs = fig.add_gridspec(2, 4, left=.085, right=.92, top=.81, bottom=.14, wspace=.12, hspace=.23)
    for j, r in enumerate(selected):
        sample = saved[("disk", r)]
        q = sample["q"]
        # Rigid alignment only: reconstruction shown in the original coordinate frame.
        x = sample["aligned3"] / sample["distance_rms"]
        ax = fig.add_subplot(gs[0, j], projection="3d")
        points = ax.scatter(*x.T, c=q, cmap=CMAP, vmin=0, vmax=1, s=9, alpha=.9, depthshade=False)
        ax.set(xlim=(-1.5,1.5), ylim=(-1.5,1.5), zlim=(-1.5,1.5), title=f"r = {r:g}")
        ax.set_box_aspect((1,1,1))
        ax.set_xticks([-1,0,1]); ax.set_yticks([-1,0,1]); ax.set_zticks([-1,0,1])
        ax.set_xlabel("x", labelpad=-2); ax.set_ylabel("y", labelpad=-2); ax.set_zlabel("z", labelpad=-2)
        ax.tick_params(labelsize=7, pad=-1)
        ax.view_init(elev=16, azim=-55)
        ax = fig.add_subplot(gs[1, j])
        x = center(sample["stress2"])
        _, _, vt = np.linalg.svd(x, full_matrices=False)
        x = x @ vt.T / sample["distance_rms"]
        if np.corrcoef(x[:,0], q)[0,1] < 0:
            x[:,0] *= -1
        ax.scatter(*x.T, c=q, cmap=CMAP, vmin=0, vmax=1, s=9, alpha=.9)
        ax.set(xlim=(-1.5,1.5), ylim=(-1.5,1.5), xlabel="Embedding PC1")
        ax.set_aspect("equal")
        ax.set_xticks([-1,0,1]); ax.set_yticks([-1,0,1])
        ax.grid(alpha=.16)
        if j == 0:
            ax.set_ylabel("PC2")
    fig.text(.012,.64,"Exact 3D MDS",rotation=90,va="center",fontsize=12)
    fig.text(.012,.31,"Raw-stress 2D MDS",rotation=90,va="center",fontsize=12)
    cax = fig.add_axes([.94,.23,.014,.42])
    fig.colorbar(points, cax=cax, label="Squared base radius / r²")
    fig.suptitle("MDS point clouds across disk radii", y=.97, fontsize=17)
    fig.text(.5,.88,f"Euclidean Δ; n = {n}; each cloud divided by its input RMS pair distance\n"
             "All panels use equal spatial units and identical limits; 3D reconstruction has zero stress to rounding",
             ha="center", fontsize=11)
    fig.text(.5,.035,"Color tracks radial position. At large r, vertical separation dominates angular separation.",ha="center",fontsize=10)
    fig.savefig(out/"embedding_snapshots.png",dpi=180)
    fig.savefig(out/"embedding_snapshots.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1,2,figsize=(11,4.6),layout="constrained")
    for sampling,color,marker,ls in [("disk",BLUE,"o","-"),("surface_area",GOLD,"s","--")]:
        subset=[a for a in rows if a["sampling"]==sampling and a["method"]=="stress_2d"]
        for ax,field in zip(axes,["relative_distance_rmse","leading_variance_fraction"]):
            ax.plot([a["radius"] for a in subset],[a[field] for a in subset],color=color,marker=marker,ls=ls,
                    label="Uniform base disk" if sampling=="disk" else "Uniform surface area")
            ax.set_xscale("log"); ax.set_xlabel("Disk radius r"); ax.grid(alpha=.18)
    axes[0].set(yscale="log",ylabel="Relative distance RMSE",title="Distance approximation")
    axes[1].set(ylabel="Variance in first principal component",ylim=(.45,1.02),title="Concentration along one dimension")
    axes[1].legend(loc="lower right")
    fig.suptitle(f"Sampling-measure comparison · raw-stress 2D MDS · Euclidean Δ · n = {n}",fontsize=13)
    fig.savefig(out/"sampling_comparison.png",dpi=180); fig.savefig(out/"sampling_comparison.pdf")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n",type=int,default=240)
    parser.add_argument("--seed",type=int,default=20260905)
    parser.add_argument("--out",type=Path,default=Path(__file__).resolve().parents[3]/"output/paraboloid-mmds-radius")
    args = parser.parse_args()
    out=args.out; out.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(args.seed)
    uniform=rng.uniform(size=args.n); theta=rng.uniform(0,2*np.pi,args.n)
    rows=[]; fits=[]; saved={}; arrays={}; checks={}
    # Check a nontrivial analytic gradient independently by finite differences.
    x0=rng.normal(size=(7,2)); d0=squareform(pdist(rng.normal(size=(7,3))))
    checks["gradient_l2_error"]=float(check_grad(lambda x: objective(x,d0,2)[0],lambda x: objective(x,d0,2)[1],x0.ravel()))
    assert checks["gradient_l2_error"] < 1e-6
    for sampling in ["disk", "surface_area"]:
        for r in [.1,.25,.5,1.,2.,4.,8.,16.,32.,64.]:
            if sampling == "disk":
                rho = r * np.sqrt(uniform)
            else:
                # Exact inverse surface-area CDF: F(rho) = ((1+4rho²)^1.5-1)/((1+4r²)^1.5-1).
                rho = .5*np.sqrt(np.expm1((2/3)*np.log1p(uniform*np.expm1(1.5*np.log1p(4*r*r)))))
            truth = np.column_stack([rho*np.cos(theta), rho*np.sin(theta), rho*rho])
            target = squareform(pdist(truth)); rms = np.sqrt(np.mean(pdist(truth)**2)); d=target/rms
            classical3 = classical(d)*rms; classical2 = classical3[:,:2]
            rotation,_=orthogonal_procrustes(center(classical3),center(truth))
            aligned3=center(classical3)@rotation
            relative_recovery=np.linalg.norm(aligned3-center(truth))/np.linalg.norm(center(truth))
            assert relative_recovery < 1e-9
            starts=[("classical",classical2/rms), ("xy",truth[:,:2]/rms),
                    ("xz",truth[:,[0,2]]/rms), ("yz",truth[:,[1,2]]/rms)]
            starts += [(f"random_{j}",rng.normal(size=(args.n,2))) for j in range(2)]
            x, diagnostics, best=stress_fit(d,starts); stress2=x*rms
            for diagnostic in diagnostics:
                fits.append(dict(sampling=sampling,radius=r,selected=diagnostic is diagnostics[best],**diagnostic))
            assert objective(x.ravel(),d,2)[0] <= objective((classical2/rms).ravel(),d,2)[0]+1e-10
            z=center(truth)[:,2]
            flat=truth[:,:2]
            # Permit the flat xy candidate its optimal overall scale.
            flat=flat*(np.dot(pdist(flat),pdist(truth))/np.dot(pdist(flat),pdist(flat)))
            for method, embedding in [("classical_3d",classical3),("classical_2d",classical2),
                                      ("stress_2d",stress2),("height_line",z[:,None]),("flat_xy_scaled",flat)]:
                rows.append(dict(sampling=sampling,radius=r,n=args.n,method=method,**metrics(embedding,target),
                                 recovery_error_3d=relative_recovery if method=="classical_3d" else ""))
            assert rows[-5]["relative_distance_rmse"] < 1e-10
            saved[(sampling,r)] = dict(q=rho*rho/(r*r),truth=truth,aligned3=aligned3,stress2=stress2,distance_rms=rms)
            key=f"{sampling}_r{r:g}"
            arrays.update({key+"_truth":truth,key+"_classical3":classical3,key+"_stress2":stress2})
            print(f"{sampling:12s} r={r:5g} 2D relative RMSE={rows[-3]['relative_distance_rmse']:.6g} "
                  f"PC1={rows[-3]['leading_variance_fraction']:.6f} best={diagnostics[best]['start']}",flush=True)
    write_csv(out/"metrics.csv",rows); write_csv(out/"optimizer_runs.csv",fits)
    np.savez_compressed(out/"embeddings.npz",uniform=uniform,theta=theta,**arrays)
    figures(out,rows,saved,args.n)
    checks.update(max_3d_relative_recovery=max(float(a["recovery_error_3d"]) for a in rows if a["method"]=="classical_3d"),
                  selected_runs_success=all(a["success"] for a in fits if a["selected"]),
                  failed_runs=sum(not a["success"] for a in fits),
                  max_selected_gradient=max(a["gradient_max"] for a in fits if a["selected"]))
    manifest=dict(generated_at=datetime.now(ZoneInfo("America/New_York")).strftime("%Y-%m-%d %H:%M:%S %Z"),
                  seed=args.seed,n=args.n,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                  matplotlib=matplotlib.__version__,input_distances="3D Euclidean chords",checks=checks,
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir()) if p.is_file() and p.name!="manifest.json"})
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps(checks,indent=2))


if __name__ == "__main__":
    main()
