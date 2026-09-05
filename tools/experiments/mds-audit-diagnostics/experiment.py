#!/usr/bin/env python3
"""Reproduce targeted accuracy, optimization, and spectral checks on saved MDS fits.

This is a separate diagnostic experiment. It never changes the original
distance matrices, embeddings, or selected minima. See README.md.
"""
import argparse
import csv
from datetime import datetime
from zoneinfo import ZoneInfo
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_bvp, solve_ivp
from scipy.linalg import eigh
from scipy.spatial.distance import pdist, squareform

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]
OUT = REPO / 'output/mds-audit-diagnostics'
SHARED = ROOT.parent / 'paraboloid-mmds-radius/experiment.py'
sys.path.insert(0, str(SHARED.parent))
from experiment import stress_fit, write_csv


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def csv_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def optimization_checks():
    spread, instability, probe = [], [], []
    rng = np.random.default_rng(640905)
    for surface in ['paraboloid', 'saddle']:
        folder = REPO / f'output/{surface}-mmds-geodesic'
        records = csv_rows(folder / 'optimizer_runs.csv')
        arrays = np.load(folder / 'embeddings.npz')
        for sampling in ['disk', 'surface_area']:
            for radius in [.1, .25, .5, 1, 2, 4, 8, 16, 32, 64]:
                runs = [a for a in records if a['sampling'] == sampling and float(a['radius']) == radius and a['dimension'] == '3']
                loss = np.array([float(a['loss']) for a in runs])
                spread.append(dict(surface=surface, sampling=sampling, radius=radius, starts=len(runs),
                                   best_loss=loss.min(), worst_loss=loss.max(), worst_over_best=loss.max()/loss.min(),
                                   relative_spread=(loss.max()-loss.min())/loss.min()))
            key = sampling + '_r64'
            target = arrays[key + '_geodesic']
            rms = np.sqrt(np.mean(squareform(target)**2))
            d = target / rms
            if surface == 'saddle':
                x = arrays[key + '_stress2'] / rms
                distances = squareform(pdist(x))
                mask = ~np.eye(len(x), dtype=bool)
                assert distances[mask].min() > 1e-10
                weights = np.zeros_like(d)
                weights[mask] = (distances[mask] - d[mask]) / distances[mask]
                pairs = len(x)*(len(x)-1)/2
                hessian = 2*(np.diag(weights.sum(axis=1))-weights)/pairs
                values, vectors = eigh(hessian, subset_by_index=[0, 0])
                eigenvalue = values[0]
                direction = vectors[:, 0]
                loss0 = np.mean((pdist(x)-squareform(d))**2)
                epsilon = .001
                loss1 = np.mean((pdist(np.c_[x, epsilon*direction])-squareform(d))**2)
                curvature_epsilon = 1e-5
                curvature_loss = np.mean((pdist(np.c_[x, curvature_epsilon*direction])-squareform(d))**2)
                finite_difference = 2*(curvature_loss-loss0)/curvature_epsilon**2
                assert eigenvalue < 0 and loss1 < loss0
                assert abs(finite_difference-eigenvalue) < 2e-6
                instability.append(dict(sampling=sampling, radius=64, smallest_hessian_eigenvalue=eigenvalue,
                                        epsilon=epsilon, baseline_mean_stress=loss0, perturbed_mean_stress=loss1,
                                        curvature_epsilon=curvature_epsilon, finite_difference_curvature=finite_difference))
            else:
                x = arrays[key + '_stress3'] / rms
                starts = [('selected_refit', x.copy())]
                starts += [(f'perturb_{scale}', x+scale*rng.normal(size=x.shape)) for scale in [.005, .02, .1]]
                starts += [(f'new_random_{i}', rng.normal(size=x.shape)) for i in range(4)]
                _, diagnostics, best = stress_fit(d, starts)
                old_loss = np.mean((pdist(x)-squareform(d))**2)
                for i, a in enumerate(diagnostics):
                    probe.append(dict(sampling=sampling, radius=64, selected=i == best,
                                      original_selected_loss=old_loss, relative_improvement=(old_loss-a['loss'])/old_loss, **a))
    return spread, instability, probe


def compare_geodesic(cache, pair_index, rtol, atol, tolerance, nodes, mesh_count):
    u = cache['u']
    ii, jj = np.triu_indices(len(u), 1)
    a, b = u[ii[pair_index]], u[jj[pair_index]]
    k = 4*64**2
    def ode(t, state):
        x, y, v, w = state
        accel = k*(v*v-w*w)/(1+k*(x*x+y*y))
        return np.array([v, w, -accel*x, accel*y])
    ivp = solve_ivp(ode, [0, 1], np.r_[a, cache['velocity'][pair_index]],
                    method='DOP853', rtol=rtol, atol=atol, dense_output=True)
    assert ivp.success
    mesh = np.unique(np.r_[ivp.t, np.linspace(0, 1, mesh_count)])
    bvp = solve_bvp(ode, lambda ya, yb: np.r_[ya[:2]-a, yb[:2]-b], mesh,
                    ivp.sol(mesh), tol=tolerance, max_nodes=nodes)
    v = bvp.sol(0)[2:]
    length = 64*np.sqrt(v@v+k*(a[0]*v[0]-a[1]*v[1])**2)
    saved = cache['distances'][pair_index]
    return dict(i=int(ii[pair_index]), j=int(jj[pair_index]), pair_index=int(pair_index),
                tolerance=tolerance, success=bool(bvp.success), bvp_nodes=len(bvp.x),
                saved_length=float(saved), comparison_length=float(length),
                absolute_disagreement=float(abs(length-saved)),
                relative_disagreement=float(abs(length-saved)/saved),
                ivp_endpoint_over_radius=float(np.linalg.norm(ivp.y[:2, -1]-b)))


def geodesic_checks():
    rng = np.random.default_rng(20260906)
    records, attempts, refined = [], [], []
    caches = {}
    for sampling in ['disk', 'surface_area']:
        cache = np.load(REPO / f'output/saddle-mmds-geodesic/distance_cache/{sampling}_r64_n240.npz')
        caches[sampling] = cache
        u, length = cache['u'], cache['distances']
        ii, jj = np.triu_indices(len(u), 1)
        height = np.abs(u[:, 0]**2-u[:, 1]**2)
        indices = np.unique(np.r_[np.argsort(length)[-4:], np.argsort(length)[:4],
                                 np.argsort(np.minimum(height[ii], height[jj]))[:4],
                                 rng.choice(len(ii), 12, replace=False)])
        assert len(indices) == 24
        for ix in indices:
            for tolerance in [1e-8, 1e-7, 1e-6]:
                result = dict(sampling=sampling, **compare_geodesic(cache, ix, 3e-13, 3e-14, tolerance, 50000, 100))
                attempts.append(result)
                if result['success']:
                    records.append(result)
                    break
            assert result['success'], result
        print(f'Additional geodesic comparisons: {sampling}, 24 pairs completed.', flush=True)
    for row in sorted(records, key=lambda a: a['relative_disagreement'], reverse=True)[:3]:
        for tolerance in [1e-8, 1e-9, 1e-10]:
            result = dict(sampling=row['sampling'], **compare_geodesic(caches[row['sampling']], row['pair_index'],
                                                                     3e-14, 1e-14, tolerance, 60000, 200))
            assert result['success'], result
            refined.append(result)
    return records, attempts, refined


def spectral_checks():
    rows = []
    for sampling in ['disk', 'surface_area']:
        for n in [256, 512]:
            x, weights = leggauss(n)
            q, weights = (x+1)/2, weights/2
            if sampling == 'surface_area':
                weights *= 1.5*np.sqrt(q)
            diff, logs = q[:, None]-q[None, :], np.log(q[:, None])-np.log(q[None, :])
            logmean = np.divide(diff, logs, out=np.broadcast_to(q[:, None], diff.shape).copy(), where=logs != 0)
            l = logmean*np.sqrt(weights[:, None]*weights[None, :])
            eigenvalues = eigh(l, eigvals_only=True)
            a = np.sqrt(weights); a /= np.linalg.norm(a)
            b = np.sqrt(weights)*(q-np.dot(weights, q)/weights.sum()); b /= np.linalg.norm(b)
            h = np.eye(n)-np.outer(a, a)-np.outer(b, b)
            radial_log = .25*diff*logs*np.sqrt(weights[:, None]*weights[None, :])
            projection_error = np.linalg.norm(h@radial_log@h)
            radial_max = eigh(-np.pi**2*h@l@h/6, eigvals_only=True)[-1]
            even_max, odd_max = max(0, -eigenvalues[0]/4), eigenvalues[-1]/9
            assert projection_error < 1e-12
            assert eigenvalues[-1] > max(eigenvalues[-2], radial_max, even_max, odd_max)
            rows.append(dict(sampling=sampling, quadrature_nodes=n, leading_kernel_eigenvalue=eigenvalues[-1],
                             second_kernel_eigenvalue=eigenvalues[-2], radial_sector_max=radial_max,
                             even_sector_max=even_max, remaining_odd_sector_max=odd_max,
                             projected_radial_log_norm=projection_error))
    return rows


def verify():
    manifest = json.loads((OUT / 'manifest.json').read_text())
    for name, checksum in manifest['inputs'].items():
        assert sha(REPO / name) == checksum, name
    for name, checksum in manifest['outputs'].items():
        assert sha(OUT / name) == checksum, name
    assert len(csv_rows(OUT / 'optimizer_spread.csv')) == 40
    assert len(csv_rows(OUT / 'additional_paraboloid_starts.csv')) == 16
    pairs, refined = csv_rows(OUT / 'additional_saddle_pairs.csv'), csv_rows(OUT / 'refined_saddle_pairs.csv')
    assert len(pairs) == 48 and len(refined) == 9
    assert all(a['success'] == 'True' for a in pairs+refined)
    assert max(float(a['relative_disagreement']) for a in refined) < 1e-7
    assert len(csv_rows(OUT / 'planar_instability.csv')) == 2
    assert len(csv_rows(OUT / 'population_spectrum.csv')) == 4
    print('Verified diagnostic source/data checksums, 48 extra pairs, 9 refinements, 16 starts, planar instability, and population sectors.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    if parser.parse_args().verify:
        verify()
        return
    OUT.mkdir(parents=True, exist_ok=True)
    spread, instability, probe = optimization_checks()
    for name, rows in [('optimizer_spread', spread), ('planar_instability', instability), ('additional_paraboloid_starts', probe),
                       ('population_spectrum', spectral_checks())]:
        write_csv(OUT / (name+'.csv'), rows)
    print('Optimization and spectral diagnostics completed.', flush=True)
    records, attempts, refined = geodesic_checks()
    for name, rows in [('additional_saddle_pairs', records), ('saddle_bvp_attempts', attempts), ('refined_saddle_pairs', refined)]:
        write_csv(OUT / (name+'.csv'), rows)
    inputs = [ROOT / 'experiment.py', SHARED]
    for surface in ['paraboloid', 'saddle']:
        folder = REPO / f'output/{surface}-mmds-geodesic'
        inputs += [folder/'optimizer_runs.csv', folder/'embeddings.npz']
    inputs += [REPO/f'output/saddle-mmds-geodesic/distance_cache/{s}_r64_n240.npz' for s in ['disk', 'surface_area']]
    manifest = dict(generated_at=datetime.now(ZoneInfo('America/New_York')).strftime('%Y-%m-%d %H:%M:%S %Z'),
                    inputs={str(p.relative_to(REPO)): sha(p) for p in inputs},
                    outputs={p.name: sha(p) for p in sorted(OUT.glob('*.csv'))})
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    verify()


if __name__ == '__main__':
    main()
