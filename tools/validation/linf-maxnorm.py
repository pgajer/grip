"""Validate against a preserved comb-V3V4 max-normalization experiment.

Usage: experiment/.venv/bin/python tools/validation/linf-maxnorm.py EXPERIMENT OUT
Requires numpy/scipy and an optimized, compiled local grip source tree. OUT must
not exist. Original inputs/results are read-only. No cohort data enter Git.
"""
import hashlib, json, os, subprocess, sys
from pathlib import Path
import numpy as np

experiment = Path(sys.argv[1]).resolve()
out = Path(sys.argv[2]).resolve()
repo = Path(__file__).resolve().parents[2]
original = experiment / 'results/2026-10-25-v1'
sys.path.insert(0, str(experiment / 'scripts'))
from geometry import distances, neighbors, constraints, score

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(8*1024**2), b''): h.update(b)
    return h.hexdigest()

def verify(path, expected):
    assert sha(path) == expected, f'Hash mismatch: {path}'

out.mkdir(parents=True, exist_ok=False)
manifest = json.loads((original/'input_manifest.json').read_text())
for name in ('maxnorm.npy', 'pilot_indices.npy'):
    verify(original/name, manifest['files'][name])
X = np.load(original/'maxnorm.npy', mmap_mode='r')
indices = np.load(original/'pilot_indices.npy')
for scope, data, h in [('pilot',X[indices],200), ('full',X,400)]:
    folder = out/scope; folder.mkdir()
    np.asarray(data).ravel(order='F').tofile(folder/'input.bin')
    (folder/'task.json').write_text(json.dumps(dict(n=len(data),p=data.shape[1],h=h)))
    old = original/scope/'linf'/f'h{h}'
    m = json.loads((old/'manifest.json').read_text())
    for name, expected in m['files'].items(): verify(old/name, expected)

# Preserve source identity independently of the historical implementation hashes.
identity = {str(p.relative_to(repo)):sha(p) for base in ('R','src')
            for p in (repo/base).glob('*') if p.suffix in ('.R','.cpp','.h')}
(out/'source_identity.json').write_text(json.dumps(identity,indent=2))
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
           MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', OMP_THREAD_LIMIT='1')
subprocess.run(['Rscript',str(Path(__file__).with_suffix('.R')),str(repo),str(original),str(out)],
               check=True,env=env)
reports = {}
for scope, data, h in [('pilot',X[indices],200), ('full',X,400)]:
    folder = out/scope; old = original/scope/'linf'/f'h{h}'
    report = json.loads((folder/'report.json').read_text()); n=len(data)
    m = json.loads((old/'manifest.json').read_text())
    assert report['fixed_pairs'] == m['pairs']
    # Exact equality is required for canonical pair IDs and endpoint counts.
    # L-infinity differences on max-normalized inputs are also exactly reproducible.
    for name in ('pairs','targets','count_i','count_j'):
        dtype = np.int32 if name=='pairs' else np.float64
        actual=np.fromfile(folder/f'fixed_{name}.bin',dtype=dtype)
        expected=np.fromfile(old/f'{name}.bin',dtype=dtype)
        np.testing.assert_array_equal(actual,expected)
    report['preserved_constraints_exact'] = True
    nn = np.load(original/scope/'linf/neighbors.npy')
    aa=np.repeat(np.arange(n),10); bb=nn.ravel()
    keys=np.unique(np.minimum(aa,bb).astype(np.int64)*n+np.maximum(aa,bb))
    edges=np.fromfile(folder/'edges.bin',dtype=np.int32).reshape((-1,2),order='F')-1
    np.testing.assert_array_equal(edges[:,0].astype(np.int64)*n+edges[:,1],keys)
    lengths=np.fromfile(folder/'lengths.bin')
    expected_lengths=np.max(np.abs(data[edges[:,0]]-data[edges[:,1]]),axis=1)
    np.testing.assert_array_equal(lengths,expected_lengths)
    report['preserved_neighbor_union_exact'] = True
    if scope=='pilot':
        pp=np.fromfile(folder/'pivots.bin',dtype=np.int32)-1
        nn_check,nd=neighbors(data,'linf',10)
        np.testing.assert_array_equal(nn,nn_check)
        oracle=constraints(distances(data[pp],data,'linf'),pp,nn,nd)
        for name,expected in zip(('pairs','targets','count_i','count_j'),oracle):
            actual=np.fromfile(folder/f'public_{name}.bin',dtype=expected.dtype)
            np.testing.assert_array_equal(actual,expected.ravel(order='F'))
        report['internally_selected_constraints_exact'] = True
        Z=np.fromfile(folder/'coordinates.bin').reshape((n,3),order='F')
        saved=np.loadtxt(original/'pilot/linf/sparse-h200-seed17-e200/coordinates.tsv',skiprows=1)
        np.testing.assert_allclose(Z,saved,rtol=0,atol=1e-12)
        report['fit']['max_coordinate_difference']=float(np.max(abs(Z-saved)))
        pairs=np.load(original/'pilot/linf/validation_pairs.npy')
        targets=np.load(original/'pilot/linf/validation_targets.npy')
        report['fit']['evaluation']=score(Z,pairs,targets)
    reports[scope]=report
    print(scope, 'preserved constraints and neighbor union match exactly',flush=True)
(out/'summary.json').write_text(json.dumps(reports,indent=2))
print('Validation complete:',out/'summary.json')
