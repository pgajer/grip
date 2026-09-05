#!/usr/bin/env python3
"""Audit the default reproduction, citation evidence, and tree inertia theorem."""
import csv
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re

import numpy as np

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parents[2] / 'output/saddle-mmds-geodesic'


class Citations(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows = {}
        self.active = None

    def handle_starttag(self, tag, attributes):
        attrs = dict(attributes)
        if tag == 'tr' and 'data-citation-key' in attrs:
            self.active = attrs['data-citation-key']
            assert self.active not in self.rows, 'Duplicate citation evidence'
            self.rows[self.active] = dict(status=attrs['data-status'], linked=False)
        if self.active and tag == 'a' and 'data-source-link' in attrs:
            assert attrs.get('href', '').startswith('https://')
            self.rows[self.active]['linked'] = True

    def handle_endtag(self, tag):
        if tag == 'tr':
            self.active = None


def read_csv(name):
    with (OUT / name).open() as stream:
        return list(csv.DictReader(stream))


def main():
    used = set(re.findall(r'\[([a-z]+\d{4})(?:\]|,)', (ROOT / 'README.md').read_text()))
    bib = set(re.findall(r'@\w+\{([\w-]+),', (ROOT / 'references.bib').read_text()))
    evidence = Citations()
    evidence.feed((ROOT / 'citation_verification.html').read_text())
    assert used == bib == set(evidence.rows)
    assert all(row['status'] == 'verified' and row['linked'] for row in evidence.rows.values())

    manifest = json.loads((OUT / 'manifest.json').read_text())
    assert manifest['n'] == 240 and manifest['seed'] == 20260905
    for name, checksum in manifest['source_hashes'].items():
        assert hashlib.sha256((ROOT.parent / name).read_bytes()).hexdigest() == checksum, name
    for name, checksum in manifest['files'].items():
        assert hashlib.sha256((OUT / name).read_bytes()).hexdigest() == checksum, name
    if manifest['paraboloid_comparison_sha256']:
        prior = OUT.parent / 'paraboloid-mmds-geodesic'
        assert hashlib.sha256((prior / 'metrics.csv').read_bytes()).hexdigest() == manifest['paraboloid_comparison_sha256']
        paired = json.loads((prior / 'manifest.json').read_text())
        assert paired['n'] == manifest['n'] and paired['seed'] == manifest['seed']

    rows = read_csv('metrics.csv')
    runs = read_csv('optimizer_runs.csv')
    checks = read_csv('distance_checks.csv')
    validation = read_csv('geodesic_validation.csv')
    assert len(rows) == 120 and len(runs) == 240
    assert len(checks) == 20 and len(validation) == 60
    assert len({(row['sampling'], row['radius'], row['method']) for row in rows}) == 120
    assert sum(row['selected'] == 'True' for row in runs) == 40
    assert all(row['success'] == 'True' for row in runs)
    assert max(float(row['gradient_max']) for row in runs if row['selected'] == 'True') < 1e-7
    assert max(float(row['endpoint_error_over_radius']) for row in validation) < 2e-5
    assert max(float(row['relative_distance_error_bvp']) for row in validation) < 1e-7
    assert max(float(row['max_triangle_violation_over_rms']) for row in checks) < 1e-7
    assert all(float(row['tree_relative_rmse']) < .002 for row in checks if float(row['radius']) == 64)

    arrays = np.load(OUT / 'embeddings.npz')
    n = manifest['n']
    j = np.eye(n) - np.ones((n, n)) / n
    for sampling in ['disk', 'surface_area']:
        radial = arrays['uniform'] ** (.5 if sampling == 'disk' else 1 / 3)
        u = radial[:, None] * np.c_[np.cos(arrays['theta']), np.sin(arrays['theta'])]
        h = u[:, 0] ** 2 - u[:, 1] ** 2
        arm = np.where(h >= 0, np.where(u[:, 0] >= 0, 0, 1), np.where(u[:, 1] >= 0, 2, 3))
        t = np.abs(h)
        f = np.eye(4)[arm] * t[:, None]
        tree = np.where(arm[:, None] == arm[None, :], np.abs(t[:, None] - t[None, :]), t[:, None] + t[None, :])
        b = -.5 * j @ (tree ** 2) @ j
        assert np.linalg.matrix_rank(j @ f) == 4
        assert np.allclose(b, j @ f @ (2 * np.eye(4) - np.ones((4, 4))) @ f.T @ j, atol=1e-12)
        values = np.linalg.eigvalsh(b)
        assert np.sum(values > 1e-9) == 3 and np.sum(values < -1e-9) == 1
        print(f'{sampling}: limiting classical ratios = {np.sqrt(values[-2]/values[-1]):.6f}, {np.sqrt(values[-3]/values[-2]):.6f}')
        for radius in [.1, .25, .5, 1, 2, 4, 8, 16, 32, 64]:
            key = f'{sampling}_r{radius:g}'
            target = arrays[key + '_geodesic']
            assert np.all(np.isfinite(target)) and np.allclose(target, target.T)
            assert np.all(np.diag(target) == 0)
            stats = {row['method']: row for row in rows if row['sampling'] == sampling and float(row['radius']) == radius}
            assert float(stats['stress_3d']['relative_distance_rmse']) <= min(float(stats[method]['relative_distance_rmse']) for method in ['stress_2d', 'classical_3d']) + 1e-6
        for method in ['stress_3d', 'classical_3d']:
            row = next(row for row in rows if row['sampling'] == sampling and float(row['radius']) == 64 and row['method'] == method)
            assert float(row['second_over_first']) > .85 and float(row['third_over_second']) > .85

    print('Verified: two supported citations, source/output checksums, 20 distance matrices, 60 independent geodesic checks, 240 optimizer runs, and tree inertia.')


if __name__ == '__main__':
    main()
