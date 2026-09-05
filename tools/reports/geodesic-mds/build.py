#!/usr/bin/env python3
"""Build the combined report from verified existing experiment artifacts."""
import argparse
import csv
from datetime import datetime
from zoneinfo import ZoneInfo
import hashlib
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]
OUT = REPO / 'output/pdf/geodesic-mds'
EXPERIMENTS = REPO / 'tools/experiments'
FIGURES = {
    'paraboloid-mmds-geodesic': ['geodesic_3d_snapshots', 'geodesic_diagnostics',
                               'geodesic_2d_snapshots', 'geodesic_sampling', 'geodesic_shape'],
    'saddle-mmds-geodesic': ['saddle_3d_snapshots', 'saddle_diagnostics', 'saddle_2d_snapshots',
                             'saddle_sampling', 'saddle_shape', 'saddle_vs_paraboloid'],
    'paraboloid-mmds-radius': ['embedding_snapshots', 'radius_diagnostics', 'sampling_comparison'],
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command, cwd=ROOT, capture=False):
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', OMP_NUM_THREADS='1')
    return subprocess.run(command, cwd=cwd, env=env, check=True, text=True,
                          stdout=subprocess.PIPE if capture else None,
                          stderr=subprocess.STDOUT if capture else None)


def read_csv(bundle, name='metrics.csv'):
    with (REPO / 'output' / bundle / name).open() as stream:
        return list(csv.DictReader(stream))


def table(name, columns, rows, caption, label):
    text = '\\par\n\\begin{table}[H]\n\\centering\n\\caption{' + caption + '}\\label{' + label + '}\n'
    text += '\\begin{tabular}{' + 'r' * len(columns) + '}\\toprule\n'
    text += ' & '.join(columns) + ' \\\\\n\\midrule\n'
    text += ''.join(' & '.join(row) + ' \\\\\n' for row in rows)
    text += '\\bottomrule\\end{tabular}\n\\end{table}\n\\par\n'
    (OUT / 'tables' / (name + '.tex')).write_text(text)


def make_tables():
    def sci(value, digits=5):
        mantissa, exponent = f'{float(value):.{digits}e}'.split('e')
        return '$' + mantissa + r'\times10^{' + str(int(exponent)) + '}$'

    for surface, bundle in [('paraboloid', 'paraboloid-mmds-geodesic'), ('saddle', 'saddle-mmds-geodesic')]:
        data = [a for a in read_csv(bundle) if a['sampling'] == 'disk' and a['method'] == 'stress_3d']
        rows = [[f"{float(a['radius']):g}", f"{float(a['relative_distance_rmse']):.6f}",
                 f"{float(a['second_over_first']):.6f}", f"{float(a['third_over_second']):.6f}"] for a in data]
        table(surface + '_sweep', ['$r$', '$E_{\\mathrm{rel}}$', '$s_2/s_1$', '$s_3/s_2$'], rows,
              surface.capitalize() + ': raw-stress MDS in three dimensions, uniform base-disk sampling, '
              '$n=240$. Rows vary the base-disk radius. Relative error uses input geodesic distances '
              'in its denominator; the singular-value ratios are scale-free.', 'tab:' + surface)
    data = read_csv('paraboloid-mmds-geodesic')
    rows = []
    for radius in [.25, 1, 4, 16, 64]:
        a = {row['method']: row for row in data if row['sampling'] == 'disk' and float(row['radius']) == radius}
        rows.append([f'{radius:g}'] + [f"{float(a[method]['quadratic_coefficient']):.4f}" for method in ['classical_3d', 'stress_3d']])
    table('paraboloid_coefficients', ['$r$', 'Classical coefficient $a$', 'Raw-stress coefficient $a$'], rows,
          'Paraboloid quadratic coefficients in original units after rigid alignment. Both methods use '
          'three-dimensional output and base-disk sampling. The original coefficient is one.', 'tab:coefficients')
    data = read_csv('paraboloid-mmds-antipodal', 'fixed_sample.csv')
    rows = [[f"{float(a['radius']):g}", f"{float(a['geodesic_over_r']):.6f}",
             f"{float(a['classical_3d_over_r']):.6f}"] for a in data]
    table('antipodal_finite', ['$r$', 'Geodesic distance$/r$', 'Classical 3D distance$/r$'], rows,
          'Boundary-antipodal distances for 240 fixed base-disk interior points plus the boundary pair '
          '($n=242$). Ratios use the base radius $r$, not RMS input distance. The classical limit '
          'is 2.92190551; the geodesic limit is $\\pi$.', 'tab:antipodal-finite')
    data = read_csv('paraboloid-mmds-antipodal', 'population_quadrature.csv')
    rows = []
    for nodes in [64, 128, 256, 512]:
        a = {row['sampling']: row for row in data if int(row['quadrature_nodes']) == nodes}
        rows.append([str(nodes)] + [f"{float(a[measure]['classical_3d_limit']):.10f}" for measure in ['disk', 'surface_area']])
    table('antipodal_population', ['Quadrature nodes', '$\\kappa_{\\mathrm{disk}}$', '$\\kappa_{\\mathrm{area}}$'], rows,
          'Classical 3D boundary-antipodal coefficient from the population second-order operator. '
          'The two columns use uniform base-disk and limiting uniform surface-area measures. '
          'Refinement differences are numerical convergence evidence, not certified error bounds.', 'tab:antipodal-population')
    data = [a for a in read_csv('paraboloid-mmds-radius') if a['sampling'] == 'disk' and a['method'] == 'stress_2d']
    rows = [[f"{float(a['radius']):g}", f"{float(a['relative_distance_rmse']):.6f}",
             f"{float(a['leading_variance_fraction']):.6f}"] for a in data]
    table('ambient_sweep', ['$r$', '$E_{\\mathrm{rel}}$', 'Leading variance fraction'], rows,
          'Ambient chord-distance comparison: raw-stress MDS in two dimensions with uniform base-disk '
          'sampling. The last column is $s_1^2/\\sum_a s_a^2$ for centered fitted coordinates. '
          'These are not geodesic results.', 'tab:ambient')
    data = read_csv('mds-audit-diagnostics', 'population_spectrum.csv')
    rows = []
    for a in data:
        if a['quadrature_nodes'] != '256':
            continue
        measure = 'Base disk' if a['sampling'] == 'disk' else 'Surface area'
        rows.append([measure, f"{float(a['leading_kernel_eigenvalue']):.6f}"] +
                    [sci(a[key], 3) for key in ['second_kernel_eigenvalue', 'radial_sector_max',
                                                'even_sector_max', 'remaining_odd_sector_max']])
    table('population_sectors', ['Measure', '$\\lambda_1(\\mathcal L)$', '$\\lambda_2(\\mathcal L)$',
                                  '$m=0$', 'Even $m$', 'Odd $m\\ge3$'], rows,
          'Numerical population mode-selection checks with 256 quadrature nodes. The first two '
          'numeric columns are eigenvalues of the logarithmic-mean operator; the remaining columns '
          'are the largest competing eigenvalues in the indicated angular sectors after the '
          'constant and height projections. Each nonzero angular radial eigenfunction has a sine/cosine pair.',
          'tab:population-sectors')
    data = read_csv('mds-audit-diagnostics', 'optimizer_spread.csv')
    rows = [[a['surface'].capitalize(), 'Base disk' if a['sampling'] == 'disk' else 'Surface area',
             f"{float(a['worst_over_best']):.6f}", sci(a['relative_spread'], 3)]
            for a in data if float(a['radius']) == 64]
    table('optimizer_spread', ['Surface', 'Measure', '$S_{\\max}/S_{\\min}$', 'Relative spread'], rows,
          'Spread of the six original 3D raw-stress starts at $r=64$. Relative spread is '
          '$(S_{\\max}-S_{\\min})/S_{\\min}$. A displayed ratio of 1.000000 for the saddle '
          'is rounding to six decimals; the final column preserves the small observed differences.',
          'tab:optimizer-spread')


class Evidence(HTMLParser):
    def __init__(self):
        super().__init__()
        self.entries = {}
        self.active = None

    def handle_starttag(self, tag, attributes):
        attrs = dict(attributes)
        if 'data-citation-key' in attrs:
            self.active = attrs['data-citation-key']
            assert self.active not in self.entries, 'Duplicate citation evidence'
            self.entries[self.active] = [attrs.get('data-status'), False]
        if self.active and tag == 'a' and 'data-source-link' in attrs:
            assert attrs.get('href', '').startswith('https://')
            self.entries[self.active][1] = True

    def handle_endtag(self, tag):
        if tag == 'tr':
            self.active = None


def verify_report():
    source = (OUT / 'report.tex').read_text()
    keys = set()
    for group in re.findall(r'\\cite\w*(?:\[[^]]*\])*\{([^}]+)\}', source):
        keys.update(group.split(','))
    bib = set(re.findall(r'@\w+\{([^,]+),', (OUT / 'references.bib').read_text()))
    evidence = Evidence()
    evidence.feed((OUT / 'citation_verification.html').read_text())
    assert keys == bib == set(evidence.entries)
    assert all(status == 'verified' and linked for status, linked in evidence.entries.values())
    log = (OUT / 'report.log').read_text()
    assert not re.search(r'undefined|multiply defined|Overfull \\[hv]box|Float too large', log, re.I), 'LaTeX warning requires review'
    plates = re.findall(r'\\plate\{([^}]+)\}', source)
    assert len(plates) == 15 and len(set(plates)) == 15
    assert all((OUT / 'figures' / name).is_file() for name in plates)
    text = run(['pdftotext', 'report.pdf', '-'], cwd=OUT, capture=True).stdout
    assert 'unknown' not in text and '\ufffd' not in text
    assert all(f'Figure {i}:' in text for i in range(1, 16)), 'Missing numbered figure'
    (OUT / 'report.txt').write_text(text)
    print('Report verified: four supported citations, resolved cross-references, 15 figures, and no overfull boxes.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    if args.verify_only:
        verify_report()
        manifest = json.loads((OUT / 'manifest.json').read_text())
        for path, checksum in manifest['inputs'].items():
            assert sha(REPO / path) == checksum, path
        for path, checksum in manifest['outputs'].items():
            assert sha(OUT / path) == checksum, path
        print('Report input and output checksums verified.')
        return
    for name in ['paraboloid-mmds-radius', 'saddle-mmds-radius']:
        run([sys.executable, 'verify.py'], cwd=EXPERIMENTS / name)
    run([sys.executable, 'experiment.py', '--verify'], cwd=EXPERIMENTS / 'mds-audit-diagnostics')
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'figures').mkdir(exist_ok=True)
    (OUT / 'tables').mkdir(exist_ok=True)
    (OUT / 'data').mkdir(exist_ok=True)
    inputs = {}
    for name in ['report.tex', 'references.bib', 'citation_verification.html', 'AUDIT_RESPONSE.md']:
        shutil.copy2(ROOT / name, OUT / name)
    for file in ROOT.iterdir():
        if file.is_file():
            inputs[str(file.relative_to(REPO))] = sha(file)
    for bundle, names in FIGURES.items():
        for name in names:
            path = REPO / 'output' / bundle / (name + '.pdf')
            shutil.copy2(path, OUT / 'figures' / path.name)
            inputs[str(path.relative_to(REPO))] = sha(path)
    for bundle in list(FIGURES) + ['paraboloid-mmds-antipodal', 'mds-audit-diagnostics']:
        for path in (REPO / 'output' / bundle).glob('*.csv'):
            inputs[str(path.relative_to(REPO))] = sha(path)
            target = OUT / 'data' / bundle
            target.mkdir(exist_ok=True)
            shutil.copy2(path, target / path.name)
        path = REPO / 'output' / bundle / 'manifest.json'
        inputs[str(path.relative_to(REPO))] = sha(path)
    from axis_scaling import render
    render(OUT / 'figures/saddle_axis_scaling.pdf')
    make_tables()
    stamp = datetime.now(ZoneInfo('America/New_York')).strftime('%Y-%m-%d %H:%M:%S %Z')
    (OUT / 'build_info.tex').write_text('\\renewcommand{\\reportbuilddatetime}{' + stamp + '}\n')
    try:
        result = run(['latexmk', '-pdf', '-interaction=nonstopmode', '-halt-on-error', 'report.tex'], cwd=OUT, capture=True)
    except subprocess.CalledProcessError as error:
        (OUT / 'build_stdout.log').write_text(error.stdout or '')
        print((error.stdout or '')[-6000:])
        raise
    (OUT / 'build_stdout.log').write_text(result.stdout)
    verify_report()
    deliverables = ['report.pdf', 'report.tex', 'references.bib', 'report.bbl', 'build_info.tex', 'citation_verification.html', 'AUDIT_RESPONSE.md']
    deliverables += [str(path.relative_to(OUT)) for folder in ['figures', 'tables'] for path in sorted((OUT / folder).iterdir())]
    deliverables += [str(path.relative_to(OUT)) for path in sorted((OUT / 'data').rglob('*.csv'))]
    (OUT / 'README.txt').write_text('Combined geodesic MDS report\n\nOpen report.pdf. This source bundle compiles with latexmk -pdf report.tex.\nThe build timestamp records the originating build. Figure plates use their native page sizes.\nCanonical source: tools/reports/geodesic-mds/ in the grip repository.\nNumerical inputs and their checksums are recorded in manifest.json; complete CSV tables are in data/.\n')
    deliverables.append('README.txt')
    manifest = dict(generated_at=stamp, inputs=inputs, outputs={name: sha(OUT / name) for name in deliverables})
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    with zipfile.ZipFile(OUT / 'geodesic-mds-report.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
        for name in deliverables + ['manifest.json']:
            archive.write(OUT / name, arcname='geodesic-mds-report/' + name)
    print('PDF:', OUT / 'report.pdf')
    print('Self-contained LaTeX bundle:', OUT / 'geodesic-mds-report.zip')


if __name__ == '__main__':
    main()
