#!/usr/bin/env python3
"""Render the methods note from saved validation results; never run fits."""
import argparse, csv, hashlib, json, re, subprocess, sys
from datetime import datetime
from pathlib import Path
from statistics import median
from zoneinfo import ZoneInfo
here = Path(__file__).resolve().parent
root = here.parents[2]
out = root / 'output/sparse-mds-methods'
run = out / 'validation_2026-10-07_final'
args = argparse.ArgumentParser()
args.add_argument('--pdf', action='store_true')
args = args.parse_args()
rows = list(csv.DictReader((run/'fits.csv').open()))
assert len(rows) == 288 and all(r['status']=='ok' for r in rows)
documentation_revisions = []
for item in csv.DictReader((run/'source_manifest.csv').open()):
    current = Path(item['path'])
    current_md5 = hashlib.md5(current.read_bytes()).hexdigest()
    if current_md5 == item['md5']:
        continue
    # The schedule-help revision changed comments only. Require the exact
    # benchmark-era source and independently compare parsed R expressions.
    archived = out/'documentation_revision_2026-10-07/metric_mds.before_schedule_documentation.R'
    assert current == root/'R/metric_mds.R' and archived.exists(), item['path']
    assert hashlib.md5(archived.read_bytes()).hexdigest() == item['md5'], archived
    expression = 'a<-commandArgs(TRUE); stopifnot(identical(parse(a[1],keep.source=FALSE),parse(a[2],keep.source=FALSE)))'
    subprocess.run(['Rscript','--vanilla','-e',expression,str(current),str(archived)],check=True)
    documentation_revisions.append({'path':str(current),'archived':str(archived),
        'benchmark_md5':item['md5'],'current_md5':current_md5,
        'verification':'identical parsed R expressions; documentation-only change'})
stamp = datetime.now(ZoneInfo('America/New_York')).strftime('%Y-%m-%d %H:%M:%S %Z')
generated = subprocess.check_output(['Rscript','--vanilla','-e',f'cat(readRDS("{run}/validation.rds")$generated_at)'],text=True)
def score(d,i,k,m):
    return median(float(r['rms_relative_percent']) for r in rows if r['case']==f'quadform_d{d}_index{i}' and int(r['output_dimension'])==k and r['method']==m)
def fmt(x): return '$<0.001$' if x < .001 else (f'{x:.3f}' if x < .01 else f'{x:.2f}')
pc = {(int(re.search(r'_d(\d)',r['case'])[1]),int(r['case'][-1]),int(r['output_dimension'])):float(r['rms_relative_percent']) for r in csv.DictReader((out/'pca-baselines.csv').open())}
e=[]; g=[]
for d in (2,3):
    for i in (0,1):
        for k in (2,3):
            e.append(' & '.join(map(str,(d,i,k)))+' & '+' & '.join(fmt(x) for x in [pc[d,i,k]]+[score(d,i,k,m) for m in ('dense_euclidean','uniform','euclidean_q10','euclidean_farthest_q10')])+r' \\')
            g.append(' & '.join(map(str,(d,i,k)))+' & '+' & '.join(fmt(score(d,i,k,m)) for m in ('dense_geodesic_q10','geodesic_q10'))+r' \\')
results = r'''
\clearpage\section{Measured results}
The final experiment completed all 288 attempted fits: zero failed, zero skipped,
and zero excluded. Each method--case--output combination has three preparation
and fitting seeds, with one start per fit. Each layout was evaluated on all
7,140 unordered pairs in its own target metric. These small dense evaluations
are validation calculations, outside the production sparse path.

\subsection{Euclidean approximation}
Table~\ref{tab:euclidean} shows the median RMS relative distance error, $R$, in
percent. The randomized regional method at $q=10$ improved on uniform landmarks
in 23 of 24 paired comparisons (one worse, no ties); farthest-first improved in
24 of 24. The much larger errors of the uniform baseline on the 3D forms embedded
in two dimensions illustrate a limitation of using landmark interactions alone.
No method wins universally from these results: the experiment uses four fixed
point clouds and only three seeds, without tuning pivot counts or epoch budgets.

\begin{table}[ht]\centering\small
\caption{Euclidean targets: median $R$ (percent; lower is better), over three
seeds. $d$ is latent dimension, $i$ is quadratic-form index, and $k$ is output
dimension. PCA is the deterministic unscaled projection. Dense is the all-pair
SGD reference; Uniform uses 16 random landmarks; Random and Farthest use 16
regional pivots and $q=10$. Values below 0.001\% are shown as such. All fits
are included; there are no missing values.}\label{tab:euclidean}
\begin{tabular}{rrr rrrrr}\toprule
$d$ & $i$ & $k$ & PCA & Dense & Uniform & Random & Farthest\\\midrule
'''+'\n'.join(e)+r'''
\bottomrule\end{tabular}\end{table}

Against the matched dense SGD reference, the randomized regional method's median
paired excess error was 0.30 percentage points, ranging from approximately zero
to 2.97. Farthest-first had median excess 0.24 points and maximum 4.14. These are
pooled descriptions of 24 paired runs; Table~\ref{tab:euclidean} and the linked
figures retain the individual cases and output dimensions. No sparse run beat
its dense reference on this evaluation. The reference is a finite-budget fit,
not a demonstrated optimum or a lower bound.

The index-0 and index-1 2D surfaces are represented in three ambient coordinates;
three output dimensions therefore admit exact Euclidean reconstruction. Dense
SGD and PCA reached numerical precision there. Regional fits approached that
limit but were not identically exact. In contrast, the 3D forms occupy four
ambient coordinates, so their three-coordinate Euclidean embeddings remain
approximate. These observations distinguish latent dimension from the number
of Euclidean coordinates required to reproduce all chord distances.

\clearpage\subsection{Graph geometry and neighborhood sensitivity}
Graph-geodesic stress has a different target. In Table~\ref{tab:geodesic}, both
columns use shortest paths on the same $q=10$ graph. The regional method's median
paired excess error over its dense reference was 0.64 percentage points, with
range 0.21--4.01. Its nonzero errors for 2D surfaces embedded in three dimensions
are compatible with the altered target geometry; they are not evidence that
the Euclidean implementation lost an exact reconstruction.

\begin{table}[ht]\centering
\caption{Graph-geodesic targets at $q=10$: median $R$ (percent) over three seeds.
$d$, $i$ and $k$ have the meanings in Table~\ref{tab:euclidean}. Dense fits all
7,140 pairs; Regional uses 16 internally selected pivots and local edges.
Each row uses the same graph target in both columns.}\label{tab:geodesic}
\begin{tabular}{rrr rr}\toprule
$d$ & $i$ & $k$ & Dense & Regional\\\midrule
'''+'\n'.join(g)+r'''
\bottomrule\end{tabular}\end{table}

Changing $q$ from 5 to 10 to 20 gave pooled median paired excess errors of
0.33, 0.30 and 0.13 percentage points for Euclidean targets, and 1.05, 0.64 and
0.30 for graph targets. The graph reference changes with $q$; those latter
numbers compare approximation error within each graph, not the biological or
geometric merits of different graphs. Improvement is not asserted for every seed.
The \href{../../../output/sparse-mds-methods/neighborhood-sensitivity.pdf}
{neighborhood-sensitivity supplement (Figures 5--8)} shows all four quadratic
forms, both output dimensions, and all three seeds at every $q$.
The \href{../../../output/sparse-mds-methods/error-comparison.pdf}
{primary comparison supplement (Figures 1--4)} shows seed-level errors at $q=10$.
Gray points are individual seeds; red markers and segments show medians and
observed ranges. These figures are descriptive. Three seeds and one sampled
cloud per form provide insufficient replication for population uncertainty
claims; ranges are not confidence intervals.

All quadratic-form graphs were connected before repair. All retained their
full three or four PCA coordinates at the 90\% threshold, with complete neighbor
agreement. Consequently these fits do not validate PCA truncation or repair
performance. Separate unit fixtures cover a 50-coordinate lifted cloud with
actual reduction, projected coincidences with positive original distances,
and three separated clusters requiring bridges. They establish implementation
properties, not the adequacy of 90\% retained variance for microbiome profiles.

\clearpage\section{Implementation verification and limits}
The complete \texttt{grip} test suite passed 6,665 checks, with no failures,
warnings or skips. The 112 targeted \texttt{dgraphs} checks also passed without failures, warnings
or skips. These regression tests cover exact
streamed component-MST cost and bridge count against a small dense reference,
tie reproducibility, fallback repair, PCA search coordinates, and original-space
edge lengths. Independent \texttt{grip} checks compare saved distance rows,
region counts and targets with direct calculations; test the all-pivot limit,
2D/3D fits, published learning rates, scaling, RNG preservation and memory refusal.

The first synthetic dispatch stopped before fitting because the new function
was not yet exported. Regenerating package documentation and exports resolved
that setup failure. A later edge-case test exposed rejection of coincident PCA
scores; the search-coordinate validator was corrected to retain those distinct
original observations. The final 288-fit run was repeated after that correction,
and its source and compiled-library fingerprints were unchanged throughout.
Earlier attempts and their logs remain available.

The final R process took 41.12 seconds of serial wall time and had maximum
resident memory 443,318,272 bytes (422.8 MiB). This includes package loading,
preparation, dense validation references, fitting, scoring and result saving;
it excludes native compilation and report rendering. It is not a sparse-only
memory benchmark. Fit-only timings and preparation timings are recorded
separately in the saved table. No large-$n$ runtime or memory-scaling claim is
made. The implicit-component repair avoids quadratic storage but can still
perform quadratic distance work.

Package documentation and API inventories were regenerated and checked.
The package-wide \texttt{make check-fast} command stopped at its repository
hygiene gate because an existing, unrelated serialized dataset contains a
personal filesystem path. That object was preserved. Thus the new tests pass,
but this work does not establish a clean package build or CRAN readiness.

\section{Reproduction and next use}
The results were generated at \textbf{GENERATEDSTAMP}. This report was rebuilt
at \reportbuilddatetime\ from saved artifacts; rendering does not rerun fits.
The authoritative result bundle is
\href{../../../output/sparse-mds-methods/validation_2026-10-07_final/validation.rds}
{validation\_2026-10-07\_final/validation.rds}, with
\href{../../../output/sparse-mds-methods/validation_2026-10-07_final/fits.csv}{all 288 scores},
\href{../../../output/sparse-mds-methods/validation_2026-10-07_final/source_manifest.csv}
{source and library fingerprints}, saved graphs, preparations, coordinates and
session information. The adjacent \href{README.md}{reproduction record}
identifies commands and logs; \href{build_manifest.json}{the build manifest}
fingerprints the report and evidence.

The implementation now supports the requested three variants without constructing
a full distance matrix during sparse preparation or fitting. For the planned
combined V4 analysis, use Euclidean regional preparation as agreed, retain
$q=10$ with 5/20 sensitivity, and measure actual PCA neighbor retention on the
unique abundance profiles. These synthetic tests do not select an adequate
pivot count, epoch budget or PCA threshold for that dataset. No combined V4
embedding was computed in this methods-validation run.
'''.replace('GENERATEDSTAMP',generated)
source = here/'methods.tex'
s = source.read_text()
s = re.sub(r'% BEGIN VALIDATION RESULTS.*?% END VALIDATION RESULTS',lambda _: '% BEGIN VALIDATION RESULTS\n'+results+'\n% END VALIDATION RESULTS',s,flags=re.S)
metadata = '\\renewcommand{\\reportbuilddatetime}{'+stamp+'}\n'
s = re.sub(r'% BEGIN BUILD METADATA.*?% END BUILD METADATA',lambda _: '% BEGIN BUILD METADATA\n% Generated by build_report.py.\n\\begin{filecontents*}[overwrite]{methods_build_info.tex}\n'+metadata+'\\end{filecontents*}\n% END BUILD METADATA',s,flags=re.S)
source.write_text(s)
(here/'methods_build_info.tex').write_text(metadata)
if args.pdf:
    for i in (1,2):
        with (here/f'build_pass_{i}.log').open('w') as stream:
            subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error',source.name],cwd=here,stdout=stream,stderr=subprocess.STDOUT,check=True)
checker = Path('/Users/pgajer/.codex/notes/agent_instructions/reports/scripts/check_citation_verification.py')
cmd=[sys.executable,str(checker),'--tex',str(source),'--bib',str(here/'references.bib'),'--html',str(here/'citation_verification.html')]
if args.pdf: cmd += ['--log',str(here/'methods.log')]
subprocess.run(cmd,check=True)
files=[source,Path(__file__).resolve(),here/'README.md',here/'Makefile',out/'dgraphs_tests_final_corrected.log',here/'references.bib',here/'citation_verification.html',run/'source_manifest.csv',run/'fits.csv',run/'graphs.csv',run/'validation.rds',root/'tools/pkg/summarize-sparse-methods.R',out/'error-comparison.pdf',out/'neighborhood-sensitivity.pdf',out/'paired_comparisons.csv',out/'pca-baselines.csv',out/'grip_full_tests_final.log']
files += [Path(item['archived']) for item in documentation_revisions]
if args.pdf: files += [here/'methods.pdf']
(here/'build_manifest.json').write_text(json.dumps({'build_time':stamp,'result_generation_time':generated,'fits':288,'pdf_exported':args.pdf,'documentation_only_source_revisions':documentation_revisions,'files':[{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files]},indent=2)+'\n')
print(stamp)
