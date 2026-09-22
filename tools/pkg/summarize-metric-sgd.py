#!/usr/bin/env python3
"""Summarize a frozen benchmark directory without discarding failed fits."""
import csv
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
from statistics import median
import sys

root=Path(sys.argv[1])
rows=list(csv.DictReader((root/"results.csv").open()))
protocol=json.loads((root/"protocol.json").read_text())
cases=list(csv.DictReader((root/"cases.csv").open()))
jobs=json.loads((root/"jobs.json").read_text())
key=lambda item:(item["case"],item["backend"],str(item["budget"]))
if len(rows)!=len(jobs) or Counter(map(key,rows))!=Counter(map(key,jobs)):
    raise SystemExit("The complete job manifest must be accounted for before generating the final report")
by_case=defaultdict(list)
for row in rows: by_case[row["case"]].append(row)

def best(case,backend,threshold):
    eligible=[r for r in by_case[case] if r["backend"]==backend and r["status"]=="ok"
              and math.isfinite(float(r["rmse"])) and float(r["fit_seconds"])<=threshold]
    return min(eligible,key=lambda r:float(r["rmse"])) if eligible else None

statuses=Counter(r["status"] for r in rows)
lines=["# Optional SGD backend: bounded pilot", "",
"This pilot asks whether the new native SGD optimizer is a useful alternative to SMACOF for the same metric-MDS objective. It does not test geometry recovery or validate a choice of graph distances.", "",
"## Methods and limits", "",
"The 50 paired cases use Euclidean distances from fixed three-dimensional Gaussian clouds (64 and 256 points), shortest-path distances on fixed weighted cycle-and-chord graphs (64 and 256 vertices), and the package's karate-club social graph (34 vertices). Each dataset is embedded in two and three dimensions from five saved random starts shared by both backends. Each backend is restarted independently with allowances of 10, 30 and 100 passes/iterations, giving 300 planned fits. Classical MDS is recorded separately as a reference.", "",
"All tuning values are provisional calibration choices. The native backend uses the hybrid schedule, initial rate 0.5, boundary target rate 0.01, switch fraction 0.4, and scoring after every pass. SMACOF uses the unchanged wrapper tolerance. The objective and target scaling agree across backends; one SGD pass is not equated with one SMACOF iteration.", "",
"Each worker has a 30-second wall-clock allowance and a sampled 1-GiB resident-memory allowance; the campaign allowance is 1,200 seconds. Fit times include wrapper work and checkpoint scoring, but exclude process startup, common preloading of SMACOF in both workers, and separately recorded distance preparation. Resident memory includes the entire R worker and is sampled every 50 ms, so it is not an exact peak or an optimizer-only allocation measure.", "",
"The comparisons below are retrospective best-completed-fit envelopes over three independently restarted schedules. They are not anytime trajectories, and the time thresholds exclude the effort spent trying the other two schedules. Their total effort is reported separately. A missing eligible fit is unavailable, not a loss or a silently omitted result.", "",
"## Results", "",
f"Recorded {len(rows)} of 300 planned rows: "+", ".join(f"{count} {status}" for status,count in sorted(statuses.items()))+".", "",
"Each comparison counts one saved dataset/dimension/start pair. A tie means the two normalized residual errors differ by at most 1e-10. The residual error is the square root of summed squared distance error divided by summed squared target distance; lower is better.", "",
"| Common fit-time threshold (seconds) | Both available | SGD lower | SMACOF lower | Ties | At least one unavailable |",
"|---:|---:|---:|---:|---:|---:|"]
comparisons=[]
for threshold in protocol["elapsed_thresholds_seconds"]:
    counts=Counter()
    for case in cases:
        a=best(case["case"],"sgd",threshold); b=best(case["case"],"smacof",threshold)
        record=dict(case=case["case"],threshold=threshold,
                    sgd_available=a is not None,smacof_available=b is not None,
                    sgd_rmse=float(a["rmse"]) if a else "",
                    smacof_rmse=float(b["rmse"]) if b else "",status="unavailable")
        if a is None or b is None: counts["unavailable"]+=1
        else:
            av=float(a["rmse"]); bv=float(b["rmse"])
            outcome="tie" if abs(av-bv)<=1e-10 else ("sgd" if av<bv else "smacof")
            counts[outcome]+=1; counts["both"]+=1
            record.update(sgd_rmse=av,smacof_rmse=bv,status=outcome)
        comparisons.append(record)
    lines.append(f'| {threshold:g} | {counts["both"]} | {counts["sgd"]} | {counts["smacof"]} | {counts["tie"]} | {counts["unavailable"]} |')
lines += ["", "At the predeclared 0.3-second threshold, the following medians use only paired available seeds; availability is shown explicitly. Classical MDS uses its own spectral initialization and is a separate reference, not a matched-start optimizer comparison.", "",
"| Target distances | Points | Dimension | Paired seeds / 5 | SGD median error | SMACOF median error | Classical reference error |",
"|---|---:|---:|---:|---:|---:|---:|"]
groups=defaultdict(list)
for case in cases: groups[(case["family"],case["n"],case["dimension"])].append(case)
labels={"euclidean":"Euclidean cloud","weighted_graph":"Weighted graph","karate_club":"Karate-club graph"}
for key,group in groups.items():
    left=[]; right=[]
    for case in group:
        a=best(case["case"],"sgd",.3); b=best(case["case"],"smacof",.3)
        if a and b: left.append(float(a["rmse"])); right.append(float(b["rmse"]))
    fmt=lambda xs:f"{median(xs):.5g}" if xs else "unavailable"
    lines.append(f'| {labels[key[0]]} | {key[1]} | {key[2]} | {len(left)} | {fmt(left)} | {fmt(right)} | {float(group[0]["classical_rmse"]):.5g} |')
lines += ["", "Total measured effort over all three scheduled fits, rather than the selected fit alone:", "",
"| Backend | Successful fits | Summed recorded fit seconds | Summed worker seconds | Largest sampled worker RSS (MiB) |",
"|---|---:|---:|---:|---:|"]
for backend in ["sgd","smacof"]:
    group=[r for r in rows if r["backend"]==backend]
    ok=[r for r in group if r["status"]=="ok"]
    fit_sum=sum(float(r["fit_seconds"]) for r in group if r["fit_seconds"])
    process_sum=sum(float(r["process_seconds"]) for r in group if r["process_seconds"])
    peak=max((int(r["sampled_peak_rss_bytes"]) for r in group if r["sampled_peak_rss_bytes"]),default=0)
    lines.append(f'| {backend} | {len(ok)} | {fit_sum:.3f} | {process_sum:.3f} | {peak/1024**2:.1f} |')
lines += ["", "A forcibly stopped worker may lack a fit-time measurement; its worker time remains included. Thus recorded fit seconds alone are not a complete effort total when workers are stopped.", "", "## Interpretation", "",
"These small problems and five starts per condition support a local pilot assessment. They do not establish behavior on large retinal datasets, high-dimensional feature data, different graph constructions, or other operating systems. The first implementation keeps quadratic storage and full-pair work. No backend default or tuning value is changed on the basis of this pilot.", "",
"All attempts, warnings, stopping reasons, fitted coordinates, and SGD checkpoint histories remain in this directory. Schedule completion is reported as an iteration limit; it is not counted as numerical convergence. Build, supervisor, worker and fixture hashes are recorded alongside R session information. See results.csv, protocol.json, fixture-sha256.json, matched-time-comparisons.csv and the per-fit logs/RDS files.", ""]
(root/"report.md").write_text("\n".join(lines))
with (root/"matched-time-comparisons.csv").open("w",newline="") as stream:
    writer=csv.DictWriter(stream,fieldnames=list(comparisons[0]));writer.writeheader();writer.writerows(comparisons)
print(root/"report.md")
